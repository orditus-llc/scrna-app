"""The analysis pipeline: runs the eight stages in order and keeps a record of everything.

`AnalysisRun` is the lab notebook for one analysis of one dataset:

- It runs the stages in a fixed order (load -> QC metrics -> filter -> preprocess ->
  reduce dimensions -> cluster -> markers -> annotate). A stage can't run before the
  stages it depends on. Optional parts (e.g. regression, upper QC limits) are settings,
  never skipped stages.
- It keeps the result of every stage (a checkpoint), so a stage can be redone with new
  settings without repeating earlier work. Redoing a stage clears the stages after it,
  because their results belonged to the old version; they are marked "needs re-running".
- It records the settings, package versions, random seeds, input file checksum, and a
  history of every stage run (with time and a short summary) in `uns["orditus"]`.
- It keeps the analyst's label review (`review_cluster`) apart from the computed results. The
  review belongs to one clustering: recomputing the clusters clears it.
- It keeps the dataset's metadata: what is known about the sample (species, tissues, a short
  description), either supplied with the dataset or confirmed by the user (`confirm_metadata`).
  Species and tissues decide which cell types Cell Type Annotation compares. An app can require
  the metadata to be confirmed before anything after loading runs (`require_metadata`).

`run_analysis(path)` runs all eight stages with the reference settings in one call.
All computation happens in the stage modules; this file only orders and records.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Callable

import anndata as ad
import pandas as pd

from scrna_analysis import (annotation, census_markers, clustering, dim_reduction, load_data, marker_db, markers,
                            preprocessing, qc, review)
from scrna_analysis.params import (
    DEFAULT_PARAMS,
    AnalysisParams,
    AnnotationParams,
    FilterParams,
    HVGParams,
    LeidenParams,
    MarkerParams,
    NeighborsParams,
    NormalizeParams,
    PCAParams,
    QCParams,
    RegressOutParams,
    ScaleParams,
    SuggestionParams,
    UMAPParams,
)

# The eight stages, in the only order they can run.
STAGES = (
    "load", "qc_metrics", "filter", "preprocess",
    "reduce_dimensions", "cluster", "markers", "annotate",
)

# Which settings (fields of AnalysisParams) each stage uses.
STAGE_SETTINGS = {
    "load": (),
    "qc_metrics": ("qc", "suggest"),
    "filter": ("filter", "suggest"),
    "preprocess": ("normalize", "hvg", "regress_out", "scale"),
    "reduce_dimensions": ("pca", "neighbors", "umap"),
    "cluster": ("leiden",),
    "markers": ("markers",),
    "annotate": ("annotation",),
}

RECORD_KEY = "orditus"  # adata.uns entry holding the record (JSON text, h5ad-safe)
# Bump if the record's structure changes (2: "review", 2026-10-29; 3: input "description", 2026-11-01;
# 4: input "metadata" replaces "description", with who confirmed it and when, 2026-11-03;
# 5: input "samples" (and "sample_column") when there are several samples, 2026-11-09).
RECORD_FORMAT = 6
MAX_ABOUT_LENGTH = 300  # the metadata's free-text description of the sample
PACKAGES = (
    "scanpy", "anndata", "numpy", "scipy", "pandas", "scikit-learn", "statsmodels",
    "umap-learn", "pynndescent", "leidenalg", "igraph", "h5py", "matplotlib",
)


class PipelineError(ValueError):
    """A stage was requested out of order, or a result doesn't exist yet."""


@dataclass(frozen=True)
class StageResult:
    """The saved result (checkpoint) of one stage."""

    stage: str
    adata: ad.AnnData  # the data after this stage
    summary: dict[str, Any]  # short, plain summary (also written to the record)
    details: Any  # the stage modules' full summaries (for display)
    settings: dict[str, Any]  # the settings this stage used
    ran_at: str  # date and time (UTC)


class AnalysisRun:
    """One analysis of one dataset: runs stages in order, keeps checkpoints, records all.

    Example
    -------
    >>> run = AnalysisRun("my_data.h5ad")            # general defaults (DEFAULT_PARAMS)
    >>> run.load(); run.compute_qc_metrics()
    >>> run.preview_filter()                          # look before filtering
    >>> run.apply_filter(); run.preprocess(); run.reduce_dimensions()
    >>> run.cluster(LeidenParams(resolution=0.5))     # any stage can take new settings
    >>> run.status()                                  # what's done / needs re-running
    """

    def __init__(self, data_path: str | Path | None, params: AnalysisParams = DEFAULT_PARAMS,
                 file_name: str | None = None, marker_db_path: str | Path = marker_db.PANGLAODB_FILE,
                 metadata: dict[str, Any] | None = None, require_metadata: bool = False,
                 sample_files: list[load_data.SampleFile] | None = None,
                 sample_columns: dict[str, str | None] | None = None,
                 census_dir: str | Path = census_markers.CENSUS_DIR,
                 ontology_path: str | Path = census_markers.ONTOLOGY_FILE):
        # None = no file chosen yet (an app starts empty); load() then explains that.
        self.data_path = Path(data_path) if data_path is not None else None
        # Several samples, one file each (then `data_path` is None), or None for one file.
        self.sample_files = list(sample_files) if sample_files else None
        # One file with several samples: the columns naming them ({"sample": ..., "condition": ...}),
        # chosen after loading (`define_samples`); None = the file is one dataset.
        self.sample_columns: dict[str, str | None] | None = dict(sample_columns) if sample_columns else None
        # The file's original name (an uploaded file is stored under a temporary name).
        self.file_name = file_name or (self.data_path.name if self.data_path else None)
        if self.sample_files and not self.file_name:
            self.file_name = f"{len(self.sample_files)} samples"
        self.marker_db_path = Path(marker_db_path)  # the cell-type marker database (PanglaoDB)
        self.census_dir = Path(census_dir)  # the Census marker tables (built ahead of time)
        self.ontology_path = Path(ontology_path)  # the Cell Ontology (which Census types are related)
        # What is known about the sample: {"species": "human", "tissues": [...], "about": "..."}.
        # Supplied with the dataset (e.g. demo data: confirmed by the dataset) or confirmed later by
        # the user (`confirm_metadata`). Input, not a result: species and tissues set Cell Type
        # Annotation's settings, and the metadata is recorded with the run.
        self.metadata: dict[str, Any] = _clean_metadata(metadata) if metadata else {}
        if self.metadata:
            self.metadata.update(confirmed_by="dataset", confirmed_at=None)
        self.require_metadata = require_metadata  # stages after load wait for confirmed metadata
        self.params = _with_metadata(params, self.metadata)  # the settings in use (updated by stages)
        self.history: list[dict[str, Any]] = []  # one entry per stage run, in order
        self._results: dict[str, StageResult] = {}  # checkpoints of completed stages
        self._needs_rerun: set[str] = set()  # stages cleared because an earlier one was redone
        # The analyst's review: cluster -> decision (never mixed into the computed results).
        self.review: dict[str, review.ClusterReview] = {}
        self.review_cleared: str | None = None  # why an earlier review was cleared, if it was

    # --- Where things stand ---------------------------------------------------------------

    @property
    def metadata_confirmed(self) -> bool:
        """Whether the species (and tissues) are known: supplied with the dataset or confirmed."""
        return bool(self.metadata.get("confirmed_by"))

    def status(self) -> dict[str, str]:
        """Each stage: "done", "needs re-running" or "not run"."""
        return {
            stage: "done" if stage in self._results
            else "needs re-running" if stage in self._needs_rerun
            else "not run"
            for stage in STAGES
        }

    def result(self, stage: str | None = None) -> ad.AnnData:
        """The data after `stage` (default: the latest completed stage)."""
        return self.stage_result(stage).adata

    def stage_result(self, stage: str | None = None) -> StageResult:
        """The full checkpoint of `stage` (default: the latest completed stage)."""
        if stage is None:
            completed = [s for s in STAGES if s in self._results]
            if not completed:
                raise PipelineError("Nothing has been run yet. Start with load().")
            stage = completed[-1]
        if stage not in STAGES:
            raise PipelineError(f"Unknown stage '{stage}'. Stages: {', '.join(STAGES)}.")
        if stage not in self._results:
            raise PipelineError(f"Stage '{stage}' has no result ({self.status()[stage]}).")
        return self._results[stage]

    def filtered_counts(self) -> ad.AnnData:
        """Raw counts for all genes that passed filtering (the Filter checkpoint).

        Later stages keep only highly variable genes in `layers["counts"]`; this keeps the
        full-gene counts available (CLAUDE.md, AnnData conventions, option A).
        """
        return self.result("filter")

    def record(self) -> dict[str, Any]:
        """The record ("methods section") of this run, as plain data."""
        load_summary = self._results["load"].summary if "load" in self._results else {}
        return {
            "record_format": RECORD_FORMAT,
            "input": {"path": str(self.data_path) if self.data_path else None, "name": self.file_name,
                      "sha256": load_summary.get("input_sha256"), "metadata": self.metadata,
                      **({"samples": load_summary["samples"]} if load_summary.get("samples") else {}),
                      **({"sample_column": self.sample_columns} if self.sample_columns else {})},
            "params": self.params.to_dict(),
            "random_states": {
                "pca": self.params.pca.random_state,
                "neighbors": self.params.neighbors.random_state,
                "leiden": self.params.leiden.random_state,
                "umap": self.params.umap.random_state,
            },
            "versions": package_versions(),
            "stages": self.status(),
            "history": self.history,
            "review": [asdict(r) for r in self.review.values()],
        }

    # --- The routine every stage goes through ----------------------------------------------

    def _run(
        self,
        stage: str,
        new_settings: dict[str, Any],
        compute: Callable[[ad.AnnData | None, AnalysisParams], tuple[ad.AnnData, dict, tuple]],
    ) -> StageResult:
        """Check the order, compute, clear later stages, record, and save the checkpoint.

        `new_settings` holds any settings passed for this stage (None = keep current ones).
        Nothing changes unless the computation succeeds.
        """
        position = STAGES.index(stage)
        if position > 0 and self.require_metadata and not self.metadata_confirmed:
            raise PipelineError("Confirm the species and tissue on Data Overview first.")
        missing = [s for s in STAGES[:position] if s not in self._results]
        if missing:
            raise PipelineError(
                f"Can't run '{stage}' yet: run '{missing[0]}' first. "
                f"Stages run in this order: {' -> '.join(STAGES)}."
            )

        params = replace(self.params, **{k: v for k, v in new_settings.items() if v is not None})
        previous = self._results[STAGES[position - 1]].adata if position else None
        adata, summary, details = compute(previous, params)  # may raise; nothing changed yet

        # Success: commit the settings, clear later stages, and record the run.
        self.params = params
        cleared = [s for s in STAGES[position + 1:] if s in self._results]
        for later in cleared:
            del self._results[later]
            self._needs_rerun.add(later)
        self._needs_rerun.discard(stage)
        # The review names clusters; once the clusters are recomputed, they may be different ones.
        cleared_review = len(self.review) if position <= STAGES.index("cluster") else 0
        if cleared_review:
            self.review = {}
            self.review_cleared = (f"The clusters were recomputed, so the label review of {cleared_review} "
                                   f"{'cluster was' if cleared_review == 1 else 'clusters was'} cleared.")

        all_settings = params.to_dict()
        settings = {name: all_settings[name] for name in STAGE_SETTINGS[stage]}
        ran_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
        self.history.append({
            "stage": stage,
            "ran_at": ran_at,
            "settings": settings,
            "summary": summary,
            "cleared_later_stages": cleared,
            **({"cleared_review": cleared_review} if cleared_review else {}),
        })

        result = StageResult(stage, adata, summary, details, settings, ran_at)
        self._results[stage] = result
        adata.uns[RECORD_KEY] = json.dumps(self.record())  # adata is this stage's own new copy
        return result

    # --- The eight stages (each calls an already-validated module) -------------------------

    def load(self) -> StageResult:
        """Stage 1: open the data file and check it holds raw counts."""
        def compute(_, params):
            samples, sample_notes = (), ()
            if self.sample_files:
                adata, s = load_data.load_samples(self.sample_files)
                # Each sample with its file's fingerprint (its stored path isn't recorded: uploads are temporary).
                samples = tuple({**{k: v for k, v in x.items() if k != "path"}, "sha256": file_sha256(x["path"])}
                                for x in s.samples)
                sha256 = hashlib.sha256("".join(x["sha256"] for x in samples).encode()).hexdigest()
            elif self.data_path is None:
                raise PipelineError("No data file has been chosen yet.")
            else:
                adata, s = load_data.load_dataset(self.data_path, self.file_name)
                sha256 = file_sha256(self.data_path)
                if self.sample_columns:
                    adata, samples, sample_notes = load_data.samples_from_column(
                        adata, self.sample_columns["sample"], self.sample_columns.get("condition"), self.file_name)
            summary = {
                "n_cells": s.n_cells,
                "n_genes": s.n_genes,
                "n_duplicate_gene_names": s.n_duplicate_gene_names,
                "file_format": s.file_format,
                "notes": list(s.notes) + list(sample_notes),
                "input_sha256": sha256,  # several files: a fingerprint of all of them, in order
            }
            if samples:  # only with several samples, so a single dataset's record is as before
                summary["samples"] = [dict(x) for x in samples]
                summary["n_genes_not_in_all_samples"] = s.n_genes_not_in_all_samples
            return adata, summary, (s,)
        return self._run("load", {}, compute)

    def compute_qc_metrics(self, qc_params: QCParams | None = None,
                           suggest_params: SuggestionParams | None = None) -> StageResult:
        """Stage 2: flag mitochondrial genes and calculate QC metrics.

        With suggestion settings (`params.suggest`, as in the app), the filter thresholds are
        then set to the ones suggested from this data (`qc.suggest_filter`); the suggestion and
        how it was found are recorded in the stage summary ("suggested_filter").

        Without them (the reference analysis), the given thresholds are kept, except: if no
        genes match the mitochondrial prefix, % mitochondrial can't be measured, so the
        % mitochondrial filter is switched off (max_pct_counts_mt=None). Recorded in the
        stage summary, never done silently.
        """
        new_settings: dict[str, Any] = {"qc": qc_params, "suggest": suggest_params}
        suggesting = (suggest_params or self.params.suggest) is not None
        switched_off = False
        if "load" in self._results and not suggesting:
            prefix = (qc_params or self.params.qc).mito_prefix
            no_mito = not qc.genes_with_prefix(self._results["load"].adata.var_names, prefix)
            if no_mito and self.params.filter.max_pct_counts_mt is not None:
                new_settings["filter"] = replace(self.params.filter, max_pct_counts_mt=None)
                switched_off = True

        suggested: list[qc.FilterSuggestion] = []

        def compute(adata, params):
            out, s = qc.calculate_qc_metrics(adata, params.qc)
            summary = {
                "n_mito_genes": s.n_mito_genes,
                "mito_qc_available": s.mito_qc_available,
                "warnings": list(s.warnings),
                "mito_filter_switched_off": switched_off,
            }
            if params.suggest is not None:
                suggestion = qc.suggest_filter(out, params.suggest)
                # No mitochondrial genes: the suggestion has no % mitochondrial limit.
                summary["mito_filter_switched_off"] = not s.mito_qc_available
                suggested.append(suggestion)
                summary["suggested_filter"] = {"filter": asdict(suggestion.filter), "details": suggestion.details}
            return out, summary, (s,)
        result = self._run("qc_metrics", new_settings, compute)
        if suggested:  # the stage succeeded: its suggestion becomes the filter to apply
            self.params = replace(self.params, filter=suggested[0].filter)
        return result

    def suggest_filter(self, suggest_params: SuggestionParams) -> qc.FilterSuggestion:
        """Filter thresholds suggested from the QC metrics with other suggestion settings (e.g.
        a different number of MADs). Changes nothing (not a stage; nothing is recorded)."""
        if "qc_metrics" not in self._results:
            raise PipelineError("Can't suggest thresholds yet: run 'qc_metrics' first.")
        return qc.suggest_filter(self.result("qc_metrics"), suggest_params)

    def preview_filter(self, filter_params: FilterParams | None = None) -> qc.FilterPreview:
        """Look before filtering: which cells/genes these thresholds would keep.

        Changes nothing (not a stage; nothing is recorded).
        """
        if "qc_metrics" not in self._results:
            raise PipelineError("Can't preview filtering yet: run 'qc_metrics' first.")
        return qc.preview_filter(self.result("qc_metrics"), filter_params or self.params.filter)

    def apply_filter(self, filter_params: FilterParams | None = None,
                     suggest_params: SuggestionParams | None = None) -> StageResult:
        """Stage 3: remove low-quality cells and rarely detected genes.

        `suggest_params`: the suggestion settings the thresholds were suggested with, if they
        changed (recorded with the stage; they don't change the thresholds themselves)."""
        def compute(adata, params):
            out, preview = qc.apply_filter(adata, params.filter)
            return out, preview.counts(), (preview,)
        return self._run("filter", {"filter": filter_params, "suggest": suggest_params}, compute)

    def preprocess(
        self,
        normalize: NormalizeParams | None = None,
        hvg: HVGParams | None = None,
        regress_out: RegressOutParams | None = None,
        scale: ScaleParams | None = None,
    ) -> StageResult:
        """Stage 4: normalize, log1p, HVGs, keep raw + subset, optional regression, scale."""
        def compute(adata, params):
            a, norm = preprocessing.normalize_and_log(adata, params.normalize)
            a, hvgs = preprocessing.find_hvgs(a, params.hvg)
            a, subset = preprocessing.keep_raw_and_subset(a)
            a, scaled = preprocessing.regress_and_scale(a, params.regress_out, params.scale)
            summary = {
                "median_counts_before": norm.median_counts_before,
                "target_sum": norm.target_sum,  # counts per cell after normalization (the value used)
                "target_is_median": norm.target_is_median,
                "n_hvg": hvgs.n_hvg,
                "raw_shape": list(subset.raw_shape),
                "working_shape": list(subset.working_shape),
                "regressed": scaled.regressed,
                "regress_keys": list(scaled.regress_keys),
                "n_values_at_clip_limit": scaled.n_values_at_clip_limit,
            }
            return a, summary, (norm, hvgs, subset, scaled)
        new = {"normalize": normalize, "hvg": hvg, "regress_out": regress_out, "scale": scale}
        return self._run("preprocess", new, compute)

    def reduce_dimensions(
        self,
        pca: PCAParams | None = None,
        neighbors: NeighborsParams | None = None,
        umap: UMAPParams | None = None,
    ) -> StageResult:
        """Stage 5: PCA, neighbor graph, UMAP."""
        def compute(adata, params):
            a, pcs = dim_reduction.run_pca(adata, params.pca)
            a, graph = dim_reduction.compute_neighbors(a, params.neighbors)
            a, layout = dim_reduction.run_umap(a, params.umap)
            summary = {
                "n_comps": pcs.n_comps,
                "total_variance_ratio": pcs.total_variance_ratio,
                "elbow_pc": pcs.elbow_pc,  # guidance only (shown on the elbow plot)
                "n_neighbors": graph.n_neighbors,
                "n_connections": graph.n_connections,
            }
            return a, summary, (pcs, graph, layout)
        return self._run("reduce_dimensions", {"pca": pca, "neighbors": neighbors, "umap": umap}, compute)

    def cluster(self, leiden: LeidenParams | None = None) -> StageResult:
        """Stage 6: Leiden clustering."""
        def compute(adata, params):
            out, s = clustering.cluster_leiden(adata, params.leiden)
            summary = {"n_clusters": s.n_clusters, "cluster_sizes": s.cluster_sizes, "notes": list(s.notes)}
            return out, summary, (s,)
        return self._run("cluster", {"leiden": leiden}, compute)

    def find_markers(self, marker_params: MarkerParams | None = None) -> StageResult:
        """Stage 7: marker genes per cluster."""
        def compute(adata, params):
            out, s = markers.rank_marker_genes(adata, params.markers)
            summary = {
                "groups": list(s.groups),
                "n_genes_tested": s.n_genes_tested,
                "top10_markers": {g: list(genes[:10]) for g, genes in s.top_markers.items()},
                "notes": list(s.notes),
            }
            return out, summary, (s,)
        return self._run("markers", {"markers": marker_params}, compute)

    def annotate(self, annotation_params: AnnotationParams | None = None) -> StageResult:
        """Stage 8: predict cell types by comparing each cluster's top marker genes with reference
        cell types. CELLxGENE Census first (its marker table for the species and Census tissue);
        PanglaoDB (the species and its tissues) for the clusters Census leaves Unresolved. A
        confident Census label is never replaced. Without a Census tissue, or with Census switched
        off, PanglaoDB alone. Every reference's answer is kept for every cluster."""
        def compute(adata, params):
            chosen = params.annotation.markers
            census = params.annotation.census
            test = params.annotation.test
            if chosen.species is not None and chosen.species not in marker_db.SPECIES:
                raise annotation.AnnotationError(
                    f"The marker references have cell types for {' and '.join(marker_db.SPECIES)} only, so cell "
                    "types can't be predicted for this species. Review the clusters with their marker genes.")
            use_census = bool(census.use and census.tissue)
            if not use_census and not chosen.tissues:
                raise annotation.AnnotationError("Choose at least one tissue, so that only cell types "
                                                 "from that tissue are compared.")
            genes = (adata.raw if params.markers.use_raw and adata.raw is not None else adata).var_names
            notes: list[str] = []
            results: dict[str, tuple] = {}  # reference -> (adata out, summary, match, source)

            if use_census:
                path = census_markers.table_path(chosen.species, census.tissue, census.release, self.census_dir)
                if not path.is_file():
                    raise annotation.AnnotationError(
                        f"There is no CELLxGENE Census marker table for {chosen.species} {census.tissue} yet, so "
                        "Census can't be used. Choose \"PanglaoDB only\" in the settings, or build the table "
                        "(scripts/build_census_markers.py).")
                try:
                    sets, provenance = census_markers.load_marker_table(path)
                    ontology = _cell_ontology(str(self.ontology_path))
                except census_markers.CensusMarkerError as error:
                    raise annotation.AnnotationError(str(error)) from error
                terms = census_markers.cell_type_terms(provenance)

                def related(a: str, b: str) -> bool:
                    return census_markers.related(ontology, terms.get(a, a), terms.get(b, b))

                match = marker_db.match_to_genes(sets, genes, chosen.min_markers)
                out, s = annotation.predict_cell_types(adata, match.sets, test, related=related)
                source = {"name": "CELLxGENE Census", "version": provenance["census_version"],
                          "file": path.name, "sha256": provenance["table_sha256"], "tissues": [census.tissue],
                          "cell_ontology_version": provenance.get("cell_ontology_version"),
                          "citation": provenance.get("citation"), "role": "first"}
                results["Census"] = (out, s, match, source)

            # PanglaoDB: on its own, or the automatic second check for Census's Unresolved clusters.
            db = None
            if chosen.tissues:
                try:
                    db = marker_db.load_panglaodb(self.marker_db_path)
                except marker_db.MarkerDatabaseError:
                    if not use_census:
                        raise
                    notes.append("PanglaoDB isn't available, so Census's Unresolved clusters weren't checked again.")
            elif use_census:
                notes.append(f"No PanglaoDB tissue matches {census.tissue}, so Census's Unresolved clusters "
                             "weren't checked again.")
            if db is not None:
                match = marker_db.match_to_genes(marker_db.select_marker_sets(db, chosen), genes, chosen.min_markers)
                out, s = annotation.predict_cell_types(adata, match.sets, test)
                source = {"name": db.name, "version": db.version, "file": db.file_name, "sha256": db.sha256,
                          "tissues": list(chosen.tissues), "role": "second" if use_census else "first"}
                results["PanglaoDB"] = (out, s, match, source)

            first = "Census" if use_census else "PanglaoDB"
            out, s, match, source = results[first]
            clusters = s.clusters.copy()
            clusters["label_source"] = [first if label != annotation.UNRESOLVED else None
                                        for label in clusters["predicted_label"]]
            if first == "Census" and "PanglaoDB" in results:
                second = results["PanglaoDB"][1].clusters.set_index("cluster")
                for i, row in clusters.iterrows():
                    other = second.loc[row["cluster"]]
                    if row["predicted_label"] == annotation.UNRESOLVED and other["predicted_label"] != annotation.UNRESOLVED:
                        for column in second.columns:
                            if column in clusters.columns:
                                clusters.at[i, column] = other[column]
                        clusters.at[i, "label_source"] = "PanglaoDB"
            label_of = dict(zip(clusters["cluster"], clusters["predicted_label"]))
            out = out.copy()
            out.obs[annotation.PREDICTION_KEY] = pd.Categorical(out.obs[s.groupby].astype(str).map(label_of))
            label_counts = {str(k): int(v) for k, v in out.obs[annotation.PREDICTION_KEY].value_counts().items()}
            combined = replace(s, clusters=clusters, label_counts=label_counts, notes=tuple(s.notes) + tuple(notes))
            record = json.loads(out.uns[annotation.RECORD_KEY])
            record.update(clusters=records(clusters), references={name: r[3] for name, r in results.items()},
                          by_reference={name: records(r[1].clusters) for name, r in results.items()})
            out.uns[annotation.RECORD_KEY] = json.dumps(record)

            unresolved = clusters["predicted_label"] == annotation.UNRESOLVED
            nothing_matched = bool(unresolved.all())
            # Unresolved clusters are also tested against every PanglaoDB tissue's cell types (a hint only).
            anywhere = (_test_all_tissues(adata, db, chosen, genes, test)
                        if db is not None and unresolved.any() else None)
            tested = set(results["PanglaoDB"][2].sets) if "PanglaoDB" in results else set()
            columns = ["cluster", "n_cells", "predicted_label", "status", "best_candidate", "overlap",
                       "n_markers", "adj_p", "label_source"]
            by_columns = ["cluster", "predicted_label", "status", "best_candidate", "overlap", "n_markers", "adj_p"]
            summary = {
                "database": source,  # the first reference
                "references": [r[3] for r in results.values()],
                "n_cell_types_tested": len(s.tests["cell_type"].unique()),
                "markers_matched_ignoring_case": match.case_matched,
                "label_counts": label_counts,
                "clusters": records(clusters[columns]),
                "by_reference": {name: records(r[1].clusters[by_columns]) for name, r in results.items()},
                "notes": list(combined.notes),
                "nothing_matched": nothing_matched,
                # Only a hint for the user; it never changes a label.
                "tissue_hint": _tissue_hint(anywhere, test) if nothing_matched else [],
                # Unresolved cluster -> its best match among every PanglaoDB tissue's cell types, when that is
                # a clear match to a cell type outside the chosen tissues. Never changes a label.
                "best_anywhere": _best_anywhere(anywhere, clusters[unresolved], tested, test),
            }
            details = {"summary": combined, "match": match,
                       "references": {name: {"summary": r[1], "match": r[2], "source": r[3]} for name, r in results.items()}}
            return out, summary, details
        return self._run("annotate", {"annotation": annotation_params}, compute)

    @property
    def samples(self) -> list[dict[str, Any]]:
        """The loaded samples (name, condition, cells, file, …); empty for a single dataset."""
        return list(self._results["load"].summary.get("samples", [])) if "load" in self._results else []

    def define_samples(self, sample_column: str | None, condition_column: str | None = None) -> StageResult:
        """For a single loaded file: use one of its columns as the samples (and optionally another as
        the conditions), or None to treat the file as one dataset again. The file is loaded again with
        that choice, so any later steps need re-running. Raises with a plain-language message."""
        if self.sample_files or self.data_path is None:
            raise PipelineError("Samples can be chosen from a column only for a single loaded file.")
        previous = self.sample_columns
        self.sample_columns = {"sample": sample_column, "condition": condition_column} if sample_column else None
        try:
            return self.load()
        except Exception:
            self.sample_columns = previous  # nothing changed
            raise

    # --- The dataset's metadata (what is known about the sample; not a stage) ---------------

    def confirm_metadata(self, species: str, tissues: tuple[str, ...] | list[str] = (),
                         about: str | None = None, tissue: str | None = None) -> dict[str, Any]:
        """Record the species and tissue of the sample (and an optional short description), as
        confirmed by the user. `tissue` is the Census tissue (e.g. "blood"); `tissues` the PanglaoDB
        tissues (filled in from the Census tissue when not given). They set Cell Type Annotation's
        settings; if it had already run with others, it needs re-running (nothing else depends on
        them). Recorded in the history. Raises `PipelineError` with a plain-language message."""
        if "load" not in self._results:
            raise PipelineError("Load the data before confirming its species and tissue.")
        if self.metadata.get("confirmed_by") == "dataset":
            raise PipelineError("This dataset's species and tissue come with its metadata and can't be changed.")
        metadata = _clean_metadata({"species": species, "tissue": tissue, "tissues": tissues, "about": about})
        if not metadata["species"]:
            raise PipelineError("Choose the species of your sample.")
        has_markers = metadata["species"] in marker_db.SPECIES
        if has_markers and not (metadata["tissue"] or metadata["tissues"]):
            raise PipelineError("Choose the tissue.")
        listed = census_markers.census_tissues(metadata["species"], self.params.annotation.census.release,
                                               self.census_dir) if metadata["tissue"] else []
        if listed and metadata["tissue"] not in {t["tissue"] for t in listed}:
            raise PipelineError(f"Unknown tissue: {metadata['tissue']}.")
        if has_markers:
            try:
                known = marker_db.tissues(marker_db.load_panglaodb(self.marker_db_path))
            except marker_db.MarkerDatabaseError:
                known = None  # without the database the names can't be checked; annotation explains
            unknown = [t for t in metadata["tissues"] if known is not None and t not in known]
            if unknown:
                raise PipelineError(f"Unknown tissue: {', '.join(unknown)}.")
        if len(metadata["about"] or "") > MAX_ABOUT_LENGTH:
            raise PipelineError(f"Keep the description under {MAX_ABOUT_LENGTH} characters.")

        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        # Compared with the settings in use (tissues in any order): the same choice clears nothing.
        markers = self.params.annotation.markers
        changed = ((metadata["species"], metadata["tissue"], set(metadata["tissues"]))
                   != (markers.species, self.params.annotation.census.tissue, set(markers.tissues)))
        cleared = ["annotate"] if changed and "annotate" in self._results else []
        for stage in cleared:
            del self._results[stage]
            self._needs_rerun.add(stage)
        self.metadata = {**metadata, "confirmed_by": "user", "confirmed_at": now}
        self.params = _with_metadata(self.params, self.metadata)
        self.history.append({
            "stage": "metadata",
            "ran_at": now,
            "settings": {},
            "summary": {k: metadata[k] for k in ("species", "tissue", "tissues", "about")},
            "cleared_later_stages": cleared,
        })
        return self.metadata

    # --- The analyst's label review (a human decision; not a stage) ------------------------

    def review_cluster(self, cluster: str, decision: str, label: str | None = None,
                       note: str | None = None) -> review.ClusterReview | None:
        """Record the analyst's decision for one cluster: "confirm" the predicted cell type,
        "rename" it (with `label`), or "clear" the decision. Recorded in the history. Raises
        `review.ReviewError` or `PipelineError` with a plain-language message."""
        if "annotate" not in self._results:
            raise PipelineError("Run Cell Type Annotation before reviewing the labels.")
        predictions = self._results["annotate"].details["summary"].clusters
        predicted = dict(zip(predictions["cluster"].astype(str), predictions["predicted_label"]))
        cluster = str(cluster)
        if cluster not in predicted:
            raise review.ReviewError(f"There is no cluster '{cluster}'.")
        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        if decision == "clear":
            if cluster not in self.review:
                return None
            del self.review[cluster]
            decided = None
        else:
            decided = review.decide(cluster, decision, predicted[cluster], label, note, now)
            self.review[cluster] = decided
        self.review_cleared = None
        self.history.append({
            "stage": "review",
            "ran_at": now,
            "settings": {},
            "summary": {"cluster": cluster, "decision": decided.decision if decided else "cleared",
                        "label": decided.label if decided else None, "predicted_label": predicted[cluster],
                        **review.review_summary(self.review, list(predicted))},
            "cleared_later_stages": [],
        })
        return decided

    def reviewed_labels(self) -> pd.Categorical | None:
        """Each cell's reviewed label ("Not reviewed" for clusters without a decision), in the
        order of the clustered cells; None before clustering."""
        if "cluster" not in self._results:
            return None
        groupby = self.params.leiden.key_added
        return review.reviewed_labels(self._results["cluster"].adata.obs[groupby], self.review)

    def run_all(self) -> AnalysisRun:
        """Run every stage in order with the current settings."""
        self.load()
        self.compute_qc_metrics()
        self.apply_filter()
        self.preprocess()
        self.reduce_dimensions()
        self.cluster()
        self.find_markers()
        self.annotate()
        return self


def run_analysis(data_path: str | Path, params: AnalysisParams = DEFAULT_PARAMS) -> AnalysisRun:
    """Run the whole analysis (all eight stages) on `data_path` and return the run."""
    return AnalysisRun(data_path, params).run_all()


# --- Small helpers ------------------------------------------------------------------------

TISSUE_HINTS = 3  # tissues suggested when no cluster matched the chosen ones


@lru_cache(maxsize=2)
def _cell_ontology(path: str) -> dict[str, dict]:
    """The Cell Ontology, read once per file (it takes about a second)."""
    return census_markers.read_cell_ontology(path)


def _clean_metadata(metadata: dict[str, Any]) -> dict[str, Any]:
    """Species, the Census tissue (or None), PanglaoDB tissues (a sorted list of unique names) and
    "about" (or None), tidied. Without PanglaoDB tissues, those matching the Census tissue are used."""
    species = str(metadata.get("species") or "").strip().lower() or None
    tissue = str(metadata.get("tissue") or "").strip().lower() or None
    tissues = sorted({str(t).strip() for t in (metadata.get("tissues") or ()) if str(t).strip()})
    if tissue and not tissues:
        tissues = sorted(census_markers.panglaodb_tissues_for(tissue))
    about = " ".join(str(metadata.get("about") or "").split()) or None
    return {"species": species, "tissue": tissue, "tissues": tissues, "about": about}


def _with_metadata(params: AnalysisParams, metadata: dict[str, Any]) -> AnalysisParams:
    """Settings with the sample's species and tissues as Cell Type Annotation's (the metadata
    decides them). Without metadata, the settings are kept as they are."""
    if not metadata.get("species"):
        return params
    markers = replace(params.annotation.markers, species=metadata["species"], tissues=tuple(metadata["tissues"]))
    census = replace(params.annotation.census, tissue=metadata.get("tissue"))
    return replace(params, annotation=replace(params.annotation, markers=markers, census=census))


def _test_all_tissues(adata: ad.AnnData, db: marker_db.MarkerDatabase, chosen, genes, test):
    """Every cluster tested against every tissue's cell types (same species, selection rule and
    test). Returns (summary, cell type -> its tissues), or None if no cell type can be tested."""
    every = replace(chosen, tissues=())
    match = marker_db.match_to_genes(marker_db.select_marker_sets(db, every), genes, chosen.min_markers)
    try:
        _, s = annotation.predict_cell_types(adata, match.sets, test)
    except annotation.AnnotationError:  # no cell type has enough markers in this dataset
        return None
    return s, marker_db.tissues_of_cell_types(db, every)


def _clear_match(row, test) -> bool:
    """The best candidate is significant and has enough matching markers (ambiguity aside)."""
    return row.adj_p < test.max_adj_p and row.overlap >= test.min_overlap


def _tissue_hint(anywhere, test) -> list[dict[str, Any]]:
    """When no cluster matched the chosen tissues: which tissues would fit better? A tissue counts
    the clusters whose best candidate (tested against every tissue) is a clear match listed under
    it. The top tissues, most clusters first; empty if nothing fits anywhere. Only a hint."""
    if anywhere is None:
        return []
    s, tissue_of = anywhere
    counts: dict[str, int] = {}
    for row in s.clusters.itertuples():
        if _clear_match(row, test):
            for tissue in tissue_of.get(row.best_candidate, ()):
                counts[tissue] = counts.get(tissue, 0) + 1
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:TISSUE_HINTS]
    return [{"tissue": tissue, "n_clusters": n} for tissue, n in ranked]


def _best_anywhere(anywhere, unresolved: pd.DataFrame, tested: set[str], test) -> dict[str, dict[str, Any]]:
    """For each Unresolved cluster: its best candidate among every tissue's cell types, when that is
    a clear match and a cell type that wasn't tested (outside the chosen tissues). Only a hint."""
    if anywhere is None:
        return {}
    s, tissue_of = anywhere
    rows = s.clusters.set_index(s.clusters["cluster"].astype(str))
    found = {}
    for cluster in unresolved["cluster"].astype(str):
        row = rows.loc[cluster]
        if _clear_match(row, test) and row.best_candidate not in tested:
            found[cluster] = {"cell_type": row.best_candidate, "tissues": list(tissue_of.get(row.best_candidate, ())),
                              "overlap": int(row.overlap), "n_markers": int(row.n_markers), "adj_p": float(row.adj_p)}
    return found

@lru_cache(maxsize=1)
def package_versions() -> dict[str, str]:
    """Python, platform and the versions of the analysis packages (looked up once)."""
    versions = {"python": sys.version.split()[0], "platform": platform.platform()}
    for package in PACKAGES:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "not installed"
    return versions


def file_sha256(path: Path) -> str:
    """Checksum ("fingerprint") of a file: identifies exactly which input was analyzed."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def records(table: pd.DataFrame) -> list[dict[str, Any]]:
    """A table as a list of plain rows (JSON-safe numbers and text; missing values -> None).
    Numbers keep full precision (pandas' to_json rounds to 10 decimals, so p-values like
    1e-22 would become 0)."""
    rows = table.astype(object).where(table.notna(), None).to_dict(orient="records")
    return [{k: v.item() if hasattr(v, "item") else v for k, v in row.items()} for row in rows]
