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
- It keeps the dataset's metadata: what is known about the sample (species, tissue, a short
  description), either supplied with the dataset or confirmed by the user (`confirm_metadata`).
  Species and tissue decide which cell types Cell Type Annotation compares. An app can require
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

from scrna_analysis import (annotation, broad_types, celltypist_models, census_markers, clustering, conditions,
                            consensus, dim_reduction, integration, load_data, marker_db, markers, preprocessing, qc,
                            review)
from scrna_analysis.params import (
    DEFAULT_PARAMS,
    AnalysisParams,
    AnnotationParams,
    ConditionParams,
    FilterParams,
    HVGParams,
    IntegrationParams,
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
    _to_plain,
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
    "reduce_dimensions": ("pca", "integration", "neighbors", "umap"),
    "cluster": ("leiden",),
    "markers": ("markers",),
    "annotate": ("annotation",),
}

RECORD_KEY = "orditus"  # adata.uns entry holding the record (JSON text, h5ad-safe)
# Bump if the record's structure changes (2: "review", 2026-10-29; 3: input "description", 2026-11-01;
# 4: input "metadata" replaces "description", with who confirmed it and when, 2026-11-03;
# 5: input "samples" (and "sample_column") when there are several samples, 2026-11-09;
# 6: metadata "tissue" and annotation "census", 2026-11-09; 7: "integration" settings, 2026-11-10;
# 8: annotation "celltypist", 2026-11-12; 9: PanglaoDB removed: metadata has no "tissues", annotation
# settings have no PanglaoDB selection or census "use", 2026-11-12; 10: annotation "consensus" settings
# and the final label per cluster (Step B); metadata "stage" (adult / developing), 2026-11-13; 11: "conditions"
# settings and "comparison" (Compare Conditions), 2026-11-16).
RECORD_FORMAT = 11
MAX_ABOUT_LENGTH = 300  # the metadata's free-text description of the sample
PACKAGES = (
    "scanpy", "anndata", "numpy", "scipy", "pandas", "scikit-learn", "statsmodels",
    "umap-learn", "pynndescent", "leidenalg", "igraph", "h5py", "matplotlib", "harmonypy", "celltypist",
    "pydeseq2",
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
                 file_name: str | None = None,
                 metadata: dict[str, Any] | None = None, require_metadata: bool = False,
                 sample_files: list[load_data.SampleFile] | None = None,
                 sample_columns: dict[str, str | None] | None = None,
                 census_dir: str | Path = census_markers.CENSUS_DIR,
                 ontology_path: str | Path = census_markers.ONTOLOGY_FILE,
                 celltypist_dir: str | Path = celltypist_models.MODEL_DIR):
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
        self.census_dir = Path(census_dir)  # the Census marker tables (built ahead of time)
        self.ontology_path = Path(ontology_path)  # the Cell Ontology (which Census types are related)
        self.celltypist_dir = Path(celltypist_dir)  # CellTypist's downloaded models
        # What is known about the sample: {"species": "human", "tissue": "blood", "about": "..."}.
        # Supplied with the dataset (e.g. demo data: confirmed by the dataset) or confirmed later by
        # the user (`confirm_metadata`). Input, not a result: species and tissue set Cell Type
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
        # Compare Conditions (several samples with conditions): the latest comparison, or None.
        self.comparison: dict[str, Any] | None = None

    # --- Where things stand ---------------------------------------------------------------

    @property
    def metadata_confirmed(self) -> bool:
        """Whether the species (and tissue) are known: supplied with the dataset or confirmed."""
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
            **({"comparison": {"status": self.comparison_status(), "settings": self.comparison["settings"],
                               "groups_to_test": self.comparison["groups_to_test"], "ran_at": self.comparison["ran_at"],
                               **conditions.record_summary(self.comparison["summary"], self.params.conditions.max_adj_p)}}
               if self.comparison else {}),
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
            # With several samples (and the setting on), the genes are picked within each sample.
            several = load_data.SAMPLE_KEY in a.obs and a.obs[load_data.SAMPLE_KEY].astype(str).nunique() > 1
            batch_key = load_data.SAMPLE_KEY if params.hvg.per_sample and several else None
            a, hvgs = preprocessing.find_hvgs(a, params.hvg, batch_key=batch_key)
            a, subset = preprocessing.keep_raw_and_subset(a)
            a, scaled = preprocessing.regress_and_scale(a, params.regress_out, params.scale)
            summary = {
                "median_counts_before": norm.median_counts_before,
                "target_sum": norm.target_sum,  # counts per cell after normalization (the value used)
                "target_is_median": norm.target_is_median,
                "n_hvg": hvgs.n_hvg,
                "hvg_per_sample": hvgs.batch_key is not None,  # picked within each sample (several samples only)
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
        integration_params: IntegrationParams | None = None,
    ) -> StageResult:
        """Stage 5: PCA, batch correction (Harmony, if switched on), neighbor graph, UMAP. With
        several batches, each cell's sample mixing score is computed on the uncorrected PCA and,
        with correction, on the corrected coordinates too (`obs["sample_mixing"]` is always the
        coordinates the graph used)."""
        def compute(adata, params):
            a, pcs = dim_reduction.run_pca(adata, params.pca)
            integ, nb = params.integration, params.neighbors
            n_pcs = nb.n_pcs
            has_batches = integ.batch_key in a.obs and a.obs[integ.batch_key].astype(str).nunique() > 1
            if integ.method != "none" and not has_batches:
                integration.check_batch_key(a, integ.batch_key)  # raises with the reason
            correction = None
            before = after = None
            if has_batches:
                labels = a.obs[integ.batch_key].astype(str).to_numpy()
                before = integration.sample_mixing(a.obsm[nb.use_rep][:, :n_pcs], labels)
            if integ.method == "harmony":
                a, correction = integration.run_harmony(a, integ, n_pcs, use_rep=nb.use_rep)
                nb = replace(nb, use_rep=integration.CORRECTED_KEY)
                after = integration.sample_mixing(a.obsm[integration.CORRECTED_KEY], labels)
            a, graph = dim_reduction.compute_neighbors(a, nb)
            a, layout = dim_reduction.run_umap(a, params.umap)
            if has_batches:
                a.obs[integration.MIXING_KEY] = after if after is not None else before
                if after is not None:
                    a.obs[integration.MIXING_BEFORE_KEY] = before
            summary = {
                "n_comps": pcs.n_comps,
                "total_variance_ratio": pcs.total_variance_ratio,
                "elbow_pc": pcs.elbow_pc,  # guidance only (shown on the elbow plot)
                "n_neighbors": graph.n_neighbors,
                "n_connections": graph.n_connections,
            }
            if has_batches:
                summary["integration"] = {
                    "method": integ.method,
                    "batch_key": integ.batch_key,
                    "n_batches": int(a.obs[integ.batch_key].astype(str).nunique()),
                    "graph_from": nb.use_rep,
                    "mixing_before": float(before.mean()),
                    "mixing_after": float(after.mean()) if after is not None else None,
                    "warnings": list(correction.warnings) if correction else [],
                }
            return a, summary, (pcs, graph, layout, correction)
        return self._run("reduce_dimensions", {"pca": pca, "neighbors": neighbors, "umap": umap,
                                                "integration": integration_params}, compute)

    def cluster(self, leiden: LeidenParams | None = None) -> StageResult:
        """Stage 6: Leiden clustering."""
        def compute(adata, params):
            out, s = clustering.cluster_leiden(adata, params.leiden)
            summary = {"n_clusters": s.n_clusters, "cluster_sizes": s.cluster_sizes, "notes": list(s.notes)}
            batch_key = params.integration.batch_key
            if integration.MIXING_KEY in out.obs and batch_key in out.obs:
                # Per cluster: how well its samples mix; clusters below the threshold are flagged.
                table = integration.mixing_by_group(out, params.leiden.key_added, batch_key)
                summary["mixing"] = records(table)
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
        """Stage 8: predict cell types. The marker check compares each cluster's top marker genes
        with the CELLxGENE Census marker table for the sample's species and tissue (over-representation
        test). Without a table for that tissue the marker check doesn't run: every cluster is
        Unresolved with status "no_reference", and a note says so. CellTypist's answer (a second,
        independent opinion) is added either way when a model is chosen. The two answers are then
        combined per cluster (`consensus.py`, Step B): `obs["cell_type_predicted"]` holds the final
        label (or "Unknown"), `obs["cell_type_marker_check"]` the marker check's own label."""
        def compute(adata, params):
            chosen = params.annotation.markers
            census = params.annotation.census
            test = params.annotation.test
            if chosen.species is not None and chosen.species not in marker_db.SPECIES:
                raise annotation.AnnotationError(
                    f"The cell-type references cover {' and '.join(marker_db.SPECIES)} only, so cell "
                    "types can't be predicted for this species. Review the clusters with their marker genes.")
            if not census.tissue:
                raise annotation.AnnotationError("Choose the tissue on Data Overview, so that only cell types "
                                                 "from that tissue are compared.")
            genes = (adata.raw if params.markers.use_raw and adata.raw is not None else adata).var_names
            path = census_markers.table_path(chosen.species, census.tissue, census.release, self.census_dir)

            if path.is_file():
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
                marker_broad = broad_types.census_broad_types(terms, ontology)
                source = {"name": "CELLxGENE Census", "version": provenance["census_version"],
                          "file": path.name, "sha256": provenance["table_sha256"], "tissues": [census.tissue],
                          "cell_ontology_version": provenance.get("cell_ontology_version"),
                          "citation": provenance.get("citation")}
                record = json.loads(out.uns[annotation.RECORD_KEY])
                record["reference"] = source
                out.uns[annotation.RECORD_KEY] = json.dumps(record)
            else:
                match, source, marker_broad = None, None, {}
                out, s = annotation.no_reference(
                    adata, f"There's no CELLxGENE Census reference for {census.tissue} yet, so the marker "
                           "check didn't run. More Census references are being added.")

            clusters = s.clusters
            columns = ["cluster", "n_cells", "predicted_label", "status", "best_candidate", "overlap",
                       "n_markers", "adj_p"]
            summary = {
                "database": source,  # None: no Census reference for this tissue yet
                "no_reference": source is None,
                "n_cell_types_tested": int(s.tests["cell_type"].nunique()),
                "markers_matched_ignoring_case": match.case_matched if match else 0,
                "label_counts": dict(s.label_counts),
                "clusters": records(clusters[columns]),
                "notes": list(s.notes),
                "nothing_matched": bool(source is not None and (clusters["predicted_label"] == annotation.UNRESOLVED).all()),
            }
            # CellTypist's second opinion (shown next to the marker check, never changes its labels).
            ct = None
            model_name = params.annotation.celltypist.model
            if model_name:
                try:
                    counts = self.result("filter")[out.obs_names]
                    per_cell, ct = celltypist_models.predict(counts, out.obs[s.groupby], model_name, self.celltypist_dir)
                    for column in per_cell.columns:
                        out.obs[column] = pd.Categorical(per_cell[column])
                    summary["celltypist"] = {**ct.model, "n_model_genes": ct.n_model_genes,
                                             "n_genes_matched": ct.n_genes_matched, "notes": list(ct.notes),
                                             "clusters": records(ct.clusters)}
                except celltypist_models.CellTypistError as problem:
                    summary["celltypist"] = {"model": model_name, "error": str(problem)}

            # Step B: one final label per cluster from both answers, compared as broad cell types.
            final = consensus.combine(
                clusters, marker_broad, ct.clusters if ct else None,
                broad_types.celltypist_broad_types(model_name) if ct else None,
                ct.gene_match if ct else None, params.annotation.consensus,
                celltypist_problem=None if ct else (consensus.CT_ERROR if model_name else consensus.CT_NOT_RUN))
            final_of = dict(zip(final["cluster"], final["final_label"]))
            out.obs[consensus.MARKER_CHECK_KEY] = out.obs[annotation.PREDICTION_KEY].astype(str).astype("category")
            out.obs[annotation.PREDICTION_KEY] = pd.Categorical(
                out.obs[s.groupby].astype(str).map(final_of), categories=list(dict.fromkeys(final["final_label"])))
            record = json.loads(out.uns[annotation.RECORD_KEY])
            record["consensus"] = {"settings": asdict(params.annotation.consensus), "clusters": records(final)}
            out.uns[annotation.RECORD_KEY] = json.dumps(record)
            summary["final"] = records(final)
            summary["marker_label_counts"] = summary["label_counts"]
            summary["label_counts"] = consensus.label_counts(final)
            summary["agreement_counts"] = {a: int((final["agreement"] == a).sum()) for a in consensus.AGREEMENTS}
            details = {"summary": s, "match": match, "celltypist": ct, "final": final}
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

    def confirm_metadata(self, species: str, tissue: str | None = None,
                         about: str | None = None, stage: str | None = None) -> dict[str, Any]:
        """Record the species and tissue of the sample (and an optional short description), as
        confirmed by the user. `tissue` is a CELLxGENE Census tissue name (e.g. "blood"); `stage` is
        optional ("adult" or "developing"; it picks the suggested CellTypist model). They set
        Cell Type Annotation's settings; if it had already run with others, it needs re-running
        (nothing else depends on them). Recorded in the history. Raises `PipelineError` with a
        plain-language message."""
        if "load" not in self._results:
            raise PipelineError("Load the data before confirming its species and tissue.")
        if self.metadata.get("confirmed_by") == "dataset":
            raise PipelineError("This dataset's species and tissue come with its metadata and can't be changed.")
        metadata = _clean_metadata({"species": species, "tissue": tissue, "about": about, "stage": stage})
        if stage and metadata["stage"] is None:
            raise PipelineError("Choose adult or developing for the stage, or leave it empty.")
        if not metadata["species"]:
            raise PipelineError("Choose the species of your sample.")
        if metadata["species"] in marker_db.SPECIES and not metadata["tissue"]:
            raise PipelineError("Choose the tissue.")
        listed = census_markers.census_tissues(metadata["species"], self.params.annotation.census.release,
                                               self.census_dir) if metadata["tissue"] else []
        if listed and metadata["tissue"] not in {t["tissue"] for t in listed}:
            raise PipelineError(f"Unknown tissue: {metadata['tissue']}.")
        if len(metadata["about"] or "") > MAX_ABOUT_LENGTH:
            raise PipelineError(f"Keep the description under {MAX_ABOUT_LENGTH} characters.")

        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        # Compared with the settings in use: the same choice clears nothing.
        changed = ((metadata["species"], metadata["tissue"], metadata["stage"])
                   != (self.params.annotation.markers.species, self.params.annotation.census.tissue,
                       self.metadata.get("stage")))
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
            "summary": {k: metadata[k] for k in ("species", "tissue", "stage", "about")},
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
        predictions = self._results["annotate"].details["final"]
        predicted = dict(zip(predictions["cluster"].astype(str), predictions["final_label"]))
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

    # --- Compare Conditions (several samples with conditions; not one of the eight stages) ----

    # How the cells can be grouped: obs column -> the stage it needs.
    CONDITION_GROUPINGS = {"leiden": "cluster", annotation.PREDICTION_KEY: "annotate", review.REVIEWED_KEY: "annotate"}

    def condition_groups(self, group_by: str) -> pd.Series:
        """Each clustered cell's group for Compare Conditions: its cluster, its suggested cell type, or
        your label ("Not reviewed" for clusters without a decision). Raises `PipelineError`."""
        if group_by not in self.CONDITION_GROUPINGS:
            raise PipelineError("Group the cells by clusters, suggested cell types or your labels.")
        needed = self.CONDITION_GROUPINGS[group_by]
        if needed not in self._results:
            raise PipelineError("Run Clustering first." if needed == "cluster" else "Run Cell Type Annotation first.")
        obs = self._results["cluster"].adata.obs
        if group_by == "leiden":
            values = obs[self.params.leiden.key_added]
        elif group_by == annotation.PREDICTION_KEY:
            values = self._results["annotate"].adata.obs[annotation.PREDICTION_KEY].reindex(obs.index)
        else:
            values = pd.Series(self.reviewed_labels(), index=obs.index)
        return values.astype(str).rename(group_by)

    def _comparison_fingerprint(self, group_by: str) -> str:
        """What a comparison depends on: the counts (the Filter checkpoint) and each cell's group."""
        groups = self.condition_groups(group_by)
        text = self._results["filter"].ran_at + "\n" + "\n".join(groups.index.astype(str) + "\t" + groups.to_numpy())
        return hashlib.sha256(text.encode()).hexdigest()

    def compare_conditions(self, params: ConditionParams | None = None,
                           groups: list[str] | None = None) -> conditions.ComparisonSummary:
        """Compare two conditions across samples: each group's share of the cells, and (for `groups`,
        None = all) each group's genes, by pseudobulk (`conditions.py`). Needs several samples with
        conditions and Clustering (Cell Type Annotation when grouping by cell types). Recorded in the
        history; never changes any other result. Raises `PipelineError` / `ConditionError`."""
        params = params or self.params.conditions
        if "cluster" not in self._results:
            raise PipelineError("Run Clustering first.")
        obs = self._results["cluster"].adata.obs.copy()
        if load_data.CONDITION_KEY not in obs or obs[load_data.SAMPLE_KEY].astype(str).nunique() < 2:
            raise PipelineError("Compare Conditions needs several samples, each with a condition (set on Data Overview).")
        obs[params.group_by] = self.condition_groups(params.group_by)
        counts = self._results["filter"].adata[obs.index]
        summary = conditions.compare_conditions(counts, obs, load_data.SAMPLE_KEY, load_data.CONDITION_KEY,
                                                params, groups)
        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        self.params = replace(self.params, conditions=params)
        settings = _to_plain(params)
        self.comparison = {"settings": settings, "groups_to_test": list(groups) if groups is not None else None,
                           "ran_at": now, "summary": summary, "fingerprint": self._comparison_fingerprint(params.group_by)}
        self.history.append({
            "stage": "compare_conditions",
            "ran_at": now,
            "settings": {"conditions": settings},
            "summary": conditions.short_summary(summary, params.max_adj_p),
            "cleared_later_stages": [],
        })
        return summary

    def comparison_status(self) -> str:
        """Compare Conditions: "done", "needs re-running" (the cells, counts or groups changed since)
        or "not run"."""
        if self.comparison is None:
            return "not run"
        try:
            current = self._comparison_fingerprint(self.comparison["settings"]["group_by"])
        except PipelineError:
            return "needs re-running"
        return "done" if current == self.comparison["fingerprint"] else "needs re-running"

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



@lru_cache(maxsize=2)
def _cell_ontology(path: str) -> dict[str, dict]:
    """The Cell Ontology, read once per file (it takes about a second)."""
    return census_markers.read_cell_ontology(path)


def _clean_metadata(metadata: dict[str, Any]) -> dict[str, Any]:
    """Species, the Census tissue (or None), the stage ("adult", "developing" or None) and "about"
    (or None), tidied. (Metadata recorded before 2026-11-12 may also hold PanglaoDB "tissues"; they
    are no longer used.)"""
    species = str(metadata.get("species") or "").strip().lower() or None
    tissue = str(metadata.get("tissue") or "").strip().lower() or None
    stage = str(metadata.get("stage") or "").strip().lower() or None
    about = " ".join(str(metadata.get("about") or "").split()) or None
    return {"species": species, "tissue": tissue,
            "stage": stage if stage in celltypist_models.STAGES else None, "about": about}


def _with_metadata(params: AnalysisParams, metadata: dict[str, Any]) -> AnalysisParams:
    """Settings with the sample's species and tissue as Cell Type Annotation's (the metadata
    decides them). Without metadata, the settings are kept as they are."""
    if not metadata.get("species"):
        return params
    markers = replace(params.annotation.markers, species=metadata["species"])
    census = replace(params.annotation.census, tissue=metadata.get("tissue"))
    # CellTypist's model follows the species and tissue (suggested; the user can change it afterwards).
    celltypist = replace(params.annotation.celltypist,
                         model=celltypist_models.suggest_model(metadata["species"], metadata.get("tissue"),
                                                               metadata.get("stage")))
    return replace(params, annotation=replace(params.annotation, markers=markers, census=census, celltypist=celltypist))


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
