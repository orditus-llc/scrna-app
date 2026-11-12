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

`run_analysis(path)` runs all eight stages with the reference settings in one call.
All computation happens in the stage modules; this file only orders and records.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Callable

import anndata as ad
import pandas as pd

from pbmc_analysis import annotation, clustering, dim_reduction, load_data, markers, preprocessing, qc
from pbmc_analysis.params import (
    REFERENCE_PRESET,
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
    "qc_metrics": ("qc",),
    "filter": ("filter",),
    "preprocess": ("normalize", "hvg", "regress_out", "scale"),
    "reduce_dimensions": ("pca", "neighbors", "umap"),
    "cluster": ("leiden",),
    "markers": ("markers",),
    "annotate": ("annotation",),
}

RECORD_KEY = "orditus"  # adata.uns entry holding the record (JSON text, h5ad-safe)
RECORD_FORMAT = 1  # bump if the record's structure changes
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
    details: tuple[Any, ...]  # the stage modules' full summaries (for display)
    settings: dict[str, Any]  # the settings this stage used
    ran_at: str  # date and time (UTC)


class AnalysisRun:
    """One analysis of one dataset: runs stages in order, keeps checkpoints, records all.

    Example
    -------
    >>> run = AnalysisRun("data/pbmc3k.h5ad")       # reference settings by default
    >>> run.load(); run.compute_qc_metrics()
    >>> run.preview_filter()                          # look before filtering
    >>> run.apply_filter(); run.preprocess(); run.reduce_dimensions()
    >>> run.cluster(LeidenParams(resolution=0.5))     # any stage can take new settings
    >>> run.status()                                  # what's done / needs re-running
    """

    def __init__(self, data_path: str | Path, params: AnalysisParams = REFERENCE_PRESET):
        self.data_path = Path(data_path)
        self.params = params  # the settings currently in use (updated when a stage changes them)
        self.history: list[dict[str, Any]] = []  # one entry per stage run, in order
        self._results: dict[str, StageResult] = {}  # checkpoints of completed stages
        self._needs_rerun: set[str] = set()  # stages cleared because an earlier one was redone

    # --- Where things stand ---------------------------------------------------------------

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
            "input": {"path": str(self.data_path), "sha256": load_summary.get("input_sha256")},
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

        all_settings = params.to_dict()
        settings = {name: all_settings[name] for name in STAGE_SETTINGS[stage]}
        ran_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
        self.history.append({
            "stage": stage,
            "ran_at": ran_at,
            "settings": settings,
            "summary": summary,
            "cleared_later_stages": cleared,
        })

        result = StageResult(stage, adata, summary, details, settings, ran_at)
        self._results[stage] = result
        adata.uns[RECORD_KEY] = json.dumps(self.record())  # adata is this stage's own new copy
        return result

    # --- The eight stages (each calls an already-validated module) -------------------------

    def load(self) -> StageResult:
        """Stage 1: open the data file and check it holds raw counts."""
        def compute(_, params):
            adata, s = load_data.load_h5ad(self.data_path)
            summary = {
                "n_cells": s.n_cells,
                "n_genes": s.n_genes,
                "n_duplicate_gene_names": s.n_duplicate_gene_names,
                "input_sha256": file_sha256(self.data_path),
            }
            return adata, summary, (s,)
        return self._run("load", {}, compute)

    def compute_qc_metrics(self, qc_params: QCParams | None = None) -> StageResult:
        """Stage 2: flag mitochondrial genes and calculate QC metrics."""
        def compute(adata, params):
            out, s = qc.calculate_qc_metrics(adata, params.qc)
            summary = {
                "n_mito_genes": s.n_mito_genes,
                "mito_qc_available": s.mito_qc_available,
                "warnings": list(s.warnings),
            }
            return out, summary, (s,)
        return self._run("qc_metrics", {"qc": qc_params}, compute)

    def preview_filter(self, filter_params: FilterParams | None = None) -> qc.FilterPreview:
        """Look before filtering: which cells/genes these thresholds would keep.

        Changes nothing (not a stage; nothing is recorded).
        """
        if "qc_metrics" not in self._results:
            raise PipelineError("Can't preview filtering yet: run 'qc_metrics' first.")
        return qc.preview_filter(self.result("qc_metrics"), filter_params or self.params.filter)

    def apply_filter(self, filter_params: FilterParams | None = None) -> StageResult:
        """Stage 3: remove low-quality cells and rarely detected genes."""
        def compute(adata, params):
            out, preview = qc.apply_filter(adata, params.filter)
            return out, preview.counts(), (preview,)
        return self._run("filter", {"filter": filter_params}, compute)

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
        """Stage 8: predict cell types from marker evidence."""
        def compute(adata, params):
            out, s = annotation.predict_cell_types(adata, params.annotation)
            columns = ["cluster", "n_cells", "predicted_label", "status", "best_candidate",
                       "score", "margin", "conflicting_markers"]
            summary = {
                "label_counts": s.label_counts,
                "clusters": _records(s.clusters[columns]),
                "missing_genes": {ct: list(gs) for ct, gs in s.missing_genes.items()},
                "notes": list(s.notes),
            }
            return out, summary, (s,)
        return self._run("annotate", {"annotation": annotation_params}, compute)

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


def run_analysis(data_path: str | Path, params: AnalysisParams = REFERENCE_PRESET) -> AnalysisRun:
    """Run the whole analysis (all eight stages) on `data_path` and return the run."""
    return AnalysisRun(data_path, params).run_all()


# --- Small helpers ------------------------------------------------------------------------

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


def _records(table: pd.DataFrame) -> list[dict[str, Any]]:
    """A table as a list of plain rows (JSON-safe numbers and text)."""
    return json.loads(table.to_json(orient="records"))
