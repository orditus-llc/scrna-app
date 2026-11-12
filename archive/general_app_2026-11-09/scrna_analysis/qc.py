"""Quality control: QC metrics, a filter preview, and filtering.

Order (same as the validated reference analysis):

1. `calculate_qc_metrics`: flag mitochondrial genes and calculate per-cell and per-gene
   QC metrics, once, on all genes.
2. `preview_filter`: predict exactly which cells and genes a set of thresholds keeps,
   without changing anything (cheap enough for live UI updates).
3. `apply_filter`: filter a copy of the data with Scanpy, then check the result is
   exactly what the preview predicted.

QC metrics are **not** recalculated after filtering. The upper thresholds and the later
`regress_out` step use the values calculated in step 1, as in the reference.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

from scrna_analysis.params import MITO_PATTERNS, FilterParams, QCParams, SuggestionParams

MITO_COLUMN = "mt"  # var column flagging mitochondrial genes; gives obs["pct_counts_mt"]
CELL_METRICS = ("n_genes_by_counts", "total_counts", "pct_counts_mt")
STATISTICS = ("count", "mean", "std", "min", "25%", "50%", "75%", "max")


class QCError(ValueError):
    """QC or filtering can't be carried out as requested. The message says why."""


@dataclass(frozen=True)
class QCSummary:
    """Result of `calculate_qc_metrics`."""

    n_cells: int
    n_genes: int
    mito_prefix: str
    mito_genes: tuple[str, ...]
    # False when no mitochondrial genes were found: pct_counts_mt is then unavailable,
    # which is different from cells genuinely having 0% mitochondrial counts.
    mito_qc_available: bool
    # metric -> {"count", "mean", "std", "min", "25%", "50%", "75%", "max"}
    metric_stats: dict[str, dict[str, float]]
    warnings: tuple[str, ...] = ()

    @property
    def n_mito_genes(self) -> int:
        return len(self.mito_genes)


@dataclass(frozen=True)
class FilterSuggestion:
    """Filter thresholds suggested from a dataset's own QC metrics, with how each was found.

    `details` has one entry per threshold ("min_genes", "max_n_genes_by_counts",
    "max_pct_counts_mt"): the value, the median and MAD it came from (on the scale used), the
    number of MADs, and a note when no limit could be suggested.
    """

    filter: FilterParams
    details: dict[str, dict]


def suggest_filter(adata: ad.AnnData, params: SuggestionParams) -> FilterSuggestion:
    """Suggest filter thresholds from the QC metrics of `adata` (after `calculate_qc_metrics`).

    A cell is an outlier when its metric is more than N median absolute deviations (MADs) from
    the median (Heumos et al. 2023):
    - genes per cell, on a log1p scale: below median - N MADs -> `min_genes` (rounded up to a
      whole number, inclusive), above median + N MADs -> `max_n_genes_by_counts` (a whole
      number; kept cells have fewer genes than it);
    - % mitochondrial, above median + N MADs -> `max_pct_counts_mt` (rounded up to 0.1%).
      No limit if no mitochondrial genes are flagged.
    A metric with no spread (MAD = 0) gets no limit, with a note. `min_cells` is not
    data-driven; it comes from `params`.
    """
    for column in ("n_genes_by_counts", "pct_counts_mt"):
        if column not in adata.obs:
            raise QCError("QC metrics are missing: calculate QC metrics before suggesting thresholds.")

    details: dict[str, dict] = {}
    log_genes = np.log1p(adata.obs["n_genes_by_counts"].to_numpy(dtype=float))
    median, mad = float(np.median(log_genes)), float(median_abs_deviation(log_genes))
    k = params.n_mads_genes
    min_genes, max_genes = 0, None
    if mad > 0:
        min_genes = max(0, int(np.ceil(np.expm1(median - k * mad))))
        max_genes = int(np.floor(np.expm1(median + k * mad))) + 1
        note = None
    else:
        note = "Genes per cell barely vary, so no limit is suggested."
    common = {"median": median, "mad": mad, "n_mads": k, "scale": "log1p", "note": note}
    details["min_genes"] = {"value": min_genes, **common}
    details["max_n_genes_by_counts"] = {"value": max_genes, **common}

    max_mt = None
    if MITO_COLUMN in adata.var and adata.var[MITO_COLUMN].any():
        pct = adata.obs["pct_counts_mt"].to_numpy(dtype=float)
        median, mad = float(np.median(pct)), float(median_abs_deviation(pct))
        if mad > 0:
            max_mt = float(np.ceil((median + params.n_mads_mito * mad) * 10) / 10)
            note = None
        else:
            note = "% mitochondrial barely varies, so no limit is suggested."
        details["max_pct_counts_mt"] = {"value": max_mt, "median": median, "mad": mad,
                                        "n_mads": params.n_mads_mito, "scale": "percent", "note": note}
    else:
        details["max_pct_counts_mt"] = {"value": None, "median": None, "mad": None, "n_mads": params.n_mads_mito,
                                        "scale": "percent",
                                        "note": "No mitochondrial genes found, so no limit is suggested."}

    return FilterSuggestion(
        filter=FilterParams(min_genes=min_genes, min_cells=params.min_cells,
                            max_n_genes_by_counts=max_genes, max_pct_counts_mt=max_mt),
        details=details,
    )


@dataclass(frozen=True)
class FilterPreview:
    """Which cells and genes a set of thresholds keeps. Computed without changing data."""

    n_cells_before: int
    n_genes_before: int
    kept_cells: pd.Index  # cell barcodes, in their original order
    kept_genes: pd.Index  # gene names, in their original order
    # Removed at each step, in the order the filters are applied (no double counting):
    # "cells: min_genes", "genes: min_cells", "cells: max_n_genes_by_counts",
    # "cells: max_pct_counts_mt". Disabled filters are left out.
    removed_in_order: dict[str, int] = field(default_factory=dict)
    # Cells failing each cell criterion on its own (these can overlap).
    cells_failing_each: dict[str, int] = field(default_factory=dict)

    @property
    def n_cells_kept(self) -> int:
        return len(self.kept_cells)

    @property
    def n_genes_kept(self) -> int:
        return len(self.kept_genes)

    def counts(self) -> dict:
        """The numbers only (no cell/gene names), JSON-serializable for recording."""
        return {
            "n_cells_before": self.n_cells_before,
            "n_genes_before": self.n_genes_before,
            "n_cells_kept": self.n_cells_kept,
            "n_genes_kept": self.n_genes_kept,
            "removed_in_order": dict(self.removed_in_order),
            "cells_failing_each": dict(self.cells_failing_each),
        }


def genes_with_prefix(var_names, prefix: str) -> tuple[str, ...]:
    """The gene names that start with `prefix` (case-sensitive: human "MT-" isn't mouse "mt-")."""
    if not prefix:
        return ()
    names = pd.Index(var_names).astype(str)
    return tuple(names[names.str.startswith(prefix)])


def mito_prefix_matches(var_names, patterns=MITO_PATTERNS) -> list[dict]:
    """How many genes match each known mitochondrial prefix, most matches first.

    Returns one entry per prefix with at least one match: {"prefix", "species", "n_genes",
    "examples" (up to 5 names)}. Empty if no known prefix matches (e.g. the mitochondrial
    genes were removed, or the species names them differently).
    """
    found = []
    for prefix, species in patterns.items():
        genes = genes_with_prefix(var_names, prefix)
        if genes:
            found.append({"prefix": prefix, "species": species, "n_genes": len(genes), "examples": list(genes[:5])})
    return sorted(found, key=lambda m: -m["n_genes"])


def suggest_mito_prefix(var_names, patterns=MITO_PATTERNS) -> str | None:
    """The known prefix that matches the most genes, or None if none matches."""
    matches = mito_prefix_matches(var_names, patterns)
    return matches[0]["prefix"] if matches else None


def calculate_qc_metrics(adata: ad.AnnData, params: QCParams) -> tuple[ad.AnnData, QCSummary]:
    """Flag mitochondrial genes and calculate QC metrics on a copy of `adata`.

    Adds `var["mt"]`, the per-cell metrics (`n_genes_by_counts`, `total_counts`,
    `total_counts_mt`, `pct_counts_mt`) and the per-gene metrics (`n_cells_by_counts`,
    `mean_counts`, `pct_dropout_by_counts`, `total_counts`).

    If no gene names start with `params.mito_prefix`, a warning is issued and the
    summary marks mitochondrial QC as unavailable. The calculation itself still runs.
    """
    if set(params.qc_vars) != {MITO_COLUMN}:
        raise QCError(
            f"Only the mitochondrial QC variable ('{MITO_COLUMN}') is supported for now; "
            f"got qc_vars={list(params.qc_vars)}."
        )

    adata = adata.copy()
    adata.var[MITO_COLUMN] = adata.var_names.str.startswith(params.mito_prefix)
    mito_genes = tuple(adata.var_names[adata.var[MITO_COLUMN]])

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=list(params.qc_vars),
        percent_top=list(params.percent_top) if params.percent_top is not None else None,
        log1p=params.log1p,
        inplace=True,
    )

    notes = []
    metrics = list(CELL_METRICS)
    if not mito_genes:
        message = (
            f"No mitochondrial genes found: no gene names start with "
            f"'{params.mito_prefix}', so % mitochondrial can't be measured (it isn't a real 0%). "
            f"Check the species or the gene prefix, or switch off the % mitochondrial filter."
        )
        warnings.warn(message, UserWarning, stacklevel=2)
        notes.append(message)
        metrics.remove("pct_counts_mt")

    summary = QCSummary(
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        mito_prefix=params.mito_prefix,
        mito_genes=mito_genes,
        mito_qc_available=bool(mito_genes),
        metric_stats={m: _describe(adata.obs[m]) for m in metrics},
        warnings=tuple(notes),
    )
    return adata, summary


def preview_filter(adata: ad.AnnData, params: FilterParams) -> FilterPreview:
    """Predict exactly which cells and genes `apply_filter` keeps. Changes nothing.

    Follows the reference order: (1) cells with fewer than `min_genes` detected genes
    are removed; (2) genes detected in fewer than `min_cells` of the *remaining* cells
    are removed; (3) cells with `n_genes_by_counts >= max_n_genes_by_counts` and
    (4) cells with `pct_counts_mt >= max_pct_counts_mt` are removed. Steps 3–4 use the
    QC metrics from `calculate_qc_metrics` and are skipped when their limit is `None`.

    Raises `QCError` if QC metrics are missing, if the mitochondrial filter is on but
    no mitochondrial genes are flagged, or if no cells or genes would remain.
    """
    _check_ready_to_filter(adata, params)

    pass_min_genes, below_max_genes, below_max_mt = _cell_checks(adata, params)
    # Step 2: the same count Scanpy's filter_genes uses (detected = > 0), on the remaining cells.
    cells_per_gene = np.asarray((adata.X[pass_min_genes] > 0).sum(axis=0)).ravel()
    keep_genes = cells_per_gene >= params.min_cells
    keep_cells = pass_min_genes & below_max_genes & below_max_mt

    removed = {
        "cells: min_genes": int((~pass_min_genes).sum()),
        "genes: min_cells": int((~keep_genes).sum()),
    }
    failing = {"min_genes": int((~pass_min_genes).sum())}
    if params.max_n_genes_by_counts is not None:
        removed["cells: max_n_genes_by_counts"] = int((pass_min_genes & ~below_max_genes).sum())
        failing["max_n_genes_by_counts"] = int((~below_max_genes).sum())
    if params.max_pct_counts_mt is not None:
        removed["cells: max_pct_counts_mt"] = int(
            (pass_min_genes & below_max_genes & ~below_max_mt).sum()
        )
        failing["max_pct_counts_mt"] = int((~below_max_mt).sum())

    preview = FilterPreview(
        n_cells_before=adata.n_obs,
        n_genes_before=adata.n_vars,
        kept_cells=adata.obs_names[keep_cells],
        kept_genes=adata.var_names[keep_genes],
        removed_in_order=removed,
        cells_failing_each=failing,
    )
    if preview.n_cells_kept == 0 or preview.n_genes_kept == 0:
        raise QCError(
            f"These settings would remove all {'cells' if preview.n_cells_kept == 0 else 'genes'} "
            f"({preview.n_cells_kept:,} cells × {preview.n_genes_kept:,} genes would remain)."
        )
    return preview


def cell_filter_failures(adata: ad.AnnData, params: FilterParams) -> pd.DataFrame:
    """Which cell filters each cell fails, for inspecting cells. Changes nothing.

    One row per cell (index = cell barcodes, in their original order) and one True/False
    column per cell filter that is switched on: `min_genes`, and `max_n_genes_by_counts` /
    `max_pct_counts_mt` unless their limit is None. Each column is tested on its own (a cell
    can fail several), with exactly the rules `preview_filter` uses, so a cell is kept by
    the filter if and only if its row has no True. The gene filter (`min_cells`) removes
    genes, not cells, so it has no column. Raises `QCError` like `preview_filter`.
    """
    _check_ready_to_filter(adata, params)
    pass_min_genes, below_max_genes, below_max_mt = _cell_checks(adata, params)
    failures = {"min_genes": ~pass_min_genes}
    if params.max_n_genes_by_counts is not None:
        failures["max_n_genes_by_counts"] = ~below_max_genes
    if params.max_pct_counts_mt is not None:
        failures["max_pct_counts_mt"] = ~below_max_mt
    return pd.DataFrame(failures, index=adata.obs_names)


def _cell_checks(adata: ad.AnnData, params: FilterParams) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The three cell tests, shared by `preview_filter` and `cell_filter_failures`.

    Step 1 counts detected genes (> 0) like Scanpy's filter_cells; steps 3-4 compare the
    stored QC metrics with the upper limits (strict "<", as in the reference; a None
    limit keeps every cell).
    """
    detected = adata.X > 0
    genes_per_cell = np.asarray(detected.sum(axis=1)).ravel()
    pass_min_genes = genes_per_cell >= params.min_genes
    below_max_genes = _below(adata.obs["n_genes_by_counts"], params.max_n_genes_by_counts)
    below_max_mt = _below(adata.obs["pct_counts_mt"], params.max_pct_counts_mt)
    return pass_min_genes, below_max_genes, below_max_mt


def apply_filter(adata: ad.AnnData, params: FilterParams) -> tuple[ad.AnnData, FilterPreview]:
    """Filter a copy of `adata` and return it with a summary of what was removed.

    Uses the same Scanpy calls as the reference: `filter_cells(min_genes)`,
    `filter_genes(min_cells)`, then the upper thresholds on the stored QC metrics.
    The result is checked against `preview_filter`: the kept cell barcodes and gene
    names must be identical, in the same order.
    """
    preview = preview_filter(adata, params)  # also runs all input checks

    filtered = adata.copy()
    sc.pp.filter_cells(filtered, min_genes=params.min_genes)
    sc.pp.filter_genes(filtered, min_cells=params.min_cells)
    keep = _below(filtered.obs["n_genes_by_counts"], params.max_n_genes_by_counts) & _below(
        filtered.obs["pct_counts_mt"], params.max_pct_counts_mt
    )
    filtered = filtered[keep, :].copy()

    if not (filtered.obs_names.equals(preview.kept_cells)
            and filtered.var_names.equals(preview.kept_genes)):
        raise RuntimeError(
            "Internal error: filtered cells/genes differ from preview_filter's prediction."
        )
    return filtered, preview


def _check_ready_to_filter(adata: ad.AnnData, params: FilterParams) -> None:
    missing = [m for m in ("n_genes_by_counts", "pct_counts_mt") if m not in adata.obs]
    if missing or MITO_COLUMN not in adata.var:
        raise QCError("QC metrics are missing. Run calculate_qc_metrics before filtering.")
    if params.max_pct_counts_mt is not None and not adata.var[MITO_COLUMN].any():
        raise QCError(
            f"The % mitochondrial filter is on ({params.max_pct_counts_mt:g}%), but no "
            "mitochondrial genes were found, so % mitochondrial can't be measured (it isn't a "
            "real 0%). Check the species or the gene prefix and run QC again, or switch off "
            "the % mitochondrial filter."
        )


def _below(values: pd.Series, limit: float | None) -> np.ndarray:
    """True where value < limit; all True when the limit is switched off (None)."""
    if limit is None:
        return np.ones(len(values), dtype=bool)
    return values.to_numpy() < limit


def _describe(values: pd.Series) -> dict[str, float]:
    """count/mean/std/min/quartiles/max, as pandas `describe()` (same as the notebook)."""
    stats = values.describe()
    return {name: float(stats[name]) for name in STATISTICS}
