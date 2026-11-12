"""Translators between the screen and the backend.

Turns what the user typed into backend settings, and backend results into tables and
labels for display. No calculations and no Shiny code, so it can be tested on its own.
All default values come from `pbmc_analysis.params`; nothing here is PBMC3k-specific.
"""

from __future__ import annotations

import pandas as pd

from dataclasses import replace

import anndata as ad

from pbmc_analysis.markers import marker_table
from pbmc_analysis.params import REFERENCE_PRESET, AnalysisParams, FilterParams
from pbmc_analysis.qc import FilterPreview, QCSummary

P = REFERENCE_PRESET  # every default shown on screen comes from here

# How the eight stages and their statuses are named on screen.
STAGE_LABELS = {
    "load": "Load",
    "qc_metrics": "QC metrics",
    "filter": "Filter",
    "preprocess": "Preprocess",
    "reduce_dimensions": "PCA / UMAP",
    "cluster": "Cluster",
    "markers": "Markers",
    "annotate": "Annotate",
}

METRIC_LABELS = {
    "n_genes_by_counts": "Genes per cell",
    "total_counts": "Counts per cell",
    "pct_counts_mt": "% mitochondrial",
}

STATISTIC_LABELS = {
    "count": "Cells", "mean": "Mean", "std": "Std. dev.", "min": "Min",
    "25%": "25th percentile", "50%": "Median", "75%": "75th percentile", "max": "Max",
}

REMOVAL_LABELS = {
    "cells: min_genes": ("Cells", "too few genes detected"),
    "genes: min_cells": ("Genes", "detected in too few cells"),
    "cells: max_n_genes_by_counts": ("Cells", "too many genes (possible doublets)"),
    "cells: max_pct_counts_mt": ("Cells", "high % mitochondrial counts"),
}


class InputError(ValueError):
    """A value on screen can't be used. The message is written for the user."""


def default_filter_inputs(params: FilterParams = REFERENCE_PRESET.filter) -> dict:
    """The filter boxes' starting values, taken from the settings."""
    return {
        "min_genes": params.min_genes,
        "min_cells": params.min_cells,
        "max_genes": params.max_n_genes_by_counts,
        "max_genes_off": params.max_n_genes_by_counts is None,
        "max_mt": params.max_pct_counts_mt,
        "max_mt_off": params.max_pct_counts_mt is None,
    }


def filter_params_from_inputs(
    min_genes, min_cells, max_genes, max_genes_off: bool, max_mt, max_mt_off: bool
) -> FilterParams:
    """Turn the filter boxes into backend settings. A ticked "no limit" box becomes None.

    Raises `InputError` (with a message for the user) for empty or negative values.
    """
    def required_number(value, label, whole=True):
        if value is None:
            raise InputError(f"Enter a value for '{label}'.")
        if value < 0:
            raise InputError(f"'{label}' can't be negative.")
        if whole and value != int(value):
            raise InputError(f"'{label}' must be a whole number.")
        return int(value) if whole else float(value)

    return FilterParams(
        min_genes=required_number(min_genes, "Min genes per cell"),
        min_cells=required_number(min_cells, "Min cells per gene"),
        max_n_genes_by_counts=None if max_genes_off else required_number(max_genes, "Max genes per cell"),
        max_pct_counts_mt=None if max_mt_off else required_number(max_mt, "Max % mitochondrial", whole=False),
    )


def qc_stats_table(summary: QCSummary) -> pd.DataFrame:
    """The QC summary statistics (as in the notebook's describe() table), for display."""
    table = pd.DataFrame(summary.metric_stats).rename(columns=METRIC_LABELS)
    table.index = [STATISTIC_LABELS[name] for name in table.index]
    table = table.round(2)
    return table.reset_index(names="Statistic")


def removal_table(preview: FilterPreview) -> pd.DataFrame:
    """What each filter step removes, in the order the steps are applied."""
    rows = [
        {"Step": index + 1, "What": REMOVAL_LABELS[step][0], "Reason": REMOVAL_LABELS[step][1],
         "Removed": removed}
        for index, (step, removed) in enumerate(preview.removed_in_order.items())
    ]
    return pd.DataFrame(rows)


def html_table(table: pd.DataFrame, decimals: int = 2) -> str:
    """A table as HTML with each header aligned like its column.

    Text columns are left-aligned and number columns right-aligned, in both the header
    and the body, so values always sit under their own header. Whole numbers are shown
    without decimals (e.g. 2,700 instead of 2700.00).
    """
    def cell(value, numeric):
        if not numeric:
            return str(value)
        if float(value).is_integer():
            return f"{int(value):,}"
        return f"{value:,.{decimals}f}"

    numeric = {col: pd.api.types.is_numeric_dtype(table[col]) for col in table.columns}
    align = {col: "text-end" if numeric[col] else "text-start" for col in table.columns}
    header = "".join(f'<th class="{align[c]}">{c}</th>' for c in table.columns)
    rows = "".join(
        "<tr>" + "".join(f'<td class="{align[c]}">{cell(row[c], numeric[c])}</td>'
                         for c in table.columns) + "</tr>"
        for _, row in table.iterrows()
    )
    return (f'<table class="table table-sm table-striped">'
            f"<thead><tr>{header}</tr></thead><tbody>{rows}</tbody></table>")


# --- Later stages: settings boxes -> backend settings ------------------------------------
#
# Each function starts from the current settings and replaces only what the screen shows,
# so anything not on screen (e.g. seeds) keeps its reference value.

def _number(value, label, *, positive=False, whole=False, at_most=None):
    if value is None:
        raise InputError(f"Enter a value for '{label}'.")
    if value < 0 or (positive and value == 0):
        raise InputError(f"'{label}' must be {'greater than 0' if positive else 'at least 0'}.")
    if whole and value != int(value):
        raise InputError(f"'{label}' must be a whole number.")
    if at_most is not None and value > at_most:
        raise InputError(f"'{label}' can be at most {at_most}.")
    return int(value) if whole else float(value)


def default_stage_inputs(params: AnalysisParams = P) -> dict:
    """Starting values of every later-stage box, taken from the settings."""
    return {
        "target_sum": params.normalize.target_sum,
        "hvg_min_mean": params.hvg.min_mean,
        "hvg_max_mean": params.hvg.max_mean,
        "hvg_min_disp": params.hvg.min_disp,
        "regress": params.regress_out.enabled,
        "scale_max": params.scale.max_value,
        "n_comps": params.pca.n_comps,
        "n_neighbors": params.neighbors.n_neighbors,
        "n_pcs": params.neighbors.n_pcs,
        "resolution": params.leiden.resolution,
        "marker_method": params.markers.method,
        "ambiguity_margin": params.annotation.ambiguity_margin,
        "max_de_rank": params.annotation.max_de_rank,
        "conflict_z": params.annotation.conflict_z,
    }


def preprocess_settings(current: AnalysisParams, target_sum, min_mean, max_mean, min_disp,
                        regress: bool, scale_max) -> dict:
    """Settings for `run.preprocess(...)` from the Preprocessing boxes."""
    low, high = _number(min_mean, "HVG min mean"), _number(max_mean, "HVG max mean")
    if low >= high:
        raise InputError("'HVG min mean' must be smaller than 'HVG max mean'.")
    return {
        "normalize": replace(current.normalize, target_sum=_number(target_sum, "Target sum", positive=True)),
        "hvg": replace(current.hvg, min_mean=low, max_mean=high,
                       min_disp=_number(min_disp, "HVG min dispersion")),
        "regress_out": replace(current.regress_out, enabled=bool(regress)),
        "scale": replace(current.scale, max_value=_number(scale_max, "Scaling limit", positive=True)),
    }


def reduce_settings(current: AnalysisParams, n_comps, n_neighbors, n_pcs) -> dict:
    """Settings for `run.reduce_dimensions(...)` from the PCA & UMAP boxes."""
    comps = _number(n_comps, "Number of PCs", positive=True, whole=True)
    used = _number(n_pcs, "PCs used for the neighbor graph", positive=True, whole=True, at_most=comps)
    return {
        "pca": replace(current.pca, n_comps=comps),
        "neighbors": replace(current.neighbors, n_neighbors=_number(
            n_neighbors, "Number of neighbors", positive=True, whole=True), n_pcs=used),
    }


def cluster_settings(current: AnalysisParams, resolution) -> dict:
    return {"leiden": replace(current.leiden, resolution=_number(resolution, "Resolution", positive=True))}


def marker_settings(current: AnalysisParams, method: str) -> dict:
    return {"marker_params": replace(current.markers, method=method)}


def annotation_settings(current: AnalysisParams, margin, max_rank, conflict_z) -> dict:
    return {"annotation_params": replace(
        current.annotation,
        ambiguity_margin=_number(margin, "Ambiguity margin"),
        max_de_rank=_number(max_rank, "Top-N marker rank", positive=True, whole=True),
        conflict_z=_number(conflict_z, "Conflict threshold"),
    )}


# --- Later stages: results -> display tables ----------------------------------------------

def cluster_sizes_table(sizes: dict[str, int]) -> pd.DataFrame:
    total = sum(sizes.values())
    return pd.DataFrame({
        "Cluster": list(sizes),
        "Cells": list(sizes.values()),
        "% of cells": [round(100 * n / total, 1) for n in sizes.values()],
    })


def marker_genes_table(adata: ad.AnnData, cluster: str, n_top: int = 10) -> pd.DataFrame:
    """Top marker genes of one cluster, in plain terms (percentages, not fractions)."""
    table = marker_table(adata, group=cluster, n_top=n_top)
    out = pd.DataFrame({
        "Rank": table["rank"],
        "Gene": table["gene"],
        "Score": table["score"].round(2),
        "Log fold change": table["logfoldchange"].round(2),
    })
    if "frac_in_group" in table:
        out["% in cluster"] = (100 * table["frac_in_group"]).round(1)
        out["% in other cells"] = (100 * table["frac_in_rest"]).round(1)
    out["Adjusted p (ranking only)"] = table["pval_adj"].map(lambda p: f"{p:.1e}")
    return out


ANNOTATION_COLUMNS = {
    "cluster": "Cluster",
    "n_cells": "Cells",
    "predicted_label": "Predicted cell type",
    "status": "Status",
    "best_candidate": "Best candidate",
    "score": "Score",
    "margin": "Margin to 2nd",
    "supporting_markers": "Supporting markers",
    "conflicting_markers": "Conflicting markers",
}


def annotation_table(clusters: pd.DataFrame) -> pd.DataFrame:
    """Per-cluster prediction with its evidence, for display."""
    table = clusters[list(ANNOTATION_COLUMNS)].rename(columns=ANNOTATION_COLUMNS).copy()
    table["Score"] = table["Score"].round(2)
    table["Margin to 2nd"] = table["Margin to 2nd"].round(2)
    table["Conflicting markers"] = table["Conflicting markers"].replace("", "none")
    return table


def uses_reference_settings(params: AnalysisParams) -> bool:
    """True while every setting equals the validated reference preset."""
    return params == REFERENCE_PRESET


def short_fingerprint(sha256: str | None) -> str:
    """First 12 characters of a file checksum: enough to tell files apart on screen."""
    return sha256[:12] if sha256 else "unknown"
