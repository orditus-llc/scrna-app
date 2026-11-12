"""Cell-type prediction from marker-gene evidence: an over-representation test.

For each cluster and each cell type of a marker database (see `marker_db.py`): how many of
the cell type's markers are among the cluster's top marker genes, and how unlikely is that
by chance? One-sided Fisher's exact test (the hypergeometric test), as in decoupler's
PanglaoDB tutorial; p-values are adjusted across the cell types tested for each cluster
(Benjamini-Hochberg). The cell type with the smallest adjusted p-value is proposed only if
it passes three checks (`EnrichmentParams`); otherwise the cluster stays "Unresolved" (its
best candidate is still shown):

- too_few_markers: fewer than `min_overlap` of its markers among the top genes
- not_significant: adjusted p-value not below `max_adj_p`
- ambiguous:       a significant runner-up's adjusted p-value is less than `ambiguity_fold`
                   times the best one's

Labels come only from the data and the marker sets. Cluster numbers are never used to
assign labels, and analyst-reviewed labels are never read here (they live in
`validation/` and are used only by tests).

Returns a new AnnData with `obs["cell_type_predicted"]`; the input is never changed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Callable, Iterable, Mapping

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from scrna_analysis.params import EnrichmentParams

UNRESOLVED = "Unresolved"
PREDICTION_KEY = "cell_type_predicted"  # obs column with the predicted labels
RECORD_KEY = "cell_type_prediction"  # uns entry: JSON text with settings and evidence
MARKER_RESULT_KEY = "rank_genes_groups"


class AnnotationError(ValueError):
    """Cell types can't be predicted with this input. The message says why."""


@dataclass(frozen=True)
class AnnotationSummary:
    groupby: str
    clusters: pd.DataFrame  # one row per cluster: prediction, status, best candidate and runner-up
    tests: pd.DataFrame  # one row per cluster per cell type tested (`enrichment_table`)
    scores: pd.DataFrame  # cluster x cell type: -log10(adjusted p-value)
    label_counts: dict[str, int]  # predicted label -> number of cells
    missing_genes: dict[str, tuple[str, ...]]  # cell type -> marker genes not in the data
    skipped_cell_types: tuple[str, ...]  # cell types not tested (too few markers in the data)
    notes: tuple[str, ...] = field(default_factory=tuple)


def enrichment_table(
    adata: ad.AnnData, marker_sets: Mapping[str, Iterable[str]], params: EnrichmentParams
) -> tuple[pd.DataFrame, dict[str, tuple[str, ...]], tuple[str, ...]]:
    """The over-representation test for every cluster and cell type.

    Returns (table, missing genes per cell type, cell types not tested). The table has one row
    per cluster per cell type tested: cluster, cell_type, n_markers (markers among the genes
    tested), overlap (markers among the cluster's top `top_n_genes` marker genes),
    matching_markers (in marker-gene order), p_value and adj_p. Cell types with fewer than
    `min_overlap` markers among the genes tested can never pass, so they aren't tested.
    Doesn't change `adata`.
    """
    _check_enrichment_params(params)
    _marker_settings(adata)
    names = pd.DataFrame(adata.uns[MARKER_RESULT_KEY]["names"])
    universe = set(names.iloc[:, 0])  # every gene the marker test ranked (same for each cluster)
    n_top = min(params.top_n_genes, len(names))

    present, missing, untested = {}, {}, []
    for cell_type, genes in marker_sets.items():
        genes = list(dict.fromkeys(genes))
        absent = tuple(g for g in genes if g not in universe)
        if absent:
            missing[cell_type] = absent
        found = [g for g in genes if g in universe]
        if len(found) >= params.min_overlap:
            present[cell_type] = set(found)
        else:
            untested.append(cell_type)
    if not present:
        raise AnnotationError(
            f"No cell type has at least {params.min_overlap} of its marker genes in this dataset, so no "
            "cell types can be predicted. The marker sets may be for another species or tissue."
        )

    rows = []
    for cluster in names.columns:
        top = list(names[cluster].head(n_top))
        top_set = set(top)
        cluster_rows = []
        for cell_type, markers in present.items():
            overlap = [g for g in top if g in markers]
            # P(at least this many markers among n_top genes drawn from the universe by chance)
            p = float(stats.hypergeom.sf(len(overlap) - 1, len(universe), len(markers), n_top))
            cluster_rows.append({"cluster": str(cluster), "cell_type": cell_type, "n_markers": len(markers),
                                 "overlap": len(top_set & markers), "matching_markers": ", ".join(overlap),
                                 "p_value": p})
        adjusted = stats.false_discovery_control([r["p_value"] for r in cluster_rows], method="bh")
        for row, adj in zip(cluster_rows, adjusted):
            row["adj_p"] = float(adj)
        rows.extend(cluster_rows)
    return pd.DataFrame(rows), missing, tuple(untested)


def predict_cell_types(
    adata: ad.AnnData, marker_sets: Mapping[str, Iterable[str]], params: EnrichmentParams,
    related: Callable[[str, str], bool] | None = None,
) -> tuple[ad.AnnData, AnnotationSummary]:
    """Predict a cell type per cluster with the over-representation test (module docstring).

    `related(a, b)`: for references whose cell types are nested (Census: "naive B cell" is a kind
    of "B cell"), whether two cell types are the same lineage. The runner-up is then the best cell
    type *not* related to the best one, so a parent or child of the best match never makes a
    cluster ambiguous. None (PanglaoDB): every other cell type counts.

    Writes `obs["cell_type_predicted"]` (may contain "Unresolved") and
    `uns["cell_type_prediction"]` (JSON text with the method, settings and per-cluster result).
    `summary.scores` holds -log10(adjusted p-value) per cluster and cell type.
    """
    groupby, _ = _marker_settings(adata)
    table, missing, untested = enrichment_table(adata, marker_sets, params)
    sizes = adata.obs[groupby].astype(str).value_counts()

    rows = []
    for cluster, group in table.groupby("cluster", sort=False):
        ordered = group.sort_values(["adj_p", "overlap", "cell_type"], ascending=[True, False, True])
        best = ordered.iloc[0]
        others = ordered.iloc[1:]
        if related is not None:
            others = others[[not related(best["cell_type"], other) for other in others["cell_type"]]]
        second = others.iloc[0] if len(others) else None
        reasons = []
        if best["overlap"] < params.min_overlap:
            reasons.append("too_few_markers")
        if not best["adj_p"] < params.max_adj_p:
            reasons.append("not_significant")
        if (second is not None and second["adj_p"] < params.max_adj_p
                and second["adj_p"] < best["adj_p"] * params.ambiguity_fold):
            reasons.append("ambiguous")
        rows.append({
            "cluster": cluster,
            "n_cells": int(sizes[cluster]),
            "predicted_label": best["cell_type"] if not reasons else UNRESOLVED,
            "status": reasons[0] if reasons else "confident",
            "reasons": ", ".join(reasons),
            "best_candidate": best["cell_type"],
            "overlap": int(best["overlap"]),
            "n_markers": int(best["n_markers"]),
            "adj_p": float(best["adj_p"]),
            "matching_markers": best["matching_markers"],
            "second_candidate": second["cell_type"] if second is not None else None,
            "second_overlap": int(second["overlap"]) if second is not None else None,
            "second_adj_p": float(second["adj_p"]) if second is not None else None,
        })
    clusters = pd.DataFrame(rows)

    adata = adata.copy()
    label_of = dict(zip(clusters["cluster"], clusters["predicted_label"]))
    adata.obs[PREDICTION_KEY] = pd.Categorical(adata.obs[groupby].astype(str).map(label_of))

    notes = []
    if untested:
        notes.append(f"{len(untested)} cell types not tested: fewer than {params.min_overlap} of their "
                     f"markers are in the data ({', '.join(untested)}).")
    n_missing = len({g for gs in missing.values() for g in gs})
    if n_missing:
        notes.append(f"{n_missing} marker genes aren't in the data.")

    scores = table.pivot(index="cluster", columns="cell_type", values="adj_p")
    scores = -np.log10(scores.clip(lower=np.finfo(float).tiny))
    summary = AnnotationSummary(
        groupby=groupby,
        clusters=clusters,
        tests=table,
        scores=scores.loc[clusters["cluster"], list(dict.fromkeys(table["cell_type"]))],
        label_counts={str(k): int(v) for k, v in adata.obs[PREDICTION_KEY].value_counts().items()},
        missing_genes=missing,
        skipped_cell_types=untested,
        notes=tuple(notes),
    )
    adata.uns[RECORD_KEY] = json.dumps({
        "groupby": groupby,
        "method": "over-representation test",
        "settings": {
            "top_n_genes": params.top_n_genes, "max_adj_p": params.max_adj_p,
            "min_overlap": params.min_overlap, "ambiguity_fold": params.ambiguity_fold,
            "marker_sets": {ct: list(gs) for ct, gs in marker_sets.items()},
        },
        "clusters": _records(clusters),
        "missing_genes": {ct: list(gs) for ct, gs in missing.items()},
        "notes": notes,
    })
    return adata, summary


def _check_enrichment_params(params: EnrichmentParams) -> None:
    if params.top_n_genes < 1:
        raise AnnotationError("The number of top marker genes must be at least 1.")
    if params.min_overlap < 1:
        raise AnnotationError("The minimum number of matching markers must be at least 1.")
    if not 0 < params.max_adj_p <= 1:
        raise AnnotationError("The significance level must be above 0 and at most 1.")
    if params.ambiguity_fold < 1:
        raise AnnotationError("The ambiguity factor must be at least 1.")


def _marker_settings(adata: ad.AnnData) -> tuple[str, bool]:
    """The grouping and expression source used by the marker step (keeps them in sync)."""
    if MARKER_RESULT_KEY not in adata.uns:
        raise AnnotationError("No marker results found. Run rank_marker_genes before predicting cell types.")
    recorded = adata.uns[MARKER_RESULT_KEY]["params"]
    groupby, use_raw = str(recorded["groupby"]), bool(recorded["use_raw"])
    if groupby not in adata.obs:
        raise AnnotationError(f"The marker results refer to obs['{groupby}'], which is missing.")
    if use_raw and adata.raw is None:
        raise AnnotationError("The marker results used adata.raw, but adata.raw is empty.")
    return groupby, use_raw


def _records(table: pd.DataFrame) -> list[dict]:
    """A table as plain rows, numbers at full precision (pandas' to_json rounds them)."""
    rows = table.astype(object).where(table.notna(), None).to_dict(orient="records")
    return [{k: v.item() if hasattr(v, "item") else v for k, v in row.items()} for row in rows]
