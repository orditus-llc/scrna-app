"""Cell-type prediction from marker-gene evidence.

For each cluster, the expression of every marker gene in the configured marker sets is
compared across clusters (z-score of the cluster means), and each cell type is scored as
the mean z-score of its markers. The best-scoring cell type is proposed only if it passes
three checks; otherwise the cluster stays "Unresolved" (its best candidate is still shown):

- ambiguous:   best score - second-best score < `ambiguity_margin`
- unsupported: none of the cell type's markers is in the cluster's top `max_de_rank`
               marker genes (from `rank_marker_genes`)
- low_score:   best score < `min_score` (only if `min_score` is set)

Labels come only from the data and the marker sets. Cluster numbers are never used to
assign labels, and analyst-reviewed labels are never read here (they live in
`validation/` and are used only by tests).

Returns a new AnnData with `obs["cell_type_predicted"]`; the input is never changed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from pbmc_analysis.params import AnnotationParams

UNRESOLVED = "Unresolved"
CONTEXT = "context"  # cell_type value for context markers (shown, not scored)
PREDICTION_KEY = "cell_type_predicted"  # obs column with the predicted labels
RECORD_KEY = "cell_type_prediction"  # uns entry: JSON text with settings and evidence
MARKER_RESULT_KEY = "rank_genes_groups"


class AnnotationError(ValueError):
    """Cell types can't be predicted with this input. The message says why."""


@dataclass(frozen=True)
class AnnotationSummary:
    groupby: str
    clusters: pd.DataFrame  # one row per cluster: prediction, status, evidence
    scores: pd.DataFrame  # cluster x cell type scores
    label_counts: dict[str, int]  # predicted label -> number of cells
    missing_genes: dict[str, tuple[str, ...]]  # cell type -> marker genes not in the data
    skipped_cell_types: tuple[str, ...]  # cell types with none of their markers in the data
    notes: tuple[str, ...] = field(default_factory=tuple)


def cluster_marker_evidence(adata: ad.AnnData, params: AnnotationParams) -> pd.DataFrame:
    """Tidy table of marker evidence: one row per cluster per marker gene.

    Columns: cluster, cell_type ("context" for context markers), gene, mean_expr
    (mean log-normalized expression), frac_expressing, z_score (cluster mean compared
    across clusters), de_rank (1 = top marker gene of the cluster) and logfoldchange.
    Only genes present in the data are included. Doesn't change `adata`.
    """
    groupby, use_raw = _marker_settings(adata)
    genes, _, _ = _present_genes(adata, params, use_raw)
    expr = sc.get.obs_df(adata, keys=genes, use_raw=use_raw)
    clusters = adata.obs[groupby]

    mean_expr = expr.groupby(clusters, observed=True).mean()
    frac = (expr > 0).groupby(clusters, observed=True).mean()
    z = _zscore(mean_expr)
    ranks, fold_changes = _de_lookup(adata)

    owner = _gene_owner(params)
    rows = []
    for cluster in mean_expr.index:
        for gene in genes:
            for cell_type in owner.get(gene, [CONTEXT]):
                rows.append({
                    "cluster": str(cluster),
                    "cell_type": cell_type,
                    "gene": gene,
                    "mean_expr": float(mean_expr.loc[cluster, gene]),
                    "frac_expressing": float(frac.loc[cluster, gene]),
                    "z_score": float(z.loc[cluster, gene]),
                    "de_rank": ranks[str(cluster)].get(gene),
                    "logfoldchange": fold_changes[str(cluster)].get(gene),
                })
    return pd.DataFrame(rows)


def predict_cell_types(
    adata: ad.AnnData, params: AnnotationParams
) -> tuple[ad.AnnData, AnnotationSummary]:
    """Predict a cell type per cluster from marker evidence (see module docstring).

    Writes `obs["cell_type_predicted"]` (categorical, may contain "Unresolved") and
    `uns["cell_type_prediction"]` (JSON text with settings and the per-cluster evidence).

    Raises `AnnotationError` if marker results are missing or no marker set can be scored.
    """
    groupby, use_raw = _marker_settings(adata)
    genes, missing, skipped = _present_genes(adata, params, use_raw)
    scored_sets = {ct: [g for g in gs if g not in missing.get(ct, ())]
                   for ct, gs in params.marker_sets.items() if ct not in skipped}
    if not scored_sets:
        raise AnnotationError("None of the marker genes are in the data; nothing can be scored.")

    # Same calculation as the validated notebook (section 13), kept in the same order.
    expr = sc.get.obs_df(adata, keys=genes, use_raw=use_raw)
    clusters = adata.obs[groupby]
    mean_expr = expr.groupby(clusters, observed=True).mean()
    frac = (expr > 0).groupby(clusters, observed=True).mean()
    z = _zscore(mean_expr)
    scores = pd.DataFrame({ct: z[gs].mean(axis=1) for ct, gs in scored_sets.items()})
    ranks, _ = _de_lookup(adata)
    sizes = clusters.value_counts()

    rows = []
    for cluster, cluster_scores in scores.iterrows():
        rows.append(_assess_cluster(
            str(cluster), cluster_scores, scored_sets, z.loc[cluster], frac.loc[cluster],
            ranks[str(cluster)], int(sizes[cluster]), params,
        ))
    table = pd.DataFrame(rows)

    adata = adata.copy()
    label_of = dict(zip(table["cluster"], table["predicted_label"]))
    adata.obs[PREDICTION_KEY] = pd.Categorical(adata.obs[groupby].astype(str).map(label_of))

    notes = []
    for cell_type, absent in missing.items():
        if cell_type in skipped:
            notes.append(f"{cell_type} not scored: none of its markers ({', '.join(absent)}) are in the data.")
        else:
            notes.append(f"{cell_type} scored without its missing markers ({', '.join(absent)}).")
    zero_variation = [g for g in z.columns if mean_expr[g].std(ddof=0) == 0]
    if zero_variation:
        notes.append(f"Markers with the same mean in every cluster (z-score set to 0): {', '.join(zero_variation)}.")

    summary = AnnotationSummary(
        groupby=groupby,
        clusters=table,
        scores=scores,
        label_counts={str(k): int(v) for k, v in adata.obs[PREDICTION_KEY].value_counts().items()},
        missing_genes={ct: tuple(gs) for ct, gs in missing.items()},
        skipped_cell_types=tuple(skipped),
        notes=tuple(notes),
    )
    adata.uns[RECORD_KEY] = json.dumps({
        "groupby": groupby,
        "settings": _plain_settings(params),
        "clusters": json.loads(table.to_json(orient="records")),
        "missing_genes": {ct: list(gs) for ct, gs in missing.items()},
        "notes": notes,
    })
    return adata, summary


def _assess_cluster(cluster, cluster_scores, scored_sets, z_row, frac_row, ranks, n_cells, params):
    """Best candidate, checks, and supporting/conflicting markers for one cluster."""
    ordered = cluster_scores.sort_values(ascending=False)  # as in the validated notebook
    best = ordered.index[0]
    second = ordered.index[1] if len(ordered) > 1 else None
    score = float(ordered.iloc[0])
    margin = float(ordered.iloc[0] - ordered.iloc[1]) if second is not None else None

    marker_ranks = {g: ranks.get(g) for g in scored_sets[best]}
    found = [r for r in marker_ranks.values() if r is not None]
    best_rank = min(found) if found else None

    reasons = []
    if margin is not None and margin < params.ambiguity_margin:
        reasons.append("ambiguous")
    if best_rank is None or best_rank > params.max_de_rank:
        reasons.append("unsupported")
    if params.min_score is not None and score < params.min_score:
        reasons.append("low_score")

    own_genes = set(scored_sets[best])
    conflicting = [
        f"{g} ({ct}, z={z_row[g]:.1f}, {frac_row[g]:.0%})"
        for ct, gs in scored_sets.items() if ct != best
        for g in gs if g not in own_genes and z_row[g] >= params.conflict_z
    ]
    return {
        "cluster": cluster,
        "n_cells": n_cells,
        "predicted_label": best if not reasons else UNRESOLVED,
        "status": reasons[0] if reasons else "confident",
        "reasons": ", ".join(reasons),
        "best_candidate": best,
        "score": score,
        "second_candidate": second,
        "margin": margin,
        "marker_DE_ranks": ", ".join(f"{g}:{r}" for g, r in marker_ranks.items()),
        "best_marker_DE_rank": best_rank,
        "supporting_markers": ", ".join(
            f"{g} ({frac_row[g]:.0%}, rank {marker_ranks[g]})" for g in scored_sets[best]
        ),
        "conflicting_markers": "; ".join(conflicting),
    }


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


def _present_genes(adata, params, use_raw):
    """Marker + context genes found in the data, plus what is missing per cell type."""
    available = set((adata.raw if use_raw else adata).var_names)
    missing = {ct: tuple(g for g in gs if g not in available)
               for ct, gs in params.marker_sets.items()}
    missing = {ct: gs for ct, gs in missing.items() if gs}
    skipped = [ct for ct, gs in params.marker_sets.items() if all(g not in available for g in gs)]
    genes = []
    for g in [*(g for gs in params.marker_sets.values() for g in gs), *params.context_markers]:
        if g in available and g not in genes:
            genes.append(g)
    return genes, missing, skipped


def _zscore(mean_expr: pd.DataFrame) -> pd.DataFrame:
    """Each gene's cluster means as z-scores across clusters (0 where a gene doesn't vary)."""
    z = (mean_expr - mean_expr.mean()) / mean_expr.std(ddof=0)
    return z.fillna(0.0)


def _de_lookup(adata: ad.AnnData) -> tuple[dict[str, dict[str, int]], dict[str, dict[str, float]]]:
    """Per cluster: gene -> rank (1 = top) and gene -> log fold change, from the marker step."""
    result = adata.uns[MARKER_RESULT_KEY]
    names = pd.DataFrame(result["names"])
    fold = pd.DataFrame(result["logfoldchanges"]) if "logfoldchanges" in result else None
    ranks, fold_changes = {}, {}
    for cluster in names.columns:
        genes = names[cluster].to_numpy()
        ranks[str(cluster)] = {g: i + 1 for i, g in enumerate(genes)}
        fold_changes[str(cluster)] = (
            {g: float(v) for g, v in zip(genes, fold[cluster].to_numpy())} if fold is not None else {}
        )
    return ranks, fold_changes


def _gene_owner(params: AnnotationParams) -> dict[str, list[str]]:
    owner: dict[str, list[str]] = {}
    for cell_type, genes in params.marker_sets.items():
        for g in genes:
            owner.setdefault(g, []).append(cell_type)
    return owner


def _plain_settings(params: AnnotationParams) -> dict:
    return {
        "marker_sets": {ct: list(gs) for ct, gs in params.marker_sets.items()},
        "context_markers": list(params.context_markers),
        "ambiguity_margin": params.ambiguity_margin,
        "max_de_rank": params.max_de_rank,
        "min_score": params.min_score,
        "conflict_z": params.conflict_z,
    }
