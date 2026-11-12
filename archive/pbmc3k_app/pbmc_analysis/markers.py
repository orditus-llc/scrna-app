"""Marker genes: which genes distinguish each cluster from all other cells.

`rank_marker_genes` tests every gene in every cluster against the rest of the cells
(Wilcoxon rank-sum test in the reference) and ranks the genes, using `adata.raw`
(log-normalized, all genes) so markers aren't limited to the highly variable genes.
`marker_table` turns the result into one tidy table for display and annotation.
`expression_by_group` gives, for chosen genes, the mean expression and fraction of cells
expressing in each group (what a marker dot plot shows).

Statistical caveat: the clusters were defined from the same cells, so p-values are
strongly inflated. Use them to rank genes, not as evidence of significance.
"""

from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import pandas as pd
import scanpy as sc

from pbmc_analysis.params import MarkerParams

RESULT_KEY = "rank_genes_groups"  # where Scanpy stores the results in adata.uns
N_TOP_IN_SUMMARY = 25
DOUBLE_DIPPING_NOTE = (
    "Clusters were defined from the same cells that were tested, so p-values are strongly "
    "inflated. Use the ranking to find candidate markers, not as evidence of significance."
)

# Scanpy's column names -> clearer names for the tidy table.
_COLUMN_NAMES = {
    "names": "gene",
    "scores": "score",
    "logfoldchanges": "logfoldchange",
    "pvals": "pval",
    "pvals_adj": "pval_adj",
    "pct_nz_group": "frac_in_group",  # fraction (0-1) of the cluster's cells expressing it
    "pct_nz_reference": "frac_in_rest",  # same, among all other cells
}


class MarkerError(ValueError):
    """Marker genes can't be computed or shown as requested. The message says why."""


@dataclass(frozen=True)
class MarkerSummary:
    groupby: str
    groups: tuple[str, ...]  # the clusters compared, in their category order
    n_genes_tested: int
    method: str
    use_raw: bool
    reference: str
    corr_method: str
    tie_correct: bool
    top_markers: dict[str, tuple[str, ...]]  # cluster -> its top-ranked genes
    notes: tuple[str, ...] = ()


def rank_marker_genes(adata: ad.AnnData, params: MarkerParams) -> tuple[ad.AnnData, MarkerSummary]:
    """Rank all genes for each group in `obs[params.groupby]` against the rest.

    Adds Scanpy's results to `uns["rank_genes_groups"]`: for every group, all genes in
    rank order with scores, log fold changes, p-values, adjusted p-values and (with
    `pts=True`) the fraction of cells expressing each gene inside and outside the group.

    Raises `MarkerError` if the grouping column is missing, there are fewer than two
    groups, a group has fewer than two cells, or `use_raw=True` without `adata.raw`.
    """
    _check_inputs(adata, params)

    adata = adata.copy()
    # Scanpy needs the grouping as a categorical column; Leiden clusters already are.
    if not isinstance(adata.obs[params.groupby].dtype, pd.CategoricalDtype):
        adata.obs[params.groupby] = adata.obs[params.groupby].astype("category")

    sc.tl.rank_genes_groups(
        adata,
        params.groupby,
        method=params.method,
        use_raw=params.use_raw,
        pts=params.pts,
        reference=params.reference,
        corr_method=params.corr_method,
        tie_correct=params.tie_correct,
    )

    groups = tuple(str(g) for g in adata.obs[params.groupby].cat.categories)
    ranked_names = pd.DataFrame(adata.uns[RESULT_KEY]["names"])
    summary = MarkerSummary(
        groupby=params.groupby,
        groups=groups,
        n_genes_tested=len(ranked_names),
        method=params.method,
        use_raw=params.use_raw,
        reference=params.reference,
        corr_method=params.corr_method,
        tie_correct=params.tie_correct,
        top_markers={g: tuple(ranked_names[g].head(N_TOP_IN_SUMMARY)) for g in ranked_names.columns},
        notes=(DOUBLE_DIPPING_NOTE,),
    )
    return adata, summary


def marker_table(
    adata: ad.AnnData, group: str | None = None, n_top: int | None = None
) -> pd.DataFrame:
    """Marker results as one tidy table (doesn't change `adata`).

    One row per gene per group, in rank order, with columns: group, rank (1 = most
    distinctive), gene, score, logfoldchange, pval, pval_adj, and, if computed,
    frac_in_group and frac_in_rest (fractions from 0 to 1).

    Parameters
    ----------
    group
        Only this group; `None` for all groups.
    n_top
        Only the top `n_top` genes per group; `None` for all genes.
    """
    if RESULT_KEY not in adata.uns:
        raise MarkerError("No marker results found. Run rank_marker_genes first.")
    available = list(pd.DataFrame(adata.uns[RESULT_KEY]["names"]).columns)
    if group is not None and group not in available:
        raise MarkerError(f"Unknown group '{group}'. Available groups: {available}.")
    if n_top is not None and n_top < 1:
        raise MarkerError(f"n_top must be at least 1, got {n_top}.")

    selected = available if group is None else [group]
    table = sc.get.rank_genes_groups_df(adata, group=selected)
    if "group" not in table.columns:  # Scanpy leaves it out for a single group
        table.insert(0, "group", group)
    table["group"] = table["group"].astype(str)
    table = table.rename(columns=_COLUMN_NAMES)
    table.insert(1, "rank", table.groupby("group", sort=False).cumcount() + 1)

    if n_top is not None:
        table = table[table["rank"] <= n_top]
    return table.reset_index(drop=True)


def expression_by_group(
    adata: ad.AnnData, genes: list[str], groupby: str, use_raw: bool
) -> pd.DataFrame:
    """Mean expression and fraction of cells expressing, per group, for `genes`.

    One row per group per gene, groups in their category order, with columns: group, gene,
    mean_expr (mean of the values in `adata.raw` if `use_raw`, else `adata.X`) and
    frac_expressing (fraction of the group's cells with a value above 0). Doesn't change
    `adata`. Raises `MarkerError` if the grouping or a gene is missing.
    """
    if groupby not in adata.obs:
        raise MarkerError(f"obs['{groupby}'] not found.")
    if use_raw and adata.raw is None:
        raise MarkerError("use_raw=True, but adata.raw is empty.")
    available = set((adata.raw if use_raw else adata).var_names)
    missing = [g for g in genes if g not in available]
    if missing:
        raise MarkerError(f"Genes not found in the data: {missing}.")
    if not genes:
        return pd.DataFrame(columns=["group", "gene", "mean_expr", "frac_expressing"])

    expr = sc.get.obs_df(adata, keys=list(dict.fromkeys(genes)), use_raw=use_raw)
    groups = adata.obs[groupby]
    mean_expr = expr.groupby(groups, observed=True).mean()
    frac = (expr > 0).groupby(groups, observed=True).mean()
    rows = [
        {"group": str(group), "gene": gene,
         "mean_expr": float(mean_expr.loc[group, gene]), "frac_expressing": float(frac.loc[group, gene])}
        for group in mean_expr.index for gene in expr.columns
    ]
    return pd.DataFrame(rows)


def _check_inputs(adata: ad.AnnData, params: MarkerParams) -> None:
    if params.groupby not in adata.obs:
        raise MarkerError(
            f"obs['{params.groupby}'] not found. Run clustering first, or set groupby to "
            "an existing column."
        )
    sizes = adata.obs[params.groupby].astype(str).value_counts()
    if len(sizes) < 2:
        raise MarkerError(
            f"Marker genes need at least 2 groups to compare; obs['{params.groupby}'] has {len(sizes)}."
        )
    too_small = sorted(sizes[sizes < 2].index)
    if too_small:
        raise MarkerError(
            f"{'Group' if len(too_small) == 1 else 'Groups'} {', '.join(map(str, too_small))} "
            f"{'has' if len(too_small) == 1 else 'have'} fewer than 2 cells; the test needs at least 2 per group."
        )
    if params.reference != "rest" and params.reference not in sizes.index:
        raise MarkerError(
            f"reference='{params.reference}' is not a group in obs['{params.groupby}']."
        )
    if params.use_raw and adata.raw is None:
        raise MarkerError(
            "use_raw=True but adata.raw is empty. Run keep_raw_and_subset first, or set use_raw=False."
        )
