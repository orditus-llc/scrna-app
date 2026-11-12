"""QC plots for the app (matplotlib). Each function takes data and returns a figure.

Same plots as the reference notebook (violins of the three QC metrics; counts vs %
mitochondrial and counts vs genes), plus the current filter thresholds as lines and the
cells the filter would remove in grey. No Shiny code and no calculations.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # draw to images (no desktop windows); the app shows the images

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.figure import Figure

from pbmc_analysis.params import FilterParams

KEPT_COLOR = "#2a7f8a"
REMOVED_COLOR = "#c3cbcf"
THRESHOLD_COLOR = "#c0392b"
LABELS = {
    "n_genes_by_counts": "Genes per cell",
    "total_counts": "Counts per cell",
    "pct_counts_mt": "% mitochondrial counts",
}


def qc_violins(obs: pd.DataFrame, mito_available: bool = True, seed: int = 0) -> Figure:
    """Violin + jittered points for genes per cell, counts per cell and % mitochondrial."""
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    rng = np.random.default_rng(seed)  # jitter only; doesn't affect any value
    for ax, metric in zip(axes, metrics):
        ax.set_title(LABELS[metric], fontsize=11)
        ax.set_xticks([])
        if metric == "pct_counts_mt" and not mito_available:
            ax.text(0.5, 0.5, "Unavailable:\nno mitochondrial genes found",
                    ha="center", va="center", transform=ax.transAxes, color="#555555")
            continue
        values = obs[metric].to_numpy(dtype=float)
        parts = ax.violinplot(values, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(KEPT_COLOR)
            body.set_alpha(0.35)
        ax.scatter(1 + rng.uniform(-0.2, 0.2, len(values)), values, s=2, color=KEPT_COLOR, alpha=0.5)
    fig.tight_layout()
    return fig


def qc_scatters(
    obs: pd.DataFrame,
    kept_cells: pd.Index | None = None,
    thresholds: FilterParams | None = None,
) -> Figure:
    """Counts vs % mitochondrial, and counts vs genes per cell.

    Cells not in `kept_cells` are drawn in grey (removed by the current filter settings);
    the active thresholds from `thresholds` are drawn as dashed lines.
    """
    kept = np.ones(len(obs), dtype=bool) if kept_cells is None else obs.index.isin(kept_cells)
    panels = [("pct_counts_mt", "max_pct_counts_mt", None),
              ("n_genes_by_counts", "max_n_genes_by_counts", "min_genes")]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, (y_metric, upper_name, lower_name) in zip(axes, panels):
        x, y = obs["total_counts"].to_numpy(float), obs[y_metric].to_numpy(float)
        ax.scatter(x[~kept], y[~kept], s=4, color=REMOVED_COLOR, label=f"removed ({(~kept).sum()})")
        ax.scatter(x[kept], y[kept], s=4, color=KEPT_COLOR, label=f"kept ({kept.sum()})")
        if thresholds is not None:
            for name, style in ((upper_name, "--"), (lower_name, ":")):
                limit = getattr(thresholds, name) if name else None
                if limit is not None:
                    ax.axhline(limit, color=THRESHOLD_COLOR, linestyle=style, linewidth=1.2)
        ax.set_xlabel(LABELS["total_counts"])
        ax.set_ylabel(LABELS[y_metric])
        ax.legend(loc="upper right", fontsize=8, markerscale=3, frameon=False)
    fig.tight_layout()
    return fig


# --- Later stages: Scanpy's own plots ------------------------------------------------------
#
# Some Scanpy plotting functions quietly store things in the data they are given (e.g. sc.pl.umap
# saves cluster colors; the marker dot plot computes a dendrogram). These functions therefore
# always plot from a copy, so the pipeline's checkpoints are never changed by drawing.

def hvg_plot(adata: ad.AnnData) -> Figure:
    """Mean expression vs dispersion, highly variable genes highlighted (all genes)."""
    plt.close("all")
    genes = adata.raw.var if adata.raw is not None else adata.var  # X keeps only the HVGs
    columns = ["means", "dispersions", "dispersions_norm", "highly_variable"]
    sc.pl.highly_variable_genes(genes[columns].copy(), show=False)
    fig = plt.gcf()
    fig.set_size_inches(10, 4)
    return fig


def pca_variance_plot(adata: ad.AnnData) -> Figure:
    """Variance explained by each principal component (the "elbow" plot), log scale."""
    plt.close("all")
    n_pcs = adata.obsm["X_pca"].shape[1]
    sc.pl.pca_variance_ratio(adata.copy(), n_pcs=n_pcs, log=True, show=False)
    fig = plt.gcf()
    fig.set_size_inches(6, 4)
    return fig


def umap_plot(adata: ad.AnnData, color: str | None = None, title: str | None = None,
              labels_on_plot: bool = False) -> Figure:
    """UMAP, optionally colored by an obs column (e.g. clusters or predicted cell types)."""
    plt.close("all")
    fig = sc.pl.umap(
        adata.copy(), color=color, title=title, show=False, return_fig=True,
        legend_loc="on data" if labels_on_plot else "right margin", frameon=False,
    )
    fig.set_size_inches(7, 5.5)
    return fig


def marker_dotplot(adata: ad.AnnData, n_genes: int = 3) -> Figure:
    """Dot plot of the top `n_genes` marker genes of every cluster (size = % of cells
    expressing, color = mean expression). No dendrogram: clusters keep their own order."""
    plt.close("all")
    dotplot = sc.pl.rank_genes_groups_dotplot(
        adata.copy(), n_genes=n_genes, dendrogram=False, show=False, return_fig=True,
    )
    dotplot.make_figure()
    return dotplot.fig
