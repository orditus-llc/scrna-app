"""Dimensionality reduction: PCA, the nearest-neighbor graph, and UMAP.

Order (same as the validated reference analysis), applied to scaled data:

1. `run_pca`: summarize the genes as a smaller number of principal components (PCs).
2. `compute_neighbors`: connect each cell to its most similar cells in PC space.
   Clustering and UMAP are both built on this graph.
3. `run_umap`: a 2-D layout of the neighbor graph, for visualization only.

UMAP uses only the neighbor graph (spectral initialization), not the clusters. However,
the UMAP library prunes the weakest edges of the graph it is given *in place*, which
would silently change the graph that Leiden clusters (found 2026-10-13: about 21 PBMC3k
cells moved between T/NK clusters when UMAP ran first). `run_umap` therefore restores
the true neighbor graph afterwards, so clustering gives the same result whether it runs
before or after UMAP.

Every function returns a new AnnData; the input is never changed.
"""

from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse

from pbmc_analysis.params import NeighborsParams, PCAParams, UMAPParams


class DimReductionError(ValueError):
    """A dimensionality-reduction step can't run on this input. The message says why."""


@dataclass(frozen=True)
class PCASummary:
    n_comps: int
    svd_solver: str
    variance_ratio: tuple[float, ...]  # fraction of variance explained by each PC
    total_variance_ratio: float  # sum over all computed PCs


@dataclass(frozen=True)
class NeighborsSummary:
    n_neighbors: int
    n_pcs: int
    use_rep: str
    method: str
    metric: str
    n_connections: int  # non-zero entries in the connectivity matrix


@dataclass(frozen=True)
class UMAPSummary:
    n_cells: int
    n_dims: int


def run_pca(adata: ad.AnnData, params: PCAParams) -> tuple[ad.AnnData, PCASummary]:
    """Principal component analysis on `adata.X`. Adds `obsm["X_pca"]`, `varm["PCs"]`
    and `uns["pca"]` (variance and variance ratio per PC).

    Raises `DimReductionError` if `X` contains NaN/infinite values or more components
    are requested than the data supports.
    """
    _check_finite(adata.X)
    limit = min(adata.n_obs, adata.n_vars)
    if not 0 < params.n_comps < limit:
        raise DimReductionError(
            f"The number of principal components ({params.n_comps:,}) must be between 1 and "
            f"{limit - 1:,} (fewer than the {adata.n_obs:,} cells and {adata.n_vars:,} genes)."
        )

    adata = adata.copy()
    sc.pp.pca(
        adata,
        n_comps=params.n_comps,
        svd_solver=params.svd_solver,
        random_state=params.random_state,
    )
    ratios = tuple(float(r) for r in adata.uns["pca"]["variance_ratio"])
    summary = PCASummary(
        n_comps=params.n_comps,
        svd_solver=params.svd_solver,
        variance_ratio=ratios,
        total_variance_ratio=float(sum(ratios)),
    )
    return adata, summary


def compute_neighbors(
    adata: ad.AnnData, params: NeighborsParams
) -> tuple[ad.AnnData, NeighborsSummary]:
    """Nearest-neighbor graph on the first `n_pcs` columns of `obsm[use_rep]`.

    Adds `obsp["connectivities"]`, `obsp["distances"]` and `uns["neighbors"]`.
    Raises `DimReductionError` if the representation is missing, has too few
    components, or there are too few cells for the requested number of neighbors.
    """
    if params.use_rep not in adata.obsm:
        raise DimReductionError(
            f"obsm['{params.use_rep}'] not found. Run run_pca before compute_neighbors."
        )
    n_available = adata.obsm[params.use_rep].shape[1]
    if params.n_pcs > n_available:
        raise DimReductionError(
            f"n_pcs={params.n_pcs} but obsm['{params.use_rep}'] has only {n_available} "
            "components. Lower n_pcs or re-run PCA with more components."
        )
    if not 1 < params.n_neighbors <= adata.n_obs:
        raise DimReductionError(
            f"n_neighbors={params.n_neighbors} must be between 2 and the number of cells "
            f"({adata.n_obs})."
        )

    adata = adata.copy()
    sc.pp.neighbors(
        adata,
        n_neighbors=params.n_neighbors,
        n_pcs=params.n_pcs,
        use_rep=params.use_rep,
        random_state=params.random_state,
    )
    recorded = adata.uns["neighbors"]["params"]
    summary = NeighborsSummary(
        n_neighbors=params.n_neighbors,
        n_pcs=params.n_pcs,
        use_rep=params.use_rep,
        method=str(recorded["method"]),
        metric=str(recorded["metric"]),
        n_connections=int(adata.obsp["connectivities"].nnz),
    )
    return adata, summary


def run_umap(adata: ad.AnnData, params: UMAPParams) -> tuple[ad.AnnData, UMAPSummary]:
    """2-D UMAP layout of the neighbor graph. Adds `obsm["X_umap"]`.

    For visualization only: distances on a UMAP aren't reliable measures of similarity.
    The stored neighbor graph is left exactly as `compute_neighbors` made it (see the
    module docstring). Raises `DimReductionError` if the neighbor graph hasn't been computed.
    """
    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        raise DimReductionError("No neighbor graph found. Run compute_neighbors before run_umap.")

    adata = adata.copy()
    true_graph = adata.obsp["connectivities"].copy()
    sc.tl.umap(adata, random_state=params.random_state)
    # UMAP pruned weak edges of the stored graph as a side effect; put the true graph back.
    adata.obsp["connectivities"] = true_graph
    n_cells, n_dims = adata.obsm["X_umap"].shape
    return adata, UMAPSummary(n_cells=n_cells, n_dims=n_dims)


def _check_finite(matrix) -> None:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all():
        raise DimReductionError("X contains NaN or infinite values; PCA needs finite data.")
