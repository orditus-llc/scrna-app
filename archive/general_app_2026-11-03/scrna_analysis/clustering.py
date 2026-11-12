"""Clustering: group cells with the Leiden algorithm on the neighbor graph.

The number of clusters is not chosen; it comes out of the data and the `resolution`
setting (higher resolution -> more, smaller clusters). Leiden numbers clusters by size
("0" is the largest). The numbers are labels only and carry no biological meaning.

Returns a new AnnData; the input is never changed.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import anndata as ad
import scanpy as sc

from scrna_analysis.params import SMALL_CLUSTER_CELLS, LeidenParams

# Scanpy warns on every leidenalg run that its *default* flavor will change to igraph.
# We pin the flavor explicitly, so the warning doesn't apply; it becomes a summary note.
_FUTURE_DEFAULT_WARNING = "In the future, the default backend for leiden will be igraph"
PINNED_FLAVOR_NOTE = (
    "Leiden flavor pinned to 'leidenalg'. Scanpy's future default ('igraph') is not used, "
    "because switching implementation can change the clusters."
)


class ClusteringError(ValueError):
    """Clustering can't run with this input or these settings. The message says why."""


@dataclass(frozen=True)
class ClusterSummary:
    n_clusters: int
    cluster_sizes: dict[str, int]  # cluster label -> number of cells, in label order
    smallest_cluster_size: int
    key_added: str
    flavor: str
    resolution: float
    n_iterations: int
    random_state: int
    directed: bool
    replaced_existing: bool  # True if the input already had a column named key_added
    notes: tuple[str, ...] = ()  # about the method (e.g. the pinned flavor)
    warnings: tuple[str, ...] = ()  # about these clusters, for the user (one cluster, tiny clusters)


def cluster_leiden(adata: ad.AnnData, params: LeidenParams) -> tuple[ad.AnnData, ClusterSummary]:
    """Leiden clustering on the neighbor graph. Adds `obs[params.key_added]` (categorical
    cluster labels "0", "1", ...) and `uns[params.key_added]` (Scanpy's record).

    Raises `ClusteringError` if the neighbor graph is missing, the resolution isn't
    positive, or the igraph flavor is combined with a directed graph.
    """
    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        raise ClusteringError("No neighbor graph found. Run compute_neighbors before clustering.")
    if not params.resolution > 0:
        raise ClusteringError(f"The Leiden resolution must be greater than 0 (got {params.resolution:g}).")
    if params.flavor == "igraph" and params.directed:
        raise ClusteringError(
            "The igraph Leiden flavor only works on an undirected graph; set directed=False."
        )

    replaced_existing = params.key_added in adata.obs
    adata = adata.copy()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=_FUTURE_DEFAULT_WARNING, category=FutureWarning)
        sc.tl.leiden(
            adata,
            flavor=params.flavor,
            resolution=params.resolution,
            n_iterations=params.n_iterations,
            random_state=params.random_state,
            directed=params.directed,
            key_added=params.key_added,
        )

    labels = adata.obs[params.key_added]
    sizes = labels.value_counts()
    ordered = sorted(sizes.index, key=_label_sort_key)
    cluster_sizes = {str(label): int(sizes[label]) for label in ordered}

    summary = ClusterSummary(
        n_clusters=len(cluster_sizes),
        cluster_sizes=cluster_sizes,
        smallest_cluster_size=min(cluster_sizes.values()),
        key_added=params.key_added,
        flavor=params.flavor,
        resolution=params.resolution,
        n_iterations=params.n_iterations,
        random_state=params.random_state,
        directed=params.directed,
        replaced_existing=replaced_existing,
        notes=(PINNED_FLAVOR_NOTE,) if params.flavor == "leidenalg" else (),
        warnings=cluster_warnings(cluster_sizes),
    )
    return adata, summary


def cluster_warnings(cluster_sizes: dict[str, int], small: int = SMALL_CLUSTER_CELLS) -> tuple[str, ...]:
    """Plain-language warnings about a clustering: only one cluster (nothing to compare, so
    marker genes can't run), or clusters with fewer than `small` cells."""
    warnings_ = []
    if len(cluster_sizes) == 1:
        warnings_.append("Only 1 cluster was found, so there is nothing to compare. "
                         "Raise the Leiden resolution to split the cells.")
    tiny = [label for label, n in cluster_sizes.items() if n < small]
    if tiny and len(cluster_sizes) > 1:
        which = "cluster" if len(tiny) == 1 else "clusters"
        warnings_.append(f"{len(tiny)} {which} ({', '.join(tiny)}) {'has' if len(tiny) == 1 else 'have'} fewer than "
                         f"{small} cells; {'its' if len(tiny) == 1 else 'their'} marker genes are less reliable. "
                         "A lower resolution merges small clusters.")
    return tuple(warnings_)


def _label_sort_key(label: str) -> tuple[int, int | str]:
    """Numeric labels in numeric order ("2" before "10"); any others after, alphabetically."""
    text = str(label)
    return (0, int(text)) if text.isdigit() else (1, text)
