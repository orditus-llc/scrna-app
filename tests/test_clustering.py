"""Tests for `scrna_analysis.clustering`.

- Reference tests: Leiden on the validated PBMC3k neighbor graph with the reference
  preset reproduces the validated notebook's clusters (read from `validation/expected/`).
- Generality tests: a PBMC3k subsample at other resolutions compared against direct
  Scanpy calls, the igraph flavor, warning handling, and edge cases.
"""

from __future__ import annotations

import warnings
from dataclasses import replace

import anndata as ad
import numpy as np
import pytest
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from scrna_analysis.clustering import PINNED_FLAVOR_NOTE, ClusteringError, cluster_leiden, cluster_warnings
from scrna_analysis.dim_reduction import compute_neighbors, run_pca
from scrna_analysis.params import LeidenParams, NeighborsParams, PCAParams
from tests.reference_preset import REFERENCE_PRESET

# regress_out (used to build the shared fixtures) triggers harmless Numba performance hints.
pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")

P = REFERENCE_PRESET


def labels(adata: ad.AnnData, key: str = "leiden") -> np.ndarray:
    return adata.obs[key].astype(str).to_numpy()


def direct_leiden(adata: ad.AnnData, **kwargs) -> ad.AnnData:
    """Call Scanpy's Leiden directly (silencing its known future-default warning)."""
    direct = adata.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        sc.tl.leiden(direct, **kwargs)
    return direct


# --- Reference tests (PBMC3k) ----------------------------------------------------------

@pytest.fixture(scope="module")
def clustered(graph_pbmc3k):
    return cluster_leiden(graph_pbmc3k, P.leiden)


def test_every_cell_has_reference_cluster(clustered, reference_results):
    adata, _ = clustered
    assert adata.obs_names.equals(reference_results.obs_names)
    np.testing.assert_array_equal(labels(adata), labels(reference_results))


def test_adjusted_rand_index_is_one(clustered, reference_results):
    adata, _ = clustered
    assert adjusted_rand_score(labels(reference_results), labels(adata)) == 1.0


def test_cluster_count_and_sizes_match_reference(clustered, reference_summary):
    _, summary = clustered
    assert summary.n_clusters == reference_summary["n_clusters"]
    assert summary.cluster_sizes == reference_summary["cluster_sizes"]
    assert list(summary.cluster_sizes) == list(reference_summary["cluster_sizes"])  # same order
    assert summary.smallest_cluster_size == min(reference_summary["cluster_sizes"].values())


def test_cluster_categories_match_reference(clustered, reference_results):
    adata, _ = clustered
    assert list(adata.obs["leiden"].cat.categories) == list(reference_results.obs["leiden"].cat.categories)


def test_clusters_unchanged_when_umap_runs_first(graph_pbmc3k, reference_results):
    """Guards the 2026-10-13 finding: UMAP must not alter the graph Leiden clusters."""
    from scrna_analysis.dim_reduction import run_umap

    after_umap, _ = run_umap(graph_pbmc3k, P.umap)
    clustered_after_umap, _ = cluster_leiden(after_umap, P.leiden)
    np.testing.assert_array_equal(labels(clustered_after_umap), labels(reference_results))


def test_scanpy_record_matches_reference(clustered, reference_results):
    adata, _ = clustered
    assert dict(adata.uns["leiden"]["params"]) == dict(reference_results.uns["leiden"]["params"])


def test_summary_records_all_settings(clustered):
    _, summary = clustered
    assert (summary.flavor, summary.resolution, summary.n_iterations, summary.random_state,
            summary.directed, summary.key_added) == ("leidenalg", 1.0, -1, 0, True, "leiden")
    assert summary.notes == (PINNED_FLAVOR_NOTE,)
    assert not summary.replaced_existing


# --- Generality tests: subsample with other settings ------------------------------------

@pytest.fixture(scope="module")
def small_graph(scaled_pbmc3k):
    """600 random cells with their own PCA and neighbor graph: a different input."""
    rng = np.random.default_rng(5)
    small = scaled_pbmc3k[np.sort(rng.choice(scaled_pbmc3k.n_obs, 600, replace=False))].copy()
    small, _ = run_pca(small, PCAParams(n_comps=20))
    small, _ = compute_neighbors(small, NeighborsParams(n_neighbors=15, n_pcs=20))
    return small


@pytest.mark.parametrize("resolution", [0.5, 2.0])
def test_other_resolutions_match_direct_scanpy(small_graph, resolution):
    params = LeidenParams(resolution=resolution, random_state=4)
    ours, summary = cluster_leiden(small_graph, params)
    direct = direct_leiden(small_graph, flavor="leidenalg", resolution=resolution,
                           n_iterations=-1, random_state=4, directed=True)
    np.testing.assert_array_equal(labels(ours), labels(direct))
    # The number of clusters is whatever the data gives, never assumed.
    assert summary.n_clusters == direct.obs["leiden"].nunique()
    assert sum(summary.cluster_sizes.values()) == small_graph.n_obs


def test_higher_resolution_gives_more_clusters(small_graph):
    _, coarse = cluster_leiden(small_graph, LeidenParams(resolution=0.3))
    _, fine = cluster_leiden(small_graph, LeidenParams(resolution=2.0))
    assert fine.n_clusters > coarse.n_clusters


def test_same_seed_gives_identical_clusters(small_graph):
    first, _ = cluster_leiden(small_graph, LeidenParams(random_state=9))
    again, _ = cluster_leiden(small_graph, LeidenParams(random_state=9))
    np.testing.assert_array_equal(labels(first), labels(again))


def test_igraph_flavor_matches_direct_scanpy(small_graph):
    params = LeidenParams(flavor="igraph", directed=False, n_iterations=2, random_state=0)
    ours, summary = cluster_leiden(small_graph, params)
    direct = direct_leiden(small_graph, flavor="igraph", directed=False, n_iterations=2,
                           random_state=0)
    np.testing.assert_array_equal(labels(ours), labels(direct))
    assert summary.flavor == "igraph" and summary.notes == ()


def test_custom_key_and_existing_clusters_are_kept_separate(small_graph):
    first, _ = cluster_leiden(small_graph, LeidenParams(resolution=0.5))
    second, summary = cluster_leiden(first, LeidenParams(resolution=2.0, key_added="leiden_r2"))
    np.testing.assert_array_equal(labels(second, "leiden"), labels(first, "leiden"))
    assert "leiden_r2" in second.obs and not summary.replaced_existing

    labels_before = labels(first).copy()
    _, rerun = cluster_leiden(first, LeidenParams(resolution=2.0))
    assert rerun.replaced_existing  # reusing the same key is recorded in the summary...
    np.testing.assert_array_equal(labels(first), labels_before)  # ...and the input keeps its labels


# --- Warning handling ------------------------------------------------------------------

def test_future_default_warning_becomes_a_note(small_graph, monkeypatch):
    """The irrelevant 'default will change' warning is silenced and noted in the summary."""
    original = sc.tl.leiden

    def leiden_that_warns(*args, **kwargs):  # emit the warning every time, to test reliably
        warnings.warn("In the future, the default backend for leiden will be igraph instead "
                      "of leidenalg.", FutureWarning)
        return original(*args, **kwargs)

    monkeypatch.setattr(sc.tl, "leiden", leiden_that_warns)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _, summary = cluster_leiden(small_graph, LeidenParams())
    assert not [w for w in caught if "default backend for leiden" in str(w.message)]
    assert summary.notes == (PINNED_FLAVOR_NOTE,)


def test_other_warnings_still_come_through(small_graph, monkeypatch):
    original = sc.tl.leiden

    def leiden_with_other_warning(*args, **kwargs):
        warnings.warn("something else worth knowing", UserWarning)
        return original(*args, **kwargs)

    monkeypatch.setattr(sc.tl, "leiden", leiden_with_other_warning)
    with pytest.warns(UserWarning, match="something else worth knowing"):
        cluster_leiden(small_graph, LeidenParams())


# --- Edge cases ------------------------------------------------------------------------

def test_clustering_before_neighbors_raises():
    adata = ad.AnnData(np.random.default_rng(0).normal(size=(20, 5)).astype(np.float32))
    with pytest.raises(ClusteringError, match="Run compute_neighbors before clustering"):
        cluster_leiden(adata, LeidenParams())


@pytest.mark.parametrize("resolution", [0, -1.0])
def test_non_positive_resolution_raises(small_graph, resolution):
    with pytest.raises(ClusteringError, match="resolution must be greater than 0"):
        cluster_leiden(small_graph, LeidenParams(resolution=resolution))


def test_igraph_with_directed_graph_raises(small_graph):
    with pytest.raises(ClusteringError, match="set directed=False"):
        cluster_leiden(small_graph, LeidenParams(flavor="igraph", directed=True))


def test_input_is_not_modified(small_graph):
    obs_before = small_graph.obs.copy()
    uns_before = set(small_graph.uns)
    cluster_leiden(small_graph, replace(P.leiden, key_added="check_input"))
    assert "check_input" not in small_graph.obs
    assert set(small_graph.uns) == uns_before
    assert small_graph.obs.equals(obs_before)


# --- Warnings about the clusters (shown above the results) --------------------------------------

def test_warnings_for_one_cluster_and_tiny_clusters():
    one = cluster_warnings({"0": 500})
    assert len(one) == 1 and one[0].startswith("Only 1 cluster was found") and "Raise the Leiden resolution" in one[0]
    tiny = cluster_warnings({"0": 400, "1": 9, "2": 3, "3": 10})
    assert len(tiny) == 1 and tiny[0].startswith("2 clusters (1, 2) have fewer than 10 cells")
    assert cluster_warnings({"0": 400, "1": 4})[0].startswith("1 cluster (1) has fewer than 10 cells; its marker genes")
    assert cluster_warnings({"0": 400, "1": 10}) == ()


def test_very_low_resolution_gives_one_cluster_with_a_warning(graph_pbmc3k):
    adata, summary = cluster_leiden(graph_pbmc3k, LeidenParams(resolution=0.001))
    assert summary.n_clusters == 1 and summary.warnings[0].startswith("Only 1 cluster")
