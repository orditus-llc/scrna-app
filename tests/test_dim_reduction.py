"""Tests for `scrna_analysis.dim_reduction`.

- Reference tests: PCA, neighbors and UMAP on preprocessed PBMC3k with the reference
  preset reproduce the validated notebook (values read from `validation/expected/`).
- Generality tests: a PBMC3k subsample with different settings compared against direct
  Scanpy calls, seed behavior, and edge cases.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
import scanpy as sc

from scrna_analysis.dim_reduction import (
    find_elbow,
    DimReductionError,
    compute_neighbors,
    run_pca,
    run_umap,
)
from scrna_analysis.params import NeighborsParams, PCAParams, UMAPParams
from tests.reference_preset import REFERENCE_PRESET

# regress_out (used to build the subsample) triggers harmless Numba performance hints.
pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")

P = REFERENCE_PRESET
TOL = 1e-6  # strict tolerance agreed for floats (CLAUDE.md, architecture decision 2)


def max_diff_allowing_sign_flips(ours: np.ndarray, expected: np.ndarray) -> float:
    """Largest difference after flipping any PC whose sign is opposite to the reference.

    A PC and its negative describe the same direction, so a flip is not a difference.
    """
    signs = np.sign(np.sum(ours * expected, axis=0))
    signs[signs == 0] = 1
    return float(np.max(np.abs(ours * signs - expected)))


def sparse_max_diff(a, b) -> float:
    diff = abs(a - b)
    return float(diff.max()) if diff.nnz else 0.0


# --- Reference tests (PBMC3k) ----------------------------------------------------------

@pytest.fixture(scope="module")
def reduced(scaled_pbmc3k):
    with_pca, pca_summary = run_pca(scaled_pbmc3k, P.pca)
    with_graph, neighbors_summary = compute_neighbors(with_pca, P.neighbors)
    with_umap, umap_summary = run_umap(with_graph, P.umap)
    return with_umap, pca_summary, neighbors_summary, umap_summary


def test_pca_coordinates_match_reference(reduced, reference_results):
    adata = reduced[0]
    assert adata.obsm["X_pca"].shape == reference_results.obsm["X_pca"].shape
    assert max_diff_allowing_sign_flips(adata.obsm["X_pca"], reference_results.obsm["X_pca"]) <= TOL


def test_pca_loadings_match_reference(reduced, reference_results):
    adata = reduced[0]
    assert max_diff_allowing_sign_flips(adata.varm["PCs"], reference_results.varm["PCs"]) <= TOL


def test_pca_variance_ratio_matches_reference(reduced, reference_summary):
    summary = reduced[1]
    np.testing.assert_allclose(
        summary.variance_ratio, reference_summary["pca_variance_ratio"], rtol=0, atol=TOL
    )
    assert summary.n_comps == len(reference_summary["pca_variance_ratio"])


def test_neighbor_distances_match_reference(reduced, reference_results):
    ours, expected = reduced[0].obsp["distances"], reference_results.obsp["distances"]
    assert ours.shape == expected.shape and ours.nnz == expected.nnz
    assert ((ours != 0) != (expected != 0)).nnz == 0  # same neighbors
    assert sparse_max_diff(ours, expected) <= TOL


def test_connectivities_correspond_to_reference(reduced, reference_results):
    """The notebook saved its graph *after* UMAP had pruned it in place (a UMAP side effect).
    We keep the true graph; applying UMAP's pruning to it must give the notebook's graph."""
    pruned = reduced[0].copy()
    sc.tl.umap(pruned, random_state=P.umap.random_state)  # direct call: prunes the stored graph
    expected = reference_results.obsp["connectivities"]
    assert pruned.obsp["connectivities"].nnz == expected.nnz
    assert sparse_max_diff(pruned.obsp["connectivities"], expected) <= TOL


def test_run_umap_leaves_the_neighbor_graph_unchanged(scaled_pbmc3k):
    with_pca, _ = run_pca(scaled_pbmc3k, P.pca)
    with_graph, _ = compute_neighbors(with_pca, P.neighbors)
    with_umap, _ = run_umap(with_graph, P.umap)
    for matrix in ("connectivities", "distances"):
        assert sparse_max_diff(with_umap.obsp[matrix], with_graph.obsp[matrix]) == 0, matrix


def test_neighbor_settings_match_reference(reduced, reference_results):
    ours = reduced[0].uns["neighbors"]["params"]
    expected = reference_results.uns["neighbors"]["params"]
    for key in ("n_neighbors", "method", "metric", "random_state", "n_pcs", "use_rep"):
        assert ours[key] == expected[key], key


def test_umap_matches_reference(reduced, reference_results):
    ours, expected = reduced[0].obsm["X_umap"], reference_results.obsm["X_umap"]
    assert ours.shape == expected.shape
    assert float(np.max(np.abs(ours - expected))) <= TOL


def test_reference_summaries(reduced, reference_results):
    _, pca_summary, neighbors_summary, umap_summary = reduced
    assert pca_summary.svd_solver == "arpack"
    assert neighbors_summary.n_connections == reference_results.obsp["connectivities"].nnz
    assert (umap_summary.n_cells, umap_summary.n_dims) == reference_results.obsm["X_umap"].shape


# --- Generality tests ------------------------------------------------------------------

@pytest.fixture(scope="module")
def small_scaled(scaled_pbmc3k):
    """600 random cells of the scaled data: a different input from the reference."""
    rng = np.random.default_rng(11)
    return scaled_pbmc3k[np.sort(rng.choice(scaled_pbmc3k.n_obs, 600, replace=False))].copy()


def test_other_settings_match_direct_scanpy(small_scaled):
    pca = PCAParams(n_comps=20, random_state=3)
    neighbors = NeighborsParams(n_neighbors=15, n_pcs=15, random_state=3)
    umap = UMAPParams(random_state=3)

    a, _ = run_pca(small_scaled, pca)
    b, _ = compute_neighbors(a, neighbors)
    c, _ = run_umap(b, umap)

    direct = small_scaled.copy()
    sc.pp.pca(direct, n_comps=20, svd_solver="arpack", random_state=3)
    sc.pp.neighbors(direct, n_neighbors=15, n_pcs=15, use_rep="X_pca", random_state=3)
    sc.tl.umap(direct, random_state=3)

    np.testing.assert_array_equal(c.obsm["X_pca"], direct.obsm["X_pca"])
    # Compare the graph before direct UMAP pruned it (ours keeps the true graph).
    assert sparse_max_diff(b.obsp["connectivities"], c.obsp["connectivities"]) == 0
    np.testing.assert_array_equal(c.obsm["X_umap"], direct.obsm["X_umap"])


def test_variance_ratios_are_sensible(small_scaled):
    _, summary = run_pca(small_scaled, PCAParams(n_comps=20))
    ratios = np.array(summary.variance_ratio)
    assert np.all(np.diff(ratios) <= 1e-12)  # each PC explains no more than the previous
    assert 0 < summary.total_variance_ratio <= 1


def test_umap_seed_is_used(small_scaled):
    graph, _ = compute_neighbors(run_pca(small_scaled, PCAParams(n_comps=20))[0],
                                 NeighborsParams(n_neighbors=10, n_pcs=20))
    first, _ = run_umap(graph, UMAPParams(random_state=0))
    again, _ = run_umap(graph, UMAPParams(random_state=0))
    other, _ = run_umap(graph, UMAPParams(random_state=1))
    np.testing.assert_array_equal(first.obsm["X_umap"], again.obsm["X_umap"])
    assert not np.array_equal(first.obsm["X_umap"], other.obsm["X_umap"])


# --- Edge cases ------------------------------------------------------------------------

def tiny_scaled(n_cells: int = 30, n_genes: int = 12) -> ad.AnnData:
    rng = np.random.default_rng(0)
    return ad.AnnData(rng.normal(size=(n_cells, n_genes)).astype(np.float32))


def test_too_many_components_raises():
    with pytest.raises(DimReductionError, match=r"principal components \(12\) must be between 1 and 11"):
        run_pca(tiny_scaled(), PCAParams(n_comps=12))


def test_non_finite_data_raises():
    adata = tiny_scaled()
    adata.X[0, 0] = np.nan
    with pytest.raises(DimReductionError, match="NaN or infinite"):
        run_pca(adata, PCAParams(n_comps=5))


def test_neighbors_before_pca_raises():
    with pytest.raises(DimReductionError, match="Run run_pca before compute_neighbors"):
        compute_neighbors(tiny_scaled(), NeighborsParams())


def test_more_pcs_than_computed_raises():
    adata, _ = run_pca(tiny_scaled(), PCAParams(n_comps=5))
    with pytest.raises(DimReductionError, match="has only 5 components"):
        compute_neighbors(adata, NeighborsParams(n_pcs=8))


def test_more_neighbors_than_cells_raises():
    adata, _ = run_pca(tiny_scaled(), PCAParams(n_comps=5))
    with pytest.raises(DimReductionError, match="n_neighbors=50"):
        compute_neighbors(adata, NeighborsParams(n_neighbors=50, n_pcs=5))


def test_umap_before_neighbors_raises():
    with pytest.raises(DimReductionError, match="Run compute_neighbors before run_umap"):
        run_umap(tiny_scaled(), UMAPParams())


def test_inputs_are_not_modified():
    original = tiny_scaled()
    x_before = original.X.copy()
    with_pca, _ = run_pca(original, PCAParams(n_comps=5))
    assert "X_pca" not in original.obsm and np.array_equal(original.X, x_before)

    with_graph, _ = compute_neighbors(with_pca, NeighborsParams(n_neighbors=5, n_pcs=5))
    assert "neighbors" not in with_pca.uns and not with_pca.obsp

    run_umap(with_graph, UMAPParams())
    assert "X_umap" not in with_graph.obsm


# --- The elbow of the variance curve (guidance only) ---------------------------------------------

def test_elbow_on_curves_with_a_known_bend():
    # Steep drop over the first 4 PCs, then flat: the bend is at PC 4.
    assert find_elbow([10, 7, 4, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4]) == 4
    # A straight line has no bend: every point is on the line, so the first is returned.
    assert find_elbow([5, 4, 3, 2, 1]) == 1
    assert find_elbow([3, 1]) == 2 and find_elbow([]) == 0  # too short to bend


def test_elbow_is_reported_and_never_applied(scaled_pbmc3k):
    adata, summary = run_pca(scaled_pbmc3k, PCAParams(n_comps=30))
    assert 1 <= summary.elbow_pc <= 30
    assert summary.elbow_pc == find_elbow(adata.uns["pca"]["variance_ratio"])
    assert adata.obsm["X_pca"].shape[1] == 30  # all requested PCs are kept
