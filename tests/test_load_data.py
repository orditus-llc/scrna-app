"""Tests for `scrna_analysis.load_data`.

- Reference tests: loading `data/pbmc3k.h5ad` gives the input the validated notebook used.
- Generality tests: small synthetic files, to check behavior that doesn't depend on PBMC3k.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import scanpy as sc

from scrna_analysis.load_data import (InputValidationError, load_dataset, load_h5ad, looks_like_gene_ids,
                                     validate_raw_counts)


# --- Reference tests (PBMC3k) ----------------------------------------------------------

@pytest.fixture(scope="module")
def loaded(raw_data_path):
    return load_h5ad(raw_data_path)


def test_shape_matches_reference(loaded, reference_summary):
    adata, summary = loaded
    assert list(adata.shape) == reference_summary["raw_shape"]
    assert [summary.n_cells, summary.n_genes] == reference_summary["raw_shape"]


def test_reference_cells_and_genes_are_present_in_order(loaded, reference_results):
    """The reference's filtered cells and genes are a subset of the loaded data, in order."""
    adata, _ = loaded
    for loaded_names, reference_names in [
        (adata.obs_names, reference_results.obs_names),
        (adata.var_names, reference_results.raw.var_names),
    ]:
        positions = loaded_names.get_indexer(reference_names)
        assert (positions >= 0).all(), "reference names missing from loaded data"
        assert (np.diff(positions) > 0).all(), "reference names in a different order"


def test_mito_genes_match_reference(loaded, reference_summary):
    adata, _ = loaded
    mito = adata.var_names[adata.var_names.str.startswith("MT-")].tolist()
    assert mito == reference_summary["mito_genes"]


def test_names_already_unique(loaded):
    _, summary = loaded
    assert summary.n_duplicate_gene_names == 0


def test_raw_counts_match_reference_counts_layer(loaded, reference_results):
    """Spot-check: loaded counts equal the reference's saved raw counts for the same cells/genes."""
    adata, _ = loaded
    ref_counts = reference_results.layers["counts"]  # HVG genes only (as in the reference)
    loaded_subset = adata[reference_results.obs_names, reference_results.var_names].X
    assert (loaded_subset != ref_counts).nnz == 0


# --- Generality tests (synthetic data) -------------------------------------------------

def _write(tmp_path, X, gene_names=None, cell_names=None, name="data.h5ad"):
    n_cells, n_genes = X.shape
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=cell_names or [f"cell{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=gene_names or [f"gene{j}" for j in range(n_genes)]),
    )
    path = tmp_path / name
    adata.write_h5ad(path)
    return path


@pytest.mark.parametrize("to_matrix", [np.asarray, sparse.csr_matrix], ids=["dense", "sparse"])
def test_loads_any_valid_count_matrix(tmp_path, to_matrix):
    counts = np.random.default_rng(1).poisson(1.0, size=(7, 5)).astype(np.float32)
    adata, summary = load_h5ad(_write(tmp_path, to_matrix(counts)))
    assert adata.shape == (7, 5)
    assert (summary.n_cells, summary.n_genes) == (7, 5)
    loaded_X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    np.testing.assert_array_equal(loaded_X, counts)


# AnnData warns when duplicate names are created; here they are created on purpose.
DUPLICATE_NAMES_WARNING = "ignore:.* names are not unique:UserWarning"


@pytest.mark.filterwarnings(DUPLICATE_NAMES_WARNING)
def test_duplicate_gene_names_are_made_unique_and_reported(tmp_path):
    path = _write(tmp_path, np.ones((2, 3)), gene_names=["A", "B", "A"])
    adata, summary = load_h5ad(path)
    assert summary.n_duplicate_gene_names == 1
    assert adata.var_names.tolist() == ["A", "B", "A-1"]


@pytest.mark.filterwarnings(DUPLICATE_NAMES_WARNING)
def test_input_file_is_not_modified(tmp_path):
    path = _write(tmp_path, np.ones((2, 3)), gene_names=["A", "B", "A"])
    before = path.read_bytes()
    load_h5ad(path)
    assert path.read_bytes() == before


@pytest.mark.parametrize(
    "X, message",
    [
        (np.array([[1.0, 0.5], [2.0, 3.0]]), "non-integer"),
        (np.array([[1.0, -1.0], [2.0, 3.0]]), "negative"),
        (np.array([[1.0, np.nan], [2.0, 3.0]]), "NaN"),
        (np.zeros((2, 2)), "All counts are zero"),
    ],
    ids=["normalized", "negative", "nan", "all-zero"],
)
def test_rejects_non_count_data(tmp_path, X, message):
    with pytest.raises(InputValidationError, match=message):
        load_h5ad(_write(tmp_path, X))


@pytest.mark.filterwarnings(DUPLICATE_NAMES_WARNING)
def test_rejects_duplicate_cell_names():
    adata = ad.AnnData(np.ones((2, 2)), obs=pd.DataFrame(index=["c1", "c1"]))
    with pytest.raises(InputValidationError, match="cell names are duplicated"):
        validate_raw_counts(adata)


def test_rejects_empty_matrix():
    with pytest.raises(InputValidationError, match="empty"):
        validate_raw_counts(ad.AnnData(np.zeros((0, 3))))


def test_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_h5ad(tmp_path / "does_not_exist.h5ad")


def test_wrong_file_type(tmp_path):
    path = tmp_path / "counts.csv"
    path.write_text("a,b\n1,2\n")
    with pytest.raises(InputValidationError, match=".h5ad"):
        load_h5ad(path)


# --- Other formats and datasets (generality) -------------------------------------------------

def test_10x_h5_loads_like_scanpy(mouse_10x_path):
    """A 10x .h5 file (mouse brain) loads exactly as Scanpy's own reader gives it."""
    adata, summary = load_dataset(mouse_10x_path)
    direct = sc.read_10x_h5(mouse_10x_path, gex_only=True)
    n_dup = int(direct.var_names.duplicated().sum())
    direct.var_names_make_unique()
    assert adata.shape == direct.shape
    assert adata.obs_names.equals(direct.obs_names) and adata.var_names.equals(direct.var_names)
    assert (adata.X != direct.X).nnz == 0
    assert summary.file_format == "10x Genomics (.h5)"
    assert summary.n_duplicate_gene_names == n_dup > 0  # this file has repeated gene names
    assert summary.notes == ()  # gene symbols, not IDs


def test_format_comes_from_the_original_file_name(mouse_10x_path, tmp_path):
    """An upload is stored under a temporary name; the original name decides the format."""
    stored = tmp_path / "0.h5"
    stored.write_bytes(mouse_10x_path.read_bytes())
    adata, summary = load_dataset(stored, file_name="my_sample.h5")
    assert summary.file_format == "10x Genomics (.h5)" and adata.n_obs > 0


def test_unsupported_format_is_explained(tmp_path):
    path = tmp_path / "counts.csv"
    path.write_text("cell,gene\n")
    with pytest.raises(InputValidationError, match=r"'counts.csv' isn't a supported format\. Use AnnData"):
        load_dataset(path)


def test_unreadable_file_is_explained(tmp_path):
    path = tmp_path / "broken.h5"
    path.write_bytes(b"not an hdf5 file")
    with pytest.raises(InputValidationError, match=r"'broken.h5' couldn't be read as 10x Genomics \(\.h5\)"):
        load_dataset(path)


def test_gene_ids_instead_of_symbols_get_a_note(tmp_path):
    counts = np.array([[1, 0, 3], [0, 2, 1]], dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = ["c1", "c2"]
    adata.var_names = ["ENSMUSG00000051951", "ENSMUSG00000089699", "ENSMUSG00000102343"]
    path = tmp_path / "ids.h5ad"
    adata.write_h5ad(path)
    _, summary = load_dataset(path)
    assert len(summary.notes) == 1 and "Ensembl IDs" in summary.notes[0]
    assert looks_like_gene_ids(["CD3D", "ENSG00000167286"]) is False  # half or fewer: symbols
