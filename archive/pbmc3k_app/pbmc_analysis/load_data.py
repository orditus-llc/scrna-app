"""Load an AnnData file and check that it is a usable raw-count matrix.

The checks are generic (they apply to any dataset): the matrix is non-empty, contains
raw counts (non-negative whole numbers, no NaN), and cell names are unique. Duplicate
gene names are made unique, as in the reference analysis, and the number changed is
reported.

Dataset-specific expectations (e.g. the PBMC3k shape) belong in tests, not here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


class InputValidationError(ValueError):
    """The input data can't be analyzed as raw counts. The message says why."""


@dataclass(frozen=True)
class LoadSummary:
    """What was loaded and what (if anything) was changed on loading."""

    path: str
    n_cells: int
    n_genes: int
    n_duplicate_gene_names: int  # gene names that were made unique (e.g. "TBCE" -> "TBCE-1")
    matrix_format: str  # e.g. "csr_matrix" or "ndarray"


def load_h5ad(path: str | Path) -> tuple[ad.AnnData, LoadSummary]:
    """Read an `.h5ad` file, make gene names unique, and check it holds raw counts.

    Parameters
    ----------
    path
        Path to an `.h5ad` file with cells as rows and genes as columns.

    Returns
    -------
    The loaded AnnData (a new object; the file is not changed) and a `LoadSummary`.

    Raises
    ------
    FileNotFoundError
        If the file doesn't exist.
    InputValidationError
        If the file isn't `.h5ad`, or the matrix fails a raw-count check.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    if path.suffix != ".h5ad":
        raise InputValidationError(f"Expected an .h5ad file, got '{path.name}'.")

    adata = ad.read_h5ad(path)

    n_duplicate_gene_names = int(adata.var_names.duplicated().sum())
    if n_duplicate_gene_names:
        adata.var_names_make_unique()

    validate_raw_counts(adata)

    summary = LoadSummary(
        path=str(path),
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        n_duplicate_gene_names=n_duplicate_gene_names,
        matrix_format=type(adata.X).__name__,
    )
    return adata, summary


def validate_raw_counts(adata: ad.AnnData) -> None:
    """Raise `InputValidationError` if `adata.X` isn't a usable raw-count matrix.

    Checks, in order: there are cells and genes; cell names are unique; there are no
    NaN values; values are non-negative; values are whole numbers (normalized or
    log-transformed data has fractions); and at least one count is non-zero.
    """
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise InputValidationError(
            f"The matrix is empty ({adata.n_obs:,} cells × {adata.n_vars:,} genes)."
        )
    if not adata.obs_names.is_unique:
        n_dup = int(adata.obs_names.duplicated().sum())
        raise InputValidationError(
            f"{n_dup} cell names are duplicated. Cell barcodes must be unique."
        )

    values = _stored_values(adata.X)
    if np.isnan(values).any():
        raise InputValidationError("The matrix contains NaN values; raw counts can't be NaN.")
    if (values < 0).any():
        raise InputValidationError(
            "The matrix contains negative values, so it isn't raw counts "
            "(it may be scaled or batch-corrected data)."
        )
    if not np.array_equal(values, np.round(values)):
        raise InputValidationError(
            "The matrix contains non-integer values, so it isn't raw counts "
            "(it may already be normalized or log-transformed)."
        )
    if not values.any():
        raise InputValidationError("All counts are zero.")


def _stored_values(matrix) -> np.ndarray:
    """The stored values of a dense or sparse matrix, as a flat array.

    For sparse matrices only the stored (non-zero) entries are returned; the implicit
    zeros are valid counts and don't need checking.
    """
    if sparse.issparse(matrix):
        return np.asarray(matrix.data)
    return np.asarray(matrix).ravel()
