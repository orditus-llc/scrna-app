"""Load a single-cell data file and check that it is a usable raw-count matrix.

Formats: AnnData `.h5ad`, and the 10x Genomics `.h5` file from Cell Ranger
(filtered_feature_bc_matrix.h5; gene expression features only, named by gene symbol).
The format is recognised from the file name.

The checks are generic (they apply to any dataset): the matrix is non-empty, contains
raw counts (non-negative whole numbers, no NaN), and cell names are unique. Duplicate
gene names are made unique, as in the reference analysis, and the number changed is
reported.

Dataset-specific expectations (e.g. the PBMC3k shape) belong in tests, not here.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse

# File-name ending -> format. Shown to users as the accepted formats.
FORMATS = {".h5ad": "AnnData (.h5ad)", ".h5": "10x Genomics (.h5)"}

# Ensembl gene IDs (ENSG00000…, ENSMUSG00000…): a dataset named by these instead of gene
# symbols won't match mitochondrial-gene patterns or marker genes.
ENSEMBL_ID = re.compile(r"^ENS[A-Z]*G\d{6,}")


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
    file_format: str = FORMATS[".h5ad"]
    notes: tuple[str, ...] = field(default_factory=tuple)  # things the user should know (not errors)


def load_dataset(path: str | Path, file_name: str | None = None) -> tuple[ad.AnnData, LoadSummary]:
    """Read a `.h5ad` or 10x `.h5` file, make gene names unique, and check it holds raw counts.

    Parameters
    ----------
    path
        Where the file is stored.
    file_name
        The file's original name, if `path` is a stored copy (e.g. an upload saved as "0.h5").
        Used to recognise the format and in messages. Defaults to the name in `path`.

    Returns
    -------
    The loaded AnnData (a new object; the file is not changed) and a `LoadSummary`.

    Raises
    ------
    FileNotFoundError
        If the file doesn't exist.
    InputValidationError
        If the format isn't supported, the file can't be read, or the matrix fails a
        raw-count check.
    """
    path = Path(path)
    name = file_name or path.name
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {name}")
    suffix = Path(name).suffix.lower()
    if suffix not in FORMATS:
        raise InputValidationError(
            f"'{name}' isn't a supported format. Use {' or '.join(FORMATS.values())}."
        )

    try:
        with warnings.catch_warnings():
            # Repeated gene names are expected and handled below (made unique and reported).
            warnings.filterwarnings("ignore", message="Variable names are not unique")
            if suffix == ".h5ad":
                adata = ad.read_h5ad(path)
            else:  # 10x .h5: gene expression only (antibody or CRISPR features are left out)
                adata = sc.read_10x_h5(path, gex_only=True)
    except Exception as error:  # a damaged or mislabelled file
        raise InputValidationError(
            f"'{name}' couldn't be read as {FORMATS[suffix]}: {error}"
        ) from error

    n_duplicate_gene_names = int(adata.var_names.duplicated().sum())
    if n_duplicate_gene_names:
        adata.var_names_make_unique()

    validate_raw_counts(adata)

    notes = []
    if looks_like_gene_ids(adata.var_names):
        notes.append(
            "Genes are named by Ensembl IDs (e.g. ENSG…), not gene symbols. Mitochondrial genes "
            "and marker genes are found by symbol, so they may not be found."
        )

    summary = LoadSummary(
        path=str(path),
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        n_duplicate_gene_names=n_duplicate_gene_names,
        matrix_format=type(adata.X).__name__,
        file_format=FORMATS[suffix],
        notes=tuple(notes),
    )
    return adata, summary


def looks_like_gene_ids(names) -> bool:
    """True if most gene names are Ensembl IDs rather than gene symbols."""
    names = list(names)
    if not names:
        return False
    return sum(bool(ENSEMBL_ID.match(str(n))) for n in names) > len(names) / 2


def load_h5ad(path: str | Path) -> tuple[ad.AnnData, LoadSummary]:
    """Read an `.h5ad` file (see `load_dataset`). Raises `InputValidationError` for other formats."""
    path = Path(path)
    if path.suffix != ".h5ad":
        raise InputValidationError(f"Expected an .h5ad file, got '{path.name}'.")
    return load_dataset(path)


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
