"""Load single-cell data and check that it is a usable raw-count matrix: one file, or several
samples (one file per sample, or one file with a sample column).

Formats: AnnData `.h5ad`, and the 10x Genomics `.h5` file from Cell Ranger
(filtered_feature_bc_matrix.h5; gene expression features only, named by gene symbol).
The format is recognized from the file name.

The checks are generic (they apply to any dataset): the matrix is non-empty, contains
raw counts (non-negative whole numbers, no NaN), and cell names are unique. Duplicate
gene names are made unique, as in the reference analysis, and the number changed is
reported.

Several samples are combined into one dataset with a `sample` column (and a `condition` column
when conditions are given), keeping the genes measured in every sample; each cell's name gets
its sample's name in front, since different files often reuse the same barcodes. A single
dataset is loaded as it always was (no sample column), so single-sample analyses don't change.

Dataset-specific expectations (e.g. the PBMC3k shape) belong in tests, not here.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from scrna_analysis.params import CONDITION_COLUMN_HINTS, MAX_SAMPLES, SAMPLE_COLUMN_HINTS

# File-name ending -> format. Shown to users as the accepted formats.
FORMATS = {".h5ad": "AnnData (.h5ad)", ".h5": "10x Genomics (.h5)"}

# Ensembl gene IDs (ENSG00000…, ENSMUSG00000…): a dataset named by these instead of gene
# symbols won't match mitochondrial-gene patterns or marker genes.
ENSEMBL_ID = re.compile(r"^ENS[A-Z]*G\d{6,}")

SAMPLE_KEY = "sample"  # obs column naming each cell's sample (only when there are several)
CONDITION_KEY = "condition"  # obs column with each cell's sample's condition (when given)
MAX_NAME_LENGTH = 40  # sample and condition names


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
    # Several samples: one entry per sample (name, condition, cells, file); empty for one dataset.
    samples: tuple[dict, ...] = field(default_factory=tuple)
    n_genes_not_in_all_samples: int = 0  # genes left out because some sample didn't measure them


def load_dataset(path: str | Path, file_name: str | None = None) -> tuple[ad.AnnData, LoadSummary]:
    """Read a `.h5ad` or 10x `.h5` file, make gene names unique, and check it holds raw counts.

    Parameters
    ----------
    path
        Where the file is stored.
    file_name
        The file's original name, if `path` is a stored copy (e.g. an upload saved as "0.h5").
        Used to recognize the format and in messages. Defaults to the name in `path`.

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
    except Exception as error:  # a damaged or mislabeled file
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


# --- Several samples -----------------------------------------------------------------------

@dataclass(frozen=True)
class SampleFile:
    """One sample's file, as the user named it (`condition` is optional, e.g. "control")."""

    path: str
    name: str
    file_name: str | None = None  # the original name, if `path` is a stored copy (an upload)
    condition: str | None = None


def clean_name(value, what: str) -> str:
    """A sample or condition name, tidied (spaces). Raises `InputValidationError` if unusable."""
    name = " ".join(str(value).split())
    if not name:
        raise InputValidationError(f"Every {what} needs a name.")
    if len(name) > MAX_NAME_LENGTH:
        raise InputValidationError(f"The {what} name '{name[:30]}…' is too long; keep it under {MAX_NAME_LENGTH} characters.")
    return name


def _check_sample_count(names: list[str]) -> None:
    if len(names) < 2:
        raise InputValidationError("Several samples need at least two; load a single file as one dataset instead.")
    if len(names) > MAX_SAMPLES:
        raise InputValidationError(f"{len(names)} samples is more than the app handles for now ({MAX_SAMPLES}).")
    repeated = sorted({n for n in names if names.count(n) > 1})
    if repeated:
        raise InputValidationError(f"Sample names must be different; repeated: {', '.join(repeated)}.")


def load_samples(files: list[SampleFile]) -> tuple[ad.AnnData, LoadSummary]:
    """Load one file per sample (each checked like `load_dataset`) and combine them: a `sample`
    column, a `condition` column if any condition is given, each cell named "<sample>:<barcode>",
    and only the genes measured in every sample (how many were left out is reported)."""
    names = [clean_name(f.name, "sample") for f in files]
    _check_sample_count(names)
    parts, samples, notes, formats = [], [], [], set()
    id_named = []
    for f, name in zip(files, names):
        adata, s = load_dataset(f.path, f.file_name)
        condition = clean_name(f.condition, "condition") if f.condition and str(f.condition).strip() else None
        adata.obs_names = [f"{name}:{barcode}" for barcode in adata.obs_names]
        adata.obs[SAMPLE_KEY] = name
        if condition is not None:
            adata.obs[CONDITION_KEY] = condition
        id_named.append(looks_like_gene_ids(adata.var_names))
        formats.add(s.file_format)
        samples.append({"name": name, "condition": condition, "n_cells": s.n_cells, "n_genes": s.n_genes,
                        "file": f.file_name or Path(f.path).name, "path": str(f.path)})
        parts.append(adata)
    if len(set(id_named)) > 1:
        raise InputValidationError("Some samples name genes by Ensembl ID and others by gene symbol, so their genes "
                                   "can't be matched. Use the same gene naming for every sample.")
    if any(s["condition"] for s in samples):
        for adata in parts:
            if CONDITION_KEY not in adata.obs:
                adata.obs[CONDITION_KEY] = ""
    all_genes = set().union(*(set(a.var_names) for a in parts))
    combined = ad.concat(parts, join="inner", merge="same")
    combined.obs[SAMPLE_KEY] = pd.Categorical(combined.obs[SAMPLE_KEY], categories=names)
    if CONDITION_KEY in combined.obs:
        combined.obs[CONDITION_KEY] = combined.obs[CONDITION_KEY].astype("category")
    if combined.n_vars == 0:
        raise InputValidationError("The samples have no genes in common, so they can't be combined.")
    dropped = len(all_genes) - combined.n_vars
    if dropped:
        notes.append(f"{dropped:,} genes measured in only some samples were left out; {combined.n_vars:,} genes "
                     "are in every sample.")
    if any(id_named):
        notes.append("Genes are named by Ensembl IDs (e.g. ENSG…), not gene symbols. Mitochondrial genes "
                     "and marker genes are found by symbol, so they may not be found.")
    validate_raw_counts(combined)
    return combined, LoadSummary(
        path=str(files[0].path), n_cells=combined.n_obs, n_genes=combined.n_vars, n_duplicate_gene_names=0,
        matrix_format=type(combined.X).__name__, file_format=" + ".join(sorted(formats)), notes=tuple(notes),
        samples=tuple(samples), n_genes_not_in_all_samples=dropped)


def sample_column_choices(adata: ad.AnnData) -> list[dict]:
    """Columns of a loaded file that could name samples (or conditions): text or category columns
    with 2 to `MAX_SAMPLES` values and none missing. Each: {column, n_values, values (up to 6),
    hint: "sample" | "condition" | None}, likely sample columns first. Only suggestions."""
    def hint(column: str) -> str | None:
        low = column.lower()
        if low in SAMPLE_COLUMN_HINTS:
            return "sample"
        if low in CONDITION_COLUMN_HINTS:
            return "condition"
        return None

    found = []
    for column in adata.obs.columns:
        values = adata.obs[column]
        if not (isinstance(values.dtype, pd.CategoricalDtype) or values.dtype == object or pd.api.types.is_string_dtype(values)):
            continue
        if values.isna().any():
            continue
        unique = pd.unique(values.astype(str))
        if 2 <= len(unique) <= MAX_SAMPLES:
            found.append({"column": str(column), "n_values": len(unique), "values": [str(v) for v in unique[:6]],
                          "hint": hint(str(column))})
    order = {"sample": 0, None: 1, "condition": 2}
    return sorted(found, key=lambda c: (order[c["hint"]], SAMPLE_COLUMN_HINTS.index(c["column"].lower())
                                        if c["hint"] == "sample" else 0, c["column"]))


def samples_from_column(adata: ad.AnnData, sample_column: str, condition_column: str | None = None,
                        file_name: str | None = None) -> tuple[ad.AnnData, tuple[dict, ...], tuple[str, ...]]:
    """Use a column of a loaded file as its samples (and optionally another as conditions). Returns
    a new AnnData with `sample` (and `condition`) columns, the samples, and notes. A sample can have
    only one condition: when a sample's cells fall under several conditions, each combination is
    its own sample ("patient_1 · ctrl"), and a note says so (the app's own rule)."""
    for column in (sample_column, condition_column):
        if column is not None and column not in adata.obs:
            raise InputValidationError(f"The file has no column named '{column}'.")
    if condition_column == sample_column:
        raise InputValidationError("Choose different columns for the samples and the conditions.")
    if adata.obs[sample_column].isna().any():
        raise InputValidationError(f"Some cells have no value in '{sample_column}', so their sample is unknown.")
    out = adata.copy()
    sample = out.obs[sample_column].astype(str).map(lambda v: " ".join(v.split()))
    notes = []
    if condition_column is not None:
        if out.obs[condition_column].isna().any():
            raise InputValidationError(f"Some cells have no value in '{condition_column}', so their condition is unknown.")
        condition = out.obs[condition_column].astype(str).map(lambda v: " ".join(v.split()))
        mixed = pd.crosstab(sample, condition).gt(0).sum(axis=1)
        if (mixed > 1).any():
            sample = sample + " · " + condition
            notes.append(f"Each value of '{sample_column}' has cells under several conditions, so each combination of "
                         f"'{sample_column}' and '{condition_column}' is treated as its own sample.")
        out.obs[CONDITION_KEY] = pd.Categorical(condition)
    names = list(pd.unique(sample))
    for name in names:
        clean_name(name, "sample")
    _check_sample_count(names)
    out.obs[SAMPLE_KEY] = pd.Categorical(sample, categories=names)
    counts = out.obs[SAMPLE_KEY].value_counts()
    first_condition = (out.obs.groupby(SAMPLE_KEY, observed=True)[CONDITION_KEY].first()
                       if CONDITION_KEY in out.obs else {})
    samples = tuple({"name": n, "condition": (str(first_condition[n]) if CONDITION_KEY in out.obs else None),
                     "n_cells": int(counts[n]), "file": file_name, "column": sample_column} for n in names)
    return out, samples, tuple(notes)
