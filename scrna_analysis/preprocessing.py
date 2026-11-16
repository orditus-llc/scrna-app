"""Preprocessing: normalization, highly variable genes, and scaling.

Order (same as the validated reference analysis), applied to QC-filtered raw counts:

1. `normalize_and_log`: keep raw counts in `layers["counts"]`, normalize each cell to the
   same total, then log1p.
2. `find_hvgs`: flag highly variable genes (adds columns only; removes nothing).
3. `keep_raw_and_subset`: store the log-normalized data for all genes in `adata.raw`,
   then keep only the HVGs in the working matrix.
4. `regress_and_scale`: optionally regress out technical covariates, then scale.

Every function returns a new AnnData; the input is never changed. In particular the
QC-filtered object passed to step 1 still holds full-gene raw counts afterwards. The
pipeline keeps it, because `layers["counts"]` here is subset to HVGs (as in the reference).
"""

from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import numpy as np
import scanpy as sc

from scrna_analysis.load_data import validate_raw_counts
from scrna_analysis.params import HVGParams, NormalizeParams, RegressOutParams, ScaleParams

SUPPORTED_HVG_FLAVORS = ("seurat",)  # flavors that expect log-normalized data


class PreprocessingError(ValueError):
    """A preprocessing step can't run on this input. The message says why."""


@dataclass(frozen=True)
class NormalizeSummary:
    n_cells: int
    n_genes: int
    target_sum: float  # counts per cell after normalization (the median when params.target_sum is None)
    median_counts_before: float  # median total counts per cell before normalization
    target_is_median: bool = False  # True when target_sum came from the data (params.target_sum None)


@dataclass(frozen=True)
class HVGSummary:
    n_genes: int
    n_hvg: int
    hvg_genes: tuple[str, ...]  # in the data's gene order
    batch_key: str | None = None  # the obs column of samples when the genes were picked within each sample
    n_batches: int = 1


@dataclass(frozen=True)
class SubsetSummary:
    raw_shape: tuple[int, int]  # all genes, log-normalized (adata.raw)
    working_shape: tuple[int, int]  # HVGs only (adata.X)


@dataclass(frozen=True)
class ScaleSummary:
    regressed: bool
    regress_keys: tuple[str, ...]  # empty if regression was switched off
    max_value: float | None
    n_values_at_clip_limit: int  # values equal to +/- max_value after scaling


def normalize_and_log(
    adata: ad.AnnData, params: NormalizeParams
) -> tuple[ad.AnnData, NormalizeSummary]:
    """Save raw counts to `layers["counts"]`, normalize to `target_sum` per cell, log1p.

    `params.target_sum=None` normalizes every cell to the median total counts of the cells
    (Scanpy's default); the value used is reported in the summary.

    Raises `PreprocessingError` if the data has already been log-transformed, and
    `InputValidationError` if `X` doesn't hold raw counts (e.g. already normalized).
    """
    if "log1p" in adata.uns:
        raise PreprocessingError("The data is already log-transformed (uns['log1p'] exists).")
    validate_raw_counts(adata)

    adata = adata.copy()
    median_before = float(np.median(np.asarray(adata.X.sum(axis=1)).ravel()))
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=params.target_sum)
    sc.pp.log1p(adata)

    totals = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    summary = NormalizeSummary(
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        # Scanpy's target when None: the median over cells with any counts.
        target_sum=params.target_sum if params.target_sum is not None else float(np.median(totals[totals > 0])),
        median_counts_before=median_before,
        target_is_median=params.target_sum is None,
    )
    return adata, summary


def find_hvgs(adata: ad.AnnData, params: HVGParams, batch_key: str | None = None) -> tuple[ad.AnnData, HVGSummary]:
    """Flag highly variable genes. Adds `var` columns only; no genes are removed.

    Adds `highly_variable`, `means`, `dispersions` and `dispersions_norm`.
    With `params.n_top_genes`, the N genes with the highest normalized dispersion are kept
    and the cut-offs are not used.
    With `batch_key` (an obs column of samples), Scanpy picks the genes within each sample: the
    cut-offs apply to the mean and normalized dispersion averaged over the samples, and top N ranks
    genes by how many samples find them variable, then by dispersion. Adds
    `highly_variable_nbatches` (in how many samples a gene is variable).
    Raises `PreprocessingError` if the data isn't log-normalized, the flavor isn't
    supported, N is not between 1 and the number of genes, or no gene passes the cut-offs.
    """
    if params.flavor not in SUPPORTED_HVG_FLAVORS:
        raise PreprocessingError(
            f"HVG flavor '{params.flavor}' isn't supported yet; use one of {SUPPORTED_HVG_FLAVORS}."
        )
    if "log1p" not in adata.uns:
        raise PreprocessingError(
            f"The '{params.flavor}' HVG flavor expects log-normalized data. "
            "Run normalize_and_log first."
        )

    n_top = params.n_top_genes
    if n_top is not None and not 1 <= n_top <= adata.n_vars:
        raise PreprocessingError(
            f"The number of highly variable genes ({n_top:,}) must be between 1 and the "
            f"{adata.n_vars:,} genes in the data."
        )

    n_batches = 1
    if batch_key is not None:
        if batch_key not in adata.obs:
            raise PreprocessingError(f"There is no '{batch_key}' column to pick genes within each sample.")
        n_batches = int(adata.obs[batch_key].astype(str).nunique())

    adata = adata.copy()
    if n_top is None:
        sc.pp.highly_variable_genes(
            adata,
            flavor=params.flavor,
            min_mean=params.min_mean,
            max_mean=params.max_mean,
            min_disp=params.min_disp,
            batch_key=batch_key,
        )
    else:
        sc.pp.highly_variable_genes(adata, flavor=params.flavor, n_top_genes=n_top, batch_key=batch_key)
    hvg_genes = tuple(adata.var_names[adata.var["highly_variable"]])
    if not hvg_genes:
        raise PreprocessingError(
            f"No genes pass the highly variable gene cut-offs (min mean {params.min_mean:g}, "
            f"max mean {params.max_mean:g}, min dispersion {params.min_disp:g}). Loosen the cut-offs."
        )
    return adata, HVGSummary(n_genes=adata.n_vars, n_hvg=len(hvg_genes), hvg_genes=hvg_genes,
                             batch_key=batch_key, n_batches=n_batches)


def keep_raw_and_subset(adata: ad.AnnData) -> tuple[ad.AnnData, SubsetSummary]:
    """Store all genes (log-normalized) in `adata.raw`, then keep only the HVGs in `X`."""
    if "highly_variable" not in adata.var:
        raise PreprocessingError("No HVG flags found. Run find_hvgs first.")

    adata = adata.copy()
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    summary = SubsetSummary(raw_shape=tuple(adata.raw.shape), working_shape=tuple(adata.shape))
    return adata, summary


def regress_and_scale(
    adata: ad.AnnData, regress: RegressOutParams, scale: ScaleParams
) -> tuple[ad.AnnData, ScaleSummary]:
    """Optionally regress out `regress.keys` (if enabled), then scale each gene.

    Raises `PreprocessingError` if a regression key is missing from `obs`, or if
    `pct_counts_mt` is requested while mitochondrial QC is unavailable.
    """
    if regress.enabled:
        _check_regress_keys(adata, regress.keys)

    adata = adata.copy()
    if regress.enabled:
        sc.pp.regress_out(adata, list(regress.keys))
    sc.pp.scale(adata, max_value=scale.max_value)

    n_clipped = 0
    if scale.max_value is not None:
        n_clipped = int((np.abs(np.asarray(adata.X)) == scale.max_value).sum())
    summary = ScaleSummary(
        regressed=regress.enabled,
        regress_keys=tuple(regress.keys) if regress.enabled else (),
        max_value=scale.max_value,
        n_values_at_clip_limit=n_clipped,
    )
    return adata, summary


def _check_regress_keys(adata: ad.AnnData, keys: tuple[str, ...]) -> None:
    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise PreprocessingError(
            f"Can't regress out {missing}: not found in obs. QC metrics such as "
            "total_counts and pct_counts_mt come from calculate_qc_metrics."
        )
    # Mito genes may not be HVGs, so look in the full-gene table when it exists.
    gene_table = adata.raw.var if adata.raw is not None else adata.var
    if "pct_counts_mt" in keys and "mt" in gene_table and not gene_table["mt"].any():
        raise PreprocessingError(
            "Can't regress out pct_counts_mt: no mitochondrial genes were found, so it is "
            "unavailable rather than 0%. Remove it from the regression keys."
        )
