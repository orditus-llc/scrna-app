"""Tests for `scrna_analysis.preprocessing`.

- Reference tests: QC + preprocessing on `data/pbmc3k.h5ad` with the reference preset
  reproduces the validated notebook (values read from `validation/expected/`).
- Generality tests: a PBMC3k subsample with different settings compared against direct
  Scanpy calls, a tiny hand-checkable case, and edge cases.
"""

from __future__ import annotations

from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse

from scrna_analysis.load_data import InputValidationError, load_h5ad
from scrna_analysis.params import HVGParams, NormalizeParams, QCParams, RegressOutParams, ScaleParams
from tests.reference_preset import REFERENCE_PRESET
from scrna_analysis.preprocessing import (
    PreprocessingError,
    find_hvgs,
    keep_raw_and_subset,
    normalize_and_log,
    regress_and_scale,
)
from scrna_analysis.qc import apply_filter, calculate_qc_metrics

# regress_out triggers Numba performance hints (also seen in the notebook); harmless.
pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")

# Without regression, X is still sparse when scaled; Scanpy densifies it (same as calling
# sc.pp.scale directly) and says so. Expected in the tests that switch regression off.
DENSIFY_WARNING = "ignore:zero-centering a sparse array/matrix densifies it:UserWarning"

P = REFERENCE_PRESET
TOL = 1e-6  # strict tolerance agreed for floats (CLAUDE.md, architecture decision 2)


def to_dense(matrix) -> np.ndarray:
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


def max_abs_diff(a, b) -> float:
    return float(np.max(np.abs(to_dense(a).astype(np.float64) - to_dense(b).astype(np.float64))))


# --- Reference tests (PBMC3k) ----------------------------------------------------------

@pytest.fixture(scope="module")
def raw(raw_data_path):
    adata, _ = load_h5ad(raw_data_path)
    return adata


@pytest.fixture(scope="module")
def steps(raw):
    """Run QC and every preprocessing step once with the reference preset."""
    with_qc, _ = calculate_qc_metrics(raw, P.qc)
    filtered, _ = apply_filter(with_qc, P.filter)
    normalized, norm_summary = normalize_and_log(filtered, P.normalize)
    flagged, hvg_summary = find_hvgs(normalized, P.hvg)
    subset, subset_summary = keep_raw_and_subset(flagged)
    scaled, scale_summary = regress_and_scale(subset, P.regress_out, P.scale)
    return {
        "filtered": filtered, "normalized": normalized, "flagged": flagged,
        "subset": subset, "scaled": scaled,
        "norm_summary": norm_summary, "hvg_summary": hvg_summary,
        "subset_summary": subset_summary, "scale_summary": scale_summary,
    }


def test_hvgs_match_reference(steps, reference_summary):
    summary = steps["hvg_summary"]
    assert summary.n_hvg == reference_summary["n_hvg"]
    assert list(summary.hvg_genes) == reference_summary["hvg_genes"]
    assert steps["subset"].var_names.tolist() == reference_summary["hvg_genes"]


@pytest.mark.parametrize("column", ["highly_variable", "means", "dispersions", "dispersions_norm"])
def test_hvg_statistics_match_reference(steps, reference_results, column):
    ours = steps["subset"].raw.var[column].to_numpy()
    expected = reference_results.raw.var[column].to_numpy()
    if column == "highly_variable":
        np.testing.assert_array_equal(ours, expected)
    else:
        np.testing.assert_allclose(ours, expected, rtol=0, atol=TOL, equal_nan=True)


def test_raw_log_normalized_matches_reference(steps, reference_results):
    ours, expected = steps["subset"].raw, reference_results.raw
    assert ours.shape == expected.shape
    assert ours.var_names.equals(expected.var_names)
    assert max_abs_diff(ours.X, expected.X) <= TOL


def test_counts_layer_matches_reference(steps, reference_results):
    ours = steps["scaled"].layers["counts"]
    expected = reference_results.layers["counts"]
    assert max_abs_diff(ours, expected) == 0


def test_scaled_matrix_matches_reference(steps, reference_results):
    ours, expected = steps["scaled"], reference_results
    assert ours.shape == expected.shape
    assert ours.obs_names.equals(expected.obs_names) and ours.var_names.equals(expected.var_names)
    assert max_abs_diff(ours.X, expected.X) <= TOL


@pytest.mark.parametrize(
    "column", ["n_genes_by_counts", "total_counts", "total_counts_mt", "pct_counts_mt", "n_genes"]
)
def test_obs_unchanged_by_preprocessing(steps, reference_results, column):
    np.testing.assert_array_equal(
        steps["scaled"].obs[column].to_numpy(), reference_results.obs[column].to_numpy()
    )


def test_reference_summaries(steps, reference_summary):
    n_cells, n_genes = reference_summary["filtered_shape"]
    assert steps["subset_summary"].raw_shape == (n_cells, n_genes)
    assert steps["subset_summary"].working_shape == (n_cells, reference_summary["n_hvg"])
    assert steps["scale_summary"].regressed
    assert steps["scale_summary"].regress_keys == ("total_counts", "pct_counts_mt")


def test_full_gene_counts_still_available_in_filtered_object(steps):
    """Option A: the QC-filtered input keeps full-gene raw counts after preprocessing."""
    filtered, scaled = steps["filtered"], steps["scaled"]
    assert "counts" not in filtered.layers and "log1p" not in filtered.uns  # input untouched
    assert filtered.shape == scaled.raw.shape  # all genes
    hvg_counts = filtered[:, scaled.var_names].X
    assert max_abs_diff(hvg_counts, scaled.layers["counts"]) == 0


# --- Generality tests: subsample with different settings vs direct Scanpy -------------

@pytest.mark.filterwarnings(DENSIFY_WARNING)
def test_subsample_with_other_settings_matches_direct_scanpy(raw):
    rng = np.random.default_rng(7)
    subset = raw[np.sort(rng.choice(raw.n_obs, size=600, replace=False))].copy()
    with_qc, _ = calculate_qc_metrics(subset, P.qc)
    filtered, _ = apply_filter(with_qc, P.filter)

    norm = NormalizeParams(target_sum=1e6)
    hvg = HVGParams(min_mean=0.02, max_mean=4, min_disp=0.3)
    no_regress = RegressOutParams(enabled=False)
    scale = ScaleParams(max_value=5)

    a, _ = normalize_and_log(filtered, norm)
    b, hvg_summary = find_hvgs(a, hvg)
    c, _ = keep_raw_and_subset(b)
    d, scale_summary = regress_and_scale(c, no_regress, scale)

    direct = filtered.copy()
    sc.pp.normalize_total(direct, target_sum=1e6)
    sc.pp.log1p(direct)
    assert max_abs_diff(a.X, direct.X) == 0
    sc.pp.highly_variable_genes(direct, flavor="seurat", min_mean=0.02, max_mean=4, min_disp=0.3)
    assert list(hvg_summary.hvg_genes) == direct.var_names[direct.var.highly_variable].tolist()
    direct.raw = direct
    direct = direct[:, direct.var.highly_variable].copy()
    sc.pp.scale(direct, max_value=5)

    assert d.var_names.equals(direct.var_names)
    assert max_abs_diff(d.X, direct.X) == 0
    assert max_abs_diff(d.raw.X, direct.raw.X) == 0
    assert not scale_summary.regressed and scale_summary.regress_keys == ()
    assert np.abs(d.X).max() <= 5


# --- Generality tests: tiny hand-checkable case and edge cases --------------------------

def tiny_counts() -> ad.AnnData:
    X = np.array([[1, 3, 0, 6], [2, 2, 4, 0], [5, 0, 1, 4]], dtype=np.float32)
    return ad.AnnData(
        X=sparse.csr_matrix(X),
        obs=pd.DataFrame(index=["c0", "c1", "c2"]),
        var=pd.DataFrame(index=["A", "B", "MT-1", "C"]),
    )


def test_each_cell_normalized_to_target_sum():
    adata, summary = normalize_and_log(tiny_counts(), NormalizeParams(target_sum=100))
    totals = np.expm1(to_dense(adata.X)).sum(axis=1)  # undo log1p
    np.testing.assert_allclose(totals, [100, 100, 100], rtol=1e-6)
    assert summary.median_counts_before == 10  # totals 10, 8, 10
    np.testing.assert_array_equal(to_dense(adata.layers["counts"]), to_dense(tiny_counts().X))


def test_already_normalized_input_is_refused():
    adata = tiny_counts()
    sc.pp.normalize_total(adata, target_sum=7)  # gives non-integer values
    with pytest.raises(InputValidationError, match="non-integer"):
        normalize_and_log(adata, NormalizeParams())


def test_double_log_is_refused():
    adata, _ = normalize_and_log(tiny_counts(), NormalizeParams())
    with pytest.raises(PreprocessingError, match="already log-transformed"):
        normalize_and_log(adata, NormalizeParams())


def test_hvgs_need_log_normalized_data():
    with pytest.raises(PreprocessingError, match="expects log-normalized data"):
        find_hvgs(tiny_counts(), HVGParams())


def test_unsupported_hvg_flavor_is_refused():
    adata, _ = normalize_and_log(tiny_counts(), NormalizeParams())
    with pytest.raises(PreprocessingError, match="isn't supported"):
        find_hvgs(adata, HVGParams(flavor="seurat_v3"))  # type: ignore[arg-type]


def test_no_hvgs_raises(raw):
    adata, _ = normalize_and_log(raw[:300].copy(), NormalizeParams())
    with pytest.raises(PreprocessingError, match="No genes pass the highly variable gene cut-offs"):
        find_hvgs(adata, HVGParams(min_mean=1000, max_mean=2000))


def test_subset_needs_hvg_flags():
    with pytest.raises(PreprocessingError, match="Run find_hvgs first"):
        keep_raw_and_subset(tiny_counts())


def test_missing_regression_key_raises():
    with pytest.raises(PreprocessingError, match="not found in obs"):
        regress_and_scale(tiny_counts(), RegressOutParams(enabled=True), ScaleParams())


@pytest.mark.filterwarnings("ignore:No mitochondrial genes found:UserWarning")
def test_regressing_unavailable_mito_percentage_is_refused():
    adata, _ = calculate_qc_metrics(tiny_counts(), QCParams(mito_prefix="mt-"))
    with pytest.raises(PreprocessingError, match="unavailable rather than 0%"):
        regress_and_scale(adata, RegressOutParams(enabled=True), ScaleParams())
    # Works when only available covariates are regressed out.
    _, summary = regress_and_scale(
        adata, RegressOutParams(enabled=True, keys=("total_counts",)), ScaleParams()
    )
    assert summary.regress_keys == ("total_counts",)


@pytest.mark.filterwarnings(DENSIFY_WARNING)
def test_inputs_are_not_modified(raw):
    small = raw[:300].copy()
    before = small.X.copy()
    normalized, _ = normalize_and_log(small, NormalizeParams())
    assert max_abs_diff(small.X, before) == 0 and not small.layers and "log1p" not in small.uns

    flagged, _ = find_hvgs(normalized, HVGParams())
    assert "highly_variable" not in normalized.var

    subset, _ = keep_raw_and_subset(flagged)
    assert flagged.raw is None and flagged.n_vars == small.n_vars

    x_before = subset.X.copy()
    regress_and_scale(subset, RegressOutParams(enabled=False), ScaleParams())
    assert max_abs_diff(subset.X, x_before) == 0


# --- Top N most variable genes (instead of the cut-offs) -----------------------------------------

@pytest.mark.parametrize("dataset", ["pbmc3k", "mouse"])
def test_top_n_genes_match_scanpy(dataset, raw_data_path, mouse_10x_path):
    from scrna_analysis.load_data import load_dataset
    path, prefix = (raw_data_path, "MT-") if dataset == "pbmc3k" else (mouse_10x_path, "mt-")
    adata, _ = calculate_qc_metrics(load_dataset(path)[0], QCParams(mito_prefix=prefix))
    adata, _ = apply_filter(adata, replace(REFERENCE_PRESET.filter, max_pct_counts_mt=None))
    normalized, _ = normalize_and_log(adata, NormalizeParams(target_sum=None))
    for n in (500, 2000):
        ours, summary = find_hvgs(normalized, HVGParams(n_top_genes=n))
        direct = normalized.copy()
        sc.pp.highly_variable_genes(direct, flavor="seurat", n_top_genes=n)
        assert summary.n_hvg == n == int(direct.var["highly_variable"].sum())
        assert summary.hvg_genes == tuple(direct.var_names[direct.var["highly_variable"]])


def test_top_n_genes_must_fit_the_data(steps):
    normalized = steps["normalized"]
    with pytest.raises(PreprocessingError, match=r"must be between 1 and the [\d,]+ genes"):
        find_hvgs(normalized, HVGParams(n_top_genes=normalized.n_vars + 1))
    with pytest.raises(PreprocessingError, match="must be between 1"):
        find_hvgs(normalized, HVGParams(n_top_genes=0))


# --- Highly variable genes within each sample (several samples) ------------------------------------

def test_genes_within_each_sample_match_scanpy(steps):
    """Not PBMC3k-specific: PBMC3k's cells split at random into three made-up samples. With a batch
    column the genes equal Scanpy's own per-batch selection, for the cut-offs and for top N."""
    normalized = steps["normalized"].copy()
    rng = np.random.default_rng(1)
    normalized.obs["sample"] = pd.Categorical(rng.choice(["s1", "s2", "s3"], normalized.n_obs))
    for params in (HVGParams(), HVGParams(n_top_genes=1000)):
        ours, summary = find_hvgs(normalized, params, batch_key="sample")
        direct = normalized.copy()
        kwargs = ({"n_top_genes": 1000} if params.n_top_genes else
                  {"min_mean": params.min_mean, "max_mean": params.max_mean, "min_disp": params.min_disp})
        sc.pp.highly_variable_genes(direct, flavor="seurat", batch_key="sample", **kwargs)
        assert summary.hvg_genes == tuple(direct.var_names[direct.var["highly_variable"]])
        assert (summary.batch_key, summary.n_batches) == ("sample", 3)
        assert ours.var["highly_variable_nbatches"].between(0, 3).all()
    pooled, summary = find_hvgs(normalized, HVGParams())
    assert summary.batch_key is None and summary.n_batches == 1
    assert "highly_variable_nbatches" not in pooled.var


def test_genes_within_each_sample_need_the_sample_column(steps):
    with pytest.raises(PreprocessingError, match="no 'sample' column"):
        find_hvgs(steps["normalized"], HVGParams(), batch_key="sample")


def test_a_single_dataset_ignores_the_per_sample_setting(raw_data_path):
    """The default (per_sample on) changes nothing without samples: same genes as switched off."""
    from scrna_analysis.params import DEFAULT_PARAMS
    from scrna_analysis.pipeline import AnalysisRun

    genes = {}
    for per_sample in (True, False):
        params = replace(DEFAULT_PARAMS, hvg=replace(DEFAULT_PARAMS.hvg, per_sample=per_sample))
        run = AnalysisRun(raw_data_path, params)
        for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess"):
            getattr(run, stage)()
        assert run.stage_result("preprocess").summary["hvg_per_sample"] is False
        genes[per_sample] = list(run.result("preprocess").var_names)
    assert genes[True] == genes[False]
