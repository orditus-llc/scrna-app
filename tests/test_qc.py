"""Tests for `scrna_analysis.qc`.

- Reference tests: QC on `data/pbmc3k.h5ad` with the reference preset reproduces the
  validated notebook (values read from `validation/expected/`).
- Generality tests: a tiny hand-calculated dataset, a PBMC3k subsample with different
  thresholds compared against direct Scanpy calls, and edge cases.
"""

from __future__ import annotations

from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse

from scrna_analysis.load_data import load_h5ad
from scrna_analysis.params import FilterParams, QCParams, SuggestionParams
from tests.reference_preset import REFERENCE_PRESET
from scrna_analysis.qc import (
    genes_with_prefix,
    mito_prefix_matches,
    suggest_filter,
    suggest_mito_prefix,
    QCError,
    apply_filter,
    calculate_qc_metrics,
    cell_filter_failures,
    preview_filter,
)

QC = REFERENCE_PRESET.qc
FILTER = REFERENCE_PRESET.filter
OBS_QC_COLUMNS = ["n_genes_by_counts", "total_counts", "total_counts_mt", "pct_counts_mt", "n_genes"]
VAR_QC_COLUMNS = ["mt", "n_cells_by_counts", "mean_counts", "pct_dropout_by_counts",
                  "total_counts", "n_cells"]


# --- Reference tests (PBMC3k) ----------------------------------------------------------

@pytest.fixture(scope="module")
def raw(raw_data_path):
    adata, _ = load_h5ad(raw_data_path)
    return adata


@pytest.fixture(scope="module")
def with_qc(raw):
    return calculate_qc_metrics(raw, QC)


@pytest.fixture(scope="module")
def filtered(with_qc):
    adata, _ = with_qc
    return apply_filter(adata, FILTER)


def test_mito_genes_match_reference(with_qc, reference_summary):
    _, summary = with_qc
    assert list(summary.mito_genes) == reference_summary["mito_genes"]
    assert summary.n_mito_genes == reference_summary["n_mito_genes"]
    assert summary.mito_qc_available
    assert summary.warnings == ()


def test_qc_statistics_match_reference(with_qc, reference_summary):
    _, summary = with_qc
    expected = reference_summary["qc_describe"]
    assert set(summary.metric_stats) == set(expected)
    for metric, stats in expected.items():
        for name, value in stats.items():
            assert summary.metric_stats[metric][name] == pytest.approx(value, rel=0, abs=1e-6), (
                f"{metric} {name}"
            )


def test_filtered_shape_matches_reference(filtered, reference_summary):
    adata, preview = filtered
    assert list(adata.shape) == reference_summary["filtered_shape"]
    assert [preview.n_cells_kept, preview.n_genes_kept] == reference_summary["filtered_shape"]


def test_kept_cells_and_genes_match_reference(filtered, reference_results):
    adata, _ = filtered
    assert adata.obs_names.equals(reference_results.obs_names)
    assert adata.var_names.equals(reference_results.raw.var_names)


@pytest.mark.parametrize("column", OBS_QC_COLUMNS)
def test_cell_qc_values_match_reference(filtered, reference_results, column):
    """Exact equality also proves the metrics were not recalculated after filtering."""
    adata, _ = filtered
    expected = reference_results.obs[column]
    assert adata.obs[column].dtype == expected.dtype
    np.testing.assert_array_equal(adata.obs[column].to_numpy(), expected.to_numpy())


@pytest.mark.parametrize("column", VAR_QC_COLUMNS)
def test_gene_qc_values_match_reference(filtered, reference_results, column):
    adata, _ = filtered
    expected = reference_results.raw.var[column]
    assert adata.var[column].dtype == expected.dtype
    np.testing.assert_array_equal(adata.var[column].to_numpy(), expected.to_numpy())


def test_preview_and_apply_keep_identical_cells_and_genes(with_qc, filtered):
    adata, _ = with_qc
    preview = preview_filter(adata, FILTER)
    result, _ = filtered
    assert result.obs_names.equals(preview.kept_cells)
    assert result.var_names.equals(preview.kept_genes)


def test_reference_removal_breakdown_adds_up(filtered, reference_summary):
    _, preview = filtered
    removed = preview.removed_in_order
    n_cells_removed = sum(v for k, v in removed.items() if k.startswith("cells"))
    assert preview.n_cells_before - n_cells_removed == reference_summary["filtered_shape"][0]
    assert preview.n_genes_before - removed["genes: min_cells"] == reference_summary["filtered_shape"][1]


# --- Generality tests: tiny hand-calculated dataset -------------------------------------
#
# Thresholds: min_genes=2, min_cells=2, max_n_genes_by_counts=5, max_pct_counts_mt=20.
#
#        G0 G1 G2 G3 G4 G5 MT-1 MT-2 | genes total  %mt  | outcome
#   c0    1  1  1  0  0  0   0    0  |   3     3     0.0 | keep
#   c1    2  1  1  0  0  0   1    0  |   4     5    20.0 | removed: max_pct_counts_mt (20 is not < 20)
#   c2    0  0  0  5  0  0   0    0  |   1     5     0.0 | removed: min_genes (1 < 2)
#   c3    1  1  1  1  1  1   0    0  |   6     6     0.0 | removed: max_n_genes_by_counts (6 is not < 5)
#   c4    4  4  0  0  0  0   0    1  |   3     9    11.1 | keep
#   c5    0  1  0  0  1  0   0    0  |   2     2     0.0 | keep
#
# Genes are counted among cells that pass min_genes (c0, c1, c3, c4, c5):
#   G0 4, G1 5, G2 3, G4 2 -> keep;  G3 1, G5 1, MT-1 1, MT-2 1 -> removed.
#   G3 is detected in 2 cells overall, so it would wrongly be kept if genes were counted
#   before the cell filter. This checks the order of the steps.

TINY_GENES = ["G0", "G1", "G2", "G3", "G4", "G5", "MT-1", "MT-2"]
TINY_COUNTS = np.array([
    [1, 1, 1, 0, 0, 0, 0, 0],
    [2, 1, 1, 0, 0, 0, 1, 0],
    [0, 0, 0, 5, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 1, 0, 0],
    [4, 4, 0, 0, 0, 0, 0, 1],
    [0, 1, 0, 0, 1, 0, 0, 0],
], dtype=np.float32)
TINY_FILTER = FilterParams(min_genes=2, min_cells=2, max_n_genes_by_counts=5, max_pct_counts_mt=20)


def tiny_adata(as_sparse: bool = True) -> ad.AnnData:
    X = sparse.csr_matrix(TINY_COUNTS) if as_sparse else TINY_COUNTS.copy()
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(6)]),
        var=pd.DataFrame(index=TINY_GENES),
    )


@pytest.mark.parametrize("as_sparse", [True, False], ids=["sparse", "dense"])
def test_tiny_qc_metrics_by_hand(as_sparse):
    adata, summary = calculate_qc_metrics(tiny_adata(as_sparse), QC)
    assert summary.mito_genes == ("MT-1", "MT-2")
    np.testing.assert_array_equal(adata.obs["n_genes_by_counts"], [3, 4, 1, 6, 3, 2])
    np.testing.assert_array_equal(adata.obs["total_counts"], [3, 5, 5, 6, 9, 2])
    np.testing.assert_allclose(adata.obs["pct_counts_mt"], [0, 20, 0, 0, 100 / 9, 0], rtol=1e-6)


@pytest.mark.parametrize("as_sparse", [True, False], ids=["sparse", "dense"])
def test_tiny_filter_by_hand(as_sparse):
    adata, _ = calculate_qc_metrics(tiny_adata(as_sparse), QC)
    filtered, preview = apply_filter(adata, TINY_FILTER)
    assert filtered.obs_names.tolist() == ["c0", "c4", "c5"]
    assert filtered.var_names.tolist() == ["G0", "G1", "G2", "G4"]
    assert preview.removed_in_order == {
        "cells: min_genes": 1,  # c2
        "genes: min_cells": 4,  # G3, G5, MT-1, MT-2
        "cells: max_n_genes_by_counts": 1,  # c3
        "cells: max_pct_counts_mt": 1,  # c1 (exactly at the limit)
    }
    assert preview.cells_failing_each == {
        "min_genes": 1, "max_n_genes_by_counts": 1, "max_pct_counts_mt": 1,
    }


def test_tiny_disabled_upper_limits_keep_those_cells():
    adata, _ = calculate_qc_metrics(tiny_adata(), QC)
    params = replace(TINY_FILTER, max_n_genes_by_counts=None, max_pct_counts_mt=None)
    filtered, preview = apply_filter(adata, params)
    assert filtered.obs_names.tolist() == ["c0", "c1", "c3", "c4", "c5"]  # only c2 removed
    assert set(preview.removed_in_order) == {"cells: min_genes", "genes: min_cells"}


# --- Generality tests: PBMC3k subsample vs direct Scanpy calls --------------------------

def test_subsample_with_other_thresholds_matches_direct_scanpy(raw):
    rng = np.random.default_rng(42)
    subset = raw[np.sort(rng.choice(raw.n_obs, size=500, replace=False))].copy()
    params = FilterParams(min_genes=300, min_cells=5, max_n_genes_by_counts=2000, max_pct_counts_mt=5.0)

    adata, _ = calculate_qc_metrics(subset, QC)
    ours, preview = apply_filter(adata, params)

    direct = subset.copy()
    direct.var["mt"] = direct.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(direct, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(direct, min_genes=300)
    sc.pp.filter_genes(direct, min_cells=5)
    direct = direct[(direct.obs.n_genes_by_counts < 2000) & (direct.obs.pct_counts_mt < 5.0)].copy()

    assert ours.obs_names.equals(direct.obs_names)
    assert ours.var_names.equals(direct.var_names)
    assert preview.kept_cells.equals(direct.obs_names)
    assert preview.kept_genes.equals(direct.var_names)
    assert (ours.X != direct.X).nnz == 0


# --- Edge cases ------------------------------------------------------------------------

def test_no_mito_genes_warns_and_marks_mito_qc_unavailable():
    with pytest.warns(UserWarning, match="No mitochondrial genes found"):
        adata, summary = calculate_qc_metrics(tiny_adata(), QCParams(mito_prefix="mt-"))
    assert summary.mito_genes == ()
    assert not summary.mito_qc_available
    assert "pct_counts_mt" not in summary.metric_stats  # not reported as a fake 0%
    assert summary.warnings and "switch off the % mitochondrial filter" in summary.warnings[0]


@pytest.mark.filterwarnings("ignore:No mitochondrial genes found:UserWarning")
def test_mito_filter_refused_when_mito_qc_unavailable():
    adata, _ = calculate_qc_metrics(tiny_adata(), QCParams(mito_prefix="mt-"))
    for function in (preview_filter, apply_filter):
        with pytest.raises(QCError, match="no mitochondrial genes were found"):
            function(adata, TINY_FILTER)
    # Works once the mitochondrial filter is explicitly switched off.
    filtered, _ = apply_filter(adata, replace(TINY_FILTER, max_pct_counts_mt=None))
    assert "c1" in filtered.obs_names  # c1 is no longer removed by a mito filter


def test_genuine_zero_percent_mito_is_not_unavailable():
    """Cells with 0% mito counts are real measurements when mito genes exist."""
    adata, summary = calculate_qc_metrics(tiny_adata(), QC)
    assert summary.mito_qc_available
    assert summary.metric_stats["pct_counts_mt"]["min"] == 0.0


def test_filter_removing_all_cells_raises():
    adata, _ = calculate_qc_metrics(tiny_adata(), QC)
    with pytest.raises(QCError, match="would remove all cells"):
        apply_filter(adata, replace(TINY_FILTER, min_genes=100))


def test_filter_without_qc_metrics_raises():
    with pytest.raises(QCError, match="Run calculate_qc_metrics"):
        preview_filter(tiny_adata(), TINY_FILTER)


def test_unsupported_qc_vars_raise():
    with pytest.raises(QCError, match="Only the mitochondrial QC variable"):
        calculate_qc_metrics(tiny_adata(), QCParams(qc_vars=("mt", "ribo")))


def test_inputs_are_not_modified():
    original = tiny_adata()
    adata, _ = calculate_qc_metrics(original, QC)
    assert list(original.obs.columns) == [] and list(original.var.columns) == []

    obs_before, var_before = adata.obs.copy(), adata.var.copy()
    preview_filter(adata, TINY_FILTER)
    apply_filter(adata, TINY_FILTER)
    assert adata.shape == (6, 8)
    pd.testing.assert_frame_equal(adata.obs, obs_before)
    pd.testing.assert_frame_equal(adata.var, var_before)


# --- cell_filter_failures (inspecting cells) -------------------------------------------------

def test_tiny_cell_failures_by_hand():
    """genes per cell [3, 4, 1, 6, 3, 2]; % mito [0, 20, 0, 0, 11.1, 0]."""
    adata, _ = calculate_qc_metrics(tiny_adata(), QC)
    failures = cell_filter_failures(adata, replace(TINY_FILTER, max_n_genes_by_counts=4))
    assert failures.index.tolist() == ["c0", "c1", "c2", "c3", "c4", "c5"]
    assert failures.to_dict("list") == {
        "min_genes": [False, False, True, False, False, False],  # c2: 1 gene < 2
        "max_n_genes_by_counts": [False, True, False, True, False, False],  # c1: 4, c3: 6 (>= 4)
        "max_pct_counts_mt": [False, True, False, False, False, False],  # c1: exactly 20%
    }
    # c1 fails two filters at once; every cell with no failure is kept by the filter.
    kept = failures.index[~failures.any(axis=1)]
    assert kept.tolist() == preview_filter(adata, replace(TINY_FILTER, max_n_genes_by_counts=4)).kept_cells.tolist()


def test_switched_off_filters_have_no_column():
    adata, _ = calculate_qc_metrics(tiny_adata(), QC)
    failures = cell_filter_failures(adata, replace(TINY_FILTER, max_n_genes_by_counts=None, max_pct_counts_mt=None))
    assert failures.columns.tolist() == ["min_genes"]


@pytest.mark.parametrize("params", [
    FILTER,
    replace(FILTER, min_genes=300, max_pct_counts_mt=5.0),
    replace(FILTER, max_n_genes_by_counts=None, max_pct_counts_mt=3.5),
], ids=["reference", "stricter", "no-max-genes"])
def test_cells_without_failures_are_exactly_the_kept_cells(with_qc, params):
    adata, _ = with_qc
    failures = cell_filter_failures(adata, params)
    preview = preview_filter(adata, params)
    assert failures.index.equals(adata.obs_names)
    assert failures.index[~failures.any(axis=1)].equals(preview.kept_cells)
    # The per-filter counts agree with the preview's own count of cells failing each one.
    assert {name: int(col.sum()) for name, col in failures.items()} == preview.cells_failing_each


def test_cell_failures_leave_input_unchanged(with_qc):
    adata, _ = with_qc
    before = adata.obs.copy()
    cell_filter_failures(adata, FILTER)
    pd.testing.assert_frame_equal(adata.obs, before)


# --- Mitochondrial gene patterns for any species ---------------------------------------------

def test_pbmc3k_suggests_the_human_pattern(raw_data_path, reference_summary):
    adata, _ = load_h5ad(raw_data_path)
    matches = mito_prefix_matches(adata.var_names)
    assert matches[0]["prefix"] == "MT-" and matches[0]["n_genes"] == reference_summary["n_mito_genes"]
    assert suggest_mito_prefix(adata.var_names) == "MT-" == REFERENCE_PRESET.qc.mito_prefix


def test_mouse_10x_suggests_the_mouse_pattern(mouse_10x_path):
    from scrna_analysis.load_data import load_dataset
    adata, _ = load_dataset(mouse_10x_path)
    expected = int(adata.var_names.str.startswith("mt-").sum())
    assert expected > 0
    assert suggest_mito_prefix(adata.var_names) == "mt-"
    assert mito_prefix_matches(adata.var_names)[0]["n_genes"] == expected


def test_patterns_are_case_sensitive_and_custom_ones_work():
    names = ["mt:CoI", "mt:ND1", "Act5C", "MT-like", "mtfoo"]
    assert genes_with_prefix(names, "mt:") == ("mt:CoI", "mt:ND1")
    assert genes_with_prefix(names, "MT-") == ("MT-like",)
    assert genes_with_prefix(names, "mt-") == ()
    assert genes_with_prefix(names, "") == ()
    assert suggest_mito_prefix(names) == "mt:"  # fruit fly: 2 matches beats MT- with 1
    assert suggest_mito_prefix(["GAPDH", "ACTB"]) is None
    assert mito_prefix_matches(["GAPDH", "ACTB"]) == []


# --- Thresholds suggested from the data (median ± MADs) -----------------------------------------

def expected_thresholds(obs, n_mads_genes, n_mads_mito):
    """The same rule written out independently: log1p genes ± N MADs, % mito + N MADs."""
    from scipy.stats import median_abs_deviation
    g = np.log1p(obs["n_genes_by_counts"].to_numpy(dtype=float))
    lo, hi = np.median(g) - n_mads_genes * median_abs_deviation(g), np.median(g) + n_mads_genes * median_abs_deviation(g)
    m = obs["pct_counts_mt"].to_numpy(dtype=float)
    mt = np.median(m) + n_mads_mito * median_abs_deviation(m)
    return int(np.ceil(np.expm1(lo))), int(np.floor(np.expm1(hi))) + 1, float(np.ceil(mt * 10) / 10)


@pytest.mark.parametrize("dataset", ["pbmc3k", "mouse"])
def test_suggested_thresholds_follow_the_mad_rule(dataset, raw_data_path, mouse_10x_path):
    from scrna_analysis.load_data import load_dataset
    path, prefix = (raw_data_path, "MT-") if dataset == "pbmc3k" else (mouse_10x_path, "mt-")
    adata, _ = calculate_qc_metrics(load_dataset(path)[0], QCParams(mito_prefix=prefix))
    for genes_k, mito_k in [(5, 3), (3, 5)]:
        s = suggest_filter(adata, SuggestionParams(n_mads_genes=genes_k, n_mads_mito=mito_k))
        assert (s.filter.min_genes, s.filter.max_n_genes_by_counts, s.filter.max_pct_counts_mt) == \
            expected_thresholds(adata.obs, genes_k, mito_k)
        assert s.filter.min_cells == 3
        assert s.details["max_pct_counts_mt"]["n_mads"] == mito_k and s.details["min_genes"]["scale"] == "log1p"
        apply_filter(adata, s.filter)  # usable as-is


def test_suggestion_adapts_to_the_data(raw_data_path, mouse_10x_path):
    """Brain cells have a higher mitochondrial share than blood cells, so their suggested
    limit is higher; a fixed blood-based limit would remove most brain cells."""
    from scrna_analysis.load_data import load_dataset
    blood, _ = calculate_qc_metrics(load_dataset(raw_data_path)[0], QCParams(mito_prefix="MT-"))
    brain, _ = calculate_qc_metrics(load_dataset(mouse_10x_path)[0], QCParams(mito_prefix="mt-"))
    sp = SuggestionParams()
    assert suggest_filter(brain, sp).filter.max_pct_counts_mt > suggest_filter(blood, sp).filter.max_pct_counts_mt
    kept_fixed = preview_filter(brain, REFERENCE_PRESET.filter).n_cells_kept
    kept_suggested = preview_filter(brain, suggest_filter(brain, sp).filter).n_cells_kept
    assert kept_suggested > 2 * kept_fixed


def test_no_mito_genes_or_no_spread_gives_no_limit():
    adata, _ = calculate_qc_metrics(tiny_adata(), QCParams(mito_prefix="mt-"))  # no "mt-" genes
    s = suggest_filter(adata, SuggestionParams())
    assert s.filter.max_pct_counts_mt is None and "No mitochondrial genes" in s.details["max_pct_counts_mt"]["note"]
    # Every cell with the same number of genes: MAD = 0, so no gene limits either.
    same = ad.AnnData(np.ones((5, 4), dtype=np.float32))
    same.var_names = ["MT-1", "A", "B", "C"]
    same, _ = calculate_qc_metrics(same, QCParams())
    s = suggest_filter(same, SuggestionParams())
    assert s.filter.min_genes == 0 and s.filter.max_n_genes_by_counts is None
    assert "barely vary" in s.details["min_genes"]["note"]


def test_suggestion_needs_qc_metrics():
    with pytest.raises(QCError, match="calculate QC metrics before suggesting"):
        suggest_filter(tiny_adata(), SuggestionParams())
