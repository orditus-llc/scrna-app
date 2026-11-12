"""Tests for `scrna_analysis.markers`.

- Reference tests: marker genes for the validated PBMC3k clusters with the reference
  preset reproduce the validated notebook (values read from `validation/expected/`).
- Generality tests: a tiny dataset with an obvious answer, a PBMC3k subsample with a
  different method compared against direct Scanpy calls, the tidy table, and edge cases.
"""

from __future__ import annotations

from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from scrna_analysis.markers import (
    DOUBLE_DIPPING_NOTE,
    MarkerError,
    marker_table,
    rank_marker_genes,
    top_genes_per_group,
)
from scrna_analysis.params import MarkerParams
from tests.reference_preset import REFERENCE_PRESET

# regress_out (used to build the shared fixtures) triggers harmless Numba performance hints.
pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")

P = REFERENCE_PRESET
TOL = 1e-6  # strict tolerance agreed for floats (CLAUDE.md, architecture decision 2)
FIELDS = ["scores", "logfoldchanges", "pvals", "pvals_adj"]


def as_frame(result: dict, field: str) -> pd.DataFrame:
    """One field of Scanpy's results as a table: columns = clusters, rows = rank."""
    return pd.DataFrame(result[field])


# --- Reference tests (PBMC3k) ----------------------------------------------------------

@pytest.fixture(scope="module")
def ranked(clustered_pbmc3k):
    return rank_marker_genes(clustered_pbmc3k, P.markers)


def test_full_gene_ranking_matches_reference(ranked, reference_results):
    ours = as_frame(ranked[0].uns["rank_genes_groups"], "names")
    expected = as_frame(reference_results.uns["rank_genes_groups"], "names")
    assert ours.shape == expected.shape  # all genes, every cluster
    pd.testing.assert_frame_equal(ours, expected)


@pytest.mark.parametrize("field", ["scores", "logfoldchanges"])
def test_scores_and_fold_changes_match_reference(ranked, reference_results, field):
    ours = as_frame(ranked[0].uns["rank_genes_groups"], field).to_numpy(float)
    expected = as_frame(reference_results.uns["rank_genes_groups"], field).to_numpy(float)
    np.testing.assert_allclose(ours, expected, rtol=0, atol=TOL)


@pytest.mark.parametrize("field", ["pvals", "pvals_adj"])
def test_p_values_match_reference(ranked, reference_results, field):
    """p-values can be tiny (e.g. 1e-200), so compare them relative to their size."""
    ours = as_frame(ranked[0].uns["rank_genes_groups"], field).to_numpy(float)
    expected = as_frame(reference_results.uns["rank_genes_groups"], field).to_numpy(float)
    np.testing.assert_allclose(ours, expected, rtol=TOL, atol=0)


@pytest.mark.parametrize("field", ["pts", "pts_rest"])
def test_fraction_expressing_matches_reference(ranked, reference_results, field):
    ours = ranked[0].uns["rank_genes_groups"][field]
    expected = reference_results.uns["rank_genes_groups"][field]
    pd.testing.assert_frame_equal(ours.sort_index(), expected.sort_index(), check_names=False)


def test_top25_markers_match_reference_summary(ranked, reference_summary):
    _, summary = ranked
    expected = reference_summary["top25_markers"]
    assert {group: list(genes) for group, genes in summary.top_markers.items()} == expected


def test_annotation_marker_ranks_match_reference_evidence(ranked, reference_summary):
    """The ranks annotation.py will use (e.g. IL7R in cluster 0) equal the notebook's."""
    table = marker_table(ranked[0])
    rank_of = {(row.group, row.gene): row.rank for row in table.itertuples()}
    for evidence in reference_summary["annotation"]["evidence"]:
        for item in evidence["marker_DE_ranks"].split(", "):
            gene, rank = item.split(":")
            assert rank_of[(evidence["cluster"], gene)] == int(rank), (evidence["cluster"], gene)


def test_settings_recorded_like_reference(ranked, reference_results):
    ours = dict(ranked[0].uns["rank_genes_groups"]["params"])
    expected = dict(reference_results.uns["rank_genes_groups"]["params"])
    assert ours == expected


def test_summary(ranked, reference_summary):
    _, summary = ranked
    assert list(summary.groups) == list(reference_summary["cluster_sizes"])
    assert summary.n_genes_tested == reference_summary["filtered_shape"][1]  # all genes (raw)
    assert (summary.method, summary.use_raw, summary.reference) == ("wilcoxon", True, "rest")
    assert summary.notes == (DOUBLE_DIPPING_NOTE,)


# --- Generality tests: a tiny dataset with an obvious answer ----------------------------
#
# 8 cells in two groups (A, B) x 4 genes, already "log-normalized" (used as is):
#   ONLY_A is expressed in every A cell and no B cell -> must be A's #1 marker,
#   ONLY_B the reverse; the two background genes are the same everywhere.

def tiny_grouped() -> ad.AnnData:
    X = np.array([
        # ONLY_A ONLY_B BG1 BG2
        [3.0, 0.0, 1.0, 2.0],
        [2.5, 0.0, 1.0, 2.0],
        [3.5, 0.0, 1.0, 2.0],
        [2.8, 0.0, 1.0, 2.0],
        [0.0, 3.1, 1.0, 2.0],
        [0.0, 2.2, 1.0, 2.0],
        [0.0, 2.9, 1.0, 2.0],
        [0.0, 3.3, 1.0, 2.0],
    ], dtype=np.float32)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"cluster": pd.Categorical(list("AAAABBBB"))},
                         index=[f"c{i}" for i in range(8)]),
        var=pd.DataFrame(index=["ONLY_A", "ONLY_B", "BG1", "BG2"]),
    )
    adata.uns["log1p"] = {"base": None}  # mark as log data, as after normalize_and_log
    return adata


TINY_PARAMS = MarkerParams(groupby="cluster", use_raw=False)


def test_obvious_markers_rank_first():
    adata, summary = rank_marker_genes(tiny_grouped(), TINY_PARAMS)
    assert summary.top_markers["A"][0] == "ONLY_A"
    assert summary.top_markers["B"][0] == "ONLY_B"
    table = marker_table(adata)
    top_a = table[(table.group == "A") & (table["rank"] == 1)].iloc[0]
    assert (top_a.gene, top_a.frac_in_group, top_a.frac_in_rest) == ("ONLY_A", 1.0, 0.0)


def test_non_categorical_grouping_is_accepted():
    adata = tiny_grouped()
    adata.obs["cluster"] = adata.obs["cluster"].astype(str)  # plain text column
    _, summary = rank_marker_genes(adata, TINY_PARAMS)
    assert summary.top_markers["A"][0] == "ONLY_A"


# --- Generality tests: subsample with a different method vs direct Scanpy ---------------

def test_other_method_matches_direct_scanpy(clustered_pbmc3k):
    rng = np.random.default_rng(13)
    small = clustered_pbmc3k[np.sort(rng.choice(clustered_pbmc3k.n_obs, 600, replace=False))].copy()
    small.obs["leiden"] = small.obs["leiden"].cat.remove_unused_categories()
    params = replace(P.markers, method="t-test", pts=False)

    ours, _ = rank_marker_genes(small, params)
    direct = small.copy()
    sc.tl.rank_genes_groups(direct, "leiden", method="t-test", use_raw=True, pts=False,
                            reference="rest", corr_method="benjamini-hochberg", tie_correct=False)

    for field in ["names", *FIELDS]:
        pd.testing.assert_frame_equal(as_frame(ours.uns["rank_genes_groups"], field),
                                      as_frame(direct.uns["rank_genes_groups"], field))


# --- marker_table ----------------------------------------------------------------------

def test_marker_table_layout(ranked):
    table = marker_table(ranked[0], n_top=5)
    assert list(table.columns) == ["group", "rank", "gene", "score", "logfoldchange", "pval",
                                   "pval_adj", "frac_in_group", "frac_in_rest"]
    assert len(table) == 5 * len(ranked[1].groups)
    for _, rows in table.groupby("group"):
        assert rows["rank"].tolist() == [1, 2, 3, 4, 5]
        assert rows["score"].is_monotonic_decreasing  # rank order = score order


def test_marker_table_one_group(ranked):
    table = marker_table(ranked[0], group="2", n_top=3)
    assert set(table["group"]) == {"2"} and table["rank"].tolist() == [1, 2, 3]
    assert table["gene"].tolist() == list(ranked[1].top_markers["2"][:3])


def test_marker_table_errors(ranked):
    with pytest.raises(MarkerError, match="Unknown group 'nope'"):
        marker_table(ranked[0], group="nope")
    with pytest.raises(MarkerError, match="n_top must be at least 1"):
        marker_table(ranked[0], n_top=0)
    with pytest.raises(MarkerError, match="Run rank_marker_genes first"):
        marker_table(tiny_grouped())


# --- Edge cases ------------------------------------------------------------------------

def test_missing_grouping_column_raises():
    with pytest.raises(MarkerError, match="No clusters found"):
        rank_marker_genes(tiny_grouped(), MarkerParams(use_raw=False))


def test_single_group_raises():
    adata = tiny_grouped()
    adata.obs["cluster"] = pd.Categorical(["A"] * 8)
    with pytest.raises(MarkerError, match="at least 2 clusters"):
        rank_marker_genes(adata, TINY_PARAMS)


def test_one_cell_group_raises():
    adata = tiny_grouped()
    adata.obs["cluster"] = pd.Categorical(list("AAAAAAAB"))
    with pytest.raises(MarkerError, match="fewer than 2 cells"):
        rank_marker_genes(adata, TINY_PARAMS)


def test_use_raw_without_raw_raises():
    with pytest.raises(MarkerError, match=r"adata\.raw\) is missing"):
        rank_marker_genes(tiny_grouped(), replace(TINY_PARAMS, use_raw=True))


def test_unknown_reference_group_raises():
    with pytest.raises(MarkerError, match="reference='Z' is not a group"):
        rank_marker_genes(tiny_grouped(), replace(TINY_PARAMS, reference="Z"))


def test_input_is_not_modified():
    original = tiny_grouped()
    original.obs["cluster"] = original.obs["cluster"].astype(str)
    rank_marker_genes(original, TINY_PARAMS)
    assert "rank_genes_groups" not in original.uns
    assert original.obs["cluster"].dtype == object  # not converted in the input


# --- expression_by_group (dot plot values) -----------------------------------------------------

def test_expression_by_group_on_a_small_example():
    """Hand-checkable: 2 groups x 2 genes."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    from scrna_analysis.markers import expression_by_group

    X = np.array([[0.0, 2.0], [4.0, 0.0], [1.0, 1.0], [0.0, 0.0], [3.0, 0.0]])
    adata = ad.AnnData(X, obs=pd.DataFrame({"grp": pd.Categorical(["a", "a", "b", "b", "b"])},
                                           index=[f"c{i}" for i in range(5)]),
                       var=pd.DataFrame(index=["G1", "G2"]))
    table = expression_by_group(adata, ["G1", "G2"], "grp", use_raw=False).set_index(["group", "gene"])
    assert table.loc[("a", "G1"), "mean_expr"] == pytest.approx(2.0)
    assert table.loc[("a", "G1"), "frac_expressing"] == pytest.approx(0.5)
    assert table.loc[("b", "G1"), "mean_expr"] == pytest.approx(4 / 3)
    assert table.loc[("b", "G2"), "frac_expressing"] == pytest.approx(1 / 3)
    assert len(table) == 4


def test_expression_by_group_matches_a_direct_calculation(marked_pbmc3k):
    """Same numbers as averaging adata.raw per cluster with pandas."""
    import scanpy as sc
    from scrna_analysis.markers import expression_by_group

    genes = list(marked_pbmc3k.raw.var_names[:: max(1, marked_pbmc3k.raw.n_vars // 25)])
    expr = sc.get.obs_df(marked_pbmc3k, keys=genes, use_raw=True)
    clusters = marked_pbmc3k.obs["leiden"]
    mean = expr.groupby(clusters, observed=True).mean()
    frac = (expr > 0).groupby(clusters, observed=True).mean()
    table = expression_by_group(marked_pbmc3k, genes, "leiden", use_raw=True).set_index(["group", "gene"])
    assert len(table) == len(genes) * len(mean)
    for (group, gene), row in table.iterrows():
        assert row["mean_expr"] == pytest.approx(mean.loc[group, gene], abs=1e-12)
        assert row["frac_expressing"] == pytest.approx(frac.loc[group, gene], abs=1e-12)


def test_expression_by_group_reports_missing_genes(marked_pbmc3k):
    from scrna_analysis.markers import MarkerError, expression_by_group

    with pytest.raises(MarkerError, match="NOT_A_GENE"):
        expression_by_group(marked_pbmc3k, ["CD14", "NOT_A_GENE"], "leiden", use_raw=True)


# --- Heatmap of top markers (top_genes_per_group, scaled_group_means) ----------------------------

def test_top_genes_per_group_on_a_small_example():
    adata, _ = rank_marker_genes(tiny_grouped(), TINY_PARAMS)
    pairs = top_genes_per_group(adata, 1)
    assert pairs == [("A", "ONLY_A"), ("B", "ONLY_B")]
    two = top_genes_per_group(adata, 2)
    assert [g for _, g in two].count("ONLY_A") == 1 and len(two) == 4  # 4 genes, each listed once
    with pytest.raises(MarkerError, match="at least 1"):
        top_genes_per_group(adata, 0)


@pytest.mark.parametrize("n", [3, 5, 10])
def test_top_genes_per_group_takes_each_groups_best_genes_once(marked_pbmc3k, n):
    pairs = top_genes_per_group(marked_pbmc3k, n)
    genes = [g for _, g in pairs]
    assert len(genes) == len(set(genes))  # each gene once
    table = marker_table(marked_pbmc3k)
    groups = list(dict.fromkeys(table["group"]))
    assert [grp for grp, _ in pairs] == [grp for grp in groups for _ in range(n)]  # n per group, in order
    taken: set[str] = set()
    for group in groups:
        mine = [g for grp, g in pairs if grp == group]
        ranked = table[table["group"] == group].sort_values("rank")["gene"].tolist()
        # Every gene ranked above this group's last pick is either picked here or taken earlier.
        above = ranked[: ranked.index(mine[-1]) + 1]
        assert set(above) <= set(mine) | taken
        taken |= set(mine)


def test_scaled_group_means_on_a_small_example():
    """3 groups: G1 means 1, 2, 3 -> 0, 0.5, 1; G2 the same everywhere -> 0."""
    from scrna_analysis.markers import scaled_group_means

    X = np.array([[1.0, 5.0], [1.0, 5.0], [2.0, 5.0], [2.0, 5.0], [3.0, 5.0], [3.0, 5.0]])
    adata = ad.AnnData(X, obs=pd.DataFrame({"grp": pd.Categorical(list("aabbcc"))},
                                           index=[f"c{i}" for i in range(6)]),
                       var=pd.DataFrame(index=["G1", "G2"]))
    table = scaled_group_means(adata, ["G1", "G2"], "grp", use_raw=False).set_index(["group", "gene"])
    assert [table.loc[(g, "G1"), "scaled"] for g in "abc"] == pytest.approx([0.0, 0.5, 1.0])
    assert [table.loc[(g, "G2"), "scaled"] for g in "abc"] == [0.0, 0.0, 0.0]
    assert table.loc[("c", "G1"), "mean_expr"] == pytest.approx(3.0)


def test_scaled_group_means_equal_scanpys_matrixplot(marked_pbmc3k):
    """The same values as Scanpy's matrix plot of group means with standard_scale="var"."""
    from scrna_analysis.markers import scaled_group_means

    genes = [g for _, g in top_genes_per_group(marked_pbmc3k, 5)]
    ours = scaled_group_means(marked_pbmc3k, genes, "leiden", use_raw=True)
    scanpy = sc.pl.MatrixPlot(marked_pbmc3k, genes, "leiden", use_raw=True, standard_scale="var").values_df
    assert len(ours) == scanpy.size
    for row in ours.itertuples(index=False):
        assert row.scaled == pytest.approx(scanpy.loc[row.group, row.gene], abs=1e-6)


def test_one_cluster_explains_what_to_do(clustered_pbmc3k):
    one = clustered_pbmc3k.copy()
    one.obs["leiden"] = pd.Categorical(["0"] * one.n_obs)
    with pytest.raises(MarkerError, match=r"need at least 2 clusters to compare, but there is 1\. Raise the Leiden resolution"):
        rank_marker_genes(one, MarkerParams())
