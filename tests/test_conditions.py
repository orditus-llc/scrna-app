"""Compare Conditions (`scrna_analysis/conditions.py`): composition test and pseudobulk expression test.

Generality tests use hand-made data with known answers; the Kang test checks the biology recorded by
the independent Kang notebook (interferon genes higher in stimulated cells), read from its saved
summary, never retyped.
"""

from __future__ import annotations

import json
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse, stats

from scrna_analysis import conditions as cond
from scrna_analysis.params import DEFAULT_PARAMS, ConditionParams
from scrna_analysis.pipeline import AnalysisRun
from tests.conftest import PROJECT_ROOT


def samples_table(n_per_condition=3, paired=True):
    rows = []
    for i in range(n_per_condition):
        for c in ("ctrl", "treated"):
            rows.append({"sample": f"s{i}_{c}", "condition": c, "pair": f"donor{i}" if paired else f"x{i}{c}"})
    return pd.DataFrame(rows).set_index("sample")


# --- The design ---------------------------------------------------------------------------------------

def test_sample_conditions_and_pairing_choices():
    obs = pd.DataFrame({
        "sample": ["a", "a", "b", "c", "d"],
        "condition": ["x", "x", "y", "x", "y"],
        "donor": ["d1", "d1", "d1", "d2", "d2"],
        "barcode_ish": ["1", "2", "3", "4", "5"],  # one value per cell: not one per sample
        "batch": ["b1", "b1", "b2", "b3", "b4"],  # one value per sample: pairs nothing
    })
    table = cond.sample_conditions(obs, "sample", "condition", "donor")
    assert table.to_dict("index") == {"a": {"condition": "x", "pair": "d1"}, "b": {"condition": "y", "pair": "d1"},
                                      "c": {"condition": "x", "pair": "d2"}, "d": {"condition": "y", "pair": "d2"}}
    assert cond.pairing_choices(obs, "sample", "condition") == ["donor"]
    mixed = obs.assign(condition=["x", "y", "y", "x", "y"])
    with pytest.raises(cond.ConditionError, match="Sample a has cells with more than one condition"):
        cond.sample_conditions(mixed, "sample", "condition")


def test_design_refusals_and_exclusions():
    s = samples_table(3)
    d = cond.make_design(s, "ctrl", "treated", "pair")
    assert d.paired and list(d.samples["side"].value_counts().sort_index()) == [3, 3] and d.excluded == ()
    with pytest.raises(cond.ConditionError, match="Choose the two conditions"):
        cond.make_design(s, None, "treated")
    with pytest.raises(cond.ConditionError, match="There is no condition 'other'"):
        cond.make_design(s, "ctrl", "other")
    with pytest.raises(cond.ConditionError, match="two different conditions"):
        cond.make_design(s, "ctrl", "ctrl")
    with pytest.raises(cond.ConditionError, match="treated has 1 sample; .*at least 2"):
        cond.make_design(s.drop(["s1_treated", "s2_treated"]), "ctrl", "treated")
    # A third condition and an incomplete pair are left out, and said so.
    more = pd.concat([s, pd.DataFrame({"condition": ["other"], "pair": ["donor0"]}, index=["s0_other"])])
    d = cond.make_design(more.drop("s2_treated"), "ctrl", "treated", "pair")
    assert set(d.excluded) == {"s0_other", "s2_ctrl"} and len(d.samples) == 4
    with pytest.raises(cond.ConditionError, match="Only 1 subject has a sample in both"):
        cond.make_design(s.drop(["s1_treated", "s2_treated"]), "ctrl", "treated", "pair")


# --- 1. Composition (known answers) ---------------------------------------------------------------

def test_composition_test_is_a_t_test_on_logit_shares():
    s = samples_table(3)
    # Donors differ a lot (A is 30%, 60% or 45% in their controls); treatment adds ~10 points in each.
    counts = pd.DataFrame({"A": [30, 40, 60, 70, 45, 55], "B": [70, 60, 40, 30, 55, 45], "C": [0, 0, 0, 0, 0, 0]},
                          index=s.index)
    for paired in (True, False):
        d = cond.make_design(s, "ctrl", "treated", "pair" if paired else None)
        out = cond.composition_test(counts, d).set_index("group")
        assert list(out.index) == ["A", "B"]  # C has no cells
        present = counts[["A", "B"]]  # groups without cells in the compared samples are left out first
        logit = np.log(((present + 0.5).div((present + 0.5).sum(axis=1), axis=0)).pipe(lambda x: x / (1 - x)))
        ref, test = ["s0_ctrl", "s1_ctrl", "s2_ctrl"], ["s0_treated", "s1_treated", "s2_treated"]
        t = (stats.ttest_rel if paired else stats.ttest_ind)(logit.loc[test, "A"], logit.loc[ref, "A"])
        assert out.loc["A", "t"] == pytest.approx(t.statistic) and out.loc["A", "p"] == pytest.approx(t.pvalue)
        assert out.loc["A", "mean_share_reference"] == pytest.approx(np.mean([0.30, 0.60, 0.45]))
        assert out.loc["A", "mean_share_test"] == pytest.approx(np.mean([0.40, 0.70, 0.55]))
        assert out.loc["A", "t"] > 0 > out.loc["B", "t"]
        assert np.allclose(out["adj_p"], stats.false_discovery_control(out["p"], method="bh"))
    # Pairing removes the donor-to-donor differences, so the paired test is the more sensitive here.
    p_paired = cond.composition_test(counts, cond.make_design(s, "ctrl", "treated", "pair")).set_index("group")["p"]
    p_unpaired = cond.composition_test(counts, cond.make_design(s, "ctrl", "treated")).set_index("group")["p"]
    assert p_paired["A"] < p_unpaired["A"]


# --- 2. Pseudobulk expression (known answers) -----------------------------------------------------------

def synthetic(seed=0, fold=4.0, n_cells=30):
    """6 samples (3 donors × ctrl/treated), two groups; gene g0 is `fold`× higher in treated cells of
    group G1 only; the other genes have the same expected counts everywhere."""
    rng = np.random.default_rng(seed)
    s = samples_table(3)
    n_genes = 200
    base = rng.uniform(1, 10, n_genes)
    donor_effect = {f"donor{i}": rng.uniform(0.8, 1.25) for i in range(3)}
    blocks, obs = [], []
    for name, row in s.iterrows():
        for group in ("G1", "G2"):
            mean = base * donor_effect[row["pair"]]
            if group == "G1" and row["condition"] == "treated":
                mean = mean.copy(); mean[0] *= fold
            blocks.append(rng.negative_binomial(5, 5 / (5 + np.tile(mean, (n_cells, 1)))))
            obs += [{"sample": name, "condition": row["condition"], "donor": row["pair"], "group": group}] * n_cells
    obs = pd.DataFrame(obs, index=[f"c{i}" for i in range(len(obs))])
    adata = ad.AnnData(X=sparse.csr_matrix(np.vstack(blocks).astype(np.float32)), obs=obs,
                       var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))
    return adata


def test_pseudobulk_sums_each_samples_cells():
    a = synthetic()
    table, n_cells = cond.pseudobulk(a, a.obs["group"], a.obs["sample"], "G1", ["s0_ctrl", "s1_treated"])
    for s in ("s0_ctrl", "s1_treated"):
        mask = ((a.obs["group"] == "G1") & (a.obs["sample"] == s)).to_numpy()
        assert np.array_equal(table.loc[s].to_numpy(), np.asarray(a.X[mask].sum(axis=0)).ravel().astype(int))
        assert n_cells[s] == mask.sum()


def test_expression_test_finds_the_planted_change():
    a = synthetic()
    params = ConditionParams(group_by="group", reference="ctrl", test="treated", paired_by="donor")
    summary = cond.compare_conditions(a, a.obs, "sample", "condition", params)
    g1, g2 = summary.expression["G1"].results.set_index("gene"), summary.expression["G2"].results.set_index("gene")
    assert summary.expression["G1"].status == cond.TESTED and summary.expression["G1"].n_genes_tested == 200
    assert g1.index[0] == "g0" and g1.loc["g0", "adj_p"] < 1e-6 and g1.loc["g0", "log2_fold_change"] == pytest.approx(2, abs=0.3)
    assert g2.loc["g0", "adj_p"] > 0.05  # unchanged in the other group
    assert (g1.drop("g0")["adj_p"] < 0.05).sum() <= 2  # few false positives among the unchanged genes
    # The same counts give the same answer again (deterministic).
    again = cond.compare_conditions(a, a.obs, "sample", "condition", params)
    pd.testing.assert_frame_equal(summary.expression["G1"].results, again.expression["G1"].results)


def test_groups_without_enough_samples_are_skipped_with_the_reason():
    a = synthetic(n_cells=30)
    # Leave group G2 with 5 cells in two of the treated samples: too few for its pseudobulk.
    drop = []
    for s in ("s1_treated", "s2_treated"):
        drop += list(a.obs.index[(a.obs["group"] == "G2") & (a.obs["sample"] == s)][5:])
    a = a[~a.obs.index.isin(drop)].copy()
    params = ConditionParams(group_by="group", reference="ctrl", test="treated", paired_by="donor")
    summary = cond.compare_conditions(a, a.obs, "sample", "condition", params)
    g2 = summary.expression["G2"]
    assert g2.status == cond.TOO_FEW_SAMPLES and "Only 1 complete pairs" in g2.reason
    assert set(g2.samples_left_out) == {"s1_treated", "s2_treated", "s1_ctrl", "s2_ctrl"}  # partners too
    assert summary.expression["G1"].status == cond.TESTED
    only = cond.compare_conditions(a, a.obs, "sample", "condition", params, groups_to_test=["G1"])
    assert list(only.expression) == ["G1"]


# --- Kang: the interferon response (biology recorded by the independent notebook) -----------------------

@pytest.fixture(scope="module")
def kang_filtered(kang_path):
    run = AnalysisRun(kang_path, DEFAULT_PARAMS, sample_columns={"sample": "replicate", "condition": "label"},
                      metadata={"species": "human", "tissue": "blood"})
    for stage in ("load", "compute_qc_metrics", "apply_filter"):
        getattr(run, stage)()
    return run.result("filter")


def test_kang_interferon_genes_go_up_in_stimulated_cells(kang_filtered):
    """Grouped by the authors' cell types and paired by donor: every interferon gene × cell type pair
    the notebook found higher in stimulated cells comes out up, and significant, in the pseudobulk test."""
    summary = json.loads((PROJECT_ROOT / "validation" / "expected_multisample" / "kang_summary.json").read_text())
    table = pd.DataFrame(summary["interferon"]["table"])
    a = kang_filtered
    assert cond.pairing_choices(a.obs, "sample", "condition")[0] == "replicate"
    params = ConditionParams(group_by="cell_type", reference="ctrl", test="stim", paired_by="replicate")
    out = cond.compare_conditions(a, a.obs, "sample", "condition", params,
                                  groups_to_test=sorted(table["cell_type"].unique()))
    assert len(out.design.samples) == 16 and out.design.paired
    for row in table.itertuples():
        test = out.expression[row.cell_type]
        assert test.status == cond.TESTED, (row.cell_type, test.reason)
        r = test.results.set_index("gene").loc[row.gene]
        assert r["log2_fold_change"] > 0 and r["adj_p"] < 0.05, (row.cell_type, row.gene)
    comp = out.composition
    assert np.isclose(comp["mean_share_reference"].sum(), 1) and np.isclose(comp["mean_share_test"].sum(), 1)


# --- In the analysis (pipeline) ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def piece_run(kang_path, tmp_path_factory):
    """4 of Kang's samples (2 donors × ctrl/stim, 400 cells each), analyzed through Cell Type Annotation."""
    full = ad.read_h5ad(kang_path)
    donors = sorted(full.obs["replicate"].astype(str).unique())[:2]
    part = full[full.obs["replicate"].astype(str).isin(donors)]
    rng = np.random.default_rng(0)
    chosen = np.concatenate([rng.choice(np.where((part.obs["replicate"].astype(str) == d) & (part.obs["label"].astype(str) == c))[0],
                                        400, replace=False) for d in donors for c in ("ctrl", "stim")])
    path = tmp_path_factory.mktemp("kang") / "kang_piece.h5ad"
    part[np.sort(chosen)].copy().write_h5ad(path)
    run = AnalysisRun(path, DEFAULT_PARAMS, sample_columns={"sample": "replicate", "condition": "label"},
                      metadata={"species": "human", "tissue": "blood"})
    run.run_all()
    return run


def test_comparison_in_the_analysis_is_recorded_and_goes_stale(piece_run):
    run = piece_run
    params = ConditionParams(reference="ctrl", test="stim", paired_by="replicate")
    biggest = run.result("cluster").obs["leiden"].value_counts().index[0]
    summary = run.compare_conditions(params, groups=[str(biggest)])
    assert list(summary.expression) == [str(biggest)] and summary.expression[str(biggest)].status == cond.TESTED
    assert len(summary.design.samples) == 4 and summary.design.paired
    assert run.comparison_status() == "done" and run.params.conditions == params
    record = run.record()
    assert record["comparison"]["design"]["test"] == "stim" and record["comparison"]["status"] == "done"
    assert record["comparison"]["expression"][str(biggest)]["n_up"] > 0
    assert run.history[-1]["stage"] == "compare_conditions"
    json.dumps(record)  # plain data
    # Grouped by your labels: changing a label makes the comparison out of date.
    run.compare_conditions(replace(params, group_by="cell_type_reviewed"), groups=[])
    assert run.comparison_status() == "done"
    run.review_cluster(str(biggest), "rename", "my label")
    assert run.comparison_status() == "needs re-running"
    run.review_cluster(str(biggest), "clear")
    assert run.comparison_status() == "done"


def test_comparison_needs_clusters_and_conditions(piece_run, raw_data_path):
    single = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    with pytest.raises(Exception, match="Run Clustering first"):
        single.compare_conditions(ConditionParams(reference="a", test="b"))
    from scrna_analysis.pipeline import PipelineError
    with pytest.raises(PipelineError, match="Group the cells by"):
        piece_run.condition_groups("something_else")


# --- Page data ------------------------------------------------------------------------------------------

def test_page_data_and_settings_from_the_page(piece_run, raw_data_path):
    import sys
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import payloads
    import report

    run = piece_run
    page = payloads.conditions_payload(run)
    assert page["visible"] and page["enabled"]
    assert [c["name"] for c in page["options"]["conditions"]] == ["ctrl", "stim"]
    assert page["options"]["pairing"] == ["replicate"]
    assert {g["value"] for g in page["options"]["groupings"]} == {"leiden", "cell_type_predicted", "cell_type_reviewed"}
    s = page["settings"]["current"]
    assert (s["reference"], s["test"], s["paired_by"]) == ("ctrl", "stim", "replicate")  # "ctrl" looks like a control
    assert payloads.suggest_conditions(["treated", "untreated"]) == ("untreated", "treated")
    assert payloads.suggest_conditions(["a", "b"]) == ("a", "b") and payloads.suggest_conditions(["a"]) == (None, None)

    biggest = page["options"]["groupings"][0]["groups"]
    largest = max(biggest, key=lambda g: g["n_cells"])["group"]
    boxes = {**s, "group_by": "leiden", "groups": [largest]}  # an earlier test may have grouped by your labels
    params, groups = payloads.condition_params_from_browser(boxes, run.params)
    assert groups == [largest] and params.paired_by == "replicate"
    run.compare_conditions(params, groups)
    page = payloads.conditions_payload(run)
    r = page["results"]
    assert r["selected"] == largest and len(r["genes"]) == r["groups"][0]["n_genes_tested"] > 0
    assert len(r["composition"]) == len(biggest) and r["design"]["paired_by"] == "replicate"
    assert "Compared stim (2 samples) with ctrl (2 samples), paired by replicate" in page["confirmation"]
    assert any(e["label"] == "Compare Conditions" for e in payloads.run_log(run))
    assert payloads.workflow_payload(run)["conditions"]["status"] == "done"
    json.dumps(page)
    assert "propeller" in report.methods_text(run) and "PyDESeq2" in report.methods_text(run)
    html = report.run_report_html(run, "piece")
    assert "<h2>Compare Conditions</h2>" in html and "paired by replicate" in html
    for bad, message in (({**s, "test": ""}, "Choose the two conditions"), ({**s, "max_adj_p": 1}, "smaller than 1"),
                         ({**s, "min_cells": 0}, "at least 1"), ({**s, "group_by": "x"}, "Group the cells")):
        with pytest.raises(payloads.InputError, match=message):
            payloads.condition_params_from_browser(bad, run.params)
    assert payloads.condition_params_from_browser({**s, "groups": "all"}, run.params)[1] is None

    single = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    single.load()
    assert payloads.conditions_payload(single) == {"visible": False, "enabled": False, "reason": None,
                                                   "status": "not_run", "options": None, "settings": None,
                                                   "results": None}


def test_groups_are_listed_in_natural_order():
    s = samples_table(2)
    obs = pd.DataFrame({"sample": np.repeat(s.index, 3), "group": ["10", "2", "b"] * len(s)})
    d = cond.make_design(s, "ctrl", "treated")
    assert list(cond.cell_counts(obs, "sample", "group", d).columns) == ["2", "10", "b"]
