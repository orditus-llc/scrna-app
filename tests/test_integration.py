"""Tests for batch correction and the sample mixing score (scrna_analysis/integration.py) and how the
pipeline uses them. Known answers on made-up coordinates; Harmony compared with harmonypy called
directly; a small piece of the Kang et al. 2018 data (4 samples) for the whole flow; and a check that a
single dataset is analyzed exactly as before."""

from __future__ import annotations

import logging
from dataclasses import replace

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
import pytest

from scrna_analysis import integration
from scrna_analysis.params import DEFAULT_PARAMS, IntegrationParams
from scrna_analysis.pipeline import AnalysisRun, PipelineError


# --- The mixing score (known answers) -------------------------------------------------------------

def test_mixing_is_zero_when_samples_sit_apart_and_near_one_when_mixed():
    rng = np.random.default_rng(0)
    apart = np.vstack([rng.normal(0, 1, (200, 5)), rng.normal(50, 1, (200, 5))])
    labels = ["a"] * 200 + ["b"] * 200
    assert integration.sample_mixing(apart, labels).max() == 0.0
    mixed = rng.normal(0, 1, (400, 5))
    shuffled = rng.permutation(labels)
    assert integration.sample_mixing(mixed, shuffled).mean() > 0.9
    assert (integration.sample_mixing(mixed, ["a"] * 400) == 1).all()  # one sample: nothing to mix


def test_mixing_by_group_summarizes_each_cluster():
    a = ad.AnnData(np.zeros((6, 1)))
    a.obs["leiden"] = ["0", "0", "0", "1", "1", "1"]
    a.obs["sample"] = ["x", "x", "y", "x", "x", "x"]
    a.obs[integration.MIXING_KEY] = [0.9, 0.8, 0.7, 0.1, 0.2, 0.3]
    table = integration.mixing_by_group(a, "leiden", "sample").set_index("group")
    assert table.loc["0", "n_samples"] == 2 and table.loc["0", "top_sample"] == "x"
    assert table.loc["0", "top_share"] == pytest.approx(2 / 3) and not table.loc["0", "low"]
    assert table.loc["1", "low"] and table.loc["1", "mixing"] == pytest.approx(0.2)


# --- Harmony --------------------------------------------------------------------------------------

def two_batches(n=50, seed=0):
    """Made-up PCA coordinates like real data (centered): three cell types in different directions,
    each present in both batches, with batch b shifted a little. Batch a = ctrl, b = stim."""
    rng = np.random.default_rng(seed)
    centers = rng.normal(0, 1, (3, 8))
    centers = 10 * centers / np.linalg.norm(centers, axis=1, keepdims=True)
    shift = rng.normal(0, 1, 8)
    shift = 3 * shift / np.linalg.norm(shift)
    parts, batch = [], []
    for center in centers:
        for name, k in (("a", 0), ("b", 1)):
            parts.append(rng.normal(0, 1, (n, 8)) + center + k * shift)
            batch += [name] * n
    coords = np.vstack(parts)
    a = ad.AnnData(np.zeros((len(batch), 1)))
    a.obsm["X_pca"] = coords - coords.mean(axis=0)
    a.obs["sample"] = batch
    a.obs["condition"] = ["ctrl" if b == "a" else "stim" for b in batch]
    return a


def test_harmony_equals_harmonypy_called_directly_and_mixes_the_batches():
    a = two_batches()
    out, summary = integration.run_harmony(a, IntegrationParams(method="harmony"), n_components=8)
    logging.getLogger("harmonypy").setLevel(logging.WARNING)
    direct = harmonypy.run_harmony(a.obsm["X_pca"], pd.DataFrame({"batch": a.obs["sample"].to_numpy()}), ["batch"],
                                   theta=2.0, max_iter_harmony=10, random_state=0, verbose=False)
    assert np.array_equal(out.obsm[integration.CORRECTED_KEY], np.asarray(direct.Z_corr))
    assert "X_pca_harmony" not in a.obsm  # the input is unchanged
    before = integration.sample_mixing(a.obsm["X_pca"], a.obs["sample"]).mean()
    after = integration.sample_mixing(out.obsm[integration.CORRECTED_KEY], a.obs["sample"]).mean()
    assert after > before + 0.25 and after > 0.9
    assert (summary.batch_key, summary.n_batches, summary.n_components) == ("sample", 2, 8)
    # Each sample here is exactly one condition, so batches and conditions can't be told apart.
    assert summary.warnings and "can't be told apart" in summary.warnings[0]


def test_bad_batch_choices_are_explained():
    a = two_batches()
    for key, message in [("nope", "no column named"), ("one", "only one value")]:
        if key == "one":
            a.obs["one"] = "x"
        with pytest.raises(integration.IntegrationError, match=message):
            integration.run_harmony(a, IntegrationParams(method="harmony", batch_key=key), n_components=8)
    a.obs["gaps"] = pd.Series(["x", None] * (a.n_obs // 2), index=a.obs_names, dtype=object)
    with pytest.raises(integration.IntegrationError, match="no value"):
        integration.check_batch_key(a, "gaps")


def test_condition_warning():
    a = two_batches()
    assert "condition itself" in integration.condition_warning(a, "condition")
    a.obs["donor"] = ["d1", "d2"] * (a.n_obs // 2)  # every donor has both conditions: fine
    assert integration.condition_warning(a, "donor") is None
    b = two_batches()
    del b.obs["condition"]
    assert integration.condition_warning(b, "sample") is None
    assert integration.batch_column_choices(a) == ["sample", "condition", "donor"]


# --- In the analysis ------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def kang_piece(kang_path, tmp_path_factory):
    """4 of Kang's samples (2 donors × ctrl/stim), 400 cells each, as one file with sample columns."""
    full = ad.read_h5ad(kang_path)
    donors = sorted(full.obs["replicate"].astype(str).unique())[:2]
    keep = full.obs["replicate"].astype(str).isin(donors)
    part = full[keep]
    rng = np.random.default_rng(0)
    chosen = np.concatenate([rng.choice(np.where((part.obs["replicate"].astype(str) == d) & (part.obs["label"].astype(str) == c))[0],
                                        400, replace=False)
                             for d in donors for c in sorted(part.obs["label"].astype(str).unique())])
    path = tmp_path_factory.mktemp("kang") / "kang_piece.h5ad"
    part[np.sort(chosen)].copy().write_h5ad(path)
    return path


def analyzed(path, integration_params=None):
    run = AnalysisRun(path, DEFAULT_PARAMS, sample_columns={"sample": "replicate", "condition": "label"},
                      metadata={"species": "human", "tissues": ["Blood"]})
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess"):
        getattr(run, stage)()
    run.reduce_dimensions(integration_params=integration_params)
    run.cluster()
    return run


def test_harmony_in_the_analysis_mixes_samples_and_is_recorded(kang_piece):
    plain = analyzed(kang_piece)
    corrected = analyzed(kang_piece, IntegrationParams(method="harmony"))
    s0 = plain.stage_result("reduce_dimensions").summary["integration"]
    s1 = corrected.stage_result("reduce_dimensions").summary["integration"]
    assert s0["method"] == "none" and s0["mixing_after"] is None and s0["graph_from"] == "X_pca"
    assert s1["method"] == "harmony" and s1["graph_from"] == "X_pca_harmony" and s1["n_batches"] == 4
    assert s1["mixing_before"] == pytest.approx(s0["mixing_before"])  # same PCA, same starting point
    assert s1["mixing_after"] > s1["mixing_before"]
    out = corrected.result("reduce_dimensions")
    assert {integration.MIXING_KEY, integration.MIXING_BEFORE_KEY} <= set(out.obs)
    assert np.allclose(out.obs[integration.MIXING_KEY].mean(), s1["mixing_after"])
    mixing = corrected.stage_result("cluster").summary["mixing"]
    assert sum(r["n_cells"] for r in mixing) == out.n_obs
    record = corrected.record()
    assert record["params"]["integration"]["method"] == "harmony"
    assert record["history"][-2]["settings"]["integration"]["batch_key"] == "sample"
    # Expression is untouched: the same normalized values with and without correction.
    assert (corrected.result("reduce_dimensions").raw.X != plain.result("reduce_dimensions").raw.X).nnz == 0


def test_a_single_dataset_is_analyzed_as_before(raw_data_path):
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster"):
        getattr(run, stage)()
    assert "integration" not in run.stage_result("reduce_dimensions").summary
    assert "mixing" not in run.stage_result("cluster").summary
    assert integration.MIXING_KEY not in run.result().obs
    with pytest.raises(integration.IntegrationError, match="no column named"):
        run.reduce_dimensions(integration_params=IntegrationParams(method="harmony"))
    assert run.status()["cluster"] == "done"  # the failed attempt changed nothing


def test_page_data_for_batch_correction(kang_piece):
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import payloads
    import report

    run = analyzed(kang_piece)
    red = payloads.reduction_payload(run)
    choices = [c["value"] for c in red["batch"]["choices"]]
    assert choices[:2] == ["sample", "condition"]
    by_value = {c["value"]: c for c in red["batch"]["choices"]}
    assert by_value["sample"]["warning"] is None and "condition itself" in by_value["condition"]["warning"]
    integ = red["results"]["integration"]
    assert integ["method"] == "none" and len(integ["batches"]) == 4 and len(integ["cell_batch"]) == run.result().n_obs
    assert integ["suggest"] == (integ["mixing_before"] < payloads.LOW_MIXING)
    mixing = payloads.cluster_payload(run)["results"]["mixing"]
    assert mixing["batch_key"] == "sample" and not mixing["corrected"]
    assert (mixing["note"] is not None) == any(r["low"] for r in mixing["rows"])
    assert "No batch correction was applied" in report.methods_text(run)
    assert {"label": "Batch correction", "value": "off"} in payloads.settings_record(run)

    run.reduce_dimensions(integration_params=IntegrationParams(method="harmony"))
    run.cluster()
    text = report.methods_text(run)
    assert "corrected with Harmony" in text and "4 batches from \"sample\"" in text
    assert "Harmony by sample" in payloads._batch_text(run.stage_result("reduce_dimensions").summary["integration"])
    assert payloads.cluster_payload(run)["results"]["mixing"]["corrected"]


def test_single_dataset_pages_have_no_batch_parts(raw_data_path):
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import payloads

    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster"):
        getattr(run, stage)()
    red = payloads.reduction_payload(run)
    assert red["batch"] is None and red["results"]["integration"] is None
    assert payloads.cluster_payload(run)["results"]["mixing"] is None
    assert all(item["label"] != "Batch correction" for item in payloads.settings_record(run))


def test_composition_counts_each_samples_cells_per_group():
    a = ad.AnnData(np.zeros((6, 1)))
    a.obs["sample"] = ["s1", "s1", "s1", "s2", "s2", "s2"]
    a.obs["leiden"] = ["0", "0", "1", "1", "1", "1"]
    table = integration.composition(a, "sample", "leiden").set_index(["sample", "group"])
    assert table.loc[("s1", "0"), "n_cells"] == 2 and table.loc[("s1", "0"), "share"] == pytest.approx(2 / 3)
    assert table.loc[("s2", "1"), "share"] == 1.0 and ("s2", "0") not in table.index  # only groups with cells
    assert table.groupby(level="sample")["share"].sum().tolist() == pytest.approx([1.0, 1.0])
    with pytest.raises(integration.IntegrationError, match="no column"):
        integration.composition(a, "sample", "cell_type")


def test_explore_results_get_samples_and_composition(kang_piece):
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import payloads

    run = analyzed(kang_piece)
    explore = payloads.explore_payload(run)["results"]["samples"]
    adata = run.result("cluster")
    assert explore["names"] == [s["name"] for s in run.samples]
    assert len(explore["cell_sample"]) == adata.n_obs
    assert [explore["names"][i] for i in explore["cell_sample"][:5]] == adata.obs["sample"].astype(str).tolist()[:5]
    assert explore["conditions"] == ["ctrl", "stim"]
    assert {explore["condition_of"][n] for n in explore["names"]} == {"ctrl", "stim"}
    rows = explore["composition"]["cluster"]
    assert sum(r["n_cells"] for r in rows) == adata.n_obs
    assert "cell_type" not in explore["composition"]  # annotation hasn't run yet
    run.find_markers()
    run.annotate()
    explore = payloads.explore_payload(run)["results"]["samples"]
    types = explore["composition"]["cell_type"]
    assert sum(r["n_cells"] for r in types) == adata.n_obs
    by_sample = {}
    for r in types:
        by_sample[r["sample"]] = by_sample.get(r["sample"], 0) + r["share"]
    assert all(abs(v - 1) < 1e-9 for v in by_sample.values())


def test_single_dataset_explore_has_no_samples(raw_data_path):
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import payloads

    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster"):
        getattr(run, stage)()
    assert payloads.explore_payload(run)["results"]["samples"] is None


def test_report_and_export_describe_the_samples(kang_piece):
    import json
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import report
    from scrna_analysis.export import analyzed_anndata

    run = analyzed(kang_piece, IntegrationParams(method="harmony"))
    run.find_markers()
    run.annotate()
    html = report.run_report_html(run, "test")
    assert "<h2>Samples</h2>" in html and "Largest predicted cell types" in html
    kept = run.result("filter").obs["sample"].astype(str).value_counts()
    for s in run.samples:
        assert s["name"] in html and f"{kept[s['name']]:,}" in html
    out = analyzed_anndata(run)
    assert {"sample", "condition", "sample_mixing"} <= set(out.obs) and "X_pca_harmony" in out.obsm
    notes = json.loads(out.uns["orditus"])["export"]
    assert notes["obs_columns"]["sample"] == "sample" and notes["obsm"]["X_pca_harmony"]


def test_single_dataset_report_has_no_samples_section(raw_data_path):
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import report

    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess"):
        getattr(run, stage)()
    assert "<h2>Samples</h2>" not in report.run_report_html(run, "test")
