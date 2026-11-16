"""Tests for CellTypist's second opinion (scrna_analysis/celltypist_models.py) and how the annotation
stage uses it: the same answers as CellTypist called directly, its own normalization, the suggested
model, clear messages when it can't run, and the marker check's labels never changed by it."""

from __future__ import annotations

import json
import shutil

import celltypist
import numpy as np
import pytest
import scanpy as sc

from scrna_analysis import celltypist_models as ct
from scrna_analysis.params import DEFAULT_PARAMS
from scrna_analysis.pipeline import AnalysisRun

MODEL = "Immune_All_Low.pkl"


@pytest.fixture(scope="module")
def models_ready():
    for name in (MODEL, "Mouse_Whole_Brain.pkl", "Developing_Mouse_Brain.pkl"):
        if not (ct.MODEL_DIR / name).is_file():
            pytest.fail(f"Missing {ct.MODEL_DIR / name}. Run scripts/download_celltypist_models.py.")


@pytest.fixture(scope="module")
def demo_run(raw_data_path, models_ready):
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS, metadata={"species": "human", "tissue": "blood"})
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster", "find_markers"):
        getattr(run, stage)()
    return run


def test_suggested_models():
    assert ct.suggest_model("human", "blood") == MODEL and ct.suggest_model("Human", "Spleen") == MODEL
    assert ct.suggest_model("mouse", "brain") == "Mouse_Whole_Brain.pkl"  # adult suggested (analyst's choice)
    assert ct.suggest_model("human", "kidney") is None and ct.suggest_model(None, None) is None
    mouse = [c["model"] for c in ct.model_choices("mouse")]
    assert mouse == ["Mouse_Whole_Brain.pkl", "Developing_Mouse_Brain.pkl"]
    assert ct.model_choices("other") == []


def test_its_own_normalization_is_10000_counts_then_log1p(demo_run):
    x = ct.normalized_for_celltypist(demo_run.result("filter"))
    totals = np.asarray(np.expm1(x.X).sum(axis=1)).ravel()
    assert np.allclose(totals, 1e4, rtol=1e-4)
    assert demo_run.result("filter").X.max() > 100  # the counts themselves are untouched


def test_same_answers_as_celltypist_called_directly(demo_run):
    counts, clusters = demo_run.result("filter"), demo_run.result("cluster").obs["leiden"]
    per_cell, summary = ct.predict(counts, clusters, MODEL)
    x = counts.copy()
    sc.pp.normalize_total(x, target_sum=1e4)
    sc.pp.log1p(x)
    direct = celltypist.annotate(x, model=str(ct.MODEL_DIR / MODEL), majority_voting=True,
                                 over_clustering=clusters.astype(str).to_numpy()).predicted_labels
    assert (per_cell[ct.CELL_KEY].to_numpy() == direct["predicted_labels"].astype(str).to_numpy()).all()
    assert (per_cell[ct.PREDICTION_KEY].to_numpy() == direct["majority_voting"].astype(str).to_numpy()).all()
    # Per cluster: the majority label, its share of the cluster's cells, and the runner-up.
    for row in summary.clusters.itertuples():
        in_cluster = clusters.astype(str).to_numpy() == row.cluster
        shares = direct.loc[in_cluster, "predicted_labels"].astype(str).value_counts(normalize=True)
        assert row.label == direct.loc[in_cluster, "majority_voting"].astype(str).iloc[0]
        assert row.share == pytest.approx(shares[row.label])
        assert 0 <= row.mean_probability <= 1
    assert list(summary.clusters["cluster"]) == list(clusters.cat.categories)
    assert summary.model["version"] and len(summary.model["sha256"]) == 64 and "Domínguez Conde" in summary.model["citation"]
    assert 0 < summary.n_genes_matched <= summary.n_model_genes


def test_annotation_records_celltypist_and_keeps_the_marker_check_labels(demo_run):
    from dataclasses import replace
    with_ct = demo_run.params.annotation
    assert with_ct.celltypist.model == MODEL  # suggested from human blood
    demo_run.annotate(replace(with_ct, celltypist=replace(with_ct.celltypist, model=None)))
    assert "celltypist" not in demo_run.stage_result("annotate").summary  # no model chosen
    marker_labels = demo_run.result().obs["cell_type_marker_check"].astype(str)
    demo_run.annotate(replace(with_ct, celltypist=replace(with_ct.celltypist, model=MODEL)))
    s = demo_run.stage_result("annotate").summary["celltypist"]
    assert s["model"] == MODEL and len(s["clusters"]) == demo_run.result().obs["leiden"].nunique()
    out = demo_run.result()
    # The marker check's own labels are unchanged; the final label (Step B) may differ.
    assert (out.obs["cell_type_marker_check"].astype(str) == marker_labels).all()
    assert {ct.CELL_KEY, ct.PREDICTION_KEY} <= set(out.obs)
    assert demo_run.record()["params"]["annotation"]["celltypist"]["model"] == MODEL


def test_a_missing_or_changed_model_is_explained(demo_run, tmp_path):
    from dataclasses import replace
    a = demo_run.params.annotation
    demo_run.annotate(replace(a, celltypist=replace(a.celltypist, model="Not_A_Model.pkl")))
    s = demo_run.stage_result("annotate").summary
    assert "isn't downloaded" in s["celltypist"]["error"] and s["clusters"]  # the marker check still ran
    shutil.copy(ct.MODEL_DIR / MODEL, tmp_path / MODEL)
    (tmp_path / ct.LISTING).write_text(json.dumps({MODEL: {"sha256": "0" * 64}}))
    with pytest.raises(ct.CellTypistError, match="fingerprint"):
        ct.load_model(MODEL, tmp_path)


def test_mouse_data_with_a_human_model_warns_or_stops(mouse_10x_path, models_ready):
    run = AnalysisRun(mouse_10x_path, DEFAULT_PARAMS, metadata={"species": "mouse", "tissue": "brain"})
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster"):
        getattr(run, stage)()
    try:
        _, summary = ct.predict(run.result("filter"), run.result("cluster").obs["leiden"], MODEL)
        assert summary.gene_match < ct.LOW_GENE_MATCH and "less reliable" in summary.notes[0]
    except ct.CellTypistError as problem:
        assert "None of the genes" in str(problem)
    _, summary = ct.predict(run.result("filter"), run.result("cluster").obs["leiden"], "Developing_Mouse_Brain.pkl")
    assert summary.gene_match >= ct.LOW_GENE_MATCH and not summary.notes


def test_page_data_and_report(raw_data_path, models_ready):
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import payloads
    import report

    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS, metadata={"species": "human", "tissue": "blood"})
    info = payloads.annotation_payload(run)["celltypist"]
    assert info["suggested"] == MODEL and info["current"] == MODEL
    assert [c["model"] for c in info["choices"]] == [MODEL] and info["choices"][0]["downloaded"]
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster",
                  "find_markers", "annotate"):
        getattr(run, stage)()
    results = payloads.annotation_payload(run)["results"]
    assert results["celltypist"]["model"] == MODEL and len(results["celltypist"]["clusters"]) == len(results["clusters"])
    assert "CellTypist" in report.methods_text(run) and "broad cell types" in report.methods_text(run)
    html = report.run_report_html(run, "test")
    assert "CellTypist (share of votes)" in html and "Domínguez Conde" in html
    # The page's choice: none, or one of the offered models.
    (a,) = payloads.annotation_params_from_browser({**_boxes(run), "celltypist_model": ""}, run.params)
    assert a.celltypist.model is None
    with pytest.raises(payloads.InputError, match="CellTypist models offered"):
        payloads.annotation_params_from_browser({**_boxes(run), "celltypist_model": "Mouse_Whole_Brain.pkl"}, run.params)


def _boxes(run):
    import payloads
    s = payloads.annotation_settings(run.params)
    return {k: str(s[k]) for k in ("top_n_genes", "max_adj_p", "min_overlap", "ambiguity_fold")}
