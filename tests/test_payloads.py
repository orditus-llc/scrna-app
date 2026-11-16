"""Tests for `app/payloads.py`: the JSON the React interface receives.

Checks that every number shown on the Workflow page comes from the run (compared with the
reference files), that step status follows the pipeline (including the 8 stages -> 7 steps
mapping and "needs re-running" after a redo), and that the payload is plain JSON.
"""

from __future__ import annotations

import dataclasses
from dataclasses import replace
import json
import re
import sys

import pytest

from scrna_analysis import marker_db
from scrna_analysis.params import LeidenParams
from scrna_analysis.params import DEFAULT_PARAMS
from tests.reference_preset import REFERENCE_PRESET
from scrna_analysis.pipeline import AnalysisRun, file_sha256, run_analysis
from tests.conftest import PROJECT_ROOT

APP_DIR = PROJECT_ROOT / "app"
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import payloads  # noqa: E402

# regress_out triggers harmless Numba performance hints (also seen in the notebook).
pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")

STEP_KEYS = [key for key, _, _ in payloads.STEPS]


def statuses(payload) -> dict[str, str]:
    return {s["key"]: s["status"] for s in payload["steps"]}


def step(payload, key) -> dict:
    return next(s for s in payload["steps"] if s["key"] == key)


# --- Before anything runs -----------------------------------------------------------------

def test_fresh_run(raw_data_path):
    payload = payloads.workflow_payload(AnalysisRun(raw_data_path, REFERENCE_PRESET))
    assert [s["key"] for s in payload["steps"]] == STEP_KEYS
    assert [s["number"] for s in payload["steps"]] == list(range(1, 8))
    assert set(statuses(payload).values()) == {"not_run"}
    assert all(s["summary"] == payloads.DESCRIPTIONS[s["key"]] for s in payload["steps"])
    assert payload["run_log"] == [] and payload["notice"] is None
    assert payload["dataset"] == {"file": raw_data_path.name, "label": raw_data_path.name, "n_cells": None,
                                  "n_genes": None, "fingerprint": None, "file_format": None, "notes": [],
                                  "mito_matches": [], "metadata": None, "samples": None}


def test_settings_record_comes_from_the_settings(raw_data_path):
    record = {item["label"]: item["value"] for item in payloads.workflow_payload(AnalysisRun(raw_data_path, REFERENCE_PRESET))["settings_record"]}
    f, P = REFERENCE_PRESET.filter, REFERENCE_PRESET
    assert record["Min genes per cell"] == f"{f.min_genes:,} (inclusive)"
    assert record["Max genes per cell"] == f"< {f.max_n_genes_by_counts:,}"
    assert record["Max % mitochondrial"] == f"< {f.max_pct_counts_mt:g}%"
    assert record["Principal components"] == f"{P.pca.n_comps} computed, {P.neighbors.n_pcs} used"
    assert record["Random seed"] == str(P.pca.random_state)
    assert record["Leiden resolution"] == repr(float(P.leiden.resolution))  # keeps its decimal, e.g. "1.0"


# --- After Load ---------------------------------------------------------------------------

def test_after_load(raw_data_path, reference_summary):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load()
    payload = payloads.workflow_payload(run)
    n_cells, n_genes = reference_summary["raw_shape"]
    fingerprint = file_sha256(raw_data_path)[:12]

    assert statuses(payload)["data_overview"] == "done"
    assert [v for k, v in statuses(payload).items() if k != "data_overview"] == ["not_run"] * 6
    assert step(payload, "data_overview")["summary"] == f"{n_cells:,} cells × {n_genes:,} genes · raw counts"
    assert payload["dataset"] == {"file": raw_data_path.name, "label": raw_data_path.name, "n_cells": n_cells,
                                  "n_genes": n_genes, "fingerprint": fingerprint,
                                  "file_format": "AnnData (.h5ad)", "notes": [],
                                  "mito_matches": payload["dataset"]["mito_matches"],
                                  "metadata": payload["dataset"]["metadata"], "samples": payload["dataset"]["samples"]}
    # A single file: no samples, and (PBMC3k has no text columns) nothing to choose a sample column from.
    assert payload["dataset"]["samples"]["list"] == [] and payload["dataset"]["samples"]["from_files"] is False
    meta = payload["dataset"]["metadata"]
    assert meta["confirmed"] is False and meta["suggested_species"] == "human"  # from the "MT-" genes
    top = payload["dataset"]["mito_matches"][0]
    assert top["prefix"] == "MT-" and top["n_genes"] == reference_summary["n_mito_genes"]
    assert [entry["label"] for entry in payload["run_log"]] == ["Data Overview"]


def test_quality_control_needs_both_stages(raw_data_path, reference_summary):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load()
    run.compute_qc_metrics()
    qc = step(payloads.workflow_payload(run), "quality_control")
    assert qc["status"] == "not_run"  # metrics alone aren't the whole step
    assert qc["summary"].startswith(f"QC metrics calculated ({reference_summary['n_mito_genes']} mitochondrial genes)")


# --- Full run and redo --------------------------------------------------------------------

@pytest.fixture(scope="module")
def full_run(raw_data_path):
    return run_analysis(raw_data_path, REFERENCE_PRESET)


def test_full_run_summaries_come_from_results(full_run, reference_summary):
    payload = payloads.workflow_payload(full_run)
    assert set(statuses(payload).values()) == {"done"}
    cells, genes = reference_summary["filtered_shape"]
    raw_cells, raw_genes = reference_summary["raw_shape"]
    assert step(payload, "quality_control")["summary"] == (
        f"{cells:,} cells × {genes:,} genes kept ({raw_cells - cells:,} cells removed)")
    assert step(payload, "preprocessing")["summary"].startswith(f"{reference_summary['n_hvg']:,} highly variable genes")
    assert step(payload, "clustering")["summary"].startswith(f"{reference_summary['n_clusters']} clusters")
    n = reference_summary["n_clusters"]
    assert re.match(rf"\d+ of {n} clusters named · \d+ Unknown$", step(payload, "cell_type_annotation")["summary"])
    assert [e["label"] for e in payload["run_log"]] == [payloads.STAGE_LABELS[s] for s in
                                                       ("load", "qc_metrics", "filter", "preprocess",
                                                        "reduce_dimensions", "cluster", "markers", "annotate")]
    assert json.loads(json.dumps(payload)) == payload  # plain JSON


def test_redo_marks_later_steps_with_a_reason(raw_data_path):
    run = run_analysis(raw_data_path, REFERENCE_PRESET)
    run.cluster(LeidenParams(resolution=0.5))
    payload = payloads.workflow_payload(run)
    assert statuses(payload)["clustering"] == "done"
    for key in ("marker_genes", "cell_type_annotation"):
        assert step(payload, key)["status"] == "needs_rerun"
        assert step(payload, key)["reason"].startswith("Clustering was redone")
    record = {item["label"]: item["value"] for item in payload["settings_record"]}
    assert record["Leiden resolution"] == "0.5"
    assert payload["run_log"][-1]["text"].endswith("cleared Marker Genes, Cell Type Annotation")


# --- Guidance: next step, confirmations, prerequisites ------------------------------------

def test_next_step_moves_along_the_workflow(raw_data_path, full_run):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    assert payloads.workflow_payload(run)["next_step"]["key"] == "data_overview"
    run.load()
    assert payloads.workflow_payload(run)["next_step"] == {
        "key": "quality_control", "title": "Quality Control", "text": payloads.NEXT_STEP_TEXT["quality_control"]}
    run.compute_qc_metrics()  # QC metrics alone don't finish the Quality Control step
    assert payloads.workflow_payload(run)["next_step"]["key"] == "quality_control"
    assert payloads.workflow_payload(full_run)["next_step"] is None


def test_next_step_points_at_a_step_that_needs_rerunning(raw_data_path):
    run = run_analysis(raw_data_path, REFERENCE_PRESET)
    run.cluster(LeidenParams(resolution=0.5))
    assert payloads.workflow_payload(run)["next_step"]["key"] == "marker_genes"


def test_load_confirmation_uses_the_loaded_numbers(raw_data_path, reference_summary):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    assert step(payloads.workflow_payload(run), "data_overview")["confirmation"] is None
    run.load()
    n_cells, n_genes = reference_summary["raw_shape"]
    assert step(payloads.workflow_payload(run), "data_overview")["confirmation"] == (
        f"Data loaded: {n_cells:,} cells × {n_genes:,} genes of raw counts. Ready for Quality Control.")
    with_metadata = AnalysisRun(raw_data_path, REFERENCE_PRESET, metadata={"species": "human", "tissue": "blood"})
    with_metadata.load()
    assert "genes of raw counts (Human · Blood)." in step(payloads.workflow_payload(with_metadata), "data_overview")["confirmation"]


def test_prerequisite_names_the_step_that_comes_first(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    payload = payloads.workflow_payload(run)
    assert step(payload, "data_overview")["prerequisite"] is None
    assert step(payload, "preprocessing")["prerequisite"] == "Data Overview"
    run.load()
    payload = payloads.workflow_payload(run)
    assert step(payload, "data_overview")["prerequisite"] is None  # done
    assert step(payload, "quality_control")["prerequisite"] is None  # can be started now
    assert step(payload, "clustering")["prerequisite"] == "Quality Control"


# --- Quality Control page -----------------------------------------------------------------

@pytest.fixture(scope="module")
def qc_run(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load()
    run.compute_qc_metrics()
    return run


REFERENCE_BOXES = {"min_genes": "200", "min_cells": "3", "max_genes": "2500", "max_genes_off": False,
                   "max_mt": "6", "max_mt_off": False}


def test_qc_payload_before_qc(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    payload = payloads.qc_payload(run)
    assert payload["metrics_done"] is False and payload["can_run"] is False and payload["applied"] is None
    run.load()
    assert payloads.qc_payload(run)["can_run"] is True
    assert payloads.preview_payload(run, None) is None  # nothing to preview before QC


def test_qc_payload_matches_reference(qc_run, reference_summary, reference_results):
    payload = payloads.qc_payload(qc_run)
    assert payload["metrics_done"] is True
    assert [payload["n_cells"], payload["n_genes"]] == reference_summary["raw_shape"]
    assert payload["mito"]["genes"] == reference_summary["mito_genes"]
    # The summary table equals the notebook's describe() table.
    names = dict(payloads.QC_STATISTICS)
    for row in payload["stats"]:
        statistic = next(k for k, v in names.items() if v == row["statistic"])
        for metric, stats in reference_summary["qc_describe"].items():
            assert row[metric] == pytest.approx(stats[statistic], rel=0, abs=1e-6), (metric, statistic)
    # Per-cell values (sent to the plots) equal the reference for every kept cell.
    cells = payload["cells"]
    position = {barcode: i for i, barcode in enumerate(cells["barcode"])}
    obs = reference_results.obs
    idx = [position[b] for b in obs.index]
    assert [cells["n_genes"][i] for i in idx] == obs["n_genes_by_counts"].astype(int).tolist()
    assert [cells["total_counts"][i] for i in idx] == [float(v) for v in obs["total_counts"]]
    assert [cells["pct_mt"][i] for i in idx] == [float(v) for v in obs["pct_counts_mt"]]
    # A reference run keeps the notebook's fixed thresholds; the page also gets the data's suggestion.
    assert payload["filter"]["current"] == dataclasses.asdict(REFERENCE_PRESET.filter)
    assert payload["filter"]["suggested"] is not None and payload["filter"]["details"]["min_genes"]["n_mads"] == 5


def test_preview_removes_exactly_the_reference_cells(qc_run, reference_results, reference_summary):
    preview = payloads.preview_payload(qc_run, None)  # nothing sent yet = current settings
    barcodes = payloads.qc_payload(qc_run)["cells"]["barcode"]
    assert set(preview["removed_barcodes"]) == set(barcodes) - set(reference_results.obs_names)
    assert set(preview["failures"]) == set(preview["removed_barcodes"])
    assert all(preview["failures"][b] for b in preview["removed_barcodes"])  # each fails at least one filter
    assert [preview["n_cells_kept"], preview["n_genes_kept"]] == reference_summary["filtered_shape"]
    assert payloads.preview_payload(qc_run, REFERENCE_BOXES) == preview  # the boxes say the same


def test_preview_follows_other_thresholds(qc_run):
    preview = payloads.preview_payload(qc_run, {**REFERENCE_BOXES, "max_mt": "5"})
    expected = qc_run.preview_filter(dataclasses.replace(REFERENCE_PRESET.filter, max_pct_counts_mt=5.0))
    assert preview["n_cells_kept"] == expected.n_cells_kept
    assert preview["settings"]["max_pct_counts_mt"] == 5.0
    unlimited = payloads.preview_payload(qc_run, {**REFERENCE_BOXES, "max_mt_off": True, "max_genes_off": True})
    assert unlimited["settings"]["max_pct_counts_mt"] is None and unlimited["n_cells_kept"] > expected.n_cells_kept


@pytest.mark.parametrize("change, message", [
    ({"min_genes": ""}, "Enter a value for 'Min genes per cell'"),
    ({"min_cells": "-1"}, "can't be negative"),
    ({"max_genes": "2500.5"}, "must be a whole number"),
    ({"min_genes": "100000"}, "would remove all cells"),
])
def test_preview_explains_unusable_settings(qc_run, change, message):
    preview = payloads.preview_payload(qc_run, {**REFERENCE_BOXES, **change})
    assert preview["ok"] is False and message in preview["message"]


def test_after_applying_the_filter(raw_data_path, reference_summary):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load()
    run.compute_qc_metrics()
    run.apply_filter(payloads.filter_params_from_browser(REFERENCE_BOXES, run.params.filter))
    cells, genes = reference_summary["filtered_shape"]
    raw_cells, raw_genes = reference_summary["raw_shape"]
    confirmation = f"Filtering applied: kept {cells:,} of {raw_cells:,} cells and {genes:,} of {raw_genes:,} genes."
    workflow = payloads.workflow_payload(run)
    assert step(workflow, "quality_control")["status"] == "done"
    assert step(workflow, "quality_control")["confirmation"] == confirmation
    assert workflow["next_step"]["key"] == "preprocessing"
    assert payloads.qc_payload(run)["applied"]["confirmation"] == confirmation


def test_notice_is_passed_through(raw_data_path):
    notice = {"id": 3, "kind": "error", "text": "Something to fix"}
    assert payloads.workflow_payload(AnalysisRun(raw_data_path, REFERENCE_PRESET), notice)["notice"] == notice


# --- Preprocessing and Dimensionality Reduction cards -------------------------------------------

def as_boxes(settings: dict) -> dict:
    """Settings as the page's boxes send them (numbers as typed text, switches as booleans)."""
    return {k: v if isinstance(v, bool) else str(v) for k, v in settings.items()}


def test_cards_before_running(qc_run):
    pre, red = payloads.preprocess_payload(qc_run), payloads.reduction_payload(qc_run)
    assert pre["status"] == red["status"] == "not_run"
    assert pre["results"] is None and red["results"] is None
    # The run uses the notebook's settings; the page compares them with the app's defaults.
    assert pre["settings"]["current"] == payloads.preprocess_settings(REFERENCE_PRESET)
    assert red["settings"]["current"] == payloads.reduction_settings(REFERENCE_PRESET)
    assert pre["settings"]["reference"] == payloads.preprocess_settings(DEFAULT_PARAMS)
    assert red["settings"]["reference"] == payloads.reduction_settings(DEFAULT_PARAMS)
    assert pre["fixed"] == {"hvg_flavor": REFERENCE_PRESET.hvg.flavor,
                            "regress_keys": list(REFERENCE_PRESET.regress_out.keys)}


def test_preprocess_card_matches_reference(full_run, reference_summary):
    payload = payloads.preprocess_payload(full_run)
    json.dumps(payload, allow_nan=False)  # plain JSON (NaN dispersions become null)
    results = payload["results"]
    genes = results["genes"]
    assert results["n_hvg"] == reference_summary["n_hvg"]
    assert results["n_genes"] == len(genes["name"]) == reference_summary["filtered_shape"][1]
    selected = {name for name, hvg in zip(genes["name"], genes["hvg"]) if hvg}
    assert selected == set(reference_summary["hvg_genes"])
    assert f"{reference_summary['n_hvg']:,} of {reference_summary['filtered_shape'][1]:,} genes" in \
        step(payloads.workflow_payload(full_run), "preprocessing")["confirmation"]


def test_reduction_card_matches_reference(full_run, reference_summary, reference_results):
    payload = payloads.reduction_payload(full_run)
    json.dumps(payload, allow_nan=False)
    results = payload["results"]
    assert results["variance_ratio"] == pytest.approx(reference_summary["pca_variance_ratio"], abs=1e-6)
    assert results["cells"]["barcode"] == reference_results.obs_names.tolist()
    expected = reference_results.obsm["X_umap"]
    assert results["umap"]["x"] == pytest.approx(expected[:, 0].tolist(), abs=1e-6)
    assert results["umap"]["y"] == pytest.approx(expected[:, 1].tolist(), abs=1e-6)
    assert results["n_pcs_used"] == REFERENCE_PRESET.neighbors.n_pcs
    confirmation = step(payloads.workflow_payload(full_run), "dimensionality_reduction")["confirmation"]
    assert f"UMAP drawn for {reference_summary['filtered_shape'][0]:,} cells" in confirmation


def test_cards_say_what_rerunning_would_clear(full_run, qc_run):
    titles = [title for _, title, _ in payloads.STEPS]
    assert payloads.preprocess_payload(full_run)["will_clear"] == titles[3:]
    assert payloads.reduction_payload(full_run)["will_clear"] == titles[4:]
    assert payloads.preprocess_payload(qc_run)["will_clear"] == []


def test_reference_boxes_give_the_reference_settings():
    P = REFERENCE_PRESET
    pre = payloads.preprocess_params_from_browser(as_boxes(payloads.preprocess_settings(P)), P)
    assert pre == (P.normalize, P.hvg, P.regress_out, P.scale)
    red = payloads.reduction_params_from_browser(as_boxes(payloads.reduction_settings(P)), P)
    assert red == (P.pca, P.neighbors, P.umap, P.integration)  # batch correction off, as in the preset


def test_switching_regression_off_keeps_everything_else():
    P = REFERENCE_PRESET
    boxes = {**as_boxes(payloads.preprocess_settings(P)), "regress": False}
    normalize, hvg, regress, scale = payloads.preprocess_params_from_browser(boxes, P)
    assert regress.enabled is False and regress.keys == P.regress_out.keys
    assert (normalize, hvg, scale) == (P.normalize, P.hvg, P.scale)


def test_one_seed_for_pca_neighbors_and_umap():
    P = REFERENCE_PRESET
    boxes = {**as_boxes(payloads.reduction_settings(P)), "seed": "7"}
    pca, neighbors, umap, integ = payloads.reduction_params_from_browser(boxes, P)
    assert pca.random_state == neighbors.random_state == umap.random_state == integ.random_state == 7
    new_params = dataclasses.replace(P, pca=pca, neighbors=neighbors, umap=umap)
    record = {i["label"]: i["value"] for i in payloads.settings_record(AnalysisRun("unused.h5ad", params=new_params))}
    assert record["Random seed"] == f"PCA 7, neighbors 7, UMAP 7, Leiden {P.leiden.random_state}"


@pytest.mark.parametrize("reader, settings, change, message", [
    ("preprocess", payloads.preprocess_settings, {"target_sum": ""}, "Enter a value for 'Counts per cell after normalization'."),
    ("preprocess", payloads.preprocess_settings, {"target_sum": "0"}, "'Counts per cell after normalization' must be greater than 0."),
    ("preprocess", payloads.preprocess_settings, {"min_mean": "3", "max_mean": "3"}, "'Min mean expression' must be smaller than 'Max mean expression'."),
    ("preprocess", payloads.preprocess_settings, {"max_value": "abc"}, "'Clip scaled values at' must be a number."),
    ("reduction", payloads.reduction_settings, {"n_comps": "20", "n_pcs": "30"}, "'PCs used for the neighbor graph' can't be more than 'Principal components'."),
    ("reduction", payloads.reduction_settings, {"n_neighbors": "1"}, "'Neighbors per cell' must be at least 2."),
    ("reduction", payloads.reduction_settings, {"seed": "1.5"}, "'Random seed' must be a whole number."),
    ("reduction", payloads.reduction_settings, {"seed": "-1"}, "'Random seed' can't be negative."),
])
def test_unusable_card_values_are_explained(reader, settings, change, message):
    read = {"preprocess": payloads.preprocess_params_from_browser,
            "reduction": payloads.reduction_params_from_browser}[reader]
    with pytest.raises(payloads.InputError) as problem:
        read({**as_boxes(settings(REFERENCE_PRESET)), **change}, REFERENCE_PRESET)
    assert str(problem.value) == message


def test_negative_dispersion_cutoff_is_allowed():
    """Normalized dispersion is a z-score, so a negative minimum is a valid choice."""
    boxes = {**as_boxes(payloads.preprocess_settings(REFERENCE_PRESET)), "min_disp": "-0.5"}
    assert payloads.preprocess_params_from_browser(boxes, REFERENCE_PRESET)[1].min_disp == -0.5


def test_qc_page_says_what_reapplying_would_clear(full_run, qc_run):
    titles = [title for _, title, _ in payloads.STEPS]
    assert payloads.qc_payload(full_run)["will_clear"] == titles[2:]
    assert payloads.qc_payload(qc_run)["will_clear"] == []


# --- Quality Control: comparison with the reference -------------------------------------------

REFERENCE_FILTER_BOXES = {"min_genes": "200", "min_cells": "3", "max_genes": "2500", "max_genes_off": False,
                          "max_mt": "6", "max_mt_off": False}


def suggested_boxes(run) -> dict:
    s = payloads.suggestion_payload(run, None)["suggested"]
    return {"min_genes": str(s["min_genes"]), "min_cells": str(s["min_cells"]),
            "max_genes": str(s["max_n_genes_by_counts"]), "max_genes_off": False,
            "max_mt": str(s["max_pct_counts_mt"]), "max_mt_off": False}


def test_preview_compares_with_the_suggestion(qc_run):
    """At the suggested thresholds the preview equals the suggestion's own preview."""
    preview = payloads.preview_payload(qc_run, suggested_boxes(qc_run))
    json.dumps(preview, allow_nan=False)
    assert preview["suggested"] == {"n_cells_kept": preview["n_cells_kept"], "n_genes_kept": preview["n_genes_kept"]}
    assert not any(row["differs"] for row in preview["filters"])
    assert all(row["removed"] == row["suggested_removed"] for row in preview["filters"])


def test_preview_at_the_reference_thresholds_adds_up(qc_run, reference_summary):
    """At the notebook's thresholds the preview keeps the reference cells and genes; the
    per-filter table adds up to them."""
    cells, genes = reference_summary["filtered_shape"]
    preview = payloads.preview_payload(qc_run, {**REFERENCE_FILTER_BOXES,
                                                "max_mt": str(REFERENCE_PRESET.filter.max_pct_counts_mt)})
    assert [preview["n_cells_kept"], preview["n_genes_kept"]] == [cells, genes]
    removed = {what: sum(r["removed"] for r in preview["filters"] if r["what"] == what) for what in ("cells", "genes")}
    assert preview["n_cells_before"] - removed["cells"] == cells
    assert preview["n_genes_before"] - removed["genes"] == genes


def test_filter_table_marks_the_setting_that_differs(qc_run):
    boxes = suggested_boxes(qc_run)
    suggested_mt = float(boxes["max_mt"])
    stricter = suggested_mt - 1
    preview = payloads.preview_payload(qc_run, {**boxes, "max_mt": str(stricter)})
    rows = {row["key"]: row for row in preview["filters"]}
    assert [key for key, row in rows.items() if row["differs"]] == ["max_pct_counts_mt"]
    mito = rows["max_pct_counts_mt"]
    assert mito["setting"] == stricter and mito["suggested_setting"] == suggested_mt
    assert mito["removed"] > mito["suggested_removed"]
    # The suggested numbers don't move when your settings change.
    at_suggestion = payloads.preview_payload(qc_run, boxes)["n_cells_kept"]
    assert preview["suggested"]["n_cells_kept"] == at_suggestion
    assert preview["n_cells_kept"] == at_suggestion - (mito["removed"] - mito["suggested_removed"])


def test_switched_off_filter_is_shown_as_not_run(qc_run):
    preview = payloads.preview_payload(qc_run, {**REFERENCE_FILTER_BOXES, "max_genes_off": True})
    row = next(r for r in preview["filters"] if r["key"] == "max_n_genes_by_counts")
    assert row["setting"] is None and row["removed"] is None and row["differs"] is True
    assert row["suggested_removed"] is not None


# --- Clustering page ---------------------------------------------------------------------------

def test_clustering_page_matches_reference(full_run, reference_summary, reference_results):
    from sklearn.metrics import adjusted_rand_score

    payload = payloads.cluster_payload(full_run)
    json.dumps(payload, allow_nan=False)
    r = payload["results"]
    assert payload["status"] == "done" and r["n_clusters"] == reference_summary["n_clusters"]
    assert sorted(c["n_cells"] for c in r["clusters"]) == sorted(reference_summary["cluster_sizes"].values())
    assert r["cells"]["barcode"] == reference_results.obs_names.tolist()
    assert adjusted_rand_score(reference_results.obs["leiden"].astype(str), r["labels"]) == 1.0
    assert sum(c["n_cells"] for c in r["clusters"]) == len(r["labels"]) == len(r["umap"]["x"])
    assert r["settings"] == payloads.cluster_settings(REFERENCE_PRESET)
    confirmation = step(payloads.workflow_payload(full_run), "clustering")["confirmation"]
    assert confirmation.startswith(f"Clustering done: {reference_summary['n_clusters']} clusters")


def test_cluster_boxes_give_the_settings():
    P = REFERENCE_PRESET
    (leiden,) = payloads.cluster_params_from_browser(as_boxes(payloads.cluster_settings(P)), P)
    assert leiden == P.leiden
    (leiden,) = payloads.cluster_params_from_browser({"resolution": "0.5", "seed": "3"}, P)
    assert (leiden.resolution, leiden.random_state) == (0.5, 3)
    assert (leiden.flavor, leiden.n_iterations, leiden.directed) == (P.leiden.flavor, P.leiden.n_iterations, P.leiden.directed)


@pytest.mark.parametrize("change, message", [
    ({"resolution": "0"}, "'Leiden resolution' must be greater than 0."),
    ({"resolution": ""}, "Enter a value for 'Leiden resolution'."),
    ({"seed": "2.5"}, "'Leiden random seed' must be a whole number."),
])
def test_unusable_cluster_values_are_explained(change, message):
    with pytest.raises(payloads.InputError) as problem:
        payloads.cluster_params_from_browser({**as_boxes(payloads.cluster_settings(REFERENCE_PRESET)), **change},
                                             REFERENCE_PRESET)
    assert str(problem.value) == message


def test_lower_resolution_gives_fewer_clusters(raw_data_path, reference_summary):
    """The number of clusters comes from the data at the chosen resolution, not from the reference."""
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    for stage in (run.load, run.compute_qc_metrics, run.apply_filter, run.preprocess, run.reduce_dimensions):
        stage()
    run.cluster(LeidenParams(resolution=0.3))
    payload = payloads.cluster_payload(run)
    assert payload["results"]["n_clusters"] < reference_summary["n_clusters"]
    assert payload["results"]["settings"]["resolution"] == 0.3


# --- Marker Genes page -------------------------------------------------------------------------

def test_marker_page_matches_reference(full_run, reference_summary):
    payload = payloads.markers_payload(full_run)
    json.dumps(payload, allow_nan=False)
    r = payload["results"]
    assert r["n_genes_tested"] == reference_summary["filtered_shape"][1]
    assert [g["cluster"] for g in r["groups"]] == list(reference_summary["top25_markers"])
    for cluster, genes in reference_summary["top25_markers"].items():
        assert [row["gene"] for row in r["top"][cluster][:25]] == genes
        assert [row["rank"] for row in r["top"][cluster][:25]] == list(range(1, 26))
    assert sum(g["n_cells"] for g in r["groups"]) == reference_summary["filtered_shape"][0]


def test_marker_dotplot_uses_each_clusters_top_genes(full_run, reference_summary):
    dot = payloads.markers_payload(full_run)["results"]["dotplot"]
    clusters = list(reference_summary["top25_markers"])
    assert len(dot["genes"]) == len(set(dot["genes"]))  # no gene twice
    for cluster, genes in reference_summary["top25_markers"].items():
        assert genes[0] in dot["genes"]  # every cluster's top marker is shown
    assert len(dot["values"]) == len(dot["genes"]) * len(clusters)  # every gene in every cluster
    assert all(0 <= v["frac_expressing"] <= 1 and v["mean_expr"] >= 0 for v in dot["values"])


def test_marker_heatmap_has_every_option_and_its_values(full_run, reference_summary):
    heat = payloads.markers_payload(full_run)["results"]["heatmap"]
    clusters = list(reference_summary["top25_markers"])
    assert heat["default"] in heat["options"]
    values = {(v["group"], v["gene"]): v for v in heat["values"]}
    for n in heat["options"]:
        rows = heat["rows"][str(n)]
        assert [r["cluster"] for r in rows] == [c for c in clusters for _ in range(n)]  # n per cluster, in order
        assert len({r["gene"] for r in rows}) == len(rows)  # no gene twice
        for cluster in clusters:  # every cluster's top marker heads its rows (if no earlier cluster took it)
            top = reference_summary["top25_markers"][cluster][0]
            assert any(r["gene"] == top for r in rows)
        # Every row's gene has a value in every cluster.
        assert all((c, r["gene"]) in values for r in rows for c in clusters)
    for gene in {g for _, g in values}:
        scaled = [values[(c, gene)]["scaled"] for c in clusters]
        assert min(scaled) == 0 and max(scaled) in (0, 1)  # scaled within each gene


def test_marker_test_choice():
    P = REFERENCE_PRESET
    assert payloads.marker_params_from_browser({"method": "wilcoxon"}, P) == (P.markers,)
    (m,) = payloads.marker_params_from_browser({"method": "t-test"}, P)
    assert m.method == "t-test" and (m.use_raw, m.reference, m.corr_method) == (P.markers.use_raw, P.markers.reference, P.markers.corr_method)
    with pytest.raises(payloads.InputError, match="Choose one of the tests"):
        payloads.marker_params_from_browser({"method": "magic"}, P)


# --- Cell Type Annotation page -------------------------------------------------------------------

def test_annotation_page_shows_the_evidence(full_run, reference_results):
    payload = payloads.annotation_payload(full_run)
    json.dumps(payload, allow_nan=False)
    assert "database" not in payload and payload["selection"] == {"min_markers": 3}  # no PanglaoDB
    assert payload["census"] == {"tissue": "Blood", "release": "2025-11-08", "has_table": True, "available": True}
    assert payload["metadata"] == {"species": "Human", "tissue": "Blood", "confirmed_by": None, "has_markers": True}
    r = payload["results"]
    assert r["database"]["name"] == "CELLxGENE Census" and r["database"]["tissue"] == "Blood"
    assert r["no_reference"] is False
    assert not {"references", "first_reference", "answers", "candidates_by_reference"} & set(r)
    assert r["cells"]["barcode"] == reference_results.obs_names.tolist()
    assert len(r["cell_label"]) == len(r["cell_cluster"]) == len(r["umap"]["x"])
    assert sum(r["label_counts"].values()) == len(r["cell_label"])
    clusters = {c["cluster"]: c for c in r["clusters"]}
    assert set(r["candidates"]) == set(clusters)
    for cluster, candidates in r["candidates"].items():
        assert 1 <= len(candidates) <= payloads.CANDIDATES_SHOWN
        assert [c["adj_p"] for c in candidates] == sorted(c["adj_p"] for c in candidates)
        assert all(c["adj_p"] > 0 for c in candidates)  # tiny p-values aren't rounded to 0
        assert candidates[0]["cell_type"] == clusters[cluster]["best_candidate"]  # best first
        assert clusters[cluster]["predicted_label"] in (clusters[cluster]["best_candidate"], "Unresolved")
    # The dot plot shows matching markers of each cluster's best candidate.
    matching = {g for c in r["clusters"] for g in (c["matching_markers"] or "").split(", ") if g}
    assert {g["gene"] for g in r["dotplot"]["genes"]} <= matching
    assert len(r["dotplot"]["values"]) == len(r["dotplot"]["genes"]) * len(clusters)
    confirmation = step(payloads.workflow_payload(full_run), "cell_type_annotation")["confirmation"]
    assert "CELLxGENE Census" in confirmation and "predictions for you to review" in confirmation


def test_annotation_settings_from_the_page():
    """The page sends only the checks (and the CellTypist model); species and tissue come from the
    dataset's metadata, so any the browser sends are ignored."""
    P = REFERENCE_PRESET
    boxes = {"top_n_genes": "100", "max_adj_p": "0.05", "min_overlap": "3", "ambiguity_fold": "10"}
    assert payloads.annotation_params_from_browser(boxes, P) == (P.annotation,)
    (a,) = payloads.annotation_params_from_browser({**boxes, "species": "mouse", "tissue": "brain",
                                                    "use_census": False, "min_overlap": "4"}, P)
    assert (a.markers, a.census, a.test.min_overlap) == (P.annotation.markers, P.annotation.census, 4)
    for change, message in [({"top_n_genes": "0"}, "must be at least 1"), ({"max_adj_p": "2"}, "can't be more than 1"),
                            ({"max_adj_p": "0"}, "must be greater than 0"), ({"ambiguity_fold": "0.5"}, "at least 1")]:
        with pytest.raises(payloads.InputError, match=message):
            payloads.annotation_params_from_browser({**boxes, **change}, P)


@pytest.fixture(scope="module")
def spare_run(raw_data_path):
    """A second full analysis that tests may re-annotate (so `full_run` stays as it is)."""
    return run_analysis(raw_data_path, REFERENCE_PRESET)


def test_stricter_checks_can_leave_clusters_unresolved(spare_run):
    """Labels come from the evidence and the checks: an extreme significance level leaves every
    cluster Unresolved instead of forcing a label."""
    run = spare_run
    run.annotate(replace(run.params.annotation, test=replace(run.params.annotation.test, max_adj_p=1e-300)))
    r = payloads.annotation_payload(run)["results"]
    assert all(c["predicted_label"] == "Unresolved" and "not_significant" in c["reasons"] for c in r["clusters"])
    assert set(r["cell_label"]) == {"Unknown"}  # the final label: no CellTypist model in this preset either
    assert "0 of" in step(payloads.workflow_payload(run), "cell_type_annotation")["confirmation"]


def test_a_tissue_without_a_census_reference_is_shown_on_the_page(raw_data_path, tmp_path):
    """No Census table for the tissue: the page says the marker check won't run, before and after."""
    from scrna_analysis.pipeline import AnalysisRun
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET, census_dir=tmp_path)  # no tables in this folder
    assert payloads.annotation_payload(run)["census"]["has_table"] is False
    for stage in payloads.remaining_stages(run):
        payloads.run_stage(run, stage)
    r = payloads.annotation_payload(run)["results"]
    assert r["database"] is None and r["no_reference"] is True and r["candidates"] == {}
    assert {c["status"] for c in r["clusters"]} == {"no_reference"} and r["unresolved_note"] is None
    assert all("no CELLxGENE Census reference for blood yet" in g for g in r["guidance"].values())
    assert "Review the clusters with their marker genes" in r["guidance"]["0"] or "top marker genes" in r["guidance"]["0"]
    confirmation = step(payloads.workflow_payload(run), "cell_type_annotation")["confirmation"]
    assert confirmation.startswith("Cell Type Annotation done: there's no CELLxGENE Census reference for blood yet")
    assert "predictions for you to review" not in confirmation  # nothing was predicted


def test_annotation_page_has_no_marker_list_upload(full_run):
    """The user's own marker list was removed (supervisor, 2026-11-03), so the page data has no
    upload, template or "Markers from" setting."""
    page = payloads.annotation_payload(full_run)
    assert not {"marker_list", "marker_list_template"} & set(page)
    assert "source" not in page["settings"]["current"] and "source" not in page["results"]["settings"]


def test_nothing_matched_message():
    assert payloads.nothing_matched_text({"nothing_matched": False}, "brain") is None
    assert payloads.nothing_matched_text({"nothing_matched": True}, "skin of body") == (
        "No cluster matched the cell types in the Census reference for skin of body. Check that the species and "
        "tissue on Data Overview fit your sample, compare with CellTypist's answer, or review the clusters yourself.")


def test_label_review_on_the_pages(spare_run):
    """The review shows on Cell Type Annotation (per cluster, with counts), in Explore Results
    and in the run log; confirming keeps the prediction's name, renaming uses the analyst's."""
    run = spare_run
    run.annotate(replace(run.params.annotation, test=replace(run.params.annotation.test, max_adj_p=0.05)))
    clusters = payloads.annotation_payload(run)["results"]["clusters"]
    first, second = clusters[0], clusters[1]
    run.review_cluster(*payloads.review_params_from_browser({"cluster": first["cluster"], "decision": "confirm"}))
    run.review_cluster(*payloads.review_params_from_browser(
        {"cluster": second["cluster"], "decision": "rename", "label": "CD14+ monocytes", "note": "LYZ high"}))

    review = payloads.annotation_payload(run)["results"]["review"]
    json.dumps(review, allow_nan=False)
    final = {f["cluster"]: f["final_label"] for f in payloads.annotation_payload(run)["results"]["final"]}
    assert review["clusters"][first["cluster"]]["label"] == final[first["cluster"]]
    assert review["clusters"][first["cluster"]]["prediction_changed"] is False
    assert review["clusters"][second["cluster"]] | {"reviewed_at": None} == {
        "cluster": second["cluster"], "decision": "renamed", "label": "CD14+ monocytes", "note": "LYZ high",
        "reviewed_at": None, "prediction_changed": False}
    assert (review["n_confirmed"], review["n_renamed"], review["n_not_reviewed"]) == (1, 1, len(clusters) - 2)
    assert payloads.explore_payload(run)["results"]["reviewed"] == {
        first["cluster"]: final[first["cluster"]], second["cluster"]: "CD14+ monocytes"}
    log = payloads.run_log(run)[-1]
    assert log["label"] == "Label review"
    assert log["text"].startswith(f"Cluster {second['cluster']}: renamed to CD14+ monocytes · 1 confirmed, 1 renamed")
    assert "You reviewed 2 of" in step(payloads.workflow_payload(run), "cell_type_annotation")["confirmation"]

    # Annotation re-run with other settings: a confirmed label whose prediction changed is flagged.
    run.annotate(replace(run.params.annotation, test=replace(run.params.annotation.test, max_adj_p=1e-300)))
    again = payloads.annotation_payload(run)["results"]["review"]["clusters"]
    assert again[first["cluster"]]["prediction_changed"] is (final[first["cluster"]] != "Unknown")

    for bad, message in [(None, "Choose a cluster"), ({"cluster": "0", "decision": "approve"}, "confirm the prediction")]:
        with pytest.raises(payloads.InputError, match=message):
            payloads.review_params_from_browser(bad)


# --- Explore Results page -------------------------------------------------------------------------

def test_explore_needs_clusters(qc_run):
    payload = payloads.explore_payload(qc_run)
    assert payload["available"] is False and payload["results"] is None
    assert payloads.gene_payload(qc_run, "CD14") is None


def test_explore_shows_what_the_steps_computed(full_run, reference_results, reference_summary):
    payload = payloads.explore_payload(full_run)
    json.dumps(payload, allow_nan=False)
    r = payload["results"]
    clusters = payloads.cluster_payload(full_run)["results"]
    assert r["cell_cluster"] == clusters["labels"] and r["umap"] == clusters["umap"]
    annotated = payloads.annotation_payload(full_run)["results"]
    assert r["predictions"]["cell_label"] == annotated["cell_label"]
    assert r["predictions"]["label_order"] == list(dict.fromkeys(c["final_label"] for c in annotated["final"]))
    for cluster, genes in reference_summary["top25_markers"].items():
        assert r["top_markers"][cluster] == genes[:payloads.EXPLORE_TOP_MARKERS]
    assert len(r["genes"]) == reference_summary["filtered_shape"][1] and r["genes"] == sorted(r["genes"])
    assert r["unavailable"] == []


def test_gene_values_come_straight_from_the_data(full_run):
    import scanpy as sc

    from scrna_analysis.markers import expression_by_group

    adata = full_run.result("cluster")
    gene = "MS4A1"
    payload = payloads.gene_payload(full_run, gene.lower())  # case doesn't matter
    assert payload["ok"] and payload["gene"] == gene
    direct = sc.get.obs_df(adata, keys=[gene], use_raw=True)[gene].to_numpy()
    assert payload["values"] == pytest.approx(direct.tolist(), abs=1e-6)
    assert payload["n_expressing"] == int((direct > 0).sum())
    expected = expression_by_group(adata, [gene], "leiden", use_raw=True)
    assert [b["mean_expr"] for b in payload["by_cluster"]] == pytest.approx(expected["mean_expr"].tolist())
    assert [b["frac_expressing"] for b in payload["by_cluster"]] == pytest.approx(expected["frac_expressing"].tolist())


def test_unknown_gene_is_explained_with_close_matches(full_run):
    payload = payloads.gene_payload(full_run, "CD144")
    assert payload["ok"] is False and payload["message"] == "Gene 'CD144' isn't in this data."
    assert "CD14" in payload["suggestions"]
    assert payloads.gene_payload(full_run, "") is None


def test_explore_leaves_out_steps_that_need_rerunning(spare_run):
    spare_run.cluster(LeidenParams(resolution=0.5))  # markers and cell types are now out of date
    r = payloads.explore_payload(spare_run)["results"]
    assert r["top_markers"] is None and r["predictions"] is None
    assert r["unavailable"] == ["Marker Genes needs re-running", "Cell Type Annotation needs re-running"]



# --- Inspecting cells on the QC page: removal by barcode, with the filters each cell fails --------

def test_failures_agree_with_the_filter_counts(qc_run):
    preview = payloads.preview_payload(qc_run, {**REFERENCE_BOXES, "max_mt": "5"})
    counted = {}
    for names in preview["failures"].values():
        for name in names:
            counted[name] = counted.get(name, 0) + 1
    expected = qc_run.preview_filter(dataclasses.replace(REFERENCE_PRESET.filter, max_pct_counts_mt=5.0))
    assert counted == {k: v for k, v in expected.cells_failing_each.items() if v}
    assert len(preview["removed_barcodes"]) == preview["n_cells_before"] - preview["n_cells_kept"]


def test_changing_a_threshold_changes_status_not_identity(qc_run):
    """The same barcodes are reported; only which ones are removed (and why) changes."""
    cells = payloads.qc_payload(qc_run)["cells"]
    reference = payloads.preview_payload(qc_run, REFERENCE_BOXES)
    stricter = payloads.preview_payload(qc_run, {**REFERENCE_BOXES, "max_mt": "5"})
    pct = dict(zip(cells["barcode"], cells["pct_mt"]))
    assert set(reference["removed_barcodes"]) <= set(stricter["removed_barcodes"]) <= set(cells["barcode"])
    newly = set(stricter["removed_barcodes"]) - set(reference["removed_barcodes"])
    assert newly and all(5 <= pct[b] < 6 for b in newly)  # exactly the cells between the two limits
    assert all("max_pct_counts_mt" in stricter["failures"][b] for b in newly)


def test_switched_off_filter_is_never_a_reason(qc_run):
    preview = payloads.preview_payload(qc_run, {**REFERENCE_BOXES, "max_mt_off": True, "max_genes_off": True})
    reasons = {name for names in preview["failures"].values() for name in names}
    assert reasons <= {"min_genes"}


# --- Run all remaining steps ------------------------------------------------------------------

def test_remaining_stages_are_the_ones_not_up_to_date(raw_data_path):
    """Every stage, in order, until something is run; stage names map to real methods."""
    assert payloads.remaining_stages(AnalysisRun(raw_data_path, REFERENCE_PRESET)) == list(payloads.STAGES_IN_ORDER)
    assert set(payloads.STAGE_METHODS) == set(payloads.STAGES_IN_ORDER)
    assert all(callable(getattr(AnalysisRun, method)) for method in payloads.STAGE_METHODS.values())


def test_run_remaining_payload_reports_progress_by_step():
    assert payloads.run_remaining_payload(None) is None
    stages = ("qc_metrics", "filter", "preprocess", "reduce_dimensions")
    running = {"id": 3, "stages": stages, "left": ("filter", "preprocess", "reduce_dimensions"), "status": "running"}
    # QC is two stages but one step, so the total counts steps, not stages.
    assert payloads.run_remaining_payload(running) == {
        "id": 3, "status": "running", "current": "Quality Control", "index": 1, "total": 3, "failed_step": None,
        "stopped_step": None}
    done = {**running, "left": (), "status": "done"}
    assert payloads.run_remaining_payload(done)["current"] is None


def test_run_remaining_stops_where_a_choice_is_needed(raw_data_path):
    """Cell Type Annotation has no default tissue, so the job stops there (not a failure)."""
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    assert payloads.choice_needed(run, "annotate") == (
        "Confirm the species and tissue on Data Overview, then run Cell Type Annotation.")
    other = AnalysisRun(raw_data_path, DEFAULT_PARAMS, metadata={"species": "other"})
    assert payloads.choice_needed(other, "annotate").startswith("Cell types can't be predicted for this species")
    assert all(payloads.choice_needed(run, stage) is None for stage in payloads.STAGES_IN_ORDER if stage != "annotate")
    chosen = replace(DEFAULT_PARAMS, annotation=replace(
        DEFAULT_PARAMS.annotation, markers=replace(DEFAULT_PARAMS.annotation.markers, species="human"),
        census=replace(DEFAULT_PARAMS.annotation.census, tissue="blood")))
    assert payloads.choice_needed(AnalysisRun(raw_data_path, chosen), "annotate") is None
    job = {"stages": ("qc_metrics", "filter", "preprocess", "markers", "annotate")}
    assert payloads.steps_done_before(job, "annotate") == 3  # Quality Control, Preprocessing, Marker Genes


# --- Mitochondrial genes setting (any species) ------------------------------------------------

def test_qc_setting_starts_at_the_pattern_found_in_the_data(mouse_10x_path):
    run = AnalysisRun(mouse_10x_path, REFERENCE_PRESET)
    run.load()
    settings = payloads.qc_settings(run)
    assert settings["suggested"] == "mt-" and settings["mito_prefix"] == "mt-"  # not the preset's "MT-"
    assert {"prefix": "mt:", "species": "fruit fly"} in settings["patterns"]
    run.compute_qc_metrics(*payloads.qc_params_from_browser({"mito_prefix": " MT- "}, run.params))
    assert payloads.qc_settings(run)["mito_prefix"] == "MT-"  # after QC: what it ran with


def test_run_remaining_uses_the_mito_pattern_found_in_the_data(mouse_10x_path):
    """"Run the remaining steps" runs QC as the page would: with the data's own pattern ("mt-"),
    so mouse data gets a % mitochondrial filter and Cell Type Annotation suggests mouse."""
    run = AnalysisRun(mouse_10x_path, DEFAULT_PARAMS)
    run.load()
    payloads.run_stage(run, "qc_metrics")
    assert run.params.qc.mito_prefix == "mt-"
    assert run.params.filter.max_pct_counts_mt is not None
    assert payloads.workflow_payload(run)["dataset"]["metadata"]["suggested_species"] == "mouse"
    # Once QC has run, its prefix is kept (a re-run uses what the user chose).
    run.compute_qc_metrics(replace(run.params.qc, mito_prefix="MT-"))
    payloads.run_stage(run, "qc_metrics")
    assert run.params.qc.mito_prefix == "MT-"


def test_qc_settings_from_the_browser():
    assert payloads.qc_params_from_browser(None, REFERENCE_PRESET) == (REFERENCE_PRESET.qc,)
    assert payloads.qc_params_from_browser({"mito_prefix": "mt:"}, REFERENCE_PRESET)[0].mito_prefix == "mt:"
    with pytest.raises(payloads.InputError, match="Enter the start of the mitochondrial gene names"):
        payloads.qc_params_from_browser({"mito_prefix": "  "}, REFERENCE_PRESET)


def test_live_mito_count(mouse_10x_path):
    run = AnalysisRun(mouse_10x_path, REFERENCE_PRESET)
    assert payloads.mito_check_payload(run, "mt-") is None  # nothing loaded yet
    run.load()
    check = payloads.mito_check_payload(run, "mt-")
    assert check["n_genes"] > 0 and all(g.startswith("mt-") for g in check["examples"])
    assert payloads.mito_check_payload(run, "MT-")["n_genes"] == 0


# --- Suggested thresholds on the page ---------------------------------------------------------

def test_suggestion_for_other_mads_and_unusable_ones(qc_run):
    base = payloads.suggestion_payload(qc_run, None)
    looser = payloads.suggestion_payload(qc_run, {"n_mads_genes": "8", "n_mads_mito": "6"})
    assert looser["suggest_settings"]["n_mads_mito"] == 6 and looser["problem"] is None
    assert looser["suggested"]["max_pct_counts_mt"] > base["suggested"]["max_pct_counts_mt"]
    bad = payloads.suggestion_payload(qc_run, {"n_mads_genes": "0", "n_mads_mito": "3"})
    assert bad["suggested"] is None and "'MADs for genes per cell' must be greater than 0" in bad["problem"]


def test_no_suggestion_before_qc(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load()
    s = payloads.suggestion_payload(run, None)
    assert s["suggested"] is None and s["suggest_settings"]["n_mads_genes"] == 5


def test_normalize_to_the_median_of_the_cells():
    values = {"target_sum": "", "target_median": True, "min_mean": "0.0125", "max_mean": "3", "min_disp": "0.5",
              "max_value": "10", "regress": False}
    normalize, *_ = payloads.preprocess_params_from_browser(values, DEFAULT_PARAMS)
    assert normalize.target_sum is None
    assert payloads.preprocess_settings(DEFAULT_PARAMS)["target_sum"] is None
    with pytest.raises(payloads.InputError, match="Counts per cell after normalization"):
        payloads.preprocess_params_from_browser({**values, "target_median": False}, DEFAULT_PARAMS)


def test_hvg_method_from_the_browser():
    base = {"target_sum": "", "target_median": True, "min_mean": "0.0125", "max_mean": "3", "min_disp": "0.5",
            "max_value": "10", "regress": False}
    _, hvg, *_ = payloads.preprocess_params_from_browser({**base, "hvg_method": "top", "n_top_genes": "3000"},
                                                         DEFAULT_PARAMS)
    assert hvg.n_top_genes == 3000
    _, hvg, *_ = payloads.preprocess_params_from_browser({**base, "hvg_method": "cutoffs"}, replace_hvg(3000))
    assert hvg.n_top_genes is None  # back to the cut-offs
    with pytest.raises(payloads.InputError, match="'Number of genes' must be a whole number"):
        payloads.preprocess_params_from_browser({**base, "hvg_method": "top", "n_top_genes": "2.5"}, DEFAULT_PARAMS)
    settings = payloads.preprocess_settings(DEFAULT_PARAMS)
    assert settings["hvg_method"] == "cutoffs" and settings["n_top_genes"] == 2000


def replace_hvg(n):
    return dataclasses.replace(DEFAULT_PARAMS, hvg=dataclasses.replace(DEFAULT_PARAMS.hvg, n_top_genes=n))


# --- Dataset metadata: confirming species and tissue on Data Overview -------------------------

def test_later_steps_are_locked_until_species_and_tissue_are_confirmed(mouse_10x_path):
    """In the app, uploaded data has no metadata: Data Overview waits for the confirmation, every
    later step is locked with the reason, and Quality Control can't run. Confirming unlocks them."""
    run = AnalysisRun(mouse_10x_path, DEFAULT_PARAMS, require_metadata=True)
    before = payloads.workflow_payload(run)
    assert all(s["locked"] for s in before["steps"][1:]) and before["dataset"]["metadata"] is None  # nothing loaded
    run.load()
    payload = payloads.workflow_payload(run)
    overview = step(payload, "data_overview")
    assert overview["status"] == "not_run" and overview["waiting"] == "Confirm species and tissue"
    for s in payload["steps"][1:]:
        assert s["locked"] and s["prerequisite"] == "Data Overview" and s["prerequisite_text"] == payloads.CONFIRM_FIRST
    assert payloads.qc_payload(run)["can_run"] is False
    meta = payload["dataset"]["metadata"]
    assert (meta["confirmed"], meta["suggested_species"], meta["mito_prefix"]) == (False, "mouse", "mt-")
    assert meta["species_options"][-1] == {"value": "other", "label": "Other species"}
    # Census's mouse tissues, those with a marker reference marked (brain), each with its CellTypist model.
    brain = next(t for t in meta["census_tissues"]["mouse"] if t["tissue"] == "brain")
    assert brain == {"tissue": "brain", "label": "Brain", "n_cells": brain["n_cells"], "has_table": True,
                     "celltypist_model": "Mouse_Whole_Brain.pkl"}
    assert not {"tissues", "tissue_cell_types", "tissues_available", "tissue_problem"} & set(meta)  # no PanglaoDB

    run.confirm_metadata("mouse", "brain")
    payload = payloads.workflow_payload(run)
    assert step(payload, "data_overview")["status"] == "done" and step(payload, "data_overview")["waiting"] is None
    assert not any(s["locked"] for s in payload["steps"]) and payloads.qc_payload(run)["can_run"] is True
    assert "(Mouse · Brain)" in step(payload, "data_overview")["confirmation"]
    assert payload["run_log"][-1] == {"label": "Data Overview", "text": "Species and tissue confirmed: Mouse · Brain",
                                      "ran_at": run.history[-1]["ran_at"]}


def test_species_that_disagrees_with_the_gene_names_is_flagged(mouse_10x_path):
    run = AnalysisRun(mouse_10x_path, DEFAULT_PARAMS, require_metadata=True)
    run.load()
    meta = payloads.metadata_payload(run)
    assert meta["suggested_species"] == "mouse" and meta["suggestion_reason"] == (
        "Mouse: Ensembl IDs start with ENSMUSG, gene names are written like Cd3d and mitochondrial genes are named mt-…")
    run.confirm_metadata("human", "brain")
    assert payloads.metadata_payload(run)["mismatch"] == (
        "Check the species: you chose Human, but Ensembl IDs start with ENSMUSG (as in mouse data), gene names are "
        "written like Cd3d (as in mouse data) and mitochondrial genes are named mt-… (as in mouse data).")
    run.confirm_metadata("mouse", "brain")
    assert payloads.metadata_payload(run)["mismatch"] is None
    run.confirm_metadata("other")
    assert payloads.metadata_payload(run)["mismatch"].startswith("The gene names look like mouse data")


def test_species_suggestion_for_species_the_app_does_not_list():
    clues = [{"species": "zebrafish", "clue": "Ensembl IDs start with ENSDARG"}]
    assert payloads._species_suggestion(clues) == ("other", "Zebrafish: Ensembl IDs start with ENSDARG")
    assert payloads._species_suggestion([]) == (None, None)
    assert payloads._species_mismatch("other", clues) is None  # nothing the app could predict instead


def test_metadata_from_the_browser():
    assert payloads.metadata_from_browser({"species": "human", "tissue": "blood", "about": "x"}) == (
        "human", "blood", "x", None)
    assert payloads.metadata_from_browser({"species": "human", "tissue": "blood", "tissues": ["Blood"]}) == (
        "human", "blood", None, None)  # PanglaoDB tissues sent by an old page are ignored
    assert payloads.metadata_from_browser({"species": "other"}) == ("other", None, None, None)
    for values, message in [({}, "Choose the species"), ({"species": "zebrafish"}, "Choose the species"),
                            ({"species": "mouse"}, "Choose the tissue"),
                            ({"species": "mouse", "tissues": ["Brain"]}, "Choose the tissue"),
                            ({"species": "mouse", "tissue": 3}, "Choose the tissue"),
                            ({"species": "other", "about": 3}, "must be text")]:
        with pytest.raises(payloads.InputError, match=message):
            payloads.metadata_from_browser(values)


def test_changed_tissue_explains_why_annotation_needs_rerunning(raw_data_path):
    run = run_analysis(raw_data_path, REFERENCE_PRESET)
    run.confirm_metadata("human", "spleen")
    annotation = step(payloads.workflow_payload(run), "cell_type_annotation")
    assert annotation["status"] == "needs_rerun"
    assert annotation["reason"] == "The species or tissue changed on Data Overview, so run this step again."
    assert payloads.annotation_payload(run)["metadata"]["tissue"] == "Spleen"



# --- Unresolved clusters: why, and what to check ---------------------------------------------

def test_unresolved_note_mentions_species_and_tissue_only_when_most_clusters_fail():
    def summary(labels, nothing=False):
        return {"nothing_matched": nothing, "final": [{"final_label": label} for label in labels]}
    assert payloads.unresolved_note(summary(["B cells", "T cells"])) is None
    assert payloads.unresolved_note(summary(["Unknown"] * 2, nothing=True)) is None  # the other box says it
    some = payloads.unresolved_note(summary(["Unknown", "B cells", "T cells"]))
    assert some.startswith("1 cluster is Unknown: the checks disagree or neither is confident") and "species" not in some
    assert "Open the To review tab and click a cluster to see why" in some
    most = payloads.unresolved_note(summary(["Unknown", "Unknown", "B cells"]))
    assert most.startswith("2 of 3 clusters are Unknown. When most clusters fail, check that the species")


def test_unresolved_guidance_follows_the_reason():
    from types import SimpleNamespace
    row = lambda status: SimpleNamespace(cluster=3, status=status, best_candidate="astrocyte",  # noqa: E731
                                         second_candidate="glioblast")
    assert payloads.unresolved_guidance(row("ambiguous"), "brain").startswith(
        "Two unrelated cell types fit almost equally: astrocyte and glioblast.")
    assert "rule out chance" in payloads.unresolved_guidance(row("not_significant"), "brain")
    assert "Census reference for brain" in payloads.unresolved_guidance(row("too_few_markers"), "brain")
    assert payloads.unresolved_guidance(row("no_reference"), "lung") == (
        "There's no CELLxGENE Census reference for lung yet, so the marker check didn't run. Use CellTypist's "
        "answer and the cluster's top marker genes in your review.")
    # Without CellTypist's answer the sentences don't point to it.
    assert all("CellTypist" not in payloads.unresolved_guidance(row(status), "lung", celltypist=False)
               for status in ("no_reference", "ambiguous", "not_significant", "too_few_markers"))
    assert "Closest match in another tissue" not in payloads.unresolved_guidance(row("ambiguous"), "brain")  # hints gone


def test_every_unresolved_cluster_gets_guidance(mouse_10x_path):
    """Mouse brain with the Census brain reference: each Unresolved cluster gets a reason sentence."""
    run = AnalysisRun(mouse_10x_path, DEFAULT_PARAMS, metadata={"species": "mouse", "tissue": "brain"})
    run.load()
    for stage in payloads.remaining_stages(run):
        payloads.run_stage(run, stage)
    s = run.stage_result("annotate").summary
    unresolved = {str(c["cluster"]) for c in s["clusters"] if c["predicted_label"] == "Unresolved"}
    results = payloads.annotation_payload(run)["results"]
    assert unresolved and set(results["guidance"]) == unresolved
    assert not {"best_anywhere", "tissue_hint"} & set(s)
    unknown = sum(f["final_label"] == "Unknown" for f in s["final"])
    if unknown:  # the note is about the final labels (CellTypist may name clusters the marker check couldn't)
        assert results["unresolved_note"].startswith(
            f"{unknown} of {len(s['final'])} clusters are Unknown" if unknown > len(s["final"]) / 2
            else f"{unknown} cluster" + ("s are" if unknown > 1 else " is") + " Unknown")


def test_batch_correction_from_the_page():
    P = REFERENCE_PRESET
    boxes = {**as_boxes(payloads.reduction_settings(P)), "batch_correction": True, "batch_key": "condition"}
    *_, integ = payloads.reduction_params_from_browser(boxes, P)
    assert (integ.method, integ.batch_key) == ("harmony", "condition")
    for bad, message in [({"batch_correction": "yes"}, "whether to correct"), ({"batch_correction": True, "batch_key": ""}, "column")]:
        with pytest.raises(payloads.InputError, match=message):
            payloads.reduction_params_from_browser({**as_boxes(payloads.reduction_settings(P)), **bad}, P)
    assert payloads.reduction_settings(P)["batch_correction"] is False  # off by default
