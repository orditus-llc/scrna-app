"""Server tests for the React app (app/app.py), with no browser.

Uses Shiny's built-in test server (`shiny.testserver.test_server`): it runs the real server
function, lets a test set inputs (e.g. a button click) and returns exactly the JSON the
page would receive.
"""

from __future__ import annotations

import asyncio
import sys
import time

import pytest
from shiny.testserver import TestServerSession, test_server, test_server_async

from tests.reference_preset import REFERENCE_PRESET
from scrna_analysis.pipeline import AnalysisRun, PipelineError
from tests.conftest import PROJECT_ROOT

sys.path.insert(0, str(PROJECT_ROOT / "app"))
import payloads  # noqa: E402

APP = PROJECT_ROOT / "app" / "app.py"


def refresh(ts, max_rounds: int = 50) -> None:
    """Let the server finish before a value is read.

    set_inputs can return after a flush in which the outputs were refreshed before the step
    ran, and the step's own refresh follows in a later pass. Each flush() gives the server a
    turn and re-reads the values; once two reads in a row are equal, nothing is pending (a step
    runs to the end within one turn). Async sessions (the "Run the remaining steps" tests)
    await their own flush, so theirs isn't started here."""
    if not isinstance(ts, TestServerSession):
        return
    previous = ts.flush().to_values()
    for _ in range(max_rounds):
        current = ts.flush().to_values()
        if current == previous:
            return
        previous = current
    raise AssertionError("The server's outputs kept changing; nothing settled to read.")


def workflow(ts) -> dict:
    refresh(ts)
    output = ts.get_output("workflow")
    assert output.status == "ok", output.error
    return output.value


def test_app_starts_with_nothing_run():
    with test_server(APP, timeout_secs=30) as ts:
        data = workflow(ts)
        assert [s["status"] for s in data["steps"]] == ["not_run"] * 7
        assert data["notice"] is None


def test_initial_zero_is_not_a_click():
    """The page sends load_demo=0 when it opens; that must not load anything."""
    with test_server(APP, timeout_secs=30) as ts:
        ts.set_inputs(load_demo=0)
        assert workflow(ts)["steps"][0]["status"] == "not_run"


def test_clicking_load_runs_the_load_stage(reference_summary):
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(load_demo=1)
        data = workflow(ts)
        n_cells, n_genes = reference_summary["raw_shape"]
        assert data["steps"][0]["status"] == "done"
        assert data["steps"][0]["summary"].startswith(f"{n_cells:,} cells × {n_genes:,} genes")
        assert data["dataset"]["n_cells"] == n_cells
        assert [entry["label"] for entry in data["run_log"]] == ["Data Overview"]
        assert data["notice"] is None


def test_a_failing_stage_becomes_a_message(monkeypatch):
    def broken_load(self):
        raise PipelineError("The data file could not be read.")

    monkeypatch.setattr(AnalysisRun, "load", broken_load)
    with test_server(APP, timeout_secs=30) as ts:
        ts.set_inputs(load_demo=1)
        data = workflow(ts)
        assert data["notice"] == {"id": 1, "kind": "error", "text": "The data file could not be read.",
                                  "step": "data_overview"}  # shown in that step's card too
        assert data["steps"][0]["status"] == "not_run"  # nothing changed

        ts.set_inputs(load_demo=2)  # a second failure gets a new id, so the page shows it again
        assert workflow(ts)["notice"]["id"] == 2


# --- Loading data: upload or demo -------------------------------------------------------------

def uploaded(path, name: str) -> list[dict]:
    """What Shiny's upload gives the server: the original name and where the file was stored."""
    return [{"name": name, "size": path.stat().st_size, "type": "", "datapath": str(path)}]


def test_uploading_a_10x_file_starts_its_analysis(mouse_10x_path):
    """An upload is loaded and checked, then waits for its species and tissue: Quality Control
    and "Run the remaining steps" refuse until they are confirmed."""
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(upload_data=uploaded(mouse_10x_path, "my_brain_cells.h5"))
        data = workflow(ts)
        assert data["notice"] is None
        assert data["steps"][0]["status"] == "not_run" and data["steps"][0]["waiting"] == "Confirm species and tissue"
        assert all(s["locked"] for s in data["steps"][1:])
        ts.set_inputs(run_qc={"n": 1, "settings": {"mito_prefix": "mt-"}})
        assert workflow(ts)["notice"]["text"] == payloads.CONFIRM_FIRST
        ts.set_inputs(run_remaining=1)
        assert workflow(ts)["notice"]["text"] == payloads.CONFIRM_FIRST and workflow(ts)["run_remaining"] is None
        ts.set_inputs(confirm_metadata={"n": 1, "settings": {"species": "mouse", "tissues": []}})
        assert workflow(ts)["notice"]["text"] == "Choose the tissue."  # nothing changed
        ts.set_inputs(confirm_metadata={"n": 2, "settings": {"species": "mouse", "tissues": ["Brain"], "about": "E18"}})
        data = workflow(ts)
        assert data["steps"][0]["status"] == "done" and not any(s["locked"] for s in data["steps"])
        assert data["dataset"]["metadata"]["text"] == "Mouse · Brain" and data["dataset"]["metadata"]["about"] == "E18"
        assert data["dataset"]["label"] == data["dataset"]["file"] == "my_brain_cells.h5"
        assert data["dataset"]["file_format"] == "10x Genomics (.h5)"
        import scanpy as sc
        assert data["dataset"]["n_cells"] == sc.read_10x_h5(mouse_10x_path).n_obs


def test_demo_data_is_named_as_demo(reference_summary):
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(load_demo=1)
        data = workflow(ts)
        dataset = data["dataset"]
        assert dataset["label"].startswith("Demo data:")
        # The demo's metadata arrives confirmed, so nothing is locked.
        assert dataset["metadata"]["confirmed_by"] == "dataset" and dataset["metadata"]["text"] == "Human · Blood"
        assert data["steps"][0]["status"] == "done" and not any(s["locked"] for s in data["steps"])
        assert dataset["n_cells"] == reference_summary["raw_shape"][0]


def test_a_bad_upload_leaves_the_current_analysis_alone(tmp_path):
    bad = tmp_path / "0.csv"
    bad.write_text("not single-cell data")
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(upload_data=uploaded(bad, "counts.csv"))
        data = workflow(ts)
        assert data["notice"]["step"] == "data_overview"
        assert "'counts.csv' isn't a supported format" in data["notice"]["text"]
        assert data["dataset"]["label"].startswith("Demo data:")  # still the demo
        assert [e["label"] for e in data["run_log"]] == ["Data Overview", "QC metrics"]  # nothing lost


def test_new_data_replaces_the_previous_analysis(mouse_10x_path):
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(upload_data=uploaded(mouse_10x_path, "brain.h5"))
        data = workflow(ts)
        assert data["dataset"]["label"] == "brain.h5"
        assert [s["status"] for s in data["steps"]] == ["not_run"] * 7  # a fresh start, waiting for its species
        assert [e["label"] for e in data["run_log"]] == ["Data Overview"]


def test_mouse_upload_uses_the_mouse_mito_pattern(mouse_10x_path):
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(upload_data=uploaded(mouse_10x_path, "brain.h5"))
        ts.set_inputs(confirm_metadata=CONFIRM_MOUSE)
        assert output(ts, "qc")["settings"]["mito_prefix"] == "mt-"
        assert workflow(ts)["dataset"]["mito_matches"][0]["prefix"] == "mt-"
        ts.set_inputs(mito_prefix="mt-")
        assert output(ts, "mito_check")["n_genes"] > 0
        ts.set_inputs(run_qc={"n": 1, "settings": {"mito_prefix": "mt-"}})
        qc = output(ts, "qc")
        assert qc["mito"]["prefix"] == "mt-" and qc["mito"]["available"] and not qc["mito"]["filter_switched_off"]
        assert qc["filter"]["current"]["max_pct_counts_mt"] is not None


def test_no_matching_mito_genes_switches_the_filter_off(mouse_10x_path):
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(upload_data=uploaded(mouse_10x_path, "brain.h5"))
        ts.set_inputs(confirm_metadata=CONFIRM_MOUSE)
        ts.set_inputs(run_qc={"n": 1, "settings": {"mito_prefix": "MT-"}})  # the human pattern
        qc = output(ts, "qc")
        assert qc["mito"]["available"] is False and qc["mito"]["filter_switched_off"] is True
        assert qc["filter"]["current"]["max_pct_counts_mt"] is None
        assert "No mitochondrial genes found" in qc["filter"]["details"]["max_pct_counts_mt"]["note"]
        log = workflow(ts)["run_log"][-1]
        assert log["label"] == "QC metrics" and "% mitochondrial filter switched off" in log["text"]


REFERENCE_BOXES = {"min_genes": "200", "min_cells": "3", "max_genes": "2500", "max_genes_off": False,
                   "max_mt": "6", "max_mt_off": False}


def output(ts, name):
    refresh(ts)
    value = ts.get_output(name)
    assert value.status == "ok", value.error
    return value.value


def test_quality_control_flow(reference_summary):
    """Run QC -> change a threshold (live preview) -> Apply Filtering, as the page does it."""
    with test_server(APP, timeout_secs=90) as ts:
        ts.set_inputs(load_demo=1)
        assert output(ts, "qc_preview") is None  # nothing to preview before QC

        ts.set_inputs(run_qc=1)
        qc = output(ts, "qc")
        assert qc["metrics_done"] is True
        assert qc["filter"]["current"] == qc["filter"]["suggested"]  # the boxes start at the suggestion
        preview = output(ts, "qc_preview")  # before the page sends its boxes: the suggestion
        assert preview["n_cells_kept"] == preview["suggested"]["n_cells_kept"]
        ts.set_inputs(filter_settings=REFERENCE_BOXES)  # the notebook's thresholds
        assert output(ts, "qc_preview")["n_cells_kept"] == reference_summary["filtered_shape"][0]

        ts.set_inputs(filter_settings={**REFERENCE_BOXES, "max_mt": "5"})
        stricter = output(ts, "qc_preview")
        assert stricter["n_cells_kept"] < reference_summary["filtered_shape"][0]
        assert output(ts, "workflow")["steps"][1]["status"] == "not_run"  # preview changes nothing

        ts.set_inputs(filter_settings=REFERENCE_BOXES, apply_filter=1)
        workflow = output(ts, "workflow")
        assert workflow["steps"][1]["status"] == "done"
        assert workflow["next_step"]["key"] == "preprocessing"
        assert output(ts, "qc")["applied"]["n_cells_kept"] == reference_summary["filtered_shape"][0]


def test_unusable_filter_values_become_a_message():
    with test_server(APP, timeout_secs=90) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(filter_settings={**REFERENCE_BOXES, "min_genes": ""}, apply_filter=1)
        workflow = output(ts, "workflow")
        assert workflow["notice"]["text"] == "Enter a value for 'Min genes per cell'."
        assert workflow["notice"]["step"] == "quality_control"
        assert workflow["steps"][1]["status"] == "not_run"  # nothing was applied
        assert output(ts, "qc_preview")["ok"] is False


def test_apply_uses_the_boxes_sent_with_the_click(reference_summary):
    """A value typed just before clicking can't be lost: Apply sends the boxes with the click."""
    with test_server(APP, timeout_secs=90) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(filter_settings={**REFERENCE_BOXES, "max_mt": "5"})  # preview still on 5%
        ts.set_inputs(apply_filter={"n": 1, "settings": REFERENCE_BOXES})
        assert output(ts, "qc")["applied"]["n_cells_kept"] == reference_summary["filtered_shape"][0]


def boxes(settings: dict) -> dict:
    return {k: v if isinstance(v, bool) else str(v) for k, v in settings.items()}


PREPROCESS_BOXES = boxes(payloads.preprocess_settings(REFERENCE_PRESET))
REDUCTION_BOXES = boxes(payloads.reduction_settings(REFERENCE_PRESET))
# Cell Type Annotation as the page sends it: the test's checks (defaults). Species and tissues
# come from the dataset's metadata (the demo's: human; Immune system and Blood).
ANNOTATION_BOXES = {"top_n_genes": "100", "max_adj_p": "0.05", "min_overlap": "3", "ambiguity_fold": "10"}
DEMO_TISSUES = ["Blood", "Immune system"]
# Confirming an uploaded mouse brain file's species and tissue on Data Overview.
CONFIRM_MOUSE = {"n": 1, "settings": {"species": "mouse", "tissues": ["Brain"], "about": ""}}


# With regression off, scaling a sparse matrix makes it dense; Scanpy says so (expected).
@pytest.mark.filterwarnings("ignore:zero-centering a sparse array:UserWarning")
def test_preprocessing_and_pca_umap_from_their_cards(reference_summary):
    """Run both cards as the page does, then redo Preprocessing with regression off."""
    with test_server(APP, timeout_secs=180) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(apply_filter={"n": 1, "settings": REFERENCE_BOXES})
        assert output(ts, "preprocessing")["status"] == "not_run"

        ts.set_inputs(run_preprocess={"n": 1, "settings": PREPROCESS_BOXES})
        assert output(ts, "preprocessing")["results"]["n_hvg"] == reference_summary["n_hvg"]

        ts.set_inputs(run_reduction={"n": 1, "settings": REDUCTION_BOXES})
        reduction = output(ts, "reduction")["results"]
        assert reduction["variance_ratio"] == pytest.approx(reference_summary["pca_variance_ratio"], abs=1e-6)
        data = workflow(ts)
        assert [s["status"] for s in data["steps"][:4]] == ["done"] * 4
        assert data["next_step"]["key"] == "clustering"

        ts.set_inputs(run_preprocess={"n": 2, "settings": {**PREPROCESS_BOXES, "regress": False}})
        data = workflow(ts)
        assert data["steps"][2]["status"] == "done" and "regression off" in data["steps"][2]["summary"]
        assert data["steps"][3]["status"] == "needs_rerun"  # Dimensionality Reduction was cleared
        assert output(ts, "preprocessing")["settings"]["current"]["regress"] is False


def test_unusable_card_values_change_nothing():
    with test_server(APP, timeout_secs=120) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(apply_filter={"n": 1, "settings": REFERENCE_BOXES})
        ts.set_inputs(run_preprocess={"n": 1, "settings": PREPROCESS_BOXES})
        ts.set_inputs(run_reduction={"n": 1, "settings": {**REDUCTION_BOXES, "n_comps": "10", "n_pcs": "20"}})
        data = workflow(ts)
        assert data["notice"]["text"] == "'PCs used for the neighbor graph' can't be more than 'Principal components'."
        assert data["notice"]["step"] == "dimensionality_reduction"
        assert data["steps"][3]["status"] == "not_run"


def test_clustering_from_its_page(reference_summary):
    with test_server(APP, timeout_secs=180) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(apply_filter={"n": 1, "settings": REFERENCE_BOXES})
        ts.set_inputs(run_preprocess={"n": 1, "settings": PREPROCESS_BOXES})
        ts.set_inputs(run_reduction={"n": 1, "settings": REDUCTION_BOXES})
        assert output(ts, "clustering")["status"] == "not_run"
        ts.set_inputs(run_cluster={"n": 1, "settings": boxes(payloads.cluster_settings(REFERENCE_PRESET))})
        clusters = output(ts, "clustering")["results"]["clusters"]
        assert sorted(c["n_cells"] for c in clusters) == sorted(reference_summary["cluster_sizes"].values())
        assert workflow(ts)["next_step"]["key"] == "marker_genes"

        ts.set_inputs(run_cluster={"n": 2, "settings": {"resolution": "-1", "seed": "0"}})
        assert workflow(ts)["notice"]["text"] == "'Leiden resolution' must be greater than 0."
        assert workflow(ts)["steps"][4]["status"] == "done"  # the earlier result is untouched


def test_marker_genes_from_their_page(reference_summary):
    with test_server(APP, timeout_secs=180) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(apply_filter={"n": 1, "settings": REFERENCE_BOXES})
        ts.set_inputs(run_preprocess={"n": 1, "settings": PREPROCESS_BOXES})
        ts.set_inputs(run_reduction={"n": 1, "settings": REDUCTION_BOXES})
        ts.set_inputs(run_cluster={"n": 1, "settings": boxes(payloads.cluster_settings(REFERENCE_PRESET))})
        ts.set_inputs(run_markers={"n": 1, "settings": {"method": "wilcoxon"}})
        top = output(ts, "markers")["results"]["top"]
        for cluster, genes in reference_summary["top25_markers"].items():
            assert [row["gene"] for row in top[cluster][:25]] == genes

        ts.set_inputs(run_markers={"n": 2, "settings": {"method": "t-test"}})
        data = workflow(ts)
        assert "t-test" in data["steps"][5]["summary"]


def test_whole_analysis_from_the_pages():
    """All seven steps the way the pages run them, with cell types from PanglaoDB."""
    with test_server(APP, timeout_secs=240) as ts:
        ts.set_inputs(load_demo=1)
        ts.set_inputs(run_qc=1)
        ts.set_inputs(apply_filter={"n": 1, "settings": REFERENCE_BOXES})
        ts.set_inputs(run_preprocess={"n": 1, "settings": PREPROCESS_BOXES})
        ts.set_inputs(run_reduction={"n": 1, "settings": REDUCTION_BOXES})
        ts.set_inputs(run_cluster={"n": 1, "settings": boxes(payloads.cluster_settings(REFERENCE_PRESET))})
        ts.set_inputs(run_markers={"n": 1, "settings": {"method": "wilcoxon"}})
        ts.set_inputs(run_annotation={"n": 1, "settings": ANNOTATION_BOXES})
        result = output(ts, "annotation")["results"]
        # The demo's tissue is Census's "blood": Census first, PanglaoDB second (its tissues mapped from blood).
        assert result["database"]["name"] == "CELLxGENE Census" and result["settings"]["census_tissue"] == "blood"
        assert [r["key"] for r in result["references"]] == ["Census", "PanglaoDB"]
        assert result["settings"]["tissues"] == DEMO_TISSUES
        for c in result["clusters"]:
            assert c["predicted_label"] in (c["best_candidate"], "Unresolved")
        data = workflow(ts)
        assert [s["status"] for s in data["steps"]] == ["done"] * 7
        assert data["next_step"] is None

        # Label review, as the page sends it: one cluster renamed, then a missing label explained.
        cluster = result["clusters"][0]["cluster"]
        ts.set_inputs(review_cluster={"n": 1, "settings": {"cluster": cluster, "decision": "rename", "label": "My cells"}})
        review = output(ts, "annotation")["results"]["review"]
        assert review["clusters"][cluster]["label"] == "My cells" and review["n_renamed"] == 1
        assert output(ts, "explore")["results"]["reviewed"] == {cluster: "My cells"}
        ts.set_inputs(review_cluster={"n": 2, "settings": {"cluster": cluster, "decision": "rename", "label": " "}})
        assert workflow(ts)["notice"]["text"] == f"Type a label for cluster {cluster}, or confirm the prediction."


# --- Run all remaining steps ------------------------------------------------------------------
# The job runs one stage per timer tick, so these tests use the async test server: it keeps
# one event loop running, which the timer needs (the synchronous one only runs its loop during
# each call).

@pytest.fixture
def fresh_reactive_lock():
    """Shiny keeps one lock for its reactive graph per Python process, tied to the first event
    loop that waits on it. Each async test below runs its own loop (asyncio.run), so the lock is
    reset first. Test-only: the real app runs a single loop. (Private Shiny attribute; check it
    still exists if Shiny is updated.)"""
    from shiny.reactive import _core
    _core._reactive_environment._lock = None
    yield
    _core._reactive_environment._lock = None


async def wait_for_job(ts, timeout: float = 240) -> dict:
    """Let the server work through a "Run all remaining steps" job; return how it ended,
    with the progress the page would have seen on the way."""
    deadline = time.monotonic() + timeout
    seen_running = []
    while time.monotonic() < deadline:
        await asyncio.sleep(0.2)
        await ts.flush()
        job = workflow(ts)["run_remaining"]
        if job and job["status"] == "running":
            seen_running.append((job["index"], job["current"]))
        elif job:
            return {**job, "seen_running": seen_running}
    raise AssertionError("the job didn't finish in time")


def test_run_remaining_from_the_start_uses_the_suggestions(fresh_reactive_lock):
    """One click after loading the demo runs every later step, filtering with the thresholds
    suggested from the data; progress is reported step by step. The demo comes with a description
    of the sample (human blood), so Cell Type Annotation runs too, with those tissues."""
    async def check():
        async with test_server_async(APP, timeout_secs=300) as ts:
            await ts.set_inputs(load_demo=1)
            assert workflow(ts)["remaining_steps"] == 6
            await ts.set_inputs(run_remaining=1)
            job = await wait_for_job(ts)
            assert job["status"] == "done" and job["total"] == 6 and job["failed_step"] is None
            steps_seen = [i for i, _ in job["seen_running"]]
            assert steps_seen == sorted(steps_seen) and len(set(steps_seen)) > 1  # in order, moving on
            data = workflow(ts)
            assert [s["status"] for s in data["steps"]] == ["done"] * 7 and data["remaining_steps"] == 0
            assert data["notice"]["kind"] == "info" and data["notice"]["text"] == "Ran 6 steps with the current settings."
            assert [e["label"] for e in data["run_log"]] == [payloads.STAGE_LABELS[s] for s in payloads.STAGES_IN_ORDER]
            qc = output(ts, "qc")
            assert qc["filter"]["current"] == qc["filter"]["suggested"]  # filtered at the suggestion
            annotation = output(ts, "annotation")
            assert annotation["settings"]["current"]["species"] == "human"
            assert annotation["settings"]["current"]["tissues"] == DEMO_TISSUES
            assert annotation["metadata"]["confirmed_by"] == "dataset"
    asyncio.run(check())


def test_run_remaining_stops_where_a_choice_is_needed(mouse_10x_path, fresh_reactive_lock):
    """Data from a species without cell types in the marker database ("Other species" on Data
    Overview): the job stops at Cell Type Annotation (a planned stop, not a failure) and says why."""
    async def check():
        async with test_server_async(APP, timeout_secs=300) as ts:
            await ts.set_inputs(upload_data=uploaded(mouse_10x_path, "brain.h5"))
            await ts.set_inputs(confirm_metadata={"n": 1, "settings": {"species": "other", "tissues": []}})
            await ts.set_inputs(run_remaining=1)
            job = await wait_for_job(ts)
            assert job["status"] == "stopped" and job["total"] == 6 and job["stopped_step"] == "cell_type_annotation"
            assert job["failed_step"] is None
            data = workflow(ts)
            assert [s["status"] for s in data["steps"]] == ["done"] * 6 + ["not_run"]
            assert data["notice"]["kind"] == "info"
            assert data["notice"]["text"].startswith("Ran 5 steps. Cell types can't be predicted for this species")
            annotation = output(ts, "annotation")
            assert annotation["metadata"]["has_markers"] is False and annotation["metadata"]["species"] == "Other species"
    asyncio.run(check())


def test_run_remaining_after_the_notebook_thresholds_matches_the_reference(reference_summary, fresh_reactive_lock):
    """With the notebook's settings for filtering, preprocessing and PCA/neighbors/UMAP, one
    click runs the rest (whose defaults equal the notebook's): the marker genes equal the
    reference's. (It stops at Cell Type Annotation: no tissue chosen yet.)"""
    async def check():
        async with test_server_async(APP, timeout_secs=300) as ts:
            await ts.set_inputs(load_demo=1)
            await ts.set_inputs(run_qc={"n": 1, "settings": {"mito_prefix": "MT-"}})
            await ts.set_inputs(apply_filter={"n": 1, "settings": REFERENCE_BOXES})
            await ts.set_inputs(run_preprocess={"n": 1, "settings": PREPROCESS_BOXES})
            await ts.set_inputs(run_reduction={"n": 1, "settings": REDUCTION_BOXES})
            await ts.set_inputs(run_remaining=1)
            job = await wait_for_job(ts)
            assert job["status"] == "done" and job["total"] == 3  # clustering, markers, annotation (demo tissues)
            top = output(ts, "markers")["results"]["top"]
            for cluster, genes in reference_summary["top25_markers"].items():
                assert [row["gene"] for row in top[cluster][:25]] == genes
    asyncio.run(check())


def test_run_remaining_keeps_finished_steps_and_uses_changed_settings(fresh_reactive_lock):
    """Steps already up to date aren't run again, and a changed setting is used by the job."""
    async def check():
        async with test_server_async(APP, timeout_secs=300) as ts:
            await ts.set_inputs(load_demo=1)
            await ts.set_inputs(run_qc=1)
            await ts.set_inputs(apply_filter={"n": 1, "settings": REFERENCE_BOXES})
            await ts.set_inputs(run_preprocess={"n": 1, "settings": PREPROCESS_BOXES})
            await ts.set_inputs(run_reduction={"n": 1, "settings": REDUCTION_BOXES})
            cluster_boxes = boxes(payloads.cluster_settings(REFERENCE_PRESET))
            await ts.set_inputs(run_cluster={"n": 1, "settings": cluster_boxes})
            await ts.set_inputs(run_markers={"n": 1, "settings": {"method": "wilcoxon"}})
            await ts.set_inputs(run_annotation={"n": 1, "settings": ANNOTATION_BOXES})
            await ts.set_inputs(run_cluster={"n": 2, "settings": {**cluster_boxes, "resolution": "0.5"}})
            runs_before = len(workflow(ts)["run_log"])
            await ts.set_inputs(run_remaining=1)
            job = await wait_for_job(ts)
            assert job["status"] == "done" and job["total"] == 2  # markers + annotation only
            data = workflow(ts)
            assert [e["label"] for e in data["run_log"][runs_before:]] == ["Marker Genes", "Cell Type Annotation"]
            record = {item["label"]: item["value"] for item in data["settings_record"]}
            assert record["Leiden resolution"] == "0.5"
    asyncio.run(check())


def test_run_remaining_stops_at_a_failing_step(monkeypatch, fresh_reactive_lock):
    def broken_clustering(self, leiden=None):
        raise PipelineError("Clustering failed: example problem.")

    monkeypatch.setattr(AnalysisRun, "cluster", broken_clustering)

    async def check():
        async with test_server_async(APP, timeout_secs=300) as ts:
            await ts.set_inputs(load_demo=1)
            await ts.set_inputs(run_remaining=1)
            job = await wait_for_job(ts)
            assert job["status"] == "failed" and job["failed_step"] == "clustering"
            data = workflow(ts)
            assert [s["status"] for s in data["steps"]] == ["done"] * 4 + ["not_run"] * 3  # earlier steps kept
            assert data["notice"] == {"id": 1, "kind": "error", "text": "Clustering failed: example problem.",
                                      "step": "clustering"}
    asyncio.run(check())


def test_run_remaining_initial_zero_is_not_a_click():
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(run_remaining=0)  # the page's initial value is not a click
        assert workflow(ts)["run_remaining"] is None


@pytest.mark.parametrize("asset", ["ui.js", "ui.css"])
def test_built_interface_exists(asset):
    """The page needs the built bundle (run `npm run build` in app/ if this fails)."""
    path = PROJECT_ROOT / "app" / "www" / asset
    assert path.is_file() and path.stat().st_size > 0, f"Missing {path}. Build it: cd app && npm run build"


# --- Several samples ---------------------------------------------------------------------------

def test_several_sample_files_load_together(tmp_path):
    """Two files uploaded together, named on the page, become one dataset with two samples; the
    names travel with the click (in the files' order)."""
    import anndata as ad
    import numpy as np
    from scipy import sparse
    paths = []
    for i in range(2):
        rng = np.random.default_rng(i)
        a = ad.AnnData(sparse.csr_matrix(rng.poisson(1.0, (30, 50)).astype(np.float32)))
        a.var_names = [f"G{g}" for g in range(50)]
        a.obs_names = [f"CELL{c}" for c in range(30)]  # the same barcodes in both files
        paths.append(tmp_path / f"{i}.h5ad")
        a.write_h5ad(paths[-1])
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(upload_samples=[uploaded(paths[0], "donor_a.h5ad")[0], uploaded(paths[1], "donor_b.h5ad")[0]])
        ts.set_inputs(load_samples={"n": 1, "settings": {"samples": [{"name": "Donor A", "condition": "ctrl"},
                                                                     {"name": "Donor B", "condition": "stim"}]}})
        data = workflow(ts)
        assert data["notice"] is None, data["notice"]
        assert data["dataset"]["label"] == "2 samples" and data["dataset"]["n_cells"] == 60
        samples = data["dataset"]["samples"]
        assert samples["from_files"] is True and samples["choices"] == []
        assert [(s["name"], s["condition"], s["n_cells"], s["file"]) for s in samples["list"]] == [
            ("Donor A", "ctrl", 30, "donor_a.h5ad"), ("Donor B", "stim", 30, "donor_b.h5ad")]
        ts.set_inputs(load_samples={"n": 2, "settings": {"samples": [{"name": "same"}, {"name": "same"}]}})
        assert "repeated: same" in workflow(ts)["notice"]["text"]
        assert workflow(ts)["dataset"]["samples"]["list"][0]["name"] == "Donor A"  # a refused load changes nothing


def test_a_file_with_a_sample_column(kang_path):
    """One uploaded file whose columns name samples and conditions: the page offers the columns,
    the user picks them, and can go back to one dataset."""
    with test_server(APP, timeout_secs=120) as ts:
        ts.set_inputs(upload_data=uploaded(kang_path, "kang.h5ad"))
        samples = workflow(ts)["dataset"]["samples"]
        assert samples["list"] == [] and samples["choices"][0] == {
            "column": "replicate", "n_values": 8, "values": samples["choices"][0]["values"], "hint": "sample"}
        ts.set_inputs(define_samples={"n": 1, "settings": {"sample_column": "replicate", "condition_column": "label"}})
        samples = workflow(ts)["dataset"]["samples"]
        assert len(samples["list"]) == 16 and samples["sample_column"] == "replicate"
        assert "sample" not in {c["column"] for c in samples["choices"]}  # the app's own columns aren't offered
        ts.set_inputs(define_samples={"n": 2, "settings": {"sample_column": "", "condition_column": ""}})
        assert workflow(ts)["dataset"]["samples"]["list"] == []


def test_a_single_file_offers_no_sample_choice():
    """The demo (one sample, no text columns) shows nothing about samples."""
    with test_server(APP, timeout_secs=60) as ts:
        ts.set_inputs(load_demo=1)
        samples = workflow(ts)["dataset"]["samples"]
        assert samples["list"] == [] and samples["choices"] == [] and samples["sample_column"] is None


def test_quality_control_per_sample_from_the_page(kang_path):
    """Kang with its sample columns: QC suggests limits per sample; the page sends one stricter
    sample limit with Apply Filtering; only that sample loses more cells."""
    with test_server(APP, timeout_secs=240) as ts:
        ts.set_inputs(upload_data=uploaded(kang_path, "kang.h5ad"))
        ts.set_inputs(define_samples={"n": 1, "settings": {"sample_column": "replicate", "condition_column": "label"}})
        ts.set_inputs(confirm_metadata={"n": 1, "settings": {"species": "human", "tissues": ["Blood"]}})
        ts.set_inputs(run_qc={"n": 1, "settings": {"mito_prefix": "MT-"}})
        qc = output(ts, "qc")
        limits = qc["filter"]["current"]["sample_limits"]
        assert len(limits) == 16 and qc["samples"] == [s["sample"] for s in limits]
        rows = [{"sample": s["sample"], "min_genes": str(s["min_genes"]), "max_genes": str(s["max_n_genes_by_counts"] or ""),
                 "max_genes_off": s["max_n_genes_by_counts"] is None, "max_mt": "", "max_mt_off": True} for s in limits]
        rows[0]["min_genes"] = "700"
        boxes = {"min_genes": "0", "min_cells": "3", "max_genes": "", "max_genes_off": True, "max_mt": "", "max_mt_off": True,
                 "per_sample": True, "samples": rows}
        ts.set_inputs(apply_filter={"n": 1, "settings": boxes})
        assert workflow(ts)["notice"] is None, workflow(ts)["notice"]
        applied = output(ts, "qc")["applied"]["by_sample"]
        first = limits[0]["sample"]
        assert applied[first]["n_cells_kept"] < applied[first]["n_cells_before"] * 0.9  # 700 genes is strict
        assert output(ts, "qc")["filter"]["current"]["sample_limits"][0]["min_genes"] == 700


def test_several_samples_demo_runs_every_step(kang_path, fresh_reactive_lock):
    """The several-samples demo arrives as 16 samples with its species and tissue confirmed, and
    one click runs every step on the combined samples (QC limits per sample)."""
    async def check():
        async with test_server_async(APP, timeout_secs=600) as ts:
            await ts.set_inputs(load_demo_samples=1)
            data = workflow(ts)
            assert data["dataset"]["label"].startswith("Demo data: 16 samples")
            assert len(data["dataset"]["samples"]["list"]) == 16 and data["dataset"]["metadata"]["confirmed_by"] == "dataset"
            await ts.set_inputs(run_remaining=1)
            job = await wait_for_job(ts, timeout=500)
            assert job["status"] == "done", job
            data = workflow(ts)
            assert [s["status"] for s in data["steps"]] == ["done"] * 7
            assert data["steps"][0]["summary"].endswith("· 16 samples")
            assert len(output(ts, "qc")["filter"]["current"]["sample_limits"]) == 16  # limits per sample
            assert "clusters confident" in data["steps"][6]["summary"]  # Cell Type Annotation ran
    asyncio.run(check())
