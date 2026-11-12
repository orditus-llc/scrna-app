"""Single-cell analysis app: React interface (shinyreact) on top of the validated pipeline.

Build the interface once (and after every change to app/srcts/):
    conda activate scrna-ui
    cd app && npm run build
Then run it from the project folder:
    conda activate scrna
    shiny run app/app.py
and open http://127.0.0.1:8000.

The React code (app/srcts/) draws everything. This server only runs pipeline stages when
buttons are clicked and publishes the results as JSON (see payloads.py). All analysis
happens in the backend (scrna_analysis/).
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = APP_DIR.parents[1]  # snapshot: archive/general_app_2026-11-03/ -> the project folder (for data/)
if str(APP_DIR) not in sys.path:  # this snapshot's own `scrna_analysis` and `payloads`
    sys.path.insert(0, str(APP_DIR))

from shiny import App, reactive, render  # noqa: E402
from shinyreact import page_react, reactive_output  # noqa: E402

import payloads  # noqa: E402
import report  # noqa: E402
from scrna_analysis import export as export_tools  # noqa: E402
from scrna_analysis import marker_db as marker_db_tools  # noqa: E402
from scrna_analysis.params import DEFAULT_PARAMS  # noqa: E402
from scrna_analysis.pipeline import AnalysisRun  # noqa: E402

# The demo dataset ("Use demo data" on Data Overview). It is analyzed like any uploaded file.
DEMO_PATH = PROJECT_ROOT / "data" / "pbmc3k.h5ad"
DEMO_LABEL = "Demo data: human blood cells (PBMC3k, 10x Genomics)"
# What is known about the demo sample (it comes with the data, like a sample sheet): Cell Type
# Annotation starts with these, so a visitor can run the whole analysis without knowing the tissue.
# Uploaded files have no description, so users choose for their own data.
DEMO_DESCRIPTION = {"species": "human", "tissues": ["Immune system", "Blood"],
                    "about": "human blood mononuclear cells"}

app_ui = page_react(src_dir=APP_DIR / "www", title="Single-cell analysis")


def server(input, output, session):
    # This browser session's analysis; empty until data is loaded. DEFAULT_PARAMS: standard
    # starting values, with the QC thresholds suggested from each dataset's own cells.
    run = AnalysisRun(None, DEFAULT_PARAMS)
    dataset_label = None  # how the loaded data is named on screen
    changed = reactive.value(0)  # bumped after every action, so the page refreshes
    notice = reactive.value(None)  # latest message for the user (shown as a pop-up)
    job = reactive.value(None)  # the latest "Run all remaining steps" job (payloads.run_remaining_payload)
    marker_list = reactive.value(None)  # the user's own marker list (marker_db.MarkerList), if uploaded
    paused = {"before_next_stage": False}  # see _work_through_remaining

    def tell_user(text: str, step: str, kind: str = "error") -> None:
        """A message for the user: shown as a pop-up and in that step's card."""
        previous = notice.get()
        notice.set({"id": (previous["id"] + 1) if previous else 1, "kind": kind, "text": text, "step": step})

    def attempt(action, step: str) -> None:
        """Run one pipeline action for a Workflow step; turn any problem into a message
        instead of a crash."""
        try:
            action()
        except Exception as error:  # backend errors carry a plain-language message
            tell_user(str(error), step)
        changed.set(changed.get() + 1)

    def clicked(value) -> bool:
        """Buttons send a click counter, either alone or as {"n": count, "settings": {...}}.
        The page sends 0 when it first opens, so only counts above 0 are real clicks."""
        count = value.get("n", 0) if isinstance(value, dict) else value
        return (count or 0) > 0

    def sent_settings(value):
        """The settings a button sent with its click (the boxes as they were when clicked)."""
        return value.get("settings") if isinstance(value, dict) else None

    def run_with_settings(value, read_settings, action, step: str) -> None:
        """Check the settings a button sent, then run its stage (or explain the problem)."""
        try:
            settings = read_settings(sent_settings(value), run.params)
        except payloads.InputError as problem:
            tell_user(str(problem), step)
            return
        attempt(lambda: action(*settings), step=step)

    def filter_boxes():
        """The QC page's filter boxes (a dict), or None before the page has sent them."""
        return input.filter_settings() if "filter_settings" in input else None

    # --- Data Overview --------------------------------------------------------------------
    def open_dataset(path, file_name: str, label: str, description: dict | None = None) -> None:
        """Start a new analysis of a file and run its checks. The current analysis is replaced
        only if the new file passes; otherwise it stays as it was and the problem is shown."""
        nonlocal run, dataset_label
        current = job.get()
        if current and current["status"] == "running":
            tell_user("Wait for the running steps to finish before loading new data.", "data_overview")
            return
        new_run = AnalysisRun(path, DEFAULT_PARAMS, file_name=file_name, description=description)
        try:
            new_run.load()
        except Exception as error:
            tell_user(str(error), "data_overview")
        else:
            run, dataset_label = new_run, label
            job.set(None)
            marker_list.set(None)  # a new analysis starts without the previous dataset's list
        changed.set(changed.get() + 1)

    @reactive.effect
    @reactive.event(input.upload_data)
    def _upload_data():
        """A file uploaded on Data Overview (Shiny's upload: [{"name", "datapath", ...}])."""
        files = input.upload_data()
        if files:
            open_dataset(files[0]["datapath"], files[0]["name"], files[0]["name"])

    @reactive.effect
    @reactive.event(input.load_demo)
    def _load_demo():
        if clicked(input.load_demo()):
            open_dataset(DEMO_PATH, DEMO_PATH.name, DEMO_LABEL, DEMO_DESCRIPTION)

    # --- Quality Control ------------------------------------------------------------------
    @reactive.effect
    @reactive.event(input.run_qc)
    def _run_qc():
        if clicked(input.run_qc()):
            run_with_settings(input.run_qc(), payloads.qc_params_from_browser,
                              run.compute_qc_metrics, step="quality_control")

    @reactive.effect
    @reactive.event(input.apply_filter)
    def _apply_filter():
        if not clicked(input.apply_filter()):
            return
        value = input.apply_filter()
        boxes = sent_settings(value) if isinstance(value, dict) else filter_boxes()
        try:
            params = payloads.filter_params_from_browser(boxes, run.params.filter)
            # The MADs the suggestions were made with, if the page sent them (recorded only).
            sent_suggest = boxes.get("suggest") if isinstance(boxes, dict) else None
            suggest = payloads.suggest_params_from_browser(sent_suggest, run.params.suggest) if sent_suggest else None
        except payloads.InputError as problem:
            tell_user(str(problem), "quality_control")
            return
        attempt(lambda: run.apply_filter(params, suggest), step="quality_control")

    # --- Preprocessing and Dimensionality Reduction ----------------------------------------
    @reactive.effect
    @reactive.event(input.run_preprocess)
    def _run_preprocess():
        if clicked(input.run_preprocess()):
            run_with_settings(input.run_preprocess(), payloads.preprocess_params_from_browser,
                              run.preprocess, step="preprocessing")

    @reactive.effect
    @reactive.event(input.run_cluster)
    def _run_cluster():
        if clicked(input.run_cluster()):
            run_with_settings(input.run_cluster(), payloads.cluster_params_from_browser,
                              run.cluster, step="clustering")

    @reactive.effect
    @reactive.event(input.run_markers)
    def _run_markers():
        if clicked(input.run_markers()):
            run_with_settings(input.run_markers(), payloads.marker_params_from_browser,
                              run.find_markers, step="marker_genes")

    @reactive.effect
    @reactive.event(input.upload_markers)
    def _upload_markers():
        """The user's own marker list, uploaded on Cell Type Annotation (checked, kept for this session)."""
        files = input.upload_markers()
        if not files:
            return
        try:
            marker_list.set(marker_db_tools.read_marker_list(files[0]["datapath"], files[0]["name"]))
        except marker_db_tools.MarkerDatabaseError as problem:
            tell_user(str(problem), "cell_type_annotation")
        changed.set(changed.get() + 1)

    @reactive.effect
    @reactive.event(input.review_cluster)
    def _review_cluster():
        """The analyst confirms, renames or clears one cluster's label (Cell Type Annotation)."""
        value = input.review_cluster()
        if not clicked(value):
            return
        try:
            cluster, decision, label, note = payloads.review_params_from_browser(sent_settings(value))
        except payloads.InputError as problem:
            tell_user(str(problem), "cell_type_annotation")
            return
        attempt(lambda: run.review_cluster(cluster, decision, label, note), step="cell_type_annotation")

    @reactive.effect
    @reactive.event(input.run_annotation)
    def _run_annotation():
        if clicked(input.run_annotation()):
            run_with_settings(input.run_annotation(),
                              lambda values, params: payloads.annotation_params_from_browser(values, params, marker_list.get()),
                              run.annotate, step="cell_type_annotation")

    @reactive.effect
    @reactive.event(input.run_reduction)
    def _run_reduction():
        if clicked(input.run_reduction()):
            run_with_settings(input.run_reduction(), payloads.reduction_params_from_browser,
                              run.reduce_dimensions, step="dimensionality_reduction")

    # --- Run all remaining steps -----------------------------------------------------------
    @reactive.effect
    @reactive.event(input.run_remaining)
    def _start_run_remaining():
        if not clicked(input.run_remaining()):
            return
        current = job.get()
        if current and current["status"] == "running":
            return  # already running
        stages = payloads.remaining_stages(run)
        if not stages:
            return
        job.set({"id": (current["id"] + 1) if current else 1, "stages": tuple(stages),
                 "left": tuple(stages), "status": "running"})

    @reactive.effect
    def _work_through_remaining():
        """Run the job's stages one at a time. Before each stage it pauses briefly, so the
        page receives the progress so far (Shiny sends updates only between runs of the
        server's code)."""
        current = job.get()
        if not current or current["status"] != "running":
            return
        if not paused["before_next_stage"]:
            paused["before_next_stage"] = True
            reactive.invalidate_later(0.05)
            return
        paused["before_next_stage"] = False

        stage, left = current["left"][0], current["left"][1:]
        needed = payloads.choice_needed(run, stage)
        if needed:  # a planned stop: this step waits for the user's choice
            n = payloads.steps_done_before(current, stage)
            done_text = f"Ran {n} {'step' if n == 1 else 'steps'}. " if n else ""
            tell_user(done_text + needed, "run_remaining", kind="info")
            job.set({**current, "status": "stopped", "stopped_step": payloads.step_of_stage(stage)[0]})
            changed.set(changed.get() + 1)
            return
        try:
            payloads.run_stage(run, stage)
        except Exception as error:  # stop here; the finished stages are kept
            step_key = payloads.step_of_stage(stage)[0]
            tell_user(str(error), step_key)
            job.set({**current, "status": "failed", "failed_step": step_key})
        else:
            if left:
                job.set({**current, "left": left})
            else:
                job.set({**current, "left": (), "status": "done"})
                n = payloads.run_remaining_payload(current)["total"]  # steps run
                tell_user(f"Ran {n} {'step' if n == 1 else 'steps'} with the current settings.",
                          "run_remaining", kind="info")
        changed.set(changed.get() + 1)

    # --- Downloads: the analysed data (.h5ad) and the run report (HTML) ----------------------
    # Served by Shiny at session/<id>/download/<name>; the page links there (lib/serverDownload.ts)
    # once `export.available` says the analysis can be exported.
    @render.download_button(filename=lambda: payloads.export_file_name(run.file_name, "analysed", "h5ad"),
                            media_type="application/x-hdf5")
    def download_h5ad():
        data = export_tools.analysed_anndata(run)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "analysed.h5ad"
            data.write_h5ad(path, compression="gzip")
            yield path.read_bytes()

    @render.download_button(filename=lambda: payloads.export_file_name(run.file_name, "run report", "html"),
                            media_type="text/html")
    def download_report():
        yield report.run_report_html(run, dataset_label).encode("utf-8")

    # --- Outputs (JSON for the React interface) -------------------------------------------
    @reactive_output
    def workflow():
        changed()
        return payloads.workflow_payload(run, notice.get(), job.get(), dataset_label)

    @reactive_output
    def qc():
        changed()
        return payloads.qc_payload(run)

    @reactive_output
    def preprocessing():
        changed()
        return payloads.preprocess_payload(run)

    @reactive_output
    def reduction():
        changed()
        return payloads.reduction_payload(run)

    @reactive_output
    def clustering():
        changed()
        return payloads.cluster_payload(run)

    @reactive_output
    def markers():
        changed()
        return payloads.markers_payload(run)

    @reactive_output
    def annotation():
        changed()
        return payloads.annotation_payload(run, marker_list.get())

    @reactive_output
    def explore():
        changed()
        return payloads.explore_payload(run)

    @reactive_output
    def explore_gene_values():
        changed()
        gene = input.explore_gene() if "explore_gene" in input else None
        return payloads.gene_payload(run, gene)

    @reactive_output
    def qc_preview():
        changed()
        return payloads.preview_payload(run, filter_boxes(), suggest_boxes())

    def suggest_boxes():
        """The "MADs" boxes under More settings on the QC page, or None before they're sent."""
        return input.suggest_settings() if "suggest_settings" in input else None

    @reactive_output
    def qc_suggestion():
        """Thresholds suggested for the MADs in the boxes (the page fills the filter boxes with them)."""
        changed()
        return payloads.suggestion_payload(run, suggest_boxes())

    @reactive_output
    def mito_check():
        changed()
        prefix = input.mito_prefix() if "mito_prefix" in input else None
        return payloads.mito_check_payload(run, prefix)


app = App(app_ui, server)
