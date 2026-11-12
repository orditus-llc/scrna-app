"""FROZEN SNAPSHOT (2026-10-22) of the PBMC3k app, kept for demos. See README.md.

Run it from the project folder (any port that's free):
    conda activate scrna
    shiny run archive/pbmc3k_app/app.py --port 8002
It uses its own copy of the analysis package (./pbmc_analysis) and the project's
data/pbmc3k.h5ad. It is not updated or tested; the main app is app/app.py.

Original description:
PBMC3k app: React interface (shinyreact) on top of the validated pipeline.

Build the interface once (and after every change to app/srcts/):
    conda activate scrna-ui
    cd app && npm run build
Then run it from the project folder:
    conda activate scrna
    shiny run app/app.py
and open http://127.0.0.1:8000.

The React code (app/srcts/) draws everything. This server only runs pipeline stages when
buttons are clicked and publishes the results as JSON (see payloads.py). All analysis
happens in the backend (pbmc_analysis/).
"""

from __future__ import annotations

import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = APP_DIR.parents[1]  # archive/pbmc3k_app -> project folder (for data/)
if str(APP_DIR) not in sys.path:  # this snapshot's own `pbmc_analysis` and `payloads`
    sys.path.insert(0, str(APP_DIR))

from shiny import App, reactive  # noqa: E402
from shinyreact import page_react, reactive_output  # noqa: E402

import payloads  # noqa: E402
from pbmc_analysis.pipeline import AnalysisRun  # noqa: E402

DATA_PATH = PROJECT_ROOT / "data" / "pbmc3k.h5ad"

app_ui = page_react(src_dir=APP_DIR / "www", title="PBMC3k analysis")


def server(input, output, session):
    run = AnalysisRun(DATA_PATH)  # this browser session's analysis (reference settings)
    changed = reactive.value(0)  # bumped after every action, so the page refreshes
    notice = reactive.value(None)  # latest message for the user (shown as a pop-up)
    job = reactive.value(None)  # the latest "Run all remaining steps" job (payloads.run_remaining_payload)
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
    @reactive.effect
    @reactive.event(input.load_data)
    def _load_data():
        if clicked(input.load_data()):
            attempt(run.load, step="data_overview")

    # --- Quality Control ------------------------------------------------------------------
    @reactive.effect
    @reactive.event(input.run_qc)
    def _run_qc():
        if clicked(input.run_qc()):
            attempt(run.compute_qc_metrics, step="quality_control")

    @reactive.effect
    @reactive.event(input.apply_filter)
    def _apply_filter():
        if not clicked(input.apply_filter()):
            return
        value = input.apply_filter()
        boxes = sent_settings(value) if isinstance(value, dict) else filter_boxes()
        try:
            params = payloads.filter_params_from_browser(boxes, run.params.filter)
        except payloads.InputError as problem:
            tell_user(str(problem), "quality_control")
            return
        attempt(lambda: run.apply_filter(params), step="quality_control")

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
    @reactive.event(input.run_annotation)
    def _run_annotation():
        if clicked(input.run_annotation()):
            run_with_settings(input.run_annotation(), payloads.annotation_params_from_browser,
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

    # --- Outputs (JSON for the React interface) -------------------------------------------
    @reactive_output
    def workflow():
        changed()
        return payloads.workflow_payload(run, notice.get(), job.get())

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
        return payloads.annotation_payload(run)

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
        return payloads.preview_payload(run, filter_boxes())


app = App(app_ui, server)
