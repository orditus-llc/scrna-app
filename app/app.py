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
PROJECT_ROOT = APP_DIR.parent
for folder in (PROJECT_ROOT, APP_DIR):  # so `scrna_analysis` and `payloads` can be found
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))

from shiny import App, reactive, render  # noqa: E402
from shinyreact import page_react, reactive_output  # noqa: E402

import payloads  # noqa: E402
import report  # noqa: E402
from scrna_analysis import export as export_tools  # noqa: E402
from scrna_analysis.params import DEFAULT_PARAMS  # noqa: E402
from scrna_analysis.pipeline import AnalysisRun  # noqa: E402

# The demo dataset ("Use demo data" on Data Overview), with its metadata: what is known about the
# sample, like a sample sheet that comes with the data. It is analyzed like any uploaded file; its
# metadata simply arrives confirmed. Uploaded files have none, so users confirm species and tissue.
DEMO_METADATA = {
    "file": PROJECT_ROOT / "data" / "pbmc3k.h5ad",
    "label": "Demo data: human blood cells (PBMC3k, 10x Genomics)",
    "species": "human",
    "tissue": "blood",  # the Census tissue
    "tissues": ["Immune system", "Blood"],  # PanglaoDB's tissue names (as the Census tissue maps them)
    "about": "Peripheral blood mononuclear cells from a healthy donor (PBMC3k, 10x Genomics).",
}
DEMO_PATH = DEMO_METADATA["file"]
# The several-samples demo: Kang et al. 2018 (16 samples: 8 patients, control vs interferon-beta).
# Downloaded by scripts/download_kang_data.py; its sample and condition columns come with it.
DEMO_SAMPLES = {
    "file": PROJECT_ROOT / "data" / "kang_2018.h5ad",
    "label": "Demo data: 16 samples of human blood cells, control vs interferon-beta (Kang et al. 2018)",
    "species": "human",
    "tissue": "blood",
    "tissues": ["Immune system", "Blood"],
    "about": "Peripheral blood mononuclear cells from 8 lupus patients, each cultured 6 hours with or without "
             "interferon-beta (Kang et al. 2018, Nature Biotechnology).",
    "sample_columns": {"sample": "replicate", "condition": "label"},
}

app_ui = page_react(src_dir=APP_DIR / "www", title="Single-cell analysis")


def server(input, output, session):
    # This browser session's analysis; empty until data is loaded. DEFAULT_PARAMS: standard
    # starting values, with the QC thresholds suggested from each dataset's own cells.
    run = AnalysisRun(None, DEFAULT_PARAMS, require_metadata=True)
    dataset_label = None  # how the loaded data is named on screen
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
    def open_dataset(path, file_name: str | None, label: str, metadata: dict | None = None,
                     sample_files: list | None = None, sample_columns: dict | None = None) -> None:
        """Start a new analysis of a file (or of several samples' files) and run its checks. The
        current analysis is replaced only if the new data passes; otherwise it stays as it was and
        the problem is shown."""
        nonlocal run, dataset_label
        current = job.get()
        if current and current["status"] == "running":
            tell_user("Wait for the running steps to finish before loading new data.", "data_overview")
            return
        new_run = AnalysisRun(path, DEFAULT_PARAMS, file_name=file_name, metadata=metadata, require_metadata=True,
                              sample_files=sample_files, sample_columns=sample_columns)
        try:
            new_run.load()
        except Exception as error:
            tell_user(str(error), "data_overview")
        else:
            run, dataset_label = new_run, label
            job.set(None)
        changed.set(changed.get() + 1)

    @reactive.effect
    @reactive.event(input.upload_data)
    def _upload_data():
        """A file uploaded on Data Overview (Shiny's upload: [{"name", "datapath", ...}])."""
        files = input.upload_data()
        if files:
            open_dataset(files[0]["datapath"], files[0]["name"], files[0]["name"])

    @reactive.effect
    @reactive.event(input.load_demo_samples)
    def _load_demo_samples():
        if not clicked(input.load_demo_samples()):
            return
        if not DEMO_SAMPLES["file"].is_file():
            tell_user("The several-samples demo isn't downloaded on this computer. Run scripts/download_kang_data.py.",
                      "data_overview")
            return
        open_dataset(DEMO_SAMPLES["file"], DEMO_SAMPLES["file"].name, DEMO_SAMPLES["label"],
                     {k: DEMO_SAMPLES[k] for k in ("species", "tissue", "tissues", "about")},
                     sample_columns=DEMO_SAMPLES["sample_columns"])

    @reactive.effect
    @reactive.event(input.load_samples)
    def _load_samples():
        """Several samples: the page uploads their files (input `upload_samples`, in order), then
        sends each one's name and condition with this click."""
        value = input.load_samples()
        if not clicked(value):
            return
        files = input.upload_samples() if "upload_samples" in input else None
        try:
            sample_files = payloads.sample_files_from_browser(sent_settings(value), files)
        except payloads.InputError as problem:
            tell_user(str(problem), "data_overview")
            return
        open_dataset(None, None, f"{len(sample_files)} samples", sample_files=sample_files)

    @reactive.effect
    @reactive.event(input.define_samples)
    def _define_samples():
        """One loaded file: use one of its columns as the samples (and another as conditions), or not."""
        value = input.define_samples()
        if not clicked(value):
            return
        settings = sent_settings(value) or {}
        attempt(lambda: run.define_samples(settings.get("sample_column") or None, settings.get("condition_column") or None),
                step="data_overview")

    @reactive.effect
    @reactive.event(input.load_demo)
    def _load_demo():
        if clicked(input.load_demo()):
            open_dataset(DEMO_PATH, DEMO_PATH.name, DEMO_METADATA["label"],
                         {k: DEMO_METADATA[k] for k in ("species", "tissue", "tissues", "about")})

    @reactive.effect
    @reactive.event(input.confirm_metadata)
    def _confirm_metadata():
        """Species and tissue (and an optional description) confirmed on Data Overview."""
        value = input.confirm_metadata()
        if not clicked(value):
            return
        try:
            species, tissue, tissues, about = payloads.metadata_from_browser(sent_settings(value))
        except payloads.InputError as problem:
            tell_user(str(problem), "data_overview")
            return
        attempt(lambda: run.confirm_metadata(species, tissues, about, tissue=tissue), step="data_overview")

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
        if payloads.needs_confirmation(run):
            tell_user(payloads.CONFIRM_FIRST, "data_overview")
            return
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

    # --- Downloads: the analyzed data (.h5ad) and the run report (HTML) ----------------------
    # Served by Shiny at session/<id>/download/<name>; the page links there (lib/serverDownload.ts)
    # once `export.available` says the analysis can be exported.
    @render.download_button(filename=lambda: payloads.export_file_name(run.file_name, "analyzed", "h5ad"),
                            media_type="application/x-hdf5")
    def download_h5ad():
        data = export_tools.analyzed_anndata(run)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "analyzed.h5ad"
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
