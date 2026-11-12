"""Classic PBMC3k app (Milestone 4): every stage on its own tab, built-in Shiny UI.

Kept as a fallback while the React interface (app/app.py) is built; retire in step 5g.

Start it from the project folder:
    conda activate scrna
    shiny run app/classic/app.py
then open the address it prints (http://127.0.0.1:8000) in a browser.

The layout lives in `stage_panels.py`. This file is the "server": it connects each button
to a pipeline stage (`pbmc_analysis.pipeline.AnalysisRun`) and fills the results areas.
All analysis happens in the backend; this file only calls it and displays the results.
"""

from __future__ import annotations

import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = APP_DIR.parents[1]  # app/classic/ -> project folder
for folder in (PROJECT_ROOT, APP_DIR):  # so `pbmc_analysis`, `display`, `plots`, ... can be found
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))

from shiny import App, reactive, render, ui  # noqa: E402

import display  # noqa: E402
import plots  # noqa: E402
import stage_panels  # noqa: E402
from pbmc_analysis.pipeline import STAGES, AnalysisRun  # noqa: E402

DATA_PATH = PROJECT_ROOT / "data" / "pbmc3k.h5ad"
FILTER_DEFAULTS = display.default_filter_inputs()

# Which button runs each stage (used in "click ... to run this stage" messages).
RUN_BUTTON = {
    "load": "Load PBMC3k", "qc_metrics": "Run QC", "filter": "Apply Filtering",
    "preprocess": "Run Preprocessing", "reduce_dimensions": "Run PCA & UMAP",
    "cluster": "Run Clustering", "markers": "Find Marker Genes", "annotate": "Predict Cell Types",
}


# --- Layout -------------------------------------------------------------------------------

app_ui = ui.page_fluid(
    ui.h2("PBMC3k analysis", class_="mt-3"),
    ui.p("A guided Scanpy workflow using the validated analysis backend.", class_="text-muted"),
    ui.output_ui("progress"),
    ui.navset_card_underline(*[tab() for tab in stage_panels.ALL_TABS], id="stage_tabs"),
)


# --- Server -------------------------------------------------------------------------------

def server(input, output, session):
    run = AnalysisRun(DATA_PATH)  # this browser session's analysis (reference settings)
    changed = reactive.value(0)  # bumped after every stage, so displays refresh

    def attempt(action, message: str = "Running…") -> None:
        """Run one pipeline action with a progress indicator; show problems as messages."""
        with ui.Progress(min=0, max=1) as progress:
            progress.set(message=message)
            try:
                action()
            except Exception as error:  # backend errors carry a message for the user
                ui.notification_show(str(error), type="error", duration=12)
        changed.set(changed() + 1)

    def run_with_settings(build_settings, stage_method, message: str) -> None:
        """Translate the boxes into settings (or explain what's wrong), then run the stage."""
        try:
            settings = build_settings()
        except display.InputError as problem:
            ui.notification_show(str(problem), type="error", duration=12)
            return
        attempt(lambda: stage_method(**settings), message)

    def not_ready(stage: str):
        """A note explaining why a stage has no results, or None if it's done."""
        status = run.status()[stage]
        if status == "done":
            return None
        if status == "needs re-running":
            text = (f"An earlier stage was redone, so these results were cleared. "
                    f"Click {RUN_BUTTON[stage]} to run this stage again.")
            return ui.div(text, class_="alert alert-warning")
        return ui.p(f"Not run yet. Click {RUN_BUTTON[stage]}.", class_="text-muted")

    def is_done(stage: str) -> bool:
        changed()
        return run.status()[stage] == "done"

    # --- Buttons -> pipeline stages ------------------------------------------------------

    @reactive.effect
    @reactive.event(input.load)
    def _load():
        attempt(run.load, "Loading data…")

    @reactive.effect
    @reactive.event(input.run_qc)
    def _run_qc():
        attempt(run.compute_qc_metrics, "Calculating QC metrics…")

    @reactive.effect
    @reactive.event(input.apply_filter)
    def _apply_filter():
        params, problem = current_filter_settings()
        if problem:
            ui.notification_show(problem, type="error", duration=12)
            return
        attempt(lambda: run.apply_filter(params), "Filtering…")

    @reactive.effect
    @reactive.event(input.reset_filters)
    def _reset_filters():
        for box in ("min_genes", "min_cells", "max_genes", "max_mt"):
            ui.update_numeric(box, value=FILTER_DEFAULTS[box])
        for box in ("max_genes_off", "max_mt_off"):
            ui.update_checkbox(box, value=FILTER_DEFAULTS[box])

    @reactive.effect
    @reactive.event(input.run_preprocess)
    def _run_preprocess():
        run_with_settings(lambda: display.preprocess_settings(
            run.params, input.target_sum(), input.hvg_min_mean(), input.hvg_max_mean(),
            input.hvg_min_disp(), input.regress(), input.scale_max()),
            run.preprocess, "Normalizing, finding HVGs, regressing and scaling…")

    @reactive.effect
    @reactive.event(input.run_reduce)
    def _run_reduce():
        run_with_settings(lambda: display.reduce_settings(
            run.params, input.n_comps(), input.n_neighbors(), input.n_pcs()),
            run.reduce_dimensions, "Running PCA, neighbor graph and UMAP…")

    @reactive.effect
    @reactive.event(input.run_cluster)
    def _run_cluster():
        run_with_settings(lambda: display.cluster_settings(run.params, input.resolution()),
                          run.cluster, "Clustering…")

    @reactive.effect
    @reactive.event(input.run_markers)
    def _run_markers():
        run_with_settings(lambda: display.marker_settings(run.params, input.marker_method()),
                          run.find_markers, "Testing marker genes…")

    @reactive.effect
    @reactive.event(input.run_annotate)
    def _run_annotate():
        run_with_settings(lambda: display.annotation_settings(
            run.params, input.ambiguity_margin(), input.max_de_rank(), input.conflict_z()),
            run.annotate, "Predicting cell types…")

    @reactive.effect
    def _marker_cluster_choices():
        """Offer the clusters that actually exist once marker genes are available."""
        if is_done("markers"):
            groups = run.stage_result("markers").summary["groups"]
            ui.update_select("marker_cluster", choices=groups, selected=groups[0])

    # --- Shared: filter settings and live preview ----------------------------------------

    @reactive.calc
    def current_filter_settings():
        try:
            params = display.filter_params_from_inputs(
                input.min_genes(), input.min_cells(), input.max_genes(), input.max_genes_off(),
                input.max_mt(), input.max_mt_off(),
            )
            return params, None
        except display.InputError as problem:
            return None, str(problem)

    @reactive.calc
    def filter_preview():
        if not is_done("qc_metrics"):
            return None, None
        params, problem = current_filter_settings()
        if problem:
            return None, problem
        try:
            return run.preview_filter(params), None
        except Exception as error:
            return None, str(error)

    # --- Top of the page -----------------------------------------------------------------

    @render.ui
    def progress():
        changed()
        badge = {"done": "bg-success", "needs re-running": "bg-warning text-dark",
                 "not run": "bg-light text-dark border"}
        status = run.status()
        settings_badge = (
            ui.span("Using the validated PBMC3k settings", class_="badge bg-info text-dark ms-2")
            if display.uses_reference_settings(run.params) else
            ui.span("Settings differ from the validated reference", class_="badge bg-warning text-dark ms-2")
        )
        return ui.div(
            *[ui.span(f"{display.STAGE_LABELS[s]}: {status[s]}", class_=f"badge {badge[status[s]]} me-1 mb-1")
              for s in STAGES],
            settings_badge,
            class_="mb-3",
        )

    # --- Data & QC tab -------------------------------------------------------------------

    @render.ui
    def load_info():
        if not is_done("load"):
            return not_ready("load")
        summary = run.stage_result("load").summary
        return ui.div(
            ui.layout_columns(
                ui.value_box("Cells", f"{summary['n_cells']:,}"),
                ui.value_box("Genes", f"{summary['n_genes']:,}"),
                ui.value_box("Duplicate gene names fixed", f"{summary['n_duplicate_gene_names']:,}"),
            ),
            ui.p(f"File: {DATA_PATH.name} (fingerprint {display.short_fingerprint(summary['input_sha256'])})",
                 class_="small text-muted"),
        )

    @render.ui
    def qc_info():
        if not is_done("qc_metrics"):
            return not_ready("qc_metrics")
        summary = run.stage_result("qc_metrics").details[0]
        parts = [ui.div(w, class_="alert alert-warning") for w in summary.warnings]
        parts.append(ui.p(f"{summary.n_mito_genes} mitochondrial genes: {', '.join(summary.mito_genes)}"))
        parts.append(ui.HTML(display.html_table(display.qc_stats_table(summary))))
        return ui.div(*parts)

    @render.plot
    def violins():
        if not is_done("qc_metrics"):
            return None
        summary = run.stage_result("qc_metrics").details[0]
        return plots.qc_violins(run.result("qc_metrics").obs, summary.mito_qc_available)

    @render.ui
    def preview_info():
        preview, problem = filter_preview()
        if problem:
            return ui.div(problem, class_="alert alert-danger")
        if preview is None:
            return ui.p("Run QC to preview filtering.", class_="text-muted")
        return ui.div(
            ui.p(ui.strong(f"Would keep {preview.n_cells_kept:,} of {preview.n_cells_before:,} cells "
                           f"and {preview.n_genes_kept:,} of {preview.n_genes_before:,} genes.")),
            ui.HTML(display.html_table(display.removal_table(preview))),
        )

    @render.plot
    def scatters():
        preview, _ = filter_preview()
        if not is_done("qc_metrics"):
            return None
        params, _ = current_filter_settings()
        kept = preview.kept_cells if preview is not None else None
        return plots.qc_scatters(run.result("qc_metrics").obs, kept, params)

    @render.ui
    def filter_result():
        if not is_done("filter"):
            return not_ready("filter")
        summary = run.stage_result("filter").summary
        return ui.div(
            ui.strong(f"Filtered: {summary['n_cells_kept']:,} cells x {summary['n_genes_kept']:,} genes."),
            ui.span(" Recorded in the analysis history. Next: the Preprocessing tab.", class_="text-muted"),
            class_="alert alert-success mt-2",
        )

    # --- Preprocessing tab ---------------------------------------------------------------

    @render.ui
    def preprocess_info():
        if not is_done("preprocess"):
            return not_ready("preprocess")
        s = run.stage_result("preprocess").summary
        regression = ("Regressed out: " + ", ".join(s["regress_keys"])) if s["regressed"] else "No regression"
        return ui.div(
            ui.layout_columns(
                ui.value_box("Highly variable genes", f"{s['n_hvg']:,}"),
                ui.value_box("Genes kept for markers (raw)", f"{s['raw_shape'][1]:,}"),
                ui.value_box("Median counts per cell before normalizing", f"{s['median_counts_before']:,.0f}"),
            ),
            ui.p(f"{regression}. Values clipped at the scaling limit: {s['n_values_at_clip_limit']:,}.",
                 class_="small text-muted"),
        )

    @render.plot
    def hvg_plot():
        if not is_done("preprocess"):
            return None
        return plots.hvg_plot(run.result("preprocess"))

    # --- PCA & UMAP tab ------------------------------------------------------------------

    @render.ui
    def reduce_info():
        if not is_done("reduce_dimensions"):
            return not_ready("reduce_dimensions")
        s = run.stage_result("reduce_dimensions").summary
        return ui.layout_columns(
            ui.value_box(f"Variance explained by {s['n_comps']} PCs", f"{100 * s['total_variance_ratio']:.1f}%"),
            ui.value_box("Neighbors per cell", f"{s['n_neighbors']}"),
            ui.value_box("Graph connections", f"{s['n_connections']:,}"),
        )

    @render.plot
    def elbow_plot():
        if not is_done("reduce_dimensions"):
            return None
        return plots.pca_variance_plot(run.result("reduce_dimensions"))

    @render.plot
    def umap_plain():
        if not is_done("reduce_dimensions"):
            return None
        return plots.umap_plot(run.result("reduce_dimensions"), title="UMAP")

    # --- Clustering tab ------------------------------------------------------------------

    @render.ui
    def cluster_info():
        if not is_done("cluster"):
            return not_ready("cluster")
        s = run.stage_result("cluster").summary
        return ui.div(
            ui.p(ui.strong(f"{s['n_clusters']} clusters"),
                 f" at resolution {run.params.leiden.resolution}."),
            ui.HTML(display.html_table(display.cluster_sizes_table(s["cluster_sizes"]))),
        )

    @render.plot
    def umap_clusters():
        if not is_done("cluster"):
            return None
        return plots.umap_plot(run.result("cluster"), color="leiden", title="Leiden clusters",
                               labels_on_plot=True)

    # --- Marker genes tab ----------------------------------------------------------------

    @render.ui
    def markers_info():
        if not is_done("markers"):
            return not_ready("markers")
        s = run.stage_result("markers").summary
        return ui.p(f"{s['n_genes_tested']:,} genes tested in each of {len(s['groups'])} clusters "
                    f"({run.params.markers.method}). Below: top 3 markers per cluster.",
                    class_="small text-muted")

    @render.plot
    def marker_dots():
        if not is_done("markers"):
            return None
        return plots.marker_dotplot(run.result("markers"), n_genes=3)

    @render.ui
    def marker_table():
        if not is_done("markers"):
            return None
        cluster = input.marker_cluster()
        if not cluster or cluster not in run.stage_result("markers").summary["groups"]:
            return None
        table = display.marker_genes_table(run.result("markers"), cluster, n_top=10)
        return ui.HTML(display.html_table(table))

    # --- Cell types tab ------------------------------------------------------------------

    @render.ui
    def annotate_info():
        if not is_done("annotate"):
            return not_ready("annotate")
        details = run.stage_result("annotate").details[0]
        notes = [ui.div(note, class_="alert alert-warning") for note in details.notes]
        return ui.div(
            *notes,
            ui.HTML(display.html_table(display.annotation_table(details.clusters))),
            ui.p("Predicted from marker evidence. Not analyst-reviewed and not a claim of "
                 "biological correctness.", class_="small text-muted"),
        )

    @render.plot
    def umap_types():
        if not is_done("annotate"):
            return None
        return plots.umap_plot(run.result("annotate"), color="cell_type_predicted",
                               title="Predicted cell types")


app = App(app_ui, server)
