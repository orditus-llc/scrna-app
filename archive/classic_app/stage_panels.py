"""Page layout: one tab per analysis stage.

Each function returns the layout of one tab: a short explanation, the stage's settings
(filled with the reference values from `params.py`), one action button, and empty spots
(`output_...`) where the server puts results. No logic here; see `app.py` for what the
buttons do.
"""

from __future__ import annotations

from shiny import ui

import display

F = display.default_filter_inputs()
S = display.default_stage_inputs()


def _advanced(*children):
    """A collapsed "Advanced" section for less common settings."""
    return ui.accordion(ui.accordion_panel("Advanced settings", *children), open=False)


def _results(*outputs):
    return ui.div(*outputs, class_="mt-3")


def qc_tab():
    filter_controls = ui.div(
        ui.input_numeric("min_genes", "Min genes per cell", F["min_genes"], min=0, step=10),
        ui.help_text("Removes near-empty droplets."),
        ui.input_numeric("min_cells", "Min cells per gene", F["min_cells"], min=0, step=1),
        ui.help_text("Removes genes detected in almost no cells."),
        ui.input_numeric("max_genes", "Max genes per cell", F["max_genes"], min=0, step=50),
        ui.input_checkbox("max_genes_off", "No limit", F["max_genes_off"]),
        ui.help_text("Unusually many genes suggests two cells in one droplet (doublet)."),
        ui.input_numeric("max_mt", "Max % mitochondrial", F["max_mt"], min=0, step=0.5),
        ui.input_checkbox("max_mt_off", "No limit", F["max_mt_off"]),
        ui.help_text("A high mitochondrial fraction often indicates damaged or dying cells."),
        ui.input_action_button("reset_filters", "Reset to reference settings",
                               class_="btn-outline-secondary btn-sm"),
    )
    return ui.nav_panel(
        "Data & QC",
        ui.card(
            ui.card_header("Load data"),
            ui.p("Opens the raw PBMC3k counts and checks they are raw counts."),
            ui.div(ui.input_action_button("load", "Load PBMC3k", class_="btn-primary")),
            ui.output_ui("load_info"),
        ),
        ui.card(
            ui.card_header("Quality control"),
            ui.p("Flags mitochondrial genes and calculates, for every cell, the number of genes "
                 "detected, total counts and % mitochondrial counts."),
            ui.div(ui.input_action_button("run_qc", "Run QC", class_="btn-primary")),
            ui.output_ui("qc_info"),
            ui.output_plot("violins", height="320px"),
        ),
        ui.card(
            ui.card_header("Filtering"),
            ui.p("Change a threshold to preview its effect. Nothing is removed until you click "
                 "Apply Filtering."),
            ui.layout_columns(
                filter_controls,
                ui.div(ui.output_ui("preview_info"), ui.output_plot("scatters", height="360px")),
                col_widths=(3, 9),
            ),
            ui.div(ui.input_action_button("apply_filter", "Apply Filtering", class_="btn-primary")),
            ui.output_ui("filter_result"),
        ),
    )


def preprocess_tab():
    return ui.nav_panel(
        "Preprocessing",
        ui.card(
            ui.card_header("Normalize, find highly variable genes, scale"),
            ui.p("Each cell is scaled to the same total counts and log-transformed. Highly "
                 "variable genes (HVGs) carry most of the biological signal; the log-normalized "
                 "data for all genes is kept for marker testing. Then each gene is scaled."),
            ui.layout_columns(
                ui.div(
                    ui.input_numeric("target_sum", "Target sum per cell", S["target_sum"], min=1, step=1000),
                    ui.input_numeric("hvg_min_mean", "HVG min mean", S["hvg_min_mean"], min=0, step=0.0025),
                    ui.input_numeric("hvg_max_mean", "HVG max mean", S["hvg_max_mean"], min=0, step=0.5),
                    ui.input_numeric("hvg_min_disp", "HVG min dispersion", S["hvg_min_disp"], min=0, step=0.1),
                    ui.input_checkbox("regress", "Regress out total counts and % mitochondrial", S["regress"]),
                    ui.help_text("Removes these technical effects from each gene. It can also remove "
                                 "real biology; your reference analysis uses it."),
                    _advanced(
                        ui.input_numeric("scale_max", "Scaling limit (clip at ± this)", S["scale_max"], min=0.1, step=1),
                    ),
                ),
                ui.div(ui.output_ui("preprocess_info"), ui.output_plot("hvg_plot", height="360px")),
                col_widths=(3, 9),
            ),
            ui.div(ui.input_action_button("run_preprocess", "Run Preprocessing", class_="btn-primary")),
        ),
    )


def reduce_tab():
    return ui.nav_panel(
        "PCA & UMAP",
        ui.card(
            ui.card_header("Principal components, neighbor graph and UMAP"),
            ui.p("PCA summarizes the HVGs as a few components. Each cell is then linked to its most "
                 "similar cells (the neighbor graph), which clustering and UMAP both use. UMAP is a "
                 "2-D picture for viewing only; distances on it aren't reliable measurements."),
            ui.layout_columns(
                ui.div(
                    ui.input_numeric("n_comps", "Number of PCs", S["n_comps"], min=2, step=5),
                    ui.input_numeric("n_neighbors", "Number of neighbors", S["n_neighbors"], min=2, step=1),
                    _advanced(
                        ui.input_numeric("n_pcs", "PCs used for the neighbor graph", S["n_pcs"], min=2, step=5),
                    ),
                ),
                ui.div(
                    ui.output_ui("reduce_info"),
                    ui.layout_columns(ui.output_plot("elbow_plot", height="320px"),
                                      ui.output_plot("umap_plain", height="380px")),
                ),
                col_widths=(3, 9),
            ),
            ui.div(ui.input_action_button("run_reduce", "Run PCA & UMAP", class_="btn-primary")),
        ),
    )


def cluster_tab():
    return ui.nav_panel(
        "Clustering",
        ui.card(
            ui.card_header("Leiden clustering"),
            ui.p("Finds groups of cells that are more connected to each other than to the rest. "
                 "The number of clusters comes from the data and the resolution: higher resolution "
                 "gives more, smaller clusters. Cluster numbers are labels only."),
            ui.layout_columns(
                ui.div(ui.input_numeric("resolution", "Resolution", S["resolution"], min=0.05, step=0.1)),
                ui.div(
                    ui.output_ui("cluster_info"),
                    ui.output_plot("umap_clusters", height="460px"),
                ),
                col_widths=(3, 9),
            ),
            ui.div(ui.input_action_button("run_cluster", "Run Clustering", class_="btn-primary")),
        ),
    )


def markers_tab():
    return ui.nav_panel(
        "Marker genes",
        ui.card(
            ui.card_header("Marker genes per cluster"),
            ui.p("Every gene in every cluster is compared with all other cells. Clusters were defined "
                 "from these same cells, so p-values are strongly inflated: use them to rank genes, "
                 "not as evidence of significance."),
            ui.layout_columns(
                ui.div(
                    ui.input_select("marker_method", "Test", {"wilcoxon": "Wilcoxon rank-sum",
                                                             "t-test": "t-test"}, selected=S["marker_method"]),
                ),
                ui.div(
                    ui.output_ui("markers_info"),
                    ui.output_plot("marker_dots", height="420px"),
                    ui.input_select("marker_cluster", "Show top genes for cluster", choices=[]),
                    ui.output_ui("marker_table"),
                ),
                col_widths=(3, 9),
            ),
            ui.div(ui.input_action_button("run_markers", "Find Marker Genes", class_="btn-primary")),
        ),
    )


def annotate_tab():
    marker_sets = ", ".join(
        f"{cell_type} ({' / '.join(genes)})"
        for cell_type, genes in display.P.annotation.marker_sets.items()
    )
    return ui.nav_panel(
        "Cell types",
        ui.card(
            ui.card_header("Predicted cell types"),
            ui.p("Each cluster is scored against marker gene sets. A label is proposed only if it "
                 "clearly beats the next candidate and is backed by the marker genes; otherwise the "
                 "cluster stays 'Unresolved'. These are predictions from marker evidence, not "
                 "analyst-reviewed labels."),
            ui.layout_columns(
                ui.div(
                    ui.p(ui.strong("Marker sets: "), marker_sets, class_="small"),
                    _advanced(
                        ui.input_numeric("ambiguity_margin", "Ambiguity margin", S["ambiguity_margin"], min=0, step=0.1),
                        ui.input_numeric("max_de_rank", "Top-N marker rank", S["max_de_rank"], min=1, step=5),
                        ui.input_numeric("conflict_z", "Conflict threshold (z)", S["conflict_z"], min=0, step=0.25),
                    ),
                ),
                ui.div(
                    ui.output_ui("annotate_info"),
                    ui.output_plot("umap_types", height="460px"),
                ),
                col_widths=(3, 9),
            ),
            ui.div(ui.input_action_button("run_annotate", "Predict Cell Types", class_="btn-primary")),
        ),
    )


ALL_TABS = (qc_tab, preprocess_tab, reduce_tab, cluster_tab, markers_tab, annotate_tab)
