"""Tests for the Milestone 4 app's helpers, plots and wiring.

The analysis itself is tested in the backend tests; here we check that the app
translates inputs and results correctly, draws the right plots, and loads.
"""

from __future__ import annotations

import importlib.util
import sys

import numpy as np
import pandas as pd
import pytest

from tests.conftest import PROJECT_ROOT

APP_DIR = PROJECT_ROOT / "app" / "classic"  # Milestone 4 app, kept until step 5g
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import display  # noqa: E402
import plots  # noqa: E402
from pbmc_analysis.params import REFERENCE_PRESET, FilterParams  # noqa: E402
from pbmc_analysis.pipeline import AnalysisRun  # noqa: E402


@pytest.fixture(scope="module")
def qc_run(raw_data_path) -> AnalysisRun:
    run = AnalysisRun(raw_data_path)
    run.load()
    run.compute_qc_metrics()
    return run


# --- Translating the filter boxes -------------------------------------------------------

def test_default_boxes_are_the_reference_settings():
    d = display.default_filter_inputs()
    params = display.filter_params_from_inputs(
        d["min_genes"], d["min_cells"], d["max_genes"], d["max_genes_off"], d["max_mt"], d["max_mt_off"]
    )
    assert params == REFERENCE_PRESET.filter


def test_no_limit_boxes_switch_limits_off():
    params = display.filter_params_from_inputs(200, 3, 2500, True, 6, True)
    assert params.max_n_genes_by_counts is None and params.max_pct_counts_mt is None


@pytest.mark.parametrize("values, message", [
    ((None, 3, 2500, False, 6, False), "Enter a value for 'Min genes per cell'"),
    ((200, -1, 2500, False, 6, False), "'Min cells per gene' can't be negative"),
    ((200.5, 3, 2500, False, 6, False), "must be a whole number"),
    ((200, 3, 2500, False, None, False), "Enter a value for 'Max % mitochondrial'"),
])
def test_unusable_boxes_give_a_message(values, message):
    with pytest.raises(display.InputError, match=message):
        display.filter_params_from_inputs(*values)


def test_empty_box_is_fine_when_no_limit_is_ticked():
    params = display.filter_params_from_inputs(200, 3, None, True, 6, False)
    assert params.max_n_genes_by_counts is None


# --- Translating results for display ----------------------------------------------------

def test_qc_table_shows_the_backend_numbers(qc_run):
    summary = qc_run.stage_result("qc_metrics").details[0]
    table = display.qc_stats_table(summary).set_index("Statistic")
    for metric, label in display.METRIC_LABELS.items():
        assert table.loc["Median", label] == round(summary.metric_stats[metric]["50%"], 2)
        assert table.loc["Max", label] == round(summary.metric_stats[metric]["max"], 2)


def test_removal_table_follows_the_preview(qc_run):
    preview = qc_run.preview_filter()
    table = display.removal_table(preview)
    assert table["Removed"].tolist() == list(preview.removed_in_order.values())
    assert table["Step"].tolist() == list(range(1, len(table) + 1))


def test_html_table_aligns_each_header_with_its_column(qc_run):
    """Regression test for the 2026-10-14 bug: headers right-aligned, values left-aligned,
    so in Safari every value appeared under the wrong header."""
    import re
    summary = qc_run.stage_result("qc_metrics").details[0]
    table = display.qc_stats_table(summary)
    html = display.html_table(table)
    header_classes = re.findall(r'<th class="([^"]+)">', html)
    first_row = re.search(r"<tbody><tr>(.*?)</tr>", html).group(1)
    cell_classes = re.findall(r'<td class="([^"]+)">', first_row)
    assert header_classes == cell_classes  # header and column aligned the same way
    assert header_classes == ["text-start", "text-end", "text-end", "text-end"]
    assert "text-align: right" not in html  # no pandas header styling


def test_html_table_number_format():
    table = pd.DataFrame({"Statistic": ["Cells", "Median"], "Value": [2700.0, 2.0296]})
    html = display.html_table(table)
    assert ">2,700<" in html and ">2.03<" in html and "2700.00" not in html


def test_settings_badge(qc_run):
    assert display.uses_reference_settings(REFERENCE_PRESET)
    changed = REFERENCE_PRESET.__class__(filter=FilterParams(min_genes=300))
    assert not display.uses_reference_settings(changed)
    assert display.short_fingerprint("abcdef0123456789") == "abcdef012345"
    assert display.short_fingerprint(None) == "unknown"


# --- Plots ------------------------------------------------------------------------------

def test_violins_have_three_panels(qc_run):
    fig = plots.qc_violins(qc_run.result("qc_metrics").obs)
    assert [ax.get_title() for ax in fig.axes] == [
        "Genes per cell", "Counts per cell", "% mitochondrial counts"]


def test_violins_say_when_mito_is_unavailable(qc_run):
    fig = plots.qc_violins(qc_run.result("qc_metrics").obs, mito_available=False)
    assert "Unavailable" in fig.axes[2].texts[0].get_text()


def test_scatters_grey_out_exactly_the_cells_the_preview_removes(qc_run):
    obs = qc_run.result("qc_metrics").obs
    preview = qc_run.preview_filter()
    fig = plots.qc_scatters(obs, preview.kept_cells, REFERENCE_PRESET.filter)
    removed = preview.n_cells_before - preview.n_cells_kept
    for ax in fig.axes:
        removed_points, kept_points = ax.collections[0], ax.collections[1]
        assert len(removed_points.get_offsets()) == removed
        assert len(kept_points.get_offsets()) == preview.n_cells_kept


def test_threshold_lines_follow_the_settings(qc_run):
    obs = qc_run.result("qc_metrics").obs
    reference = plots.qc_scatters(obs, None, REFERENCE_PRESET.filter)
    mito_ax, genes_ax = reference.axes
    assert [line.get_ydata()[0] for line in mito_ax.get_lines()] == [REFERENCE_PRESET.filter.max_pct_counts_mt]
    assert sorted(line.get_ydata()[0] for line in genes_ax.get_lines()) == sorted(
        [REFERENCE_PRESET.filter.min_genes, REFERENCE_PRESET.filter.max_n_genes_by_counts])

    no_limits = FilterParams(max_n_genes_by_counts=None, max_pct_counts_mt=None)
    unlimited = plots.qc_scatters(obs, None, no_limits)
    assert unlimited.axes[0].get_lines() == []  # no mito line
    assert len(unlimited.axes[1].get_lines()) == 1  # only the lower genes line


# --- Later stages: settings boxes -> backend settings ------------------------------------

def test_default_stage_boxes_reproduce_the_reference_settings():
    s = display.default_stage_inputs()
    P = REFERENCE_PRESET
    pre = display.preprocess_settings(P, s["target_sum"], s["hvg_min_mean"], s["hvg_max_mean"],
                                      s["hvg_min_disp"], s["regress"], s["scale_max"])
    assert (pre["normalize"], pre["hvg"], pre["regress_out"], pre["scale"]) == (
        P.normalize, P.hvg, P.regress_out, P.scale)
    red = display.reduce_settings(P, s["n_comps"], s["n_neighbors"], s["n_pcs"])
    assert (red["pca"], red["neighbors"]) == (P.pca, P.neighbors)
    assert display.cluster_settings(P, s["resolution"])["leiden"] == P.leiden
    assert display.marker_settings(P, s["marker_method"])["marker_params"] == P.markers
    ann = display.annotation_settings(P, s["ambiguity_margin"], s["max_de_rank"], s["conflict_z"])
    assert ann["annotation_params"] == P.annotation


def test_changed_boxes_change_only_their_settings():
    P = REFERENCE_PRESET
    leiden = display.cluster_settings(P, 0.5)["leiden"]
    assert leiden.resolution == 0.5 and leiden.random_state == P.leiden.random_state
    off = display.preprocess_settings(P, 1e4, 0.0125, 3, 0.5, False, 10)["regress_out"]
    assert off.enabled is False and off.keys == P.regress_out.keys


@pytest.mark.parametrize("call, message", [
    (lambda P: display.reduce_settings(P, 40, 10, 50), "can be at most 40"),
    (lambda P: display.reduce_settings(P, 0, 10, 0), "'Number of PCs' must be greater than 0"),
    (lambda P: display.cluster_settings(P, None), "Enter a value for 'Resolution'"),
    (lambda P: display.preprocess_settings(P, 1e4, 3, 3, 0.5, True, 10), "must be smaller than"),
    (lambda P: display.annotation_settings(P, 0.5, 2.5, 1.0), "must be a whole number"),
])
def test_unusable_stage_boxes_give_a_message(call, message):
    with pytest.raises(display.InputError, match=message):
        call(REFERENCE_PRESET)


# --- Later stages: tables and plots on a full run ---------------------------------------

@pytest.fixture(scope="module")
def full_run(raw_data_path):
    from pbmc_analysis.pipeline import run_analysis
    return run_analysis(raw_data_path)


def test_cluster_sizes_table(full_run):
    sizes = full_run.stage_result("cluster").summary["cluster_sizes"]
    table = display.cluster_sizes_table(sizes)
    assert table["Cells"].tolist() == list(sizes.values())
    assert abs(table["% of cells"].sum() - 100) < 0.5


def test_marker_genes_table(full_run):
    first = full_run.stage_result("markers").summary["groups"][0]
    table = display.marker_genes_table(full_run.result("markers"), first, n_top=5)
    expected_top = full_run.stage_result("markers").summary["top10_markers"][first][:5]
    assert table["Gene"].tolist() == expected_top and table["Rank"].tolist() == [1, 2, 3, 4, 5]
    assert table["% in cluster"].between(0, 100).all()


def test_annotation_table(full_run):
    details = full_run.stage_result("annotate").details[0]
    table = display.annotation_table(details.clusters)
    assert list(table.columns) == list(display.ANNOTATION_COLUMNS.values())
    assert table["Predicted cell type"].tolist() == details.clusters["predicted_label"].tolist()


PLOTS = {
    "hvg": ("preprocess", lambda a: plots.hvg_plot(a)),
    "elbow": ("reduce_dimensions", lambda a: plots.pca_variance_plot(a)),
    "umap": ("reduce_dimensions", lambda a: plots.umap_plot(a)),
    "clusters": ("cluster", lambda a: plots.umap_plot(a, "leiden", "Leiden clusters", True)),
    "dotplot": ("markers", lambda a: plots.marker_dotplot(a)),
    "cell types": ("annotate", lambda a: plots.umap_plot(a, "cell_type_predicted")),
}


@pytest.mark.parametrize("name", PLOTS)
def test_plots_draw_without_changing_the_data(full_run, name):
    stage, draw = PLOTS[name]
    adata = full_run.result(stage)
    uns_before, obs_before = set(adata.uns), list(adata.obs.columns)
    fig = draw(adata)
    assert fig.axes
    assert set(adata.uns) == uns_before and list(adata.obs.columns) == obs_before


def test_hvg_plot_shows_all_genes_not_only_hvgs(full_run):
    adata = full_run.result("preprocess")
    fig = plots.hvg_plot(adata)
    points = sum(len(c.get_offsets()) for c in fig.axes[0].collections)
    assert points == adata.raw.n_vars  # every gene, HVGs and non-HVGs


def test_cluster_plot_labels_every_cluster(full_run):
    adata = full_run.result("cluster")
    fig = plots.umap_plot(adata, "leiden", "Leiden clusters", labels_on_plot=True)
    labels = {t.get_text() for t in fig.axes[0].texts}
    assert labels == set(adata.obs["leiden"].astype(str))  # whatever clusters the data gave


def test_cell_type_plot_legend_uses_predicted_labels(full_run):
    adata = full_run.result("annotate")
    fig = plots.umap_plot(adata, "cell_type_predicted")
    legend = {t.get_text() for t in fig.axes[0].get_legend().get_texts()}
    assert legend == set(adata.obs["cell_type_predicted"].astype(str))


def test_dotplot_has_no_dendrogram(full_run):
    adata = full_run.result("markers")
    plots.marker_dotplot(adata)
    assert not any(key.startswith("dendrogram") for key in adata.uns)


# --- The app itself ---------------------------------------------------------------------

def test_app_loads():
    spec = importlib.util.spec_from_file_location("pbmc_shiny_app", APP_DIR / "app.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from shiny import App
    assert isinstance(module.app, App)
    assert module.DATA_PATH == PROJECT_ROOT / "data" / "pbmc3k.h5ad"
