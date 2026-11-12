"""Tests for the downloads: scrna_analysis/export.py (the analyzed .h5ad) and app/report.py
(the run report and its Methods paragraph), plus the page's export information.

The file must hold exactly what the run computed (compared value by value with the checkpoints),
read back intact, keep the analyst's labels apart from the predictions, and leave out steps that
aren't up to date.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace

import anndata as ad
import numpy as np
import pytest

from scrna_analysis.export import ExportError, analyzed_anndata
from scrna_analysis.pipeline import AnalysisRun, run_analysis
from tests.conftest import PROJECT_ROOT
from tests.reference_preset import REFERENCE_PRESET

sys.path.insert(0, str(PROJECT_ROOT / "app"))
import payloads  # noqa: E402
import report  # noqa: E402

pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")


@pytest.fixture(scope="module")
def reviewed_run(raw_data_path):
    """A whole analysis with one cluster renamed (the label contains characters HTML must escape)."""
    run = run_analysis(raw_data_path, REFERENCE_PRESET)
    cluster = run.stage_result("annotate").details["summary"].clusters["cluster"].iloc[1]
    run.review_cluster(cluster, "rename", label="Cells <b>& co", note="note with <i>")
    return run


def test_the_file_holds_what_the_run_computed_and_reads_back(reviewed_run, tmp_path):
    run = reviewed_run
    out = analyzed_anndata(run)
    path = tmp_path / "analyzed.h5ad"
    out.write_h5ad(path, compression="gzip")
    back = ad.read_h5ad(path)

    annotated, counts = run.result("annotate"), run.result("filter")
    assert back.shape == annotated.raw.shape  # all genes kept after QC
    assert list(back.var_names) == list(counts.var_names) and list(back.obs_names) == list(counts.obs_names)
    assert np.array_equal(back.X.toarray(), annotated.raw.X.toarray())  # log-normalized
    assert np.array_equal(back.layers["counts"].toarray(), counts.X.toarray())  # raw counts
    assert back.var["highly_variable"].sum() == annotated.n_vars
    assert np.allclose(back.obsm["X_umap"], annotated.obsm["X_umap"])
    assert (back.obsp["connectivities"] != annotated.obsp["connectivities"]).nnz == 0
    assert {"rank_genes_groups", "pca", "neighbors", "cell_type_prediction", "orditus"} <= set(back.uns)

    # Predictions untouched; the analyst's labels in their own column.
    assert list(back.obs["cell_type_predicted"]) == list(annotated.obs["cell_type_predicted"])
    assert list(back.obs["cell_type_reviewed"]) == list(run.reviewed_labels())
    assert "Cells <b>& co" in set(back.obs["cell_type_reviewed"])

    record = json.loads(back.uns["orditus"])
    assert record["export"]["steps_included"] == ["preprocess", "reduce_dimensions", "cluster", "markers", "annotate"]
    assert record["export"]["obs_columns"]["reviewed_labels"] == "cell_type_reviewed"
    assert [r["label"] for r in record["review"]] == ["Cells <b>& co"]
    assert record["params"] == run.params.to_dict()


def test_only_steps_that_are_up_to_date_are_exported(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load(); run.compute_qc_metrics(); run.apply_filter()
    with pytest.raises(ExportError, match="Run Preprocessing first"):
        analyzed_anndata(run)
    run.preprocess()
    out = analyzed_anndata(run)
    assert json.loads(out.uns["orditus"])["export"]["steps_included"] == ["preprocess"]
    assert "leiden" not in out.obs and "X_umap" not in out.obsm and "cell_type_reviewed" not in out.obs


def test_methods_paragraph_follows_the_settings(reviewed_run):
    run = reviewed_run
    p, s = run.params, run.stage_result("filter").summary
    text = report.methods_text(run)
    for expected in (f"Scanpy {run.record()['versions']['scanpy']}", f"fewer than {p.filter.min_genes:,} genes",
                     f"leaving {s['n_cells_kept']:,} cells", f"resolution {p.leiden.resolution:g}",
                     f"top {p.annotation.test.top_n_genes} marker genes", "1 cluster confirmed or renamed"):
        assert expected in text
    assert "regressed out" in text  # the notebook's settings regress; the defaults don't
    # Another resolution changes the text (it is written from the settings, not fixed).
    other = AnalysisRun(None, replace(p, leiden=replace(p.leiden, resolution=0.4)))
    other._results, other.history = run._results, run.history  # same results, other settings
    assert "resolution 0.4" in report.methods_text(other)


def test_run_report_page(reviewed_run):
    page = report.run_report_html(reviewed_run, "My dataset")
    assert page.startswith("<!doctype html>") and page.endswith("</html>")
    for expected in ("Single-cell RNA-seq run report", "My dataset", "Methods", reviewed_run.file_name,
                     reviewed_run.record()["versions"]["scanpy"], "Run log", "Your label"):
        assert expected in page
    assert "Cells &lt;b&gt;&amp; co" in page and "<b>& co" not in page  # typed text is shown, not run
    assert "note with &lt;i&gt;" in page
    assert "<th>Species and tissue</th><td>not confirmed</td>" in page  # this run had no metadata
    reviewed_run.metadata = {"species": "human", "tissues": ["Blood"], "about": "x", "confirmed_by": "user",
                             "confirmed_at": None}
    page = report.run_report_html(reviewed_run, "My dataset")
    assert "<th>Species and tissue</th><td>Human · Blood</td>" in page and "<td>the analyst</td>" in page
    reviewed_run.metadata = {}


def test_export_information_for_the_page(raw_data_path, reviewed_run):
    empty = payloads.export_payload(AnalysisRun(None, REFERENCE_PRESET))
    assert empty["available"] is False and "Run Preprocessing first" in empty["reason"]
    full = payloads.export_payload(reviewed_run)
    assert full["available"] and full["reviewed"]
    assert full["includes"] == ["Preprocessing", "Dimensionality Reduction", "Clustering", "Marker Genes",
                                "Cell Type Annotation"]
    assert full["h5ad_name"] == "pbmc3k_analyzed.h5ad" and full["report_name"] == "pbmc3k_run-report.html"
    assert payloads.export_file_name("My Brain Cells.h5", "analyzed", "h5ad") == "my-brain-cells_analyzed.h5ad"
    assert payloads.export_file_name(None, "run report", "html") == "run-report.html"
