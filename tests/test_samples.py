"""Tests for loading several samples (scrna_analysis/load_data.py: `load_samples`,
`sample_column_choices`, `samples_from_column`; AnalysisRun with `sample_files` / `define_samples`).

Two ways in: one file per sample, or one file with a sample column. Checked on small made-up data
with known answers, and on the Kang et al. 2018 dataset (16 samples: 8 patients x control /
interferon-beta), where both ways must give the same cells and genes. A single dataset must load
exactly as before (no sample column).
"""

from __future__ import annotations

import json
import re

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scrna_analysis import load_data
from scrna_analysis.load_data import InputValidationError, SampleFile, load_samples, sample_column_choices, samples_from_column
from scrna_analysis.params import DEFAULT_PARAMS, MAX_SAMPLES
from scrna_analysis.pipeline import AnalysisRun, PipelineError


def tiny(genes, barcodes, seed=0, **obs) -> ad.AnnData:
    """A small raw-count dataset with the given gene names and cell barcodes."""
    rng = np.random.default_rng(seed)
    a = ad.AnnData(sparse.csr_matrix(rng.poisson(1.0, (len(barcodes), len(genes))).astype(np.float32)))
    a.var_names, a.obs_names = list(genes), list(barcodes)
    for key, values in obs.items():
        a.obs[key] = values
    return a


def write(tmp_path, name, adata):
    path = tmp_path / name
    adata.write_h5ad(path)
    return str(path)


# --- One file per sample -----------------------------------------------------------------------

def test_files_are_combined_with_sample_names_and_shared_genes(tmp_path):
    a = write(tmp_path, "a.h5ad", tiny(["G1", "G2", "G3"], ["AAA", "CCC"]))
    b = write(tmp_path, "b.h5ad", tiny(["G2", "G3", "G4"], ["AAA", "GGG", "TTT"], seed=1))  # same barcode "AAA"
    adata, s = load_samples([SampleFile(a, " donor 1 ", condition="ctrl"), SampleFile(b, "donor 2")])
    assert list(adata.obs_names) == ["donor 1:AAA", "donor 1:CCC", "donor 2:AAA", "donor 2:GGG", "donor 2:TTT"]
    assert list(adata.var_names) == ["G2", "G3"] and s.n_genes_not_in_all_samples == 2  # G1, G4 left out
    assert list(adata.obs["sample"]) == ["donor 1"] * 2 + ["donor 2"] * 3
    assert list(adata.obs["condition"]) == ["ctrl", "ctrl", "", "", ""]
    assert [(x["name"], x["condition"], x["n_cells"]) for x in s.samples] == [("donor 1", "ctrl", 2), ("donor 2", None, 3)]
    assert "2 genes measured in only some samples were left out" in s.notes[0]
    # The counts are the files' own, unchanged.
    original = ad.read_h5ad(b)
    assert np.array_equal(adata[adata.obs["sample"] == "donor 2", ["G2", "G3"]].X.toarray(), original[:, ["G2", "G3"]].X.toarray())


def test_no_condition_column_when_none_is_given(tmp_path):
    a = write(tmp_path, "a.h5ad", tiny(["G1", "G2"], ["A1"]))
    b = write(tmp_path, "b.h5ad", tiny(["G1", "G2"], ["B1"]))
    adata, _ = load_samples([SampleFile(a, "x"), SampleFile(b, "y")])
    assert "condition" not in adata.obs


@pytest.mark.parametrize("names, message", [
    (["only one"], "at least two"),
    (["same", "same"], "repeated: same"),
    (["ok", "   "], "needs a name"),
    (["ok", "x" * 41], "too long"),
])
def test_unusable_sample_lists_are_explained(tmp_path, names, message):
    files = [SampleFile(write(tmp_path, f"{i}.h5ad", tiny(["G1"], [f"C{i}"])), n) for i, n in enumerate(names)]
    with pytest.raises(InputValidationError, match=message):
        load_samples(files)


def test_too_many_samples(tmp_path):
    path = write(tmp_path, "a.h5ad", tiny(["G1"], ["C"]))
    with pytest.raises(InputValidationError, match=f"more than the app handles for now \\({MAX_SAMPLES}\\)"):
        load_samples([SampleFile(path, f"s{i}") for i in range(MAX_SAMPLES + 1)])


def test_samples_must_name_genes_the_same_way(tmp_path):
    a = write(tmp_path, "a.h5ad", tiny(["CD3D", "MS4A1"], ["A"]))
    b = write(tmp_path, "b.h5ad", tiny(["ENSG00000167286", "ENSG00000156738"], ["B"]))
    with pytest.raises(InputValidationError, match="Ensembl ID and others by gene symbol"):
        load_samples([SampleFile(a, "a"), SampleFile(b, "b")])


def test_each_file_is_checked_like_a_single_file(tmp_path):
    good = write(tmp_path, "good.h5ad", tiny(["G1"], ["A"]))
    bad_data = tiny(["G1"], ["B"])
    bad_data.X = bad_data.X.toarray() * 0.5 + 0.25  # not whole numbers
    bad = write(tmp_path, "bad.h5ad", bad_data)
    with pytest.raises(InputValidationError, match="non-integer values"):
        load_samples([SampleFile(good, "good"), SampleFile(bad, "bad")])


# --- One file with a sample column ---------------------------------------------------------------

def test_sample_column_choices_suggest_likely_columns():
    a = tiny(["G1"], [f"c{i}" for i in range(6)], donor=["d1", "d1", "d2", "d2", "d3", "d3"],
             treatment=["ctrl", "stim"] * 3, n_counts=[1, 2, 3, 4, 5, 6], unique_id=[f"u{i}" for i in range(6)],
             one=["x"] * 6)
    choices = sample_column_choices(a)
    assert [c["column"] for c in choices] == ["donor", "unique_id", "treatment"]  # numbers and 1-value columns left out
    assert choices[0]["hint"] == "sample" and choices[-1]["hint"] == "condition" and choices[0]["n_values"] == 3


def test_samples_from_a_column_and_conditions():
    a = tiny(["G1", "G2"], ["c1", "c2", "c3", "c4"], donor=["d1", "d1", "d2", "d2"], group=["ctrl", "ctrl", "stim", "stim"])
    out, samples, notes = samples_from_column(a, "donor", "group")
    assert list(out.obs["sample"]) == ["d1", "d1", "d2", "d2"] and list(out.obs["condition"]) == ["ctrl", "ctrl", "stim", "stim"]
    assert [(s["name"], s["condition"], s["n_cells"]) for s in samples] == [("d1", "ctrl", 2), ("d2", "stim", 2)]
    assert notes == () and "sample" not in a.obs  # the input is unchanged


def test_a_sample_under_several_conditions_becomes_one_sample_per_condition():
    a = tiny(["G1"], ["c1", "c2", "c3", "c4"], donor=["d1", "d1", "d2", "d2"], group=["ctrl", "stim", "ctrl", "stim"])
    out, samples, notes = samples_from_column(a, "donor", "group")
    assert [s["name"] for s in samples] == ["d1 · ctrl", "d1 · stim", "d2 · ctrl", "d2 · stim"]
    assert "treated as its own sample" in notes[0]


@pytest.mark.parametrize("column, condition, message", [
    ("nope", None, "no column named 'nope'"),
    ("donor", "donor", "different columns"),
    ("same", None, "at least two"),
])
def test_unusable_sample_columns_are_explained(column, condition, message):
    a = tiny(["G1"], ["c1", "c2"], donor=["d1", "d2"], same=["x", "x"])
    with pytest.raises(InputValidationError, match=message):
        samples_from_column(a, column, condition)


def test_missing_sample_values_are_refused():
    a = tiny(["G1"], ["c1", "c2", "c3"], donor=pd.Categorical(["d1", None, "d2"]))
    with pytest.raises(InputValidationError, match="have no value in 'donor'"):
        samples_from_column(a, "donor")


# --- In the analysis run ------------------------------------------------------------------------

def test_a_single_dataset_loads_exactly_as_before(raw_data_path):
    """No sample column, no samples in the summary or the record: single-sample analyses don't change."""
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    run.load()
    assert "sample" not in run.result("load").obs and "condition" not in run.result("load").obs
    assert run.samples == [] and "samples" not in run.stage_result("load").summary
    assert not {"samples", "sample_column"} & set(run.record()["input"])


def test_kang_from_its_columns_equals_kang_from_one_file_per_sample(kang_path, tmp_path):
    """The Kang data: 8 patients, each under control and interferon-beta, so 16 samples. Choosing
    its columns, or splitting it into 16 files, gives the same cells, genes and counts."""
    run = AnalysisRun(kang_path, DEFAULT_PARAMS)
    run.load()
    run.define_samples("replicate", "label")
    by_column = run.result("load")
    assert len(run.samples) == 16 and by_column.obs["sample"].nunique() == 16
    assert set(by_column.obs["condition"]) == {"ctrl", "stim"}
    assert run.record()["input"]["sample_column"] == {"sample": "replicate", "condition": "label"}

    original = ad.read_h5ad(kang_path)
    files = []
    for (patient, label), cells in original.obs.groupby(["replicate", "label"], observed=True).groups.items():
        path = tmp_path / f"{patient}_{label}.h5ad"
        original[cells].copy().write_h5ad(path)
        files.append(SampleFile(str(path), f"{patient} · {label}", condition=label))
    run2 = AnalysisRun(None, DEFAULT_PARAMS, sample_files=files)
    run2.load()
    by_files = run2.result("load")
    assert by_files.shape == by_column.shape and list(by_files.var_names) == list(by_column.var_names)
    barcode = by_files.obs_names.str.split(":").str[-1]
    reordered = by_column[barcode]
    assert np.array_equal(by_files.X.toarray(), reordered.X.toarray())
    assert list(by_files.obs["sample"].astype(str)) == list(reordered.obs["sample"].astype(str))
    record = run2.record()["input"]
    assert len(record["samples"]) == 16 and all(len(s["sha256"]) == 64 and "path" not in s for s in record["samples"])
    json.dumps(run2.record())  # storable


def test_choosing_samples_again_or_undoing_reloads_and_clears_later_steps(kang_path):
    run = AnalysisRun(kang_path, DEFAULT_PARAMS)
    run.load()
    run.compute_qc_metrics()
    run.define_samples("replicate")  # samples = patients only
    assert len(run.samples) == 8 and run.status()["qc_metrics"] == "needs re-running"
    with pytest.raises(load_data.InputValidationError, match="no column named"):
        run.define_samples("nope")
    assert len(run.samples) == 8  # a failed choice changes nothing
    run.define_samples(None)
    assert run.samples == [] and "sample" not in run.result("load").obs


def test_sample_columns_need_a_single_file(tmp_path):
    a = write(tmp_path, "a.h5ad", tiny(["G1"], ["A"]))
    b = write(tmp_path, "b.h5ad", tiny(["G1"], ["B"]))
    run = AnalysisRun(None, DEFAULT_PARAMS, sample_files=[SampleFile(a, "a"), SampleFile(b, "b")])
    run.load()
    with pytest.raises(PipelineError, match="only for a single loaded file"):
        run.define_samples("sample")


# --- Quality Control per sample ------------------------------------------------------------------

def kang_with_samples(kang_path):
    run = AnalysisRun(kang_path, DEFAULT_PARAMS)
    run.load()
    run.define_samples("replicate", "label")
    return run


def test_each_sample_gets_limits_from_its_own_cells(kang_path):
    """Known answer: each sample's suggested limits equal the single-dataset suggestion computed on
    that sample alone; the other fields keep the limits for all cells together."""
    from scrna_analysis import qc
    from scrna_analysis.params import SuggestionParams
    run = kang_with_samples(kang_path)
    run.compute_qc_metrics()
    adata = run.result("qc_metrics")
    suggestion = qc.suggest_filter(adata, SuggestionParams())
    assert [s.sample for s in suggestion.filter.sample_limits] == [s["name"] for s in run.samples]
    for limits in suggestion.filter.sample_limits:
        alone = qc.suggest_filter(adata[adata.obs["sample"] == limits.sample].copy(), SuggestionParams(per_sample=False))
        assert (limits.min_genes, limits.max_n_genes_by_counts, limits.max_pct_counts_mt) == (
            alone.filter.min_genes, alone.filter.max_n_genes_by_counts, alone.filter.max_pct_counts_mt)
    pooled = qc.suggest_filter(adata, SuggestionParams(per_sample=False))
    assert pooled.filter.sample_limits == () and pooled.filter.min_genes == suggestion.filter.min_genes
    assert len({s.max_n_genes_by_counts for s in suggestion.filter.sample_limits}) > 1  # the samples really differ
    assert run.params.filter == suggestion.filter  # what the app starts with


def test_each_cell_is_judged_by_its_own_samples_limits(kang_path):
    """Known answer: the cells kept per sample are exactly those a single-dataset filter keeps on
    that sample alone with the same limits; the preview and the filtering agree."""
    from scrna_analysis import qc
    from scrna_analysis.params import FilterParams
    run = kang_with_samples(kang_path)
    run.compute_qc_metrics()
    adata, params = run.result("qc_metrics"), run.params.filter
    preview = qc.preview_filter(adata, params)
    for limits in params.sample_limits:
        alone = adata[adata.obs["sample"] == limits.sample].copy()
        single = FilterParams(min_genes=limits.min_genes, min_cells=0, max_n_genes_by_counts=limits.max_n_genes_by_counts,
                              max_pct_counts_mt=limits.max_pct_counts_mt)
        expected = set(qc.preview_filter(alone, single).kept_cells)
        assert {c for c in preview.kept_cells if c in set(alone.obs_names)} == expected
        assert preview.by_sample[limits.sample]["n_cells_kept"] == len(expected)
    filtered = run.apply_filter().adata  # apply_filter itself checks it equals the preview
    assert filtered.n_obs == preview.n_cells_kept and "n_genes" in filtered.obs
    assert sum(v["n_cells_kept"] for v in run.stage_result("filter").summary["by_sample"].values()) == filtered.n_obs


def test_one_set_of_limits_for_all_samples(kang_path):
    from scrna_analysis.params import SuggestionParams
    run = kang_with_samples(kang_path)
    run.compute_qc_metrics(suggest_params=SuggestionParams(per_sample=False))
    assert run.params.filter.sample_limits == ()
    run.apply_filter()
    assert "by_sample" in run.stage_result("filter").summary  # kept cells are still counted per sample


def test_limits_per_sample_need_every_sample_and_a_sample_column(kang_path, raw_data_path):
    from dataclasses import replace
    from scrna_analysis import qc
    run = kang_with_samples(kang_path)
    run.compute_qc_metrics()
    params = run.params.filter
    with pytest.raises(qc.QCError, match="No limits were given for sample"):
        qc.preview_filter(run.result("qc_metrics"), replace(params, sample_limits=params.sample_limits[1:]))
    single = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    single.load()
    single.compute_qc_metrics()
    assert single.params.filter.sample_limits == ()
    assert "sample_limits" not in single.record()["params"]["filter"]  # a single dataset's record is unchanged
    with pytest.raises(qc.QCError, match="need a sample column"):
        qc.preview_filter(single.result("qc_metrics"), params)


def test_quality_control_page_data_per_sample(kang_path, raw_data_path):
    """The QC page gets the samples and each cell's sample; its boxes can send limits per sample;
    the live preview counts cells per sample, next to what the suggestion keeps."""
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import payloads

    run = kang_with_samples(kang_path)
    run.compute_qc_metrics()
    page = payloads.qc_payload(run)
    names = [s["name"] for s in run.samples]
    assert page["samples"] == names and len(page["cells"]["sample"]) == page["n_cells"]
    assert {page["samples"][i] for i in page["cells"]["sample"]} == set(names)
    assert len(page["filter"]["suggested"]["sample_limits"]) == 16

    first = run.params.filter.sample_limits[0]
    boxes = {"min_genes": "200", "min_cells": "3", "max_genes": "2500", "max_mt": "", "max_mt_off": True, "per_sample": True,
             "samples": [{"sample": s.sample, "min_genes": str(s.min_genes), "max_genes": str(s.max_n_genes_by_counts),
                          "max_mt": "", "max_mt_off": True} for s in run.params.filter.sample_limits]}
    boxes["samples"][0]["min_genes"] = "600"  # stricter for one sample only
    params = payloads.filter_params_from_browser(boxes, run.params.filter)
    assert params.sample_limits[0].min_genes == 600 and params.sample_limits[1] == run.params.filter.sample_limits[1]
    preview = payloads.preview_payload(run, boxes)
    rows = {r["sample"]: r for r in preview["by_sample"]}
    assert rows[first.sample]["n_cells_kept"] < rows[first.sample]["suggested_kept"]  # only that sample loses cells
    assert all(r["n_cells_kept"] == r["suggested_kept"] for name, r in rows.items() if name != first.sample)
    assert next(r for r in preview["filters"] if r["key"] == "min_genes")["setting"][0] == 600

    bad = {**boxes, "samples": [{**boxes["samples"][0], "min_genes": "x"}] + boxes["samples"][1:]}
    with pytest.raises(payloads.InputError, match=re.escape(f"'Min genes per cell ({first.sample})' must be a number")):
        payloads.filter_params_from_browser(bad, run.params.filter)
    one = payloads.suggestion_payload(run, {"n_mads_genes": "5", "n_mads_mito": "3", "per_sample": False})
    assert one["suggested"]["sample_limits"] == () and one["suggest_settings"]["per_sample"] is False

    single = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    single.load()
    single.compute_qc_metrics()
    assert payloads.qc_payload(single)["samples"] == [] and payloads.qc_payload(single)["cells"]["sample"] is None
    assert payloads.preview_payload(single, None)["by_sample"] == []


def test_methods_paragraph_describes_the_samples(kang_path):
    """The run report's Methods paragraph says how many samples were combined, with which
    conditions, and that QC limits were set per sample."""
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import report

    run = kang_with_samples(kang_path)
    run.compute_qc_metrics()
    run.apply_filter()
    text = report.methods_text(run)
    assert "16 samples (ctrl and stim conditions) were combined, keeping genes measured in every sample." in text
    assert "Quality control limits were set for each sample from its own cells" in text
    assert f"leaving {run.stage_result('filter').summary['n_cells_kept']:,} cells" in text
    assert "Cells with fewer than" not in text  # the pooled limits weren't used
