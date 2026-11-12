"""Tests for `scrna_analysis.pipeline`.

- Reference: the whole analysis run through the pipeline reproduces the validated
  notebook through Marker Genes, and records everything. (Cell types now come from
  PanglaoDB, so they aren't compared with the notebook's.)
- Behavior: fixed stage order, redoing a stage clears later stages, failures change
  nothing, results are repeatable, and the pipeline equals running the modules by hand
  (also on a different input with different settings).
"""

from __future__ import annotations

import json
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import adjusted_rand_score

from scrna_analysis import annotation, clustering, dim_reduction, markers, preprocessing, qc
from scrna_analysis.load_data import load_h5ad
from scrna_analysis.params import FilterParams, LeidenParams, NeighborsParams, PCAParams, QCParams
from tests.reference_preset import REFERENCE_PRESET
from scrna_analysis.pipeline import (
    RECORD_KEY,
    STAGES,
    AnalysisRun,
    PipelineError,
    file_sha256,
    run_analysis,
)
from scrna_analysis.qc import QCError

# regress_out triggers harmless Numba performance hints (also seen in the notebook).
pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")

P = REFERENCE_PRESET


def labels(adata: ad.AnnData, key: str) -> np.ndarray:
    return adata.obs[key].astype(str).to_numpy()


def without_times(history: list[dict]) -> list[dict]:
    return [{k: v for k, v in entry.items() if k != "ran_at"} for entry in history]


# --- Reference: the whole analysis through the pipeline ---------------------------------

@pytest.fixture(scope="module")
def reference_run(raw_data_path) -> AnalysisRun:
    return run_analysis(raw_data_path, REFERENCE_PRESET)


def test_all_stages_done_in_order(reference_run):
    assert set(reference_run.status().values()) == {"done"}
    assert [entry["stage"] for entry in reference_run.history] == list(STAGES)


def test_end_result_matches_validated_notebook(reference_run, reference_results):
    ours, ref = reference_run.result(), reference_results
    assert ours.obs_names.equals(ref.obs_names)
    assert ours.var_names.equals(ref.var_names)  # the same HVGs
    assert ours.raw.var_names.equals(ref.raw.var_names)  # the same genes kept
    assert float(np.max(np.abs(np.asarray(ours.X) - np.asarray(ref.X)))) == 0
    assert float(np.max(np.abs(ours.obsm["X_pca"] - ref.obsm["X_pca"]))) == 0
    assert float(np.max(np.abs(ours.obsm["X_umap"] - ref.obsm["X_umap"]))) == 0
    np.testing.assert_array_equal(labels(ours, "leiden"), labels(ref, "leiden"))
    assert adjusted_rand_score(labels(ref, "leiden"), labels(ours, "leiden")) == 1.0
    pd.testing.assert_frame_equal(pd.DataFrame(ours.uns["rank_genes_groups"]["names"]),
                                  pd.DataFrame(ref.uns["rank_genes_groups"]["names"]))


def test_stage_summaries_match_reference(reference_run, reference_summary):
    summary = {stage: reference_run.stage_result(stage).summary for stage in STAGES}
    assert [summary["load"]["n_cells"], summary["load"]["n_genes"]] == reference_summary["raw_shape"]
    assert summary["qc_metrics"]["n_mito_genes"] == reference_summary["n_mito_genes"]
    assert [summary["filter"]["n_cells_kept"], summary["filter"]["n_genes_kept"]] == reference_summary["filtered_shape"]
    assert summary["preprocess"]["n_hvg"] == reference_summary["n_hvg"]
    assert summary["preprocess"]["regressed"] is True
    assert summary["cluster"]["cluster_sizes"] == reference_summary["cluster_sizes"]


def test_annotation_records_the_database_and_choices(reference_run, panglaodb_path):
    from scrna_analysis.pipeline import file_sha256
    summary = reference_run.stage_result("annotate").summary
    # No Census tissue in the preset, so PanglaoDB is the only (first) reference.
    assert summary["database"] == {"name": "PanglaoDB", "version": "27 Mar 2020", "file": panglaodb_path.name,
                                   "sha256": file_sha256(panglaodb_path), "tissues": ["Immune system", "Blood"],
                                   "role": "first"}
    assert summary["references"] == [summary["database"]] and set(summary["by_reference"]) == {"PanglaoDB"}
    assert {c["label_source"] for c in summary["clusters"] if c["predicted_label"] != "Unresolved"} == {"PanglaoDB"}
    settings = reference_run.record()["params"]["annotation"]
    assert settings["markers"]["species"] == "human" and settings["markers"]["tissues"] == ["Immune system", "Blood"]
    assert sum(summary["label_counts"].values()) == reference_run.result().n_obs


# --- The record ("methods section") ------------------------------------------------------

def test_record_contents(reference_run, raw_data_path, reference_summary):
    record = reference_run.record()
    assert record["params"] == P.to_dict()
    assert set(record["random_states"].values()) == {reference_summary["seed"]}
    assert record["input"]["sha256"] == file_sha256(raw_data_path)
    assert record["stages"] == {stage: "done" for stage in STAGES}
    assert len(record["history"]) == len(STAGES)


def test_environment_matches_reference(reference_run, reference_summary):
    """Same package versions as when the reference was made (a mismatch means: re-validate)."""
    ours = reference_run.record()["versions"]
    for package in ("scanpy", "anndata", "numpy", "scipy", "leidenalg", "igraph", "umap-learn"):
        assert ours[package] == reference_summary["versions"][package], package


def test_record_travels_with_the_data_and_survives_saving(reference_run, tmp_path):
    adata = reference_run.result()
    assert json.loads(adata.uns[RECORD_KEY]) == json.loads(json.dumps(reference_run.record()))
    path = tmp_path / "analysis.h5ad"
    adata.write_h5ad(path)
    reloaded = ad.read_h5ad(path)
    assert json.loads(reloaded.uns[RECORD_KEY])["params"] == P.to_dict()
    np.testing.assert_array_equal(labels(reloaded, "cell_type_predicted"),
                                  labels(adata, "cell_type_predicted"))


def test_each_history_entry_records_its_settings(reference_run):
    entry = {e["stage"]: e for e in reference_run.history}
    assert entry["cluster"]["settings"] == {"leiden": P.to_dict()["leiden"]}
    assert set(entry["preprocess"]["settings"]) == {"normalize", "hvg", "regress_out", "scale"}
    assert entry["load"]["settings"] == {}


def test_full_gene_counts_available(reference_run, reference_summary):
    counts = reference_run.filtered_counts()
    assert list(counts.shape) == reference_summary["filtered_shape"]  # all genes, raw counts
    values = counts.X.data
    assert np.array_equal(values, np.round(values)) and "log1p" not in counts.uns


def test_preview_filter_changes_nothing(reference_run):
    history_before = len(reference_run.history)
    preview = reference_run.preview_filter(FilterParams(min_genes=500))
    assert preview.n_cells_kept < reference_run.stage_result("filter").summary["n_cells_kept"]
    assert len(reference_run.history) == history_before
    assert reference_run.params == P


# --- Stage order and failures ------------------------------------------------------------

def test_stages_cannot_be_skipped(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    with pytest.raises(PipelineError, match="Can't run 'cluster' yet: run 'load' first"):
        run.cluster()
    run.load()
    with pytest.raises(PipelineError, match="Can't run 'filter' yet: run 'qc_metrics' first"):
        run.apply_filter()
    with pytest.raises(PipelineError, match="run 'qc_metrics' first"):
        run.preview_filter()
    assert run.status()["load"] == "done" and run.status()["filter"] == "not run"


def test_asking_for_missing_results(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    with pytest.raises(PipelineError, match="Nothing has been run yet"):
        run.result()
    run.load()
    with pytest.raises(PipelineError, match="Stage 'cluster' has no result"):
        run.result("cluster")
    with pytest.raises(PipelineError, match="Unknown stage 'clustering'"):
        run.result("clustering")


def test_a_failed_stage_changes_nothing(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load()
    run.compute_qc_metrics()
    history_before = list(run.history)
    with pytest.raises(QCError, match="would remove all cells"):
        run.apply_filter(FilterParams(min_genes=100_000))
    assert run.status()["filter"] == "not run"
    assert run.params == P  # the failed settings were not adopted
    assert run.history == history_before


# --- Redoing a stage, and repeatability --------------------------------------------------

@pytest.fixture(scope="module")
def redone(raw_data_path):
    """A second full run; a snapshot is taken, then clustering is redone at resolution 0.5."""
    run = run_analysis(raw_data_path, REFERENCE_PRESET)
    snapshot = {
        "history": without_times(run.history),
        "leiden": labels(run.result(), "leiden"),
        "cell_type": labels(run.result(), "cell_type_predicted"),
        "umap": run.result().obsm["X_umap"].copy(),
    }
    checkpoints_before = {stage: run.result(stage) for stage in STAGES}
    run.cluster(LeidenParams(resolution=0.5))
    return run, snapshot, checkpoints_before


def test_two_runs_give_identical_results(reference_run, redone):
    _, snapshot, _ = redone
    final = reference_run.result()
    assert snapshot["history"] == without_times(reference_run.history)
    np.testing.assert_array_equal(snapshot["leiden"], labels(final, "leiden"))
    np.testing.assert_array_equal(snapshot["cell_type"], labels(final, "cell_type_predicted"))
    np.testing.assert_array_equal(snapshot["umap"], final.obsm["X_umap"])


def test_redoing_a_stage_clears_later_stages(redone):
    run, _, checkpoints_before = redone
    status = run.status()
    assert status["cluster"] == "done"
    assert status["markers"] == status["annotate"] == "needs re-running"
    assert run.history[-1]["stage"] == "cluster"
    assert run.history[-1]["cleared_later_stages"] == ["markers", "annotate"]
    # Earlier checkpoints are kept as they were (not recomputed).
    for stage in STAGES[:STAGES.index("cluster")]:
        assert run.result(stage) is checkpoints_before[stage]


def test_redo_uses_and_records_the_new_settings(redone, reference_summary):
    run, _, _ = redone
    assert run.params.leiden.resolution == 0.5
    assert run.record()["params"]["leiden"]["resolution"] == 0.5
    n_clusters = run.stage_result("cluster").summary["n_clusters"]
    assert n_clusters < reference_summary["n_clusters"]  # coarser resolution, fewer clusters


def test_later_stages_can_be_rerun_after_a_redo(redone):
    run, _, _ = redone
    run.find_markers()
    run.annotate()
    assert set(run.status().values()) == {"done"}
    groups = run.stage_result("markers").summary["groups"]
    assert len(groups) == run.stage_result("cluster").summary["n_clusters"]


# --- Same as running the modules by hand -------------------------------------------------

def run_modules_by_hand(path, params):
    a, _ = load_h5ad(path)
    a, _ = qc.calculate_qc_metrics(a, params.qc)
    a, _ = qc.apply_filter(a, params.filter)
    a, _ = preprocessing.normalize_and_log(a, params.normalize)
    a, _ = preprocessing.find_hvgs(a, params.hvg)
    a, _ = preprocessing.keep_raw_and_subset(a)
    a, _ = preprocessing.regress_and_scale(a, params.regress_out, params.scale)
    a, _ = dim_reduction.run_pca(a, params.pca)
    a, _ = dim_reduction.compute_neighbors(a, params.neighbors)
    a, _ = dim_reduction.run_umap(a, params.umap)
    a, _ = clustering.cluster_leiden(a, params.leiden)
    a, _ = markers.rank_marker_genes(a, params.markers)
    from scrna_analysis import marker_db
    sets = marker_db.select_marker_sets(marker_db.load_panglaodb(), params.annotation.markers)
    found = marker_db.match_to_genes(sets, a.raw.var_names, params.annotation.markers.min_markers)
    a, _ = annotation.predict_cell_types(a, found.sets, params.annotation.test)
    return a


def test_other_input_and_settings_match_modules_by_hand(raw_data_path, tmp_path):
    raw, _ = load_h5ad(raw_data_path)
    rng = np.random.default_rng(21)
    subset_path = tmp_path / "subset.h5ad"
    raw[np.sort(rng.choice(raw.n_obs, 800, replace=False))].copy().write_h5ad(subset_path)
    params = replace(
        P,
        filter=FilterParams(min_genes=300, max_n_genes_by_counts=2000, max_pct_counts_mt=5.0),
        pca=PCAParams(n_comps=20),
        neighbors=NeighborsParams(n_neighbors=15, n_pcs=20),
        leiden=LeidenParams(resolution=0.8),
    )

    ours = run_analysis(subset_path, params).result()
    by_hand = run_modules_by_hand(subset_path, params)

    assert ours.obs_names.equals(by_hand.obs_names) and ours.var_names.equals(by_hand.var_names)
    np.testing.assert_array_equal(ours.obsm["X_umap"], by_hand.obsm["X_umap"])
    np.testing.assert_array_equal(labels(ours, "leiden"), labels(by_hand, "leiden"))
    np.testing.assert_array_equal(labels(ours, "cell_type_predicted"),
                                  labels(by_hand, "cell_type_predicted"))
    assert json.loads(ours.uns[RECORD_KEY])["params"]["leiden"]["resolution"] == 0.8


# --- Mitochondrial genes for any species -----------------------------------------------------

def test_no_mito_genes_switches_off_the_mito_filter_and_records_it(mouse_10x_path):
    """With the human pattern, mouse data has no matching genes: the % mitochondrial filter is
    switched off (recorded), so QC isn't blocked."""
    run = AnalysisRun(mouse_10x_path, REFERENCE_PRESET)
    run.load()
    assert run.params.filter.max_pct_counts_mt is not None
    result = run.compute_qc_metrics(replace(REFERENCE_PRESET.qc, mito_prefix="MT-"))
    assert result.summary["mito_filter_switched_off"] is True and not result.summary["mito_qc_available"]
    assert run.params.filter.max_pct_counts_mt is None
    assert run.history[-1]["summary"]["mito_filter_switched_off"] is True  # in the record
    run.apply_filter()  # no longer blocked


def test_mouse_data_runs_to_the_end_with_the_mouse_pattern(mouse_10x_path):
    """A second dataset (mouse brain, 10x .h5) runs through every stage once the prefix is "mt-"
    and the species and tissue are mouse brain."""
    from scrna_analysis.params import AnnotationParams, MarkerDatabaseParams
    params = replace(REFERENCE_PRESET, qc=replace(REFERENCE_PRESET.qc, mito_prefix="mt-"),
                     annotation=AnnotationParams(markers=MarkerDatabaseParams(species="mouse", tissues=("Brain",))))
    run = AnalysisRun(mouse_10x_path, params)
    run.run_all()
    qc_summary = run.stage_result("qc_metrics").summary
    assert qc_summary["n_mito_genes"] > 0 and qc_summary["mito_filter_switched_off"] is False
    assert run.result("filter").n_obs < run.result("load").n_obs  # the mito filter did something
    assert run.result("cluster").obs["leiden"].nunique() > 1
    assert run.stage_result("annotate").summary["markers_matched_ignoring_case"] > 0  # CD3D -> Cd3d style


def test_annotation_needs_a_tissue_and_a_failure_changes_nothing(raw_data_path):
    from scrna_analysis.params import AnnotationParams, MarkerDatabaseParams
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster",
                  "find_markers"):
        getattr(run, stage)()
    with pytest.raises(annotation.AnnotationError, match="Choose at least one tissue"):
        run.annotate(AnnotationParams(markers=MarkerDatabaseParams(species="human")))
    assert run.status()["annotate"] == "not run" and run.params == REFERENCE_PRESET  # nothing changed


def test_a_missing_marker_database_is_explained(raw_data_path, tmp_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET, marker_db_path=tmp_path / "absent.tsv.gz")
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster",
                  "find_markers"):
        getattr(run, stage)()
    from scrna_analysis.marker_db import MarkerDatabaseError
    with pytest.raises(MarkerDatabaseError, match="download_marker_db.py"):
        run.annotate()


def test_mito_filter_stays_on_when_genes_are_found(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load()
    result = run.compute_qc_metrics()
    assert result.summary["mito_filter_switched_off"] is False
    assert run.params.filter == REFERENCE_PRESET.filter


# --- The app's default settings: thresholds suggested from the data ------------------------------

def test_default_params_suggest_the_filter_and_record_it(raw_data_path):
    from scrna_analysis.params import DEFAULT_PARAMS, SuggestionParams
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    run.load()
    result = run.compute_qc_metrics()
    suggested = result.summary["suggested_filter"]
    assert run.params.filter == FilterParams(**suggested["filter"])  # the filter to apply is the suggestion
    assert run.params.filter != REFERENCE_PRESET.filter
    assert run.history[-1]["settings"]["suggest"] == {"n_mads_genes": 5.0, "n_mads_mito": 3.0, "min_cells": 3,
                                                        "per_sample": True}  # no effect on a single dataset
    # Other MADs: a look only, nothing changes.
    looser = run.suggest_filter(SuggestionParams(n_mads_genes=8, n_mads_mito=6))
    assert looser.filter.max_pct_counts_mt > run.params.filter.max_pct_counts_mt
    assert run.params.filter == FilterParams(**suggested["filter"])
    # Applying records the MADs the thresholds were suggested with.
    run.apply_filter(looser.filter, SuggestionParams(n_mads_genes=8, n_mads_mito=6))
    assert run.history[-1]["settings"]["suggest"]["n_mads_mito"] == 6


def test_reference_preset_keeps_the_notebook_thresholds(raw_data_path):
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load()
    result = run.compute_qc_metrics()
    assert "suggested_filter" not in result.summary and run.params.filter == REFERENCE_PRESET.filter


# --- The whole analysis with the app's general defaults (no PBMC3k settings) ---------------------

def test_demo_runs_to_the_end_with_the_general_defaults(raw_data_path):
    """The general defaults, plus the two choices the user makes on Cell Type Annotation."""
    from scrna_analysis.params import DEFAULT_PARAMS, AnnotationParams, MarkerDatabaseParams
    chosen = AnnotationParams(markers=MarkerDatabaseParams(species="human", tissues=("Immune system", "Blood")))
    run = AnalysisRun(raw_data_path, replace(DEFAULT_PARAMS, annotation=chosen))
    run.run_all()
    assert set(run.status().values()) == {"done"}
    pre = run.stage_result("preprocess").summary
    assert pre["target_is_median"] and pre["target_sum"] == pre["median_counts_before"]
    assert pre["regressed"] is False
    assert run.result("reduce_dimensions").obsm["X_pca"].shape[1] == DEFAULT_PARAMS.pca.n_comps
    assert run.result("annotate").obs["cell_type_predicted"].notna().all()


def test_mouse_runs_through_marker_genes_with_the_general_defaults(mouse_10x_path):
    from scrna_analysis.params import DEFAULT_PARAMS
    run = AnalysisRun(mouse_10x_path, replace(DEFAULT_PARAMS, qc=QCParams(mito_prefix="mt-")))
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster",
                  "find_markers"):
        getattr(run, stage)()
    assert run.params.filter.max_pct_counts_mt is not None  # suggested from the data
    assert run.result("cluster").obs["leiden"].nunique() > 1


def test_nothing_matched_gives_a_tissue_hint(raw_data_path):
    """Blood cells annotated with brain cell types: every cluster stays Unresolved (nothing is
    forced), and the run suggests the tissues whose cell types do match. With fitting tissues
    there is no hint."""
    run = AnalysisRun(raw_data_path, REFERENCE_PRESET)
    run.load(); run.compute_qc_metrics(); run.apply_filter(); run.preprocess(); run.reduce_dimensions()
    run.cluster(); run.find_markers()
    markers = run.params.annotation.markers
    brain = run.annotate(replace(run.params.annotation, markers=replace(markers, tissues=("Brain",))))
    assert brain.summary["nothing_matched"] is True
    assert set(brain.summary["label_counts"]) == {annotation.UNRESOLVED}
    hint = brain.summary["tissue_hint"]
    assert hint[0]["tissue"] == "Immune system" and 1 <= len(hint) <= 3
    assert [h["n_clusters"] for h in hint] == sorted((h["n_clusters"] for h in hint), reverse=True)
    fitting = run.annotate(replace(run.params.annotation, markers=replace(markers, tissues=("Immune system", "Blood"))))
    assert fitting.summary["nothing_matched"] is False and fitting.summary["tissue_hint"] == []


def test_dataset_metadata_sets_the_annotation_settings(raw_data_path):
    """What is known about a sample (species, tissues) decides Cell Type Annotation's species and
    tissues, and is recorded as confirmed by the dataset."""
    from scrna_analysis.params import DEFAULT_PARAMS, MarkerDatabaseParams

    about = {"species": "Human", "tissues": ["Immune system", "Blood", "Blood"], "about": " blood   cells "}
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS, metadata=about)
    assert run.metadata_confirmed
    assert (run.params.annotation.markers.species, run.params.annotation.markers.tissues) == ("human", ("Blood", "Immune system"))
    assert run.record()["input"]["metadata"] == {"species": "human", "tissue": None, "tissues": ["Blood", "Immune system"],
                                                  "about": "blood cells", "confirmed_by": "dataset", "confirmed_at": None}
    chosen = replace(DEFAULT_PARAMS, annotation=replace(DEFAULT_PARAMS.annotation,
                                                        markers=MarkerDatabaseParams(species="mouse", tissues=("Brain",))))
    decided = AnalysisRun(raw_data_path, chosen, metadata=about).params.annotation.markers
    assert (decided.species, decided.tissues) == ("human", ("Blood", "Immune system"))  # the metadata decides
    run.load()
    with pytest.raises(PipelineError, match="come with its metadata and can't be changed"):
        run.confirm_metadata("mouse", ["Brain"])
    assert run.metadata["species"] == "human"
    no_metadata = AnalysisRun(raw_data_path, DEFAULT_PARAMS)
    assert no_metadata.params == DEFAULT_PARAMS and not no_metadata.metadata_confirmed and no_metadata.metadata == {}


def test_metadata_must_be_confirmed_before_the_analysis_continues(mouse_10x_path):
    """With `require_metadata` (as in the app), only loading runs until the species and tissue are
    confirmed; then the analysis continues and the confirmation is recorded."""
    from scrna_analysis.params import DEFAULT_PARAMS

    run = AnalysisRun(mouse_10x_path, DEFAULT_PARAMS, require_metadata=True)
    with pytest.raises(PipelineError, match="Load the data before confirming"):
        run.confirm_metadata("mouse", ["Brain"])
    run.load()
    with pytest.raises(PipelineError, match="Confirm the species and tissue on Data Overview first"):
        run.compute_qc_metrics()
    for species, tissues, about, message in [
        ("", ["Brain"], None, "Choose the species"),
        ("mouse", [], None, "Choose the tissue"),
        ("mouse", ["Not a tissue"], None, "Unknown tissue: Not a tissue"),
        ("mouse", ["Brain"], "x" * 301, "under 300 characters"),
    ]:
        with pytest.raises(PipelineError, match=message):
            run.confirm_metadata(species, tissues, about)
    assert not run.metadata_confirmed and run.history[-1]["stage"] == "load"  # nothing recorded

    run.confirm_metadata("mouse", ["Brain", "Immune system"], "E18 mouse brain")
    assert run.metadata_confirmed and run.metadata["confirmed_by"] == "user" and run.metadata["confirmed_at"]
    assert run.params.annotation.markers.tissues == ("Brain", "Immune system")
    assert run.history[-1]["stage"] == "metadata" and run.history[-1]["summary"]["about"] == "E18 mouse brain"
    run.compute_qc_metrics()  # now it runs
    assert run.record()["input"]["metadata"]["species"] == "mouse"


def test_other_species_need_no_tissue_and_skip_cell_types(mouse_10x_path):
    from scrna_analysis.params import DEFAULT_PARAMS

    run = AnalysisRun(mouse_10x_path, DEFAULT_PARAMS, require_metadata=True)
    run.load()
    run.confirm_metadata("other", [])
    for stage in ("compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster", "find_markers"):
        getattr(run, stage)()
    with pytest.raises(annotation.AnnotationError, match="can't be predicted for this species"):
        run.annotate()


def test_changing_the_metadata_only_reruns_cell_type_annotation(raw_data_path):
    """Species and tissues matter only to Cell Type Annotation: changing them marks it as needing
    re-running and keeps every other step; confirming the same ones changes nothing."""
    run = run_analysis(raw_data_path, REFERENCE_PRESET)
    markers = run.params.annotation.markers
    run.confirm_metadata(markers.species, list(markers.tissues))
    assert set(run.status().values()) == {"done"} and run.history[-1]["cleared_later_stages"] == []
    run.confirm_metadata("human", ["Immune system"])
    assert run.status()["annotate"] == "needs re-running"
    assert all(state == "done" for stage, state in run.status().items() if stage != "annotate")
    assert run.history[-1]["cleared_later_stages"] == ["annotate"]
    assert run.annotate().settings["annotation"]["markers"]["tissues"] == ["Immune system"]
