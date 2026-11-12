"""Tests for scrna_analysis/review.py and the label review in pipeline.py.

The review is the analyst's decision, kept apart from what's computed: confirming or renaming a
cluster never changes the predictions or the clusters. Hand-made examples check the rules; a
whole analysis of PBMC3k checks that the review is recorded, kept and cleared at the right times.
"""

from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from scrna_analysis.annotation import PREDICTION_KEY
from scrna_analysis.pipeline import AnalysisRun, PipelineError, run_analysis
from scrna_analysis.review import (
    MAX_LABEL_LENGTH,
    NOT_REVIEWED,
    ClusterReview,
    ReviewError,
    decide,
    review_summary,
    reviewed_labels,
)
from tests.reference_preset import REFERENCE_PRESET

pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")


# --- The rules (hand-made examples) ----------------------------------------------------------

def test_confirm_keeps_the_prediction():
    r = decide("3", "confirm", "B cells", label="ignored", note="  CD79A high ")
    assert (r.decision, r.label, r.note) == ("confirmed", "B cells", "CD79A high")


def test_rename_takes_the_analysts_label_tidied():
    r = decide("1", "rename", "Dendritic cells", label="  CD14+   monocytes ")
    assert (r.decision, r.label) == ("renamed", "CD14+ monocytes")


@pytest.mark.parametrize("decision, label, note, message", [
    ("rename", "   ", None, "Type a label for cluster 1"),
    ("rename", "x" * (MAX_LABEL_LENGTH + 1), None, "under 60 characters"),
    ("rename", NOT_REVIEWED, None, "can't be used as a label"),
    ("confirm", None, "y" * 301, "note for cluster 1 under 300"),
    ("approve", "A", None, "Unknown review decision"),
])
def test_unusable_decisions_are_explained(decision, label, note, message):
    with pytest.raises(ReviewError, match=message):
        decide("1", decision, "B cells", label=label, note=note)


def test_reviewed_labels_per_cell():
    clusters = pd.Series(["0", "1", "0", "2", "1"])
    reviews = {"0": ClusterReview("0", "confirmed", "T cells", "", ""),
               "2": ClusterReview("2", "renamed", "Platelets", "", "")}
    labels = reviewed_labels(clusters, reviews)
    assert list(labels) == ["T cells", NOT_REVIEWED, "T cells", "Platelets", NOT_REVIEWED]
    assert list(labels.categories) == ["T cells", "Platelets", NOT_REVIEWED]  # cluster order, then not reviewed
    assert review_summary(reviews, ["0", "1", "2"]) == {"n_clusters": 3, "n_confirmed": 1, "n_renamed": 1,
                                                        "n_not_reviewed": 1}


# --- In a whole analysis -------------------------------------------------------------------------

@pytest.fixture
def annotated(raw_data_path):
    """A fresh, fully run analysis (each test changes its review)."""
    return run_analysis(raw_data_path, REFERENCE_PRESET)


def test_review_is_recorded_and_never_changes_the_computed_results(annotated):
    run = annotated
    predictions = run.result("annotate").obs[PREDICTION_KEY].copy()
    clusters = run.result("cluster").obs["leiden"].copy()
    table = run.stage_result("annotate").details["summary"].clusters
    first, second = table["cluster"].iloc[0], table["cluster"].iloc[1]

    run.review_cluster(first, "confirm", note="markers agree")
    run.review_cluster(second, "rename", label="My own label")
    assert run.review[first].label == table["predicted_label"].iloc[0]
    assert run.review[second].label == "My own label"

    # Predictions and clusters are untouched; reviewed labels are built separately.
    assert run.result("annotate").obs[PREDICTION_KEY].equals(predictions)
    assert run.result("cluster").obs["leiden"].equals(clusters)
    labels = run.reviewed_labels()
    in_first = (clusters.astype(str) == first).to_numpy()
    assert set(labels[in_first]) == {run.review[first].label}
    assert NOT_REVIEWED in set(labels)  # the other clusters

    # Recorded: in the history (with the counts) and in the record.
    entry = run.history[-1]
    assert entry["stage"] == "review" and entry["summary"]["decision"] == "renamed"
    assert entry["summary"]["n_confirmed"] == 1 and entry["summary"]["n_renamed"] == 1
    assert [r["cluster"] for r in run.record()["review"]] == [first, second]

    run.review_cluster(second, "clear")
    assert second not in run.review and run.history[-1]["summary"]["decision"] == "cleared"


def test_review_needs_annotation_and_a_real_cluster(raw_data_path, annotated):
    with pytest.raises(PipelineError, match="Run Cell Type Annotation before reviewing"):
        AnalysisRun(raw_data_path, REFERENCE_PRESET).review_cluster("0", "confirm")
    with pytest.raises(ReviewError, match="There is no cluster 'nope'"):
        annotated.review_cluster("nope", "confirm")


def test_review_survives_later_steps_but_not_new_clusters(annotated):
    run = annotated
    cluster = run.stage_result("annotate").details["summary"].clusters["cluster"].iloc[0]
    run.review_cluster(cluster, "confirm")

    run.find_markers()  # same clusters: the review stays
    run.annotate()
    assert cluster in run.review and run.review_cleared is None

    run.cluster(replace(run.params.leiden, resolution=0.5))  # new clusters: the review goes
    assert run.review == {}
    assert run.review_cleared == "The clusters were recomputed, so the label review of 1 cluster was cleared."
    assert run.history[-1]["cleared_review"] == 1
