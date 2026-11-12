"""The analyst's review of the predicted cell types: confirm or rename each cluster's label.

A review is a human decision, not a calculation, so it is kept apart from everything computed:
the predictions (`obs["cell_type_predicted"]`, `annotation.py`) and the clusters
(`obs["leiden"]`) are never changed. Reviewed labels are built from the review when they are
needed (`reviewed_labels`), and the review itself is recorded with each decision's date.

Each reviewed cluster is either:

- confirmed: the analyst agrees with the prediction (its label is the prediction at that time)
- renamed:   the analyst gives a label of their own (free text, e.g. "CD14+ monocytes")

Clusters without a decision are "Not reviewed". Nothing here knows about a particular dataset.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping

import pandas as pd

NOT_REVIEWED = "Not reviewed"  # label of cells in clusters without a decision
REVIEWED_KEY = "cell_type_reviewed"  # obs column for the reviewed labels (in exports)
MAX_LABEL_LENGTH = 60
MAX_NOTE_LENGTH = 300


class ReviewError(ValueError):
    """A review decision can't be saved. The message says why."""


@dataclass(frozen=True)
class ClusterReview:
    cluster: str
    decision: Literal["confirmed", "renamed"]
    label: str  # the confirmed prediction, or the analyst's own label
    note: str  # the analyst's reason (may be empty)
    reviewed_at: str  # date and time (UTC)


def decide(
    cluster: str,
    decision: str,
    predicted_label: str,
    label: str | None = None,
    note: str | None = None,
    reviewed_at: str = "",
) -> ClusterReview:
    """One checked decision. `decision` is "confirm" (keeps `predicted_label`) or "rename"
    (needs `label`). Raises `ReviewError` with a plain-language message."""
    note = (note or "").strip()
    if len(note) > MAX_NOTE_LENGTH:
        raise ReviewError(f"Keep the note for cluster {cluster} under {MAX_NOTE_LENGTH} characters.")
    if decision == "confirm":
        return ClusterReview(cluster, "confirmed", predicted_label, note, reviewed_at)
    if decision != "rename":
        raise ReviewError(f"Unknown review decision '{decision}'; choose confirm or rename.")
    label = " ".join((label or "").split())  # trimmed, single spaces
    if not label:
        raise ReviewError(f"Type a label for cluster {cluster}, or confirm the prediction.")
    if len(label) > MAX_LABEL_LENGTH:
        raise ReviewError(f"Keep the label for cluster {cluster} under {MAX_LABEL_LENGTH} characters.")
    if label == NOT_REVIEWED:
        raise ReviewError(f"'{NOT_REVIEWED}' can't be used as a label.")
    return ClusterReview(cluster, "renamed", label, note, reviewed_at)


def reviewed_labels(cell_clusters: pd.Series, reviews: Mapping[str, ClusterReview]) -> pd.Categorical:
    """Each cell's reviewed label (its cluster's label, or "Not reviewed"), in cell order.
    Categories: the reviewed labels in cluster order, then "Not reviewed"."""
    clusters = cell_clusters.astype(str)
    labels = clusters.map(lambda c: reviews[c].label if c in reviews else NOT_REVIEWED)
    order = list(dict.fromkeys(reviews[c].label for c in pd.unique(clusters) if c in reviews))
    return pd.Categorical(labels, categories=order + ([NOT_REVIEWED] if NOT_REVIEWED not in order else []))


def review_summary(reviews: Mapping[str, ClusterReview], clusters: list[str]) -> dict[str, int]:
    """How many clusters are confirmed, renamed and not reviewed."""
    decisions = [reviews[c].decision for c in clusters if c in reviews]
    return {"n_clusters": len(clusters), "n_confirmed": decisions.count("confirmed"),
            "n_renamed": decisions.count("renamed"), "n_not_reviewed": len(clusters) - len(decisions)}
