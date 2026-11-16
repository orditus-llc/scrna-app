"""One final cell type per cluster from the two checks.

The two checks are the marker check (`annotation.py`, CELLxGENE Census markers) and CellTypist
(`celltypist_models.py`). Each gives a usable answer for a cluster only when it is confident:

- marker check: status "confident" (passed its three checks);
- CellTypist: its majority label is voted for by at least `min_celltypist_share` of the cluster's
  cells, at least `min_celltypist_gene_match` of the model's genes are in the data, and the label
  has a broad type (not "none", e.g. "Undefined").

Both answers are compared at the level of broad cell types (`broad_types.py`):

| marker check | CellTypist | final label                 | agreement           |
|--------------|------------|-----------------------------|---------------------|
| confident    | confident, same broad type      | the marker check's label | agree        |
| confident    | confident, other broad type     | Unknown                  | disagree     |
| confident    | no usable answer                | the marker check's label | marker_check_only |
| no answer    | confident                       | CellTypist's label       | celltypist_only   |
| no answer    | no usable answer                | Unknown                  | neither      |

Strict on purpose: related but different broad types (e.g. B cell vs plasma cell) count as a
disagreement, and the analyst decides. "Unknown" is a valid answer; a label is never forced.
The checks' own answers are kept unchanged next to the final label.
"""

from __future__ import annotations

from typing import Mapping

import pandas as pd

from scrna_analysis.broad_types import NO_BROAD_TYPE, OWN_LABEL
from scrna_analysis.params import ConsensusParams

UNKNOWN = "Unknown"
MARKER_CHECK_KEY = "cell_type_marker_check"  # obs: the marker check's own label (the final one is cell_type_predicted)

AGREE = "agree"
DISAGREE = "disagree"
MARKER_CHECK_ONLY = "marker_check_only"
CELLTYPIST_ONLY = "celltypist_only"
NEITHER = "neither"
AGREEMENTS = (AGREE, MARKER_CHECK_ONLY, CELLTYPIST_ONLY, DISAGREE, NEITHER)

# How sure the suggestion is, from how the two checks compare (supervisor, 2026-11-13; the app's own
# rule). It only orders the review: the user confirms or renames every cluster.
HIGH, MEDIUM, LOW = "high", "medium", "low"
CONFIDENCE = {AGREE: HIGH, MARKER_CHECK_ONLY: MEDIUM, CELLTYPIST_ONLY: MEDIUM, DISAGREE: LOW, NEITHER: LOW}

# Why CellTypist gave no usable answer for a cluster.
CT_NOT_RUN = "not_run"  # no model chosen for this species and tissue
CT_ERROR = "error"  # it couldn't run (e.g. none of the model's genes in the data)
CT_NO_TABLE = "no_table"  # the model has no broad-type table, so it can't be compared
CT_LOW_GENE_MATCH = "low_gene_match"
CT_LOW_SHARE = "low_share"
CT_NO_BROAD_TYPE = "no_broad_type"  # a vague label ("Undefined", "Mono-mac")


def _key(label: str, broad: str) -> str:
    """What two answers are compared by: the broad type, or the label itself when it is outside the list."""
    return f"label:{label.strip().lower()}" if broad == OWN_LABEL else broad


def combine(marker_clusters: pd.DataFrame, marker_broad: Mapping[str, str],
            celltypist_clusters: pd.DataFrame | None, celltypist_broad: Mapping[str, str] | None,
            gene_match: float | None, params: ConsensusParams,
            celltypist_problem: str | None = None) -> pd.DataFrame:
    """The final label of every cluster.

    `marker_clusters`: the marker check's rows (cluster, n_cells, predicted_label, status);
    `marker_broad`: its label -> broad type; `celltypist_clusters`: CellTypist's rows (cluster,
    label, share) or None; `celltypist_broad`: the model's label -> broad type (None = no table);
    `gene_match`: share of the model's genes in the data; `celltypist_problem`: CT_NOT_RUN or
    CT_ERROR when CellTypist has no rows.

    One row per cluster, in the marker check's order: cluster, n_cells, final_label,
    final_broad_type, agreement, marker_label / marker_broad_type (None unless confident),
    marker_status, celltypist_label, celltypist_share, celltypist_broad_type, celltypist_confident,
    celltypist_issue (None when confident).
    """
    _check(params)
    ct_rows = ({str(r["cluster"]): r for r in celltypist_clusters.to_dict("records")}
               if celltypist_clusters is not None else {})
    rows = []
    for m in marker_clusters.to_dict("records"):
        cluster = str(m["cluster"])
        marker_ok = m["status"] == "confident"
        marker_label = m["predicted_label"] if marker_ok else None
        marker_broad_type = marker_broad.get(marker_label, OWN_LABEL) if marker_ok else None

        ct = ct_rows.get(cluster)
        ct_label = ct["label"] if ct else None
        ct_share = float(ct["share"]) if ct else None
        ct_broad_type = celltypist_broad.get(ct_label) if (ct and celltypist_broad is not None) else None
        if ct is None:
            issue = celltypist_problem or CT_NOT_RUN
        elif celltypist_broad is None:
            issue = CT_NO_TABLE
        elif gene_match is None or gene_match < params.min_celltypist_gene_match:
            issue = CT_LOW_GENE_MATCH
        elif ct_share < params.min_celltypist_share:
            issue = CT_LOW_SHARE
        elif ct_broad_type in (None, NO_BROAD_TYPE):
            issue = CT_NO_BROAD_TYPE
        else:
            issue = None
        ct_ok = issue is None

        if marker_ok and ct_ok:
            same = _key(marker_label, marker_broad_type) == _key(ct_label, ct_broad_type)
            agreement = AGREE if same else DISAGREE
            final, final_broad = (marker_label, marker_broad_type) if same else (UNKNOWN, None)
        elif marker_ok:
            agreement, final, final_broad = MARKER_CHECK_ONLY, marker_label, marker_broad_type
        elif ct_ok:
            agreement, final, final_broad = CELLTYPIST_ONLY, ct_label, ct_broad_type
        else:
            agreement, final, final_broad = NEITHER, UNKNOWN, None
        rows.append({"cluster": cluster, "n_cells": int(m["n_cells"]), "final_label": final,
                     "final_broad_type": final_broad if final_broad != OWN_LABEL else None, "agreement": agreement,
                     "marker_label": marker_label, "marker_broad_type": marker_broad_type, "marker_status": m["status"],
                     "celltypist_label": ct_label, "celltypist_share": ct_share, "celltypist_broad_type": ct_broad_type,
                     "celltypist_confident": ct_ok, "celltypist_issue": issue})
    return pd.DataFrame(rows, columns=["cluster", "n_cells", "final_label", "final_broad_type", "agreement",
                                       "marker_label", "marker_broad_type", "marker_status", "celltypist_label",
                                       "celltypist_share", "celltypist_broad_type", "celltypist_confident",
                                       "celltypist_issue"])


def label_counts(final: pd.DataFrame) -> dict[str, int]:
    """Final label -> number of cells (largest first)."""
    counts = final.groupby("final_label", sort=False)["n_cells"].sum().sort_values(ascending=False, kind="stable")
    return {str(k): int(v) for k, v in counts.items()}


def _check(params: ConsensusParams) -> None:
    for name in ("min_celltypist_share", "min_celltypist_gene_match"):
        value = getattr(params, name)
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must be between 0 and 1 (got {value}).")
