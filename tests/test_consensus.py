"""Tests for the final label per cluster (scrna_analysis/consensus.py, Step B): the rule on made-up
answers, then the whole Cell Type Annotation stage on the demo data."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from scrna_analysis import annotation, celltypist_models as ct, consensus as cs
from scrna_analysis.broad_types import NO_BROAD_TYPE, OWN_LABEL
from scrna_analysis.params import DEFAULT_PARAMS, ConsensusParams
from scrna_analysis.pipeline import AnalysisRun

APP_DIR = Path(__file__).resolve().parent.parent / "app"
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

P = ConsensusParams()


def marker(*rows):
    return pd.DataFrame([{"cluster": c, "n_cells": 10, "predicted_label": label if status == "confident" else "Unresolved",
                          "status": status} for c, label, status in rows])


def celltypist(*rows):
    return pd.DataFrame([{"cluster": c, "label": label, "share": share} for c, label, share in rows])


MARKER_BROAD = {"naive B cell": "B cell", "IgG plasma cell": "plasma cell", "classical monocyte": "monocyte",
                "innate lymphoid cell": OWN_LABEL}
CT_BROAD = {"B cells": "B cell", "Monocytes": "monocyte", "Mono-mac": NO_BROAD_TYPE, "ILC": OWN_LABEL,
            "Innate lymphoid cell": OWN_LABEL}


def test_the_five_cases():
    final = cs.combine(
        marker(("0", "naive B cell", "confident"), ("1", "IgG plasma cell", "confident"),
               ("2", "classical monocyte", "confident"), ("3", "x", "ambiguous"), ("4", "x", "too_few_markers")),
        MARKER_BROAD,
        celltypist(("0", "B cells", 0.9), ("1", "B cells", 0.95), ("2", "Monocytes", 0.3), ("3", "Monocytes", 0.8),
                   ("4", "Monocytes", 0.2)),
        CT_BROAD, gene_match=0.7, params=P).set_index("cluster")
    assert final.loc["0", ["final_label", "agreement", "final_broad_type"]].tolist() == ["naive B cell", cs.AGREE, "B cell"]
    # Related but different broad types count as a disagreement (strict).
    assert final.loc["1", ["final_label", "agreement"]].tolist() == [cs.UNKNOWN, cs.DISAGREE]
    assert final.loc["2", ["final_label", "agreement", "celltypist_issue"]].tolist() == [
        "classical monocyte", cs.MARKER_CHECK_ONLY, cs.CT_LOW_SHARE]
    assert final.loc["3", ["final_label", "agreement"]].tolist() == ["Monocytes", cs.CELLTYPIST_ONLY]
    assert final.loc["4", ["final_label", "agreement"]].tolist() == [cs.UNKNOWN, cs.NEITHER]
    # The checks' own answers are kept next to the final label.
    assert final.loc["3", "marker_status"] == "ambiguous" and final.loc["3", "marker_label"] is None
    assert final.loc["1", "celltypist_label"] == "B cells" and final.loc["1", "marker_label"] == "IgG plasma cell"


def test_the_share_and_gene_match_limits_are_inclusive():
    m = marker(("0", "naive B cell", "ambiguous"))
    at = cs.combine(m, MARKER_BROAD, celltypist(("0", "B cells", 0.5)), CT_BROAD, 0.5, P)
    assert at.loc[0, "agreement"] == cs.CELLTYPIST_ONLY
    below = cs.combine(m, MARKER_BROAD, celltypist(("0", "B cells", 0.5)), CT_BROAD, 0.49, P)
    assert below.loc[0, ["agreement", "celltypist_issue"]].tolist() == [cs.NEITHER, cs.CT_LOW_GENE_MATCH]
    stricter = cs.combine(m, MARKER_BROAD, celltypist(("0", "B cells", 0.6)), CT_BROAD, 0.9, ConsensusParams(0.7, 0.5))
    assert stricter.loc[0, "celltypist_issue"] == cs.CT_LOW_SHARE


def test_vague_labels_missing_tables_and_no_celltypist():
    m = marker(("0", "classical monocyte", "confident"))
    vague = cs.combine(m, MARKER_BROAD, celltypist(("0", "Mono-mac", 0.9)), CT_BROAD, 0.9, P)
    assert vague.loc[0, ["agreement", "celltypist_issue"]].tolist() == [cs.MARKER_CHECK_ONLY, cs.CT_NO_BROAD_TYPE]
    no_table = cs.combine(m, MARKER_BROAD, celltypist(("0", "Monocytes", 0.9)), None, 0.9, P)
    assert no_table.loc[0, "celltypist_issue"] == cs.CT_NO_TABLE
    not_run = cs.combine(m, MARKER_BROAD, None, None, None, P)
    assert not_run.loc[0, ["agreement", "celltypist_issue"]].tolist() == [cs.MARKER_CHECK_ONLY, cs.CT_NOT_RUN]
    failed = cs.combine(m, MARKER_BROAD, None, None, None, P, celltypist_problem=cs.CT_ERROR)
    assert failed.loc[0, "celltypist_issue"] == cs.CT_ERROR


def test_types_outside_the_list_are_compared_by_name():
    m = marker(("0", "innate lymphoid cell", "confident"), ("1", "innate lymphoid cell", "confident"))
    final = cs.combine(m, MARKER_BROAD, celltypist(("0", "Innate lymphoid cell", 0.9), ("1", "ILC", 0.9)),
                       CT_BROAD, 0.9, P).set_index("cluster")
    assert final.loc["0", "agreement"] == cs.AGREE and final.loc["0", "final_broad_type"] is None
    assert final.loc["1", "agreement"] == cs.DISAGREE


def test_no_reference_and_label_counts():
    m = pd.DataFrame([{"cluster": "0", "n_cells": 30, "predicted_label": "Unresolved", "status": "no_reference"},
                      {"cluster": "1", "n_cells": 5, "predicted_label": "Unresolved", "status": "no_reference"}])
    final = cs.combine(m, {}, celltypist(("0", "B cells", 0.9), ("1", "B cells", 0.1)), CT_BROAD, 0.9, P)
    assert final["final_label"].tolist() == ["B cells", cs.UNKNOWN]
    assert cs.label_counts(final) == {"B cells": 30, cs.UNKNOWN: 5}
    with pytest.raises(ValueError, match="between 0 and 1"):
        cs.combine(m, {}, None, None, None, ConsensusParams(min_celltypist_share=2))


# --- The whole stage on the demo data --------------------------------------------------------------

@pytest.fixture(scope="module")
def annotated(raw_data_path):
    if not (ct.MODEL_DIR / "Immune_All_Low.pkl").is_file():
        pytest.fail("Missing the CellTypist model. Run scripts/download_celltypist_models.py.")
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS, metadata={"species": "human", "tissue": "blood"})
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster",
                  "find_markers", "annotate"):
        getattr(run, stage)()
    return run


def test_the_stage_follows_the_rule(annotated):
    result = annotated.stage_result("annotate")
    final = result.details["final"]
    marker_rows = result.details["summary"].clusters.set_index("cluster")
    ct_rows = result.details["celltypist"].clusters.set_index("cluster")
    assert final["cluster"].tolist() == marker_rows.index.tolist()
    for row in final.itertuples():
        m_ok = marker_rows.loc[row.cluster, "status"] == "confident"
        c_ok = ct_rows.loc[row.cluster, "share"] >= 0.5 and row.celltypist_broad_type not in (None, NO_BROAD_TYPE)
        if m_ok and c_ok:
            assert row.agreement in (cs.AGREE, cs.DISAGREE)
        else:
            assert row.agreement == {(True, False): cs.MARKER_CHECK_ONLY, (False, True): cs.CELLTYPIST_ONLY,
                                     (False, False): cs.NEITHER}[(m_ok, c_ok)]
        expected = {cs.AGREE: row.marker_label, cs.MARKER_CHECK_ONLY: row.marker_label,
                    cs.CELLTYPIST_ONLY: row.celltypist_label}.get(row.agreement, cs.UNKNOWN)
        assert row.final_label == expected


def test_the_stage_keeps_both_labels_on_the_cells(annotated):
    result = annotated.stage_result("annotate")
    obs = annotated.result("annotate").obs
    final = result.details["final"].set_index("cluster")["final_label"]
    marker_labels = result.details["summary"].clusters.set_index("cluster")["predicted_label"]
    clusters = obs[result.details["summary"].groupby].astype(str)
    assert (obs[annotation.PREDICTION_KEY].astype(str) == clusters.map(final)).all()
    assert (obs[cs.MARKER_CHECK_KEY].astype(str) == clusters.map(marker_labels)).all()
    summary = result.summary
    assert sum(summary["label_counts"].values()) == obs.shape[0] == sum(summary["marker_label_counts"].values())
    assert sum(summary["agreement_counts"].values()) == len(final)
    record = json.loads(annotated.result("annotate").uns[annotation.RECORD_KEY])
    assert record["consensus"]["settings"] == {"min_celltypist_share": 0.5, "min_celltypist_gene_match": 0.5}
    assert [r["final_label"] for r in record["consensus"]["clusters"]] == final.tolist()


def test_the_review_confirms_the_final_label(annotated):
    final = annotated.stage_result("annotate").details["final"]
    row = final.iloc[0]
    decided = annotated.review_cluster(row["cluster"], "confirm")
    assert decided.label == row["final_label"]
    annotated.review_cluster(row["cluster"], "clear")


# --- Page data, report and the developmental stage ------------------------------------------------

def test_page_data_shows_the_final_label(annotated):
    import payloads
    from report import run_report_html

    page = payloads.annotation_payload(annotated)["results"]
    final = annotated.stage_result("annotate").details["final"]
    assert [f["final_label"] for f in page["final"]] == final["final_label"].tolist()
    assert all(f["agreement_text"] == payloads.AGREEMENT_TEXT[f["agreement"]] and f["reason"] for f in page["final"])
    assert page["label_counts"] == cs.label_counts(final) and page["unknown_label"] == cs.UNKNOWN
    # The marker check's own rows are still there for the evidence.
    assert [c["cluster"] for c in page["clusters"]] == final["cluster"].tolist()
    explore = payloads.explore_payload(annotated)["results"]["predictions"]
    assert [c["predicted_label"] for c in explore["clusters"]] == final["final_label"].tolist()
    assert explore["unresolved_label"] == cs.UNKNOWN
    html = run_report_html(annotated, "demo")
    assert "Agreement" in html and "broad cell types" in html
    assert "Suggested cell type" in html and "Confidence" in html


def test_confidence_and_review_groups():
    """Confidence follows the agreement (both agree = High, one check = Medium, disagree or neither =
    Low), and every agreement belongs to exactly one "To review" group, in the agreed order."""
    import payloads

    assert cs.CONFIDENCE == {cs.AGREE: cs.HIGH, cs.MARKER_CHECK_ONLY: cs.MEDIUM, cs.CELLTYPIST_ONLY: cs.MEDIUM,
                             cs.DISAGREE: cs.LOW, cs.NEITHER: cs.LOW}
    assert [g["key"] for g in payloads.REVIEW_GROUPS] == ["disagree", "one_check", "agree", "no_suggestion"]
    grouped = [a for g in payloads.REVIEW_GROUPS for a in g["agreements"]]
    assert sorted(grouped) == sorted(cs.AGREEMENTS)
    final = cs.combine(
        marker(("0", "naive B cell", "confident"), ("1", "IgG plasma cell", "confident"), ("2", "classical monocyte", "confident"),
               ("3", "x", "too_few_markers"), ("4", "y", "too_few_markers")),
        MARKER_BROAD,
        celltypist(("0", "B cells", 0.9), ("1", "B cells", 0.95), ("2", "Monocytes", 0.3), ("3", "B cells", 0.9),
                   ("4", "Mono-mac", 0.9)),
        CT_BROAD, 0.9, P)
    rows = payloads.final_rows(final)
    assert [(r["agreement"], r["confidence_text"], r["review_group"]) for r in rows] == [
        (cs.AGREE, "High", "agree"), (cs.DISAGREE, "Low", "disagree"), (cs.MARKER_CHECK_ONLY, "Medium", "one_check"),
        (cs.CELLTYPIST_ONLY, "Medium", "one_check"), (cs.NEITHER, "Low", "no_suggestion")]


def test_page_sends_the_review_groups(annotated):
    import payloads

    page = payloads.annotation_payload(annotated)["results"]
    assert [g["key"] for g in page["review_groups"]] == [g["key"] for g in payloads.REVIEW_GROUPS]
    assert all(g["label"] and g["description"] for g in page["review_groups"])
    keys = {g["key"] for g in page["review_groups"]}
    assert all(f["review_group"] in keys and f["confidence"] in ("high", "medium", "low") for f in page["final"])


def test_reasons_say_why():
    import payloads

    final = cs.combine(
        marker(("0", "naive B cell", "confident"), ("1", "IgG plasma cell", "confident"), ("2", "classical monocyte", "confident"),
               ("3", "x", "too_few_markers")),
        MARKER_BROAD,
        celltypist(("0", "B cells", 0.9), ("1", "B cells", 0.95), ("2", "Monocytes", 0.3), ("3", "Mono-mac", 0.9)),
        CT_BROAD, 0.9, P)
    reasons = [r["reason"] for r in payloads.final_rows(final)]
    assert "same broad type: B cell" in reasons[0]
    assert "IgG plasma cell (plasma cell)" in reasons[1] and "B cell" in reasons[1]
    assert "only 30% of the cells" in reasons[2]
    assert reasons[3].startswith("Neither check is confident") and "too vague" in reasons[3]


def test_the_stage_picks_the_celltypist_model(raw_data_path):
    import payloads

    assert ct.suggest_model("mouse", "brain") == "Mouse_Whole_Brain.pkl"
    assert ct.suggest_model("mouse", "brain", "developing") == "Developing_Mouse_Brain.pkl"
    assert ct.suggest_model("human", "blood", "developing") == "Immune_All_Low.pkl"  # no developing model: the usual one
    assert payloads.metadata_from_browser({"species": "mouse", "tissue": "brain", "stage": "developing"})[3] == "developing"
    with pytest.raises(payloads.InputError, match="adult or developing"):
        payloads.metadata_from_browser({"species": "mouse", "tissue": "brain", "stage": "teenage"})
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS, require_metadata=True)
    run.load()
    run.confirm_metadata("mouse", "brain", stage="developing")
    assert run.params.annotation.celltypist.model == "Developing_Mouse_Brain.pkl"
    assert run.metadata["stage"] == "developing" and "Developing" in payloads.metadata_text(run.metadata)
    run.confirm_metadata("mouse", "brain")
    assert run.params.annotation.celltypist.model == "Mouse_Whole_Brain.pkl" and run.metadata["stage"] is None
