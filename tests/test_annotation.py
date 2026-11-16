"""Tests for `scrna_analysis.annotation` (the over-representation test).

- Known answers on a tiny dataset, unresolved outcomes, settings and edge cases.
- The calculation against independent implementations: the p-value against SciPy's Fisher's
  exact test, the overlap counts against Scanpy's `sc.tl.marker_gene_overlap`.
- The CELLxGENE Census marker tables on PBMC3k and the mouse brain data: the rules are applied (no
  expected labels; the app's predictions are never tuned to the analyst-reviewed PBMC3k labels).
- No reference for a tissue: every cluster Unresolved with status "no_reference".
- The analyst-reviewed PBMC3k labels stay a consistent record (`validation/`), checked here.

None of this is a claim of biological correctness.
"""

from __future__ import annotations

import json
import re
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scrna_analysis.annotation import (
    PREDICTION_KEY,
    RECORD_KEY,
    UNRESOLVED,
    AnnotationError,
    enrichment_table,
    predict_cell_types,
)
from scrna_analysis.markers import rank_marker_genes
from scrna_analysis.params import EnrichmentParams, MarkerParams
from tests.conftest import PROJECT_ROOT

# regress_out (used to build the shared fixtures) triggers harmless Numba performance hints.
pytestmark = pytest.mark.filterwarnings("ignore:.*is faster on contiguous arrays:Warning")


def outcome(summary) -> dict[str, tuple[str, str]]:
    return {r.cluster: (r.predicted_label, r.status) for r in summary.clusters.itertuples()}



# --- Known answers on a tiny dataset -----------------------------------------------------
#
# 3 clusters x 10 cells, 60 genes. Each cluster has 5 genes of its own that are high only there
# (a: A1-A5, b: B1-B5, c: C1-C5), so with the top 5 marker genes each cluster's top list is
# exactly its own 5 genes.

ENRICH = EnrichmentParams(top_n_genes=5)
DB_SETS = {
    "TypeA": ("A1", "A2", "A3", "A4"),          # a: 4 of 4 among the top genes
    "TypeB": ("B1", "B2", "B3", "G1"),          # b: 3 of 4 (G1 is a background gene)
    "TypeMix": ("A5", "B5", "C1"),              # at most 1 in any cluster: too few
}


def tiny_enrichment_data() -> ad.AnnData:
    own = {"a": "A", "b": "B", "c": "C"}
    genes = [f"{p}{i}" for p in "ABC" for i in range(1, 6)] + [f"G{i}" for i in range(1, 46)]
    rng = np.random.default_rng(0)
    rows, groups = [], []
    for cluster, prefix in own.items():
        for _ in range(10):
            values = rng.uniform(0, 0.1, size=len(genes))
            values[[i for i, g in enumerate(genes) if g[0] == prefix and g[1:] in "12345"]] += 4
            rows.append(values)
            groups.append(cluster)
    adata = ad.AnnData(
        X=np.array(rows, dtype=np.float32),
        obs=pd.DataFrame({"cluster": pd.Categorical(groups)}, index=[f"cell{i}" for i in range(30)]),
        var=pd.DataFrame(index=genes),
    )
    adata.uns["log1p"] = {"base": None}
    adata.raw = adata
    marked, _ = rank_marker_genes(adata, MarkerParams(groupby="cluster", use_raw=True, pts=True))
    return marked


def test_enrichment_known_answers():
    adata, summary = predict_cell_types(tiny_enrichment_data(), DB_SETS, ENRICH)
    result = outcome(summary)
    assert result == {"a": ("TypeA", "confident"), "b": ("TypeB", "confident"),
                      "c": (UNRESOLVED, "too_few_markers")}
    rows = summary.clusters.set_index("cluster")
    assert (rows.loc["a", "overlap"], rows.loc["a", "n_markers"]) == (4, 4)
    assert set(rows.loc["a", "matching_markers"].split(", ")) == {"A1", "A2", "A3", "A4"}
    assert (rows.loc["b", "overlap"], rows.loc["b", "n_markers"]) == (3, 4)
    assert set(adata.obs[PREDICTION_KEY]) == {"TypeA", "TypeB", UNRESOLVED}


def test_enrichment_p_value_is_fishers_exact_test():
    from scipy.stats import fisher_exact
    table, _, _ = enrichment_table(tiny_enrichment_data(), DB_SETS, ENRICH)
    row = table[(table.cluster == "b") & (table.cell_type == "TypeB")].iloc[0]
    n_genes, top, markers, overlap = 60, 5, 4, 3
    expected = fisher_exact([[overlap, top - overlap], [markers - overlap, n_genes - top - markers + overlap]],
                            alternative="greater").pvalue
    assert row.p_value == pytest.approx(expected, rel=1e-9)
    # Adjusted across the cell types tested for the cluster (Benjamini-Hochberg).
    from scipy.stats import false_discovery_control
    in_b = table[table.cluster == "b"]
    assert list(in_b.adj_p) == pytest.approx(list(false_discovery_control(in_b.p_value, method="bh")))


def test_enrichment_two_equally_good_cell_types_are_ambiguous():
    sets = {**DB_SETS, "TypeA2": ("A2", "A3", "A4", "A5")}
    _, summary = predict_cell_types(tiny_enrichment_data(), sets, ENRICH)
    row = summary.clusters.set_index("cluster").loc["a"]
    assert (row.predicted_label, row.status) == (UNRESOLVED, "ambiguous")
    assert {row.best_candidate, row.second_candidate} == {"TypeA", "TypeA2"}


def test_enrichment_rules_are_settings():
    data = tiny_enrichment_data()
    _, strict = predict_cell_types(data, DB_SETS, replace(ENRICH, max_adj_p=1e-12))
    assert set(strict.clusters.status) == {"not_significant", "too_few_markers"}
    _, more = predict_cell_types(data, DB_SETS, replace(ENRICH, min_overlap=4))
    assert outcome(more)["b"] == (UNRESOLVED, "too_few_markers")
    assert "TypeMix" in more.skipped_cell_types  # only 3 markers: can never reach 4


def test_enrichment_reports_missing_genes_and_untested_cell_types():
    sets = {**DB_SETS, "TypeZ": ("Z1", "Z2", "A1"), "TypeA": ("A1", "A2", "A3", "A4", "NOT_A_GENE")}
    _, summary = predict_cell_types(tiny_enrichment_data(), sets, ENRICH)
    assert summary.missing_genes == {"TypeZ": ("Z1", "Z2"), "TypeA": ("NOT_A_GENE",)}
    assert summary.skipped_cell_types == ("TypeZ",)
    assert "TypeZ" not in summary.scores.columns
    assert outcome(summary)["a"] == ("TypeA", "confident")


def test_enrichment_without_usable_marker_sets_raises():
    with pytest.raises(AnnotationError, match="No cell type has at least 3"):
        predict_cell_types(tiny_enrichment_data(), {"T": ("Z1", "Z2", "Z3")}, EnrichmentParams())


@pytest.mark.parametrize("params", [
    EnrichmentParams(top_n_genes=0), EnrichmentParams(min_overlap=0),
    EnrichmentParams(max_adj_p=0), EnrichmentParams(ambiguity_fold=0.5),
])
def test_enrichment_unusable_settings_raise(params):
    with pytest.raises(AnnotationError):
        predict_cell_types(tiny_enrichment_data(), DB_SETS, params)


def test_enrichment_input_is_not_modified_and_record_is_h5ad_safe(tmp_path):
    original = tiny_enrichment_data()
    adata, _ = predict_cell_types(original, DB_SETS, ENRICH)
    assert PREDICTION_KEY not in original.obs and RECORD_KEY not in original.uns
    record = json.loads(adata.uns[RECORD_KEY])
    assert record["method"] == "over-representation test" and record["settings"]["top_n_genes"] == 5
    adata.write_h5ad(tmp_path / "x.h5ad")
    assert json.loads(ad.read_h5ad(tmp_path / "x.h5ad").uns[RECORD_KEY]) == record


# --- Real data with the Census marker tables ------------------------------------------------

def _census_sets(species: str, tissue: str) -> dict[str, tuple[str, ...]]:
    from scrna_analysis import census_markers
    path = census_markers.table_path(species, tissue, census_markers.CensusMarkerParams.release)
    if not path.is_file():
        pytest.fail(f"Missing Census marker table: {path}. Build it with scripts/build_census_markers.py.")
    return census_markers.load_marker_table(path)[0]


def _rules_hold(summary, params: EnrichmentParams) -> None:
    for row in summary.clusters.itertuples():
        if row.status == "confident":
            assert row.predicted_label == row.best_candidate
            assert row.overlap >= params.min_overlap and row.adj_p < params.max_adj_p
            assert row.second_adj_p is None or not (
                row.second_adj_p < params.max_adj_p and row.second_adj_p < row.adj_p * params.ambiguity_fold)
        else:
            assert row.predicted_label == UNRESOLVED


def test_enrichment_overlap_equals_scanpys_marker_gene_overlap(marked_pbmc3k):
    """Scanpy's own overlap count (top 100 marker genes, its default) gives the same numbers."""
    import scanpy as sc
    from scrna_analysis.marker_db import match_to_genes
    sets = match_to_genes(_census_sets("human", "blood"), marked_pbmc3k.raw.var_names, 3).sets
    table, _, _ = enrichment_table(marked_pbmc3k, sets, EnrichmentParams())
    scanpy_counts = sc.tl.marker_gene_overlap(marked_pbmc3k, {ct: set(gs) for ct, gs in sets.items()})
    for row in table.itertuples():
        assert row.overlap == scanpy_counts.loc[row.cell_type, row.cluster], (row.cluster, row.cell_type)


def test_enrichment_on_human_blood_with_the_census_table(marked_pbmc3k):
    from scrna_analysis.marker_db import match_to_genes
    matched = match_to_genes(_census_sets("human", "blood"), marked_pbmc3k.raw.var_names, 3)
    adata, summary = predict_cell_types(marked_pbmc3k, matched.sets, EnrichmentParams())
    _rules_hold(summary, EnrichmentParams())
    assert (summary.clusters.status == "confident").any()
    assert adata.obs[PREDICTION_KEY].notna().all()


def test_enrichment_on_mouse_brain_with_the_census_table(mouse_10x_path):
    from scrna_analysis.marker_db import match_to_genes
    from scrna_analysis.params import DEFAULT_PARAMS, QCParams
    from scrna_analysis.pipeline import AnalysisRun
    run = AnalysisRun(mouse_10x_path, DEFAULT_PARAMS)
    run.load()
    run.compute_qc_metrics(QCParams(mito_prefix="mt-"), DEFAULT_PARAMS.suggest)
    run.apply_filter(); run.preprocess(); run.reduce_dimensions(); run.cluster(); run.find_markers()
    marked = run.result("markers")
    matched = match_to_genes(_census_sets("mouse", "brain"), marked.raw.var_names, 3)
    _, summary = predict_cell_types(marked, matched.sets, EnrichmentParams())
    _rules_hold(summary, EnrichmentParams())
    assert (summary.clusters.status == "confident").any()


def test_no_reference_leaves_every_cluster_unresolved(marked_pbmc3k):
    """Without a marker table for the tissue: the same outputs, every cluster "no_reference"."""
    from scrna_analysis.annotation import RECORD_KEY, no_reference
    adata, summary = no_reference(marked_pbmc3k, "No reference for this tissue.")
    groups = marked_pbmc3k.obs["leiden"].astype(str)
    assert list(summary.clusters["cluster"]) == list(marked_pbmc3k.obs["leiden"].cat.categories)
    assert set(summary.clusters["status"]) == {"no_reference"}
    assert set(summary.clusters["predicted_label"]) == {UNRESOLVED}
    assert dict(zip(summary.clusters["cluster"], summary.clusters["n_cells"])) == groups.value_counts().to_dict()
    assert (adata.obs[PREDICTION_KEY] == UNRESOLVED).all() and summary.tests.empty
    assert summary.notes == ("No reference for this tissue.",)
    assert json.loads(adata.uns[RECORD_KEY])["clusters"][0]["status"] == "no_reference"
    assert PREDICTION_KEY not in marked_pbmc3k.obs  # the input is unchanged


# --- The analyst-reviewed PBMC3k labels (validation record) ------------------------------

def test_reviewed_json_matches_review_record(reviewed_annotations):
    """The machine-readable labels must agree with the human-readable review."""
    text = (PROJECT_ROOT / "validation" / "annotation_review.md").read_text()
    comparison = text.split("## Original PBMC3k labels")[1]
    md_rows = {}
    for line in comparison.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 4 and cells[0].isdigit():
            md_rows[cells[0]] = {"original_label": cells[1], "working_label": cells[2],
                                 "change": "none" if cells[3] == "no" else cells[3]}
    decisions = dict(re.findall(r"^\| (\d+) \|.*\| \*\*(.+?)\*\*", text, flags=re.MULTILINE))

    json_rows = {c["reference_cluster"]: c for c in reviewed_annotations["clusters"]}
    assert set(md_rows) == set(json_rows) == set(decisions)
    for cluster, row in md_rows.items():
        for key in ("original_label", "working_label", "change"):
            assert json_rows[cluster][key] == row[key], (cluster, key)
        assert decisions[cluster] == json_rows[cluster]["working_label"]


def test_reviewed_labels_refer_to_current_reference(reviewed_annotations):
    import hashlib
    path = PROJECT_ROOT / "validation" / "expected" / "reference_summary.json"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert reviewed_annotations["reference"]["summary_sha256"] == digest, (
        "The reference was regenerated: the review applies to the old clusters and must be repeated."
    )


def test_prediction_before_markers_raises():
    adata = tiny_enrichment_data()
    del adata.uns["rank_genes_groups"]
    with pytest.raises(AnnotationError, match="Run rank_marker_genes before predicting"):
        predict_cell_types(adata, DB_SETS, ENRICH)
