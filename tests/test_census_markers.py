"""Tests for the Census marker tables (scrna_analysis/census_markers.py), without the internet.

The Census fetch itself isn't run here (it needs the network and takes minutes); everything after
it is: the Cell Ontology family tree, choosing cell types, sampling, computing markers (known
answers on made-up cells), and saving/loading tables with their fingerprint.
"""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scrna_analysis import census_markers as cm
from scrna_analysis.params import CensusMarkerParams


# --- The Cell Ontology ------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ontology():
    if not cm.ONTOLOGY_FILE.is_file():
        pytest.fail(f"Missing {cm.ONTOLOGY_FILE}. Run scripts/download_cell_ontology.py.")
    return cm.read_cell_ontology()


def term(ontology, name):
    return next(t for t, info in ontology.items() if info["name"] == name)


def test_the_ontology_knows_which_cell_types_are_related(ontology):
    b, naive_b, t_cell, cd4_t = (term(ontology, n) for n in ("B cell", "naive B cell", "T cell", "CD4-positive, alpha-beta T cell"))
    assert cm.related(ontology, b, naive_b) and cm.related(ontology, naive_b, b)  # parent and child
    assert cm.related(ontology, t_cell, cd4_t)
    assert not cm.related(ontology, b, t_cell) and not cm.related(ontology, naive_b, cd4_t)
    assert term(ontology, "lymphocyte") in cm.ancestors(ontology, naive_b)
    assert cm.ontology_version().startswith("cl/releases/")


def test_a_tiny_ontology_file(tmp_path):
    path = tmp_path / "tiny.obo"
    path.write_text("format-version: 1.2\ndata-version: test/1\n\n[Term]\nid: X:1\nname: cell\n\n"
                    "[Term]\nid: X:2\nname: blood cell\nis_a: X:1 ! cell\n\n[Term]\nid: X:3\nname: B\nis_a: X:2 ! blood cell\n\n"
                    "[Typedef]\nid: part_of\nname: part of\n")
    onto = cm.read_cell_ontology(path)
    assert onto["X:3"] == {"name": "B", "parents": {"X:2"}} and cm.ancestors(onto, "X:3") == {"X:1", "X:2"}
    assert "part_of" not in onto and cm.ontology_version(path) == "test/1"
    with pytest.raises(cm.CensusMarkerError, match="download_cell_ontology.py"):
        cm.read_cell_ontology(tmp_path / "missing.obo")


# --- Choosing and sampling cell types -------------------------------------------------------------

def reference_rows():
    rows = []
    for cell_type, term_id, n, studies in [("B cell", "CL:0000236", 150, 3), ("rare cell", "CL:9", 50, 3),
                                           ("one-study cell", "CL:8", 400, 1), ("unknown", "CL:7", 500, 5)]:
        rows += [{"soma_joinid": len(rows) + i, "cell_type": cell_type, "cell_type_ontology_term_id": term_id,
                  "dataset_id": f"d{i % studies}"} for i in range(n)]
    return pd.DataFrame(rows)


def test_cell_types_need_enough_cells_from_enough_studies_and_a_real_label():
    types = cm.choose_cell_types(reference_rows(), CensusMarkerParams())
    assert types.to_dict(orient="records") == [{"cell_type": "B cell", "ontology_id": "CL:0000236", "n_cells": 150, "n_studies": 3}]


def test_sampling_is_fixed_and_capped():
    cells = reference_rows()
    types = cm.choose_cell_types(cells, CensusMarkerParams())
    params = CensusMarkerParams(cells_per_type=100)
    first, again = cm.sample_cells(cells, types, params), cm.sample_cells(cells, types, params)
    assert len(first) == 100 and np.array_equal(first, again)  # the same cells every time
    assert set(first) <= set(cells.loc[cells["cell_type"] == "B cell", "soma_joinid"])


# --- Computing markers (known answers) ------------------------------------------------------------

def made_up_reference(ontology):
    """Four cell types with known marker genes: B cell and its child naive B cell share BGENE (and
    naive B also has NAIVE), T cell has TGENE, NK has NKGENE; every cell has background genes."""
    names = ["B cell", "naive B cell", "T cell", "natural killer cell"]
    marker = {"B cell": ["BGENE"], "naive B cell": ["BGENE", "NAIVE"], "T cell": ["TGENE"], "natural killer cell": ["NKGENE"]}
    genes = ["BGENE", "NAIVE", "TGENE", "NKGENE"] + [f"BG{i}" for i in range(40)]
    rng = np.random.default_rng(0)
    rows, obs = [], []
    for name in names:
        for _ in range(60):
            counts = rng.poisson(1.0, len(genes)).astype(float)
            for g in marker[name]:
                counts[genes.index(g)] += 30
            rows.append(counts)
            obs.append({"cell_type": name, "cell_type_ontology_term_id": term(ontology, name)})
    adata = ad.AnnData(sparse.csr_matrix(np.array(rows)), obs=pd.DataFrame(obs, index=[f"c{i}" for i in range(len(obs))]))
    adata.var_names = genes
    return adata


def test_markers_are_found_and_related_types_are_not_compared(ontology):
    table = cm.compute_marker_table(made_up_reference(ontology), ontology, CensusMarkerParams(top_n_genes=3))
    top = table[table["rank"] == 1].set_index("cell_type")["gene"].to_dict()
    assert top["T cell"] == "TGENE" and top["natural killer cell"] == "NKGENE"
    # B cell and naive B cell are compared only with T and NK cells, so both keep the shared BGENE:
    # compared with each other, it would cancel out.
    for name in ("B cell", "naive B cell"):
        assert "BGENE" in set(table.loc[table["cell_type"] == name, "gene"])
        assert set(table.loc[table["cell_type"] == name, "n_compared"]) == {120}  # T + NK cells only
    assert set(table.loc[table["cell_type"] == "T cell", "n_compared"]) == {180}  # the three others
    assert (table["pct_in"] >= 0.25).all() and (table["logfoldchange"] > 0).all()
    assert table.groupby("cell_type")["rank"].max().max() <= 3


# --- Saving and loading ---------------------------------------------------------------------------

def test_tables_save_and_load_with_their_fingerprint(tmp_path, ontology):
    table = cm.compute_marker_table(made_up_reference(ontology), ontology, CensusMarkerParams(top_n_genes=2))
    path = cm.table_path("human", "peripheral blood", "2025-11-08", tmp_path)
    assert path.name == "human_peripheral_blood.tsv.gz"
    provenance = cm.save_marker_table(table, {"source": "test"}, path)
    assert provenance["n_cell_types"] == 4 and len(provenance["table_sha256"]) == 64
    again = cm.save_marker_table(table, {"source": "test"}, path)
    assert again["table_sha256"] == provenance["table_sha256"]  # the same table gives the same file
    sets, meta = cm.load_marker_table(path)
    assert sets["T cell"][0] == "TGENE" and meta["source"] == "test"
    path.with_name("human_peripheral_blood.json").write_text(json.dumps({**meta, "table_sha256": "0" * 64}))
    with pytest.raises(cm.CensusMarkerError, match="doesn't match its recorded fingerprint"):
        cm.load_marker_table(path)
    with pytest.raises(cm.CensusMarkerError, match="build_census_markers.py"):
        cm.load_marker_table(tmp_path / "none.tsv.gz")


# --- Census first, PanglaoDB second, in the analysis (made-up Census tables, known answers) -------

from dataclasses import replace  # noqa: E402

from scrna_analysis import markers as marker_tools  # noqa: E402
from scrna_analysis.params import DEFAULT_PARAMS, CensusReferenceParams  # noqa: E402
from scrna_analysis.pipeline import AnalysisRun, PipelineError  # noqa: E402

RELEASE = CensusReferenceParams().release


@pytest.fixture(scope="module")
def marked_run(raw_data_path, tmp_path_factory):
    """The demo analyzed up to Marker Genes, with an empty folder for made-up Census tables."""
    folder = tmp_path_factory.mktemp("census")
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS, census_dir=folder,
                      metadata={"species": "human", "tissue": "blood"})
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess", "reduce_dimensions", "cluster", "find_markers"):
        getattr(run, stage)()
    return run, folder


def top_genes(run, n=30):
    table = marker_tools.marker_table(run.result("markers"), n_top=n)
    return {str(c): list(rows["gene"]) for c, rows in table.groupby("group", sort=False)}


def write_table(folder, types):
    """A Census table for human blood: [(cell type, ontology id, genes), …]."""
    rows = [{"cell_type": name, "ontology_id": term, "rank": i + 1, "gene": g, "score": 1.0, "logfoldchange": 1.0,
             "pct_in": 0.5, "pct_out": 0.1, "n_cells": 100, "n_compared": 100}
            for name, term, genes in types for i, g in enumerate(genes)]
    provenance = {"census_version": RELEASE, "citation": "test",
                  "cell_types": [{"cell_type": n, "ontology_id": t} for n, t, _ in types]}
    cm.save_marker_table(pd.DataFrame(rows), provenance, cm.table_path("human", "blood", RELEASE, folder))


def annotate(run, **census):
    census = {"tissue": "blood", "use": True, **census}  # each test starts from Census + blood
    params = replace(run.params.annotation, census=replace(run.params.annotation.census, **census))
    run.annotate(params)
    return run.stage_result("annotate").summary


def test_census_labels_come_first(marked_run):
    """A Census table made from each cluster's own top genes gives every cluster its own type."""
    run, folder = marked_run
    genes = top_genes(run)
    write_table(folder, [(f"made-up type {c}", f"X:{c}", g) for c, g in genes.items()])
    summary = annotate(run)
    assert {c["cluster"]: c["predicted_label"] for c in summary["clusters"]} == {c: f"made-up type {c}" for c in genes}
    assert {c["label_source"] for c in summary["clusters"]} == {"Census"}
    assert [r["name"] for r in summary["references"]] == ["CELLxGENE Census", "PanglaoDB"]
    assert summary["database"]["role"] == "first" and summary["database"]["tissues"] == ["blood"]
    assert set(summary["by_reference"]) == {"Census", "PanglaoDB"}  # PanglaoDB's answer is kept for every cluster


def test_panglaodb_fills_only_the_clusters_census_leaves_unresolved(marked_run):
    run, folder = marked_run
    genes = top_genes(run)
    covered = list(genes)[:1]  # one type only: unrelated clusters are left Unresolved by Census
    write_table(folder, [(f"made-up type {c}", f"X:{c}", genes[c]) for c in covered])
    summary = annotate(run)
    census = {c["cluster"]: c for c in summary["by_reference"]["Census"]}
    panglao = {c["cluster"]: c for c in summary["by_reference"]["PanglaoDB"]}
    for c in summary["clusters"]:
        cluster = c["cluster"]
        if census[cluster]["predicted_label"] != "Unresolved":  # a confident Census label is never replaced
            assert (c["predicted_label"], c["label_source"]) == (census[cluster]["predicted_label"], "Census")
        elif panglao[cluster]["predicted_label"] != "Unresolved":
            assert (c["predicted_label"], c["label_source"]) == (panglao[cluster]["predicted_label"], "PanglaoDB")
            assert c["best_candidate"] == panglao[cluster]["best_candidate"]
        else:
            assert (c["predicted_label"], c["label_source"]) == ("Unresolved", None)
    assert any(c["label_source"] == "PanglaoDB" for c in summary["clusters"])  # the fallback was used
    cells = run.result().obs["cell_type_predicted"].astype(str)
    assert cells.value_counts().to_dict() == summary["label_counts"]


def test_related_cell_types_dont_make_a_cluster_ambiguous(marked_run, ontology):
    """Two Census types with the same markers: ambiguous if unrelated, not if one is a kind of the other."""
    run, folder = marked_run
    genes = top_genes(run)
    first = list(genes)[0]
    for (a, a_term), (b, b_term), expected in [
        (("B cell", term(ontology, "B cell")), ("naive B cell", term(ontology, "naive B cell")), "confident"),
        (("B cell", term(ontology, "B cell")), ("T cell", term(ontology, "T cell")), "ambiguous"),
    ]:
        write_table(folder, [(a, a_term, genes[first]), (b, b_term, genes[first])])
        census = {c["cluster"]: c for c in annotate(run)["by_reference"]["Census"]}
        assert census[first]["status"] == expected


def test_without_a_census_table_the_user_chooses_panglaodb_only(marked_run):
    run, folder = marked_run
    with pytest.raises(Exception, match="PanglaoDB only"):
        annotate(run, tissue="spleen")  # no table for spleen
    summary = annotate(run, tissue="spleen", use=False)
    assert [r["name"] for r in summary["references"]] == ["PanglaoDB"]
    assert {c["label_source"] for c in summary["clusters"]} <= {"PanglaoDB", None}


def test_the_census_tissue_fills_in_panglaodb_tissues(raw_data_path, tmp_path):
    (tmp_path / RELEASE).mkdir()
    (tmp_path / RELEASE / cm.TISSUE_LIST).write_text(json.dumps(
        {"census_version": RELEASE, "tissues": {"human": [{"tissue": "blood", "n_cells": 10}], "mouse": []}}))
    run = AnalysisRun(raw_data_path, DEFAULT_PARAMS, require_metadata=True, census_dir=tmp_path)
    run.load()
    with pytest.raises(PipelineError, match="Unknown tissue: kidney"):
        run.confirm_metadata("human", tissue="kidney")
    metadata = run.confirm_metadata("human", tissue="Blood")
    assert (metadata["tissue"], metadata["tissues"]) == ("blood", ["Blood", "Immune system"])
    assert run.params.annotation.census.tissue == "blood"
    assert cm.census_tissues("human", RELEASE, tmp_path) == [{"tissue": "blood", "n_cells": 10, "has_table": False}]
    assert cm.panglaodb_tissues_for("Not a tissue") == ()


def test_methods_and_report_name_both_references_and_cite_census(marked_run):
    import sys
    from tests.conftest import PROJECT_ROOT
    sys.path.insert(0, str(PROJECT_ROOT / "app"))
    import payloads
    import report

    run, folder = marked_run
    genes = top_genes(run)
    write_table(folder, [(f"made-up type {c}", f"X:{c}", genes[c]) for c in list(genes)[:1]])
    summary = annotate(run)
    n = sum(c["label_source"] == "PanglaoDB" for c in summary["clusters"])
    text = report.methods_text(run)
    assert "CZ CELLxGENE Discover Census (release 2025-11-08" in text and "not related to it in the Cell Ontology" in text
    assert f"Clusters left unresolved were tested in the same way against PanglaoDB" in text
    assert f"{n} clusters were labeled this way" in text or f"{n} cluster was labeled this way" in text
    html = report.run_report_html(run, "test")
    assert "Cite: test" in html and "<th>From</th>" in html
    # The page: each reference's answer per cluster, and candidates from both.
    results = payloads.annotation_payload(run)["results"]
    assert results["first_reference"] == "Census" and [r["key"] for r in results["references"]] == ["Census", "PanglaoDB"]
    assert set(results["answers"]) == set(results["candidates_by_reference"]) == {"Census", "PanglaoDB"}
    assert {c["label_source"] for c in results["clusters"]} <= {"Census", "PanglaoDB", None}
    assert results["candidates"] == results["candidates_by_reference"]["Census"]


def test_ribosomal_and_mitochondrial_genes_are_never_markers(ontology):
    names = ["RPL13A", "RPS4X", "RPLP0", "RPSA", "MRPS18B", "MRPL12", "MT-CO1", "Rpl13a", "mt-Nd1",
             "RPS6KA1", "RPL", "MTOR", "CD3E", "Mrpl15b"]
    out = cm.excluded_genes(names, CensusMarkerParams())
    assert [n for n, x in zip(names, out) if x] == ["RPL13A", "RPS4X", "RPLP0", "RPSA", "MRPS18B", "MRPL12", "MT-CO1",
                                                    "Rpl13a", "mt-Nd1", "Mrpl15b"]
    # A made-up reference where T cells' strongest difference is a ribosomal gene: it isn't a marker.
    adata = made_up_reference(ontology)
    adata.var_names = ["BGENE", "NAIVE", "RPL13A", "NKGENE"] + [f"BG{i}" for i in range(40)]  # TGENE renamed
    table = cm.compute_marker_table(adata, ontology, CensusMarkerParams(top_n_genes=3))
    assert "RPL13A" not in set(table["gene"])
