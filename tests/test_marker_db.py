"""Tests for scrna_analysis/marker_db.py: matching marker genes to a dataset's gene names, and the
species clues read from gene names. Also checks that PanglaoDB is gone from the app (removed
2026-11-12; the frozen copies in archive/ keep it).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from scrna_analysis import census_markers
from scrna_analysis.load_data import load_dataset
from scrna_analysis.marker_db import MarkerDatabaseError, match_to_genes

PROJECT = Path(__file__).resolve().parent.parent


# --- Matching to the dataset's gene names ----------------------------------------------------

def test_markers_match_exactly_or_ignoring_case():
    match = match_to_genes({"T cells": ("CD3D", "CD3E", "IL7R")}, ["Cd3d", "CD3E", "Il7r", "Actb"], min_markers=3)
    assert match.sets == {"T cells": ("Cd3d", "CD3E", "Il7r")}
    assert match.case_matched == 2
    assert match.missing == {} and match.too_few == {}


def test_missing_markers_are_reported_and_small_sets_left_out():
    match = match_to_genes({"T cells": ("CD3D", "CD3E", "XYZ1"), "B cells": ("MS4A1", "CD79A")},
                           ["CD3D", "CD3E", "MS4A1"], min_markers=2)
    assert match.sets == {"T cells": ("CD3D", "CD3E")}
    assert match.missing == {"T cells": ("XYZ1",), "B cells": ("CD79A",)}
    assert match.too_few == {"B cells": 1}


def test_a_name_matching_several_genes_ignoring_case_is_not_guessed():
    match = match_to_genes({"T cells": ("CD3D",)}, ["Cd3d", "cd3d"], min_markers=1)
    assert match.sets == {} and match.missing == {"T cells": ("CD3D",)}


def test_min_markers_must_be_positive():
    with pytest.raises(MarkerDatabaseError, match="at least 1"):
        match_to_genes({}, [], min_markers=0)


def test_a_human_census_table_matches_mouse_genes_by_ignoring_case(mouse_10x_path):
    """Census's human blood markers (CD3D, …) are found in mouse data (Cd3d, …) by ignoring case."""
    path = census_markers.table_path("human", "blood", census_markers.CensusMarkerParams.release)
    if not path.is_file():
        pytest.fail(f"Missing Census marker table: {path}. Build it with scripts/build_census_markers.py.")
    sets, _ = census_markers.load_marker_table(path)
    adata, _ = load_dataset(mouse_10x_path)
    match = match_to_genes(sets, adata.var_names, min_markers=3)
    assert match.case_matched > 0 and len(match.sets) > 0
    for genes in match.sets.values():
        assert all(g in adata.var_names for g in genes)


def test_panglaodb_is_gone_from_the_app():
    """No code path reads PanglaoDB any more (analyst's decision, 2026-11-12)."""
    files = [*PROJECT.glob("scrna_analysis/*.py"), *PROJECT.glob("app/*.py"), *PROJECT.glob("app/srcts/**/*.ts*"),
             *PROJECT.glob("scripts/*.py")]
    pattern = ("read_panglaodb", "load_panglaodb", "PANGLAODB_FILE", "panglaodb_tissues_for", "PanglaoDB only",
               "download_marker_db")
    found = [(f.relative_to(PROJECT), word) for f in files for word in pattern if word in f.read_text()]
    assert found == []
    assert not (PROJECT / "scripts" / "download_marker_db.py").exists()


# --- What the gene names say about the species ----------------------------------------------

def test_species_clues_from_made_up_gene_names():
    from scrna_analysis.marker_db import species_clues
    human = ["CD3D", "MS4A1", "HLA-DRA", "RP11-34P13.7", "MT-CO1"]
    mouse = ["Cd3d", "Ms4a1", "Gm1992", "Rp1", "mt-Co1"]
    assert species_clues(human) == [{"species": "human", "clue": "gene names are written in capitals, like CD3D"}]
    assert species_clues(mouse, mito_prefix="mt-") == [
        {"species": "mouse", "clue": "gene names are written like Cd3d"},
        {"species": "mouse", "clue": "mitochondrial genes are named mt-…"}]
    zebrafish_ids = [f"ENSDARG{i:011d}" for i in range(5)]
    assert species_clues(["cd3d", "ms4a1", "rag1", "lck", "cd4-1"], zebrafish_ids) == [
        {"species": "zebrafish", "clue": "Ensembl IDs start with ENSDARG"}]  # lower-case symbols: no style clue
    # Data named by IDs only: the IDs speak; a mixed or unclear set says nothing.
    assert species_clues([f"ENSMUSG{i:011d}" for i in range(5)])[0]["clue"] == "Ensembl IDs start with ENSMUSG"
    assert species_clues(["CD3D", "Cd3d", "MS4A1", "Ms4a1"]) == []
    assert species_clues([]) == [] and species_clues(["CD3D"], mito_prefix="Mt-") == [
        {"species": "human", "clue": "gene names are written in capitals, like CD3D"}]  # rat prefix: no app species


def test_species_clues_agree_on_the_two_test_datasets(raw_data_path, mouse_10x_path):
    from scrna_analysis.marker_db import species_clues
    for path, species in [(raw_data_path, "human"), (mouse_10x_path, "mouse")]:
        adata, _ = load_dataset(path, path.name)
        clues = species_clues(adata.var_names, adata.var["gene_ids"], "MT-" if species == "human" else "mt-")
        assert len(clues) == 3 and {c["species"] for c in clues} == {species}


def test_cell_types_by_gene_lists_every_set_that_has_the_gene():
    from scrna_analysis.marker_db import cell_types_by_gene
    sets = {"type A": ("G1", "G2", "G2"), "type B": ("G2", "G3"), "type C": ()}
    assert cell_types_by_gene(sets) == {"G1": ["type A"], "G2": ["type A", "type B"], "G3": ["type B"]}
    assert cell_types_by_gene({}) == {}
