"""Tests for scrna_analysis/marker_db.py: reading PanglaoDB and choosing marker sets.

The selection rules are tested on a small hand-made file (every row's fate is known), then the
real PanglaoDB file is checked against the same filter written out independently with pandas.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from scrna_analysis.load_data import load_dataset
from scrna_analysis.marker_db import (
    MarkerDatabaseError, cell_types_by_tissue, tissues_of_cell_types, match_to_genes, read_panglaodb,
    select_marker_sets, tissues,
)
from scrna_analysis.params import MarkerDatabaseParams

HEADER = ["species", "official gene symbol", "cell type", "nicknames", "ubiquitousness index",
          "product description", "gene type", "canonical marker", "germ layer", "organ",
          "sensitivity_human", "sensitivity_mouse", "specificity_human", "specificity_mouse"]

# species, gene, cell type, canonical, organ, sensitivity human, sensitivity mouse
ROWS = [
    ("Mm Hs", "CD3D", "T cells", "1", "Immune system", "0.9", "0.8"),
    ("Mm Hs", "CD3E", "T cells", "1", "Immune system", "0.7", "0.6"),
    ("Mm Hs", "IL7R", "T cells", "1", "Immune system", "0.95", "0.2"),   # mouse: too low
    ("Mm Hs", "LCK", "T cells", "", "Immune system", "0.99", "0.99"),    # not canonical
    ("Mm Hs", "CD2", "T cells", "1", "Immune system", "0.5", "0.5"),     # not above 0.5
    ("Hs", "TRAC", "T cells", "1", "Immune system", "0.8", "0"),         # human only
    ("Mm", "Thy1", "T cells", "1", "Immune system", "0", "0.9"),         # mouse only
    ("Mm Hs", "MS4A1", "B cells", "1", "Immune system", "0.6", "0.7"),
    ("Mm Hs", "GFAP", "Astrocytes", "1", "Brain", "1", "1"),
    ("Mm Hs", "AQP4", "Astrocytes", "1", "Brain", "0.8", "0.9"),
]


@pytest.fixture
def tiny_db(tmp_path: Path):
    rows = [dict(zip(HEADER, [s, g, ct, "", "0.01", "", "protein-coding gene", c, "", o, sh, sm, "0.01", "0.01"]))
            for s, g, ct, c, o, sh, sm in ROWS]
    path = tmp_path / "PanglaoDB_markers_27_Mar_2020.tsv"
    pd.DataFrame(rows, columns=HEADER).to_csv(path, sep="\t", index=False)
    return read_panglaodb(path)


# --- Reading ---------------------------------------------------------------------------------

def test_reads_the_file_and_its_version(tiny_db):
    assert tiny_db.name == "PanglaoDB"
    assert tiny_db.version == "27 Mar 2020"
    assert len(tiny_db.sha256) == 64
    assert len(tiny_db.table) == len(ROWS)
    assert tissues(tiny_db) == ("Brain", "Immune system")


def test_a_missing_file_says_how_to_get_it(tmp_path):
    with pytest.raises(MarkerDatabaseError, match="download_marker_db.py"):
        read_panglaodb(tmp_path / "absent.tsv.gz")


def test_a_file_without_the_expected_columns_is_refused(tmp_path):
    path = tmp_path / "other.tsv"
    path.write_text("gene\tcell\nCD3D\tT\n")
    with pytest.raises(MarkerDatabaseError, match="isn't a PanglaoDB marker file"):
        read_panglaodb(path)


# --- Selecting marker sets -------------------------------------------------------------------

def test_human_selection_keeps_canonical_markers_above_the_sensitivity(tiny_db):
    sets = select_marker_sets(tiny_db, MarkerDatabaseParams(species="human"))
    # Ordered by sensitivity; LCK (not canonical), CD2 (= 0.5) and Thy1 (mouse only) left out.
    assert sets == {
        "Astrocytes": ("GFAP", "AQP4"),
        "B cells": ("MS4A1",),
        "T cells": ("IL7R", "CD3D", "TRAC", "CD3E"),
    }


def test_mouse_selection_uses_the_mouse_sensitivity_and_species(tiny_db):
    sets = select_marker_sets(tiny_db, MarkerDatabaseParams(species="mouse"))
    assert sets["T cells"] == ("Thy1", "CD3D", "CD3E")  # IL7R too low in mouse; TRAC human only


def test_tissues_limit_the_cell_types(tiny_db):
    sets = select_marker_sets(tiny_db, MarkerDatabaseParams(species="human", tissues=("Brain",)))
    assert list(sets) == ["Astrocytes"]


def test_cell_types_by_tissue_follow_the_selection_rule(tiny_db):
    # Every tissue is listed; the chosen tissues are ignored (the page shows them all).
    human = cell_types_by_tissue(tiny_db, MarkerDatabaseParams(species="human", tissues=("Brain",)))
    assert human == {"Brain": ("Astrocytes",), "Immune system": ("B cells", "T cells")}
    strict = MarkerDatabaseParams(species="mouse", min_sensitivity=0.85)  # only Thy1, GFAP, AQP4 pass
    assert cell_types_by_tissue(tiny_db, strict) == {"Brain": ("Astrocytes",), "Immune system": ("T cells",)}
    with pytest.raises(MarkerDatabaseError, match="Choose a species"):
        cell_types_by_tissue(tiny_db, MarkerDatabaseParams())


def test_tissues_of_each_cell_type(tiny_db):
    assert tissues_of_cell_types(tiny_db, MarkerDatabaseParams(species="human", tissues=("Brain",))) == {
        "Astrocytes": ("Brain",), "B cells": ("Immune system",), "T cells": ("Immune system",)}


def test_the_rules_are_settings(tiny_db):
    loose = MarkerDatabaseParams(species="human", canonical_only=False, min_sensitivity=0.0)
    assert "LCK" in select_marker_sets(tiny_db, loose)["T cells"]
    assert "CD2" in select_marker_sets(tiny_db, loose)["T cells"]


@pytest.mark.parametrize("params, message", [
    (MarkerDatabaseParams(), "Choose a species"),
    (MarkerDatabaseParams(species="zebrafish"), "human and mouse only"),
    (MarkerDatabaseParams(species="human", tissues=("Moon",)), "no tissue named Moon"),
])
def test_unusable_settings_are_explained(tiny_db, params, message):
    with pytest.raises(MarkerDatabaseError, match=message):
        select_marker_sets(tiny_db, params)


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


# --- The real PanglaoDB file -----------------------------------------------------------------

def test_real_file_selection_equals_the_filter_written_out(panglaodb_path):
    """The decoupler tutorial's filter, written independently with pandas on the raw file."""
    db = read_panglaodb(panglaodb_path)
    assert db.version == "27 Mar 2020"
    assert {"Immune system", "Brain"} <= set(tissues(db))
    raw = pd.read_csv(panglaodb_path, sep="\t")
    for species, code in [("human", "Hs"), ("mouse", "Mm")]:
        sensitivity = raw[f"sensitivity_{species}"]
        expected = raw[raw["species"].fillna("").str.contains(code) & (raw["canonical marker"] == 1)
                       & (sensitivity > 0.5) & (raw["organ"] == "Immune system")]
        sets = select_marker_sets(db, MarkerDatabaseParams(species=species, tissues=("Immune system",)))
        assert {(ct, g) for ct, gs in sets.items() for g in gs} == set(
            zip(expected["cell type"], expected["official gene symbol"]))


def test_mouse_data_finds_database_markers_by_ignoring_case(panglaodb_path, mouse_10x_path):
    db = read_panglaodb(panglaodb_path)
    adata, _ = load_dataset(mouse_10x_path)
    sets = select_marker_sets(db, MarkerDatabaseParams(species="mouse", tissues=("Brain",)))
    match = match_to_genes(sets, adata.var_names, min_markers=3)
    assert match.case_matched > 0
    assert len(match.sets) > 0
    for cell_type, genes in match.sets.items():
        assert all(g in adata.var_names for g in genes)


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
