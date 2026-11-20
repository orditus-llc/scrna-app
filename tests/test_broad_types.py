"""Tests for the broad cell types (scrna_analysis/broad_types.py): rolling Census labels up through the
Cell Ontology, and the CellTypist label tables."""

from __future__ import annotations

import csv

import pytest

from scrna_analysis import broad_types as bt
from scrna_analysis import celltypist_models as ct
from scrna_analysis import census_markers as cm
from scrna_analysis.params import CensusMarkerParams

# Every model with a broad-type table (55 drafted 2026-11-16 + the 3 written earlier).
MODELS = tuple(sorted(f"{p.stem}.pkl" for p in bt.TABLE_DIR.glob("*.csv")))
HOW = {"ontology", "pattern", "hand", "ask"}  # how a row was found; "review" must not be left


@pytest.fixture(scope="module")
def ontology():
    if not cm.ONTOLOGY_FILE.is_file():
        pytest.fail(f"Missing {cm.ONTOLOGY_FILE}. Run scripts/download_cell_ontology.py.")
    return cm.read_cell_ontology()


def term(ontology, name):
    return next(t for t, info in ontology.items() if info["name"] == name)


def test_every_broad_type_id_is_in_the_ontology(ontology):
    for name, ids in bt.BROAD_TYPES.items():
        for i in ids:
            assert i in ontology, f"{name}: {i} isn't in the Cell Ontology"


@pytest.mark.parametrize("label, expected", [
    ("naive B cell", "B cell"),
    ("IgG plasma cell", "plasma cell"),
    ("plasmablast", "plasma cell"),  # under both B cell and antibody secreting cell: the list order decides
    ("CD14-positive, CD16-negative classical monocyte", "monocyte"),
    ("CD8-positive, alpha-beta cytotoxic T cell", "T cell"),
    ("CD16-positive, CD56-dim natural killer cell, human", "NK cell"),
    ("microglial cell", "microglia"),  # microglia are macrophages in the ontology: the more specific wins
    ("hypendymal cell", "ependymal cell"),
    ("glioblast", "neural progenitor"),
    ("Purkinje cell", "neuron"),
    ("innate lymphoid cell", bt.OWN_LABEL),  # outside the list
])
def test_census_labels_roll_up_through_the_ontology(ontology, label, expected):
    assert bt.broad_type_of_term(ontology, term(ontology, label)) == expected


def test_an_unknown_term_is_its_own_label(ontology):
    assert bt.broad_type_of_term(ontology, "CL:not-a-term") == bt.OWN_LABEL


def test_the_most_specific_broad_type_wins_on_a_tiny_ontology(monkeypatch):
    onto = {"X:1": {"name": "cell", "parents": set()}, "X:2": {"name": "immune", "parents": {"X:1"}},
            "X:3": {"name": "special immune", "parents": {"X:2"}}, "X:4": {"name": "leaf", "parents": {"X:3"}}}
    monkeypatch.setattr(bt, "BROAD_TYPES", {"immune": ("X:2",), "special": ("X:3",)})
    assert bt.broad_type_of_term(onto, "X:4") == "special"
    assert bt.broad_type_of_term(onto, "X:2") == "immune"
    assert bt.broad_type_of_term(onto, "X:1") == bt.OWN_LABEL


@pytest.mark.parametrize("species, tissue", [("human", "blood"), ("mouse", "brain")])
def test_every_census_cell_type_gets_a_value(ontology, species, tissue):
    path = cm.table_path(species, tissue, CensusMarkerParams.release)
    if not path.is_file():
        pytest.fail(f"Missing {path}. Run scripts/build_census_markers.py {species} {tissue}.")
    _, provenance = cm.load_marker_table(path)
    rolled = bt.census_broad_types(cm.cell_type_terms(provenance), ontology)
    assert set(rolled.values()) <= set(bt.VALUES)
    # Nearly all of a tissue's Census types fall into the list; only a few are compared by their own name.
    assert sum(v == bt.OWN_LABEL for v in rolled.values()) <= 0.05 * len(rolled)


@pytest.mark.parametrize("model", MODELS)
def test_every_celltypist_label_is_in_its_table_once(model):
    path = ct.MODEL_DIR / model
    if not path.is_file():
        pytest.fail(f"Missing {path}. Run scripts/download_celltypist_models.py.")
    loaded, _ = ct.load_model(model)
    table = bt.celltypist_broad_types(model)
    assert table is not None, f"No broad-type table for {model}"
    assert set(table) == set(loaded.cell_types)
    assert set(table.values()) <= set(bt.VALUES)
    with open(bt.table_path(model), newline="") as f:
        rows = list(csv.DictReader(f))
    assert {r.get("how") or "hand" for r in rows} <= HOW, f"{model} has rows still to review"


def test_every_offered_model_has_a_table_and_the_drafts_are_complete():
    assert len(MODELS) >= 58
    for species in ("human", "mouse"):
        for choice in ct.model_choices(species):
            assert f"{choice['model'].removesuffix('.pkl')}.csv" in {p.name for p in bt.TABLE_DIR.glob("*.csv")}


def test_draft_spot_checks_from_the_review():
    tonsil = bt.celltypist_broad_types("Cells_Human_Tonsil.pkl")
    assert tonsil["Mature IgG+ PC"] == "plasma cell"  # "PC" = plasma cell here, not pericyte
    assert tonsil["Crypt"] == "epithelial cell"  # crypt epithelium, not the olfactory crypt neuron
    assert bt.celltypist_broad_types("Adult_Human_Vascular.pkl")["cap_pc"] == "pericyte"  # here it is
    assert bt.celltypist_broad_types("Healthy_COVID19_PBMC.pkl")["CD83_CD14_mono"] == "monocyte"  # CD83 isn't CD8
    lung = bt.celltypist_broad_types("Human_Lung_Atlas.pkl")
    assert lung["Plasmacytoid DCs"] == "dendritic cell" and lung["Monocyte-derived Mph"] == "macrophage"
    gut = bt.celltypist_broad_types("Cells_Intestinal_Tract.pkl")
    assert gut["EC cells (TAC1+)"] == "epithelial cell"  # enterochromaffin, not endothelial
    assert bt.celltypist_broad_types("Healthy_Adult_Heart.pkl")["Purkinje_cell"] == bt.OWN_LABEL  # heart muscle
    assert bt.celltypist_broad_types("Mouse_Dentate_Gyrus.pkl")["MiCajal-Retziusoglia"] == "microglia"  # garbled label


def test_analyst_decisions_on_the_unclear_labels(ontology):
    """The analyst's choices (2026-11-17) for the labels the draft couldn't settle."""
    # Myelocytes follow the ontology, so they agree with a Census "myelocyte" answer.
    census = bt.broad_type_of_term(ontology, term(ontology, "myelocyte"))
    for model, label in (("Immune_All_Low.pkl", "Myelocytes"), ("Immune_All_High.pkl", "Promyelocytes"),
                         ("Pan_Fetal_Human.pkl", "MYELOCYTE"), ("Cells_Fetal_Lung.pkl", "Myelocyte-like")):
        assert bt.celltypist_broad_types(model)[label] == census
    assert bt.celltypist_broad_types("Cells_Human_Tonsil.pkl")["FDC"] == bt.OWN_LABEL  # stromal, not immune DCs
    assert bt.celltypist_broad_types("Immune_All_High.pkl")["ILC"] == "NK cell"  # the model has no NK label
    assert bt.celltypist_broad_types("Mouse_Isocortex_Hippocampus.pkl")["SMC-Peri"] == "smooth muscle cell"


def test_celltypist_tables_spot_checks():
    immune = bt.celltypist_broad_types("Immune_All_Low.pkl")
    assert immune["Classical monocytes"] == "monocyte" and immune["Plasma cells"] == "plasma cell"
    assert immune["Mono-mac"] == bt.NO_BROAD_TYPE and immune["ILC3"] == bt.OWN_LABEL
    brain = bt.celltypist_broad_types("Mouse_Whole_Brain.pkl")
    assert brain["334 Microglia NN"] == "microglia" and brain["327 Oligo NN"] == "oligodendrocyte"
    assert brain["037 DG Glut"] == "neuron" and brain["338 Lymphoid NN"] == bt.NO_BROAD_TYPE


def test_a_model_without_a_table_and_a_broken_table(tmp_path):
    assert bt.celltypist_broad_types("No_Such_Model.pkl") is None
    assert bt.broad_type("B cells", None) is None and bt.broad_type(None, {"B cells": "B cell"}) is None
    (tmp_path / "Bad.csv").write_text("label,broad_type\nX,wizard\n")
    with pytest.raises(bt.BroadTypeError, match="isn't a broad type"):
        bt.celltypist_broad_types("Bad.pkl", tmp_path)
    (tmp_path / "Twice.csv").write_text("label,broad_type\nX,neuron\nX,neuron\n")
    with pytest.raises(bt.BroadTypeError, match="listed twice"):
        bt.celltypist_broad_types("Twice.pkl", tmp_path)
