"""Broad cell types: the common level at which the two cell-type checks are compared.

The marker check (CELLxGENE Census) and CellTypist name cell types at different levels of detail and
in different vocabularies ("CD14-positive, CD16-negative classical monocyte" vs "Classical
monocytes"). Both answers are rolled up to one of `BROAD_TYPES` before they are compared:

- **Census labels** carry a Cell Ontology ID, so they are rolled up through the ontology: the broad
  type is the most specific broad type at or above the label ("naive B cell" -> B cell).
- **CellTypist labels** have no ontology IDs, so each model has a table (`broad_type_tables/
  <model>.csv`: label, broad_type). The tables were drafted from the labels' names and reviewed by
  hand; a model without a table can't be compared.

An answer can also be outside the list:
- `OWN_LABEL` ("own label"): a specific cell type the list doesn't cover (e.g. innate lymphoid cell);
  it is compared by its own name.
- `NO_BROAD_TYPE` ("none"): too vague or mixed to compare ("Undefined", "Mono-mac"); it counts as no
  usable answer from that check.

Nothing here is dataset-specific: the list, the ontology and the tables are reference data.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

from scrna_analysis.census_markers import ancestors

TABLE_DIR = Path(__file__).resolve().parent / "broad_type_tables"
OWN_LABEL = "own label"
NO_BROAD_TYPE = "none"

# Broad type -> its Cell Ontology IDs. The order breaks ties when a label sits under two broad types
# that aren't above one another (the ontology allows several parents): the earlier one wins
# (plasmablast is under both "antibody secreting cell" and "B cell" -> plasma cell).
BROAD_TYPES: dict[str, tuple[str, ...]] = {
    # Blood and immune cells
    "T cell": ("CL:0000084",),
    "NK cell": ("CL:0000623",),
    "plasma cell": ("CL:0000946",),  # antibody secreting cell: plasma cells and plasmablasts
    "B cell": ("CL:0000236",),
    "monocyte": ("CL:0000576",),
    "macrophage": ("CL:0000235",),
    "dendritic cell": ("CL:0000451",),
    "mast cell": ("CL:0000097",),
    "granulocyte": ("CL:0000094",),
    "megakaryocyte/platelet": ("CL:0000556", "CL:0000233"),
    "erythroid cell": ("CL:0000764", "CL:0000232"),
    "blood stem/progenitor cell": ("CL:0000037", "CL:0000837", "CL:0008001"),
    # Brain
    "neuron": ("CL:0000540",),
    "neural progenitor": ("CL:0011020", "CL:0000681", "CL:0000031", "CL:0000030"),  # + radial glia, neuroblast, glioblast
    "astrocyte": ("CL:0000127",),
    "oligodendrocyte": ("CL:0000128",),
    "oligodendrocyte precursor cell": ("CL:0002453",),
    "microglia": ("CL:0000129",),
    "ependymal cell": ("CL:0000065", "CL:4023181"),  # + hypendymal cell (the ontology puts it under epithelial cell)
    # Any tissue
    "endothelial cell": ("CL:0000115",),
    "fibroblast": ("CL:0000057",),
    "smooth muscle cell": ("CL:0000192",),
    "pericyte": ("CL:0000669",),
    "epithelial cell": ("CL:0000066",),
}
VALUES = (*BROAD_TYPES, OWN_LABEL, NO_BROAD_TYPE)  # what a label can roll up to


class BroadTypeError(ValueError):
    """A broad-type table can't be read. The message says why."""


def broad_type_of_term(ontology: dict[str, dict], term_id: str) -> str:
    """The broad type of a Cell Ontology term: the most specific broad type at or above it, or
    `OWN_LABEL` if none is (or the term isn't in the ontology)."""
    above = {term_id} | ancestors(ontology, term_id)
    candidates = [(name, i) for name, ids in BROAD_TYPES.items() for i in ids if i in above]
    # Keep the most specific: drop a candidate when another candidate lies below it.
    specific = [(name, i) for name, i in candidates
                if not any(i in ancestors(ontology, j) for _, j in candidates if j != i)]
    return specific[0][0] if specific else OWN_LABEL


def census_broad_types(terms: dict[str, str], ontology: dict[str, dict]) -> dict[str, str]:
    """Census cell type -> broad type, from each type's Cell Ontology ID (`census_markers.cell_type_terms`)."""
    return {label: broad_type_of_term(ontology, term) for label, term in terms.items()}


def table_path(model: str, folder: str | Path = TABLE_DIR) -> Path:
    return Path(folder) / f"{Path(model).stem}.csv"


@lru_cache(maxsize=8)
def _read(path: str) -> dict[str, str]:
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    table: dict[str, str] = {}
    for row in rows:
        label, value = row.get("label"), row.get("broad_type")
        if label is None or value not in VALUES:
            raise BroadTypeError(f"{Path(path).name}: '{value}' for '{label}' isn't a broad type.")
        if label in table:
            raise BroadTypeError(f"{Path(path).name}: '{label}' is listed twice.")
        table[label] = value
    return table


def celltypist_broad_types(model: str, folder: str | Path = TABLE_DIR) -> dict[str, str] | None:
    """A CellTypist model's label -> broad type, or None if the model has no table."""
    path = table_path(model, folder)
    return dict(_read(str(path))) if path.is_file() else None


def broad_type(label: str | None, table: dict[str, str] | None) -> str | None:
    """A label's broad type from a table; None when there is no label, no table, or the label isn't in it."""
    if label is None or table is None:
        return None
    return table.get(label)
