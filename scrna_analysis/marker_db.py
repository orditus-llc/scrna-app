"""Marker genes matched to a dataset's genes, and clues to the dataset's species.

Used by Cell Type Annotation with the CELLxGENE Census marker tables (`census_markers.py`):

- `match_to_genes`: find each cell type's marker genes among the dataset's own gene names. A
  reference may write gene symbols one way (human style, CD3D) and a dataset another (mouse
  style, Cd3d), so a name is matched ignoring case when there is no exact match. Markers that
  aren't found are reported, never silently dropped.
- `species_clues` / `suggest_species`: what the gene names say about the species (Data Overview
  suggests the species from them and warns when the user's choice disagrees).

PanglaoDB, the marker database used here until 2026-11-12, was removed (analyst's decision); the
frozen app copies in `archive/` keep it.

Nothing here knows about a particular dataset: the result depends only on the marker lists,
the settings and the gene names given.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping

import pandas as pd

# Species with cell-type references (Census marker tables, CellTypist models). Any other species
# runs every step except Cell Type Annotation.
SPECIES = ("human", "mouse")

# The species to suggest from how mitochondrial genes are named (params.MITO_PATTERNS). Only
# a suggestion: "MT-" is also used by other mammals, and "mt-" by zebrafish.
SPECIES_BY_MITO_PREFIX = {"MT-": "human", "mt-": "mouse"}


class MarkerDatabaseError(ValueError):
    """Marker genes can't be matched with these settings. The message says why."""


@dataclass(frozen=True)
class MarkerMatch:
    sets: dict[str, tuple[str, ...]]  # cell type -> markers as named in the dataset
    missing: dict[str, tuple[str, ...]]  # cell type -> markers not found in the dataset
    too_few: dict[str, int]  # cell types left out: markers found, below min_markers
    case_matched: int  # markers found only when ignoring upper/lower case


def suggest_species(mito_prefix: str | None) -> str | None:
    """The species suggested by the mitochondrial gene prefix, if any."""
    return SPECIES_BY_MITO_PREFIX.get(mito_prefix or "")


# Ensembl gene ID prefixes -> species (Ensembl's naming; ENSG is human).
ENSEMBL_PREFIXES = {"ENSMUSG": "mouse", "ENSRNOG": "rat", "ENSDARG": "zebrafish", "ENSG": "human"}
# How gene symbols are written: HGNC (human) uses capitals (CD3D); MGI (mouse, also rat) capitalizes
# only the first letter (Cd3d).
SPECIES_BY_NAME_STYLE = {"capitals": "human", "first_capital": "mouse"}
CLEAR_SHARE = 0.8  # the app's own rule: a style or ID prefix counts when it covers this share of genes


def species_clues(gene_names: Iterable[str], gene_ids: Iterable[str] | None = None,
                  mito_prefix: str | None = None) -> list[dict[str, str]]:
    """What the gene names say about the species: Ensembl ID prefixes (if the data has IDs), the
    way symbols are written, and the mitochondrial gene prefix. Each clue: {"species", "clue"}
    (plain words). Naming only: data renamed to another species' genes would pass unnoticed."""
    names = pd.Series(list(gene_names), dtype=str)
    clues: list[dict[str, str]] = []

    ids = pd.Series(list(gene_ids), dtype=str) if gene_ids is not None else names
    ids = ids[ids.str.startswith("ENS")]
    if len(ids) and len(ids) >= 0.5 * len(names):
        for prefix, species in ENSEMBL_PREFIXES.items():
            if ids.str.startswith(prefix).mean() >= CLEAR_SHARE:
                clues.append({"species": species, "clue": f"Ensembl IDs start with {prefix}"})
                break

    symbols = names[~names.str.startswith("ENS") & names.str.contains(r"[A-Za-z]")]
    if len(symbols):
        capitals = symbols.str.fullmatch(r"[A-Z0-9][A-Z0-9.\-]*").mean()
        first_capital = symbols.str.fullmatch(r"[A-Z][a-z0-9.\-]*[a-z][a-z0-9.\-]*").mean()
        if capitals >= CLEAR_SHARE:
            clues.append({"species": SPECIES_BY_NAME_STYLE["capitals"], "clue": "gene names are written in capitals, like CD3D"})
        elif first_capital >= CLEAR_SHARE:
            clues.append({"species": SPECIES_BY_NAME_STYLE["first_capital"], "clue": "gene names are written like Cd3d"})

    if suggest_species(mito_prefix):
        clues.append({"species": suggest_species(mito_prefix), "clue": f"mitochondrial genes are named {mito_prefix}…"})
    return clues


def match_to_genes(
    marker_sets: Mapping[str, Iterable[str]], gene_names: Iterable[str], min_markers: int
) -> MarkerMatch:
    """Find each marker among `gene_names` (exact name first, otherwise ignoring case).

    A name that matches several genes when case is ignored (and none exactly) counts as not
    found, since the right one can't be told apart. Cell types with fewer than `min_markers`
    markers found are left out (listed in `too_few`).
    """
    if min_markers < 1:
        raise MarkerDatabaseError("The minimum number of marker genes must be at least 1.")
    names = list(gene_names)
    exact = set(names)
    by_lower: dict[str, list[str]] = {}
    for name in names:
        by_lower.setdefault(name.lower(), []).append(name)

    sets, missing, too_few, case_matched = {}, {}, {}, set()
    for cell_type, markers in marker_sets.items():
        found, absent = [], []
        for marker in markers:
            if marker in exact:
                gene = marker
            elif len(by_lower.get(marker.lower(), [])) == 1:
                gene = by_lower[marker.lower()][0]
                case_matched.add(marker)
            else:
                absent.append(marker)
                continue
            if gene not in found:
                found.append(gene)
        if absent:
            missing[cell_type] = tuple(absent)
        if len(found) >= min_markers:
            sets[cell_type] = tuple(found)
        else:
            too_few[cell_type] = len(found)
    return MarkerMatch(sets=sets, missing=missing, too_few=too_few, case_matched=len(case_matched))
