"""Cell-type marker genes from a marker database (PanglaoDB), for the species and tissues chosen.

PanglaoDB (Franzén et al. 2019, Database; https://panglaodb.se) lists marker genes for the cell
types of 29 organs, for human and mouse. The file is downloaded by
`scripts/download_marker_db.py` into `data/markers/`; it isn't bundled with the app, because
the database states no licence.

Three steps:

1. `read_panglaodb`: read the file and check that it has the expected columns.
2. `select_marker_sets`: keep the markers for one species and the chosen tissues
   (`MarkerDatabaseParams`: canonical markers with sensitivity > 0.5 in that species).
3. `match_to_genes`: find each marker among the dataset's own gene names. PanglaoDB writes
   gene symbols in upper case (human style), and a dataset may write them otherwise (mouse:
   Cd3d), so a name is matched ignoring case when there is no exact match. Markers that aren't
   found are reported, never silently dropped.

Users can bring their own markers instead (`read_marker_list`: a CSV or TSV with a cell type
and a gene column); they are matched and tested the same way.

Nothing here knows about a particular dataset: the result depends only on the database file
(or the user's list), the settings and the gene names given.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

from scrna_analysis.params import MarkerDatabaseParams

# Species the database covers: the app's name -> PanglaoDB's code in its "species" column.
SPECIES = {"human": "Hs", "mouse": "Mm"}

# The species to suggest from how mitochondrial genes are named (params.MITO_PATTERNS). Only
# a suggestion: "MT-" is also used by other mammals, and "mt-" by zebrafish.
SPECIES_BY_MITO_PREFIX = {"MT-": "human", "mt-": "mouse"}

# Where scripts/download_marker_db.py saves the database.
# Snapshot: the project's data/ folder (three levels up from this file).
PANGLAODB_FILE = Path(__file__).resolve().parents[3] / "data" / "markers" / "PanglaoDB_markers_27_Mar_2020.tsv.gz"

# PanglaoDB column -> the name used here.
_COLUMNS = {
    "species": "species",
    "official gene symbol": "gene",
    "cell type": "cell_type",
    "organ": "organ",
    "canonical marker": "canonical",
    "sensitivity_human": "sensitivity_human",
    "sensitivity_mouse": "sensitivity_mouse",
}


class MarkerDatabaseError(ValueError):
    """The marker database can't be read or used with these settings. The message says why."""


@dataclass(frozen=True)
class MarkerDatabase:
    name: str  # "PanglaoDB"
    version: str  # the release date in the file name, e.g. "27 Mar 2020"
    file_name: str
    sha256: str
    table: pd.DataFrame  # one row per cell type per gene, with the columns in _COLUMNS


@dataclass(frozen=True)
class MarkerMatch:
    sets: dict[str, tuple[str, ...]]  # cell type -> markers as named in the dataset
    missing: dict[str, tuple[str, ...]]  # cell type -> markers not found in the dataset
    too_few: dict[str, int]  # cell types left out: markers found, below min_markers
    case_matched: int  # markers found only when ignoring upper/lower case


def read_panglaodb(path: str | Path) -> MarkerDatabase:
    """Read a PanglaoDB marker file (.tsv or .tsv.gz)."""
    path = Path(path)
    if not path.is_file():
        raise MarkerDatabaseError(
            f"The marker database file {path.name} isn't there. Download it with "
            "scripts/download_marker_db.py."
        )
    try:
        raw = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as error:  # an unreadable file, whatever the reason
        raise MarkerDatabaseError(f"{path.name} can't be read as a PanglaoDB marker file ({error}).") from error
    absent = [c for c in _COLUMNS if c not in raw.columns]
    if absent:
        raise MarkerDatabaseError(f"{path.name} isn't a PanglaoDB marker file: it has no {', '.join(absent)} column.")

    table = raw[list(_COLUMNS)].rename(columns=_COLUMNS)
    table["canonical"] = pd.to_numeric(table["canonical"], errors="coerce") == 1
    for column in ("sensitivity_human", "sensitivity_mouse"):
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table["species"] = table["species"].fillna("")
    found = re.search(r"(\d{1,2})_([A-Za-z]{3})_(\d{4})", path.name)
    return MarkerDatabase(
        name="PanglaoDB",
        version=" ".join(found.groups()) if found else "unknown version",
        file_name=path.name,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        table=table,
    )


@lru_cache(maxsize=4)
def _read_cached(path: str, modified: float) -> MarkerDatabase:
    return read_panglaodb(path)


def load_panglaodb(path: str | Path = PANGLAODB_FILE) -> MarkerDatabase:
    """`read_panglaodb`, read once and reused until the file changes. Treat it as read-only."""
    path = Path(path)
    if not path.is_file():
        return read_panglaodb(path)  # raises with the message saying how to get it
    return _read_cached(str(path), path.stat().st_mtime)


def suggest_species(mito_prefix: str | None) -> str | None:
    """The database species suggested by the mitochondrial gene prefix, if any."""
    return SPECIES_BY_MITO_PREFIX.get(mito_prefix or "")


def tissues(db: MarkerDatabase) -> tuple[str, ...]:
    """The database's tissue (organ) names, sorted."""
    return tuple(sorted(db.table["organ"].dropna().unique()))


def select_marker_sets(db: MarkerDatabase, params: MarkerDatabaseParams) -> dict[str, tuple[str, ...]]:
    """Cell type -> marker genes (database names) for the species and tissues in `params`.

    Markers are ordered by sensitivity (most widely expressed in the cell type first), then
    by name; cell types are in alphabetical order.
    """
    chosen = _selected_rows(db, params)  # checks the species first
    unknown = [t for t in params.tissues if t not in tissues(db)]
    if unknown:
        raise MarkerDatabaseError(f"{db.name} has no tissue named {', '.join(unknown)}.")
    if params.tissues:
        chosen = chosen[chosen["organ"].isin(params.tissues)]
    chosen = chosen.sort_values([f"sensitivity_{params.species}", "gene"], ascending=[False, True], kind="stable")
    return {
        cell_type: tuple(dict.fromkeys(group["gene"]))
        for cell_type, group in sorted(chosen.groupby("cell_type", sort=False), key=lambda item: item[0])
    }


def cell_types_by_tissue(db: MarkerDatabase, params: MarkerDatabaseParams) -> dict[str, tuple[str, ...]]:
    """Tissue -> the cell types (sorted) with markers that pass the selection rule for the
    species in `params` (its tissues are ignored). Every tissue is listed, even with none."""
    chosen = _selected_rows(db, params)
    return {tissue: tuple(sorted(chosen.loc[chosen["organ"] == tissue, "cell_type"].unique()))
            for tissue in tissues(db)}


def tissues_of_cell_types(db: MarkerDatabase, params: MarkerDatabaseParams) -> dict[str, tuple[str, ...]]:
    """Cell type -> the tissues (sorted) it is listed under, for the species in `params` and the
    markers that pass the selection rule (its tissues are ignored)."""
    chosen = _selected_rows(db, params)
    return {cell_type: tuple(sorted(group["organ"].dropna().unique()))
            for cell_type, group in chosen.groupby("cell_type")}


def _selected_rows(db: MarkerDatabase, params: MarkerDatabaseParams) -> pd.DataFrame:
    """The database rows for the species in `params` that pass the selection rule (any tissue)."""
    if params.species is None:
        raise MarkerDatabaseError("Choose a species before selecting marker genes.")
    if params.species not in SPECIES:
        raise MarkerDatabaseError(
            f"{db.name} has markers for {' and '.join(SPECIES)} only, not '{params.species}'."
        )
    t = db.table
    keep = t["species"].str.split().apply(lambda codes: SPECIES[params.species] in codes)
    if params.canonical_only:
        keep &= t["canonical"]
    keep &= t[f"sensitivity_{params.species}"] > params.min_sensitivity
    return t[keep]


@dataclass(frozen=True)
class MarkerList:
    """A user's own marker list, read and checked (`read_marker_list`)."""

    file_name: str
    sha256: str
    sets: dict[str, tuple[str, ...]]  # cell type -> marker genes, in the file's order
    notes: tuple[str, ...]  # what was tidied (duplicates, empty rows), in plain words


# Column names accepted for the two columns (compared ignoring case, spaces, "_" and "-").
CELL_TYPE_COLUMNS = ("celltype", "cell", "type")
GENE_COLUMNS = ("gene", "genes", "marker", "markers", "symbol")
MAX_CELL_TYPES = 500
MAX_LABEL_LENGTH = 60

# The example file offered on the page ("Download a template"): long format, one gene per row.
MARKER_LIST_TEMPLATE = (
    "cell_type,gene\n"
    "Hepatocytes,ALB\nHepatocytes,APOA1\nHepatocytes,TTR\n"
    "Endothelial cells,PECAM1\nEndothelial cells,VWF\nEndothelial cells,CDH5\n"
    "Fibroblasts,COL1A1\nFibroblasts,DCN\nFibroblasts,LUM\n"
)


def read_marker_list(path: str | Path, file_name: str | None = None) -> MarkerList:
    """Read a user's marker list: a CSV or TSV file with a cell type column and a gene column,
    one gene per row (a cell may also hold several genes separated by commas or semicolons).

    Tidies and reports: surrounding spaces, empty rows, repeated genes. Raises
    `MarkerDatabaseError` with a plain-language message when the file can't be used.
    """
    path = Path(path)
    name = file_name or path.name
    if Path(name).suffix.lower() not in (".csv", ".tsv", ".txt"):
        raise MarkerDatabaseError(f"'{name}' isn't a CSV or TSV file. Save the marker list as .csv or .tsv.")
    try:
        raw = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    except Exception as error:  # any unreadable file
        raise MarkerDatabaseError(f"'{name}' can't be read as a table ({error}).") from error

    def simple(column: str) -> str:
        return re.sub(r"[\s_\-]+", "", str(column)).lower()

    columns = {simple(c): c for c in raw.columns}
    cell_col = next((columns[c] for c in CELL_TYPE_COLUMNS if c in columns), None)
    gene_col = next((columns[c] for c in GENE_COLUMNS if c in columns), None)
    if cell_col is None or gene_col is None:
        raise MarkerDatabaseError(
            f"'{name}' needs two columns, 'cell_type' and 'gene' (it has: {', '.join(map(str, raw.columns))}). "
            "Download the template to see the layout."
        )

    sets: dict[str, list[str]] = {}
    empty = repeated = 0
    for cell_type, genes in zip(raw[cell_col], raw[gene_col]):
        cell_type = " ".join(str(cell_type).split())
        found = [g.strip() for g in re.split(r"[,;]", str(genes)) if g.strip()]
        if not cell_type or not found:
            empty += 1
            continue
        if len(cell_type) > MAX_LABEL_LENGTH:
            raise MarkerDatabaseError(f"The cell type name '{cell_type[:40]}…' is too long; keep names under "
                                      f"{MAX_LABEL_LENGTH} characters.")
        listed = sets.setdefault(cell_type, [])
        for gene in found:
            if gene in listed:
                repeated += 1
            else:
                listed.append(gene)
    if not sets:
        raise MarkerDatabaseError(f"'{name}' has no rows with both a cell type and a gene.")
    if len(sets) > MAX_CELL_TYPES:
        raise MarkerDatabaseError(f"'{name}' has {len(sets):,} cell types; the limit is {MAX_CELL_TYPES}.")

    notes = []
    if empty:
        notes.append(f"{empty:,} {'row was' if empty == 1 else 'rows were'} left out (no cell type or no gene).")
    if repeated:
        notes.append(f"{repeated:,} repeated {'gene was' if repeated == 1 else 'genes were'} listed once.")
    return MarkerList(file_name=name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      sets={ct: tuple(genes) for ct, genes in sets.items()}, notes=tuple(notes))


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
