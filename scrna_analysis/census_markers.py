"""Cell-type marker tables built from CELLxGENE Census reference cells. 
Census doesn't publish marker genes: it holds millions of labeled cells. So a marker table is
built once per Census release, species and tissue, ahead of time (`scripts/build_census_markers.py`),
and the app reads the small saved table. Four steps:

1. `fetch_reference` (needs the internet): the tissue's healthy, primary cells; cell types with
   enough cells from enough studies; a fixed-seed sample of each; their raw counts.
2. `compute_marker_table`: each cell type's marker genes, from a Wilcoxon test of the type against
   the cell types it is **not related to** in the Cell Ontology. Census labels mix levels of detail
   ("B cell", "naive B cell"), and comparing a type with its own parent or children would cancel
   their shared genes. Vague labels ("unknown", "lymphocyte", …) are dropped first.
3. `save_marker_table`: the table (TSV) plus a provenance file (JSON: release, query, settings,
   cells and studies per type, ontology version, the table's SHA-256).
4. `load_marker_table`: back to cell type -> genes, with the provenance.

Only steps 2–4 run without the internet, so the app never waits for Census.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from scrna_analysis.params import CensusMarkerParams

CENSUS_DIR = Path(__file__).resolve().parent.parent / "data" / "markers" / "census"
ONTOLOGY_FILE = Path(__file__).resolve().parent.parent / "data" / "markers" / "cl-basic.obo"
CENSUS_SPECIES = {"human": "homo_sapiens", "mouse": "mus_musculus"}  # app name -> Census name


class CensusMarkerError(ValueError):
    """A marker table can't be built, read or used. The message says why."""


# --- The Cell Ontology ---------------------------------------------------------------------

def read_cell_ontology(path: str | Path = ONTOLOGY_FILE) -> dict[str, dict]:
    """The Cell Ontology's terms from an OBO file: id -> {"name", "parents"} (is_a only)."""
    path = Path(path)
    if not path.is_file():
        raise CensusMarkerError(f"The Cell Ontology file {path.name} isn't there. Run scripts/download_cell_ontology.py.")
    terms: dict[str, dict] = {}
    current: dict | None = None
    for line in path.read_text().splitlines():
        if line == "[Term]":
            current = {"id": None, "name": None, "parents": set()}
        elif line.startswith("[") and line.endswith("]"):
            current = None  # a [Typedef] or other stanza
        elif current is not None and line.startswith("id: "):
            current["id"] = line[4:].strip()
            terms[current["id"]] = current
        elif current is not None and line.startswith("name: "):
            current["name"] = line[6:].strip()
        elif current is not None and line.startswith("is_a: "):
            current["parents"].add(line[6:].split("!")[0].split("{")[0].strip())  # drop "{is_inferred=…}" tags
    return {term_id: {"name": t["name"], "parents": t["parents"]} for term_id, t in terms.items()}


def ontology_version(path: str | Path = ONTOLOGY_FILE) -> str:
    for line in Path(path).read_text().splitlines()[:30]:
        if line.startswith("data-version:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def ancestors(ontology: dict[str, dict], term_id: str) -> set[str]:
    """Every term above `term_id` (its parents, their parents, …), not itself."""
    seen: set[str] = set()
    stack = list(ontology.get(term_id, {}).get("parents", ()))
    while stack:
        parent = stack.pop()
        if parent not in seen:
            seen.add(parent)
            stack.extend(ontology.get(parent, {}).get("parents", ()))
    return seen


def related(ontology: dict[str, dict], a: str, b: str) -> bool:
    """True if one term is the other, or one is above the other ("B cell" and "naive B cell")."""
    return a == b or a in ancestors(ontology, b) or b in ancestors(ontology, a)


# --- Computing a table from reference cells ----------------------------------------------------

def choose_cell_types(cells: pd.DataFrame, params: CensusMarkerParams) -> pd.DataFrame:
    """From one row per reference cell (columns cell_type, cell_type_ontology_term_id, dataset_id),
    the cell types to use: not vague, at least `min_cells` cells from at least `min_studies` studies.
    One row per type: cell_type, ontology_id, n_cells, n_studies (most cells first)."""
    vague = {v.lower() for v in params.vague_labels}
    usable = cells[~cells["cell_type"].astype(str).str.lower().isin(vague)]
    per = (usable.groupby(["cell_type", "cell_type_ontology_term_id"], observed=True)
           .agg(n_cells=("dataset_id", "size"), n_studies=("dataset_id", "nunique")).reset_index()
           .rename(columns={"cell_type_ontology_term_id": "ontology_id"}))
    keep = per[(per["n_cells"] >= params.min_cells) & (per["n_studies"] >= params.min_studies)]
    return keep.sort_values(["n_cells", "cell_type"], ascending=[False, True]).reset_index(drop=True)


def sample_cells(cells: pd.DataFrame, types: pd.DataFrame, params: CensusMarkerParams) -> np.ndarray:
    """Up to `cells_per_type` cells of each chosen type, picked with a fixed seed (the same cells
    every time). Returns the `soma_joinid`s, sorted."""
    rng = np.random.default_rng(params.seed)
    chosen = []
    for cell_type in types["cell_type"]:
        ids = np.sort(cells.loc[cells["cell_type"] == cell_type, "soma_joinid"].to_numpy())
        chosen.append(rng.choice(ids, min(params.cells_per_type, len(ids)), replace=False))
    return np.sort(np.concatenate(chosen)) if chosen else np.array([], dtype=np.int64)


def excluded_genes(genes, params: CensusMarkerParams) -> np.ndarray:
    """True for genes never used as markers (`params.excluded_genes`: ribosomal and mitochondrial)."""
    return np.asarray(pd.Series(list(genes), dtype=str).str.fullmatch(params.excluded_genes, case=False))


def compute_marker_table(reference: ad.AnnData, ontology: dict[str, dict], params: CensusMarkerParams) -> pd.DataFrame:
    """Each cell type's marker genes from reference cells (raw counts in X; obs `cell_type` and
    `cell_type_ontology_term_id`; gene names in var_names). Counts are normalized per cell and
    log-transformed (Scanpy), then each type is compared with the cells of every type it isn't
    related to (Wilcoxon, Scanpy's `rank_genes_groups`). Ribosomal and mitochondrial genes are left
    out (`excluded_genes`). Kept: genes higher in the type, detected in
    at least `min_fraction` of its cells, best `top_n_genes` by score. One row per type and gene:
    cell_type, ontology_id, rank, gene, score, logfoldchange, pct_in, pct_out, n_cells, n_compared."""
    data = reference.copy()
    sc.pp.normalize_total(data)  # on all genes, so each cell's scaling isn't changed by the exclusion
    sc.pp.log1p(data)
    data = data[:, ~excluded_genes(data.var_names, params)].copy()
    data.obs["cell_type"] = data.obs["cell_type"].astype(str)
    types = (data.obs[["cell_type", "cell_type_ontology_term_id"]].drop_duplicates()
             .set_index("cell_type")["cell_type_ontology_term_id"].to_dict())
    rows = []
    for cell_type, term in sorted(types.items()):
        others = [t for t, other_term in types.items() if not related(ontology, term, other_term)]
        if not others:
            continue  # nothing unrelated to compare with
        subset = data[data.obs["cell_type"].isin([cell_type, *others])].copy()
        subset.obs["group"] = np.where(subset.obs["cell_type"] == cell_type, "type", "others")
        sc.tl.rank_genes_groups(subset, "group", groups=["type"], reference="others", method="wilcoxon",
                                use_raw=False)
        result = sc.get.rank_genes_groups_df(subset, group="type")
        # Share of cells detecting each gene, in the type and in the cells it was compared with.
        in_type = subset.obs["group"].to_numpy() == "type"
        detected = subset.X > 0
        pct_in = pd.Series(np.asarray(detected[in_type].mean(axis=0)).ravel(), index=subset.var_names)
        pct_out = pd.Series(np.asarray(detected[~in_type].mean(axis=0)).ravel(), index=subset.var_names)
        result["pct_nz_group"] = pct_in.reindex(result["names"]).to_numpy()
        result["pct_nz_reference"] = pct_out.reindex(result["names"]).to_numpy()
        result = result[(result["logfoldchanges"] > 0) & (result["pct_nz_group"] >= params.min_fraction)]
        top = result.sort_values("scores", ascending=False, kind="stable").head(params.top_n_genes)
        n_in = int((subset.obs["group"] == "type").sum())
        for rank, row in enumerate(top.itertuples(), start=1):
            rows.append({"cell_type": cell_type, "ontology_id": term, "rank": rank, "gene": row.names,
                         "score": float(row.scores), "logfoldchange": float(row.logfoldchanges),
                         "pct_in": float(row.pct_nz_group), "pct_out": float(row.pct_nz_reference),
                         "n_cells": n_in, "n_compared": int(subset.n_obs - n_in)})
    return pd.DataFrame(rows, columns=["cell_type", "ontology_id", "rank", "gene", "score", "logfoldchange",
                                       "pct_in", "pct_out", "n_cells", "n_compared"])


# --- Saving and loading -------------------------------------------------------------------------

def table_path(species: str, tissue: str, release: str, folder: str | Path = CENSUS_DIR) -> Path:
    """Where the table for a species and Census tissue is kept, e.g.
    data/markers/census/2025-11-08/human_blood.tsv.gz."""
    return Path(folder) / release / f"{species}_{tissue.replace(' ', '_')}.tsv.gz"


def save_marker_table(table: pd.DataFrame, provenance: dict, path: str | Path) -> dict:
    """Write the table (gzip TSV) and its provenance (same name, .json, with the table's SHA-256)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(path, "wb", mtime=0) as handle:  # mtime=0: the same table gives the same bytes
        handle.write(table.to_csv(sep="\t", index=False).encode())
    provenance = {**provenance, "table_file": path.name, "table_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "n_cell_types": int(table["cell_type"].nunique()) if len(table) else 0}
    path.with_name(path.name.replace(".tsv.gz", ".json")).write_text(json.dumps(provenance, indent=2))
    return provenance


def load_marker_table(path: str | Path) -> tuple[dict[str, tuple[str, ...]], dict]:
    """A saved table as cell type -> marker genes (best first), and its provenance. Checks the
    table against the SHA-256 in its provenance, so an edited table is noticed."""
    path = Path(path)
    meta_path = path.with_name(path.name.replace(".tsv.gz", ".json"))
    if not path.is_file() or not meta_path.is_file():
        raise CensusMarkerError(f"No Census marker table at {path.name}. Build it with scripts/build_census_markers.py.")
    provenance = json.loads(meta_path.read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest() != provenance.get("table_sha256"):
        raise CensusMarkerError(f"{path.name} doesn't match its recorded fingerprint; rebuild it.")
    table = pd.read_csv(path, sep="\t").sort_values(["cell_type", "rank"])
    sets = {cell_type: tuple(group["gene"].astype(str)) for cell_type, group in table.groupby("cell_type", sort=True)}
    return sets, provenance


# --- Census tissues ------------------------------------------------------------------------------

TISSUE_LIST = "tissues.json"  # per release folder: species -> Census tissues with their cell counts


def fetch_tissue_list(release: str = CensusMarkerParams.release, folder: str | Path = CENSUS_DIR) -> Path:
    """Save Census's tissues for each species with their number of cells (from Census's own
    summary table; needs the internet, ~2 s)."""
    import cellxgene_census  # only needed when building tables

    with cellxgene_census.open_soma(census_version=release) as census:
        summary = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()
    rows = summary[summary["category"] == "tissue_general"]
    tissues = {species: [{"tissue": str(r.label), "n_cells": int(r.unique_cell_count)}
                         for r in rows[rows["organism"] == census_name]
                         .sort_values("unique_cell_count", ascending=False).itertuples()]
               for species, census_name in CENSUS_SPECIES.items()}
    path = Path(folder) / release / TISSUE_LIST
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"census_version": release, "tissues": tissues}, indent=2))
    return path


def census_tissues(species: str, release: str = CensusMarkerParams.release,
                   folder: str | Path = CENSUS_DIR) -> list[dict]:
    """Census's tissues for a species (most cells first), each with `n_cells` and `has_table`
    (whether its marker table has been built). Empty if the species isn't in Census or the
    tissue list hasn't been saved (scripts/build_census_markers.py --tissues)."""
    path = Path(folder) / release / TISSUE_LIST
    if species not in CENSUS_SPECIES or not path.is_file():
        return []
    listed = json.loads(path.read_text())["tissues"].get(species, [])
    return [{**t, "has_table": table_path(species, t["tissue"], release, folder).is_file()} for t in listed]


def cell_type_terms(provenance: dict) -> dict[str, str]:
    """Cell type -> its Cell Ontology id, from a table's provenance."""
    return {t["cell_type"]: t["ontology_id"] for t in provenance.get("cell_types", [])}


# --- Fetching reference cells from Census (needs the internet) ----------------------------------

def fetch_reference(species: str, tissue: str, params: CensusMarkerParams, log=print) -> tuple[ad.AnnData, pd.DataFrame, dict]:
    """The reference cells for one species and Census tissue (`tissue_general`, e.g. "blood"):
    healthy (`disease == 'normal'`), primary (each cell once) cells; types chosen by
    `choose_cell_types`; a sample of each (`sample_cells`); raw counts and gene names. Returns the
    cells, the chosen types and the query (for the provenance)."""
    import cellxgene_census  # only needed when building tables

    if species not in CENSUS_SPECIES:
        raise CensusMarkerError(f"Census marker tables cover {' and '.join(CENSUS_SPECIES)} only.")
    query = (f"is_primary_data == True and disease == 'normal' and tissue_general == '{tissue}'")
    with cellxgene_census.open_soma(census_version=params.release) as census:
        obs = census["census_data"][CENSUS_SPECIES[species]].obs
        cells = obs.read(value_filter=query, column_names=["soma_joinid", "cell_type", "cell_type_ontology_term_id",
                                                           "dataset_id"]).concat().to_pandas()
        log(f"{len(cells):,} healthy {tissue} cells in Census {params.release}")
        types = choose_cell_types(cells, params)
        ids = sample_cells(cells, types, params)
        log(f"{len(types)} cell types pass the rule; fetching {len(ids):,} cells…")
        reference = cellxgene_census.get_anndata(
            census, organism={"human": "Homo sapiens", "mouse": "Mus musculus"}[species], obs_coords=ids,
            obs_column_names=["cell_type", "cell_type_ontology_term_id", "dataset_id"],
            var_column_names=["feature_id", "feature_name"])
    reference.var_names = reference.var["feature_name"].astype(str)
    reference.var_names_make_unique()
    return reference, types, {"value_filter": query, "census_version": params.release}


def build_marker_table(species: str, tissue: str, params: CensusMarkerParams = CensusMarkerParams(),
                       folder: str | Path = CENSUS_DIR, ontology_path: str | Path = ONTOLOGY_FILE, log=print) -> Path:
    """Fetch, compute and save one table; returns its path."""
    ontology = read_cell_ontology(ontology_path)
    reference, types, query = fetch_reference(species, tissue, params, log)
    log("Computing marker genes…")
    table = compute_marker_table(reference, ontology, params)
    path = table_path(species, tissue, params.release, folder)
    provenance = {
        "source": "CZ CELLxGENE Census", "census_version": params.release, "species": species, "tissue": tissue,
        "query": query, "settings": asdict(params), "cell_types": types.to_dict(orient="records"),
        "n_reference_cells": int(reference.n_obs), "cell_ontology_version": ontology_version(ontology_path),
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "citation": "CZ CELLxGENE Discover, Census (CZI Single-Cell Biology Program et al. 2023); CC BY 4.0",
    }
    save_marker_table(table, provenance, path)
    log(f"Saved {table['cell_type'].nunique()} cell types to {path}")
    return path
