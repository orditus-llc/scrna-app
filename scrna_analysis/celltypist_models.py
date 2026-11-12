"""CellTypist: a second, independent opinion on each cluster's cell type (Phase 0 of the supervisor's
plan, docs/Census_reference_build_plan.md; analyst choices 2026-11-12).

CellTypist (Domínguez Conde et al. 2022, Science; Xu et al. 2023, Cell) is a classifier (logistic
regression) trained on labeled reference cells. It looks at every gene of every cell and predicts a
cell type per cell; each cluster then takes the majority vote (CellTypist's own `majority_voting`
with the app's clusters as the over-clustering, so no clustering is added).

- **Models:** ready-made ones published by CellTypist, downloaded by
  `scripts/download_celltypist_models.py` into `data/markers/celltypist/` with their version and
  SHA-256 (`models.json`). Never bundled: the models state no licence.
- **Its own normalization:** CellTypist's models expect each cell normalized to 10,000 counts and
  log1p-transformed, so it gets its own copy made from the raw counts (the app normalizes to the
  median by default).
- **Suggested model:** from the species and Census tissue (`SUGGESTED_MODELS`, the app's own
  mapping, started small on purpose; analyst 2026-11-12). The user can choose another or none.
- **Step A:** the answer is shown next to the marker check. Nothing is combined or overruled yet:
  rolling both answers up to broad types and the agree / "Unknown" rule wait for the supervisor.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

MODEL_DIR = Path(__file__).resolve().parent.parent / "data" / "markers" / "celltypist"
LISTING = "models.json"  # written by the download script: version, date, source, SHA-256 per model
TARGET_SUM = 1e4  # CellTypist's models expect counts per cell normalized to 10,000, then log1p
LOW_GENE_MATCH = 0.5  # the app's own rule: warn when fewer than half of the model's genes are in the data
PREDICTION_KEY = "celltypist_majority"  # obs: each cell's cluster's CellTypist answer
CELL_KEY = "celltypist_cell"  # obs: CellTypist's answer for the cell itself
CITATION = ("CellTypist: Domínguez Conde et al. 2022, Science 376, eabl5197; "
            "Xu et al. 2023, Cell 186, 5876–5891")

# (species, Census tissue) -> the suggested model. The app's own mapping: only tissues a ready-made
# model clearly covers; everything else gets no suggestion (the marker check alone).
SUGGESTED_MODELS: dict[tuple[str, str], str] = {
    ("human", "blood"): "Immune_All_Low.pkl",
    ("human", "bone marrow"): "Immune_All_Low.pkl",
    ("human", "spleen"): "Immune_All_Low.pkl",
    ("human", "lymph node"): "Immune_All_Low.pkl",
    ("human", "immune system"): "Immune_All_Low.pkl",
    ("mouse", "brain"): "Mouse_Whole_Brain.pkl",  # most mouse brain data is adult (analyst's choice)
}
# Other models offered for a species (e.g. embryonic brain), besides its suggestions.
ALTERNATIVES: dict[str, tuple[str, ...]] = {"human": (), "mouse": ("Developing_Mouse_Brain.pkl",)}


class CellTypistError(ValueError):
    """CellTypist can't run with these settings. The message says why."""


def suggest_model(species: str | None, tissue: str | None) -> str | None:
    """The suggested model for a species and Census tissue, or None."""
    return SUGGESTED_MODELS.get(((species or "").lower(), (tissue or "").lower()))


def model_choices(species: str | None, folder: str | Path = MODEL_DIR) -> list[dict]:
    """The models offered for a species (suggestions first), each with its description from the
    download and whether it is downloaded."""
    species = (species or "").lower()
    names = list(dict.fromkeys([m for (sp, _), m in SUGGESTED_MODELS.items() if sp == species]
                               + list(ALTERNATIVES.get(species, ()))))
    listing = read_listing(folder)
    return [{"model": n, "details": (listing.get(n) or {}).get("details"), "version": (listing.get(n) or {}).get("version"),
             "downloaded": (Path(folder) / n).is_file()} for n in names]


def read_listing(folder: str | Path = MODEL_DIR) -> dict:
    path = Path(folder) / LISTING
    return json.loads(path.read_text()) if path.is_file() else {}


@lru_cache(maxsize=4)
def _load(path: str):
    from celltypist import models
    return models.Model.load(path)


def load_model(name: str, folder: str | Path = MODEL_DIR):
    """The model and its record (name, version, date, SHA-256 checked against the download)."""
    path = Path(folder) / name
    if not path.is_file():
        raise CellTypistError(f"The CellTypist model {name} isn't downloaded. Run scripts/download_celltypist_models.py "
                              "(or choose another model).")
    entry = read_listing(folder).get(name, {})
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if entry.get("sha256") and entry["sha256"] != sha:
        raise CellTypistError(f"{name} doesn't match its recorded fingerprint; download it again.")
    record = {"name": "CellTypist", "model": name, "version": entry.get("version"), "date": entry.get("date"),
              "source": entry.get("source"), "sha256": sha, "citation": CITATION}
    return _load(str(path)), record


@dataclass(frozen=True)
class CellTypistSummary:
    model: dict  # name, model, version, date, source, sha256, citation
    n_model_genes: int
    n_genes_matched: int
    clusters: pd.DataFrame  # cluster, n_cells, label, share, mean_probability, runner_up, runner_up_share
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def gene_match(self) -> float:
        return self.n_genes_matched / self.n_model_genes if self.n_model_genes else 0.0


def normalized_for_celltypist(counts: ad.AnnData) -> ad.AnnData:
    """A copy of raw counts normalized to 10,000 per cell and log1p-transformed (what the models expect)."""
    x = ad.AnnData(X=counts.X.copy(), obs=pd.DataFrame(index=counts.obs_names), var=pd.DataFrame(index=counts.var_names))
    sc.pp.normalize_total(x, target_sum=TARGET_SUM)
    sc.pp.log1p(x)
    return x


def predict(counts: ad.AnnData, clusters: pd.Series, model_name: str,
            folder: str | Path = MODEL_DIR) -> tuple[pd.DataFrame, CellTypistSummary]:
    """CellTypist's answer per cell and per cluster (majority vote over the given clusters).

    `counts`: raw counts of the clustered cells (all genes kept after QC); `clusters`: each cell's
    cluster, same cells in the same order. Returns a per-cell table (cell label, the cluster's
    majority label) and the per-cluster summary: the majority label, its share of the cluster's cells,
    the mean probability CellTypist gave that label in the cluster, and the runner-up."""
    import logging

    import celltypist

    logging.getLogger("celltypist").setLevel(logging.WARNING)  # its progress messages aren't needed
    if list(counts.obs_names) != list(clusters.index.astype(str)):
        raise CellTypistError("The raw counts and the clusters don't line up; re-run the steps from Filtering.")
    model, record = load_model(model_name, folder)
    data = normalized_for_celltypist(counts)
    features = np.asarray(model.features)
    matched = int(np.isin(features, data.var_names).sum())
    if matched == 0:
        raise CellTypistError(f"None of the genes of the CellTypist model {model_name} are in this dataset (check the "
                              "species and how genes are named), so it can't be used.")
    groups = clusters.astype(str).to_numpy()
    out = celltypist.annotate(data, model=model, majority_voting=True, over_clustering=groups)
    labels = out.predicted_labels
    probabilities = out.probability_matrix
    per_cell = pd.DataFrame({CELL_KEY: labels["predicted_labels"].astype(str).to_numpy(),
                             PREDICTION_KEY: labels["majority_voting"].astype(str).to_numpy()}, index=counts.obs_names)

    rows = []
    order = (list(clusters.cat.categories) if isinstance(clusters.dtype, pd.CategoricalDtype)
             else list(pd.unique(groups)))
    for cluster in [str(c) for c in order if (groups == str(c)).any()]:
        in_cluster = groups == cluster
        shares = per_cell.loc[in_cluster, CELL_KEY].value_counts(normalize=True)
        label = per_cell.loc[in_cluster, PREDICTION_KEY].iloc[0]
        runner = [lab for lab in shares.index if lab != label]
        rows.append({"cluster": str(cluster), "n_cells": int(in_cluster.sum()), "label": label,
                     "share": float(shares.get(label, 0.0)),
                     "mean_probability": float(probabilities.loc[in_cluster, label].mean()),
                     "runner_up": runner[0] if runner else None,
                     "runner_up_share": float(shares[runner[0]]) if runner else None})
    notes = []
    share = matched / len(features)
    if share < LOW_GENE_MATCH:
        notes.append(f"Only {matched:,} of the model's {len(features):,} genes ({100 * share:.0f}%) are in this dataset, "
                     "so CellTypist's answers are less reliable. Check the species and how genes are named.")
    summary = CellTypistSummary(model=record, n_model_genes=len(features), n_genes_matched=matched,
                                clusters=pd.DataFrame(rows), notes=tuple(notes))
    return per_cell, summary
