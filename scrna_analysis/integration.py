"""Batch correction for several samples, and how well the samples mix (Milestone 6b, 2026-11-10).

Two parts:

1. `run_harmony`: Harmony (Korsunsky et al. 2019, Nature Methods; the `harmonypy` package) adjusts
   the PCA coordinates so that cells of the same type line up across batches. It is called directly
   (not through Scanpy's `harmony_integrate`, whose Scanpy 1.11.5 version expects the older
   harmonypy 0.0.x output layout). The corrected coordinates go to `obsm["X_pca_harmony"]`; gene
   expression is not changed, so marker genes still use the measured data.
2. `sample_mixing`: for each cell, how varied the samples are among its nearest neighbors (the
   Shannon entropy of their sample labels), divided by the most that is possible: the entropy of
   the samples' overall shares, or log(k) if that is smaller. 1 = as mixed as the samples' sizes
   allow; 0 = its neighbors all come from one sample. A simplified form of the LISI score used in
   the Harmony paper.

Everything is computed from the data and the given settings; nothing here is dataset-specific.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from scrna_analysis.params import LOW_MIXING, MIXING_NEIGHBORS, IntegrationParams

CORRECTED_KEY = "X_pca_harmony"  # obsm: Harmony's corrected PCA coordinates
MIXING_KEY = "sample_mixing"  # obs: each cell's mixing score in the coordinates used for clustering
MIXING_BEFORE_KEY = "sample_mixing_uncorrected"  # obs: the same on the uncorrected PCA (with correction)


class IntegrationError(ValueError):
    """Batch correction can't run with these settings. The message says why."""


@dataclass(frozen=True)
class HarmonySummary:
    batch_key: str
    n_batches: int
    n_components: int
    theta: float
    max_iter: int
    random_state: int
    warnings: tuple[str, ...] = field(default_factory=tuple)


def batch_column_choices(adata: ad.AnnData, max_values: int = 100) -> list[str]:
    """obs columns that could hold batches: text or category, 2 to `max_values` values, none missing."""
    out = []
    for column in adata.obs.columns:
        values = adata.obs[column]
        if not (pd.api.types.is_object_dtype(values) or isinstance(values.dtype, pd.CategoricalDtype)
                or pd.api.types.is_string_dtype(values)):
            continue
        if values.isna().any():
            continue
        if 2 <= values.astype(str).nunique() <= max_values:
            out.append(str(column))
    return out


def check_batch_key(adata: ad.AnnData, batch_key: str) -> pd.Series:
    """The batch labels, checked: the column exists, has no missing values and at least two batches."""
    if batch_key not in adata.obs:
        raise IntegrationError(f"There is no column named \"{batch_key}\" to use as the batch.")
    values = adata.obs[batch_key]
    if values.isna().any():
        raise IntegrationError(f"Some cells have no value in \"{batch_key}\", so it can't be used as the batch.")
    labels = values.astype(str)
    if labels.nunique() < 2:
        raise IntegrationError(f"\"{batch_key}\" has only one value; batch correction needs at least two batches.")
    return labels


def condition_warning(adata: ad.AnnData, batch_key: str, condition_key: str = "condition") -> str | None:
    """A warning when batches and conditions can't be told apart (the batch is the condition, or each
    batch value is exactly one condition): correcting the batches would then also remove the
    conditions' real difference from the map. Gene counts are never changed."""
    if batch_key == condition_key:
        return ("The batch is the condition itself, so correction will also remove the difference between "
                "conditions on the map (the gene counts are unchanged).")
    if condition_key not in adata.obs:
        return None
    pairs = adata.obs[[batch_key, condition_key]].astype(str).drop_duplicates()
    if (pairs[batch_key].is_unique and pairs[condition_key].nunique() > 1
            and pairs[batch_key].nunique() == pairs[condition_key].nunique()):
        return (f"Each value of \"{batch_key}\" is exactly one condition, so batches and conditions can't be told "
                "apart: correction will also remove the difference between conditions on the map (the gene "
                "counts are unchanged).")
    return None


def run_harmony(adata: ad.AnnData, params: IntegrationParams, n_components: int,
                use_rep: str = "X_pca") -> tuple[ad.AnnData, HarmonySummary]:
    """Harmony on the first `n_components` columns of `obsm[use_rep]` (the components the neighbor
    graph uses). Returns a copy with `obsm["X_pca_harmony"]`."""
    import harmonypy

    if params.method != "harmony":
        raise IntegrationError("Batch correction is switched off in these settings.")
    labels = check_batch_key(adata, params.batch_key)
    if use_rep not in adata.obsm:
        raise IntegrationError("Run PCA before batch correction.")
    coords = np.asarray(adata.obsm[use_rep][:, :n_components], dtype=np.float64)
    logging.getLogger("harmonypy").setLevel(logging.WARNING)  # its progress messages aren't needed
    result = harmonypy.run_harmony(coords, pd.DataFrame({"batch": labels.to_numpy()}), ["batch"],
                                   theta=params.theta, max_iter_harmony=params.max_iter,
                                   random_state=params.random_state, verbose=False)
    corrected = np.asarray(result.Z_corr)
    if corrected.shape != coords.shape:  # older harmonypy returns components x cells
        corrected = corrected.T
    out = adata.copy()
    out.obsm[CORRECTED_KEY] = corrected
    warning = condition_warning(adata, params.batch_key)
    summary = HarmonySummary(batch_key=params.batch_key, n_batches=int(labels.nunique()),
                             n_components=int(n_components), theta=params.theta, max_iter=params.max_iter,
                             random_state=params.random_state, warnings=(warning,) if warning else ())
    return out, summary


def sample_mixing(coords: np.ndarray, labels, n_neighbors: int = MIXING_NEIGHBORS) -> np.ndarray:
    """Each cell's mixing score (module docstring): entropy of the sample labels among its
    `n_neighbors` nearest neighbors (itself excluded), divided by the highest entropy possible."""
    labels = pd.Series(np.asarray(labels)).astype(str)
    codes, uniques = pd.factorize(labels)
    n = len(codes)
    if len(uniques) < 2 or n < 3:
        return np.ones(n)
    k = min(n_neighbors, n - 1)
    neighbors = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(coords, return_distance=False)[:, 1:]
    counts = np.zeros((n, len(uniques)))
    np.add.at(counts, (np.repeat(np.arange(n), k), codes[neighbors].ravel()), 1)
    shares = counts / k
    with np.errstate(divide="ignore", invalid="ignore"):
        local = -np.nansum(np.where(shares > 0, shares * np.log(shares), 0.0), axis=1)
    overall = np.bincount(codes) / n
    best = min(float(-(overall * np.log(overall)).sum()), float(np.log(k)))
    return np.clip(local / best, 0.0, 1.0)


def mixing_by_group(adata: ad.AnnData, groupby: str, batch_key: str, key: str = MIXING_KEY,
                    low: float = LOW_MIXING) -> pd.DataFrame:
    """Per group (e.g. cluster): its cells, mean mixing score, the sample with the most of its cells
    and that sample's share, and `low` (mean score below the threshold)."""
    obs = adata.obs[[groupby, batch_key, key]].copy()
    obs[groupby] = obs[groupby].astype(str)
    obs[batch_key] = obs[batch_key].astype(str)
    rows = []
    for group, cells in obs.groupby(groupby, sort=False, observed=True):
        shares = cells[batch_key].value_counts(normalize=True)
        rows.append({"group": group, "n_cells": int(len(cells)), "mixing": float(cells[key].mean()),
                     "top_sample": str(shares.index[0]), "top_share": float(shares.iloc[0]),
                     "n_samples": int((cells[batch_key].value_counts() > 0).sum()),
                     "low": bool(cells[key].mean() < low)})
    return pd.DataFrame(rows)


def composition(adata: ad.AnnData, sample_key: str, group_key: str) -> pd.DataFrame:
    """How each sample's cells divide among groups (clusters or cell types): one row per sample and
    group that has cells, with `n_cells` and `share` (of that sample's cells, 0–1). Descriptive only:
    with one sample per condition, differences can't be tested (supervisor's plan, 2026-11-09)."""
    for key in (sample_key, group_key):
        if key not in adata.obs:
            raise IntegrationError(f"There is no column named \"{key}\".")
    obs = adata.obs[[sample_key, group_key]].astype(str)
    counts = obs.groupby([sample_key, group_key], observed=True).size().rename("n_cells").reset_index()
    counts["share"] = counts["n_cells"] / counts.groupby(sample_key)["n_cells"].transform("sum")
    return counts.rename(columns={sample_key: "sample", group_key: "group"})
