"""The analyzed data as one AnnData file (.h5ad) that other tools can open.

One object with everything computed so far, laid out for reuse (Scanpy, CELLxGENE, Galaxy):

- cells: those kept after QC; genes: all genes kept after QC (not only the highly variable ones)
- `X`: log-normalized values (what Scanpy kept in `adata.raw`); `layers["counts"]`: raw counts
  of the same cells and genes (the Filter checkpoint)
- `var["highly_variable"]` and the dispersion columns; `obs`: QC metrics, the clusters, the
  predicted cell types, and the analyst's reviewed labels in their own column
- `obsm` (PCA, UMAP), `obsp` (neighbor graph) and Scanpy's own `uns` entries (marker genes,
  PCA variance, settings of each tool)
- `uns["orditus"]`: the full record (settings, seeds, versions, history, review) plus what the
  export contains, as JSON text

Steps that aren't up to date are left out. Nothing is computed here: every value comes from the
run's checkpoints, so the file holds exactly what the app showed.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import anndata as ad
import numpy as np

from scrna_analysis import review

if TYPE_CHECKING:  # only for type hints (pipeline imports nothing from here)
    from scrna_analysis.pipeline import AnalysisRun

# The later stages in order; the export is built from the latest one that is up to date.
EXPORT_STAGES = ("preprocess", "reduce_dimensions", "cluster", "markers", "annotate")
# Scanpy's own records kept in `uns` (whichever exist).
SCANPY_UNS_KEYS = ("log1p", "hvg", "pca", "neighbors", "umap", "leiden", "rank_genes_groups",
                   "cell_type_prediction")


class ExportError(ValueError):
    """The analysis can't be exported yet. The message says why."""


def analyzed_anndata(run: AnalysisRun) -> ad.AnnData:
    """The analysis so far as one AnnData (module docstring). Needs Preprocessing to be up to
    date (the file holds normalized values). Doesn't change the run."""
    status = run.status()
    if status["preprocess"] != "done":
        raise ExportError("Run Preprocessing first: the file holds the normalized data.")
    included = []
    for stage in EXPORT_STAGES:  # up to date, without gaps
        if status[stage] != "done":
            break
        included.append(stage)
    base = run.result(included[-1])
    counts = run.result("filter")
    if list(base.raw.var_names) != list(counts.var_names) or list(base.obs_names) != list(counts.obs_names):
        raise ExportError("The normalized data and the raw counts don't line up; re-run the steps from Filtering.")

    out = ad.AnnData(X=base.raw.X.copy(), obs=base.obs.copy(), var=base.raw.var.copy())
    out.layers["counts"] = counts.X.copy()
    for key, value in base.obsm.items():
        out.obsm[key] = np.asarray(value).copy()
    for key, value in base.obsp.items():
        out.obsp[key] = value.copy()
    for key in SCANPY_UNS_KEYS:
        if key in base.uns:
            out.uns[key] = base.uns[key]

    reviewed = run.review if "cluster" in included else {}
    if reviewed:
        out.obs[review.REVIEWED_KEY] = run.reviewed_labels()

    out.uns["orditus"] = json.dumps({
        **run.record(),
        "export": {
            "exported_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "steps_included": included,
            "X": "log-normalized values, all genes kept after QC",
            "layers": {"counts": "raw counts of the same cells and genes"},
            "obs_columns": {"clusters": "leiden" if "cluster" in included else None,
                            "predicted_cell_types": "cell_type_predicted" if "annotate" in included else None,
                            "reviewed_labels": review.REVIEWED_KEY if reviewed else None},
        },
    })
    return out
