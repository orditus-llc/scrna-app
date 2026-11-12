"""The data the React interface receives: pipeline state -> plain JSON.

This is the "data contract" between the Shiny server and the browser. Each function takes
an `AnalysisRun` and returns JSON-serializable dicts and lists only. All numbers come from
the run's results and settings; nothing here is PBMC3k-specific. No Shiny code.
"""

from __future__ import annotations

import json
from dataclasses import asdict, replace
from typing import Any

import numpy as np

from pbmc_analysis.params import (
    REFERENCE_PRESET,
    AnalysisParams,
    AnnotationParams,
    FilterParams,
    HVGParams,
    LeidenParams,
    MarkerParams,
    NeighborsParams,
    NormalizeParams,
    PCAParams,
    RegressOutParams,
    ScaleParams,
    UMAPParams,
)
from pbmc_analysis import annotation as annotation_tools
from pbmc_analysis import markers as marker_tools
from pbmc_analysis import qc as qc_tools
from pbmc_analysis.pipeline import AnalysisRun, package_versions

# The 7 steps shown on the Workflow page, and the pipeline stages behind each.
# "Quality Control" covers two stages: QC metrics and filtering.
STEPS = (
    ("data_overview", "Data Overview", ("load",)),
    ("quality_control", "Quality Control", ("qc_metrics", "filter")),
    ("preprocessing", "Preprocessing", ("preprocess",)),
    ("dimensionality_reduction", "Dimensionality Reduction", ("reduce_dimensions",)),
    ("clustering", "Clustering", ("cluster",)),
    ("marker_genes", "Marker Genes", ("markers",)),
    ("cell_type_annotation", "Cell Type Annotation", ("annotate",)),
)

# What each step does, shown until it has results.
DESCRIPTIONS = {
    "data_overview": "Load the data and check it holds raw counts",
    "quality_control": "Measure every cell and remove poor-quality ones",
    "preprocessing": "Normalize and pick the most informative genes",
    "dimensionality_reduction": "Summarize cells with PCA, a neighbor graph and UMAP",
    "clustering": "Group similar cells into clusters",
    "marker_genes": "Find the genes that mark each cluster",
    "cell_type_annotation": "Predict a cell type for each cluster",
}

# Pipeline stage -> label used in the run log (the step names; Quality Control has two stages).
STAGE_LABELS = {
    "load": "Data Overview", "qc_metrics": "QC metrics", "filter": "Filtering", "preprocess": "Preprocessing",
    "reduce_dimensions": "Dimensionality Reduction", "cluster": "Clustering", "markers": "Marker Genes",
    "annotate": "Cell Type Annotation",
}
# The marker tests offered on the Marker Genes page: (Scanpy name, label shown).
MARKER_METHODS = [("wilcoxon", "Wilcoxon rank-sum"), ("t-test", "t-test"),
                  ("t-test_overestim_var", "t-test (overestimated variance)"), ("logreg", "Logistic regression")]
TEST_NAMES = dict(MARKER_METHODS)


def workflow_payload(run: AnalysisRun, notice: dict[str, Any] | None = None,
                     job: dict[str, Any] | None = None) -> dict[str, Any]:
    """Everything the Workflow page shows. `notice` is the latest message for the user;
    `job` is the latest "Run all remaining steps" job (see `run_remaining_payload`)."""
    load = run.stage_result("load").summary if run.status()["load"] == "done" else None
    steps = [_step(run, number, key, title, stages)
             for number, (key, title, stages) in enumerate(STEPS, start=1)]
    return {
        "dataset": {
            "file": run.data_path.name,
            "n_cells": load["n_cells"] if load else None,
            "n_genes": load["n_genes"] if load else None,
            "fingerprint": load["input_sha256"][:12] if load else None,
        },
        "reference_settings": run.params == REFERENCE_PRESET,
        "steps": steps,
        "next_step": next_step(steps),
        "settings_record": settings_record(run),
        "versions": {name: package_versions()[name] for name in ("scanpy", "leidenalg", "umap-learn")},
        "run_log": run_log(run),
        "notice": notice,
        "remaining_steps": sum(step["status"] != "done" for step in steps),
        "run_remaining": run_remaining_payload(job),
    }


def _step(run: AnalysisRun, number: int, key: str, title: str, stages: tuple[str, ...]) -> dict:
    status = run.status()
    states = [status[s] for s in stages]
    if all(state == "done" for state in states):
        step_status, summary = "done", _summary(run, key)
    elif any(state == "needs re-running" for state in states):
        step_status, summary = "needs_rerun", DESCRIPTIONS[key]
    else:
        step_status = "not_run"
        partial = key == "quality_control" and status["qc_metrics"] == "done"
        summary = _qc_metrics_only(run) if partial else DESCRIPTIONS[key]
    return {
        "key": key,
        "number": number,
        "title": title,
        "status": step_status,
        "summary": summary,
        "reason": _cleared_reason(run, stages) if step_status == "needs_rerun" else None,
        "confirmation": _confirmation(run, key) if step_status == "done" else None,
        "prerequisite": _prerequisite(run, number) if step_status != "done" else None,
    }


# What to do next, per step: the action to take on that step's page.
NEXT_STEP_TEXT = {
    "data_overview": "Open Data Overview and click Load PBMC3k to open the file and check it holds raw counts.",
    "quality_control": ("Open Quality Control and click Run QC. Adjust the thresholds while watching the preview, "
                        "then click Apply Filtering."),
    "preprocessing": ("Open Preprocessing, review the settings if you like (the reference values are filled in), "
                      "then click Run Preprocessing. It takes a few seconds."),
    "dimensionality_reduction": ("Open Dimensionality Reduction, review the settings if you like, then click "
                                 "Run Dimensionality Reduction. It takes about 10 seconds."),
    "clustering": ("Open Clustering, choose a Leiden resolution if you like (the reference value is filled in), "
                   "then click Run Clustering."),
    "marker_genes": ("Open Marker Genes, keep the Wilcoxon test (the reference) or choose another, "
                     "then click Run Marker Genes."),
    "cell_type_annotation": ("Open Cell Type Annotation, review the thresholds if you like (the reference values "
                             "are filled in), then click Run Cell Type Annotation."),
}


def next_step(steps: list[dict]) -> dict | None:
    """The first step that isn't up to date, with a short instruction; None when all are done."""
    for step in steps:
        if step["status"] != "done":
            return {"key": step["key"], "title": step["title"], "text": NEXT_STEP_TEXT[step["key"]]}
    return None


def _confirmation(run: AnalysisRun, key: str) -> str | None:
    """A plain confirmation once a step has finished (more steps get one as they are built)."""
    if key == "data_overview":
        s = run.stage_result("load").summary
        return (f"Data loaded successfully: {s['n_cells']:,} cells × {s['n_genes']:,} genes. "
                "The file contains raw counts, so it's ready for Quality Control.")
    if key == "quality_control":
        return _filter_confirmation(run.stage_result("filter").summary)
    if key == "preprocessing":
        return _preprocess_confirmation(run.stage_result("preprocess"))
    if key == "dimensionality_reduction":
        return _reduction_confirmation(run.stage_result("reduce_dimensions"))
    if key == "cell_type_annotation":
        clusters = run.stage_result("annotate").summary["clusters"]
        unresolved = sum(1 for c in clusters if c["predicted_label"] == annotation_tools.UNRESOLVED)
        return (f"Cell types predicted from marker evidence: {len(clusters) - unresolved} of {len(clusters)} "
                f"clusters confident, {unresolved} Unresolved. These are predictions for you to review.")
    if key == "marker_genes":
        s, p = run.stage_result("markers").summary, run.stage_result("markers").settings
        test = TEST_NAMES.get(p["markers"]["method"], p["markers"]["method"])
        return (f"Marker Genes done: all {s['n_genes_tested']:,} genes ranked for each of the "
                f"{len(s['groups'])} clusters ({test}).")
    if key == "clustering":
        s = run.stage_result("cluster").summary
        sizes = s["cluster_sizes"].values()
        return (f"Clustering done: {s['n_clusters']} clusters, from {max(sizes):,} cells down to "
                f"{min(sizes):,} cells.")
    return None


def _filter_confirmation(s: dict) -> str:
    return (f"Filtering applied: kept {s['n_cells_kept']:,} of {s['n_cells_before']:,} cells and "
            f"{s['n_genes_kept']:,} of {s['n_genes_before']:,} genes.")


def _prerequisite(run: AnalysisRun, number: int) -> str | None:
    """Title of the first earlier step that isn't up to date (it has to come first)."""
    status = run.status()
    for _, title, stages in STEPS[: number - 1]:
        if any(status[s] != "done" for s in stages):
            return title
    return None


def _summary(run: AnalysisRun, key: str) -> str:
    """The step's key result in a few words, built from its real results and settings.
    (The full settings are on the step's page and in the settings record.)"""
    result = run.stage_result
    if key == "data_overview":
        s = result("load").summary
        return f"{s['n_cells']:,} cells × {s['n_genes']:,} genes · raw counts"
    if key == "quality_control":
        s = result("filter").summary
        removed_cells = s["n_cells_before"] - s["n_cells_kept"]
        return f"{s['n_cells_kept']:,} cells × {s['n_genes_kept']:,} genes kept ({removed_cells:,} cells removed)"
    if key == "preprocessing":
        s = result("preprocess").summary
        return f"{s['n_hvg']:,} highly variable genes · regression {'on' if s['regressed'] else 'off'}"
    if key == "dimensionality_reduction":
        s, p = result("reduce_dimensions").summary, result("reduce_dimensions").settings
        return f"{s['n_comps']} PCs ({p['neighbors']['n_pcs']} used) · {s['n_neighbors']} neighbors · UMAP"
    if key == "clustering":
        s, p = result("cluster").summary, result("cluster").settings
        return f"{s['n_clusters']} clusters (resolution {_decimal(p['leiden']['resolution'])})"
    if key == "marker_genes":
        s, p = result("markers").summary, result("markers").settings
        test = TEST_NAMES.get(p["markers"]["method"], p["markers"]["method"])
        return f"{s['n_genes_tested']:,} genes ranked per cluster · {test}"
    if key == "cell_type_annotation":
        clusters = result("annotate").summary["clusters"]
        unresolved = sum(1 for c in clusters if c["predicted_label"] == "Unresolved")
        return f"{len(clusters) - unresolved} of {len(clusters)} clusters confident · {unresolved} Unresolved"
    raise ValueError(f"Unknown step '{key}'")


def _qc_metrics_only(run: AnalysisRun) -> str:
    s = run.stage_result("qc_metrics").summary
    text = f"QC metrics calculated ({s['n_mito_genes']} mitochondrial genes); filtering not applied yet."
    return text + (" " + s["warnings"][0] if s["warnings"] else "")


def _cleared_reason(run: AnalysisRun, stages: tuple[str, ...]) -> str:
    """Which earlier step was redone, from the most recent history entry that cleared these."""
    for entry in reversed(run.history):
        if any(stage in entry["cleared_later_stages"] for stage in stages):
            redone = next(title for _, title, st in STEPS if entry["stage"] in st)
            return f"{redone} was redone, so run this step again."
    return "An earlier step was redone, so run this step again."


def settings_record(run: AnalysisRun) -> list[dict[str, str]]:
    """The current settings that most affect results, in plain words."""
    p = run.params
    f = p.filter
    return [
        {"label": "Mitochondrial gene prefix", "value": p.qc.mito_prefix},
        {"label": "Min genes per cell", "value": f"{f.min_genes:,} (inclusive)"},
        {"label": "Min cells per gene", "value": f"{f.min_cells:,} (inclusive)"},
        {"label": "Max genes per cell",
         "value": "no limit" if f.max_n_genes_by_counts is None else f"< {f.max_n_genes_by_counts:,}"},
        {"label": "Max % mitochondrial",
         "value": "no limit" if f.max_pct_counts_mt is None else f"< {f.max_pct_counts_mt:g}%"},
        {"label": "Principal components", "value": f"{p.pca.n_comps} computed, {p.neighbors.n_pcs} used"},
        {"label": "Neighbors per cell", "value": str(p.neighbors.n_neighbors)},
        {"label": "Leiden resolution", "value": _decimal(p.leiden.resolution)},
        {"label": "Random seed", "value": _seeds(p)},
    ]


def _seeds(p: AnalysisParams) -> str:
    """One number if every step uses the same seed, otherwise each step's seed."""
    seeds = {"PCA": p.pca.random_state, "neighbors": p.neighbors.random_state,
             "UMAP": p.umap.random_state, "Leiden": p.leiden.random_state}
    if len(set(seeds.values())) == 1:
        return str(p.pca.random_state)
    return ", ".join(f"{name} {seed}" for name, seed in seeds.items())


def run_log(run: AnalysisRun) -> list[dict[str, Any]]:
    """One entry per stage run, oldest first, from the pipeline's history."""
    entries = []
    for entry in run.history:
        text = _log_text(entry["stage"], entry["summary"])
        if entry["cleared_later_stages"]:
            cleared = ", ".join(STAGE_LABELS[s] for s in entry["cleared_later_stages"])
            text += f" · cleared {cleared}"
        entries.append({"label": STAGE_LABELS[entry["stage"]], "text": text, "ran_at": entry["ran_at"]})
    return entries


def _log_text(stage: str, s: dict) -> str:
    if stage == "load":
        return f"{s['n_cells']:,} cells × {s['n_genes']:,} genes"
    if stage == "qc_metrics":
        return f"{s['n_mito_genes']} mitochondrial genes found"
    if stage == "filter":
        return f"{s['n_cells_kept']:,} cells × {s['n_genes_kept']:,} genes kept"
    if stage == "preprocess":
        return f"{s['n_hvg']:,} highly variable genes · regression {'on' if s['regressed'] else 'off'}"
    if stage == "reduce_dimensions":
        return f"{s['n_comps']} PCs · {s['n_connections']:,} graph connections"
    if stage == "cluster":
        return f"{s['n_clusters']} clusters"
    if stage == "markers":
        return f"{s['n_genes_tested']:,} genes ranked per cluster"
    if stage == "annotate":
        return f"{sum(s['label_counts'].values()):,} cells annotated"
    return ""


# --- Quality Control page -------------------------------------------------------------------

class InputError(ValueError):
    """A value typed on screen can't be used. The message is written for the user."""


# Statistic rows of the QC summary table (same as the notebook's describe() table).
QC_STATISTICS = (("count", "Cells"), ("mean", "Mean"), ("std", "Std. dev."), ("min", "Min"),
                 ("25%", "25th percentile"), ("50%", "Median"), ("75%", "75th percentile"), ("max", "Max"))
REMOVAL_STEPS = {
    "cells: min_genes": ("cells", "Too few genes detected"),
    "genes: min_cells": ("genes", "Detected in too few cells"),
    "cells: max_n_genes_by_counts": ("cells", "Too many genes (possible doublets)"),
    "cells: max_pct_counts_mt": ("cells", "Too high % mitochondrial"),
}


def filter_params_from_browser(values: dict[str, Any] | None, current: FilterParams) -> FilterParams:
    """The filter boxes on the QC page -> backend settings. Empty boxes are None.

    `values` holds min_genes, min_cells, max_genes, max_mt and the two "no limit" flags
    max_genes_off and max_mt_off. Missing values keep `current`. Raises `InputError`.
    """
    if not values:
        return current

    def number(key, label, *, whole):
        return _number(values, key, label, whole=whole)

    return FilterParams(
        min_genes=number("min_genes", "Min genes per cell", whole=True),
        min_cells=number("min_cells", "Min cells per gene", whole=True),
        max_n_genes_by_counts=None if values.get("max_genes_off") else number("max_genes", "Max genes per cell", whole=True),
        max_pct_counts_mt=None if values.get("max_mt_off") else number("max_mt", "Max % mitochondrial", whole=False),
    )


def _number(values: dict[str, Any], key: str, label: str, *, whole: bool, minimum: float = 0,
            above: bool = False) -> float | int:
    """One value typed on screen, checked. `minimum` is the smallest allowed value
    (`above=True`: the value must be larger than it). Raises `InputError`."""
    value = values.get(key)
    if value is None or value == "":
        raise InputError(f"Enter a value for '{label}'.")
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise InputError(f"'{label}' must be a number.") from None
    if not np.isfinite(value):
        raise InputError(f"'{label}' must be a number.")
    if above and value <= minimum:
        raise InputError(f"'{label}' must be greater than {minimum:g}.")
    if not above and value < minimum:
        raise InputError(f"'{label}' can't be negative." if minimum == 0
                         else f"'{label}' must be at least {minimum:g}.")
    if whole and value != int(value):
        raise InputError(f"'{label}' must be a whole number.")
    return int(value) if whole else value


def qc_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Quality Control page shows, except the live filter preview."""
    status = run.status()
    payload: dict[str, Any] = {
        "metrics_done": status["qc_metrics"] == "done",
        "filter_status": {"done": "done", "needs re-running": "needs_rerun"}.get(status["filter"], "not_run"),
        "can_run": status["load"] == "done",
        "filter": {"current": asdict(run.params.filter), "reference": asdict(REFERENCE_PRESET.filter)},
        "applied": None,
        "will_clear": _will_clear(run, "filter"),  # later steps that re-applying would clear
    }
    if status["filter"] == "done":
        s = run.stage_result("filter").summary
        payload["applied"] = {**s, "confirmation": _filter_confirmation(s)}
    if not payload["metrics_done"]:
        return payload

    summary = run.stage_result("qc_metrics").details[0]
    obs = run.result("qc_metrics").obs
    payload.update({
        "n_cells": int(obs.shape[0]),
        "n_genes": int(run.result("qc_metrics").n_vars),
        "mito": {"genes": list(summary.mito_genes), "available": summary.mito_qc_available,
                 "prefix": summary.mito_prefix, "warnings": list(summary.warnings)},
        "stats": [
            {"statistic": label, **{m: float(stats[name]) for m, stats in summary.metric_stats.items()}}
            for name, label in QC_STATISTICS
        ],
        "cells": {
            "barcode": obs.index.astype(str).tolist(),
            "n_genes": obs["n_genes_by_counts"].astype(int).tolist(),
            "total_counts": [float(v) for v in obs["total_counts"]],
            "pct_mt": [float(v) for v in obs["pct_counts_mt"]],
        },
    })
    return payload


def preview_payload(run: AnalysisRun, values: dict[str, Any] | None) -> dict[str, Any] | None:
    """What the filter boxes would keep, computed without changing anything (None before QC)."""
    if run.status()["qc_metrics"] != "done":
        return None
    try:
        params = filter_params_from_browser(values, run.params.filter)
        preview = run.preview_filter(params)
    except Exception as problem:  # InputError, or the backend explaining why (e.g. removes all cells)
        return {"ok": False, "message": str(problem)}
    # Which cell filters each removed cell fails (same rules as the preview; a cell can fail
    # several). Keyed by barcode, so the page matches cells by identity, never by position.
    failures = qc_tools.cell_filter_failures(run.result("qc_metrics"), params)
    failing = failures[failures.any(axis=1)]
    # The same preview at the validated reference settings, for "vs. reference" on the page.
    # Computed from this data with the preset's settings (nothing is typed in).
    reference = run.preview_filter(REFERENCE_PRESET.filter)
    return {
        "ok": True,
        "settings": asdict(params),
        "n_cells_before": preview.n_cells_before,
        "n_cells_kept": preview.n_cells_kept,
        "n_genes_before": preview.n_genes_before,
        "n_genes_kept": preview.n_genes_kept,
        "reference": {"n_cells_kept": reference.n_cells_kept, "n_genes_kept": reference.n_genes_kept},
        "filters": filter_table(params, preview, REFERENCE_PRESET.filter, reference),
        "removed_steps": [
            {"what": REMOVAL_STEPS[key][0], "reason": REMOVAL_STEPS[key][1], "removed": removed}
            for key, removed in preview.removed_in_order.items()
        ],
        # The cells these settings remove, by barcode, each with the filters it fails.
        "removed_barcodes": failing.index.astype(str).tolist(),
        "failures": {str(barcode): [name for name, failed in row.items() if failed]
                     for barcode, row in failing.iterrows()},
    }


# --- Preprocessing and Dimensionality Reduction cards ------------------------------------------

def _stage_status(run: AnalysisRun, stage: str) -> str:
    return {"done": "done", "needs re-running": "needs_rerun"}.get(run.status()[stage], "not_run")


def _will_clear(run: AnalysisRun, stage: str) -> list[str]:
    """Titles of the later Workflow steps that running `stage` again would clear."""
    status = run.status()
    later = STAGES_IN_ORDER[STAGES_IN_ORDER.index(stage) + 1:]
    return [title for _, title, stages in STEPS
            if any(s in later and status[s] == "done" for s in stages)]


STAGES_IN_ORDER = tuple(stage for _, _, stages in STEPS for stage in stages)


# --- "Run all remaining steps" -------------------------------------------------------------
# Runs every stage that isn't up to date, in order, with the settings currently in use
# (run.params: the validated ones unless the user changed a step). The server runs one stage
# at a time (app.py), so the page can show progress in between.

# The AnalysisRun method behind each stage. Called without arguments, each one uses the
# settings in run.params.
STAGE_METHODS = {
    "load": "load", "qc_metrics": "compute_qc_metrics", "filter": "apply_filter",
    "preprocess": "preprocess", "reduce_dimensions": "reduce_dimensions", "cluster": "cluster",
    "markers": "find_markers", "annotate": "annotate",
}


def remaining_stages(run: AnalysisRun) -> list[str]:
    """The stages that aren't up to date, in the order they run."""
    status = run.status()
    return [stage for stage in STAGES_IN_ORDER if status[stage] != "done"]


def run_stage(run: AnalysisRun, stage: str) -> None:
    """Run one stage with the current settings."""
    getattr(run, STAGE_METHODS[stage])()


def step_of_stage(stage: str) -> tuple[str, str]:
    """The Workflow step (key, title) a pipeline stage belongs to."""
    return next((key, title) for key, title, stages in STEPS if stage in stages)


def run_remaining_payload(job: dict[str, Any] | None) -> dict[str, Any] | None:
    """What the page shows about the latest job: which step is running ("step 3 of 5"),
    or how it ended. `job` holds the stages it started with and the ones still to run."""
    if job is None:
        return None
    steps = list(dict.fromkeys(step_of_stage(stage) for stage in job["stages"]))  # in order, once each
    current = step_of_stage(job["left"][0]) if job["left"] and job["status"] == "running" else None
    return {
        "id": job["id"],
        "status": job["status"],  # "running", "done" or "failed"
        "current": current[1] if current else None,
        "index": steps.index(current) + 1 if current else None,
        "total": len(steps),
        "failed_step": job.get("failed_step"),
    }


def preprocess_settings(params: AnalysisParams) -> dict[str, Any]:
    """The Preprocessing settings shown on screen."""
    return {
        "target_sum": params.normalize.target_sum,
        "min_mean": params.hvg.min_mean,
        "max_mean": params.hvg.max_mean,
        "min_disp": params.hvg.min_disp,
        "max_value": params.scale.max_value,
        "regress": params.regress_out.enabled,
    }


def preprocess_params_from_browser(
    values: dict[str, Any] | None, current: AnalysisParams
) -> tuple[NormalizeParams, HVGParams, RegressOutParams, ScaleParams]:
    """The Preprocessing boxes -> backend settings. Settings not on screen (HVG flavor,
    regression keys) keep their current values. Raises `InputError`."""
    if not values:
        return current.normalize, current.hvg, current.regress_out, current.scale
    target_sum = _number(values, "target_sum", "Counts per cell after normalization", whole=False, above=True)
    min_mean = _number(values, "min_mean", "Min mean expression", whole=False)
    max_mean = _number(values, "max_mean", "Max mean expression", whole=False, above=True)
    if min_mean >= max_mean:
        raise InputError("'Min mean expression' must be smaller than 'Max mean expression'.")
    min_disp = _number(values, "min_disp", "Min normalized dispersion", whole=False, minimum=-np.inf)
    max_value = _number(values, "max_value", "Clip scaled values at", whole=False, above=True)
    return (
        replace(current.normalize, target_sum=target_sum),
        replace(current.hvg, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp),
        replace(current.regress_out, enabled=bool(values.get("regress"))),
        replace(current.scale, max_value=max_value),
    )


def preprocess_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Preprocessing card shows: settings, and results once it has run."""
    payload: dict[str, Any] = {
        "status": _stage_status(run, "preprocess"),
        "settings": {"current": preprocess_settings(run.params), "reference": preprocess_settings(REFERENCE_PRESET)},
        "fixed": {"hvg_flavor": run.params.hvg.flavor, "regress_keys": list(run.params.regress_out.keys)},
        "will_clear": _will_clear(run, "preprocess"),
        "results": None,
    }
    if run.status()["preprocess"] != "done":
        return payload
    result = run.stage_result("preprocess")
    genes = result.adata.raw.var  # all genes, with the HVG statistics
    payload["results"] = {
        "n_hvg": result.summary["n_hvg"],
        "n_genes": int(genes.shape[0]),
        "settings": preprocess_settings(run.params),
        "genes": {
            "name": genes.index.astype(str).tolist(),
            "mean": _floats(genes["means"]),
            "dispersion_norm": _floats(genes["dispersions_norm"]),
            "hvg": genes["highly_variable"].astype(bool).tolist(),
        },
    }
    return payload


def _preprocess_confirmation(result) -> str:
    s, p = result.summary, result.settings
    regression = ("regressed out " + " and ".join(s["regress_keys"])) if s["regressed"] else "regression off"
    return (f"Preprocessing done: {s['n_hvg']:,} of {s['raw_shape'][1]:,} genes selected as highly variable; "
            f"{regression}; values scaled and capped at ±{p['scale']['max_value']:g}.")


def reduction_settings(params: AnalysisParams) -> dict[str, Any]:
    """The Dimensionality Reduction settings shown on screen (one seed for PCA, neighbors and UMAP)."""
    return {
        "n_comps": params.pca.n_comps,
        "n_neighbors": params.neighbors.n_neighbors,
        "n_pcs": params.neighbors.n_pcs,
        "seed": params.pca.random_state,
    }


def reduction_params_from_browser(
    values: dict[str, Any] | None, current: AnalysisParams
) -> tuple[PCAParams, NeighborsParams, UMAPParams]:
    """The Dimensionality Reduction boxes -> backend settings. Raises `InputError`."""
    if not values:
        return current.pca, current.neighbors, current.umap
    n_comps = _number(values, "n_comps", "Principal components", whole=True, minimum=1)
    n_neighbors = _number(values, "n_neighbors", "Neighbors per cell", whole=True, minimum=2)
    n_pcs = _number(values, "n_pcs", "PCs used for the neighbor graph", whole=True, minimum=1)
    if n_pcs > n_comps:
        raise InputError("'PCs used for the neighbor graph' can't be more than 'Principal components'.")
    seed = _number(values, "seed", "Random seed", whole=True)
    return (
        replace(current.pca, n_comps=n_comps, random_state=seed),
        replace(current.neighbors, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed),
        replace(current.umap, random_state=seed),
    )


def reduction_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Dimensionality Reduction card shows."""
    payload: dict[str, Any] = {
        "status": _stage_status(run, "reduce_dimensions"),
        "settings": {"current": reduction_settings(run.params), "reference": reduction_settings(REFERENCE_PRESET)},
        "will_clear": _will_clear(run, "reduce_dimensions"),
        "results": None,
    }
    if run.status()["reduce_dimensions"] != "done":
        return payload
    result = run.stage_result("reduce_dimensions")
    adata, (pcs, graph, _) = result.adata, result.details
    obs, umap = adata.obs, adata.obsm["X_umap"]
    payload["results"] = {
        "variance_ratio": list(pcs.variance_ratio),
        "total_variance_ratio": pcs.total_variance_ratio,
        "n_pcs_used": graph.n_pcs,
        "umap": {"x": _floats(umap[:, 0]), "y": _floats(umap[:, 1])},
        "cells": {
            "barcode": obs.index.astype(str).tolist(),
            "n_genes": obs["n_genes_by_counts"].astype(int).tolist(),
            "total_counts": _floats(obs["total_counts"]),
            "pct_mt": _floats(obs["pct_counts_mt"]) if "pct_counts_mt" in obs else None,
        },
    }
    return payload


def _reduction_confirmation(result) -> str:
    s, p = result.summary, result.settings
    n_cells = result.adata.n_obs
    return (f"Dimensionality Reduction done: {s['n_comps']} principal components explain "
            f"{100 * s['total_variance_ratio']:.1f}% of the variance; neighbor graph with "
            f"{s['n_neighbors']} neighbors per cell from the first {p['neighbors']['n_pcs']} PCs; "
            f"UMAP drawn for {n_cells:,} cells.")


# --- Clustering page ---------------------------------------------------------------------------

def cluster_settings(params: AnalysisParams) -> dict[str, Any]:
    """The Clustering settings shown on screen (the editable ones)."""
    return {"resolution": params.leiden.resolution, "seed": params.leiden.random_state}


def cluster_params_from_browser(values: dict[str, Any] | None, current: AnalysisParams) -> tuple[LeidenParams]:
    """The Clustering boxes -> backend settings. Algorithm settings (flavor, iterations,
    directed graph) are not on screen and keep their current values. Raises `InputError`."""
    if not values:
        return (current.leiden,)
    resolution = _number(values, "resolution", "Leiden resolution", whole=False, above=True)
    seed = _number(values, "seed", "Leiden random seed", whole=True)
    return (replace(current.leiden, resolution=resolution, random_state=seed),)


def cluster_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Clustering page shows."""
    leiden = run.params.leiden
    payload: dict[str, Any] = {
        "status": _stage_status(run, "cluster"),
        "settings": {"current": cluster_settings(run.params), "reference": cluster_settings(REFERENCE_PRESET)},
        "fixed": {"flavor": leiden.flavor, "n_iterations": leiden.n_iterations, "directed": leiden.directed},
        "will_clear": _will_clear(run, "cluster"),
        "results": None,
    }
    if run.status()["cluster"] != "done":
        return payload
    result = run.stage_result("cluster")
    adata, summary = result.adata, result.details[0]
    obs, umap = adata.obs, adata.obsm["X_umap"]
    payload["results"] = {
        "settings": cluster_settings(run.params),
        "n_clusters": summary.n_clusters,
        "clusters": [{"cluster": label, "n_cells": n} for label, n in summary.cluster_sizes.items()],
        "notes": list(summary.notes),
        # One entry per cell, in the same order as `umap` and `cells`.
        "labels": obs[summary.key_added].astype(str).tolist(),
        "umap": {"x": _floats(umap[:, 0]), "y": _floats(umap[:, 1])},
        "cells": {
            "barcode": obs.index.astype(str).tolist(),
            "n_genes": obs["n_genes_by_counts"].astype(int).tolist(),
            "total_counts": _floats(obs["total_counts"]),
            "pct_mt": _floats(obs["pct_counts_mt"]) if "pct_counts_mt" in obs else None,
        },
    }
    return payload


# --- Marker Genes page ------------------------------------------------------------------------

TOP_MARKERS_SENT = 50  # top genes per cluster sent to the page (the page shows fewer by default)
DOTPLOT_GENES_PER_CLUSTER = 3


def marker_settings(params: AnalysisParams) -> dict[str, Any]:
    return {"method": params.markers.method}


def marker_params_from_browser(values: dict[str, Any] | None, current: AnalysisParams) -> tuple[MarkerParams]:
    """The Marker Genes settings -> backend settings (only the test can change). Raises `InputError`."""
    if not values:
        return (current.markers,)
    method = values.get("method")
    if method not in dict(MARKER_METHODS):
        raise InputError(f"Choose one of the tests: {', '.join(label for _, label in MARKER_METHODS)}.")
    return (replace(current.markers, method=method),)


def markers_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Marker Genes page shows."""
    m = run.params.markers
    payload: dict[str, Any] = {
        "status": _stage_status(run, "markers"),
        "settings": {"current": marker_settings(run.params), "reference": marker_settings(REFERENCE_PRESET)},
        "methods": [{"value": v, "label": label} for v, label in MARKER_METHODS],
        "fixed": {"use_raw": m.use_raw, "reference": m.reference, "corr_method": m.corr_method,
                  "tie_correct": m.tie_correct, "groupby": m.groupby},
        "will_clear": _will_clear(run, "markers"),
        "results": None,
    }
    if run.status()["markers"] != "done":
        return payload
    result = run.stage_result("markers")
    adata, summary = result.adata, result.details[0]
    table = marker_tools.marker_table(adata)
    groups = [str(g) for g in summary.groups]
    sizes = adata.obs[summary.groupby].astype(str).value_counts()

    def rows(frame):
        return [{"rank": int(r.rank), "gene": str(r.gene), "score": _num(r.score),
                 "logfoldchange": _num(r.logfoldchange), "pval_adj": _num(r.pval_adj),
                 "frac_in_group": _num(getattr(r, "frac_in_group", None)),
                 "frac_in_rest": _num(getattr(r, "frac_in_rest", None))}
                for r in frame.itertuples(index=False)]

    top = {g: rows(table[(table["group"] == g) & (table["rank"] <= TOP_MARKERS_SENT)]) for g in groups}
    # Dot plot: for each cluster, its first top-ranked genes not already in the plot.
    dot_genes: list[str] = []
    for g in groups:
        added = 0
        for gene in table[table["group"] == g].sort_values("rank")["gene"].astype(str):
            if added == DOTPLOT_GENES_PER_CLUSTER:
                break
            if gene not in dot_genes:
                dot_genes.append(gene)
                added += 1
    expr = marker_tools.expression_by_group(adata, dot_genes, summary.groupby, use_raw=summary.use_raw)
    ranks = table[table["gene"].isin(dot_genes)].set_index(["group", "gene"])
    dots = [{"group": r.group, "gene": r.gene, "mean_expr": r.mean_expr, "frac_expressing": r.frac_expressing,
             "rank": int(ranks.loc[(r.group, r.gene), "rank"]),
             "logfoldchange": _num(ranks.loc[(r.group, r.gene), "logfoldchange"])}
            for r in expr.itertuples(index=False)]
    payload["results"] = {
        "settings": marker_settings(run.params),
        "method_label": dict(MARKER_METHODS).get(summary.method, summary.method),
        "n_genes_tested": summary.n_genes_tested,
        "groups": [{"cluster": g, "n_cells": int(sizes.get(g, 0))} for g in groups],
        "top": top,
        "dotplot": {"genes": dot_genes, "values": dots},
        "notes": list(summary.notes),
    }
    return payload


# --- Cell Type Annotation page ------------------------------------------------------------------

def annotation_settings(params: AnalysisParams) -> dict[str, Any]:
    """The editable annotation thresholds (marker sets are shown read-only)."""
    a = params.annotation
    return {"ambiguity_margin": a.ambiguity_margin, "max_de_rank": a.max_de_rank,
            "min_score": a.min_score, "conflict_z": a.conflict_z}


def annotation_params_from_browser(values: dict[str, Any] | None, current: AnalysisParams) -> tuple[AnnotationParams]:
    """The threshold boxes -> backend settings; marker sets keep their current values.
    `min_score_off` = no minimum-score check. Raises `InputError`."""
    if not values:
        return (current.annotation,)
    margin = _number(values, "ambiguity_margin", "Ambiguity margin", whole=False)
    max_rank = _number(values, "max_de_rank", "Max DE rank", whole=True, minimum=1)
    min_score = None if values.get("min_score_off") else _number(values, "min_score", "Min score",
                                                                  whole=False, minimum=-np.inf)
    conflict_z = _number(values, "conflict_z", "Conflict z-score", whole=False, minimum=-np.inf)
    return (replace(current.annotation, ambiguity_margin=margin, max_de_rank=max_rank,
                    min_score=min_score, conflict_z=conflict_z),)


def _marker_sets(params: AnnotationParams, available: set[str] | None) -> dict[str, Any]:
    """The marker sets for display; `present` is None until the genes in the data are known."""
    def genes(gs):
        return [{"gene": g, "present": None if available is None else g in available} for g in gs]
    return {
        "sets": [{"cell_type": ct, "markers": genes(gs)} for ct, gs in params.marker_sets.items()],
        "context": genes(params.context_markers),
    }


def annotation_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Cell Type Annotation page shows. Predictions only: reviewed labels are
    never part of the app."""
    status = run.status()
    available = None
    if status["markers"] == "done":
        source = run.result("markers")
        available = set((source.raw if source.raw is not None else source).var_names)
    payload: dict[str, Any] = {
        "status": _stage_status(run, "annotate"),
        "settings": {"current": annotation_settings(run.params), "reference": annotation_settings(REFERENCE_PRESET)},
        "marker_sets": _marker_sets(run.params.annotation, available),
        "unresolved_label": annotation_tools.UNRESOLVED,
        "will_clear": _will_clear(run, "annotate"),
        "results": None,
    }
    if status["annotate"] != "done":
        return payload
    result = run.stage_result("annotate")
    adata, summary = result.adata, result.details[0]
    clusters = summary.clusters
    columns = ["cluster", "n_cells", "predicted_label", "status", "reasons", "best_candidate", "second_candidate",
               "score", "margin", "best_marker_DE_rank", "supporting_markers", "conflicting_markers"]
    evidence = annotation_tools.cluster_marker_evidence(adata, run.params.annotation)
    obs, umap = adata.obs, adata.obsm["X_umap"]
    payload["results"] = {
        "settings": annotation_settings(run.params),
        "clusters": json.loads(clusters[columns].to_json(orient="records")),
        "scores": {str(c): {ct: _num(v) for ct, v in row.items()} for c, row in summary.scores.iterrows()},
        "evidence": json.loads(evidence.to_json(orient="records")),
        "label_counts": summary.label_counts,
        "missing_genes": {ct: list(gs) for ct, gs in summary.missing_genes.items()},
        "notes": list(summary.notes),
        # One entry per cell, in the same order as `umap` and `cells`.
        "cell_cluster": obs[summary.groupby].astype(str).tolist(),
        "cell_label": obs[annotation_tools.PREDICTION_KEY].astype(str).tolist(),
        "umap": {"x": _floats(umap[:, 0]), "y": _floats(umap[:, 1])},
        "cells": {
            "barcode": obs.index.astype(str).tolist(),
            "n_genes": obs["n_genes_by_counts"].astype(int).tolist(),
            "total_counts": _floats(obs["total_counts"]),
            "pct_mt": _floats(obs["pct_counts_mt"]) if "pct_counts_mt" in obs else None,
        },
    }
    return payload


# --- Explore Results page ---------------------------------------------------------------------

EXPLORE_TOP_MARKERS = 5


def explore_payload(run: AnalysisRun) -> dict[str, Any]:
    """What Explore Results shows: the clustered cells on the UMAP, plus cell types and top
    markers when those steps are up to date. Nothing is computed here that the steps
    didn't compute; results of steps that need re-running are left out (with the reason)."""
    status = run.status()
    if status["cluster"] != "done":
        return {"available": False, "reason": "Run Clustering first: Explore Results shows clustered cells.",
                "results": None}
    clustered = run.stage_result("cluster")
    adata, summary = clustered.adata, clustered.details[0]
    obs, umap = adata.obs, adata.obsm["X_umap"]
    raw = adata.raw if adata.raw is not None else adata
    results: dict[str, Any] = {
        "clusters": [{"cluster": c, "n_cells": n} for c, n in summary.cluster_sizes.items()],
        "cell_cluster": obs[summary.key_added].astype(str).tolist(),
        "umap": {"x": _floats(umap[:, 0]), "y": _floats(umap[:, 1])},
        "cells": {
            "barcode": obs.index.astype(str).tolist(),
            "n_genes": obs["n_genes_by_counts"].astype(int).tolist(),
            "total_counts": _floats(obs["total_counts"]),
            "pct_mt": _floats(obs["pct_counts_mt"]) if "pct_counts_mt" in obs else None,
        },
        "genes": sorted(map(str, raw.var_names)),  # for the gene search
        "top_markers": None,
        "predictions": None,
        "unavailable": [],
    }
    if status["markers"] == "done":
        table = marker_tools.marker_table(run.result("markers"), n_top=EXPLORE_TOP_MARKERS)
        results["top_markers"] = {str(g): list(map(str, rows["gene"])) for g, rows in table.groupby("group", sort=False)}
    else:
        results["unavailable"].append("Marker Genes " + ("needs re-running" if status["markers"] == "needs re-running" else "hasn't been run"))
    if status["annotate"] == "done":
        result = run.stage_result("annotate")
        clusters = result.details[0].clusters
        results["predictions"] = {
            "unresolved_label": annotation_tools.UNRESOLVED,
            "clusters": json.loads(clusters[["cluster", "predicted_label", "status"]].to_json(orient="records")),
            "cell_label": result.adata.obs[annotation_tools.PREDICTION_KEY].astype(str).tolist(),
            "marker_set_order": list(run.params.annotation.marker_sets),
        }
    else:
        results["unavailable"].append("Cell Type Annotation " + ("needs re-running" if status["annotate"] == "needs re-running" else "hasn't been run"))
    return {"available": True, "reason": None, "results": results}


def gene_payload(run: AnalysisRun, gene: str | None) -> dict[str, Any] | None:
    """One gene's log-normalized expression in every clustered cell, plus mean and fraction
    expressing per cluster (the same helper as the dot plots). None when no gene is asked."""
    if not gene or run.status()["cluster"] != "done":
        return None
    clustered = run.stage_result("cluster")
    adata, summary = clustered.adata, clustered.details[0]
    use_raw = adata.raw is not None
    names = list(map(str, (adata.raw if use_raw else adata).var_names))
    wanted = gene.strip()
    if wanted not in names:
        by_upper = {n.upper(): n for n in names}
        if wanted.upper() in by_upper:
            wanted = by_upper[wanted.upper()]
        else:
            import difflib
            close = difflib.get_close_matches(wanted.upper(), list(by_upper), n=5, cutoff=0.6)
            return {"ok": False, "gene": gene, "message": f"Gene '{gene}' isn't in this data.",
                    "suggestions": [by_upper[c] for c in close]}
    values = sc_obs_values(adata, wanted, use_raw)
    groups = marker_tools.expression_by_group(adata, [wanted], summary.key_added, use_raw=use_raw)
    return {
        "ok": True,
        "gene": wanted,
        "values": _floats(values),
        "max": float(np.max(values)) if len(values) else 0.0,
        "n_expressing": int(np.sum(values > 0)),
        "by_cluster": [{"cluster": r.group, "mean_expr": r.mean_expr, "frac_expressing": r.frac_expressing}
                       for r in groups.itertuples(index=False)],
    }


def sc_obs_values(adata, gene: str, use_raw: bool) -> np.ndarray:
    """A gene's values for every cell (from adata.raw when use_raw), as a plain array."""
    import scanpy as sc
    return sc.get.obs_df(adata, keys=[gene], use_raw=use_raw)[gene].to_numpy(dtype=float)


def _num(value) -> float | None:
    """One number for JSON (None for missing or NaN)."""
    if value is None:
        return None
    value = float(value)
    return value if np.isfinite(value) else None


def _floats(values) -> list[float | None]:
    """Numbers for JSON (NaN, which JSON can't hold, becomes None)."""
    return [float(v) if np.isfinite(v) else None for v in np.asarray(values, dtype=float)]


def filter_table(params: FilterParams, preview, reference_params: FilterParams, reference) -> list[dict[str, Any]]:
    """One row per filter, in the order they run: your setting and what it removes, next to the
    reference setting and what that removes. A filter switched off has setting None and
    removes None (it is not run)."""
    rows = []
    for key, (what, reason) in REMOVAL_STEPS.items():
        name = key.split(": ")[1]  # "cells: max_pct_counts_mt" -> "max_pct_counts_mt"
        setting, reference_setting = getattr(params, name), getattr(reference_params, name)
        rows.append({
            "key": name,
            "what": what,
            "reason": reason,
            "setting": setting,
            "reference_setting": reference_setting,
            "removed": preview.removed_in_order.get(key),
            "reference_removed": reference.removed_in_order.get(key),
            "differs": setting != reference_setting,
        })
    return rows


def _decimal(value: float) -> str:
    """A setting that is a decimal, always shown with its decimal point (1.0, 0.5, 0.75)."""
    return repr(float(value))


def _join(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]
