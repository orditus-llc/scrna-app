"""The data the React interface receives: pipeline state -> plain JSON.

This is the "data contract" between the Shiny server and the browser. Each function takes
an `AnalysisRun` and returns JSON-serializable dicts and lists only. All numbers come from
the run's results and settings; nothing here is PBMC3k-specific. No Shiny code.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, replace
from typing import Any

import numpy as np

from scrna_analysis.params import (
    DEFAULT_N_TOP_GENES,
    DEFAULT_PARAMS,
    LOW_MIXING,
    MAX_SAMPLES,
    MITO_PATTERNS,
    AnalysisParams,
    AnnotationParams,
    FilterParams,
    HVGParams,
    IntegrationParams,
    LeidenParams,
    MarkerParams,
    NeighborsParams,
    NormalizeParams,
    PCAParams,
    QCParams,
    RegressOutParams,
    SampleLimits,
    ScaleParams,
    SuggestionParams,
    UMAPParams,
    filter_from_dict,
)
from scrna_analysis import annotation as annotation_tools
from scrna_analysis import celltypist_models as celltypist_tools
from scrna_analysis import census_markers as census_tools
from scrna_analysis import integration as integration_tools
from scrna_analysis import load_data as load_tools
from scrna_analysis import marker_db as marker_db_tools
from scrna_analysis import markers as marker_tools
from scrna_analysis import qc as qc_tools
from scrna_analysis import review as review_tools
from scrna_analysis.pipeline import MAX_ABOUT_LENGTH, AnalysisRun, package_versions, records

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
    "annotate": "Cell Type Annotation", "review": "Label review", "metadata": "Data Overview",
}
# The marker tests offered on the Marker Genes page: (Scanpy name, label shown).
MARKER_METHODS = [("wilcoxon", "Wilcoxon rank-sum"), ("t-test", "t-test"),
                  ("t-test_overestim_var", "t-test (overestimated variance)"), ("logreg", "Logistic regression")]
TEST_NAMES = dict(MARKER_METHODS)


def workflow_payload(run: AnalysisRun, notice: dict[str, Any] | None = None,
                     job: dict[str, Any] | None = None, dataset_label: str | None = None) -> dict[str, Any]:
    """Everything the Workflow page shows. `notice` is the latest message for the user;
    `job` is the latest "Run all remaining steps" job (see `run_remaining_payload`);
    `dataset_label` is how the loaded data is named on screen (default: its file name)."""
    load = run.stage_result("load").summary if run.status()["load"] == "done" else None
    steps = [_step(run, number, key, title, stages)
             for number, (key, title, stages) in enumerate(STEPS, start=1)]
    return {
        "dataset": {
            "file": run.file_name,
            "label": dataset_label or run.file_name,
            "file_format": load["file_format"] if load else None,
            "notes": load["notes"] if load else [],
            # Which known mitochondrial gene-name patterns the data uses (most matches first).
            "mito_matches": qc_tools.mito_prefix_matches(run.result("load").var_names) if load else [],
            "n_cells": load["n_cells"] if load else None,
            "n_genes": load["n_genes"] if load else None,
            "fingerprint": load["input_sha256"][:12] if load else None,
            "metadata": metadata_payload(run) if load else None,
            "samples": samples_payload(run) if load else None,
        },
        "steps": steps,
        "next_step": next_step(steps),
        "settings_record": settings_record(run),
        "versions": {name: package_versions()[name] for name in ("scanpy", "leidenalg", "umap-learn")},
        "run_log": run_log(run),
        "notice": notice,
        "remaining_steps": sum(step["status"] != "done" for step in steps),
        "run_remaining": run_remaining_payload(job),
        "export": export_payload(run),
    }


def export_payload(run: AnalysisRun) -> dict[str, Any]:
    """Whether the analysis can be downloaded yet, what the files would hold, and their names."""
    status = run.status()
    included = []
    for key, title, stages in STEPS[2:]:  # Preprocessing onward, up to date, without gaps
        if not all(status[s] == "done" for s in stages):
            break
        included.append(title)
    available = status["preprocess"] == "done"
    return {
        "available": available,
        "reason": None if available else "Run Preprocessing first: the file holds the normalized data.",
        "includes": included,
        "reviewed": bool(run.review),
        "h5ad_name": export_file_name(run.file_name, "analyzed", "h5ad"),
        "report_name": export_file_name(run.file_name, "run report", "html"),
    }


def export_file_name(data_file: str | None, what: str, extension: str) -> str:
    """"<dataset>_<what>.<ext>", e.g. "pbmc3k_analyzed.h5ad" (as the page names chart downloads)."""
    def slug(text: str) -> str:
        text = re.sub(r"\.(h5ad|h5)$", "", text.lower())
        return re.sub(r"[^a-z0-9]+", "-", text).strip("-") or "analysis"
    return f"{slug(data_file) + '_' if data_file else ''}{slug(what)}.{extension}"


def _step(run: AnalysisRun, number: int, key: str, title: str, stages: tuple[str, ...]) -> dict:
    status = run.status()
    states = [status[s] for s in stages]
    unconfirmed = needs_confirmation(run)
    if key == "data_overview" and unconfirmed:
        step_status, summary = "not_run", "Data loaded · confirm the species and tissue"
    elif all(state == "done" for state in states):
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
        # Why the prerequisite isn't met, when it's the confirmation rather than a run.
        "prerequisite_text": CONFIRM_FIRST if unconfirmed and number > 1 else None,
        # Later steps can't be opened until the data is loaded and its metadata confirmed.
        "locked": number > 1 and (status["load"] != "done" or unconfirmed),
        # What the step is waiting for, shown instead of its status (sidebar).
        "waiting": "Confirm species and tissue" if key == "data_overview" and unconfirmed else None,
    }


CONFIRM_FIRST = "Confirm the species and tissue on Data Overview first."


def needs_confirmation(run: AnalysisRun) -> bool:
    """Data is loaded, but its species and tissue still need the user's confirmation."""
    return run.require_metadata and run.status()["load"] == "done" and not run.metadata_confirmed


# What to do next, per step: the action to take on that step's page.
NEXT_STEP_TEXT = {
    "data_overview": "Open Data Overview and upload a file (or use the demo data); it's checked for raw counts.",
    "quality_control": ("Open Quality Control and click Run QC. Adjust the thresholds while watching the preview, "
                        "then click Apply Filtering."),
    "preprocessing": ("Open Preprocessing, review the settings if you like (the reference values are filled in), "
                      "then click Run Preprocessing. It takes a few seconds."),
    "dimensionality_reduction": ("Open Dimensionality Reduction, review the settings if you like, then click "
                                 "Run Dimensionality Reduction. It takes about 10 seconds."),
    "clustering": ("Open Clustering, choose a Leiden resolution if you like (the default is filled in), "
                   "then click Run Clustering."),
    "marker_genes": ("Open Marker Genes, keep the Wilcoxon test (the default) or choose another, "
                     "then click Run Marker Genes."),
    "cell_type_annotation": ("Open Cell Type Annotation, check the species, choose the tissues, "
                             "then click Run Cell Type Annotation."),
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
        about = f" ({metadata_text(run.metadata)})" if run.metadata else ""
        samples = f" in {len(s['samples'])} samples" if s.get("samples") else ""
        return (f"Data loaded: {s['n_cells']:,} cells × {s['n_genes']:,} genes of raw counts{samples}{about}. "
                "Ready for Quality Control.")
    if key == "quality_control":
        return _filter_confirmation(run.stage_result("filter").summary)
    if key == "preprocessing":
        return _preprocess_confirmation(run.stage_result("preprocess"))
    if key == "dimensionality_reduction":
        return _reduction_confirmation(run.stage_result("reduce_dimensions"))
    if key == "cell_type_annotation":
        s = run.stage_result("annotate").summary
        clusters = s["clusters"]
        unresolved = sum(1 for c in clusters if c["predicted_label"] == annotation_tools.UNRESOLVED)
        reviewed = sum(1 for c in clusters if str(c["cluster"]) in run.review)
        review_text = (f"You reviewed {reviewed} of {len(clusters)} clusters." if reviewed
                       else "These are predictions for you to review.")
        n_second = sum(1 for c in clusters if c.get("label_source") == "PanglaoDB" and s["database"]["name"] != "PanglaoDB")
        second = f" ({n_second} from PanglaoDB, where Census was Unresolved)" if n_second else ""
        return (f"Cell types predicted from {s['database']['name']} markers: {len(clusters) - unresolved} of "
                f"{len(clusters)} clusters confident{second}, {unresolved} Unresolved. {review_text}")
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
    if number > 1 and needs_confirmation(run):
        return STEPS[0][1]
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
        samples = f" · {len(s['samples'])} samples" if s.get("samples") else ""
        return f"{s['n_cells']:,} cells × {s['n_genes']:,} genes · raw counts{samples}"
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
            if entry["stage"] == "metadata":
                return "The species or tissue changed on Data Overview, so run this step again."
            redone = next(title for _, title, st in STEPS if entry["stage"] in st)
            return f"{redone} was redone, so run this step again."
    return "An earlier step was redone, so run this step again."


def settings_record(run: AnalysisRun) -> list[dict[str, str]]:
    """The current settings that most affect results, in plain words."""
    p = run.params
    f = p.filter
    batch = ([{"label": "Batch correction",
               "value": f"Harmony by {p.integration.batch_key}" if p.integration.method == "harmony" else "off"}]
             if run.samples else [])  # only with several samples
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
        *batch,
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
        return f"{s['n_cells']:,} cells × {s['n_genes']:,} genes" + (f" · {len(s['samples'])} samples" if s.get("samples") else "")
    if stage == "qc_metrics":
        text = f"{s['n_mito_genes']} mitochondrial genes found"
        return text + (" · % mitochondrial filter switched off" if s.get("mito_filter_switched_off") else "")
    if stage == "filter":
        return f"{s['n_cells_kept']:,} cells × {s['n_genes_kept']:,} genes kept"
    if stage == "preprocess":
        return f"{s['n_hvg']:,} highly variable genes · regression {'on' if s['regressed'] else 'off'}"
    if stage == "reduce_dimensions":
        bend = f" · curve bends at PC {s['elbow_pc']}" if s.get("elbow_pc") else ""
        return f"{s['n_comps']} PCs · {s['n_connections']:,} graph connections{bend}"
    if stage == "cluster":
        return f"{s['n_clusters']} clusters"
    if stage == "markers":
        return f"{s['n_genes_tested']:,} genes ranked per cluster"
    if stage == "annotate":
        names = " + ".join(r["name"] for r in s.get("references") or [s["database"]])
        return f"{sum(s['label_counts'].values()):,} cells annotated · {names}"
    if stage == "metadata":
        return f"Species and tissue confirmed: {metadata_text(s)}"
    if stage == "review":
        did = {"confirmed": f"confirmed {s['label']}", "renamed": f"renamed to {s['label']}",
               "cleared": "review cleared"}[s["decision"]]
        return (f"Cluster {s['cluster']}: {did} · {s['n_confirmed']} confirmed, {s['n_renamed']} renamed, "
                f"{s['n_not_reviewed']} not reviewed")
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
    max_genes_off and max_mt_off. Several samples: `per_sample` True and `samples`, one entry per
    sample with the same boxes (sample, min_genes, max_genes, max_genes_off, max_mt, max_mt_off).
    Missing values keep `current`. Raises `InputError`.
    """
    if not values:
        return current

    def limits(box: dict[str, Any], where: str = "") -> tuple[int, int | None, float | None]:
        return (_number(box, "min_genes", f"Min genes per cell{where}", whole=True),
                None if box.get("max_genes_off") else _number(box, "max_genes", f"Max genes per cell{where}", whole=True),
                None if box.get("max_mt_off") else _number(box, "max_mt", f"Max % mitochondrial{where}", whole=False))

    min_genes, max_genes, max_mt = limits(values)
    sample_limits = ()
    if values.get("per_sample"):
        boxes = values.get("samples") or []
        if not isinstance(boxes, (list, tuple)) or not all(isinstance(b, dict) for b in boxes):
            raise InputError("The limits per sample are missing.")
        sample_limits = tuple(SampleLimits(str(b.get("sample")), *limits(b, f" ({b.get('sample')})")) for b in boxes)
    return FilterParams(
        min_genes=min_genes,
        min_cells=_number(values, "min_cells", "Min cells per gene", whole=True),
        max_n_genes_by_counts=max_genes,
        max_pct_counts_mt=max_mt,
        sample_limits=sample_limits,
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


def qc_settings(run: AnalysisRun) -> dict[str, Any]:
    """The "Mitochondrial genes" setting: the prefix to show in the box, the known patterns,
    and the pattern suggested by the data. The box shows the prefix QC last ran with; before
    QC has run, the data's suggestion (or the current setting if nothing matches)."""
    status = run.status()
    var_names = run.result("load").var_names if status["load"] == "done" else []
    suggested = qc_tools.suggest_mito_prefix(var_names)
    ran = status["qc_metrics"] != "not run"  # done, or needs re-running
    return {
        "mito_prefix": run.params.qc.mito_prefix if ran or suggested is None else suggested,
        "suggested": suggested,
        "patterns": [{"prefix": prefix, "species": species} for prefix, species in MITO_PATTERNS.items()],
    }


def qc_params_from_browser(values: dict[str, Any] | None, current: AnalysisParams) -> tuple[QCParams]:
    """The QC settings sent with Run QC -> backend settings. Raises `InputError`."""
    if not values or "mito_prefix" not in values:
        return (current.qc,)
    prefix = str(values["mito_prefix"]).strip()
    if not prefix:
        raise InputError("Enter the start of the mitochondrial gene names, e.g. 'MT-' or 'mt-'.")
    return (replace(current.qc, mito_prefix=prefix),)


def mito_check_payload(run: AnalysisRun, prefix: str | None) -> dict[str, Any] | None:
    """Live count under the "Mitochondrial genes" box: how many genes start with `prefix`."""
    if run.status()["load"] != "done" or prefix is None:
        return None
    genes = qc_tools.genes_with_prefix(run.result("load").var_names, str(prefix).strip())
    return {"prefix": prefix, "n_genes": len(genes), "examples": list(genes[:6])}


def qc_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Quality Control page shows, except the live filter preview."""
    status = run.status()
    payload: dict[str, Any] = {
        "settings": qc_settings(run),
        "metrics_done": status["qc_metrics"] == "done",
        "filter_status": {"done": "done", "needs re-running": "needs_rerun"}.get(status["filter"], "not_run"),
        "can_run": status["load"] == "done" and not needs_confirmation(run),
        "filter": {"current": asdict(run.params.filter), **suggestion_payload(run, None)},
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
    samples = qc_tools.sample_names(run.result("qc_metrics"))  # empty for one dataset
    payload.update({
        "n_cells": int(obs.shape[0]),
        "n_genes": int(run.result("qc_metrics").n_vars),
        "mito": {"genes": list(summary.mito_genes), "available": summary.mito_qc_available,
                 "prefix": summary.mito_prefix, "warnings": list(summary.warnings),
                 "filter_switched_off": run.stage_result("qc_metrics").summary.get("mito_filter_switched_off", False)},
        "stats": [
            {"statistic": label, **{m: float(stats[name]) for m, stats in summary.metric_stats.items()}}
            for name, label in QC_STATISTICS
        ],
        "cells": {
            "barcode": obs.index.astype(str).tolist(),
            "n_genes": obs["n_genes_by_counts"].astype(int).tolist(),
            "total_counts": [float(v) for v in obs["total_counts"]],
            "pct_mt": [float(v) for v in obs["pct_counts_mt"]],
            # Several samples: each cell's sample, as a position in `samples` (None for one dataset).
            "sample": ([samples.index(str(s)) for s in obs[load_tools.SAMPLE_KEY]] if samples else None),
        },
        "samples": samples,
    })
    return payload


def suggest_params_from_browser(values: dict[str, Any] | None, current: SuggestionParams | None) -> SuggestionParams:
    """The "MADs" boxes under More settings -> suggestion settings. Raises `InputError`."""
    base = current or SuggestionParams()
    if not values:
        return base
    return replace(
        base,
        n_mads_genes=_number(values, "n_mads_genes", "MADs for genes per cell", whole=False, above=True),
        n_mads_mito=_number(values, "n_mads_mito", "MADs for % mitochondrial", whole=False, above=True),
        per_sample=bool(values.get("per_sample", base.per_sample)),
    )


def suggestion_payload(run: AnalysisRun, values: dict[str, Any] | None) -> dict[str, Any]:
    """The thresholds suggested from this data (for the MADs in `values`, or the current ones),
    with how each was found. Before QC has run there is no suggestion yet."""
    try:
        sp = suggest_params_from_browser(values, run.params.suggest)
    except InputError as problem:
        return {"suggested": None, "details": None, "suggest_settings": None, "problem": str(problem)}
    settings = {"n_mads_genes": sp.n_mads_genes, "n_mads_mito": sp.n_mads_mito, "min_cells": sp.min_cells,
                "per_sample": sp.per_sample}
    if run.status()["qc_metrics"] == "not run":
        return {"suggested": None, "details": None, "suggest_settings": settings, "problem": None}
    suggestion = run.suggest_filter(sp)
    return {"suggested": asdict(suggestion.filter), "details": suggestion.details,
            "suggest_settings": settings, "problem": None}


def preview_payload(run: AnalysisRun, values: dict[str, Any] | None,
                    suggest_values: dict[str, Any] | None = None) -> dict[str, Any] | None:
    """What the filter boxes would keep, computed without changing anything (None before QC).
    Compared with what the suggested thresholds (for the MADs in `suggest_values`) would keep."""
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
    # The same preview at the suggested thresholds, for "vs suggested" on the page.
    suggestion = suggestion_payload(run, suggest_values)["suggested"]
    suggested_params = filter_from_dict(suggestion) if suggestion else params
    try:
        suggested = run.preview_filter(suggested_params)
    except Exception:  # e.g. the suggestion itself would remove everything: compare with nothing
        suggested_params, suggested = params, preview
    return {
        "ok": True,
        "settings": asdict(params),
        "n_cells_before": preview.n_cells_before,
        "n_cells_kept": preview.n_cells_kept,
        "n_genes_before": preview.n_genes_before,
        "n_genes_kept": preview.n_genes_kept,
        "suggested": {"n_cells_kept": suggested.n_cells_kept, "n_genes_kept": suggested.n_genes_kept},
        "filters": filter_table(params, preview, suggested_params, suggested),
        "removed_steps": [
            {"what": REMOVAL_STEPS[key][0], "reason": REMOVAL_STEPS[key][1], "removed": removed}
            for key, removed in preview.removed_in_order.items()
        ],
        # The cells these settings remove, by barcode, each with the filters it fails.
        # Several samples: cells before and kept per sample, with the suggestion's count.
        "by_sample": [{"sample": name, **counts, "suggested_kept": suggested.by_sample.get(name, {}).get("n_cells_kept")}
                      for name, counts in preview.by_sample.items()],
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
    """Run one stage with the current settings. QC that has never run uses the mitochondrial
    prefix the Quality Control page would show (the one suggested by the data)."""
    if stage == "qc_metrics":
        run.compute_qc_metrics(replace(run.params.qc, mito_prefix=qc_settings(run)["mito_prefix"]))
        return
    getattr(run, STAGE_METHODS[stage])()


def choice_needed(run: AnalysisRun, stage: str) -> str | None:
    """What the user must choose before `stage` can run, or None. "Run the remaining steps"
    stops there (a planned stop, not a failure): Cell Type Annotation has no default tissue."""
    markers = run.params.annotation.markers
    if stage == "annotate" and markers.species is not None and markers.species not in marker_db_tools.SPECIES:
        return ("Cell types can't be predicted for this species (the marker database covers human and mouse), so "
                "Cell Type Annotation was left out. Review the clusters with their marker genes.")
    if stage == "annotate" and not markers.tissues:
        return "Confirm the species and tissue on Data Overview, then run Cell Type Annotation."
    return None


def steps_done_before(job: dict[str, Any], stage: str) -> int:
    """How many of the job's steps come before the step that `stage` belongs to."""
    steps = list(dict.fromkeys(step_of_stage(s) for s in job["stages"]))
    return steps.index(step_of_stage(stage))


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
        "status": job["status"],  # "running", "done", "stopped" (a choice is needed) or "failed"
        "current": current[1] if current else None,
        "index": steps.index(current) + 1 if current else None,
        "total": len(steps),
        "failed_step": job.get("failed_step"),
        "stopped_step": job.get("stopped_step"),  # the step waiting for the user's choice
    }


def preprocess_settings(params: AnalysisParams) -> dict[str, Any]:
    """The Preprocessing settings shown on screen."""
    return {
        "target_sum": params.normalize.target_sum,
        # How highly variable genes are picked: "cutoffs" (Scanpy's default) or "top" (N genes).
        "hvg_method": "cutoffs" if params.hvg.n_top_genes is None else "top",
        "n_top_genes": params.hvg.n_top_genes if params.hvg.n_top_genes is not None else DEFAULT_N_TOP_GENES,
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
    # "Median of your cells" ticked = None (computed from the data).
    target_sum = (None if values.get("target_median")
                  else _number(values, "target_sum", "Counts per cell after normalization", whole=False, above=True))
    if values.get("hvg_method") == "top":
        hvg = replace(current.hvg, n_top_genes=_number(values, "n_top_genes", "Number of genes", whole=True, minimum=1))
    else:
        min_mean = _number(values, "min_mean", "Min mean expression", whole=False)
        max_mean = _number(values, "max_mean", "Max mean expression", whole=False, above=True)
        if min_mean >= max_mean:
            raise InputError("'Min mean expression' must be smaller than 'Max mean expression'.")
        min_disp = _number(values, "min_disp", "Min normalized dispersion", whole=False, minimum=-np.inf)
        hvg = replace(current.hvg, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, n_top_genes=None)
    max_value = _number(values, "max_value", "Clip scaled values at", whole=False, above=True)
    return (
        replace(current.normalize, target_sum=target_sum),
        hvg,
        replace(current.regress_out, enabled=bool(values.get("regress"))),
        replace(current.scale, max_value=max_value),
    )


def preprocess_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Preprocessing card shows: settings, and results once it has run."""
    payload: dict[str, Any] = {
        "status": _stage_status(run, "preprocess"),
        "settings": {"current": preprocess_settings(run.params), "reference": preprocess_settings(DEFAULT_PARAMS)},
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
        "target_sum_used": result.summary.get("target_sum"),  # the counts per cell actually used (exact)
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
    target = s.get("target_sum")
    norm = (f"normalized to {_whole(target):,} counts per cell" + (" (the median)" if s.get("target_is_median") else "")
            if target is not None else "normalized")
    how = "the most variable" if p["hvg"].get("n_top_genes") is not None else "by the cut-offs"
    return (f"Preprocessing done: {norm}; {s['n_hvg']:,} of {s['raw_shape'][1]:,} genes selected as highly "
            f"variable ({how}); {regression}; values scaled and capped at ±{p['scale']['max_value']:g}.")


def reduction_settings(params: AnalysisParams) -> dict[str, Any]:
    """The Dimensionality Reduction settings shown on screen (one seed for PCA, neighbors and UMAP)."""
    return {
        "n_comps": params.pca.n_comps,
        "n_neighbors": params.neighbors.n_neighbors,
        "n_pcs": params.neighbors.n_pcs,
        "seed": params.pca.random_state,
        "batch_correction": params.integration.method == "harmony",  # off by default (analyst, 2026-11-10)
        "batch_key": params.integration.batch_key,
    }


def reduction_params_from_browser(
    values: dict[str, Any] | None, current: AnalysisParams
) -> tuple[PCAParams, NeighborsParams, UMAPParams, IntegrationParams]:
    """The Dimensionality Reduction boxes -> backend settings (batch correction too, when sent).
    Raises `InputError`."""
    if not values:
        return current.pca, current.neighbors, current.umap, current.integration
    n_comps = _number(values, "n_comps", "Principal components", whole=True, minimum=1)
    n_neighbors = _number(values, "n_neighbors", "Neighbors per cell", whole=True, minimum=2)
    n_pcs = _number(values, "n_pcs", "PCs used for the neighbor graph", whole=True, minimum=1)
    if n_pcs > n_comps:
        raise InputError("'PCs used for the neighbor graph' can't be more than 'Principal components'.")
    seed = _number(values, "seed", "Random seed", whole=True)
    integ = current.integration
    if "batch_correction" in values:
        on = values.get("batch_correction")
        if not isinstance(on, bool):
            raise InputError("Choose whether to correct batch effects.")
        batch_key = values.get("batch_key", integ.batch_key)
        if not isinstance(batch_key, str) or not batch_key:
            raise InputError("Choose the column that holds the batches.")
        integ = replace(integ, method="harmony" if on else "none", batch_key=batch_key, random_state=seed)
    return (
        replace(current.pca, n_comps=n_comps, random_state=seed),
        replace(current.neighbors, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed),
        replace(current.umap, random_state=seed),
        integ,
    )


def batch_info(run: AnalysisRun) -> dict[str, Any] | None:
    """Batch correction choices for Dimensionality Reduction: only with several samples (None
    otherwise). The columns that could hold batches (sample first, then condition), and a warning
    for each column that can't be told apart from the condition."""
    if not run.samples or run.status()["preprocess"] != "done":
        return None
    adata = run.result("preprocess")
    columns = integration_tools.batch_column_choices(adata)
    order = {load_tools.SAMPLE_KEY: 0, load_tools.CONDITION_KEY: 1}
    columns.sort(key=lambda c: order.get(c, 2))
    labels = {load_tools.SAMPLE_KEY: "Sample", load_tools.CONDITION_KEY: "Condition"}
    return {
        "choices": [{"value": c, "label": labels.get(c, c), "n_values": int(adata.obs[c].astype(str).nunique()),
                     "warning": integration_tools.condition_warning(adata, c, load_tools.CONDITION_KEY)}
                    for c in columns],
        "low_mixing": LOW_MIXING,
    }


def reduction_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Dimensionality Reduction card shows."""
    payload: dict[str, Any] = {
        "status": _stage_status(run, "reduce_dimensions"),
        "settings": {"current": reduction_settings(run.params), "reference": reduction_settings(DEFAULT_PARAMS)},
        "batch": batch_info(run),
        "will_clear": _will_clear(run, "reduce_dimensions"),
        "results": None,
    }
    if run.status()["reduce_dimensions"] != "done":
        return payload
    result = run.stage_result("reduce_dimensions")
    adata, (pcs, graph, _, _) = result.adata, result.details
    obs, umap = adata.obs, adata.obsm["X_umap"]
    payload["results"] = {
        "variance_ratio": list(pcs.variance_ratio),
        "total_variance_ratio": pcs.total_variance_ratio,
        "n_pcs_used": graph.n_pcs,
        "elbow_pc": pcs.elbow_pc,  # where the variance curve bends (guidance only)
        # With several samples: correction used, mixing before/after, each cell's sample (for coloring).
        "integration": mixing_results(adata, result.summary.get("integration")),
        "umap": {"x": _floats(umap[:, 0]), "y": _floats(umap[:, 1])},
        "cells": {
            "barcode": obs.index.astype(str).tolist(),
            "n_genes": obs["n_genes_by_counts"].astype(int).tolist(),
            "total_counts": _floats(obs["total_counts"]),
            "pct_mt": _floats(obs["pct_counts_mt"]) if "pct_counts_mt" in obs else None,
        },
    }
    return payload


def mixing_results(adata, summary: dict[str, Any] | None) -> dict[str, Any] | None:
    """The batch part of the Dimensionality Reduction results (None for a single dataset)."""
    if not summary:
        return None
    labels = adata.obs[summary["batch_key"]].astype(str)
    names = sorted(labels.unique())
    index = {n: i for i, n in enumerate(names)}
    return {**summary, "batches": names, "cell_batch": [index[x] for x in labels],
            "suggest": summary["method"] == "none" and summary["mixing_before"] < LOW_MIXING}


def _reduction_confirmation(result) -> str:
    s, p = result.summary, result.settings
    n_cells = result.adata.n_obs
    return (f"Dimensionality Reduction done: {s['n_comps']} principal components explain "
            f"{100 * s['total_variance_ratio']:.1f}% of the variance; neighbor graph with "
            f"{s['n_neighbors']} neighbors per cell from the first {p['neighbors']['n_pcs']} PCs; "
            f"UMAP drawn for {n_cells:,} cells." + _batch_text(s.get("integration")))


def _batch_text(integ: dict[str, Any] | None) -> str:
    """" Batch correction: Harmony by sample; sample mixing 0.65 -> 0.82." or the mixing alone."""
    if not integ:
        return ""
    if integ["method"] == "harmony":
        return (f" Batch correction: Harmony by {integ['batch_key']}; sample mixing {integ['mixing_before']:.2f} → "
                f"{integ['mixing_after']:.2f}.")
    return f" No batch correction; sample mixing {integ['mixing_before']:.2f}."


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
        "settings": {"current": cluster_settings(run.params), "reference": cluster_settings(DEFAULT_PARAMS)},
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
        "warnings": list(summary.warnings),  # shown above the results (one cluster, tiny clusters)
        # With several samples: per cluster, how well its samples mix (flagged below LOW_MIXING).
        "mixing": cluster_mixing(run, result),
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


def cluster_mixing(run: AnalysisRun, result) -> dict[str, Any] | None:
    """Clustering's sample mixing table and its note (None for a single dataset)."""
    rows = result.summary.get("mixing")
    if not rows:
        return None
    integ = run.stage_result("reduce_dimensions").summary.get("integration") or {}
    low = [r for r in rows if r["low"]]
    corrected = integ.get("method") == "harmony"
    note = None
    if low:
        n = len(low)
        what = f"{n} {'cluster is' if n == 1 else 'clusters are'} made mostly of one or a few samples"
        note = (f"{what}. That can be a batch effect, or a cell state found only in some samples (for example after "
                "a treatment). " + ("Batch correction is on; check these clusters' marker genes."
                                    if corrected else "Try batch correction on Dimensionality Reduction and compare."))
    return {"rows": rows, "low_mixing": LOW_MIXING, "batch_key": integ.get("batch_key"), "corrected": corrected,
            "note": note}


# --- Marker Genes page ------------------------------------------------------------------------

TOP_MARKERS_SENT = 50  # top genes per cluster sent to the page (the page shows fewer by default)
DOTPLOT_GENES_PER_CLUSTER = 3
HEATMAP_GENE_OPTIONS = (3, 5, 10)  # top genes per cluster the heatmap offers
HEATMAP_GENES_DEFAULT = 5  # the app's own choice (analyst, 2026-10-29): fits one screen for 8-16 clusters


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
        "settings": {"current": marker_settings(run.params), "reference": marker_settings(DEFAULT_PARAMS)},
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
    dot_genes = [gene for _, gene in marker_tools.top_genes_per_group(adata, DOTPLOT_GENES_PER_CLUSTER)]
    expr = marker_tools.expression_by_group(adata, dot_genes, summary.groupby, use_raw=summary.use_raw)
    ranks = table[table["gene"].isin(dot_genes)].set_index(["group", "gene"])
    dots = [{"group": r.group, "gene": r.gene, "mean_expr": r.mean_expr, "frac_expressing": r.frac_expressing,
             "rank": int(ranks.loc[(r.group, r.gene), "rank"]),
             "logfoldchange": _num(ranks.loc[(r.group, r.gene), "logfoldchange"])}
            for r in expr.itertuples(index=False)]
    # Heatmap: each cluster's top genes (listed once) for every option, and the scaled cluster
    # means of all of them, so the page can switch between 3, 5 and 10 genes without asking.
    heat_rows = {str(n): [{"cluster": c, "gene": g} for c, g in marker_tools.top_genes_per_group(adata, n)]
                 for n in HEATMAP_GENE_OPTIONS}
    heat_genes = list(dict.fromkeys(r["gene"] for rows_n in heat_rows.values() for r in rows_n))
    heat = marker_tools.scaled_group_means(adata, heat_genes, summary.groupby, use_raw=summary.use_raw)
    payload["results"] = {
        "settings": marker_settings(run.params),
        "method_label": dict(MARKER_METHODS).get(summary.method, summary.method),
        "n_genes_tested": summary.n_genes_tested,
        "groups": [{"cluster": g, "n_cells": int(sizes.get(g, 0))} for g in groups],
        "top": top,
        "dotplot": {"genes": dot_genes, "values": dots},
        "heatmap": {
            "options": list(HEATMAP_GENE_OPTIONS),
            "default": HEATMAP_GENES_DEFAULT,
            "rows": heat_rows,
            "values": [{"group": r.group, "gene": r.gene, "mean_expr": r.mean_expr, "scaled": r.scaled}
                       for r in heat.itertuples(index=False)],
        },
        "notes": list(summary.notes),
    }
    return payload


# --- Cell Type Annotation page ------------------------------------------------------------------

CANDIDATES_SHOWN = 5  # candidate cell types listed in a cluster's evidence
DOTPLOT_MARKERS_PER_CLUSTER = 3  # matching markers of each cluster's best candidate in the dot plot


def annotation_settings(params: AnalysisParams) -> dict[str, Any]:
    """The Cell Type Annotation settings: species, the Census tissue and PanglaoDB tissues (which
    cell types are compared), whether Census is used, and the test's checks."""
    m, t, c = params.annotation.markers, params.annotation.test, params.annotation.census
    return {"species": m.species, "tissues": list(m.tissues), "census_tissue": c.tissue, "use_census": c.use,
            "celltypist_model": params.annotation.celltypist.model,
            "top_n_genes": t.top_n_genes, "max_adj_p": t.max_adj_p, "min_overlap": t.min_overlap,
            "ambiguity_fold": t.ambiguity_fold}


def _annotation_checks(values: dict[str, Any]) -> dict[str, Any]:
    """The four checks typed on the page, checked. Raises `InputError`."""
    max_adj_p = _number(values, "max_adj_p", "Max adjusted p-value", whole=False, above=True)
    if max_adj_p > 1:
        raise InputError("'Max adjusted p-value' can't be more than 1.")
    return {"top_n_genes": _number(values, "top_n_genes", "Top marker genes per cluster", whole=True, minimum=1),
            "max_adj_p": max_adj_p,
            "min_overlap": _number(values, "min_overlap", "Min matching markers", whole=True, minimum=1),
            "ambiguity_fold": _number(values, "ambiguity_fold", "Lead over the runner-up", whole=False, minimum=1)}


def annotation_params_from_browser(values: dict[str, Any] | None, current: AnalysisParams) -> tuple[AnnotationParams]:
    """The page's checks and reference choice ("PanglaoDB only" = `use_census` False) -> backend
    settings. Species and tissues aren't set here: they come from the dataset's metadata (Data
    Overview). Raises `InputError`."""
    if not values:
        return (current.annotation,)
    a = current.annotation
    use_census = values.get("use_census", a.census.use)
    if not isinstance(use_census, bool):
        raise InputError("Choose which references to use.")
    model = values.get("celltypist_model", a.celltypist.model) or None
    if model is not None:
        offered = {c["model"] for c in celltypist_tools.model_choices(a.markers.species)}
        if model not in offered:
            raise InputError("Choose one of the CellTypist models offered, or none.")
    return (replace(a, test=replace(a.test, **_annotation_checks(values)), census=replace(a.census, use=use_census),
                    celltypist=replace(a.celltypist, model=model)),)


# --- Dataset metadata (Data Overview) ---------------------------------------------------------

# Species offered on Data Overview. Cell types can be predicted for those in the marker database
# (marker_db.SPECIES); any other species still runs every other step.
SPECIES_LABELS = {"human": "Human", "mouse": "Mouse", "other": "Other species"}


def metadata_text(m: dict[str, Any]) -> str:
    """"Human · Blood and Immune system", in words."""
    species = SPECIES_LABELS.get(m.get("species"), m.get("species") or "species not set")
    if m.get("tissue"):
        return f"{species} · {tissue_label(m['tissue'])}"
    return f"{species} · {_join(list(m['tissues']))}" if m.get("tissues") else species


def tissue_label(tissue: str | None) -> str | None:
    """A Census tissue name for the page: "blood" -> "Blood", "skin of body" -> "Skin of body"."""
    return tissue[:1].upper() + tissue[1:] if tissue else tissue


def samples_payload(run: AnalysisRun) -> dict[str, Any]:
    """Data Overview's samples: the loaded samples (empty for one dataset), whether they came from
    several files, the column choice for a single file, and the columns that could name samples."""
    from_files = bool(run.sample_files)
    columns = run.sample_columns or {}
    choices = []
    if not from_files:
        loaded = run.result("load")
        added = {load_tools.SAMPLE_KEY, load_tools.CONDITION_KEY} if columns else set()
        choices = [c for c in load_tools.sample_column_choices(loaded) if c["column"] not in added]
    return {
        "list": [{"name": s["name"], "condition": s.get("condition"), "n_cells": s["n_cells"], "file": s.get("file")}
                 for s in run.samples],
        "from_files": from_files,
        "sample_column": columns.get("sample"),
        "condition_column": columns.get("condition"),
        "choices": choices,  # likely sample columns first; each has "hint": "sample" | "condition" | None
        "max_samples": MAX_SAMPLES,
    }


def sample_files_from_browser(values: dict[str, Any] | None, files: list[dict] | None) -> list:
    """The uploaded files (Shiny's list, in the order chosen) with the names and conditions typed
    for them -> sample files for loading. Raises `InputError`."""
    details = (values or {}).get("samples") or []
    if not files or len(files) != len(details):
        raise InputError("The sample files didn't arrive. Choose them again.")
    out = []
    for file, d in zip(files, details):
        if not isinstance(d, dict):
            raise InputError("Each sample needs a name.")
        condition = d.get("condition")
        out.append(load_tools.SampleFile(path=file["datapath"], name=str(d.get("name") or ""),
                                         file_name=file["name"], condition=str(condition) if condition else None))
    return out


def metadata_from_browser(values: dict[str, Any] | None) -> tuple[str, str | None, list[str], str | None]:
    """The Confirm section's values -> (species, Census tissue, PanglaoDB tissues, about). The
    PanglaoDB tissues may be empty (then they follow the Census tissue). Raises `InputError`."""
    values = values or {}
    species = values.get("species")
    if species not in SPECIES_LABELS:
        raise InputError("Choose the species of your sample.")
    tissue = values.get("tissue") or None
    if tissue is not None and not isinstance(tissue, str):
        raise InputError("Choose the tissue.")
    tissues = values.get("tissues") or []
    if not isinstance(tissues, (list, tuple)) or not all(isinstance(x, str) for x in tissues):
        raise InputError("The tissues must be a list of names.")
    if species in marker_db_tools.SPECIES and not (tissue or tissues):
        raise InputError("Choose the tissue.")
    about = values.get("about")
    if about is not None and not isinstance(about, str):
        raise InputError("The description must be text.")
    return species, tissue, list(tissues), about


def _species_label(species: str) -> str:
    """"Human", "Mouse", or the species' name for those the app files under "Other species"."""
    return SPECIES_LABELS.get(species, species.capitalize())


def _species_suggestion(clues: list[dict[str, str]]) -> tuple[str | None, str | None]:
    """The species to suggest (the strongest clue first: Ensembl IDs, then the name style, then
    the mitochondrial prefix) and why, in words. Species the app doesn't list become "other"."""
    if not clues:
        return None, None
    species = clues[0]["species"]
    agreeing = [c["clue"] for c in clues if c["species"] == species]
    return (species if species in SPECIES_LABELS else "other"), f"{_species_label(species)}: {_join(agreeing)}"


def _species_mismatch(species: str | None, clues: list[dict[str, str]]) -> str | None:
    """A warning when the chosen species doesn't fit the gene names (Ensembl ID prefix, how symbols
    are written, mitochondrial prefix). Only a warning; it never blocks."""
    if not species or not clues:
        return None
    if species == "other":
        others = [c for c in clues if c["species"] in marker_db_tools.SPECIES]
        if others and len({c["species"] for c in others}) == 1:
            return (f"The gene names look like {_species_label(others[0]['species']).lower()} data "
                    f"({_join([c['clue'] for c in others])}). Choose it to get cell types predicted.")
        return None
    disagree = [c for c in clues if c["species"] != species]
    if not disagree:
        return None
    reasons = [f"{c['clue']} (as in {_species_label(c['species']).lower()} data)" for c in disagree]
    return f"Check the species: you chose {SPECIES_LABELS[species]}, but {_join(reasons)}."


def metadata_payload(run: AnalysisRun) -> dict[str, Any]:
    """The Confirm section on Data Overview: what is set, the choices, and the suggestion."""
    loaded = run.result("load")
    mito_matches = qc_tools.mito_prefix_matches(loaded.var_names)
    top = mito_matches[0]["prefix"] if mito_matches else None
    ids = loaded.var["gene_ids"] if "gene_ids" in loaded.var else None
    clues = marker_db_tools.species_clues(loaded.var_names, ids, top)
    suggested, reason = _species_suggestion(clues)
    m = run.metadata
    db = _database_info(run)
    return {
        "confirmed": run.metadata_confirmed,
        "confirmed_by": m.get("confirmed_by"),  # "dataset" (e.g. demo data), "user" or None
        "species": m.get("species"),
        "tissue": m.get("tissue"),  # the Census tissue
        "tissues": list(m.get("tissues") or []),  # PanglaoDB's tissues
        "about": m.get("about"),
        "text": metadata_text(m) if m else None,
        # species -> Census tissues (most cells first): name, label, cells, whether its marker table is
        # built, and the matching PanglaoDB tissues. Empty when Census's tissue list isn't saved; the
        # page then offers PanglaoDB's tissues only.
        "census_tissues": census_tissue_options(run),
        "census_version": run.params.annotation.census.release,
        "species_options": [{"value": v, "label": label} for v, label in SPECIES_LABELS.items()],
        "marker_species": list(marker_db_tools.SPECIES),  # species with cell types (tissue required)
        "suggested_species": suggested,
        "suggestion_reason": reason,  # e.g. "Mouse: Ensembl IDs start with ENSMUSG and gene names …"
        "species_clues": clues,
        "mito_prefix": top,
        "tissues_available": db["available"],
        "tissue_problem": db["problem"],
        # species -> tissue -> number of cell types (tissues with none for a species are left out)
        "tissue_cell_types": {sp: {tissue: len(types) for tissue, types in by_tissue.items() if types}
                              for sp, by_tissue in db["cell_types"].items()},
        "mismatch": _species_mismatch(m.get("species"), clues),
        "max_about_length": MAX_ABOUT_LENGTH,
    }


def census_tissue_options(run: AnalysisRun) -> dict[str, list[dict[str, Any]]]:
    """Per species with cell types: Census's tissues for the tissue picker (see `metadata_payload`)."""
    release = run.params.annotation.census.release
    return {species: [{"tissue": t["tissue"], "label": tissue_label(t["tissue"]), "n_cells": t["n_cells"],
                       "has_table": t["has_table"], "panglaodb": list(census_tools.panglaodb_tissues_for(t["tissue"]))}
                      for t in census_tools.census_tissues(species, release, run.census_dir)]
            for species in marker_db_tools.SPECIES}


def _database_info(run: AnalysisRun) -> dict[str, Any]:
    """The marker database: name, version, tissues and, per species, each tissue's cell types
    (those whose markers pass the selection rule), or why the database can't be read."""
    try:
        db = marker_db_tools.load_panglaodb(run.marker_db_path)
    except marker_db_tools.MarkerDatabaseError as problem:
        return {"available": False, "problem": str(problem), "name": None, "version": None, "tissues": [],
                "cell_types": {}}
    rule = run.params.annotation.markers
    return {"available": True, "problem": None, "name": db.name, "version": db.version,
            "tissues": list(marker_db_tools.tissues(db)),
            "cell_types": {species: {tissue: list(types) for tissue, types in
                                     marker_db_tools.cell_types_by_tissue(db, replace(rule, species=species)).items()}
                           for species in marker_db_tools.SPECIES}}


def annotation_payload(run: AnalysisRun) -> dict[str, Any]:
    """Everything the Cell Type Annotation page shows. Predictions only: reviewed labels are
    never part of the app."""
    status = run.status()
    m = run.params.annotation.markers
    payload: dict[str, Any] = {
        "status": _stage_status(run, "annotate"),
        "settings": {"current": annotation_settings(run.params), "reference": annotation_settings(DEFAULT_PARAMS)},
        "metadata": {"species": SPECIES_LABELS.get(m.species, m.species), "tissues": list(m.tissues),
                     "tissue": tissue_label(run.params.annotation.census.tissue),
                     "confirmed_by": run.metadata.get("confirmed_by"),
                     "has_markers": m.species in marker_db_tools.SPECIES},
        "database": _database_info(run),
        "census": census_reference_info(run),
        "celltypist": celltypist_info(run),
        "selection": {"canonical_only": m.canonical_only, "min_sensitivity": m.min_sensitivity,
                      "min_markers": m.min_markers},
        "unresolved_label": annotation_tools.UNRESOLVED,
        "will_clear": _will_clear(run, "annotate"),
        "review_cleared": run.review_cleared,  # why an earlier review was cleared, if it was
        "results": None,
    }
    if status["annotate"] != "done":
        return payload
    result = run.stage_result("annotate")
    adata = result.adata
    summary, match = result.details["summary"], result.details["match"]
    by_reference = result.details["references"]
    candidates = {name: _candidates(ref["summary"]) for name, ref in by_reference.items()}
    answers = {name: {str(r.cluster): {"predicted_label": r.predicted_label, "status": r.status,
                                       "best_candidate": r.best_candidate, "overlap": int(r.overlap),
                                       "n_markers": int(r.n_markers), "adj_p": float(r.adj_p),
                                       "second_candidate": r.second_candidate}
                      for r in ref["summary"].clusters.itertuples()}
               for name, ref in by_reference.items()}
    # Dot plot: each cluster's first matching markers of its best candidate, once each.
    dot_genes: list[str] = []
    for row in summary.clusters.itertuples():
        found = [g for g in str(row.matching_markers or "").split(", ") if g]
        for gene in [g for g in found if g not in dot_genes][:DOTPLOT_MARKERS_PER_CLUSTER]:
            dot_genes.append(gene)
    use_raw = run.params.markers.use_raw and adata.raw is not None
    expr = (marker_tools.expression_by_group(adata, dot_genes, summary.groupby, use_raw=use_raw)
            if dot_genes else None)
    owner = {g: row.best_candidate for row in summary.clusters.itertuples()
             for g in str(row.matching_markers or "").split(", ") if g in dot_genes}
    obs, umap = adata.obs, adata.obsm["X_umap"]
    columns = ["cluster", "n_cells", "predicted_label", "status", "reasons", "best_candidate", "overlap",
               "n_markers", "adj_p", "matching_markers", "second_candidate", "second_overlap", "second_adj_p",
               "label_source"]
    db = result.summary["database"]
    first = "Census" if db["name"] == "CELLxGENE Census" else "PanglaoDB"
    payload["results"] = {
        "settings": annotation_settings(run.params),
        "database": {"name": db["name"], "version": db["version"]},
        # Each reference used: name, version, tissues, role ("first" / "second"), citation.
        "references": [{"key": name, "name": ref["source"]["name"], "version": ref["source"]["version"],
                        "tissues": [tissue_label(t) if name == "Census" else t for t in ref["source"]["tissues"]],
                        "role": ref["source"]["role"], "citation": ref["source"].get("citation")}
                       for name, ref in by_reference.items()],
        "first_reference": first,
        # reference -> cluster -> its answer (label, status, best candidate, …): for "How was this identified?"
        "answers": answers,
        "n_cell_types_tested": result.summary["n_cell_types_tested"],
        "n_cell_types_not_tested": len(summary.skipped_cell_types),
        "n_markers_not_found": len({g for gs in match.missing.values() for g in gs}),
        "clusters": records(summary.clusters[columns]),
        # CellTypist's answer per cluster (Step A: shown beside the marker check), or why it didn't run.
        "celltypist": result.summary.get("celltypist"),
        "candidates": candidates[first],  # the first reference's best candidates per cluster
        "candidates_by_reference": candidates,
        "dotplot": {
            "genes": [{"gene": g, "cell_type": owner.get(g)} for g in dot_genes],
            "values": [] if expr is None else [
                {"group": r.group, "gene": r.gene, "mean_expr": r.mean_expr, "frac_expressing": r.frac_expressing}
                for r in expr.itertuples(index=False)],
        },
        "label_counts": summary.label_counts,
        "nothing_matched": nothing_matched_text(result.summary, run.params.annotation.markers),
        # Why some clusters are Unresolved: a note for the page, and one per Unresolved cluster.
        "unresolved_note": unresolved_note(result.summary),
        "guidance": {str(row.cluster): unresolved_guidance(row, result.summary, run.params.annotation.markers)
                     for row in summary.clusters.itertuples() if row.predicted_label == annotation_tools.UNRESOLVED},
        "review": _review_payload(run, summary.clusters),
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


def _candidates(summary) -> dict[str, list[dict[str, Any]]]:
    """Cluster -> its best candidates (most significant first) in one reference's test."""
    tests = summary.tests.assign(score=-np.log10(summary.tests["adj_p"].clip(lower=np.finfo(float).tiny)))
    return {str(c): records(rows.sort_values(["adj_p", "overlap", "cell_type"], ascending=[True, False, True])
                            .head(CANDIDATES_SHOWN)[["cell_type", "overlap", "n_markers", "adj_p", "score",
                                                     "matching_markers"]])
            for c, rows in tests.groupby("cluster", sort=False)}


def celltypist_info(run: AnalysisRun) -> dict[str, Any]:
    """CellTypist on the Cell Type Annotation page: the models offered for this species (whether each is
    downloaded), the one suggested from the species and tissue, and the one in use."""
    m = run.params.annotation.markers
    return {"choices": celltypist_tools.model_choices(m.species, run.celltypist_dir),
            "suggested": celltypist_tools.suggest_model(m.species, run.params.annotation.census.tissue),
            "current": run.params.annotation.celltypist.model,
            "low_gene_match": celltypist_tools.LOW_GENE_MATCH,
            "citation": celltypist_tools.CITATION}


def census_reference_info(run: AnalysisRun) -> dict[str, Any]:
    """Census for this dataset: its tissue, whether it is used, and whether its marker table exists
    (if not, the page offers "PanglaoDB only")."""
    c, species = run.params.annotation.census, run.params.annotation.markers.species
    has_table = bool(c.tissue and species in census_tools.CENSUS_SPECIES
                     and census_tools.table_path(species, c.tissue, c.release, run.census_dir).is_file())
    return {"tissue": tissue_label(c.tissue), "use": c.use, "release": c.release, "has_table": has_table,
            "available": bool(c.tissue)}


MOSTLY_UNRESOLVED = 0.5  # the app's own rule: above this share, the species/tissue may be the cause


def unresolved_note(summary: dict[str, Any]) -> str | None:
    """A note when some (not all) clusters are Unresolved. A wrong species or tissue shows as most
    clusters failing at once, so it is mentioned only then; otherwise the usual causes are listed."""
    if summary.get("nothing_matched"):
        return None  # the "No cell type matched" box says it
    clusters = summary["clusters"]
    n = sum(c["predicted_label"] == annotation_tools.UNRESOLVED for c in clusters)
    if n == 0:
        return None
    if n / len(clusters) > MOSTLY_UNRESOLVED:
        return (f"{n} of {len(clusters)} clusters are Unresolved. When most clusters fail, check that the species "
                "and tissue on Data Overview fit your sample; otherwise the marker database may not list these cell "
                "types for these tissues. Open the Clusters tab and click a cluster to see why it is Unresolved.")
    are = "cluster is" if n == 1 else "clusters are"
    return (f"{n} {are} Unresolved. This often means a cell type the marker database doesn't list for these "
            "tissues, a rare or mixed population, or two similar cell types. Open the Clusters tab and click a "
            "cluster to see why it is Unresolved, then review it yourself.")


def unresolved_guidance(row, summary: dict[str, Any], markers) -> str:
    """One Unresolved cluster's likely cause, from the first check it failed, plus its closest
    match among other tissues' cell types when there is a clear one (a hint to check)."""
    tissues = _join(list(markers.tissues))
    if row.status == "ambiguous":
        text = (f"Two cell types fit almost equally: {row.best_candidate} and {row.second_candidate}. They are often "
                "closely related; decide between them in your review, or name the type they share.")
    elif row.status == "not_significant":
        text = ("Some of its top genes match a cell type, but too few to rule out chance. It may be a rare or mixed "
                "population (for example doublets), or a cell type the database describes only weakly. Check its "
                "top marker genes.")
    else:
        text = (f"Too few of its top genes are markers of any cell type listed for {tissues}. It may be a cell type "
                "the database doesn't list there (a rare one, or one from another tissue), or lower-quality cells. "
                "Check its top marker genes.")
    elsewhere = (summary.get("best_anywhere") or {}).get(str(row.cluster))
    if elsewhere:
        where = _join(elsewhere["tissues"]) or "another tissue"
        text += (f" Closest match in another tissue: {elsewhere['cell_type']} ({where}), {elsewhere['overlap']} of "
                 f"{elsewhere['n_markers']} markers. A hint to check, not a prediction.")
    return text


def nothing_matched_text(summary: dict[str, Any], markers) -> str | None:
    """When every cluster is Unresolved: say so, and point to tissues that would fit better (the
    run's hint: tissues whose cell types matched clusters when all tissues were tried)."""
    if not summary.get("nothing_matched"):
        return None
    ticked = _join(list(markers.tissues))
    hint = summary.get("tissue_hint") or []
    if not hint:
        return (f"No cluster matched the cell types of the ticked tissues ({ticked}), and no other tissue in "
                f"{summary['database']['name']} fits clearly. Check the species, or review the clusters yourself.")
    closest = _join([f"{h['tissue']} ({h['n_clusters']} {'cluster' if h['n_clusters'] == 1 else 'clusters'})"
                     for h in hint])
    return (f"No cluster matched the cell types of the ticked tissues ({ticked}). Check that the tissues fit "
            f"your data: the closest matches come from {closest}.")


def _review_payload(run: AnalysisRun, clusters) -> dict[str, Any]:
    """The analyst's decision for each cluster, and how many are reviewed. A confirmed label
    whose prediction has changed since (annotation re-run with other settings) is flagged."""
    predicted = dict(zip(clusters["cluster"].astype(str), clusters["predicted_label"]))
    decisions = {
        c: {**asdict(r), "prediction_changed": r.decision == "confirmed" and r.label != predicted.get(c)}
        for c, r in run.review.items() if c in predicted
    }
    return {"clusters": decisions, "not_reviewed_label": review_tools.NOT_REVIEWED,
            **review_tools.review_summary(run.review, list(predicted))}


def review_map(run: AnalysisRun) -> dict[str, str] | None:
    """Cluster -> reviewed label (only reviewed clusters), or None when nothing is reviewed."""
    return {c: r.label for c, r in run.review.items()} or None


def review_params_from_browser(values: dict[str, Any] | None) -> tuple[str, str, str | None, str | None]:
    """A review decision sent by the page -> (cluster, decision, label, note). Raises `InputError`."""
    if not isinstance(values, dict) or not values.get("cluster"):
        raise InputError("Choose a cluster to review.")
    decision = values.get("decision")
    if decision not in ("confirm", "rename", "clear"):
        raise InputError("Choose to confirm the prediction, rename it, or clear the review.")
    return str(values["cluster"]), decision, values.get("label"), values.get("note")


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
        "reviewed": review_map(run),
        "not_reviewed_label": review_tools.NOT_REVIEWED,
        "samples": explore_samples(run, adata, summary.key_added),
        "unavailable": [],
    }
    if status["markers"] == "done":
        table = marker_tools.marker_table(run.result("markers"), n_top=EXPLORE_TOP_MARKERS)
        results["top_markers"] = {str(g): list(map(str, rows["gene"])) for g, rows in table.groupby("group", sort=False)}
    else:
        results["unavailable"].append("Marker Genes " + ("needs re-running" if status["markers"] == "needs re-running" else "hasn't been run"))
    if status["annotate"] == "done":
        result = run.stage_result("annotate")
        clusters = result.details["summary"].clusters
        results["predictions"] = {
            "unresolved_label": annotation_tools.UNRESOLVED,
            "clusters": json.loads(clusters[["cluster", "predicted_label", "status"]].to_json(orient="records")),
            "cell_label": result.adata.obs[annotation_tools.PREDICTION_KEY].astype(str).tolist(),
            # Cell types in the order of their first cluster.
            "label_order": list(dict.fromkeys(clusters["predicted_label"])),
        }
        if results["samples"]:  # composition by predicted cell type too
            results["samples"]["composition"]["cell_type"] = records(integration_tools.composition(
                result.adata, load_tools.SAMPLE_KEY, annotation_tools.PREDICTION_KEY))
    else:
        results["unavailable"].append("Cell Type Annotation " + ("needs re-running" if status["annotate"] == "needs re-running" else "hasn't been run"))
    return {"available": True, "reason": None, "results": results}


def explore_samples(run: AnalysisRun, adata, cluster_key: str) -> dict[str, Any] | None:
    """Explore Results with several samples: each cell's sample and condition (indices into the
    lists, for coloring) and the composition (cells per sample per cluster; per cell type is added
    when annotation is up to date). None for a single dataset."""
    if not run.samples or load_tools.SAMPLE_KEY not in adata.obs:
        return None
    listed = [s["name"] for s in run.samples]
    names = [n for n in listed if n in set(adata.obs[load_tools.SAMPLE_KEY].astype(str))]
    index = {n: i for i, n in enumerate(names)}
    condition_of = {s["name"]: s.get("condition") for s in run.samples}
    conditions = sorted({c for c in condition_of.values() if c}) or None
    out: dict[str, Any] = {
        "names": names,  # in the order they were loaded
        "condition_of": {n: condition_of.get(n) for n in names},
        "conditions": conditions,
        "cell_sample": [index[x] for x in adata.obs[load_tools.SAMPLE_KEY].astype(str)],
        "cell_condition": None,
        "composition": {"cluster": records(integration_tools.composition(adata, load_tools.SAMPLE_KEY, cluster_key))},
    }
    if conditions and load_tools.CONDITION_KEY in adata.obs:
        cindex = {c: i for i, c in enumerate(conditions)}
        out["cell_condition"] = [cindex[x] for x in adata.obs[load_tools.CONDITION_KEY].astype(str)]
    return out


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


def filter_table(params: FilterParams, preview, suggested_params: FilterParams, suggested) -> list[dict[str, Any]]:
    """One row per filter, in the order they run: your setting and what it removes, next to the
    suggested setting and what that removes. A filter switched off has setting None and
    removes None (it is not run)."""
    rows = []
    for key, (what, reason) in REMOVAL_STEPS.items():
        name = key.split(": ")[1]  # "cells: max_pct_counts_mt" -> "max_pct_counts_mt"
        setting, suggested_setting = getattr(params, name), getattr(suggested_params, name)
        if params.sample_limits and name != "min_cells":  # each sample has its own; the rows say so
            setting = [getattr(s, name) for s in params.sample_limits]
        if suggested_params.sample_limits and name != "min_cells":
            suggested_setting = [getattr(s, name) for s in suggested_params.sample_limits]
        rows.append({
            "key": name,
            "what": what,
            "reason": reason,
            "setting": setting,
            "suggested_setting": suggested_setting,
            "removed": preview.removed_in_order.get(key),
            "suggested_removed": suggested.removed_in_order.get(key),
            "differs": setting != suggested_setting,
        })
    return rows


def _whole(value: float) -> int:
    """Rounded to a whole number, halves up (as the page rounds), e.g. 2,242.5 -> 2,243."""
    return int(np.floor(value + 0.5))


def _decimal(value: float) -> str:
    """A setting that is a decimal, always shown with its decimal point (1.0, 0.5, 0.75)."""
    return repr(float(value))


def _join(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]
