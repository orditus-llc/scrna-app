"""The run report: one self-contained web page (HTML) describing an analysis, for the user to
keep, share or print to PDF from the browser.

It is built only from the run's record and step summaries (nothing is recomputed): the dataset,
each step's settings and result, the cell types and the analyst's review, a Methods paragraph
written from the settings actually used, package versions, seeds and the run log.
"""

from __future__ import annotations

import html
from datetime import datetime, timezone
from typing import Any

from scrna_analysis import integration
from scrna_analysis.params import CensusMarkerParams
from scrna_analysis.pipeline import AnalysisRun

import payloads

CSS = """
body { font-family: -apple-system, "Helvetica Neue", Arial, sans-serif; color: #1c2529; max-width: 860px;
       margin: 32px auto; padding: 0 20px; line-height: 1.5; font-size: 14px; }
h1 { font-size: 26px; margin: 0 0 4px; } h2 { font-size: 18px; border-bottom: 2px solid #1c2529;
     padding-bottom: 4px; margin-top: 30px; } .meta { color: #56656b; font-size: 13px; }
table { border-collapse: collapse; width: 100%; margin: 8px 0 14px; font-size: 13px; }
th, td { border: 1px solid #d3dcdf; padding: 5px 7px; text-align: left; vertical-align: top; }
th { background: #eef2f3; } td.num { text-align: right; font-variant-numeric: tabular-nums; }
.methods { background: #f4f6f7; border-left: 3px solid #0d9488; padding: 10px 14px; }
.note { color: #56656b; font-size: 12.5px; } code { font-family: Menlo, monospace; font-size: 12px; }
@media print { body { margin: 0; } h2 { break-after: avoid; } tr { break-inside: avoid; } }
"""


def _e(value: Any) -> str:
    return html.escape(str(value))


def _n(value: float | int) -> str:
    return f"{value:,}" if isinstance(value, int) else f"{value:,.4g}"


def methods_text(run: AnalysisRun) -> str:
    """A Methods paragraph written from the settings and results actually recorded (plain text)."""
    status, p, versions = run.status(), run.params, run.record()["versions"]
    summary = {e["stage"]: e["summary"] for e in run.history if e["stage"] != "review"}  # latest of each
    parts = [f"Single-cell RNA-seq data were analyzed with Scanpy {versions.get('scanpy')} "
             f"(AnnData {versions.get('anndata')}, Python {versions.get('python')})."]
    samples = run.samples
    if samples:
        conditions = sorted({s["condition"] for s in samples if s.get("condition")})
        dropped = run.stage_result("load").summary.get("n_genes_not_in_all_samples", 0)
        parts.append(f"{len(samples)} samples" + (f" ({_join_or(conditions).replace(' or ', ' and ')} conditions)" if conditions else "")
                     + " were combined, keeping genes measured in every sample"
                     + (f" ({dropped:,} genes measured in only some samples were left out)" if dropped else "") + ".")
    if status["filter"] == "done":
        f, s = p.filter, summary["filter"]
        if f.sample_limits:
            parts.append(f"Quality control limits were set for each sample from its own cells (median ± median "
                         f"absolute deviations: {p.suggest.n_mads_genes:g} for genes per cell on a log scale, "
                         f"{p.suggest.n_mads_mito:g} above the median for mitochondrial share; Heumos et al. 2023)"
                         if p.suggest is not None else "Quality control limits were set for each sample")
            parts[-1] += (f"; genes detected in fewer than {f.min_cells} cells were removed, leaving "
                          f"{s['n_cells_kept']:,} cells and {s['n_genes_kept']:,} genes.")
            f = None
        limits = [f"fewer than {f.min_genes:,} genes"] if f else []
        if f and f.max_n_genes_by_counts is not None:
            limits.append(f"{f.max_n_genes_by_counts:,} or more genes")
        if f and f.max_pct_counts_mt is not None:
            limits.append(f"{f.max_pct_counts_mt:g}% or more mitochondrial counts (genes starting with "
                          f"'{p.qc.mito_prefix}')")
        how = (f" (limits set at {p.suggest.n_mads_genes:g} median absolute deviations from the median for genes "
               f"per cell, on a log scale, and {p.suggest.n_mads_mito:g} above it for mitochondrial share)"
               if p.suggest is not None else "")
        if f:
            parts.append(f"Cells with {_join_or(limits)} were removed{how}, and genes detected in fewer than "
                         f"{f.min_cells} cells were removed, leaving {s['n_cells_kept']:,} cells and "
                         f"{s['n_genes_kept']:,} genes.")
    if status["preprocess"] == "done":
        s = summary["preprocess"]
        target = (f"the median total count ({_amount(s['target_sum'])})" if s.get("target_is_median")
                  else f"{_amount(s['target_sum'])} counts")
        hvg = (f"the {p.hvg.n_top_genes:,} most variable genes" if p.hvg.n_top_genes
               else f"highly variable genes (flavor '{p.hvg.flavor}', min mean {p.hvg.min_mean:g}, max mean "
                    f"{p.hvg.max_mean:g}, min normalized dispersion {p.hvg.min_disp:g}; {s['n_hvg']:,} genes)")
        regress = (f" Effects of {', '.join(p.regress_out.keys)} were regressed out." if s.get("regressed") else "")
        parts.append(f"Counts were normalized to {target} per cell and log-transformed (log1p); {hvg} were "
                     f"selected.{regress} Data were scaled to unit variance and clipped at {p.scale.max_value:g}.")
    if status["reduce_dimensions"] == "done":
        parts.append(f"Principal component analysis ({p.pca.n_comps} components, solver '{p.pca.svd_solver}') was "
                     f"followed by a neighbor graph ({p.neighbors.n_neighbors} neighbors, {p.neighbors.n_pcs} PCs) "
                     f"and UMAP (umap-learn {versions.get('umap-learn')}) for visualization.")
        integ = summary.get("reduce_dimensions", {}).get("integration")
        if integ and integ["method"] == "harmony":
            parts.append(f"Batch effects were corrected with Harmony (harmonypy {versions.get('harmonypy', 'version not recorded')}; "
                         f"Korsunsky et al. 2019) on the {p.neighbors.n_pcs} PCs used for the graph, with "
                         f"{integ['n_batches']} batches from \"{integ['batch_key']}\" (theta {p.integration.theta:g}, "
                         f"up to {p.integration.max_iter} rounds); the neighbor graph, clustering and UMAP used the "
                         f"corrected coordinates, and expression values were not changed. The mean sample mixing score "
                         f"(neighbor-label entropy, 0–1) rose from {integ['mixing_before']:.2f} to {integ['mixing_after']:.2f}.")
        elif integ:
            parts.append(f"No batch correction was applied; the mean sample mixing score (neighbor-label entropy, 0–1) "
                         f"across the {integ['n_batches']} batches was {integ['mixing_before']:.2f}.")
    if status["cluster"] == "done":
        parts.append(f"Cells were clustered with the Leiden algorithm (leidenalg {versions.get('leidenalg')}, "
                     f"resolution {p.leiden.resolution:g}), giving {summary['cluster']['n_clusters']} clusters.")
    if status["markers"] == "done":
        test = payloads.TEST_NAMES.get(p.markers.method, p.markers.method)
        parts.append(f"Marker genes were ranked for each cluster against all other cells ({test} test, "
                     f"Benjamini–Hochberg correction); because clusters were defined from the same cells, "
                     f"p-values were used only for ranking.")
    if status["annotate"] == "done":
        m, t = p.annotation.markers, p.annotation.test
        refs = {r["name"]: r for r in summary["annotate"].get("references") or [summary["annotate"]["database"]]}
        panglao = refs.get("PanglaoDB")
        panglao_text = (f"{panglao['name']} ({panglao['version']}; {m.species}; tissues: {', '.join(m.tissues)}; "
                        f"canonical markers with sensitivity above {m.min_sensitivity:g})") if panglao else None
        census = refs.get("CELLxGENE Census")
        if census:
            c = CensusMarkerParams()
            source = (f"marker tables built from CZ CELLxGENE Discover Census (release {census['version']}; "
                      f"healthy, primary {m.species} {census['tissues'][0]} cells; cell types with at least "
                      f"{c.min_cells} cells from {c.min_studies} or more studies, up to {c.cells_per_type} cells each; "
                      f"each type's top {c.top_n_genes} genes by Wilcoxon test against the cell types it is not "
                      f"related to in the Cell Ontology, {census.get('cell_ontology_version') or 'version not recorded'})")
        else:
            source = panglao_text
        parts.append(f"Cell types were predicted with an over-representation test (one-sided Fisher's exact test of "
                     f"each cluster's top {t.top_n_genes} marker genes against each cell type's markers from "
                     f"{source}; Benjamini–Hochberg adjusted). A cell type was assigned when it had at least "
                     f"{t.min_overlap} matching markers, an adjusted p-value below {t.max_adj_p:g}, and an adjusted "
                     f"p-value at least {t.ambiguity_fold:g} times smaller than any other significant candidate"
                     + (" not related to it in the Cell Ontology" if census else "")
                     + "; otherwise the cluster was left unresolved.")
        if census and panglao:
            n = sum(c["label_source"] == "PanglaoDB" for c in summary["annotate"]["clusters"])
            parts.append(f"Clusters left unresolved were tested in the same way against {panglao_text}; "
                         f"{n} {'cluster was' if n == 1 else 'clusters were'} labeled this way. A confident Census "
                         f"label was never replaced.")
        ct = summary["annotate"].get("celltypist")
        if ct and not ct.get("error"):
            parts.append(f"As an independent second opinion, CellTypist ({versions.get('celltypist', 'version not recorded')}; "
                         f"Domínguez Conde et al. 2022) with the ready-made model {ct['model'].removesuffix('.pkl')} "
                         f"({ct['version']}) labeled every cell from its log-normalized expression (10,000 counts per cell, "
                         f"as the model expects; {ct['n_genes_matched']:,} of the model's {ct['n_model_genes']:,} genes "
                         f"present), and each cluster took the majority vote. Its answers were reported beside the "
                         f"marker-based prediction and did not change it.")
        if run.review:
            n = len(run.review)
            parts.append(f"Predicted labels were reviewed by the analyst ({n} "
                         f"{'cluster' if n == 1 else 'clusters'} confirmed or renamed).")
    parts.append(f"Random seeds were fixed (PCA, neighbors, Leiden and UMAP: {p.pca.random_state}).")
    return " ".join(parts)


def samples_section(run: AnalysisRun) -> list[str]:
    """With several samples: one row per sample (condition, cells loaded, cells kept after QC) and,
    when Cell Type Annotation is up to date, its three largest predicted cell types. Descriptive."""
    samples = run.samples
    if not samples:
        return []
    status = run.status()
    kept = (run.result("filter").obs["sample"].astype(str).value_counts().to_dict()
            if status["filter"] == "done" else {})
    top: dict[str, str] = {}
    if status["annotate"] == "done":
        table = integration.composition(run.result("annotate"), "sample", "cell_type_predicted")
        for name, rows in table.groupby("sample"):
            best = rows.sort_values("share", ascending=False).head(3)
            top[name] = ", ".join(f"{r.group} {100 * r.share:.0f}%" for r in best.itertuples())
    out = ["<h2>Samples</h2><table><tr><th>Sample</th><th>Condition</th><th class='num'>Cells loaded</th>"
           "<th class='num'>Kept after QC</th>" + ("<th>Largest predicted cell types</th>" if top else "") + "</tr>"]
    for s in samples:
        out.append(f"<tr><td>{_e(s['name'])}</td><td>{_e(s.get('condition') or '–')}</td>"
                   f"<td class='num'>{s['n_cells']:,}</td><td class='num'>"
                   f"{_count(kept.get(s['name']))}</td>"
                   + (f"<td>{_e(top.get(s['name'], '–'))}</td>" if top else "") + "</tr>")
    out.append("</table>")
    if top:
        out.append("<p class='note'>A description of each sample, not a test: with one sample per donor and condition, "
                   "differences in composition aren't called significant.</p>")
    return out


def _count(n: int | None) -> str:
    return "–" if n is None else f"{n:,}"


def _join_or(items: list[str]) -> str:
    """"a", "a or b", "a, b, or c"."""
    if len(items) <= 2:
        return " or ".join(items)
    return ", ".join(items[:-1]) + ", or " + items[-1]


def _amount(value: float) -> str:
    """A count as written: 10,000 or 2,242.5 (no rounding away a half)."""
    return f"{value:,.0f}" if float(value).is_integer() else f"{value:,.1f}"


def run_report_html(run: AnalysisRun, dataset_label: str | None) -> str:
    """The whole report as one HTML page (no outside files, so it can be emailed or printed)."""
    record = run.record()
    status = run.status()
    now = datetime.now(timezone.utc).strftime("%d %b %Y, %H:%M UTC")
    load = next((e["summary"] for e in run.history if e["stage"] == "load"), {})
    out = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Run report: "
           f"{_e(run.file_name or 'analysis')}</title><style>{CSS}</style></head><body>",
           "<h1>Single-cell RNA-seq run report</h1>",
           f"<p class='meta'>{_e(dataset_label or run.file_name or '')} · created {now} · Orditus scRNA-seq "
           "analysis app</p>"]

    out.append("<h2>Dataset</h2><table>")
    rows = [("File", run.file_name), ("Format", load.get("file_format")),
            ("Fingerprint (SHA-256)", f"<code>{_e(load.get('input_sha256', ''))}</code>"),
            ("Cells × genes as loaded", f"{load['n_cells']:,} × {load['n_genes']:,}" if load else "–"),
            ("Species and tissue", payloads.metadata_text(run.metadata) if run.metadata else "not confirmed"),
            ("Confirmed by", {"dataset": "the dataset's metadata", "user": "the analyst"}.get(
                run.metadata.get("confirmed_by"), "–")),
            ("About the sample", run.metadata.get("about") or "–")]
    out += [f"<tr><th>{_e(k)}</th><td>{v if k.startswith('Fingerprint') else _e(v)}</td></tr>" for k, v in rows]
    out.append("</table>")
    out += samples_section(run)

    out.append("<h2>Methods</h2><p class='note'>Written from the settings used; edit as needed for a paper "
               "or thesis.</p>")
    out.append(f"<p class='methods'>{_e(methods_text(run))}</p>")

    out.append("<h2>Steps</h2><table><tr><th>Step</th><th>Status</th><th>Result</th></tr>")
    for key, title, stages in payloads.STEPS:
        state = "up to date" if all(status[s] == "done" for s in stages) else (
            "needs re-running" if any(status[s] == "needs re-running" for s in stages) else "not run")
        result = payloads._summary(run, key) if state == "up to date" else ""
        out.append(f"<tr><td>{_e(title)}</td><td>{_e(state)}</td><td>{_e(result)}</td></tr>")
    out.append("</table>")

    out.append("<h2>Settings</h2><table><tr><th>Setting</th><th>Value</th></tr>")
    out += [f"<tr><td>{_e(item['label'])}</td><td>{_e(item['value'])}</td></tr>"
            for item in payloads.settings_record(run)]
    out.append("</table>")

    if status["annotate"] == "done":
        summary = run.stage_result("annotate").details["summary"]
        stage = run.stage_result("annotate").summary
        refs = stage.get("references") or [stage["database"]]
        used = "; ".join(f"{r['name']} {r['version']}" for r in refs)
        out.append(f"<h2>Cell types</h2><p class='note'>Predicted from marker evidence ({_e(used)}); predictions to "
                   "review, not confirmed identities. Your labels are your own decisions, kept separately.</p>")
        cites = [r["citation"] for r in refs if r.get("citation")]
        ct_summary = stage.get("celltypist")
        if ct_summary and not ct_summary.get("error"):
            cites.append(ct_summary["citation"])
        if cites:
            out.append(f"<p class='note'>Cite: {_e('; '.join(cites))}.</p>")
        ct_rows = ({c["cluster"]: c for c in ct_summary["clusters"]}
                   if ct_summary and not ct_summary.get("error") else None)
        out.append("<table><tr><th>Cluster</th><th class='num'>Cells</th><th>Predicted cell type</th><th>From</th>"
                   "<th>Status</th><th>Matching markers</th><th class='num'>Adjusted p-value</th>"
                   + ("<th>CellTypist (share of votes)</th>" if ct_rows else "") + "<th>Your label</th><th>Your note</th></tr>")
        for row in summary.clusters.itertuples():
            r = run.review.get(str(row.cluster))
            mine = f"{_e(r.label)} ({'confirmed' if r.decision == 'confirmed' else 'renamed'})" if r else "Not reviewed"
            source = {"Census": "CELLxGENE Census"}.get(row.label_source, row.label_source) if row.label_source else "–"
            out.append(f"<tr><td>{_e(row.cluster)}</td><td class='num'>{row.n_cells:,}</td><td>{_e(row.predicted_label)}"
                       f"</td><td>{_e(source)}</td><td>{_e(row.status.replace('_', ' '))}</td><td>{_e(row.overlap)} of {_e(row.n_markers)}: "
                       f"{_e(row.matching_markers or '–')}</td><td class='num'>{row.adj_p:.2g}</td>"
                       + (f"<td>{_e(ct_rows[str(row.cluster)]['label'])} ({100 * ct_rows[str(row.cluster)]['share']:.0f}%)</td>"
                          if ct_rows and str(row.cluster) in ct_rows else ("<td>–</td>" if ct_rows else ""))
                       + f"<td>{mine}</td>"
                       f"<td>{_e(r.note if r else '')}</td></tr>")
        out.append("</table>")

    out.append("<h2>Reproducibility</h2><table><tr><th>Package</th><th>Version</th></tr>")
    out += [f"<tr><td>{_e(k)}</td><td>{_e(v)}</td></tr>" for k, v in record["versions"].items()]
    out.append("</table>")
    seeds = ", ".join(f"{k} {v}" for k, v in record["random_states"].items())
    out.append(f"<p>Random seeds: {_e(seeds)}.</p>")

    out.append("<h2>Run log</h2><table><tr><th>Time (UTC)</th><th>Step</th><th>What happened</th></tr>")
    out += [f"<tr><td>{_e(e['ran_at'].replace('T', ' ').replace('+00:00', ''))}</td><td>{_e(e['label'])}</td>"
            f"<td>{_e(e['text'])}</td></tr>" for e in payloads.run_log(run)]
    out.append("</table>")
    out.append("<p class='note'>The full machine-readable record (every setting, the history and the review) is "
               "stored in the exported .h5ad file under <code>uns[\"orditus\"]</code>.</p></body></html>")
    return "".join(out)
