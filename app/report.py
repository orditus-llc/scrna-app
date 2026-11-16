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


def _p(value) -> str:
    """An adjusted p-value for the report ("–" when there is none, e.g. no reference)."""
    return "–" if value is None or value != value else f"{value:.2g}"


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
        within = (" within each sample (Scanpy's batch_key = sample; genes ranked by the number of samples in "
                  "which they are variable)" if s.get("hvg_per_sample") else "")
        parts.append(f"Counts were normalized to {target} per cell and log-transformed (log1p); {hvg} were "
                     f"selected{within}.{regress} Data were scaled to unit variance and clipped at {p.scale.max_value:g}.")
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
        census = summary["annotate"].get("database")
        tissue = p.annotation.census.tissue
        if census:
            c = CensusMarkerParams()
            source = (f"marker tables built from CZ CELLxGENE Discover Census (release {census['version']}; "
                      f"healthy, primary {m.species} {census['tissues'][0]} cells; cell types with at least "
                      f"{c.min_cells} cells from {c.min_studies} or more studies, up to {c.cells_per_type} cells each; "
                      f"each type's top {c.top_n_genes} genes by Wilcoxon test against the cell types it is not "
                      f"related to in the Cell Ontology, {census.get('cell_ontology_version') or 'version not recorded'})")
            parts.append(f"Cell types were predicted with an over-representation test (one-sided Fisher's exact test of "
                         f"each cluster's top {t.top_n_genes} marker genes against each cell type's markers from "
                         f"{source}; Benjamini–Hochberg adjusted). A cell type was assigned when it had at least "
                         f"{t.min_overlap} matching markers, an adjusted p-value below {t.max_adj_p:g}, and an adjusted "
                         f"p-value at least {t.ambiguity_fold:g} times smaller than any other significant candidate "
                         "not related to it in the Cell Ontology; otherwise the cluster was left unresolved.")
        else:
            parts.append(f"No CELLxGENE Census marker reference was available for {m.species} {tissue}, so no "
                         "marker-based cell type prediction was made.")
        ct = summary["annotate"].get("celltypist")
        if ct and not ct.get("error"):
            parts.append(f"As an independent second opinion, CellTypist ({versions.get('celltypist', 'version not recorded')}; "
                         f"Domínguez Conde et al. 2022) with the ready-made model {ct['model'].removesuffix('.pkl')} "
                         f"({ct['version']}) labeled every cell from its log-normalized expression (10,000 counts per cell, "
                         f"as the model expects; {ct['n_genes_matched']:,} of the model's {ct['n_model_genes']:,} genes "
                         f"present), and each cluster took the majority vote.")
        final = summary["annotate"].get("final")
        if final:
            cs = run.params.annotation.consensus
            named = sum(1 for c in final if c["final_label"] != "Unknown")
            parts.append("Both answers were mapped to broad cell types (Census labels through the Cell Ontology, "
                         "CellTypist labels through the app's label table) and combined per cluster: a cluster was named "
                         "when both were confident and agreed, or when only one was confident (marked as supported by one "
                         "check only); it was labeled Unknown when they disagreed or neither was confident. CellTypist "
                         f"counted as confident when at least {cs.min_celltypist_share:.0%} of the cluster's cells voted "
                         f"for its label and at least {cs.min_celltypist_gene_match:.0%} of the model's genes were present. "
                         f"{named} of {len(final)} clusters were named.")
        if run.review:
            n = len(run.review)
            parts.append(f"Predicted labels were reviewed by the analyst ({n} "
                         f"{'cluster' if n == 1 else 'clusters'} confirmed or renamed).")
    if run.comparison and run.comparison_status() == "done":
        c = run.comparison
        cs, d = c["settings"], c["summary"].design
        grouping = dict(payloads.CONDITION_GROUPINGS)[cs["group_by"]].lower()
        paired = (f", paired by {d.paired_by} (subjects with a sample in both conditions)" if d.paired_by else "")
        n_ref = int((d.samples["side"] == "reference").sum())
        n_test = int((d.samples["side"] == "test").sum())
        parts.append(
            f"Conditions were compared with every sample counted once: {d.test} ({n_test} samples) against "
            f"{d.reference} ({n_ref} samples){paired}, with cells grouped by {grouping}. Changes in each group's share "
            "of cells were tested on logit-transformed per-sample shares (pseudo-count 0.5, as in propeller; Phipson et al. "
            f"2022) with {'a paired' if d.paired_by else 'an'} ordinary t-test, Benjamini–Hochberg-adjusted across groups. "
            "Changes in gene expression were tested per group on pseudobulk profiles (raw counts summed per sample; samples "
            f"with fewer than {cs['min_cells']} of the group's cells left out) with PyDESeq2 "
            f"{versions.get('pydeseq2', '(version not recorded)')} (design {'~ pair + condition' if d.paired_by else '~ condition'}; "
            f"genes with fewer than {cs['min_gene_counts']} counts left out); significant = adjusted p-value below "
            f"{cs['max_adj_p']:g}.")
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
        final = run.stage_result("annotate").details["final"].set_index("cluster")
        stage = run.stage_result("annotate").summary
        db = stage.get("database")
        used = (f"{db['name']} {db['version']}" if db
                else f"no CELLxGENE Census reference for {run.params.annotation.census.tissue} yet")
        out.append(f"<h2>Cell types</h2><p class='note'>Predicted from two checks, Census marker genes ({_e(used)}) and "
                   "CellTypist, compared as broad cell types; predictions to review, not confirmed identities. Your "
                   "labels are your own decisions, kept separately.</p>")
        cites = [db["citation"]] if db and db.get("citation") else []
        ct_summary = stage.get("celltypist")
        if ct_summary and not ct_summary.get("error"):
            cites.append(ct_summary["citation"])
        if cites:
            out.append(f"<p class='note'>Cite: {_e('; '.join(cites))}.</p>")
        ct_rows = ({c["cluster"]: c for c in ct_summary["clusters"]}
                   if ct_summary and not ct_summary.get("error") else None)
        out.append("<table><tr><th>Cluster</th><th class='num'>Cells</th><th>Suggested cell type</th><th>Confidence</th><th>Agreement</th>"
                   "<th>Marker check</th><th>Marker check status</th><th>Matching markers</th><th class='num'>Adjusted p-value</th>"
                   + ("<th>CellTypist (share of votes)</th>" if ct_rows else "") + "<th>Your label</th><th>Your note</th></tr>")
        for row in summary.clusters.itertuples():
            r = run.review.get(str(row.cluster))
            mine = f"{_e(r.label)} ({'confirmed' if r.decision == 'confirmed' else 'renamed'})" if r else "Not reviewed"
            f = final.loc[str(row.cluster)]
            out.append(f"<tr><td>{_e(row.cluster)}</td><td class='num'>{row.n_cells:,}</td><td>{_e(f['final_label'])}</td>"
                       f"<td>{_e(payloads.CONFIDENCE_TEXT[payloads.consensus_tools.CONFIDENCE[f['agreement']]])}</td>"
                       f"<td>{_e(payloads.AGREEMENT_TEXT[f['agreement']])}</td><td>{_e(row.predicted_label)}"
                       f"</td><td>{_e(row.status.replace('_', ' '))}</td><td>{_e(row.overlap)} of {_e(row.n_markers)}: "
                       f"{_e(row.matching_markers or '–')}</td><td class='num'>{_p(row.adj_p)}</td>"
                       + (f"<td>{_e(ct_rows[str(row.cluster)]['label'])} ({100 * ct_rows[str(row.cluster)]['share']:.0f}%)</td>"
                          if ct_rows and str(row.cluster) in ct_rows else ("<td>–</td>" if ct_rows else ""))
                       + f"<td>{mine}</td>"
                       f"<td>{_e(r.note if r else '')}</td></tr>")
        out.append("</table>")

    if run.comparison and run.comparison_status() == "done":
        c = run.comparison
        summary, max_adj_p = c["summary"], c["settings"]["max_adj_p"]
        d = summary.design
        out.append(f"<h2>Compare Conditions</h2><p class='note'>{_e(d.test)} against {_e(d.reference)}"
                   + (f", paired by {_e(d.paired_by)}" if d.paired_by else "") + ". Each sample counted once; "
                   f"significant = adjusted p-value below {max_adj_p:g}.</p>")
        out.append("<table><tr><th>Group</th><th class='num'>Share in " + _e(d.reference) + "</th><th class='num'>Share in "
                   + _e(d.test) + "</th><th class='num'>Adjusted p-value (share)</th><th>Genes up / down</th>"
                   "<th>Top genes (by p-value)</th></tr>")
        for row in summary.composition.itertuples():
            t = summary.expression.get(row.group)
            if t is None:
                genes, top = "not tested", ""
            elif t.status != "tested":
                genes, top = _e(t.reason), ""
            else:
                up, down = payloads.condition_tools.significant_counts(t.results, max_adj_p)
                genes = f"{up:,} / {down:,}"
                top = ", ".join(f"{g} ({'+' if lfc > 0 else ''}{lfc:.1f})"
                                for g, lfc in t.results.head(5)[["gene", "log2_fold_change"]].itertuples(index=False))
            out.append(f"<tr><td>{_e(row.group)}</td><td class='num'>{row.mean_share_reference:.1%}</td>"
                       f"<td class='num'>{row.mean_share_test:.1%}</td><td class='num'>{_p(row.adj_p)}</td>"
                       f"<td>{genes}</td><td>{_e(top)}</td></tr>")
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
