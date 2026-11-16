// Step 6: marker genes. Settings on the left; on the right a cluster picker, that cluster's
// top markers, a dot plot of every cluster's top markers across all clusters, and a heatmap
// of those markers' mean expression per cluster.
import type { Data, Datum, Layout, Shape } from "plotly.js-dist-min";
import { useMemo, useState } from "react";

import { CsvDownload } from "../components/DownloadMenu";
import { InfoBox, Note, RuledSelect } from "../components/fields";
import { MoreSection } from "../components/MoreSection";
import { PlotlyChart } from "../components/PlotlyChart";
import { Tabs } from "../components/Tabs";
import { Prerequisite, SettingsAndResults, StepHeader } from "../components/StepPage";
import { clusterColors } from "../lib/clusterColors";
import { cn } from "../lib/utils";
import { goTo } from "../navigation";
import type { MarkerRow, MarkersPayload, Step } from "../types";
import { Empty, fmt, Loading, prerequisiteKey, SettingsReset, StepBar, SectionTitle, stepOf, useRunButton, type PageProps } from "./stepShared";

// Mean expression: light to dark teal (the accent's family).
const TEAL_SCALE: [number, string][] = [
  [0, "#f0fdfa"],
  [0.35, "#5eead4"],
  [0.7, "#0d9488"],
  [1, "#134e4a"],
];
const SHOWN_OPTIONS = [5, 10, 25, 50].map((n) => ({ value: String(n), label: `${n} genes` }));

export function MarkersPage({ workflow, data, onNavigate }: PageProps & { data: MarkersPayload | undefined }) {
  const step = stepOf(workflow, "marker_genes");
  if (!data || !step) return <Loading />;
  return <MarkersContent workflow={workflow} data={data} step={step} onNavigate={onNavigate} />;
}

function MarkersContent({ workflow, data, step, onNavigate }: PageProps & { data: MarkersPayload; step: Step }) {
  const reference = data.settings.reference;
  const [method, setMethod] = useState(data.settings.current.method);
  const [shown, setShown] = useState("10"); // display only: genes listed per cluster
  const run = useRunButton<{ method: string }>("run_markers", "marker_genes", data, workflow.notice);
  const r = data.results;
  const [selected, setSelected] = useState<string | null>(null);
  const cluster = r ? (selected ?? r.groups[0]?.cluster ?? null) : null;
  const changedSinceRun = r !== null && method !== data.settings.current.method;

  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <StepHeader step={step} />
      {step.prerequisite && <Prerequisite title={step.prerequisite} text={step.prerequisite_text} onOpen={() => goTo(onNavigate, prerequisiteKey(workflow, step))} />}
      {step.reason && <Note dot="bg-accent">{step.reason}</Note>}

      <SettingsAndResults
        settings={
          <>
            <SectionTitle>Settings</SectionTitle>
            <RuledSelect
              id="mk_method"
              label="Statistical test"
              help="How each cluster is compared with all other cells. Default: Wilcoxon, the test Scanpy's documentation recommends for marker genes."
              value={method}
              options={data.methods}
              onChange={setMethod}
              reference={reference.method}
            />
            <MoreSection title="Fixed settings (validated method)">
            <InfoBox label="Fixed settings" className="my-2">
            <dl className="m-0 grid grid-cols-[max-content_minmax(0,1fr)] gap-x-6 gap-y-[6px] text-[14px]">
              <dt className="text-muted">Genes</dt>
              <dd className="m-0">{data.fixed.use_raw ? "all genes, log-normalized (adata.raw)" : "highly variable genes"}</dd>
              <dt className="text-muted">Compared with</dt>
              <dd className="m-0">{data.fixed.reference === "rest" ? "all other cells" : data.fixed.reference}</dd>
              <dt className="text-muted">Correction</dt>
              <dd className="m-0">{data.fixed.corr_method === "benjamini-hochberg" ? "Benjamini–Hochberg" : data.fixed.corr_method}</dd>
              <dt className="text-muted">Tie correction</dt>
              <dd className="m-0">{data.fixed.tie_correct ? "on" : "off"}</dd>
            </dl>
            </InfoBox>
            </MoreSection>
            <SettingsReset onReset={() => setMethod(reference.method)} changedSinceRun={changedSinceRun} />
          </>
        }
        results={
          r && cluster ? (
            <MarkerResults results={r} cluster={cluster} onSelect={setSelected} shown={shown} setShown={setShown} />
          ) : (
            <>
              <SectionTitle>Marker genes</SectionTitle>
              <Empty>Run Marker Genes to rank the genes of every cluster.</Empty>
            </>
          )
        }
      />
      <StepBar
        label="Marker Genes"
        nextTitle={stepOf(workflow, "cell_type_annotation")?.title}
        status={step.status}
        blocked={step.prerequisite !== null}
        busy={run.busy}
        busyText="Running Marker Genes…"
        error={run.error}
        changedSinceRun={changedSinceRun}
        willClear={data.will_clear}
        onRun={() => run.click({ method })}
        onContinue={() => goTo(onNavigate, "cell_type_annotation")}
        prerequisite={step.prerequisite}
        prerequisiteText={step.prerequisite_text}
        confirmation={step.confirmation}
      />
    </div>
  );
}

type Results = NonNullable<MarkersPayload["results"]>;

function MarkerResults({
  results,
  cluster,
  onSelect,
  shown,
  setShown,
}: {
  results: Results;
  cluster: string;
  onSelect: (cluster: string) => void;
  shown: string;
  setShown: (n: string) => void;
}) {
  const [view, setView] = useState<"table" | "dots" | "heatmap">("table");
  const [heatGenes, setHeatGenes] = useState(String(results.heatmap.default)); // display only
  const order = results.groups.map((g) => g.cluster);
  const colors = useMemo(() => clusterColors(order), [results]); // eslint-disable-line react-hooks/exhaustive-deps
  const size = results.groups.find((g) => g.cluster === cluster)?.n_cells ?? 0;
  const rows = (results.top[cluster] ?? []).slice(0, Number(shown));

  return (
    <>
      <p className="m-0 text-[14px] text-body">
        <b className="text-ink">Read p-values as a ranking:</b> the same cells made the clusters, so they are far too small.
      </p>
      <Tabs
        label="Marker gene views"
        tabs={[
          { key: "table", label: "Top markers" },
          { key: "dots", label: "Dot plot" },
          { key: "heatmap", label: "Heatmap" },
        ]}
        active={view}
        onChange={setView}
        right={
          view === "table" ? (
            <label className="flex items-center gap-2 text-[13px]">
              Show
              <select
                id="mk_shown"
                value={shown}
                onChange={(e) => setShown(e.target.value)}
                className="rounded-none border border-field bg-white px-2 py-1 text-[13px]"
              >
                {SHOWN_OPTIONS.map((o) => (
                  <option key={o.value} value={o.value}>
                    {o.label}
                  </option>
                ))}
              </select>
            </label>
          ) : view === "heatmap" ? (
            <label className="flex items-center gap-2 text-[13px]">
              Genes per cluster
              <select
                id="mk_heat_genes"
                value={heatGenes}
                onChange={(e) => setHeatGenes(e.target.value)}
                className="rounded-none border border-field bg-white px-2 py-1 text-[13px]"
              >
                {results.heatmap.options.map((n) => (
                  <option key={n} value={String(n)}>
                    {n}
                  </option>
                ))}
              </select>
            </label>
          ) : undefined
        }
      >
      {view === "table" ? (
      <>
      <p className="m-0 text-[13px] text-body">
        Cluster <b className="text-ink">{cluster}</b> · <span className="font-mono">{fmt(size)}</span> cells. Pick another cluster:
      </p>
      <div role="group" aria-label="Cluster" className="flex flex-wrap gap-2">
        {order.map((c) => (
          <button
            key={c}
            type="button"
            aria-pressed={c === cluster}
            onClick={() => onSelect(c)}
            className={cn(
              "inline-flex items-center gap-2 rounded-none border px-3 py-[6px] text-[14px] font-bold",
              c === cluster ? "border-ink bg-ink text-white" : "border-field bg-white text-ink hover:bg-chip",
            )}
          >
            <span aria-hidden="true" className="size-[10px] rounded-full ring-1 ring-black/25" style={{ background: colors[c] }} />
            {c}
          </button>
        ))}
      </div>
      <div className="flex justify-end">
        <CsvDownload
          name="Marker Genes top markers all clusters"
          columns={[
            { key: "cluster", header: "Cluster" },
            { key: "rank", header: "Rank" },
            { key: "gene", header: "Gene" },
            { key: "logfoldchange", header: "Log fold change" },
            { key: "frac_in_group", header: "Fraction expressing in cluster" },
            { key: "frac_in_rest", header: "Fraction expressing in other cells" },
            { key: "score", header: "Score" },
            { key: "pval_adj", header: "Adjusted p-value (ranking only)" },
          ]}
          rows={() => order.flatMap((c) => (results.top[c] ?? []).map((row) => ({ cluster: c, ...row })))}
        />
      </div>
      <MarkerTable rows={rows} />
      </>
      ) : view === "dots" ? (
        <DotPlot
          results={results}
          order={order}
          cluster={cluster}
          onSelect={(c) => {
            onSelect(c);
            setView("table");
          }}
        />
      ) : (
        <Heatmap
          results={results}
          order={order}
          colors={colors}
          genesPerCluster={heatGenes}
          onSelect={(c) => {
            onSelect(c);
            setView("table");
          }}
        />
      )}
      </Tabs>
    </>
  );
}

function MarkerTable({ rows }: { rows: MarkerRow[] }) {
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  return (
    <div className="overflow-x-auto">
      <table className="w-full min-w-[620px] border-collapse text-[14px]">
        <thead>
          <tr className="border-b border-ink text-left">
            <th scope="col" className={cn(head, "w-12 pr-2")}>Rank</th>
            <th scope="col" className={cn(head, "pr-3")}>Gene</th>
            <th scope="col" className={cn(head, "px-3 text-right")}>Log fold change</th>
            <th scope="col" className={cn(head, "px-3")}>% expressing: cluster vs rest</th>
            <th scope="col" className={cn(head, "px-3 text-right")}>Score</th>
            <th scope="col" className={cn(head, "pl-3 text-right")}>Adj. p</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((row) => (
            <tr key={row.gene} className="border-b border-line">
              <td className="py-[10px] pr-2 font-mono text-muted">{row.rank}</td>
              <td className="py-[10px] pr-3 font-mono font-bold">{row.gene}</td>
              <td className="px-3 py-[10px] text-right font-mono">{num(row.logfoldchange, 2)}</td>
              <td className="px-3 py-[10px]">
                <div className="flex items-center gap-3">
                  <span className="w-[86px] font-mono text-[13px] whitespace-nowrap">
                    {percent(row.frac_in_group)} <span className="text-muted">/ {percent(row.frac_in_rest)}</span>
                  </span>
                  <span aria-hidden="true" className="flex w-[110px] flex-col gap-[3px]">
                    <span className="block h-[5px] bg-line">
                      <span className="block h-[5px] bg-accent" style={{ width: `${100 * (row.frac_in_group ?? 0)}%` }} />
                    </span>
                    <span className="block h-[5px] bg-line">
                      <span className="block h-[5px] bg-[#a1a1aa]" style={{ width: `${100 * (row.frac_in_rest ?? 0)}%` }} />
                    </span>
                  </span>
                </div>
              </td>
              <td className="px-3 py-[10px] text-right font-mono">{num(row.score, 1)}</td>
              <td className="py-[10px] pl-3 text-right font-mono text-muted">{pvalue(row.pval_adj)}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <p className="mt-2 mb-0 text-[13px] text-muted">
        Teal = this cluster, grey = all other cells.
      </p>
    </div>
  );
}

/** Every cluster's top markers (rows) across all clusters (columns). Size = % of cells
 * expressing; color = mean expression. Hover a dot for its numbers; click to pick a cluster. */
function DotPlot({
  results,
  order,
  cluster,
  onSelect,
}: {
  results: Results;
  order: string[];
  cluster: string;
  onSelect: (cluster: string) => void;
}) {
  const [hovered, setHovered] = useState<string | null>(null);
  const genes = results.dotplot.genes;
  const chart = useMemo(() => {
    const v = results.dotplot.values;
    const maxMean = Math.max(...v.map((d) => d.mean_expr), 1e-9);
    const data: Data[] = [
      {
        type: "scatter",
        mode: "markers",
        x: v.map((d) => `Cluster ${d.group}`),
        y: v.map((d) => d.gene),
        customdata: v.map((d) => [d.group, d.mean_expr, 100 * d.frac_expressing, d.rank, d.logfoldchange ?? NaN]),
        marker: {
          size: v.map((d) => 3 + 19 * Math.sqrt(d.frac_expressing)),
          color: v.map((d) => d.mean_expr),
          cmin: 0,
          cmax: maxMean,
          colorscale: TEAL_SCALE,
          showscale: true,
          colorbar: { title: { text: "Mean expression", side: "right" }, thickness: 12, outlinewidth: 0 },
          line: { color: "rgba(17,17,17,0.35)", width: 0.6 },
        },
        hovertemplate:
          "<b>%{y}</b> in cluster %{customdata[0]}<br>Mean expression %{customdata[1]:.2f}<br>" +
          "%{customdata[2]:.0f}% of cells express it<br>DE rank in this cluster: %{customdata[3]:,}<br>" +
          "Log fold change %{customdata[4]:.2f}<extra></extra>",
      },
    ];
    const band = (c: string, fill: string): Partial<Shape> => ({
      type: "rect", xref: "x", yref: "paper", x0: order.indexOf(c) - 0.5, x1: order.indexOf(c) + 0.5, y0: 0, y1: 1,
      fillcolor: fill, line: { width: 0 }, layer: "below",
    });
    const shapes = [band(cluster, "rgba(17,17,17,0.06)")];
    if (hovered && hovered !== cluster) shapes.push(band(hovered, "rgba(13,148,136,0.08)"));
    const layout: Partial<Layout> = {
      margin: { l: 90, r: 16, t: 16, b: 70 },
      xaxis: {
        type: "category", categoryorder: "array", categoryarray: order.map((c) => `Cluster ${c}`),
        tickangle: -40, showgrid: false, zeroline: false,
      },
      yaxis: { type: "category", categoryorder: "array", categoryarray: [...genes].reverse(), gridcolor: "#f4f4f5", tickfont: { family: "IBM Plex Mono, monospace", size: 12 } },
      shapes,
      hoverlabel: { bgcolor: "#ffffff", bordercolor: "#d4d4d8", font: { color: "#111111" } },
    };
    return { data, layout };
  }, [results, order, genes, cluster, hovered]);

  const groupOf = (pointIndex: number) => results.dotplot.values[pointIndex]?.group ?? null;
  return (
    <>
      <p className="m-0 text-[13px] text-body">
        Each cluster's top 3 markers. Dot size = share of cells expressing; color = mean expression. Click a dot to see
        that cluster's markers.
      </p>
      <PlotlyChart
        data={chart.data}
        layout={chart.layout}
        height={Math.max(320, 22 * genes.length + 110)}
        label="Dot plot of each cluster's top marker genes across all clusters"
        download={{
          name: "Marker Genes dot plot",
          title: "Top 3 marker genes of each cluster (size = % of cells expressing; color = mean expression)",
        }}
        onHover={(p) => setHovered(groupOf(p.pointIndex))}
        onUnhover={() => setHovered(null)}
        onClick={(p) => {
          const g = groupOf(p.pointIndex);
          if (g) onSelect(g);
        }}
      />
    </>
  );
}

/** Each cluster's top markers (rows, grouped by the cluster they mark) by cluster (columns):
 * mean expression, scaled from 0 to 1 within each gene. Click a square to see that cluster's
 * markers. */
function Heatmap({
  results,
  order,
  colors,
  genesPerCluster,
  onSelect,
}: {
  results: Results;
  order: string[];
  colors: Record<string, string>;
  genesPerCluster: string;
  onSelect: (cluster: string) => void;
}) {
  const rows = results.heatmap.rows[genesPerCluster] ?? [];
  const chart = useMemo(() => {
    const value = new Map(results.heatmap.values.map((v) => [`${v.group}|${v.gene}`, v]));
    const cell = (gene: string, c: string) => value.get(`${c}|${gene}`);
    const data: Data[] = [
      {
        type: "heatmap",
        x: order.map((c) => `Cluster ${c}`),
        y: rows.map((r) => r.gene),
        z: rows.map((r) => order.map((c) => cell(r.gene, c)?.scaled ?? null)),
        customdata: rows.map((r) => order.map((c) => [r.cluster, c, cell(r.gene, c)?.mean_expr ?? NaN])) as unknown as Datum[][],
        zmin: 0,
        zmax: 1,
        colorscale: TEAL_SCALE,
        xgap: 1,
        ygap: 1,
        colorbar: {
          title: { text: "Scaled mean expression", side: "right" },
          thickness: 12,
          outlinewidth: 0,
          tickvals: [0, 1],
          ticktext: ["lowest cluster", "highest"],
        },
        hovertemplate:
          "<b>%{y}</b> (top marker of cluster %{customdata[0]}) in cluster %{customdata[1]}<br>" +
          "Mean expression %{customdata[2]:.2f}<br>Scaled %{z:.2f} (1 = highest cluster)<extra></extra>",
      },
    ];
    // A colored square and a rule mark where each cluster's group of genes starts.
    const shapes: Partial<Shape>[] = [];
    rows.forEach((r, i) => {
      if (i > 0 && rows[i - 1].cluster === r.cluster) return;
      if (i > 0) shapes.push({ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: i - 0.5, y1: i - 0.5, line: { color: "#111111", width: 1 } });
      const end = rows.findIndex((s, j) => j > i && s.cluster !== r.cluster);
      shapes.push({
        type: "rect", xref: "paper", yref: "y", x0: -0.012, x1: -0.004, y0: i - 0.5, y1: (end === -1 ? rows.length : end) - 0.5,
        fillcolor: colors[r.cluster], line: { width: 0 },
      });
    });
    const layout: Partial<Layout> = {
      margin: { l: 96, r: 16, t: 16, b: 70 },
      xaxis: { type: "category", tickangle: -40, showgrid: false, zeroline: false },
      yaxis: { type: "category", autorange: "reversed", showgrid: false, tickfont: { family: "IBM Plex Mono, monospace", size: 12 } },
      shapes,
      hoverlabel: { bgcolor: "#ffffff", bordercolor: "#d4d4d8", font: { color: "#111111" } },
    };
    return { data, layout };
  }, [results, order, colors, rows]);

  return (
    <>
      <p className="m-0 text-[13px] text-body">
        Each cluster's top {genesPerCluster} markers, grouped by the cluster they mark. Color = mean expression in each
        cluster, scaled within each gene (darkest = the cluster where it is highest). Click a square to see that
        cluster's markers.
      </p>
      <PlotlyChart
        data={chart.data}
        layout={chart.layout}
        height={Math.max(320, 18 * rows.length + 110)}
        label="Heatmap of each cluster's top marker genes, mean expression per cluster"
        download={{
          name: `Marker Genes heatmap top ${genesPerCluster}`,
          title: `Top ${genesPerCluster} marker genes per cluster: mean expression per cluster, scaled within each gene`,
        }}
        onClick={(p) => {
          // On a heatmap Plotly reports the square as [row, column].
          const at = p.pointIndex as unknown as number | [number, number];
          const c = Array.isArray(at) ? order[at[1]] : undefined;
          if (c) onSelect(c);
        }}
      />
    </>
  );
}

function num(v: number | null, digits: number): string {
  return v === null ? "–" : v.toFixed(digits);
}

function percent(v: number | null): string {
  return v === null ? "–" : `${Math.round(100 * v)}%`;
}

function pvalue(v: number | null): string {
  if (v === null) return "–";
  if (v === 0) return "< 1e-300";
  return v < 0.001 ? v.toExponential(1) : v.toFixed(3);
}
