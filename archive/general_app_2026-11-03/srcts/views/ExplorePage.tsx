// Explore Results: for looking, not running. One large UMAP colored by cluster, predicted
// cell type, a QC value or any gene; groups highlight on hover, with their details beside.
import { useShinyInput } from "@posit-dev/shinyreact";
import { ArrowRight, Search } from "lucide-react";
import type { Data, Layout } from "plotly.js-dist-min";
import { useMemo, useState, type ReactNode } from "react";

import { Note } from "../components/fields";
import { ExportButtons } from "../components/ExportButtons";
import { GroupedUmap, percent, type UmapGroup } from "../components/GroupedUmap";
import { PlotlyChart } from "../components/PlotlyChart";
import { Button } from "../components/ui/button";
import { clusterColors } from "../lib/clusterColors";
import { reviewedGroups } from "../lib/reviewedGroups";
import { cn } from "../lib/utils";
import { goTo, type Page } from "../navigation";
import type { ExportInfo, ExplorePayload, GenePayload } from "../types";
import { STATUS_TEXT } from "./AnnotationPage";

type Mode = "cluster" | "cell_type" | "reviewed" | "n_genes" | "total_counts" | "pct_mt" | "gene";
type Results = NonNullable<ExplorePayload["results"]>;

const UNRESOLVED_COLOR = "#a1a1aa";
const NOT_EXPRESSED = "#e4e4e7";
// Expression: light to dark teal (the app's accent family).
const TEAL_SCALE: [number, string][] = [
  [0, "#ccfbf1"],
  [0.35, "#2dd4bf"],
  [0.7, "#0d9488"],
  [1, "#134e4a"],
];
const QC_LABELS = { n_genes: "Genes per cell", total_counts: "Counts per cell", pct_mt: "% mitochondrial" } as const;

export function ExplorePage({
  data,
  gene,
  exportInfo,
  onNavigate,
}: {
  data: ExplorePayload | undefined;
  gene: GenePayload | null | undefined;
  exportInfo: ExportInfo | undefined;
  onNavigate: (page: Page) => void;
}) {
  if (!data) return <p className="m-0 px-8 py-10 text-muted">Loading…</p>;
  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <div className="flex flex-col gap-[6px]">
        <h1 className="m-0 text-[34px] font-extrabold tracking-[-0.02em]">Explore Results</h1>
        <p className="m-0 max-w-[780px] text-[15px] leading-normal text-muted">
          Color the UMAP by cluster, cell type, a QC value or any gene. Nothing is run here.
        </p>
      </div>
      {exportInfo && <ExportButtons info={exportInfo} />}
      {!data.available || !data.results ? (
        <div className="flex flex-wrap items-center gap-3 rounded-[4px] border border-line bg-panel px-4 py-3 text-[14px]">
          <span>{data.reason}</span>
          <Button variant="secondary" size="small" onClick={() => goTo(onNavigate, "clustering")}>
            Open Clustering <ArrowRight aria-hidden="true" size={14} />
          </Button>
        </div>
      ) : (
        <Explorer results={data.results} gene={gene ?? null} onNavigate={onNavigate} />
      )}
    </div>
  );
}

function Explorer({ results, gene, onNavigate }: { results: Results; gene: GenePayload | null; onNavigate: (page: Page) => void }) {
  const [mode, setMode] = useState<Mode>("cluster");
  const [, sendGene] = useShinyInput<string | null>("explore_gene", null, { debounceMs: 0 });
  const [query, setQuery] = useState("");

  const clusterOrder = results.clusters.map((c) => c.cluster);
  const cColors = useMemo(() => clusterColors(clusterOrder), [results]); // eslint-disable-line react-hooks/exhaustive-deps
  const predictions = results.predictions;
  const labelOf = useMemo(
    () => Object.fromEntries((predictions?.clusters ?? []).map((c) => [c.cluster, c])),
    [predictions],
  );

  const pickGene = (name: string) => {
    setQuery(name);
    sendGene(name);
    setMode("gene");
  };

  const suggestions = useMemo(() => {
    const q = query.trim().toUpperCase();
    if (!q) return [];
    const starts = results.genes.filter((g) => g.toUpperCase().startsWith(q));
    return starts.slice(0, 8);
  }, [query, results.genes]);

  const modeButton = (key: Mode, label: string, enabled = true, note?: string) => (
    <button
      key={key}
      type="button"
      aria-pressed={mode === key}
      disabled={!enabled}
      onClick={() => setMode(key)}
      className={cn(
        "flex w-full flex-col rounded-[4px] px-3 py-2 text-left text-[14px] disabled:cursor-not-allowed",
        mode === key ? "bg-ink font-bold text-white" : "text-ink enabled:hover:bg-chip disabled:text-off",
      )}
    >
      {label}
      {note && <span className={cn("text-[12px] font-normal", mode === key ? "text-white/80" : "text-muted")}>{note}</span>}
    </button>
  );

  return (
    <div className="flex flex-wrap gap-8">
      <aside aria-label="Color cells by" className="flex min-w-0 flex-[1_1_240px] flex-col gap-4 sm:max-w-[280px]">
        <section className="flex flex-col gap-1">
          <h2 className="m-0 border-b border-ink pb-[10px] text-lg font-extrabold">Color cells by</h2>
          <div className="flex flex-col gap-1 pt-2">
            {modeButton("cluster", "Cluster", true, `${clusterOrder.length} clusters`)}
            {modeButton("cell_type", "Predicted cell type", predictions !== null, predictions ? "from marker evidence" : "run Cell Type Annotation first")}
            {modeButton(
              "reviewed",
              "Your labels",
              results.reviewed !== null,
              results.reviewed ? `${Object.keys(results.reviewed).length} of ${clusterOrder.length} clusters reviewed` : "review clusters in Cell Type Annotation",
            )}
            {modeButton("n_genes", QC_LABELS.n_genes)}
            {modeButton("total_counts", QC_LABELS.total_counts)}
            {modeButton("pct_mt", QC_LABELS.pct_mt, results.cells.pct_mt !== null)}
            {modeButton("gene", gene && gene.ok ? `Gene: ${gene.gene}` : "Gene", gene !== null && gene.ok)}
          </div>
        </section>

        <section className="flex flex-col gap-2">
          <h2 className="m-0 border-b border-ink pb-[10px] text-lg font-extrabold">Find a gene</h2>
          <form
            role="search"
            onSubmit={(e) => {
              e.preventDefault();
              const exact = results.genes.find((g) => g.toUpperCase() === query.trim().toUpperCase());
              if (exact || suggestions[0]) pickGene(exact ?? suggestions[0]);
              else if (query.trim()) pickGene(query.trim()); // the server explains and suggests close matches
            }}
            className="flex flex-col gap-2 pt-2"
          >
            <label htmlFor="explore_gene" className="text-[13px] text-muted">
              Any of the {results.genes.length.toLocaleString("en-US")} genes (log-normalized expression)
            </label>
            <div className="flex gap-2">
              <input
                id="explore_gene"
                type="search"
                autoComplete="off"
                spellCheck={false}
                value={query}
                placeholder="e.g. MS4A1"
                onChange={(e) => setQuery(e.target.value)}
                className="min-w-0 flex-auto rounded-[4px] border border-field bg-white px-[10px] py-2 font-mono text-[14px]"
              />
              <Button type="submit" size="small" aria-label="Show gene">
                <Search aria-hidden="true" size={15} />
              </Button>
            </div>
          </form>
          {suggestions.length > 0 && !(gene?.ok && gene.gene === query.trim()) && (
            <ul aria-label="Matching genes" className="m-0 flex list-none flex-wrap gap-1 p-0">
              {suggestions.map((g) => (
                <li key={g}>
                  <button type="button" onClick={() => pickGene(g)} className="rounded-[3px] bg-chip px-2 py-[2px] font-mono text-[12.5px] hover:bg-line">
                    {g}
                  </button>
                </li>
              ))}
            </ul>
          )}
          {gene && !gene.ok && (
            <Note dot="bg-accent">
              {gene.message}
              {gene.suggestions.length > 0 && (
                <span className="mt-1 flex flex-wrap gap-1">
                  Did you mean:
                  {gene.suggestions.map((s) => (
                    <button key={s} type="button" onClick={() => pickGene(s)} className="font-mono font-bold underline underline-offset-2">
                      {s}
                    </button>
                  ))}
                </span>
              )}
            </Note>
          )}
        </section>

        {results.unavailable.map((u) => (
          <Note key={u} dot="bg-off">
            {u}, so its results aren't shown here.
          </Note>
        ))}
      </aside>

      <section aria-label="UMAP" className="flex min-w-0 flex-[999_1_640px] flex-col gap-3">
        {mode === "cluster" && (
          <GroupedUmap
            title="Clusters on the UMAP"
            description="Hover or click a cluster to highlight it and see its details."
            listHeading="Cluster"
            groups={results.clusters.map((c) => ({ key: c.cluster, label: `Cluster ${c.cluster}`, n: c.n_cells }))}
            colors={cColors}
            cellGroup={results.cell_cluster}
            umap={results.umap}
            cells={results.cells}
            extraHover={predictions ? (i) => labelOf[results.cell_cluster[i]]?.predicted_label ?? "" : undefined}
            details={(key) => (
              <ClusterDetails
                cluster={key}
                results={results}
                prediction={key ? labelOf[key] : undefined}
                onGene={pickGene}
                onNavigate={onNavigate}
              />
            )}
          />
        )}
        {mode === "reviewed" && results.reviewed && <ReviewedView results={results} cColors={cColors} clusterOrder={clusterOrder} />}
        {mode === "cell_type" && predictions && (
          <CellTypeView results={results} cColors={cColors} onGene={pickGene} onNavigate={onNavigate} />
        )}
        {(mode === "n_genes" || mode === "total_counts" || mode === "pct_mt") && (
          <ContinuousView
            title={QC_LABELS[mode]}
            description="Colored by this QC value (before filtering)."
            values={mode === "pct_mt" ? results.cells.pct_mt! : mode === "n_genes" ? results.cells.n_genes : results.cells.total_counts}
            results={results}
            greyZeros={false}
            colorTitle={QC_LABELS[mode]}
          />
        )}
        {mode === "gene" && gene && gene.ok && (
          <>
            <ContinuousView
              title={`${gene.gene} expression`}
              description={`Each dot is a cell. Grey = ${gene.gene} not detected; light to dark teal = higher log-normalized expression. ${gene.n_expressing.toLocaleString("en-US")} of ${results.cells.barcode.length.toLocaleString("en-US")} cells (${percent(gene.n_expressing, results.cells.barcode.length)}) express it.`}
              values={gene.values}
              results={results}
              greyZeros
              colorTitle="log-normalized"
            />
            <GeneByCluster gene={gene} cColors={cColors} />
          </>
        )}
      </section>
    </div>
  );
}

/** Details of the highlighted cluster: size, prediction, top markers (click to show a gene). */
function ClusterDetails({
  cluster,
  results,
  prediction,
  onGene,
  onNavigate,
}: {
  cluster: string | null;
  results: Results;
  prediction: { predicted_label: string; status: string } | undefined;
  onGene: (gene: string) => void;
  onNavigate: (page: Page) => void;
}) {
  if (!cluster) return <DetailsBox>Hover or click a cluster to see its details here.</DetailsBox>;
  const n = results.clusters.find((c) => c.cluster === cluster)?.n_cells ?? 0;
  const markers = results.top_markers?.[cluster];
  return (
    <DetailsBox>
      <span className="text-[15px] font-extrabold text-ink">Cluster {cluster}</span>
      <span className="font-mono text-[13px] text-muted">
        {n.toLocaleString("en-US")} cells · {percent(n, results.cell_cluster.length)}
      </span>
      {prediction && (
        <span className="text-[13.5px]">
          Predicted: <b>{prediction.predicted_label}</b>{" "}
          <span className={prediction.status === "confident" ? "text-muted" : "font-bold text-accent-text"}>
            ({prediction.status === "confident" ? STATUS_TEXT.confident : `▲ ${STATUS_TEXT[prediction.status as keyof typeof STATUS_TEXT] ?? prediction.status}`})
          </span>
        </span>
      )}
      {markers && <MarkerButtons genes={markers} onGene={onGene} />}
      <div className="flex flex-wrap gap-2 pt-1">
        <Button variant="secondary" size="small" onClick={() => goTo(onNavigate, "marker_genes")}>
          Open Marker Genes <ArrowRight aria-hidden="true" size={13} />
        </Button>
        {prediction && (
          <Button variant="secondary" size="small" onClick={() => goTo(onNavigate, "cell_type_annotation")}>
            Open Cell Type Annotation <ArrowRight aria-hidden="true" size={13} />
          </Button>
        )}
      </div>
    </DetailsBox>
  );
}

function CellTypeView({
  results,
  cColors,
  onGene,
  onNavigate,
}: {
  results: Results;
  cColors: Record<string, string>;
  onGene: (gene: string) => void;
  onNavigate: (page: Page) => void;
}) {
  const p = results.predictions!;
  const { groups, colors, clustersOf } = useMemo(() => {
    const counts: Record<string, number> = {};
    p.cell_label.forEach((l) => (counts[l] = (counts[l] ?? 0) + 1));
    const labels = Object.keys(counts).sort((a, b) => {
      const ia = a === p.unresolved_label ? Infinity : p.label_order.indexOf(a);
      const ib = b === p.unresolved_label ? Infinity : p.label_order.indexOf(b);
      return ia - ib;
    });
    const of: Record<string, string[]> = {};
    p.clusters.forEach((c) => (of[c.predicted_label] = [...(of[c.predicted_label] ?? []), c.cluster]));
    const colorOf = Object.fromEntries(
      labels.map((l) => [l, l === p.unresolved_label ? UNRESOLVED_COLOR : cColors[of[l][0]]]),
    );
    return { groups: labels.map((l): UmapGroup => ({ key: l, label: l, n: counts[l] })), colors: colorOf, clustersOf: of };
  }, [p, cColors]);

  return (
    <GroupedUmap
      title="Predicted cell types on the UMAP"
      description="Colored by predicted cell type (Unresolved in grey). Predictions to review."
      listHeading="Cell type"
      groups={groups}
      colors={colors}
      cellGroup={p.cell_label}
      umap={results.umap}
      cells={results.cells}
      extraHover={(i) => `Cluster ${results.cell_cluster[i]}`}
      details={(label) =>
        !label ? (
          <DetailsBox>Hover or click a cell type to see its clusters and markers here.</DetailsBox>
        ) : (
          <DetailsBox>
            <span className="text-[15px] font-extrabold text-ink">{label}</span>
            <span className="text-[13.5px]">
              {clustersOf[label]?.length === 1 ? "Cluster" : "Clusters"} {clustersOf[label]?.join(", ")}
            </span>
            {clustersOf[label]?.map((c) =>
              results.top_markers?.[c] ? (
                <div key={c} className="flex flex-col gap-1">
                  <span className="text-[12.5px] text-muted">Top markers of cluster {c}</span>
                  <MarkerButtons genes={results.top_markers[c]} onGene={onGene} />
                </div>
              ) : null,
            )}
            <div className="pt-1">
              <Button variant="secondary" size="small" onClick={() => goTo(onNavigate, "cell_type_annotation")}>
                Open Cell Type Annotation <ArrowRight aria-hidden="true" size={13} />
              </Button>
            </div>
          </DetailsBox>
        )
      }
    />
  );
}

/** One continuous value per cell (a gene or a QC value) on the UMAP. */
function ContinuousView({
  title,
  description,
  values,
  results,
  greyZeros,
  colorTitle,
}: {
  title: string;
  description: string;
  values: (number | null)[];
  results: Results;
  greyZeros: boolean;
  colorTitle: string;
}) {
  const chart = useMemo(() => {
    const c = results.cells;
    const idx = values.map((_, i) => i);
    const zero = greyZeros ? idx.filter((i) => !values[i]) : [];
    // Expressing cells drawn on top, highest last, so they aren't hidden under others.
    const shown = (greyZeros ? idx.filter((i) => values[i]) : idx).sort((a, b) => (values[a] ?? 0) - (values[b] ?? 0));
    const custom = (i: number) => [c.barcode[i], results.cell_cluster[i], values[i] ?? 0];
    const hover = `<b>%{customdata[0]}</b><br>Cluster %{customdata[1]}<br>${title.replace(" expression", "")}: %{customdata[2]:.2f}<extra></extra>`;
    const data: Data[] = [];
    if (zero.length) {
      data.push({
        type: "scattergl", mode: "markers", name: "not detected",
        x: zero.map((i) => results.umap.x[i]), y: zero.map((i) => results.umap.y[i]),
        customdata: zero.map(custom), marker: { size: 4, color: NOT_EXPRESSED }, hovertemplate: hover,
      });
    }
    data.push({
      type: "scattergl", mode: "markers", name: title,
      x: shown.map((i) => results.umap.x[i]), y: shown.map((i) => results.umap.y[i]),
      customdata: shown.map(custom),
      marker: {
        size: 4.5, color: shown.map((i) => values[i] ?? 0), colorscale: TEAL_SCALE, showscale: true,
        colorbar: { title: { text: colorTitle, side: "right" }, thickness: 12, outlinewidth: 0 },
        line: { color: "rgba(17,17,17,0.25)", width: 0.3 },
      },
      hovertemplate: hover,
    });
    const layout: Partial<Layout> = {
      xaxis: { title: { text: "UMAP 1" }, showticklabels: false, showgrid: false, zeroline: false },
      yaxis: { title: { text: "UMAP 2" }, showticklabels: false, showgrid: false, zeroline: false, scaleanchor: "x" },
      margin: { l: 40, r: 16, t: 16, b: 40 },
      hoverlabel: { bgcolor: "#ffffff", bordercolor: "#d4d4d8", font: { color: "#111111" } },
    };
    return { data, layout };
  }, [values, results, greyZeros, title, colorTitle]);

  return (
    <>
      <div className="border-b border-ink pb-[10px]">
        <h2 className="m-0 text-lg font-extrabold">{title}</h2>
      </div>
      <p className="m-0 text-[13px] text-body">{description}</p>
      <PlotlyChart
        data={chart.data}
        layout={chart.layout}
        height={540}
        label={`UMAP of all cells colored by ${title}`}
        download={{
          name: `Explore Results UMAP ${title}`,
          title: `UMAP of all cells colored by ${title}${greyZeros ? " (grey = not detected)" : ""}`,
        }}
      />
    </>
  );
}

/** Mean expression and share of cells expressing the gene, per cluster. */
function GeneByCluster({ gene, cColors }: { gene: Extract<GenePayload, { ok: true }>; cColors: Record<string, string> }) {
  const chart = useMemo(() => {
    const rows = gene.by_cluster;
    const data: Data[] = [
      {
        type: "bar",
        x: rows.map((r) => `Cluster ${r.cluster}`),
        y: rows.map((r) => r.mean_expr),
        text: rows.map((r) => `${Math.round(100 * r.frac_expressing)}%`),
        textposition: "outside",
        marker: { color: rows.map((r) => cColors[r.cluster]), line: { color: "rgba(17,17,17,0.35)", width: 0.6 } },
        hovertemplate: "%{x}<br>Mean expression %{y:.2f}<br>%{text} of cells express it<extra></extra>",
      },
    ];
    const layout: Partial<Layout> = {
      margin: { l: 56, r: 16, t: 24, b: 60 },
      xaxis: { tickangle: -40 },
      yaxis: { title: { text: "Mean expression" }, gridcolor: "#ececef", rangemode: "tozero" },
    };
    return { data, layout };
  }, [gene, cColors]);
  return (
    <>
      <div className="mt-2 border-b border-ink pb-[10px]">
        <h2 className="m-0 text-lg font-extrabold">{gene.gene} by cluster</h2>
      </div>
      <p className="m-0 text-[13px] text-body">Bar = mean expression; label = share of cells expressing it.</p>
      <PlotlyChart
        data={chart.data}
        layout={chart.layout}
        height={300}
        label={`${gene.gene} mean expression per cluster`}
        download={{ name: `${gene.gene} by cluster`, title: `${gene.gene}: mean expression per cluster (label = % of cells expressing it)` }}
      />
    </>
  );
}

function MarkerButtons({ genes, onGene }: { genes: string[]; onGene: (gene: string) => void }) {
  return (
    <span className="flex flex-wrap gap-1">
      {genes.map((g) => (
        <button
          key={g}
          type="button"
          title={`Color the UMAP by ${g}`}
          onClick={() => onGene(g)}
          className="rounded-[3px] bg-chip px-[6px] py-[1px] font-mono text-[12.5px] text-ink hover:bg-line"
        >
          {g}
        </button>
      ))}
    </span>
  );
}

function DetailsBox({ children }: { children: ReactNode }) {
  return <div className="flex flex-col gap-2 rounded-[4px] border border-line bg-panel p-3 text-[13px] text-muted">{children}</div>;
}

/** The analyst's reviewed labels on the UMAP (clusters not reviewed in grey). */
function ReviewedView({ results, cColors, clusterOrder }: { results: Results; cColors: Record<string, string>; clusterOrder: string[] }) {
  const r = useMemo(
    () => reviewedGroups(results.cell_cluster, clusterOrder, results.reviewed ?? {}, cColors, results.not_reviewed_label),
    [results, cColors, clusterOrder],
  );
  const predicted = useMemo(
    () => Object.fromEntries((results.predictions?.clusters ?? []).map((c) => [c.cluster, c.predicted_label])),
    [results],
  );
  return (
    <GroupedUmap
      title="Your labels on the UMAP"
      description="Your reviewed labels (Not reviewed in grey). Predictions are kept separately."
      listHeading="Your label"
      groups={r.groups}
      colors={r.colors}
      cellGroup={r.cellLabel}
      umap={results.umap}
      cells={results.cells}
      extraHover={(i) => `Cluster ${results.cell_cluster[i]} · predicted: ${predicted[results.cell_cluster[i]] ?? "–"}`}
    />
  );
}
