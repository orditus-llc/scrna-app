// Step 7: cell types predicted from marker evidence. Settings (thresholds, read-only marker
// sets) on the left; on the right the UMAP by predicted cell type, one row per cluster with
// its evidence, and a dot plot of the marker genes. Predictions only: analyst-reviewed labels
// are never part of the app.
import type { Data, Layout } from "plotly.js-dist-min";
import { useMemo, useState, type ReactNode } from "react";

import { Note, RuledField } from "../components/fields";
import { GroupedUmap, percent } from "../components/GroupedUmap";
import { MoreSection } from "../components/MoreSection";
import { SidePanel } from "../components/SidePanel";
import { PlotlyChart } from "../components/PlotlyChart";
import { Tabs } from "../components/Tabs";
import { Prerequisite, SettingsAndResults, StepHeader, SummaryStrip } from "../components/StepPage";
import { clusterColors } from "../lib/clusterColors";
import { cn } from "../lib/utils";
import { goTo } from "../navigation";
import type {
  AnnotationBoxes,
  AnnotationPayload,
  AnnotationSettings,
  ClusterPrediction,
  EvidenceRow,
  Step,
} from "../types";
import { Empty, fmt, Loading, prerequisiteKey, SettingsReset, StepBar, sameAs, SectionTitle, stepOf, useRunButton, type PageProps } from "./stepShared";

const UNRESOLVED_COLOR = "#a1a1aa";
const TEAL_SCALE: [number, string][] = [
  [0, "#f0fdfa"],
  [0.35, "#5eead4"],
  [0.7, "#0d9488"],
  [1, "#134e4a"],
];
export const STATUS_TEXT: Record<ClusterPrediction["status"], string> = {
  confident: "Confident",
  ambiguous: "Ambiguous",
  unsupported: "Unsupported",
  low_score: "Low score",
};

/** The server's reasons for "Unresolved" ("ambiguous, low_score") in words ("Ambiguous, Low score"). */
function reasonsText(reasons: string): string {
  return reasons
    .split(/,\s*/)
    .map((r) => STATUS_TEXT[r as ClusterPrediction["status"]] ?? r)
    .join(", ");
}

export function AnnotationPage({ workflow, data, onNavigate }: PageProps & { data: AnnotationPayload | undefined }) {
  const step = stepOf(workflow, "cell_type_annotation");
  if (!data || !step) return <Loading />;
  return <AnnotationContent workflow={workflow} data={data} step={step} onNavigate={onNavigate} />;
}

function AnnotationContent({ workflow, data, step, onNavigate }: PageProps & { data: AnnotationPayload; step: Step }) {
  const reference = data.settings.reference;
  const [boxes, setBoxes] = useState<AnnotationBoxes>(() => toBoxes(data.settings.current));
  const set = (change: Partial<AnnotationBoxes>) => setBoxes({ ...boxes, ...change });
  const run = useRunButton<AnnotationBoxes>("run_annotation", "cell_type_annotation", data, workflow.notice);
  const changedSinceRun = data.results !== null && !sameBoxes(boxes, toBoxes(data.settings.current));
  const r = data.results;

  let summary = null;
  if (r) {
    const confident = r.clusters.filter((c) => c.status === "confident");
    const unresolved = r.clusters.length - confident.length;
    const types = new Set(confident.map((c) => c.predicted_label));
    const cellsConfident = confident.reduce((sum, c) => sum + c.n_cells, 0);
    summary = (
      <SummaryStrip
        items={[
          { label: "Confident clusters", value: `${confident.length} of ${r.clusters.length}`, note: `${types.size} cell types found` },
          {
            label: "Unresolved clusters",
            value: String(unresolved),
            note: unresolved === 0 ? "none" : "need your review",
            differs: unresolved > 0,
          },
          {
            label: "Cells with a confident label",
            value: fmt(cellsConfident),
            note: `${percent(cellsConfident, r.cell_label.length)} of ${fmt(r.cell_label.length)} cells`,
          },
        ]}
      />
    );
  }

  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <StepHeader step={step}>
        Proposes a cell type for each cluster from its marker genes. These are predictions to review; unclear clusters stay Unresolved.
      </StepHeader>
      {step.prerequisite && <Prerequisite title={step.prerequisite} onOpen={() => goTo(onNavigate, prerequisiteKey(workflow, step))} />}
      {step.reason && <Note dot="bg-accent">{step.reason}</Note>}
      {summary}

      <SettingsAndResults
        settings={
          <>
            <SectionTitle>Settings</SectionTitle>
            <p className="m-0 pt-3 text-[13.5px] text-body">The validated checks and PBMC marker sets are used. Open to review or change them.</p>
            <MoreSection title="Checks" note={sameBoxes(boxes, toBoxes(reference)) ? undefined : "changed"}>
            <RuledField
              id="an_margin"
              label="Ambiguity margin"
              help="How far the best cell type must lead the second-best. Closer calls are marked Ambiguous."
              value={boxes.ambiguity_margin}
              onChange={(v) => set({ ambiguity_margin: v })}
              step={0.1}
              reference={reference.ambiguity_margin}
            />
            <RuledField
              id="an_rank"
              label="Max DE rank"
              help="A marker of the best cell type must rank in the cluster's top this-many genes. Otherwise it's Unsupported."
              value={boxes.max_de_rank}
              onChange={(v) => set({ max_de_rank: v })}
              step={5}
              min={1}
              reference={reference.max_de_rank}
            />
            <RuledField
              id="an_min_score"
              label="Min score"
              help="Clusters scoring below this are marked Low score. Off in the reference."
              value={boxes.min_score}
              onChange={(v) => set({ min_score: v })}
              step={0.1}
              min={-10}
              reference={reference.min_score}
              off={boxes.min_score_off}
              onOff={(off) => set({ min_score_off: off })}
            />
            <RuledField
              id="an_conflict"
              label="Conflict z-score"
              help="Other cell types' markers at or above this z-score are listed as conflicting. It never changes a label."
              value={boxes.conflict_z}
              onChange={(v) => set({ conflict_z: v })}
              step={0.1}
              min={-10}
              reference={reference.conflict_z}
            />
            </MoreSection>
            <MoreSection title="Marker sets (read-only)">
              <MarkerSets data={data} />
            </MoreSection>
            <SettingsReset onReset={() => setBoxes(toBoxes(reference))} changedSinceRun={changedSinceRun} />
          </>
        }
        results={
          r ? (
            <AnnotationResults data={data} results={r} />
          ) : (
            <>
              <SectionTitle>Predicted cell types</SectionTitle>
              <Empty>Run Cell Type Annotation to predict a cell type for each cluster.</Empty>
            </>
          )
        }
      />
      <StepBar
        label="Cell Type Annotation"
        nextTitle="Explore Results"
        status={step.status}
        blocked={step.prerequisite !== null}
        busy={run.busy}
        busyText="Running Cell Type Annotation…"
        error={run.error}
        changedSinceRun={changedSinceRun}
        willClear={data.will_clear}
        onRun={() => run.click(boxes)}
        onContinue={() => goTo(onNavigate, "explorer")}
        prerequisite={step.prerequisite}
        confirmation={step.confirmation}
      />
    </div>
  );
}

function MarkerSets({ data }: { data: AnnotationPayload }) {
  const gene = (g: { gene: string; present: boolean | null }) => (
    <span
      key={g.gene}
      title={g.present === false ? "Not in the data: left out of the score" : undefined}
      className={cn(
        "rounded-[3px] px-[6px] py-[1px] font-mono text-[12.5px]",
        g.present === false ? "bg-accent-tint text-accent-text line-through" : "bg-chip text-ink",
      )}
    >
      {g.gene}
    </span>
  );
  return (
    <div className="flex flex-col border-b border-line py-2">
      <ul className="m-0 flex list-none flex-col gap-[7px] p-0 py-1">
        {data.marker_sets.sets.map((s) => (
          <li key={s.cell_type} className="flex flex-wrap items-baseline justify-between gap-x-3 gap-y-1 text-[13.5px]">
            <span className="font-bold">{s.cell_type}</span>
            <span className="flex flex-wrap gap-1">{s.markers.map(gene)}</span>
          </li>
        ))}
        {data.marker_sets.context.length > 0 && (
          <li className="flex flex-wrap items-baseline justify-between gap-x-3 gap-y-1 border-t border-line pt-[7px] text-[13.5px]">
            <span className="text-muted">Context (shown, not scored)</span>
            <span className="flex flex-wrap gap-1">{data.marker_sets.context.map(gene)}</span>
          </li>
        )}
      </ul>
      <span className="text-[12.5px] leading-snug text-muted">
        Editing marker sets comes later. Genes not found in the data would be crossed out and left out of the score.
      </span>
    </div>
  );
}

type Results = NonNullable<AnnotationPayload["results"]>;

function AnnotationResults({ data, results }: { data: AnnotationPayload; results: Results }) {
  const clusterOrder = results.clusters.map((c) => c.cluster);
  const cColors = useMemo(() => clusterColors(clusterOrder), [results]); // eslint-disable-line react-hooks/exhaustive-deps
  const [selected, setSelected] = useState<string>(clusterOrder[0]);

  // Cell types in marker-set order, Unresolved last. Each takes the color of its first
  // cluster; Unresolved is grey.
  const { groups, colors } = useMemo(() => {
    const setOrder = data.marker_sets.sets.map((s) => s.cell_type);
    const labels = [...new Set(results.clusters.map((c) => c.predicted_label))].sort((a, b) => {
      const ia = a === data.unresolved_label ? Infinity : setOrder.indexOf(a);
      const ib = b === data.unresolved_label ? Infinity : setOrder.indexOf(b);
      return ia - ib;
    });
    const colorOf: Record<string, string> = {};
    for (const label of labels) {
      const first = results.clusters.find((c) => c.predicted_label === label)!;
      colorOf[label] = label === data.unresolved_label ? UNRESOLVED_COLOR : cColors[first.cluster];
    }
    return {
      groups: labels.map((label) => ({ key: label, label, n: results.label_counts[label] ?? 0 })),
      colors: colorOf,
    };
  }, [results, data, cColors]);

  const [view, setView] = useState<"types" | "clusters" | "dots">("types");
  const [panel, setPanel] = useState<string | null>(null); // the cluster whose evidence is open
  const open = (cluster: string) => {
    setSelected(cluster);
    setPanel(cluster);
  };

  return (
    <>
      <p className="m-0 text-[13.5px] text-body">
        <b className="text-ink">Predicted from marker evidence:</b> proposals to review, not confirmed identities.
      </p>
      {results.notes.map((n) => (
        <Note key={n} dot="bg-accent">
          {n}
        </Note>
      ))}
      <Tabs
        label="Cell type views"
        tabs={[
          { key: "types", label: "Cell types" },
          { key: "clusters", label: "Clusters" },
          { key: "dots", label: "Marker dot plot" },
        ]}
        active={view}
        onChange={setView}
      >
        {view === "types" && (
          <GroupedUmap
            title="Predicted cell types on the UMAP"
            description="Hover or click a cell type to highlight it."
            listHeading="Cell type"
            groups={groups}
            colors={colors}
            cellGroup={results.cell_label}
            umap={results.umap}
            cells={results.cells}
            extraHover={(i) => `Cluster ${results.cell_cluster[i]}`}
          />
        )}
        {view === "clusters" && (
          <>
            <p className="m-0 text-[13px] text-body">Click a cluster to see the evidence behind its prediction.</p>
            <PredictionTable clusters={results.clusters} colors={cColors} selected={selected} onSelect={open} />
            <p className="m-0 text-[12.5px] text-muted">
              Limitation: scores use cluster averages, so a small group of different cells can be missed.
            </p>
          </>
        )}
        {view === "dots" && <MarkerDotPlot results={results} clusterOrder={clusterOrder} selected={selected} onSelect={open} />}
      </Tabs>
      {panel !== null && (
        <SidePanel title={`Cluster ${panel}: evidence`} onClose={() => setPanel(null)}>
          <Evidence results={results} cluster={panel} color={cColors[panel]} />
        </SidePanel>
      )}
    </>
  );
}

function StatusMark({ status }: { status: ClusterPrediction["status"] }) {
  return status === "confident" ? (
    <span className="inline-flex items-center gap-[6px] text-[13px] font-bold whitespace-nowrap text-ink">
      <span aria-hidden="true" className="size-2 rounded-full bg-ok" />
      {STATUS_TEXT[status]}
    </span>
  ) : (
    <span className="inline-flex items-center gap-[6px] text-[13px] font-bold whitespace-nowrap text-accent-text">
      <span aria-hidden="true" className="text-[10px]">▲</span>
      {STATUS_TEXT[status]}
    </span>
  );
}

function PredictionTable({
  clusters,
  colors,
  selected,
  onSelect,
}: {
  clusters: ClusterPrediction[];
  colors: Record<string, string>;
  selected: string;
  onSelect: (cluster: string) => void;
}) {
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  return (
    <div className="overflow-x-auto">
      <table className="w-full min-w-[720px] border-collapse text-[14px]">
        <thead>
          <tr className="border-b border-ink text-left">
            <th scope="col" className={cn(head, "pr-3")}>Cluster</th>
            <th scope="col" className={cn(head, "px-3 text-right")}>Cells</th>
            <th scope="col" className={cn(head, "px-3")}>Predicted cell type</th>
            <th scope="col" className={cn(head, "px-3")}>Status</th>
            <th scope="col" className={cn(head, "px-3")}>Runner-up</th>
            <th scope="col" className={cn(head, "px-3 text-right")}>Score</th>
            <th scope="col" className={cn(head, "pl-3 text-right")}>Margin</th>
          </tr>
        </thead>
        <tbody>
          {clusters.map((c) => (
            <tr
              key={c.cluster}
              onClick={() => onSelect(c.cluster)}
              className={cn(
                "cursor-pointer border-b border-line",
                c.cluster === selected ? "bg-chip" : "hover:bg-panel",
                c.status !== "confident" && "bg-accent-tint",
              )}
            >
              <td className="py-[11px] pr-3">
                <button
                  type="button"
                  aria-pressed={c.cluster === selected}
                  onClick={(e) => {
                    e.stopPropagation();
                    onSelect(c.cluster);
                  }}
                  className="inline-flex items-center gap-2 font-bold"
                >
                  <span aria-hidden="true" className="size-3 rounded-full ring-1 ring-black/25" style={{ background: colors[c.cluster] }} />
                  {c.cluster}
                </button>
              </td>
              <td className="px-3 py-[11px] text-right font-mono">{fmt(c.n_cells)}</td>
              <td className="px-3 py-[11px] font-bold">
                {c.predicted_label}
                {c.status !== "confident" && <span className="block text-[12.5px] font-normal text-muted">best candidate: {c.best_candidate}</span>}
              </td>
              <td className="px-3 py-[11px]">
                <StatusMark status={c.status} />
              </td>
              <td className="px-3 py-[11px] text-body">{c.second_candidate ?? "–"}</td>
              <td className="px-3 py-[11px] text-right font-mono">{c.score.toFixed(2)}</td>
              <td className="py-[11px] pl-3 text-right font-mono">{c.margin === null ? "–" : c.margin.toFixed(2)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

/** Why one cluster got its prediction: every cell type's score, and each marker's numbers. */
function Evidence({ results, cluster, color }: { results: Results; cluster: string; color: string }) {
  const row = results.clusters.find((c) => c.cluster === cluster);
  const scores = results.scores[cluster] ?? {};
  const ordered = Object.entries(scores)
    .filter((e): e is [string, number] => e[1] !== null)
    .sort((a, b) => b[1] - a[1]);
  const chart = useMemo(() => {
    const best = ordered[0]?.[0];
    const data: Data[] = [
      {
        type: "bar",
        orientation: "h",
        x: ordered.map(([, s]) => s),
        y: ordered.map(([ct]) => ct),
        marker: { color: ordered.map(([ct]) => (ct === best ? "#111111" : "#d4d4d8")) },
        hovertemplate: "%{y}: score %{x:.2f}<extra></extra>",
      },
    ];
    const layout: Partial<Layout> = {
      margin: { l: 150, r: 16, t: 8, b: 40 },
      xaxis: { title: { text: "Score (mean z-score of the cell type's markers)" }, gridcolor: "#ececef", zeroline: true, zerolinecolor: "#a1a1aa" },
      yaxis: { autorange: "reversed" },
      bargap: 0.35,
    };
    return { data, layout };
  }, [cluster, results]); // eslint-disable-line react-hooks/exhaustive-deps

  if (!row) return null;
  // The best cell type's markers first, then context markers, then the rest.
  const rank = (e: EvidenceRow) => (e.cell_type === row.best_candidate ? 0 : e.cell_type === "context" ? 1 : 2);
  const evidence = results.evidence.filter((e) => e.cluster === cluster).sort((a, b) => rank(a) - rank(b));
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";

  return (
    <section aria-label={`Evidence for cluster ${cluster}`} className="flex flex-col gap-3">
      <div className="flex flex-wrap items-center gap-3">
        <span aria-hidden="true" className="size-3 rounded-full ring-1 ring-black/25" style={{ background: color }} />
        <h3 className="m-0 text-[16px] font-extrabold">
          Cluster {cluster}: {row.predicted_label}
        </h3>
        <StatusMark status={row.status} />
        <span className="font-mono text-[13px] text-muted">{fmt(row.n_cells)} cells</span>
      </div>
      <dl className="m-0 grid grid-cols-[max-content_minmax(0,1fr)] gap-x-6 gap-y-[6px] text-[13.5px]">
        <dt className="text-muted">Supporting markers</dt>
        <dd className="m-0">{row.supporting_markers || "–"}</dd>
        <dt className="text-muted">Best marker's DE rank</dt>
        <dd className="m-0 font-mono">{row.best_marker_DE_rank ?? "not in the ranking"}</dd>
        <dt className="text-muted">Conflicting markers</dt>
        <dd className="m-0">{row.conflicting_markers || "none"}</dd>
        {row.reasons && (
          <>
            <dt className="text-muted">Why Unresolved</dt>
            <dd className="m-0 font-bold text-accent-text">{reasonsText(row.reasons)}</dd>
          </>
        )}
      </dl>
      <h4 className="m-0 mt-1 text-[13px] font-bold tracking-[0.04em] text-muted uppercase">Score for every cell type</h4>
      <PlotlyChart data={chart.data} layout={chart.layout} height={Math.max(180, 30 * ordered.length + 60)} label={`Cell type scores for cluster ${cluster}`} />
      <h4 className="m-0 mt-1 text-[13px] font-bold tracking-[0.04em] text-muted uppercase">Marker genes in this cluster</h4>
      <div className="overflow-x-auto">
        <table className="w-full min-w-[620px] border-collapse text-[13.5px]">
          <thead>
            <tr className="border-b border-ink text-left">
              <th scope="col" className={cn(head, "pr-3")}>Gene</th>
              <th scope="col" className={cn(head, "px-3")}>Marker of</th>
              <th scope="col" className={cn(head, "px-3 text-right")}>Mean expression</th>
              <th scope="col" className={cn(head, "px-3 text-right")}>% expressing</th>
              <th scope="col" className={cn(head, "px-3 text-right")}>z-score</th>
              <th scope="col" className={cn(head, "px-3 text-right")}>DE rank</th>
              <th scope="col" className={cn(head, "pl-3 text-right")}>Log fold change</th>
            </tr>
          </thead>
          <tbody>
            {evidence.map((e) => (
              <tr key={`${e.gene}-${e.cell_type}`} className={cn("border-b border-line", e.cell_type === row.best_candidate && "font-bold")}>
                <td className="py-[8px] pr-3 font-mono">{e.gene}</td>
                <td className="px-3 py-[8px] font-normal text-body">{e.cell_type === "context" ? "context (not scored)" : e.cell_type}</td>
                <td className="px-3 py-[8px] text-right font-mono">{e.mean_expr.toFixed(2)}</td>
                <td className="px-3 py-[8px] text-right font-mono">{Math.round(100 * e.frac_expressing)}%</td>
                <td className="px-3 py-[8px] text-right font-mono">{e.z_score.toFixed(2)}</td>
                <td className="px-3 py-[8px] text-right font-mono">{e.de_rank === null ? "–" : fmt(e.de_rank)}</td>
                <td className="py-[8px] pl-3 text-right font-mono">{e.logfoldchange === null ? "–" : e.logfoldchange.toFixed(2)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <span className="text-[12.5px] text-muted">Bold = markers of the predicted cell type.</span>
    </section>
  );
}

/** Marker-set genes (rows, grouped by cell type) across clusters (columns). */
function MarkerDotPlot({
  results,
  clusterOrder,
  selected,
  onSelect,
}: {
  results: Results;
  clusterOrder: string[];
  selected: string;
  onSelect: (cluster: string) => void;
}) {
  const values = useMemo(() => {
    const seen = new Set<string>();
    return results.evidence.filter((e) => {
      const key = `${e.cluster}|${e.gene}`;
      if (seen.has(key)) return false;
      seen.add(key);
      return true;
    });
  }, [results]);
  const genes = useMemo(() => [...new Set(values.map((v) => v.gene))], [values]);
  const ownerOf = useMemo(() => Object.fromEntries(values.map((v) => [v.gene, v.cell_type])), [values]);

  const chart = useMemo(() => {
    const maxMean = Math.max(...values.map((v) => v.mean_expr), 1e-9);
    const data: Data[] = [
      {
        type: "scatter",
        mode: "markers",
        x: values.map((v) => `Cluster ${v.cluster}`),
        y: values.map((v) => v.gene),
        customdata: values.map((v) => [v.cluster, ownerOf[v.gene] === "context" ? "context marker" : ownerOf[v.gene], v.mean_expr, 100 * v.frac_expressing, v.de_rank ?? NaN]),
        marker: {
          size: values.map((v) => 3 + 19 * Math.sqrt(v.frac_expressing)),
          color: values.map((v) => v.mean_expr),
          cmin: 0,
          cmax: maxMean,
          colorscale: TEAL_SCALE,
          showscale: true,
          colorbar: { title: { text: "Mean expression", side: "right" }, thickness: 12, outlinewidth: 0 },
          line: { color: "rgba(17,17,17,0.35)", width: 0.6 },
        },
        hovertemplate:
          "<b>%{y}</b> (%{customdata[1]}) in cluster %{customdata[0]}<br>Mean expression %{customdata[2]:.2f}<br>" +
          "%{customdata[3]:.0f}% of cells express it<br>DE rank in this cluster: %{customdata[4]:,}<extra></extra>",
      },
    ];
    const layout: Partial<Layout> = {
      margin: { l: 90, r: 16, t: 16, b: 70 },
      xaxis: { type: "category", categoryorder: "array", categoryarray: clusterOrder.map((c) => `Cluster ${c}`), tickangle: -40, showgrid: false },
      yaxis: { type: "category", categoryorder: "array", categoryarray: [...genes].reverse(), gridcolor: "#f4f4f5", tickfont: { family: "IBM Plex Mono, monospace", size: 12 } },
      shapes: [
        {
          type: "rect", xref: "x", yref: "paper", x0: clusterOrder.indexOf(selected) - 0.5, x1: clusterOrder.indexOf(selected) + 0.5,
          y0: 0, y1: 1, fillcolor: "rgba(17,17,17,0.06)", line: { width: 0 }, layer: "below",
        },
      ],
      hoverlabel: { bgcolor: "#ffffff", bordercolor: "#d4d4d8", font: { color: "#111111" } },
    };
    return { data, layout };
  }, [values, genes, ownerOf, clusterOrder, selected]);

  return (
    <>
      <div className="mt-4 border-b border-ink pb-[10px]">
        <h2 className="m-0 text-lg font-extrabold">Marker genes across clusters</h2>
      </div>
      <Caption>
        Dot size = share of cells expressing; color = mean expression. Click a dot to see that cluster.
      </Caption>
      <PlotlyChart
        data={chart.data}
        layout={chart.layout}
        height={Math.max(320, 22 * genes.length + 110)}
        label="Dot plot of the marker-set genes across clusters"
        onClick={(p) => {
          const c = values[p.pointIndex]?.cluster;
          if (c) onSelect(c);
        }}
      />
    </>
  );
}

function Caption({ children }: { children: ReactNode }) {
  return <p className="m-0 text-[13px] text-body">{children}</p>;
}

function toBoxes(s: AnnotationSettings): AnnotationBoxes {
  return {
    ambiguity_margin: decimal(s.ambiguity_margin),
    max_de_rank: String(s.max_de_rank),
    min_score: s.min_score === null ? "" : String(s.min_score),
    min_score_off: s.min_score === null,
    conflict_z: decimal(s.conflict_z),
  };
}

/** A decimal setting keeps its decimal point ("1.0", not "1"). */
function decimal(v: number): string {
  return Number.isInteger(v) ? v.toFixed(1) : String(v);
}

function sameBoxes(a: AnnotationBoxes, b: AnnotationBoxes): boolean {
  const { min_score: ma, min_score_off: oa, ...restA } = a;
  const { min_score: mb, min_score_off: ob, ...restB } = b;
  return oa === ob && (oa || Number(ma) === Number(mb)) && sameAs(restA, restB);
}
