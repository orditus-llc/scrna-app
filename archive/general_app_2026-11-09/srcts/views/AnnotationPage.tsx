// Step 7: cell types predicted from a marker database (PanglaoDB). Settings on the left:
// species, tissues, the database line and the test's checks. On the right the UMAP by
// predicted cell type, one row per cluster (click for its evidence: the best candidates and
// their matching markers), and a dot plot of the matching markers. Predictions only:
// analyst-reviewed labels are never part of the app.
import { useShinyInput } from "@posit-dev/shinyreact";
import type { Data, Layout } from "plotly.js-dist-min";
import { useMemo, useState, type ReactNode } from "react";

import { CsvDownload } from "../components/DownloadMenu";
import { InfoBox, Note, RuledField } from "../components/fields";
import { GroupedUmap } from "../components/GroupedUmap";
import { InfoTip } from "../components/InfoTip";
import { MoreSection } from "../components/MoreSection";
import { SidePanel } from "../components/SidePanel";
import { PlotlyChart } from "../components/PlotlyChart";
import { Tabs } from "../components/Tabs";
import { Button } from "../components/ui/button";
import { Prerequisite, SettingsAndResults, StepHeader } from "../components/StepPage";
import { clusterColors } from "../lib/clusterColors";
import { reviewedGroups } from "../lib/reviewedGroups";
import { cn } from "../lib/utils";
import { goTo } from "../navigation";
import {
  clickCount,
  type AnnotationBoxes,
  type AnnotationPayload,
  type AnnotationSettings,
  type AnnotationStatus,
  type Candidate,
  type Click,
  type ClusterPrediction,
  type ClusterReview,
  type Step,
} from "../types";
import { Empty, fmt, Loading, prerequisiteKey, SettingsReset, StepBar, sameAs, SectionTitle, stepOf, useRunButton, type PageProps } from "./stepShared";

const UNRESOLVED_COLOR = "#a1a1aa";
const TEAL_SCALE: [number, string][] = [
  [0, "#f0fdfa"],
  [0.35, "#5eead4"],
  [0.7, "#0d9488"],
  [1, "#134e4a"],
];
export const STATUS_TEXT: Record<AnnotationStatus, string> = {
  confident: "Confident",
  too_few_markers: "Too few matching markers",
  not_significant: "Not significant",
  ambiguous: "Ambiguous",
};

/** The server's reasons for "Unresolved" ("not_significant, ambiguous") in words. */
function reasonsText(reasons: string): string {
  return reasons
    .split(/,\s*/)
    .map((r) => STATUS_TEXT[r as AnnotationStatus] ?? r)
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
  const changedSinceRun = data.results !== null && !sameAs(boxes, toBoxes(data.settings.current));
  const r = data.results;
  const db = data.database;
  const meta = data.metadata;
  const speciesKey = data.settings.current.species;
  const byTissue = (speciesKey && db.cell_types[speciesKey]) || {};
  const chosenTypes = new Set(meta.tissues.flatMap((t) => byTissue[t] ?? []));

  let missing: string | null = null;
  if (!db.available) missing = "The marker database isn't available, so cell types can't be predicted.";
  else if (!meta.species) missing = "Confirm the species and tissue on Data Overview, then run.";
  else if (!meta.has_markers) missing = "Cell types can't be predicted for this species (the marker database covers human and mouse).";
  else if (meta.tissues.length === 0) missing = "Choose a tissue on Data Overview, then run.";


  const checks: (keyof AnnotationBoxes)[] = ["top_n_genes", "max_adj_p", "min_overlap", "ambiguity_fold"];
  const checksOf = (b: AnnotationBoxes) => Object.fromEntries(checks.map((k) => [k, b[k]]));
  const referenceBoxes = toBoxes(reference);

  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <StepHeader step={step} />
      {step.prerequisite && <Prerequisite title={step.prerequisite} text={step.prerequisite_text} onOpen={() => goTo(onNavigate, prerequisiteKey(workflow, step))} />}
      {step.reason && <Note dot="bg-accent">{step.reason}</Note>}
      {data.review_cleared && <Note dot="bg-accent">{data.review_cleared}</Note>}

      <SettingsAndResults
        settings={
          <>
            <SectionTitle>Settings</SectionTitle>
            <div className="border-b border-line py-[14px]">
            <InfoBox label="Species, tissue and marker source" className="flex flex-col gap-[6px]">
              <span className="text-[14px] leading-snug">
                <span className="font-bold">Species and tissue</span>{" "}
                <InfoTip label="About species and tissue">
                  Set on Data Overview{meta.confirmed_by === "dataset" ? " (from the dataset's metadata)" : ""}. Only the
                  cell types of these tissues are compared with your clusters.
                </InfoTip>
              </span>
              <span className="text-[14px] text-ink">
                {meta.species ?? "Not confirmed"}
                {meta.tissues.length > 0 && <span className="text-body"> · {meta.tissues.join(", ")}</span>}
              </span>
              <button
                type="button"
                className="self-start text-[13px] font-bold text-ink underline underline-offset-2"
                onClick={() => goTo(onNavigate, "data_overview")}
              >
                Change in Data Overview
              </button>
            {db.available ? (
              <p className="m-0 mt-1 border-t border-field pt-2 text-[14px] text-body">
                Markers: <b className="text-ink">{db.name}</b> ({db.version}) ·{" "}
                {meta.has_markers ? `${fmt(chosenTypes.size)} cell types from these tissues` : "no cell types for this species"}{" "}
                <InfoTip label="About the marker selection">
                  {data.selection.canonical_only ? "Canonical markers" : "Markers"} expressed in more than{" "}
                  {Math.round(100 * data.selection.min_sensitivity)}% of the cell type's cells in this species, as in
                  decoupler's {db.name} tutorial. A cell type is tested only if at least {data.selection.min_markers} of
                  its markers are in your data.
                </InfoTip>
              </p>
            ) : (
              <div className="pt-2">
                <Note dot="bg-accent">{db.problem}</Note>
              </div>
            )}
            </InfoBox>
            </div>
            <MoreSection title="Checks" note={sameAs(checksOf(boxes), checksOf(referenceBoxes)) ? undefined : "changed"}>
              <RuledField
                id="an_top_n"
                label="Top marker genes per cluster"
                help="How many of each cluster's top-ranked marker genes are compared with each cell type's markers. Default: 100, Scanpy's default for comparing marker lists."
                value={boxes.top_n_genes}
                onChange={(v) => set({ top_n_genes: v })}
                step={10}
                min={1}
                reference={reference.top_n_genes}
              />
              <RuledField
                id="an_max_p"
                label="Max adjusted p-value"
                help="How unlikely by chance the best cell type's overlap must be (Fisher's exact test, Benjamini–Hochberg adjusted). Otherwise: Not significant. Default: 0.05."
                value={boxes.max_adj_p}
                onChange={(v) => set({ max_adj_p: v })}
                step={0.01}
                reference={reference.max_adj_p}
              />
              <RuledField
                id="an_min_overlap"
                label="Min matching markers"
                help="The best cell type needs at least this many of its markers among the cluster's top genes. Otherwise: Too few matching markers. Default: the app's own rule."
                value={boxes.min_overlap}
                onChange={(v) => set({ min_overlap: v })}
                step={1}
                min={1}
                reference={reference.min_overlap}
              />
              <RuledField
                id="an_fold"
                label="Lead over the runner-up"
                help="The best cell type's adjusted p-value must be at least this many times smaller than a significant runner-up's. Otherwise: Ambiguous. Default: the app's own rule."
                value={boxes.ambiguity_fold}
                onChange={(v) => set({ ambiguity_fold: v })}
                step={1}
                min={1}
                unit="×"
                reference={reference.ambiguity_fold}
              />
            </MoreSection>
            <SettingsReset
              label="Reset checks to default"
              onReset={() => set(checksOf(referenceBoxes) as Partial<AnnotationBoxes>)}
              changedSinceRun={changedSinceRun}
            />
          </>
        }
        results={
          r ? (
            <AnnotationResults data={data} results={r} />
          ) : (
            <>
              <SectionTitle>Predicted cell types</SectionTitle>
              <Empty>
                Run Cell Type Annotation to predict a cell type for each cluster.
              </Empty>
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
        missing={missing}
        onRun={() => run.click(boxes)}
        onContinue={() => goTo(onNavigate, "explorer")}
        prerequisite={step.prerequisite}
        prerequisiteText={step.prerequisite_text}
        confirmation={step.confirmation}
      />
    </div>
  );
}

type Results = NonNullable<AnnotationPayload["results"]>;

function AnnotationResults({ data, results }: { data: AnnotationPayload; results: Results }) {
  const clusterOrder = results.clusters.map((c) => c.cluster);
  const cColors = useMemo(() => clusterColors(clusterOrder), [results]); // eslint-disable-line react-hooks/exhaustive-deps

  // Cell types in the order of their first cluster, Unresolved last. Each takes the color of
  // its first cluster; Unresolved is grey.
  const { groups, colors } = useMemo(() => {
    const labels = [...new Set(results.clusters.map((c) => c.predicted_label))].sort(
      (a, b) => Number(a === data.unresolved_label) - Number(b === data.unresolved_label),
    );
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

  // The analyst's review: each decision is sent to the server, which records it.
  const reviews = results.review.clusters;
  const [reviewClick, setReviewClick] = useShinyInput<Click<ReviewSend>>("review_cluster", 0);
  const sendReview = (send: ReviewSend) => setReviewClick({ n: clickCount(reviewClick) + 1, settings: send });
  const anyReview = Object.keys(reviews).length > 0;
  const [show, setShow] = useState<"predicted" | "reviewed">("predicted");
  const reviewed = useMemo(
    () =>
      reviewedGroups(
        results.cell_cluster,
        clusterOrder,
        Object.fromEntries(Object.values(reviews).map((r) => [r.cluster, r.label])),
        cColors,
        results.review.not_reviewed_label,
      ),
    [results, cColors], // eslint-disable-line react-hooks/exhaustive-deps
  );
  const showReviewed = show === "reviewed" && anyReview;
  // When nothing matched, one grey "Unresolved" group hides the structure: color by cluster instead.
  const allUnresolved = results.clusters.every((c) => c.predicted_label === data.unresolved_label);

  return (
    <>
      <InfoBox label="About these results" className="flex flex-col gap-1">
        <p className="m-0 text-[14px] text-body">
          <b className="text-ink">Predicted from marker evidence:</b> proposals to review, not confirmed identities.
          Confident means a cluster passed the checks, not that its label is right.
        </p>
        <p className="m-0 text-[13px] text-muted">
          {fmt(results.n_cell_types_tested)} cell types tested, from {results.database.name} ({results.database.version}).
        </p>
      </InfoBox>
      {results.nothing_matched && (
        <div role="status" className="flex items-start gap-[9px] rounded-[4px] border border-accent bg-accent-tint px-3 py-[10px] text-[14px] leading-snug text-ink">
          <span aria-hidden="true" className="mt-[1px] text-[11px] text-accent-text">▲</span>
          <span>
            <b>No cell type matched.</b> {results.nothing_matched}
          </span>
        </div>
      )}
      {results.unresolved_note && (
        <div className="flex flex-wrap items-center gap-x-3 gap-y-2">
          <Note dot="bg-accent">{results.unresolved_note}</Note>
          {view !== "clusters" && (
            <Button variant="secondary" size="small" onClick={() => setView("clusters")}>
              Open the Clusters tab
            </Button>
          )}
        </div>
      )}
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
        right={
          view === "types" && anyReview ? (
            <label className="flex items-center gap-2 text-[13px]">
              Show
              <select
                id="an_show"
                value={show}
                onChange={(e) => setShow(e.target.value as "predicted" | "reviewed")}
                className="rounded-[4px] border border-field bg-white px-2 py-1 text-[13px]"
              >
                <option value="predicted">Predicted cell types</option>
                <option value="reviewed">Your labels</option>
              </select>
            </label>
          ) : undefined
        }
      >
        {view === "types" && !showReviewed && allUnresolved && (
          <GroupedUmap
            title="Clusters on the UMAP"
            description="No cluster has a confident cell type, so cells are colored by cluster. Hover or click a cluster to highlight it."
            listHeading="Cluster (all Unresolved)"
            groups={results.clusters.map((c) => ({ key: c.cluster, label: `Cluster ${c.cluster}`, n: c.n_cells }))}
            colors={cColors}
            cellGroup={results.cell_cluster}
            umap={results.umap}
            cells={results.cells}
            extraHover={(i) => {
              const row = results.clusters.find((x) => x.cluster === results.cell_cluster[i]);
              return row && row.overlap > 0 ? `Unresolved · best candidate: ${row.best_candidate}` : "Unresolved";
            }}
          />
        )}
        {view === "types" && !showReviewed && !allUnresolved && (
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
        {view === "types" && showReviewed && (
          <GroupedUmap
            title="Your labels on the UMAP"
            description="Your reviewed labels (Not reviewed in grey). Review clusters in the Clusters tab."
            listHeading="Your label"
            groups={reviewed.groups}
            colors={reviewed.colors}
            cellGroup={reviewed.cellLabel}
            umap={results.umap}
            cells={results.cells}
            extraHover={(i) => {
              const c = results.cell_cluster[i];
              const predicted = results.clusters.find((x) => x.cluster === c)?.predicted_label;
              return `Cluster ${c} · predicted: ${predicted ?? "–"}`;
            }}
          />
        )}
        {view === "clusters" && (
          <>
            <div className="flex flex-wrap items-center justify-between gap-3">
              <p className="m-0 text-[13px] text-body">
                Click a cluster to see its evidence and to confirm or rename its label.
              </p>
              <CsvDownload
                name="Cell Type Annotation predictions"
                columns={[
                  { key: "cluster", header: "Cluster" },
                  { key: "n_cells", header: "Cells" },
                  { key: "predicted_label", header: "Predicted cell type" },
                  { key: "status", header: "Status" },
                  { key: "best_candidate", header: "Best candidate" },
                  { key: "overlap", header: "Matching markers" },
                  { key: "n_markers", header: "Markers tested" },
                  { key: "adj_p", header: "Adjusted p-value" },
                  { key: "matching_markers", header: "Matching marker genes" },
                  { key: "second_candidate", header: "Runner-up" },
                  { key: "second_overlap", header: "Runner-up matching markers" },
                  { key: "second_adj_p", header: "Runner-up adjusted p-value" },
                  { key: "database", header: "Marker database" },
                  { key: "your_label", header: "Your label" },
                  { key: "review", header: "Your review" },
                  { key: "review_note", header: "Your note" },
                  { key: "reviewed_at", header: "Reviewed at (UTC)" },
                ]}
                rows={() =>
                  results.clusters.map((c) => {
                    const r = reviews[c.cluster];
                    return {
                      ...c,
                      status: STATUS_TEXT[c.status],
                      best_candidate: c.overlap === 0 ? "" : c.best_candidate,
                      database: `${results.database.name} (${results.database.version})`,
                      your_label: r?.label ?? "",
                      review: r ? REVIEW_TEXT[r.decision] : results.review.not_reviewed_label,
                      review_note: r?.note ?? "",
                      reviewed_at: r?.reviewed_at ?? "",
                    };
                  })
                }
              />
            </div>
            <PredictionTable
              clusters={results.clusters}
              colors={cColors}
              selected={panel}
              onSelect={setPanel}
              reviews={reviews}
              onConfirm={(cluster) => sendReview({ cluster, decision: "confirm" })}
            />
            <p className="m-0 text-[13px] text-muted">
              Limitation: each cluster is tested as a whole, so a small group of different cells inside it can be missed.
            </p>
          </>
        )}
        {view === "dots" && <MarkerDotPlot results={results} clusterOrder={clusterOrder} selected={panel} onSelect={setPanel} />}
      </Tabs>
      {panel !== null && (
        <SidePanel title={`Cluster ${panel}: evidence`} onClose={() => setPanel(null)}>
          <Evidence key={panel} results={results} cluster={panel} color={cColors[panel]} review={reviews[panel]} onReview={sendReview} />
        </SidePanel>
      )}
    </>
  );
}

/** What the page sends when the analyst decides about a cluster (app.py, input review_cluster). */
interface ReviewSend {
  cluster: string;
  decision: "confirm" | "rename" | "clear";
  label?: string;
  note?: string;
}

const REVIEW_TEXT: Record<ClusterReview["decision"], string> = { confirmed: "Confirmed", renamed: "Renamed" };

/** The "Your label" cell: the reviewed label and how it was decided, or a quick Confirm. */
function ReviewCell({ review, onConfirm }: { review: ClusterReview | undefined; onConfirm: () => void }) {
  if (!review)
    return (
      <Button
        variant="secondary"
        size="small"
        title="Confirm the predicted cell type (open the cluster to rename it instead)"
        onClick={(e) => {
          e.stopPropagation();
          onConfirm();
        }}
      >
        Confirm
      </Button>
    );
  return (
    <span className="flex flex-col">
      <span className="font-bold text-ink">{review.label}</span>
      <span className={cn("text-[13px]", review.prediction_changed ? "font-bold text-accent-text" : "text-muted")}>
        {REVIEW_TEXT[review.decision]}
        {review.prediction_changed && " · ▲ prediction changed"}
      </span>
    </span>
  );
}

/** The analyst's decision for one cluster, at the top of its evidence: confirm the prediction,
 * or rename it (free text, with the candidate cell types as suggestions), with an optional note.
 * Saved on the server and recorded with the date; never changes the prediction itself. */
function ReviewBox({
  row,
  review,
  candidates,
  onReview,
}: {
  row: ClusterPrediction;
  review: ClusterReview | undefined;
  candidates: Candidate[];
  onReview: (send: ReviewSend) => void;
}) {
  const [label, setLabel] = useState(review?.decision === "renamed" ? review.label : "");
  const [note, setNote] = useState(review?.note ?? "");
  const listId = `review_suggestions_${row.cluster}`;
  const suggestions = [...new Set([row.best_candidate, ...candidates.map((c) => c.cell_type)])];
  const confirmedNow = review?.decision === "confirmed" && !review.prediction_changed;
  return (
    <section aria-label="Your review" className="flex flex-col gap-3 rounded-[4px] border border-line bg-panel p-4">
      <div className="flex flex-wrap items-baseline justify-between gap-2">
        <h4 className="m-0 text-[13px] font-bold tracking-[0.04em] text-muted uppercase">Your review</h4>
        <span className="text-[13px]">
          {review ? (
            <>
              <b className="text-ink">{review.label}</b>{" "}
              <span className="text-muted">
                · {REVIEW_TEXT[review.decision]} {review.reviewed_at.replace("T", " ").replace("+00:00", " UTC")}
              </span>
            </>
          ) : (
            <span className="text-muted">Not reviewed</span>
          )}
        </span>
      </div>
      {review?.prediction_changed && (
        <Note dot="bg-accent">
          You confirmed “{review.label}”, but the prediction is now “{row.predicted_label}”. Confirm it again or rename it.
        </Note>
      )}
      <div className="flex flex-wrap items-center gap-2">
        <Button variant="secondary" size="small" disabled={confirmedNow} onClick={() => onReview({ cluster: row.cluster, decision: "confirm", note })}>
          Confirm “{row.predicted_label}”
        </Button>
        <span className="text-[13px] text-muted">or give your own label:</span>
      </div>
      <div className="flex flex-wrap items-end gap-2">
        <label className="flex min-w-0 flex-[1_1_220px] flex-col gap-1 text-[13px] font-bold">
          Your label
          <input
            id={`review_label_${row.cluster}`}
            list={listId}
            value={label}
            maxLength={60}
            placeholder="e.g. CD14+ monocytes"
            onChange={(e) => setLabel(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Enter" && label.trim()) onReview({ cluster: row.cluster, decision: "rename", label, note });
            }}
            className="rounded-[4px] border border-field bg-white px-[10px] py-2 text-[14px] font-normal"
          />
          <datalist id={listId}>
            {suggestions.map((s) => (
              <option key={s} value={s} />
            ))}
          </datalist>
        </label>
        <Button size="small" disabled={!label.trim()} onClick={() => onReview({ cluster: row.cluster, decision: "rename", label, note })}>
          Save label
        </Button>
      </div>
      <label className="flex flex-col gap-1 text-[13px] font-bold">
        <span>
          Note <span className="font-normal text-muted">(optional, saved with your decision)</span>
        </span>
        <textarea
          id={`review_note_${row.cluster}`}
          value={note}
          maxLength={300}
          rows={2}
          placeholder="Why, e.g. LYZ and S100A8/9 high"
          onChange={(e) => setNote(e.target.value)}
          className="rounded-[4px] border border-field bg-white px-[10px] py-2 text-[14px] font-normal"
        />
      </label>
      {review && (
        <div>
          <Button variant="secondary" size="small" onClick={() => onReview({ cluster: row.cluster, decision: "clear" })}>
            Clear review
          </Button>
        </div>
      )}
    </section>
  );
}

function StatusMark({ status }: { status: AnnotationStatus }) {
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

/** A p-value as a number a reader can take in: 0.012, or 2.3 × 10⁻¹² when it's tiny. */
function PValue({ p }: { p: number }) {
  if (p >= 0.001) return <>{p.toFixed(3)}</>;
  const [mantissa, exponent] = p.toExponential(1).split("e");
  return (
    <>
      {mantissa} × 10<sup>{Number(exponent)}</sup>
    </>
  );
}

/** The same, as plain text (for chart hover boxes). */
function pText(p: number): string {
  return p >= 0.001 ? p.toFixed(3) : p.toExponential(1).replace("e", " × 10^");
}

function PredictionTable({
  clusters,
  colors,
  selected,
  onSelect,
  reviews,
  onConfirm,
}: {
  clusters: ClusterPrediction[];
  colors: Record<string, string>;
  selected: string | null;
  onSelect: (cluster: string) => void;
  reviews: Record<string, ClusterReview>;
  onConfirm: (cluster: string) => void;
}) {
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  return (
    <div className="overflow-x-auto">
      <table className="w-full min-w-[760px] border-collapse text-[14px]">
        <thead>
          <tr className="border-b border-ink text-left">
            <th scope="col" className={cn(head, "pr-3")}>Cluster</th>
            <th scope="col" className={cn(head, "px-3 text-right")}>Cells</th>
            <th scope="col" className={cn(head, "px-3")}>Predicted cell type</th>
            <th scope="col" className={cn(head, "px-3")}>Your label</th>
            <th scope="col" className={cn(head, "px-3")}>Status</th>
            <th scope="col" className={cn(head, "px-3 text-right")}>Matching markers</th>
            <th scope="col" className={cn(head, "px-3 text-right")}>Adjusted p-value</th>
            <th scope="col" className={cn(head, "pl-3")}>Runner-up</th>
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
                {c.status !== "confident" && (
                  <span className="block text-[13px] font-normal text-muted">
                    {c.overlap === 0 ? "no matching markers" : `best candidate: ${c.best_candidate}`}
                  </span>
                )}
              </td>
              <td className="px-3 py-[11px]">
                <ReviewCell review={reviews[c.cluster]} onConfirm={() => onConfirm(c.cluster)} />
              </td>
              <td className="px-3 py-[11px]">
                <StatusMark status={c.status} />
              </td>
              <td className="px-3 py-[11px] text-right font-mono whitespace-nowrap">
                {c.overlap} of {c.n_markers}
              </td>
              <td className="px-3 py-[11px] text-right font-mono whitespace-nowrap">
                <PValue p={c.adj_p} />
              </td>
              <td className="py-[11px] pl-3 text-body">{c.second_candidate ?? "–"}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

/** Why one cluster got its prediction: its best candidate cell types and their matching markers. */
function Evidence({
  results,
  cluster,
  color,
  review,
  onReview,
}: {
  results: Results;
  cluster: string;
  color: string;
  review: ClusterReview | undefined;
  onReview: (send: ReviewSend) => void;
}) {
  const row = results.clusters.find((c) => c.cluster === cluster);
  const candidates = results.candidates[cluster] ?? [];
  const settings = results.settings;
  const chart = useMemo(() => {
    const data: Data[] = [
      {
        type: "bar",
        orientation: "h",
        x: candidates.map((c) => c.score),
        y: candidates.map((c) => c.cell_type),
        customdata: candidates.map((c) => [c.overlap, c.n_markers, pText(c.adj_p)]),
        marker: { color: candidates.map((_, i) => (i === 0 ? "#111111" : "#d4d4d8")) },
        hovertemplate: "<b>%{y}</b><br>%{customdata[0]} of %{customdata[1]} markers match<br>Adjusted p-value %{customdata[2]}<extra></extra>",
      },
    ];
    const cut = -Math.log10(settings.max_adj_p);
    const layout: Partial<Layout> = {
      margin: { l: 170, r: 16, t: 24, b: 44 },
      xaxis: { title: { text: "−log₁₀ adjusted p-value (longer = stronger)" }, gridcolor: "#ececef", rangemode: "tozero" },
      yaxis: { autorange: "reversed" },
      bargap: 0.35,
      shapes: [{ type: "line", xref: "x", yref: "paper", x0: cut, x1: cut, y0: 0, y1: 1, line: { color: "#0d9488", width: 1.5, dash: "dash" } }],
      annotations: [
        {
          xref: "x", yref: "paper", x: cut, y: 1, yanchor: "bottom", showarrow: false,
          text: `p = ${settings.max_adj_p}`, font: { size: 11, color: "#0f766e" },
        },
      ],
    };
    return { data, layout };
  }, [cluster, results]); // eslint-disable-line react-hooks/exhaustive-deps

  if (!row) return null;
  const markers = (row.matching_markers || "").split(", ").filter(Boolean);
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";

  return (
    <section aria-label={`Evidence for cluster ${cluster}`} className="flex flex-col gap-3">
      <div className="flex flex-wrap items-center gap-3">
        <span aria-hidden="true" className="size-3 rounded-full ring-1 ring-black/25" style={{ background: color }} />
        <h3 className="m-0 text-[16px] font-bold">
          Cluster {cluster}: {row.predicted_label}
        </h3>
        <StatusMark status={row.status} />
        <span className="font-mono text-[13px] text-muted">{fmt(row.n_cells)} cells</span>
      </div>
      <ReviewBox row={row} review={review} candidates={candidates} onReview={onReview} />
      <dl className="m-0 grid grid-cols-[max-content_minmax(0,1fr)] gap-x-6 gap-y-[6px] text-[14px]">
        <dt className="text-muted">Best candidate</dt>
        <dd className="m-0 font-bold">{row.overlap === 0 ? "none: no cell type has a matching marker" : row.best_candidate}</dd>
        <dt className="text-muted">Matching markers</dt>
        <dd className="m-0">
          {row.overlap} of its {row.n_markers} markers are among the cluster's top {fmt(settings.top_n_genes)} genes
          {markers.length > 0 && (
            <span className="mt-1 flex flex-wrap gap-1">
              {markers.map((g) => (
                <span key={g} className="rounded-[3px] bg-chip px-[6px] py-[1px] font-mono text-[13px]">
                  {g}
                </span>
              ))}
            </span>
          )}
        </dd>
        <dt className="text-muted">Runner-up</dt>
        <dd className="m-0">
          {row.second_candidate === null ? (
            "none"
          ) : (
            <>
              {row.second_candidate} ({row.second_overlap} matching, adjusted p-value <PValue p={row.second_adj_p!} />)
            </>
          )}
        </dd>
        {row.reasons && (
          <>
            <dt className="text-muted">Why Unresolved</dt>
            <dd className="m-0 flex flex-col gap-1">
              <span className="font-bold text-accent-text">{reasonsText(row.reasons)}</span>
              {results.guidance[cluster] && <span className="text-[13px] text-body">{results.guidance[cluster]}</span>}
            </dd>
          </>
        )}
      </dl>
      <h4 className="m-0 mt-1 text-[13px] font-bold tracking-[0.04em] text-muted uppercase">
        Best {candidates.length} of {fmt(results.n_cell_types_tested)} cell types tested
      </h4>
      <PlotlyChart
        data={chart.data}
        layout={chart.layout}
        height={Math.max(180, 34 * candidates.length + 80)}
        label={`Best candidate cell types for cluster ${cluster}`}
        download={{
          name: `Cell Type Annotation cluster ${cluster} candidates`,
          title: `Cluster ${cluster}: best candidate cell types (−log₁₀ adjusted p-value; dashed line = significance level)`,
        }}
      />
      <div className="flex flex-wrap items-center justify-between gap-3">
        <span className="text-[13px] text-muted">Bars past the dashed line pass the significance check.</span>
        <CsvDownload
          name={`Cell Type Annotation cluster ${cluster} candidates`}
          columns={[
            { key: "cell_type", header: "Cell type" },
            { key: "overlap", header: "Matching markers" },
            { key: "n_markers", header: "Markers tested" },
            { key: "adj_p", header: "Adjusted p-value" },
            { key: "matching_markers", header: "Matching marker genes" },
          ]}
          rows={() => candidates.map((c) => ({ ...c }))}
        />
      </div>
      <div className="overflow-x-auto">
        <table className="w-full min-w-[560px] border-collapse text-[14px]">
          <thead>
            <tr className="border-b border-ink text-left">
              <th scope="col" className={cn(head, "pr-3")}>Cell type</th>
              <th scope="col" className={cn(head, "px-3 text-right")}>Matching</th>
              <th scope="col" className={cn(head, "px-3 text-right")}>Adjusted p-value</th>
              <th scope="col" className={cn(head, "pl-3")}>Matching markers</th>
            </tr>
          </thead>
          <tbody>
            {candidates.map((c, i) => (
              <tr key={c.cell_type} className={cn("border-b border-line", i === 0 && "font-bold")}>
                <td className="py-[8px] pr-3">{c.cell_type}</td>
                <td className="px-3 py-[8px] text-right font-mono whitespace-nowrap">
                  {c.overlap} of {c.n_markers}
                </td>
                <td className="px-3 py-[8px] text-right font-mono whitespace-nowrap">
                  <PValue p={c.adj_p} />
                </td>
                <td className="py-[8px] pl-3 font-mono text-[13px] font-normal text-body">{c.matching_markers || "–"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  );
}

/** Matching markers of each cluster's best candidate (rows) across clusters (columns). */
function MarkerDotPlot({
  results,
  clusterOrder,
  selected,
  onSelect,
}: {
  results: Results;
  clusterOrder: string[];
  selected: string | null;
  onSelect: (cluster: string) => void;
}) {
  const { genes, values } = results.dotplot;
  const ownerOf = useMemo(() => Object.fromEntries(genes.map((g) => [g.gene, g.cell_type ?? ""])), [genes]);

  const chart = useMemo(() => {
    const maxMean = Math.max(...values.map((v) => v.mean_expr), 1e-9);
    const data: Data[] = [
      {
        type: "scatter",
        mode: "markers",
        x: values.map((v) => `Cluster ${v.group}`),
        y: values.map((v) => v.gene),
        customdata: values.map((v) => [v.group, ownerOf[v.gene], v.mean_expr, 100 * v.frac_expressing]),
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
          "<b>%{y}</b> (marker of %{customdata[1]}) in cluster %{customdata[0]}<br>Mean expression %{customdata[2]:.2f}<br>" +
          "%{customdata[3]:.0f}% of cells express it<extra></extra>",
      },
    ];
    const at = selected === null ? -1 : clusterOrder.indexOf(selected);
    const layout: Partial<Layout> = {
      margin: { l: 90, r: 16, t: 16, b: 70 },
      xaxis: { type: "category", categoryorder: "array", categoryarray: clusterOrder.map((c) => `Cluster ${c}`), tickangle: -40, showgrid: false },
      yaxis: {
        type: "category",
        categoryorder: "array",
        categoryarray: genes.map((g) => g.gene).reverse(),
        gridcolor: "#f4f4f5",
        tickfont: { family: "IBM Plex Mono, monospace", size: 12 },
      },
      shapes:
        at < 0
          ? []
          : [{ type: "rect", xref: "x", yref: "paper", x0: at - 0.5, x1: at + 0.5, y0: 0, y1: 1, fillcolor: "rgba(17,17,17,0.06)", line: { width: 0 }, layer: "below" }],
      hoverlabel: { bgcolor: "#ffffff", bordercolor: "#d4d4d8", font: { color: "#111111" } },
    };
    return { data, layout };
  }, [values, genes, ownerOf, clusterOrder, selected]);

  if (genes.length === 0) return <Empty>No cluster has matching markers to show.</Empty>;
  return (
    <>
      <div className="mt-4 border-b border-field pb-[10px]">
        <h2 className="m-0 text-[16px] font-bold">Matching markers across clusters</h2>
      </div>
      <Caption>
        Up to 3 matching markers of each cluster's best candidate. Dot size = share of cells expressing; color = mean
        expression. Click a dot to see that cluster's evidence.
      </Caption>
      <PlotlyChart
        data={chart.data}
        layout={chart.layout}
        height={Math.max(320, 22 * genes.length + 110)}
        label="Dot plot of the matching marker genes across clusters"
        download={{
          name: "Cell Type Annotation matching markers dot plot",
          title: "Matching markers of each cluster's best candidate cell type (size = % of cells expressing; color = mean expression)",
        }}
        onClick={(p) => {
          const c = values[p.pointIndex]?.group;
          if (c) onSelect(c);
        }}
      />
    </>
  );
}

function Caption({ children }: { children: ReactNode }) {
  return <p className="m-0 text-[13px] text-body">{children}</p>;
}

/** The checks as boxes (species and tissues come from Data Overview). */
function toBoxes(s: AnnotationSettings): AnnotationBoxes {
  return {
    top_n_genes: String(s.top_n_genes),
    max_adj_p: String(s.max_adj_p),
    min_overlap: String(s.min_overlap),
    ambiguity_fold: String(s.ambiguity_fold),
  };
}
