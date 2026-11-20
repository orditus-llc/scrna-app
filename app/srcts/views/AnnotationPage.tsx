// Step 7: cell types predicted from the CELLxGENE Census marker reference for the sample's tissue,
// with CellTypist's answer beside it. Settings on the left: species and tissue (from Data Overview),
// the reference line, the CellTypist model and the test's checks. On the right the UMAP by
// predicted cell type, one row per cluster (click for its evidence: the best candidates and
// their matching markers), and a dot plot of the matching markers. Predictions only:
// analyst-reviewed labels are never part of the app.
import { useShinyInput } from "@posit-dev/shinyreact";
import type { Data, Layout } from "plotly.js-dist-min";
import { useMemo, useState, type ReactNode } from "react";

import { CsvDownload } from "../components/DownloadMenu";
import { InfoBox, Note, RuledField, RuledSelect } from "../components/fields";
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
import { nextToReview, splitForReview } from "../lib/reviewSplit";
import { cn } from "../lib/utils";
import { goTo } from "../navigation";
import {
  clickCount,
  type AnnotationBoxes,
  type AnnotationPayload,
  type AnnotationSettings,
  type AnnotationStatus,
  type Candidate,
  type CellTypistCluster,
  type Click,
  type ClusterPrediction,
  type ClusterReview,
  type FinalCluster,
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
  no_reference: "No reference yet",
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
  const meta = data.metadata;
  const census = data.census;
  let missing: string | null = null;
  if (!meta.species) missing = "Confirm the species and tissue on Data Overview, then run.";
  else if (!meta.has_markers) missing = "Cell types can't be predicted for this species (the references cover human and mouse).";
  else if (!census.available) missing = "Choose the tissue on Data Overview, then run.";

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
            <InfoBox label="Species, tissue and reference" className="flex flex-col gap-[6px]">
              <span className="text-[14px] leading-snug">
                <span className="font-bold">Species and tissue</span>{" "}
                <InfoTip label="About species and tissue">
                  Set on Data Overview{meta.confirmed_by === "dataset" ? " (from the dataset's metadata)" : ""}. Only the
                  cell types of this tissue are compared with your clusters.
                </InfoTip>
              </span>
              <span className="text-[14px] text-ink">
                {meta.species ?? "Not confirmed"}
                {meta.tissue && <span className="text-body"> · {meta.tissue}</span>}
              </span>
              <button
                type="button"
                className="self-start text-[13px] font-bold text-ink underline underline-offset-2"
                onClick={() => goTo(onNavigate, "data_overview")}
              >
                Change in Data Overview
              </button>
            {census.available && meta.has_markers && (
              <p className="m-0 mt-1 border-t border-field pt-2 text-[14px] text-body">
                Marker check: <b className="text-ink">CELLxGENE Census</b> ({census.release}) ·{" "}
                {census.has_table
                  ? `reference cells from ${census.tissue?.toLowerCase()}`
                  : `no reference for ${census.tissue?.toLowerCase()} yet, so the marker check won't run`}{" "}
                <InfoTip label="About the Census reference">
                  Marker genes computed ahead of time from healthy reference cells in CELLxGENE Census (CC BY 4.0): each cell
                  type with at least 100 cells from 2 or more studies, compared with the cell types it isn't related to in
                  the Cell Ontology. A cell type is tested only if at least {data.selection.min_markers} of its markers are in
                  your data. More tissues are being added.
                </InfoTip>
              </p>
            )}
            </InfoBox>
            </div>
            {meta.has_markers && data.celltypist.choices.length > 0 && (
              <RuledSelect
                id="an_celltypist"
                label="CellTypist model"
                help="A second, independent opinion: a classifier trained on labeled reference cells predicts every cell's type, and each cluster takes the majority vote. Shown beside the marker check; it doesn't change its labels. Cite Domínguez Conde et al. 2022."
                value={boxes.celltypist_model}
                options={[
                  ...data.celltypist.choices.map((c) => ({
                    value: c.model,
                    label: `${c.model.replace(/\.pkl$/, "").replaceAll("_", " ")}${c.downloaded ? "" : " (not downloaded)"}`,
                    group: c.group,
                  })),
                  { value: "", label: census.has_table ? "None (marker check only)" : "None" },
                ]}
                onChange={(v) => set({ celltypist_model: v })}
                reference={data.celltypist.suggested ?? ""}
                note={
                  data.celltypist.choices.find((c) => c.model === boxes.celltypist_model)?.details ??
                  (data.celltypist.suggested ? undefined : "No ready-made model is suggested for this tissue yet.")
                }
              />
            )}
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
            <AnnotationResults results={r} />
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

function AnnotationResults({ results }: { results: Results }) {
  const clusterOrder = results.clusters.map((c) => c.cluster);
  const cColors = useMemo(() => clusterColors(clusterOrder), [results]); // eslint-disable-line react-hooks/exhaustive-deps

  // Final cell types in the order of their first cluster, Unknown last. Each takes the color of
  // its first cluster; Unknown is grey.
  const unknown = results.unknown_label;
  const finalOf = useMemo(() => Object.fromEntries(results.final.map((f) => [f.cluster, f])), [results]);
  const { groups, colors } = useMemo(() => {
    const labels = [...new Set(results.final.map((c) => c.final_label))].sort(
      (a, b) => Number(a === unknown) - Number(b === unknown),
    );
    const colorOf: Record<string, string> = {};
    for (const label of labels) {
      const first = results.final.find((c) => c.final_label === label)!;
      colorOf[label] = label === unknown ? UNRESOLVED_COLOR : cColors[first.cluster];
    }
    return {
      groups: labels.map((label) => ({ key: label, label, n: results.label_counts[label] ?? 0 })),
      colors: colorOf,
    };
  }, [results, unknown, cColors]);

  const [view, setView] = useState<"types" | "review" | "labeled" | "dots">("types");
  const [panel, setPanel] = useState<string | null>(null); // the cluster whose evidence is open

  // The analyst's review: each decision is sent to the server, which records it.
  const reviews = results.review.clusters;
  const [reviewClick, setReviewClick] = useShinyInput<Click<ReviewSend>>("review_cluster", 0);
  const sendReview = (send: ReviewSend) => setReviewClick({ n: clickCount(reviewClick) + 1, settings: send });
  const anyReview = Object.keys(reviews).length > 0;
  // The map shows your label for each reviewed cluster and the suggestion for the rest ("current",
  // analyst 2026-11-16); the suggestions alone or your labels alone are a choice away. The computed
  // suggestions themselves never change.
  const [show, setShow] = useState<"current" | "predicted" | "reviewed">("current");
  const current = useMemo(() => {
    const labels: Record<string, string> = {};
    for (const f of results.final) {
      const mine = reviews[f.cluster]?.label;
      if (mine !== undefined) labels[f.cluster] = mine;
      else if (f.final_label !== unknown) labels[f.cluster] = f.final_label;
    }
    return reviewedGroups(results.cell_cluster, clusterOrder, labels, cColors, unknown);
  }, [results, cColors, unknown]); // eslint-disable-line react-hooks/exhaustive-deps
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
  const showCurrent = show === "current" && anyReview;
  // When no cluster is named, one grey "Unknown" group hides the structure: color by cluster instead.
  const allUnresolved = results.final.every((c) => c.final_label === unknown);
  // The two review tabs: clusters still to label (grouped by how the checks compare) and labeled ones.
  const split = useMemo(() => splitForReview(results.final, reviews, results.review_groups), [results, reviews]);
  const markerOf = useMemo(() => Object.fromEntries(results.clusters.map((c) => [c.cluster, c])), [results]);
  const nLabeled = split.labeled.length;
  const nClusters = results.final.length;
  const csv = (
    <CsvDownload
      name="Cell Type Annotation predictions"
      columns={[
        { key: "cluster", header: "Cluster" },
        { key: "n_cells", header: "Cells" },
        { key: "final_label", header: "Suggested cell type" },
        { key: "confidence", header: "Confidence" },
        { key: "agreement_text", header: "Agreement" },
        { key: "reason", header: "Why" },
        { key: "predicted_label", header: "Marker check" },
        { key: "status", header: "Marker check status" },
        { key: "best_candidate", header: "Best candidate" },
        { key: "overlap", header: "Matching markers" },
        { key: "n_markers", header: "Markers tested" },
        { key: "adj_p", header: "Adjusted p-value" },
        { key: "matching_markers", header: "Matching marker genes" },
        { key: "second_candidate", header: "Runner-up" },
        { key: "second_overlap", header: "Runner-up matching markers" },
        { key: "second_adj_p", header: "Runner-up adjusted p-value" },
        { key: "database", header: "Marker reference" },
        { key: "celltypist", header: "CellTypist (majority vote)" },
        { key: "celltypist_share", header: "CellTypist share of cells" },
        { key: "your_label", header: "Your label" },
        { key: "review", header: "Your review" },
        { key: "review_note", header: "Your note" },
        { key: "reviewed_at", header: "Reviewed at (UTC)" },
      ]}
      rows={() =>
        results.clusters.map((c) => {
          const r = reviews[c.cluster];
          const f = finalOf[c.cluster];
          return {
            ...c,
            final_label: f?.final_label ?? "",
            confidence: f?.confidence_text ?? "",
            agreement_text: f?.agreement_text ?? "",
            reason: f?.reason ?? "",
            status: STATUS_TEXT[c.status],
            best_candidate: c.overlap === 0 ? "" : c.best_candidate,
            database: referenceName(results),
            celltypist: celltypistByCluster(results)?.[c.cluster]?.label ?? "",
            celltypist_share: celltypistByCluster(results)?.[c.cluster]?.share ?? "",
            your_label: r?.label ?? "",
            review: r ? REVIEW_TEXT[r.decision] : results.review.not_reviewed_label,
            review_note: r?.note ?? "",
            reviewed_at: r?.reviewed_at ?? "",
          };
        })
      }
    />
  );
  const tableProps = {
    markers: markerOf,
    celltypist: celltypistByCluster(results),
    colors: cColors,
    selected: panel,
    onSelect: setPanel,
    reviews,
    unknown,
    onConfirm: (cluster: string) => sendReview({ cluster, decision: "confirm" }),
  };

  return (
    <>
      <InfoBox label="About these results" className="flex flex-col gap-1">
        <p className="m-0 flex flex-wrap items-center gap-2 text-[14px] text-body">
          <span>
            <b className="text-ink">Suggestions from two checks.</b> You confirm or rename every cluster.
          </span>
          <InfoTip label="About the suggestions">
            Each cluster is checked twice: its top genes against {results.database ? `${fmt(results.n_cell_types_tested)} cell types' markers in ${results.database.name} (${results.database.version})` : "a marker reference"}, and with
            CellTypist. They are suggestions, not confirmed identities. Each cluster is tested as a whole, so a small group of
            different cells inside it can be missed.
          </InfoTip>
        </p>
        {!results.database && (
          <p className="m-0 text-[13px] text-muted">
            No CELLxGENE Census reference for {(results.settings.census_tissue ?? "this tissue").toLowerCase()} yet, so the marker check didn't run.
            {celltypistByCluster(results) ? " CellTypist's answer is shown for each cluster." : ""}
          </p>
        )}
      </InfoBox>
      {results.nothing_matched && (
        <div role="status" className="flex items-start gap-[9px] rounded-none border border-accent bg-accent-tint px-3 py-[10px] text-[14px] leading-snug text-ink">
          <span aria-hidden="true" className="mt-[1px] text-[11px] text-accent-text">▲</span>
          <span>
            <b>No cell type matched.</b> {results.nothing_matched}
          </span>
        </div>
      )}
      {results.unresolved_note && (
        <div className="flex flex-wrap items-center gap-x-3 gap-y-2">
          <Note dot="bg-accent">{results.unresolved_note}</Note>
        </div>
      )}
      {results.notes.map((n) => (
        <Note key={n} dot="bg-accent">
          {n}
        </Note>
      ))}
      {results.celltypist && "error" in results.celltypist && (
        <Note dot="bg-accent">CellTypist didn't run: {results.celltypist.error}</Note>
      )}
      {results.celltypist && !("error" in results.celltypist) && results.celltypist.notes.map((n) => (
        <Note key={n} dot="bg-accent">
          CellTypist: {n}
        </Note>
      ))}
      <ReviewProgress labeled={nLabeled} total={nClusters} onOpen={split.toReview.length > 0 ? () => { setView("review"); setPanel(split.toReview[0].clusters[0].cluster); } : undefined} />
      <Tabs
        label="Cell type views"
        tabs={[
          { key: "types", label: "Cell types" },
          { key: "review", label: `To review (${nClusters - nLabeled})` },
          { key: "labeled", label: `Labeled clusters (${nLabeled})` },
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
                onChange={(e) => setShow(e.target.value as "current" | "predicted" | "reviewed")}
                className="rounded-none border border-field bg-white px-2 py-1 text-[13px]"
              >
                <option value="current">Your labels, then suggestions</option>
                <option value="predicted">Suggested cell types only</option>
                <option value="reviewed">Your labels only</option>
              </select>
            </label>
          ) : undefined
        }
      >
        {view === "types" && showCurrent && (
          <GroupedUmap
            title="Cell types on the UMAP"
            description="Your label for each cluster you've reviewed; the suggested cell type for the others (Unknown in grey). Hover or click a cell type to highlight it."
            listHeading="Cell type"
            groups={current.groups}
            colors={current.colors}
            cellGroup={current.cellLabel}
            umap={results.umap}
            cells={results.cells}
            extraHover={(i) => {
              const c = results.cell_cluster[i];
              return reviews[c]
                ? `Cluster ${c} · your label (suggested: ${finalOf[c]?.final_label ?? "–"})`
                : `Cluster ${c} · suggestion, not reviewed yet`;
            }}
          />
        )}
        {view === "types" && !showReviewed && !showCurrent && allUnresolved && (
          <GroupedUmap
            title="Clusters on the UMAP"
            description="No cluster could be named, so cells are colored by cluster. Hover or click a cluster to highlight it."
            listHeading={`Cluster (all ${unknown})`}
            groups={results.clusters.map((c) => ({ key: c.cluster, label: `Cluster ${c.cluster}`, n: c.n_cells }))}
            colors={cColors}
            cellGroup={results.cell_cluster}
            umap={results.umap}
            cells={results.cells}
            extraHover={(i) => `${unknown} · ${finalOf[results.cell_cluster[i]]?.agreement_text ?? ""}`}
          />
        )}
        {view === "types" && !showReviewed && !showCurrent && !allUnresolved && (
          <GroupedUmap
            title="Suggested cell types on the UMAP"
            description="Hover or click a cell type to highlight it."
            listHeading="Cell type"
            groups={groups}
            colors={colors}
            cellGroup={results.cell_label}
            umap={results.umap}
            cells={results.cells}
            extraHover={(i) => `Cluster ${results.cell_cluster[i]} · ${finalOf[results.cell_cluster[i]]?.agreement_text ?? ""}`}
          />
        )}
        {view === "types" && showReviewed && (
          <GroupedUmap
            title="Your labels on the UMAP"
            description="Your labels (Not reviewed in grey). Label clusters in the To review tab."
            listHeading="Your label"
            groups={reviewed.groups}
            colors={reviewed.colors}
            cellGroup={reviewed.cellLabel}
            umap={results.umap}
            cells={results.cells}
            extraHover={(i) => {
              const c = results.cell_cluster[i];
              const predicted = finalOf[c]?.final_label;
              return `Cluster ${c} · predicted: ${predicted ?? "–"}`;
            }}
          />
        )}
        {view === "review" && (
          <div className="flex flex-col gap-9">
            {split.toReview.length === 0 ? (
              <Note dot="bg-ok">
                Every cluster is labeled.{" "}
                <button type="button" className="font-bold underline" onClick={() => setView("labeled")}>
                  Open Labeled clusters
                </button>
              </Note>
            ) : (
              split.toReview.map((g) => (
                <section key={g.key} aria-label={g.label} className="flex flex-col gap-3">
                  <h3 className="m-0 flex items-center gap-2 border-b-2 border-ink pb-1.5 text-[16px] font-extrabold text-ink">
                    <span>
                      {g.label} <span className="text-[14px] font-semibold text-muted">· {g.clusters.length} {g.clusters.length === 1 ? "cluster" : "clusters"}</span>
                    </span>
                    {g.description && <InfoTip label={`About ${g.label}`}>{g.description}</InfoTip>}
                  </h3>
                  <ClusterTable mode="review" rows={g.clusters} {...tableProps} />
                </section>
              ))
            )}
          </div>
        )}
        {view === "labeled" && (
          <div className="flex flex-col gap-3">
            <div className="flex flex-wrap items-center justify-end gap-3">{csv}</div>
            {nLabeled === 0 ? (
              <Note dot="bg-accent">No cluster is labeled yet. Confirm or rename clusters in To review.</Note>
            ) : (
              <ClusterTable mode="labeled" rows={split.labeled} {...tableProps} />
            )}
          </div>
        )}
        {view === "dots" && <MarkerDotPlot results={results} clusterOrder={clusterOrder} selected={panel} onSelect={setPanel} />}
      </Tabs>
      {panel !== null && (
        <SidePanel title={`Cluster ${panel}: review`} onClose={() => setPanel(null)}>
          <Evidence
            key={panel}
            results={results}
            final={finalOf[panel]}
            cluster={panel}
            color={cColors[panel]}
            review={reviews[panel]}
            onReview={sendReview}
            next={nextToReview(split.toReview.flatMap((g) => g.clusters.map((c) => c.cluster)), panel)}
            onOpen={setPanel}
          />
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
function ReviewCell({ review, canConfirm, onConfirm }: { review: ClusterReview | undefined; canConfirm: boolean; onConfirm: () => void }) {
  if (!review && !canConfirm) return <span className="text-[13px] text-muted">Open to label</span>;
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

/** The analyst's decision for one cluster: confirm the suggestion, or give your own label (free text,
 * with the checks' answers and the candidate cell types as suggestions), with an optional note. Saved on
 * the server and recorded with the date; never changes the suggestion itself. "Next cluster" moves on in
 * To review order. */
function ReviewBox({
  row,
  finalLabel,
  review,
  candidates,
  onReview,
  extraSuggestions,
  next,
  onOpen,
}: {
  row: ClusterPrediction;
  finalLabel: string;
  review: ClusterReview | undefined;
  candidates: Candidate[];
  onReview: (send: ReviewSend) => void;
  extraSuggestions: string[]; // e.g. CellTypist's answer
  next: string | null;
  onOpen: (cluster: string) => void;
}) {
  const [label, setLabel] = useState(review?.decision === "renamed" ? review.label : "");
  const [note, setNote] = useState(review?.note ?? "");
  const listId = `review_suggestions_${row.cluster}`;
  const suggestions = [...new Set([finalLabel, row.best_candidate, ...extraSuggestions, ...candidates.map((c) => c.cell_type)].filter((x): x is string => !!x && x !== "Unknown"))];
  const confirmedNow = review?.decision === "confirmed" && !review.prediction_changed;
  return (
    <section aria-label="Your review" className="flex flex-col gap-3 rounded-none border border-line bg-panel p-4">
      <div className="flex flex-wrap items-baseline justify-between gap-2">
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
          You confirmed “{review.label}”, but the prediction is now “{finalLabel}”. Confirm it again or rename it.
        </Note>
      )}
      {finalLabel === "Unknown" ? (
        <span className="text-[13px] text-muted">No cell type was predicted. Give your own label:</span>
      ) : (
        <div className="flex flex-wrap items-center gap-2">
          <Button variant="secondary" size="small" disabled={confirmedNow} onClick={() => onReview({ cluster: row.cluster, decision: "confirm", note })}>
            Confirm “{finalLabel}”
          </Button>
          <span className="text-[13px] text-muted">or give your own label:</span>
        </div>
      )}
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
            className="rounded-none border border-field bg-white px-[10px] py-2 text-[14px] font-normal"
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
          className="rounded-none border border-field bg-white px-[10px] py-2 text-[14px] font-normal"
        />
      </label>
      {(review || next) && (
        <div className="flex flex-wrap items-center gap-2">
          {review && (
            <Button variant="secondary" size="small" onClick={() => onReview({ cluster: row.cluster, decision: "clear" })}>
              Clear review
            </Button>
          )}
          {next && (
            <Button variant="secondary" size="small" className="ml-auto" onClick={() => onOpen(next)}>
              Next cluster →
            </Button>
          )}
        </div>
      )}
    </section>
  );
}

/** A p-value as a number a reader can take in: 0.012, or 2.3 × 10⁻¹² when it's tiny. */
function PValue({ p }: { p: number | null }) {
  if (p === null) return <>–</>;
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

/** CellTypist's answer per cluster, or null when it didn't run. */
function celltypistByCluster(results: Results): Record<string, CellTypistCluster> | null {
  const ct = results.celltypist;
  if (!ct || "error" in ct) return null;
  return Object.fromEntries(ct.clusters.map((c) => [c.cluster, c]));
}

/** The final label's agreement: green dot when both checks agree, teal ▲ otherwise. */
function AgreementMark({ final }: { final: FinalCluster }) {
  const agree = final.agreement === "agree";
  return (
    <span className={cn("inline-flex items-center gap-[6px] text-[13px] font-bold whitespace-nowrap", agree ? "text-ink" : "text-accent-text")}>
      {agree ? <span aria-hidden="true" className="size-2 rounded-full bg-ok" /> : <span aria-hidden="true" className="text-[10px]">▲</span>}
      {final.agreement_text}
    </span>
  );
}

/** "Labeled 3 of 8 clusters" with a bar: how far the user's review has come. */
function ReviewProgress({ labeled, total, onOpen }: { labeled: number; total: number; onOpen?: () => void }) {
  const done = total > 0 && labeled === total;
  return (
    <div className="flex flex-wrap items-center gap-x-4 gap-y-2" role="status">
      <span className="text-[14px] text-ink">
        <b>
          Labeled {fmt(labeled)} of {fmt(total)} clusters
        </b>
        {done && <span className="text-muted"> · review complete</span>}
      </span>
      <span
        aria-hidden="true"
        className="h-[6px] min-w-[120px] flex-[1_1_160px] overflow-hidden rounded-none bg-chip"
      >
        <span className={cn("block h-full", done ? "bg-ok" : "bg-ink")} style={{ width: `${total ? (100 * labeled) / total : 0}%` }} />
      </span>
      {onOpen && (
        <Button variant="secondary" size="small" onClick={onOpen}>
          Review next cluster
        </Button>
      )}
    </div>
  );
}

/** How sure the suggestion is: green dot for High, teal ▲ for Medium and Low. */
function ConfidenceMark({ final }: { final: FinalCluster }) {
  const high = final.confidence === "high";
  return (
    <span
      className={cn("inline-flex items-center gap-[6px] text-[13px] font-bold whitespace-nowrap", high ? "text-ink" : "text-accent-text")}
      title={final.agreement_text}
    >
      {high ? <span aria-hidden="true" className="size-2 rounded-full bg-ok" /> : <span aria-hidden="true" className="text-[10px]">▲</span>}
      {final.confidence_text}
    </span>
  );
}

/** One table of clusters. "review": still to label (Confirm + Review). "labeled": the user's labels
 * with the suggestion beside them (Review to change). Both keep the two checks' columns. */
function ClusterTable({
  mode,
  rows,
  markers,
  celltypist,
  colors,
  selected,
  onSelect,
  reviews,
  unknown,
  onConfirm,
}: {
  mode: "review" | "labeled";
  rows: FinalCluster[];
  markers: Record<string, ClusterPrediction>;
  celltypist: Record<string, CellTypistCluster> | null;
  colors: Record<string, string>;
  selected: string | null;
  onSelect: (cluster: string | null) => void;
  reviews: Record<string, ClusterReview>;
  unknown: string;
  onConfirm: (cluster: string) => void;
}) {
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  const review = mode === "review";
  return (
    <div className="overflow-x-auto">
      <table className="w-full min-w-[640px] border-collapse text-[14px]">
        <thead>
          <tr className="border-b border-ink text-left">
            <th scope="col" className={cn(head, "pr-3")}>Cluster</th>
            <th scope="col" className={cn(head, "px-3 text-right")}>Cells</th>
            {!review && <th scope="col" className={cn(head, "px-3")}>Your label</th>}
            <th scope="col" className={cn(head, "px-3")}>Suggested cell type</th>
            <th scope="col" className={cn(head, "px-3")}>Confidence</th>
            <th scope="col" className={cn(head, "px-3")}>Marker check</th>
            {celltypist && <th scope="col" className={cn(head, "px-3")}>CellTypist</th>}
            {!review && <th scope="col" className={cn(head, "px-3")}>Your note</th>}
            <th scope="col" className={cn(head, "pl-3 text-right")}>
              <span className="sr-only">Actions</span>
            </th>
          </tr>
        </thead>
        <tbody>
          {rows.map((f) => {
            const c = markers[f.cluster];
            const r = reviews[f.cluster];
            const suggested = f.final_label !== unknown;
            return (
              <tr
                key={f.cluster}
                onClick={() => onSelect(f.cluster)}
                className={cn("cursor-pointer border-b border-line", f.cluster === selected ? "bg-chip" : "hover:bg-panel")}
              >
                <td className="py-[11px] pr-3">
                  <span className="inline-flex items-center gap-2 font-bold">
                    <span aria-hidden="true" className="size-3 rounded-full ring-1 ring-black/25" style={{ background: colors[f.cluster] }} />
                    {f.cluster}
                  </span>
                </td>
                <td className="px-3 py-[11px] text-right font-mono">{fmt(f.n_cells)}</td>
                {!review && (
                  <td className="px-3 py-[11px]">
                    <ReviewCell review={r} canConfirm={false} onConfirm={() => undefined} />
                  </td>
                )}
                <td className="px-3 py-[11px]">
                  {suggested ? <span className="font-bold">{f.final_label}</span> : <span className="text-muted">No suggestion</span>}
                </td>
                <td className="px-3 py-[11px]">
                  <ConfidenceMark final={f} />
                  {!review && <span className="block text-[13px] text-muted">{f.agreement_text}</span>}
                </td>
                <td className="px-3 py-[11px]">
                  {c?.predicted_label ?? "–"}
                  {c && c.status !== "confident" && (
                    <span className="block text-[13px] font-normal text-muted">
                      {c.status === "no_reference"
                        ? "no reference for this tissue yet"
                        : c.overlap === 0
                          ? "no matching markers"
                          : `best candidate: ${c.best_candidate}`}
                    </span>
                  )}
                </td>
                {celltypist && (
                  <td className="px-3 py-[11px]">
                    {celltypist[f.cluster] ? (
                      <>
                        {celltypist[f.cluster].label}
                        <span className="block text-[13px] whitespace-nowrap text-muted">{Math.round(100 * celltypist[f.cluster].share)}% of cells</span>
                      </>
                    ) : (
                      "–"
                    )}
                  </td>
                )}
                {!review && <td className="max-w-[240px] px-3 py-[11px] text-[13px] text-body">{r?.note || "–"}</td>}
                <td className="py-[11px] pl-3 text-right whitespace-nowrap">
                  <span className="inline-flex items-center gap-2">
                    {review && suggested && (
                      <Button
                        variant="secondary"
                        size="small"
                        title={`Confirm “${f.final_label}” (click Review to rename it instead)`}
                        onClick={(e) => {
                          e.stopPropagation();
                          onConfirm(f.cluster);
                        }}
                      >
                        Confirm
                      </Button>
                    )}
                    <Button
                      size="small"
                      variant={review ? "primary" : "secondary"}
                      aria-pressed={f.cluster === selected}
                      onClick={(e) => {
                        e.stopPropagation();
                        onSelect(f.cluster);
                      }}
                    >
                      Review →
                    </Button>
                  </span>
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

/** The marker check's reference for the page, e.g. "CELLxGENE Census (2025-11-08)". */
function referenceName(results: Results, withVersion = true): string {
  const db = results.database;
  if (!db) return "no reference";
  return withVersion ? `${db.name} (${db.version})` : db.name;
}

/** Why CellTypist's answer can't be used for a cluster, in words (consensus.py's issue codes). */
const CELLTYPIST_ISSUE: Record<string, string> = {
  not_run: "no model for this species and tissue",
  error: "it couldn't run on this data",
  no_table: "this model's labels can't be compared with the marker check yet",
  low_gene_match: "too few of the model's genes are in the data",
  low_share: "fewer than half of the cluster's cells agree",
  no_broad_type: "the label is too vague to compare",
};

/** "How was this identified?": the two checks side by side, each in a few words. */
function ChecksSideBySide({ results, row, final }: { results: Results; row: ClusterPrediction; final: FinalCluster | undefined }) {
  const ct = celltypistByCluster(results)?.[row.cluster];
  const model = results.celltypist && !("error" in results.celltypist) ? results.celltypist.model : null;
  const card = "flex flex-col gap-1 rounded-none border border-line bg-white p-3";
  const title = "text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  return (
    <div className="grid grid-cols-1 gap-2 sm:grid-cols-2">
        <div className={card}>
          <span className={title}>Marker check</span>
          {row.status === "confident" ? (
            <>
              <b className="text-[14px]">{row.predicted_label}</b>
              <span className="text-[13px] text-muted">
                {row.overlap} of {row.n_markers} markers match
              </span>
            </>
          ) : row.status === "no_reference" ? (
            <span className="text-[14px] text-muted">Not run: no reference for this tissue yet</span>
          ) : (
            <>
              <b className="text-[14px]">Unresolved</b>
              <span className="text-[13px] text-muted">
                {STATUS_TEXT[row.status]}
                {row.overlap > 0 && `; closest: ${row.best_candidate}`}
              </span>
            </>
          )}
        </div>
        <div className={card}>
          <span className={title}>CellTypist</span>
          {ct && model ? (
            <>
              <b className="text-[14px]">{ct.label}</b>
              <span className="text-[13px] text-muted">
                {Math.round(100 * ct.share)}% of cells vote for it
                {final && !final.celltypist_confident && final.celltypist_issue && ` · not used: ${CELLTYPIST_ISSUE[final.celltypist_issue] ?? final.celltypist_issue}`}
              </span>
            </>
          ) : (
            <span className="text-[14px] text-muted">
              Not run{final?.celltypist_issue ? `: ${CELLTYPIST_ISSUE[final.celltypist_issue] ?? final.celltypist_issue}` : ""}
            </span>
          )}
        </div>
      </div>
  );
}

/** The cluster's top marker genes, each with the reference cell types it is a known marker of. The
 * cell types this cluster's checks point to are listed first. */
function TopGenes({ results, cluster, final, candidates }: { results: Results; cluster: string; final: FinalCluster | undefined; candidates: Candidate[] }) {
  const genes = results.top_genes[cluster] ?? [];
  if (genes.length === 0) return null;
  const first = new Set([final?.marker_label, ...candidates.map((c) => c.cell_type)].filter(Boolean) as string[]);
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  const pct = (x: number | null) => (x === null ? "–" : `${Math.round(100 * x)}%`);
  const hasReference = !!results.database && !results.no_reference;
  return (
    <PanelSection
      title="Top genes in this cluster"
      info={
        <>
          The genes most specific to this cluster (from Marker Genes). “Known marker of” lists the cell types in{" "}
          {hasReference ? referenceName(results) : "the marker reference"} that have the gene among their markers.
        </>
      }
    >
      <div className="overflow-x-auto">
        <table className="w-full min-w-[460px] border-collapse text-[14px]">
          <thead>
            <tr className="border-b border-ink text-left">
              <th scope="col" className={cn(head, "pr-3")}>Gene</th>
              <th scope="col" className={cn(head, "px-3 text-right")}>Log fold change</th>
              <th scope="col" className={cn(head, "px-3 text-right")}>% in cluster vs rest</th>
              <th scope="col" className={cn(head, "pl-3")}>Known marker of</th>
            </tr>
          </thead>
          <tbody>
            {genes.map((g) => {
              const types = [...g.marker_of.filter((t) => first.has(t)), ...g.marker_of.filter((t) => !first.has(t))];
              const shown = types.slice(0, 2);
              return (
                <tr key={g.gene} className="border-b border-line">
                  <td className="py-[7px] pr-3 font-mono text-[13px] font-bold">{g.gene}</td>
                  <td className="px-3 py-[7px] text-right font-mono text-[13px] whitespace-nowrap">
                    {g.log_fold_change === null ? "–" : g.log_fold_change.toFixed(1)}
                  </td>
                  <td className="px-3 py-[7px] text-right font-mono text-[13px] whitespace-nowrap">
                    {pct(g.frac_in_cluster)} <span className="text-muted">vs {pct(g.frac_in_rest)}</span>
                  </td>
                  <td className="py-[7px] pl-3 text-[13px]" title={types.length > 2 ? types.join("\n") : undefined}>
                    {shown.length === 0 ? (
                      <span className="text-muted">{hasReference ? "–" : "no reference"}</span>
                    ) : (
                      <>
                        {shown.map((t, i) => (
                          <span key={t} className={first.has(t) ? "font-bold text-ink" : "text-body"}>
                            {i > 0 && ", "}
                            {t}
                          </span>
                        ))}
                        {types.length > 2 && <span className="text-muted"> +{types.length - 2} more</span>}
                      </>
                    )}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </PanelSection>
  );
}

/** One section of the review panel: a bold heading over a line, so the parts stand apart. */
function PanelSection({ title, info, aside, children }: { title: string; info?: ReactNode; aside?: ReactNode; children: ReactNode }) {
  return (
    <section aria-label={title} className="flex flex-col gap-3">
      <div className="flex flex-wrap items-center justify-between gap-2 border-b border-ink pb-[6px]">
        <h4 className="m-0 flex items-center gap-2 text-[16px] font-bold text-ink">
          {title}
          {info && <InfoTip label={`About ${title.toLowerCase()}`}>{info}</InfoTip>}
        </h4>
        {aside}
      </div>
      {children}
    </section>
  );
}

/** One cluster's review panel, in sections with bold headings: the suggestion, your review, the two
 * checks side by side, the likely cell types (bar chart) and the top genes; the remaining statistics
 * are folded under "More details". */
function Evidence({
  results,
  final,
  cluster,
  color,
  review,
  onReview,
  next,
  onOpen,
}: {
  results: Results;
  final: FinalCluster | undefined;
  cluster: string;
  color: string;
  review: ClusterReview | undefined;
  onReview: (send: ReviewSend) => void;
  next: string | null;
  onOpen: (cluster: string) => void;
}) {
  const row = results.clusters.find((c) => c.cluster === cluster);
  const candidates = results.candidates[cluster] ?? [];
  const nTested = results.n_cell_types_tested;
  const settings = results.settings;
  const ct = celltypistByCluster(results)?.[cluster];
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
    <section aria-label={`Evidence for cluster ${cluster}`} className="flex flex-col gap-6">
      <div className="flex flex-wrap items-center gap-3">
        <span aria-hidden="true" className="size-3 rounded-full ring-1 ring-black/25" style={{ background: color }} />
        <h3 className="m-0 text-[16px] font-bold">Cluster {cluster}</h3>
        <span className="font-mono text-[13px] text-muted">{fmt(row.n_cells)} cells</span>
        {final && (
          <span className="ml-auto">
            <ConfidenceMark final={final} /> <span className="text-[13px] text-muted">confidence</span>
          </span>
        )}
      </div>

      <PanelSection title="Suggested cell type">
        <div className="flex flex-col gap-1">
          <span className="flex flex-wrap items-center gap-x-3 gap-y-1">
            <b className="text-[16px]">{final ? (final.final_label === "Unknown" ? "No suggestion" : final.final_label) : row.predicted_label}</b>
            {final && <AgreementMark final={final} />}
          </span>
          {final && <span className="text-[13px] text-body">{final.reason}</span>}
        </div>
      </PanelSection>

      <PanelSection title="Your review">
        <ReviewBox
          key={`${cluster}-${review?.reviewed_at ?? ""}`}
          row={row}
          finalLabel={final?.final_label ?? row.predicted_label}
          review={review}
          candidates={candidates}
          onReview={onReview}
          extraSuggestions={[final?.marker_label, final?.celltypist_label].filter((x): x is string => !!x)}
          next={next}
          onOpen={onOpen}
        />
      </PanelSection>

      <PanelSection title="How was this identified?">
        <ChecksSideBySide results={results} row={row} final={final} />
      </PanelSection>

      {candidates.length > 0 && (
        <PanelSection
          title="Marker check: best matches"
          info={
            <>
              The marker check's best {candidates.length} of {fmt(nTested)} cell types in {referenceName(results)}: how strongly
              the cluster's top genes match each one's markers. Bars past the dashed line pass the significance check.
            </>
          }
        >
          <PlotlyChart
            data={chart.data}
            layout={chart.layout}
            height={Math.max(180, 34 * candidates.length + 80)}
            label={`Best candidate cell types for cluster ${cluster}`}
            download={{
              name: `Cell Type Annotation cluster ${cluster} candidates`,
              title: `Cluster ${cluster}: marker check's best matches in ${referenceName(results)} (−log₁₀ adjusted p-value; dashed line = significance level)`,
            }}
          />
        </PanelSection>
      )}

      <PanelSection
        title="Best candidates"
        info={
          <>
            Both answers are compared as broad cell types (for example T cell, monocyte, neuron). The suggestion is a cell
            type when they agree, or when one is confident and the other gives no confident answer.
          </>
        }
      >
          <dl className="m-0 grid grid-cols-[max-content_minmax(0,1fr)] gap-x-6 gap-y-[6px] text-[14px]">
            <dt className="text-muted">Marker check, best candidate</dt>
            <dd className="m-0 font-bold">
              {row.status === "no_reference"
                ? "none: no reference for this tissue yet"
                : row.overlap === 0
                  ? "none: no cell type has a matching marker"
                  : row.best_candidate}
            </dd>
            <dt className="text-muted">Matching markers</dt>
            <dd className="m-0">
              {row.overlap} of its {row.n_markers} markers are among the cluster's top {fmt(settings.top_n_genes)} genes
              {row.adj_p !== null && (
                <>
                  {" "}
                  (adjusted p-value <PValue p={row.adj_p} />)
                </>
              )}
              {markers.length > 0 && (
                <span className="mt-1 flex flex-wrap gap-1">
                  {markers.map((g) => (
                    <span key={g} className="rounded-none bg-chip px-[6px] py-[1px] font-mono text-[13px]">
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
                <dt className="text-muted">Why the marker check is Unresolved</dt>
                <dd className="m-0 flex flex-col gap-1">
                  <span className="font-bold text-accent-text">{reasonsText(row.reasons)}</span>
                  {results.guidance[cluster] && <span className="text-[13px] text-body">{results.guidance[cluster]}</span>}
                </dd>
              </>
            )}
            {ct && (
              <>
                <dt className="text-muted">CellTypist</dt>
                <dd className="m-0">
                  {ct.label}: {Math.round(100 * ct.share)}% of cells, mean probability {ct.mean_probability.toFixed(2)}
                  {ct.runner_up ? `; next: ${ct.runner_up} (${Math.round(100 * (ct.runner_up_share ?? 0))}%)` : ""}
                </dd>
              </>
            )}
          </dl>
      </PanelSection>

      <TopGenes results={results} cluster={cluster} final={final} candidates={candidates} />

      {candidates.length > 0 && (
        // Open at first where the checks disagree or neither is confident (Low confidence), so the full table is in view.
        <MoreSection title="More details" defaultOpen={final?.review_group === "disagree" || final?.review_group === "no_suggestion"}>
          <div className="flex flex-col gap-3">
            <>
              <div className="flex flex-wrap items-center justify-between gap-3">
                <span className="text-[13px] font-bold">Marker check: best matches as a table</span>
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
            </>
          </div>
        </MoreSection>
      )}
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

/** The checks as boxes (species and tissue come from Data Overview). */
function toBoxes(s: AnnotationSettings): AnnotationBoxes {
  return {
    celltypist_model: s.celltypist_model ?? "",
    top_n_genes: String(s.top_n_genes),
    max_adj_p: String(s.max_adj_p),
    min_overlap: String(s.min_overlap),
    ambiguity_fold: String(s.ambiguity_fold),
  };
}
