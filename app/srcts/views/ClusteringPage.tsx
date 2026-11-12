// Step 5: Leiden clustering. Settings on the left; the UMAP colored by cluster and the
// cluster list on the right. Hovering a cell or a cluster highlights that cluster and fades
// the rest; clicking keeps the highlight (Esc or "Clear" removes it).
import { useMemo, useState } from "react";

import { CsvDownload } from "../components/DownloadMenu";
import { InfoBox, Note, RuledField } from "../components/fields";
import { InfoTip } from "../components/InfoTip";
import { GroupedUmap } from "../components/GroupedUmap";
import { MoreSection } from "../components/MoreSection";
import { Tabs } from "../components/Tabs";
import { Prerequisite, SettingsAndResults, StepHeader } from "../components/StepPage";
import { clusterColors } from "../lib/clusterColors";
import { goTo } from "../navigation";
import type { ClusterBoxes, ClusterMixing, ClusterPayload, ClusterSettings, Step } from "../types";
import {
  Empty,
  Loading,
  prerequisiteKey,
  SettingsReset,
  StepBar,
  sameAs,
  SectionTitle,
  stepOf,
  SubTitle,
  useRunButton,
  type PageProps,
} from "./stepShared";

export function ClusteringPage({ workflow, data, onNavigate }: PageProps & { data: ClusterPayload | undefined }) {
  const step = stepOf(workflow, "clustering");
  if (!data || !step) return <Loading />;
  return <ClusteringContent workflow={workflow} data={data} step={step} onNavigate={onNavigate} />;
}

function ClusteringContent({ workflow, data, step, onNavigate }: PageProps & { data: ClusterPayload; step: Step }) {
  const reference = data.settings.reference;
  const [boxes, setBoxes] = useState<ClusterBoxes>(() => clusterBoxes(data.settings.current));
  const set = (change: Partial<ClusterBoxes>) => setBoxes({ ...boxes, ...change });
  const run = useRunButton<ClusterBoxes>("run_cluster", "clustering", data, workflow.notice);
  const changedSinceRun = data.results !== null && !sameAs(boxes, clusterBoxes(data.settings.current));
  const r = data.results;


  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <StepHeader step={step} />
      {step.prerequisite && <Prerequisite title={step.prerequisite} text={step.prerequisite_text} onOpen={() => goTo(onNavigate, prerequisiteKey(workflow, step))} />}
      {step.reason && <Note dot="bg-accent">{step.reason}</Note>}

      <SettingsAndResults
        settings={
          <>
            <SectionTitle>Settings</SectionTitle>
            <RuledField
              id="cl_resolution"
              label="Leiden resolution"
              help="How finely to split the cells. Higher gives more, smaller clusters. Default: Scanpy's default."
              value={boxes.resolution}
              onChange={(v) => set({ resolution: v })}
              step={0.1}
              reference={reference.resolution}
            />
            <MoreSection title="More settings" note={Number(boxes.seed) === reference.seed ? undefined : "changed"}>
            <RuledField
              id="cl_seed"
              label="Leiden random seed"
              help="The same seed gives the same clusters every time. Default: Scanpy's default."
              value={boxes.seed}
              onChange={(v) => set({ seed: v })}
              step={1}
              reference={reference.seed}
            />
            <SubTitle>Fixed settings (validated method)</SubTitle>
            <InfoBox label="Fixed settings" className="my-3">
            <dl className="m-0 grid grid-cols-[max-content_minmax(0,1fr)] gap-x-6 gap-y-[6px] text-[14px]">
              <dt className="text-muted">Algorithm</dt>
              <dd className="m-0 font-mono">{data.fixed.flavor}</dd>
              <dt className="text-muted">Iterations</dt>
              <dd className="m-0">{data.fixed.n_iterations < 0 ? "until the clusters stop changing" : data.fixed.n_iterations}</dd>
              <dt className="text-muted">Graph</dt>
              <dd className="m-0">{data.fixed.directed ? "directed" : "undirected"}</dd>
            </dl>
            </InfoBox>
            {r?.notes.map((n) => (
              <p key={n} className="m-0 pb-2 text-[13px] leading-snug text-muted">
                {n}
              </p>
            ))}
            </MoreSection>
            <SettingsReset onReset={() => setBoxes(clusterBoxes(reference))} changedSinceRun={changedSinceRun} />
          </>
        }
        results={
          r ? (
            <>
              {r.warnings.map((w) => (
                <Note key={w} dot="bg-accent">
                  {w}
                </Note>
              ))}
              {r.mixing?.note && <Note dot="bg-accent">{r.mixing.note}</Note>}
              {r.mixing ? <ClusterTabs results={r} mixing={r.mixing} /> : <ClusterView results={r} />}
            </>
          ) : (
            <>
              <SectionTitle>Clusters</SectionTitle>
              <Empty>Run Clustering to see the clusters on the UMAP.</Empty>
            </>
          )
        }
      />
      <StepBar
        label="Clustering"
        nextTitle={stepOf(workflow, "marker_genes")?.title}
        status={step.status}
        blocked={step.prerequisite !== null}
        busy={run.busy}
        busyText="Running Clustering…"
        error={run.error}
        changedSinceRun={changedSinceRun}
        willClear={data.will_clear}
        onRun={() => run.click(boxes)}
        onContinue={() => goTo(onNavigate, "marker_genes")}
        prerequisite={step.prerequisite}
        prerequisiteText={step.prerequisite_text}
        confirmation={step.confirmation}
      />
    </div>
  );
}

type Results = NonNullable<ClusterPayload["results"]>;

/** The UMAP colored by cluster plus the cluster list, highlighting each other. */
function ClusterView({ results }: { results: Results }) {
  const groups = useMemo(
    () => results.clusters.map((c) => ({ key: c.cluster, label: `Cluster ${c.cluster}`, n: c.n_cells })),
    [results],
  );
  const colors = useMemo(() => clusterColors(groups.map((g) => g.key)), [groups]);
  return (
    <GroupedUmap
      title="Clusters on the UMAP"
      description="Hover or click a cluster to highlight it."
      listHeading="Cluster"
      groups={groups}
      colors={colors}
      cellGroup={results.labels}
      umap={results.umap}
      cells={results.cells}
    />
  );
}

/** With several samples: the clusters, and how well each cluster's samples mix. */
function ClusterTabs({ results, mixing }: { results: Results; mixing: ClusterMixing }) {
  const [view, setView] = useState<"clusters" | "mixing">("clusters");
  return (
    <Tabs
      label="Clustering views"
      tabs={[
        { key: "clusters", label: "Clusters" },
        { key: "mixing", label: "Sample mixing" },
      ]}
      active={view}
      onChange={setView}
    >
      {view === "clusters" ? <ClusterView results={results} /> : <MixingTable mixing={mixing} />}
    </Tabs>
  );
}

/** One row per cluster: cells, how many batches it holds, its mixing score, and the batch with the
 * most of its cells. Clusters below the threshold are teal. */
function MixingTable({ mixing }: { mixing: ClusterMixing }) {
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  const pct = (x: number) => `${Math.round(100 * x)}%`;
  return (
    <section aria-label="Sample mixing per cluster" className="flex flex-col gap-3">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <p className="m-0 text-[14px] text-body">
          How varied the {mixing.batch_key === "sample" ? "samples" : `“${mixing.batch_key}” batches`} are around each
          cluster's cells{mixing.corrected ? ", after batch correction" : ""}.{" "}
          <InfoTip label="About the mixing score">
            0 = a cluster's cells have neighbors from their own batch only; 1 = as mixed as the batch sizes allow.
            Clusters below {mixing.low_mixing} are marked: a batch effect, or a cell state found only in some samples.
          </InfoTip>
        </p>
        <CsvDownload
          name="Clustering sample mixing"
          columns={[
            { key: "group", header: "Cluster" },
            { key: "n_cells", header: "Cells" },
            { key: "n_samples", header: "Batches present" },
            { key: "mixing", header: "Mixing score" },
            { key: "top_sample", header: "Largest batch" },
            { key: "top_share", header: "Largest batch share" },
          ]}
          rows={() => mixing.rows.map((r) => ({ ...r }))}
        />
      </div>
      <div className="overflow-x-auto">
        <table className="w-full min-w-[560px] border-collapse text-[14px]">
          <thead>
            <tr className="border-b border-field text-left">
              <th scope="col" className={`${head} pr-3`}>Cluster</th>
              <th scope="col" className={`${head} px-3 text-right`}>Cells</th>
              <th scope="col" className={`${head} px-3 text-right`}>Batches</th>
              <th scope="col" className={`${head} px-3 text-right`}>Mixing</th>
              <th scope="col" className={`${head} pl-3`}>Largest batch</th>
            </tr>
          </thead>
          <tbody>
            {mixing.rows.map((r) => (
              <tr key={r.group} className={`border-b border-line ${r.low ? "bg-accent-tint" : ""}`}>
                <td className={`py-2 pr-3 font-bold ${r.low ? "pl-[10px] text-accent-text" : "text-ink"}`}>
                  {r.low ? "▲ " : ""}Cluster {r.group}
                </td>
                <td className="px-3 py-2 text-right font-mono">{r.n_cells.toLocaleString("en-US")}</td>
                <td className="px-3 py-2 text-right font-mono">{r.n_samples}</td>
                <td className="px-3 py-2 text-right font-mono">{r.mixing.toFixed(2)}</td>
                <td className="py-2 pl-3 text-body">
                  {r.top_sample} <span className="font-mono text-[13px] text-muted">({pct(r.top_share)})</span>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  );
}

function clusterBoxes(s: ClusterSettings): ClusterBoxes {
  return { resolution: decimal(s.resolution), seed: String(s.seed) }; // "1.0", not "1"
}

function decimal(v: number): string {
  return Number.isInteger(v) ? v.toFixed(1) : String(v);
}

