// Step 5: Leiden clustering. Settings on the left; the UMAP colored by cluster and the
// cluster list on the right. Hovering a cell or a cluster highlights that cluster and fades
// the rest; clicking keeps the highlight (Esc or "Clear" removes it).
import { useMemo, useState } from "react";

import { InfoBox, Note, RuledField } from "../components/fields";
import { GroupedUmap } from "../components/GroupedUmap";
import { MoreSection } from "../components/MoreSection";
import { Prerequisite, SettingsAndResults, StepHeader } from "../components/StepPage";
import { clusterColors } from "../lib/clusterColors";
import { goTo } from "../navigation";
import type { ClusterBoxes, ClusterPayload, ClusterSettings, Step } from "../types";
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
              <ClusterView results={r} />
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

function clusterBoxes(s: ClusterSettings): ClusterBoxes {
  return { resolution: decimal(s.resolution), seed: String(s.seed) }; // "1.0", not "1"
}

function decimal(v: number): string {
  return Number.isInteger(v) ? v.toFixed(1) : String(v);
}

