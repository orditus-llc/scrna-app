// The Preprocessing and Dimensionality Reduction pages: settings on the left, results on the
// right, and one button that runs the step and opens the next page.
import type { Data, Layout, Shape } from "plotly.js-dist-min";
import { useMemo, useState } from "react";

import { Note, RuledCheck, RuledField, RuledSelect } from "../components/fields";
import { InfoTip } from "../components/InfoTip";
import { MoreSection } from "../components/MoreSection";
import { Tabs } from "../components/Tabs";
import { PlotlyChart } from "../components/PlotlyChart";
import { Prerequisite, SettingsAndResults, StepHeader } from "../components/StepPage";
import { clusterColors } from "../lib/clusterColors";
import { goTo, type StepKey } from "../navigation";
import {
  type BatchInfo,
  type IntegrationResults,
  type PreprocessBoxes,
  type PreprocessPayload,
  type PreprocessSettings,
  type ReductionBoxes,
  type ReductionPayload,
  type ReductionSettings,
  type Step,
} from "../types";
import {
  Empty,
  Loading,
  prerequisiteKey,
  SettingsReset,
  StepBar,
  SectionTitle,
  stepOf,
  SubTitle,
  sameAs,
  useRunButton,
  type PageProps,
} from "./stepShared";

const INK = "#111111";
const MUTED = "#c4c4cc";
const TEAL = "#0d9488"; // accent (theme.css)
const GRID = "#ececef";

// --- Preprocessing ------------------------------------------------------------------------

export function PreprocessingPage({ workflow, data, onNavigate }: PageProps & { data: PreprocessPayload | undefined }) {
  const step = stepOf(workflow, "preprocessing");
  if (!data || !step) return <Loading />;
  return <PreprocessingContent workflow={workflow} data={data} step={step} onNavigate={onNavigate} />;
}

function PreprocessingContent({ workflow, data, step, onNavigate }: PageProps & { data: PreprocessPayload; step: Step }) {
  const reference = data.settings.reference;
  const [boxes, setBoxes] = useState<PreprocessBoxes>(() => preprocessBoxes(data.settings.current));
  const set = (change: Partial<PreprocessBoxes>) => setBoxes({ ...boxes, ...change });
  const run = useRunButton<PreprocessBoxes>("run_preprocess", "preprocessing", data, workflow.notice);
  const changedSinceRun = data.results !== null && !sameAs(boxes, preprocessBoxes(data.settings.current));
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
              id="pp_target_sum"
              label="Counts per cell after normalization"
              help="Each cell is scaled to this total, then log-transformed, so cells are comparable. Default: the median total counts of your cells (Scanpy's default)."
              value={boxes.target_sum}
              onChange={(v) => set({ target_sum: v })}
              step={1000}
              reference={reference.target_sum}
              off={boxes.target_median}
              onOff={(target_median) => set({ target_median })}
              offLabel="Median of your cells"
              offText="median of your cells"
            />
            <SubTitle>Highly variable genes · {data.fixed.hvg_flavor === "seurat" ? "Seurat" : data.fixed.hvg_flavor} method</SubTitle>
            <RuledSelect
              id="pp_hvg_method"
              label="How to pick them"
              help="Cut-offs keep every gene above a variability threshold, so the number depends on your data (Scanpy's default). Top N keeps a fixed number of the most variable genes."
              value={boxes.hvg_method}
              options={[
                { value: "cutoffs", label: "Cut-offs (number depends on the data)" },
                { value: "top", label: "Top N most variable genes" },
              ]}
              onChange={(v) => set({ hvg_method: v as PreprocessBoxes["hvg_method"] })}
              reference={reference.hvg_method}
            />
            {data.n_samples > 1 && (
              <RuledSelect
                id="pp_hvg_per_sample"
                label="Pick them"
                help="Within each sample: genes are judged sample by sample and ranked by how many samples find them variable, so genes that only differ between samples aren't picked (recommended by the single-cell best-practices guide). Across all cells: all samples judged together."
                value={boxes.per_sample ? "per_sample" : "pooled"}
                options={[
                  { value: "per_sample", label: `Within each sample (${data.n_samples} samples)` },
                  { value: "pooled", label: "Across all cells together" },
                ]}
                onChange={(v) => set({ per_sample: v === "per_sample" })}
                reference={reference.per_sample ? "per_sample" : "pooled"}
              />
            )}
            {boxes.hvg_method === "top" ? (
              <RuledField
                id="pp_n_top_genes"
                label="Number of genes"
                help="How many of the most variable genes to keep. Default: 2,000 (Seurat's default); results usually change little between about 1,000 and 5,000."
                value={boxes.n_top_genes}
                onChange={(v) => set({ n_top_genes: v })}
                step={500}
                min={1}
                reference={reference.n_top_genes}
              />
            ) : (
              <>
              <RuledField
                id="pp_min_mean"
                label="Min mean expression"
                help="Skips genes expressed too rarely to judge. Default: Scanpy's default."
                value={boxes.min_mean}
                onChange={(v) => set({ min_mean: v })}
                step={0.0025}
                reference={reference.min_mean}
              />
              <RuledField
                id="pp_max_mean"
                label="Max mean expression"
                help="Skips very highly expressed genes, often housekeeping genes. Default: Scanpy's default."
                value={boxes.max_mean}
                onChange={(v) => set({ max_mean: v })}
                step={0.5}
                reference={reference.max_mean}
              />
              <RuledField
                id="pp_min_disp"
                label="Min normalized dispersion"
                help="How much more variable than similar genes a gene must be. Default: Scanpy's default."
                value={boxes.min_disp}
                onChange={(v) => set({ min_disp: v })}
                step={0.1}
                min={-10}
                reference={reference.min_disp}
              />
              </>
            )}
            <RuledCheck
              id="pp_regress"
              label={`Regress out ${data.fixed.regress_keys.join(" and ")}`}
              help="Removes variation explained by sequencing depth and mitochondrial share. Off by default: the single-cell best-practices guide advises it only when needed, as it can remove real signal."
              checked={boxes.regress}
              onChange={(regress) => set({ regress })}
              reference={reference.regress}
            />
            <MoreSection title="More settings" note={Number(boxes.max_value) === reference.max_value ? undefined : "changed"}>
              <RuledField
                id="pp_max_value"
                label="Clip scaled values at"
                help="Caps scaled values at ± this, so a few extreme cells can't dominate. Default: Seurat's default."
                value={boxes.max_value}
                onChange={(v) => set({ max_value: v })}
                step={1}
                reference={reference.max_value}
              />
            </MoreSection>
            <SettingsReset onReset={() => setBoxes(preprocessBoxes(reference))} changedSinceRun={changedSinceRun} />
          </>
        }
        results={
          <>
            <SectionTitle>Highly variable genes</SectionTitle>
            {r ? <HvgPlot results={r} nSamples={data.n_samples} /> : <Empty>Run Preprocessing to see which genes are selected.</Empty>}
          </>
        }
      />
      <StepBar
        label="Preprocessing"
        nextTitle={stepOf(workflow, "dimensionality_reduction")?.title}
        status={step.status}
        blocked={step.prerequisite !== null}
        busy={run.busy}
        busyText="Running Preprocessing…"
        error={run.error}
        changedSinceRun={changedSinceRun}
        willClear={data.will_clear}
        onRun={() => run.click(boxes)}
        onContinue={() => goTo(onNavigate, "dimensionality_reduction")}
        prerequisite={step.prerequisite}
        prerequisiteText={step.prerequisite_text}
        confirmation={step.confirmation}
      />
    </div>
  );
}

function HvgPlot({ results, nSamples }: { results: NonNullable<PreprocessPayload["results"]>; nSamples: number }) {
  const plot = useMemo(() => {
    const { genes, settings } = results;
    const selected: number[] = [];
    const other: number[] = [];
    genes.name.forEach((_, i) => {
      if (genes.mean[i] === null || genes.dispersion_norm[i] === null) return; // no dispersion (e.g. constant gene)
      (genes.hvg[i] ? selected : other).push(i);
    });
    const trace = (idx: number[], name: string, color: string, size: number): Data => ({
      type: "scattergl",
      mode: "markers",
      name,
      x: idx.map((i) => genes.mean[i]),
      y: idx.map((i) => genes.dispersion_norm[i]),
      text: idx.map((i) => genes.name[i]),
      customdata: genes.n_samples_variable ? idx.map((i) => genes.n_samples_variable![i]) : undefined,
      marker: { size, color, opacity: 0.7 },
      hovertemplate:
        `<b>%{text}</b><br>Mean %{x:.3f}<br>Normalized dispersion %{y:.2f}` +
        (genes.n_samples_variable ? `<br>Variable in %{customdata} of ${nSamples} samples` : "") +
        `<extra>${name}</extra>`,
    });
    const vline = (x: number): Partial<Shape> => ({
      type: "line", xref: "x", yref: "paper", x0: x, x1: x, y0: 0, y1: 1, line: { color: TEAL, width: 1.5, dash: "dash" },
    });
    return {
      data: [trace(other, "not selected", MUTED, 3), trace(selected, "highly variable", INK, 3.5)],
      layout: {
        xaxis: { title: { text: "Mean expression (log)" }, gridcolor: GRID },
        yaxis: { title: { text: "Normalized dispersion" }, gridcolor: GRID },
        // Cut-off lines only when the cut-offs picked the genes (top N has no threshold).
        shapes:
          settings.hvg_method === "cutoffs"
            ? [
                vline(settings.min_mean),
                vline(settings.max_mean),
                { type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: settings.min_disp, y1: settings.min_disp, line: { color: TEAL, width: 1.5, dash: "dash" } },
              ]
            : [],
      } as Partial<Layout>,
    };
  }, [results, nSamples]);
  const s = results.settings;
  return (
    <>
      {results.per_sample_used && (
        <p className="m-0 text-[13px] text-body">
          Picked within each of the {nSamples} samples; the dots show each gene's values averaged over the samples.
        </p>
      )}
      {s.hvg_method === "cutoffs" ? (
        <p className="m-0 text-[13px] text-body">
          Each dot is a gene: dark = highly variable, grey = not. Teal lines are the cut-offs.{" "}
          <InfoTip label="About the cut-offs">
            Selected genes have a mean between {s.min_mean} and {s.max_mean} and a normalized dispersion above {s.min_disp}.
            Only these genes are used for PCA; all genes are kept for marker genes.
          </InfoTip>
        </p>
      ) : (
        <p className="m-0 text-[13px] text-body">
          Each dot is a gene: dark = the {s.n_top_genes.toLocaleString("en-US")} most variable (highest normalized
          dispersion), grey = not.{" "}
          <InfoTip label="About the top genes">
            Genes are ranked by normalized dispersion within bins of similar mean, and the top ones are kept. Only these
            genes are used for PCA; all genes are kept for marker genes.
          </InfoTip>
        </p>
      )}
      <PlotlyChart
        data={plot.data}
        layout={plot.layout}
        height={420}
        label="Mean expression against normalized dispersion for every gene"
        download={{
          name: "Preprocessing highly variable genes",
          title: "Highly variable genes: mean expression against normalized dispersion (dark = kept)",
        }}
      />
    </>
  );
}

function preprocessBoxes(s: PreprocessSettings): PreprocessBoxes {
  return {
    target_sum: s.target_sum === null ? "" : String(s.target_sum),
    target_median: s.target_sum === null,
    hvg_method: s.hvg_method,
    n_top_genes: String(s.n_top_genes),
    min_mean: String(s.min_mean),
    max_mean: String(s.max_mean),
    min_disp: String(s.min_disp),
    per_sample: s.per_sample,
    max_value: String(s.max_value),
    regress: s.regress,
  };
}

// --- Dimensionality Reduction ----------------------------------------------------------------

export function ReductionPage({ workflow, data, onNavigate }: PageProps & { data: ReductionPayload | undefined }) {
  const step = stepOf(workflow, "dimensionality_reduction");
  if (!data || !step) return <Loading />;
  return <ReductionContent workflow={workflow} data={data} step={step} onNavigate={onNavigate} />;
}

function ReductionContent({ workflow, data, step, onNavigate }: PageProps & { data: ReductionPayload; step: Step }) {
  const reference = data.settings.reference;
  const [boxes, setBoxes] = useState<ReductionBoxes>(() => reductionBoxes(data.settings.current));
  const set = (change: Partial<ReductionBoxes>) => setBoxes({ ...boxes, ...change });
  const run = useRunButton<ReductionBoxes>("run_reduction", "dimensionality_reduction", data, workflow.notice);
  const changedSinceRun = data.results !== null && !sameAs(boxes, reductionBoxes(data.settings.current));
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
              id="dr_n_comps"
              label="Principal components"
              help="How many principal components to compute. The elbow plot shows what each one adds. Default: Scanpy's default."
              value={boxes.n_comps}
              onChange={(v) => set({ n_comps: v })}
              step={1}
              min={1}
              reference={reference.n_comps}
            />
            <RuledField
              id="dr_n_pcs"
              label="PCs used for the neighbor graph"
              help="How many PCs are used to find each cell's neighbors. Default: all computed PCs (Scanpy's default)."
              value={boxes.n_pcs}
              onChange={(v) => set({ n_pcs: v })}
              step={1}
              min={1}
              reference={reference.n_pcs}
            />
            <RuledField
              id="dr_n_neighbors"
              label="Neighbors per cell"
              help="How many similar cells each cell is linked to. Fewer keep small groups separate. Default: Scanpy's default."
              value={boxes.n_neighbors}
              onChange={(v) => set({ n_neighbors: v })}
              step={1}
              min={2}
              reference={reference.n_neighbors}
            />
            {data.batch && (
              <BatchSettings
                info={data.batch}
                boxes={boxes}
                reference={reference}
                suggest={r?.integration?.suggest ?? false}
                onChange={set}
              />
            )}
            <MoreSection title="More settings" note={Number(boxes.seed) === reference.seed ? undefined : "changed"}>
              <RuledField
                id="dr_seed"
                label="Random seed"
                help="The same seed gives the same PCA, graph and UMAP every time. Default: Scanpy's default."
                value={boxes.seed}
                onChange={(v) => set({ seed: v })}
                step={1}
                reference={reference.seed}
              />
            </MoreSection>
            <SettingsReset onReset={() => setBoxes(reductionBoxes(reference))} changedSinceRun={changedSinceRun} />
          </>
        }
        results={r ? <ReductionPlots results={r} /> : (
          <>
            <SectionTitle>Results</SectionTitle>
            <Empty>Run Dimensionality Reduction to see the elbow plot and the UMAP.</Empty>
          </>
        )}
      />
      <StepBar
        label="Dimensionality Reduction"
        nextTitle={stepOf(workflow, "clustering")?.title}
        status={step.status}
        blocked={step.prerequisite !== null}
        busy={run.busy}
        busyText="Running Dimensionality Reduction… (about 10 seconds)"
        error={run.error}
        changedSinceRun={changedSinceRun}
        willClear={data.will_clear}
        onRun={() => run.click(boxes)}
        onContinue={() => goTo(onNavigate, "clustering")}
        prerequisite={step.prerequisite}
        prerequisiteText={step.prerequisite_text}
        confirmation={step.confirmation}
      />
    </div>
  );
}

type Metric = "n_genes" | "total_counts" | "pct_mt" | "batch";
const METRICS: Record<Metric, string> = {
  n_genes: "Genes per cell",
  total_counts: "Counts per cell",
  pct_mt: "% mitochondrial",
  batch: "Batch",
};

function ReductionPlots({ results }: { results: NonNullable<ReductionPayload["results"]> }) {
  const integ = results.integration;
  const [metric, setMetric] = useState<Metric>(integ ? "batch" : "total_counts");
  const [view, setView] = useState<"umap" | "elbow">("umap");
  const cells = results.cells;
  const available = (Object.keys(METRICS) as Metric[]).filter(
    (m) => (m !== "pct_mt" || cells.pct_mt !== null) && (m !== "batch" || integ !== null),
  );

  const elbow = useMemo(() => {
    const pcs = results.variance_ratio.map((_, i) => i + 1);
    const used = results.n_pcs_used;
    const knee = results.elbow_pc;
    return {
      data: [
        {
          type: "scatter",
          mode: "lines+markers",
          x: pcs,
          y: results.variance_ratio.map((v) => 100 * v),
          marker: { size: 6, color: pcs.map((pc) => (pc <= used ? INK : MUTED)) },
          line: { color: "#a1a1aa", width: 1 },
          hovertemplate: "PC %{x}: %{y:.2f}% of variance<extra></extra>",
        },
        {
          // The elbow: a ring around its point (guidance only).
          type: "scatter",
          mode: "markers",
          x: [knee],
          y: [100 * results.variance_ratio[knee - 1]],
          marker: { size: 16, color: "rgba(0,0,0,0)", line: { color: TEAL, width: 2 } },
          hovertemplate: `The curve bends here (PC ${knee})<extra></extra>`,
        },
      ] as Data[],
      layout: {
        xaxis: { title: { text: "Principal component" }, gridcolor: GRID },
        yaxis: { title: { text: "% of variance explained" }, gridcolor: GRID, rangemode: "tozero" },
        showlegend: false,
        annotations: [
          {
            x: knee, y: 100 * results.variance_ratio[knee - 1], xref: "x", yref: "y", ax: 40, ay: -30,
            text: `Bends at PC ${knee}`, showarrow: true, arrowcolor: TEAL, font: { size: 12, color: "#0F766E" },
            bgcolor: "rgba(255,255,255,0.85)",
          },
        ],
        shapes:
          used < pcs.length
            ? [{ type: "line", xref: "x", yref: "paper", x0: used + 0.5, x1: used + 0.5, y0: 0, y1: 1, line: { color: TEAL, width: 1.5, dash: "dash" } }]
            : [],
      } as Partial<Layout>,
    };
  }, [results]);

  const umap = useMemo(() => {
    const pct = cells.pct_mt;
    if (metric === "batch" && integ) {
      // One trace per batch, so the legend names them.
      const colors = clusterColors(integ.batches);
      return {
        data: integ.batches.map((name, b) => {
          const idx = integ.cell_batch.flatMap((x, i) => (x === b ? [i] : []));
          return {
            type: "scattergl",
            mode: "markers",
            name,
            x: idx.map((i) => results.umap.x[i]),
            y: idx.map((i) => results.umap.y[i]),
            text: idx.map((i) => cells.barcode[i]),
            marker: { size: 4, color: colors[name], opacity: 0.8 },
            hovertemplate: `<b>%{text}</b><br>${name}<extra></extra>`,
          };
        }) as Data[],
        layout: {
          xaxis: { title: { text: "UMAP 1" }, showticklabels: false, showgrid: false, zeroline: false },
          yaxis: { title: { text: "UMAP 2" }, showticklabels: false, showgrid: false, zeroline: false, scaleanchor: "x" },
          margin: { l: 40, r: 16, t: 16, b: 40 },
          showlegend: true,
          legend: { itemsizing: "constant", font: { size: 12 }, title: { text: integ.batch_key } },
        } as Partial<Layout>,
      };
    }
    const values = metric === "pct_mt" ? cells.pct_mt! : metric === "n_genes" ? cells.n_genes : cells.total_counts;
    return {
      data: [
        {
          type: "scattergl",
          mode: "markers",
          x: results.umap.x,
          y: results.umap.y,
          text: cells.barcode,
          customdata: cells.barcode.map((_, i) => [cells.n_genes[i], cells.total_counts[i], pct ? pct[i] : null]),
          marker: {
            size: 4,
            color: values,
            colorscale: "Viridis",
            showscale: true,
            colorbar: { title: { text: METRICS[metric], side: "right" }, thickness: 12, outlinewidth: 0 },
          },
          hovertemplate:
            "<b>%{text}</b><br>Genes %{customdata[0]:,}<br>Counts %{customdata[1]:,.0f}" +
            (pct ? "<br>% mito %{customdata[2]:.2f}" : "") +
            "<extra></extra>",
        },
      ] as Data[],
      layout: {
        xaxis: { title: { text: "UMAP 1" }, showticklabels: false, showgrid: false, zeroline: false },
        yaxis: { title: { text: "UMAP 2" }, showticklabels: false, showgrid: false, zeroline: false, scaleanchor: "x" },
        margin: { l: 40, r: 16, t: 16, b: 40 },
      } as Partial<Layout>,
    };
  }, [results, metric, cells, integ]);

  return (
    <Tabs
      label="Dimensionality Reduction views"
      tabs={[
        { key: "umap", label: "UMAP" },
        { key: "elbow", label: "Elbow plot" },
      ]}
      active={view}
      onChange={setView}
      right={
        view === "umap" ? (
          <label className="flex items-center gap-2 text-[13px]">
            Color cells by
            <select
              value={metric}
              onChange={(e) => setMetric(e.target.value as Metric)}
              className="rounded-[4px] border border-field bg-white px-2 py-1 text-[13px]"
            >
              {available.map((m) => (
                <option key={m} value={m}>
                  {METRICS[m]}
                </option>
              ))}
            </select>
          </label>
        ) : undefined
      }
    >
      {view === "umap" ? (
        <>
          {integ && <MixingSummary integ={integ} />}
          <p className="m-0 text-[13px] text-body">
            {metric === "batch" ? "Each dot is a cell, colored by its batch: well-corrected batches overlap." : "Each dot is a cell; hover for its QC values."}
          </p>
          <PlotlyChart
            data={umap.data}
            layout={umap.layout}
            height={460}
            label={`UMAP of all cells, colored by ${METRICS[metric]}`}
            download={{ name: `Dimensionality Reduction UMAP ${METRICS[metric]}`, title: `UMAP of all cells, colored by ${METRICS[metric]}` }}
          />
        </>
      ) : (
        <>
          <p className="m-0 text-[13px] text-body">
            Variance explained by each PC; dark ones are used for the graph. The curve bends around PC {results.elbow_pc}.{" "}
            <InfoTip label="About the bend">
              Found from your data: the point farthest from a straight line between the first and last PC. Treat it as a
              lower bound: later PCs can still separate small cell populations, so using more is usually safe. All
              computed PCs are used by default (Scanpy's default).
            </InfoTip>
          </p>
          <PlotlyChart
            data={elbow.data}
            layout={elbow.layout}
            height={380}
            label="Variance explained by each principal component"
            download={{
              name: "Dimensionality Reduction elbow plot",
              title: `Variance explained by each principal component (curve bends at PC ${results.elbow_pc})`,
            }}
          />
        </>
      )}
    </Tabs>
  );
}

/** The sample mixing score, before and (with Harmony) after correction, plus any warning. */
function MixingSummary({ integ }: { integ: IntegrationResults }) {
  const fmt2 = (x: number) => x.toFixed(2);
  return (
    <div className="flex flex-col gap-2">
      <p className="m-0 text-[14px] text-body">
        Sample mixing{" "}
        <InfoTip label="About sample mixing">
          For each cell: how varied the batches are among its 30 nearest neighbors, from 0 (all from its own batch) to 1
          (as mixed as the batch sizes allow), averaged over all cells. Low mixing can be a batch effect, or real
          differences between conditions.
        </InfoTip>
        :{" "}
        {integ.method === "harmony" ? (
          <>
            <b className="font-mono text-ink">{fmt2(integ.mixing_before)}</b> before →{" "}
            <b className="font-mono text-ink">{fmt2(integ.mixing_after!)}</b> after Harmony ({integ.n_batches} batches from
            “{integ.batch_key}”)
          </>
        ) : (
          <>
            <b className="font-mono text-ink">{fmt2(integ.mixing_before)}</b> across {integ.n_batches} batches, no
            correction
          </>
        )}
      </p>
      {integ.warnings.map((w) => (
        <Note key={w} dot="bg-accent">{w}</Note>
      ))}
      {integ.suggest && (
        <Note dot="bg-accent">The batches mix poorly. Try “Correct batch effects (Harmony)” in the settings and compare.</Note>
      )}
    </div>
  );
}

function reductionBoxes(s: ReductionSettings): ReductionBoxes {
  return {
    n_comps: String(s.n_comps),
    n_neighbors: String(s.n_neighbors),
    n_pcs: String(s.n_pcs),
    seed: String(s.seed),
    batch_correction: s.batch_correction,
    batch_key: s.batch_key,
  };
}

/** Batch correction (several samples only): Harmony on/off and which column holds the batches. */
function BatchSettings({
  info,
  boxes,
  reference,
  suggest,
  onChange,
}: {
  info: BatchInfo;
  boxes: ReductionBoxes;
  reference: ReductionSettings;
  suggest: boolean;
  onChange: (change: Partial<ReductionBoxes>) => void;
}) {
  const chosen = info.choices.find((c) => c.value === boxes.batch_key);
  return (
    <>
      <RuledCheck
        id="dr_batch"
        label="Correct batch effects (Harmony)"
        help="Lines up cells of the same type across batches, so clusters follow cell types rather than samples. Only the map used for the neighbor graph, clustering and UMAP changes; gene counts and marker genes stay as measured. Off by default: turn it on when samples don't mix."
        checked={boxes.batch_correction}
        onChange={(on) => onChange({ batch_correction: on })}
        reference={reference.batch_correction}
      />
      {suggest && !boxes.batch_correction && (
        <p className="m-0 border-b border-line py-[10px] text-[13px] font-bold text-accent-text">
          Samples mix poorly on the last run: try batch correction.
        </p>
      )}
      {boxes.batch_correction && (
        <RuledSelect
          id="dr_batch_key"
          label="Batch"
          help="The column whose values are the batches: usually each sample, or the sequencing run if several samples were run together."
          value={boxes.batch_key}
          options={info.choices.map((c) => ({ value: c.value, label: `${c.label} (${c.n_values})` }))}
          onChange={(v) => onChange({ batch_key: v })}
          reference={reference.batch_key}
          note={chosen?.warning ?? undefined}
        />
      )}
    </>
  );
}

// --- A step that isn't built yet ----------------------------------------------------------

export function ComingSoonPage({ workflow, stepKey, onNavigate }: PageProps & { stepKey: StepKey }) {
  const step = stepOf(workflow, stepKey);
  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <StepHeader step={step} />
      {step?.prerequisite && <Prerequisite title={step.prerequisite} text={step.prerequisite_text} onOpen={() => goTo(onNavigate, prerequisiteKey(workflow, step))} />}
      <Note dot="bg-off">This page is coming soon: the analysis behind it is built and validated, and its settings and results will appear here next.</Note>
    </div>
  );
}
