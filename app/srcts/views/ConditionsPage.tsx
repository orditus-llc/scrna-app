// Compare Conditions (several samples with conditions): what changed between two conditions, with
// every sample counted once. Settings on the left (the two conditions, how cells are grouped, how
// samples are paired, which groups' genes to test); on the right the design, then two tabs: each
// group's share of the cells, and each group's genes (pseudobulk, volcano plot + table).
import { useShinyInput } from "@posit-dev/shinyreact";
import type { Data, Layout } from "plotly.js-dist-min";
import { useMemo, useState } from "react";

import { CsvDownload } from "../components/DownloadMenu";
import { InfoBox, Note, RuledField, RuledSelect } from "../components/fields";
import { InfoTip } from "../components/InfoTip";
import { MoreSection } from "../components/MoreSection";
import { PlotlyChart } from "../components/PlotlyChart";
import { Prerequisite, SettingsAndResults } from "../components/StepPage";
import { Tabs } from "../components/Tabs";
import { cn } from "../lib/utils";
import { goTo, type Page } from "../navigation";
import type { ConditionSettings, ConditionsPayload, GeneChange, WorkflowPayload } from "../types";
import { Empty, fmt, Loading, SectionTitle, SettingsReset, StepBar, SubTitle, sameAs, useRunButton } from "./stepShared";

const INK = "#111111";
const MUTED = "#c4c4cc";
const TEAL = "#0d9488";
const GRID = "#ececef";
const TOP_ROWS = 50; // genes listed in the table (all of them in the CSV)

/** What the page sends with Run (strings while typing). */
interface ConditionBoxes {
  reference: string;
  test: string;
  group_by: string;
  paired_by: string; // "" = don't pair
  groups: string[] | "all";
  min_cells: string;
  min_gene_counts: string;
  max_adj_p: string;
}

function toBoxes(s: ConditionSettings, groups: string[] | null): ConditionBoxes {
  return {
    reference: s.reference ?? "",
    test: s.test ?? "",
    group_by: s.group_by,
    paired_by: s.paired_by ?? "",
    groups: groups ?? "all",
    min_cells: String(s.min_cells),
    min_gene_counts: String(s.min_gene_counts),
    max_adj_p: String(s.max_adj_p),
  };
}

export function ConditionsPage({
  workflow,
  data,
  onNavigate,
}: {
  workflow: WorkflowPayload;
  data: ConditionsPayload | undefined;
  onNavigate: (page: Page) => void;
}) {
  if (!data) return <Loading />;
  if (!data.visible || !data.options || !data.settings)
    return (
      <div className="px-8 pt-[30px] pb-10">
        <Empty>Compare Conditions needs several samples, each with a condition.</Empty>
      </div>
    );
  return <ConditionsContent workflow={workflow} data={data} onNavigate={onNavigate} />;
}

function ConditionsContent({
  workflow,
  data,
  onNavigate,
}: {
  workflow: WorkflowPayload;
  data: ConditionsPayload;
  onNavigate: (page: Page) => void;
}) {
  const options = data.options!;
  const settings = data.settings!;
  const [boxes, setBoxes] = useState<ConditionBoxes>(() => toBoxes(settings.current, settings.groups));
  const set = (change: Partial<ConditionBoxes>) => setBoxes({ ...boxes, ...change });
  const run = useRunButton<ConditionBoxes>("run_conditions", "compare_conditions", data, workflow.notice);
  const r = data.results;
  // The group list is compared as text (sameAs compares single values only).
  const comparable = (b: ConditionBoxes) => ({ ...b, groups: b.groups === "all" ? "all" : b.groups.join("\n") });
  const changedSinceRun = r !== null && !sameAs(comparable(boxes), comparable(toBoxes(settings.current, settings.groups)));
  const grouping = options.groupings.find((g) => g.value === boxes.group_by) ?? options.groupings[0];
  const nChosen = boxes.groups === "all" ? grouping?.groups.length ?? 0 : boxes.groups.length;
  const conditionLabel = (name: string) => {
    const c = options.conditions.find((x) => x.name === name);
    return c ? `${name} (${c.n_samples} ${c.n_samples === 1 ? "sample" : "samples"})` : name;
  };
  const missing =
    !boxes.reference || !boxes.test
      ? "Choose the two conditions to compare."
      : boxes.reference === boxes.test
        ? "Choose two different conditions."
        : null;
  const prerequisite = data.enabled ? null : data.reason;

  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <h1 className="sr-only">Compare Conditions</h1>
      {prerequisite && <Prerequisite title="Clustering" text={prerequisite} onOpen={() => goTo(onNavigate, "clustering")} />}
      {data.status === "needs_rerun" && (
        <Note dot="bg-accent">The cells or their groups changed since the last comparison, so run it again.</Note>
      )}
      <SettingsAndResults
        settings={
          <>
            <SectionTitle>Settings</SectionTitle>
            <RuledSelect
              id="cc_test"
              label="Compare"
              help="The condition you want to know about, for example the treated or stimulated samples."
              value={boxes.test}
              options={options.conditions.map((c) => ({ value: c.name, label: conditionLabel(c.name) }))}
              onChange={(test) => set({ test })}
              reference={settings.default.test ?? undefined}
            />
            <RuledSelect
              id="cc_reference"
              label="Against (baseline)"
              help="The condition to compare with, usually the controls. A condition named like a control is suggested."
              value={boxes.reference}
              options={options.conditions.map((c) => ({ value: c.name, label: conditionLabel(c.name) }))}
              onChange={(reference) => set({ reference })}
              reference={settings.default.reference ?? undefined}
            />
            <RuledSelect
              id="cc_group_by"
              label="Group cells by"
              help="Which groups to compare: clusters, the suggested cell types, or your own labels from Cell Type Annotation."
              value={boxes.group_by}
              options={options.groupings.map((g) => ({ value: g.value, label: g.label }))}
              onChange={(group_by) => set({ group_by, groups: "all" })}
              reference={settings.default.group_by}
            />
            <RuledSelect
              id="cc_paired_by"
              label="Pair samples by"
              help="When the same person (or animal) gave a sample in each condition, pairing compares each person with themselves, which removes person-to-person differences. Suggested from a column with one value per sample that appears under both conditions."
              value={boxes.paired_by}
              options={[{ value: "", label: "Don't pair" }, ...options.pairing.map((c) => ({ value: c, label: c }))]}
              onChange={(paired_by) => set({ paired_by })}
              reference={settings.default.paired_by ?? ""}
              note={options.pairing.length === 0 ? "No column pairs the samples, so they are compared unpaired." : undefined}
            />
            <SubTitle>Groups to test for gene changes</SubTitle>
            <GroupPicker
              groups={grouping?.groups ?? []}
              chosen={boxes.groups}
              onChange={(groups) => set({ groups })}
            />
            <p className="m-0 text-[13px] text-muted">
              {nChosen} {nChosen === 1 ? "group" : "groups"} · about {Math.max(1, Math.round((nChosen * 8) / 60))} min. Shares are
              compared for every group.
            </p>
            <MoreSection
              title="More settings"
              note={
                sameAs(
                  { a: boxes.min_cells, b: boxes.min_gene_counts, c: boxes.max_adj_p },
                  { a: String(settings.default.min_cells), b: String(settings.default.min_gene_counts), c: String(settings.default.max_adj_p) },
                )
                  ? undefined
                  : "changed"
              }
            >
              <RuledField
                id="cc_min_cells"
                label="Min cells per sample"
                help="A sample counts for a group only with at least this many of its cells; fewer give unreliable sums. The app's own rule."
                value={boxes.min_cells}
                onChange={(min_cells) => set({ min_cells })}
                step={5}
                min={1}
                reference={settings.default.min_cells}
              />
              <RuledField
                id="cc_min_gene_counts"
                label="Min counts per gene"
                help="Genes with fewer counts in total over a group's samples are left out (DESeq2's usual pre-filter)."
                value={boxes.min_gene_counts}
                onChange={(min_gene_counts) => set({ min_gene_counts })}
                step={5}
                min={0}
                reference={settings.default.min_gene_counts}
              />
              <RuledField
                id="cc_max_adj_p"
                label="Max adjusted p-value"
                help="What counts as a significant change in the summaries (DESeq2's default)."
                value={boxes.max_adj_p}
                onChange={(max_adj_p) => set({ max_adj_p })}
                step={0.01}
                reference={settings.default.max_adj_p}
              />
            </MoreSection>
            <SettingsReset onReset={() => setBoxes(toBoxes(settings.default, null))} changedSinceRun={changedSinceRun} />
          </>
        }
        results={r ? <ConditionResults results={r} /> : <Empty>Run Compare Conditions to see what changed between the conditions.</Empty>}
      />
      <StepBar
        label="Compare Conditions"
        nextTitle={undefined}
        status={data.status}
        blocked={!data.enabled}
        busy={run.busy}
        busyText={`Running Compare Conditions (about ${Math.max(1, Math.round((nChosen * 8) / 60))} min)…`}
        error={run.error}
        changedSinceRun={changedSinceRun}
        willClear={[]}
        missing={missing}
        onRun={() => run.click(boxes)}
        onContinue={() => undefined}
        prerequisite={prerequisite ? "Clustering" : null}
        prerequisiteText={prerequisite}
        confirmation={data.confirmation ?? null}
      />
    </div>
  );
}

/** Tick the groups whose genes to test (all by default). */
function GroupPicker({
  groups,
  chosen,
  onChange,
}: {
  groups: { group: string; n_cells: number }[];
  chosen: string[] | "all";
  onChange: (chosen: string[] | "all") => void;
}) {
  const isOn = (g: string) => chosen === "all" || chosen.includes(g);
  const toggle = (g: string) => {
    const current = chosen === "all" ? groups.map((x) => x.group) : chosen;
    const next = current.includes(g) ? current.filter((x) => x !== g) : [...current, g];
    onChange(next.length === groups.length ? "all" : groups.map((x) => x.group).filter((x) => next.includes(x)));
  };
  return (
    <div className="flex flex-col gap-2 py-2">
      <div className="flex gap-3 text-[13px]">
        <button type="button" className="font-bold underline" onClick={() => onChange("all")}>
          All
        </button>
        <button type="button" className="font-bold underline" onClick={() => onChange([])}>
          None
        </button>
      </div>
      <div className="flex max-h-[220px] flex-col gap-1 overflow-y-auto">
        {groups.map((g) => (
          <label key={g.group} className="flex items-center gap-2 text-[14px]">
            <input type="checkbox" checked={isOn(g.group)} onChange={() => toggle(g.group)} className="size-4 accent-ink" />
            <span className="min-w-0 flex-1 break-words">{g.group}</span>
            <span className="font-mono text-[13px] text-muted">{fmt(g.n_cells)}</span>
          </label>
        ))}
      </div>
    </div>
  );
}

type Results = NonNullable<ConditionsPayload["results"]>;

function ConditionResults({ results }: { results: Results }) {
  const [view, setView] = useState<"shares" | "genes">("shares");
  const d = results.design;
  const nRef = d.samples.filter((s) => s.side === "reference").length;
  const nTest = d.samples.filter((s) => s.side === "test").length;
  return (
    <>
      <InfoBox label="What was compared" className="flex flex-col gap-1">
        <p className="m-0 text-[14px] text-body">
          <b className="text-ink">{d.test}</b> ({nTest} samples) against <b className="text-ink">{d.reference}</b> ({nRef} samples)
          {d.paired_by ? `, paired by ${d.paired_by} (${nTest} pairs)` : ", unpaired"} · cells grouped by{" "}
          {results.grouping_label.toLowerCase()}.
        </p>
        <p className="m-0 text-[13px] text-muted">
          Each sample counts once, so the tests ask whether a change is larger than the differences between samples. With few
          samples, real changes can be missed.
          {d.excluded.length > 0 && ` Left out: ${d.excluded.join(", ")}.`}
        </p>
      </InfoBox>
      <Tabs
        label="Comparison views"
        tabs={[
          { key: "shares", label: "Cell shares" },
          { key: "genes", label: "Gene changes" },
        ]}
        active={view}
        onChange={setView}
      >
        {view === "shares" ? <SharesView results={results} /> : <GenesView results={results} />}
      </Tabs>
    </>
  );
}

function pText(p: number | null): string {
  if (p === null) return "–";
  if (p === 0) return "< 1 × 10⁻³⁰⁰";
  return p >= 0.001 ? p.toFixed(3) : p.toExponential(1).replace("e", " × 10^");
}

function pct(x: number): string {
  return `${(100 * x).toFixed(1)}%`;
}

function SharesView({ results }: { results: Results }) {
  const d = results.design;
  const alpha = results.settings.max_adj_p;
  const rows = results.composition;
  const sig = (adj: number | null) => adj !== null && adj < alpha;
  const plot = useMemo(() => {
    const y = rows.map((r) => r.group);
    const lines: Data[] = rows.map((r) => ({
      type: "scatter",
      mode: "lines",
      x: [100 * r.mean_share_reference, 100 * r.mean_share_test],
      y: [r.group, r.group],
      line: { color: sig(r.adj_p) ? TEAL : MUTED, width: 2 },
      hoverinfo: "skip",
      showlegend: false,
    }));
    const dots = (name: string, values: number[], color: string): Data => ({
      type: "scatter",
      mode: "markers",
      name,
      x: values.map((v) => 100 * v),
      y,
      marker: { size: 9, color },
      customdata: rows.map((r) => [pText(r.adj_p)]),
      hovertemplate: `<b>%{y}</b><br>${name}: %{x:.1f}% of cells<br>Adjusted p-value %{customdata[0]}<extra></extra>`,
    });
    return {
      data: [
        ...lines,
        dots(d.reference, rows.map((r) => r.mean_share_reference), "#a1a1aa"),
        dots(d.test, rows.map((r) => r.mean_share_test), INK),
      ],
      layout: {
        xaxis: { title: { text: "Mean share of the sample's cells (%)" }, gridcolor: GRID, rangemode: "tozero" },
        yaxis: { automargin: true, autorange: "reversed", type: "category" },
        legend: { orientation: "h", y: 1.08 },
        margin: { l: 10 },
      } as Partial<Layout>,
    };
  }, [rows, d.reference, d.test, alpha]); // eslint-disable-line react-hooks/exhaustive-deps
  const nSig = rows.filter((r) => sig(r.adj_p)).length;
  return (
    <div className="flex flex-col gap-3">
      <p className="m-0 text-[13px] text-body">
        Each group's mean share of the cells in {d.reference} (grey) and {d.test} (dark); teal lines changed significantly
        (adjusted p-value below {alpha}). {nSig} of {rows.length} groups changed.
        {results.settings.group_by === "leiden" &&
          " When a condition changes expression strongly, clusters can split by condition, so their shares change without any change in cell numbers; compare cell types to see changes in numbers."}{" "}
        <InfoTip label="About the share test">
          Each sample's share of the group is logit-transformed and compared between the conditions with a t-test (paired when
          samples are paired), as in the propeller method. p-values are adjusted for testing many groups.
        </InfoTip>
      </p>
      <PlotlyChart
        data={plot.data}
        layout={plot.layout}
        height={Math.max(260, 28 * rows.length + 90)}
        label="Each group's mean share of cells in the two conditions"
        download={{ name: "Compare Conditions cell shares", title: `Share of cells per group: ${d.test} vs ${d.reference}` }}
      />
      <div className="flex justify-end">
        <CsvDownload
          name="Compare Conditions cell shares"
          columns={[
            { key: "group", header: "Group" },
            { key: "mean_share_reference", header: `Mean share in ${d.reference}` },
            { key: "mean_share_test", header: `Mean share in ${d.test}` },
            { key: "log2_ratio", header: "log2 ratio of shares" },
            { key: "t", header: "t statistic" },
            { key: "p", header: "p-value" },
            { key: "adj_p", header: "Adjusted p-value" },
            { key: "cells_reference", header: `Cells in ${d.reference}` },
            { key: "cells_test", header: `Cells in ${d.test}` },
          ]}
          rows={() => rows.map((r) => ({ ...r }))}
        />
      </div>
      <table className="w-full border-collapse text-[14px]">
        <thead>
          <tr className="border-b border-ink text-left text-[12px] font-bold tracking-[0.04em] text-muted uppercase">
            <th scope="col" className="py-2 pr-3">Group</th>
            <th scope="col" className="px-3 py-2 text-right">{d.reference}</th>
            <th scope="col" className="px-3 py-2 text-right">{d.test}</th>
            <th scope="col" className="px-3 py-2 text-right">Change</th>
            <th scope="col" className="py-2 pl-3 text-right">Adjusted p-value</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((r) => (
            <tr key={r.group} className={cn("border-b border-line", sig(r.adj_p) && "bg-accent-tint")}>
              <td className="py-[9px] pr-3 font-bold">{r.group}</td>
              <td className="px-3 py-[9px] text-right font-mono">{pct(r.mean_share_reference)}</td>
              <td className="px-3 py-[9px] text-right font-mono">{pct(r.mean_share_test)}</td>
              <td className="px-3 py-[9px] text-right font-mono">
                {r.log2_ratio > 0 ? "▲" : "▼"} {Math.abs(2 ** r.log2_ratio).toFixed(2)}×
              </td>
              <td className={cn("py-[9px] pl-3 text-right font-mono", sig(r.adj_p) && "font-bold text-accent-text")}>
                {pText(r.adj_p)}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function GenesView({ results }: { results: Results }) {
  const [, setGroup] = useShinyInput<string | null>("conditions_group", results.selected);
  const d = results.design;
  const alpha = results.settings.max_adj_p;
  const tested = results.groups.filter((g) => g.status === "tested");
  const skipped = results.groups.filter((g) => g.status !== "tested");
  const genes = results.genes;
  const sig = (g: GeneChange) => g.adj_p !== null && g.adj_p < alpha && g.log2_fold_change !== null;
  const plot = useMemo(() => {
    const finite = genes.filter((g) => g.p !== null && g.log2_fold_change !== null);
    const smallest = Math.min(...finite.filter((g) => (g.p ?? 0) > 0).map((g) => g.p as number), 1);
    const y = (g: GeneChange) => -Math.log10(Math.max(g.p as number, smallest / 10)); // p = 0 drawn just above the rest
    const trace = (list: GeneChange[], name: string, color: string): Data => ({
      type: "scattergl",
      mode: "markers",
      name,
      x: list.map((g) => g.log2_fold_change as number),
      y: list.map(y),
      text: list.map((g) => g.gene),
      customdata: list.map((g) => [pText(g.adj_p)]),
      marker: { size: 5, color, opacity: 0.8 },
      hovertemplate: "<b>%{text}</b><br>log2 fold change %{x:.2f}<br>Adjusted p-value %{customdata[0]}<extra></extra>",
    });
    const top = finite.filter(sig).slice(0, 10);
    return {
      data: [trace(finite.filter((g) => !sig(g)), "not significant", MUTED), trace(finite.filter(sig), "significant", TEAL)],
      layout: {
        xaxis: { title: { text: `log2 fold change (${d.test} vs ${d.reference})` }, gridcolor: GRID, zeroline: true },
        yaxis: { title: { text: "−log10 p-value" }, gridcolor: GRID },
        annotations: top.map((g) => ({
          x: g.log2_fold_change, y: y(g), text: g.gene, showarrow: false, yshift: 10, font: { size: 11, color: INK },
        })),
        legend: { orientation: "h", y: 1.08 },
      } as Partial<Layout>,
    };
  }, [genes, alpha, d.test, d.reference]); // eslint-disable-line react-hooks/exhaustive-deps
  const current = tested.find((g) => g.group === results.selected);
  if (tested.length === 0)
    return (
      <div className="flex flex-col gap-2">
        <Empty>No group could be tested for gene changes.</Empty>
        <SkippedList skipped={skipped} />
      </div>
    );
  return (
    <div className="flex flex-col gap-3">
      <label className="flex flex-wrap items-center gap-2 text-[14px] font-bold">
        Group
        <select
          id="cc_gene_group"
          value={results.selected ?? ""}
          onChange={(e) => setGroup(e.target.value)}
          className="rounded-[4px] border border-field bg-white px-2 py-1 text-[14px] font-normal"
        >
          {tested.map((g) => (
            <option key={g.group} value={g.group}>
              {g.group} · {fmt(g.n_up)} up, {fmt(g.n_down)} down
            </option>
          ))}
        </select>
      </label>
      {current && (
        <p className="m-0 text-[13px] text-body">
          {fmt(current.n_up)} genes up and {fmt(current.n_down)} down in {d.test} (adjusted p-value below {alpha}), of{" "}
          {fmt(current.n_genes_tested)} tested in {current.n_samples_used} samples
          {current.samples_left_out.length > 0 && ` (left out, too few cells: ${current.samples_left_out.join(", ")})`}.{" "}
          <InfoTip label="About the gene test">
            Each sample's cells of this group are summed into one profile (pseudobulk), and the profiles are compared with DESeq2
            {d.paired_by ? `, pairing samples by ${d.paired_by}` : ""}. A fold change is only meaningful with its adjusted p-value.
          </InfoTip>
        </p>
      )}
      <PlotlyChart
        data={plot.data}
        layout={plot.layout}
        height={420}
        label="Volcano plot: log2 fold change against significance for each gene"
        download={{ name: `Compare Conditions genes ${results.selected ?? ""}`, title: `Gene changes in ${results.selected}: ${d.test} vs ${d.reference}` }}
      />
      <div className="flex items-center justify-between gap-3">
        <span className="text-[13px] text-muted">The {Math.min(TOP_ROWS, genes.length)} genes with the smallest p-values; the CSV has all of them.</span>
        <CsvDownload
          name={`Compare Conditions genes ${results.selected ?? ""}`}
          columns={[
            { key: "gene", header: "Gene" },
            { key: "log2_fold_change", header: `log2 fold change (${d.test} vs ${d.reference})` },
            { key: "lfc_se", header: "Standard error" },
            { key: "base_mean", header: "Mean normalized count" },
            { key: "stat", header: "Wald statistic" },
            { key: "p", header: "p-value" },
            { key: "adj_p", header: "Adjusted p-value" },
          ]}
          rows={() => genes.map((g) => ({ ...g }))}
        />
      </div>
      <table className="w-full border-collapse text-[14px]">
        <thead>
          <tr className="border-b border-ink text-left text-[12px] font-bold tracking-[0.04em] text-muted uppercase">
            <th scope="col" className="py-2 pr-3">Gene</th>
            <th scope="col" className="px-3 py-2 text-right">log2 fold change</th>
            <th scope="col" className="px-3 py-2 text-right">Mean count</th>
            <th scope="col" className="py-2 pl-3 text-right">Adjusted p-value</th>
          </tr>
        </thead>
        <tbody>
          {genes.slice(0, TOP_ROWS).map((g) => (
            <tr key={g.gene} className="border-b border-line">
              <td className="py-[8px] pr-3 font-bold">{g.gene}</td>
              <td className="px-3 py-[8px] text-right font-mono">
                {g.log2_fold_change === null ? "–" : `${g.log2_fold_change > 0 ? "+" : ""}${g.log2_fold_change.toFixed(2)}`}
              </td>
              <td className="px-3 py-[8px] text-right font-mono">{g.base_mean === null ? "–" : fmt(Math.round(g.base_mean))}</td>
              <td className={cn("py-[8px] pl-3 text-right font-mono", sig(g) && "font-bold text-accent-text")}>{pText(g.adj_p)}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <SkippedList skipped={skipped} />
    </div>
  );
}

function SkippedList({ skipped }: { skipped: Results["groups"] }) {
  if (skipped.length === 0) return null;
  return (
    <div className="flex flex-col gap-1">
      <span className="text-[13px] font-bold text-ink">Not tested ({skipped.length})</span>
      {skipped.map((g) => (
        <Note key={g.group} dot="bg-accent">
          {g.group}: {g.reason}
        </Note>
      ))}
    </div>
  );
}
