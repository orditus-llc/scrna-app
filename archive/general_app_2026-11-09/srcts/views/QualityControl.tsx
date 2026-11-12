import { useShinyInput } from "@posit-dev/shinyreact";
import { ChevronDown, MoveHorizontal, RotateCcw, X } from "lucide-react";
import type { Annotation, Data, Layout, Shape } from "plotly.js-dist-min";
import { useEffect, useMemo, useRef, useState, type ReactNode } from "react";

import { Note, RuledField } from "../components/fields";
import { CsvDownload } from "../components/DownloadMenu";
import { MitoSetting } from "../components/MitoSetting";
import { MoreSection } from "../components/MoreSection";
import { StepHeader, SummaryStrip, type SummaryItem } from "../components/StepPage";
import { Tabs } from "../components/Tabs";
import { goTo, type Page } from "../navigation";
import { ActionBar } from "./stepShared";
import { PlotlyChart } from "../components/PlotlyChart";
import { Button } from "../components/ui/button";
import {
  groupByLook,
  METRIC_LABELS,
  parseRange,
  selectCells,
  summarize,
  type QcCells,
  type QcMetric,
} from "../lib/qcSelection";
import { cn } from "../lib/utils";
import {
  clickCount,
  type Click,
  type FilterBoxes,
  type FilterRow,
  type FilterSettings,
  type MitoCheck,
  type Notice,
  type QcSuggestion,
  type SuggestionDetail,
  type PreviewPayload,
  type QcPayload,
  type Step,
} from "../types";

// Colors (theme.css): kept cells dark grey, removed cells and your limits teal, suggested limits grey.
const INK = "#111111";
const KEPT = "#3f3f46";
const TEAL = "#0d9488";
const TEAL_TEXT = "#0f766e";
const SUGGESTED = "#71717a";

type OkPreview = Extract<PreviewPayload, { ok: true }>;

export function QualityControl({
  qc,
  preview,
  step,
  notice,
  mitoCheck,
  qcSuggestion,
  onNavigate,
}: {
  qc: QcPayload | undefined;
  preview: PreviewPayload | null | undefined;
  step: Step | undefined;
  notice: Notice | null;
  mitoCheck: MitoCheck | null;
  qcSuggestion: QcSuggestion | null; // thresholds for the MADs in "More settings"
  onNavigate: (page: Page) => void;
}) {
  // Click counters (the server ignores the initial 0) and the live filter boxes.
  const [runClick, setRunClick] = useShinyInput<Click<{ mito_prefix: string }>>("run_qc", 0);
  // The "Mitochondrial genes" setting; each change is sent so the page can show how many genes match.
  const [mitoPrefix, setMitoPrefix] = useState<string | null>(null);
  const [, sendMitoPrefix] = useShinyInput<string | null>("mito_prefix", null, { debounceMs: 300 });
  // How strict the suggestions are (MADs, under "More settings"). A change asks the server for
  // new suggestions, which then fill the filter boxes.
  const [mads, setMads] = useState<MadBoxes | null>(null);
  const [, sendMads] = useShinyInput<MadBoxes | null>("suggest_settings", null, { debounceMs: 400 });
  const madsEdited = useRef(false);
  // Apply sends the boxes with the click, so a value typed just before clicking is used.
  const [applyClick, setApplyClick] = useShinyInput<Click<FilterBoxes & { suggest?: MadBoxes }>>("apply_filter", 0);
  const [, sendBoxes] = useShinyInput<FilterBoxes | null>("filter_settings", null, { debounceMs: 400 });

  const [boxes, setBoxes] = useState<FilterBoxes | null>(null);
  const [busy, setBusy] = useState<"run" | "apply" | null>(null);
  // Inspecting a group of cells (looking only: it never changes the filter settings).
  const [inspect, setInspect] = useState<InspectBoxes>({ metric: "n_genes", from: "", to: "" });
  const [focus, setFocus] = useState<string | null>(null);
  const [inspectOpen, setInspectOpen] = useState(false); // Inspect is closed until used
  const [view, setView] = useState<"preview" | "distributions" | "numbers">("preview");
  const parsed = parseRange(inspect.metric, inspect.from, inspect.to);
  const inspectError = parsed && "error" in parsed ? parsed.error : null;
  const qcCells = qc?.cells as QcCells | undefined;
  const clearInspection = () => {
    setInspect((current) => ({ ...current, from: "", to: "" }));
    setFocus(null);
  };
  const inspecting = inspect.from !== "" || inspect.to !== "" || focus !== null;
  // Esc removes the highlight while one is shown.
  useEffect(() => {
    if (!inspecting) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && clearInspection();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [inspecting]); // eslint-disable-line react-hooks/exhaustive-deps -- clearInspection only uses setters
  const selected = useMemo(
    () => (qcCells && parsed && !("error" in parsed) ? selectCells(qcCells, parsed) : new Set<string>()),
    [qcCells, inspect.metric, inspect.from, inspect.to], // eslint-disable-line react-hooks/exhaustive-deps
  );

  // Fill the boxes from the current settings once, then send every change (debounced).
  useEffect(() => {
    if (qc && boxes === null) setBoxes(toBoxes(qc.filter.current));
  }, [qc, boxes]);
  useEffect(() => {
    if (qc && mitoPrefix === null) setMitoPrefix(qc.settings.mito_prefix);
  }, [qc, mitoPrefix]);
  useEffect(() => {
    if (mitoPrefix !== null) sendMitoPrefix(mitoPrefix);
  }, [mitoPrefix]); // eslint-disable-line react-hooks/exhaustive-deps -- sendMitoPrefix is stable
  // When the settings in use change (QC ran and suggested thresholds), re-read the filter boxes.
  const currentKey = qc ? JSON.stringify(qc.filter.current) : null;
  useEffect(() => {
    if (qc) setBoxes(toBoxes(qc.filter.current));
  }, [currentKey]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => {
    const s = qc?.filter.suggest_settings;
    if (s && mads === null) setMads({ n_mads_genes: String(s.n_mads_genes), n_mads_mito: String(s.n_mads_mito) });
  }, [qc, mads]);
  useEffect(() => {
    if (mads && madsEdited.current) sendMads(mads);
  }, [mads]); // eslint-disable-line react-hooks/exhaustive-deps -- sendMads is stable
  // New suggestions for edited MADs fill the filter boxes.
  useEffect(() => {
    if (madsEdited.current && qcSuggestion?.suggested) setBoxes(toBoxes(qcSuggestion.suggested));
  }, [qcSuggestion]); // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => {
    if (boxes) sendBoxes(boxes);
  }, [boxes]); // eslint-disable-line react-hooks/exhaustive-deps -- sendBoxes is stable
  // Stop "busy" when the server answers (new results or a message).
  useEffect(() => setBusy(null), [qc, notice?.id]);

  if (!qc) return <p className="m-0 px-8 py-10 text-muted">Loading…</p>;
  // The thresholds suggested from the data: for the edited MADs if any, else from QC.
  const suggestion: QcSuggestion = madsEdited.current && qcSuggestion?.suggested ? qcSuggestion : qc.filter;
  const suggested = suggestion.suggested ?? qc.filter.current;
  const runQc = () => {
    setBusy("run");
    setRunClick({ n: clickCount(runClick) + 1, settings: { mito_prefix: mitoPrefix ?? qc.settings.mito_prefix } });
  };
  const mitoSetting = mitoPrefix !== null && (
    <MitoSetting settings={qc.settings} value={mitoPrefix} onChange={setMitoPrefix} check={mitoCheck} />
  );
  // Filtering is applied with exactly the settings in the boxes: nothing to apply.
  const appliedAsShown = qc.filter_status === "done" && boxes !== null && sameFilter(boxes, toBoxes(qc.filter.current));
  const error = notice?.step === "quality_control" ? notice.text : null;
  const ok = preview && preview.ok ? preview : null;

  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <StepHeader step={step} />

      {!qc.metrics_done ? (
        <section className="flex flex-col gap-3 rounded-[4px] border-2 border-ink p-5">
          <h2 className="m-0 text-[16px] font-bold">Measure every cell</h2>
          <p className="m-0 max-w-[760px] text-[14px] text-body">
            Flags mitochondrial genes and calculates, for every cell, the number of genes detected, the total counts and
            the percentage of counts from mitochondrial genes. Damaged cells tend to have a high mitochondrial share; empty
            droplets have very few genes; doublets have unusually many.
          </p>
          {qc.can_run && <div className="max-w-[420px]">{mitoSetting}</div>}
          {error && <Note dot="bg-accent">{error}</Note>}
          <div className="flex flex-wrap items-center gap-3">
            <Button disabled={!qc.can_run || busy === "run" || !mitoPrefix?.trim()} onClick={runQc}>
              {busy === "run" ? "Running QC…" : "Run QC"}
            </Button>
            {!qc.can_run && (
              <span className="text-[13px] text-muted">{step?.prerequisite_text ?? "Run Data Overview first."}</span>
            )}
          </div>
        </section>
      ) : (
        <>
          {/* 1. The key numbers. */}
          <Summary qc={qc} preview={ok} />

          {/* 2. Settings (left) and one view at a time (right). Side by side only when the
              content area is wide enough (a container query, not the window width); the
              settings stick while scrolling only then, so they never float over the charts. */}
          <div className="@container">
          <div className="flex flex-col gap-10 @min-[900px]:flex-row">
            <FilterPanel
              boxes={boxes}
              setBoxes={setBoxes}
              suggested={suggested}
              details={suggestion.details}
              mitoAvailable={qc.mito?.available ?? true}
              mads={mads}
              madsProblem={madsEdited.current ? qcSuggestion?.problem ?? null : null}
              onMads={(next) => {
                madsEdited.current = true;
                setMads(next);
              }}
              mito={
                <>
                  {mitoSetting}
                  {mitoPrefix !== null && mitoPrefix !== qc.mito?.prefix && (
                    <div className="flex flex-col items-start gap-2 border-b border-line py-3">
                      <Button size="small" disabled={busy === "run" || !mitoPrefix.trim()} onClick={runQc}>
                        {busy === "run" ? "Running QC…" : "Re-run QC with this pattern"}
                      </Button>
                      <span className="text-[13px] text-muted">Re-running QC clears filtering and the later steps.</span>
                    </div>
                  )}
                </>
              }
            />
            <div className="flex min-w-0 flex-1 flex-col">
              <Tabs
                label="Quality control views"
                tabs={[
                  { key: "preview", label: "Preview" },
                  { key: "distributions", label: "Distributions" },
                  { key: "numbers", label: "Numbers" },
                ]}
                active={view}
                onChange={setView}
              >
                {view === "preview" && (
                  <>
                    <PreviewCharts
                      qc={qc}
                      suggested={suggested}
                      preview={preview ?? null}
                      selected={selected}
                      focus={focus}
                      onClear={clearInspection}
                      onDragRange={(metric, from, to) => {
                        const round = (v: number) => (metric === "pct_mt" ? String(Math.round(v * 100) / 100) : String(Math.max(0, Math.round(v))));
                        setInspect({ metric, from: round(from), to: round(to) });
                        setInspectOpen(true);
                      }}
                    />
                    <InspectSection
                      qc={qc}
                      preview={preview ?? null}
                      boxes={inspect}
                      setBoxes={setInspect}
                      selected={selected}
                      error={inspectError}
                      focus={focus}
                      setFocus={setFocus}
                      onClear={clearInspection}
                      open={inspectOpen || inspecting}
                      onToggle={() => {
                        if (inspectOpen || inspecting) {
                          clearInspection();
                          setInspectOpen(false);
                        } else setInspectOpen(true);
                      }}
                    />
                  </>
                )}
                {view === "distributions" && <MetricsSection qc={qc} />}
                {view === "numbers" && (
                  <div className="flex flex-col gap-6">
                    {ok && <FilterTable preview={ok} />}
                    <StatsSection qc={qc} />
                    <MitoLine qc={qc} />
                  </div>
                )}
              </Tabs>
            </div>
          </div>
          </div>

          {/* 3. The actions, always on screen. */}
          <QcActionBar
            preview={ok}
            previewProblem={preview && !preview.ok ? preview.message : null}
            error={error}
            confirmation={qc.applied && qc.filter_status === "done" ? qc.applied.confirmation : null}
            needsRerun={qc.filter_status === "needs_rerun"}
            runLabel={busy === "apply" ? "Applying…" : `${qc.filter_status === "not_run" ? "Apply" : "Re-apply"} Filtering`}
            runDisabled={busy === "apply" || !ok || appliedAsShown}
            upToDate={appliedAsShown}
            willClear={!appliedAsShown && busy !== "apply" ? qc.will_clear : []}
            onRun={() => {
              if (!boxes) return;
              setBusy("apply");
              setApplyClick({ n: clickCount(applyClick) + 1, settings: { ...boxes, suggest: mads ?? undefined } });
            }}
            onContinue={() => goTo(onNavigate, "preprocessing")}
          />
        </>
      )}
    </div>
  );
}

// --- Pieces of the page -------------------------------------------------------------------

/** The headline numbers between two rules: what these settings keep, compared with what the
 * suggested thresholds keep. */
function Summary({ qc, preview }: { qc: QcPayload; preview: OkPreview | null }) {
  const vsSuggested = (kept: number, suggested: number) =>
    kept === suggested
      ? { text: "= suggested", differs: false }
      : { text: `${kept > suggested ? "+" : "−"}${fmt(Math.abs(kept - suggested))} vs suggested (${fmt(suggested)})`, differs: true };

  const items: SummaryItem[] = [];
  if (preview) {
    const cells = vsSuggested(preview.n_cells_kept, preview.suggested.n_cells_kept);
    const genes = vsSuggested(preview.n_genes_kept, preview.suggested.n_genes_kept);
    const removed = preview.n_cells_before - preview.n_cells_kept;
    items.push(
      { label: "Cells kept", value: fmt(preview.n_cells_kept), note: `of ${fmt(preview.n_cells_before)} · ${cells.text}`, differs: cells.differs },
      { label: "Genes kept", value: fmt(preview.n_genes_kept), note: `of ${fmt(preview.n_genes_before)} · ${genes.text}`, differs: genes.differs },
      { label: "Cells removed", value: fmt(removed), note: `${((100 * removed) / preview.n_cells_before).toFixed(1)}% of cells` },
    );
  } else {
    items.push(
      { label: "Cells kept", value: "–", note: `of ${fmt(qc.n_cells)} · check the settings` },
      { label: "Genes kept", value: "–", note: `of ${fmt(qc.n_genes)} · check the settings` },
      { label: "Cells removed", value: "–", note: "check the settings" },
    );
  }
  return <SummaryStrip items={items} />;
}

function MitoLine({ qc }: { qc: QcPayload }) {
  if (!qc.mito) return null;
  return (
    <div className="flex flex-col gap-2">
      {qc.mito.warnings.map((w) => (
        <Note key={w} dot="bg-accent">
          {w}
        </Note>
      ))}
      {qc.mito.available && (
        <p className="m-0 text-[14px] text-body">
          <strong>{qc.mito.genes.length} mitochondrial genes</strong> (names starting with "{qc.mito.prefix}"):{" "}
          <span className="font-mono text-[13px]">{qc.mito.genes.join(", ")}</span>
        </p>
      )}
    </div>
  );
}

/** The four filters as a ruled list: label and explanation on the left, the box on the right. */
function FilterPanel({
  boxes,
  setBoxes,
  suggested,
  details,
  mitoAvailable,
  mito,
  mads,
  madsProblem,
  onMads,
}: {
  boxes: FilterBoxes | null;
  setBoxes: (b: FilterBoxes) => void;
  suggested: FilterSettings; // the thresholds suggested from the data
  details: QcSuggestion["details"]; // how each was found (for the ⓘ)
  mitoAvailable: boolean;
  mito: ReactNode; // the "Mitochondrial genes" setting, first in the list
  mads: MadBoxes | null;
  madsProblem: string | null;
  onMads: (mads: MadBoxes) => void;
}) {
  if (!boxes) return null;
  const reference = suggested;
  const set = (patch: Partial<FilterBoxes>) => setBoxes({ ...boxes, ...patch });
  return (
    <section aria-label="Filter settings" className="flex w-full min-w-0 max-w-[360px] flex-none flex-col self-start @min-[900px]:sticky @min-[900px]:top-4 @min-[900px]:w-[320px]">
      <h2 className="m-0 border-b border-field pb-[10px] text-[16px] font-bold">Filter settings</h2>
      {mito}
      <RuledField
        id="min_genes"
        label="Min genes per cell"
        help={"Removes cells with fewer genes than this, which are often empty droplets." + explain(details?.min_genes, "below")}
        value={boxes.min_genes}
        onChange={(v) => set({ min_genes: v })}
        step={10}
        reference={reference.min_genes}
        compareLabel="Suggested"
      />
      <RuledField
        id="min_cells"
        label="Min cells per gene"
        help="Removes genes found in fewer cells than this. Suggested: 3, a common convention (Seurat tutorial), not calculated from the data."
        value={boxes.min_cells}
        onChange={(v) => set({ min_cells: v })}
        step={1}
        reference={reference.min_cells}
        compareLabel="Suggested"
      />
      <RuledField
        id="max_genes"
        label="Max genes per cell"
        help={"Removes cells with this many genes or more, which are often doublets." + explain(details?.max_n_genes_by_counts, "above")}
        value={boxes.max_genes}
        onChange={(v) => set({ max_genes: v })}
        step={50}
        reference={reference.max_n_genes_by_counts}
        compareLabel="Suggested"
        off={boxes.max_genes_off}
        onOff={(off) => set({ max_genes_off: off })}
      />
      <RuledField
        id="max_mt"
        label="Max % mitochondrial"
        help={"Removes cells at or above this share of mitochondrial counts, which are often damaged." + explain(details?.max_pct_counts_mt, "above")}
        note={mitoAvailable ? undefined : "Unavailable: no mitochondrial genes found, so this filter is off."}
        value={boxes.max_mt}
        onChange={(v) => set({ max_mt: v })}
        step={0.5}
        unit="%"
        reference={reference.max_pct_counts_mt}
        compareLabel="Suggested"
        off={boxes.max_mt_off}
        onOff={(off) => set({ max_mt_off: off })}
      />
      {mads && (
        <MoreSection title="More settings" note={madsProblem ? "check" : undefined}>
          <p className="m-0 pt-2 text-[13px] text-body">
            How far from your data's median a cell must be to be suggested for removal, in median absolute deviations
            (MADs). Changing these updates the suggested thresholds above.
          </p>
          <RuledField
            id="n_mads_genes"
            label="MADs for genes per cell"
            help="Distance from the median number of genes per cell (on a log scale), below and above. Larger = keeps more cells. Default: 5, as in the single-cell best-practices guide."
            value={mads.n_mads_genes}
            onChange={(v) => onMads({ ...mads, n_mads_genes: v })}
            step={0.5}
            reference={5}
          />
          <RuledField
            id="n_mads_mito"
            label="MADs for % mitochondrial"
            help="Distance above the median % mitochondrial. Larger = keeps more cells. Default: 3, as in the single-cell best-practices guide."
            value={mads.n_mads_mito}
            onChange={(v) => onMads({ ...mads, n_mads_mito: v })}
            step={0.5}
            reference={3}
          />
          {madsProblem && <Note dot="bg-accent">{madsProblem}</Note>}
        </MoreSection>
      )}
      <div className="flex flex-col items-start gap-3 pt-[18px]">
        <Button variant="secondary" size="small" onClick={() => setBoxes(toBoxes(reference))}>
          <RotateCcw aria-hidden="true" size={14} /> Reset to suggested settings
        </Button>
      </div>
    </section>
  );
}

/** The MAD boxes under "More settings", as typed. */
export interface MadBoxes {
  n_mads_genes: string;
  n_mads_mito: string;
}

/** One sentence for a field's ⓘ: how its suggested value was found in this data. */
function explain(d: SuggestionDetail | undefined, side: "below" | "above"): string {
  if (!d) return "";
  if (d.value === null || d.median === null || d.mad === null) return d.note ? ` ${d.note}` : "";
  if (d.scale === "log1p") {
    const median = Math.round(Math.expm1(d.median));
    return ` Suggested from your data: ${d.n_mads} MADs ${side} the median of ${fmt(median)} genes per cell (log scale).`;
  }
  return ` Suggested from your data: ${d.n_mads} MADs ${side} the median of ${d.median.toFixed(1)}% (MAD ${d.mad.toFixed(1)}%), rounded up.`;
}


// Filters as people read them (keys come from the server's filter table).
const FILTER_NAMES: Record<string, string> = {
  min_genes: "min genes",
  max_n_genes_by_counts: "max genes",
  max_pct_counts_mt: "max % mitochondrial",
};

/** The inspection a user is making: which measure and the typed range ("" = open end). */
export interface InspectBoxes {
  metric: QcMetric;
  from: string;
  to: string;
}

/** Both QC measures as dot charts, one above the other, every cell in each. Dark = kept,
 * teal = removed by these settings. An inspected group of cells (chosen by a range, or by
 * dragging across either chart) is drawn on top with a dark ring in BOTH charts and the
 * rest fade; cells are matched by barcode. Lines: your limits teal, the reference dashed grey. */
function PreviewCharts({
  qc,
  suggested,
  preview,
  selected,
  focus,
  onClear,
  onDragRange,
}: {
  qc: QcPayload;
  suggested: FilterSettings;
  preview: PreviewPayload | null;
  selected: Set<string>;
  focus: string | null;
  onClear: () => void;
  onDragRange: (metric: QcMetric, from: number, to: number) => void;
}) {
  const mitoAvailable = qc.mito?.available ?? true;
  const ok = preview && preview.ok ? preview : null;
  const removed = useMemo(() => new Set(ok ? ok.removed_barcodes : []), [ok]);
  const cells = qc.cells as QcCells;
  const nSelected = selected.size;
  const nRemovedSelected = useMemo(() => [...selected].filter((b) => removed.has(b)).length, [selected, removed]);

  return (
    <section aria-label="Preview" className="flex min-w-0 flex-[999_1_560px] flex-col">
      <p className="m-0 flex items-start gap-[10px] rounded-[4px] bg-accent-tint px-3 py-[9px] text-[14px] leading-snug text-ink">
        <MoveHorizontal aria-hidden="true" size={18} className="mt-[1px] flex-none text-accent-text" />
        <span>
          <b>Drag sideways across a chart</b> to highlight those cells in both charts. Double-click or press Esc to clear.
        </span>
      </p>
      {preview && !preview.ok && (
        <div className="mt-3">
          <Note dot="bg-accent">{preview.message}</Note>
        </div>
      )}
      <p className="mt-[10px] mb-0 text-[13px] text-body">
        Each dot is a cell: dark = kept, teal = removed.
      </p>
      {(nSelected > 0 || focus !== null) && (
        <div className="mt-3 flex flex-wrap items-center gap-3 text-[14px]">
          <span className="font-bold text-ink">
            {nSelected > 0 ? `${fmt(nSelected)} cells highlighted` : "1 cell marked"}
          </span>
          <Button variant="secondary" size="small" onClick={onClear}>
            <X aria-hidden="true" size={14} /> Clear highlight
          </Button>
          <span className="text-[13px] text-muted">or double-click a chart, or press Esc</span>
        </div>
      )}
      {mitoAvailable && (
        <QcStrip metric="pct_mt" qc={qc} settings={ok ? ok.settings : qc.filter.current} suggested={suggested} removed={removed} selected={selected} focus={focus} onDragRange={onDragRange} onClear={onClear} />
      )}
      <QcStrip metric="n_genes" qc={qc} settings={ok ? ok.settings : qc.filter.current} suggested={suggested} removed={removed} selected={selected} focus={focus} onDragRange={onDragRange} onClear={onClear} />
      <ul aria-label="Legend" className="m-0 flex list-none flex-wrap gap-x-[20px] gap-y-2 p-0 text-[13px] text-body">
        <li className="inline-flex items-center gap-2">
          <span aria-hidden="true" className="size-[9px] rounded-full bg-[#3f3f46]" />
          kept · {fmt(cells.barcode.length - removed.size)}
        </li>
        <li className="inline-flex items-center gap-2 font-bold text-accent-text">
          <span aria-hidden="true" className="size-[9px] rounded-full bg-accent" />
          removed by these settings · {fmt(removed.size)}
        </li>
        {nSelected > 0 && (
          <>
            <li className="inline-flex items-center gap-2 font-bold text-ink">
              <span aria-hidden="true" className="size-[11px] rounded-full border-2 border-ink bg-[#3f3f46]" />
              inspected, kept · {fmt(nSelected - nRemovedSelected)}
            </li>
            <li className="inline-flex items-center gap-2 font-bold text-ink">
              <span aria-hidden="true" className="size-[11px] rounded-full border-2 border-ink bg-accent" />
              inspected, would be removed · {fmt(nRemovedSelected)}
            </li>
          </>
        )}
      </ul>
      {nSelected > 0 && (
        <p className="mt-2 mb-0 text-[13px] text-muted">
          Being inspected doesn't mark a cell as low quality.
        </p>
      )}
    </section>
  );
}

function QcStrip({
  metric,
  qc,
  settings,
  suggested,
  removed,
  selected,
  focus,
  onDragRange,
  onClear,
}: {
  metric: "pct_mt" | "n_genes";
  qc: QcPayload;
  settings: FilterSettings;
  suggested: FilterSettings; // the thresholds suggested from the data (dashed lines)
  removed: Set<string>;
  selected: Set<string>;
  focus: string | null;
  onDragRange: (metric: QcMetric, from: number, to: number) => void;
  onClear: () => void;
}) {
  const cells = qc.cells as QcCells;
  const mitoAvailable = qc.mito?.available ?? true;

  const chart = useMemo(() => {
    const values = metric === "pct_mt" ? cells.pct_mt : cells.n_genes;
    const groups = groupByLook(cells, selected, removed);
    const inspecting = selected.size > 0;
    const trace = (idx: number[], name: string, marker: object): Data => ({
      type: "scattergl",
      mode: "markers",
      name,
      x: idx.map((i) => values[i]),
      y: idx.map((i) => spread(i)), // the same height for a cell in both charts
      customdata: idx.map((i) => [cells.barcode[i], cells.n_genes[i], cells.total_counts[i], cells.pct_mt[i]]),
      marker,
      hovertemplate:
        "<b>%{customdata[0]}</b><br>Genes %{customdata[1]:,}<br>Counts %{customdata[2]:,.0f}" +
        (mitoAvailable ? "<br>% mito %{customdata[3]:.2f}" : "") +
        `<extra>${name}</extra>`,
    });
    const ring = { color: INK, width: 1.6 };
    const data: Data[] = [
      trace(groups.kept, "kept", { size: 5, color: KEPT, opacity: inspecting ? 0.14 : 0.4 }),
      trace(groups.removed, "removed", { size: 7.5, color: TEAL, opacity: inspecting ? 0.3 : 1, line: { color: "#ffffff", width: 1 } }),
      trace(groups.selected_kept, "inspected, kept", { size: 8, color: KEPT, opacity: 0.9, line: ring }),
      trace(groups.selected_removed, "inspected, would be removed", { size: 9, color: TEAL, opacity: 1, line: ring }),
    ];
    const focusIndex = focus ? cells.barcode.indexOf(focus) : -1;
    if (focusIndex >= 0) {
      data.push(trace([focusIndex], "this cell", { size: 18, color: "rgba(0,0,0,0)", line: { color: INK, width: 2.5 } }));
    }

    // Lines: your limits (teal) and the suggested ones (dashed grey); one label when they're equal.
    const lines: { x: number; label: string; yours: boolean }[] = [];
    const add = (yours: number | null, ref: number | null, name: string, unit: string) => {
      if (yours !== null) lines.push({ x: yours, label: `${name}: ${fmt(yours)}${unit}${yours === ref ? " = suggested" : ""}`, yours: true });
      if (ref !== null && ref !== yours) lines.push({ x: ref, label: `Suggested: ${fmt(ref)}${unit}`, yours: false });
    };
    if (metric === "pct_mt") add(settings.max_pct_counts_mt, suggested.max_pct_counts_mt, "Max % mitochondrial", "%");
    else {
      add(settings.min_genes, suggested.min_genes, "Min genes", "");
      add(settings.max_n_genes_by_counts, suggested.max_n_genes_by_counts, "Max genes", "");
    }
    const shapes: Partial<Shape>[] = lines.map((l) => ({
      type: "line", xref: "x", yref: "paper", x0: l.x, x1: l.x, y0: 0, y1: 1,
      line: l.yours ? { color: TEAL, width: 2 } : { color: SUGGESTED, width: 1.5, dash: "dash" },
    }));
    const annotations: Partial<Annotation>[] = lines.map((l) => ({
      x: l.x, xref: "x", yref: "paper", y: l.yours ? 1.02 : 0.96, yanchor: l.yours ? "bottom" : "top",
      xanchor: "left", xshift: 4, showarrow: false, text: l.label,
      font: { size: 12, color: l.yours ? TEAL_TEXT : "#52525b" }, bgcolor: "rgba(255,255,255,0.85)",
    }));
    return {
      data,
      layout: {
        margin: { l: 16, r: 16, t: 34, b: 48 },
        dragmode: "select",
        selectdirection: "h",
        xaxis: { title: { text: METRIC_AXIS[metric] }, gridcolor: "#ececef", zeroline: false, rangemode: "tozero" },
        yaxis: { visible: false, range: [-0.04, 1.04], fixedrange: true },
        shapes,
        annotations,
        showlegend: false,
      } as Partial<Layout>,
    };
  }, [metric, cells, removed, selected, focus, settings, suggested, mitoAvailable]);

  return (
    <div className="mt-3 flex flex-col">
      <h3 className="m-0 text-[14px] font-extrabold">{METRIC_LABELS[metric]}</h3>
      <PlotlyChart
        data={chart.data}
        layout={chart.layout}
        height={250}
        label={`Every cell by ${METRIC_LABELS[metric]}; removed cells in teal; inspected cells ringed`}
        download={{
          name: `Quality Control preview ${METRIC_LABELS[metric]}`,
          title: `Every cell by ${METRIC_LABELS[metric]}: dark = kept, teal = removed by the current filter settings`,
        }}
        onSelectRange={(x0, x1) => onDragRange(metric, x0, x1)}
        onDeselect={onClear}
      />
    </div>
  );
}

const METRIC_AXIS: Record<"pct_mt" | "n_genes", string> = {
  pct_mt: "% of each cell's counts from mitochondrial genes",
  n_genes: "Genes detected per cell",
};

/** Choose a group of cells by a range on one QC measure, and see what's in it. */
function InspectSection({
  qc,
  preview,
  boxes,
  setBoxes,
  selected,
  error,
  focus,
  setFocus,
  onClear,
  open,
  onToggle,
}: {
  qc: QcPayload;
  preview: PreviewPayload | null;
  boxes: InspectBoxes;
  setBoxes: (b: InspectBoxes) => void;
  selected: Set<string>;
  error: string | null;
  focus: string | null;
  setFocus: (barcode: string | null) => void;
  onClear: () => void;
  open: boolean;
  onToggle: () => void; // open, or close (closing also clears the highlight)
}) {
  const cells = qc.cells as QcCells;
  const ok = preview && preview.ok ? preview : null;
  const failures = ok ? ok.failures : {};
  const summary = useMemo(() => summarize(cells, selected, failures), [cells, selected, failures]);
  const mitoAvailable = qc.mito?.available ?? true;
  const metrics = (Object.keys(METRIC_LABELS) as QcMetric[]).filter((m) => m !== "pct_mt" || mitoAvailable);
  const [sortBy, setSortBy] = useState<QcMetric>(boxes.metric);
  const shown = useMemo(() => {
    const idx = cells.barcode.flatMap((b, i) => (selected.has(b) ? [i] : []));
    return idx.sort((a, b) => cells[sortBy][a] - cells[sortBy][b]);
  }, [cells, selected, sortBy]);
  const LIST = 50;

  const spreadRow = (m: QcMetric) => {
    const s = summary.spread[m];
    const f = (v: number) => (m === "pct_mt" ? v.toFixed(2) : fmt(Math.round(v)));
    return (
      <tr key={m} className="border-b border-line">
        <td className="py-[7px] pr-3">{METRIC_LABELS[m]}</td>
        <td className="px-3 py-[7px] text-right font-mono">{s ? f(s.min) : "–"}</td>
        <td className="px-3 py-[7px] text-right font-mono">{s ? f(s.median) : "–"}</td>
        <td className="py-[7px] pl-3 text-right font-mono">{s ? f(s.max) : "–"}</td>
      </tr>
    );
  };
  const head = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";

  return (
    <section aria-label="Inspect a group of cells" className="flex flex-col gap-3 rounded-[4px] border border-line bg-panel px-4 py-3">
      <button
        type="button"
        aria-expanded={open}
        onClick={onToggle}
        className="flex w-full flex-wrap items-center justify-between gap-2 text-left"
      >
        <span className="flex items-center gap-2">
          <ChevronDown aria-hidden="true" size={16} className={cn("transition-transform", !open && "-rotate-90")} />
          <span className="text-[14px] font-extrabold text-ink">Inspect cells</span>
          {!open && <span className="text-[13px] text-muted">Type a range or drag across a chart to look at a group of cells</span>}
        </span>
        <span className="text-[13px] font-bold text-muted">{open ? "Close inspection" : "Inspection only: nothing is filtered"}</span>
      </button>
      {open && (
      <>
      <p className="m-0 max-w-[780px] text-[14px] text-body">
        Type a range or drag across a chart to look at a group of cells, e.g. 200–300 genes. Nothing is filtered.
      </p>
      <form
        className="flex flex-wrap items-end gap-3"
        onSubmit={(e) => e.preventDefault()}
        aria-describedby={error ? "inspect_error" : undefined}
      >
        <label className="flex flex-col gap-1 text-[13px] font-bold">
          Measure
          <select
            id="inspect_metric"
            value={boxes.metric}
            onChange={(e) => {
              setBoxes({ metric: e.target.value as QcMetric, from: "", to: "" });
              setSortBy(e.target.value as QcMetric);
            }}
            className="rounded-[4px] border border-field bg-white px-[10px] py-2 text-[14px] font-normal"
          >
            {metrics.map((m) => (
              <option key={m} value={m}>
                {METRIC_LABELS[m]}
              </option>
            ))}
          </select>
        </label>
        <label className="flex flex-col gap-1 text-[13px] font-bold">
          From (included)
          <input
            id="inspect_from"
            type="number"
            value={boxes.from}
            onChange={(e) => setBoxes({ ...boxes, from: e.target.value })}
            placeholder="no lower end"
            className="w-[150px] rounded-[4px] border border-field bg-white px-[10px] py-2 text-right font-mono text-[14px]"
          />
        </label>
        <label className="flex flex-col gap-1 text-[13px] font-bold">
          To (not included)
          <input
            id="inspect_to"
            type="number"
            value={boxes.to}
            onChange={(e) => setBoxes({ ...boxes, to: e.target.value })}
            placeholder="no upper end"
            className="w-[150px] rounded-[4px] border border-field bg-white px-[10px] py-2 text-right font-mono text-[14px]"
          />
        </label>
        <Button variant="secondary" size="small" disabled={boxes.from === "" && boxes.to === "" && focus === null} onClick={onClear}>
          <X aria-hidden="true" size={14} /> Clear highlight
        </Button>
      </form>
      {error && (
        <p id="inspect_error" className="m-0 text-[13px] font-bold text-accent-text">
          {error}
        </p>
      )}

      {selected.size > 0 && (
        <div className="flex flex-wrap gap-8 pt-1">
          <div className="flex min-w-0 flex-[1_1_300px] flex-col gap-3">
            <p className="m-0 text-[14px]">
              <b className="font-mono text-[22px]">{fmt(summary.n)}</b>{" "}
              <span className="text-body">
                of {fmt(cells.barcode.length)} cells inspected ({((100 * summary.n) / cells.barcode.length).toFixed(1)}%)
              </span>
            </p>
            <table className="w-full border-collapse text-[14px]">
              <thead>
                <tr className="border-b border-ink text-left">
                  <th scope="col" className={cn(head, "pr-3")}>Measure</th>
                  <th scope="col" className={cn(head, "px-3 text-right")}>Min</th>
                  <th scope="col" className={cn(head, "px-3 text-right")}>Median</th>
                  <th scope="col" className={cn(head, "pl-3 text-right")}>Max</th>
                </tr>
              </thead>
              <tbody>{metrics.map(spreadRow)}</tbody>
            </table>
            <div className="flex flex-col gap-1 text-[14px]">
              <span className="font-bold text-ink">Under the current settings</span>
              {!ok ? (
                <span className="text-muted">Fix the filter settings to see which of these cells would be removed.</span>
              ) : summary.nRemoved === 0 ? (
                <span className="text-body">None of these cells would be removed.</span>
              ) : (
                <>
                  <span className="font-bold text-accent-text">
                    {fmt(summary.nRemoved)} of {fmt(summary.n)} would be removed
                  </span>
                  <ul className="m-0 list-none p-0 text-body">
                    {Object.entries(summary.failing).map(([name, n]) => (
                      <li key={name}>
                        fail {FILTER_NAMES[name] ?? name}: <span className="font-mono">{fmt(n)}</span>
                      </li>
                    ))}
                  </ul>
                  {Object.keys(summary.failing).length > 1 && (
                    <span className="text-[13px] text-muted">A cell can fail more than one filter, so these can add up to more.</span>
                  )}
                </>
              )}
              <span className="text-[13px] text-muted">{fmt(summary.n - summary.nRemoved)} would be kept.</span>
            </div>
          </div>

          <div className="flex min-w-0 flex-[1_1_100%] flex-col gap-2">
            <div className="overflow-x-auto">
              <table className="w-full min-w-[520px] border-collapse text-[14px]">
                <thead>
                  <tr className="border-b border-ink text-left">
                    <th scope="col" className={cn(head, "pr-3")}>Cell</th>
                    {metrics.map((m) => (
                      <th key={m} scope="col" className={cn(head, "px-3 text-right")} aria-sort={sortBy === m ? "ascending" : "none"}>
                        <button type="button" onClick={() => setSortBy(m)} className={cn("uppercase", sortBy === m && "text-ink underline underline-offset-4")}>
                          {m === "n_genes" ? "Genes" : m === "total_counts" ? "Counts" : "% mitochondrial"}
                        </button>
                      </th>
                    ))}
                    <th scope="col" className={cn(head, "pl-3")}>Under the current settings</th>
                  </tr>
                </thead>
                <tbody>
                  {shown.slice(0, LIST).map((i) => {
                    const barcode = cells.barcode[i];
                    const fails = failures[barcode] ?? [];
                    return (
                      <tr
                        key={barcode}
                        onClick={() => setFocus(focus === barcode ? null : barcode)}
                        className={cn("cursor-pointer border-b border-line", focus === barcode ? "bg-chip" : "hover:bg-panel")}
                      >
                        <td className="py-[7px] pr-3">
                          <button
                            type="button"
                            aria-pressed={focus === barcode}
                            onClick={(e) => {
                              e.stopPropagation();
                              setFocus(focus === barcode ? null : barcode);
                            }}
                            className="font-mono text-[13px] whitespace-nowrap"
                          >
                            {barcode}
                          </button>
                        </td>
                        {metrics.map((m) => (
                          <td key={m} className="px-3 py-[7px] text-right font-mono">
                            {m === "pct_mt" ? cells.pct_mt[i].toFixed(2) : fmt(cells[m][i])}
                          </td>
                        ))}
                        <td className={cn("py-[7px] pl-3", fails.length ? "font-bold text-accent-text" : "text-body")}>
                          {fails.length ? `removed: ${fails.map((f) => FILTER_NAMES[f] ?? f).join(", ")}` : "kept"}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
            <div className="flex flex-wrap items-start justify-between gap-3">
              <span className="text-[13px] text-muted">
                {shown.length > LIST ? `Showing the first ${LIST} of ${fmt(shown.length)}, sorted by ${METRIC_LABELS[sortBy].toLowerCase()}. ` : ""}
                Click a cell to mark it in both charts; click it again to unmark it.
              </span>
              <CsvDownload
                name="Quality Control inspected cells"
                columns={[
                  { key: "barcode", header: "Cell" },
                  { key: "n_genes", header: "Genes" },
                  { key: "total_counts", header: "Counts" },
                  { key: "pct_mt", header: "% mitochondrial" },
                  { key: "status", header: "Under the current settings" },
                ]}
                rows={() =>
                  shown.map((i) => {
                    const fails = failures[cells.barcode[i]] ?? [];
                    return {
                      barcode: cells.barcode[i],
                      n_genes: cells.n_genes[i],
                      total_counts: cells.total_counts[i],
                      pct_mt: cells.pct_mt[i],
                      status: fails.length ? `removed: ${fails.map((f) => FILTER_NAMES[f] ?? f).join(", ")}` : "kept",
                    };
                  })
                }
              />
            </div>
          </div>
        </div>
      )}
      {selected.size === 0 && !error && (boxes.from !== "" || boxes.to !== "") && (
        <p className="m-0 text-[14px] text-muted">No cells in this range.</p>
      )}
      {selected.size === 0 && !error && boxes.from === "" && boxes.to === "" && (
        <p className="m-0 rounded-[4px] bg-white px-4 py-5 text-center text-[14px] text-muted">
          No cells selected. Type a range or drag across a chart.
        </p>
      )}
      </>
      )}
    </section>
  );
}

/** What each filter removes, in the order they run, next to the reference. */
function FilterTable({ preview }: { preview: OkPreview }) {
  const setting = (row: FilterRow, value: number | null) =>
    value === null ? "no limit" : `${fmt(value)}${row.key === "max_pct_counts_mt" ? "%" : ""}`;
  const total = (row: FilterRow) => (row.what === "cells" ? preview.n_cells_before : preview.n_genes_before);
  const header = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  return (
    <section aria-label="Removed by each filter" className="flex flex-col gap-3">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <h2 className="m-0 text-[16px] font-bold">Removed by each filter</h2>
        <CsvDownload
          name="Quality Control removed by each filter"
          columns={[
            { key: "order", header: "Order" },
            { key: "filter", header: "Filter" },
            { key: "what", header: "Removes" },
            { key: "setting", header: "Setting" },
            { key: "suggested", header: "Suggested" },
            { key: "removed", header: "Removed" },
            { key: "suggested_removed", header: "Removed at the suggestion" },
            { key: "out_of", header: "Out of" },
          ]}
          rows={() =>
            preview.filters.map((row, i) => ({
              order: i + 1,
              filter: row.reason,
              what: row.what,
              setting: row.setting === null ? "no limit" : row.setting,
              suggested: row.suggested_setting === null ? "no limit" : row.suggested_setting,
              removed: row.removed,
              suggested_removed: row.suggested_removed,
              out_of: total(row),
            }))
          }
        />
      </div>
      <div className="overflow-x-auto">
        <table className="w-full min-w-[760px] border-collapse text-[14px]">
          <thead>
            <tr className="border-b border-ink text-left">
              <th scope="col" className={cn(header, "pr-3")}>Filter</th>
              <th scope="col" className={cn(header, "px-3")}>Removes</th>
              <th scope="col" className={cn(header, "px-3 text-right")}>Setting</th>
              <th scope="col" className={cn(header, "px-3 text-right")}>Suggested</th>
              <th scope="col" className={cn(header, "px-3 text-right")}>Removed</th>
              <th scope="col" className={cn(header, "w-[190px] pl-3")}>Share</th>
            </tr>
          </thead>
          <tbody>
            {preview.filters.map((row, i) => {
              const share = row.removed === null ? 0 : row.removed / total(row);
              const delta = row.removed !== null && row.suggested_removed !== null ? row.removed - row.suggested_removed : 0;
              return (
                <tr key={row.key} className={cn("border-b border-line", row.differs && "bg-accent-tint")}>
                  <td className={cn("py-[13px] pr-3 font-bold", row.differs ? "pl-[10px] text-accent-text" : "text-ink")}>
                    {i + 1} · {row.reason}
                  </td>
                  <td className="px-3 py-[13px] text-body">{row.what}</td>
                  <td className={cn("px-3 py-[13px] text-right font-mono", row.differs && "font-bold text-accent-text")}>
                    {setting(row, row.setting)}
                  </td>
                  <td className="px-3 py-[13px] text-right font-mono text-muted">{setting(row, row.suggested_setting)}</td>
                  <td className={cn("px-3 py-[13px] text-right font-mono", delta !== 0 ? "font-bold text-accent-text" : "text-ink")}>
                    {row.removed === null ? "not run" : `${fmt(row.removed)} ${row.what}`}
                    {delta !== 0 && <span className="font-normal"> ({delta > 0 ? "+" : "−"}{fmt(Math.abs(delta))})</span>}
                  </td>
                  <td className="py-[13px] pl-3">
                    <div className="flex items-center gap-2">
                      <div className="h-1 flex-auto bg-line">
                        <div
                          className={cn("h-1", row.differs ? "bg-accent" : "bg-[#3f3f46]")}
                          style={{ width: `${(100 * share).toFixed(2)}%`, minWidth: row.removed ? 2 : 0 }}
                        />
                      </div>
                      <span className="w-[44px] text-right font-mono text-[13px] text-muted">{formatShare(share)}</span>
                    </div>
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
      <p className="m-0 text-[13px] text-muted">
        Filters run in order; teal marks a difference from the suggestion. Share: of all{" "}
        {fmt(preview.n_cells_before)} cells, or all {fmt(preview.n_genes_before)} genes.
      </p>
    </section>
  );
}

function MetricsSection({ qc }: { qc: QcPayload }) {
  const cells = qc.cells!;
  const violins = useMemo(() => {
    const metrics: [string, number[] | null][] = [
      ["Genes per cell", cells.n_genes],
      ["Counts per cell", cells.total_counts],
      ["% mitochondrial", qc.mito?.available ? cells.pct_mt : null],
    ];
    return metrics.map(([title, values]) => ({
      title,
      data: values
        ? ([
            {
              type: "violin",
              y: values,
              text: cells.barcode,
              name: title,
              points: "all",
              jitter: 0.45,
              pointpos: 0,
              marker: { size: 2.5, color: INK, opacity: 0.35 },
              line: { color: "#52525b", width: 1 },
              fillcolor: "rgba(17,17,17,0.06)",
              box: { visible: true, width: 0.12 },
              hoveron: "points",
              hovertemplate: `%{text}<br>${title} %{y}<extra></extra>`,
            },
          ] as Data[])
        : null,
    }));
  }, [cells, qc.mito]);

  return (
    <section aria-label="QC distributions" className="flex flex-col gap-3">
      <p className="m-0 text-[14px] text-muted">
        Calculated once for every cell, before any filtering (the same values the filters use).
      </p>
      <div className="grid grid-cols-[repeat(auto-fit,minmax(240px,1fr))] gap-4">
        {violins.map((v) =>
          v.data ? (
            <PlotlyChart
              key={v.title}
              data={v.data}
              height={300}
              label={`Violin plot of ${v.title} for every cell`}
              download={{ name: `Quality Control distribution ${v.title}`, title: `${v.title}, every cell before filtering` }}
              layout={{ title: { text: v.title, font: { size: 14 } }, xaxis: { visible: false }, yaxis: { gridcolor: "#ececef" } }}
            />
          ) : (
            <div key={v.title} className="flex h-[300px] items-center justify-center rounded-[4px] bg-chip text-center text-[14px] text-muted">
              {v.title}. Unavailable: no mitochondrial genes found.
            </div>
          ),
        )}
      </div>
    </section>
  );
}

/** The per-cell QC statistics (before filtering), next to the other numbers. */
function StatsSection({ qc }: { qc: QcPayload }) {
  return (
    <section aria-label="QC statistics per cell" className="flex flex-col gap-2">
      <h2 className="m-0 text-[16px] font-bold">QC statistics per cell</h2>
      <p className="m-0 text-[13px] text-muted">Measured once for every cell, before filtering.</p>
      <StatsTable qc={qc} />
    </section>
  );
}

/** Fixed to the bottom of the screen: what the settings would keep, Apply and Continue
 * (the shared step bar, with a status built from the live filter preview). */
function QcActionBar(props: {
  preview: OkPreview | null;
  previewProblem: string | null;
  error: string | null;
  confirmation: string | null;
  needsRerun: boolean;
  runLabel: string;
  runDisabled: boolean;
  upToDate: boolean;
  willClear: string[];
  onRun: () => void;
  onContinue: () => void;
}) {
  const p = props.preview;
  let status: ReactNode = null;
  if (props.error) status = <span className="font-bold text-accent-text">{props.error}</span>;
  else if (props.previewProblem) status = <span className="font-bold text-accent-text">{props.previewProblem}</span>;
  else if (props.upToDate && props.confirmation) status = <span className="text-body">{props.confirmation}</span>;
  else if (p)
    status = (
      <span className="text-body">
        Would keep <b className="font-mono text-ink">{fmt(p.n_cells_kept)}</b> of {fmt(p.n_cells_before)} cells and{" "}
        <b className="font-mono text-ink">{fmt(p.n_genes_kept)}</b> of {fmt(p.n_genes_before)} genes
        {props.needsRerun && <span className="font-bold text-accent-text"> · Quality Control was re-run. Apply Filtering again.</span>}
      </span>
    );
  return (
    <ActionBar
      label="Apply Filtering"
      status={status}
      note={props.willClear.length > 0 ? `Re-applying will clear: ${props.willClear.join(", ")}` : undefined}
      runLabel={props.runLabel}
      runDisabled={props.runDisabled}
      upToDate={props.upToDate}
      nextTitle="Preprocessing"
      onRun={props.onRun}
      onContinue={props.onContinue}
    />
  );
}

function StatsTable({ qc }: { qc: QcPayload }) {
  if (!qc.stats) return null;
  const hasMito = qc.stats.some((row) => row.pct_counts_mt !== undefined);
  const cell = (row: (typeof qc.stats)[number], value: number | undefined) =>
    value === undefined ? "–" : row.statistic === "Cells" ? fmt(value) : value.toLocaleString("en-US", { maximumFractionDigits: 2 });
  return (
    <div className="flex flex-col gap-2">
    <div className="flex justify-end">
      <CsvDownload
        name="Quality Control statistics per cell"
        columns={[
          { key: "statistic", header: "Statistic" },
          { key: "n_genes_by_counts", header: "Genes per cell" },
          { key: "total_counts", header: "Counts per cell" },
          ...(hasMito ? [{ key: "pct_counts_mt", header: "% mitochondrial" }] : []),
        ]}
        rows={() => qc.stats!.map((row) => ({ ...row }))}
      />
    </div>
    <div className="overflow-x-auto">
      <table className="w-full border-collapse text-[14px] tabular-nums">
        <thead>
          <tr className="border-b border-ink text-left text-[12px] tracking-[0.04em] text-muted uppercase">
            <th className="py-2 pr-4 font-bold">Statistic</th>
            <th className="py-2 pr-4 text-right font-bold">Genes per cell</th>
            <th className="py-2 pr-4 text-right font-bold">Counts per cell</th>
            {hasMito && <th className="py-2 text-right font-bold">% mitochondrial</th>}
          </tr>
        </thead>
        <tbody>
          {qc.stats.map((row) => (
            <tr key={row.statistic} className="border-b border-line last:border-0">
              <td className="py-2 pr-4">{row.statistic}</td>
              <td className="py-2 pr-4 text-right font-mono">{cell(row, row.n_genes_by_counts)}</td>
              <td className="py-2 pr-4 text-right font-mono">{cell(row, row.total_counts)}</td>
              {hasMito && <td className="py-2 text-right font-mono">{cell(row, row.pct_counts_mt)}</td>}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
    </div>
  );
}

// --- Helpers ------------------------------------------------------------------------------

function toBoxes(s: FilterSettings): FilterBoxes {
  return {
    min_genes: String(s.min_genes),
    min_cells: String(s.min_cells),
    max_genes: s.max_n_genes_by_counts === null ? "" : String(s.max_n_genes_by_counts),
    max_genes_off: s.max_n_genes_by_counts === null,
    max_mt: s.max_pct_counts_mt === null ? "" : String(s.max_pct_counts_mt),
    max_mt_off: s.max_pct_counts_mt === null,
  };
}

/** True when two sets of filter boxes mean the same settings ("No limit" ignores the value). */
function sameFilter(a: FilterBoxes, b: FilterBoxes): boolean {
  const num = (x: string) => (x.trim() === "" ? NaN : Number(x));
  return (
    num(a.min_genes) === num(b.min_genes) &&
    num(a.min_cells) === num(b.min_cells) &&
    a.max_genes_off === b.max_genes_off &&
    (a.max_genes_off || num(a.max_genes) === num(b.max_genes)) &&
    a.max_mt_off === b.max_mt_off &&
    (a.max_mt_off || num(a.max_mt) === num(b.max_mt))
  );
}

/** A fixed vertical position for cell i between 0 and 1, bunched towards the middle, so the
 * dots spread out instead of sitting on one line. The same cell always lands in the same place. */
function spread(i: number): number {
  const a = ((i * 2654435761) % 4294967296) / 4294967296;
  const b = ((i * 40503 + 12345) % 65536) / 65536;
  return 0.5 + ((a + b) / 2 - 0.5) * 0.92; // stays within 0.04–0.96
}

function formatShare(share: number): string {
  if (share === 0) return "0%";
  return share < 0.01 ? `${(100 * share).toFixed(1)}%` : `${Math.round(100 * share)}%`;
}

function fmt(value: number | undefined | null): string {
  return value === undefined || value === null ? "–" : value.toLocaleString("en-US", { maximumFractionDigits: 2 });
}
