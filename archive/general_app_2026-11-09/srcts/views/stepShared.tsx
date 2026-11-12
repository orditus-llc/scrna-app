// Pieces shared by the step pages: the run button with "& Continue", section titles,
// empty/loading states and small helpers.
import { useShinyInput } from "@posit-dev/shinyreact";
import { ArrowRight, FastForward, RotateCcw } from "lucide-react";
import { useContext, useEffect, useRef, useState, type ReactNode } from "react";

import { Note } from "../components/fields";
import { RunRemainingContext } from "../components/RunRemaining";
import { Button } from "../components/ui/button";
import type { Page, StepKey } from "../navigation";
import { clickCount, type Click, type Notice, type Step, type WorkflowPayload } from "../types";

export interface PageProps {
  workflow: WorkflowPayload;
  onNavigate: (page: Page) => void;
}

export function stepOf(workflow: WorkflowPayload, key: StepKey): Step | undefined {
  return workflow.steps.find((s) => s.key === key);
}

// --- Shared pieces --------------------------------------------------------------------------

/** The key of the first earlier step that isn't done (the one to open first). */
export function prerequisiteKey(workflow: WorkflowPayload, step: Step): Page {
  return (workflow.steps.find((s) => s.title === step.prerequisite)?.key ?? "data_overview") as Page;
}

/** A run button that sends the boxes with the click, shows "Running…" until the server
 * answers, and keeps the error that answered this click (if any). The page stays open, so
 * the new results can be looked at; "Continue" is a separate button. */
export function useRunButton<B>(
  inputId: string,
  step: string,
  answer: { status: string },
  notice: Notice | null,
) {
  const [click, setClick] = useShinyInput<Click<B>>(inputId, 0);
  const [busy, setBusy] = useState(false);
  const noticeAtClick = useRef<number | null>(null);
  const clickedOnce = useRef(false);

  // The server answers with a new payload for this page (even a re-run whose status stays
  // "done") or, if something went wrong, with a new notice.
  useEffect(() => setBusy(false), [answer, notice?.id]);

  const error =
    clickedOnce.current && notice?.step === step && notice.id !== noticeAtClick.current ? notice.text : null;

  return {
    busy,
    error,
    click(boxes: B) {
      noticeAtClick.current = notice?.id ?? null;
      clickedOnce.current = true;
      setBusy(true);
      setClick({ n: clickCount(click) + 1, settings: boxes });
    },
  };
}

/** Reset, the run button ("Run X & Continue", "Re-run …", or "Continue to …" when nothing
 * changed), and what re-running would clear. */
/** Under a page's settings: Reset, and a note when the settings were changed after the last
 * run. Running and continuing are in the bar at the bottom of the screen (`StepBar`). */
export function SettingsReset({
  onReset,
  changedSinceRun,
  label = "Reset to default settings",
}: {
  onReset: () => void;
  changedSinceRun: boolean;
  label?: string;
}) {
  return (
    <div className="flex flex-col items-start gap-3 pt-[18px]">
      <Button variant="secondary" size="small" onClick={onReset}>
        <RotateCcw aria-hidden="true" size={14} /> {label}
      </Button>
      {changedSinceRun && <Note dot="bg-accent">Settings changed since the last run. Re-run to apply them.</Note>}
    </div>
  );
}

/** The bar fixed to the bottom of every step page: where the step stands (left) and the Run
 * and Continue buttons (right). The button to press next is the dark one. Once the step is up
 * to date, "Run the remaining N steps" runs the rest in one go; while that runs, the bar shows
 * its progress and the other buttons wait. When it ends, the page stays where it is and the bar
 * says how it went, with "Open …" for the step that needs attention (or Explore Results). */
export function ActionBar(props: {
  label: string; // e.g. "Run Preprocessing" (read by screen readers as the bar's name)
  status: ReactNode;
  note?: ReactNode; // a second, quieter line (e.g. "Re-running will clear: …")
  runLabel?: string; // no Run button when left out (Data Overview loads data in the page itself)
  runDisabled?: boolean;
  upToDate: boolean;
  nextTitle: string | undefined;
  onRun?: () => void;
  onContinue: () => void;
}) {
  const rest = useContext(RunRemainingContext);
  const running = rest?.running ?? false;
  const outcome = rest?.outcome ?? null; // how "Run the remaining steps" ended, started on this page
  const offerRest = !running && !outcome && props.upToDate && (rest?.remaining ?? 0) > 0;
  let status: ReactNode = props.status;
  if (running) status = <span role="status" className="font-bold text-ink">{rest?.progress}</span>;
  else if (outcome)
    status = (
      <span role={outcome.kind === "error" ? "alert" : "status"} className={outcome.kind === "error" ? "font-bold text-accent-text" : "text-ink"}>
        {outcome.text}
      </span>
    );
  return (
    <div
      role="region"
      aria-label={props.label}
      className="sticky bottom-0 z-40 -mx-8 -mb-10 flex flex-wrap items-center justify-between gap-x-6 gap-y-2 border-t border-ink bg-white/95 px-8 py-3 backdrop-blur-sm"
    >
      <div className="flex min-w-0 flex-[1_1_320px] flex-col gap-[2px] text-[14px]">
        {status}
        {!running && !outcome && props.note && <span className="text-[13px] text-muted">{props.note}</span>}
      </div>
      <div className="flex flex-wrap gap-[10px]">
        {props.runLabel && (
          <Button variant={props.upToDate ? "secondary" : "primary"} disabled={props.runDisabled || running} onClick={props.onRun}>
            {props.runLabel}
          </Button>
        )}
        {offerRest && rest && (
          <Button variant="secondary" onClick={rest.start} title="Runs every later step with the current settings">
            <FastForward aria-hidden="true" size={15} />
            Run the remaining {rest.remaining} {rest.remaining === 1 ? "step" : "steps"}
          </Button>
        )}
        {outcome && (
          <Button variant="primary" onClick={() => (window.location.hash = outcome.open.page)}>
            Open {outcome.open.title} <ArrowRight aria-hidden="true" size={16} />
          </Button>
        )}
        {props.nextTitle && (
          <Button
            variant={props.upToDate && !outcome ? "primary" : "secondary"}
            disabled={!props.upToDate || running}
            title={props.upToDate ? undefined : "Run this step first"}
            onClick={props.onContinue}
          >
            Continue to {props.nextTitle} <ArrowRight aria-hidden="true" size={16} />
          </Button>
        )}
      </div>
    </div>
  );
}

/** The bottom bar for a step that runs one stage (Preprocessing, Dimensionality Reduction, Clustering, …). */
export function StepBar(props: {
  label: string; // the step, as in "Run Preprocessing"
  nextTitle: string | undefined;
  status: string;
  blocked: boolean;
  prerequisite: string | null;
  prerequisiteText?: string | null; // why, when it's the species/tissue confirmation rather than a run
  confirmation: string | null;
  busy: boolean;
  busyText: string;
  error: string | null;
  changedSinceRun: boolean;
  willClear: string[];
  missing?: string | null; // a choice still to make before the step can run (Run stays off)
  onRun: () => void;
  onContinue: () => void;
}) {
  const upToDate = props.status === "done" && !props.changedSinceRun;
  let status: ReactNode;
  if (props.error) status = <span className="font-bold text-accent-text">{props.error}</span>;
  else if (props.busy) status = <span className="text-body">{props.busyText}</span>;
  else if (props.prerequisite) status = <span className="text-body">{props.prerequisiteText ?? `Run ${props.prerequisite} first.`}</span>;
  else if (props.missing) status = <span className="font-bold text-accent-text">{props.missing}</span>;
  else if (props.changedSinceRun) status = <span className="font-bold text-accent-text">Settings changed since the last run. Re-run to apply them.</span>;
  else if (props.status === "done") status = <span className="text-body">{props.confirmation ?? "Up to date."}</span>;
  else if (props.status === "needs_rerun") status = <span className="font-bold text-accent-text">An earlier step was redone, so run this step again.</span>;
  else status = <span className="text-body">Not run yet. Review the settings, then run it.</span>;
  return (
    <ActionBar
      label={`Run ${props.label}`}
      status={status}
      note={!upToDate && !props.busy && props.willClear.length > 0 ? `Re-running will clear: ${props.willClear.join(", ")}` : undefined}
      runLabel={props.busy ? props.busyText : `${props.status === "not_run" ? "Run" : "Re-run"} ${props.label}`}
      runDisabled={props.busy || props.blocked || Boolean(props.missing)}
      upToDate={upToDate}
      nextTitle={props.nextTitle}
      onRun={props.onRun}
      onContinue={props.onContinue}
    />
  );
}

export function SectionTitle({ children }: { children: ReactNode }) {
  return <h2 className="m-0 border-b border-field pb-[10px] text-[16px] font-bold">{children}</h2>;
}

export function SubTitle({ children }: { children: ReactNode }) {
  return <h3 className="mt-5 mb-0 text-[12px] font-bold tracking-[0.04em] text-muted uppercase">{children}</h3>;
}

export function Empty({ children }: { children: ReactNode }) {
  return <p className="m-0 rounded-[4px] bg-panel px-4 py-10 text-center text-[14px] text-muted">{children}</p>;
}

export function Loading() {
  return <p className="m-0 px-8 py-10 text-muted">Loading…</p>;
}

/** True when two sets of boxes hold the same values (numbers compared as numbers). */
export function sameAs<B extends object>(a: B, b: B): boolean {
  return (Object.keys(a) as (keyof B)[]).every((k) =>
    // Identical text (including two empty boxes, e.g. a number replaced by a ticked option)
    // or the same number counts as unchanged.
    typeof a[k] === "boolean" ? a[k] === b[k] : a[k] === b[k] || (Number(a[k]) === Number(b[k]) && String(a[k]).trim() !== ""),
  );
}

export function fmt(value: number): string {
  return value.toLocaleString("en-US", { maximumFractionDigits: 2 });
}
