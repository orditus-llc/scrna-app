import { Check, ClipboardList } from "lucide-react";

import { STEP_DESCRIPTIONS } from "../lib/stepDescriptions";
import { cn } from "../lib/utils";
import type { Page } from "../navigation";
import type { Step } from "../types";
import { InfoTip } from "./InfoTip";

const STATUS_TEXT = { done: "Up to date", needs_rerun: "Needs re-running", not_run: "Not run" } as const;

/** The step list on every page (on a very light grey background): the progress line (✓ and a solid line when done, teal ▲
 * when it needs re-running, hollow and dotted when not run yet), each step a link to its page;
 * then Explore Results and the Run record (settings and run log, opened in a side panel). */
export function Sidebar({
  steps,
  nextKey,
  page,
  exploreAvailable,
  onNavigate,
  onOpenRecord,
  runCount,
}: {
  steps: Step[];
  nextKey: string | null;
  page: Page;
  exploreAvailable: boolean; // clustering is done
  onNavigate: (page: Page) => void;
  onOpenRecord: () => void;
  runCount: number;
}) {
  const done = steps.filter((s) => s.status === "done").length;
  return (
    <nav aria-label="Analysis" className="flex flex-col gap-5 px-5 py-6">
      <div className="flex flex-col">
        <span className="mb-3 px-3 text-[12px] font-bold tracking-[0.04em] text-side-muted uppercase">
          Steps · {done} of {steps.length} done
        </span>
        <ol className="m-0 flex list-none flex-col p-0">
          {steps.map((step, i) => {
            const current = page === step.key;
            const isNext = step.key === nextKey;
            const last = i === steps.length - 1;
            return (
              <li key={step.key} className="relative">
                {!last && (
                  <span
                    aria-hidden="true"
                    className={cn(
                      "absolute top-[36px] -bottom-[6px] left-[25px] w-0 border-l-2",
                      step.status === "done" ? "border-solid border-ink" : "border-dotted border-side-line",
                    )}
                  />
                )}
                <button
                  type="button"
                  aria-current={current ? "page" : undefined}
                  disabled={step.locked}
                  title={step.locked ? (step.prerequisite_text ?? "Load your data on Data Overview first.") : undefined}
                  onClick={() => onNavigate(step.key as Page)}
                  className={cn(
                    "relative flex w-full items-start gap-3 rounded-[4px] py-[6px] pr-7 pl-3 text-left disabled:cursor-not-allowed",
                    step.locked && "opacity-45",
                    current ? "bg-side-raised" : "enabled:hover:bg-side-raised",
                  )}
                >
                  <span
                    aria-hidden="true"
                    className={cn(
                      "mt-[2px] flex size-7 flex-none items-center justify-center rounded-[4px] text-[13px] font-extrabold",
                      step.status === "done" && "bg-ink text-white",
                      step.status === "needs_rerun" && "bg-accent text-white",
                      step.status === "not_run" && (isNext ? "bg-ink text-white" : "border-2 border-field bg-white text-muted"),
                    )}
                  >
                    {step.status === "done" ? <Check size={15} strokeWidth={3} /> : step.number}
                  </span>
                  <span className="flex min-w-0 flex-col gap-[1px] pb-3">
                    <span className={cn("text-[14px] leading-snug", current ? "font-extrabold text-side-text" : "font-semibold text-side-body")}>
                      {step.title}
                    </span>
                    <span
                      className={cn(
                        "text-[12px]",
                        step.status === "needs_rerun" ? "font-bold text-side-accent" : isNext ? "font-bold text-side-text" : "text-side-muted",
                      )}
                    >
                      {step.status === "needs_rerun" ? "▲ " : ""}
                      {step.waiting
                        ? step.waiting
                        : isNext && step.status === "not_run"
                          ? "Next"
                          : STATUS_TEXT[step.status]}
                    </span>
                  </span>
                </button>
                {/* What the step does (outside the link: a button can't hold another button). */}
                <span className="absolute top-[8px] right-1">
                  <InfoTip label={`About ${step.title}`} onDark>
                    {STEP_DESCRIPTIONS[step.key]}
                  </InfoTip>
                </span>
              </li>
            );
          })}
        </ol>
      </div>

      <div className="relative border-t border-side-line pt-4">
        <span className="absolute top-[24px] right-1">
          <InfoTip label="About Explore Results" onDark>
            {STEP_DESCRIPTIONS.explorer}
          </InfoTip>
        </span>
        <button
          type="button"
          aria-current={page === "explorer" ? "page" : undefined}
          disabled={!exploreAvailable}
          onClick={() => onNavigate("explorer")}
          className={cn(
            "flex w-full flex-col gap-[2px] rounded-[4px] px-3 py-2 text-left disabled:cursor-not-allowed",
            page === "explorer" ? "bg-side-raised" : "enabled:hover:bg-side-raised",
          )}
        >
          <span className={cn("text-[14px] font-bold", exploreAvailable ? "text-side-text" : "text-side-off")}>Explore Results</span>
          <span className="text-[12px] text-side-muted">
            {exploreAvailable ? "Clusters, cell types and genes on one UMAP" : "Opens after Clustering"}
          </span>
        </button>
      </div>

      <div className="border-t border-side-line pt-4">
        <button
          type="button"
          onClick={onOpenRecord}
          className="flex w-full items-start gap-2 rounded-[4px] px-3 py-2 text-left hover:bg-side-raised"
        >
          <ClipboardList aria-hidden="true" size={16} className="mt-[2px] flex-none text-side-muted" />
          <span className="flex flex-col gap-[2px]">
            <span className="text-[14px] font-bold text-side-text">Run record</span>
            <span className="text-[12px] text-side-muted">
              Settings and run log{runCount > 0 ? ` · ${runCount} ${runCount === 1 ? "run" : "runs"}` : ""}
            </span>
          </span>
        </button>
      </div>
    </nav>
  );
}
