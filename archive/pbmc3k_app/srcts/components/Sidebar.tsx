import { Check, ClipboardList } from "lucide-react";

import { cn } from "../lib/utils";
import type { Page } from "../navigation";
import type { Step } from "../types";

const STATUS_TEXT = { done: "Up to date", needs_rerun: "Needs re-running", not_run: "Not run" } as const;

/** The step list on every page: the progress line (✓ and a solid line when done, teal ▲
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
  referenceSettings,
}: {
  steps: Step[];
  nextKey: string | null;
  page: Page;
  exploreAvailable: boolean; // clustering is done
  onNavigate: (page: Page) => void;
  onOpenRecord: () => void;
  runCount: number;
  referenceSettings: boolean; // false when any setting differs from the validated preset
}) {
  const done = steps.filter((s) => s.status === "done").length;
  return (
    <nav aria-label="Analysis" className="flex flex-col gap-5 px-5 py-6">
      {!referenceSettings && (
        <p className="m-0 flex items-start gap-2 rounded-[4px] border border-accent-line bg-accent-tint px-3 py-2 text-[12.5px] leading-snug font-bold text-accent-text">
          <span aria-hidden="true" className="text-[10px] leading-[18px]">▲</span>
          Settings differ from the validated reference
        </p>
      )}
      <div className="flex flex-col">
        <span className="mb-3 px-3 text-[12px] font-bold tracking-[0.04em] text-muted uppercase">
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
                      step.status === "done" ? "border-solid border-ink" : "border-dotted border-field",
                    )}
                  />
                )}
                <button
                  type="button"
                  aria-current={current ? "page" : undefined}
                  onClick={() => onNavigate(step.key as Page)}
                  className={cn(
                    "relative flex w-full items-start gap-3 rounded-[4px] px-3 py-[6px] text-left",
                    current ? "bg-chip" : "hover:bg-chip",
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
                    <span className={cn("text-[14px] leading-snug", current ? "font-extrabold text-ink" : "font-semibold text-body")}>
                      {step.title}
                    </span>
                    <span
                      className={cn(
                        "text-[12px]",
                        step.status === "needs_rerun" ? "font-bold text-accent-text" : isNext ? "font-bold text-ink" : "text-muted",
                      )}
                    >
                      {step.status === "needs_rerun" ? "▲ " : ""}
                      {isNext && step.status === "not_run" ? "Next" : STATUS_TEXT[step.status]}
                    </span>
                  </span>
                </button>
              </li>
            );
          })}
        </ol>
      </div>

      <div className="border-t border-line pt-4">
        <button
          type="button"
          aria-current={page === "explorer" ? "page" : undefined}
          disabled={!exploreAvailable}
          onClick={() => onNavigate("explorer")}
          className={cn(
            "flex w-full flex-col gap-[2px] rounded-[4px] px-3 py-2 text-left disabled:cursor-not-allowed",
            page === "explorer" ? "bg-chip" : "enabled:hover:bg-chip",
          )}
        >
          <span className={cn("text-[14.5px] font-bold", exploreAvailable ? "text-ink" : "text-off")}>Explore Results</span>
          <span className="text-[12px] text-muted">
            {exploreAvailable ? "Clusters, cell types and genes on one UMAP" : "Opens after Clustering"}
          </span>
        </button>
      </div>

      <div className="border-t border-line pt-4">
        <button
          type="button"
          onClick={onOpenRecord}
          className="flex w-full items-start gap-2 rounded-[4px] px-3 py-2 text-left hover:bg-chip"
        >
          <ClipboardList aria-hidden="true" size={16} className="mt-[2px] flex-none text-muted" />
          <span className="flex flex-col gap-[2px]">
            <span className="text-[14px] font-bold text-ink">Run record</span>
            <span className="text-[12px] text-muted">
              Settings and run log{runCount > 0 ? ` · ${runCount} ${runCount === 1 ? "run" : "runs"}` : ""}
            </span>
          </span>
        </button>
      </div>
    </nav>
  );
}
