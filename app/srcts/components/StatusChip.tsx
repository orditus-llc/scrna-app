import type { StepStatus } from "../types";

const LOOK: Record<StepStatus, { dot: string; text: string }> = {
  done: { dot: "bg-ok", text: "Up to date" },
  needs_rerun: { dot: "bg-accent", text: "Needs re-running" },
  not_run: { dot: "bg-off", text: "Not run" },
};

/** Grey chip with a status mark: readable by text, by color and by shape (▲ = needs re-running). */
export function StatusChip({ status }: { status: StepStatus }) {
  const look = LOOK[status];
  return (
    <span className="inline-flex items-center gap-[7px] rounded-none bg-chip px-[11px] py-1 text-[13px] font-bold text-[#27272a]">
      {status === "needs_rerun" ? (
        <span aria-hidden="true" className="text-[10px] leading-none text-accent">▲</span>
      ) : (
        <span aria-hidden="true" className={`inline-block size-2 rounded-full ${look.dot}`} />
      )}
      {look.text}
    </span>
  );
}
