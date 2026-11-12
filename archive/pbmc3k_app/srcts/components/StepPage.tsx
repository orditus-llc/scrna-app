// The layout every step page shares: title and explanation, a strip of headline numbers,
// then settings on the left and results on the right.
import { ArrowRight } from "lucide-react";
import type { ReactNode } from "react";

import { cn } from "../lib/utils";
import type { Step } from "../types";
import { StatusChip } from "./StatusChip";
import { Button } from "./ui/button";

export function StepHeader({ step, children }: { step: Step | undefined; children: ReactNode }) {
  return (
    <div className="flex flex-col gap-[6px]">
      <div className="flex flex-wrap items-center gap-x-4 gap-y-2">
        {step && <span className="text-[13px] font-bold text-muted">Step {step.number}</span>}
        <h1 className="m-0 text-[34px] font-extrabold tracking-[-0.02em]">{step?.title}</h1>
        {step && <StatusChip status={step.status} />}
      </div>
      <p className="m-0 max-w-[780px] text-[15px] leading-normal text-muted">{children}</p>
    </div>
  );
}

export interface SummaryItem {
  label: string;
  value: string;
  note: string;
  differs?: boolean; // shown in teal: differs from the validated reference
}

/** Headline numbers between two rules. */
export function SummaryStrip({ items }: { items: SummaryItem[] }) {
  return (
    <section aria-label="Summary" className="grid grid-cols-[repeat(auto-fit,minmax(190px,1fr))] border-t-2 border-b border-ink">
      {items.map((item, i) => (
        <div key={item.label} className={cn("flex flex-col gap-2 py-[18px] pr-5 pb-5", i > 0 && "border-l border-line pl-5")}>
          <span className="text-[12px] font-bold tracking-[0.04em] text-muted uppercase">{item.label}</span>
          <span className="text-[34px] leading-none font-extrabold tracking-[-0.02em] tabular-nums">{item.value}</span>
          <span className={cn("text-[13px] tabular-nums", item.differs ? "font-bold text-accent-text" : "text-muted")}>{item.note}</span>
        </div>
      ))}
    </section>
  );
}

/** Settings (left) and results (right); they stack on narrow screens. */
export function SettingsAndResults({ settings, results }: { settings: ReactNode; results: ReactNode }) {
  return (
    <div className="flex flex-wrap gap-10">
      <section aria-label="Settings" className="flex min-w-0 flex-[1_1_300px] flex-col sm:max-w-[380px]">
        {settings}
      </section>
      <section aria-label="Results" className="flex min-w-0 flex-[999_1_560px] flex-col gap-3">
        {results}
      </section>
    </div>
  );
}

/** Shown on a step whose earlier step isn't done yet. */
export function Prerequisite({ title, onOpen }: { title: string; onOpen: () => void }) {
  return (
    <div className="flex flex-wrap items-center gap-3 rounded-[4px] border border-line bg-panel px-4 py-3 text-[14px]">
      <span>
        <b>Run {title} first.</b> This step uses its results. You can still look at the settings here.
      </span>
      <Button variant="secondary" size="small" onClick={onOpen}>
        Open {title} <ArrowRight aria-hidden="true" size={14} />
      </Button>
    </div>
  );
}
