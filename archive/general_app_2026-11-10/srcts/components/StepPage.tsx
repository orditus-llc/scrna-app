// The layout every step page shares: a heading for screen readers, headline numbers (Quality Control),
// then settings on the left and results on the right.
import { ArrowRight } from "lucide-react";
import type { ReactNode } from "react";

import { cn } from "../lib/utils";
import type { Step } from "../types";
import { Button } from "./ui/button";

/** The page's heading, for screen readers only: the step's name, status and description are in
 * the sidebar (analyst, 2026-11-04: no title or description at the top of the pages). */
export function StepHeader({ step }: { step: Step | undefined }) {
  return <h1 className="sr-only">{step ? `Step ${step.number}: ${step.title}` : ""}</h1>;
}

export interface SummaryItem {
  label: string;
  value: string;
  note: string;
  differs?: boolean; // shown in teal: differs from the validated reference
}

/** Headline numbers in a grey bordered box (read-only results; used by Quality Control). */
export function SummaryStrip({ items }: { items: SummaryItem[] }) {
  return (
    <section aria-label="Summary" className="grid grid-cols-[repeat(auto-fit,minmax(190px,1fr))] rounded-[4px] border border-line bg-chip px-5">
      {items.map((item, i) => (
        <div key={item.label} className={cn("flex flex-col gap-2 py-[18px] pr-5 pb-5", i > 0 && "border-l border-field pl-5")}>
          <span className="text-[12px] font-bold tracking-[0.04em] text-muted uppercase">{item.label}</span>
          <span className="text-[30px] leading-none font-semibold tracking-[-0.01em] tabular-nums">{item.value}</span>
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
export function Prerequisite({ title, text, onOpen }: { title: string; text?: string | null; onOpen: () => void }) {
  return (
    <div className="flex flex-wrap items-center gap-3 rounded-[4px] border border-line bg-panel px-4 py-3 text-[14px]">
      <span>
        <b>{text ?? `Run ${title} first.`}</b> This step uses its results. You can still look at the settings here.
      </span>
      <Button variant="secondary" size="small" onClick={onOpen}>
        Open {title} <ArrowRight aria-hidden="true" size={14} />
      </Button>
    </div>
  );
}
