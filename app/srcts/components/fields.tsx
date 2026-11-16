// Settings rows and messages shared by the step pages.
import type { ReactNode } from "react";

import { cn } from "../lib/utils";
import { InfoTip } from "./InfoTip";

/** One setting in a ruled list: label (with an ⓘ for its explanation) on the left, the box
 * on the right. When the value differs from the value it is compared with (`reference`: the
 * default, or on Quality Control the suggestion), the box gets a teal border and
 * "Default: …" (or `compareLabel`) appears (null = no limit). `note` is shown under the
 * label for things that must stay visible (e.g. "unavailable"). */
export function RuledField(props: {
  id: string;
  label: string;
  help: string;
  value: string;
  onChange: (value: string) => void;
  step: number;
  min?: number;
  unit?: string;
  reference: number | null;
  compareLabel?: string; // what `reference` is called on screen (default "Default")
  off?: boolean; // "No limit" (or `offLabel`) ticked
  onOff?: (off: boolean) => void; // shows a "No limit" checkbox
  offLabel?: string; // the checkbox's text, e.g. "Median of your cells" (default "No limit")
  offText?: string; // how a null reference reads, e.g. "median of your cells" (default "no limit")
  note?: string;
}) {
  const current = props.off ? null : props.value.trim() === "" ? NaN : Number(props.value);
  const differs = current !== props.reference;
  const referenceText =
    props.reference === null ? (props.offText ?? "no limit") : `${props.reference.toLocaleString("en-US", { maximumFractionDigits: 6 })}${props.unit ?? ""}`;
  return (
    <div className="flex items-start justify-between gap-4 border-b border-line py-[14px]">
      <div className="flex min-w-0 flex-col gap-[3px]">
        <span className="text-[14px] leading-snug">
          <label htmlFor={props.id} className="font-bold">
            {props.label}
          </label>{" "}
          <InfoTip label={`About ${props.label}`}>{props.help}</InfoTip>
        </span>
        {props.note && <span className="text-[13px] leading-snug font-bold text-accent-text">{props.note}</span>}
        {props.onOff && (
          <label className="mt-1 flex items-center gap-2 text-[13px]">
            <input
              type="checkbox"
              checked={props.off}
              onChange={(e) => props.onOff?.(e.target.checked)}
              className="size-4 accent-ink"
            />
            {props.offLabel ?? "No limit"}
          </label>
        )}
        {differs && (
          <span className="text-[13px] font-bold text-accent-text">
            {props.compareLabel ?? "Default"}: {referenceText}
          </span>
        )}
      </div>
      <div className="flex flex-none items-center gap-[6px] pt-[2px]">
        <input
          id={props.id}
          type="number"
          min={props.min ?? 0}
          step={props.step}
          value={props.value}
          disabled={props.off}
          onChange={(e) => props.onChange(e.target.value)}
          className={cn(
            "w-[92px] rounded-none border bg-white px-[10px] py-2 text-right font-mono text-[14px] tabular-nums disabled:bg-chip disabled:text-off",
            differs ? "border-2 border-accent" : "border-field",
          )}
        />
        {props.unit && <span className="font-mono text-[14px] text-body">{props.unit}</span>}
      </div>
    </div>
  );
}

/** An on/off setting in a ruled list. */
export function RuledCheck(props: {
  id: string;
  label: string;
  help: string;
  checked: boolean;
  onChange: (checked: boolean) => void;
  reference: boolean;
}) {
  const differs = props.checked !== props.reference;
  return (
    <div className="flex items-start justify-between gap-4 border-b border-line py-[14px]">
      <div className="flex min-w-0 flex-col gap-[3px]">
        <span className="text-[14px] leading-snug">
          <label htmlFor={props.id} className="font-bold">
            {props.label}
          </label>{" "}
          <InfoTip label={`About ${props.label}`}>{props.help}</InfoTip>
        </span>
        {differs && <span className="text-[13px] font-bold text-accent-text">Default: {props.reference ? "on" : "off"}</span>}
      </div>
      <input
        id={props.id}
        type="checkbox"
        checked={props.checked}
        onChange={(e) => props.onChange(e.target.checked)}
        className={cn("mt-1 size-5 flex-none accent-ink", differs && "outline-2 outline-offset-2 outline-accent")}
      />
    </div>
  );
}

/** A short message with a colored status dot. */
/** Read-only information (computed results, fixed or confirmed settings) in a grey bordered
 * box, so it reads as information rather than something to change; settings stay on white. */
export function InfoBox({ children, label, className }: { children: ReactNode; label?: string; className?: string }) {
  return (
    <div role={label ? "group" : undefined} aria-label={label} className={cn("rounded-none border border-line bg-chip px-4 py-3", className)}>
      {children}
    </div>
  );
}

export function Note({ dot, children }: { dot: string; children: ReactNode }) {
  return (
    <p className="m-0 flex items-start gap-[9px] rounded-none bg-chip px-3 py-2 text-[14px] leading-snug text-body">
      <span aria-hidden="true" className={`mt-[6px] inline-block size-2 flex-none rounded-full ${dot}`} />
      <span>{children}</span>
    </p>
  );
}

/** A choice from a list in a ruled list (teal when it differs from the default). */
export function RuledSelect(props: {
  id: string;
  label: string;
  help: string;
  value: string;
  options: { value: string; label: string }[];
  onChange: (value: string) => void;
  reference?: string; // leave out for a display-only choice
  note?: string; // shown under the label, e.g. "Suggested: Mouse"
}) {
  const differs = props.reference !== undefined && props.value !== props.reference;
  const referenceLabel = props.options.find((o) => o.value === props.reference)?.label;
  return (
    <div className="flex flex-col gap-2 border-b border-line py-[14px]">
      <div className="flex min-w-0 flex-col gap-[3px]">
        <span className="text-[14px] leading-snug">
          <label htmlFor={props.id} className="font-bold">
            {props.label}
          </label>{" "}
          <InfoTip label={`About ${props.label}`}>{props.help}</InfoTip>
        </span>
        {props.note && <span className="text-[13px] leading-snug text-muted">{props.note}</span>}
        {differs && <span className="text-[13px] font-bold text-accent-text">Default: {referenceLabel}</span>}
      </div>
      <select
        id={props.id}
        value={props.value}
        onChange={(e) => props.onChange(e.target.value)}
        className={cn(
          "rounded-none border bg-white px-[10px] py-2 text-[14px]",
          differs ? "border-2 border-accent" : "border-field",
        )}
      >
        {props.options.map((o) => (
          <option key={o.value} value={o.value}>
            {o.label}
          </option>
        ))}
      </select>
    </div>
  );
}
