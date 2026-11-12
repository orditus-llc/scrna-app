// Settings rows and messages shared by the step pages.
import type { ReactNode } from "react";

import { cn } from "../lib/utils";
import { InfoTip } from "./InfoTip";

/** One setting in a ruled list: label (with an ⓘ for its explanation) on the left, the box
 * on the right. When the value differs from the validated `reference`, the box gets a teal
 * border and "Reference value: …" appears (null reference = no limit). `note` is shown
 * under the label for things that must stay visible (e.g. "unavailable"). */
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
  off?: boolean; // "No limit" ticked
  onOff?: (off: boolean) => void; // shows a "No limit" checkbox
  note?: string;
}) {
  const current = props.off ? null : props.value.trim() === "" ? NaN : Number(props.value);
  const differs = current !== props.reference;
  const referenceText =
    props.reference === null ? "no limit" : `${props.reference.toLocaleString("en-US", { maximumFractionDigits: 6 })}${props.unit ?? ""}`;
  return (
    <div className="flex items-start justify-between gap-4 border-b border-line py-[14px]">
      <div className="flex min-w-0 flex-col gap-[3px]">
        <span className="text-[14px] leading-snug">
          <label htmlFor={props.id} className="font-bold">
            {props.label}
          </label>{" "}
          <InfoTip label={`About ${props.label}`}>{props.help}</InfoTip>
        </span>
        {props.note && <span className="text-[12.5px] leading-snug font-bold text-accent-text">{props.note}</span>}
        {props.onOff && (
          <label className="mt-1 flex items-center gap-2 text-[13px]">
            <input
              type="checkbox"
              checked={props.off}
              onChange={(e) => props.onOff?.(e.target.checked)}
              className="size-4 accent-ink"
            />
            No limit
          </label>
        )}
        {differs && <span className="text-[12.5px] font-bold text-accent-text">Reference value: {referenceText}</span>}
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
            "w-[92px] rounded-[4px] border bg-white px-[10px] py-2 text-right font-mono text-[15px] tabular-nums disabled:bg-chip disabled:text-off",
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
        {differs && <span className="text-[12.5px] font-bold text-accent-text">Reference value: {props.reference ? "on" : "off"}</span>}
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
export function Note({ dot, children }: { dot: string; children: ReactNode }) {
  return (
    <p className="m-0 flex items-start gap-[9px] rounded-[4px] bg-chip px-3 py-2 text-[13.5px] leading-snug text-body">
      <span aria-hidden="true" className={`mt-[6px] inline-block size-2 flex-none rounded-full ${dot}`} />
      <span>{children}</span>
    </p>
  );
}

/** A choice from a list in a ruled list (teal when it differs from the reference). */
export function RuledSelect(props: {
  id: string;
  label: string;
  help: string;
  value: string;
  options: { value: string; label: string }[];
  onChange: (value: string) => void;
  reference?: string; // leave out for a display-only choice
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
        {differs && <span className="text-[12.5px] font-bold text-accent-text">Reference value: {referenceLabel}</span>}
      </div>
      <select
        id={props.id}
        value={props.value}
        onChange={(e) => props.onChange(e.target.value)}
        className={cn(
          "rounded-[4px] border bg-white px-[10px] py-2 text-[14px]",
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
