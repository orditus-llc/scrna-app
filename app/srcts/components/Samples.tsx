// Several samples on Data Overview (Milestone 6b). Shown only when they're needed, so loading a
// single file works exactly as before:
// - SampleTable: several files were chosen; name each sample (and optionally its condition), then load.
// - SampleColumnChoice: one file was loaded that may contain several samples; pick the column.
// - SampleList: the loaded samples (read-only, inside the Loaded data box).
import { useShinyInput } from "@posit-dev/shinyreact";
import { X } from "lucide-react";
import { useState } from "react";

import { clickCount, type Click, type SamplesInfo } from "../types";
import { InfoBox } from "./fields";
import { InfoTip } from "./InfoTip";
import { Button } from "./ui/button";

export interface PendingSample {
  file: File;
  name: string;
  condition: string;
}

/** A sample's starting name: the file name without its ending ("donor_1.h5ad" -> "donor_1"). */
export function sampleNameFromFile(fileName: string): string {
  return fileName.replace(/\.(h5ad|h5)$/i, "");
}

/** What's wrong with the typed names, if anything (the server checks again). */
export function sampleTableProblem(samples: PendingSample[], maxSamples: number): string | null {
  if (samples.length < 2) return "Choose at least two files for several samples.";
  if (samples.length > maxSamples) return `The app handles up to ${maxSamples} samples for now.`;
  const names = samples.map((s) => s.name.trim());
  if (names.some((n) => !n)) return "Every sample needs a name.";
  const repeated = names.find((n, i) => names.indexOf(n) !== i);
  if (repeated) return `Sample names must be different ("${repeated}" is used twice).`;
  return null;
}

/** Name each chosen file's sample (and optionally its condition), then load them together. */
export function SampleTable({
  samples,
  maxSamples,
  busy,
  onChange,
  onLoad,
  onCancel,
}: {
  samples: PendingSample[];
  maxSamples: number;
  busy: boolean;
  onChange: (samples: PendingSample[]) => void;
  onLoad: () => void;
  onCancel: () => void;
}) {
  const problem = sampleTableProblem(samples, maxSamples);
  const set = (i: number, change: Partial<PendingSample>) => onChange(samples.map((s, j) => (j === i ? { ...s, ...change } : s)));
  return (
    <section aria-label="Name your samples" className="flex flex-col gap-3">
      <h2 className="m-0 border-b border-field pb-[10px] text-[16px] font-bold">
        Name your samples{" "}
        <InfoTip label="About several samples">
          Each file is one sample (for example one donor or one condition). They are checked like a single file,
          then combined; only genes measured in every sample are kept.
        </InfoTip>
      </h2>
      <div className="max-w-[760px] overflow-x-auto">
        <table className="w-full border-collapse text-[14px]">
          <thead>
            <tr className="border-b border-line text-left text-[12px] font-bold tracking-[0.04em] text-muted uppercase">
              <th className="py-2 pr-3 font-bold">File</th>
              <th className="py-2 pr-3 font-bold">Sample name</th>
              <th className="py-2 pr-3 font-bold">
                Condition <span className="font-normal normal-case">(optional)</span>
              </th>
              <th className="py-2" />
            </tr>
          </thead>
          <tbody>
            {samples.map((s, i) => (
              <tr key={`${s.file.name}-${i}`} className="border-b border-line">
                <td className="py-2 pr-3 font-mono text-[13px] text-body">{s.file.name}</td>
                <td className="py-2 pr-3">
                  <input
                    aria-label={`Sample name for ${s.file.name}`}
                    value={s.name}
                    maxLength={40}
                    onChange={(e) => set(i, { name: e.target.value })}
                    className="w-full min-w-[140px] rounded-none border border-field bg-white px-2 py-[6px] text-[14px]"
                  />
                </td>
                <td className="py-2 pr-3">
                  <input
                    aria-label={`Condition for ${s.file.name}`}
                    value={s.condition}
                    maxLength={40}
                    placeholder="e.g. control"
                    onChange={(e) => set(i, { condition: e.target.value })}
                    className="w-full min-w-[120px] rounded-none border border-field bg-white px-2 py-[6px] text-[14px]"
                  />
                </td>
                <td className="py-2 text-right">
                  <button
                    type="button"
                    aria-label={`Remove ${s.file.name}`}
                    onClick={() => onChange(samples.filter((_, j) => j !== i))}
                    className="inline-flex size-7 items-center justify-center rounded-none text-muted hover:bg-chip hover:text-ink"
                  >
                    <X size={15} aria-hidden="true" />
                  </button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div className="flex flex-wrap items-center gap-3">
        <Button disabled={problem !== null || busy} onClick={onLoad}>
          Load {samples.length} samples
        </Button>
        <Button variant="secondary" disabled={busy} onClick={onCancel}>
          Cancel
        </Button>
        {problem && <span className="text-[13px] text-muted">{problem}</span>}
      </div>
    </section>
  );
}

/** One loaded file that may contain several samples: choose the column naming them (and optionally
 * the one naming conditions), or keep it as one dataset. Opens by itself only when a column looks
 * like a sample column; otherwise it is a quiet link, so single-sample users aren't asked. */
export function SampleColumnChoice({ info }: { info: SamplesInfo }) {
  const [click, setClick] = useShinyInput<Click<{ sample_column: string; condition_column: string }>>("define_samples", 0);
  const likely = info.choices.find((c) => c.hint === "sample") ?? null;
  const [open, setOpen] = useState(false);
  const [sampleColumn, setSampleColumn] = useState(info.sample_column ?? likely?.column ?? "");
  const [conditionColumn, setConditionColumn] = useState(
    info.condition_column ?? info.choices.find((c) => c.hint === "condition")?.column ?? "",
  );
  const send = (sample: string, condition: string) =>
    setClick({ n: clickCount(click) + 1, settings: { sample_column: sample, condition_column: condition } });

  if (info.choices.length === 0 && !info.sample_column) return null;

  if (info.sample_column && !open) {
    return (
      <InfoBox label="Samples in this file" className="flex flex-wrap items-center gap-x-4 gap-y-2 text-[14px]">
        <span className="min-w-0 flex-[1_1_360px]">
          <b>{info.list.length} samples</b> from the column <span className="font-mono">{info.sample_column}</span>
          {info.condition_column && (
            <>
              {" "}
              · conditions from <span className="font-mono">{info.condition_column}</span>
            </>
          )}
        </span>
        <Button variant="secondary" size="small" onClick={() => setOpen(true)}>
          Change
        </Button>
        <Button variant="secondary" size="small" onClick={() => send("", "")}>
          Treat as one dataset
        </Button>
      </InfoBox>
    );
  }

  if (!likely && !open) {
    return (
      <button
        type="button"
        className="self-start text-[13px] font-bold text-ink underline underline-offset-2"
        onClick={() => setOpen(true)}
      >
        Does this file contain several samples?
      </button>
    );
  }

  const label = (c: SamplesInfo["choices"][number]) =>
    `${c.column} · ${c.n_values} values (${c.values.slice(0, 3).join(", ")}${c.n_values > 3 ? ", …" : ""})`;
  return (
    <section aria-label="Samples in this file" className="flex max-w-[760px] flex-col gap-3">
      <p className="m-0 text-[14px] text-body">
        {likely && !info.sample_column ? (
          <>
            This file may contain several samples (column <span className="font-mono">{likely.column}</span>). Use them as
            samples, or keep it as one dataset.
          </>
        ) : (
          "Choose the column that names each cell's sample."
        )}
      </p>
      <div className="flex flex-wrap gap-4">
        <label className="flex min-w-[260px] flex-1 flex-col gap-[6px] text-[14px] font-bold">
          Sample column
          <select
            value={sampleColumn}
            onChange={(e) => setSampleColumn(e.target.value)}
            className="rounded-none border border-field bg-white px-[10px] py-2 text-[14px] font-normal"
          >
            <option value="">Choose a column</option>
            {info.choices.map((c) => (
              <option key={c.column} value={c.column}>
                {label(c)}
              </option>
            ))}
          </select>
        </label>
        <label className="flex min-w-[260px] flex-1 flex-col gap-[6px] text-[14px] font-bold">
          <span>
            Condition column <span className="font-normal text-muted">(optional)</span>
          </span>
          <select
            value={conditionColumn}
            onChange={(e) => setConditionColumn(e.target.value)}
            className="rounded-none border border-field bg-white px-[10px] py-2 text-[14px] font-normal"
          >
            <option value="">None</option>
            {info.choices
              .filter((c) => c.column !== sampleColumn)
              .map((c) => (
                <option key={c.column} value={c.column}>
                  {label(c)}
                </option>
              ))}
          </select>
        </label>
      </div>
      <div className="flex flex-wrap items-center gap-3">
        <Button disabled={!sampleColumn} onClick={() => (send(sampleColumn, conditionColumn), setOpen(false))}>
          Use as samples
        </Button>
        <Button
          variant="secondary"
          onClick={() => {
            if (info.sample_column) send("", "");
            setOpen(false);
            setSampleColumn("");
          }}
        >
          Keep as one dataset
        </Button>
      </div>
    </section>
  );
}

/** The loaded samples, for the Loaded data box. */
export function SampleList({ info }: { info: SamplesInfo }) {
  if (info.list.length === 0) return null;
  const anyCondition = info.list.some((s) => s.condition);
  return (
    <table className="w-full max-w-[720px] border-collapse text-[14px]">
      <thead>
        <tr className="text-left text-[12px] tracking-[0.04em] text-muted uppercase">
          <th className="py-1 pr-4 font-bold">Sample</th>
          {anyCondition && <th className="py-1 pr-4 font-bold">Condition</th>}
          <th className="py-1 pr-4 text-right font-bold">Cells</th>
          {info.from_files && <th className="py-1 font-bold">File</th>}
        </tr>
      </thead>
      <tbody>
        {info.list.map((s) => (
          <tr key={s.name} className="border-t border-field">
            <td className="py-1 pr-4 font-bold">{s.name}</td>
            {anyCondition && <td className="py-1 pr-4">{s.condition ?? "–"}</td>}
            <td className="py-1 pr-4 text-right tabular-nums">{s.n_cells.toLocaleString("en-US")}</td>
            {info.from_files && <td className="py-1 font-mono text-[13px] text-body">{s.file}</td>}
          </tr>
        ))}
      </tbody>
    </table>
  );
}
