import { ChevronDown } from "lucide-react";
import { useState } from "react";

import { cn } from "../lib/utils";
import type { LogEntry, SettingItem } from "../types";

/** The settings that most affect results, plus what is saved with every run. */
export function SettingsRecord({ items, versions }: { items: SettingItem[]; versions: Record<string, string> }) {
  const versionText = Object.entries(versions)
    .map(([name, version]) => `${name} ${version}`)
    .join(", ");
  return (
    <section className="flex flex-col gap-3 rounded-[4px] border border-line bg-panel p-5">
      <h2 className="m-0 text-[16px] font-bold">Settings record</h2>
      <dl className="m-0 grid grid-cols-[max-content_minmax(0,1fr)] gap-x-4 gap-y-2 text-[14px]">
        {items.map((item) => (
          <div key={item.label} className="contents">
            <dt className="whitespace-nowrap text-muted">{item.label}</dt>
            <dd className="m-0 font-mono">{item.value}</dd>
          </div>
        ))}
      </dl>
      <p className="m-0 border-t border-line pt-[10px] text-[13px] leading-snug text-muted">
        Saved with each run, with seeds and package versions ({versionText}).
      </p>
    </section>
  );
}

/** One entry per stage run, oldest first (from the pipeline's history). Closed by default,
 * so the record stays short; one click opens it. */
export function RunLog({ entries }: { entries: LogEntry[] }) {
  const [open, setOpen] = useState(false);
  return (
    <section className="flex flex-col gap-[10px] rounded-[4px] border border-line bg-panel p-5">
      <button
        type="button"
        aria-expanded={open}
        onClick={() => setOpen((o) => !o)}
        className="flex items-center justify-between gap-2 text-left"
      >
        <h2 className="m-0 text-[16px] font-bold">Run log</h2>
        <span className="flex items-center gap-1 text-[13px] text-muted">
          {entries.length === 0 ? "nothing yet" : `${entries.length} ${entries.length === 1 ? "run" : "runs"}`}
          <ChevronDown aria-hidden="true" size={16} className={cn("transition-transform", !open && "-rotate-90")} />
        </span>
      </button>
      {open &&
        (entries.length === 0 ? (
          <p className="m-0 text-[14px] text-muted">Each step you run is recorded here.</p>
        ) : (
          <ol className="m-0 flex list-none flex-col gap-2 p-0">
            {entries.map((entry, index) => (
              <li key={index} className="border-l-2 border-field py-[2px] pl-[10px] text-[13px] leading-snug">
                <div className="font-mono text-[12px] text-muted">
                  {entry.label} · {formatTime(entry.ran_at)}
                </div>
                <div>{entry.text}</div>
              </li>
            ))}
          </ol>
        ))}
    </section>
  );
}

function formatTime(iso: string): string {
  const time = new Date(iso);
  return Number.isNaN(time.getTime()) ? iso : time.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
}
