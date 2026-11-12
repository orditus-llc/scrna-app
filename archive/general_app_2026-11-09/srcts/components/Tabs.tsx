// Tabs above a plot area: one view at a time, so a page shows one main plot instead of a
// stack. Keyboard: Left/Right arrows move between tabs (the WAI-ARIA tabs pattern).
import { useId, type KeyboardEvent, type ReactNode } from "react";

import { cn } from "../lib/utils";

export interface TabItem<K extends string> {
  key: K;
  label: string;
  disabled?: boolean;
}

export function Tabs<K extends string>({
  tabs,
  active,
  onChange,
  label,
  right,
  children,
}: {
  tabs: TabItem<K>[];
  active: K;
  onChange: (key: K) => void;
  label: string; // what the tabs switch between, for screen readers
  right?: ReactNode; // controls shown at the right end of the tab row
  children: ReactNode; // the active tab's content
}) {
  const id = useId();
  const enabled = tabs.filter((t) => !t.disabled);
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== "ArrowRight" && e.key !== "ArrowLeft") return;
    const i = enabled.findIndex((t) => t.key === active);
    const next = enabled[(i + (e.key === "ArrowRight" ? 1 : enabled.length - 1)) % enabled.length];
    onChange(next.key);
    document.getElementById(`${id}-${next.key}`)?.focus();
  };
  return (
    <div className="flex min-w-0 flex-col gap-3">
      <div className="flex flex-wrap items-end justify-between gap-3 border-b border-ink">
        <div role="tablist" aria-label={label} onKeyDown={onKey} className="flex flex-wrap gap-6">
          {tabs.map((t) => (
            <button
              key={t.key}
              id={`${id}-${t.key}`}
              type="button"
              role="tab"
              aria-selected={t.key === active}
              aria-controls={`${id}-panel`}
              tabIndex={t.key === active ? 0 : -1}
              disabled={t.disabled}
              onClick={() => onChange(t.key)}
              className={cn(
                "-mb-px border-b-2 pt-1 pb-[9px] text-[14px] disabled:cursor-not-allowed disabled:text-off",
                t.key === active ? "border-ink font-extrabold text-ink" : "border-transparent font-semibold text-muted enabled:hover:text-ink",
              )}
            >
              {t.label}
            </button>
          ))}
        </div>
        {right && <div className="pb-2">{right}</div>}
      </div>
      <div id={`${id}-panel`} role="tabpanel" aria-labelledby={`${id}-${active}`} className="flex min-w-0 flex-col gap-3">
        {children}
      </div>
    </div>
  );
}
