// While one of the longer steps runs (the Run buttons of Preprocessing through Cell Type Annotation and
// Compare Conditions, or "Run the remaining steps"), the whole screen is greyed out with a small message
// in the middle: what is running and for how long (analyst request, 2026-11-16). Quick actions (Run QC,
// Apply Filtering, loading data, live previews, gene look-ups) don't use it (analyst, same day).
//
// Pages report what they are running through `useBusyOverlay(label)`; App collects the labels and shows
// one overlay. The overlay appears after a short delay so quick steps don't flash.
import { createContext, useContext, useEffect, useId, useState } from "react";

/** Set (or clear, with null) the busy message of one part of the page. Provided by App. */
export const BusyContext = createContext<(key: string, label: string | null) => void>(() => {});

/** Show the overlay with `label` while it isn't null ("Running Clustering…"). */
export function useBusyOverlay(label: string | null) {
  const set = useContext(BusyContext);
  const key = useId();
  useEffect(() => {
    set(key, label);
    return () => set(key, null);
  }, [label]); // eslint-disable-line react-hooks/exhaustive-deps
}

/** "Running Clustering…" for a step's key; the step names as on screen. */
export const STEP_NAMES: Record<string, string> = {
  data_overview: "Data Overview",
  quality_control: "Quality Control",
  preprocessing: "Preprocessing",
  dimensionality_reduction: "Dimensionality Reduction",
  clustering: "Clustering",
  marker_genes: "Marker Genes",
  cell_type_annotation: "Cell Type Annotation",
  compare_conditions: "Compare Conditions",
};

const SHOW_AFTER_MS = 300;

export function BusyOverlay({ label }: { label: string | null }) {
  const [since, setSince] = useState<number | null>(null);
  const [now, setNow] = useState(() => Date.now());
  const busy = label !== null;

  useEffect(() => {
    if (!busy) {
      setSince(null);
      return;
    }
    const start = Date.now();
    setSince(start);
    setNow(start);
    const timer = window.setInterval(() => setNow(Date.now()), 250);
    return () => window.clearInterval(timer);
  }, [busy]);

  if (!busy || since === null || now - since < SHOW_AFTER_MS) return null;
  const seconds = Math.floor((now - since) / 1000);
  return (
    <div className="fixed inset-0 z-[60] flex items-center justify-center bg-[#18181b]/25" aria-busy="true">
      <div role="status" aria-live="polite" className="flex items-center gap-3 border border-field bg-white px-5 py-3 shadow-sm">
        <span aria-hidden="true" className="size-4 flex-none animate-spin rounded-full border-2 border-line border-t-ink" />
        <span className="text-[14px] font-bold text-ink">{label}</span>
        {seconds > 0 && <span className="font-mono text-[13px] text-muted">{seconds} s</span>}
      </div>
    </div>
  );
}
