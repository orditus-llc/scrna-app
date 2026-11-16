// Quality Control with several samples (Milestone 6b): each sample's own cell limits, and how many
// cells each sample keeps. Only shown when the data has several samples.
import type { FilterSettings, PreviewPayload, SampleLimitBoxes } from "../types";
import { CsvDownload } from "./DownloadMenu";
import { InfoTip } from "./InfoTip";

/** One sample's limits as boxes. */
export function sampleLimitBoxes(s: NonNullable<FilterSettings["sample_limits"]>[number]): SampleLimitBoxes {
  return {
    sample: s.sample,
    min_genes: String(s.min_genes),
    max_genes: s.max_n_genes_by_counts === null ? "" : String(s.max_n_genes_by_counts),
    max_genes_off: s.max_n_genes_by_counts === null,
    max_mt: s.max_pct_counts_mt === null ? "" : String(s.max_pct_counts_mt),
    max_mt_off: s.max_pct_counts_mt === null,
  };
}

/** One sample's limits as settings (for the chart lines); empty or "no limit" boxes are null. */
export function sampleLimitSettings(b: SampleLimitBoxes, base: FilterSettings): FilterSettings {
  const num = (x: string) => (x.trim() === "" || Number.isNaN(Number(x)) ? null : Number(x));
  return {
    ...base,
    min_genes: num(b.min_genes) ?? 0,
    max_n_genes_by_counts: b.max_genes_off ? null : num(b.max_genes),
    max_pct_counts_mt: b.max_mt_off ? null : num(b.max_mt),
    sample_limits: [],
  };
}

/** Each sample's cell limits, one row per sample. An empty max box = no limit. */
export function SampleLimitsTable({
  rows,
  suggested,
  mitoAvailable,
  onChange,
}: {
  rows: SampleLimitBoxes[];
  suggested: FilterSettings["sample_limits"];
  mitoAvailable: boolean;
  onChange: (rows: SampleLimitBoxes[]) => void;
}) {
  const bySample = Object.fromEntries((suggested ?? []).map((s) => [s.sample, s]));
  const set = (i: number, patch: Partial<SampleLimitBoxes>) => onChange(rows.map((r, j) => (j === i ? { ...r, ...patch } : r)));
  const cell = "rounded-none border bg-white px-2 py-[5px] text-right font-mono text-[13px] w-full min-w-[64px]";
  const changed = (value: string, ref: number | null | undefined) =>
    ref !== undefined && (value.trim() === "" ? ref !== null : Number(value) !== ref);
  return (
    <div className="flex flex-col gap-2 border-b border-line py-[14px]">
      <span className="text-[14px] leading-snug">
        <span className="font-bold">Cell limits per sample</span>{" "}
        <InfoTip label="About cell limits per sample">
          Each sample's limits are suggested from its own cells, since samples differ in sequencing depth and quality
          (single-cell best-practices guide). Leave a max box empty for no limit. Teal = changed from the suggestion.
        </InfoTip>
      </span>
      <div className="max-h-[420px] overflow-auto">
        <table className="w-full border-collapse text-[13px]">
          <thead className="sticky top-0 bg-white">
            <tr className="text-left text-[12px] tracking-[0.04em] text-muted uppercase">
              <th className="py-1 pr-2 font-bold">Sample</th>
              <th className="py-1 pr-2 text-right font-bold">Min genes</th>
              <th className="py-1 pr-2 text-right font-bold">Max genes</th>
              {mitoAvailable && <th className="py-1 text-right font-bold">Max % mito</th>}
            </tr>
          </thead>
          <tbody>
            {rows.map((r, i) => {
              const ref = bySample[r.sample];
              return (
                <tr key={r.sample} className="border-t border-line">
                  <td className="max-w-[110px] py-1 pr-2 leading-tight break-words text-body">
                    {r.sample}
                  </td>
                  <td className="py-1 pr-2">
                    <input
                      aria-label={`Min genes per cell for ${r.sample}`}
                      inputMode="numeric"
                      value={r.min_genes}
                      onChange={(e) => set(i, { min_genes: e.target.value })}
                      className={`${cell} ${changed(r.min_genes, ref?.min_genes) ? "border-2 border-accent" : "border-field"}`}
                    />
                  </td>
                  <td className="py-1 pr-2">
                    <input
                      aria-label={`Max genes per cell for ${r.sample} (empty = no limit)`}
                      inputMode="numeric"
                      value={r.max_genes_off ? "" : r.max_genes}
                      placeholder="none"
                      onChange={(e) => set(i, { max_genes: e.target.value, max_genes_off: e.target.value.trim() === "" })}
                      className={`${cell} ${changed(r.max_genes_off ? "" : r.max_genes, ref?.max_n_genes_by_counts) ? "border-2 border-accent" : "border-field"}`}
                    />
                  </td>
                  {mitoAvailable && (
                    <td className="py-1">
                      <input
                        aria-label={`Max % mitochondrial for ${r.sample} (empty = no limit)`}
                        inputMode="decimal"
                        value={r.max_mt_off ? "" : r.max_mt}
                        placeholder="none"
                        onChange={(e) => set(i, { max_mt: e.target.value, max_mt_off: e.target.value.trim() === "" })}
                        className={`${cell} ${changed(r.max_mt_off ? "" : r.max_mt, ref?.max_pct_counts_mt) ? "border-2 border-accent" : "border-field"}`}
                      />
                    </td>
                  )}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}

type OkPreview = Extract<PreviewPayload, { ok: true }>;

/** How many cells each sample keeps with these settings, next to the suggestion. */
export function KeptPerSample({ preview }: { preview: OkPreview }) {
  if (preview.by_sample.length === 0) return null;
  const header = "py-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase";
  const pct = (kept: number, before: number) => `${((100 * (before - kept)) / before).toFixed(1)}%`;
  return (
    <section aria-label="Kept per sample" className="flex flex-col gap-3">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <h2 className="m-0 text-[16px] font-bold">Kept per sample</h2>
        <CsvDownload
          name="Quality Control kept per sample"
          columns={[
            { key: "sample", header: "Sample" },
            { key: "before", header: "Cells before" },
            { key: "kept", header: "Cells kept" },
            { key: "removed", header: "Cells removed" },
            { key: "suggested_kept", header: "Kept at the suggestion" },
          ]}
          rows={() =>
            preview.by_sample.map((r) => ({
              sample: r.sample,
              before: r.n_cells_before,
              kept: r.n_cells_kept,
              removed: r.n_cells_before - r.n_cells_kept,
              suggested_kept: r.suggested_kept,
            }))
          }
        />
      </div>
      <div className="overflow-x-auto">
        <table className="w-full min-w-[560px] border-collapse text-[14px]">
          <thead>
            <tr className="border-b border-field text-left">
              <th scope="col" className={`${header} pr-3`}>Sample</th>
              <th scope="col" className={`${header} px-3 text-right`}>Before</th>
              <th scope="col" className={`${header} px-3 text-right`}>Kept</th>
              <th scope="col" className={`${header} px-3 text-right`}>Removed</th>
              <th scope="col" className={`${header} pl-3 text-right`}>Kept at the suggestion</th>
            </tr>
          </thead>
          <tbody>
            {preview.by_sample.map((r) => {
              const differs = r.suggested_kept !== null && r.suggested_kept !== r.n_cells_kept;
              return (
                <tr key={r.sample} className={`border-b border-line ${differs ? "bg-accent-tint" : ""}`}>
                  <td className={`py-2 pr-3 font-bold ${differs ? "pl-[10px] text-accent-text" : "text-ink"}`}>{r.sample}</td>
                  <td className="px-3 py-2 text-right font-mono">{r.n_cells_before.toLocaleString("en-US")}</td>
                  <td className="px-3 py-2 text-right font-mono">{r.n_cells_kept.toLocaleString("en-US")}</td>
                  <td className="px-3 py-2 text-right font-mono text-muted">
                    {(r.n_cells_before - r.n_cells_kept).toLocaleString("en-US")} ({pct(r.n_cells_kept, r.n_cells_before)})
                  </td>
                  <td className="py-2 pl-3 text-right font-mono text-muted">
                    {r.suggested_kept === null ? "–" : r.suggested_kept.toLocaleString("en-US")}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
      <p className="m-0 text-[13px] text-muted">Teal = a sample keeps a different number of cells than at the suggestion.</p>
    </section>
  );
}
