// Inspecting a group of cells on the Quality Control page. Pure functions (no React, no
// Plotly) so they can be unit-tested: which cells a range selects, how each cell is drawn,
// and a summary of the selected cells. Cells are always matched by barcode.

export type QcMetric = "n_genes" | "total_counts" | "pct_mt";

export interface QcCells {
  barcode: string[];
  n_genes: number[];
  total_counts: number[];
  pct_mt: number[];
}

/** An inspection range on one QC measure: from ≤ value < to. A missing end is open. */
export interface InspectRange {
  metric: QcMetric;
  from: number | null;
  to: number | null;
}

export const METRIC_LABELS: Record<QcMetric, string> = {
  n_genes: "Genes per cell",
  total_counts: "Counts per cell",
  pct_mt: "% mitochondrial",
};

/** The barcodes of the cells whose value is in the range (from ≤ value < to). An empty
 * range (no ends) selects nothing: an inspection always needs at least one end. */
export function selectCells(cells: QcCells, range: InspectRange | null): Set<string> {
  const selected = new Set<string>();
  if (!range || (range.from === null && range.to === null)) return selected;
  const values = cells[range.metric];
  cells.barcode.forEach((barcode, i) => {
    const v = values[i];
    if (v === null || v === undefined || Number.isNaN(v)) return;
    if (range.from !== null && v < range.from) return;
    if (range.to !== null && v >= range.to) return;
    selected.add(barcode);
  });
  return selected;
}

/** How a cell is drawn: kept or would be removed (by the current filter preview), and
 * whether it is in the inspected group. Every chart uses this, so all charts agree. */
export type CellLook = "kept" | "removed" | "selected_kept" | "selected_removed";

export function cellLook(barcode: string, selected: Set<string>, removed: Set<string>): CellLook {
  const isRemoved = removed.has(barcode);
  if (selected.has(barcode)) return isRemoved ? "selected_removed" : "selected_kept";
  return isRemoved ? "removed" : "kept";
}

/** The positions (in `cells`) of the cells with each look, so a chart can draw one layer
 * per look (selected cells on top). Positions only index the same QC table the barcodes
 * came from; membership is decided by barcode. */
export function groupByLook(cells: QcCells, selected: Set<string>, removed: Set<string>): Record<CellLook, number[]> {
  const groups: Record<CellLook, number[]> = { kept: [], removed: [], selected_kept: [], selected_removed: [] };
  cells.barcode.forEach((barcode, i) => groups[cellLook(barcode, selected, removed)].push(i));
  return groups;
}

export interface Spread {
  min: number;
  median: number;
  max: number;
}

export interface SelectionSummary {
  n: number;
  spread: Record<QcMetric, Spread | null>;
  nRemoved: number;
  failing: Record<string, number>; // filter -> selected cells failing it (a cell can fail several)
}

/** A compact summary of the selected cells under the current filter preview. */
export function summarize(
  cells: QcCells,
  selected: Set<string>,
  failures: Record<string, string[]>,
): SelectionSummary {
  const idx = cells.barcode.flatMap((b, i) => (selected.has(b) ? [i] : []));
  const spread = (metric: QcMetric): Spread | null => {
    const v = idx.map((i) => cells[metric][i]).filter((x) => x !== null && !Number.isNaN(x)).sort((a, b) => a - b);
    if (v.length === 0) return null;
    const mid = Math.floor(v.length / 2);
    return { min: v[0], median: v.length % 2 ? v[mid] : (v[mid - 1] + v[mid]) / 2, max: v[v.length - 1] };
  };
  const failing: Record<string, number> = {};
  let nRemoved = 0;
  for (const i of idx) {
    const reasons = failures[cells.barcode[i]];
    if (!reasons || reasons.length === 0) continue;
    nRemoved += 1;
    for (const r of reasons) failing[r] = (failing[r] ?? 0) + 1;
  }
  return {
    n: idx.length,
    spread: { n_genes: spread("n_genes"), total_counts: spread("total_counts"), pct_mt: spread("pct_mt") },
    nRemoved,
    failing,
  };
}

/** Turn the typed range boxes into a range ("" = open end). Returns a message instead
 * when the boxes can't be used. */
export function parseRange(metric: QcMetric, from: string, to: string): InspectRange | { error: string } | null {
  const num = (s: string) => (s.trim() === "" ? null : Number(s));
  const a = num(from);
  const b = num(to);
  if (a === null && b === null) return null;
  if ((a !== null && !Number.isFinite(a)) || (b !== null && !Number.isFinite(b))) return { error: "Enter numbers for the range." };
  if (a !== null && b !== null && a >= b) return { error: "'From' must be smaller than 'To'." };
  return { metric, from: a, to: b };
}
