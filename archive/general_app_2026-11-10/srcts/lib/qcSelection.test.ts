// Unit tests for inspecting cells on the QC page (run with `npm test`). Uses a small
// hand-made table; nothing here is PBMC3k-specific.
import { describe, expect, it } from "vitest";

import { cellLook, groupByLook, parseRange, selectCells, summarize, type QcCells } from "./qcSelection";

const cells: QcCells = {
  barcode: ["AAA-1", "BBB-1", "CCC-1", "DDD-1", "EEE-1", "FFF-1"],
  n_genes: [150, 200, 250, 299, 300, 900],
  total_counts: [400, 600, 700, 800, 900, 3000],
  pct_mt: [1.0, 7.5, 2.0, 12.0, 3.0, 4.0],
};

describe("selectCells", () => {
  it("selects from ≤ value < to on the chosen measure, by barcode", () => {
    const s = selectCells(cells, { metric: "n_genes", from: 200, to: 300 });
    expect([...s]).toEqual(["BBB-1", "CCC-1", "DDD-1"]); // 200 in, 299 in, 300 out
  });

  it("allows an open end", () => {
    expect([...selectCells(cells, { metric: "pct_mt", from: 5, to: null })]).toEqual(["BBB-1", "DDD-1"]);
    expect([...selectCells(cells, { metric: "total_counts", from: null, to: 700 })]).toEqual(["AAA-1", "BBB-1"]);
  });

  it("selects nothing without a range", () => {
    expect(selectCells(cells, null).size).toBe(0);
    expect(selectCells(cells, { metric: "n_genes", from: null, to: null }).size).toBe(0);
  });

  it("matches cells by barcode, not by position", () => {
    // The same cells in a different order select the same barcodes.
    const order = [5, 3, 1, 0, 4, 2];
    const shuffled: QcCells = {
      barcode: order.map((i) => cells.barcode[i]),
      n_genes: order.map((i) => cells.n_genes[i]),
      total_counts: order.map((i) => cells.total_counts[i]),
      pct_mt: order.map((i) => cells.pct_mt[i]),
    };
    const a = selectCells(cells, { metric: "n_genes", from: 200, to: 300 });
    const b = selectCells(shuffled, { metric: "n_genes", from: 200, to: 300 });
    expect(new Set(b)).toEqual(new Set(a));
  });
});

describe("both charts draw the same cells the same way", () => {
  const selected = selectCells(cells, { metric: "n_genes", from: 200, to: 300 });
  const removed = new Set(["AAA-1", "DDD-1"]);

  it("gives every cell one look from its barcode", () => {
    expect(cellLook("BBB-1", selected, removed)).toBe("selected_kept");
    expect(cellLook("DDD-1", selected, removed)).toBe("selected_removed");
    expect(cellLook("AAA-1", selected, removed)).toBe("removed");
    expect(cellLook("FFF-1", selected, removed)).toBe("kept");
  });

  it("puts the same barcodes in the selected layers whichever measure is plotted", () => {
    // Each chart groups the same table; the genes chart and the % mito chart only differ in
    // which value they plot on the x axis, so their layers hold the same barcodes.
    const forGenesChart = groupByLook(cells, selected, removed);
    const forMitoChart = groupByLook(cells, selected, removed);
    const barcodes = (idx: number[]) => idx.map((i) => cells.barcode[i]);
    expect(barcodes(forMitoChart.selected_kept)).toEqual(barcodes(forGenesChart.selected_kept));
    expect(barcodes(forMitoChart.selected_removed)).toEqual(barcodes(forGenesChart.selected_removed));
    expect(new Set([...barcodes(forGenesChart.selected_kept), ...barcodes(forGenesChart.selected_removed)])).toEqual(selected);
  });

  it("keeps the selection when a threshold changes; only the statuses change", () => {
    const stricter = new Set([...removed, "BBB-1"]); // e.g. a lower max % mito removes BBB-1 too
    const before = groupByLook(cells, selected, removed);
    const after = groupByLook(cells, selected, stricter);
    const inGroup = (g: Record<string, number[]>) => new Set([...g.selected_kept, ...g.selected_removed]);
    expect(inGroup(after)).toEqual(inGroup(before)); // the same cells are still selected
    expect(after.selected_removed.map((i) => cells.barcode[i])).toEqual(["BBB-1", "DDD-1"]);
  });
});

describe("summarize", () => {
  it("describes the selected cells and what the current settings would remove", () => {
    const selected = selectCells(cells, { metric: "n_genes", from: 200, to: 300 });
    const failures = { "DDD-1": ["max_pct_counts_mt"], "AAA-1": ["min_genes"], "BBB-1": ["max_pct_counts_mt", "max_n_genes_by_counts"] };
    const s = summarize(cells, selected, failures);
    expect(s.n).toBe(3);
    expect(s.spread.n_genes).toEqual({ min: 200, median: 250, max: 299 });
    expect(s.spread.pct_mt).toEqual({ min: 2.0, median: 7.5, max: 12.0 });
    expect(s.nRemoved).toBe(2); // AAA-1 isn't selected
    expect(s.failing).toEqual({ max_pct_counts_mt: 2, max_n_genes_by_counts: 1 });
  });

  it("handles an empty selection", () => {
    const s = summarize(cells, new Set(), {});
    expect(s.n).toBe(0);
    expect(s.spread.n_genes).toBeNull();
  });
});

describe("parseRange", () => {
  it("reads the boxes, with empty = open end", () => {
    expect(parseRange("n_genes", "200", "300")).toEqual({ metric: "n_genes", from: 200, to: 300 });
    expect(parseRange("n_genes", "", "300")).toEqual({ metric: "n_genes", from: null, to: 300 });
    expect(parseRange("n_genes", "", "")).toBeNull();
  });

  it("explains unusable ranges", () => {
    expect(parseRange("n_genes", "300", "200")).toEqual({ error: "'From' must be smaller than 'To'." });
    expect(parseRange("n_genes", "abc", "")).toEqual({ error: "Enter numbers for the range." });
  });
});
