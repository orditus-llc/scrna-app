// Unit tests for the reviewed-label groups (run with `npm test`). Hand-made, neutral names
// (no real cluster numbers or cell types, so nothing resembles a cluster -> label answer key).
import { describe, expect, it } from "vitest";

import { NOT_REVIEWED_COLOR, reviewedGroups } from "./reviewedGroups";

describe("reviewedGroups", () => {
  const colors = { x: "#111", y: "#222", z: "#333" };
  it("labels cells by their cluster's review; clusters with the same label share one group", () => {
    const r = reviewedGroups(["x", "y", "z", "x", "z"], ["x", "y", "z"], { x: "Type A", z: "Type A" }, colors, "Not reviewed");
    expect(r.cellLabel).toEqual(["Type A", "Not reviewed", "Type A", "Type A", "Type A"]);
    expect(r.groups).toEqual([
      { key: "Type A", label: "Type A", n: 4 },
      { key: "Not reviewed", label: "Not reviewed", n: 1 },
    ]);
    expect(r.colors).toEqual({ "Type A": "#111", "Not reviewed": NOT_REVIEWED_COLOR });
  });
  it("leaves out 'Not reviewed' when every cluster is reviewed", () => {
    const r = reviewedGroups(["y", "x"], ["x", "y"], { x: "Type A", y: "Type B" }, colors, "Not reviewed");
    expect(r.groups.map((g) => g.key)).toEqual(["Type A", "Type B"]);
  });
});
