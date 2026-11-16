// Unit tests for the review tabs' split (run with `npm test`). Hand-made, neutral names.
import { describe, expect, it } from "vitest";

import { byCluster, splitForReview } from "./reviewSplit";

const groups = [
  { key: "g1", label: "First", description: "d1" },
  { key: "g2", label: "Second", description: "d2" },
  { key: "g3", label: "Third", description: "d3" },
];

describe("splitForReview", () => {
  const clusters = [
    { cluster: "0", review_group: "g2" },
    { cluster: "10", review_group: "g1" },
    { cluster: "2", review_group: "g1" },
    { cluster: "3", review_group: "g2" },
  ];
  it("groups clusters not yet labeled in the given group order, empty groups left out", () => {
    const s = splitForReview(clusters, {}, groups);
    expect(s.toReview.map((g) => [g.key, g.clusters.map((c) => c.cluster)])).toEqual([
      ["g1", ["2", "10"]],
      ["g2", ["0", "3"]],
    ]);
    expect(s.labeled).toEqual([]);
  });
  it("moves labeled clusters to the labeled list, in cluster order", () => {
    const s = splitForReview(clusters, { "10": {}, "0": {}, "2": {} }, groups);
    expect(s.toReview.map((g) => [g.key, g.clusters.map((c) => c.cluster)])).toEqual([["g2", ["3"]]]);
    expect(s.labeled.map((c) => c.cluster)).toEqual(["0", "2", "10"]);
  });
  it("never loses a cluster whose group is unknown", () => {
    const s = splitForReview([{ cluster: "5", review_group: "new" }], new Set<string>(), groups);
    expect(s.toReview).toEqual([{ key: "other", label: "Other", description: "", clusters: [{ cluster: "5", review_group: "new" }] }]);
  });
  it("sorts cluster names naturally", () => {
    expect(["10", "2", "1", "b", "a"].sort(byCluster)).toEqual(["1", "2", "10", "a", "b"]);
  });
});
