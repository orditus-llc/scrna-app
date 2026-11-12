// The analyst's reviewed labels as UMAP groups: each cell gets its cluster's reviewed label (or
// "Not reviewed"); each label takes the color of its first cluster, "Not reviewed" is grey.
// Used by Cell Type Annotation and Explore Results.

export const NOT_REVIEWED_COLOR = "#a1a1aa";

export function reviewedGroups(
  cellCluster: string[],
  clusterOrder: string[],
  reviewed: Record<string, string>,
  clusterColors: Record<string, string>,
  notReviewed: string,
): { cellLabel: string[]; groups: { key: string; label: string; n: number }[]; colors: Record<string, string> } {
  const cellLabel = cellCluster.map((c) => reviewed[c] ?? notReviewed);
  const counts: Record<string, number> = {};
  cellLabel.forEach((l) => (counts[l] = (counts[l] ?? 0) + 1));
  const colors: Record<string, string> = {};
  const order: string[] = [];
  for (const c of clusterOrder) {
    const label = reviewed[c];
    if (label !== undefined && !(label in colors)) {
      colors[label] = clusterColors[c];
      order.push(label);
    }
  }
  if (counts[notReviewed]) {
    colors[notReviewed] = NOT_REVIEWED_COLOR;
    order.push(notReviewed);
  }
  return { cellLabel, groups: order.map((l) => ({ key: l, label: l, n: counts[l] ?? 0 })), colors };
}
