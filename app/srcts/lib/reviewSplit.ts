// Splits the clusters of Cell Type Annotation into the two review tabs (analyst, 2026-11-13):
// "To review" (not yet labeled by the user, in the server's groups and order) and "Labeled clusters"
// (labeled by the user, in cluster order). The user confirms or renames every cluster.

export interface Reviewable {
  cluster: string;
  review_group: string;
}

export interface ReviewSplit<T extends Reviewable> {
  toReview: { key: string; label: string; description: string; clusters: T[] }[]; // non-empty groups only
  labeled: T[];
}

/** Cluster names in natural order: "2" before "10"; names that aren't numbers sort as text. */
export function byCluster(a: string, b: string): number {
  return a.localeCompare(b, undefined, { numeric: true });
}

export function splitForReview<T extends Reviewable>(
  clusters: T[],
  reviewed: Set<string> | Record<string, unknown>,
  groups: { key: string; label: string; description: string }[],
): ReviewSplit<T> {
  const isReviewed = (c: string) => (reviewed instanceof Set ? reviewed.has(c) : c in reviewed);
  const open = clusters.filter((c) => !isReviewed(c.cluster));
  const known = new Set(groups.map((g) => g.key));
  const toReview = groups
    .map((g) => ({ ...g, clusters: open.filter((c) => c.review_group === g.key).sort((a, b) => byCluster(a.cluster, b.cluster)) }))
    .filter((g) => g.clusters.length > 0);
  // A cluster in a group the page doesn't know must never disappear: list it at the end.
  const other = open.filter((c) => !known.has(c.review_group));
  if (other.length) toReview.push({ key: "other", label: "Other", description: "", clusters: other });
  const labeled = clusters.filter((c) => isReviewed(c.cluster)).sort((a, b) => byCluster(a.cluster, b.cluster));
  return { toReview, labeled };
}
