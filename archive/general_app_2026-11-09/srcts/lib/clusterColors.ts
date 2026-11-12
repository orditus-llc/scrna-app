// Cluster colors: CELLxGENE Explorer's scheme (d3 interpolateRainbow at i / n for cluster i
// of n), computed from the actual number of clusters, never a fixed list.
import { interpolateRainbow } from "d3-scale-chromatic";

export function clusterColors(clusters: string[]): Record<string, string> {
  const n = clusters.length;
  return Object.fromEntries(clusters.map((c, i) => [c, interpolateRainbow(i / n)]));
}
