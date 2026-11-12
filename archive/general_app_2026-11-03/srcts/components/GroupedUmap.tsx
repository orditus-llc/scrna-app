// A UMAP colored by group (cluster or cell type) with the list of groups beside it. Hovering
// a cell or a group in the list keeps that group in full color and fades the rest; clicking
// keeps the highlight (Clear or Esc removes it). Used by the Clustering and Cell Type pages.
import type { Data, Layout } from "plotly.js-dist-min";
import { useEffect, useMemo, useState, type ReactNode } from "react";

import { cn } from "../lib/utils";
import { PlotlyChart, type ChartPoint } from "./PlotlyChart";

export interface UmapGroup {
  key: string;
  label: string; // e.g. "Cluster 3" or "B cells"
  n: number; // cells in the group
}

export interface UmapCells {
  barcode: string[];
  n_genes: number[];
  total_counts: (number | null)[];
  pct_mt: (number | null)[] | null;
}

export function GroupedUmap({
  title,
  description,
  listHeading,
  groups,
  colors,
  cellGroup,
  umap,
  cells,
  extraHover,
  children,
  details,
}: {
  title: string;
  description: ReactNode;
  listHeading: string;
  groups: UmapGroup[]; // in display order
  colors: Record<string, string>;
  cellGroup: string[]; // group key of each cell, same order as umap and cells
  umap: { x: (number | null)[]; y: (number | null)[] };
  cells: UmapCells;
  extraHover?: (cellIndex: number) => string; // another line in the hover box (e.g. the cluster)
  children?: ReactNode; // notes shown under the description
  details?: (shown: string | null, pinned: boolean) => ReactNode; // a panel under the list for the highlighted group
}) {
  const order = groups.map((g) => g.key);
  const byKey = useMemo(() => Object.fromEntries(groups.map((g) => [g.key, g])), [groups]);
  const total = cellGroup.length;
  const largest = Math.max(...groups.map((g) => g.n), 1);
  const [hovered, setHovered] = useState<string | null>(null);
  const [pinned, setPinned] = useState<string | null>(null);
  const shown = pinned ?? hovered;

  // Esc clears a kept highlight.
  useEffect(() => {
    if (!pinned) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && setPinned(null);
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [pinned]);

  // One trace per group, so a whole group can be highlighted at once.
  const members = useMemo(() => {
    const m: Record<string, number[]> = Object.fromEntries(order.map((k) => [k, []]));
    cellGroup.forEach((k, i) => m[k]?.push(i));
    return m;
  }, [cellGroup, groups]); // eslint-disable-line react-hooks/exhaustive-deps

  const chart = useMemo(() => {
    const pct = cells.pct_mt;
    const data: Data[] = order.map((key) => {
      const idx = members[key];
      const faded = shown !== null && shown !== key;
      return {
        type: "scattergl",
        mode: "markers",
        name: byKey[key].label,
        x: idx.map((i) => umap.x[i]),
        y: idx.map((i) => umap.y[i]),
        customdata: idx.map((i) => [cells.barcode[i], cells.n_genes[i], cells.total_counts[i], pct ? pct[i] : null, extraHover?.(i) ?? ""]),
        marker: {
          size: shown === key ? 5.5 : 4.5,
          color: colors[key],
          opacity: faded ? 0.08 : 0.9,
          line: { color: "rgba(17,17,17,0.35)", width: faded ? 0 : 0.4 },
        },
        hovertemplate:
          `<b>%{customdata[0]}</b><br>${byKey[key].label} · ${byKey[key].n.toLocaleString("en-US")} cells<br>` +
          (extraHover ? "%{customdata[4]}<br>" : "") +
          "Genes %{customdata[1]:,}<br>Counts %{customdata[2]:,.0f}" +
          (pct ? "<br>% mito %{customdata[3]:.2f}" : "") +
          "<extra></extra>",
      };
    });
    const layout: Partial<Layout> = {
      xaxis: { title: { text: "UMAP 1" }, showticklabels: false, showgrid: false, zeroline: false },
      yaxis: { title: { text: "UMAP 2" }, showticklabels: false, showgrid: false, zeroline: false, scaleanchor: "x" },
      margin: { l: 40, r: 16, t: 16, b: 40 },
      hoverlabel: { bgcolor: "#ffffff", bordercolor: "#d4d4d8", font: { color: "#111111" } },
    };
    // For a saved image: every group at full color (no hover highlight), named with its size.
    const exportData: Data[] = data.map((trace, i) => {
      const t = trace as { marker: Record<string, unknown> };
      const g = byKey[order[i]];
      return {
        ...trace,
        name: `${g.label} (${g.n.toLocaleString("en-US")})`,
        marker: { ...t.marker, size: 4.5, opacity: 0.9, line: { color: "rgba(17,17,17,0.35)", width: 0.4 } },
      } as Data;
    });
    return { data, layout, exportData };
  }, [members, colors, shown, byKey, umap, cells]); // eslint-disable-line react-hooks/exhaustive-deps

  const fromPoint = (p: ChartPoint) => order[p.curveNumber] ?? null;
  const toggle = (key: string | null) => setPinned((current) => (current === key ? null : key));

  return (
    <>
      <div className="flex flex-wrap items-end justify-between gap-2 border-b border-ink pb-[10px]">
        <h2 className="m-0 text-lg font-extrabold">{title}</h2>
        {pinned !== null && (
          <span className="text-[13px]">
            Showing <b>{byKey[pinned]?.label}</b> ·{" "}
            <button type="button" onClick={() => setPinned(null)} className="font-bold underline underline-offset-4">
              Clear
            </button>{" "}
            <span className="text-muted">(Esc)</span>
          </span>
        )}
      </div>
      <p className="m-0 text-[13px] text-body">{description}</p>
      {children}
      <div className="flex flex-wrap gap-6">
        <div className="min-w-0 flex-[999_1_420px]">
          <PlotlyChart
            data={chart.data}
            layout={chart.layout}
            height={500}
            label={`UMAP of all cells colored by ${listHeading.toLowerCase()}`}
            download={{ name: title, title, legend: true, data: chart.exportData }}
            onHover={(p) => setHovered(fromPoint(p))}
            onUnhover={() => setHovered(null)}
            onClick={(p) => toggle(fromPoint(p))}
          />
        </div>
        <div className="flex min-w-0 flex-[1_1_230px] flex-col gap-4">
        <ol aria-label={listHeading} className="m-0 flex list-none flex-col p-0">
          <li className="flex justify-between border-b border-ink pb-2 text-[12px] font-bold tracking-[0.04em] text-muted uppercase">
            <span>{listHeading}</span>
            <span>Cells</span>
          </li>
          {groups.map((g) => {
            const active = shown === g.key;
            return (
              <li key={g.key} className="border-b border-line">
                <button
                  type="button"
                  aria-pressed={pinned === g.key}
                  onMouseEnter={() => setHovered(g.key)}
                  onMouseLeave={() => setHovered(null)}
                  onFocus={() => setHovered(g.key)}
                  onBlur={() => setHovered(null)}
                  onClick={() => toggle(g.key)}
                  className={cn(
                    "flex w-full flex-col gap-[6px] px-2 py-[9px] text-left",
                    active ? "bg-chip" : shown !== null ? "opacity-45" : "hover:bg-chip",
                  )}
                >
                  <span className="flex items-center justify-between gap-3 text-[14px]">
                    <span className="flex min-w-0 items-center gap-[10px] font-bold">
                      <span aria-hidden="true" className="size-3 flex-none rounded-full ring-1 ring-black/25" style={{ background: colors[g.key] }} />
                      <span className="truncate" title={g.label}>{g.label}</span>
                    </span>
                    <span className="font-mono whitespace-nowrap tabular-nums">
                      {g.n.toLocaleString("en-US")} <span className="text-muted">· {percent(g.n, total)}</span>
                    </span>
                  </span>
                  <span aria-hidden="true" className="block h-1 bg-line">
                    <span className="block h-1" style={{ width: `${(100 * g.n) / largest}%`, background: colors[g.key] }} />
                  </span>
                </button>
              </li>
            );
          })}
        </ol>
        {details?.(shown, pinned !== null)}
        </div>
      </div>
    </>
  );
}

export function percent(n: number, total: number): string {
  const p = (100 * n) / total;
  return p < 1 ? `${p.toFixed(1)}%` : `${Math.round(p)}%`;
}
