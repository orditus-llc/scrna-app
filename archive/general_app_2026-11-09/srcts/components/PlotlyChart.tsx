// A small wrapper that draws a Plotly chart, updates it in place when its data changes, and
// offers it for download (PNG for slides, SVG for papers) with a title and caption.
import type { Data, Layout } from "plotly.js-dist-min";
import Plotly from "plotly.js-dist-min";
import { useContext, useEffect, useRef } from "react";

import { fileName, saveFile } from "../lib/download";
import { DownloadContext, DownloadMenu } from "./DownloadMenu";

export const PLOT_FONT = { family: "Manrope, 'Helvetica Neue', sans-serif", size: 12, color: "#3f3f46" };

/** True when the browser can draw WebGL (fast scatter plots for many cells). */
const WEBGL_AVAILABLE = (() => {
  try {
    const canvas = document.createElement("canvas");
    return Boolean(canvas.getContext("webgl2") || canvas.getContext("webgl"));
  } catch {
    return false;
  }
})();

/** Use WebGL scatter plots where possible, and standard ones where WebGL is unavailable
 * (some older or remote computers). Same look either way. */
function withFallback(data: Data[]): Data[] {
  if (WEBGL_AVAILABLE) return data;
  return data.map((trace) => ((trace as { type?: string }).type === "scattergl" ? { ...trace, type: "scatter" } : trace) as Data);
}

/** The layout every chart starts from (the chart's own layout is laid over it). */
function baseLayout(height: number, layout: Partial<Layout>): Partial<Layout> {
  return {
    height,
    margin: { l: 56, r: 16, t: 36, b: 48 },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: PLOT_FONT,
    showlegend: false,
    ...layout,
  };
}

const CAPTION_SPACE = 70; // px above the chart for the title and caption in a saved image

/** The chart as an image file: a title and a caption (dataset · date) on top, white
 * background, the chart's current width. PNG at 3x resolution; SVG as vectors (fast WebGL
 * dot plots are drawn as ordinary SVG dots so they stay sharp at any size). */
async function exportChart(opts: {
  data: Data[];
  layout: Partial<Layout>;
  height: number;
  width: number;
  title: string;
  caption: string;
  format: "png" | "svg";
  legend: boolean;
}): Promise<string> {
  const escape = (t: string) => t.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  const base = baseLayout(opts.height, opts.layout);
  const margin = base.margin ?? {};
  const height = opts.height + CAPTION_SPACE;
  const layout: Partial<Layout> = {
    ...base,
    height,
    width: opts.width,
    paper_bgcolor: "#ffffff",
    plot_bgcolor: "#ffffff",
    margin: { ...margin, t: (margin.t ?? 36) + CAPTION_SPACE, r: opts.legend ? 200 : margin.r },
    // A key for charts whose colors are groups (the on-screen key is a list beside the chart).
    ...(opts.legend
      ? {
          showlegend: true,
          legend: { x: 1.02, xanchor: "left", y: 1, yanchor: "top", font: { ...PLOT_FONT, size: 12 }, itemsizing: "constant" },
        }
      : {}),
    title: {
      text: `<b>${escape(opts.title)}</b><br><span style="font-size:12px;color:#63636b">${escape(opts.caption)}</span>`,
      x: 0,
      xref: "paper",
      xanchor: "left",
      yref: "container",
      y: 1,
      yanchor: "top",
      pad: { t: 18 },
      font: { ...PLOT_FONT, size: 16, color: "#111111" },
    },
  };
  const traces =
    opts.format === "svg"
      ? opts.data.map((t) => ((t as { type?: string }).type === "scattergl" ? { ...t, type: "scatter" } : t) as Data)
      : withFallback(opts.data);
  return Plotly.toImage({ data: traces, layout }, { format: opts.format, width: opts.width, height, scale: opts.format === "png" ? 3 : 1 });
}

/** The point under the mouse: which trace (curveNumber) and which point in it (pointIndex). */
export interface ChartPoint {
  curveNumber: number;
  pointIndex: number;
}

export function PlotlyChart({
  data,
  layout,
  height,
  label,
  download,
  onHover,
  onUnhover,
  onClick,
  onSelectRange,
  onDeselect,
}: {
  data: Data[];
  layout: Partial<Layout>;
  height: number;
  label: string; // describes the chart for screen readers
  // Every chart can be saved: `name` becomes the file name (after the dataset's), `title` the
  // heading on the saved image. Make the title say what the chart shows on its own.
  // `legend`: add a key to the saved image (for charts colored by group); `data`: what to save
  // instead of what's on screen (e.g. without a hover highlight).
  download: { name: string; title: string; legend?: boolean; data?: Data[] };
  onHover?: (point: ChartPoint) => void;
  onUnhover?: () => void;
  onClick?: (point: ChartPoint) => void;
  onSelectRange?: (x0: number, x1: number) => void; // a drag across the chart (layout.dragmode "select")
  onDeselect?: () => void; // a double-click on the chart, clearing a drag selection
}) {
  const ref = useRef<HTMLDivElement>(null);
  const dataset = useContext(DownloadContext);
  // The latest handlers, so Plotly's listeners (attached once) always call the current ones.
  const handlers = useRef({ onHover, onUnhover, onClick, onSelectRange, onDeselect });
  handlers.current = { onHover, onUnhover, onClick, onSelectRange, onDeselect };

  useEffect(() => {
    if (!ref.current) return;
    void Plotly.react(
      ref.current,
      withFallback(data),
      baseLayout(height, layout),
      // Plotly's own camera button is replaced by the Download menu above the chart.
      { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["select2d", "lasso2d", "toggleSpikelines", "toImage"] },
    );
  }, [data, layout, height]);

  // Hover and click events (attached once, after the first draw).
  useEffect(() => {
    const node = ref.current as (HTMLDivElement & { on?: (event: string, fn: (e: unknown) => void) => void }) | null;
    if (!node) return;
    const first = (e: unknown): ChartPoint | null => {
      const p = (e as { points?: { curveNumber: number; pointIndex?: number; pointNumber?: number }[] }).points?.[0];
      return p ? { curveNumber: p.curveNumber, pointIndex: p.pointIndex ?? p.pointNumber ?? 0 } : null;
    };
    const attach = () => {
      node.on?.("plotly_hover", (e) => {
        const p = first(e);
        if (p) handlers.current.onHover?.(p);
      });
      node.on?.("plotly_unhover", () => handlers.current.onUnhover?.());
      node.on?.("plotly_click", (e) => {
        const p = first(e);
        if (p) handlers.current.onClick?.(p);
      });
      node.on?.("plotly_selected", (e) => {
        const x = (e as { range?: { x?: [number, number] } } | undefined)?.range?.x;
        if (x) handlers.current.onSelectRange?.(Math.min(x[0], x[1]), Math.max(x[0], x[1]));
      });
      node.on?.("plotly_deselect", () => handlers.current.onDeselect?.());
    };
    // Plotly adds .on() to the element during the first draw.
    if (node.on) attach();
    else {
      const timer = window.setInterval(() => {
        if (node.on) {
          attach();
          window.clearInterval(timer);
        }
      }, 50);
      return () => window.clearInterval(timer);
    }
  }, []);

  useEffect(() => {
    const node = ref.current;
    return () => {
      if (node) Plotly.purge(node);
    };
  }, []);

  const save = async (format: "png" | "svg") => {
    const today = new Date().toLocaleDateString("en-GB", { day: "numeric", month: "short", year: "numeric" });
    const caption = [dataset.label ?? dataset.file, today].filter(Boolean).join(" · ");
    const width = Math.max(ref.current?.offsetWidth ?? 0, 720);
    const url = await exportChart({
      data: download.data ?? data,
      layout,
      height,
      width: download.legend ? width + 200 : width,
      title: download.title,
      caption,
      format,
      legend: download.legend ?? false,
    });
    saveFile(fileName(dataset.file, download.name, format), url);
  };

  return (
    <div className="flex flex-col gap-1">
      <div className="flex justify-end">
        <DownloadMenu
          what={download.title}
          choices={[
            { label: "PNG image (for slides)", onSelect: () => save("png") },
            { label: "SVG vector (for papers)", onSelect: () => save("svg") },
          ]}
        />
      </div>
      <div ref={ref} role="img" aria-label={label} style={{ width: "100%", height }} />
    </div>
  );
}
