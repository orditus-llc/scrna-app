// A small wrapper that draws a Plotly chart and updates it in place when its data changes.
import type { Data, Layout } from "plotly.js-dist-min";
import Plotly from "plotly.js-dist-min";
import { useEffect, useRef } from "react";

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
  onHover?: (point: ChartPoint) => void;
  onUnhover?: () => void;
  onClick?: (point: ChartPoint) => void;
  onSelectRange?: (x0: number, x1: number) => void; // a drag across the chart (layout.dragmode "select")
  onDeselect?: () => void; // a double-click on the chart, clearing a drag selection
}) {
  const ref = useRef<HTMLDivElement>(null);
  // The latest handlers, so Plotly's listeners (attached once) always call the current ones.
  const handlers = useRef({ onHover, onUnhover, onClick, onSelectRange, onDeselect });
  handlers.current = { onHover, onUnhover, onClick, onSelectRange, onDeselect };

  useEffect(() => {
    if (!ref.current) return;
    void Plotly.react(
      ref.current,
      withFallback(data),
      {
        height,
        margin: { l: 56, r: 16, t: 36, b: 48 },
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        font: PLOT_FONT,
        showlegend: false,
        ...layout,
      },
      { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["select2d", "lasso2d", "toggleSpikelines"] },
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

  return <div ref={ref} role="img" aria-label={label} style={{ width: "100%", height }} />;
}
