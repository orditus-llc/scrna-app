// Step 1: open the data file and check it holds raw counts.
import { useShinyInput } from "@posit-dev/shinyreact";
import { useEffect, useState } from "react";

import { StepHeader, SummaryStrip } from "../components/StepPage";
import { goTo, type Page } from "../navigation";
import type { WorkflowPayload } from "../types";
import { ActionBar } from "./stepShared";

export function DataOverview({ workflow, onNavigate }: { workflow: WorkflowPayload; onNavigate: (page: Page) => void }) {
  const step = workflow.steps.find((s) => s.key === "data_overview");
  // Each click sends a new count; the server ignores the initial 0 (see app.py).
  const [loadClicks, setLoadClicks] = useShinyInput<number>("load_data", 0);
  const [loading, setLoading] = useState(false);
  // Stop "Loading…" when the server answers (the data, or a message).
  useEffect(() => setLoading(false), [step?.status, workflow.notice?.id]);

  if (!step) return null;
  const d = workflow.dataset;
  const error = workflow.notice?.step === "data_overview" && step.status !== "done" ? workflow.notice.text : null;

  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <StepHeader step={step}>
        Opens the data file and checks that it holds raw counts.
      </StepHeader>

      {step.status === "done" && d.n_cells !== null && d.n_genes !== null && (
        <SummaryStrip
          items={[
            { label: "Cells", value: fmt(d.n_cells), note: "barcodes in the file" },
            { label: "Genes", value: fmt(d.n_genes), note: "measured per cell" },
            { label: "Counts", value: "Raw", note: "whole numbers, not normalized" },
            { label: "File fingerprint", value: d.fingerprint ?? "–", note: "SHA-256, recorded with every run" },
          ]}
        />
      )}

      <section className="flex max-w-[760px] flex-col gap-3">
        <h2 className="m-0 border-b border-ink pb-[10px] text-lg font-extrabold">Input</h2>
        <dl className="m-0 grid grid-cols-[max-content_minmax(0,1fr)] gap-x-8 gap-y-2 text-[14px]">
          <dt className="text-muted">File</dt>
          <dd className="m-0 font-mono">{d.file}</dd>
          <dt className="text-muted">Format</dt>
          <dd className="m-0">AnnData (.h5ad)</dd>
          <dt className="text-muted">Checks</dt>
          <dd className="m-0">counts are whole numbers ≥ 0 · cell and gene names are unique</dd>
        </dl>
      </section>

      <ActionBar
        label="Load the data"
        status={
          error ? (
            <span className="font-bold text-accent-text">{error}</span>
          ) : loading ? (
            <span className="text-body">Loading…</span>
          ) : step.confirmation ? (
            <span className="text-body">{step.confirmation}</span>
          ) : (
            <span className="text-body">Not loaded yet.</span>
          )
        }
        runLabel={loading ? "Loading…" : "Load PBMC3k"}
        runDisabled={loading || step.status === "done"}
        upToDate={step.status === "done"}
        nextTitle="Quality Control"
        onRun={() => {
          setLoading(true);
          setLoadClicks(loadClicks + 1);
        }}
        onContinue={() => goTo(onNavigate, "quality_control")}
      />
    </div>
  );
}

function fmt(value: number): string {
  return value.toLocaleString("en-US");
}
