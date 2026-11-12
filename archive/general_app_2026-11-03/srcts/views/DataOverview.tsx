// Step 1: load a dataset (upload a file, or use the demo data) and check it holds raw counts.
import { useShinyInput } from "@posit-dev/shinyreact";
import { Database, Upload } from "lucide-react";
import { useEffect, useRef, useState } from "react";

import { InfoTip } from "../components/InfoTip";
import { StepHeader, SummaryStrip } from "../components/StepPage";
import { Button } from "../components/ui/button";
import { ACCEPTED_FILES, acceptedEnding, uploadToShiny } from "../lib/shinyUpload";
import { goTo, type Page } from "../navigation";
import type { WorkflowPayload } from "../types";
import { ActionBar, fmt } from "./stepShared";

/** What the user asked to load: a file from their computer, or the demo data. */
type Choice = { kind: "file"; file: File } | { kind: "demo" };

/** Where loading stands: uploading (with progress), then the server checking the file. */
type Busy = { what: string; progress: number | null } | null;

export function DataOverview({ workflow, onNavigate }: { workflow: WorkflowPayload; onNavigate: (page: Page) => void }) {
  const step = workflow.steps.find((s) => s.key === "data_overview");
  // Each click sends a new count; the server ignores the initial 0 (see app.py).
  const [demoClicks, setDemoClicks] = useShinyInput<number>("load_demo", 0);
  const fileInput = useRef<HTMLInputElement>(null);
  const [busy, setBusy] = useState<Busy>(null);
  const [problem, setProblem] = useState<string | null>(null); // found in the browser (wrong file type, upload failed)
  const [confirming, setConfirming] = useState<Choice | null>(null); // waiting for "Replace" when data is loaded
  // Only messages that arrive after the page opened (or after the latest load) are shown here.
  const noticeAtStart = useRef<number | null>(workflow.notice?.id ?? null);

  const d = workflow.dataset;
  const loaded = step?.status === "done";

  // The server has answered when the loaded data changes, or a new message arrives.
  useEffect(() => setBusy(null), [d.fingerprint, d.label, workflow.notice?.id]);

  if (!step) return null;
  const notice = workflow.notice;
  const serverError = notice?.step === "data_overview" && notice.id !== noticeAtStart.current ? notice.text : null;
  const error = problem ?? serverError;

  function choose(choice: Choice) {
    setProblem(null);
    if (choice.kind === "file" && !acceptedEnding(choice.file.name)) {
      setProblem(`'${choice.file.name}' isn't a supported format. Use an AnnData (.h5ad) or 10x Genomics (.h5) file.`);
      return;
    }
    if (loaded) setConfirming(choice);
    else void load(choice);
  }

  async function load(choice: Choice) {
    setConfirming(null);
    noticeAtStart.current = workflow.notice?.id ?? null;
    if (choice.kind === "demo") {
      setBusy({ what: "Loading the demo data…", progress: null });
      setDemoClicks(demoClicks + 1);
      return;
    }
    setBusy({ what: `Uploading ${choice.file.name}…`, progress: 0 });
    try {
      await uploadToShiny(choice.file, "upload_data", (fraction) =>
        setBusy({ what: fraction < 1 ? `Uploading ${choice.file.name}…` : `Checking ${choice.file.name}…`, progress: fraction }),
      );
    } catch (e) {
      setBusy(null);
      setProblem(e instanceof Error ? e.message : String(e));
    }
  }

  return (
    <div className="flex flex-col gap-6 px-8 pt-[30px] pb-10">
      <StepHeader step={step}>Load a single-cell dataset and check that it holds raw counts.</StepHeader>

      {loaded && d.n_cells !== null && d.n_genes !== null && (
        <SummaryStrip
          items={[
            { label: "Cells", value: fmt(d.n_cells), note: "barcodes in the file" },
            { label: "Genes", value: fmt(d.n_genes), note: "measured per cell" },
            { label: "Counts", value: "Raw", note: "whole numbers, not normalized" },
          ]}
        />
      )}

      <section aria-label="Upload data" className="flex max-w-[760px] flex-col gap-3">
        <h2 className="m-0 border-b border-ink pb-[10px] text-lg font-extrabold">Upload data</h2>
        <div className="flex flex-wrap items-center gap-x-4 gap-y-3 pt-1">
          <div className="flex flex-col gap-1">
            <Button variant={loaded ? "secondary" : "primary"} disabled={busy !== null} onClick={() => fileInput.current?.click()}>
              <Upload aria-hidden="true" size={15} /> Choose a file…
            </Button>
            <span className="text-[12.5px] text-muted">.h5ad or 10x .h5</span>
          </div>
          <span className="pb-5 text-[13px] text-muted">or</span>
          <div className="flex items-center gap-2 pb-5">
            <Button variant="secondary" disabled={busy !== null} onClick={() => choose({ kind: "demo" })}>
              <Database aria-hidden="true" size={15} /> Use demo data
            </Button>
            <InfoTip label="About the demo data">
              Demo data: human blood cells (PBMC3k, 10x Genomics). About 3,000 peripheral blood mononuclear cells from
              a healthy donor. Use it to try the app before uploading your own data.
            </InfoTip>
          </div>
          <input
            ref={fileInput}
            type="file"
            accept={ACCEPTED_FILES.join(",")}
            className="hidden"
            aria-label="Choose a data file"
            onChange={(e) => {
              const file = e.target.files?.[0];
              e.target.value = ""; // so choosing the same file again still counts
              if (file) choose({ kind: "file", file });
            }}
          />
        </div>

        {confirming && (
          <div role="alertdialog" aria-label="Replace the current analysis" className="flex flex-wrap items-center gap-3 rounded-[4px] border border-ink bg-panel px-4 py-3 text-[14px]">
            <span className="min-w-0 flex-[1_1_320px]">
              <b>Loading new data replaces the current analysis.</b>{" "}
              {confirming.kind === "file" ? `Load ${confirming.file.name}?` : "Load the demo data?"}
            </span>
            <Button size="small" onClick={() => void load(confirming)}>
              Replace
            </Button>
            <Button size="small" variant="secondary" onClick={() => setConfirming(null)}>
              Cancel
            </Button>
          </div>
        )}

        {busy && (
          <div role="status" className="flex flex-col gap-[6px] text-[13.5px] text-body">
            <span>{busy.what}</span>
            {busy.progress !== null && (
              <span className="h-[6px] w-full max-w-[420px] overflow-hidden rounded-[2px] bg-chip" aria-hidden="true">
                <span className="block h-full bg-ink transition-[width]" style={{ width: `${Math.round(busy.progress * 100)}%` }} />
              </span>
            )}
          </div>
        )}
      </section>

      {loaded && (
        <section aria-label="Loaded data" className="flex max-w-[760px] flex-col gap-3">
          <h2 className="m-0 border-b border-ink pb-[10px] text-lg font-extrabold">Loaded data</h2>
          <dl className="m-0 grid grid-cols-[max-content_minmax(0,1fr)] gap-x-8 gap-y-2 text-[14px]">
            <dt className="text-muted">Dataset</dt>
            <dd className="m-0 font-bold">{d.label}</dd>
            {d.label !== d.file && (
              <>
                <dt className="text-muted">File</dt>
                <dd className="m-0 font-mono">{d.file}</dd>
              </>
            )}
            <dt className="text-muted">Format</dt>
            <dd className="m-0">{d.file_format}</dd>
            <dt className="text-muted">Checks</dt>
            <dd className="m-0">counts are whole numbers ≥ 0 · cell and gene names are unique</dd>
            <dt className="text-muted">Mitochondrial genes</dt>
            <dd className="m-0">
              {d.mito_matches.length > 0 ? (
                <>
                  <b>{fmt(d.mito_matches[0].n_genes)} found</b>, named <span className="font-mono">{d.mito_matches[0].prefix}…</span>{" "}
                  <span className="text-muted">(as in {d.mito_matches[0].species})</span>
                </>
              ) : (
                <span className="font-bold text-accent-text">None found with the usual name patterns. You can set the pattern on Quality Control.</span>
              )}
            </dd>
            <dt className="text-muted">Fingerprint</dt>
            <dd className="m-0">
              <span className="font-mono">{d.fingerprint}</span> <span className="text-muted">(SHA-256, recorded with every run)</span>
            </dd>
          </dl>
          {d.notes.map((note) => (
            <p key={note} className="m-0 text-[13.5px] font-bold text-accent-text">
              ▲ {note}
            </p>
          ))}
        </section>
      )}

      <ActionBar
        label="Load data"
        status={
          error ? (
            <span className="font-bold text-accent-text">{error}</span>
          ) : busy ? (
            <span className="text-body">{busy.what}</span>
          ) : step.confirmation ? (
            <span className="text-body">{step.confirmation}</span>
          ) : (
            <span className="text-body">No data loaded yet. Upload a file or use the demo data.</span>
          )
        }
        upToDate={loaded}
        nextTitle="Quality Control"
        onContinue={() => goTo(onNavigate, "quality_control")}
      />
    </div>
  );
}
