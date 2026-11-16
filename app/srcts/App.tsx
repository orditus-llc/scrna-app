import { useShinyOutputValue } from "@posit-dev/shinyreact";
import { useCallback, useEffect, useState } from "react";

import { BusyContext, BusyOverlay } from "./components/BusyOverlay";
import { DownloadContext } from "./components/DownloadMenu";
import { ExportButtons } from "./components/ExportButtons";
import { ConnectionLost } from "./components/ConnectionLost";
import { RunRemainingContext, useRunRemaining } from "./components/RunRemaining";
import { SidePanel } from "./components/SidePanel";
import { RunLog, SettingsRecord } from "./components/SidePanels";
import { Sidebar } from "./components/Sidebar";
import { Toast } from "./components/Toast";
import { BUILT_STEPS, pageFromHash, type Page, type StepKey } from "./navigation";
import type { AnnotationPayload, ConditionsPayload, MitoCheck, QcSuggestion, ClusterPayload, ExplorePayload, GenePayload, MarkersPayload, PreprocessPayload, PreviewPayload, QcPayload, ReductionPayload, WorkflowPayload } from "./types";
import { AnnotationPage } from "./views/AnnotationPage";
import { ClusteringPage } from "./views/ClusteringPage";
import { ConditionsPage } from "./views/ConditionsPage";
import { DataOverview } from "./views/DataOverview";
import { ExplorePage } from "./views/ExplorePage";
import { MarkersPage } from "./views/MarkersPage";
import { QualityControl } from "./views/QualityControl";
import { ComingSoonPage, PreprocessingPage, ReductionPage } from "./views/StepPages";

export function App() {
  // Every server output is subscribed here, once, for the life of the page. Subscribing
  // inside a page would re-register the output each time the page reopens, and
  // shinyreact 0.1.0 doesn't release the old one cleanly ("Duplicate output IDs").
  const data = useShinyOutputValue<WorkflowPayload>("workflow");
  const qc = useShinyOutputValue<QcPayload>("qc");
  const qcPreview = useShinyOutputValue<PreviewPayload | null>("qc_preview");
  const mitoCheck = useShinyOutputValue<MitoCheck | null>("mito_check");
  const qcSuggestion = useShinyOutputValue<QcSuggestion>("qc_suggestion");
  const preprocessing = useShinyOutputValue<PreprocessPayload>("preprocessing");
  const reduction = useShinyOutputValue<ReductionPayload>("reduction");
  const clustering = useShinyOutputValue<ClusterPayload>("clustering");
  const markers = useShinyOutputValue<MarkersPayload>("markers");
  const annotation = useShinyOutputValue<AnnotationPayload>("annotation");
  const explore = useShinyOutputValue<ExplorePayload>("explore");
  const exploreGene = useShinyOutputValue<GenePayload | null>("explore_gene_values");
  const conditions = useShinyOutputValue<ConditionsPayload>("conditions");
  // The app always opens on the first step: a new browser session starts a new, empty
  // analysis on the server, so a later page would have nothing to show.
  const [page, setPage] = useState<Page>("data_overview");
  const [recordOpen, setRecordOpen] = useState(false);

  // Each page has its own address (#quality_control, …), so Back and Forward move between
  // pages. Opening a page sets the address; the address change then shows the page.
  useEffect(() => {
    window.history.replaceState(null, "", "#data_overview");
    const onHashChange = () => {
      setPage(pageFromHash(window.location.hash) ?? "data_overview");
      window.scrollTo({ top: 0 });
    };
    window.addEventListener("hashchange", onHashChange);
    return () => window.removeEventListener("hashchange", onHashChange);
  }, []);
  function navigate(next: Page) {
    if (window.location.hash === `#${next}`) setPage(next);
    else window.location.hash = next;
  }

  const runRemaining = useRunRemaining(data, page);
  // What is running right now, reported by the pages (components/BusyOverlay.tsx): the screen greys out.
  const [busyLabels, setBusyLabels] = useState<Record<string, string>>({});
  const setBusyLabel = useCallback((key: string, label: string | null) => {
    setBusyLabels((old) => {
      if ((old[key] ?? null) === label) return old;
      const next = { ...old };
      if (label === null) delete next[key];
      else next[key] = label;
      return next;
    });
  }, []);

  // Results live in the server's memory for this browser tab; reloading or closing it starts
  // over. Once anything has been run, the browser asks before leaving ("Leave site?").
  const hasResults = (data?.run_log.length ?? 0) > 0;
  useEffect(() => {
    if (!hasResults) return;
    const warn = (event: BeforeUnloadEvent) => {
      event.preventDefault();
      event.returnValue = ""; // older browsers need this to show the warning
    };
    window.addEventListener("beforeunload", warn);
    return () => window.removeEventListener("beforeunload", warn);
  }, [hasResults]);

  function content() {
    if (!data) return <p className="m-0 px-8 py-10 text-muted">Connecting to the analysis server…</p>;
    switch (page) {
      case "data_overview":
        return <DataOverview workflow={data} onNavigate={navigate} />;
      case "quality_control":
        return (
          <QualityControl
            qc={qc}
            preview={qcPreview}
            mitoCheck={mitoCheck ?? null}
            qcSuggestion={qcSuggestion ?? null}
            step={data.steps.find((s) => s.key === "quality_control")}
            notice={data.notice}
            onNavigate={navigate}
          />
        );
      case "preprocessing":
        return <PreprocessingPage workflow={data} data={preprocessing} onNavigate={navigate} />;
      case "dimensionality_reduction":
        return <ReductionPage workflow={data} data={reduction} onNavigate={navigate} />;
      case "clustering":
        return <ClusteringPage workflow={data} data={clustering} onNavigate={navigate} />;
      case "marker_genes":
        return <MarkersPage workflow={data} data={markers} onNavigate={navigate} />;
      case "cell_type_annotation":
        return <AnnotationPage workflow={data} data={annotation} onNavigate={navigate} />;
      case "explorer":
        return <ExplorePage data={explore} gene={exploreGene} exportInfo={data?.export} onNavigate={navigate} />;
      case "compare_conditions":
        return <ConditionsPage workflow={data} data={conditions} onNavigate={navigate} />;
      default:
        return BUILT_STEPS[page as StepKey] === false ? (
          <ComingSoonPage workflow={data} stepKey={page as StepKey} onNavigate={navigate} />
        ) : null;
    }
  }

  return (
    <BusyContext.Provider value={setBusyLabel}>
    <RunRemainingContext.Provider value={runRemaining}>
      <DownloadContext.Provider value={{ file: data?.dataset.file ?? null, label: data?.dataset.label ?? null }}>
      <div className="page-ground min-h-screen font-sans text-ink">
        <div className="mx-auto min-h-screen max-w-[1440px] border-x border-line bg-white">
          {/* The top bar (components/Header.tsx) is left out for now (analyst request); only its
              "settings differ" warning is kept, at the top of the sidebar. */}
          {/* The sidebar sits above the bottom bar (z-40) so its ⓘ tooltips are never covered, and below
              the side panels and pop-ups (z-50). */}
          <div className="flex flex-wrap items-stretch md:min-h-screen">
            <aside className="relative z-[45] min-w-0 flex-[0_0_264px] p-3 max-md:flex-[1_1_100%]">
              {/* A panel like VS Code's Explorer: light grey, thin border, rounded, set in from the edges. */}
              <div className="rounded-none border border-field bg-side md:sticky md:top-3 md:min-h-[calc(100vh-24px)]">
              {data && (
                <Sidebar
                  steps={data.steps}
                  nextKey={data.next_step?.key ?? null}
                  page={page}
                  exploreAvailable={explore?.available ?? false}
                  conditions={data.conditions}
                  onNavigate={navigate}
                  onOpenRecord={() => setRecordOpen(true)}
                  runCount={data.run_log.length}
                  />
              )}
              </div>
            </aside>
            {/* key: each page starts fresh (its settings boxes are re-read from the server) */}
            <main key={page} className="min-w-0 flex-[999_1_720px] pb-12">
              {content()}
            </main>
          </div>
        </div>
        {recordOpen && data && (
          <SidePanel title="Run record" onClose={() => setRecordOpen(false)}>
            <div className="flex flex-col gap-4">
              <ExportButtons info={data.export} />
              <SettingsRecord items={data.settings_record} versions={data.versions} />
              <RunLog entries={data.run_log} />
            </div>
          </SidePanel>
        )}
        {/* A finished "Run the remaining steps" says how it went in the page's bottom bar instead. */}
        <Toast notice={data?.notice ?? null} shownElsewhere={runRemaining.outcome?.noticeId ?? null} />
        <ConnectionLost />
        <BusyOverlay label={runRemaining.running ? runRemaining.progress : (Object.values(busyLabels)[0] ?? null)} />
      </div>
      </DownloadContext.Provider>
    </RunRemainingContext.Provider>
    </BusyContext.Provider>
  );
}
