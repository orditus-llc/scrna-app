// "Run the remaining steps": one click runs every step that isn't up to date, in order, with
// the settings in use. The server runs one step at a time and reports progress (app.py,
// _work_through_remaining). The state lives in App and reaches every page's bottom bar
// (`ActionBar`) through this context, so the pages don't change.
import { useShinyInput } from "@posit-dev/shinyreact";
import { createContext, useEffect, useRef, useState } from "react";

import type { Page } from "../navigation";
import type { WorkflowPayload } from "../types";

/** How the last job ended, shown in the bottom bar of the page it was started from. */
export interface RunOutcome {
  kind: "info" | "error";
  text: string; // e.g. "Ran 5 steps. Choose the species and tissues for Cell Type Annotation, then run it."
  noticeId: number | null; // the server message with the same text (not shown again as a pop-up)
  open: { page: Page; title: string }; // where to go next, if the user wants to
}

export interface RunRemaining {
  remaining: number; // steps not up to date
  running: boolean;
  progress: string; // e.g. "Running Clustering · step 3 of 5…"
  outcome: RunOutcome | null; // the last job's result, on the page it was started from
  start: () => void;
}

export const RunRemainingContext = createContext<RunRemaining | null>(null);

/** The state behind the button. The page stays where it is (analyst decision 2026-10-29): when
 * a job ends, its bottom bar says how it went and offers "Open …" for the step that needs
 * attention (failed, or waiting for a choice) or Explore Results. Leaving the page clears it. */
export function useRunRemaining(workflow: WorkflowPayload | undefined, page: Page): RunRemaining {
  const [clicks, setClicks] = useShinyInput<number>("run_remaining", 0);
  const job = workflow?.run_remaining ?? null;
  const [waitingAfter, setWaitingAfter] = useState<number | null>(null); // job id before our click
  const [startedOn, setStartedOn] = useState<Page | null>(null);
  const [outcome, setOutcome] = useState<RunOutcome | null>(null);
  const handled = useRef<number | null>(null); // the last finished job we reported

  // The server's answer (a new job) ends the wait for it to start.
  useEffect(() => {
    if (waitingAfter !== null && job && job.id > waitingAfter) setWaitingAfter(null);
  }, [job?.id]); // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (!workflow || !job || job.status === "running" || handled.current === job.id) return;
    handled.current = job.id;
    const titleOf = (key: string) => workflow.steps.find((s) => s.key === key)?.title ?? key;
    const notice = workflow.notice;
    const target = job.status === "failed" ? job.failed_step : job.status === "stopped" ? job.stopped_step : null;
    setOutcome({
      kind: job.status === "failed" ? "error" : "info",
      text: notice?.text ?? "",
      noticeId: notice?.id ?? null,
      open: target ? { page: target as Page, title: titleOf(target) } : { page: "explorer", title: "Explore Results" },
    });
  }, [job?.id, job?.status]); // eslint-disable-line react-hooks/exhaustive-deps

  // The outcome belongs to the page the button was pressed on; opening another page clears it.
  useEffect(() => {
    if (startedOn !== null && page !== startedOn) {
      setOutcome(null);
      setStartedOn(null);
    }
  }, [page]); // eslint-disable-line react-hooks/exhaustive-deps

  const running = waitingAfter !== null || job?.status === "running";
  return {
    remaining: workflow?.remaining_steps ?? 0,
    running,
    progress:
      job?.status === "running" && job.current
        ? `Running ${job.current} · step ${job.index} of ${job.total}…`
        : "Starting the remaining steps…",
    outcome: !running && startedOn === page ? outcome : null,
    start() {
      setOutcome(null);
      setStartedOn(page);
      setWaitingAfter(job?.id ?? 0);
      setClicks(clicks + 1);
    },
  };
}
