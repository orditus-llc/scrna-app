// "Run the remaining steps": one click runs every step that isn't up to date, in order, with
// the settings in use (the validated ones unless a step was changed). The server runs one step
// at a time and reports progress (app.py, _work_through_remaining). The state lives in App and
// reaches every page's bottom bar (`ActionBar`) through this context, so the pages don't change.
import { useShinyInput } from "@posit-dev/shinyreact";
import { createContext, useEffect, useRef, useState } from "react";

import type { Page } from "../navigation";
import type { WorkflowPayload } from "../types";

export interface RunRemaining {
  remaining: number; // steps not up to date
  running: boolean;
  progress: string; // e.g. "Running Clustering · step 3 of 5…"
  start: () => void;
}

export const RunRemainingContext = createContext<RunRemaining | null>(null);

/** The state behind the button. When a job ends it opens Explore Results, or the page of the
 * step that failed (its error is shown in the pop-up). */
export function useRunRemaining(workflow: WorkflowPayload | undefined, navigate: (page: Page) => void): RunRemaining {
  const [clicks, setClicks] = useShinyInput<number>("run_remaining", 0);
  const job = workflow?.run_remaining ?? null;
  const [waitingAfter, setWaitingAfter] = useState<number | null>(null); // job id before our click
  const handled = useRef<number | null>(null); // the last finished job we navigated for

  // The server's answer (a new job) ends the wait for it to start.
  useEffect(() => {
    if (waitingAfter !== null && job && job.id > waitingAfter) setWaitingAfter(null);
  }, [job?.id]); // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (!job || job.status === "running" || handled.current === job.id) return;
    handled.current = job.id;
    navigate(job.status === "failed" && job.failed_step ? (job.failed_step as Page) : "explorer");
  }, [job?.id, job?.status]); // eslint-disable-line react-hooks/exhaustive-deps

  const running = waitingAfter !== null || job?.status === "running";
  return {
    remaining: workflow?.remaining_steps ?? 0,
    running,
    progress:
      job?.status === "running" && job.current
        ? `Running ${job.current} · step ${job.index} of ${job.total}…`
        : "Starting the remaining steps…",
    start() {
      setWaitingAfter(job?.id ?? 0);
      setClicks(clicks + 1);
    },
  };
}
