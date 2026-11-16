import { X } from "lucide-react";
import { useEffect, useRef, useState } from "react";

import type { Notice } from "../types";

/** Shows the server's latest message (e.g. an error) in the bottom corner for 12 seconds.
 * `shownElsewhere`: a message the page already shows (the outcome of "Run the remaining steps"
 * in the bottom bar); it is hidden, and stays hidden once seen. */
export function Toast({ notice, shownElsewhere }: { notice: Notice | null; shownElsewhere: number | null }) {
  const [visible, setVisible] = useState<Notice | null>(null);
  const hidden = useRef(new Set<number>());
  if (shownElsewhere !== null) hidden.current.add(shownElsewhere);

  useEffect(() => {
    if (!notice) return;
    setVisible(notice);
    const timer = window.setTimeout(() => setVisible(null), 12_000);
    return () => window.clearTimeout(timer);
  }, [notice?.id]); // eslint-disable-line react-hooks/exhaustive-deps -- a new id = a new message

  if (!visible || hidden.current.has(visible.id)) return null;
  return (
    <div
      role={visible.kind === "error" ? "alert" : "status"}
      className="fixed right-6 bottom-6 z-50 flex max-w-[420px] items-start gap-3 rounded-none bg-ink px-4 py-3 text-sm text-white shadow-[0_2px_8px_rgba(0,0,0,0.18)]"
    >
      <span
        aria-hidden="true"
        className={`mt-[6px] inline-block size-2 flex-none rounded-full ${visible.kind === "error" ? "bg-accent" : "bg-ok"}`}
      />
      <p className="m-0 flex-1 leading-snug">{visible.text}</p>
      <button type="button" aria-label="Close message" onClick={() => setVisible(null)} className="text-white/80 hover:text-white">
        <X size={16} />
      </button>
    </div>
  );
}
