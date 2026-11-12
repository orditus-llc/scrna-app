// A message across the top when the connection to the analysis server is lost (server
// stopped, laptop slept, network dropped). Shiny greys the page out with an overlay element
// (#shiny-disconnected-overlay) when that happens; we watch for that element and explain what
// it means. The analysis lives in the server's memory, so the results of this session are gone.
import { RotateCw, Unplug } from "lucide-react";
import { useEffect, useState } from "react";

import { Button } from "./ui/button";

const OVERLAY_ID = "shiny-disconnected-overlay";

export function ConnectionLost() {
  const [lost, setLost] = useState(false);

  useEffect(() => {
    const check = () => setLost(document.getElementById(OVERLAY_ID) !== null);
    check();
    const observer = new MutationObserver(check);
    observer.observe(document.body, { childList: true });
    return () => observer.disconnect();
  }, []);

  if (!lost) return null;
  return (
    <div
      role="alert"
      className="fixed inset-x-0 top-0 z-[100000] flex flex-wrap items-center justify-center gap-x-4 gap-y-2 border-b border-ink bg-white px-6 py-3 text-[14px] text-ink shadow-[0_2px_10px_rgba(0,0,0,0.12)]"
    >
      <Unplug aria-hidden="true" size={18} className="flex-none text-accent-text" />
      <span>
        <b>Lost connection to the analysis server.</b> This session's results are gone. Reload to start again.
      </span>
      <Button size="small" onClick={() => window.location.reload()}>
        <RotateCw aria-hidden="true" size={14} /> Reload
      </Button>
    </div>
  );
}
