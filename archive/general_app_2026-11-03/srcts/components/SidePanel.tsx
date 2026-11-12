// A panel that slides in from the right over the page, closed with its Close button or Esc.
// Used for a cluster's evidence (Cell Type Annotation) and the Run record (sidebar).
import { X } from "lucide-react";
import { useEffect, type ReactNode } from "react";

import { Button } from "./ui/button";

/** A panel that slides in from the right over the page (Close or Esc to close). */
export function SidePanel({ title, onClose, children }: { title: string; onClose: () => void; children: ReactNode }) {
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && onClose();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [onClose]);
  return (
    <div
      role="dialog"
      aria-label={title}
      className="fixed top-0 right-0 z-50 flex h-full w-[min(620px,100vw)] flex-col border-l border-ink bg-white shadow-[-8px_0_24px_rgba(0,0,0,0.12)]"
    >
      <div className="flex items-center justify-between gap-3 border-b border-line px-5 py-3">
        <h2 className="m-0 text-[16px] font-extrabold">{title}</h2>
        <Button variant="secondary" size="small" onClick={onClose} autoFocus>
          <X aria-hidden="true" size={14} /> Close
        </Button>
      </div>
      <div className="flex-auto overflow-y-auto px-5 py-4">{children}</div>
    </div>
  );
}
