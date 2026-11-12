// A section that starts folded ("More settings ▸") so a page shows only its main settings.
import { ChevronDown } from "lucide-react";
import { useId, useState, type ReactNode } from "react";

import { cn } from "../lib/utils";

export function MoreSection({
  title,
  note,
  defaultOpen = false,
  children,
}: {
  title: string;
  note?: string; // shown next to the title while folded, e.g. "changed"
  defaultOpen?: boolean;
  children: ReactNode;
}) {
  const id = useId();
  const [open, setOpen] = useState(defaultOpen);
  return (
    <div className="flex flex-col border-b border-line">
      <button
        type="button"
        aria-expanded={open}
        aria-controls={id}
        onClick={() => setOpen((o) => !o)}
        className="flex items-center gap-2 py-3 text-left text-[14px] font-bold text-ink hover:text-body"
      >
        <ChevronDown aria-hidden="true" size={16} className={cn("transition-transform", !open && "-rotate-90")} />
        {title}
        {note && !open && <span className="text-[12.5px] font-bold text-accent-text">{note}</span>}
      </button>
      {open && (
        <div id={id} className="flex flex-col pb-2">
          {children}
        </div>
      )}
    </div>
  );
}
