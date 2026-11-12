// A small ⓘ button that shows a short explanation. Opens on hover, on keyboard focus and on
// tap (click toggles); Esc or moving away closes it. Built here (no pop-up library): the
// text sits in a box next to the icon, so it needs no portal.
import { Info } from "lucide-react";
import { useEffect, useId, useState, type ReactNode } from "react";

import { cn } from "../lib/utils";

export function InfoTip({
  label,
  children,
  align = "left",
}: {
  label: string; // what it explains, for screen readers: "About Max genes per cell"
  children: ReactNode;
  align?: "left" | "right"; // which edge of the icon the box lines up with
}) {
  const id = useId();
  const [hover, setHover] = useState(false);
  const [pinned, setPinned] = useState(false); // opened by a click or tap
  const open = hover || pinned;

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") {
        setHover(false);
        setPinned(false);
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open]);

  return (
    <span className="relative inline-flex align-middle" onMouseEnter={() => setHover(true)} onMouseLeave={() => setHover(false)}>
      <button
        type="button"
        aria-label={label}
        aria-describedby={open ? id : undefined}
        aria-expanded={open}
        onFocus={() => setHover(true)}
        onBlur={() => {
          setHover(false);
          setPinned(false);
        }}
        onClick={() => setPinned((p) => !p)}
        className={cn(
          "inline-flex size-[18px] items-center justify-center rounded-full text-muted hover:text-ink",
          open && "text-ink",
        )}
      >
        <Info aria-hidden="true" size={15} strokeWidth={2.2} />
      </button>
      {open && (
        <span
          id={id}
          role="tooltip"
          className={cn(
            "absolute top-[calc(100%+6px)] z-50 w-[260px] max-w-[70vw] rounded-[4px] border border-line bg-white px-3 py-2 text-left text-[13px] leading-snug font-normal text-body shadow-[0_2px_10px_rgba(0,0,0,0.12)]",
            align === "left" ? "left-0" : "right-0",
          )}
        >
          {children}
        </span>
      )}
    </span>
  );
}
