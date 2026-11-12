// "Download" next to a chart or table: one button, or a short menu when there are several
// formats (PNG / SVG). Built here without a pop-up library: the menu sits under the button,
// closes on Esc, on a click elsewhere, or after a choice.
import { Download } from "lucide-react";
import { createContext, useContext, useEffect, useId, useRef, useState } from "react";

import { fileName, saveCsv, toCsv, type CsvValue } from "../lib/download";
import { cn } from "../lib/utils";

/** The loaded data file's name (for file names) and how it's named on screen (for image
 * captions). Provided once by App. */
export const DownloadContext = createContext<{ file: string | null; label: string | null }>({ file: null, label: null });

export interface DownloadChoice {
  label: string; // e.g. "PNG image (high resolution)"
  onSelect: () => void | Promise<void>;
}

export function DownloadMenu({ what, choices }: { what: string; choices: DownloadChoice[] }) {
  const [open, setOpen] = useState(false);
  const [busy, setBusy] = useState(false);
  const box = useRef<HTMLDivElement>(null);
  const id = useId();

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && setOpen(false);
    const onClick = (e: MouseEvent) => !box.current?.contains(e.target as Node) && setOpen(false);
    document.addEventListener("keydown", onKey);
    document.addEventListener("mousedown", onClick);
    return () => {
      document.removeEventListener("keydown", onKey);
      document.removeEventListener("mousedown", onClick);
    };
  }, [open]);

  const run = async (choice: DownloadChoice) => {
    setOpen(false);
    setBusy(true);
    try {
      await choice.onSelect();
    } finally {
      setBusy(false);
    }
  };

  const button =
    "inline-flex items-center gap-[6px] rounded-[4px] border border-field bg-white px-[9px] py-[5px] text-[12.5px] font-bold text-ink enabled:hover:bg-chip disabled:text-muted";
  if (choices.length === 1) {
    return (
      <button type="button" className={button} disabled={busy} onClick={() => run(choices[0])} aria-label={`Download ${what} (${choices[0].label})`}>
        <Download aria-hidden="true" size={14} />
        {busy ? "Preparing…" : choices[0].label}
      </button>
    );
  }
  return (
    <div ref={box} className="relative">
      <button
        type="button"
        className={button}
        disabled={busy}
        aria-haspopup="menu"
        aria-expanded={open}
        aria-controls={id}
        aria-label={`Download ${what}`}
        onClick={() => setOpen((o) => !o)}
      >
        <Download aria-hidden="true" size={14} />
        {busy ? "Preparing…" : "Download"}
      </button>
      {open && (
        <ul
          id={id}
          role="menu"
          className="absolute right-0 z-30 m-0 mt-1 flex min-w-[220px] list-none flex-col rounded-[4px] border border-line bg-white p-1 shadow-[0_2px_8px_rgba(0,0,0,0.12)]"
        >
          {choices.map((c) => (
            <li key={c.label} role="none">
              <button
                type="button"
                role="menuitem"
                onClick={() => run(c)}
                className={cn("w-full rounded-[3px] px-3 py-2 text-left text-[13px] hover:bg-chip focus-visible:bg-chip")}
              >
                {c.label}
              </button>
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}

/** "CSV" next to a table: saves `rows` (all of them, not only the ones shown) as
 * "<dataset>_<name>.csv". */
export function CsvDownload({
  name,
  columns,
  rows,
}: {
  name: string; // what the table holds, e.g. "Marker Genes top markers"
  columns: { key: string; header: string }[];
  rows: () => Record<string, CsvValue>[]; // built only when clicked
}) {
  const dataset = useContext(DownloadContext);
  return (
    <DownloadMenu
      what={name}
      choices={[{ label: "CSV", onSelect: () => saveCsv(fileName(dataset.file, name, "csv"), toCsv(columns, rows())) }]}
    />
  );
}
