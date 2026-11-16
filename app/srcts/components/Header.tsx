// Not shown for now (analyst request 2026-10-21); App.tsx doesn't render it. The sidebar shows
// the "settings differ" warning instead.
// The logo is left out for now (placeholder until the official Orditus file); components/Logo.tsx is kept.

/** The bar across the top: what is being analyzed, and whether the settings are the
 * validated reference ones. Pages are chosen in the sidebar. */
export function Header({ datasetLine, referenceSettings }: { datasetLine: string; referenceSettings: boolean }) {
  return (
    <header className="flex flex-wrap items-center justify-between gap-4 border-b border-line px-7 py-4">
      <div className="flex min-w-0 items-center gap-[14px]">
        <div className="flex min-w-0 flex-col gap-[2px]">
          <div className="text-[16px] font-bold tracking-[-0.01em]">PBMC3k analysis</div>
          <div className="text-[13px] text-muted">{datasetLine}</div>
        </div>
      </div>

      {referenceSettings ? (
        <span className="inline-flex items-center gap-2 rounded-none border border-line bg-white px-3 py-[6px] text-[13px] font-bold text-[#27272a]">
          <span aria-hidden="true" className="inline-block size-2 rounded-full bg-ok" />
          Validated PBMC3k settings
        </span>
      ) : (
        <span className="inline-flex items-center gap-2 rounded-none border border-accent-line bg-accent-tint px-3 py-[6px] text-[13px] font-bold text-accent-text">
          <span aria-hidden="true" className="text-[11px]">▲</span>
          Settings differ from the validated reference
        </span>
      )}
    </header>
  );
}
