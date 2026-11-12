// The "Mitochondrial genes" setting on the Quality Control page: which gene-name prefix marks
// mitochondrial genes (human "MT-", mouse/zebrafish "mt-", rat "Mt-", fly "mt:", or typed in).
// It starts at the pattern found in the data, and a live count shows how many genes match.
import type { MitoCheck, QcSettings } from "../types";
import { cn } from "../lib/utils";
import { InfoTip } from "./InfoTip";

const OTHER = "__other__";

export function MitoSetting({
  settings,
  value,
  onChange,
  check,
}: {
  settings: QcSettings;
  value: string;
  onChange: (prefix: string) => void;
  check: MitoCheck | null; // the live count for `value` (null until the server answers)
}) {
  const known = settings.patterns.some((p) => p.prefix === value);
  const counted = check && check.prefix === value ? check : null;
  return (
    <div className="flex flex-col gap-2 border-b border-line py-[14px]">
      <span className="text-[14px] leading-snug">
        <label htmlFor="mito_pattern" className="font-bold">
          Mitochondrial genes
        </label>{" "}
        <InfoTip label="About Mitochondrial genes">
          How mitochondrial gene names start in your data. Damaged cells have a high share of counts from these genes.
          Set from your data; choose another pattern or type your own.
        </InfoTip>
      </span>
      <select
        id="mito_pattern"
        value={known ? value : OTHER}
        onChange={(e) => onChange(e.target.value === OTHER ? "" : e.target.value)}
        className="rounded-[4px] border border-field bg-white px-[10px] py-2 text-[14px]"
      >
        {settings.patterns.map((p) => (
          <option key={p.prefix} value={p.prefix}>
            {p.prefix} ({p.species}){p.prefix === settings.suggested ? " · found in your data" : ""}
          </option>
        ))}
        <option value={OTHER}>Other…</option>
      </select>
      {!known && (
        <input
          id="mito_prefix_other"
          aria-label="Start of mitochondrial gene names"
          placeholder="e.g. MT-"
          value={value}
          onChange={(e) => onChange(e.target.value)}
          className="rounded-[4px] border border-field bg-white px-[10px] py-2 font-mono text-[14px]"
        />
      )}
      {counted && (
        <span className={cn("text-[12.5px] leading-snug", counted.n_genes ? "text-muted" : "font-bold text-accent-text")}>
          {counted.n_genes
            ? `${counted.n_genes.toLocaleString("en-US")} ${counted.n_genes === 1 ? "gene matches" : "genes match"}: ${counted.examples.join(", ")}${counted.n_genes > counted.examples.length ? ", …" : ""}`
            : "No genes match. % mitochondrial can't be measured, so its filter will be switched off."}
        </span>
      )}
    </div>
  );
}
