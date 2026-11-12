// Data Overview's "Confirm your data": the species and tissue of the sample (and an optional
// description). Required before the analysis continues, because they decide which cell types Cell
// Type Annotation compares. The demo's arrive confirmed from its metadata; uploads start empty,
// with the species suggested from how the mitochondrial genes are named.
import { useShinyInput } from "@posit-dev/shinyreact";
import { Check } from "lucide-react";
import { useEffect, useRef, useState } from "react";

import { clickCount, type Click, type MetadataInfo } from "../types";
import { InfoBox } from "./fields";
import { InfoTip } from "./InfoTip";
import { Button } from "./ui/button";

interface MetadataSend {
  species: string;
  tissues: string[];
  about: string;
}

export function ConfirmMetadata({ meta }: { meta: MetadataInfo }) {
  const [click, setClick] = useShinyInput<Click<MetadataSend>>("confirm_metadata", 0);
  const [editing, setEditing] = useState(!meta.confirmed);
  const [species, setSpecies] = useState(meta.species ?? meta.suggested_species ?? "");
  const [tissues, setTissues] = useState<string[]>(meta.tissues);
  const [about, setAbout] = useState(meta.about ?? "");
  const sent = useRef(false);

  // The server has accepted the confirmation when the metadata changes: close the form.
  useEffect(() => {
    if (sent.current) {
      sent.current = false;
      setEditing(false);
    }
  }, [meta.text, meta.about, meta.confirmed_by]);

  const counts = meta.tissue_cell_types[species] ?? {};
  const needsTissue = meta.marker_species.includes(species);
  const ready = species !== "" && (!needsTissue || tissues.length > 0);

  const startEditing = () => {
    setSpecies(meta.species ?? meta.suggested_species ?? "");
    setTissues(meta.tissues);
    setAbout(meta.about ?? "");
    setEditing(true);
  };
  const confirm = () => {
    sent.current = true;
    setClick({ n: clickCount(click) + 1, settings: { species, tissues: needsTissue ? tissues : [], about } });
  };

  return (
    <section aria-label="Confirm your data" className="flex flex-col gap-3">
      <h2 className="m-0 border-b border-field pb-[10px] text-[16px] font-bold">Confirm your data</h2>

      {!editing && meta.confirmed ? (
        <InfoBox label="Confirmed species and tissue" className="flex flex-wrap items-start gap-x-4 gap-y-2 text-[14px]">
          <span className="flex min-w-0 flex-[1_1_360px] flex-col gap-[3px]">
            <span className="flex items-center gap-2 font-bold text-ink">
              <Check aria-hidden="true" size={16} strokeWidth={3} className="text-ok" /> {meta.text}
            </span>
            <span className="text-[13px] text-muted">
              {meta.confirmed_by === "dataset" ? "From the dataset's metadata." : "Confirmed by you."}
              {meta.about ? ` ${meta.about}` : ""}
            </span>
          </span>
          {/* Species and tissue supplied with the dataset (e.g. the demo) are fixed. */}
          {meta.confirmed_by !== "dataset" && (
            <Button variant="secondary" size="small" onClick={startEditing}>
              Edit
            </Button>
          )}
        </InfoBox>
      ) : (
        <div className="flex max-w-[760px] flex-col gap-4">
          <p className="m-0 text-[14px] text-body">
            Fields marked <Required /> are required before the analysis can continue. They decide which cell types are
            compared on Cell Type Annotation.
          </p>

          <div className="flex max-w-[420px] flex-col gap-[6px]">
            <span className="text-[14px] leading-snug">
              <label htmlFor="meta_species" className="font-bold">
                Species <Required />
              </label>{" "}
              <InfoTip label="About Species">
                The species your cells come from. Cell types can be predicted for human and mouse; other species run every
                other step.
              </InfoTip>
            </span>
            {meta.suggestion_reason && (
              <span className="text-[13px] text-muted">
                Suggested from the gene names: {meta.suggestion_reason}{meta.suggestion_reason.endsWith("…") ? "" : "."}{" "}
                <InfoTip label="About the suggestion">
                  Read from how the genes are named (Ensembl ID prefix, how symbols are written, mitochondrial genes),
                  not from their expression. Data renamed to another species' genes would not be noticed.
                </InfoTip>
              </span>
            )}
            <select
              id="meta_species"
              value={species}
              aria-required="true"
              onChange={(e) => {
                const next = e.target.value;
                setSpecies(next);
                setTissues(tissues.filter((t) => (meta.tissue_cell_types[next]?.[t] ?? 0) > 0));
              }}
              className="rounded-[4px] border border-field bg-white px-[10px] py-2 text-[14px]"
            >
              <option value="">Choose a species</option>
              {meta.species_options.map((o) => (
                <option key={o.value} value={o.value}>
                  {o.label}
                </option>
              ))}
            </select>
          </div>

          <fieldset className="m-0 flex flex-col gap-2 border-0 p-0">
            <legend className="float-left w-full p-0 text-[14px] leading-snug">
              <span className="font-bold">
                Tissue {needsTissue && <Required />}
              </span>{" "}
              <InfoTip label="About Tissue">
                Where your cells come from. Tick the tissue, plus Immune system if the sample holds immune cells. The number
                is how many cell types the marker database lists for it.
              </InfoTip>
            </legend>
            {species === "" ? (
              <span className="text-[13px] text-muted">Choose a species to see its tissues.</span>
            ) : !needsTissue ? (
              <span className="text-[13px] text-muted">
                Not needed: the marker database has no cell types for this species, so Cell Type Annotation will be left
                out. Every other step works.
              </span>
            ) : !meta.tissues_available ? (
              <span className="text-[13px] font-bold text-accent-text">{meta.tissue_problem}</span>
            ) : (
              <TissueChoices counts={counts} chosen={tissues} onChange={setTissues} />
            )}
          </fieldset>

          <div className="flex max-w-[560px] flex-col gap-[6px]">
            <label htmlFor="meta_about" className="text-[14px] font-bold">
              About this sample <span className="font-normal text-muted">(optional)</span>
            </label>
            <textarea
              id="meta_about"
              value={about}
              maxLength={meta.max_about_length}
              rows={2}
              placeholder="e.g. healthy donor, fresh tissue"
              onChange={(e) => setAbout(e.target.value)}
              className="rounded-[4px] border border-field bg-white px-[10px] py-2 text-[14px]"
            />
            <span className="text-[12px] text-muted">Shown in the run report; not used in the analysis.</span>
          </div>

          <div className="flex flex-wrap items-center gap-3">
            <Button disabled={!ready} onClick={confirm}>
              Confirm species and tissue
            </Button>
            {meta.confirmed && (
              <Button variant="secondary" onClick={() => setEditing(false)}>
                Cancel
              </Button>
            )}
            {!ready && (
              <span className="text-[13px] text-muted">
                {species === "" ? "Choose the species." : "Tick at least one tissue."}
              </span>
            )}
          </div>
        </div>
      )}

      {meta.mismatch && (
        <p className="m-0 text-[14px] font-bold text-accent-text" role="status">
          ▲ {meta.mismatch}
        </p>
      )}
    </section>
  );
}

/** The red asterisk of a required field. */
function Required() {
  return (
    <span className="font-extrabold text-required" aria-label="required">
      *
    </span>
  );
}

/** Tissue checkboxes, each with its number of cell types (tissues with none aren't listed). */
function TissueChoices({
  counts,
  chosen,
  onChange,
}: {
  counts: Record<string, number>;
  chosen: string[];
  onChange: (tissues: string[]) => void;
}) {
  const toggle = (t: string, on: boolean) => onChange(on ? [...chosen, t] : chosen.filter((c) => c !== t));
  return (
    <div className="grid grid-cols-[repeat(auto-fill,minmax(210px,1fr))] gap-x-4 gap-y-[6px]">
      {Object.keys(counts)
        .sort()
        .map((t) => (
          <label key={t} className="flex items-start gap-2 text-[14px] leading-snug">
            <input
              type="checkbox"
              checked={chosen.includes(t)}
              onChange={(e) => toggle(t, e.target.checked)}
              className="mt-[2px] size-4 flex-none accent-ink"
            />
            <span className="min-w-0">
              {t} <span className="font-mono text-[12px] whitespace-nowrap text-muted">{counts[t]}</span>
            </span>
          </label>
        ))}
    </div>
  );
}
