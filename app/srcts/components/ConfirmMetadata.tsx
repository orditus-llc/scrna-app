// Data Overview's "Confirm your data": the species and tissue of the sample (and an optional
// description). Required before the analysis continues, because they decide which cell types Cell
// Type Annotation compares. The demo's arrive confirmed from its metadata; uploads start empty,
// with the species suggested from how the mitochondrial genes are named. The tissue is one of
// CELLxGENE Census's tissues; those with a Census marker reference come first.
import { useShinyInput } from "@posit-dev/shinyreact";
import { Check } from "lucide-react";
import { useEffect, useRef, useState } from "react";

import { clickCount, type CensusTissue, type Click, type MetadataInfo } from "../types";
import { InfoBox } from "./fields";
import { InfoTip } from "./InfoTip";
import { Button } from "./ui/button";

interface MetadataSend {
  species: string;
  tissue: string | null; // the Census tissue
  about: string;
  stage: string | null; // "adult", "developing" or null (not given)
}

export function ConfirmMetadata({ meta }: { meta: MetadataInfo }) {
  const [click, setClick] = useShinyInput<Click<MetadataSend>>("confirm_metadata", 0);
  const [editing, setEditing] = useState(!meta.confirmed);
  const [species, setSpecies] = useState(meta.species ?? meta.suggested_species ?? "");
  const [tissue, setTissue] = useState(meta.tissue ?? "");
  const [about, setAbout] = useState(meta.about ?? "");
  const [stage, setStage] = useState(meta.stage ?? "");
  const sent = useRef(false);

  // The server has accepted the confirmation when the metadata changes: close the form.
  useEffect(() => {
    if (sent.current) {
      sent.current = false;
      setEditing(false);
    }
  }, [meta.text, meta.about, meta.confirmed_by]);

  const needsTissue = meta.marker_species.includes(species);
  const censusTissues = meta.census_tissues[species] ?? [];
  const ready = species !== "" && (!needsTissue || tissue !== "");

  const startEditing = () => {
    setSpecies(meta.species ?? meta.suggested_species ?? "");
    setTissue(meta.tissue ?? "");
    setAbout(meta.about ?? "");
    setStage(meta.stage ?? "");
    setEditing(true);
  };
  const confirm = () => {
    sent.current = true;
    setClick({
      n: clickCount(click) + 1,
      settings: {
        species,
        tissue: needsTissue ? tissue : null,
        about,
        stage: stage || null,
      },
    });
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
                setTissue("");
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
                Where your cells come from. Cell types are compared with CELLxGENE Census reference cells from this tissue,
                and CellTypist gives a second opinion where it has a model.
              </InfoTip>
            </legend>
            {species === "" ? (
              <span className="text-[13px] text-muted">Choose a species to see its tissues.</span>
            ) : !needsTissue ? (
              <span className="text-[13px] text-muted">
                Not needed: there are no cell-type references for this species, so Cell Type Annotation will be left out.
                Every other step works.
              </span>
            ) : censusTissues.length > 0 ? (
              <CensusTissuePicker options={censusTissues} tissue={tissue} onChoose={setTissue} />
            ) : (
              <span className="text-[13px] font-bold text-accent-text">
                The list of tissues isn't available on this installation, so the tissue can't be chosen yet.
              </span>
            )}
          </fieldset>

          <div className="flex max-w-[560px] flex-col gap-[6px]">
            <span className="flex items-center gap-1">
              <label htmlFor="meta_stage" className="text-[14px] font-bold">
                Stage <span className="font-normal text-muted">(optional)</span>
              </label>
              <InfoTip label="About Stage">
                Whether the sample comes from an adult or from an embryo or fetus. Used to suggest a matching CellTypist
                model, for example the developing mouse brain model for embryonic brain.
              </InfoTip>
            </span>
            <select
              id="meta_stage"
              value={stage}
              onChange={(e) => setStage(e.target.value)}
              className="max-w-[320px] rounded-[4px] border border-field bg-white px-[10px] py-2 text-[14px]"
            >
              <option value="">Not specified</option>
              {meta.stage_options.map((o) => (
                <option key={o.value} value={o.value}>
                  {o.label}
                </option>
              ))}
            </select>
          </div>

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
                {species === "" ? "Choose the species." : "Choose the tissue."}
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

/** One Census tissue: those with a Census marker reference first, then the rest. */
function CensusTissuePicker({
  options,
  tissue,
  onChoose,
}: {
  options: CensusTissue[];
  tissue: string;
  onChoose: (tissue: string) => void;
}) {
  const ready = options.filter((o) => o.has_table);
  const others = options.filter((o) => !o.has_table);
  const chosen = options.find((o) => o.tissue === tissue);
  return (
    <div className="flex flex-col gap-2">
      <select
        id="meta_tissue"
        aria-label="Tissue"
        aria-required="true"
        value={tissue}
        onChange={(e) => onChoose(e.target.value)}
        className="max-w-[420px] rounded-[4px] border border-field bg-white px-[10px] py-2 text-[14px]"
      >
        <option value="">Choose a tissue</option>
        {ready.length > 0 && (
          <optgroup label="Census reference ready">
            {ready.map((o) => (
              <option key={o.tissue} value={o.tissue}>
                {o.label}
              </option>
            ))}
          </optgroup>
        )}
        <optgroup label={ready.length > 0 ? "No Census reference yet" : "Tissues"}>
          {others.map((o) => (
            <option key={o.tissue} value={o.tissue}>
              {o.label}
            </option>
          ))}
        </optgroup>
      </select>
      {chosen && !chosen.has_table && (
        <span className="text-[13px] text-muted">
          No Census reference is built for {chosen.label.toLowerCase()} yet, so Cell Type Annotation won't run the marker
          check.{" "}
          {chosen.celltypist_model
            ? "CellTypist will still give its answer for each cluster."
            : "CellTypist has no model for it either, so you'll review the clusters with their marker genes."}
        </span>
      )}
    </div>
  );
}
