"""Draft a broad-type table for each downloaded CellTypist model (label -> broad type).

The two cell-type checks are compared at the level of broad types (`scrna_analysis/broad_types.py`).
Census labels roll up through the Cell Ontology automatically; CellTypist labels have no ontology IDs,
so each model needs a table `scrna_analysis/broad_type_tables/<model>.csv`. This script drafts one:

1. **ontology**: the label (tidied: lower case, plural "cells" -> "cell") equals a Cell Ontology term's
   name or exact synonym; its broad type is then found through the ontology, as for Census labels.
2. **pattern**: otherwise a word pattern in the label (e.g. "CD8", "Mono", "Astro") gives the broad type.
3. **review**: otherwise "own label" (compared by its own name), marked for a person to check.

Every row records how it was found (column `how`), so drafts can be checked; rows a person has
checked or changed get `how = hand` (or `ask`, waiting for a decision) and are kept as they are when
the script runs again. Tables that
already exist are left alone unless their model is named on the command line with --redo.

Run: python scripts/draft_broad_type_tables.py            (models that have no table yet)
     python scripts/draft_broad_type_tables.py --report   (how many rows of each kind, per table)
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scrna_analysis import broad_types, celltypist_models  # noqa: E402
from scrna_analysis.census_markers import ONTOLOGY_FILE, read_cell_ontology  # noqa: E402

# Models for species the app can't type cells for yet (analyst, 2026-11-16): no table drafted.
SKIPPED = {"Adult_CynomolgusMacaque_Hippocampus.pkl", "Adult_Pig_Hippocampus.pkl",
           "Adult_RhesusMacaque_Hippocampus.pkl"}

NONE, OWN = broad_types.NO_BROAD_TYPE, broad_types.OWN_LABEL
# Rows a person has decided ("hand") or that wait for a decision ("ask") are never redrafted.
KEEP = ("hand", "ask")

# Word patterns, checked in order (the first match wins); case is ignored. Specific before general:
# NKT before NK, OPC before oligodendrocyte, plasma before B, macrophage before monocyte
# ("monocyte-derived macrophage").
PATTERNS: list[tuple[str, str]] = [
    (r"^(undefined|unknown|unassigned|unclassified|ambiguous|mixed|other|others|nan|na|none|lq|low[ _-]?quality.*)$", NONE),
    (r"doublet|low[ _-]?quality|damaged|dying|apoptotic|unknown|undefined", NONE),
    (r"plasma(?!cytoid)|plasmablast|antibody[ -]secreting|\bascs?\b|\bpcs?\b|\bpbs?\b", "plasma cell"),
    (r"\bnkt\b|\bnk[ _-]?t\b|\bmait\b|gamma[ -]?delta|\bgd[ _]?t|γδ|\bt[ _-]?cells?\b|\bt[ _](reg|helper|cell|naive|memory|cytotoxic|effector)"
     r"|\btregs?\b|\bcd4(?![0-9])|\bcd8(?![0-9])|\btfh\b|\bth(1|2|17|22)\b|\btem\b|\btcm\b|\btrm\b|\btemra\b|\bctls?\b|thymocyte|\bdn[0-9(]|\bdp\b"
     r"|\btcr|\bcd3\b|\bdouble[ -](negative|positive) t|^t\b", "T cell"),
    (r"^nk\b|\bnk\b|\bnks\b|\bnk[ _-]?cells?\b|natural killer|\bnk[0-9]|\bcd56", "NK cell"),
    (r"\bb[ _-]?cells?\b|\bb[ _](naive|memory|cell|activated)|\bnaive[ _]b\b|\bmemory[ _]b\b|\bpre-?b\b|\bpro-?b\b|\bgc[ _]b\b"
     r"|germinal cent(er|re) b|\bbcells?\b|^b\b", "B cell"),
    (r"microglia|\bmicro\b|\bmg\b", "microglia"),
    (r"macrophage|\bmacro|\bmacs?\b|kupffer|\bmφ\b|\bmf\b|\bmph\b|\bmdms?\b", "macrophage"),
    (r"monocyte|\bmono\b|\bmonos\b|\bcmono|\bncmono|\bcd14\+? ?mono|\bcd16\+? ?mono|\bmo\b", "monocyte"),
    (r"dendritic|\bdcs?[0-9]*\b|\bcdc[0-9s]?\b|\bpdcs?\b|\bmigdc|langerhans|\bas-?dc", "dendritic cell"),
    (r"\bmast", "mast cell"),
    (r"neutrophil|eosinophil|basophil|granulocyte", "granulocyte"),
    (r"megakaryocyte|platelet|\bmk\b|\bmks\b|\bmegak|\bmega\b", "megakaryocyte/platelet"),
    (r"erythro|reticulocyte|\brbcs?\b|red blood|normoblast", "erythroid cell"),
    (r"\bhscs?\b|\bhspcs?\b|\bmpps?\b|\blmpps?\b|\bcmps?\b|\bgmps?\b|\bmeps?\b|\bclps?\b|haem?atopoietic stem|haem?atopoietic progenitor"
     r"|myeloid progenitor|lymphoid progenitor", "blood stem/progenitor cell"),
    (r"\bopcs?\b|oligodendrocyte (precursor|progenitor)|\bcops?\b", "oligodendrocyte precursor cell"),
    (r"oligodendrocyte|\boligos?\b|\bodcs?\b|\bmols?\b|\bmfols?\b|\bnfols?\b", "oligodendrocyte"),
    (r"astrocyte|\bastro|bergmann", "astrocyte"),
    (r"radial glia|\brgs?\b|neuroblast|glioblast|neural (progenitor|stem)|\bnpcs?\b|\bnscs?\b|\bipcs?\b|intermediate progenitor"
     r"|neuroepithel", "neural progenitor"),
    (r"ependym|hypendym|tanycyte", "ependymal cell"),
    (r"neuron|neuronal|\bexc\b|\binh\b|\bex[0-9]|\bin[0-9]|glutamatergic|gabaergic|\bglut\b|\bgaba\b|interneuron|pyramidal"
     r"|purkinje|cajal|\bcr\b|granule cell|motor neuron|^l[0-9]|\b(lamp5|sst|pvalb|vip|sncg|chandelier)\b", "neuron"),
    (r"pericyte", "pericyte"),
    (r"smooth muscle|muscle[ _]smooth|\bsmcs?[0-9]*\b|\bvsmcs?\b", "smooth muscle cell"),
    (r"endothel|\becs?[0-9]*\b|\bcap\b|capillary|arterial|venous|\bvein\b|\bartery\b|lymphatic|\blecs?\b|\bendo\b|sinusoidal"
     r"|\bhev\b", "endothelial cell"),
    (r"fibroblast|\bfib\b|\bfibro|\bfbs?[0-9]*\b", "fibroblast"),
    (r"epithel|\bepi\b|enterocyte|colonocyte|goblet|\bclub\b|ionocyte|\btuft\b|\bat[12]\b|alveolar type|pneumocyte|keratinocyte"
     r"|hepatocyte|cholangiocyte|acinar|\bductal?\b", "epithelial cell"),
]


def tidy(text: str) -> str:
    """Lower case, underscores and hyphens as spaces, plural words made singular."""
    text = re.sub(r"[_\s]+", " ", text.strip().lower())
    words = [w[:-1] if len(w) > 3 and w.endswith("s") and not w.endswith("ss") else w for w in text.split(" ")]
    return " ".join(words)


def ontology_names(path: Path = ONTOLOGY_FILE) -> dict[str, str]:
    """Tidied term name or exact synonym -> term ID (only names that point to one term)."""
    names: dict[str, set[str]] = {}
    term = None
    for line in path.read_text().splitlines():
        if line == "[Term]":
            term = None
        elif line.startswith("id: CL:"):
            term = line[4:].strip()
        elif term and line.startswith("name: "):
            names.setdefault(tidy(line[6:]), set()).add(term)
        elif term and line.startswith("synonym: ") and " EXACT" in line:
            names.setdefault(tidy(line.split('"')[1]), set()).add(term)
    return {name: next(iter(ids)) for name, ids in names.items() if len(ids) == 1}


def draft(label: str, names: dict[str, str], ontology: dict[str, dict]) -> tuple[str, str]:
    """(broad type, how) for one label."""
    key = tidy(label)
    for candidate in (key, f"{key} cell"):
        if candidate in names:
            return broad_types.broad_type_of_term(ontology, names[candidate]), "ontology"
    spaced = re.sub(r"[_\-]+", " ", label)  # "B_c01-TCL1A" -> "B c01 TCL1A", so word boundaries work
    for pattern, value in PATTERNS:
        if re.search(pattern, spaced, flags=re.IGNORECASE):
            return value, "pattern"
    return OWN, "review"


def read_existing(path: Path) -> dict[str, dict]:
    if not path.is_file():
        return {}
    with open(path, newline="") as f:
        return {row["label"]: row for row in csv.DictReader(f)}


def write_table(model: str, names: dict[str, str], ontology: dict[str, dict]) -> Path:
    from celltypist import models as ct_models

    labels = list(ct_models.Model.load(str(celltypist_models.MODEL_DIR / model)).cell_types)
    path = broad_types.table_path(model)
    kept = {label: row for label, row in read_existing(path).items() if row.get("how") in KEEP}
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["label", "broad_type", "how"])
        for label in labels:
            if label in kept:
                writer.writerow([label, kept[label]["broad_type"], kept[label]["how"]])
            else:
                writer.writerow([label, *draft(label, names, ontology)])
    return path


def report() -> None:
    for path in sorted(broad_types.TABLE_DIR.glob("*.csv")):
        rows = list(read_existing(path).values())
        counts: dict[str, int] = {}
        for row in rows:
            counts[row.get("how") or "hand"] = counts.get(row.get("how") or "hand", 0) + 1
        print(f"{path.stem:40s} {len(rows):4d} labels  " + "  ".join(f"{k} {v}" for k, v in sorted(counts.items())))


def main() -> None:
    args = sys.argv[1:]
    if "--report" in args:
        report()
        return
    redo = "--redo" in args
    named = [a for a in args if not a.startswith("--")]
    names, ontology = ontology_names(), read_cell_ontology()
    models = named or sorted(p.name for p in celltypist_models.MODEL_DIR.glob("*.pkl"))
    for model in models:
        if model in SKIPPED:
            continue
        if broad_types.table_path(model).is_file() and not (redo and model in named):
            continue
        print("Drafted", write_table(model, names, ontology).name)
    report()


if __name__ == "__main__":
    main()
