# General app: frozen snapshot (archived 2026-11-10)

The app as it was on 2026-11-10, before the supervisor's two-check annotation plan (CellTypist + a
marker check side by side, broad Cell Ontology labels) and batch correction. Kept so this version can
be run and shown at any time; the main app keeps changing, this copy does not. (This folder can't use
git, so a snapshot plays the role of a branch.)

**Run it** from the project folder:

```
conda activate scrna
shiny run archive/general_app_2026-11-10/app.py --port 8005
```

then open http://127.0.0.1:8005.

**What it has** (CLAUDE.md, Milestone 6 up to 2026-11-09): single datasets and several samples (one file
per sample, or a sample column; Kang et al. 2018 as the second demo); Quality Control limits per sample;
species and tissue confirmed on Data Overview, with a CELLxGENE Census tissue picker; Cell Type
Annotation with Census marker tables first (human blood, mouse brain; ribosomal and mitochondrial genes
left out) and PanglaoDB automatically for clusters Census leaves Unresolved, "How was this
identified?" per cluster; label review; downloads; `.h5ad` export and run report; QC preview with kept
cells teal and removed cells grey; American spelling.

**Known issues at archiving time:** two wrong Census labels (the demo's B cells as "IgG plasma cell", a
183-cell, 2-study Census type; young mouse neurons as "erythrocyte"); sibling cell types (e.g.
immature vs naive B cell) count as ambiguous. No batch correction.

**What's in the folder:** `app.py`, `payloads.py`, `report.py` (only change: `app.py` finds this
folder's own code first and the project folder for `data/`); `scrna_analysis/` (own frozen copy; only
change: `marker_db.py` and `census_markers.py` look for PanglaoDB, the Census tables and the Cell
Ontology in the project's `data/markers/`); `www/` (the built interface); `srcts/` and build files for
reference (no `node_modules`; not rebuilt).

**Not tested.** `pytest` only collects `tests/`. At archiving time the same code passed 503 pytest and
27 Vitest tests. Needs the `scrna` environment as pinned in `environment.yml` on 2026-11-10 and the
Census tables in `data/markers/census/2025-11-08/` (human_blood, mouse_brain).
