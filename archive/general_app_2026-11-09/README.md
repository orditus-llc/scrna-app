# General app: frozen snapshot (archived 2026-11-09)

The single-sample app as it was on 2026-11-09, before multiple samples (Milestone 6b). Kept so this
version can be run and shown at any time; the main app keeps changing, this copy does not. (This
folder can't use git, so a snapshot plays the role of a branch.)

**Run it** from the project folder:

```
conda activate scrna
shiny run archive/general_app_2026-11-09/app.py --port 8004
```

then open http://127.0.0.1:8004.

**What it has** (CLAUDE.md, Milestone 6 up to 2026-11-05): one dataset per analysis; species and tissue
confirmed on Data Overview (demo's fixed); QC suggested from the data; PanglaoDB annotation with
guidance for Unresolved clusters and species clues; label review; downloads; `.h5ad` export and run
report; grey info boxes, sidebar panel with ⓘ descriptions, four-size typography.

**What's in the folder:** `app.py`, `payloads.py`, `report.py` (only change: `app.py` finds this
folder's own code first and the project folder for `data/`); `scrna_analysis/` (own frozen copy; only
change: `marker_db.py` looks for PanglaoDB in the project's `data/markers/`); `www/` (the built
interface); `srcts/` and build files for reference (no `node_modules`; not rebuilt).

**Not tested.** `pytest` only collects `tests/`. At archiving time the same code passed 454 pytest and
25 Vitest tests. Needs the `scrna` environment as pinned in `environment.yml` on 2026-11-09.
