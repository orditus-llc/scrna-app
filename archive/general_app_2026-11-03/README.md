# General app: frozen snapshot (archived 2026-11-03)

The generalized app as it was on 2026-11-03, before the dataset metadata / species-and-tissue
confirmation changes. Kept so this version can be run and shown at any time; the main app
(`app/app.py`) keeps changing, this copy does not. (This folder can't use git, so a snapshot
plays the role of a branch.)

**Run it** from the project folder:

```
conda activate scrna
shiny run archive/general_app_2026-11-03/app.py --port 8003
```

then open http://127.0.0.1:8003 (a different port, so it can run next to the main app on 8000).

**What it has** (see CLAUDE.md, Milestone 6, up to "Demo pre-selection"): upload `.h5ad` / 10x `.h5`
or the demo; QC suggested from each dataset; general defaults; clustering, marker genes with a
heatmap; Cell Type Annotation with **Markers from: PanglaoDB or your own marker list (upload)**,
species and tissues chosen on that page (pre-selected for the demo); label review; downloads of
every chart and table; `.h5ad` export and HTML run report.

**What's in the folder**
- `app.py`, `payloads.py`, `report.py`: the server, as on 2026-11-03. Only change: `app.py` finds
  this folder's own code first and the project folder for `data/`.
- `scrna_analysis/`: its own frozen copy of the analysis package. Only change: `marker_db.py` looks
  for PanglaoDB in the project's `data/markers/`.
- `www/`: the built interface it serves. `srcts/` and the build files are the React source, kept
  for reference (no `node_modules`; not rebuilt).
- Data: reads the project's `data/pbmc3k.h5ad` and `data/markers/` (run the scripts in `scripts/`
  if they're missing).

**Not tested.** `pytest` only collects `tests/`. At archiving time the same code passed 448 pytest
and 25 Vitest tests. It needs the `scrna` environment as pinned in `environment.yml` on 2026-11-03.
