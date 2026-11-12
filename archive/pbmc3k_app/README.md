# PBMC3k app: frozen snapshot (archived 2026-10-22)

The app as it was when it reproduced the validated PBMC3k analysis
(`reference/scrna-seq-validated.ipynb`), just before Milestone 6 (generalizing). Kept so the
PBMC3k demo can be shown at any time. The main app (`app/app.py`) keeps changing; this copy
does not.

**Run it** from the project folder:

```
conda activate scrna
shiny run archive/pbmc3k_app/app.py --port 8002
```

then open http://127.0.0.1:8002 (a different port, so it can run next to the main app on 8000).

**What's in it**
- `app.py`, `payloads.py`: the server, as on 2026-10-22.
- `pbmc_analysis/`: its own frozen copy of the analysis package (the main code was renamed
  `scrna_analysis/` and keeps changing).
- `www/`: the built interface it serves. `srcts/` and the build files are the React source,
  kept for reference only (no `node_modules`; not rebuilt).
- Data: reads the project's `data/pbmc3k.h5ad`.

**Not tested.** `pytest` only collects `tests/`. At archiving time the same code passed 332
pytest and 14 Vitest tests and reproduced the validated notebook exactly.
It needs the `scrna` environment as pinned in `environment.yml` on 2026-10-22.
