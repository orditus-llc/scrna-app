# Classic Milestone 4 app (archived 2026-10-16)

The first Shiny "Core" app: one tab per stage, static matplotlib plots
(`app.py`, `stage_panels.py`, `display.py`, `plots.py`) and its tests (`test_app.py`).
Replaced by the React app in `app/` (Milestone 5) and kept here only as a record.

It is no longer run or tested: `pytest` only collects `tests/`, and the file paths inside
assume the old location `app/classic/`. To run it again, move the four `.py` files back to
`app/classic/` and `test_app.py` back to `tests/`.

Note (2026-10-22): the analysis package was renamed `pbmc_analysis` -> `scrna_analysis`. These
files still import the old name; running them would need that rename too (or the frozen copy in
`archive/pbmc3k_app/pbmc_analysis/`).
