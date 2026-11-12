# Census Reference Build Plan

Nov 9, 2026 · supervisor's plan (copied here 2026-11-10; adopted by the analyst on 2026-11-10)

## Goal

We will use CELLxGENE Census as the main public reference for cell type annotation. We pull it once with a build script, not live from the app. The script pins one Census release, keeps human and mouse cells from the tissues we support, and writes per-tissue reference files and marker tables. The app loads those files at startup.

Why prebuild instead of querying Census per user:

- **Reproducible:** the same upload gets the same labels until we deliberately ship a new build.
- **Fast:** Census queries take minutes and pull large data; a prebuilt per-tissue file loads in seconds.
- **Cheap to host:** the server doesn't download gigabytes for each user.

Scope for the proof of concept: human and mouse, broad cell type labels, Python only, free academic and research use.

## Phase 0: proof of concept while the build is written

Get labels showing in the app within days using ready-made tools, then replace or extend them with the Census build.

1. **Reference-based labels:** run CellTypist with a pretrained model (an immune model, or the tissue model closest to the user's data). Pass the app's own cluster labels as the over-clustering input with majority voting, so each cluster gets one label.
2. **Marker-based cross-check:** score each cluster against a marker table (CellMarker 2.0 and CELLxGENE CellGuide), using an overlap test on the cluster's top genes.
3. **Show both results side by side**, each with a confidence score. Label a cluster "Unknown" when confidence is low or the two methods disagree badly. Never force a label.

Nothing here is throwaway: Phase 1 trains CellTypist models on Census data, so the app code that calls CellTypist stays the same.

## Phase 1: the Census build script

One Python script (using the `cellxgene_census` package) runs these steps for each species and tissue. Run it by hand when we want a new data release, never at app runtime.

1. **Pin the release.** Open Census with an explicit `census_version` (a long-term-support release), never "latest". Record that version in the output.
2. **Pick tissues.** Start with a short list (for example blood, lung, brain, liver, kidney) using the `tissue_general` column. Add tissues later as users need them.
3. **Filter cells.** Keep only:
   - `is_primary_data == True` (Census contains duplicate copies of some cells across datasets; this drops them),
   - healthy tissue (`disease == 'normal'`),
   - standard single-cell assays (10x and similar; leave out unusual assays at first).
4. **Check labels before pulling expression.** Query only the cell metadata (`obs`) first. Count cells per cell type and per dataset, and drop cell types with too few cells (for example under 100) or that come from a single dataset.
5. **Downsample.** Cap each cell type at a fixed number of cells per tissue (for example 1,000–2,000), spread across datasets. This keeps files small and stops common cell types drowning out rare ones.
6. **Roll labels up to broad types.** Map each `cell_type_ontology_term_id` up the Cell Ontology to a fixed list of broad types (T cell, B cell, macrophage, fibroblast, ...). Keep the fine label too, in a separate column, for later.
7. **Pull expression** for the selected cells as an AnnData object (`get_anndata`). Keep both gene IDs from Census: `feature_id` (Ensembl) and `feature_name` (symbol).
8. **Normalize** the same way the app normalizes user data (for example total counts to 10,000, then log1p). Reference and user data must be processed identically.
9. **Build three outputs per tissue:**
   - a **CellTypist model** trained on the downsampled cells (`celltypist.train`), with broad labels;
   - a **pseudobulk profile table**: average expression per broad cell type, for simple correlation matching and for debugging;
   - a **marker table**: top genes per broad cell type from a one-vs-rest test (`scanpy.tl.rank_genes_groups`), keeping only genes that are specific (high in that type, low elsewhere).
10. **Merge public marker tables.** Add CellMarker 2.0 and CellGuide markers after mapping their names to Cell Ontology IDs. Keep a `source` column and a count of how many sources support each marker.
11. **Write a manifest** (see the next section) and run the checks before shipping.

## Output files and how the app uses them

Each build writes one versioned folder, for example `refdata_v1/`. The app reads the folder named in its config, so moving to a new build is a one-line change and the old one stays available.

| File | One per | What the app does with it |
| --- | --- | --- |
| `manifest.json` | build | Shows the Census version, build date, tissues, cell counts and source versions on an "About the reference" page |
| `models/<species>_<tissue>.pkl` | species + tissue | CellTypist model used to label clusters |
| `pseudobulk/<species>_<tissue>.parquet` | species + tissue | Average expression per broad cell type, for correlation scores |
| `markers/<species>.parquet` | species | Marker genes per cell type (Census, CellMarker, CellGuide), for the marker cross-check and for showing *why* a label was given |
| `ontology_rollup.csv` | build | Fine Cell Ontology ID to broad label lookup |
| `citations.csv` | build | Every Census dataset used, with its title, DOI and license |

Checks to run before a build ships:

- Label a known, already-annotated public dataset that was **not** in the reference. Broad-label accuracy should be high (agree on a target with the intern, for example 85%+).
- Confirm every tissue has at least the core cell types expected for it.
- Confirm gene IDs in user uploads (symbols or Ensembl) match the reference for at least most genes; warn the user when the match rate is low.

## Watch-outs

1. **Consistent labels matter more than more data.** One source says "T cell", another "T lymphocyte", another "CD4+ T cell". Map every name to a Cell Ontology ID (like `CL:0000084`) before merging. Census and CellGuide already use these IDs; CellMarker 2.0 only partly does, so unmapped CellMarker rows get dropped or mapped by hand.
2. **Broad labels come from the ontology.** Map fine labels up the hierarchy (CD8 effector memory T cell → T cell). We can train on detailed references and still report broad labels.
3. **Tissue matters.** Use the reference for the user's tissue when they tell us. A mismatched reference gives confident but wrong labels. If the tissue is unknown, use a multi-tissue reference and say so in the results.
4. **Gene IDs differ between sources.** Census uses Ensembl IDs; marker tables and many user files use gene symbols. Build one mapping layer and use it everywhere, and report how many of the user's genes matched.
5. **Licenses.** Most Census data is CC BY 4.0, which requires attribution; `citations.csv` and the About page cover that. Check the terms of CellMarker 2.0, PanglaoDB and MSigDB before bundling them. Free academic use is usually fine; redistributing a merged copy inside our app is the part to confirm.
6. **Record versions.** Results change between source releases. The manifest records the Census version, every other source's version or download date, and the build date, and the app shows them next to results.

## Comparing two samples

The two-sample feature should use the same reference build for both samples, so "T cell" means the same thing in each.

- **Correct for batch effects first.** Integrate the two samples (Harmony or scVI) before clustering, or processing differences show up as fake biology.
- **Two kinds of comparison:**
  - *differential abundance*: does the share of each cell type change between samples;
  - *differential expression within a cell type*: do the T cells themselves behave differently.
- **One sample per condition cannot be tested statistically.** Thousands of cells are still one sample. With one each, show descriptive differences only (for example 30% vs 45% T cells), with no p-values. With replicates, use pseudobulk tests (sum counts per cell type per sample, then a DESeq2-style test) and, for abundance, scCODA or Milo from the `pertpy` package.

Running per-cell statistical tests between two single samples produces thousands of "significant" genes that mean little. Avoid that.

## Open questions

- What file format do users upload (.h5ad, 10x output, CSV)?
- Which annotation method does the app use today, if any?
- Which tissues come first? The build list in step 2 is a placeholder.
- For the two-sample feature, will users have replicates, or one sample per condition?
- Where will the reference files live on the server, and how large can they be?

## Milestones

- [ ] Phase 0: CellTypist pretrained model + marker cross-check showing labels in the app
- [ ] Pick the Census release and the first tissue list
- [ ] Build script runs end to end for one tissue in human
- [ ] Ontology roll-up to broad labels working
- [ ] Validation on a held-out public dataset meets the agreed accuracy
- [ ] Extend to all first-round tissues, human and mouse
- [ ] CellMarker 2.0 and CellGuide markers merged; licenses confirmed
- [ ] Manifest, citations and "About the reference" page in the app
- [ ] Two-sample comparison uses integration and the same reference
