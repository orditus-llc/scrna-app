# PBMC3k cell-type annotation review

**Status: analyst-reviewed working annotations (reviewed 2026-10-08).**
These are the labels used for the PBMC3k prototype. They are an analyst's interpretation
of marker-gene evidence. They are **not** independently validated biological ground truth
(no protein data, reference atlas or orthogonal method has been used to confirm them).

- **Use: dataset-specific validation reference only.** The app predicts cell types from
  marker evidence (`annotation.py`) and never uses these labels to annotate. They are
  used to measure agreement between the app's predictions and this review.
- Source: `validation/expected/reference_results.h5ad` (generated 2026-10-08; see
  "Reference provenance" in CLAUDE.md, summary SHA-256 `eb568f94…5465`).
  Expression = `adata.raw` (log-normalized, all genes).
- **The decisions apply to the cluster IDs of that specific reference run.** Leiden IDs can
  change with the Scanpy version, flavor or seed. If the reference is regenerated and the
  clusters change, these decisions must be reviewed again.
- **% = percent of cells in the cluster with any detected expression.** DE rank = rank in that
  cluster's Wilcoxon marker list (1 = top; out of 13,714 genes).
- 10x v1 PBMC3k has high dropout. Low % for a lowly expressed gene (e.g. CD4, CD8B) is weak
  evidence of absence.
- The saved reference outputs were **not** changed. `obs["cell_type"]` in
  `reference_results.h5ad` still holds the automated proposal (the original-style names).

## Evidence and decisions

| Cluster | Cells | Automated proposal | Supporting markers | Conflicting / unexpected markers | Uncertainty | Decision (working label) |
|---|---|---|---|---|---|---|
| 0 | 1131 | CD4+ T cells | CD3D 87% (rank 12), IL7R 66% (rank 43), CCR7 37% (rank 87), SELL 57%; top DE genes are ribosomal genes + LDHB (typical of naive/resting T cells) | CD8A or CD8B detected in 19% of cells (CD8A in 120 cells), so naive CD8 T cells are probably mixed in. CD4 itself is only 10% (dropout, not evidence against). Mix of naive (CCR7) and memory (S100A4 79%) | **Medium.** The T-cell identity is clear. "CD4+" describes most, but not all, cells | **CD4-enriched T cells** (qualified) |
| 1 | 430 | CD14+ Monocytes | LYZ 100% (rank 2), S100A8 99% (rank 3), S100A9 (rank 1), CD14 69% (rank 30), CST3 99% | MS4A7 24% and FCGR3A 9%: the expected monocyte continuum. CD4 25% (monocytes express low CD4) | **Low** | **CD14+ Monocytes** |
| 2 | 342 | B cells | CD79A 94% (rank 2), MS4A1 87% (rank 6), CD74/HLA class II at the top | None notable | **Low** | **B cells** |
| 3 | 326 | CD8+ T cells | CD3D 91% (rank 10), CD8A 46% (rank 14), CD8B 35% (rank 31), GZMK 58% (rank 8), CCL5 93% (rank 1), NKG7 87% (rank 2) | 43% of cells have no CD8A/CD8B detected (dropout and/or cytotoxic CD4/NKT/γδ T cells). IL7R 50%, GNLY 25%, PRF1 41%: an effector/memory cytotoxic profile overlapping NK | **Medium** | **Cytotoxic/memory T cells (CD8-enriched)** (qualified) |
| 4 | 214 | FCGR3A+ Monocytes | FCGR3A 84% (rank 19), MS4A7 70% (rank 37), LST1 99% (rank 1), AIF1, FCER1G | CD14 23%, S100A8 59%: some intermediate (CD14+CD16+) monocytes, expected biologically | **Low** | **FCGR3A+ Monocytes** |
| 5 | 178 | NK cells | NKG7 100% (rank 1), GNLY 90% (rank 5), GZMB 93% (rank 2), PRF1 93% (rank 3), KLRD1 53%, FCGR3A 81% (CD56-dim NK) | CD3D in 19% (34 cells), so a few NKT or cytotoxic T cells may be included. CD3E 44% and CD8A 13% (NK cells can express CD3E mRNA and CD8A) | **Low–medium** | **NK cells** |
| 6 | 36 | Dendritic cells | FCER1A 83% (rank 8), CLEC10A 69% (rank 15), CD1C 50% (rank 64), HLA class II top, CST3 100% | Monocyte genes LYZ 97% and S100A8 42% (partly normal for DCs; possible monocyte doublets). Plasmacytoid DC markers LILRA4/IL3RA in ~5 cells, so a few pDCs are mixed in. CD79A 17%. Small cluster | **Medium** | **Conventional dendritic cells (cDC2-enriched)** (qualified) |
| 7 | 13 | Megakaryocytes | PPBP 100% (rank 4), PF4 100% (rank 3), GP9 92%, ITGA2B 85%, SDPR/GNG11 top | True megakaryocytes rarely circulate; in PBMC these are most likely platelets (or platelet–leukocyte aggregates). LYZ 46%, CST3 62% suggest some aggregates. Only 13 cells | **Low for lineage, high for the specific cell type** | **Platelet-lineage cells** (qualified) |

## Reasoning for the four qualified labels

- **Cluster 0, "CD4-enriched T cells".** The T-cell identity is well supported (CD3D 87%,
  IL7R 66%, CCR7 37%). But 19% of cells have CD8A or CD8B detected, which suggests some naive
  CD8 T cells cluster with the CD4 cells at this resolution. CD4 transcript detection (10%)
  is too low to confirm CD4 identity cell by cell. "Enriched" states the majority population
  without claiming purity.
- **Cluster 3, "Cytotoxic/memory T cells (CD8-enriched)".** The cluster is defined mainly by
  a cytotoxic/effector-memory program (CCL5, NKG7, GZMK, GZMA, PRF1) on a T-cell background
  (CD3D 91%). CD8A/CD8B are the top lineage markers, but 43% of cells have neither detected.
  That may be dropout, or cytotoxic CD4, NKT or γδ T cells. The label describes the shared
  phenotype first and CD8 as the dominant, not exclusive, lineage.
- **Cluster 6, "Conventional dendritic cells (cDC2-enriched)".** FCER1A, CLEC10A and CD1C
  support CD1C+ cDC2 identity. About 5 cells express plasmacytoid DC markers (LILRA4/IL3RA),
  and monocyte genes are partly shared. With only 36 cells, the cluster can't be split further
  reliably. The label names the dominant subtype without excluding the minor populations.
- **Cluster 7, "Platelet-lineage cells".** PPBP, PF4, GP9 and ITGA2B clearly indicate the
  megakaryocyte/platelet lineage. Mature megakaryocytes are rare in peripheral blood, so these
  13 cells are more plausibly platelets or platelet–leukocyte aggregates (supported by LYZ/CST3
  signal). "Platelet-lineage" states what the evidence supports without asserting either.

## Original PBMC3k labels (for comparison)

From `reference/scrna-seq.ipynb` (hard-coded by cluster number; same numbering as the
validated reference run).

| Cluster | Original label | Working label | Changed? |
|---|---|---|---|
| 0 | CD4+ T cells | CD4-enriched T cells | qualified |
| 1 | CD14+ Monocytes | CD14+ Monocytes | no |
| 2 | B cells | B cells | no |
| 3 | CD8+ T cells | Cytotoxic/memory T cells (CD8-enriched) | qualified |
| 4 | FCGR3A+ Monocytes | FCGR3A+ Monocytes | no |
| 5 | NK cells | NK cells | no |
| 6 | Dendritic cells | Conventional dendritic cells (cDC2-enriched) | qualified |
| 7 | Megakaryocytes | Platelet-lineage cells | renamed |

The automated proposal in `validation/expected/` uses the original label names. Both are
kept so results stay comparable with the original analysis and the PBMC3k tutorial.
