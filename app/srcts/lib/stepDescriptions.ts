// What each step does, in one or two sentences: shown behind the ⓘ next to the step in the
// sidebar (analyst, 2026-11-04: no description at the top of the pages).
export const STEP_DESCRIPTIONS: Record<string, string> = {
  data_overview: "Load a single-cell dataset, check that it holds raw counts, and confirm its species and tissue.",
  quality_control: "Choose which cells to keep. Changes are previewed; nothing is removed until you click Apply Filtering.",
  preprocessing: "Makes cells comparable and keeps the most informative genes.",
  dimensionality_reduction:
    "Summarizes cells with principal components, links each cell to its most similar cells, and draws a UMAP.",
  clustering: "Groups similar cells into clusters with the Leiden algorithm. Cluster numbers are only labels.",
  marker_genes: "Ranks the genes that set each cluster apart from all other cells. These markers are used to name the clusters.",
  cell_type_annotation:
    "Proposes a cell type for each cluster: which cell type's known markers are over-represented among the cluster's top marker genes. These are predictions to review; unclear clusters stay Unresolved.",
  explorer: "Color the UMAP by cluster, cell type, a QC value or any gene. Nothing is run here.",
  compare_conditions:
    "Tests what changed between two conditions, counting each sample once: each group's share of the cells, and each group's genes (pseudobulk).",
};
