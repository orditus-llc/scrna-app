"""The settings of the validated PBMC3k analysis (`reference/scrna-seq-validated.ipynb`).

Used only by the tests: given these settings and data/pbmc3k.h5ad, the code must reproduce
the notebook's saved results (`validation/expected/`). That checks the calculations are right.
These are not the app's defaults (see `scrna_analysis/params.py`, `DEFAULT_PARAMS`): the
notebook's own choices differ from the general defaults where shown below, and its QC
thresholds are fixed instead of suggested from the data.
"""

from scrna_analysis.params import (
    AnalysisParams,
    AnnotationParams,
    CensusReferenceParams,
    MarkerDatabaseParams,
    FilterParams,
    NeighborsParams,
    NormalizeParams,
    PCAParams,
    RegressOutParams,
)

REFERENCE_PRESET = AnalysisParams(
    # Fixed QC thresholds (the app suggests them from the data instead).
    filter=FilterParams(min_genes=200, min_cells=3, max_n_genes_by_counts=2500, max_pct_counts_mt=6.0),
    normalize=NormalizeParams(target_sum=1e4),  # app default: the median of the cells
    regress_out=RegressOutParams(enabled=True),  # app default: off
    pca=PCAParams(n_comps=40),  # app default: 50
    neighbors=NeighborsParams(n_neighbors=10, n_pcs=40),  # app default: 15 neighbors, 50 PCs
    # Cell Type Annotation: the notebook's own marker-set scoring was retired (2026-10-28); the
    # app compares clusters with the CELLxGENE Census reference. Species and tissue are set here
    # only so that whole analyses can run in the tests; their labels aren't compared with the
    # notebook's. (No CellTypist model: the reference analysis had none.)
    annotation=AnnotationParams(markers=MarkerDatabaseParams(species="human"),
                                census=CensusReferenceParams(tissue="blood")),
    # Everything else (HVG cut-offs, scaling, Leiden, UMAP, markers, seed) is the same as the
    # app's defaults. suggest=None: thresholds are used as given.
)
