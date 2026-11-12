"""Analysis parameters: one place for every setting used by the backend.

Each analysis step has a small, typed, read-only dataclass. Defaults are the values of
the validated PBMC3k reference analysis (`reference/scrna-seq-validated.ipynb`), which in
turn are the original notebook's values with Scanpy's defaults written out.

These are *inputs* to the analysis (thresholds, settings, seeds, marker sets). Expected
*results* (cell counts, number of clusters, labels) never belong here; they live in
`validation/` and `tests/`.

Use `REFERENCE_PRESET` to reproduce the reference analysis, or build an `AnalysisParams`
with changed values, e.g. `replace(REFERENCE_PRESET, pca=PCAParams(n_comps=50))`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from types import MappingProxyType
from typing import Any, Literal, Mapping

# One seed for every step that uses randomness (PCA, neighbors, Leiden, UMAP).
# 0 is Scanpy's default, which the original notebook used implicitly.
SEED = 0

# Canonical human PBMC marker sets used to score clusters (cell type -> marker genes).
# CD3D marks all T cells, so it is a context marker only (see AnnotationParams).
PBMC_MARKER_SETS: Mapping[str, tuple[str, ...]] = MappingProxyType({
    "CD4+ T cells": ("IL7R", "CCR7"),
    "CD8+ T cells": ("CD8A", "CD8B"),
    "NK cells": ("GNLY", "NKG7"),
    "B cells": ("MS4A1", "CD79A"),
    "CD14+ Monocytes": ("CD14", "LYZ"),
    "FCGR3A+ Monocytes": ("FCGR3A", "MS4A7"),
    "Dendritic cells": ("FCER1A", "CST3"),
    "Megakaryocytes": ("PPBP",),
})


@dataclass(frozen=True)
class QCParams:
    """QC metric calculation (`sc.pp.calculate_qc_metrics`)."""

    mito_prefix: str = "MT-"  # human gene symbols; mouse would be "mt-"
    qc_vars: tuple[str, ...] = ("mt",)
    percent_top: tuple[int, ...] | None = None
    log1p: bool = False


@dataclass(frozen=True)
class FilterParams:
    """Cell and gene filtering. Upper limits are exclusive (`<`), as in the reference.

    Set an upper limit to `None` to switch that filter off.
    """

    min_genes: int = 200  # sc.pp.filter_cells
    min_cells: int = 3  # sc.pp.filter_genes
    max_n_genes_by_counts: int | None = 2500  # keep cells with n_genes_by_counts < this
    max_pct_counts_mt: float | None = 6.0  # keep cells with pct_counts_mt < this


@dataclass(frozen=True)
class NormalizeParams:
    """Library-size normalization (`sc.pp.normalize_total`), followed by log1p."""

    target_sum: float = 1e4


@dataclass(frozen=True)
class HVGParams:
    """Highly variable genes (`sc.pp.highly_variable_genes`)."""

    flavor: Literal["seurat"] = "seurat"  # expects log-normalized data
    min_mean: float = 0.0125
    max_mean: float = 3
    min_disp: float = 0.5


@dataclass(frozen=True)
class RegressOutParams:
    """Optional regression of technical covariates (`sc.pp.regress_out`)."""

    enabled: bool = True
    keys: tuple[str, ...] = ("total_counts", "pct_counts_mt")


@dataclass(frozen=True)
class ScaleParams:
    """Scaling to unit variance (`sc.pp.scale`)."""

    max_value: float = 10


@dataclass(frozen=True)
class PCAParams:
    """Principal component analysis (`sc.pp.pca`)."""

    n_comps: int = 40
    # Scanpy's default for dense data, written out so a future default change can't
    # silently alter results. Verified bit-identical to the reference (2026-10-12).
    svd_solver: Literal["arpack"] = "arpack"
    random_state: int = SEED


@dataclass(frozen=True)
class NeighborsParams:
    """Nearest-neighbor graph (`sc.pp.neighbors`)."""

    n_neighbors: int = 10
    n_pcs: int = 40
    use_rep: str = "X_pca"
    random_state: int = SEED


@dataclass(frozen=True)
class LeidenParams:
    """Leiden clustering (`sc.tl.leiden`). The flavor is pinned on purpose."""

    flavor: Literal["leidenalg", "igraph"] = "leidenalg"
    resolution: float = 1.0
    n_iterations: int = -1  # -1 = run until the partition stops changing
    random_state: int = SEED
    key_added: str = "leiden"
    # Scanpy's default for leidenalg, written out. The igraph flavor requires False.
    # Verified bit-identical to the reference (2026-10-13).
    directed: bool = True


@dataclass(frozen=True)
class UMAPParams:
    """UMAP embedding (`sc.tl.umap`), for visualization only."""

    random_state: int = SEED


@dataclass(frozen=True)
class MarkerParams:
    """Marker genes per cluster (`sc.tl.rank_genes_groups`)."""

    groupby: str = "leiden"
    method: Literal["wilcoxon", "t-test", "t-test_overestim_var", "logreg"] = "wilcoxon"
    use_raw: bool = True  # test all genes, log-normalized (adata.raw)
    pts: bool = True  # also record the fraction of cells expressing each gene
    # Scanpy defaults the reference used implicitly, written out. Verified bit-identical
    # to the reference (2026-10-13).
    reference: str = "rest"  # compare each cluster with all other cells
    corr_method: Literal["benjamini-hochberg", "bonferroni"] = "benjamini-hochberg"
    tie_correct: bool = False  # no tie correction in the Wilcoxon test


@dataclass(frozen=True)
class AnnotationParams:
    """Marker-based cell-type prediction (see CLAUDE.md, `annotation.py` design).

    The marker sets are a configurable input. Replace them for other tissues or species.
    """

    marker_sets: Mapping[str, tuple[str, ...]] = field(default_factory=lambda: PBMC_MARKER_SETS)
    context_markers: tuple[str, ...] = ("CD3D",)  # shown as evidence, not scored
    ambiguity_margin: float = 0.5  # best - second-best score below this -> ambiguous
    max_de_rank: int = 50  # a supporting marker must rank within the top N DE genes
    min_score: float | None = None  # None = no minimum-score check (as in the reference)
    # App-only display setting (not in the reference): another cell type's marker is shown
    # as "conflicting" when its z-score in the cluster is at least this. Evidence only;
    # it never changes a label. Agreed 2026-10-13.
    conflict_z: float = 1.0


@dataclass(frozen=True)
class AnalysisParams:
    """All parameters for one complete analysis run."""

    qc: QCParams = field(default_factory=QCParams)
    filter: FilterParams = field(default_factory=FilterParams)
    normalize: NormalizeParams = field(default_factory=NormalizeParams)
    hvg: HVGParams = field(default_factory=HVGParams)
    regress_out: RegressOutParams = field(default_factory=RegressOutParams)
    scale: ScaleParams = field(default_factory=ScaleParams)
    pca: PCAParams = field(default_factory=PCAParams)
    neighbors: NeighborsParams = field(default_factory=NeighborsParams)
    leiden: LeidenParams = field(default_factory=LeidenParams)
    umap: UMAPParams = field(default_factory=UMAPParams)
    markers: MarkerParams = field(default_factory=MarkerParams)
    annotation: AnnotationParams = field(default_factory=AnnotationParams)

    def to_dict(self) -> dict[str, Any]:
        """Plain, JSON-serializable dictionary (tuples become lists), for recording."""
        return _to_plain(self)


def _to_plain(value: Any) -> Any:
    """Convert dataclasses, mappings and tuples recursively into dicts and lists."""
    if is_dataclass(value):
        return {f.name: _to_plain(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, Mapping):
        return {key: _to_plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_plain(item) for item in value]
    return value


# The validated PBMC3k reference analysis. Same values as the dataclass defaults; the name
# makes it explicit that a run uses the reference settings.
REFERENCE_PRESET = AnalysisParams()
