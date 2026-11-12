"""Analysis parameters: one place for every setting used by the backend.

Each analysis step has a small, typed, read-only dataclass. The defaults are general
starting points for any dataset, each with a named source written next to it (mostly
Scanpy's own defaults, the single-cell best-practices guide (Heumos et al. 2023,
Nat Rev Genet), or Seurat's defaults). Where a value is the app's own rule, it says so.
QC filter thresholds are suggested from each dataset's own cells (`SuggestionParams`).

These are *inputs* to the analysis (thresholds, settings, seeds, marker sets). Expected
*results* (cell counts, number of clusters, labels) never belong here.

`DEFAULT_PARAMS` is what the app uses. Build other settings with `dataclasses.replace`,
e.g. `replace(DEFAULT_PARAMS, pca=PCAParams(n_comps=30))`. (The settings of the validated
PBMC3k analysis, used to check the calculations, live in `tests/reference_preset.py`.)
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from types import MappingProxyType
from typing import Any, Literal, Mapping

# One seed for every step that uses randomness (PCA, neighbors, Leiden, UMAP).
# Source: Scanpy's default (random_state=0).
SEED = 0

# How mitochondrial genes are named in different species: gene-name prefix -> species that use
# it. Used to suggest the prefix for a dataset (qc.mito_prefix_matches); any other prefix can be
# entered. Mouse and zebrafish share "mt-", so the data alone can't tell them apart.
MITO_PATTERNS: Mapping[str, str] = MappingProxyType({
    "MT-": "human and many other mammals",
    "mt-": "mouse, zebrafish",
    "Mt-": "rat",
    "mt:": "fruit fly",
})


@dataclass(frozen=True)
class QCParams:
    """QC metric calculation (`sc.pp.calculate_qc_metrics`)."""

    # Human gene symbols. The app sets it from the data (qc.suggest_mito_prefix); see
    # MITO_PATTERNS for other species.
    mito_prefix: str = "MT-"
    qc_vars: tuple[str, ...] = ("mt",)
    percent_top: tuple[int, ...] | None = None
    log1p: bool = False


@dataclass(frozen=True)
class SampleLimits:
    """One sample's cell limits (several samples, Milestone 6b); same meaning as in FilterParams."""

    sample: str
    min_genes: int = 0
    max_n_genes_by_counts: int | None = None
    max_pct_counts_mt: float | None = None


@dataclass(frozen=True)
class FilterParams:
    """Cell and gene filtering. Upper limits are exclusive (`<`).

    Set an upper limit to `None` to switch that filter off. The defaults filter nothing but
    genes seen in fewer than 3 cells: in the app, the thresholds are suggested from each
    dataset's own QC metrics after QC runs (`SuggestionParams`, `qc.suggest_filter`).

    Several samples: `sample_limits` gives each sample its own cell limits (each cell is judged
    by its sample's); the three cell limits above are then not used for cells (they hold the
    limits for all samples together, offered when switching to one set of limits). `min_cells`
    always applies to all cells together. Empty = the same limits for every cell (as before).
    """

    min_genes: int = 0  # sc.pp.filter_cells (0 = keep every cell)
    min_cells: int = 3  # sc.pp.filter_genes; common convention (Seurat tutorial, min.cells = 3)
    max_n_genes_by_counts: int | None = None  # keep cells with n_genes_by_counts < this
    max_pct_counts_mt: float | None = None  # keep cells with pct_counts_mt < this
    # Left out of the record when empty, so a single dataset's record is as before.
    sample_limits: tuple[SampleLimits, ...] = field(default=(), metadata={"omit_if_empty": True})


def filter_from_dict(values: Mapping[str, Any]) -> FilterParams:
    """FilterParams from its plain form (e.g. `dataclasses.asdict`), per-sample limits included."""
    values = dict(values)
    values["sample_limits"] = tuple(SampleLimits(**s) if isinstance(s, Mapping) else s
                                    for s in values.get("sample_limits") or ())
    return FilterParams(**values)


@dataclass(frozen=True)
class NormalizeParams:
    """Library-size normalization (`sc.pp.normalize_total`), followed by log1p."""

    # Counts per cell after normalization. None = the median total counts of the cells
    # (computed from the data). Source: Scanpy's default.
    target_sum: float | None = None


@dataclass(frozen=True)
class HVGParams:
    """Highly variable genes (`sc.pp.highly_variable_genes`).

    Source: Scanpy's defaults (flavor "seurat" with these cut-offs, taken from Seurat).
    """

    flavor: Literal["seurat"] = "seurat"  # expects log-normalized data
    min_mean: float = 0.0125
    max_mean: float = 3
    min_disp: float = 0.5
    # Set = keep the N most variable genes instead of applying the cut-offs above (same method,
    # ranked by normalized dispersion). None = use the cut-offs (Scanpy's default).
    n_top_genes: int | None = None


# What the app offers when the user picks "Top N most variable genes".
# Source: Seurat's default (FindVariableFeatures, nfeatures = 2000).
DEFAULT_N_TOP_GENES = 2000


@dataclass(frozen=True)
class RegressOutParams:
    """Optional regression of technical covariates (`sc.pp.regress_out`).

    Off by default. Source: the best-practices guide (Heumos et al. 2023) advises against it
    unless there is a specific reason, as it can remove biological signal.
    """

    enabled: bool = False
    keys: tuple[str, ...] = ("total_counts", "pct_counts_mt")


@dataclass(frozen=True)
class ScaleParams:
    """Scaling to unit variance (`sc.pp.scale`)."""

    # Scaled values above this are clipped. Source: Seurat's default (ScaleData, scale.max = 10).
    max_value: float = 10


@dataclass(frozen=True)
class PCAParams:
    """Principal component analysis (`sc.pp.pca`)."""

    n_comps: int = 50  # Source: Scanpy's default (settings.N_PCS = 50)
    # Written out so a future Scanpy default change can't silently alter results.
    svd_solver: Literal["arpack"] = "arpack"
    random_state: int = SEED


@dataclass(frozen=True)
class NeighborsParams:
    """Nearest-neighbor graph (`sc.pp.neighbors`)."""

    n_neighbors: int = 15  # Source: Scanpy's default
    n_pcs: int = 50  # all computed components; Source: Scanpy's default (n_pcs=None uses them all)
    use_rep: str = "X_pca"
    random_state: int = SEED


@dataclass(frozen=True)
class LeidenParams:
    """Leiden clustering (`sc.tl.leiden`). The flavor is pinned on purpose.

    Source: Scanpy's defaults (flavor "leidenalg", resolution 1, iterations until stable).
    """

    flavor: Literal["leidenalg", "igraph"] = "leidenalg"
    resolution: float = 1.0
    n_iterations: int = -1  # -1 = run until the partition stops changing
    random_state: int = SEED
    key_added: str = "leiden"
    # Scanpy's default for leidenalg, written out. The igraph flavor requires False.
    # Verified bit-identical to the reference (2026-10-13).
    directed: bool = True


# Clusters with fewer cells than this get a note that their marker genes are less reliable.
# The app's own rule (a display threshold; it never changes the clusters).
SMALL_CLUSTER_CELLS = 10

# Several samples in one analysis (Milestone 6b). The app's own limit for now: every sample's cells
# are held in memory together, and copies are kept per step.
MAX_SAMPLES = 20
# Column names that usually hold the sample in a combined file, best guess first (compared ignoring
# case). Only used to suggest a column; the user chooses.
SAMPLE_COLUMN_HINTS = ("sample", "sample_id", "orig.ident", "batch", "donor", "replicate", "patient",
                       "library", "channel", "individual")
CONDITION_COLUMN_HINTS = ("condition", "treatment", "group", "label", "stim", "status", "disease")


@dataclass(frozen=True)
class IntegrationParams:
    """Batch correction for several samples (Harmony; analyst + supervisor, 2026-11-10).

    Harmony (Korsunsky et al. 2019, harmonypy) adjusts the PCA coordinates so that cells of the same
    type line up across batches; the neighbor graph, clustering and UMAP then use the corrected
    coordinates. Gene expression is not changed. Off by default (analyst's choice): the mixing score
    shows whether it is needed.
    """

    method: Literal["none", "harmony"] = "none"
    batch_key: str = "sample"  # the obs column whose values are the batches (analyst: sample by default)
    theta: float = 2.0  # Source: Harmony's default diversity penalty
    max_iter: int = 10  # Source: Harmony's default number of rounds
    random_state: int = SEED


# Sample mixing score (how varied each cell's nearest neighbors are across samples). The app's own
# rules: neighbors looked at (as LISI's default perplexity of 30), and the score below which a cluster
# is flagged as made mostly of some samples.
MIXING_NEIGHBORS = 30
LOW_MIXING = 0.5


@dataclass(frozen=True)
class UMAPParams:
    """UMAP embedding (`sc.tl.umap`), for visualization only."""

    random_state: int = SEED


@dataclass(frozen=True)
class MarkerParams:
    """Marker genes per cluster (`sc.tl.rank_genes_groups`)."""

    groupby: str = "leiden"
    # Source: Wilcoxon is the test Scanpy's documentation recommends for marker genes (its
    # built-in default, the t-test, is meant as a quick first look).
    method: Literal["wilcoxon", "t-test", "t-test_overestim_var", "logreg"] = "wilcoxon"
    use_raw: bool = True  # test all genes, log-normalized (adata.raw)
    pts: bool = True  # also record the fraction of cells expressing each gene
    # Scanpy defaults the reference used implicitly, written out. Verified bit-identical
    # to the reference (2026-10-13).
    reference: str = "rest"  # compare each cluster with all other cells
    corr_method: Literal["benjamini-hochberg", "bonferroni"] = "benjamini-hochberg"
    tie_correct: bool = False  # no tie correction in the Wilcoxon test


@dataclass(frozen=True)
class MarkerDatabaseParams:
    """Which cell types and marker genes to take from a marker database (PanglaoDB).

    The selection follows decoupler's tutorial "Cell type annotation from marker genes"
    (scverse; decoupler 1.9, PanglaoDB): canonical markers only, expressed in more than half
    of that cell type's cells in the chosen species (sensitivity > 0.5).
    """

    species: str | None = None  # "human" or "mouse"; None = not chosen yet
    tissues: tuple[str, ...] = ()  # the database's organ names; () = every tissue
    canonical_only: bool = True  # decoupler tutorial
    min_sensitivity: float = 0.5  # decoupler tutorial; markers must be above this
    min_markers: int = 3  # the app's own rule: fewer markers in the data -> cell type not tested


@dataclass(frozen=True)
class CensusMarkerParams:
    """How a Census marker table is built from reference cells (supervisor-approved proposal,
    2026-11-03/09; values are the app's own rules unless a source is named)."""

    release: str = "2025-11-08"  # the Census release (long-term support), pinned for reproducibility
    min_cells: int = 100  # a cell type needs this many healthy, primary cells in the tissue…
    min_studies: int = 2  # …from at least this many studies, so one study's labeling can't decide it
    cells_per_type: int = 200  # cells sampled per type (at most); see `fetch_reference`
    top_n_genes: int = 50  # marker genes kept per cell type
    min_fraction: float = 0.25  # a marker must be detected in at least this share of the type's cells
    seed: int = 0  # for sampling cells (same cells every time)
    # Genes never used as markers: ribosomal proteins (RPL13A, RPS4X, RPLP0, RPSA, MRPS18B, …) and
    # mitochondrial genes (MT-…). They track sequencing depth and sample handling, which differ between
    # studies, not cell identity (common practice; analyst's decision 2026-11-09 after 26 of 91 blood
    # types' markers were mostly these). Full-name match, case ignored: RPS6KA1 (a kinase) is kept.
    excluded_genes: str = r"^(M?RP[LS]\d+[A-Z]?\d*|RPLP\d|RPSA|MT-.*)$"
    # Labels too vague to be a cell type (dropped before computing markers).
    vague_labels: tuple[str, ...] = ("unknown", "cell", "native cell", "animal cell", "eukaryotic cell",
                                     "somatic cell", "blood cell", "lymphocyte", "leukocyte",
                                     "mononuclear cell", "hematopoietic cell", "myeloid cell")


@dataclass(frozen=True)
class EnrichmentParams:
    """Cell types from many marker sets (a marker database): an over-representation test.

    For each cluster: are more of a cell type's markers among the cluster's top marker genes
    than chance would give? One-sided Fisher's exact test (the hypergeometric test), as used by
    decoupler's `run_ora` in its PanglaoDB tutorial. P-values are adjusted across the cell
    types tested for each cluster (Benjamini-Hochberg).
    """

    top_n_genes: int = 100  # Scanpy's default in sc.tl.marker_gene_overlap (top 100 marker genes)
    max_adj_p: float = 0.05  # the usual significance level
    # The two checks below are the app's own rules (agreed 2026-10-28).
    min_overlap: int = 3  # at least this many of the cell type's markers among the top genes
    ambiguity_fold: float = 10.0  # best adjusted p-value at least this many times smaller than a significant runner-up's


@dataclass(frozen=True)
class CensusReferenceParams:
    """CELLxGENE Census as the first cell-type reference (supervisor's decision, 2026-11-03):
    its marker table for the species and tissue is tested first; PanglaoDB then runs
    automatically for the clusters Census leaves Unresolved, and never replaces a confident
    Census label."""

    use: bool = True  # False = PanglaoDB only (chosen by the user, e.g. when no Census table exists)
    tissue: str | None = None  # the Census tissue (`tissue_general`, e.g. "blood"); None = not chosen
    release: str = CensusMarkerParams.release  # which release's tables to read


@dataclass(frozen=True)
class CellTypistParams:
    """CellTypist's second opinion (supervisor's plan, Phase 0; analyst choices 2026-11-12): a
    ready-made model, majority vote per cluster. The model is suggested from the species and tissue
    (`celltypist_models.SUGGESTED_MODELS`) and can be changed; None = not used."""

    model: str | None = None  # file name, e.g. "Immune_All_Low.pkl"


@dataclass(frozen=True)
class AnnotationParams:
    """Cell Type Annotation: the references (Census, then PanglaoDB) and the test's checks."""

    markers: MarkerDatabaseParams = field(default_factory=MarkerDatabaseParams)  # PanglaoDB
    test: EnrichmentParams = field(default_factory=EnrichmentParams)
    census: CensusReferenceParams = field(default_factory=CensusReferenceParams)
    celltypist: CellTypistParams = field(default_factory=CellTypistParams)


@dataclass(frozen=True)
class SuggestionParams:
    """Filter thresholds suggested from each dataset's own QC metrics (median ± MADs).

    Following Heumos et al. 2023, "Best practices for single-cell analysis across modalities"
    (Nat Rev Genet): a cell is an outlier when a QC metric is more than N median absolute
    deviations (MADs) from the dataset's median. Genes per cell are compared on a log scale
    (log1p), both sides; % mitochondrial on its own scale, upper side only. Analyst choice
    (2026-10-22): 5 MADs for genes, 3 MADs for % mitochondrial.
    """

    n_mads_genes: float = 5.0  # genes per cell, log1p scale, below and above the median
    n_mads_mito: float = 3.0  # % mitochondrial, above the median
    min_cells: int = 3  # genes detected in fewer cells are removed (a standard starting point)
    # Several samples: suggest limits from each sample's own cells (Heumos et al. 2023 recommend QC
    # per sample, since samples differ in depth and quality); False = one set for all samples.
    per_sample: bool = True


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
    integration: IntegrationParams = field(default_factory=IntegrationParams)
    leiden: LeidenParams = field(default_factory=LeidenParams)
    umap: UMAPParams = field(default_factory=UMAPParams)
    markers: MarkerParams = field(default_factory=MarkerParams)
    annotation: AnnotationParams = field(default_factory=AnnotationParams)
    # None = use the `filter` thresholds as given (the reference analysis). Set = after QC,
    # the filter thresholds are suggested from the data (what the app does for any dataset).
    suggest: SuggestionParams | None = None

    def to_dict(self) -> dict[str, Any]:
        """Plain, JSON-serializable dictionary (tuples become lists), for recording."""
        return _to_plain(self)


def _to_plain(value: Any) -> Any:
    """Convert dataclasses, mappings and tuples recursively into dicts and lists."""
    if is_dataclass(value):
        return {f.name: _to_plain(getattr(value, f.name)) for f in fields(value)
                if not (f.metadata.get("omit_if_empty") and not getattr(value, f.name))}
    if isinstance(value, Mapping):
        return {key: _to_plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_plain(item) for item in value]
    return value


# The app's starting settings for any dataset: the general defaults above, with the QC filter
# thresholds suggested from each dataset's own cells (SuggestionParams).
DEFAULT_PARAMS = AnalysisParams(suggest=SuggestionParams())
