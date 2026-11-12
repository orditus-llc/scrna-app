// The exact shape of the data the Shiny server sends (app/payloads.py, workflow_payload).
// If the server changes the payload, update this file: the type check then shows every
// place in the interface that needs to change.

export type StepStatus = "done" | "needs_rerun" | "not_run";

export interface Step {
  key: string;
  number: number;
  title: string;
  status: StepStatus;
  summary: string;
  reason: string | null;
  confirmation: string | null; // e.g. "Data loaded successfully: …" once the step is done
  prerequisite: string | null; // title of the earlier step that has to come first
}

export interface NextStep {
  key: string;
  title: string;
  text: string;
}

export interface Notice {
  id: number;
  kind: "error" | "info";
  text: string;
  step: string | null; // the Workflow step the message belongs to
}

export interface SettingItem {
  label: string;
  value: string;
}

export interface LogEntry {
  label: string;
  text: string;
  ran_at: string;
}

// --- Quality Control page (app/payloads.py, qc_payload and preview_payload) ---------------

export interface FilterSettings {
  min_genes: number;
  min_cells: number;
  max_n_genes_by_counts: number | null; // null = no limit
  max_pct_counts_mt: number | null; // null = no limit
}

export interface QcStatRow {
  statistic: string;
  n_genes_by_counts: number;
  total_counts: number;
  pct_counts_mt?: number; // missing when mitochondrial QC is unavailable
}

/** The "Mitochondrial genes" setting (app/payloads.py, qc_settings). */
export interface QcSettings {
  mito_prefix: string; // what the box shows: the prefix QC ran with, or the data's suggestion
  suggested: string | null; // the known pattern that matches the most genes, if any
  patterns: { prefix: string; species: string }[];
}

/** How one suggested threshold was found (app/payloads.py, suggestion_payload). `median` and
 * `mad` are on the scale used: log1p(genes per cell), or % mitochondrial. */
export interface SuggestionDetail {
  value: number | null; // null = no limit suggested (see note)
  median: number | null;
  mad: number | null;
  n_mads: number;
  scale: "log1p" | "percent";
  note: string | null;
}

/** Thresholds suggested from the data for some MAD settings (also the `qc_suggestion` output). */
export interface QcSuggestion {
  suggested: FilterSettings | null; // null before QC has run
  details: Record<"min_genes" | "max_n_genes_by_counts" | "max_pct_counts_mt", SuggestionDetail> | null;
  suggest_settings: { n_mads_genes: number; n_mads_mito: number; min_cells: number } | null;
  problem: string | null; // the MAD boxes can't be used (e.g. not a number)
}

/** Live count under the setting (app/payloads.py, mito_check_payload). */
export interface MitoCheck {
  prefix: string;
  n_genes: number;
  examples: string[];
}

export interface QcPayload {
  settings: QcSettings;
  metrics_done: boolean;
  filter_status: "done" | "needs_rerun" | "not_run";
  can_run: boolean;
  filter: { current: FilterSettings } & QcSuggestion;
  applied: (Record<string, unknown> & { confirmation: string }) | null;
  will_clear: string[]; // later steps that re-applying the filter would clear
  n_cells?: number;
  n_genes?: number;
  mito?: { genes: string[]; available: boolean; prefix: string; warnings: string[]; filter_switched_off: boolean };
  stats?: QcStatRow[];
  cells?: { barcode: string[]; n_genes: number[]; total_counts: number[]; pct_mt: number[] };
}

export type PreviewPayload =
  | {
      ok: true;
      settings: FilterSettings;
      n_cells_before: number;
      n_cells_kept: number;
      n_genes_before: number;
      n_genes_kept: number;
      removed_steps: { what: "cells" | "genes"; reason: string; removed: number }[];
      removed_barcodes: string[]; // cells these settings remove, by barcode
      failures: Record<string, string[]>; // barcode -> the filters that cell fails (removed cells only)
      suggested: { n_cells_kept: number; n_genes_kept: number }; // the same preview at the suggested thresholds
      filters: FilterRow[];
    }
  | { ok: false; message: string };

/** One filter, in the order they run: your setting next to the suggested one (null = no limit). */
export interface FilterRow {
  key: keyof FilterSettings;
  what: "cells" | "genes";
  reason: string;
  setting: number | null;
  suggested_setting: number | null;
  removed: number | null; // null when the filter is switched off
  suggested_removed: number | null;
  differs: boolean;
}

/** What the filter boxes send to the server (strings while typing; "" = empty box). */
export interface FilterBoxes {
  min_genes: string;
  min_cells: string;
  max_genes: string;
  max_genes_off: boolean;
  max_mt: string;
  max_mt_off: boolean;
}

export interface WorkflowPayload {
  dataset: {
    file: string | null; // the file's original name (null until data is loaded)
    label: string | null; // how the data is named on screen (e.g. "Demo data: …" or the file name)
    n_cells: number | null;
    n_genes: number | null;
    fingerprint: string | null;
    file_format: string | null; // e.g. "10x Genomics (.h5)"
    notes: string[]; // things to know about the file (e.g. genes named by Ensembl IDs)
    mito_matches: { prefix: string; species: string; n_genes: number; examples: string[] }[]; // most first
  };
  steps: Step[];
  next_step: NextStep | null; // null when every step is up to date
  settings_record: SettingItem[];
  versions: Record<string, string>;
  run_log: LogEntry[];
  notice: Notice | null;
  remaining_steps: number; // steps not up to date
  run_remaining: RunRemainingJob | null; // the latest "Run the remaining steps" job
  export: ExportInfo; // whether the analysis can be downloaded yet (app/payloads.py, export_payload)
}

/** A "Run the remaining steps" job (app/payloads.py, run_remaining_payload). */
export interface RunRemainingJob {
  id: number;
  status: "running" | "done" | "stopped" | "failed"; // stopped = a step waits for the user's choice
  current: string | null; // the step running now
  index: number | null; // its position in the job, from 1
  total: number; // steps in the job
  failed_step: string | null; // the step that failed, if any
  stopped_step: string | null; // the step waiting for the user's choice, if any
}

// --- Buttons that send settings with the click (app.py, `clicked` and `sent_settings`) ---

/** A click counter, or the counter plus the settings on screen when the button was clicked. */
export type Click<S> = number | { n: number; settings: S };

export function clickCount<S>(click: Click<S>): number {
  return typeof click === "number" ? click : click.n;
}

// --- Preprocessing card (app/payloads.py, preprocess_payload) -------------------------------

export interface PreprocessSettings {
  target_sum: number | null; // null = the median total counts of the cells
  hvg_method: "cutoffs" | "top"; // cut-offs (Scanpy's default) or the N most variable genes
  n_top_genes: number; // N for "top" (Seurat's default 2,000 when not set)
  min_mean: number;
  max_mean: number;
  min_disp: number;
  max_value: number;
  regress: boolean;
}

/** What the Preprocessing boxes send (strings while typing). */
export interface PreprocessBoxes {
  target_sum: string;
  target_median: boolean; // "Median of your cells" ticked
  hvg_method: "cutoffs" | "top";
  n_top_genes: string;
  min_mean: string;
  max_mean: string;
  min_disp: string;
  max_value: string;
  regress: boolean;
}

export interface PreprocessPayload {
  status: StepStatus;
  settings: { current: PreprocessSettings; reference: PreprocessSettings };
  fixed: { hvg_flavor: string; regress_keys: string[] };
  will_clear: string[]; // later steps that re-running would clear
  results: {
    n_hvg: number;
    n_genes: number;
    target_sum_used: number | null; // counts per cell actually used (the median when the setting is null)
    settings: PreprocessSettings; // the settings these results were computed with
    genes: { name: string[]; mean: (number | null)[]; dispersion_norm: (number | null)[]; hvg: boolean[] };
  } | null;
}

// --- Dimensionality Reduction card (app/payloads.py, reduction_payload) ----------------------

export interface ReductionSettings {
  n_comps: number;
  n_neighbors: number;
  n_pcs: number;
  seed: number;
}

export interface ReductionBoxes {
  n_comps: string;
  n_neighbors: string;
  n_pcs: string;
  seed: string;
}

export interface ReductionPayload {
  status: StepStatus;
  settings: { current: ReductionSettings; reference: ReductionSettings };
  will_clear: string[];
  results: {
    variance_ratio: number[];
    total_variance_ratio: number;
    n_pcs_used: number;
    elbow_pc: number; // where the variance curve bends (guidance only, never applied)
    umap: { x: (number | null)[]; y: (number | null)[] };
    cells: {
      barcode: string[];
      n_genes: number[];
      total_counts: (number | null)[];
      pct_mt: (number | null)[] | null; // null when mitochondrial QC is unavailable
    };
  } | null;
}

// --- Clustering page (app/payloads.py, cluster_payload) -------------------------------------

export interface ClusterSettings {
  resolution: number;
  seed: number;
}

export interface ClusterBoxes {
  resolution: string;
  seed: string;
}

export interface ClusterPayload {
  status: StepStatus;
  settings: { current: ClusterSettings; reference: ClusterSettings };
  fixed: { flavor: string; n_iterations: number; directed: boolean };
  will_clear: string[];
  results: {
    settings: ClusterSettings; // the settings these clusters were computed with
    n_clusters: number;
    clusters: { cluster: string; n_cells: number }[];
    notes: string[];
    warnings: string[]; // one cluster, or clusters with very few cells
    labels: string[]; // cluster of each cell, same order as umap and cells
    umap: { x: (number | null)[]; y: (number | null)[] };
    cells: { barcode: string[]; n_genes: number[]; total_counts: (number | null)[]; pct_mt: (number | null)[] | null };
  } | null;
}

// --- Marker Genes page (app/payloads.py, markers_payload) -----------------------------------

export interface MarkerSettings {
  method: string;
}

export interface MarkerRow {
  rank: number;
  gene: string;
  score: number | null;
  logfoldchange: number | null;
  pval_adj: number | null;
  frac_in_group: number | null; // 0–1
  frac_in_rest: number | null; // 0–1
}

export interface DotValue {
  group: string;
  gene: string;
  mean_expr: number;
  frac_expressing: number; // 0–1
  rank: number;
  logfoldchange: number | null;
}

export interface MarkersPayload {
  status: StepStatus;
  settings: { current: MarkerSettings; reference: MarkerSettings };
  methods: { value: string; label: string }[];
  fixed: { use_raw: boolean; reference: string; corr_method: string; tie_correct: boolean; groupby: string };
  will_clear: string[];
  results: {
    settings: MarkerSettings;
    method_label: string;
    n_genes_tested: number;
    groups: { cluster: string; n_cells: number }[];
    top: Record<string, MarkerRow[]>; // cluster -> its top genes, in rank order
    dotplot: { genes: string[]; values: DotValue[] };
    heatmap: {
      options: number[]; // top genes per cluster the page offers
      default: number;
      rows: Record<string, { cluster: string; gene: string }[]>; // option -> rows, top to bottom (each gene once)
      values: { group: string; gene: string; mean_expr: number; scaled: number }[]; // scaled 0–1 within each gene
    };
    notes: string[];
  } | null;
}

// --- Cell Type Annotation page (app/payloads.py, annotation_payload) ------------------------

export interface AnnotationSettings {
  source: "panglaodb" | "custom"; // where the markers come from: the database, or the user's own list
  custom_file: string | null; // the user's list's file name (source "custom")
  custom_fingerprint: string | null; // start of its SHA-256 (tells a re-uploaded, edited file apart)
  species: string | null; // "human" | "mouse"; null = not chosen yet
  tissues: string[]; // the database's tissue names; none chosen at first
  top_n_genes: number;
  max_adj_p: number;
  min_overlap: number;
  ambiguity_fold: number;
}

export interface AnnotationBoxes {
  source: "panglaodb" | "custom";
  species: string; // "" = not chosen
  tissues: string[];
  top_n_genes: string;
  max_adj_p: string;
  min_overlap: string;
  ambiguity_fold: string;
}

export type AnnotationStatus = "confident" | "too_few_markers" | "not_significant" | "ambiguous";

export interface ClusterPrediction {
  cluster: string;
  n_cells: number;
  predicted_label: string;
  status: AnnotationStatus;
  reasons: string; // every check that failed, e.g. "not_significant, ambiguous"
  best_candidate: string;
  overlap: number; // the best candidate's markers among the cluster's top marker genes
  n_markers: number; // the best candidate's markers among the genes tested
  adj_p: number;
  matching_markers: string; // comma-separated, in marker-gene order
  second_candidate: string | null;
  second_overlap: number | null;
  second_adj_p: number | null;
}

export interface Candidate {
  cell_type: string;
  overlap: number;
  n_markers: number;
  adj_p: number;
  score: number; // -log10(adjusted p-value)
  matching_markers: string;
}

/** The user's own marker list as the page shows it (app/payloads.py, marker_list_info). */
export interface MarkerListInfo {
  file_name: string;
  fingerprint: string; // start of the file's SHA-256
  n_genes: number;
  cell_types: { cell_type: string; n_genes: number; n_found: number | null }[]; // n_found: genes in this dataset
  min_markers: number; // cell types with fewer genes in the data aren't tested
  notes: string[];
}

/** The analyst's decision for one cluster (app/payloads.py, _review_payload). */
export interface ClusterReview {
  cluster: string;
  decision: "confirmed" | "renamed";
  label: string; // the confirmed prediction, or the analyst's own label
  note: string;
  reviewed_at: string;
  prediction_changed: boolean; // confirmed, but annotation was re-run and now predicts something else
}

export interface ReviewInfo {
  clusters: Record<string, ClusterReview>; // reviewed clusters only
  not_reviewed_label: string; // "Not reviewed"
  n_clusters: number;
  n_confirmed: number;
  n_renamed: number;
  n_not_reviewed: number;
}

export interface AnnotationPayload {
  status: StepStatus;
  settings: { current: AnnotationSettings; reference: AnnotationSettings };
  species_options: { value: string; label: string }[];
  suggested_species: string | null; // from how the mitochondrial genes are named
  mito_prefix: string | null;
  database: {
    available: boolean;
    problem: string | null;
    name: string | null;
    version: string | null;
    tissues: string[];
    cell_types: Record<string, Record<string, string[]>>; // species -> tissue -> cell types passing the selection rule
  };
  selection: { canonical_only: boolean; min_sensitivity: number; min_markers: number };
  unresolved_label: string;
  will_clear: string[];
  review_cleared: string | null; // why an earlier label review was cleared (clusters recomputed)
  marker_list: MarkerListInfo | null; // the user's uploaded marker list (this session), if any
  preset: string | null; // why species/tissues are pre-selected (the dataset's description), if they are
  marker_list_template: string; // the example CSV offered as "Download a template"
  results: {
    settings: AnnotationSettings;
    database: { name: string; version: string };
    n_cell_types_tested: number;
    n_cell_types_not_tested: number;
    n_markers_not_found: number;
    clusters: ClusterPrediction[];
    candidates: Record<string, Candidate[]>; // cluster -> best candidates first
    dotplot: {
      genes: { gene: string; cell_type: string | null }[]; // the cell type it is a marker of
      values: { group: string; gene: string; mean_expr: number; frac_expressing: number }[];
    };
    label_counts: Record<string, number>;
    nothing_matched: string | null; // every cluster Unresolved: what to check, with a tissue hint
    review: ReviewInfo;
    notes: string[];
    cell_cluster: string[]; // per cell, same order as umap and cells
    cell_label: string[];
    umap: { x: (number | null)[]; y: (number | null)[] };
    cells: { barcode: string[]; n_genes: number[]; total_counts: (number | null)[]; pct_mt: (number | null)[] | null };
  } | null;
}

// --- Explore Results page (app/payloads.py, explore_payload and gene_payload) ----------------

export interface ExplorePayload {
  available: boolean;
  reason: string | null;
  results: {
    clusters: { cluster: string; n_cells: number }[];
    cell_cluster: string[]; // per cell, same order as umap and cells
    umap: { x: (number | null)[]; y: (number | null)[] };
    cells: { barcode: string[]; n_genes: number[]; total_counts: (number | null)[]; pct_mt: (number | null)[] | null };
    genes: string[]; // every gene, alphabetical (for the search)
    top_markers: Record<string, string[]> | null; // cluster -> top genes (null if markers aren't up to date)
    predictions: {
      unresolved_label: string;
      clusters: { cluster: string; predicted_label: string; status: string }[];
      cell_label: string[];
      label_order: string[]; // cell types in the order of their first cluster
    } | null;
    reviewed: Record<string, string> | null; // cluster -> the analyst's label (reviewed clusters only)
    not_reviewed_label: string; // "Not reviewed"
    unavailable: string[]; // e.g. "Cell Type Annotation needs re-running"
  } | null;
}

export type GenePayload =
  | {
      ok: true;
      gene: string;
      values: (number | null)[]; // per cell, same order as the explore payload
      max: number;
      n_expressing: number;
      by_cluster: { cluster: string; mean_expr: number; frac_expressing: number }[];
    }
  | { ok: false; gene: string; message: string; suggestions: string[] };

/** The downloads of the whole analysis (app/payloads.py, export_payload). */
export interface ExportInfo {
  available: boolean;
  reason: string | null; // why not, when not available
  includes: string[]; // the steps whose results the file holds
  reviewed: boolean; // the analyst's reviewed labels are included
  h5ad_name: string;
  report_name: string;
}
