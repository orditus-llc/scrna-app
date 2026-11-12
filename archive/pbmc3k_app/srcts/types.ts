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

export interface QcPayload {
  metrics_done: boolean;
  filter_status: "done" | "needs_rerun" | "not_run";
  can_run: boolean;
  filter: { current: FilterSettings; reference: FilterSettings };
  applied: (Record<string, unknown> & { confirmation: string }) | null;
  will_clear: string[]; // later steps that re-applying the filter would clear
  n_cells?: number;
  n_genes?: number;
  mito?: { genes: string[]; available: boolean; prefix: string; warnings: string[] };
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
      reference: { n_cells_kept: number; n_genes_kept: number }; // the same preview at the reference settings
      filters: FilterRow[];
    }
  | { ok: false; message: string };

/** One filter, in the order they run: your setting next to the reference (null = no limit). */
export interface FilterRow {
  key: keyof FilterSettings;
  what: "cells" | "genes";
  reason: string;
  setting: number | null;
  reference_setting: number | null;
  removed: number | null; // null when the filter is switched off
  reference_removed: number | null;
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
  dataset: { file: string; n_cells: number | null; n_genes: number | null; fingerprint: string | null };
  reference_settings: boolean;
  steps: Step[];
  next_step: NextStep | null; // null when every step is up to date
  settings_record: SettingItem[];
  versions: Record<string, string>;
  run_log: LogEntry[];
  notice: Notice | null;
  remaining_steps: number; // steps not up to date
  run_remaining: RunRemainingJob | null; // the latest "Run the remaining steps" job
}

/** A "Run the remaining steps" job (app/payloads.py, run_remaining_payload). */
export interface RunRemainingJob {
  id: number;
  status: "running" | "done" | "failed";
  current: string | null; // the step running now
  index: number | null; // its position in the job, from 1
  total: number; // steps in the job
  failed_step: string | null; // the step that failed, if any
}

// --- Buttons that send settings with the click (app.py, `clicked` and `sent_settings`) ---

/** A click counter, or the counter plus the settings on screen when the button was clicked. */
export type Click<S> = number | { n: number; settings: S };

export function clickCount<S>(click: Click<S>): number {
  return typeof click === "number" ? click : click.n;
}

// --- Preprocessing card (app/payloads.py, preprocess_payload) -------------------------------

export interface PreprocessSettings {
  target_sum: number;
  min_mean: number;
  max_mean: number;
  min_disp: number;
  max_value: number;
  regress: boolean;
}

/** What the Preprocessing boxes send (strings while typing). */
export interface PreprocessBoxes {
  target_sum: string;
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
    notes: string[];
  } | null;
}

// --- Cell Type Annotation page (app/payloads.py, annotation_payload) ------------------------

export interface AnnotationSettings {
  ambiguity_margin: number;
  max_de_rank: number;
  min_score: number | null; // null = no minimum-score check
  conflict_z: number;
}

export interface AnnotationBoxes {
  ambiguity_margin: string;
  max_de_rank: string;
  min_score: string;
  min_score_off: boolean;
  conflict_z: string;
}

export interface MarkerGene {
  gene: string;
  present: boolean | null; // null until the genes in the data are known
}

export interface ClusterPrediction {
  cluster: string;
  n_cells: number;
  predicted_label: string;
  status: "confident" | "ambiguous" | "unsupported" | "low_score";
  reasons: string;
  best_candidate: string;
  second_candidate: string | null;
  score: number;
  margin: number | null;
  best_marker_DE_rank: number | null;
  supporting_markers: string;
  conflicting_markers: string;
}

export interface EvidenceRow {
  cluster: string;
  cell_type: string; // the marker set the gene belongs to ("context" = shown, not scored)
  gene: string;
  mean_expr: number;
  frac_expressing: number;
  z_score: number;
  de_rank: number | null;
  logfoldchange: number | null;
}

export interface AnnotationPayload {
  status: StepStatus;
  settings: { current: AnnotationSettings; reference: AnnotationSettings };
  marker_sets: { sets: { cell_type: string; markers: MarkerGene[] }[]; context: MarkerGene[] };
  unresolved_label: string;
  will_clear: string[];
  results: {
    settings: AnnotationSettings;
    clusters: ClusterPrediction[];
    scores: Record<string, Record<string, number | null>>; // cluster -> cell type -> score
    evidence: EvidenceRow[];
    label_counts: Record<string, number>;
    missing_genes: Record<string, string[]>;
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
      marker_set_order: string[];
    } | null;
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
