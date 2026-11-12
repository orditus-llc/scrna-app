// The app's pages: one page per analysis step (in the order they run), then Explore Results.
// The app opens on the first step. The step keys match the server's (app/payloads.py, STEPS).

export const STEP_KEYS = [
  "data_overview",
  "quality_control",
  "preprocessing",
  "dimensionality_reduction",
  "clustering",
  "marker_genes",
  "cell_type_annotation",
] as const;

export type StepKey = (typeof STEP_KEYS)[number];
export type Page = StepKey | "explorer";

/** Step pages whose controls are built. The others open with a "Coming soon" note. */
export const BUILT_STEPS: Record<StepKey, boolean> = {
  data_overview: true,
  quality_control: true,
  preprocessing: true,
  dimensionality_reduction: true,
  clustering: true,
  marker_genes: true,
  cell_type_annotation: true,
};

/** The step after `key`, or undefined for the last one. */
export function nextStepKey(key: StepKey): StepKey | undefined {
  return STEP_KEYS[STEP_KEYS.indexOf(key) + 1];
}

/** Every page, for reading one from the address. */
export const PAGES: readonly Page[] = [...STEP_KEYS, "explorer"];

/** The page named in the address ("#clustering" -> "clustering"), or null if there is none.
 * Each page has its own address, so the browser's Back and Forward buttons move between
 * pages instead of leaving the app. */
export function pageFromHash(hash: string): Page | null {
  const name = hash.replace(/^#/, "");
  return (PAGES as readonly string[]).includes(name) ? (name as Page) : null;
}

/** Open a page and start at its top. */
export function goTo(navigate: (page: Page) => void, page: Page) {
  navigate(page);
  window.scrollTo({ top: 0 });
}
