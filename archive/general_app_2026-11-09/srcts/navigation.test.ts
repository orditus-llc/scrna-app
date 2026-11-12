// Unit tests for page addresses (run with `npm test`).
import { describe, expect, it } from "vitest";

import { PAGES, pageFromHash, STEP_KEYS } from "./navigation";

describe("pageFromHash", () => {
  it("reads every page from its address", () => {
    for (const page of PAGES) expect(pageFromHash(`#${page}`)).toBe(page);
  });

  it("covers all seven steps and Explore Results", () => {
    expect(PAGES).toEqual([...STEP_KEYS, "explorer"]);
  });

  it("returns null for an empty or unknown address", () => {
    expect(pageFromHash("")).toBeNull();
    expect(pageFromHash("#")).toBeNull();
    expect(pageFromHash("#overview")).toBeNull(); // the Overview page was removed
    expect(pageFromHash("#nonsense")).toBeNull();
  });
});
