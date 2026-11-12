// Unit tests for the file check before uploading (run with `npm test`).
import { describe, expect, it } from "vitest";

import { acceptedEnding } from "./shinyUpload";

describe("acceptedEnding", () => {
  it("accepts .h5ad and 10x .h5 files, in any letter case", () => {
    expect(acceptedEnding("sample.h5ad")).toBe(".h5ad");
    expect(acceptedEnding("filtered_feature_bc_matrix.h5")).toBe(".h5");
    expect(acceptedEnding("SAMPLE.H5AD")).toBe(".h5ad");
  });

  it("rejects other files", () => {
    expect(acceptedEnding("counts.csv")).toBeNull();
    expect(acceptedEnding("matrix.mtx.gz")).toBeNull();
    expect(acceptedEnding("sample.h5ad.zip")).toBeNull();
  });
});
