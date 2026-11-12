import { describe, expect, it } from "vitest";

import { sampleNameFromFile, sampleTableProblem, type PendingSample } from "./Samples";

const sample = (name: string): PendingSample => ({ file: new File([""], `${name}.h5ad`), name, condition: "" });

describe("naming several samples", () => {
  it("starts each name from the file name", () => {
    expect(sampleNameFromFile("donor_1.h5ad")).toBe("donor_1");
    expect(sampleNameFromFile("run2_filtered.H5")).toBe("run2_filtered");
  });

  it("explains what is wrong before loading", () => {
    expect(sampleTableProblem([sample("a")], 20)).toBe("Choose at least two files for several samples.");
    expect(sampleTableProblem([sample("a"), sample(" ")], 20)).toBe("Every sample needs a name.");
    expect(sampleTableProblem([sample("a"), sample("a ")], 20)).toBe('Sample names must be different ("a" is used twice).');
    expect(sampleTableProblem([sample("a"), sample("b"), sample("c")], 2)).toBe("The app handles up to 2 samples for now.");
    expect(sampleTableProblem([sample("a"), sample("b")], 20)).toBeNull();
  });
});
