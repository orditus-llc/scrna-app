// Unit tests for file names and CSV text (run with `npm test`).
import { describe, expect, it } from "vitest";

import { fileName, slug, toCsv } from "./download";

describe("slug", () => {
  it("makes names safe for files", () => {
    expect(slug("Marker Genes heatmap")).toBe("marker-genes-heatmap");
    expect(slug("UMAP: % mitochondrial")).toBe("umap-mitochondrial");
    expect(slug("Café  cells")).toBe("cafe-cells");
    expect(slug("***")).toBe("file");
  });
  it("drops the data file's extension", () => {
    expect(slug("pbmc3k.h5ad")).toBe("pbmc3k");
    expect(slug("neuron_1k_v3_filtered_feature_bc_matrix.h5")).toBe("neuron-1k-v3-filtered-feature-bc-matrix");
  });
});

describe("fileName", () => {
  it("starts with the dataset", () => {
    expect(fileName("pbmc3k.h5ad", "Marker Genes heatmap", "png")).toBe("pbmc3k_marker-genes-heatmap.png");
  });
  it("works without a dataset", () => {
    expect(fileName(null, "Run log", "csv")).toBe("run-log.csv");
  });
});

describe("toCsv", () => {
  const columns = [
    { key: "gene", header: "Gene" },
    { key: "lfc", header: "Log fold change" },
    { key: "note", header: "Note" },
  ];
  it("writes a header and one line per row, numbers at full precision", () => {
    expect(toCsv(columns, [{ gene: "CD3D", lfc: 2.123456789, note: "" }])).toBe(
      "Gene,Log fold change,Note\r\nCD3D,2.123456789,\r\n",
    );
  });
  it("quotes commas, quotes and line breaks", () => {
    expect(toCsv(columns, [{ gene: "A", lfc: 1, note: 'LTB, CD3D "top"\nnext' }])).toBe(
      'Gene,Log fold change,Note\r\nA,1,"LTB, CD3D ""top""\nnext"\r\n',
    );
  });
  it("leaves missing values blank", () => {
    expect(toCsv(columns, [{ gene: "A", lfc: null, note: undefined }, { gene: "B", lfc: NaN }])).toBe(
      "Gene,Log fold change,Note\r\nA,,\r\nB,,\r\n",
    );
  });
});
