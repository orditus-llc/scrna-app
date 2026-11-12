// Saving charts and tables as files: file names, CSV text, and the browser download itself.

/** A file-name-safe version of a name: lower case, words joined by "-" ("Marker Genes" -> "marker-genes"). */
export function slug(text: string): string {
  return (
    text
      .normalize("NFKD")
      .replace(/[̀-ͯ]/g, "") // accents
      .toLowerCase()
      .replace(/\.(h5ad|h5)$/, "") // the data file's extension
      .replace(/[^a-z0-9]+/g, "-")
      .replace(/^-+|-+$/g, "") || "file"
  );
}

/** "<dataset>_<what>.<ext>", e.g. fileName("pbmc3k.h5ad", "Marker Genes heatmap", "png") ->
 * "pbmc3k_marker-genes-heatmap.png". Without a dataset, just "<what>.<ext>". */
export function fileName(dataset: string | null, what: string, extension: string): string {
  return `${dataset ? `${slug(dataset)}_` : ""}${slug(what)}.${extension}`;
}

export type CsvValue = string | number | boolean | null | undefined;

/** Rows as CSV text (RFC 4180): a header line, then one line per row. Values with a comma,
 * quote or line break are quoted; empty values (null, undefined, NaN) are left blank. Numbers
 * keep their full precision. */
export function toCsv(columns: { key: string; header: string }[], rows: Record<string, CsvValue>[]): string {
  const cell = (v: CsvValue): string => {
    if (v === null || v === undefined || (typeof v === "number" && Number.isNaN(v))) return "";
    const text = String(v);
    return /[",\r\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
  };
  const lines = [columns.map((c) => cell(c.header)).join(",")];
  for (const row of rows) lines.push(columns.map((c) => cell(row[c.key])).join(","));
  return lines.join("\r\n") + "\r\n";
}

/** Hands a file to the browser to save (a data: URL or a Blob). */
export function saveFile(name: string, content: Blob | string): void {
  const url = typeof content === "string" ? content : URL.createObjectURL(content);
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  document.body.appendChild(link);
  link.click();
  link.remove();
  if (typeof content !== "string") window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/** Saves CSV text. A byte-order mark lets Excel read accents and "µ" correctly. */
export function saveCsv(name: string, csv: string): void {
  saveFile(name, new Blob(["﻿", csv], { type: "text/csv;charset=utf-8" }));
}
