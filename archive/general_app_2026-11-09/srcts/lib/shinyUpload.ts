// Send a file to the server with Shiny's own upload system (the one behind Shiny's standard
// file button), so Shiny stores it in the session's upload folder and deletes it when the
// session ends. Three parts, as in Shiny's file input:
//   1. "uploadInit": tell the server the file's name, size and type; it answers with an
//      upload address,
//   2. POST the file's bytes to that address (progress is reported as it goes),
//   3. "uploadEnd": the server sets the input `inputId` to [{name, size, type, datapath}].
// shinyreact has no upload component, so this uses Shiny's in-page API (window.Shiny.shinyapp),
// which Shiny's own file button uses too. A browser check covers it.

interface ShinyApp {
  makeRequest(
    method: string,
    args: unknown[],
    onSuccess: (value: unknown) => void,
    onError: (message: string) => void,
    blobs: undefined,
  ): void;
}

declare global {
  interface Window {
    Shiny?: { shinyapp?: ShinyApp };
  }
}

function request<T>(method: string, args: unknown[]): Promise<T> {
  const app = window.Shiny?.shinyapp;
  if (!app) return Promise.reject(new Error("Not connected to the analysis server."));
  return new Promise((resolve, reject) =>
    app.makeRequest(method, args, (value: unknown) => resolve(value as T), (message: string) => reject(new Error(message)), undefined),
  );
}

/** Upload `file` to the server as the value of input `inputId`. `onProgress` gets 0–1. */
export async function uploadToShiny(file: File, inputId: string, onProgress: (fraction: number) => void): Promise<void> {
  const { jobId, uploadUrl } = await request<{ jobId: string; uploadUrl: string }>("uploadInit", [
    [{ name: file.name, size: file.size, type: file.type }],
  ]);
  await new Promise<void>((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.open("POST", uploadUrl);
    xhr.setRequestHeader("Content-Type", "application/octet-stream");
    xhr.upload.onprogress = (e) => e.lengthComputable && onProgress(e.loaded / e.total);
    xhr.onload = () => (xhr.status >= 200 && xhr.status < 300 ? resolve() : reject(new Error(`Upload failed (${xhr.status}).`)));
    xhr.onerror = () => reject(new Error("Upload failed: the connection was interrupted."));
    xhr.send(file);
  });
  onProgress(1);
  await request("uploadEnd", [jobId, inputId]);
}

/** File endings the app accepts (app/payloads.py and scrna_analysis/load_data.py FORMATS). */
export const ACCEPTED_FILES = [".h5ad", ".h5"];

/** A file name's ending, if accepted; otherwise null. */
export function acceptedEnding(name: string): string | null {
  const lower = name.toLowerCase();
  return ACCEPTED_FILES.find((ending) => lower.endsWith(ending)) ?? null;
}
