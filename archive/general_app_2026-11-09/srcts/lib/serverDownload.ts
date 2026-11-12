// Files the server makes on request (the analysed .h5ad, the run report). Shiny serves each
// registered download (app.py, @render.download_button) at session/<session id>/download/<name>;
// shinyreact has no download control, so the page links to that address itself.

/** The address of the server's download `name` for this browser session, or null before the
 * connection is ready. */
export function serverDownloadUrl(name: string): string | null {
  const shiny = (window as { Shiny?: { shinyapp?: { config?: { sessionId?: string } } } }).Shiny;
  const session = shiny?.shinyapp?.config?.sessionId;
  return session ? `session/${encodeURIComponent(session)}/download/${encodeURIComponent(name)}?w=` : null;
}
