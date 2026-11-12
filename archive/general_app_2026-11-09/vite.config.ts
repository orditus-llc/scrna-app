// Builds the React interface (srcts/ui.tsx) into www/ui.js + www/ui.css, which the Shiny
// app serves (see app.py, page_react). Setup follows shinyreact's "TSX files and
// JavaScript build tools" guide.
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";
import { defineConfig } from "vite";

export default defineConfig({
  plugins: [react(), tailwindcss()],
  define: { "process.env.NODE_ENV": JSON.stringify("production") },
  build: {
    outDir: "www",
    emptyOutDir: false,
    // Embed fonts in ui.css so the app works offline and www/ holds only two files.
    assetsInlineLimit: 1_000_000,
    lib: {
      entry: "srcts/ui.tsx",
      formats: ["iife"],
      name: "PbmcAnalysisUI",
      fileName: () => "ui.js",
    },
    rollupOptions: {
      // React and the Shiny hooks come from shinyreact's bundle, which the server already
      // loads (window.shinyreact). Never ship a second copy: hooks would break. The npm
      // package @posit-dev/shinyreact is used only for its TypeScript types.
      external: ["react", "react-dom", "react-dom/client", "@posit-dev/shinyreact"],
      output: {
        globals: {
          react: "window.shinyreact.React",
          "react-dom": "window.shinyreact.ReactDOM",
          "react-dom/client": "window.shinyreact.ReactDOM",
          "@posit-dev/shinyreact": "window.shinyreact",
        },
        assetFileNames: "ui.[ext]",
      },
    },
  },
});
