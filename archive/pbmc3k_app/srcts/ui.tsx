// Entry point: mounts the React interface into the page the Shiny app serves.
import { createRoot } from "react-dom/client";

import { App } from "./App";
import "./theme.css";

const mount = document.body.appendChild(document.createElement("div"));
createRoot(mount).render(<App />);
