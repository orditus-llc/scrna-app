// "Download results": the analysed data as one .h5ad file and the run report (HTML), made by
// the server from the current analysis. Shown in Explore Results and in the Run record panel.
import { Download, FileText } from "lucide-react";

import { serverDownloadUrl } from "../lib/serverDownload";
import { cn } from "../lib/utils";
import type { ExportInfo } from "../types";

export function ExportButtons({ info }: { info: ExportInfo }) {
  const link = (name: string, file: string, label: string, icon: React.ReactNode, primary: boolean) => {
    const url = info.available ? serverDownloadUrl(name) : null;
    const style = cn(
      "inline-flex items-center gap-2 rounded-[4px] px-3 py-[7px] text-[13px] font-bold no-underline",
      primary ? "bg-ink text-white hover:bg-[#2a2a2e]" : "border border-field bg-white text-ink hover:bg-chip",
    );
    return url ? (
      <a href={url} download={file} className={style}>
        {icon} {label}
      </a>
    ) : (
      <span aria-disabled="true" className="inline-flex items-center gap-2 rounded-[4px] border border-line bg-chip px-3 py-[7px] text-[13px] font-bold text-off">
        {icon} {label}
      </span>
    );
  };
  return (
    <section aria-label="Download results" className="flex flex-col gap-2">
      <div className="flex flex-wrap gap-2">
        {link("download_h5ad", info.h5ad_name, "Download analysed data (.h5ad)", <Download aria-hidden="true" size={14} />, true)}
        {link("download_report", info.report_name, "Download run report", <FileText aria-hidden="true" size={14} />, false)}
      </div>
      <span className="text-[12.5px] leading-snug text-muted">
        {info.available
          ? `Includes: ${info.includes.join(", ")}${info.reviewed ? ", your reviewed labels" : ""}. The .h5ad opens in Scanpy, CELLxGENE and Galaxy; the report holds the settings, a Methods paragraph and the run log.`
          : info.reason}
      </span>
    </section>
  );
}
