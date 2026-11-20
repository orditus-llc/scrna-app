"""Download CellTypist's ready-made models (Phase 0 of the supervisor's plan,
docs/Census_reference_build_plan.md). Saved to data/markers/celltypist/ (never committed or bundled:
the models state no licence; celltypist.org asks users to cite Domínguez Conde et al. 2022, Science,
and Xu et al. 2023, Cell). Also saves the server's model index entry for each (version, date, source)
and each file's SHA-256 in models.json.

    conda activate scrna
    python scripts/download_celltypist_models.py                # every published model (61, ~86 MB)
    python scripts/download_celltypist_models.py Human_Lung_Atlas.pkl   # only the models named

The app offers every human and mouse model that has a broad-type table (scrna_analysis/broad_type_tables/;
scripts/draft_broad_type_tables.py drafts one for a new model).
"""

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

INDEX_URL = "https://celltypist.cog.sanger.ac.uk/models/models.json"
FOLDER = Path(__file__).resolve().parent.parent / "data" / "markers" / "celltypist"


def main() -> None:
    with urllib.request.urlopen(INDEX_URL) as response:
        index = {m["filename"]: m for m in json.load(response)["models"]}
    wanted = sys.argv[1:] or sorted(index)  # all of them by default (analyst, 2026-11-16)
    unknown = [w for w in wanted if w not in index]
    if unknown:
        sys.exit(f"Not in CellTypist's model list: {', '.join(unknown)}")
    FOLDER.mkdir(parents=True, exist_ok=True)
    listing_path = FOLDER / "models.json"
    listing = json.loads(listing_path.read_text()) if listing_path.is_file() else {}
    for name in wanted:
        path = FOLDER / name
        if not path.is_file():
            print(f"Downloading {name} ...", flush=True)
            with urllib.request.urlopen(index[name]["url"]) as response:
                path.write_bytes(response.read())
        entry = {k: index[name].get(k) for k in ("filename", "url", "version", "date", "details", "No_celltypes", "source")}
        entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        listing[name] = entry
        print(f"{name}: version {entry['version']}, {entry['No_celltypes']} cell types, {path.stat().st_size:,} bytes")
    listing_path.write_text(json.dumps(listing, indent=2))
    print(f"Saved to {FOLDER}")


if __name__ == "__main__":
    main()
