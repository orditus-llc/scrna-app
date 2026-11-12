"""Download the Cell Ontology (CL, "basic" edition): the standard family tree of cell types
("naive B cell" is a kind of "B cell", which is a kind of "lymphocyte", …) that CELLxGENE Census
uses for its cell type labels. The Census marker tables use it to compare each cell type only with
cell types it isn't related to. Licence: CC BY 4.0 (cite the Cell Ontology). Saved to data/markers/
(never committed).

    conda activate scrna
    python scripts/download_cell_ontology.py
"""

import hashlib
import urllib.request
from pathlib import Path

URL = "http://purl.obolibrary.org/obo/cl/cl-basic.obo"
OUTPUT = Path(__file__).resolve().parent.parent / "data" / "markers" / "cl-basic.obo"


def main() -> None:
    if OUTPUT.exists():
        print(f"Already downloaded: {OUTPUT}")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {URL} ...")
        request = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request) as response:
            OUTPUT.write_bytes(response.read())
        print(f"Saved to {OUTPUT}")
    version = next((line.split(":", 1)[1].strip() for line in OUTPUT.read_text().splitlines()[:30]
                    if line.startswith("data-version:")), "unknown")
    print(f"Version: {version}")
    print(f"Size: {OUTPUT.stat().st_size:,} bytes")
    print(f"SHA-256: {hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
