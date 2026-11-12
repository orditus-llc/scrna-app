"""Download the multiple-samples test dataset: Kang et al. 2018 (Nature Biotechnology 36, 89–94).

Blood cells (PBMCs) from 8 lupus patients, each sample split in two and cultured for 6 hours with
or without interferon-beta: 16 samples, 2 conditions. The standard example for analysing several
samples together (used by Seurat's and Scanpy's integration tutorials). Original data: GEO
GSE96583. This copy is the one scverse (the makers of Scanpy) provide for pertpy's
`kang_2018()`. Used only by the tests. Saved to data/ (never committed).

    conda activate scrna
    python scripts/download_kang_data.py
"""

import hashlib
import urllib.request
from pathlib import Path

URL = "https://exampledata.scverse.org/pertpy/kang_2018.h5ad"
OUTPUT = Path(__file__).resolve().parent.parent / "data" / "kang_2018.h5ad"


def main() -> None:
    if OUTPUT.exists():
        print(f"Already downloaded: {OUTPUT}")
    else:
        OUTPUT.parent.mkdir(exist_ok=True)
        print(f"Downloading {URL} ...")
        request = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request) as response:
            OUTPUT.write_bytes(response.read())
        print(f"Saved to {OUTPUT}")
    print(f"Size: {OUTPUT.stat().st_size:,} bytes")
    print(f"SHA-256: {hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
