"""Download the second test dataset: 10x Genomics "1k Brain Cells from an E18 Mouse" (v3).

Used only by the tests, to check that the app works on data other than PBMC3k: a different
species (mouse, mitochondrial genes named "mt-"), tissue and file format (10x .h5). Public
example data from 10x Genomics. Saved to data/ (never committed). Users never see it.

    conda activate scrna
    python scripts/download_test_data.py
"""

import hashlib
import urllib.request
from pathlib import Path

URL = ("https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/"
       "neuron_1k_v3_filtered_feature_bc_matrix.h5")
OUTPUT = Path(__file__).resolve().parent.parent / "data" / "neuron_1k_v3_filtered_feature_bc_matrix.h5"


def main() -> None:
    if OUTPUT.exists():
        print(f"Already downloaded: {OUTPUT}")
    else:
        OUTPUT.parent.mkdir(exist_ok=True)
        print(f"Downloading {URL} ...")
        # 10x's server refuses Python's default User-Agent (HTTP 403), so send a plain one.
        request = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request) as response:
            OUTPUT.write_bytes(response.read())
        print(f"Saved to {OUTPUT}")
    print(f"Size: {OUTPUT.stat().st_size:,} bytes")
    print(f"SHA-256: {hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
