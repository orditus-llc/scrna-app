"""Download the PanglaoDB cell-type marker database (release of 27 March 2020).

PanglaoDB (Franzén et al. 2019, Database, https://panglaodb.se) lists marker genes for cell
types of human and mouse; Cell Type Annotation uses it to suggest a cell type for each
cluster. The database states no licence, so the file is downloaded here rather than bundled
with the app. Saved to data/markers/ (never committed).

    conda activate scrna
    python scripts/download_marker_db.py
"""

import hashlib
import urllib.request
from pathlib import Path

URL = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
OUTPUT = Path(__file__).resolve().parent.parent / "data" / "markers" / "PanglaoDB_markers_27_Mar_2020.tsv.gz"
# Fingerprint of the file when it was first downloaded (2026-10-28), to notice if it changes.
EXPECTED_SHA256 = "6779952ad40aa5a124de7bd0e18975c6630bd6006d6b3ef210a916caaa6b53c9"


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
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print(f"Size: {OUTPUT.stat().st_size:,} bytes")
    print(f"SHA-256: {digest}")
    if digest != EXPECTED_SHA256:
        print(f"Warning: this isn't the file first downloaded (SHA-256 {EXPECTED_SHA256}). "
              "PanglaoDB may have changed it; check before using it.")


if __name__ == "__main__":
    main()
