"""Build CELLxGENE Census marker tables (one per species and Census tissue), ahead of time.

Pulls reference cells from the Census API (needs the internet), computes each cell type's marker
genes, and saves a small table in data/markers/census/<release>/ with its provenance. The app then
reads the table instantly. See scrna_analysis/census_markers.py for how. Census data: CC BY 4.0
(cite the Census).

    conda activate scrna
    python scripts/build_census_markers.py human blood
    python scripts/build_census_markers.py mouse brain
    python scripts/build_census_markers.py human blood --cells-per-type 300
    python scripts/build_census_markers.py --tissues     # only save Census's list of tissues
"""

import argparse
import sys
import time
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scrna_analysis.census_markers import build_marker_table, fetch_tissue_list  # noqa: E402
from scrna_analysis.params import CensusMarkerParams  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("species", nargs="?", choices=["human", "mouse"])
    parser.add_argument("tissue", nargs="?", help="a Census tissue_general name, e.g. blood, brain, spleen, lung")
    parser.add_argument("--tissues", action="store_true", help="only save Census's list of tissues")
    parser.add_argument("--cells-per-type", type=int, default=CensusMarkerParams.cells_per_type)
    args = parser.parse_args()
    params = replace(CensusMarkerParams(), cells_per_type=args.cells_per_type)
    print(f"Tissue list: {fetch_tissue_list(params.release)}", flush=True)  # small; always refreshed
    if args.tissues:
        return
    if not (args.species and args.tissue):
        parser.error("give a species and a tissue, or --tissues")
    start = time.time()
    path = build_marker_table(args.species, args.tissue, params, log=lambda m: print(m, flush=True))
    print(f"Done in {(time.time() - start) / 60:.1f} minutes: {path}")


if __name__ == "__main__":
    main()
