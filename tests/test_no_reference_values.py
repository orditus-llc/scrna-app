"""Guard: the backend must not contain expected PBMC3k results.

Backend functions compute results from the data. Expected values (cell counts, HVG
count, cluster sizes, reviewed labels, cluster-number -> label maps) belong only in
`validation/` and `tests/`. This scan is a safety net, not a replacement for review.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parent.parent / "scrna_analysis"

# Numbers that are PBMC3k *results* (not parameters): raw and filtered shapes, HVG count,
# and the larger cluster sizes. Small results (13 mito genes, 8 clusters, the 36- and
# 13-cell clusters) aren't searched for, because those numbers occur innocently in code.
REFERENCE_NUMBERS = ["2700", "32738", "2670", "13714", "1846",
                     "1131", "430", "342", "326", "214", "178"]

# Analyst-reviewed working labels (validation/annotation_review.md).
REVIEWED_LABELS = [
    "CD4-enriched T cells",
    "Cytotoxic/memory T cells",
    "Conventional dendritic cells",
    "Platelet-lineage cells",
]

# A dictionary keyed by cluster numbers, e.g. {"0": "CD4+ T cells", ...}.
CLUSTER_ID_MAP = re.compile(r"""['"]\d+['"]\s*:\s*['"][A-Za-z]""")


APP = BACKEND.parent / "app"


def _backend_files() -> list[Path]:
    """Backend modules, the app's Python files and the React source
    (app/srcts/). None of them may contain reference results."""
    app_python = [p for p in APP.rglob("*.py") if "node_modules" not in p.parts]
    react_source = list((APP / "srcts").rglob("*.ts")) + list((APP / "srcts").rglob("*.tsx"))
    files = sorted(BACKEND.glob("*.py")) + sorted(app_python) + sorted(react_source)
    assert files, f"No backend files found in {BACKEND}"
    return files


def find_reference_values(text: str) -> list[str]:
    """Describe every reference-only value found in `text`."""
    problems = []
    for number in REFERENCE_NUMBERS:
        if re.search(rf"(?<![\d.]){number}(?![\d.])", text):
            problems.append(f"reference number {number}")
    for label in REVIEWED_LABELS:
        if label in text:
            problems.append(f"reviewed label '{label}'")
    if CLUSTER_ID_MAP.search(text):
        problems.append("cluster-number -> label mapping")
    return problems


@pytest.mark.parametrize("path", _backend_files(), ids=lambda p: p.name)
def test_no_reference_results_in_backend(path):
    problems = find_reference_values(path.read_text())
    assert not problems, f"{path.name} contains reference-only values: {problems}"


def test_guard_detects_planted_values():
    """The scan must actually catch violations (otherwise passing means nothing)."""
    bad = (
        'assert adata.n_obs == 2670\n'
        'labels = {"0": "CD4+ T cells", "1": "B cells"}\n'
        'name = "Platelet-lineage cells"\n'
    )
    assert len(find_reference_values(bad)) == 3
    # Parameters and ordinary numbers are allowed.
    ok = "max_n_genes_by_counts: int = 2500\nn_comps: int = 40\ntarget_sum = 1e4\n"
    assert find_reference_values(ok) == []
