"""The multi-sample path checked against an independent reference: `reference/kang-multisample-validated.ipynb`.

That notebook was written in a separate session from a written spec (`docs/kang_validation_spec.md`),
without seeing the app's code, and saved its answers to `validation/expected_multisample/`. Here the app
runs the same analysis (Kang et al. 2018, 16 samples, default settings) and is compared step by step.
Values are read from the saved files, never retyped.

Tolerances (CLAUDE.md, "Strict test tolerances"; the reason for each looser check is written here):
- Exact: per-sample QC limits, kept cells and genes, highly variable genes, mixing scores, clusters and
  marker genes without batch correction, and every step when it is given the notebook's own input.
- PCA within 1e-6. The app's and the notebook's PCA differ by ~2e-7 (float32 rounding): the notebook
  subsets its genes without Scanpy's `highly_variable` column, so `sc.pp.pca` takes a different internal
  path. Neither is wrong; the app's path is the one that is bit-identical to the PBMC3k reference.
- That ~2e-7 difference reaches the end of the Harmony run: Harmony enlarges it to ~1e-5 and ~2% of cells
  then fall in a neighboring Leiden cluster (found 2026-11-12). So the Harmony run is checked (a) exactly,
  step by step from the notebook's own saved coordinates, and (b) end to end with ARI >= 0.95.
- UMAP is checked exactly from the notebook's PCA. End to end it is a picture that magnifies tiny graph
  differences (4 of 584,438 graph entries differ), so only its overall agreement is checked.
"""

from __future__ import annotations

import json
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from scrna_analysis import integration, markers as marker_tools
from scrna_analysis.params import DEFAULT_PARAMS, IntegrationParams
from scrna_analysis.pipeline import AnalysisRun
from tests.conftest import PROJECT_ROOT

EXPECTED = PROJECT_ROOT / "validation" / "expected_multisample"
HINT = ("Regenerate it (with approval) by executing the Kang notebook:\n"
        "  jupyter nbconvert --to notebook --execute --inplace reference/kang-multisample-validated.ipynb")


@pytest.fixture(scope="module")
def expected():
    for name in ("kang_summary.json", "kang_results.h5ad"):
        if not (EXPECTED / name).is_file():
            pytest.fail(f"Missing reference file: validation/expected_multisample/{name}\n{HINT}")
    return json.loads((EXPECTED / "kang_summary.json").read_text()), ad.read_h5ad(EXPECTED / "kang_results.h5ad")


# The app's defaults pick highly variable genes within each sample (the notebook's "per_sample"); the
# notebook's "pooled" variant is the same settings with that switched off.
POOLED = replace(DEFAULT_PARAMS, hvg=replace(DEFAULT_PARAMS.hvg, per_sample=False))


def _kang(kang_path, params):
    """The app on Kang up to Preprocessing (samples from the file's columns)."""
    run = AnalysisRun(kang_path, params, sample_columns={"sample": "replicate", "condition": "label"},
                      metadata={"species": "human", "tissue": "blood"})
    for stage in ("load", "compute_qc_metrics", "apply_filter", "preprocess"):
        getattr(run, stage)()
    return run


@pytest.fixture(scope="module")
def kang_run(kang_path):
    """Pooled highly variable genes (the notebook's "pooled" runs)."""
    return _kang(kang_path, POOLED)


@pytest.fixture(scope="module")
def kang_run_per_sample(kang_path):
    """The app's default settings: genes picked within each sample (the notebook's "per_sample" runs)."""
    return _kang(kang_path, DEFAULT_PARAMS)


def _downstream(run, integ):
    run.reduce_dimensions(integration_params=integ)
    run.cluster()
    run.find_markers()
    return {"reduced": run.stage_result("reduce_dimensions"), "clustered": run.result("cluster"),
            "markers": run.result("markers")}


@pytest.fixture(scope="module")
def runs(kang_run, kang_run_per_sample):
    """Dimensionality Reduction, Clustering and Marker Genes without and with Harmony (by sample)."""
    out = {name: _downstream(kang_run, integ)
           for name, integ in (("pooled_none", IntegrationParams()), ("pooled_harmony", IntegrationParams(method="harmony")))}
    out["per_sample_none"] = _downstream(kang_run_per_sample, IntegrationParams())
    return out


# --- Loading and Quality Control (exact) ---------------------------------------------------------

def test_samples_and_qc_limits_match(kang_run, expected):
    summary, _ = expected
    assert [s["name"] for s in kang_run.samples] == [s["name"] for s in summary["samples"]]
    assert [s["condition"] for s in kang_run.samples] == [s["condition"] for s in summary["samples"]]
    assert [s["n_cells"] for s in kang_run.samples] == [s["cells_loaded"] for s in summary["samples"]]
    app = {lim.sample: lim for lim in kang_run.params.filter.sample_limits}
    for ref in summary["qc_limits"]:
        lim = app[ref["sample"]]
        assert (lim.min_genes, lim.max_n_genes_by_counts, lim.max_pct_counts_mt) == (
            ref["min_genes"], ref["max_genes"], ref["max_pct_mt"])
    details = kang_run.stage_result("qc_metrics").summary["suggested_filter"]["details"]["by_sample"]
    for ref in summary["qc_limits"]:
        assert details[ref["sample"]]["min_genes"]["median"] == pytest.approx(ref["log1p_genes_median"], abs=1e-12)
        assert details[ref["sample"]]["min_genes"]["mad"] == pytest.approx(ref["log1p_genes_mad"], abs=1e-12)


def test_kept_cells_and_genes_match(kang_run, expected):
    summary, ref = expected
    filtered = kang_run.result("filter")
    assert list(filtered.obs_names) == list(ref.obs_names)
    assert filtered.n_vars == summary["n_genes_kept"]
    kept = filtered.obs["sample"].astype(str).value_counts()
    for row in summary["cells_kept_per_sample"]:
        assert kept[row["sample"]] == row["cells_after"]


def test_normalization_and_highly_variable_genes_match(kang_run, expected):
    summary, _ = expected
    pre = kang_run.result("preprocess")
    assert sorted(map(str, pre.var_names)) == sorted(summary["hvg"]["pooled"])
    totals = np.asarray(kang_run.result("filter").X.sum(axis=1)).ravel()
    assert float(np.median(totals)) == summary["normalize_target_sum_used"]
    assert kang_run.stage_result("preprocess").summary["hvg_per_sample"] is False


def test_highly_variable_genes_within_each_sample_match(kang_run_per_sample, expected):
    """The app's default with several samples: genes picked within each sample, exactly the notebook's list."""
    summary, _ = expected
    pre = kang_run_per_sample.result("preprocess")
    assert sorted(map(str, pre.var_names)) == sorted(summary["hvg"]["per_sample"])
    assert len(summary["hvg"]["per_sample"]) == summary["hvg"]["n_per_sample"]
    assert kang_run_per_sample.stage_result("preprocess").summary["hvg_per_sample"] is True
    assert pre.var["highly_variable_nbatches"].max() <= len(kang_run_per_sample.samples)
    assert kang_run_per_sample.record()["params"]["hvg"]["per_sample"] is True


# --- Without batch correction ----------------------------------------------------------------------

def aligned(a, b):
    """`a` with each column's sign flipped to match `b` (PC signs are arbitrary)."""
    return a * np.sign((a * b).sum(axis=0))


def test_pca_and_mixing_match(runs, expected):
    summary, ref = expected
    for name in ("pooled_none", "pooled_harmony", "per_sample_none"):
        reduced = runs[name]["reduced"]
        variance = np.asarray(reduced.details[0].variance_ratio)
        assert np.allclose(variance, summary["runs"][name]["pca_variance_ratio"], rtol=0, atol=1e-6)
        pca = reduced.adata.obsm["X_pca"]
        assert np.abs(aligned(pca, ref.obsm[f"X_pca_{name}"]) - ref.obsm[f"X_pca_{name}"]).max() < 1e-6
        before = reduced.adata.obs.get(integration.MIXING_BEFORE_KEY, reduced.adata.obs[integration.MIXING_KEY])
        assert np.allclose(before.to_numpy(), ref.obs[f"mixing_before_{name}"].to_numpy(), rtol=0, atol=1e-9)
        assert reduced.summary["integration"]["mixing_before"] == pytest.approx(
            summary["runs"][name]["mixing_before_mean"], abs=1e-9)


@pytest.mark.parametrize("name", ["pooled_none", "per_sample_none"])
def test_clusters_markers_and_composition_match_without_correction(runs, expected, name):
    summary, ref = expected
    r = summary["runs"][name]
    clustered = runs[name]["clustered"]
    assert adjusted_rand_score(clustered.obs["leiden"], ref.obs[f"leiden_{name}"]) == 1.0
    assert {str(k): v for k, v in clustered.obs["leiden"].value_counts().items()} == r["cluster_sizes"]
    table = marker_tools.marker_table(runs[name]["markers"], n_top=25)
    assert {str(g): list(rows["gene"]) for g, rows in table.groupby("group")} == r["top25_markers"]
    comp = integration.composition(clustered, "sample", "leiden").set_index(["sample", "group"])
    for row in r["composition"]:
        if row["n_cells"]:
            assert comp.loc[(row["sample"], row["cluster"]), "n_cells"] == row["n_cells"]
            assert comp.loc[(row["sample"], row["cluster"]), "share"] == pytest.approx(row["share"], abs=1e-12)


# --- With Harmony ------------------------------------------------------------------------------------

@pytest.mark.parametrize("variant", ["pooled", "per_sample"])
def test_harmony_step_by_step_from_the_notebooks_own_coordinates(kang_run, kang_run_per_sample, expected, variant):
    """Each step exactly, given the previous step's saved answer: Harmony from the notebook's PCA,
    the mixing score from its corrected coordinates, and clustering from those coordinates."""
    _, ref = expected
    name = f"{variant}_harmony"
    a = (kang_run if variant == "pooled" else kang_run_per_sample).result("preprocess").copy()
    a.obsm["X_pca"] = ref.obsm[f"X_pca_{name}"].copy()
    out, _ = integration.run_harmony(a, IntegrationParams(method="harmony"), n_components=50)
    assert np.array_equal(out.obsm[integration.CORRECTED_KEY], ref.obsm[f"X_pca_harmony_{name}"])
    after = integration.sample_mixing(out.obsm[integration.CORRECTED_KEY], out.obs["sample"])
    assert np.allclose(after, ref.obs[f"mixing_after_{name}"].to_numpy(), rtol=0, atol=1e-9)
    sc.pp.neighbors(out, n_neighbors=15, n_pcs=50, use_rep=integration.CORRECTED_KEY, random_state=0)
    from scrna_analysis.clustering import cluster_leiden
    clustered, _ = cluster_leiden(out, DEFAULT_PARAMS.leiden)
    assert adjusted_rand_score(clustered.obs["leiden"], ref.obs[f"leiden_{name}"]) == 1.0


def test_harmony_end_to_end_agrees_closely(runs, expected):
    """End to end the ~2e-7 PCA difference (module docstring) moves a few boundary cells."""
    summary, ref = expected
    reduced, clustered = runs["pooled_harmony"]["reduced"], runs["pooled_harmony"]["clustered"]
    corrected = reduced.adata.obsm[integration.CORRECTED_KEY]
    assert np.abs(corrected - ref.obsm["X_pca_harmony_pooled_harmony"]).max() < 1e-4
    assert reduced.summary["integration"]["mixing_after"] == pytest.approx(
        summary["runs"]["pooled_harmony"]["mixing_after_mean"], abs=1e-3)
    assert adjusted_rand_score(clustered.obs["leiden"], ref.obs["leiden_pooled_harmony"]) >= 0.95
    assert clustered.obs["leiden"].nunique() == summary["runs"]["pooled_harmony"]["n_clusters"]


# --- UMAP --------------------------------------------------------------------------------------------

def test_umap_from_the_notebooks_pca_is_identical(kang_run, expected):
    _, ref = expected
    from scrna_analysis.dim_reduction import compute_neighbors, run_umap
    a = kang_run.result("preprocess").copy()
    a.obsm["X_pca"] = ref.obsm["X_pca_pooled_none"].copy()
    a, _ = compute_neighbors(a, DEFAULT_PARAMS.neighbors)
    a, _ = run_umap(a, DEFAULT_PARAMS.umap)
    assert np.array_equal(a.obsm["X_umap"], ref.obsm["X_umap_pooled_none"])


def test_umap_end_to_end_has_the_same_shape(runs, expected):
    _, ref = expected
    mine, theirs = runs["pooled_none"]["reduced"].adata.obsm["X_umap"], ref.obsm["X_umap_pooled_none"]
    for axis in (0, 1):
        assert abs(np.corrcoef(mine[:, axis], theirs[:, axis])[0, 1]) > 0.95


# --- The notebook's biological checks (recorded there; read back here) ---------------------------------

def test_the_reference_passed_its_biological_checks(expected):
    summary, _ = expected
    table = pd.DataFrame(summary["interferon"]["table"])
    assert len(table) > 0 and (table["stim_minus_ctrl"] > 0).all()  # interferon genes up in every cell type
    assert summary["warnings"] == []
