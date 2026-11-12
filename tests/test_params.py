"""Tests for `scrna_analysis.params`.

Reference test: the preset must equal the parameters recorded by the validated notebook
in `reference_summary.json`, section by section.
"""

from __future__ import annotations

import dataclasses
import json

import pytest

import inspect

from scrna_analysis.params import (
    DEFAULT_PARAMS,
    SEED,
    AnalysisParams,
    FilterParams,
    PCAParams,
    SuggestionParams,
)
from tests.reference_preset import REFERENCE_PRESET

# Sections of the notebook's PARAMS that map one-to-one onto AnalysisParams fields.
# ("load" holds the input file path, which is an input, not an analysis parameter.)
PARAM_SECTIONS = [
    "qc", "filter", "normalize", "hvg", "regress_out", "scale",
    "pca", "neighbors", "leiden", "umap", "markers",
]

# Scanpy defaults the notebook used implicitly and the app now writes out. They aren't in
# the reference's recorded params, so they're checked separately (each was verified to
# give bit-identical results to the reference).
WRITTEN_OUT_DEFAULTS = {
    "hvg": {"n_top_genes": None},  # Scanpy's default (use the cut-offs); added 2026-10-26
    "pca": {"svd_solver": "arpack"},
    "leiden": {"directed": True},
    "markers": {"reference": "rest", "corr_method": "benjamini-hochberg", "tie_correct": False},
}


# --- Reference tests (PBMC3k) ----------------------------------------------------------

@pytest.mark.parametrize("section", PARAM_SECTIONS)
def test_preset_matches_reference_params(section, reference_summary):
    preset = REFERENCE_PRESET.to_dict()[section]
    for name in WRITTEN_OUT_DEFAULTS.get(section, {}):
        preset.pop(name)
    assert preset == reference_summary["params"][section]


@pytest.mark.parametrize("section", WRITTEN_OUT_DEFAULTS)
def test_written_out_defaults(section):
    preset = REFERENCE_PRESET.to_dict()[section]
    for name, value in WRITTEN_OUT_DEFAULTS[section].items():
        assert preset[name] == value


def test_preset_covers_all_reference_sections(reference_summary):
    assert set(reference_summary["params"]) - {"load"} == set(PARAM_SECTIONS)


def test_seed_matches_reference(reference_summary):
    assert SEED == reference_summary["seed"]
    params = REFERENCE_PRESET
    seeds = [params.pca.random_state, params.neighbors.random_state,
             params.leiden.random_state, params.umap.random_state]
    assert seeds == [SEED] * 4


# --- General behavior ------------------------------------------------------------------

def test_annotation_has_no_species_or_tissue_by_default():
    """Nothing about the dataset is assumed: species and tissues are chosen per dataset."""
    m = DEFAULT_PARAMS.annotation.markers
    assert m.species is None and m.tissues == ()


def test_to_dict_is_json_serializable():
    as_dict = REFERENCE_PRESET.to_dict()
    assert json.loads(json.dumps(as_dict)) == as_dict


def test_params_are_read_only():
    with pytest.raises(dataclasses.FrozenInstanceError):
        REFERENCE_PRESET.pca.n_comps = 50  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        REFERENCE_PRESET.annotation.markers.tissues = ("Brain",)  # type: ignore[misc]


def test_override_one_value_leaves_preset_unchanged():
    custom = dataclasses.replace(REFERENCE_PRESET, pca=PCAParams(n_comps=50))
    assert custom.pca.n_comps == 50
    assert REFERENCE_PRESET.pca.n_comps == 40
    assert custom.leiden == REFERENCE_PRESET.leiden


# --- The app's general defaults and their named sources -------------------------------------

def _default(function, name):
    return inspect.signature(function).parameters[name].default


def test_defaults_are_scanpys_own_where_params_py_says_so():
    """The defaults documented as "Source: Scanpy's default" are Scanpy's actual defaults in
    the installed version, so a Scanpy update that changes one is noticed."""
    import scanpy as sc
    d = AnalysisParams()
    assert d.normalize.target_sum == _default(sc.pp.normalize_total, "target_sum") is None
    hvg = sc.pp.highly_variable_genes
    assert (d.hvg.flavor, d.hvg.min_mean, d.hvg.max_mean, d.hvg.min_disp) == (
        _default(hvg, "flavor"), _default(hvg, "min_mean"), _default(hvg, "max_mean"), _default(hvg, "min_disp"))
    assert d.pca.n_comps == sc.settings.N_PCS and _default(sc.pp.pca, "n_comps") is None
    assert d.neighbors.n_neighbors == _default(sc.pp.neighbors, "n_neighbors")
    assert d.neighbors.n_pcs == d.pca.n_comps  # n_pcs=None: all computed components
    leiden = sc.tl.leiden
    assert (d.leiden.flavor, d.leiden.resolution, d.leiden.n_iterations) == (
        _default(leiden, "flavor"), _default(leiden, "resolution"), _default(leiden, "n_iterations"))
    assert SEED == _default(sc.pp.pca, "random_state") == _default(sc.tl.umap, "random_state")
    assert d.regress_out.enabled is False  # best-practices guide: off unless needed
    assert d.scale.max_value == 10  # Seurat's ScaleData default (scale.max = 10)


def test_default_params_suggest_the_qc_thresholds():
    assert DEFAULT_PARAMS == AnalysisParams(suggest=SuggestionParams())
    assert AnalysisParams().filter == FilterParams(min_genes=0, min_cells=3)  # filters nothing before QC


def test_reference_preset_differs_from_the_defaults_only_where_documented():
    """tests/reference_preset.py: the notebook's own choices."""
    differs = {f.name for f in dataclasses.fields(AnalysisParams)
               if getattr(REFERENCE_PRESET, f.name) != getattr(AnalysisParams(), f.name)}
    assert differs == {"filter", "normalize", "regress_out", "pca", "neighbors", "annotation"}
