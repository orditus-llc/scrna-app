"""Compare Conditions: what changed between two conditions, with every sample counted once
(several samples with conditions; analyst + supervisor, 2026-11-16).

Two questions, two tests:

1. **Did the share of a group of cells change?** (composition). Each sample's share of each group,
   logit-transformed with a pseudo-count of 0.5 (as in propeller, Phipson et al. 2022, Bioinformatics),
   then a t-test per group between the conditions: paired when the same subject (e.g. donor) has a
   sample in each condition. Unlike propeller, an ordinary t-test is used (not limma's moderated
   one), so it is slightly more conservative with few samples. p-values are BH-adjusted across groups.
2. **Did a group's genes change?** (pseudobulk differential expression). For each group, the raw
   counts of its cells are summed per sample (one "pseudobulk" profile per sample), then tested with
   PyDESeq2 (Muzellec et al. 2023, the Python version of DESeq2, Love et al. 2014), design
   `~ pair + condition` when paired, else `~ condition`. Comparing samples instead of cells is the
   current best practice (Squair et al. 2021; Heumos et al. 2023): cells from one sample aren't
   independent, so cell-level tests give misleadingly small p-values.

Both need replicates: at least 2 samples per condition (2 complete pairs when paired). A group's
pseudobulk test uses only the samples with at least `min_cells` of its cells; a group without enough
samples is skipped with the reason. Nothing here knows about any particular dataset.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from scrna_analysis.params import ConditionParams

MIN_SAMPLES_PER_CONDITION = 2  # sample-to-sample variation can't be estimated from fewer
PSEUDO_COUNT = 0.5  # propeller's pseudo-count before the logit transform
MAX_GROUP_VALUES = 200  # a pairing or group column with more values than this isn't offered

# Statuses of a group's pseudobulk test.
TESTED = "tested"
TOO_FEW_SAMPLES = "too_few_samples"


class ConditionError(ValueError):
    """The comparison can't run with this data or these settings. The message says why."""


@dataclass(frozen=True)
class Design:
    """The samples compared: one row per sample with its condition ("reference"/"test") and pair."""

    samples: pd.DataFrame  # index sample; columns condition (original name), side ("reference"/"test"), pair
    reference: str
    test: str
    paired_by: str | None
    excluded: tuple[str, ...] = ()  # samples of other conditions (or incomplete pairs), left out

    @property
    def paired(self) -> bool:
        return self.paired_by is not None


@dataclass(frozen=True)
class GroupTest:
    """One group's pseudobulk test."""

    group: str
    status: str  # TESTED or TOO_FEW_SAMPLES
    samples_used: tuple[str, ...]
    samples_left_out: tuple[str, ...]  # too few of this group's cells
    n_genes_tested: int = 0
    results: pd.DataFrame | None = None  # gene, base_mean, log2_fold_change, lfc_se, stat, p, adj_p
    reason: str | None = None


@dataclass(frozen=True)
class ComparisonSummary:
    design: Design
    composition: pd.DataFrame
    expression: dict[str, GroupTest] = field(default_factory=dict)


# --- The design --------------------------------------------------------------------------------------

def sample_conditions(obs: pd.DataFrame, sample_key: str, condition_key: str,
                      paired_by: str | None = None) -> pd.DataFrame:
    """One row per sample: its condition (and its pair, if `paired_by`). Raises when a sample has
    several conditions or several pair values."""
    for column in (sample_key, condition_key, paired_by):
        if column is not None and column not in obs:
            raise ConditionError(f"There is no '{column}' column in the data.")
    cols = [condition_key] + ([paired_by] if paired_by else [])
    per = obs[[sample_key] + cols].astype(str).drop_duplicates()
    many = per[sample_key].value_counts()
    if (many > 1).any():
        name = many.index[many > 1][0]
        what = "condition" if per.groupby(sample_key)[condition_key].nunique().get(name, 1) > 1 else "subject"
        raise ConditionError(f"Sample {name} has cells with more than one {what}, so it can't be compared.")
    out = per.set_index(sample_key).rename(columns={condition_key: "condition"})
    if paired_by:
        out = out.rename(columns={paired_by: "pair"})
    return out.sort_index()


def pairing_choices(obs: pd.DataFrame, sample_key: str, condition_key: str) -> list[str]:
    """Columns that can pair samples: one value per sample, and at least two values that each have
    samples under two or more conditions (e.g. a donor column)."""
    out = []
    for column in obs.columns:
        if column in (sample_key, condition_key):
            continue
        values = obs[column]
        if not (isinstance(values.dtype, pd.CategoricalDtype) or values.dtype == object) or values.isna().any():
            continue
        if values.astype(str).nunique() > MAX_GROUP_VALUES:
            continue
        per = obs[[sample_key, condition_key, column]].astype(str).drop_duplicates()
        if per[sample_key].duplicated().any():
            continue  # not one value per sample
        if per[column].nunique() == per[sample_key].nunique():
            continue  # one value per sample: pairs nothing
        conditions_per_value = per.groupby(column)[condition_key].nunique()
        if (conditions_per_value >= 2).sum() >= MIN_SAMPLES_PER_CONDITION:
            out.append(str(column))
    return out


def make_design(samples: pd.DataFrame, reference: str | None, test: str | None,
                paired_by: str | None = None) -> Design:
    """Choose the samples to compare. Samples of other conditions are left out; when paired, only
    subjects with a sample in both conditions are kept. Raises with the reason when there aren't
    enough samples."""
    conditions = sorted(samples["condition"].unique())
    if reference is None or test is None:
        raise ConditionError("Choose the two conditions to compare.")
    for name in (reference, test):
        if name not in conditions:
            raise ConditionError(f"There is no condition '{name}'. The conditions are: {', '.join(conditions)}.")
    if reference == test:
        raise ConditionError("Choose two different conditions to compare.")
    chosen = samples[samples["condition"].isin([reference, test])].copy()
    excluded = [s for s in samples.index if s not in chosen.index]
    if paired_by:
        if "pair" not in chosen:
            raise ConditionError("The samples have no subject column to pair them.")
        both = chosen.groupby("pair")["condition"].nunique()
        complete = set(both.index[both == 2])
        per_side = chosen[chosen["pair"].isin(complete)].groupby(["pair", "condition"]).size()
        if (per_side > 1).any():
            subject = per_side.index[per_side > 1][0][0]
            raise ConditionError(f"Subject {subject} has more than one sample in the same condition, so samples "
                                 "can't be paired one to one. Compare without pairing.")
        excluded += [s for s in chosen.index if chosen.loc[s, "pair"] not in complete]
        chosen = chosen[chosen["pair"].isin(complete)]
        if len(complete) < MIN_SAMPLES_PER_CONDITION:
            raise ConditionError(f"Only {len(complete)} subject{'s' if len(complete) != 1 else ''} "
                                 f"{'have' if len(complete) != 1 else 'has'} a sample in both {reference} and {test}; "
                                 f"pairing needs at least {MIN_SAMPLES_PER_CONDITION}. Compare without pairing.")
    counts = chosen["condition"].value_counts()
    for name in (reference, test):
        if counts.get(name, 0) < MIN_SAMPLES_PER_CONDITION:
            raise ConditionError(f"{name} has {counts.get(name, 0)} sample{'s' if counts.get(name, 0) != 1 else ''}; "
                                 f"comparing conditions needs at least {MIN_SAMPLES_PER_CONDITION} samples in each "
                                 "(replicates), so sample-to-sample variation can be measured.")
    chosen["side"] = np.where(chosen["condition"] == reference, "reference", "test")
    return Design(samples=chosen, reference=reference, test=test, paired_by=paired_by, excluded=tuple(sorted(excluded)))


def natural_key(name: str) -> tuple:
    """Sort key putting "2" before "10" (cluster numbers), names that aren't numbers after them."""
    return (0, int(name), "") if str(name).isdigit() else (1, 0, str(name))


# --- 1. Composition -----------------------------------------------------------------------------------

def cell_counts(obs: pd.DataFrame, sample_key: str, group_key: str, design: Design) -> pd.DataFrame:
    """Cells of each group (columns) in each compared sample (rows), zeros included."""
    kept = obs[obs[sample_key].astype(str).isin(design.samples.index)]
    groups = (list(obs[group_key].cat.categories) if isinstance(obs[group_key].dtype, pd.CategoricalDtype)
              else sorted(obs[group_key].astype(str).unique(), key=natural_key))
    table = pd.crosstab(kept[sample_key].astype(str), kept[group_key].astype(str))
    return table.reindex(index=design.samples.index, columns=[str(g) for g in groups], fill_value=0)


def logit_shares(counts: pd.DataFrame) -> pd.DataFrame:
    """Each sample's share of each group with propeller's pseudo-count, logit-transformed."""
    pseudo = counts + PSEUDO_COUNT
    shares = pseudo.div(pseudo.sum(axis=1), axis=0)
    return np.log(shares / (1 - shares))


def composition_test(counts: pd.DataFrame, design: Design) -> pd.DataFrame:
    """Per group: mean share in each condition, the change, t statistic, p-value and BH-adjusted
    p-value (across groups). Groups with no cells in the compared samples are left out."""
    counts = counts.loc[:, counts.sum(axis=0) > 0]
    shares = counts.div(counts.sum(axis=1), axis=0)
    transformed = logit_shares(counts)
    side = design.samples["side"]
    ref, test = side.index[side == "reference"], side.index[side == "test"]
    if design.paired:
        pairs = design.samples["pair"]
        ref = sorted(ref, key=lambda s: pairs[s])
        test = sorted(test, key=lambda s: pairs[s])
    rows = []
    for group in counts.columns:
        a, b = transformed.loc[ref, group].to_numpy(), transformed.loc[test, group].to_numpy()
        t = stats.ttest_rel(b, a) if design.paired else stats.ttest_ind(b, a, equal_var=True)
        mean_ref, mean_test = float(shares.loc[ref, group].mean()), float(shares.loc[test, group].mean())
        rows.append({"group": group, "mean_share_reference": mean_ref, "mean_share_test": mean_test,
                     "log2_ratio": float(np.log2((mean_test + 1e-12) / (mean_ref + 1e-12))),
                     "t": float(t.statistic), "p": float(t.pvalue),
                     "cells_reference": int(counts.loc[ref, group].sum()), "cells_test": int(counts.loc[test, group].sum())})
    out = pd.DataFrame(rows)
    p = out["p"].to_numpy()
    ok = ~np.isnan(p)
    adj = np.full(len(p), np.nan)
    if ok.any():
        adj[ok] = stats.false_discovery_control(p[ok], method="bh")
    out["adj_p"] = adj
    return out


# --- 2. Pseudobulk expression ---------------------------------------------------------------------------

def pseudobulk(counts: ad.AnnData, groups: pd.Series, samples: pd.Series, group: str,
               keep_samples: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    """Raw counts of one group's cells summed per sample (samples × genes), and its cells per sample.
    `groups` and `samples` give each cell's group and sample, in `counts`' cell order."""
    if list(groups.index.astype(str)) != list(counts.obs_names) or list(samples.index.astype(str)) != list(counts.obs_names):
        raise ConditionError("The counts and the groups don't line up; re-run the steps from Filtering.")
    in_group = (groups.astype(str) == group).to_numpy()
    sample_of = samples.astype(str).to_numpy()
    rows, n_cells = [], {}
    for s in keep_samples:
        mask = in_group & (sample_of == s)
        n_cells[s] = int(mask.sum())
        x = counts.X[mask]
        rows.append(np.asarray(x.sum(axis=0)).ravel() if mask.any() else np.zeros(counts.n_vars))
    table = pd.DataFrame(np.vstack(rows) if rows else np.zeros((0, counts.n_vars)), index=keep_samples,
                         columns=counts.var_names.astype(str))
    return table.round().astype(np.int64), pd.Series(n_cells, dtype=int)


def _codes(design: Design, samples: list[str]) -> pd.DataFrame:
    """The design table for PyDESeq2 with safe level names (condition A = reference, B = test; pairs
    p0, p1, …), so any condition or subject name works in its formula."""
    d = design.samples.loc[samples]
    meta = pd.DataFrame({"condition": np.where(d["side"] == "reference", "A", "B")}, index=samples)
    if design.paired:
        pairs = sorted(d["pair"].unique())
        meta["pair"] = d["pair"].map({p: f"p{i}" for i, p in enumerate(pairs)}).to_numpy()
    return meta


def expression_test(counts: ad.AnnData, groups: pd.Series, samples: pd.Series, group: str, design: Design,
                    params: ConditionParams) -> GroupTest:
    """One group's pseudobulk test (PyDESeq2). Samples with fewer than `params.min_cells` of the
    group's cells are left out (and, when paired, their partner too); too few samples -> skipped."""
    names = list(design.samples.index)
    _, n_cells = pseudobulk(counts, groups, samples, group, names)
    keep = [s for s in names if n_cells[s] >= params.min_cells]
    if design.paired:
        pair = design.samples["pair"]
        sides = design.samples.loc[keep].groupby("pair")["side"].nunique()
        whole = set(sides.index[sides == 2])
        keep = [s for s in keep if pair[s] in whole]
    left_out = tuple(s for s in names if s not in keep)
    side = design.samples.loc[keep, "side"].value_counts()
    enough = (min(side.get("reference", 0), side.get("test", 0)) >= MIN_SAMPLES_PER_CONDITION)
    if not enough:
        unit = "complete pairs" if design.paired else "samples in each condition"
        have = (len(keep) // 2) if design.paired else min(side.get("reference", 0), side.get("test", 0))
        return GroupTest(group=group, status=TOO_FEW_SAMPLES, samples_used=tuple(keep), samples_left_out=left_out,
                         reason=(f"Only {have} {unit} with at least {params.min_cells} cells of this group; the test "
                                 f"needs {MIN_SAMPLES_PER_CONDITION}."))
    table, _ = pseudobulk(counts, groups, samples, group, keep)
    table = table.loc[:, table.sum(axis=0) >= params.min_gene_counts]
    if table.shape[1] == 0:
        return GroupTest(group=group, status=TOO_FEW_SAMPLES, samples_used=tuple(keep), samples_left_out=left_out,
                         reason=f"No gene has at least {params.min_gene_counts} counts in this group's samples.")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    logging.getLogger("pydeseq2").setLevel(logging.WARNING)
    meta = _codes(design, keep)
    formula = "~pair + condition" if design.paired else "~condition"
    dds = DeseqDataSet(counts=table, metadata=meta, design=formula, quiet=True, n_cpus=1)
    dds.deseq2()
    res = DeseqStats(dds, contrast=["condition", "B", "A"], alpha=params.max_adj_p, quiet=True, n_cpus=1)
    res.summary()
    df = res.results_df.rename(columns={"baseMean": "base_mean", "log2FoldChange": "log2_fold_change",
                                        "lfcSE": "lfc_se", "pvalue": "p", "padj": "adj_p"})
    df = df.rename_axis("gene").reset_index()[["gene", "base_mean", "log2_fold_change", "lfc_se", "stat", "p", "adj_p"]]
    df = df.sort_values(["p", "gene"], na_position="last", kind="mergesort").reset_index(drop=True)
    return GroupTest(group=group, status=TESTED, samples_used=tuple(keep), samples_left_out=left_out,
                     n_genes_tested=int(table.shape[1]), results=df)


# --- Both --------------------------------------------------------------------------------------------------

def compare_conditions(counts: ad.AnnData, obs: pd.DataFrame, sample_key: str, condition_key: str,
                       params: ConditionParams, groups_to_test: list[str] | None = None) -> ComparisonSummary:
    """Both tests. `counts`: raw counts of the analyzed cells (all genes kept after QC); `obs`: the
    same cells with the sample, condition, group and (if paired) subject columns. `groups_to_test`:
    the groups for the expression test (None = all)."""
    if params.group_by not in obs:
        raise ConditionError(f"There is no '{params.group_by}' column to group the cells by.")
    if params.min_cells < 1 or params.min_gene_counts < 0 or not 0 < params.max_adj_p < 1:
        raise ConditionError("Check the settings: min cells at least 1, min gene counts at least 0, and the adjusted "
                             "p-value between 0 and 1.")
    obs = obs.loc[counts.obs_names]
    samples = sample_conditions(obs, sample_key, condition_key, params.paired_by)
    design = make_design(samples, params.reference, params.test, params.paired_by)
    composition = composition_test(cell_counts(obs, sample_key, params.group_by, design), design)
    all_groups = list(composition["group"])
    chosen = all_groups if groups_to_test is None else [g for g in all_groups if g in set(groups_to_test)]
    expression = {g: expression_test(counts, obs[params.group_by], obs[sample_key], g, design, params) for g in chosen}
    return ComparisonSummary(design=design, composition=composition, expression=expression)


# --- For the record ---------------------------------------------------------------------------------------

def significant_counts(result: pd.DataFrame | None, max_adj_p: float) -> tuple[int, int]:
    """Genes significantly up and down (adjusted p-value below `max_adj_p`)."""
    if result is None:
        return 0, 0
    sig = result[result["adj_p"] < max_adj_p]
    return int((sig["log2_fold_change"] > 0).sum()), int((sig["log2_fold_change"] < 0).sum())


def short_summary(summary: ComparisonSummary, max_adj_p: float) -> dict:
    """A few numbers for the run log."""
    d = summary.design
    tested = [t for t in summary.expression.values() if t.status == TESTED]
    return {"reference": d.reference, "test": d.test, "paired_by": d.paired_by,
            "n_samples": {d.reference: int((d.samples["side"] == "reference").sum()),
                          d.test: int((d.samples["side"] == "test").sum())},
            "n_groups": int(len(summary.composition)),
            "n_groups_changed": int((summary.composition["adj_p"] < max_adj_p).sum()),
            "n_groups_tested": len(tested), "n_groups_skipped": len(summary.expression) - len(tested),
            "n_genes_changed": int(sum(sum(significant_counts(t.results, max_adj_p)) for t in tested))}


def record_summary(summary: ComparisonSummary, max_adj_p: float) -> dict:
    """The comparison for the record: the design, the composition table and, per group, what was
    tested and how many genes changed (the full gene tables are downloaded from the page)."""
    from scrna_analysis.pipeline import records
    d = summary.design
    samples = d.samples.reset_index().rename(columns={d.samples.index.name or "index": "sample"})
    return {
        "design": {"reference": d.reference, "test": d.test, "paired_by": d.paired_by,
                   "samples": records(samples), "excluded": list(d.excluded)},
        "composition": records(summary.composition),
        "expression": {g: {"status": t.status, "reason": t.reason, "samples_used": list(t.samples_used),
                           "samples_left_out": list(t.samples_left_out), "n_genes_tested": t.n_genes_tested,
                           "n_up": significant_counts(t.results, max_adj_p)[0],
                           "n_down": significant_counts(t.results, max_adj_p)[1]}
                       for g, t in summary.expression.items()},
    }
