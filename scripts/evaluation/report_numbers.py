"""Recompute the report's derived numbers and check them against ``report/results_numbers.tex``.

**Why this exists.** ``results_numbers.tex`` is transcribed by hand (TODO item 12), and on 12.10.2026 a
survey found three classes of macro with no committed artifact behind the value itself (TODO, *Numbers in
the report without a committed artifact*):

* **B** -- printed by a committed script but never written: the global label-quality figures of
  ``label_quality.py``.
* **C** -- computed ad hoc from the data files: the CTRPv2 line count, IC50 availability, the target
  density and the cells per line.
* **D** -- derived from committed CSVs by a reduction stated only in a macro comment.

This script performs each reduction in code, writes every value to
``notebooks/outputs/report/report_numbers.csv`` with its source and reduction, and compares it with the
macro, rounded to the macro's own decimals. It **does not write** ``results_numbers.tex``: the macros stay
hand-maintained, and this is the check that they still match their sources.

Class A (no artifact at all: the input-scale test, the CPU-vs-MPS embedding check, the per-compound
fragility values) is out of scope -- those need a run, not a reduction.

Exit status 1 on any mismatch, and also when nothing was compared, because "0 mismatches" and
"0 compared" are otherwise indistinguishable (the rule of ``scripts/gate/``).

Run:  .venv/bin/python scripts/evaluation/report_numbers.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.layout import PipelinePaths  # noqa: E402

OUT = ROOT / "notebooks/outputs"
TEX = ROOT / "report/results_numbers.tex"
CSV = OUT / "report/report_numbers.csv"

#: The highest-scoring configuration -- what \SpBest, the per-compound range and the flat-curve
#: correlation are read on.
BEST = dict(rep="X_pca", alpha=0.5, loss="mae")


def label_quality_globals() -> list[dict]:
    """B: the fit statistics of the curves behind the panel, as ``label_quality.py`` section 1 prints
    them -- same filter (panel compounds, split lines, a defined AUC), written instead of printed."""
    from scripts.preprocessing.ctrp_to_h5ad import CTRP_CELL_LINE_ALIASES, _normalize_cell_line, _normalize_drug

    p = PipelinePaths.build(None, "hvg5000", "auc_cc")
    # Exact duplicate rows dropped first, over the full row, as the pipeline drops them
    # (ctrp_to_h5ad._deduplicate_measurements; Selin, 11.08.2026). Added 22.10.2026: without it 64 curves
    # behind the panel were counted twice (1,699 rows for 1,635 pairs).
    resp = pd.read_csv(p.ctrp_response_csv).drop_duplicates()
    # Same filter as label_quality.py, kept in step with it (21.10.2026, Selin): the pipeline's alias
    # table on the CTRP side, and the 153 CV lines only -- not the 28 test lines.
    resp["line"] = _normalize_cell_line(resp.ccl_name).replace(CTRP_CELL_LINE_ALIASES)
    resp["drug"] = _normalize_drug(resp.drug_name)
    panel = pd.read_csv(OUT / "panel/panel.csv")
    split = pd.read_csv(ROOT / "splits/split_ctrp.csv")
    split = split[split.split.isin(["train", "val"])]
    drugs = set(_normalize_drug(panel.drug_ctrp))
    lines = set(_normalize_cell_line(split.Cell_line.str.split("_").str[0]))
    r = resp[resp.drug.isin(drugs) & resp.line.isin(lines)].dropna(subset=["AUC_curvecurator"])
    assert not r.duplicated(["line", "drug"]).any(), "a (line, compound) pair still carries two curves"
    src, how = "DrEval CTRPv2 response table", "curves of panel compounds x CV lines (label_quality.py filter)"
    return [
        row("NCurves", len(r), src, how + "; count"),
        row("CurveMedRsq", r.R2.median(), src, how + "; median R2"),
        row("CurveFracLowRsq", 100 * (r.R2 < 0.5).mean(), src, how + "; % R2 < 0.5"),
        row("CurveFracNs", 100 * (r.pValue > 0.05).mean(), src, how + "; % pValue > 0.05 (flat)"),
    ]


def data_counts() -> list[dict]:
    """C: counts read from the data files rather than from a model output."""
    import anndata as ad

    p = PipelinePaths.build(None, "hvg5000", "auc_cc")
    resp = pd.read_csv(p.ctrp_response_csv, usecols=["ccl_name", "IC50_curvecurator"])
    a = ad.read_h5ad(p.targets_h5ad, backed="r")
    m = np.asarray(a.obsm["M_ctrp"])
    labelled = m.sum(axis=1) > 0
    line = a.obs["Cell_line"].to_numpy()
    line_mask = pd.DataFrame(m[labelled]).groupby(line[labelled]).max()
    first = pd.Series(np.arange(len(line)))[labelled].groupby(line[labelled]).first().to_numpy()
    y_line, m_line = np.asarray(a.obsm["Y_ctrp"])[first], m[first].astype(bool)
    per_line = pd.Series(line[labelled]).value_counts()
    # The panel's share of screened pairs over the CV lines -- the columns and lines the objective sees
    # (cv.py: `kcol` from panel.csv, `grouped_folds` eligible splits train + val).
    panel = list(pd.read_csv(OUT / "panel/panel.csv").drug_key)
    kcol = [list(a.uns["ctrp_drugs"]).index(d) for d in panel]
    cv = a.obs["split_ctrp"].isin(("train", "val")).to_numpy()
    cv_mask = pd.DataFrame(m[cv][:, kcol]).groupby(line[cv]).max()
    src_r, src_t = "DrEval CTRPv2 response table", "targets h5ad (hvg5000, auc_cc)"
    return [
        row("NCtrpLines", resp.ccl_name.nunique(), src_r, "distinct ccl_name"),
        row("IcDefinedPct", 100 * resp.IC50_curvecurator.notna().mean(), src_r, "% rows with an IC50"),
        row("TargetDensity", 100 * line_mask.to_numpy().mean(), src_t,
            "% of (labelled line, compound) pairs with M_ctrp = 1"),
        row("PanelDensity", 100 * cv_mask.to_numpy().mean(), src_t,
            "% of (CV line, panel compound) pairs with M_ctrp = 1; CV = split_ctrp train + val"),
        row("NCancerTypes", a.obs["Cancer_type"].nunique(), src_t, "distinct Cancer_type"),
        row("AucMedian", np.median(y_line[m_line]), src_t,
            "median auc_cc over screened (labelled line, compound) pairs, all compounds"),
        row("CellsMin", per_line.min(), src_t, "fewest cells in a labelled line"),
        row("CellsMax", per_line.max(), src_t, "most cells in a labelled line"),
        row("CellsPerLineMean", per_line.mean(), src_t, "mean cells per labelled line"),
        # Added 21.10.2026: NCells (53,513) is the whole atlas, 198 lines, but the report set it beside
        # the 181 labelled lines and the 153 CV lines.
        row("NCellsLabelled", int(labelled.sum()), src_t, "cells of the labelled lines (any M_ctrp = 1)"),
        row("NCellsCV", int(cv.sum()), src_t, "cells of the CV lines (split_ctrp train + val)"),
    ]


def sweep() -> list[dict]:
    """D: the twelve configurations of Table 1, the per-compound range and the flat-curve correlation."""
    src = "panel/panel_per_drug_correlation.csv"
    corr = pd.read_csv(OUT / src)
    rows = []
    names = {("mse", 0.0): "MseZero", ("mse", 0.5): "MseHalf", ("mse", 1.0): "MseOne",
             ("mae", 0.0): "MaeZero", ("mae", 0.5): "MaeHalf", ("mae", 1.0): "MaeOne"}
    for rep, tag in (("X_pca", "Pca"), ("X_scGPT", "Scgpt")):
        for (loss, alpha), name in names.items():
            g = corr[(corr.rep == rep) & (corr.loss == loss) & (corr.alpha == alpha)]
            rows.append(row(f"Sw{tag}{name}", g.groupby("seed").spearman.mean().mean(), src,
                            f"{rep} {loss} alpha={alpha}: mean over compounds per seed, then over seeds"))
    best = corr[(corr.rep == BEST["rep"]) & (corr.alpha == BEST["alpha"]) & (corr.loss == BEST["loss"])]
    per_drug = best.groupby("drug").spearman.mean()
    rows += [row("SpDrugLo", per_drug.min(), src, f"best configuration, lowest compound ({per_drug.idxmin()})"),
             row("SpDrugHi", per_drug.max(), src, f"best configuration, highest compound ({per_drug.idxmax()})")]

    lq_src = "diagnostics/label_quality_vs_performance.csv"
    lq = pd.read_csv(OUT / lq_src).set_index("drug")
    j = lq.join(per_drug.rename("rho")).dropna(subset=["frac_ns", "rho"])
    rho, pval = spearmanr(j.frac_ns, j.rho)
    rows += [row("FlatRho", rho, f"{lq_src} + {src}", "Spearman(frac_ns, best-configuration score) over compounds"),
             row("FlatP", pval, f"{lq_src} + {src}", "its p-value"),
             row("PlatinFracNs", 100 * lq.loc["platin", "frac_ns"], lq_src, "cisplatin frac_ns, %"),
             row("PlatinMedRsq", lq.loc["platin", "median_R2"], lq_src, "cisplatin median_R2")]

    band_src = "panel/panel_execution_band.csv"
    band = pd.read_csv(OUT / band_src)
    alt = band[(band.arm == "MLP mse alpha=0") & (band.rep == "X_pca")]["exec_1_committed_13.08"].iloc[0]
    rows.append(row("SwPcaMseZeroAlt", alt, band_src, "X_pca mse alpha=0, first execution"))
    # Added 21.10.2026 (Selin): the band holds three executions since the audit rerun, and the report
    # names all three values. The third is SwPcaMseZero itself, read from panel_leaderboard.csv.
    alt2 = band[(band.arm == "MLP mse alpha=0") & (band.rep == "X_pca")]["exec_2_full_rerun_13.08"].iloc[0]
    rows.append(row("SwPcaMseZeroAltB", alt2, band_src, "X_pca mse alpha=0, second execution"))
    return rows


def controls() -> list[dict]:
    """D: ridge, the unselected catalogue, DrEval LCO and LPO."""
    rows = []
    src = "panel/panel_ridge_baseline.csv"
    ridge = pd.read_csv(OUT / src).groupby("rep").spearman.mean()
    rows += [row("RidgePca", ridge["X_pca"], src, "mean spearman over compounds"),
             row("RidgeScgpt", ridge["X_scGPT"], src, "mean spearman over compounds")]
    # Where the penalty grid still binds (21.10.2026): per fit, from ridge_grid_refit.py --edges, which
    # asserts that its refit reproduces panel_ridge_baseline.csv before writing.
    src_e = "panel/ridge_alpha_per_fit.csv"
    pf = pd.read_csv(OUT / src_e)
    n = pf.groupby("rep").size()
    assert n.nunique() == 1, "fit counts differ between representations"
    edges = pf.groupby("rep")[["at_floor", "at_ceiling"]].sum()
    rows += [row("RidgeFitsPerRep", n.iloc[0], src_e, "ridge fits per representation (compounds x folds)"),
             row("RidgeFloorPca", edges.loc["X_pca", "at_floor"], src_e, "fits at the grid floor, X_pca"),
             row("RidgeCeilPca", edges.loc["X_pca", "at_ceiling"], src_e, "fits at the grid ceiling, X_pca"),
             row("RidgeFloorScgpt", edges.loc["X_scGPT", "at_floor"], src_e, "fits at the grid floor, X_scGPT"),
             row("RidgeCeilScgpt", edges.loc["X_scGPT", "at_ceiling"], src_e, "fits at the grid ceiling, X_scGPT")]

    src = "panel/panel_heads_summary.csv"
    heads = pd.read_csv(OUT / src)
    allc = heads[heads.scored_on.str.startswith("all ")]
    rows += [row("SpAllDrugs", allc[(allc.arch == "linear") & (allc.rep == "X_pca")].spearman.iloc[0], src,
                 "all compounds, linear head, X_pca")]
    # Both architectures at both scorings, so the (capacity x head count) factorial of 4a sections C
    # and D is reported whole (Selin, 21.10.2026, option ii). Until then only the linear PCA row was.
    # NullGainMax (largest vs_null of the four rows, in MSE units, 0.0004) was retired the same day:
    # it sat beside a Spearman and came from a different row than the score it qualified.
    pan = heads[~heads.scored_on.str.startswith("all ")]
    for arch, a_name in (("linear", "Lin"), ("trunk", "Trunk")):
        for rep, r_name in (("X_pca", "Pca"), ("X_scGPT", "Scgpt")):
            sel = lambda df: df[(df.arch == arch) & (df.rep == rep)].spearman.iloc[0]
            if (arch, rep) != ("linear", "X_pca"):   # that one is SpAllDrugs
                rows.append(row(f"HdAll{a_name}{r_name}", sel(allc), src, f"all compounds, {arch}, {rep}"))
            rows.append(row(f"HdPan{a_name}{r_name}", sel(pan), src, f"534-head model, panel compounds, {arch}, {rep}"))
    rows += heads_vs_constant()

    src = "dreval/dreval_lco_results.csv"
    lco = pd.read_csv(OUT / src).groupby("algorithm")
    norm, raw = lco["Spearman: normalized"], lco["Spearman"]
    rows += [row("DrMlpNorm", norm.mean()["OncoMLP (X_pca)"], src, "normalized Spearman, mean over folds"),
             row("DrRfNorm", norm.mean()["SingleDrugRF (scgpt)"], src, "normalized Spearman, mean over folds"),
             row("DrMlpSd", norm.std()["OncoMLP (X_pca)"], src, "normalized Spearman, SD over folds"),
             row("DrRfSd", norm.std()["SingleDrugRF (scgpt)"], src, "normalized Spearman, SD over folds"),
             # The forest on the MLP's own representation, so the comparison holds the features fixed
             # (Selin, 21.10.2026, option b). DrRfNorm stays the scGPT forest.
             row("DrRfPcaNorm", norm.mean()["SingleDrugRF (pca)"], src, "normalized Spearman, mean over folds"),
             row("DrRfPcaSd", norm.std()["SingleDrugRF (pca)"], src, "normalized Spearman, SD over folds")]
    # The same forest on the first 10 PCA components only (22.10.2026), which the Discussion compares with.
    src_rf = "dreval/dreval_rf_checks.csv"
    rf = pd.read_csv(OUT / src_rf)
    rows.append(row("DrRfPcaTen", rf[(rf.view == "pca") & (rf.k == 10) & (rf.max_samples == 0.2)]
                    ["Spearman: normalized"].mean(), src_rf, "normalized Spearman, k=10, default forest, mean over folds"))

    src = "dreval/dreval_lpo_results.csv"
    lpo = pd.read_csv(OUT / src).groupby("algorithm").per_drug_Spearman.mean()
    rows += [row("LpoLineMean", lpo["NaiveCellLineMeanPredictor"], src, "per_drug_Spearman, mean over folds"),
             row("LpoTissueMean", lpo["NaiveTissueMeanPredictor"], src, "per_drug_Spearman, mean over folds")]
    return rows


def heads_vs_constant() -> list[dict]:
    """The 534-head model against a per-compound constant, over all compounds, per (arch, rep).

    Two quantities (Selin, 21.10.2026, options b and a-plus-share):

    * the **relative error reduction**, (constant MSE - model MSE) / constant MSE, each first averaged
      over compounds -- a ratio of means, so compounds with more response variance weigh more, as they
      already do in ``vs_null``;
    * the **share of compounds** on which the model's MSE is below the constant's.

    Both per seed, then averaged over the three seeds. The constant is 4a section D2's ``_null_mse``,
    restated here: the compound's mean label on the other folds' lines, never on the held-out truth.
    Compounds with fewer than 5 lines or no label variance are skipped, as section D2 skips them.
    """
    src = "panel/panel_heads_oof.csv"
    d = pd.read_csv(OUT / src)

    def null_mse(g):
        err = [(held.y_true.to_numpy() - g.loc[g.fold != f, "y_true"].mean()) ** 2
               for f, held in g.groupby("fold", sort=False)]
        return float(np.concatenate(err).mean())

    recs = []
    for (arch, rep, seed, _drug), g in d.groupby(["arch", "rep", "seed", "drug"]):
        if len(g) < 5 or g.y_true.std(ddof=0) == 0:
            continue
        recs.append(dict(arch=arch, rep=rep, seed=seed, null=null_mse(g),
                         model=float(((g.y_pred - g.y_true) ** 2).mean())))
    x = pd.DataFrame(recs)
    per_seed = x.groupby(["arch", "rep", "seed"]).agg(null=("null", "mean"), model=("model", "mean"),
                                                      beat=("model", "size"))
    per_seed["beat"] = (x.assign(b=x.model < x.null).groupby(["arch", "rep", "seed"]).b.mean())
    per_seed["gain"] = 100 * (per_seed.null - per_seed.model) / per_seed.null
    m = per_seed.groupby(["arch", "rep"])[["gain", "beat"]].mean()
    rows = []
    for arch, a_name in (("linear", "Lin"), ("trunk", "Trunk")):
        for rep, r_name in (("X_pca", "Pca"), ("X_scGPT", "Scgpt")):
            rows += [row(f"HdGain{a_name}{r_name}", m.loc[(arch, rep), "gain"], src,
                         "% MSE reduction vs per-compound constant, ratio of means, mean over seeds"),
                     row(f"HdBeat{a_name}{r_name}", 100 * m.loc[(arch, rep), "beat"], src,
                         "% of compounds with model MSE below the constant's, mean over seeds")]
    rows.append(row("HdRiseLinPca", -m.loc[("linear", "X_pca"), "gain"], src,
                    "% by which the linear head's MSE exceeds the constant's (X_pca); stated as a rise in the text"))
    return rows


def fragility() -> list[dict]:
    """D: the fragility accounting, mean over the six configurations per representation."""
    src = "dreval/bias_accounting.csv"
    b = pd.read_csv(OUT / src)
    m = b.groupby("rep")[["rho_pred_fragility", "spearman_raw_per_drug", "partial_rho"]].mean()
    best = b[(b.alpha == BEST["alpha"]) & (b.loss == BEST["loss"])].set_index("rep")
    rows = []
    for rep, tag in (("X_pca", "Pca"), ("X_scGPT", "Scgpt")):
        raw, part = m.loc[rep, "spearman_raw_per_drug"], m.loc[rep, "partial_rho"]
        rows += [row(f"FragRho{tag}", m.loc[rep, "rho_pred_fragility"], src, "mean over six configurations"),
                 row(f"FragRaw{tag}", raw, src, "mean spearman_raw_per_drug over six configurations"),
                 row(f"FragPartial{tag}", part, src, "mean partial_rho over six configurations"),
                 row(f"FragCost{tag}", 100 * (1 - part / raw), src, "1 - partial/raw, %"),
                 row(f"FragBest{tag}", best.loc[rep, "partial_rho"], src, "partial_rho, best configuration")]
    pooled = b[(b.rep == BEST["rep"]) & (b.alpha == BEST["alpha"]) & (b.loss == BEST["loss"])]
    rows.append(row("SpPooled", pooled.spearman_raw_pooled.iloc[0], src, "pooled raw Spearman, best configuration"))
    # Same folds, pairs and pooling as SpPooled (Selin, 21.10.2026, option a); replaces DrDrugMean, which
    # came from DrEval's own LCO splits and was pooled per fold.
    rows.append(row("SpDrugMeanPooled", pooled.spearman_drug_mean_pooled.iloc[0], src,
                    "pooled raw Spearman of DrEval's NaiveDrugMeanPredictor on the same out-of-fold pairs"))
    return rows


def representations() -> list[dict]:
    """D: the input ratio, the four alternative explanations and the input-dropout sweep."""
    rows = []
    src = "mil/stage0_input_ceiling.csv"
    w = pd.read_csv(OUT / src).set_index("rep").within_line_share
    for rep, tag in (("X_pca", "Pca"), ("X_scGPT", "Scgpt")):
        rows += [row(f"BtwShare{tag}", 1 - w[rep], src, "1 - within_line_share"),
                 row(f"BtwRatio{tag}", (1 - w[rep]) / w[rep], src, "(1 - w) / w")]
    # The paired line bootstrap of that ratio (4b §3.2b, added 19.10.2026).
    src = "mil/stage0_bootstrap.csv"
    bs = pd.read_csv(OUT / src).set_index("quantity")
    d = bs.loc["ratio_diff_pca_minus_scgpt"]
    rows += [row("BtwRatioDiff", d.point, src, "PCA ratio minus scGPT ratio"),
             row("BtwRatioDiffLo", d.ci_lo, src, "2.5th percentile of the paired difference"),
             row("BtwRatioDiffHi", d.ci_hi, src, "97.5th percentile of the paired difference"),
             row("BtwBootPcaAheadPct", 100 * d.frac_resamples_pca_ahead, src, "% of resamples, PCA ratio higher")]

    src = "panel/panel_curve_summary.csv"
    c = pd.read_csv(OUT / src).pivot(index="budget", columns="rep", values="mean")
    full = 0  # budget 0 is the full label supply (labelled_lines = 103)
    rows += [row("CurveGainPca", c.loc[full, "X_pca"] - c.loc[25, "X_pca"], src, "full minus 25 lines"),
             row("CurveGainScgpt", c.loc[full, "X_scGPT"] - c.loc[25, "X_scGPT"], src, "full minus 25 lines"),
             row("CurveMarginLo", c.loc[25, "X_pca"] - c.loc[25, "X_scGPT"], src, "PCA minus scGPT at 25 lines"),
             row("CurveMarginHi", c.loc[full, "X_pca"] - c.loc[full, "X_scGPT"], src, "PCA minus scGPT, full")]

    src = "panel/panel_arch_summary.csv"
    a = pd.read_csv(OUT / src).set_index(["arch", "rep"])["mean"]
    rows += [row("ArchLinPca", a["linear", "X_pca"], src, "linear head"),
             row("ArchTrunkPca", a["trunk", "X_pca"], src, "(128,64) trunk"),
             row("ArchLinScgpt", a["linear", "X_scGPT"], src, "linear head"),
             row("ArchTrunkScgpt", a["trunk", "X_scGPT"], src, "(128,64) trunk")]

    src = "diagnostics/input_scale.csv"
    s = pd.read_csv(OUT / src).set_index("rep").std_median
    rows.append(row("InputScaleRatio", s["X_pca"] / s["X_scGPT"], src, "ratio of median per-dimension SD"))
    # Checked since 21.10.2026: the rescaling test's epochs, until then hand-copied (and the report said
    # the output was not committed, though it was).
    src = "diagnostics/input_scale_test.csv"
    e = pd.read_csv(OUT / src).groupby("rep").best_epoch.mean()
    rows += [row("ScaleEpochNative", e["X_scGPT"], src, "mean best epoch over seeds x folds, X_scGPT"),
             row("ScaleEpochScaled", e["X_scGPT_scaled"], src, "mean best epoch, X_scGPT rescaled to X_pca")]

    srcs = ["diagnostics/input_dropout_test.csv", "diagnostics/input_dropout_sweep_extra.csv"]
    d = pd.concat([pd.read_csv(OUT / f) for f in srcs]).groupby(["rep", "input_dropout"]).order.mean()
    src = " + ".join(srcs)
    rows += [row("DropRangePca", d["X_pca"].max() - d["X_pca"].min(), src, "range over the six rates"),
             row("DropRangeScgpt", d["X_scGPT"].max() - d["X_scGPT"].min(), src, "range over the six rates"),
             row("DropBestPca", d["X_pca", 0.02], src, "X_pca at rate 0.02"),
             row("DropShippedPca", d["X_pca", 0.1], src, "X_pca at rate 0.1")]
    return rows


def heterogeneity() -> list[dict]:
    """D: the MIL criterion, one row per representation of the verdict table."""
    src = "mil/q2_verdict.csv"
    v = pd.read_csv(OUT / src).set_index("rep")
    rows = []
    for rep, tag in (("X_pca", "Pca"), ("X_scGPT", "Scgpt")):
        repr_ = v.loc[rep, "stage2_median_rho"] ** 2
        conf = v.loc[rep, "stage6_median_r2_adj"]
        rows += [row(f"QSpread{tag}", 100 * v.loc[rep, "stage1_frac_mil_exceeds"], src, "% pairs, stage 1"),
                 row(f"QSeedRho{tag}", v.loc[rep, "stage2_median_rho"], src, "median cross-seed rho"),
                 row(f"QAuroc{tag}", v.loc[rep, "stage7_median_auroc"], src, "median within-bag AUROC"),
                 row(f"QConf{tag}", conf, src, "median adjusted R2 of the confounds"),
                 row(f"QRepr{tag}", repr_, src, "squared stage2_median_rho (median over seed pairs of each pair's median)"),
                 row(f"QBar{tag}", 100 * conf / repr_, src, "confound share of the veto bar, %"),
                 row(f"QSeedSig{tag}", 100 * v.loc[rep, "stage2_frac_significant"], src,
                     "% of (compound, line) pairs beating the within-line permutation null, stage 2"),
                 row(f"QConfRaw{tag}", v.loc[rep, "stage6_median_r2"], src, "median unadjusted R2 of the confounds"),
                 row(f"QAurocSig{tag}", 100 * v.loc[rep, "stage7_frac_significant"], src,
                     "% of positive-control pairs beating the permutation null, stage 7")]
    # Stage 1 per seed (19.10.2026): the range over the three seeds of each seed's median within-line SD,
    # MIL and per-cell model. A range, not a pooled value, so no aggregation over seeds is chosen here.
    src = "mil/stage1_spread_vs_percell.csv"
    s1 = pd.read_csv(OUT / src)
    for rep, tag in (("X_pca", "Pca"), ("X_scGPT", "Scgpt")):
        s = s1[s1.rep == rep]
        rows += [row(f"QSdMil{tag}Lo", s.median_sd_mil.min(), src, "min over seeds of median within-line SD, MIL"),
                 row(f"QSdMil{tag}Hi", s.median_sd_mil.max(), src, "max over seeds, MIL"),
                 row(f"QSdCell{tag}Lo", s.median_sd_4a.min(), src, "min over seeds, per-cell model"),
                 row(f"QSdCell{tag}Hi", s.median_sd_4a.max(), src, "max over seeds, per-cell model")]
    rows.append(row("QPairs", s1.n_pairs.iloc[0], src, "(compound, line) pairs per seed"))
    # Stage 7's own unit (21.10.2026, Selin, item 10 option a): mixed two-line bags, not (compound, line)
    # pairs, so the report no longer gives QPairs as the count of every stage.
    s7 = pd.read_csv(OUT / "mil/stage7_positive_control.csv").groupby(["rep", "seed"]).size()
    assert s7.nunique() == 1, "stage 7 bag count differs between arms or seeds"
    rows.append(row("QBags", s7.iloc[0], "mil/stage7_positive_control.csv", "mixed two-line bags per seed and arm"))
    # The (compound, line) combinations the label-free stages include although CTRPv2 never screened
    # them (21.10.2026, Selin): all combinations less the screened ones, which are the pairs scored out of
    # fold in panel_oof_predictions.csv.
    oof = pd.read_csv(OUT / "panel/panel_oof_predictions.csv", usecols=["drug", "cell_line"]).drop_duplicates()
    rows.append(row("QUnscreened", s1.n_pairs.iloc[0] - len(oof), f"{src} + panel/panel_oof_predictions.csv",
                    "QPairs less the screened (compound, line) pairs"))
    # Closing description of Q3 (4b §4, cells 29 and 31; added 19.10.2026): the scGPT arm's within-line
    # association of per-cell predictions with the proteasomal-degradation program, with and without
    # mitochondrial fraction controlled, and that program's own within-line relation to depth.
    src = "mil/closing_program_partial_mito.csv"
    pp = pd.read_csv(OUT / src)
    d = pp[(pp.rep == "X_scGPT") & (pp.program == "ProtDegra_score")]
    rows += [row("QProgDegraScgpt", d.rho.median(), src, "median within-line rho, prediction vs program"),
             row("QProgDegraScgptPartial", d.rho_partial_mito.median(), src,
                 "median within-line partial rho, mitochondrial fraction controlled")]
    # E2, the smaller-study test (scripts/evaluation/section_e2_smaller_study.py; added to the report
    # 19.10.2026): PCA minus scGPT, mean over seeds, at each number of lines kept.
    src = "panel/panel_curve_e2.csv"
    e2 = pd.read_csv(OUT / src).groupby(["n_lines_kept", "rep"]).order.mean().unstack("rep")
    m2 = e2["X_pca"] - e2["X_scGPT"]
    ks = sorted(e2.index)
    rows += [row("EtwoKeepMin", ks[0], src, "fewest lines kept"),
             row("EtwoKeepA", ks[1], src, "second budget"), row("EtwoKeepB", ks[2], src, "third budget"),
             row("EtwoPcaMin", e2.loc[ks[0], "X_pca"], src, "X_pca at the fewest lines"),
             row("EtwoScgptMin", e2.loc[ks[0], "X_scGPT"], src, "X_scGPT at the fewest lines"),
             row("EtwoMarginA", m2[ks[1]], src, "PCA minus scGPT, second budget"),
             row("EtwoMarginB", m2[ks[2]], src, "PCA minus scGPT, third budget"),
             row("EtwoMarginAll", m2[ks[3]], src, "PCA minus scGPT, all lines")]
    # The MIL model's own line-level score against the per-cell model at the same alpha and loss.
    src = "panel/panel_metrics.csv"
    pm = pd.read_csv(OUT / src)
    pm = pm[(pm.alpha == 0.5) & (pm.loss == "mse")].groupby(["model", "rep"]).order.mean()
    rows += [row("MilScorePca", pm["mil", "X_pca"], src, "MIL, alpha 0.5, mse, mean over seeds"),
             row("MilScoreScgpt", pm["mil", "X_scGPT"], src, "MIL, alpha 0.5, mse, mean over seeds"),
             row("MilMarginPercell", pm["mlp", "X_pca"] - pm["mlp", "X_scGPT"], src, "PCA minus scGPT, per-cell"),
             row("MilMarginBag", pm["mil", "X_pca"] - pm["mil", "X_scGPT"], src, "PCA minus scGPT, MIL")]
    src = "diagnostics/mil_by_within_line_spread.csv"
    sp = pd.read_csv(OUT / src).set_index(["half", "rep"]).bag_minus_percell
    rows += [row("MilGapLowPca", sp["low", "X_pca"], src, "MIL minus per-cell, low-spread half"),
             row("MilGapHighPca", sp["high", "X_pca"], src, "high-spread half"),
             row("MilGapLowScgpt", sp["low", "X_scGPT"], src, "low-spread half"),
             row("MilGapHighScgpt", sp["high", "X_scGPT"], src, "high-spread half")]
    src = "mil/closing_program_depth.csv"
    pd_ = pd.read_csv(OUT / src)
    rows.append(row("QProgDegraDepth", pd_[(pd_.program == "ProtDegra_score")
                                           & (pd_.covariate == "total_counts")].rho.median(), src,
                    "median within-line rho, program score vs total UMI count"))
    return rows


def row(macro: str, value, source: str, reduction: str) -> dict:
    return {"macro": macro, "value": float(value), "source": source, "reduction": reduction}


def read_macros(path: Path) -> dict[str, str]:
    """``\\newcommand{\\Name}{value}`` -> {Name: value}, for numeric values only (spelled words and
    macro references such as ``\\SpBest{\\SwPcaMaeHalf}`` are skipped)."""
    out = {}
    # Anything after a macro's closing brace must be whitespace, a comment or another definition. Added
    # 21.10.2026: an edit left `\NCurves{1{,}699}971}`, the pattern below read `1{,}699` and the check
    # passed on a macro that typesets as "1,699971}".
    text = path.read_text()
    for m in re.finditer(r"\\newcommand\{\\(\w+)\}\{((?:[^{}]|\{,\})*)\}([^\n]*)", text):
        rest = m.group(3).lstrip()
        if rest and not rest.startswith(("%", "\\newcommand")):
            raise SystemExit(f"\\{m.group(1)}: unexpected text after the definition: {rest[:40]!r}")
    for name, val in re.findall(r"\\newcommand\{\\(\w+)\}\{((?:[^{}]|\{,\})*)\}", text):
        val = val.replace("{,}", "")
        if re.fullmatch(r"-?\d+(\.\d+)?", val):
            out[name] = val
    return out


def compare(values: pd.DataFrame, macros: dict[str, str]) -> pd.DataFrame:
    """Round each computed value to the macro's own number of decimals and compare the strings."""
    recs = []
    for r in values.itertuples():
        if r.macro not in macros:
            recs.append((r.macro, r.value, None, "NO MACRO"))
            continue
        m = macros[r.macro]
        dec = len(m.split(".")[1]) if "." in m else 0
        got = f"{r.value:.{dec}f}"
        recs.append((r.macro, got, m, "ok" if float(got) == float(m) else "MISMATCH"))
    return pd.DataFrame(recs, columns=["macro", "computed", "macro_value", "status"])


def main() -> int:
    values = pd.DataFrame(label_quality_globals() + data_counts() + sweep() + controls()
                          + fragility() + representations() + heterogeneity())
    CSV.parent.mkdir(parents=True, exist_ok=True)
    values.to_csv(CSV, index=False)
    print(f"wrote {CSV.relative_to(ROOT)}  ({len(values)} values)")

    result = compare(values, read_macros(TEX))
    checked = result[result.status != "NO MACRO"]
    bad = result[result.status != "ok"]
    print(f"compared {len(checked)} of {len(values)} values with {TEX.relative_to(ROOT)}")
    if len(bad):
        print(bad.to_string(index=False))
    if len(checked) == 0:
        print("FAIL: nothing was compared")
        return 1
    print("FAIL" if len(bad) else "OK: every computed value matches its macro")
    return 1 if len(bad) else 0


if __name__ == "__main__":
    raise SystemExit(main())
