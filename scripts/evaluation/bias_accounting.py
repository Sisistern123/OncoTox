"""Bias accounting: DrEval's normalized metric against ours, and how fragility-driven each arm is.

**Why this file exists at all.** Both tables below were already stated in `docs/steps/05` on
14.08.2026 — and both had been computed in an *uncommitted* shell session, which `CLAUDE.md` says is
not a result. Re-deriving them here found that neither reproduced to the quoted decimals (the
metric-equivalence row read `0.2421 -> 0.2692` and reads **0.2473 -> 0.2750**; the fragility means
read `0.1660 / 0.2224` and `0.1200 / 0.1778` and read **0.1602 / 0.2138** and **0.1172 / 0.1733**).
**Every conclusion drawn from them is unchanged** — see the two sections below — but the numbers
moved in the third decimal, which is exactly the failure mode an uncommitted computation produces.

---

## Table 1 — is our metric a way of dodging DrEval's?

**The objection.** DrEval exists because most published drug-response models score well through mean
effects alone, and their answer is a *normalized* correlation: subtract
``overall mean + cell-line effect + drug effect`` (their ``NaiveMeanEffectsPredictor``) from truth and
from prediction, then correlate the residuals. This project reports a **mean per-drug Spearman**
instead. Reporting one's own metric in place of the field's bias metric is exactly what a DrEval
reader should be suspicious of.

**The answer, measured rather than argued.** Under **leave-cell-line-out** the two nearly coincide,
because their normalization has no line effect available: a held-out line was never seen, so
``cell_line_effects.get(cl, 0)`` returns 0. What their baseline actually removes here is the drug
effect — which is what correlating *within* a drug removes by construction. This is verified, not
assumed: the ``n_distinct_naive`` column below counts how many distinct fitted baseline values occur
inside a single drug, and it is **5** — one per fold, none varying by cell line.

So the per-drug figure barely moves under their normalization, while the *pooled* figure collapses
from ~0.77 to ~0.31. **The pooled number is the flattering one and this project does not quote it.**

## Table 2 — which representation leans harder on line fragility?

**Why it is a Q1 question, not a footnote.** If ``X_pca`` leads because it exploits general
fragility — some lines die to everything — that is a different claim from "PCA carries more of what
the label varies over", and it would undercut Q1 rather than support it.

**Method (fixed 14.08.2026, unchanged here).** For each drug, the fragility proxy of a cell line is
its mean ``y_true`` over the **other ten** panel drugs — leave-one-drug-out, so a drug never enters
its own proxy. Two quantities per arm: ``rho_pred_fragility``, how fragility-like the predictions
are; and ``partial_rho``, the rank-partial correlation between truth and prediction with the proxy
controlled out. Averaged over drugs, then over seeds, then reported per arm.

⚠️ **Limits, and they are not small.** Eleven drugs; the proxy is built from those same eleven, so it
approximates "general fragility" rather than measuring it; rank-partial correlation controls only
linearly in the ranks. Nothing here decomposes the score into a fragility share — that would need a
different design, and `dreval_normalize.py` records why this project does not have one.

---

**Aggregation, stated because it is a choice.** Every score is computed **per seed and then averaged
over the three seeds**, which is the convention the leaderboard uses; the ``spearman_raw_per_drug``
column therefore reproduces ``panel_leaderboard.csv`` exactly, and that is worth reading as a
coherence check on this file rather than as a separate result.

⚠️ **The committed ``dreval_normalized.csv`` uses a different convention** — ``dreval_normalize.py``
filters on ``rep`` and ``alpha`` only, so it pools two losses x three seeds (9,810 rows) and scores
that stack as one. That is not wrong, but it is not comparable to a per-arm number, and the two files
will disagree in the second decimal. Nothing here changes that file.

Run:  .venv/bin/python scripts/evaluation/bias_accounting.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.evaluation.dreval_normalize import naive_predictions  # noqa: E402

OOF = ROOT / "notebooks/outputs/panel/panel_oof_predictions.csv"
OUT = ROOT / "notebooks/outputs/dreval/bias_accounting.csv"

#: The arm docs/steps/05 quotes the metric-equivalence row for. Recorded, not re-decided here.
QUOTED_ARM = ("X_pca", 0.0, "mse")


def fragility_proxies(oof: pd.DataFrame) -> dict[str, pd.Series]:
    """Per drug, each cell line's mean ``y_true`` over the **other** panel drugs.

    Built from the truth alone, so it is identical across arms and seeds — a property of the labels,
    not of any model. Leave-one-drug-out is what keeps a drug out of its own proxy; without it the
    proxy would contain the very column it is used to control for.
    """
    truth = oof[["drug", "cell_line", "y_true"]].drop_duplicates()
    tot = truth.groupby("cell_line")["y_true"].agg(["sum", "count"])
    out = {}
    for drug in truth["drug"].unique():
        own = truth[truth.drug == drug].set_index("cell_line")["y_true"].reindex(tot.index)
        out[drug] = (tot["sum"] - own.fillna(0)) / (tot["count"] - own.notna().astype(int))
    return out


def partial_rho(y_true, y_pred, proxy) -> float:
    """Rank-partial correlation of truth with prediction, controlling for the fragility proxy."""
    ry, rp, rf = rankdata(y_true), rankdata(y_pred), rankdata(proxy)
    typ, tyf, tpf = (np.corrcoef(a, b)[0, 1] for a, b in ((ry, rp), (ry, rf), (rp, rf)))
    return float((typ - tyf * tpf) / np.sqrt((1 - tyf ** 2) * (1 - tpf ** 2)))


def score_one_run(run: pd.DataFrame, proxies: dict[str, pd.Series]) -> dict:
    """Both tables for a single (rep, alpha, loss, seed) set of out-of-fold predictions."""
    run = run.reset_index(drop=True)
    scored = run.assign(y_naive=naive_predictions(run))

    def rho(a, b):
        return spearmanr(a, b).statistic

    per_drug, per_drug_norm, frag, partial, n_naive = [], [], [], [], []
    for drug, g in scored.groupby("drug"):
        per_drug.append(rho(g.y_true, g.y_pred))
        per_drug_norm.append(rho(g.y_true - g.y_naive, g.y_pred - g.y_naive))
        n_naive.append(g.y_naive.round(10).nunique())
        f = proxies[drug].reindex(g.cell_line).to_numpy()
        ok = np.isfinite(f)
        frag.append(rho(g.y_pred.to_numpy()[ok], f[ok]))
        partial.append(partial_rho(g.y_true.to_numpy()[ok], g.y_pred.to_numpy()[ok], f[ok]))

    return {
        "spearman_raw_pooled": rho(scored.y_true, scored.y_pred),
        "spearman_norm_pooled": rho(scored.y_true - scored.y_naive, scored.y_pred - scored.y_naive),
        "spearman_raw_per_drug": float(np.mean(per_drug)),
        "spearman_norm_per_drug": float(np.mean(per_drug_norm)),
        "rho_pred_fragility": float(np.mean(frag)),
        "partial_rho": float(np.mean(partial)),
        "n_distinct_naive": int(np.median(n_naive)),
    }


def main() -> None:
    oof = pd.read_csv(OOF)
    proxies = fragility_proxies(oof)
    keys = ["rep", "alpha", "loss", "seed"]
    rows = [{**dict(zip(keys, k)), **score_one_run(g, proxies)} for k, g in oof.groupby(keys)]
    if not rows:
        raise SystemExit(f"no runs found in {OOF} — refusing to write an empty artifact")

    per_seed = pd.DataFrame(rows)
    arms = per_seed.groupby(["rep", "alpha", "loss"]).mean(numeric_only=True).drop(columns="seed")
    arms["n_distinct_naive"] = arms["n_distinct_naive"].round().astype(int)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    arms.round(4).to_csv(OUT)
    print(f"{len(per_seed)} runs -> {len(arms)} arms | wrote {OUT.relative_to(ROOT)}\n")

    print("=== Table 1 — DrEval's normalized metric vs ours, on our own predictions (LCO) ===")
    print(arms[["spearman_raw_pooled", "spearman_norm_pooled",
                "spearman_raw_per_drug", "spearman_norm_per_drug",
                "n_distinct_naive"]].round(4).to_string())
    q = arms.loc[QUOTED_ARM]
    print(f"\n  quoted arm {QUOTED_ARM}:  pooled {q.spearman_raw_pooled:.4f} -> "
          f"{q.spearman_norm_pooled:.4f} | per-drug {q.spearman_raw_per_drug:.4f} -> "
          f"{q.spearman_norm_per_drug:.4f}")
    print(f"  distinct fitted baseline values within one drug: {int(q.n_distinct_naive)} "
          f"(one per fold — no cell-line term, which is the whole point)")

    print("\n=== Table 2 — how fragility-driven is each arm? ===")
    print(arms[["rho_pred_fragility", "partial_rho", "spearman_raw_per_drug"]].round(4).to_string())
    print("\n  mean over the six arms (share = how much of the raw score the control costs):")
    m = arms.groupby("rep")[["rho_pred_fragility", "partial_rho", "spearman_raw_per_drug"]].mean()
    m["fragility_share"] = 1 - m.partial_rho / m.spearman_raw_per_drug
    print(m.round(4).to_string())
    margin = arms.reset_index().pivot(index=["alpha", "loss"], columns="rep", values="partial_rho")
    gap = (margin["X_pca"] - margin["X_scGPT"]).round(4)
    print(f"\n  partial-rho margin X_pca - X_scGPT, per arm: {gap.min():+.4f} to {gap.max():+.4f}, "
          f"positive in {int((gap > 0).sum())}/{len(gap)}")


if __name__ == "__main__":
    main()
