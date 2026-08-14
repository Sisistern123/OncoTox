"""What the model predicts *before* it is trained — and how that compares to the label's own spread.

**Why this file exists.** `docs/steps/03` and `OncoMLP.init_head_bias_` both carry the warning that
initializing the head bias at the per-drug mean starts the model near the null predictor's *level*
but not *at* the null predictor, because the head's weight rows are still random. That warning was
quantified by a script called ``init_spread.py`` that **was never committed anywhere** — not on
`main`, not on any branch (`git log --all -- '*init_spread.py'` is empty). It is the second instance
of that defect after ``arch_facts.py``; both are recorded in `docs/TODO.md`.

Re-deriving on 14.08.2026 reproduced the head row norm and both mean predictions but **not** the
prediction spread, which the docs gave as ~0.31 and which measures ~0.37 here. Rather than swap one
unsourced number for another, the measurement is committed so every figure in that paragraph is
checkable from now on. Selin's decision, 14.08.2026.

**What it measures, and the choices in it** (stated because they are choices, not defaults):

* **Synthetic ``randn`` input, not real embeddings.** The docstring being checked says "synthetic
  unit-variance input", so that is what is used. It matters: the two real representations differ in
  input scale by ~100×, so the initial scatter on real input is arm-dependent and this number is a
  stand-in for both. Anyone quoting it should say "on unit-variance input".
* **``hidden_dims=(128, 64)``** — `DEFAULT_HIDDEN_DIMS` for both arms. `(64, 32)`, the class default,
  was tried and reproduces neither the means nor the row norm, so it was not the original setup.
* **Seeds 42 and 0**, the two the original claim quoted.
* **``.eval()``**, so dropout is off: the claim is about what the network emits at initialization,
  not about a training step.
* **A uniform requested mean of 0.90**, matching the claim. Real head biases are initialized per drug
  (`cv.per_drug_line_mean`); a single value is used here so the spread is not confounded by the
  spread *between* drug means, which is a property of the labels rather than of the initialization.

The label-side comparison is read from the panel artifacts rather than asserted, so the ratio this
paragraph turns on cannot go stale silently.

Run:  .venv/bin/python scripts/evaluation/init_spread.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[2]
import sys

sys.path.insert(0, str(ROOT))
from scripts.model.OncoMLP import DEFAULT_HIDDEN_DIMS, OncoMLP, init_head_bias_  # noqa: E402

HIDDEN = DEFAULT_HIDDEN_DIMS["X_pca"]
SEEDS = (42, 0)
N_CELLS = 20_000
REQUESTED_MEAN = 0.90


def _panel() -> list[str]:
    return pd.read_csv(ROOT / "notebooks/outputs/panel/panel.csv")["drug_key"].tolist()


def init_stats(seed: int, n_heads: int) -> dict:
    """Row norm, prediction spread and mean prediction of an untrained, bias-initialized model."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = OncoMLP(input_dim=512, hidden_dims=HIDDEN, output_dim=n_heads)
    init_head_bias_(model, np.full(n_heads, REQUESTED_MEAN))
    model.eval()
    with torch.no_grad():
        pred = model(torch.randn(N_CELLS, 512))
    head = [m for m in model.modules() if isinstance(m, torch.nn.Linear)][-1]
    return {
        "seed": seed,
        "head_row_norm": float(head.weight.detach().norm(dim=1).mean()),
        "pred_sd": float(pred.std()),
        "pred_mean": float(pred.mean()),
    }


def label_spread() -> pd.Series:
    """Per-drug spread of ``auc_cc`` **across cell lines** — what the initial scatter is compared to.

    Taken from the committed out-of-fold predictions' ``y_true`` rather than from a summary column,
    so it is the same labels every other panel number is scored against. One arm is selected purely
    to pick one row per (drug, line); ``y_true`` does not depend on the arm.
    """
    oof = pd.read_csv(ROOT / "notebooks/outputs/panel/panel_oof_predictions.csv")
    one = oof[(oof.rep == "X_scGPT") & (oof.alpha == 0.0)
              & (oof.loss == "mse") & (oof.seed == 42)]
    return one.groupby("drug").y_true.std(ddof=1).reindex(_panel())


def main() -> None:
    panel = _panel()
    rows = [init_stats(s, len(panel)) for s in SEEDS]
    df = pd.DataFrame(rows)
    print(f"OncoMLP(input_dim=512, hidden_dims={HIDDEN}, output_dim={len(panel)}), "
          f"head bias initialized to a uniform {REQUESTED_MEAN}, .eval(), "
          f"{N_CELLS:,} synthetic unit-variance cells\n")
    print(df.round(4).to_string(index=False))

    spread = label_spread()
    print(f"\nper-drug across-line spread of auc_cc (panel, n={spread.size}): "
          f"mean {spread.mean():.4f}, median {spread.median():.4f}, "
          f"range {spread.min():.4f}-{spread.max():.4f}")
    print(f"\ninitial scatter / label spread = "
          f"{df.pred_sd.mean():.4f} / {spread.mean():.4f} = {df.pred_sd.mean()/spread.mean():.2f}x")
    print("\nThe claim this replaces gave the scatter as ~0.31 and the label spread as ~0.17,")
    print("i.e. ~1.8x. Both were unsourced; see this module's docstring.")


if __name__ == "__main__":
    main()
