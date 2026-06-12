# Step 08 — Foundation model & clinical fine-tuning (overarching goal / horizon)

> **Status: ❌ NOT STARTED — long-horizon placeholder.** This is the project's **main goal** as
> stated in the plan, beyond the concrete sub-goals. It is the reason every earlier step is built
> the way it is. Documented here so the ultimate direction is visible, not because work is imminent.

*Part of [OncoTox project progress](../project_progress.md). Covers (future): turning the trained
predictor into a reusable, fine-tunable pan-cancer single-cell foundation model.*

---

## Main goal (from the project plan)

> "To develop a domain-specific pan-cancer single-cell foundation model capable of predicting
> pharmacological response (efficacy/toxicity), which can subsequently be **fine-tuned for specific
> cancer types and/or clinical datasets (binary clinical outcomes)**."

The continuous bulk labels (CTRPv2/PRISM/GDSC) are deliberately chosen as **pre-training** signal:
they capture a richer response distribution than binary cutoffs, so a model pre-trained on them
should yield a representation that **fine-tunes well onto scarce clinical data with binary
outcomes** (responder / non-responder).

## Why this sits on top of Steps 01–06

- The whole pipeline (frozen scGPT prior → supervised regression head, masked cross-database
  multi-task) exists to produce a **transferable representation**, not just a CTRPv2 predictor.
- The **inference payoff** the plan describes: run the trained model on individual cells to predict
  a **distribution of sensitivities within one sample**, computationally flagging rare
  naturally-resistant sub-clones before selection pressure — directly addressing how baseline tumor
  heterogeneity drives resistance.

## What this needs first

- A stable cross-database multi-task model ([Step 06](06-cross-database-integration.md)) — its
  `runs/<…>/best_model.pt` becomes the pre-trained **trunk**, with `OncoMLP`
  (`scripts/model/OncoMLP.py`) given a swappable clinical (binary) head.
- A new fine-tune entrypoint alongside `scripts/training/train_multitask.py`, plus a clinical
  dataset loader analogous to `scripts/model/dataset.py`.
- Access to a clinical dataset with binary outcomes (not yet in the project — the plan notes
  large-scale clinical data is the standing bottleneck).

## Open questions (long horizon)

- Fine-tuning protocol: freeze trunk + retrain head, vs. full fine-tune on clinical data?
- How to bridge in-vitro cell-line training to in-vivo clinical tumors (domain shift).
- Single-cell inference + sub-clone distribution analysis as a downstream evaluation.

## Definition of done

- Pre-trained model checkpoint packaged as reusable starting weights.
- At least one fine-tuning / transfer demonstration documented here.
