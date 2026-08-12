# Gate 5 — the rerun, as it happened

**Started 13.08.2026.** Authorised by Selin, who lifted the [03.08 freeze](TODO.md) for this run:
*"the rerun should execute the notebooks directly … if a bug arises, fix it. if a bigger decision
arises, do a decision by yourself first to make the pipeline run and document it in the gate 5 report
for me to double check against the other options that were there."*

**How it is being run.** Each notebook is executed **in place** with
`jupyter nbconvert --to notebook --execute --inplace`, so its outputs live in the notebook and are
readable there. Nothing is run through a scratch script and nothing is staged in a temporary
directory — the notebook and its stored output *are* the record.

**The order is Gate 4's**, recorded in [TODO](TODO.md): the chain first and strictly sequential,
then everything else.

---

## Decisions taken to make the run proceed

*Each is one I made on Selin's instruction to decide rather than block. **Every one lists the
alternatives it was chosen over**, so it can be overturned on the merits rather than re-derived.*

### D1 · `data/` symlinked into the worktree

**The problem.** `drug_catalog` and `2_drug_selection` read `<repo>/data/drug/` and
`<repo>/data/GDSC2_*.xlsx`. That directory is gitignored, so it exists in the main checkout and
**not** in this worktree — the same fact that made me wrongly report GDSC2 as missing during Gate 4.

**Chosen:** a symlink `data -> /Users/selin/PycharmProjects/OncoTox/data`, excluded via
`.git/info/exclude` so it can never be committed. `/data/` in `.gitignore` is anchored *and*
directory-suffixed, so it does not match a symlink — checked, not assumed.

**Over:** *(a)* running the whole chain in the main checkout, which would put a many-hour run on top
of Selin's working tree and lose the isolation the worktree exists for; *(b)* copying `data/drug/`
in, which would fork the catalog `drug_catalog` writes and leave two copies to reconcile.

**Cost:** the run writes to the main checkout's `data/drug/`, which is where those files belong and
where `2_drug_selection` reads them from — so this is the real location, not a staging one.

### D2 · PCA keeps 512 components, and what it retains is reported (Selin)

**Chosen:** keep 512 — the width scGPT's embedding has — and report the retained variance as a
number rather than as a gate.

**Over:** *(b)* setting the count from the variance curve (an elbow, or 90 % retained), which is
data-driven but makes the two arms differ in **width as well as representation**, so any measured
difference confounds the two — the one thing the design exists to separate; *(c)* keeping 512 and
saying nothing about what it retains, which is the status quo the question was raised against.

**Cost:** components far down the spectrum may be near-noise, and the input dropout deletes them at
the same rate as PC1. That is D3's subject, not this one.

### D3 · Stage 6's veto fires on a comparison, not a threshold (Selin)

**Chosen:** the confound veto fires when the confounds explain **as much or more** of the within-line
variation as the reproducible signal does — stage 6's adjusted R² against stage 2's cross-seed
agreement, both measured within line, both produced by this run. **No constant is chosen.**

**Over:** *(a)* a fixed adjusted-R² bar (0.05, 0.10, …), arbitrary, and this project has already
retracted three numbers of exactly that kind; *(c)* reporting the value and never vetoing, which
leaves every positive Q2 result permanently provisional on a check that was run but not used.

**Why it works where a null does not:** a permutation null answers *"is the confound effect
non-zero"*, which with hundreds of cells per line is always yes. The question the veto actually asks
is *"is it as large as the thing we are claiming"*, and that is a comparison between two quantities
the run measures on the same scale.

### D4 · The three test aggregations are deferred (Selin)

Stages 7, 1 and 2 run as written (Benjamini–Hochberg at FDR 0.05; the paired fraction for stage 1).
All three are computed from the **saved per-cell prediction arrays**, so they can be re-decided after
the distributions are visible without retraining anything.

### ⬜ D5 · `input_dropout` — deferred to the variance numbers (Selin)

**Not yet decided, deliberately.** `input_dropout=0.1` removes 10 % of input variance in expectation
from **both** arms, so it is matched in value. It is not matched in *distribution*: scGPT's 512
dimensions are comparable in magnitude, so each draw removes about the same small amount, while
PCA's are variance-ordered, so **one draw in ten deletes PC1 alone**.

`3_representations` writes `uns["pca_fits"]["variance_ratio"]`, which says exactly how much of the
signal that is. The choice is made against that number rather than by judgement:

* **(a) keep 0.1 on both arms** — status quo. ⚠️ It handicaps **PCA, the control**, so it flatters
  the project's own hypothesis. One argument for it is already dead: *"comparability with previous
  runs"* has no force, because every previous run is void on target and panel grounds.
* **(b) set 0 on both arms** — removes the asymmetry outright; trunk dropout 0.5 remains.
* **(c) scale per arm** to equalise the variance removed — matched in effect, but needs a constant
  nobody has sourced.

Selin's instruction, 13.08.2026: decide it once the number exists; if she is unavailable, take it
alone and document the alternatives, since **retraining costs far less than `3_representations`** and
the choice can be revisited.

---

## Execution log

| # | notebook | status | notes |
|---|---|---|---|
| 1 | `1_data` | ✅ **clean** | both variants rebuilt, 53,513 cells × 198 lines. UMI covariates attached, join check **r = 1.0000**. hvg5000 22,722→5,000 genes; 1,129 symbols renamed, 23 collisions withheld |
| 2 | `analysis/harmonization/drug_catalog` | ✅ **clean** (2nd attempt) | **7,120 catalog rows** = 7,415 − 295 GDSC, exactly as predicted. PRISM 6,575 + CTRPv2 545. First attempt failed at §6 — see **B2**, and **M1** for why it looked like success |
| 3 | `2_drug_selection` | ✅ **clean** | **panel unchanged: the same 11 drugs**, so removing GDSC did not touch selection — as predicted from `drug_annotation`'s `dataset == "CTRPv2"` filter. **181 overlapping cell lines** (the `H292` alias fix, 180 → 181). 57 of 150 FDA drugs screened by CTRPv2; coverage 91.2–98.3 % |
| 4 | `3_representations` | ⏳ running | scGPT embedding × 2 variants, then targets → splits → pca. The longest non-training step |
| 5 | `4a_percell_training` | — | |
| 6 | `4b_mil_training` | — | |
| 7 | `5_evaluation` §1 | — | |
| 8 | `analysis/qc/hvg_sweep_build` | — | |
| 9 | `analysis/qc/verify_variants` | — | |
| 10 | `analysis/qc/gene_symbol_rescue` | — | |
| 11 | `analysis/harmonization/drug_coverage` | — | |
| 12 | `analysis/harmonization/cell_line_join_verification` | — | |
| 13 | `analysis/evaluation/diagnostics` | — | |
| 14 | `analysis/evaluation/dreval_benchmark` | — | |
| 15 | `5_evaluation` §2/§3 | — | |

---

## Bugs found and fixed during the run

*Nothing in the chain had ever been executed — zero executed cells across all six numbered
notebooks — so every runtime defect below is being seen for the first time.*

| # | where | what | fix |
|---|---|---|---|
| **B1** | `4b_mil_training` §3.8 and §4.1 | **The cell-cycle columns are `G1/S_score` and `G2/M_score`, not `G1`/`G2`.** Stage 6 would have raised `KeyError` on its first run, and §4.1 appended two non-existent columns to `PROGRAMS` while double-counting the two it was trying to add. | Both corrected to the real names. |

**How B1 survived every static check, which is the part worth keeping.** I took the names from
`list(h5file['obs'].keys())` while building 4b. **A forward slash in a column name becomes a group
hierarchy in an h5 file**, so `obs/G1` is a *group* containing `S_score`, and anything reading the
file's structure sees a column called `G1` that does not exist. The authority is the `column-order`
attribute, which says `G1/S_score` — and `anndata` reads it correctly, so the data was never wrong;
only my reading of it was. No signature check, path check or name check could see this: the code was
valid, the column was a plausible string, and the file did contain a `G1` object.

It is also a live warning rather than only a past one: `anndata` emits *"Forward slashes will be
disallowed in h5 stores in the next minor release"* on every write of these columns.
| **B2** | `drug_catalog` §6 | **`auc_cc` is not a column in DrEval's `CTRPv2.csv`.** It is the *project's* name for the measure; the file's column is `AUC_curvecurator`. My §6 repoint passed `paths.score` straight to `usecols=`, so it raised `ValueError: Usecols do not match columns … ['auc_cc']`. | Imports `ctrp_to_h5ad.SCORE_COLUMNS`, the one mapping between the project's score names and the file's columns, rather than restating either. |

### M1 · A defect in how I was running the notebooks, not in the notebooks

**`jupyter nbconvert … 2>&1 | tail -30` reports the exit code of `tail`, which is always 0.** So the
first `drug_catalog` run reported **success while having failed**, and I read its *stale* stored
outputs — from a pre-rewrite execution, still showing 7,415 catalog rows with 295 GDSC ones — as if
they were the run's. Had I not noticed that those numbers were impossible under the new code, a
failed notebook would have been logged as clean and the whole execution record would have inherited
the error.

It is the same failure this review has been finding all day, committed by me while checking for it:
**something reported success without having checked anything.** Notebooks are now executed with no
pipe, so the exit code is `nbconvert`'s.

⚠️ **Consequence for reading this log:** a notebook's stored outputs are only evidence of *this* run
if its execution actually succeeded. Stale outputs from an earlier execution survive a failed
`--inplace` run untouched and look identical to fresh ones.
