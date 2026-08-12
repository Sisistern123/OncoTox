# OncoTox — TODO

> # ⛔ 28.07.2026 — THE CURRENT DRUG PANEL IS VOID
>
> **The 8-drug literature panel is discarded.** Its candidate list was ranked by `min(kill, spare)` on
> our own response values before the literature criterion was applied, so the selection inherited the
> label dependency it was built to remove — and 32 approved or clinical compounds, `nutlin-3` among them,
> were excluded for the wrong reason.
>
> **Everything computed on it is therefore provisional:** the step-1 run (`4a_percell_training`), the
> distributions and weighting design (`panel_distributions`), the dispersion figures (`diagnostics` §5), the panel rows in
> [Step 05](./steps/05-multitask-results.md), and the corresponding numbers in the report. Do not quote
> any of them, and do not build on them.
>
> ⚠️ **Release condition corrected 12.08.2026.** This used to read *"until the panel is rebuilt"*, and the
> panel **was** rebuilt on 12.08.2026 — which by the old wording would have made these numbers quotable
> again. It does not: they were computed on the *old* panel, so rebuilding it retires them rather than
> restoring them. The banner lifts when a run exists on
> [the new panel](./steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026),
> i.e. at **R4** of the sweep.
>
> **Nothing new is started before the review below is done.** The 15.07 progress report was postponed
> rather than presented on a panel we knew to be flawed; the point of that decision is lost if we patch
> around it.

> # ⛔ 03.08.2026 — NOTHING IS RE-RUN UNTIL SELIN'S REVIEW IS FINISHED
>
> **No re-embedding, no retraining, no recomputing of `X_pca` — nothing — until Selin has completed her
> own double-check of the pipeline and the repository and says so.** Code and documentation fixes
> continue; producing new numbers does not.
>
> **Why it matters now:** the 03.08 preprocessing pass changed what the code *produces* without
> re-running anything. `gen_embeds.py` seeds with 42 and runs on MPS; `add_pca.py` passes
> `random_state=42` to `sc.pp.pca`; the training `DataLoader`s take an explicit generator. So **every
> artifact on disk predates the code that now exists** — embeddings, `X_pca`, and every run under
> `runs/`. That divergence is recorded, not hidden: the affected numbers carry dated markers in
> [Step 02](./steps/02-preprocessing-and-embeddings.md),
> [Step 05](./steps/05-multitask-results.md) and
> [Corrections](./steps/corrections-and-dead-ends.md).
>
> **The order is deliberate.** Re-running first would burn hours regenerating artifacts from a pipeline
> the review may still change, and would destroy the inputs behind the committed results before they
> have been read. Decisions still open that would each force another recompute: gene scaling before PCA,
> the post-HVG renormalization, the HVG ranking scale, and the re-embedding scope.
>
> **Everything is regenerated in one clean sweep at the end of the review — decided 05.08.2026 (Selin).**
> Not piecemeal as each item closes. Two narrow exceptions were allowed during the audit, both for
> measurements whose values *cannot* depend on anything the review might change, and both read-only:
> `verify_variants.ipynb` §10a–§10c (input-length measurement) and `gene_symbol_rescue.ipynb`. Anything
> that regenerates a **committed artifact** waits for the sweep — including `verify_variants` §7/§8,
> which rebuild the UMAPs and overwrite `outputs/embeddings/umap_cancertype_pca_vs_scgpt.png`. Anything
> that **trains** waits unconditionally, `verify_variants` §9 among it.

> # ⛔ 12.08.2026 — NUMBERS AND CLAIMS CLEARED FROM THE DOCS AND THE REPORT
>
> **What this was.** A docs-vs-code audit, asked for by Selin: *are the docs and the report up to date
> and non-redundant with the code, and are there numbers or hypotheses stated without a re-run?* The
> current pipeline was confirmed with the supervising session rather than read off the docs. Its ruling,
> which governs everything below: **nothing has been re-run and nothing on disk is current** — the only
> sanctioned exceptions since the 03.08 freeze are read-only measurements (`verify_variants` §10a–§10c,
> `gene_symbol_rescue`) and `2_drug_selection`, which reads the response CSV and no pipeline artifact.
>
> **What "cleared" meant here, since it needed a rule.** Taken literally, *every* number in the docs is
> unsupported by a current run, and clearing all of them would delete the record that
> [Corrections](./steps/corrections-and-dead-ends.md) and [Step 05](./steps/05-multitask-results.md)
> exist to hold — which Selin decided on 11.08.2026 to keep. So the rule applied was: **clear a number
> where it is presented as currently holding, and leave it where it is already marked as the record of
> what was believed.** Nothing under an existing ⛔/⚠️ banner was touched. If the intended scope was
> wider, say so — the wider pass is mechanical from here.
>
> **⚠️ These are removals, not re-derivations. Every item below is now *unevidenced*, not *answered*.**
> Several were load-bearing arguments; deleting the number does not settle the question it answered, and
> R6 has to write the replacement from regenerated artifacts rather than restore what was here.
>
> ### Cleared — and therefore open again
>
> | What was stated as live | Where | Why it could not stand |
> |---|---|---|
> | *"The gene-set size is not critical"* + the sweep's heads-beating counts and val MSEs | `03_methods.tex` §Data; [Step 02](./steps/02-preprocessing-and-embeddings.md) §Why HVG-5000 is the default, reason 1 | The sweep has **no live numbers** — [Step 05](./steps/05-multitask-results.md) already banners the same table as superseded. Step 02 restated it unmarked *and* as reason 1 for the live default, so **HVG-5000 now rests on reasons 2 and 3 alone** |
> | Pearson tracks Spearman within 0.02 | `03_methods.tex` §Evaluation | Empirical, on the retired target and voided panel |
> | Synthetic mean-effects predictor scores normalized ρ = **0.98** | `03_methods.tex` §Evaluation; `dreval_normalize.py` docstring; `scripts/archive/README.md` | **No code in the repo produces it and no artifact records it** — it exists only in prose, presented as "demonstrated". The *mechanism* is verifiable and was kept |
> | MPS band = 0.313 / 0.315 / 0.317 / 0.320 | `06_limitations.tex` | Void runs. That MPS is non-deterministic stands; the **width** of the band does not, and differences are currently being called interpretable against a number that no longer applies |
> | `pred_std` 0.53 / 0.47 against "a true spread of **1.0**" | `06_limitations.tex` | The 1.0 is `auc_z`'s unit variance — stated in units the pipeline no longer produces |
> | Density weighting "does not help" | `project_progress.md` | Void panel, and the deltas fall inside the MPS band. **Superseded the same day by audit 09**, which found the sharper fault — the metrics could not have seen the effect — and re-tests rather than retires it (item 9A). Audit 09 owns the `06_limitations.tex` correction; this audit's weaker note there was dropped on rebase |
> | *"averaging a line's cells loses nothing measurable"* / ridge ties the MLP / model-side tuning closed | `06_limitations.tex`, `project_progress.md`, scorecard | Only the network arm early-stopped on the fold it was scored on. The tie is an **upper bound on the network's side** — the honest ordering may favour ridge (item 8C) |
> | ~a fifth of the signal is the cell-line effect | `06_limitations.tex` | Void run, deleted normalization, **and cited at `outputs/dreval/…`, a path that does not exist** |
> | scGPT−PCA margin "was **not** sign-consistent across seeds" | `06_limitations.tex` | ⚠️ **Direct contradiction with [Step 05](./steps/05-multitask-results.md)** — both marked unsupported, neither resolved. Own item below |
> | MIL is *"the only untested capacity lever"*; must beat ridge to be worthwhile | `06_limitations.tex`, `project_progress.md` | **Already retracted 11.08.2026** in this file and left standing in both. Controls are a floor; what counts as a positive Q2 result is still **Selin's open decision**, to be fixed *before* the run |
> | +0.077 / 0.048 / 0.011 / 82× / ~78× on the index page | `project_progress.md` | Void, **and** the index is barred from holding numbers by its own conventions |
>
> ### ⛔ Open — the seed sign-consistency conflict, to be settled by the rerun and not by choosing now
>
> **Two documents say opposite things about the same three-seed check, and neither is re-derivable under
> the freeze.** Ruled 12.08.2026 (Selin): **do not pick a side.** This is corrected *after* the rerun,
> from the regenerated artifact.
>
> | Location | What it currently claims |
> |---|---|
> | `report/sections/06_limitations.tex` | the scGPT−PCA margin *"was **not** sign-consistent across seeds, and did not survive DrEval's normalized metric"* |
> | [Step 05](./steps/05-multitask-results.md#learnability-filtered-subset--the-signal-was-there-all-along-13072026) | *"Gap = +0.075 ± 0.038, **sign-consistent** across all three seeds"* — seeds 42/1/7, gaps +0.043 / +0.066 / +0.117 |
>
> **Both are marked unsupported in the meantime** — the report's characterisation was cleared on
> 12.08.2026, and Step 05's table already sits under that page's ⛔ void banner. Neither may be quoted.
>
> **Actionable at R4/R6, in that order.** R4 produces ≥ 3 seeds on the rebuilt panel and the current
> target (already blocking there). At R6, read the sign consistency off that run and write **one**
> statement into both locations, deleting the loser rather than reconciling the wording. The old runs
> cannot arbitrate it: they were on a retired target, a voided panel, representations predating the
> preprocessing corrections, and — for the MLP arms — an early-stopping leak that was not uniform across
> the two representations, which is itself a mechanism capable of flipping a small margin's sign.
>
> ### Not cleared — surfaced instead, because fixing them means deciding something
>
> - **180 vs 181 — owned elsewhere, deliberately untouched here.** Selin is deciding this one herself
>   (12.08.2026), so this audit changed **no** line-count number in the report and files no task for it.
>
> ### Fixed outright — code contradicted, no judgement needed
>
> - `notebooks/README.md` §4 and `4a_percell_training.ipynb` §B both claimed §B's outputs go to
>   **`outputs/matrix/`**, a directory that has never existed — the documentation of the `OUT_MATRIX`
>   defect fixed in `f6cbef4`, left behind when the code was fixed. Both now name
>   `outputs/archive/training_545_mean_pv/`, which is also where the previous run's artifacts sit, so §B
>   **overwrites in place**.
> - `report/README.md` documented `pdflatex` twice and *"no `bibtex` run needed — references are a manual
>   `thebibliography`"*. `07_bibliography.tex` is `\bibliographystyle` + `\bibliography{../references}`,
>   i.e. real bibtex, and `.bbl` is gitignored. **Verified in this worktree, which had no `.bbl`:**
>   `pdflatex` alone reports `No file main.bbl`; the four-step sequence gives 0 undefined references over
>   17 pages. Same defect `cf3ad3f` fixed in `.gitignore` and did not carry here.
> - `report/README.md` named `\rhoFullScgpt` / `\ridgePca` as example macros and described "two results
>   tables" — all removed with the 12.08 withdrawal — and gave `cp` commands refreshing
>   **`fig_rescue.png` and `fig_dreval.png`, which no `.tex` file references** and whose sources are a
>   dead notebook and a void run. `fig_umap.png` is the only figure the report uses.
> - **`\NIndep` named two quantities at once and is retired** (Selin's ruling, 12.08.2026). It was
>   defined as *"independent examples after held-out splits"* but carried the value **150**, which was a
>   ridge control's line count — a different quantity, and one the report no longer reports. Replaced by
>   two macros named for what each counts: **`\NLinesCV` = 153** (lines inside the cross-validation, test
>   held out) and **`\NLinesFit` = 126** (the fixed training split alone). **Every site citing the
>   effective sample size takes `\NLinesCV` — confirmed by Selin 12.08.2026, so it is a decision and not
>   a reading.** What settles it is two sentences of `03_methods.tex` §Evaluation: it names *"the 153
>   training and validation lines"*, then cites the effective sample size for that same protocol, so the
>   two cannot differ, and one macro then serves all four sites. It also matches the quantity's role —
>   every use of it says uncertainty is governed by lines rather than cells, and the metric reported is a
>   correlation across the lines CV predicts exactly once. **126 is a training-set size and no reported
>   uncertainty is taken over it.** The *"held-out split**s**"* plural was loose: under this protocol
>   there is no standing validation holdout at all, since `inner_holdout` draws its slice from each
>   fold's own training lines. ⚠️ **The choice of quantity is settled; the value is not** — both macros
>   are derived from `\NLines`, so 153 becomes ~154 at the sweep. Re-read at R6, never hand-adjusted.
> - **The two false statements about `dreval_benchmark` are corrected** — item 6's sub-bullet and R5,
>   with the R5 hold lifted by Selin for exactly this and nothing else. `dreval_normalize.py` is live at
>   `scripts/evaluation/`; `e804f07` fixed the `'auc'` literal; the real blocker is that three imported
>   functions were deleted with the cell-line-effect diagnostic (item 11). One blocker, not two.
> - [Step 02](./steps/02-preprocessing-and-embeddings.md) §HVG-5000 pipeline outputs stated the on-disk
>   counts with no marker; it now records that the symbol repair moves 4,576 → 4,704, the `H292` alias
>   moves 180 → 181 and every split size with it, and that **`add_pca.TRAIN_SPLIT_COLS` no longer writes
>   `X_pca_train_paclitaxel`** (verified: `TRAIN_SPLIT_COLS = ("split_ctrp",)`).
>
> ### ⛔ Frozen until R4 — four sentences that are TRUE TODAY and become false when the code changes
>
> **Filed here rather than in R4 deliberately** (agreed 12.08.2026): R4 is plan-of-record and read-only,
> and a list that needs permission to be written is a list that gets lost. This is the **audit-11 class**
> — code changed, docs did not — which cost us two stale claims within a day of it happening.
>
> **The rule these come from, worth stating once:** a sentence is safe to edit when a decision is taken
> if it describes **a decision**; it is frozen until the code changes if it describes **what the code
> does**. Item 7's decision on the three fits closes the cross-validation leak, but it closes it *at R4*,
> when training changes — so every sentence below is **accurate right now** and rewriting it early would
> describe code that does not exist. **Do not sweep these when the decision lands. Sweep them when R4
> runs.**
>
> | Where | The sentence | Why it is not editable yet |
> |---|---|---|
> | `docs/steps/02` §"Still not fixed: cross-validation" | *"`resolve_rep` leaves them on the all-cells `X_pca` … every CV number still carries it"* | True until training resolves the representation differently |
> | `scripts/model/dataset.py:26` (`resolve_rep` docstring) | *"still leaky, and documented as such"* | Same, and it is the code's own description of itself |
> | `docs/steps/02`, the `X_pca` table row | *"used for: UMAPs and latent-space validation — **and every CV run**"* | ⚠️ A **table cell**, not prose: a sweep for "leaky" or "not fixed" will not find it |
> | `report/sections/06_limitations_and_outlook.tex:81` | *"those remain on the all-cells decomposition, **and gene selection remains an all-cells step for both arms**"* | ⚠️ **Half-true at R4**, the nastiest of the four — the second clause stays true *permanently* under decision 1, while the first changes. Editing it as one unit will break the half that is right |
> | ⭐ `report/sections/06_limitations_and_outlook.tex:73-74` | *"**The bias runs toward the control, so any scGPT advantage measured this way is conservative**"*, and the `\revision` after it, *"both asymmetries … accumulate rather than offset"* | **The one that matters most: this is the retracted lower-bound claim, in its other home.** Valid today, because its premise — the rotation estimated over every cell — is what the code does. `fitc` dissolves the premise at R4 and the conclusion with it. §Methods no longer asserts a lower bound; **this passage still does**, so the retraction stays half-propagated until R4. Now carries a dated marker naming this exact sentence as the one to revisit. ⚠️ The **second** asymmetry in it, the gene-symbol one, survives R4 untouched — no fitting set affects which genes reach the model — so this is another sentence that must not be edited as one unit |
>
> ✅ **One defect in that passage was wrong *today* rather than at R4, and was corrected before merge
> (12.08.2026).** Line 67 read *"The scGPT embedding is **not fitted on this data at all**"*. It is:
> scGPT is fed our CPM matrix, and shares the gene set besides. Selin corrected that phrasing herself.
>
> **Why it was fixed here rather than deferred**, since the first instinct was to leave it: the freeze
> rule protects sentences that *accurately describe current code* and would become wrong at R4. **A
> sentence that is false today is not frozen, it is just false** — the rule was never meant to shelter
> it. And the timing settles what looked like a judgement call: **on main there is no contradiction**,
> because §Methods does not yet carry the corrected account. This branch *creates* the contradiction by
> introducing the accurate version alongside the false one, so repairing it is part of landing the
> change rather than cleanup taken on the way past. The conclusion two sentences later was left standing
> and flagged instead, which is a correction plus a dated marker — the repo's own convention — not a
> half-repair.
>
> ### ⛔ A second freeze class — accurate until Huber leaves the code, false immediately after
>
> **Different trigger, so it is listed separately.** The four above turn false at **R4**, when the
> representations are regenerated. These turn false at **a specific commit**: the model session is
> removing `--loss huber`, `TrainConfig.huber_beta`, the name check and the `smooth_l1_loss` branch
> (Selin, 12.08.2026 — Huber is dropped from the loss comparison *and* from the code). Whoever lands
> that commit should sweep these in the same change; they are not R4's problem and will be stale for
> however long the gap is.
>
> **Do not pre-emptively rewrite them.** The removal is on a branch. Until it lands, every sentence below
> is an accurate description of what the code does, and this is exactly the inversion worth naming: they
> were left alone *because* they were accurate, and the same accuracy is what makes them false the moment
> the option goes.
>
> ⚠️ **Each of the four is wrong twice over once that branch lands, not once.** The branch **adds MAE**
> as well as removing Huber (`291440e`, then `f16b3ec`), and **MAE is not on `main` today** — `main`'s
> code offers `{mse, huber}`. So the fix is a *substitution*, not a deletion, and a sweep that only
> strikes "Huber" leaves all four still wrong, now by omission. Written out so it is mechanical:
>
> | Where | Reads today (accurate against `main`) | Becomes |
> |---|---|---|
> | `docs/steps/03-model-and-training-design.md:14` | optimizes a masked **MSE** or **Huber** loss | **MSE** or **MAE** |
> | `docs/steps/03-model-and-training-design.md:199` | `sq = (pred − y)²` (MSE), or `smooth_l1_loss(beta=0.05)` for `--loss huber` | MSE, or `l1_loss` for `--loss mae` — `smooth_l1_loss` and `huber_beta` are both gone |
> | `docs/steps/03-model-and-training-design.md:429` | the CLI flag list, `--loss {mse,huber}` | `--loss {mse,mae}` |
> | `docs/project_progress.md:137` | *"fully supervised regression (masked **MSE/Huber**)"* | masked **MSE/MAE** — ⚠️ an index page, so it reads as current |
>
> *(Line numbers refreshed 12.08.2026 — the CLI site moved 415 → 429 when the loss-comparison note above
> was added to the same file. Re-grep rather than trusting them if that file moves again.)*
>
> ✅ The *decision* half is already done: Step 03's loss-comparison grid is **MSE / MAE × α ∈ {off, 0.5,
> 1.0}, six arms**, with the grounds recorded. Only the code descriptions are waiting.
>
> **Why these are not corrected in advance, stated because it looks like an omission.** There is **no
> defect on `main` today** — all four match `main`'s code exactly. Correcting them now would put the docs
> ahead of the code in *both* directions at once: claiming MAE exists when it does not, and denying Huber
> when it does. That is the same self-contradiction-on-merge that `06_limitations:67` was, except
> self-inflicted. **They land in the commit that merges the R4 training branch**, which is also the only
> moment at which the replacement column above becomes true.
>
> ⏸️ **The mechanism is not final, and the flag holds either way.** Decision 2 — *how* the CV PCA is
> fitted — was reopened on cost the same day it was taken: a per-fold fit at training time needs
> `paths.raw_h5ad` (~2.15 GB), which the training path has never opened, because the targets `.X` has the
> scGPT OOV genes dropped and would give a different gene set; and that cost compounds across R4's loss
> grid. It may be re-taken, or become a precomputed-at-R2 change instead. **The leak closes either way**,
> so all four sentences stop being true either way — only the date and the replacement wording move.
>
> ### Re-based onto `58fadd7`, and two of the clearings above were overtaken by it
>
> This audit was cut against `f6cbef4`; Selin committed **`58fadd7`** while it was running, splitting
> stage 4 into `4a_percell_training` / `4b_mil_training` and **pre-registering the Q2 criterion**. The
> branch was rebased and the notebook edit followed the rename. Two consequences the audit had to take
> back, because they are now wrong in the *other* direction:
>
> - **The Q2 success criterion is no longer an open decision.** It is fixed in
>   [`4b_mil_training.ipynb`](../notebooks/4b_mil_training.ipynb) §2 as four stages with distinct roles —
>   synthetic positive control (precondition), within-line spread (necessary condition), cross-seed
>   reproducibility against a shuffled-cell control (**the test**), confound regression (**veto**). What
>   remains open is **one number**, `Q2_CONTROL_THRESHOLD`, plus the stage-1 and stage-2 fractions
>   derived from it. Wherever this audit wrote "an open decision for Selin", read "pre-registered; one
>   threshold outstanding" — corrected in the report and in `project_progress.md`.
>   ⚠️ ***Agreed plan, Step 2* below still describes the whole criterion as undecided** and was not
>   touched by `58fadd7`'s path repointing. Pointer added there; the bullet itself is Selin's to retire.
> - **"Attention weights are the readout" is now wrong.** `58fadd7` chose **instance-level MIL, not
>   attention pooling** (Ilse, Tomczak & Welling, ICML 2018): every cell gets a predicted *response*,
>   not a *weight* over a pooled embedding, trading predictive performance for readability at the
>   instance. The report and `project_progress.md` both still described attention weights as the
>   clinically interesting readout; both corrected. The subpopulation-predictivity test that would have
>   used top-k cells is deliberately **not** in the criterion — selecting cells by their predicted value
>   and scoring them against the line's true response is biased by construction.

- start bei data download, schau genauer auf drug selection, suche publications dafür raus
- data harmonization genauer anschauen
      - bulk und sc annotation merge -- wie genau wurde es gemacht, ist es valide?
      - scGPT OOV check -- sind die wirklich OOV?
- scGPT check
      - wie genau wird es trainiert
      - wie werden die embeds generiert
      - kann man finetunen?
- redundanz, staleness, file overload
      - schauen, ob es redundanten code gibt, code duplication
      - schauen, ob es veralteten stale code gibt, der nirgends genutzt wird
      - restrukturierung der notebooks, archivierung von veraltetem code
      - schauen, ob zu viele files für nichts gemacht werden


## 🔍 28.07.2026 — full pipeline review, data download to evaluation

Walk the whole thing once, in order, deciding three things at each stage: **what is settled**, **what is
open**, **what has never actually been verified**. The last column is where today's problems came from —
every one of them was a step that looked settled and had never been checked.

> ### ✅ Gate 1 sweep — all 65 unchecked boxes read against the code (13.08.2026)
>
> **Ten of the sixty-five were not work.** Every box making a checkable claim about the code was
> verified against the source rather than believed. Five described work already done — the
> `Adam → AdamW` migration, both clauses of the per-drug-mean null, `dreval_benchmark`'s epoch count,
> and the 🔴 head-bias item — and five were written against a **retired target or a retracted
> criterion** (`z-score train-only` on the withdrawn `auc_z`, and four boxes scoped to the 5-drug
> learnability subset). All ten are now closed or marked, each keeping what it claimed in the past
> tense so the record of what was believed survives the correction.
>
> **The one worth learning from is the 🔴.** It was dated the same day it was swept, read *"NOT
> fixed — Selin's call"*, and the fix was already committed — in `8b6a678`, whose subject says so.
> A red flag that outlives its defect costs the same attention as a real one, and this is the second
> instance of the same shape in two days: `notebooks/README.md` row 5 sent readers to repair a
> notebook that already worked. **The failure mode is a list that is only ever appended to.**
>
> **What the sweep leaves genuinely open before the rerun is small:** `dreval` cell 8's missing
> `DataLoader` generator (benign, fragile),
> `NaiveMeanEffects` as the default baseline, and two decisions that are Selin's — where the ridge
> baseline lives, and whether this model needs weight decay at all.
>
> **Three items are unblocked by R2 rather than blocking it.** PCA's 512 components, the
> `input_dropout` asymmetry between the arms, and item 4A's input-scale measurement all resolve from
> `uns["pca_fits"]["variance_ratio"]`, which `add_pca.py` writes. They sit **between R2 and R4**.
> Everything else is either the rerun itself (R2–R6) or post-project.

- [x] **1 · Data acquisition — walked 10.08.2026.** Release, files, date and origin now recorded per
      source in [Step 01](./steps/01-datasets-and-harmonization.md#provenance--what-was-retrieved-from-where-when),
      together with the two data roots (`~/Desktop/OncoTox/data` and `<repo>/data`), which nothing had
      stated before. Found and fixed: GDSC2's release was only ever in its filenames (**27Oct23**);
      DrugBank had no date (**06.03.2026**, from the file mtime — flagged as such); the Repurposing Hub
      file is the **20200324** release. **CTRPv2 verified byte-identical to the published release** — all
      15 files match the MD5s in its `MANIFEST.txt`, the only source that shipped any. SCP542's licence
      question closed (item F) and its GEO accession **GSE157220** recorded. The two orphan downloads
      identified. *Settled, not open:* SCP542 exposes no version identifier, so the retrieval date is its
      version reference (Selin, 10.08.2026); no checksums generated for the three sources that shipped
      none; GDSC's `LN_IC50` processing stays undocumented — `GDSC_Raw_Data_Description.pdf` came with
      the download and may close it, but GDSC is not a modelling priority.
- [x] **2 · Harmonization — walked 10.08.2026, two defects fixed in code, replicate handling settled.**
      Full record in
      [Step 01](./steps/01-datasets-and-harmonization.md#the-join-audit--what-was-checked-and-what-held-10082026).
      **Verified:** no normalized-name collisions on either side, no rows lost to the merges, and — the
      check that had never been run — **no false matches**: 189 of the 190 name matches also agree on
      CCLE primary site, the exception being a line for which CTRPv2 records no site. **Fixed:**
      (a) CTRPv2 spells `NCI-H292` as `H292`, so the name join silently dropped a **screened** line —
      an explicit sourced alias (Cellosaurus `CVCL_0455`) takes the trainable overlap **180 → 181**,
      recovering 213 cells and 454 drug labels; (b) `v20.meta.per_experiment.txt` lists an experiment
      once per calendar day it ran, so 153 experiments were double-counted in the per-(line, drug) mean,
      moving 460 of `NCIH1299`'s targets and its rank on 427 of its 469 drugs. Both take effect **at the
      sweep only**; `\NLines` in `report/results_numbers.tex` stays at 180 until an artifact supports 181.
  - [x] **Replicate handling — settled 10.08.2026: keep averaging, but measure the disagreement.**
        2,637 of 81,626 (line, drug) pairs were screened twice (never three times, so median = mean).
        They come from just **six** cell lines, each re-screened against 534 of the 545 drugs. The two
        measurements differ by a median of **0.49×** the drug's spread across cell lines, and **27.3 %
        of them differ by more than that full spread**. Quantified in
        `notebooks/archive/replicate_variation.ipynb` →
        `outputs/archive/replicate_variation.{png,csv}`; written up in
        [Step 01](./steps/01-datasets-and-harmonization.md#genuine-repeats-are-averaged-and-they-disagree-more-than-the-targets-own-spread-10082026).
        **Six of 181 lines is not a random sample**, so this bounds nothing numerically — it says only
        that a substantial share of the target is screening noise, which items 5 (target), 6 (drug
        selection) and 11 (evaluation) should each read before drawing conclusions from a modest ρ.
  - [x] ~~**B · Migrate to persistent identifiers.**~~ **Dropped 10.08.2026 (Selin)** — the premise did
        not survive the audit and the join it was meant to protect is now verified.
        [Dead ends](./steps/corrections-and-dead-ends.md#migrating-the-cell-line-join-to-persistent-identifiers).
- [x] **3 · Preprocessing — answered by the 05.08.2026 transform audit; closed 10.08.2026.** The facts
      this item asked for are all in [Step 02](./steps/02-preprocessing-and-embeddings.md): CPM arrives
      already library-size normalized and the second normalization that had crept in was removed; the
      log step is the dataset authors' own `log2(1 + CPM/10)`; HVG is selected **once**, at `convert`,
      ranked on a log copy and applied to the CPM original
      ([what it removes and what `.X` holds](./steps/02-preprocessing-and-embeddings.md#what-hvg-filtering-removes-and-what-x-holds-at-each-stage),
      including the per-step table); the count of 5,000 has its own justification; and the scGPT OOV
      set was found to be mostly a symbol-matching defect, repaired in code (item A below).
      **The question it posed — does the HVG set depend on all cells including test? — is answered:
      yes**, and so do `sc.pp.scale` and the all-cells rotation
      ([what transform PCA sees](./steps/02-preprocessing-and-embeddings.md#what-transform-pca-sees--corrected-05082026)).
      Establishing that closes item 3; **deciding what to do about it is item 7**, where it sits with
      the other two open fits.
  - [x] ~~**A · Apply the gene-symbol repair — BEFORE the clean sweep.**~~ **Done in code 05.08.2026**
        (`scripts/annotation/gene_symbols.py`; `scp542_conversion.py` annotates,
        `gen_embeds.py::resolve_gene_names` resolves, own symbol first so nothing embedded today is
        lost). Takes effect at the sweep: **4,576 → 4,704** genes for `hvg5000` and
        **20,570 → 21,332** for `all_genes`. Three decisions and what each rejects:
        [Corrections](./steps/corrections-and-dead-ends.md#scgpt-discarded-genes-that-are-in-its-vocabulary-under-their-current-symbols).
        Source table: `reference/hgnc_complete_set.txt`. (FAIRER: **I**)
    - [ ] **Open:** `gene_symbol_rescue.ipynb` and its artifact `gene_symbol_rescue.csv` predate the
          reassignment guard (decision 3), so the CSV still counts `RNU12` and `EPB41L4A-AS2` as
          rescuable — 773 rather than 775, no effect at the precision anything quotes. Re-running the
          notebook is read-only and takes seconds, but it **overwrites a committed artifact**, so it
          waits for the sweep unless released separately.
- [x] **4 · Representations — walked 10.08.2026. The one strand that cannot be answered by reading has
      been moved out.** Three strands: whether the PCA fit is recoverable at all (**B**, closed
      10.08.2026), what scGPT is fed (**C**, closed 10.08.2026), and the **~78× input-scale asymmetry
      between the two arms under one shared learning rate** (**A**). A asks what a scale difference does
      to a *trained* model, so no amount of code-reading settles it — it is scheduled with the runs that
      can, under [After the sweep](#after-the-sweep--the-one-review-item-that-needs-new-runs). Nothing
      in the review still waits on it.
  - [x] **B · The PCA fit is now stored, not just the coordinates — done in code 10.08.2026.**
        `obsm["X_pca*"]` held only where each cell landed; the loadings, the variance ratios and the
        standardization statistics were computed and thrown away, so the pipeline could answer neither
        *what fraction of variance does PCA(512) retain* nor *which genes dominate PC1*, and could not
        project a new cell into an existing space — which the cross-database and XAI stages both need.
        `add_pca.py::_pca_record` writes all of it to `uns["pca_fits"]` per key; `varm` was not an option
        because the targets file's gene axis differs from the PCA's — the `scgpt` step drops OOV genes
        from `.X` while PCA keeps the full HVG set, and *how many* that leaves moves with the
        gene-symbol repair (item 3A), which is why the record carries its own gene vector.
        Structure, the reprojection formula and how it was checked:
        [Step 02](./steps/02-preprocessing-and-embeddings.md#the-fit-is-stored-not-only-the-coordinates-10082026).
        Costs ~10 MB (`hvg5000`) / ~47 MB (`all_genes`) per fit. Takes effect at the sweep.
    - [x] **Found on the way: the two PCA fits were standardized differently.**
          `_pca_fitted_on_train` used numpy's `ddof=0`, `sc.pp.scale` uses `ddof=1`. Worth under 0.01 %
          on this atlas, so no conclusion moves — but the two **PCA fits** were not produced by
          identical code. **Harmonized to `ddof=1` (Selin, 10.08.2026)**, leaving *which
          cells are seen* as the only difference between the fits. Changes `X_pca_train_*` in its last
          digits, at the sweep.
    - [x] ~~Count how many values actually reach the ±10 clip.~~ **Not needed (Selin, 10.08.2026)** —
          the cap is Seurat's `ScaleData(scale.max = 10)` default and stays at the default; measuring
          how often it binds would not change it.
  - [x] **C · What scGPT is fed — settled 10.08.2026, no measurement needed.** This item asked
        "raw counts vs CPM — never confirmed". **It was never a choice:** SCP542 distributes only
        `CPM_data.txt`, so raw counts do not exist in this project. What remained was whether feeding
        CPM puts us off-distribution from a model pretrained on binned counts, and it does not: bin
        edges are quantiles of each cell's own non-zero values, so any strictly monotone transform
        moves values and edges together and — at a given RNG state, see the tie-breaking below — the
        bins are unchanged
        ([Step 02](./steps/02-preprocessing-and-embeddings.md#what-scgpt-is-fed-and-why-its-scale-does-not-matter)). Also
        established there: binning **breaks ties at random** (`scgpt/preprocess.py:239`), which is what
        the `np.random.seed(42)` in `gen_embeds.py` is for. The docs' earlier claim that this had been
        *measured* on 200 cells was retracted — no code behind it, and nothing to compare against:
        [Retracted claims](./steps/corrections-and-dead-ends.md#the-scgpt-binning-invariance-was-verified-on-200-cells).
- [x] **5 · Target — ANSWERED 11.08.2026 by audit 05; marked closed 13.08.2026.** The target moved to
      **DrEval's reprocessed CTRPv2 (`auc_cc`)**, re-fitted from CTRPv2's raw dose-response with
      CurveCurator and normalised per replicate against the no-drug control. Full record:
      [Step 01](./steps/01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026).
      The original AUC-vs-EC50-vs-Emax question is answered by that move rather than by a separate
      comparison: the objection was that our own AUC conflated potency with efficacy over a
      0.13–600 µM range, and the replacement is a curve fit rather than a trapezoid over sampled
      concentrations.
      ⚠️ **This item stayed unchecked for two days while its body described a pipeline that no longer
      existed**, which is why it is recorded rather than quietly ticked. Two things it said were wrong
      by the time anyone read it: the report was said to *"state raw `auc`"* — the target is `auc_cc`;
      and it named as a *"related open leak"* the per-drug target mean/std computed over every cell
      line including val and test. **That leak no longer exists.** It was the `auc_z` per-drug
      z-scoring, retired 27.07.2026, and `ctrp_to_h5ad.py` records that the `ctrp_score_center` /
      `ctrp_score_scale` keys which existed only to invert it went with it.
      **Left open deliberately, and routed:** the winsorizing threshold is retired (`DEFAULT_WINSOR`,
      11.08.2026 — the benchmark applies no clipping "to maintain comparability to previous studies and
      avoid data loss"), and "are all statistics per fold?" was answered by items 7 and 8, not here.
- [x] **6 · Drug selection — REBUILT 12.08.2026. The panel is 11 drugs.** Full record in
      [Step 01](./steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026);
      produced by `notebooks/2_drug_selection.ipynb` → `notebooks/outputs/panel/panel.csv`.
      **The criterion changed from what this item asked for (Selin, 12.08.2026):** not "coverage and
      spread", but **FDA approval + a verified published claim**, with coverage as the only property of
      our own data that enters. Spread was dropped because it is still *our* label statistic, so
      selecting on it stays label-dependent — the objection that voided both previous panels, only
      subtler ([Corrections](./steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability)).
      The replicate-noise reuse this item suggested therefore has nothing to apply to.
      **Found on the way, and both fixed in code:** the **drug**-name join lost 102 of 545 compounds, 15
      of them single-agent FDA/clinical, now joined on `master_cpd_id` (545/545); and **cisplatin** was
      invisible under CTRP's name `Platin` with DrEval's `pubchem_id` pointing at *elemental platinum*
      ([Step 01](./steps/01-datasets-and-harmonization.md#the-drug-name-join-and-the-compounds-it-hid-12082026)).
      Of the expected re-entries, `oxaliplatin` and `bortezomib` are candidates but not panel members,
      and `nutlin-3` is excluded on the criterion itself — it has never been FDA-approved, which is a
      defensible reason where the old gate's was not. **Takes effect at the sweep:** every result still
      on record was computed on the void 8-drug panel.
  - [ ] **Left open by this item, routed elsewhere.** The 28.07 panel-void banner lifts only once a run
        exists on the new panel, i.e. at R4 of the sweep, not here.
        `notebooks/analysis/evaluation/dreval_benchmark.ipynb` cannot run: **three of the functions it
        imports from `scripts/evaluation/dreval_normalize.py` were deleted on 12.08.2026** with the
        cell-line-effect diagnostic, so its import cell raises.
        ✅ **DECIDED 12.08.2026 (Selin): rewire the import cell to DrEval's own recipe — `load_oof` +
        `normalized_evaluation`.** The deleted fragility diagnostic is **not** restored; it stays
        retired and is recoverable at `bf93084` if it is ever wanted. Two consequences beyond fixing
        the import: the notebook stops **re-training 20 models of its own**, and it scores the
        line-level out-of-fold predictions `4a_percell_training` already writes — so it benchmarks the
        model this project actually produces rather than a re-fit of it, against the same predictions
        every other evaluation reads. Rejected: restoring the diagnostic under its own name outside a
        file called after DrEval, which keeps a capability nothing currently asks for and leaves the
        benchmark re-training regardless.
        ⚠️ **Check the interaction with item 10's epoch fix before writing either.** If the rewire
        removes training from this notebook entirely, then setting `epochs=50` on a config nothing
        trains with is dead code that reads as meaningful, and the epoch defect is closed *by the
        rewire* rather than by the fix — a closed-by-accident, to be labelled as one rather than
        quietly disappearing.
        **This is on the critical path**: it unblocks `5_evaluation`, which must be **authored before
        R4 runs**, because audit 09's finding is that the loss comparison and the capacity
        re-derivation both need their metric set fixed before the run rather than after seeing one.
        *(Corrected 12.08.2026 — this said the notebook "imports the now-archived `dreval_normalize.py`
        and hardcodes the removed `'auc'` score, so it is broken twice over". Both were false: the module
        is **live** at `scripts/evaluation/dreval_normalize.py`, restored paper-only the same day, and
        the `'auc'` literal was fixed in `e804f07`. One blocker, not two.)*
- [x] **7 · Splits — walked 12.08.2026; CLOSED 13.08.2026 when the three fits were decided.** 7A fixed
      in code, 7B accepted rather than patched, 7C routed to item 11, and the three fits handed over
      from item 3 — what may a fit see — are now decided by Selin (HVG stays all-cells; the CV PCA is
      fitted per fold on `fitc`; the never-training cells stay in). Nothing in this item is outstanding.
      Confirmed sound: grouping is by cell line everywhere, the
      fixed `test` set is outside CV by construction (`eligible_splits=("train","val")`) and **has never
      been used by anything**, the MLP and the ridge control share one partition through
      `cv.grouped_folds` rather than two that agree by seed, the per-fold statistics are fitted inside
      the fold, and `X_pca_train_ctrp` sees only `split_ctrp=="train"` cells. What was *not* standard is
      written up in [Step 03](./steps/03-model-and-training-design.md#what-is-and-is-not-standard-about-the-cross-validation).
  - [x] **A · Early stopping ran on the fold being scored — FIXED IN CODE 12.08.2026.** `train_model`
        restores the lowest-validation-MSE checkpoint, and both `cv.oof_predictions` and
        `train_multitask.cv_evaluate` handed it the scored fold, so every out-of-fold prediction and CV
        metric was a minimum over epochs on its own evaluation data. **The identical defect was found
        and fixed in the DrEval benchmark on 14.07.2026 (`ee07b00`) and never carried into the code
        that produces the headline numbers.** Not uniform across the arms — selected epochs `[1,1,3,1,1]`
        (PCA) vs `[10,11,2,21,4]` (scGPT) — so it biased the comparison, not only the level; and the
        ridge control has no early stopping at all, so the MLP-vs-ridge gap was flattered too.
        Fixed by `cv.py::inner_holdout`: 15 % of each fold's training lines, grouped, become the
        early-stopping set (Selin, 12.08.2026 — 15 % chosen over reusing a neighbouring fold at 20 %;
        the fraction is arbitrary and documented as such, the inner seed is deliberately separate from
        `TrainConfig.seed`). Design: [Step 03](./steps/03-model-and-training-design.md#the-early-stopping-set-is-nested-inside-the-training-lines-12082026);
        record: [Corrections](./steps/corrections-and-dead-ends.md#the-same-val-split-leak-in-the-code-that-produced-everything-else).
        **Takes effect at R4.**
  - [x] **B · `splits/split_ctrp.csv` does not exist — accepted, not patched (Selin, 12.08.2026).**
        Never committed (`git log --all -- splits/` is empty), not on disk, not gitignored, while
        [Step 02](./steps/02-preprocessing-and-embeddings.md) and `report/sections/03_methods.tex` both
        stated it *was* versioned. So `frozen_split` has taken its redraw branch on every run and the
        guard has never been in force. **Both claims corrected.** Freezing the current assignment was
        considered and rejected: it is recoverable anyway — `outputs/archive/panel_void_8drug/panel_oof_predictions.csv`
        names all 153 train+val lines, so the test set is the labelled lines it omits — and every number
        scored on it is void on target and panel grounds. R2 creates the file itself; committing it
        there is where the guard starts to protect something (added to R2).
  - [x] **C · CLOSED 13.08.2026 (Gate 1 sweep) — both clauses were fixed and the box was not.**
        It read: *"The per-drug-mean null is computed two different ways. `cv_evaluate` fits its
        constant on the fold's fitting lines (honest); `4a_percell_training.ipynb` §4 computes
        `null_mse` from the variance of the **held-out** truth, an oracle constant fitted on the rows
        it is scored against … Also `_per_drug_train_mean` averages over cells, so lines weigh by
        their cell count."* Both were true when written and neither is now: `4a` §4's `_oof_null`
        fits each fold's constant on the rows that fold did **not** hold out (audit 11, 12.08.2026,
        and the cell carries the reasoning), and `_per_drug_train_mean` calls
        `cv.per_drug_line_mean` — per **cell line**, not per cell. Verified by reading both against
        the code, not by spot-check.
  - [x] **DECIDED 12.08.2026 (Selin) — what may a fit see. Handed over from item 3 (10.08.2026):
        three fits, one question.** Decided together, as the item required, so the three answers are
        consistent. The shared grounds: all three are **unsupervised** — no fit sees a response label.
        ⚠️ **Retracted 12.08.2026, same day, by decision 2's own resolution.** This line first read that
        all three fits bias *toward* the PCA control, "so any scGPT-over-PCA margin measured under them
        is a lower bound rather than an inflated one". **That does not survive `fitc`.** The lower-bound
        claim rested specifically on the baseline's fit being estimated over cells the model never
        trained on — which is exactly what restricting it to `fitc` removes. Gene selection stays
        all-cells (decision 1), but both arms receive the identical set, so it favours neither. Kept
        here rather than deleted because it is the kind of sentence that gets quoted into an abstract,
        and it would have been quoted on grounds this decision dissolved. **What survives is narrower
        and on a different footing:** the embedding reads only its in-vocabulary subset of the selected
        genes — 4,704 of 5,000 on `hvg5000` — which no fitting-set choice touches. The restriction to
        `fitc` is therefore about **attributability**, not conservatism: it is the one knob that could
        move a difference for a reason unrelated to the representations.
        **1 — HVG stays all-cells.** Keeping one gene set keeps folds, arms and Step 05's gene-set sweep
        comparable; a train-only HVG set would be fold-dependent, so "the 5,000 HVGs" would stop being a
        single object. **2 — the cross-validated PCA is fitted per fold, at training time**, on
        **`fitc` — the cells the model's weights are actually fitted on**, i.e. the fold's training side
        *minus* the 15 % early-stopping slice. This is the one fit that changes, because every CV number
        carried it and the CV numbers are the headline.
        ⚠️ **`fitc` rather than `trc` decided 12.08.2026 (Selin), on comparability.** "That fold's
        training lines" became ambiguous once audit 07 nested the early-stopping set: `trc` is the whole
        training side, `fitc` is `trc` minus the stopping slice, and every other per-fold statistic —
        `fit_weight_fns`, the head-bias init — already uses `fitc`. The deciding argument is **what the
        two arms see**, not the leak. PCA and scGPT already receive the same cells, the same per-cell
        CPM and log transform, and the same all-cells HVG gene set (decision 1). The *only* place they
        differ is that **PCA needs a fit — mean, std, rotation — and scGPT needs none**: its weights are
        pretrained elsewhere and frozen, and its value binning digitizes each cell against its own
        distribution, so nothing is estimated across cells. Since that fit is the single asymmetry,
        keeping it as narrow as possible is what makes a measured difference attributable to the
        *representation* rather than to how much the control's fit was allowed to see. `fitc` also
        removes the exception: PCA is now fitted like every other per-fold statistic.
        **Cost, stated rather than glossed:** ~15 % fewer cells for a 512-component fit, so if scGPT
        wins, PCA can be asked whether it was undersold. The answer is that the fit still uses ~85 % of
        the training cells, and the alternative hands the control an advantage scGPT structurally cannot
        receive. Rejected: `trc`, better estimated but it widens the one place the arms are not alike.
        `X_pca` is still written and remains correct for UMAPs and
        other descriptive use; CV stops reading it. Implementation is in the training path
        (`model/dataset.py::resolve_rep` today falls through to `X_pca` whenever `split_col is None`),
        with **one helper called by both `cv.oof_predictions` and `train_multitask.cv_evaluate`** — never
        two, since the §A/§B asymmetry has already produced a real defect twice.
        ⚠️ **Re-confirmed 12.08.2026 after the cost was corrected.** This was first recorded as needing
        no preprocessing change, which was wrong: `4a` opens the targets h5ad *backed* and copies only
        `obs`/`obsm`, so the CV path holds no expression matrix, and the matrix it needs is not that file
        — `add_pca` fits on `counts_h5ad` (`paths.raw_h5ad`, **2.15 GB**), because the targets `.X` has
        the scGPT OOV genes dropped and would give a different gene set, so per-fold fits taken from it
        would not be comparable with the stored `X_pca`. **True cost:** the training path opens a file it
        has never touched, carries ~1 GB alongside the embeddings, and does five 512-component fits per
        CV run for the PCA arm — repeated across R4's grid (MSE/MAE/Huber × α ∈ {off, 0.5, 1.0}, plus
        item 8C), which is why the fits are **cached per run, keyed by the fold assignment** (Selin,
        12.08.2026). Folds are deterministic given the seed and the eligible line set, so the cache is
        sound and the projections need not be stored. The cache is **in-process, not on disk under
        `runs/`**: a cached projection outliving the code that produced it is the exact failure class
        the 03.08 freeze exists for, `runs/` is gitignored so it would be invisible in review, and it
        would fail by producing numbers that are plausible and wrong.
        ⚠️ **float32 — BOTH fits, decided 12.08.2026 (Selin), and this makes decision 2 gate R2.** The
        counts matrix is stored dense float64 (53,513 × 5,000 = 2.14 GB, the whole file), so casting on
        load halves the per-fold cost: ~1.07 GB resident, ~0.86 GB per fold subset, peak ~2.2 GB rather
        than ~4.4. Casting **only** the per-fold fits was rejected: `add_pca`'s descriptive all-cells fit
        runs float64, so the two would then differ in dtype as well as in which cells they see — undoing
        the property Selin established on 10.08.2026 when she harmonised `ddof` to 1 precisely so that
        *which cells are seen* would be the only difference between them. So `add_pca` casts too. That
        is a preprocessing change: it moves `X_pca` in its last digits, and **R2 must run after it**,
        not before. The precision cost for PCA on log-CPM is nil at the precision anything quotes.
        **Rejected: precomputing and storing five fold-keyed matrices at R2.** It would make the fold
        PCAs inspectable artifacts, consistent with `uns["pca_fits"]` (audit 04b) — but it adds ~548 MB
        per variant (512 comps × 53,513 cells × 5 folds), ~1.1 GB across `hvg5000` + `all_genes`, and it
        pushes CV configuration (`n_splits`, seed, `group_col`, `eligible_splits`) into preprocessing, so
        changing a *training* parameter would force a *preprocessing* re-run. The artifact argument is
        also weaker here than in 04b: that stored a **fit** — loadings and variance ratios, answering
        questions nothing else could — where this would store **projections**, recomputable from a
        deterministic fold assignment.
        **3 — the never-training cells stay in.** Decision 2 already removes them from the PCA the model
        reads, because a fold's training lines are eligible (labelled) by construction; what remained was
        only whether they should help choose the HVG set, and they should — no label exists to leak, more
        cells estimate gene variance better, and both arms share the result so neither is advantaged.
        ⚠️ **Consequence for the record:** the leak in `resolve_rep` is closed for CV rather than
        documented, so Step 02's "still leaky, and documented as such" and the matching docstring stop
        being true at R4 and must be rewritten there, not before — the code changes at R4, not now.
        Detail on what each fit sees:
        [Step 02](./steps/02-preprocessing-and-embeddings.md#what-transform-pca-sees--corrected-05082026).
        1. **HVG selection is all-cells**, for both arms — the one fit the two representations share.
        2. **The cross-validated PCA is all-cells.** The fixed splits were fixed 05.08.2026
           (`X_pca_train_ctrp`); CV folds are drawn at training time, so five fold-specific matrices
           cannot be stored and `resolve_rep` leaves them on the leaky `X_pca`. Every CV number carries it.
        3. **Cells that never train are in all three fits** — the 17 lines / 6,073 cells (18 / 6,286 on
           disk today) with no CTRPv2 label, 11.4 % of the atlas. Not a test leak; a separate question
           about whether the representation should be shaped by data the model never sees.
        Any change here alters the gene set and therefore every number, so it lands in the sweep.
- [ ] **8 · Model — walked 12.08.2026.** Confirmed sound: the architecture is what
      [Step 03](./steps/03-model-and-training-design.md#model-architecture--regularization-oncomlppy-25052026)
      says it is, the trunk is genuinely matched between the arms, the "K heads" are the K rows of one
      output `Linear` over a shared trunk with no per-drug sub-network, and inference runs under
      `.eval()` so dropout is off when predictions are made. **The shared trunk is 74,304 parameters at
      every panel size — only the head layer scales with K**, at 65 parameters per drug: 715 of 75,019
      at K=11 (1.0 %), 35,425 of 109,729 at K=545 (32.3 %). So the capacity the heads compete *for* is
      fixed no matter how many of them there are, which is what makes "capacity competition between
      heads" the right description of the K=545 collapse and a bigger panel no cheaper in trunk than a
      smaller one. Against ~153 independent labels.
      ⚠️ **Provenance repaired 13.08.2026 — the four counts are correct, the citation was dead.**
      `arch_facts.py` was **never committed anywhere** — not on `main`, not on any branch, not in any
      commit in this repository's history; it existed only in a session scratch directory. All four
      counts were re-derived against the real `scripts/model/OncoMLP.py` with `DEFAULT_HIDDEN_DIMS` and
      `input_dim` 512, and **every one matches** (trunk 74,304 at every K; 65 parameters per head;
      715 / 75,019 at K=11; 35,425 / 109,729 at K=545). Stated this way deliberately: after two
      dead-provenance findings on the same day the reflex is to distrust the numbers, and here the
      numbers are fine. The computation is pure architecture — no data, no training — so it is cheap to
      re-derive; **whether it earns a gated cell in `4a` beside the weight-decay one is open.** The
      weight-decay cell cleared that bar because it *decided* something; this one only *describes*
      something, which may be the distinction worth keeping.
      **The per-cell framing stands until MIL replaces it, and what MIL has to replace is the dataset and
      the loss, not `OncoMLP`** — the encoder carries over unchanged as the instance model. Q2's design,
      its controls and the open decision on what counts as a positive result stay in *Agreed plan,
      Step 2*; the objective side — that the loss penalizes the within-line variation — is item 9.
  - [x] **A · The two uncentred-target mechanics ran in one training path of three — FIXED IN CODE
        12.08.2026.** Head-bias initialization and the weight-decay exemption were applied only in
        `cv.oof_predictions`; `train_multitask.cv_evaluate` and `train_rep`, which produce the entire
        8-run matrix, had neither, so on a target centred near 0.9 the matrix trained against an offset
        the panel run did not. The flag was also not what the docs described: `exclude_output_from_decay`
        exempted the whole output `Linear`, weight matrix included, while still decaying the LayerNorm
        parameters. **Replaced with the standard grouping (Selin, 12.08.2026):**
        `TrainConfig.no_decay_bias_and_norm`, default on — every weight matrix decayed, every bias and
        normalization parameter exempt, as in HuggingFace `transformers`' `Trainer.create_optimizer` from
        the BERT reference implementation — plus `init_head_bias=True` on all three paths via
        `OncoMLP.init_head_bias_`, with the means taken per **line** (`cv.per_drug_line_mean`).
        Record: [Corrections](./steps/corrections-and-dead-ends.md#the-two-uncentred-target-mechanics-ran-in-one-training-path-of-three).
        **Takes effect at R4.** Note the target itself was left uncentred: per-drug mean-centring is the
        more standard fix and was rejected here because it is a target change and pre-empts half of S1.
    - [x] ✅ **CLOSED 13.08.2026 (Gate 1 sweep) — the fourth path was fixed, hours after this was
          written, and the box stayed red.** `dreval_benchmark.ipynb` cell 8 now calls
          `init_head_bias_(model, per_drug_line_mean(_y_lines, _obs_lines))` at line 56, landed in
          `8b6a678` — a commit whose subject is literally *"the head-bias fix reached three training
          paths of four -- this is the fourth"*. **The sub-item (i) below is NOT closed by this**; the
          missing `DataLoader` generator is a separate defect in the same cell and is still live.
          ⚠️ Worth keeping as a process finding rather than only a code one: this box was marked 🔴,
          dated the same day, and read *"NOT fixed — Selin's call"* while the fix was already in the
          working tree. A red flag that outlives its defect costs the same attention as a real one.
          What it said, in the past tense: cell 8's `run_oncomlp` **also** trained an `OncoMLP` and
          nobody had enumerated it — it constructed the model and called `train_model` with **no
          `init_head_bias_` anywhere**, while `cv.py:375` called
          `init_head_bias_(model, per_drug_line_mean(...))`. The consequence it named was real for the
          window it describes: `auc_cc` centres near 0.9, a zero-initialized head starts every drug an
          offset from the base rate, and early epochs are spent travelling there — so between 12.08 and
          the fix the benchmark was biased *against* our model, the same direction as the 25-vs-50 epoch
          defect, in the same cell, found the same way. It was not fixed unilaterally at the time because
          it changes what the benchmark measures, which is the same class of decision as the epoch count.
          **Both were then decided and applied together.** No number is affected: the benchmark has not
          been re-run, and its existing numbers were already void on three other grounds.
          Found by enumerating cell 8 against `cv.oof_predictions` rather than spot-checking — spot-checking
          is what missed it the first time.
      - Two lesser differences from the same enumeration, recorded so they are not re-found.
        **(i) No `DataLoader` generator** — `dreval`'s train loader is `shuffle=True` with no `generator`,
        which is precisely the ordering dependency `cv.py:378-386` carries a comment about having removed.
        **Benign today**, because `train_model` calls `set_seed(config.seed)` before the loader is first
        iterated, so the shuffle is in fact seeded — fragile rather than broken, and fragile in the exact
        way the project already decided against. **(ii) Validation `batch_size` 256 vs the pipeline's 128**
        — numerically nothing: dropout is off at eval and `LayerNorm` is per-sample.
        Everything else matches: dropout 0.5, input dropout 0.1, norm layer, train batch 128, lr 1e-3,
        weight decay 0.0, loss `mse`, `no_decay_bias_and_norm` on, epochs 50, seed 42.
    - [ ] **Open, not decided: the bias init does not start the model at the null predictor.** It fixes
          the level only approximately — the randomly initialized head weight rows already scatter
          predictions with sd ≈0.31 at initialization, against a true across-line spread of order 0.17,
          and the mean lands at 0.76 / 0.96 for a requested 0.90 at seeds 42 / 0 (`init_spread.py`,
          synthetic). Starting genuinely at the null needs the output layer's *weights* initialized small
          as well — Lin et al. use σ = 0.01 alongside the bias prior. **Selin's call**, and it is a
          second architecture change, so it should not ride along with A.
  - [ ] **B · `input_dropout=0.1` is matched in value but not in effect between the arms.** It zeroes
        each input coordinate independently. scGPT's 512 dimensions are entangled and comparable in
        magnitude, so the perturbation is uniformly small; PCA's are variance-ordered, so one draw in ten
        deletes PC1. Same expected variance removed, far heavier-tailed for PCA — and the "matched trunk
        ⇒ fair comparison" argument covers the trunk, not the input regularizer. Settleable from
        `uns["pca_fits"]["variance_ratio"]` once R2 writes it (item 4B), before any decision to change it.
  - [ ] **C · The evidence that closed model-side tuning is void — decide the minimal re-derivation.**
        `ablations_and_rescue.ipynb` early-stopped on the fold it scored, on retired `auc_z`, over the
        five drugs of the discredited gate, and cannot be re-run.
        [Corrections](./steps/corrections-and-dead-ends.md#the-evidence-that-closed-model-side-tuning).
        Only the MLP rows are flattered — ridge has no early stopping — and ridge *ties* the PCA MLP, so
        the honest ordering may be ridge above it. **Scope agreed (Selin, 12.08.2026): the minimal
        re-run** — trunk `(128,64)` vs a bare linear head, both representations, against `RidgeCV` on the
        same folds and the rebuilt panel; not the four-knob sweep. Scheduled at R4.
  - [ ] **Found on the way, routed elsewhere — NARROWED 13.08.2026 (Gate 1 sweep).** It read
        *"`--epochs` defaults to 50 in the CLI, 25 in `TrainConfig` **and `4a_percell_training`**, and
        `4a_percell_training` §B sets 50"*. The notebook half is no longer true: `4a` sets
        `TrainConfig(epochs=50, seed=SEED)` once, for **both** sections. ⚠️ **What remains is
        sharper than the original, not softer:** `TrainConfig.epochs` was still **25** while **every
        caller in the project passed 50** — the CLI default, `4a` §A and §B, `dreval_benchmark`
        cell 8, and `4b_mil_training`. A default that no caller uses is not a default, it is a trap
        for the next one written. ✅ **CLOSED 13.08.2026 (Selin): the default is now 50**, which
        also settles item 10's epoch question. The 36-run evidence behind 25 is kept in the field's
        own comment — it argued that 25 was *enough*, never that 50 was wrong, and it predates
        `cv.inner_holdout`, so it no longer describes where training peaks. Early stopping
        (`patience=10`) is what ends training in every recorded run. `dreval_benchmark.ipynb` builds
        `OncoMLP` by hand and so has neither mechanic from A → **item 11**, which owns that notebook.
        `ScGPTDrugDataset` has no consumers — the `train_baseline.py` / `train_scGPT.py` its docstring
        names were deleted in `090f957` — and the `norm="batch"` / `"none"` branches are never exercised
        → **item 13**. LayerNorm makes the forward pass invariant to input rescaling up to the first
        `Linear`'s bias (78× moves the output by 0.074 against a spread of 0.377; zeroing that bias drops
        it to 8e-6), so **4A's premise needs restating before 4A is run** — the naive "different
        effective step size" argument does not go through under LayerNorm plus Adam.
- [x] **9 · Loss — walked 12.08.2026.** Confirmed sound: the masked loss is `Σ(sq·M)/ΣM` over observed
      (cell, drug) entries, so unscreened pairs contribute nothing to loss, gradient or metric; sample
      weights ride in the mask, turning it into `Σ(w·sq)/Σw` exactly with no change to the training
      loop; the density is fitted per fold on the fitting lines only. **The objective stays plain masked
      MSE until MIL (Selin):** no spread term, no ranking term, no tuned weights — the reasoning, and
      what replaces each, is in
      [Step 03](./steps/03-model-and-training-design.md#the-loss-is-plain-masked-mse-and-stays-that-way-until-mil-audit-09-12082026).
  - [x] **A · The density weighting is re-tested, not retired, and `alpha` is swept instead of
        justified (Selin, 12.08.2026).** Its parameters had lost their source when
        `panel_distributions.ipynb` was archived, and the *explanation* for its null turned on a
        winsorization retired 11.08.2026 and a target since replaced. Rather than derive a number in
        advance to test a hypothesis the same run tests anyway, `alpha` becomes an arm of the loss
        comparison over **{off, 0.5, 1.0}**; `cap=3` is held fixed and **documented as arbitrary**.
        [Corrections](./steps/corrections-and-dead-ends.md#inverse-density-loss-weighting-improves-ranking).
  - [x] **B · The 82× per-line imbalance stands until MIL.** Line-balanced reweighting was already
        tested and is empty, and MIL removes the defect structurally — one bag is one line is one
        example — so a weighting built now is machinery built to be discarded.
  - [x] **D · The optimizer is AdamW, and the decay that rode on it is 0.0 (Selin, 12.08.2026).**
        Loshchilov & Hutter, *Decoupled Weight Decay Regularization*, ICLR 2019: `weight_decay`
        passed to Adam is an L2 term added to the gradient, which Adam's adaptive scaling divides
        through, so decay strength is entangled with gradient history; AdamW applies it to the
        weights directly. That argument is **unaddressed** — we are not using the standard form.
        **AdamW was implemented and reverted the same day (Selin, 12.08.2026)**, because measuring
        it first showed the migration carries an uncosted parameter decision (8 epochs, real
        `OncoMLP` and grouping, synthetic data at workload scale):
        **The numbers are in `4a_percell_training`'s `RUN_WD_CHECK` cell, not here** — run it and the
        comparison reproduces. **AdamW at the same nominal value is indistinguishable from no weight
        decay**, where Adam at 1e-3 measurably shrinks the weight matrices. AdamW's step is
        `θ -= lr·wd·θ`, a factor of `1 − 1e-6` at these settings; Adam's L2 enters the gradient and is
        amplified by `1/sqrt(v)`.
        ⛔ **Corrected 13.08.2026. Three figures stood here — 5.597 / 8.971 / 8.975, "0.04 %", "38 %"
        and "`wd ≈ 1.0`, three orders of magnitude up" — and none of them is re-derivable.** The script
        that produced them was never saved; what survives in scratch is a different measurement. When
        the comparison was rebuilt as a committed cell it gave different absolute norms, because the
        settings behind the originals are unrecoverable. **The conclusion is unchanged and slightly
        stronger** — the two AdamW arms differ by less on the reproducible measurement than on the lost
        one — but the numbers were the justification for a merged configuration choice, and they could
        not be checked. They are removed rather than restated: the cell is the source, and a number
        stated in two places is one that can drift. `wd ≈ 1.0` goes with them; nothing now supports it,
        and it was only ever the rejected alternative.
        *This is the failure item 12's closure predicted in its own words — a check written to confirm
        work, passing on first run, never entering the repo, and missing when it turned out to matter.*
        **So `weight_decay = 0.0`, chosen rather than defaulted.** Carrying 1e-3 across would have
        stated a setting the optimizer ignores. The rejected alternative was `wd ≈ 1.0`, which
        reverse-engineers a value to reproduce the shrinkage of runs that are themselves void.
        Aligned in three places so the CLI, the dataclass and the sweep notebook cannot train
        different models — the defect this item found with `--epochs`.
    - [ ] **What opens in its place: does this model need weight decay at all?** Not settled.
          It must **not** be justified by ["the model is over-regularized"](./steps/corrections-and-dead-ends.md#the-model-is-over-regularized-or-too-small)
          being refuted — that refutation rested on the weight-decay axis of an ablation which is
          itself void, so the evidence is gone. The honest justification for 0.0 is narrower: no
          sourced value exists for either optimizer, dropout 0.5 and input dropout 0.1 are already
          substantial, and a decay nobody can justify is worse than none. If regularization is ever
          *claimed* as tuned, it has to be derived first.
  - [x] **C · Two defects found and fixed in code.** In weighted runs the per-drug log printed `Σw` as
        `n=`, i.e. a weight sum labelled as a sample size — now printed as `w=` when it is not an
        integer count. And `TrainConfig.huber_beta = 0.05` is unsourced and **mis-scaled for `auc_cc`**:
        `smooth_l1` is linear above `beta`, and 0.05 sits well below the typical residual (~0.163 RMSE),
        putting roughly three quarters of residuals in the linear region — so `--loss huber` behaves
        closer to L1 than to what the docs describe. Left at its value rather than silently rescaled,
        because `beta` is a threshold on the residual scale and choosing one is an analysis decision;
        it is derived in the loss comparison if Huber is included.
  - [ ] **Left for item 11:** the four quantities the loss is *not* asked to optimize — order
        (Spearman), order at the top (NDCG@K against a random null), values (RMSE in AUC units) and
        spread (calibration slope) — get a section each in `notebooks/5_evaluation.ipynb`, which is
        written but not yet built. Item 11 owns the metric set and records it.
    - ⚠️ **The four sections exist; the table that puts them side by side does not (13.08.2026).**
      §1.4–1.7 each compute their quantity and §1.3 defines `decide()`, which applies the rule — but
      nothing assembles the four into one row per arm, and **`decide()` has no call site in the
      notebook at all**: its `def` is the only occurrence of the name. So *"`5_evaluation` is done"*
      is true of the quantities and false of the notebook. Two requirements are already recorded in
      §1.6 and have to survive into the table: each quantity's **`n_drugs`** appears beside it,
      because they are not always over the same panel (on a synthetic arm `order` averaged over 11
      drugs and `top_of_order` over 10, the thin drug having fallen below `TOP_K`); and the
      **different-drug-set check between arms**, which §1.6 flags but leaves to §1.3, where arms
      meet. **Where this sits in R1–R6 is open** — the table needs no rerun output to be *written*,
      but it cannot be *exercised* until R4 supplies a file with more than one arm in it.
      Reported by the preprocessing session on handover, 13.08.2026; the no-call-site claim verified
      against `notebooks/5_evaluation.ipynb` the same day.
  - [x] **How `order` is computed — decided by Selin, 13.08.2026.** Four sub-choices, all settled
        before any of the six arms is run, because the last comparison failed by choosing after seeing
        the numbers. **Spearman**, **per drug then averaged**, **unweighted mean across the 11 drugs**,
        **pooled across the five folds**.
    - **Spearman, not Pearson** — so that `order` and the calibration slope stay independent. Under
      Pearson a single defect (predictions compressed toward the mean) moves the primary *and* its
      spread guard together, and the guard stops being a check. Source: the primary-plus-guards rule
      itself (Selin, 12.08.2026, this file); comparability with `4a`'s own markdown, which calls
      within-drug Spearman "the quantity the whole project is judged on".
    - **Per drug, not pooled over all (line × drug) pairs** — pooling is the potency artifact the
      DrEval work is written about. Source in this repository, not the paper:
      `scripts/evaluation/dreval_normalize.py`'s docstring records a synthetic predictor with **zero**
      drug-specific signal scoring normalized Spearman **0.98** under this split design. Arm-vs-arm
      would survive pooling; the headline number would not be defensible on its own.
    - **Unweighted mean across drugs** — every drug counts equally. Chosen for continuity: it is what
      `4a` already computes and what every recorded number used. ⚠️ **Stated cost, not hidden:** panel
      coverage is *not* equal, so a thinly covered drug moves the headline as much as a well covered
      one. The count-weighted mean is the more defensible statistic and was rejected only because no
      earlier number would be comparable to it.
    - **Pooled across folds** — every held-out line appears exactly once across the five folds, which
      is what the out-of-fold table is for. Note that `analysis/qc/diagnostics.ipynb`'s
      `sd_across_folds` is the *other* combination (per fold, then averaged) and gives a different
      number; **that spread is not the seed band and must not be used as one.**
    - ⚠️ **Why the aggregation mattered more than it looks, and why it was fixed in advance.**
      `SEED_BAND` is now measured from the run's own ≥3 seeds, and the band is the seed-to-seed spread
      *of this aggregate*. A median across drugs is more stable than a mean, so it would have produced
      a **narrower band — an easier bar for an arm to clear.** The aggregation choice therefore sets
      the sensitivity of the whole loss comparison before any data exists. Surfaced by the item-11
      session rather than defaulted, which is the reason it is on the record at all.
  - [x] **Guard margins are named parameters, defaulting to each guard's own measured band —
        13.08.2026.** ±0.04 lives on Spearman's scale and does **not** transfer: `values` is in the
        target's own units and the calibration slope is centred on 1.0. This follows from the
        `SEED_BAND` decision (measure it from the run) rather than being a new choice; it is recorded
        because the alternative — one margin reused across three scales — is the kind of thing that
        reads as principled once it is in a table.
- [x] **10 · Training — walked 12.08.2026. Nothing here gates R2.** Optimizer, weight-decay groups,
      epochs, early stopping and the `mps` nondeterminism, read against the code.
  - [x] **The `mps` nondeterminism does not reproduce under current code.** Measured on the real
        classes and the real `train_model` loop, synthetic input, nothing read or written: at the real
        workload's scale — 34,000 cells, 512-d input, 11 heads, 25 epochs — three identical runs are
        **bit-identical**, max |Δ val MSE| 0.000e+00, same best epoch, same final loss. Identical on CPU,
        and identical with and without the explicit DataLoader generator, because `train_model` calls
        `set_seed` before the loader is first iterated so `RandomSampler` draws from a freshly seeded
        global RNG either way. Run at two scales deliberately — concluding "deterministic" from an
        underpowered test is the silent-zero error. **So `06_limitations`'s 0.313 / 0.315 / 0.317 /
        0.320 is not re-derivable**, and the honest replacement is stronger than a caveat: at a fixed
        seed the training loop is bit-reproducible on `mps`. ⚠️ **Scope:** this exercises the *training
        loop*, not the full CV chain on real data — five folds, the line-level aggregation, the
        Spearman. Those are deterministic numpy given deterministic input, but they have not been run,
        so the claim is the loop's, not the chain's. Certifying the chain needs real data and waits for
        R2. **This also discharges audit 12's outstanding clause** — the 28.07.2026 no-action decision
        required measuring the non-determinism where it matters, and it is now measured.
  - [x] **Adam → AdamW — DONE, and this box was stale (closed 13.08.2026, Gate 1 sweep).**
        `training_utils.py` builds `optim.AdamW` on both branches; the migration landed 12.08.2026
        together with `weight_decay = 0.0`. ⚠️ **Its first sentence was false about the code from
        that day onward** and is preserved here in the past tense, because it is what a reader
        would have acted on: *"Weight decay **is currently** added as L2 to the gradient under
        `optim.Adam` rather than decoupled."* The item's own closing parenthetical below already
        said the finding "closes here" — the box and its opening claim were simply never updated,
        which is how a done item keeps asking to be done. The argument for the switch stands: Loshchilov & Hutter, *Decoupled Weight Decay
        Regularization*, ICLR 2019, showed the two are not equivalent under Adam and that the decoupled
        form is the one that behaves as intended. It interacts directly with the grouping audit 08
        changed — *which* parameters are decayed and *how* decay is applied are one decision — so the
        grouping is unchanged and only the application moves. **R4-side; lands with the training branch,
        not before R2.**
        ⚠️ **And `weight_decay = 0.0`, written explicitly (Selin, 12.08.2026).** Not a default and not
        an oversight: **at `1e-3` under AdamW there is effectively no weight decay at all**, while Adam
        at the same nominal value measurably shrinks the weight matrices. **The measurement is the
        `RUN_WD_CHECK` cell in `4a_percell_training`** — real `OncoMLP`, real parameter grouping,
        synthetic input, settings fixed in its docstring — and the numbers live there so they can be
        re-run rather than believed. The arithmetic: AdamW's step is `θ -= lr·wd·θ`, a factor of
        `1 − 1e-6` per step at `lr=1e-3`, where Adam's L2 enters the gradient and is amplified by
        `1/sqrt(v)`.
        ⛔ **The figures that stood here — 8.971 / 8.975 / 5.597, "0.04 %", "~38 %" — were removed
        13.08.2026 as not re-derivable**; see the corrected block under item 10. The conclusion is
        unchanged; the numbers had no reproducible source, and this record now points at one.
        So the choice was between carrying a number that does nothing, reverse-engineering a decay
        value to reproduce the shrinkage of runs that are **themselves void**, or saying plainly that this model
        trains **without weight decay** — dropout 0.5 and the input dropout as the only regularizers.
        The third was taken. **The justification is narrow and deliberately so:** no sourced value exists
        for either optimizer, dropout is already substantial, and a decay nobody can justify is worse
        than none. ⛔ It is **not** justified by "the model is not regularization-limited" — that
        hypothesis's refutation rested on the weight-decay axis of the
        [void ablation](./steps/corrections-and-dead-ends.md#the-evidence-that-closed-model-side-tuning),
        so **whether this model needs weight decay is now an open question, not a settled one.**
        ⚠️ **Consequence for R4:** its model differs from every recorded run in its regularization —
        those trained under `Adam` at `weight_decay=1e-3`, which measurably shrinks the weight matrices;
        R4 has none. One more reason R4 establishes a baseline rather than measuring a delta.
        *(Item 10's optimizer finding closes here — decoupled decay is what Loshchilov & Hutter argue
        for and it is now in use. What replaces it is the open regularization question above, not
        nothing.)*
  - [x] **DONE — `dreval_benchmark` now sets `epochs=50` explicitly** (closed 13.08.2026, Gate 1
        sweep; the change itself is dated 12.08.2026 in the cell). Cell 8 builds
        `TrainConfig(epochs=50, seed=42, log_every=1000)`. What the box said, kept because the
        reasoning is the reason the fix was right: it was the only caller that passed no `epochs`
        and so fell to `TrainConfig`'s default, which every other caller overrides. **This is not a change to the
        benchmark — it is the benchmark silently running half the training the pipeline uses**, so the
        notebook answering "how strong is this by the field's standard" has been scoring a weaker model
        than the field would see. Nobody chose 25 for it. R5-side.
  - [x] **Weight-decay groups — nothing left beyond AdamW above.** Audit 08 replaced
        `exclude_output_from_decay` with the conventional grouping and verified it by fault injection.
  - [x] **Early stopping — nothing left that gates anything.** Audit 07 closed the leak with
        `inner_holdout`. Two observations, both R4-side and neither blocking: the improvement threshold
        is an absolute `1e-6` on val MSE, which against a typical `auc_cc` MSE around 0.026 is a
        0.004 % relative bar, so nearly any change registers as improvement and `patience` runs longer
        than it reads; and at 50 epochs rather than 25, `patience=10` can now actually fire, where
        before the cap was reached first.
- [ ] **11 · Evaluation** — metric, cell→line aggregation, baselines (ridge, `NaiveMeanEffects`),
      dispersion across folds *and* drugs, raw vs normalized.
  - [ ] **⚠️ A blind spot in the metric set, not in any one result (found 12.08.2026, audit 09).**
        Every intervention this project has judged was judged on **rank correlation and MSE**. Neither
        can see a change in *calibration*: widening the predictions without reordering them leaves
        Spearman identical, and MSE is minimised by the shrunken predictor to begin with. So for any
        intervention whose effect is on spread rather than order, "flat on both" never distinguished
        *no effect* from *an effect the metrics cannot see*. This is a property of the instrument, and
        it silently qualifies a class of conclusions rather than one of them:
        [inverse-density weighting](./steps/corrections-and-dead-ends.md#inverse-density-loss-weighting-improves-ranking)
        (re-tested at R4, since the mechanism demonstrably fired — predicted spread moved and the
        ranking did not), [line-balanced reweighting](./steps/corrections-and-dead-ends.md#line-balanced-reweighting-will-help)
        (not re-tested — MIL dissolves it either way), and the regularization axis of the
        [void ablation table](./steps/corrections-and-dead-ends.md#the-evidence-that-closed-model-side-tuning),
        whose own text attributes the harm of heavy regularization to *over-shrinkage* — an effect its
        metrics could not have measured. The claim that survives untouched is that shrinkage is
        MSE-optimal at low ρ, because that is an argument rather than a measurement.
        **Fix: the calibration slope has to exist before any of these is re-read**, which makes
        building `5_evaluation` a prerequisite for R4's loss comparison rather than a follow-up to it.
        The report passage that drew the strongest version of this inference is corrected in
        `report/sections/06_limitations_and_outlook.tex`.
  - [ ] **⚠️ Where the ridge baseline lives — an open choice, raised 12.08.2026 (Selin).** `RidgeCV`
        on cell-line mean embeddings is fitted *inside* `4a_percell_training` §A, which calls it "the
        baseline that actually binds". Since `4b_mil_training` is to be scored by the same
        `5_evaluation` over the shared line-level out-of-fold format both notebooks emit, that
        placement decides whether the two architectures are measured against **one** ridge or against
        two separate fits of it. The choice:
        **(a) Move it into `5_evaluation`**, beside `NaiveMeanEffects`, computed once from the shared
        line-level predictions. The comparison then holds by construction rather than by convention —
        the same argument `4b` §2 makes for one scorer — and the baseline cannot drift between arms.
        Cost: `4a` §A stops being a self-contained result, and the ridge fit needs the fold assignment
        carried in the shared format rather than taken from the training notebook's own
        `cv.grouped_folds` call.
        **(b) Leave it in `4a` §A** and fit a second one in `4b`. Each notebook stays readable end to
        end, at the cost of two fits that must be held identical by hand — the drift the shared-scorer
        decision was taken to remove.
        **What (a) rests on:** that `4a` and `4b` genuinely emit one common line-level format. If they
        diverge, (a) is unavailable and the question is moot — so this is decidable only once `4b` is
        past a stub.
        **A reading, not a decision:** (a). This matters more than a tidying choice because ridge has
        no early stopping and so was never flattered by the leak that flattered the MLP, and on PCA the
        two were equal — the honest ordering may put ridge ahead (recorded with item 8C). A baseline
        that may outrank the model is not a throwaway control and should not be computed twice.
        Sequenced **after item 8C**, which re-derives the ridge comparison at R4 on the rebuilt panel,
        and **after `5_evaluation` exists**, which item 11 already blocks.
- [x] **12 · Reproducibility — walked 12.08.2026; CLOSED 13.08.2026. Nothing here gates R2.** Seeds
      swept across 22 script modules and 14 live notebooks for every stochastic call site, and
      **independently re-verified**: every `random_state` in `scripts/` takes `seed`; both
      `DataLoader(shuffle=True)` sites pass an explicit `torch.Generator().manual_seed(config.seed)`;
      `GroupKFold` does not shuffle; `GroupShuffleSplit` passes `random_state=seed`; every scanpy
      `neighbors`/`umap` call in the notebooks passes `SEED`; and `dreval_benchmark`'s `split_dataset`
      passes `random_state=42`. **Nothing is unseeded.** Two candidates were flagged and both are false
      positives on inspection — a docstring containing the string "PCA(512)", and `RidgeCV`, which takes
      no `random_state` because its solver is closed-form.
      ⚠️ **The sweep's denominator is not recorded and could not be re-derived** (13.08.2026). The
      original pass reported "19 sites checked"; a second session re-verified every claim *about* those
      sites but could not reproduce the count, because what counted as a site was never written down.
      The list above is what is independently checkable and is stated instead. Same class as the
      weight-decay figures removed under item 10 — a number quoted from a pass whose definition did not
      survive it. Determinism: the 28.07.2026 no-action decision holds — `set_seed`
      sets random/numpy/torch/cuda/mps and forces no flags, exactly as specified — and its second clause,
      *measure the non-determinism where it matters*, was discharged by audit 10, which found the
      training loop bit-reproducible on `mps` at a fixed seed. `frozen_split` **raises** rather than
      guessing when the data holds an eligible line the frozen file does not cover, so the 180 → 181
      recovery fails loudly instead of silently reassigning.
      ⚠️ **Two things stay open and are tracked at R6, not here.** (1) The hand-transcription below —
      16 macros, all cited, but `results_numbers.tex`'s dataset block points at `notebooks/04`, removed
      by the restructure, and that dead pointer is the provenance for **13 of the file's 16 citations**.
      (2) `\NLinesCV` and `\NLinesFit` are **derived**, not transcribed — hand-arithmetic at 85 % and
      70 % of `\NLines`, with **no artifact to extract from**, so an extraction script pointed at CSVs
      cannot produce them. That distinction has to be settled before item 12's fix is scoped.
      ⚠️ **And the item's own standard is currently failed by our tooling**: `scripts/` holds only
      `layout.py` and `check_resolved_paths.py`, so every other verification number quoted during the
      audits — links, artifact references, command paths, notebook validity, module imports, report page
      counts — came from scripts in agents' scratch directories and **cannot be re-derived by anyone**.
      The asymmetry that explains it: `check_resolved_paths.py` was committed because it caught a defect
      that had already bitten; checks written to *confirm* work passed on first run and never entered the
      repo. That predicts which check is lost next — whichever currently passes.
      *Seeds are now fixed
      (`gen_embeds.py`, `add_pca.py`, the training `DataLoader`s) and determinism was decided
      28.07.2026 as no-action, so what is left is the hand-copying:* **`report/results_numbers.tex` is
      transcribed by hand from CSVs under `notebooks/outputs/` with no extraction script**, so every
      macro can drift silently from its source — and it grew again on 10.08.2026 with four replicate
      macros (`\NRepPairs`, `\NRepLines`, `\RepDiffMed`, `\RepDiffPct`). Related and already fixed:
      `report/main.pdf` was tracked and went stale whenever a `.tex` changed; it is now untracked and
      built from source.
- [x] **13 · Code redundancy, stale code, notebook restructuring — CLOSED 13.08.2026.** The **code**
      half of Selin's "redundanz, staleness, file overload" note; the documentation half was done
      05.08.2026. Done across 12.08: the notebooks restructured into the five-stage pipeline
      (`1_data` → `5_evaluation`, analysis regrouped under `analysis/{qc,harmonization,evaluation}`),
      `1_preprocessing` and the `run_preprocessing.py` CLI archived in favour of
      `scripts/preprocessing/pipeline.py`, dead outputs moved to `notebooks/outputs/archive/`, dead code
      archived (`ScGPTDrugDataset` and the single-drug chain), and `scripts/check_resolved_paths.py`
      added — a pre-merge check for paths built from variables, written after three such defects in one
      week, one of which would have silently recomputed a cross-validation while reporting that it had
      loaded committed folds.
      ⚠️ **Deliberately left, and Selin's to schedule rather than this item's to finish:** 22 hardcoded
      absolute `/Users` paths across 6 live notebooks — `drug_catalog` holds 14 of them and needs a
      rewrite rather than path edits, since it also reads CTRPv2's retired `v20.*` tables and inputs
      that are not in the repository; four cwd-dependent relative writes, two of which write **into**
      `outputs/archive/` from a live notebook; three notebooks pointing `OUT` at the flat `outputs/`
      root (latent — every current call site re-appends its subdirectory); and two orphaned report
      figures referenced by no `.tex`.

**Working agreement for the session, restated because it was broken today:** no step of the analysis gets
decided silently. If a choice affects what enters the model or how a number is computed, it is proposed
first and agreed before it is executed.

## The sweep — R1 to R6, in order (written down 11.08.2026)

*This sequence existed only in the session task list until now, which meant it lived in a per-session
cache and would not have survived the session. Nothing here may start before the review above is
finished and Selin says so (03.08 banner); R1 is a decision, not a run.*

> ### ⭐ Every unarchived notebook is re-run, in a dependency-respecting order (Selin, 13.08.2026)
>
> **This governs R2–R5 and takes precedence over their per-item notebook lists wherever the two
> disagree.** The sequence was written as a set of *targeted* re-runs — R2's two preprocessing drivers,
> R3's two refreshers, R5's named analysis notebooks — and a targeted list silently exempts every
> notebook nobody thought to name, leaving it holding results computed before the code changed. The
> scope is now **every notebook under `notebooks/` that is not in `notebooks/archive/`**.
>
> `notebooks/archive/` stays out, by the same rule that put it there: those notebooks read targets that
> no longer exist, so re-running them either raises or manufactures numbers on a retired scale.
>
> **The order — the numbered chain in sequence, then everything else (Selin, 13.08.2026).** Two groups,
> run one after the other rather than interleaved.
>
> ⚠️ **PROVISIONAL — this was produced out of sequence and is not yet the plan of record.** It was
> written while Gate 1 ("close everything still open or wrong") was still open, and Gate 4 exists
> precisely to *re-decide the order* once Gates 1–3 have reported. That dependency is not theoretical:
> on the day this was written, a Gate 1 finding moved `diagnostics` and `dreval_benchmark` from "blocked,
> cannot run at all" to "waits for R4" (`35fe0bc`), which changed where both sit. Any further Gate 1
> finding can move it again. **Read the structure and the verified dependency counts as established;
> read the sequence as a proposal until Gate 4 confirms it.**
>
> **1 · The numbered chain, strictly sequential.** `1_data` → `2_drug_selection` → `3_representations`
> → `4a_percell_training` → `4b_mil_training` → `5_evaluation`. `4b` **stopped being a stub on
> 13.08.2026**: its criterion is implemented in §3 against `scripts/training/mil.py`, and it now runs
> after `4a` rather than being outside the chain, because stage 1 reads `4a`'s within-line spread
> table. `PCA_SEED` was settled at 42 on 13.08.2026; what remains is the ⬜ aggregation of the
> per-pair tests, and the magnitude at which stage 6's veto fires. Its two missing covariates were
> recovered from the raw `UMIcount_data.txt` on 13.08.2026 and are now written by `convert`.
> **Neither gates the rerun** — both are computed from
> the saved per-cell prediction arrays and can be settled afterwards without retraining. See the
> notebook's own §3.3, §3.5 and §3.8.
> **Verified 13.08.2026: no
> numbered notebook reads any analysis notebook's output — zero, across all six.** The chain is
> self-contained and nothing outside it can gate it.
>
> **2 · Every other notebook, after the chain has finished.** Each of the seven runnable ones reads a
> pipeline artifact, so running the group *after* the whole chain satisfies every prerequisite at once
> and removes any need to interleave:
>
> | notebook | what it reads | earliest the reads allow |
> |---|---|---|
> | `analysis/qc/hvg_sweep_build` | — it **writes** `hvg1000/2000/3000`, driving all six steps itself | after `1_data` |
> | `analysis/qc/gene_symbol_rescue` | `SCP542_CCLE.h5ad`, the scGPT OOV tables | after `3_representations` |
> | `analysis/qc/verify_variants` | five variants' targets, raw h5ad | after `hvg_sweep_build` |
> | `analysis/harmonization/drug_coverage` | `hvg5000` and `all_genes` targets | after `3_representations` |
> | `analysis/harmonization/cell_line_join_verification` | targets h5ad, DrEval CTRPv2 tables | after `3_representations` |
> | `analysis/evaluation/diagnostics` | `outputs/panel/panel_oof_predictions.csv` | after `4a` |
> | `analysis/evaluation/dreval_benchmark` | `outputs/panel/panel.csv`, `panel_oof_predictions.csv` | after `4a` |
>
> The last column is what the file reads permit, **not a schedule**. Running the whole group after the
> chain is what makes it independent of the chain — the property that matters operationally: these
> notebooks never gate the pipeline, never feed back into it, and need no fixed order among themselves.
> **One exception, and it is inside the group:** `hvg_sweep_build` must precede `verify_variants`,
> because it builds three of the five variants that notebook reads. It is a pipeline driver living in
> `analysis/qc/`, which is why it does not look like a dependency until you read its cells.
>
> ⚠️ **Not in the rerun: `analysis/harmonization/drug_catalog`.** 15 absolute `/Users/` paths and 6 reads
> of CTRPv2's retired `v20.*` tables. It is simultaneously the **only** notebook with no pipeline
> dependency at all and the only one that cannot run on another machine. It needs the rewrite already
> recorded above, not a position in the order.

- [x] **R1 · DECIDED 12.08.2026 (Selin): re-embed `hvg5000` + `all_genes`.** Not all five, not
      `hvg5000` alone. This covers every number the report currently quotes, at the middle cost —
      scGPT embedding is the expensive step, which is why the scope had to be set before R2.
      Five variants exist on disk (`hvg1000/2000/3000/5000`, `all_genes` — `layout.py`,
      `VARIANT_N_TOP_GENES`); the three that are not re-embedded keep their current artifacts.
      ⚠️ **What this costs, recorded so it is not rediscovered:** [Step 05](./steps/05-multitask-results.md)'s
      gene-set sweep spans `hvg1000/2000/3000/5000`, so after R2 it **mixes re-embedded `hvg5000` with
      three variants embedded by the older code** — before the gene-symbol repair, the seeding and the
      `ddof=1` harmonization. The sweep is therefore not like-for-like and any conclusion drawn across
      its points needs that stated, or the three remaining variants re-embedded later as a top-up.
      Rejected alternatives and why: **all five** keeps the sweep like-for-like but is the longest run;
      **`hvg5000` only** is cheapest but voids the sweep entirely and leaves every `all_genes` number in
      the report stale.
- [ ] **R2 · Re-run preprocessing end to end.** Driver: the notebooks —
      `notebooks/1_data.ipynb` (`fetch`, `convert`) then `notebooks/3_representations.ipynb`
      (`scgpt`, `targets`, `splits`, `pca`), both calling `scripts/preprocessing/pipeline.py`. Needs
      `overwrite=True` on the guarded steps and the separate scGPT venv as `SCGPT_PYTHON`.
      ⚠️ **`run_preprocessing.py` was archived 12.08.2026** and this item used to name it as the
      driver; there is no CLI replacement. Every artifact under
      `data/processed/scRNAseq_SCP542/<variant>/` predates the code that now produces it:
  - `scp542_conversion.py` annotates `var["hgnc_symbol"]` (`gene_symbols.py`, 05.08.2026) → `SCP542_CCLE.h5ad`
  - `gen_embeds.py` seeds with 42, runs on MPS, resolves through `resolve_gene_names` → embeddings and
    the OOV table; **4,576 → 4,704** genes (`hvg5000`), **20,570 → 21,332** (`all_genes`)
  - `add_pca.py` passes `random_state=42`, plus gene scaling, the post-HVG renormalization fix, the HVG
    ranking scale, and `uns["pca_fits"]` (item 4B)
  - `ctrp_to_h5ad.py` (audit 02): deduplicated experiment table + the `H292` alias → **180 → 181** lines,
    213 cells and 454 labels gained, 460 of `NCIH1299`'s targets changed
  - [ ] **Commit `splits/split_ctrp.csv` once this run writes it** (item 7B). The file has never existed,
        so `frozen_split` has redrawn on every run; `create_splits` writes it when it is missing, and
        committing it is what makes the guard real from then on. The 181st line means this draw shares
        nothing with the current assignment — expected, and the reason the file is committed *after* R2
        rather than before.
- [ ] **R3 · Refresh the committed read-only artifacts.** Held back only because they overwrite tracked
      files. `gene_symbol_rescue.ipynb` → `gene_symbol_rescue.csv` (773 vs 775 — no effect at quoted
      precision; on re-run it should import `load_rename_map` from `gene_symbols.py` rather than rebuild
      the map, and it reads only the vocabulary and HGNC so it could be released early).
      `verify_variants.ipynb` §7 `variants.png`, §8a `umap_cancertype_pca_vs_scgpt.png`, §8b
      `umap_sweep_cancertype.png`, §10b `scgpt_nonzero_per_cell.npz` — all need R2 first, and §10a–c are
      stale again despite running under the audit exception. **§8b and §10b cover all five variants**, so
      if R1 regenerates only some they will mix new and old embeddings — check before rendering. §9 trains,
      so it goes to R4.
> ### ⚠️ R4 carries six simultaneous changes, and that is deliberate — read this before attributing anything
>
> The standing rule is one change at a time, so an outcome can be attributed. R4 breaks it, carrying all
> of: the **cross-validated PCA fitted per fold** on each fold's fitting set, **float32** in both PCA
> fits, **MAE** added as a second loss arm, **Huber** removed, **AdamW** with **`weight_decay = 0.0`**,
> and `dreval_benchmark`'s **epoch fix**.
>
> ⚠️ **R4's model differs from every recorded run in its regularization.** Those runs used Adam at
> 1e-3, which measurably shrank the weight matrices; R4 has no weight decay at all. That is not a side
> effect of the optimizer switch — it is the decision taken with it (item 10D), and it is one more
> reason R4 is a baseline rather than a delta.
>
> **The defence is that there is no "before" to attribute against.** Every prior result is void on two
> independent grounds — the target was replaced on 11.08.2026 and the panel rebuilt on 12.08.2026 — so
> R4 is not measuring a delta from a previous run, it is *establishing the baseline* that later runs
> will be attributed against. Bundling changes is only dangerous when something is being compared to
> what came before, and here nothing is: the comparison R4 supports is internal to itself (MSE vs MAE,
> weighting off vs on, trunk vs linear head), and every arm of it carries all six changes equally.
>
> **What this costs, stated rather than hidden.** If R4's numbers disappoint, the six cannot be
> separated after the fact — no one will be able to say which contributed. That is the price of one
> clean baseline over six sequential runs, and it is acceptable only because none of the six is a
> hypothesis under test. Each is a defect fixed or a decision recorded, with its own justification:
> the per-fold PCA closes a leak, float32 preserves a property Selin established, MAE is an arm the
> plan already specified, Huber went because its position was set entirely by an unsourced parameter,
> AdamW is the standard form of a decay the code was already applying, and the epoch fix stops a
> benchmark training for half as long as the pipeline it represents.
>
> Written down here rather than left in commit bodies because a reader finding six bundled changes with
> no stated reason reaches for carelessness first, and would be right to.

- [ ] **R4 · Retrain.** The `DataLoader`s now take an explicit generator and the inputs change under R2,
      so no run under `runs/` is reproducible from current code. Requires the target (item 5) and the
      panel (item 6) settled first, and **never both in one run**. Scope: the 8-run matrix + 5-fold CV
      (`4a_percell_training.ipynb` (§B), `train_multitask.py`) — expected to **overturn** the Steps 04–05 numbers, not
      refresh them; `4a_percell_training.ipynb` on the rebuilt panel; ridge / `NaiveMeanEffects` on the same
      folds; `verify_variants.ipynb` §9; and 4A below. Blockers already recorded: **≥ 3 seeds** before any
      scGPT − PCA margin is quoted, and train-only drug selection inside each fold.
  - [ ] **The loss comparison** (item 9A, 12.08.2026): **MSE / MAE** × density weighting
        `alpha` ∈ {off, 0.5, 1.0} — **six arms**, on the per-cell architecture only; ranking losses wait
        for MIL, where one bag is one cell line and they are well-posed. Two conditions, both from the
        failure of the last comparison: **the decision rule is fixed before the run**, and **≥3 seeds**.
        Scored on all four quantities in `5_evaluation`, which is what lets a spread effect be seen at
        all — the previous null was measured on Spearman and MSE, both blind to it.
        ✅ **The margin is MEASURED from the run's own seeds, not inherited (Selin, 13.08.2026).**
        `5_evaluation`'s `SEED_BAND` stays unset and §1.3 **raises** until R4 supplies it. What is
        pre-registered is the *rule* — an arm must beat the others by more than the seed band — not the
        number; the bar is computed from the ≥3 seeds of the run being judged, which is the same
        construction that produced ±0.04 originally, on data that is not void.
        ⛔ **Why ±0.04 was not simply adopted.** It is the sample **sd of the gap between two arms**,
        from seeds 42/1/7 on the learnability-filtered subset — a **retracted** selection, on the
        **retired `auc`** target, **before** the early-stopping leak was fixed — and it was being
        applied as "the seed band on Spearman", which is a different quantity. Whether it was an sd, a
        range or a half-range was never recorded either. Three defects and an undefined definition, so
        pre-registering it would have fixed a bar that could not be justified.
        ⚠️ **Open sub-question, and it must be settled before the rule is applied:** ±0.04 is on
        Spearman's scale and does not transfer to the guards — values is in the target's own units and
        the calibration slope is centred on 1.0. Either each guard gets its own measured band, or the
        guards need separately stated bars.
        ✅ **All six arms compete for one winner (Selin, 13.08.2026)** — not loss-fixed-sweep-alpha, and
        not best-alpha-per-loss-then-compare. Recorded with the cost accepted rather than glossed: if
        the winner differs from the runners-up on **both** axes, the result cannot be attributed to
        either, which is the "one change at a time" constraint appearing inside a single comparison.
        Selin's note: analyse the attribution in more detail after the run rather than constraining the
        grid before it.
        ⚠️ **If Huber ever returns, `beta` returns unsourced with it.** Audit 09 left
        `TrainConfig.huber_beta = 0.05` at its value deliberately, to be derived in the comparison *if*
        Huber were included. Dropping Huber made that defect stop applying rather than fixing it — a
        closed-by-accident, recorded here so it does not have to be rediscovered. `0.05` against
        ~0.163 RMSE puts roughly three quarters of residuals in the linear region.
        ⚠️ **HUBER DROPPED 12.08.2026 (Selin); this read MSE / MAE / Huber until then, and item 9C —
        derive Huber's `beta` from the residual scale — goes with it.** Two grounds. Huber's role in the
        grid was to be the point *between* L2 and L1, and its position is set entirely by `beta`: at
        `TrainConfig`'s 0.05 against ~0.163 RMSE roughly three quarters of residuals fall in the linear
        region, so it behaves close to MAE and the grid would carry two near-duplicate columns. Fixing
        that means choosing `beta`, and every non-arbitrary choice (e.g. the textbook 1.345·σ for 95 %
        asymptotic efficiency) introduces a new sourced-but-imported constant into a comparison whose
        whole purpose is to *be* the justification — the same objection that retired `DEFAULT_WINSOR`
        and that keeps `cap=3` documented as arbitrary. **MSE and MAE already bracket the robustness
        axis**, which is what the comparison is testing; what is given up is only the ability to say a
        middle ground was tried. `MAE` is implemented for this (it has no parameter, which is part of
        why it survives the same objection Huber does not).
  - [ ] **The minimal capacity re-derivation** (item 8C, 12.08.2026): trunk `(128,64)` vs a bare linear
        head (`hidden_dims=()`), both representations, against `RidgeCV` on the *same* folds via
        `cv.grouped_folds`, on the rebuilt panel. ~20 fits. It re-establishes the one claim that carries
        weight — *scGPT needs the nonlinearity and PCA does not* — whose current evidence is void. Do not
        widen it into the four-knob sweep; that question is not what the load-bearing claim rests on.
  - [ ] **Two settings change here that no run on record used** (item 8A): every training path now
        initializes head biases at the fitting-fold per-drug means, and weight decay skips the biases and
        the LayerNorm parameters while now applying to the output weight matrix, which
        `exclude_output_from_decay` had exempted. So the panel run's configuration also changed, not only
        the matrix's — `4a_percell_training` §6 is not a like-for-like predecessor of its own next execution.
  - [ ] **Expect every CV number to move down, for two reasons at once** (item 7, 12.08.2026): the
        optimistic epoch selection is gone, and each fold now fits on ~104 lines instead of 122. A drop
        is the change working, not a regression. The two causes cannot be separated after the fact, so
        if the size of the drop matters, it needs `INNER_VAL_FRACTION` varied deliberately — not
        inferred from the difference to the old numbers, which also change target and panel.
  - [ ] **Lift the 28.07 panel-void banner here** (decided 12.08.2026, Selin). Item 6 rebuilt the panel,
        but the banner's live consequence is that *numbers* computed on the old one are void, and that
        stays true until a run exists on `outputs/panel/panel.csv`. Lifting it when the panel changed
        would have declared the numbers current a run too early.
- [ ] **R5 · Re-run the analysis notebooks that read the retrained outputs.**
      `analysis/evaluation/`: `diagnostics` (§5 dispersion) and `dreval_benchmark`. Both are held by
      **data, not code** — each reads `outputs/panel/panel_oof_predictions.csv`, which R4 writes and
      which does not exist (`outputs/panel/` holds `panel.csv` and `literature_panel_candidates.csv`
      and nothing else). Neither has a code blocker left.
      *(Corrected twice, both times downward. 12.08.2026, R5 hold lifted for exactly this by Selin: this
      read "broken twice over (imports the archived `dreval_normalize`, hardcodes the removed `'auc'`)",
      and neither was true — the module is **live**, archived and restored paper-only the same day, and
      `e804f07` fixed the `'auc'` literal. **13.08.2026:** the one remaining blocker — "raises on its
      import cell, because three of the functions it takes from `dreval_normalize.py` were deleted with
      the cell-line-effect diagnostic" — is **also gone**. It imports `load_panel`, `load_oof` and
      `normalized_evaluation`; all three are present, at lines 93, 108 and 190. `diagnostics` likewise
      uses `DEFAULT_CTRP_SCORE` throughout, its one `'auc'` sitting inside a comment about that
      replacement. **A stale "raises" is the costly direction of wrong** — it sends the next reader to
      repair working code. Established by reading the cells and the module's exports; neither notebook
      has been **executed**, and neither can be until R2 writes the `auc_cc` targets file.)*
      `analysis/harmonization/drug_coverage` — **not optional**, the line count moves 180 → 181.
      `2_drug_selection` does **not** need the retrained outputs — it reads only the
      response CSV — but re-run it anyway once preprocessing has, so the panel is regenerated against
      the same artifacts as everything else and `panel.csv` cannot silently predate them.
      **Two notebooks still have their stored outputs cleared** (11.08.2026), because their score
      literal changed to `auc_cc` and the old results could not be refreshed under the freeze — the
      code is correct and only the results are missing: `4a_percell_training` (c1) and
      `analysis/qc/verify_variants` (c24).
      **Archived and not re-run:** `target_comparison`, `ablations_and_rescue`, `replicate_variation`
      (targets that no longer exist, [why](./steps/corrections-and-dead-ends.md#retired-code-paths)),
      and — since 12.08.2026 — `learnability_filter`, `learnable_subset_training` and
      `panel_distributions`, superseded by the
      [rebuilt panel](./steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026).
      ⚠️ `panel_distributions` also held the justification for the density-weighting parameters
      (`alpha=0.5`, `cap=3`); **item 9 must re-derive it or drop the weighting**, since the notebook
      that evidenced it is no longer live.
      **Render every figure and look at it before anything is reported from it.**
- [ ] **R6 · Update the docs and the report from the refreshed artifacts.** Nothing is re-run here;
      everything is re-read from what R2–R5 produced. `report/results_numbers.tex`: `\NLines` 180 → 181,
      `\NVocab` 4,576 → measured (`\NVocabRepaired` then goes), every ρ / gap / DrEval macro — ideally
      via the extraction script from item 12 rather than by hand. `report/sections/03_methods.tex`: the
      revision block still says "no representation has been regenerated", false once R2 lands.
      **The report's numbers were stripped, not bannered (12.08.2026, Selin)**, so R6 is a rewrite and
      not an edit: `04_results.tex` was emptied to a withdrawal note, and `01_abstract.tex`,
      `05_discussion.tex` and `06_limitations_and_outlook.tex` had every quantitative claim removed.
      They are written once, from the regenerated artifacts.
      **The loss and the metric set both need writing (added 12.08.2026, audits 08–09).** The report
      says only that training uses "a masked mean-squared-error loss"; it does not say that the density
      weighting is swept as an arm rather than fixed, nor that four quantities are reported — order,
      order at the top, values, spread — where earlier versions reported rank correlation and error
      alone. The correction already in `06_limitations_and_outlook.tex` explains why that mattered;
      the methods section has to state the replacement.
      **`03_methods.tex` §Representation and model describes the architecture but not the training
      configuration** — no optimizer, learning rate, weight-decay grouping, early stopping or head-bias
      initialization appears anywhere in the report, so a reader cannot reconstruct the run. Write it
      here from `TrainConfig` and [Step 03](./steps/03-model-and-training-design.md#the-uncentred-target-is-handled-the-same-way-in-every-training-path)
      (noted 12.08.2026, audit 08).
      **Define the drug-panel counts as macros** rather than the inline digits currently in
      `03_methods.tex` §Drug panel selection (150 / 120 / 57 / 44 / 11, the 90 % cut, 91.2–98.3 %
      coverage, 102 and 15 unmatched compounds, 13 parent-CID matches) — via item 12's extraction
      script, since adding them by hand is the defect that item exists to remove.
  - [ ] **Re-read every number *derived* from the cell-line count — not just `\NLines` itself**
        (added 12.08.2026 by the docs audit, on Selin's ruling). R6's list above is **instance-based**:
        it names `\NLines`, the ⚠️ 181 markers and the `project_progress.md` note, i.e. the places
        180 appears *as itself*. A second class moves with it and is named nowhere — numbers **computed
        against** the line count, which are wrong the moment it changes even though they never print it.
        **The rule, because the list below is not the whole of it: any number derived from the cell-line
        count is provisional until R2–R5 regenerate, and R6 re-reads it from the artifact rather than
        adjusting it by hand.** Hand-adjustment is exactly the defect item 12 exists to remove, and a
        derived number is the easiest place for it to hide — nothing in the text looks stale.
        Known instances, to start from and not to stop at:
    - [ ] `\NLinesCV` (**153**) and `\NLinesFit` (**126**) in `report/results_numbers.tex` — both are
          70+15 % and 70 % of `\NLines`, so both move; 153 becomes ~154. They are cited in
          `02_introduction`, `03_methods` (§Data and §Evaluation) and `05_discussion`.
    - [ ] Prose of the shape *"X of the 180 …"*, where the 180 is load-bearing rather than decorative —
          [Step 03](./steps/03-model-and-training-design.md) has *"covers 153 of the 180 labelled
          lines"*; sweep [Step 05](./steps/05-multitask-results.md) and `project_progress.md` for the
          same shape, and the split-distribution tables (126 / 27 / 27 lines, and their cell counts)
          with them.
    - [ ] The 18 unassigned lines / 6,286 cells, which is `198 − \NLines` and therefore moves to 17.
    - [ ] ⚠️ **Not in this class, and must not be swept into it:** the drug-panel funnel and coverage
          figures in `03_methods.tex` §Drug panel selection. `2_drug_selection` reads the response CSV
          rather than a pipeline artifact, so those were **already computed on 181** — they are the one
          place where the corrected count is live today. Whether that is left as is or reconciled is
          Selin's own open item; this sub-item does not touch it.
      [Step 01](./steps/01-datasets-and-harmonization.md) the ⚠️ 181 marker;
      [Step 02](./steps/02-preprocessing-and-embeddings.md) the ⛔/⚠️ blocks on truncated, unseeded,
      symbol-limited embeddings; [Step 05](./steps/05-multitask-results.md) the ⚠️ 181 marker and every
      pre-sweep dated marker; `project_progress.md`'s "181 after the next sweep" note;
      [Corrections](./steps/corrections-and-dead-ends.md) — move the "repaired in code, takes effect at
      the sweep" entries to applied. **Lift the 03.08 freeze banner here**; the 28.07 panel banner lifts
      with item 6, not here. Rebuild `main.pdf` and check the log is clean.

## After the sweep — the one review item that needs new runs

Every other item above was settled by reading code, data or a paper. This one asks what a difference in
*input scale* does to a trained model, which only training can answer. It is parked here rather than left
open above, so the review closes on what the review can decide — and so it is not mistaken for something
that can be picked up early.

- [ ] **4A · The ~78× input-scale asymmetry between the two arms, under one shared learning rate.**
      Moved out of review item 4 on 10.08.2026. It qualifies every PCA-vs-scGPT claim made so far: if one
      arm reaches the optimizer with values ~78× larger than the other, one learning rate is not one
      setting, and an arm can look worse for a reason that has nothing to do with the representation.
      Two halves with different prerequisites:
  - [ ] **Measure the asymmetry — needs the sweep's preprocessing pass, no training.** The 78× is
        stale: `sc.pp.scale` entered the PCA path 05.08.2026 and standardizing genes changes component
        magnitudes, so the number has to be re-measured, not re-used. Once PCA is recomputed it is a
        read off the stored record — `sqrt(variance[i])` **is** the standard deviation of PCA
        coordinate *i*, so the arm's input scale no longer costs a multi-GB matrix load
        ([Step 02](./steps/02-preprocessing-and-embeddings.md#the-fit-is-stored-not-only-the-coordinates-10082026)).
  - [ ] **Test whether it matters — needs the retrained baseline to compare against.** Not before the
        sweep's retraining step: a controlled comparison needs a current baseline from the current
        pipeline, and every run under `runs/` predates the code that now exists.
  - [ ] **Recorded 10.08.2026 (Selin): how many components should PCA keep?** 512 was chosen to match
        scGPT's embedding width — a fair basis for a controlled comparison, but not a statement about
        the data, and the item has never asked it. Once the sweep writes `variance_ratio`, *what does
        512 retain* and *where does the curve flatten* are a two-line read off `uns["pca_fits"]`.
        Answer it here, since changing the width changes the very asymmetry this item measures.


Action list. Scientific narrative + full numbers live in
[project_progress.md](./project_progress.md) and [`docs/steps/`](./steps/); this is the running tasks.
A standalone write-up of the current state is `../report/` (LaTeX → `main.pdf`).

> **What every item here is ultimately for.** Two questions, and each task should be traceable to one:
> **Q1 — is scGPT a viable representation for drug-response prediction at all**, against the standard
> dimensionality reduction? **Q2 — does a model trained on single cells learn cellular heterogeneity
> *implicitly*?** Q2 is the clinically consequential one, because relapse is driven by rare surviving
> subpopulations rather than by the average cell. As of 27.07.2026 Q1 has evidence (scGPT clears the
> ridge control, replicated on an independently chosen drug panel) and **Q2 has none — the current
> objective penalizes the within-line variation it would have to express**, which is why MIL is next
> and not one of the several cheaper items below.

## FAIRER — data stewardship (28.07.2026)

The project works under **FAIRER**: FAIR as defined by Wilkinson et al. 2016, plus **E**thical and
**R**eproducible. FAIRER has no canonical primary source; it is treated here as a convention in
health-AI contexts, with FAIR-Health (Holub et al. 2018) as the closest published anchor. Both are in
`references.bib`. Full six-letter assessment of where the project stands: see item H below.

**Already closed (28.07.2026):** MIT `LICENSE`; `CITATION.cff` with ORCID and both affiliations, tagged
`v0.1.0`; `references.bib` as the single bibliography, with the report generated from it; and the
source-data terms of use checked against the providers and recorded in
[Step 01](./steps/01-datasets-and-harmonization.md) — CTRPv2 and PRISM CC BY 4.0, SCP542 unrestricted
under the portal ToS, **GDSC granting no redistribution**, **DrugBank non-commercial**.

- [ ] **D · Column-level schema for the derived outputs.** *(Re-scoped 30.07.2026 — the original item
      asked for an `outputs/README.md`, which now exists.)* `notebooks/outputs/README.md` already maps each
      directory to its producing notebook and explains what the figures show. What is still missing is
      **column** documentation: a reader must open a CSV to learn what its columns mean, and the targets
      h5ads document neither their `obsm`/`uns` keys nor their units. Extend the existing README rather
      than adding a file. (FAIRER: **F**)
- [ ] **E · Data-availability path.** The raw data sits under `~/Desktop/OncoTox/data`, gitignored and
      several GB, so nothing is reconstructible by a third party. *(Partly closed 30.07.2026: which
      release of each source, retrieved when and from where, is now recorded in
      [Step 01](./steps/01-datasets-and-harmonization.md#provenance--what-was-retrieved-from-where-when)
      instead of decaying in a dated log.)* Still open: **a script that rebuilds the derived artifacts
      from those sources**, and a separate decision on whether anything gets deposited (Zenodo) at
      publication. (FAIRER: **A**)
- [x] **F · One residual on the data terms — closed 10.08.2026.** SCP542 carries **no** study-level
      licence beyond the portal ToS. Checked in Kinker et al. 2020 itself: the Data availability
      statement (p. 13) names no licence and declares no restriction, and the Reporting Summary (p. 28)
      lists none under the heading that asks for restrictions. It also yields a persistent accession we
      did not have — **GEO `GSE157220`**. Both recorded in
      [Step 01](./steps/01-datasets-and-harmonization.md#provenance--what-was-retrieved-from-where-when).
      (FAIRER: **E**)
- [ ] **G · State the FAIRER commitment in the docs, not only in the bibliography.** The sources are in
      `references.bib`, but no document says the project works under FAIRER or what the six letters mean
      here. Natural home is `project_progress.md`, next to H. (FAIRER: **R**)
- [ ] **H · Record the six-letter assessment.** Where the project stands against F, A, I, R, E, R — with
      the two live gaps named: name-based joins instead of persistent identifiers (item B), and the fact
      that `set_seed` sets no determinism flags, so the PCA arm moved across four identical runs
      (0.313 / 0.315 / 0.317 / 0.320) at a fixed seed. **Decision needed when this is written:** own
      section in `project_progress.md`, or a new `docs/steps/` file. (FAIRER: all six)

**Determinism — decided 28.07.2026, no action:** do *not* force `torch.use_deterministic_algorithms`;
several operations have no deterministic MPS implementation and it would likely fail rather than merely
slow training. Report the non-determinism instead, and measure it at the specific points where it
matters — folded into review items 10 and 12.

## Next up — prioritized (15.07.2026, from the progress-report feedback)

> **Framing that governs this list:** *more performance ≠ a bigger model.* Model-side tuning is
> demonstrably closed (see "Model-side tuning is closed" below); the levers are label-/data-side.
> The audience's "bigger MLP / more capacity" suggestion was already tested and is flat. *(This line
> previously named MIL as the one remaining capacity lever. Corrected 11.08.2026: MIL is Q2's
> instrument, not a performance lever — see S2. There is no untested capacity lever.)*
>
> **Update 25.07.2026:** the report's first next-step — *drug selection from the literature instead of my
> filter* — is **done as a definition** (8-drug panel, see "Next focus" below); the training run on it is
> pending and should be paired with S1. The panel is literature-anchored but **not yet label-blind**, so
> train-only selection remains blocking for any headline number.
>
> ⚠️ **Superseded 12.08.2026.** That 8-drug panel was voided three days later and replaced by the
> [11-drug panel](./steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026),
> which *is* label-blind — so the "not yet label-blind" caveat above is resolved, not outstanding.

1. **S1 — DrEval-aligned target (the top performance lever).** Train on the double-normalized residual
   `resid[i,j] = auc[i,j] − (μ_drug[j] + μ_line[i] − μ_global)`, means computed **train-only per fold**,
   so the objective matches DrEval's normalized metric. Motivated by `dreval/dreval_normalized.csv`
   (~20% of the signal is pure cell-line effect; `kx2-391` is entirely that artifact). New `--score`
   option in `ctrp_to_h5ad.py` (pattern: `_zscore_per_drug`, `DEFAULT_CTRP_SCORE`); compare in `dreval_benchmark`.
   Success = normalized DrEval ρ rises above the current scGPT value without inflating the raw
   correlation (the bar itself lives in [Step 05](./steps/05-multitask-results.md), not here).
2. **S2 — MIL / attention pooling over a line's cells (bag of cells → line label).** **Not a performance
   item — it is the instrument for Q2**, reframed 11.08.2026; ρ against ridge and the per-cell MLP is a
   floor, not the success criterion. Full item, with the open decision on what a positive Q2 result is:
   *Agreed plan, Step 2*. *(The control values this entry used to quote were measured on the voided
   panel and are withdrawn.)*
3. **S3 — More independent cell lines** — SCP542×CTRPv2 caps at 180; CTRPv2 has ~1,100. Attacks the
   real ceiling. *(Overlaps the scDEAL/label-side lever under "Levers / later".)*
4. **S4 — Diagnostic explainability (now, low-risk)** — where errors concentrate (drugs/tissues/lines
   with ρ<0) and how much residual error is the cell-line effect. Uses existing per-drug-ρ CSVs.
5. **Later — Discovery explainability (gated)** — gene-level XAI only once per-drug ρ is substantially
   higher and stable, else it interprets noise. *(Project-plan stretch goal; see Step 07.)*

**Communication fixes:** rebuild the single overview image (data → rep → shared trunk → K heads →
out-of-fold eval); make the rescue figure show the **ceiling** (ridge = MLP, no-reg memorizes), because
the "model-side is closed" message did not land in the talk.

**Working agreement:** for any analysis beyond the explicitly agreed fix — *especially how a plot is
computed/displayed* — confirm first, don't decide silently. (Unasked decisions produced process bugs and
undefendable slides last round.)

## Agreed plan — order of work (27.07.2026)

**Governing rule: never change the target and the architecture in the same run.** The June result took
weeks to unpick precisely because two changes landed together; if MIL and a new target arrive at once and
the number moves, the cause is unattributable.

**Step 1 — target + loss weighting, on the existing per-cell MLP.** Everything else unchanged
(architecture, splits, optimizer, batching), so the change is attributable.

**Decision (27.07.2026): `auc_z` is retired; the target is raw `auc`.** The decomposition that retired it
is in [Corrections](./steps/corrections-and-dead-ends.md#auc_z-as-the-training-target); the current target
definition and the two mechanics it forces are in
[Step 03](./steps/03-model-and-training-design.md#target-y--the-response-score-and-at-what-resolution-it-is-defined).

The per-drug variance-weighting prerequisites this plan originally carried (`σ_noise`, `w_j = 1/σ_j²`,
then `r_j/σ_j²`) were **never needed** — per-drug variance is a K=545 problem and the panel dissolved it:
[Corrections](./steps/corrections-and-dead-ends.md#per-drug-variance-weighting--dissolved-by-a-scope-change-never-needed).

### Step 1 — run 27.07.2026 (`notebooks/4a_percell_training.ipynb`), open items only

Raw `auc` winsorized at 1.1, 8-drug panel, per-sample inverse-density weights fitted per fold on training
lines only, output layer excluded from weight decay, head biases initialized to train-fold per-drug means,
one seed. ⛔ **Run on the [voided panel](./steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it), so its numbers are provisional.**
What it settled — the collapse was a head-count effect, the ridge tie replicated, density weighting is a
clean null — is written up in [Step 05](./steps/05-multitask-results.md) and
[Corrections](./steps/corrections-and-dead-ends.md#inverse-density-loss-weighting-improves-ranking).

- [ ] **Reproducibility.** The PCA-unweighted arm is not bit-reproducible on `mps` — 0.313 / 0.315 /
      0.317 / 0.320 over four identical runs, every other arm exact. **Do not report the sign of the
      weighting deltas**; they lie inside that band. Cause (PCA peaks at epoch 1) in
      [Corrections](./steps/corrections-and-dead-ends.md#inverse-density-loss-weighting-improves-ranking).
- [ ] **Seeds.** One seed against ±0.04 documented seed variation. Repeat over ≥ 3 seeds before
      scGPT − PCA (+0.061) or scGPT − ridge (+0.077) is quoted as a margin. **Blocking for any headline
      number.**
- [ ] **Report raw + normalized** (DrEval) on this panel, once seeds are in.

**Step 2 — MIL / attention pooling: the instrument for Q2, against the target fixed in Step 1.**

**Reframed 11.08.2026 (Selin): MIL is how Q2 gets answered, not a capacity lever.** The per-cell MLP
hands every cell of a line the same label, so the objective penalizes precisely the within-line variation
Q2 asks about — under that architecture Q2 is not merely unanswered, it is unanswerable. A bag of cells →
one line label is the smallest change that lets the model express heterogeneity rather than be punished
for it. What this item said before the reframing — *"if MIL beats neither control, the single-cell
resolution has again failed to justify itself"* — scored a Q2 experiment on a Q1 criterion.

- **The controls are a floor, not the criterion.** RidgeCV on cell-line mean embeddings and the per-cell
  MLP still run on the same panel and the same folds: if MIL lands far below them, its attention weights
  are not evidence of anything. But *tying* them is not a failure here — it says the heterogeneity
  signal does not raise line-level ρ, which is a narrower claim than the one this item used to make.
  Their values are pending: the numbers previously quoted here were measured on the
  [voided panel](./steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it),
  so the floor is re-measured at R4, on the panel item 6 rebuilt on 12.08.2026.
- **✅ What counts as a positive Q2 result — decided 12.08.2026 (Selin), and written down before any
  model exists.** It lives in `notebooks/4b_mil_training.ipynb` §2, next to the code it will judge,
  rather than here. Four stages with distinct roles: a **synthetic positive control** as precondition
  (bags mixed from two lines of known response — without it a negative cannot be told apart from "the
  method cannot detect heterogeneity"), **within-line spread** as necessary condition, **cross-seed
  reproducibility** as the test, and **confound regression** as a veto, since predictions that replicate
  *and* are explained by library size are a sequencing artifact. Kinker program alignment is
  characterisation and gates nothing.
  **Two consequences of the architecture choice** (instance-level, not attention pooling — every cell
  gets a predicted response rather than a weight): the subpopulation-predictivity test is *not* usable,
  because selecting cells by their predicted value and scoring that against the line's truth is biased
  by construction; and the biology readout needs no top-k, since both sides are continuous and can
  simply be correlated. **One blank remains**, and it is the only number needing judgement: the
  threshold at which the synthetic control counts as recovered. Stages 1 and 2 are then expressed as
  fractions of what it achieves, rather than invented in advance.
- MIL makes the per-line weighting problem disappear structurally (one bag = one example), so the
  82× cells×labels imbalance needs no separate fix under it.
- Open design decisions, to settle before building: fixed bag size with per-epoch subsampling (acts as
  augmentation) vs. padding + mask, given 56–1,990 cells per line; and the optimizer regime, since an
  epoch becomes ~120 bags instead of ~34,000 cells so the current epoch/LR/early-stop settings do not
  carry over.
- Attention weights are the readout for *which subpopulation drives response* — the link to the relapse
  motivation and to the annotated heterogeneity programs.

**Not in scope for either step:** the base-quantity question (AUC vs EC50/Emax, T3) and learned/adaptive
task weights. Both are deferred deliberately; adaptive weights estimate *residual* variance, which mixes
label noise with model error and risks a self-reinforcing loop, especially at 545 tasks over ~120 bags.

## Target & drug-selection defects found 27.07.2026 — all three closed

> ✅ **Closed 11–12.08.2026, and none was closed the way it is written below.** T3 was answered by
> [audit 05](./steps/01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026):
> the target moved to DrEval's CurveCurator re-fit, so the CTRP columns T3 proposed choosing between are
> no longer what we read. T2 dissolved with `auc_z`, [retired](./steps/corrections-and-dead-ends.md#auc_z-as-the-training-target)
> — there is no denominator left to fix, and between-drug scaling is now audit item 9's question. T1's
> replacement criterion (`auc_std` + coverage) was itself rejected on 12.08.2026 for selecting on our own
> labels; the panel is built on
> [FDA approval and published determinants](./steps/01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026)
> instead. Kept unedited below as the diagnosis that was correct — the gate did select on potency — with
> the prescriptions superseded.

Both found by asking why `nutlin-3` was rejected by the filter. Write-ups:
[Corrections](./steps/corrections-and-dead-ends.md#auc_z-as-the-training-target),
[Step 05](./steps/05-multitask-results.md#the-learnability-gate-measured-the-wrong-quantity-27072026).

- [x] ~~**T1 — Replace the kill/spare gate with `auc_std` + coverage.**~~ The gate filters on absolute
      potency; `auc_z` subtracts the per-drug mean and Spearman reads only the ordering, so we selected
      on a quantity the model never sees. `nutlin-3` σ = 0.147 vs `dasatinib` σ = 0.155 — same spread,
      rejected only because it is cytostatic. **116/545** drugs have zero kills but σ ≥ 0.10 and
      coverage ≥ 90 %. Everything downstream (10-drug panel, 8-drug literature panel, all K=10 numbers)
      rests on the old gate and must be re-derived.
- [x] ~~**T2 — Fix the `auc_z` denominator.**~~ Dividing by `auc_std` forces noise-floor drugs to variance 1
      and hands them full weight in the shared loss — the mirror of the June σ² bug. Use
      `sqrt(auc_std² + σ_noise²)`, or weight each drug by its reliable variance fraction.
      `σ_noise` is estimable **pooled** from the 7,708 replicated (line × compound) fits (2.0 % of
      387,130) in `v20.data.curves_post_qc.txt`; per-drug is not feasible.
- [x] ~~**T3 — Reconsider AUC as the target**~~ (raised by Selin's supervisor, DrEval co-author; DrEval lists
      inconsistent viability data as an obstacle and recommends **CurveCurator**). AUC conflates potency
      with efficacy, and CTRP's own fit already separates them in the file we parse:
      `apparent_ec50_umol`, `pred_pv_high_conc` (≈ Emax), `p3_total_decline`. Top test concentration
      spans **0.13–600 µM** across the 545 — harmless across drugs (z-scoring is within-drug) but it
      compresses spread *within* a drug, so it is a **cause of T2**, not a separate issue.
      Interacts with S1 — decide the target once, for both.

## Completed work — where it is written up

This file no longer restates finished work. Completed items live where their numbers are:

| What was done | Written up in |
|---|---|
| 8-run matrix, matched `(128,64)` trunk + 512-d width, shared `split_ctrp` | [Step 05](./steps/05-multitask-results.md) — conclusions superseded, see [Corrections](./steps/corrections-and-dead-ends.md#the-8-run-matrix-conclusions) |
| 5-fold GroupKFold CV, per-drug correlation, gene-set sweep | [Step 05](./steps/05-multitask-results.md) |
| Target distribution + per-drug coverage & learnability | [Step 05](./steps/05-multitask-results.md#the-learnability-gate-measured-the-wrong-quantity-27072026) |
| Cancer-type UMAPs and the latent-space validation | [Step 02](./steps/02-preprocessing-and-embeddings.md#latent-space-validation-umap-fig-3--fig-4) |
| 190 vs 180 cell-line overlap, source licences, the externally shared drug list | [Step 01](./steps/01-datasets-and-harmonization.md) |
| Target score → `auc_z`, and its retirement | [Step 03](./steps/03-model-and-training-design.md), [Corrections](./steps/corrections-and-dead-ends.md#auc_z-as-the-training-target) |
| Learnability filter + the 5-/10-drug results | [Step 05](./steps/05-multitask-results.md), [Corrections](./steps/corrections-and-dead-ends.md#the-1307-five-drug-numbers) |
| The 13.07 "net read" that replaced the label-ceiling conclusion | [Corrections](./steps/corrections-and-dead-ends.md#neither-representation-ranks-cell-lines--the-k545-null-result) |

## Next focus — make the 5-drug result honest (13.07.2026)

The `learnability_filter`/`learnable_subset_training` numbers are a **best-case diagnostic**: the 5 drugs were selected using all 180 lines,
val/test included, so the selection saw held-out labels. Turning it into a reportable result:

The panel that was built to fix this is **void** — see
[Corrections](./steps/corrections-and-dead-ends.md#the-8-drug-literature-panel-and-every-number-computed-on-it)
for what it was, why it failed, and the rebuild criterion. The seed checks and the three-target
comparison are written up in [Step 03](./steps/03-model-and-training-design.md) and
[Step 05](./steps/05-multitask-results.md).

- [ ] **Train-only selection** — run the drug-selection criterion *inside each CV fold* (train lines only)
      and re-measure. If the effect survives, it is real. **Still the blocking item** for any headline
      number, on any panel.
> ⚠️ **The three boxes below were written for a panel that no longer exists (13.08.2026, Gate 1
> sweep).** They target the **5-drug learnability subset**, whose selection criterion was
> [retracted](./steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability)
> — it measured potency, not rankability — and the 8-drug panel that replaced it is void. The live
> panel is the **rebuilt 11-drug literature panel** (item 6, 12.08.2026), selected on published
> evidence rather than on our own response values, which dissolves the label-dependency each of these
> was designed to remove. **They are kept rather than deleted because two of them describe procedures
> that would still be worth having on the current panel** — externalizing the spread condition, and
> asking where the signal dies as the criterion loosens — but neither is scoped to it, and neither
> should be run as written. **Re-scope or retire: Selin's call, and it is not a rerun blocker.**
>
> **`Train-only selection` above is deliberately NOT in this group** — it says "on any panel", so it
> generalises past the artifact it was written for and remains the blocking item it claims to be.

- [ ] ⬜ **Externalize the spread requirement** *(needs re-scoping — see the note above)* — re-derive the panel with the spread condition measured on
      **GDSC2** (`data/GDSC2_fitted_dose_response_27Oct23.xlsx`) or PRISM instead of on the CTRP labels we
      train on. Cheaper than fold-internal selection and would make the panel genuinely label-blind.
- [ ] ⬜ **Loosen to ~20–50 drugs** *(needs re-scoping — see the note above)* — a handful is a diagnostic, not a model. Where does the signal die as
      the criterion relaxes? (`outputs/archive/learnability/ctrp_drug_learnability_auc.csv` is already ranked for
      this, though it is ranked on the [discredited gate](./steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability).)
- [ ] ⬜ **Re-run the full 8-run matrix + CV on the current target** *(needs re-scoping — see the note above; and `4a` §B is already scheduled at R4, so check this is not a duplicate)* — for a like-for-like against the
      `mean_pv` Steps 04–05 numbers. **Expect this to overturn them, not refresh them**
      ([why](./steps/corrections-and-dead-ends.md#the-8-run-matrix-conclusions)).
      ⚠️ *This item originally said `--score auc_z`, which is retired — use the current default.*
- [ ] ⬜ **More seeds + a wider drug set** *(the seed half is already R4's; the "wider drug set" half needs re-scoping — see the note above)* before scGPT > PCA becomes a headline claim. Pair it with
      train-only selection.

## Model-side tuning is closed (13.07.2026)

> ⛔ **The evidence for this section is void (12.08.2026, audit 08) — the instruction stands only as an
> argument, not as a measurement.** The notebook behind every number below chose each model's checkpoint
> on the fold it then scored, on the retired `auc_z` target, over the five drugs of the discredited
> learnability gate, and it can no longer be run.
> [Corrections](./steps/corrections-and-dead-ends.md#the-evidence-that-closed-model-side-tuning).
> **The ridge tie is the load-bearing part and it is the part most at risk:** ridge has no early stopping
> so it was never flattered, the MLP was, and on PCA the two are equal — the honest ordering may put ridge
> ahead. Re-derived at R4 by the minimal re-run in item 8C. The argument that survives without the
> numbers is that ~153 independent labels cannot support architecture search, which is a fact about the
> data rather than a result.

Four knobs — regularization, capacity, batch size, sample reweighting — are **all flat**, and `RidgeCV` on
the 150 cell-line mean embeddings **ties the PCA MLP**. Tables, the rescue test on the broken setting, and
what each result rules out:
[Step 03](./steps/03-model-and-training-design.md#these-hyperparameters-are-not-worth-tuning-ablated-13072026)
and [Corrections](./steps/corrections-and-dead-ends.md#the-model-is-over-regularized-or-too-small).
**Don't spend more time on architecture or hyperparameters** — the remaining levers are label-side.

- [ ] **Averaging a line's cells into one vector currently loses nothing** — `RidgeCV` on the 150
      line-mean embeddings ties the per-cell MLP, so at line-level ρ the single-cell dimension is not
      earning itself. *(This item used to propose MIL as the fix. Moved 11.08.2026: MIL is Q2's
      instrument, not a performance lever — the item lives in* Agreed plan, Step 2 *and is not scored
      on beating ridge.)* What stays open here is narrower: if line-level ρ is all that is ever
      reported, the per-cell framing needs a justification that does not depend on Q2 succeeding.
- [ ] **Add ridge (line-level) to `4a_percell_training` §B's comparison tables** so every future claim is scored against it.
- [x] ⛔ **VOID — there is no z-scoring left to fix** (13.08.2026, Gate 1 sweep). It read:
      *"z-score train-only. The per-drug mean/std currently use all 180 lines, val/test included —
      mild leakage."* That leak was real while the target was `auc_z`, which was **retired
      27.07.2026**; `layout.CTRP_SCORES` is now `('auc_cc', 'ln_ic50_cc')` and no per-drug mean or
      standard deviation is estimated anywhere in the target path. The concern did not go away by
      being fixed — the target that carried it was withdrawn, which is a different thing and worth
      distinguishing: if a standardized target ever returns, this leak returns with it.
- [ ] *(Stretch)* cluster cell lines by response and **stratify train/val/test** (high/med/low) for
      lower-variance evaluation.

## DrEval alignment (14.07.2026)

Paper: Bernett, Iversen, Picciani, **Wilhelm**, Baum, List — *Critical evaluation of drug response
prediction models with DrEval*, Nat. Commun. 2026. Half of published models don't beat a naive
drug-mean + cell-line-mean predictor. **Our ridge ≈ MLP finding is the field's norm, independently
reproduced.** Our split *is* their LCO.

Benchmarked with the real package (`drevalpy` 1.5.1) in `notebooks/analysis/evaluation/dreval_benchmark.ipynb`, and
separately re-normalized to remove the cell-line effect. Results in
[Step 05](./steps/05-multitask-results.md); the first run's leak and the numbers it invalidated in
[Corrections](./steps/corrections-and-dead-ends.md#the-first-dreval-benchmark--a-val-split-leak).

- [ ] **Make `NaiveMeanEffects` the default baseline** in `train_multitask.py` (currently: per-drug mean,
      too weak).
- [ ] **Report raw + normalized** correlations everywhere from now on.
- [ ] **Run DrEval on all 545 drugs**, not just the best-case 5 — and with their LTO / LDO settings.
- [ ] *(Consider)* their other splits — LTO (leave-tissue-out) and LDO (leave-drug-out). LDO would test
      whether anything generalizes across chemical space; DrEval found **no model** beats naive there.
- [ ] *(Consider)* **CurveCurator** for standardized dose-response fitting, as they recommend.

## Open scientific questions

*Moved here 30.07.2026 from `project_progress.md`, which now indexes rather than carries these.
Unlike the items above these are questions, not scheduled work.*

- Does multi-task help or hurt paclitaxel? Single-task on `split_ctrp` exists (scGPT 0.0406, PCA 0.0372
  on `hvg5000`); compare against the paclitaxel **head** inside the K=545 run
  ([Step 05](./steps/05-multitask-results.md)).
- Which low-coverage heads should be dropped or down-weighted? The ≈16-line drugs (n_val 221) are the
  unreliable ones — quantified in `notebooks/analysis/harmonization/drug_coverage.ipynb` and
  [Step 05](./steps/05-multitask-results.md).
- Should the loss move from uniform-per-entry to per-head / uncertainty weighting? *(Note: adaptive task
  weights are deliberately deferred — they estimate residual variance, mixing label noise with model
  error.)*
- Does HVG-5000 lose signal against the full transcriptome? The gene-set sweep says no
  ([Step 05](./steps/05-multitask-results.md)), but the `all_genes` arm has never been evaluated on the
  current target.
- When should PRISM/GDSC come in as additional masked heads
  ([Step 06 · A](./steps/06-planned-work.md#a-cross-database-integration))?

## Levers / later

- [ ] **Bulk RNA-seq pretraining / scDEAL-style denoising + domain adaptation** — attacks the
      noisy-label bottleneck (the real ceiling). **Promoted by the `ablations_and_rescue` ablations:** with model-side
      tuning closed and ridge-on-150-lines matching the MLP, the only remaining levers are label-side —
      above all **more independent cell lines** (SCP542×CTRPv2 caps at 180; CTRPv2 itself has ~1,100).
- [ ] **Cross-database PRISM** (masked multi-task) — [Step 06](./steps/06-planned-work.md#a-cross-database-integration).
      (GDSC not a modelling priority; was only for the externally shared list.)
- [ ] **XAI** — feature importance → resistance drivers — [Step 06 · B](./steps/06-planned-work.md#b-xai-and-feature-interpretability).
- [x] ~~Confirm scGPT input preprocessing in `gen_embeds.py` (raw counts vs CPM) so scGPT isn't
      handicapped.~~ **Answered 03.08.2026: CPM does not handicap scGPT, and the paper sanctions it.**
      Value binning uses per-cell quantiles of that cell's own non-zero values, so it is rank-based and
      invariant to any monotone per-cell transform — CPM, raw counts, `normalize_total`, `log1p` all
      give the same bins at a given RNG state. Cui et al. present binning as the *replacement*
      for TPM-normalization and log1p, and state that `X` "represent[s] both the raw and preprocessed
      data matrices before binning". Detail in
      [Step 02](./steps/02-preprocessing-and-embeddings.md#what-scgpt-is-fed-and-why-its-scale-does-not-matter).
      **Correction 10.08.2026:** this line previously said "byte-identical … (verified on 200 cells)".
      No such check exists, and the raw counts it names are not in this project — the claim is an
      argument from `scgpt/preprocess.py`, not a measurement, and the bins are not byte-identical
      because `_digitize` breaks ties with an unseeded-at-that-point RNG
      ([Corrections](./steps/corrections-and-dead-ends.md#the-scgpt-binning-invariance-was-verified-on-200-cells)).
- [ ] Regenerate scGPT embeddings — **no longer optional and no longer identical**: `gen_embeds.py` now
      seeds with 42 and runs on MPS, so every embeddings file on disk predates the change. Scope
      (all_genes only vs all five variants) still undecided.
- [ ] *(Optional)* re-run `split_paclitaxel` single-task to fill [Step 04](./steps/04-single-task-results.md)'s
      PCA column, or retire that progression.

## Roadmap (project plan)

- [ ] Cross-database integration — PRISM then GDSC, efficacy + toxicity ([Step 06](./steps/06-planned-work.md#a-cross-database-integration)).
- [ ] XAI / feature interpretability ([Step 06 · B](./steps/06-planned-work.md#b-xai-and-feature-interpretability)).
- [ ] Foundation model + clinical fine-tuning ([Step 06 · C](./steps/06-planned-work.md#c-foundation-model-and-clinical-fine-tuning)).
