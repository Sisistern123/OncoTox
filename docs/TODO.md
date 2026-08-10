# OncoTox — TODO

> # ⛔ 28.07.2026 — THE CURRENT DRUG PANEL IS VOID
>
> **The 8-drug literature panel is discarded.** Its candidate list was ranked by `min(kill, spare)` on
> our own response values before the literature criterion was applied, so the selection inherited the
> label dependency it was built to remove — and 32 approved or clinical compounds, `nutlin-3` among them,
> were excluded for the wrong reason.
>
> **Everything computed on it is therefore provisional:** the step-1 run (`3_panel_training`), the
> distributions and weighting design (`panel_distributions`), the dispersion figures (`diagnostics` §5), the panel rows in
> [Step 05](./steps/05-multitask-results.md), and the corresponding numbers in the report. Do not quote
> any of them, and do not build on them, until the panel is rebuilt.
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
        `notebooks/data_and_harmonization/replicate_variation.ipynb` →
        `outputs/data/replicate_variation.{png,csv}`; written up in
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
        (`scripts/preprocessing/gene_symbols.py`; `scp542_conversion.py` annotates,
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
- [ ] **4 · Representations** — PCA components, scGPT embedding generation and its input format
      (raw counts vs CPM — never confirmed), and **the ~78× input-scale asymmetry between the two under
      one shared learning rate**, which is untested and qualifies every PCA-vs-scGPT claim we have made.
- [ ] **5 · Target** — AUC vs EC50 vs Emax: AUC conflates potency with efficacy and the tested
      concentration range spans 0.13–600 µM. Winsorizing threshold. Are all statistics per fold?
- [ ] **6 · Drug selection — REBUILD, this is the main deliverable.** Pool on coverage and spread only,
      **no kill counts at any point**, then apply the literature criterion to that pool. Decide the
      spread threshold explicitly. Expect `nutlin-3`, `oxaliplatin`, `bortezomib` and others to re-enter.
- [ ] **7 · Splits** — grouped by cell line, test held out, folds shared between model and baselines.
      Confirm nothing leaks through statistics computed outside the fold.
  - [ ] **Handed over from item 3 (10.08.2026): three fits, one question — what may a fit see?** All
        three are established facts, none is decided, and they should be decided together rather than
        separately, or they will get three inconsistent answers. All are unsupervised (no fit sees a
        response label), which is why standard pipelines tolerate them; and all bias **toward** the PCA
        control, since scGPT's per-cell binning draws on no other cell, so any scGPT-over-PCA margin
        measured today is conservative. Detail:
        [Step 02](./steps/02-preprocessing-and-embeddings.md#what-transform-pca-sees--corrected-05082026).
        1. **HVG selection is all-cells**, for both arms — the one fit the two representations share.
        2. **The cross-validated PCA is all-cells.** The fixed splits were fixed 05.08.2026
           (`X_pca_train_ctrp`); CV folds are drawn at training time, so five fold-specific matrices
           cannot be stored and `resolve_rep` leaves them on the leaky `X_pca`. Every CV number carries it.
        3. **Cells that never train are in all three fits** — the 17 lines / 6,073 cells (18 / 6,286 on
           disk today) with no CTRPv2 label, 11.4 % of the atlas. Not a test leak; a separate question
           about whether the representation should be shaped by data the model never sees.
        Any change here alters the gene set and therefore every number, so it lands in the sweep.
- [ ] **8 · Model** — architecture, capacity, and whether the shared-trunk multi-head design is still the
      right one for a small panel.
- [ ] **9 · Loss** — masking, per-sample weighting (the density weighting was a null — drop or keep?),
      per-line weighting (the 82× artifact), and what the objective actually rewards.
- [ ] **10 · Training** — optimizer, weight decay groups, epochs, early stopping, and the `mps`
      nondeterminism that moves the PCA arm across identical runs.
- [ ] **11 · Evaluation** — metric, cell→line aggregation, baselines (ridge, `NaiveMeanEffects`),
      dispersion across folds *and* drugs, raw vs normalized.
- [ ] **12 · Reproducibility** — seeds, determinism, what is derived in code versus typed by hand.
      Anything that exists only as a shell command is not a result.

**Working agreement for the session, restated because it was broken today:** no step of the analysis gets
decided silently. If a choice affects what enters the model or how a number is computed, it is proposed
first and agreed before it is executed.


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
> The audience's "bigger MLP / more capacity" suggestion was already tested and is flat — only MIL
> (S2) is a genuinely untested capacity lever.
>
> **Update 25.07.2026:** the report's first next-step — *drug selection from the literature instead of my
> filter* — is **done as a definition** (8-drug panel, see "Next focus" below); the training run on it is
> pending and should be paired with S1. The panel is literature-anchored but **not yet label-blind**, so
> train-only selection remains blocking for any headline number.

1. **S1 — DrEval-aligned target (the top performance lever).** Train on the double-normalized residual
   `resid[i,j] = auc[i,j] − (μ_drug[j] + μ_line[i] − μ_global)`, means computed **train-only per fold**,
   so the objective matches DrEval's normalized metric. Motivated by `dreval/dreval_normalized.csv`
   (~20% of the signal is pure cell-line effect; `kx2-391` is entirely that artifact). New `--score`
   option in `ctrp_to_h5ad.py` (pattern: `_zscore_per_drug`, `DEFAULT_CTRP_SCORE`); compare in `dreval_benchmark`.
   Success = normalized DrEval ρ rises above the current scGPT value without inflating the raw
   correlation (the bar itself lives in [Step 05](./steps/05-multitask-results.md), not here).
2. **S2 — MIL / attention pooling over a line's cells** — the only untested capacity lever (bag of cells
   → line label). Must beat the ridge baseline **and** the per-cell MLP. On the literature panel
   (27.07.2026) that is ridge 0.306 / 0.299 and MLP 0.316 / 0.377 for PCA / scGPT — the 0.342 quoted
   here previously was the old 10-drug panel. *(Detailed item under "Model-side tuning".)*
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

### Step 1 — run 27.07.2026 (`notebooks/3_panel_training.ipynb`), open items only

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

**Step 2 — MIL / attention pooling, against the target fixed in Step 1.** Controls that must both be
beaten: the per-cell MLP and RidgeCV on cell-line mean embeddings. If MIL beats neither, the single-cell
resolution has again failed to justify itself and that is the reportable result.

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

## Target & drug-selection defects found 27.07.2026 (do these before any new headline number)

Both found by asking why `nutlin-3` was rejected by the filter. Write-ups:
[Corrections](./steps/corrections-and-dead-ends.md#auc_z-as-the-training-target),
[Step 05](./steps/05-multitask-results.md#the-learnability-gate-measured-the-wrong-quantity-27072026).

- [ ] **T1 — Replace the kill/spare gate with `auc_std` + coverage.** The gate filters on absolute
      potency; `auc_z` subtracts the per-drug mean and Spearman reads only the ordering, so we selected
      on a quantity the model never sees. `nutlin-3` σ = 0.147 vs `dasatinib` σ = 0.155 — same spread,
      rejected only because it is cytostatic. **116/545** drugs have zero kills but σ ≥ 0.10 and
      coverage ≥ 90 %. Everything downstream (10-drug panel, 8-drug literature panel, all K=10 numbers)
      rests on the old gate and must be re-derived.
- [ ] **T2 — Fix the `auc_z` denominator.** Dividing by `auc_std` forces noise-floor drugs to variance 1
      and hands them full weight in the shared loss — the mirror of the June σ² bug. Use
      `sqrt(auc_std² + σ_noise²)`, or weight each drug by its reliable variance fraction.
      `σ_noise` is estimable **pooled** from the 7,708 replicated (line × compound) fits (2.0 % of
      387,130) in `v20.data.curves_post_qc.txt`; per-drug is not feasible.
- [ ] **T3 — Reconsider AUC as the target** (raised by Selin's supervisor, DrEval co-author; DrEval lists
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
- [ ] **Externalize the spread requirement** — re-derive the panel with the spread condition measured on
      **GDSC2** (`data/GDSC2_fitted_dose_response_27Oct23.xlsx`) or PRISM instead of on the CTRP labels we
      train on. Cheaper than fold-internal selection and would make the panel genuinely label-blind.
- [ ] **Loosen to ~20–50 drugs** — a handful is a diagnostic, not a model. Where does the signal die as
      the criterion relaxes? (`outputs/learnability/ctrp_drug_learnability_auc.csv` is already ranked for
      this, though it is ranked on the [discredited gate](./steps/corrections-and-dead-ends.md#the-learnability-gate-measured-potency-not-rankability).)
- [ ] **Re-run the full 8-run matrix + CV on the current target** for a like-for-like against the
      `mean_pv` Steps 04–05 numbers. **Expect this to overturn them, not refresh them**
      ([why](./steps/corrections-and-dead-ends.md#the-8-run-matrix-conclusions)).
      ⚠️ *This item originally said `--score auc_z`, which is retired — use the current default.*
- [ ] **More seeds + a wider drug set** before scGPT > PCA becomes a headline claim. Pair it with
      train-only selection.

## Model-side tuning is closed (13.07.2026)

Four knobs — regularization, capacity, batch size, sample reweighting — are **all flat**, and `RidgeCV` on
the 150 cell-line mean embeddings **ties the PCA MLP**. Tables, the rescue test on the broken setting, and
what each result rules out:
[Step 03](./steps/03-model-and-training-design.md#these-hyperparameters-are-not-worth-tuning-ablated-13072026)
and [Corrections](./steps/corrections-and-dead-ends.md#the-model-is-over-regularized-or-too-small).
**Don't spend more time on architecture or hyperparameters** — the remaining levers are label-side.

- [ ] **Make the single-cell dimension earn itself** — averaging a line's cells into one vector currently
      loses nothing. Test MIL / attention pooling over a line's cells (predict the line label from a *bag*
      of cells), which at least matches the true label resolution. If that doesn't beat ridge either, the
      per-cell framing needs a different justification.
- [ ] **Add ridge (line-level) to `2_training`'s comparison tables** so every future claim is scored against it.
- [ ] *(Optional)* **z-score train-only.** The per-drug mean/std currently use all 180 lines, val/test
      included — mild leakage. Fixing it means computing splits before the targets step.
- [ ] *(Stretch)* cluster cell lines by response and **stratify train/val/test** (high/med/low) for
      lower-variance evaluation.

## DrEval alignment (14.07.2026)

Paper: Bernett, Iversen, Picciani, **Wilhelm**, Baum, List — *Critical evaluation of drug response
prediction models with DrEval*, Nat. Commun. 2026. Half of published models don't beat a naive
drug-mean + cell-line-mean predictor. **Our ridge ≈ MLP finding is the field's norm, independently
reproduced.** Our split *is* their LCO.

Benchmarked with the real package (`drevalpy` 1.5.1) in `notebooks/result_evaluation/dreval_benchmark.ipynb`, and
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
  unreliable ones — quantified in `notebooks/data_and_harmonization/drug_coverage.ipynb` and
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
      give byte-identical bins (verified on 200 cells). Cui et al. present binning as the *replacement*
      for TPM-normalization and log1p, and state that `X` "represent[s] both the raw and preprocessed
      data matrices before binning". Detail in
      [Step 02](./steps/02-preprocessing-and-embeddings.md#compute-environment-and-its-limits-03082026).
- [ ] Regenerate scGPT embeddings — **no longer optional and no longer identical**: `gen_embeds.py` now
      seeds with 42 and runs on MPS, so every embeddings file on disk predates the change. Scope
      (all_genes only vs all five variants) still undecided.
- [ ] *(Optional)* re-run `split_paclitaxel` single-task to fill [Step 04](./steps/04-single-task-results.md)'s
      PCA column, or retire that progression.

## Roadmap (project plan)

- [ ] Cross-database integration — PRISM then GDSC, efficacy + toxicity ([Step 06](./steps/06-planned-work.md#a-cross-database-integration)).
- [ ] XAI / feature interpretability ([Step 06 · B](./steps/06-planned-work.md#b-xai-and-feature-interpretability)).
- [ ] Foundation model + clinical fine-tuning ([Step 06 · C](./steps/06-planned-work.md#c-foundation-model-and-clinical-fine-tuning)).
