# OncoTox — TODO

> # ⛔ 28.07.2026 — THE CURRENT DRUG PANEL IS VOID
>
> **The 8-drug literature panel is discarded.** Its candidate list was ranked by `min(kill, spare)` on
> our own response values before the literature criterion was applied, so the selection inherited the
> label dependency it was built to remove — and 32 approved or clinical compounds, `nutlin-3` among them,
> were excluded for the wrong reason.
>
> **Everything computed on it is therefore provisional:** the step-1 run (`4_training`), the
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
        `notebooks/archive/replicate_variation.ipynb` →
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
- [ ] **5 · Target** — AUC vs EC50 vs Emax: AUC conflates potency with efficacy and the tested
      concentration range spans 0.13–600 µM. Winsorizing threshold. Are all statistics per fold?
      Related open leak: the per-drug target mean/std are computed over every cell line, val and test
      included. **If this audit changes the target, `report/sections/03_methods.tex` §Response target
      changes with it** — it was rewritten 10.08.2026 to state raw `auc`.
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
        `notebooks/analysis/evaluation/dreval_benchmark.ipynb` imports the now-archived
        `dreval_normalize.py` and hardcodes the removed `'auc'` score, so it is broken twice over;
        untouched pending **item 11 (Evaluation)**, which also decides whether the fragility diagnostic
        returns ([why it was archived](../scripts/archive/README.md)).
- [ ] **7 · Splits — walked 12.08.2026.** Confirmed sound: grouping is by cell line everywhere, the
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
        considered and rejected: it is recoverable anyway — `outputs/panel/panel_oof_predictions.csv`
        names all 153 train+val lines, so the test set is the labelled lines it omits — and every number
        scored on it is void on target and panel grounds. R2 creates the file itself; committing it
        there is where the guard starts to protect something (added to R2).
  - [ ] **C · The per-drug-mean null is computed two different ways.** `cv_evaluate` fits its constant
        on the fold's fitting lines (honest); `4_training.ipynb` §4 computes `null_mse` from the
        variance of the **held-out** truth, an oracle constant fitted on the rows it is scored against.
        Conservative — it makes the model look worse — but they are not the same bar and one figure
        cannot be read against the other. Also `_per_drug_train_mean` averages over cells, so lines
        weigh by their cell count, where `density_weighting.line_level` is per line. Routed to
        **item 11 (Evaluation)**, which owns the baselines.
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
- [ ] **8 · Model — walked 12.08.2026.** Confirmed sound: the architecture is what
      [Step 03](./steps/03-model-and-training-design.md#model-architecture--regularization-oncomlppy-25052026)
      says it is, the trunk is genuinely matched between the arms, the "K heads" are the K rows of one
      output `Linear` over a shared trunk with no per-drug sub-network, and inference runs under
      `.eval()` so dropout is off when predictions are made. **The shared trunk is 74,304 parameters at
      every panel size — only the head layer scales with K**, at 65 parameters per drug: 715 of 75,019
      at K=11 (1.0 %), 35,425 of 109,729 at K=545 (32.3 %). So the capacity the heads compete *for* is
      fixed no matter how many of them there are, which is what makes "capacity competition between
      heads" the right description of the K=545 collapse and a bigger panel no cheaper in trunk than a
      smaller one. Against ~153 independent labels. Counts from `arch_facts.py` (audit 08, synthetic
      input, nothing trained).
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
  - [ ] **Found on the way, routed elsewhere.** `--epochs` defaults to 50 in the CLI, 25 in `TrainConfig`
        and `4_training`, and `4_training` §B sets 50 → **item 10**. `dreval_benchmark.ipynb` builds
        `OncoMLP` by hand and so has neither mechanic from A → **item 11**, which owns that notebook.
        `ScGPTDrugDataset` has no consumers — the `train_baseline.py` / `train_scGPT.py` its docstring
        names were deleted in `090f957` — and the `norm="batch"` / `"none"` branches are never exercised
        → **item 13**. LayerNorm makes the forward pass invariant to input rescaling up to the first
        `Linear`'s bias (78× moves the output by 0.074 against a spread of 0.377; zeroing that bias drops
        it to 8e-6), so **4A's premise needs restating before 4A is run** — the naive "different
        effective step size" argument does not go through under LayerNorm plus Adam.
- [ ] **9 · Loss** — masking, per-sample weighting (the density weighting was a null — drop or keep?),
      per-line weighting (the 82× artifact), and what the objective actually rewards.
- [ ] **10 · Training** — optimizer, weight decay groups, epochs, early stopping, and the `mps`
      nondeterminism that moves the PCA arm across identical runs.
- [ ] **11 · Evaluation** — metric, cell→line aggregation, baselines (ridge, `NaiveMeanEffects`),
      dispersion across folds *and* drugs, raw vs normalized.
- [ ] **12 · Reproducibility** — seeds, determinism, what is derived in code versus typed by hand.
      Anything that exists only as a shell command is not a result. *Seeds are now fixed
      (`gen_embeds.py`, `add_pca.py`, the training `DataLoader`s) and determinism was decided
      28.07.2026 as no-action, so what is left is the hand-copying:* **`report/results_numbers.tex` is
      transcribed by hand from CSVs under `notebooks/outputs/` with no extraction script**, so every
      macro can drift silently from its source — and it grew again on 10.08.2026 with four replicate
      macros (`\NRepPairs`, `\NRepLines`, `\RepDiffMed`, `\RepDiffPct`). Related and already fixed:
      `report/main.pdf` was tracked and went stale whenever a `.tex` changed; it is now untracked and
      built from source.
- [ ] **13 · Code redundancy, stale code, notebook restructuring** — not a pipeline stage, so it runs
      last. The **code** half of Selin's "redundanz, staleness, file overload" note; the documentation
      half was done 05.08.2026. Open: duplicated code across scripts and notebooks, code used nowhere,
      restructuring and archiving of outdated notebooks, and whether too many files are produced for
      nothing.

**Working agreement for the session, restated because it was broken today:** no step of the analysis gets
decided silently. If a choice affects what enters the model or how a number is computed, it is proposed
first and agreed before it is executed.

## The sweep — R1 to R6, in order (written down 11.08.2026)

*This sequence existed only in the session task list until now, which meant it lived in a per-session
cache and would not have survived the session. Nothing here may start before the review above is
finished and Selin says so (03.08 banner); R1 is a decision, not a run.*

- [ ] **R1 · Decide the re-embedding scope — which variants get regenerated.** **Selin's decision, and
      it sizes everything below**, because scGPT embedding is the expensive step. Five variants exist on
      disk (`hvg1000/2000/3000/5000`, `all_genes` — `layout.py:31`, `VARIANT_N_TOP_GENES`). What each
      option implies: **all five** keeps [Step 05](./steps/05-multitask-results.md)'s gene-set sweep
      like-for-like; **`hvg5000` + `all_genes`** covers every number the report currently quotes but
      leaves the sweep mixing old and new embeddings; **`hvg5000` only** is cheapest and voids the sweep.
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
- [ ] **R4 · Retrain.** The `DataLoader`s now take an explicit generator and the inputs change under R2,
      so no run under `runs/` is reproducible from current code. Requires the target (item 5) and the
      panel (item 6) settled first, and **never both in one run**. Scope: the 8-run matrix + 5-fold CV
      (`4_training.ipynb` (§B), `train_multitask.py`) — expected to **overturn** the Steps 04–05 numbers, not
      refresh them; `4_training.ipynb` on the rebuilt panel; ridge / `NaiveMeanEffects` on the same
      folds; `verify_variants.ipynb` §9; and 4A below. Blockers already recorded: **≥ 3 seeds** before any
      scGPT − PCA margin is quoted, and train-only drug selection inside each fold.
  - [ ] **The minimal capacity re-derivation** (item 8C, 12.08.2026): trunk `(128,64)` vs a bare linear
        head (`hidden_dims=()`), both representations, against `RidgeCV` on the *same* folds via
        `cv.grouped_folds`, on the rebuilt panel. ~20 fits. It re-establishes the one claim that carries
        weight — *scGPT needs the nonlinearity and PCA does not* — whose current evidence is void. Do not
        widen it into the four-knob sweep; that question is not what the load-bearing claim rests on.
  - [ ] **Two settings change here that no run on record used** (item 8A): every training path now
        initializes head biases at the fitting-fold per-drug means, and weight decay skips the biases and
        the LayerNorm parameters while now applying to the output weight matrix, which
        `exclude_output_from_decay` had exempted. So the panel run's configuration also changed, not only
        the matrix's — `4_training` §6 is not a like-for-like predecessor of its own next execution.
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
      `analysis/evaluation/`: `diagnostics` (§5 dispersion), `dreval_benchmark` — the latter is
      **broken twice over** (imports the archived `dreval_normalize`, hardcodes the removed `'auc'`)
      and belongs to review item 11 before it can run at all.
      `analysis/harmonization/drug_coverage` — **not optional**, the line count moves 180 → 181.
      `2_drug_selection` does **not** need the retrained outputs — it reads only the
      response CSV — but re-run it anyway once preprocessing has, so the panel is regenerated against
      the same artifacts as everything else and `panel.csv` cannot silently predate them.
      **Two notebooks still have their stored outputs cleared** (11.08.2026), because their score
      literal changed to `auc_cc` and the old results could not be refreshed under the freeze — the
      code is correct and only the results are missing: `4_training` (c1) and
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
      **`03_methods.tex` §Representation and model describes the architecture but not the training
      configuration** — no optimizer, learning rate, weight-decay grouping, early stopping or head-bias
      initialization appears anywhere in the report, so a reader cannot reconstruct the run. Write it
      here from `TrainConfig` and [Step 03](./steps/03-model-and-training-design.md#the-uncentred-target-is-handled-the-same-way-in-every-training-path)
      (noted 12.08.2026, audit 08).
      **Define the drug-panel counts as macros** rather than the inline digits currently in
      `03_methods.tex` §Drug panel selection (150 / 120 / 57 / 44 / 11, the 90 % cut, 91.2–98.3 %
      coverage, 102 and 15 unmatched compounds, 13 parent-CID matches) — via item 12's extraction
      script, since adding them by hand is the defect that item exists to remove.
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

### Step 1 — run 27.07.2026 (`notebooks/4_training.ipynb`), open items only

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
- **⚠ Open decision — what counts as a positive Q2 result, fixed before the run.** Under the Q1 framing
  the readout was ρ against the controls; under this one it is not, and nothing has replaced it. Whatever
  is chosen needs a control that says what *no* heterogeneity signal looks like (a shuffled-cell or
  uniform-attention null), or a structured-looking attention map will be read as a positive result by
  default. **Selin's call, and it has to be made before the run rather than after seeing one.**
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
- [ ] **Add ridge (line-level) to `4_training` §B's comparison tables** so every future claim is scored against it.
- [ ] *(Optional)* **z-score train-only.** The per-drug mean/std currently use all 180 lines, val/test
      included — mild leakage. Fixing it means computing splits before the targets step.
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
