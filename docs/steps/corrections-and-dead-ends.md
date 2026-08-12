# Corrections, superseded results & dead ends

*Part of [OncoTox project progress](../project_progress.md). The record of what did **not** survive:
results superseded by better measurement, claims retracted, hypotheses tested and refuted, approaches
abandoned, and working-method failures worth not repeating.*

**Steps 01–05 state only what currently holds.** Where an earlier version said something else, the step
carries a one-line pointer here. So a number in Steps 01–05 is current unless a pointer says otherwise,
and nothing in this file is a live result.

**Why keep any of it.** A refuted hypothesis is a result — it rules something out, and the ruling-out is
usually what justifies the next step. A retracted claim that is quietly deleted leaves the docs looking as
though they were always right, which makes the same error easy to repeat and the record impossible to
audit. Under FAIRER this is **R** and **E** at once: the reasoning has to be reconstructible, and a
correction that is silently overwritten is a reporting problem rather than a bookkeeping one.

**Each entry carries:** what was believed, when it was established and when it was overturned, what
overturned it (notebook, script or commit — re-runnable, not a chat message), what it affected, and what
replaced it. Where nothing has replaced it yet, that is stated.

---

## Index

Grouped by what each entry *is*, because "improved on" and "was wrong" are not the same thing.

**Refuted hypotheses — tested deliberately, answer was no.** Each rules something out.

| Date | What |
|---|---|
| 27.07.2026 | [Inverse-density loss weighting improves ranking](#inverse-density-loss-weighting-improves-ranking) |
| 27.07.2026 | [The cell-line effect is largely proliferation](#the-cell-line-effect-is-largely-proliferation) |
| 13.07.2026 | [The model is over-regularized, or too small](#the-model-is-over-regularized-or-too-small) |
| 13.07.2026 | [Line-balanced reweighting will help](#line-balanced-reweighting-will-help) |

**Superseded by better measurement or a closed confound.** The work improved; the older numbers are kept
so the improvement is auditable.

| Date | What |
|---|---|
| 27.07.2026 | [`auc_z` as the training target](#auc_z-as-the-training-target) — decomposed into what it actually did |
| 27.07.2026 | [Per-drug variance weighting (`σ_noise`, plan items 1.0–1.2)](#per-drug-variance-weighting--dissolved-by-a-scope-change-never-needed) — dissolved by a scope change |
| 14.07.2026 | [The 13.07 five-drug numbers](#the-1307-five-drug-numbers) — re-run on a wider set; neither column was ever a generalization number |
| 13.07.2026 | ["Neither representation ranks cell lines" — the K=545 null](#neither-representation-ranks-cell-lines--the-k545-null-result), and the "ceiling is the label" net read it supported |
| 13.07.2026 · 05.08.2026 | [The 8-run matrix conclusions](#the-8-run-matrix-conclusions) — design stands, conclusions do not; the `all_genes` axis struck separately (scGPT saw a capped random draw, PCA the whole gene set) |
| 13.07.2026 | [The Steps 04–05 numbers as a comparable baseline](#the-steps-0405-numbers-as-a-comparable-baseline) — scale change, not error |
| 27.06.2026 | [The ~50-d PCA baseline and the `(64,32)` trunk](#the-50-d-pca-baseline-and-the-6432-pca-trunk) — two confounds found and closed |

**Genuine errors.** These cost time or invalidated results.

| Date | What |
|---|---|
| 12.08.2026 | [The same val-split leak, in the code that produced everything else](#the-same-val-split-leak-in-the-code-that-produced-everything-else) — fixed once in the DrEval notebook 14.07.2026 and never carried into `cv.py`; **every out-of-fold prediction and CV metric on record is a minimum over epochs on its own scored data** |
| 12.08.2026 | [The evidence that closed model-side tuning](#the-evidence-that-closed-model-side-tuning) — the same leak, in the ablation notebook, on a retired target and a discredited drug set; only the MLP rows are flattered and ridge ties them |
| 12.08.2026 | [The two uncentred-target mechanics ran in one training path of three](#the-two-uncentred-target-mechanics-ran-in-one-training-path-of-three) — the 8-run matrix trained against an offset the panel run did not, and the flag did not do what the docs said |
| 11.08.2026 | [The `auc` target was divided by the wrong quantity](#the-auc-target-was-divided-by-the-wrong-quantity) — `conc_pts_fit` counts points surviving outlier censoring, not the width of the integral; **every `auc` number is void**, and the target source was replaced rather than the divisor fixed |
| 05.08.2026 | [scGPT discarded genes that are in its vocabulary under their current symbols](#scgpt-discarded-genes-that-are-in-its-vocabulary-under-their-current-symbols) — exact symbol match against an older annotation; **repaired in code 05.08.2026, takes effect at the sweep** |
| 28.07.2026 | [The 8-drug literature panel, and every number computed on it](#the-8-drug-literature-panel-and-every-number-computed-on-it) — drawn from a pre-filtered pool |
| 27.07.2026 | [The kill/spare learnability gate measured potency, not rankability](#the-learnability-gate-measured-potency-not-rankability) — survived months |
| 14.07.2026 | [The first DrEval benchmark — a val-split leak](#the-first-dreval-benchmark--a-val-split-leak) — found and fixed same day (`ee07b00`) |

**Retracted claims.** Written into the docs, then withdrawn — most within a day.

| Date | What |
|---|---|
| 10.08.2026 | [The scGPT binning invariance was "verified on 200 cells"](#the-scgpt-binning-invariance-was-verified-on-200-cells) — no code behind it, and the raw counts it compared against do not exist here; the conclusion holds as an argument |
| 27.07.2026 | ["scGPT clears the ridge control" was a first](#scgpt-clears-the-ridge-control-was-a-first) — it was a replication |
| 25.07.2026 | [The panel was chosen blind to our labels](#the-panel-was-chosen-blind-to-our-labels) |
| 14.07.2026 | [`ml210` was rejected on coverage](#ml210-was-rejected-on-coverage) |
| 14.07.2026 | [The learnability filter was validated against achieved ρ](#the-learnability-filter-was-validated-against-the-ρ-the-model-achieves) — unreproducible source |
| 13.07.2026 | [The curve fit preserves signal the dose-average destroys](#the-curve-fit-preserves-signal-the-dose-average-destroys) |
| 13.07.2026 | [The prediction shrinkage is a defect to fix](#the-prediction-shrinkage-is-a-defect-to-fix-with-lighter-regularization) |
| 28.06.2026 | [PCA prefers the full transcriptome](#pca-prefers-the-full-transcriptome) — refutation holds for PCA; narrowed 05.08.2026 for scGPT |

**Dead ends.** Checked, then not pursued.

| What |
|---|
| [scDrugAtlas and ClinTox as data sources](#scdrugatlas-and-clintox-as-data-sources) |
| [Migrating the cell-line join to persistent identifiers](#migrating-the-cell-line-join-to-persistent-identifiers) |
| [Kinker's two named associations do not transfer to this task](#kinkers-two-named-associations-do-not-transfer-to-this-task) |
| [`kx2-391` carries drug-specific signal](#kx2-391-carries-drug-specific-signal) |
| [Considered and never pursued](#considered-and-never-pursued) — a bespoke transformer/VAE, DeepInsight |
| [Retired code paths](#retired-code-paths) |

**Also here:** [the Step 1 run on the voided panel](#the-step-1-training-run-on-the-voided-panel) and
[Process failures](#process-failures) — working-method problems rather than results.

---

## Superseded results

### scGPT discarded genes that are in its vocabulary under their current symbols

**Found** 05.08.2026 while checking [TODO](../TODO.md)'s standing question *"scGPT OOV check — sind die
wirklich OOV?"*. The answer is **no**.

**What it is.** `gen_embeds.py` resolves genes by an **exact string match** of SCP542's symbols against
`scGPT_human/vocab.json`; `embed_data` then drops everything unmatched (`scgpt/tasks/cell_emb.py:220`).
But SCP542 carries an **older HGNC annotation than the vocabulary**, so every gene renamed between the
two releases fails the match and is thrown away — while being present in the vocabulary under its
current symbol. The mismatch is not a formatting problem: of the 424 / 2,152 dropped symbols, **0** are
recovered by case-folding and **5** by stripping punctuation.

The six highest-expression casualties in `all_genes`, with the vocabulary membership of each symbol
checked directly:

| dropped as OOV | expressed in | current symbol | in vocab |
|---|---|---|---|
| `GNB2L1` | 99.96 % of cells | `RACK1` | ✓ |
| `ATP5E` | 99.95 % | `ATP5F1E` | ✓ |
| `H2AFZ` | 99.47 % | `H2AZ1` | ✓ |
| `H3F3B` | 99.82 % | `H3-3B` | ✓ |
| `HIST1H4C` | 92.30 % | `H4C3` | ✓ |
| `H3F3A` | 99.65 % | `H3-3A` | ✓ |

Also `PTRF`→`CAVIN1`, `CYR61`→`CCN1`, `LINC00152`→`CYTOR`, `H2AFJ`→`H2AJ`, `C6orf48`→`SNHG32`. Every
old symbol absent from the vocabulary, every current one present.

**Scale**, counted 05.08.2026 in `notebooks/analysis/qc/gene_symbol_rescue.ipynb` against
`reference/hgnc_complete_set.txt` (HGNC approved set, published 04.08.2026, pinned by SHA-256 because
HGNC overwrites its download URL in place and publishes no dated archive —
[provenance](../../reference/README.md)). Percentages are of the **whole transcriptome**: the source
matrix is CPM-normalized, so each cell carries a fixed 10⁶ budget and a gene set's CPM-per-cell is its
share of it.

| | genes discarded | **recoverable renames** | discarded expression | **of which recoverable** |
|---|---|---|---|---|
| `hvg5000` | 424 | **129** | 0.433 % | **0.378 %** — 87 % of the loss |
| `all_genes` | 2,152 | **775** | 3.874 % | **3.608 %** — 93 % of the loss |

**Nearly the whole loss is recoverable.** In `all_genes`, 775 genes carrying **3.6 % of every cell's
transcriptome** were discarded for no reason, 392 of them expressed in over 10 % of cells. Genuinely
absent genes — clone-based names that were never HGNC entities — number 1,348 but carry only 0.21 %.
A further 26 genes were renamed into a symbol the vocabulary also lacks (0.013 %), and 3 former symbols
map to several current genes and were left unresolved rather than guessed.

Recovery is by HGNC `prev_symbol` only; `alias_symbol` is deliberately not used, since a synonym can map
two distinct genes onto one name. The counts are therefore a **lower bound**.

> ⚠️ **Correction, 05.08.2026 (same day).** This entry first gave the discarded share as **1.48 %** for
> `hvg5000` and 3.90 % for `all_genes` under one column heading. Those had different denominators: the
> `all_genes` figure was of the transcriptome (3.87 % exact, confirmed), but the `hvg5000` figure was of
> the **HVG-5000 subset's** expression. Against the transcriptome it is **0.433 %**. The table above uses
> one denominator throughout.

**Collisions — 1 in `hvg5000`, 11 in `all_genes`.** A recovered symbol that already exists as its own row
cannot simply be remapped; the two rows would have to be merged, which combines expression values and is
an analysis decision about the data. They are listed by the notebook and **left unresolved**. The largest
shows why: `HNRNPU-AS1 → HNRNPU` would fold an antisense lncRNA into its sense gene. Also
`C2orf48 → RRM2`, `C10orf12 → LCOR`, `CTAGE5 → MIA2`.

**What it affects — scGPT only, and worse than the gene count suggests.** PCA reads the `convert` file,
before the OOV drop, so it receives all 5,000 / 22,722 genes; HVG selection also runs before the drop.
The damage is confined to the scGPT arm. Within it, the loss is **not** limited to the discarded genes:
`DataCollator(do_binning=True)` bins each cell by quantiles of that cell's own non-zero values, so
removing the most-expressed genes moves every quantile edge and changes the bin assigned to *every
remaining gene*.

**Direction of the bias.** It handicaps the representation the project argues for. Every reported
scGPT-vs-PCA margin is therefore **conservative** with respect to this defect — it can only have
understated scGPT, never inflated it.

**Repair decided 05.08.2026 (Selin); applied to the code the same day, `scripts/annotation/gene_symbols.py`
plus the two callers below. Nothing is re-run** — see *Standing* at the end of this entry. Three choices
had to be settled; all were taken so that **no expression value anywhere changes**, which keeps the
eventual re-embed attributable to the recovered genes alone.

**1 — The rename is recorded as a new `var['hgnc_symbol']` column, not by rewriting `var_names`.**
`scp542_conversion.py` adds the column, `gen_embeds.py` resolves the vocabulary through it. SCP542's
identifiers stay exactly as distributed, so the matrix remains comparable to Kinker et al.'s published
one, while the modern symbol becomes available as a join key for everything downstream — the direction
[review item 2B](../TODO.md) asks for. `var` is currently empty, so the column is purely additive.
*Rejected:* rewriting `var_names` (destroys provenance for nothing the column does not give) and
resolving only inside `gen_embeds.py` (leaves stale symbols in the data, so every later join — CTRPv2
gene sets, XAI naming, PRISM/GDSC integration — meets the same problem again).

**2 — The 12 collisions are left unmapped.** They stay out-of-vocabulary exactly as today. They carry
**21.4 CPM/cell = 0.0021 %** of the transcriptome, against the rescue's 3.608 % — **0.06 % of what the
rescue is worth**. *Rejected:* merging the old row into the existing one, which would recover that
0.0021 % at the price of changing expression values and risking a biologically wrong merge
(`HNRNPU-AS1 → HNRNPU` folds an antisense lncRNA into its sense gene); and dropping the old row, which
loses the same expression while still perturbing the gene set. Recovering them case by case remains
possible **later, as a separately attributable change**.

**3 — A rename is refused where the old symbol is itself a current approved symbol.** *Added
05.08.2026 while implementing the repair, after the guard was found missing from the map the notebook
built.* HGNC's `prev_symbol` records **reassignments** as well as renames: `OSR1` is the approved symbol
of odd-skipped related 1 *and* a former symbol of `OXSR1`; likewise `NTNG1`→`NTNG2`, `ADCY3`→`ADCY8`,
`SRGAP2`→`SRGAP3` (verified directly in `reference/hgnc_complete_set.txt`). Where SCP542 uses a symbol
that is approved today, it is read as meaning the gene that holds it today. The assumption — that
SCP542's annotation postdates each reassignment — cannot be checked, since SCP542 ships no annotation
version, which is the reason for resolving the doubt toward leaving the gene alone. This refuses 228
further pairs (314 refused in total, with the 86 ambiguous ones), and costs **two** rescues in
`all_genes` — `RNU12`, which would otherwise be fed to scGPT under the pseudogene token `RNU12-2P`, and
`EPB41L4A-AS2` — worth 0.003 % of the rescued expression, and **none** in `hvg5000`.

**Consequence for the numbers above:** the recoverable counts become **773** (`all_genes`) and **129**
(`hvg5000`), which round to the same 3.6 % / 0.378 % quoted in the table and in the report. The table is
left as measured; `gene_symbol_rescue.csv` predates the guard.

Because renaming is numerically inert for the value-based steps — PCA works on values, HVG ranks by
dispersion, neither reads gene names — these choices leave `X_pca` and the HVG selection
**bit-identical** and move only scGPT's token set.

**How it resolves, and why a row's own symbol wins.** `gen_embeds.py::resolve_gene_names` consults
`hgnc_symbol` **only where the row's own symbol is not a token**. Renaming unconditionally would *lose*
**336** genes in `all_genes` and **53** in `hvg5000`: genes the vocabulary holds under the symbol SCP542
uses, which HGNC has since reassigned to a gene the vocabulary does *not* hold (`TP73-AS1`→`GFOD3P`,
`HSPB11`→`IFT25`, `C1orf127`→`CIROZ`). Preferring the row's own symbol makes the repair strictly
additive — no gene kept today is dropped by it — which is what keeps the re-embed attributable.

The annotation runs **before** HVG selection, so the collision check sees the whole transcriptome and a
symbol is never reused for a gene that HVG happened to drop. Measured cost against checking each
variant's own gene set: one recovered gene in `hvg5000`, none in the smaller variants.

**Net effect on the gene set handed to scGPT**, measured 05.08.2026 on the matrices as they stand
(read-only, `backed='r'`; nothing regenerated):

| variant | embedded today | after the repair | recovered |
|---|---|---|---|
| `hvg5000` | 4,576 | **4,704** | +128 |
| `all_genes` | 20,570 | **21,332** | +762 |

The recovered counts are below the 773 / 129 *recoverable* figures because collisions stay unmapped
(11 and 1 of them respectively) — the defect's size and the repair's yield are different quantities.

**Standing — code only, nothing recomputed.** Every embedding on disk, and therefore every scGPT number
in the docs and the report, still carries the defect in full. The repair takes effect at the
[clean sweep](../TODO.md); it landed before it so that the sweep does not have to be paid for twice.

### The learnability gate measured potency, not rankability

**Established** 13.07.2026 (`notebooks/archive/learnability_filter.ipynb`) — a drug was kept only if it both
**killed** a real population of lines (raw `auc ≤ 0.5`) and **spared** one (`auc ≥ 0.8`), on top of
coverage ≥ 90 %. The differential-response condition was the part the loose `drug_coverage` gates lacked, and it
was what took 545 drugs down to 10.

**Overturned** 27.07.2026 (`notebooks/analysis/evaluation/diagnostics.ipynb`; the question that surfaced it was why
`nutlin-3` was not in the panel). `auc ≤ 0.5` asks *does the line die* — absolute potency, which is
essentially the per-drug mean. But the target subtracts that mean and the metric is Spearman, which
reads only the ordering of lines. **The gate selected on the one quantity the model is neither given
nor scored on.**

| drug | `auc_mean` | `auc_std` | kill (`≤ 0.5`) | gate verdict |
|---|---|---|---|---|
| `dasatinib` | 0.631 | **0.155** | 35 | selected |
| `nutlin-3` | 0.874 | **0.147** | **0** | rejected |

Nutlin-3's lines differ as much as dasatinib's; the distribution simply sits higher, because nutlin-3
is **cytostatic** — p53 drives arrest and senescence, so viability never crosses 50 % however sensitive
the line. Any "does it kill" threshold is structurally blind to every cytostatic agent, which is a large
share of targeted therapy. This is not one unlucky drug: **116 of 545** have zero kills yet
`auc_std ≥ 0.10` and coverage ≥ 90 %, `oxaliplatin` among them, all silently discarded.

The cost is sharpened by which drug it was: `nutlin-3`/TP53 is the strongest association in the GDSC
pharmacogenomic screen, and it is *expression*-readable — `MDM2`, `CDKN1A`, `RPL22L1` appear in
~90–100 % of published gene sets predicting nutlin-3a sensitivity. Close to a best case for this model,
and the filter threw it out.

**What it affected.** Everything selected through it: the 10-drug panel and all K=10 numbers, the 8-drug
literature panel drawn from a pool this gate had pre-filtered (below), and the 5-drug best-case
diagnostic. Not the K=545 results, which apply no gate.

**Replaced by** — the [literature panel of 11](01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026)
(12.08.2026), which selects on **FDA approval and published determinants**, with no statistic of our own
labels in the criterion at all.

⚠️ **The replacement first agreed here — spread (`auc_std`) plus coverage, no kill counts — was itself
abandoned on 12.08.2026 (Selin), and this paragraph used to state it as the plan.** Two reasons. Spread
is still a property of *our* response values, so selecting on it keeps the selection label-dependent and
optimistic in exactly the way the entry above objects to, only less crudely. And the second half of the
rationale — that high `auc_std` is "a safe z-score denominator" — belonged to `auc_z`, which was
[retired](#auc_z-as-the-training-target); `uns["ctrp_score_scale"]` no longer exists, so that sentence
described machinery the tree does not contain.

### The 8-drug literature panel, and every number computed on it

**Established** 25.07.2026. To escape the gate above, the panel was anchored in *published* cell-line
sensitivity determinants instead of our own label statistics: CTRPv2 rows of
`data/drug/all_sources_drug_catalog.csv` restricted to single agents with
`compound_status ∈ {FDA, clinical}` (173 of 545), keeping those with an independently published
determinant. The eight: `methotrexate` (SLC19A1), `dasatinib` (six-gene signature / LYN), `paclitaxel`
and `vincristine` (ABCB1 / TUBB3), `afatinib` (EGFR+ERBB2 amplification), `topotecan` (SLFN11),
`tanespimycin` (NQO1), `selumetinib` (BRAF/RAS). Citations are in
[Step 05](05-multitask-results.md).

**Overturned** 28.07.2026. The candidate list was produced by ranking those 173 compounds by
`min(kill, spare)` — computed on our own `auc` values over all 180 lines, val and test included — and
the literature criterion was applied only to that ranked list. So the panel inherited the exact label
dependency it was built to remove, filtered by the exact quantity the entry above discredits.

Two consequences, and the second is the reason the panel is void rather than merely caveated:

- Compounds with a published determinant dropped out on **our** label statistics, not on the
  literature: `sirolimus` (6 kill / 63 spare), `neratinib` (12 / 75), `clofarabine` (15 / 68),
  `cytarabine hydrochloride` (5 / 109), `gdc-0941` (5 / 45).
- **32 of the 116 wrongly-discarded drugs are approved or in clinical trials** — `oxaliplatin`,
  `bortezomib`, `ruxolitinib`, `regorafenib`, `entinostat`, and **`nutlin-3` itself**, the compound used
  to demonstrate the defect (spread 0.147, coverage 0.96, status `clinical`, zero kills, so balance 0
  and never a candidate).

The eight compounds remain defensible *as compounds*; what was silently pre-filtered is the **pool they
were drawn from**. The honest description of the panel as executed: literature-anchored,
spread-verified, and **drawn from a kill-filtered pool**.

**What it affected.** The Step 1 training run (`notebooks/4_training.ipynb`), the distribution
and weighting design (`notebooks/archive/panel_distributions.ipynb`), the dispersion figures
(`notebooks/analysis/evaluation/diagnostics.ipynb` §5), the panel rows in [Step 05](05-multitask-results.md), and the
corresponding numbers in `report/`. The *methodological* findings from that run survive — the collapse
was a head-count effect, density weighting is a null, the ridge tie replicates — because none of them
depends on which eight drugs were chosen. The *numbers* must be re-derived.

**Replaced by** — the [literature panel of 11](01-datasets-and-harmonization.md#the-drug-panel--fda-approved-compounds-this-screen-covers-12082026),
built 12.08.2026 and label-blind by construction: FDA approval (Sun et al. 2017, Table 1) and a
verified published claim decide membership, and the only property of our own data that enters is
**coverage** — how many of the 181 overlapping lines were screened at all, which is a statement about
the experiment rather than about the response.

Two things this entry recommended were **not** done, and the reasons are worth keeping:

- *"Rebuild the pool on coverage and `auc_std`"* — dropped, because spread is still our own label
  statistic ([above](#the-learnability-gate-measured-potency-not-rankability)).
- *"Measure the spread requirement on GDSC2 or PRISM instead"* — unnecessary once the criterion stopped
  using spread at all. It would have made the selection blind to *CTRP's* labels while still selecting
  on somebody's response values; the published-determinant criterion is blind to all of them.

**Three of the eight re-enter on their own evidence** — `paclitaxel`, `dasatinib` and `afatinib`, each
a named worked example in Seashore-Ludlow 2015. The other five fall out in three different ways, which
is a fair test of the new criterion:

- `methotrexate`, `vincristine`, `topotecan` — **candidates**, clearing FDA approval and the coverage
  cut (0.956 / 0.989 / 0.978), and each carries a determinant in the general literature. They are not in
  the panel because none appears in the four papers this project is built on, and `doxorubicin` is the
  only compound admitted on an outside citation alone, to represent the anthracycline class.
- `tanespimycin` — **not a candidate: never FDA-approved.** Its development was discontinued. The old
  panel carried it on an NQO1 bioactivation determinant, which is real biology and not an approval.
- `selumetinib` — **not a candidate under this list**, approved in 2020 and therefore past Sun et al.'s
  2014 cut-off. This one is a limitation of the source rather than a judgement on the compound, and it
  is the concrete cost of using a list that stops at 2014.

The 15.07 progress report was postponed rather than presented on this panel.

**The reasoning that led here, kept so the rebuild does not repeat it.** The panel replaced the 10-drug
filtered set for two reasons, and only the second was fully solved:

- The `learnability_filter` gates were computed on all 180 lines, val and test included, so selection saw held-out labels.
- The gates were *unstable*: shifting the kill/spare thresholds from 0.5/0.8 to 0.7/0.8 yields a
  completely different ten drugs of the same quality. The filter enriched reliably, but *which* drugs it
  named was arbitrary — and that arbitrariness was invisible until it was written down.

Against those, the panel as executed **fixed** the arbitrariness of which drugs (they are named by
citation, reproducible by someone who never sees our AUCs) and the threshold instability. It **did not
fix** the optimistic component in per-drug ρ, because the panel is still enriched for drugs that happen
to separate *our* 180 lines. The honest description at the time was "literature-anchored,
spread-verified" — and after 27.07, "drawn from a kill-filtered pool" as well. Train-only selection
therefore stayed blocking throughout and still is.

Two further notes on the mechanics. All eight passed the `learnability_filter` gate unchanged, so the panel was a
**re-ranking inside the gate-passing set rather than a relaxation of it** — which is precisely why the
gate defect propagated into it. And only `dasatinib` and `methotrexate` overlapped the old ten, so six of
the eight were compounds the previous filter never named.

The selection also **promoted `data/drug/all_sources_drug_catalog.csv` from "exploratory, consumed by no
model" to a selection input.** The catalog is built in `notebooks/analysis/harmonization/drug_catalog.ipynb` from CTRP's
official `v20.meta.per_compound.txt`, mapping `gene_symbol_of_protein_target` → `target`,
`target_or_activity_of_compound` → `moa_or_pathway`, `cpd_status` → `compound_status`
([Step 01](01-datasets-and-harmonization.md)). The rebuild will use it the same way, so that promotion
stands even though the panel does not.

### The Step 1 training run on the voided panel

**Run** 27.07.2026, `notebooks/4_training.ipynb`. The first execution of the retired-`auc_z`
setup: target raw `auc` winsorized at 1.1, the 8-drug literature panel, per-sample inverse-density
weights fitted per fold on training lines only, output layer excluded from weight decay, head biases
initialized to the train-fold per-drug means. Architecture, splits, optimizer and batching unchanged, so
the change was attributable. 5-fold GroupKFold over the 153 train+val lines, **one seed (42)**.

**Provisional as of 28.07.2026** — computed on the
[voided panel](#the-8-drug-literature-panel-and-every-number-computed-on-it), so every number below must
be re-derived on the rebuilt panel.

| model | ρ `X_pca` | ρ `X_scGPT` | MSE `X_pca` | MSE `X_scGPT` |
|---|---|---|---|---|
| MLP, unweighted | 0.316 ± 0.003 | **0.377** | 0.0265 | 0.0254 |
| MLP, density-weighted | 0.308 | 0.369 | 0.0274 | 0.0254 |
| `RidgeCV` on line means | 0.306 | 0.299 | 0.0270 | 0.0268 |

Null (per-drug mean) MSE is 0.030, so these read directly: the scGPT model explains ~15 % of the variance
**in AUC units**, RMSE ≈ 0.16 viability.

**Dispersion**, computed from the stored out-of-fold predictions without retraining
(`notebooks/analysis/evaluation/diagnostics.ipynb` §5, `outputs/diagnostics/result_dispersion.csv`):

| | pooled ρ | sd across the 5 folds | sd across the 8 drugs | per-drug range |
|---|---|---|---|---|
| PCA, unweighted | 0.315 | ±0.028 | 0.111 | 0.19 – 0.53 |
| scGPT, unweighted | 0.377 | ±0.043 | 0.091 | 0.30 – 0.55 |

The **scGPT−PCA gap of +0.062 is about one fold standard deviation**, so it was consistent evidence and
never an established margin — on top of being a single seed. (The reporting rule this run established —
pooled value as the point estimate, fold spread as the dispersion, never fused into one `mean ± sd` — is
live and kept in [Step 05](05-multitask-results.md).)

**Absolute numbers were lower than the 14.07 panel (0.356 / 0.402), which was the expected direction:**
those drugs were selected using all 180 lines including val/test, and these were not chosen for spread on
our labels at all. Lower and more defensible was the deliberate trade.

**Most of the signal survived removing the cell-line effect.** Computed with **zero model fits** from the
stored out-of-fold predictions → `outputs/dreval/dreval_normalized_panel.csv`.

> ⛔ **The code that produced this table no longer exists (12.08.2026).** It was
> `scripts/evaluation/dreval_normalize.py`, whose cell-line-effect removal was a **local invention with
> no counterpart in DrEval's paper** and was deleted rather than archived, on the grounds that a
> home-grown metric inside a file named after DrEval gets read as DrEval's
> ([why](../../scripts/archive/README.md)). A script of that name is still live but computes the
> paper's normalization only, so **it will not reproduce the `normalized ρ` or `fragility alone`
> columns below**. Recoverable at `git show bf93084:scripts/evaluation/dreval_normalize.py`; whether
> the diagnostic returns is review item 11.

| | raw ρ | normalized ρ | fragility alone |
|---|---|---|---|
| scGPT | 0.377 | **0.329** | 0.491 |
| PCA | 0.315 | 0.304 | 0.491 |

Removing the cell-line effect cost scGPT 0.048 and PCA 0.011, so the bulk of the per-drug correlation was
**drug-specific** rather than "this line is fragile". Per drug, the weakest were `topotecan`
(0.296 → 0.192) and `vincristine` (0.407 → 0.280), while `dasatinib` (0.546 → 0.558) and `afatinib`
(0.325 → 0.368) *gained* — their signal is orthogonal to fragility. Nothing on this panel behaved like
`kx2-391` on the earlier one, which [collapsed to 0.006](#kx2-391-carries-drug-specific-signal).

Note the fragility baseline itself scores **0.491, higher than either model.** It is **not** a legitimate
predictor — it reads held-out labels — so it is a diagnostic ceiling rather than a competitor, but it is a
useful statement of how much of the raw ranking is line-level rather than drug-level.

**What the run settled, and does not depend on which eight drugs:** the June collapse was a head-count
effect ([above](#auc_z-as-the-training-target)); inverse-density weighting is a clean null
([below](#inverse-density-loss-weighting-improves-ranking)); the PCA ridge tie and the scGPT margin over
ridge both replicate ([retraction](#scgpt-clears-the-ridge-control-was-a-first)); and most of the signal
is drug-specific rather than cell-line fragility. **What it did not touch:** research question 2, which
remains structurally untestable under a constant-within-line label.

### The `auc` target was divided by the wrong quantity

**Overturned 11.08.2026**, during audit item 5 (Target). **Every number computed on the `auc` target is
void** — Steps 04 and 05, the report, and every figure derived from them. This was not repaired: the
target source was replaced instead (below), so the defective code path leaves the tree rather than
getting a corrected divisor.

**What the code did.** `scripts/preprocessing/ctrp_to_h5ad.py` built `obsm["Y_ctrp"]` as

```python
curves["score"] = curves["area_under_curve"] / curves["conc_pts_fit"]
```

The intent, stated in the comment above it, was correct: CTRP's `area_under_curve` is an *integrated*
area, so it grows with the extent of the concentration range and must be reduced to an average before
compounds sit on a common axis.

**Why it is wrong: `conc_pts_fit` is not the extent of the range.** CTRP's own column documentation
(`v20._COLUMNS.txt`, 17.12.2015, distributed with `CTRPv2.0_2015_ctd2_ExpandedDataset`) defines it as
*"number of concentration-response points participating in curve fit"*. The CTRPv2 publications say
where that number comes from — the design is fixed, and the count varies only because points are
discarded during fitting:

> "The effects of small molecules were measured over a **16-point concentration range (2-fold
> dilution)** in duplicate." … "we used **Cook's distance during the curve-fit procedure to censor
> outlier points**, but did not allow the highest two concentration points to be censored this way." …
> "average 15.6 concentration points each"
>
> — Seashore-Ludlow et al., *Cancer Discovery* **5**(11):1210–1223 (2015), Methods. The same 16-point
> design and "integration under the 16-point curve" in Rees et al., *Nature Chemical Biology*
> **12**:109–116 (2016), Methods.

So the divisor counts **the points that survived outlier censoring**, while the integral's width is the
concentration range — which censoring of interior points does not change. Dividing a fixed-width
integral by a quality-dependent count inflates the target for precisely those cell lines whose
measurements were noisy enough to lose points: a curve reduced from 16 points to 13 is divided by 13
instead of 16, so its value rises by 23 % with no change in biology. **Noisier measurements read as more
drug-resistant.**

That lands on the reported metric directly. Ranking cell lines *within* a compound is what the Spearman
figures in Steps 04–05 measure, and the inflation varies between the cell lines of a compound rather
than shifting the compound as a whole.

**On the size of it.** The magnitude was measured on 11.08.2026 while the defect was found — how many
curves lost points, and how far the corrected ranking departs from the current one. **Those checks were
one-off scripts and are not in the repository, so no figure from them is quoted here.** The mechanism
above is sourced and sufficient to void the numbers; the replacement target makes re-measuring the old
one pointless.

**What replaced it.** The target source moved to DrEval's reprocessed CTRPv2 (decided 11.08.2026,
Selin), pinned to Zenodo record `21807175`. See [Retired code paths](#retired-code-paths).

### `auc_z` as the training target

**Established** 13.07.2026 (`notebooks/archive/target_comparison.ipynb`) — per-drug z-scored AUC,
`auc_z[l,d] = (auc[l,d] − mean_d) / std_d`, computed by `_zscore_per_drug` in
`scripts/preprocessing/ctrp_to_h5ad.py`. It was adopted because raw `auc` collapsed at 545 heads
(scGPT −0.087, PCA +0.016) while `auc_z` held (+0.430 / +0.378) on the same drugs, model and split —
at the time the single largest improvement the project had made.

**Overturned** 27.07.2026, by decomposing it into the two transforms it actually applies:

- **The centering is inert.** Each drug is its own output row with its own bias term, so the head
  absorbs the per-drug mean either way — and per-drug Spearman is shift-invariant, so the metric never
  saw it.
- **The scaling is the defect.** Dividing by `auc_std` forces *every* drug to variance 1, including
  drugs whose true spread sits at the assay noise floor, whose noise is then amplified to variance 1 and
  enters the shared loss at full weight. This is the mirror image of the June σ² bug: that one
  over-weighted wide drugs, `auc_z` over-weights narrow ones. The correct weight is signal-to-noise, not
  equal variance.

Both are within-drug monotone transforms, so neither was ever visible to the metric — **`auc_z` was a
loss-weighting scheme in disguise**. And the collapse it was introduced to fix turned out to be a
**head-count** effect rather than a target property: the same raw `auc` scores −0.069 at K=545 and
**+0.377** at K=8 (`notebooks/4_training.ipynb`). Removing the cause works at least as well as
compensating for it, without amplifying noise-dominated drugs.

It also carried a standing leak: `center` and `scale` were computed once over all 180 overlapping lines,
val and test included, and baked into the targets h5ad.

**The three arguments originally made for it**, recorded because two of them were sound and only the
second survives as a reason for anything:

1. **It equalizes the heads in the shared loss** — the real reason at the time, and correct as far as it
   went. Per-drug `auc` spread ranges from **0.034 to 0.302** across the 545 drugs (a 9× span, ~80× in
   squared error), so under a masked MSE on raw `auc` the wide-spread drugs supply nearly all the
   gradient to the shared trunk purely because of their *units*. What was missed: equal variance is the
   wrong equalizer, because it also promotes drugs whose spread is pure noise.
2. **It makes the metric readable** — with unit variance the per-drug-mean null scores MSE = 1.0 exactly,
   so any value below 1 is real per-drug signal. True, and the one genuine loss in retiring it: on raw
   `auc` every drug has its own null value. Mitigated by the fact that raw AUC is interpretable in
   viability units instead (MSE 0.0254 against a null of 0.030 ≈ 15 % of variance, RMSE ≈ 0.16).
3. **It removes the per-drug potency offset** — was flagged as the weakest of the three even when it was
   written, and it is simply vacuous: each drug is its own output row with its own bias term, so the head
   absorbs the offset either way.

**And the known problems, written up before it was retired** — these are what the retirement acted on:

- **The target discards `auc_mean`, so nothing downstream may select on it.** Subtracting the per-drug
  mean removes *potency* entirely, so any drug filter phrased as "kills ≥ N lines" selects on the one
  quantity the target throws away. This is what invalidated the learnability gate
  ([above](#the-learnability-gate-measured-potency-not-rankability)).
- **Dividing by `auc_std` amplifies noise for narrow-spread drugs**, as described above.
  `_zscore_per_drug`'s own docstring anticipated it — zero-spread drugs would "blow up" — and relied on
  `--min-cell-lines` plus the learnability filter to remove them, a safeguard the gate did not actually
  provide.
- **Candidate fixes that were never implemented**, cheapest first: floor the scale by dividing by
  `sqrt(auc_std² + σ_noise²)`; weight each drug by its reliable variance fraction
  `(auc_std² − σ_noise²)/auc_std²`; or select on `auc_std` directly and not patch the scale at all — the
  last of which is what the rebuild criterion now does.

The map was exactly invertible (`auc = auc_z * scale + center`, with `center`/`scale` saved in
`uns["ctrp_score_center"]` / `["ctrp_score_scale"]`), so any stored `auc_z` prediction can be returned to
native units. The `scale` vector is still used — it *is* the per-drug std of `auc`, which is the quantity
the rebuilt drug criterion selects on.

**What it affected.** Every result dated 13.07–27.07. Absolute MSEs on `auc_z` are not comparable to
those on any unstandardized target — a z-scored target has unit variance, so its per-drug-mean null sits
at ≈ 1.0 rather than ≈ 0.0097.

**Replaced by** raw `auc` winsorized at 1.1, with the weighting moved into the loss where it can be
estimated per fold on training lines only — which also closed the leak. No per-drug weight is applied:
the variance imbalance is a K=545 problem, and on a comparable panel the σ range is ~2.5× rather than
81×. Two mechanics are forced by a target near 0.7 instead of 0 — weight decay is applied to the weight
matrices only (`TrainConfig.no_decay_bias_and_norm`; until 12.08.2026 the narrower and non-standard
`exclude_output_from_decay`), and each head's bias is initialized at the fold's per-drug mean. Current
definition: [Step 03](03-model-and-training-design.md).

### Per-drug variance weighting — dissolved by a scope change, never needed

**Established** 27.07.2026 as plan items 1.0–1.2: estimate a pooled `σ_noise` from the replicated curve
fits, then weight each drug by `w_j = 1/σ_j²`, then by the reliable-variance fraction
`r_j = (σ_j² − σ_noise²)/σ_j²`.

**Overturned** the same day, by scope rather than by measurement. Per-drug variance weighting is a
**K=545** problem; on an 8-drug panel the variance ratio is 2.5× and the raw target works unweighted
(raw `auc` scores −0.069 at 545 heads and +0.377 at 8). `σ_noise` was therefore never estimated.

Kept here so the reasoning is not re-derived from scratch if the project returns to a large head count.
The estimate is feasible **pooled but not per drug**: `v20.data.curves_post_qc.txt` holds 387,130
(cell line × compound) fits of which only **7,708 (2.0 %)** are replicated, at most 3-fold. Per-curve
confidence intervals (`p1_conf_int_low/high`) are an alternative route.

*The process failure underneath this — continuing to work a problem the agreed scope had already
dissolved — is recorded under [Process failures](#process-failures).*

### "Neither representation ranks cell lines" — the K=545 null result

**Established** 27.06.2026 (`notebooks/4_training.ipynb` §B3): per-drug Spearman between predicted and
true response across held-out lines, over the 461 drugs with real per-line variance, came out at
**−0.02 (PCA) / −0.05 (scGPT)**, with only ~4 % of drugs above ρ = 0.3. Read at the time as the
project's central finding — that at this label resolution the task is barely learnable beyond the mean,
for either representation.

**Overturned** 13.07.2026, for two independent reasons:

1. **It is an average over 545 drugs, and the average destroys it.** On the drugs that carry real
   signal the same architecture reaches ρ ≈ 0.43–0.49 (`notebooks/archive/learnable_subset_training.ipynb`).
2. **The multi-task loss was unstandardized.** These runs used `mean_pv`, whose per-drug variance is
   wildly heterogeneous (spreads span 9×, ~80× in squared error), so a minority of wide-spread heads
   monopolized the shared trunk's gradient **by unit size, not by learnability**.
   `notebooks/archive/target_comparison.ipynb` reproduces the failure on demand: raw `auc` at K=545 gives −0.087
   (scGPT) / +0.016 (PCA), while a standardized target on the *same* drugs, model and split reaches
   +0.430 / +0.378.

Decomposed on the same 5 drugs throughout (`learnable_subset_training`), the move from the null to the working result
splits as: **target ≈ +0.29 (PCA) / +0.64 (scGPT)**, honest out-of-fold measurement ≈ +0.1 (27 → ~150
held-out lines), drug filtering ≈ +0.06. The target term dominates, and it is a genuine improvement in
the predictions rather than in the metric — on `mean_pv` the model's ranking of those drugs was
*negative*.

**What it affected.** The null was never clean evidence about scGPT vs PCA. Anything resting on it is
suspect rather than merely stale — above all the 8-run matrix conclusions (below).

**Replaced by** the per-drug results in [Step 05](05-multitask-results.md).

**It also took the standing net read down with it** — *"the ceiling is the label; the model learns the
per-drug mean, not cross-line sensitivity; PCA ≈ scGPT"*, held through June 2026. That reading is **true
on average and false on the drugs that matter**. The label ceiling is real (~150 independent cell lines;
ridge on line means ties the MLP), but "no gene representation can help" was never established.

### The 8-run matrix conclusions

**Established** 27.06.2026 at matched trunk and matched 512-d width, over
`{hvg5000, all_genes} × {X_pca, X_scGPT} × {single-paclitaxel, all-drugs K=545}`: scGPT overfits far
less (single-task gap 0.004 vs PCA 0.033) but PCA is competitive or better on raw accuracy
(heads-beating `hvg5000` 169 vs 147, ~~`all_genes` 138 vs 131~~, val MSEs within 0.0003) — so the
representation affects generalization, not predictive power.

**Overturned** 13.07.2026 in its *conclusions only*. Every all-drugs cell rests on the K=545 null above,
which the unstandardized loss substantially produced. The matrix **design** stands — the axes, the
matched `(128,64)` trunk, the 512-d width, the shared cell-line-grouped split — and so does the
single-task overfitting comparison, which involves no multi-task coupling.

**What it affected.** The all-drugs rows and the "PCA competitive on accuracy" reading.

**Replaced by** nothing yet: a re-run on a standardized or comparable target is required, and should be
expected to **change** the conclusions rather than refresh them ([TODO](../TODO.md)).

**Also struck 05.08.2026 — the `all_genes` half of the gene-set axis, on a second and independent
ground.** scGPT is embedded at `max_length=1200`, and at `all_genes` the cap binds in every cell, so it
received a random fraction of each cell's expressed genes while PCA received the whole gene set
(counts: [Step 02](02-preprocessing-and-embeddings.md#why-hvg-5000-is-the-default-03082026)).
The `all_genes` arms of the matrix therefore differ in gene set as well as in encoding and were never a
like-for-like comparison — this holds independently of the target defect above and would survive a
re-run on a standardized target. `hvg5000` is unaffected (1 cell of 53,513 above the cap), so the
matrix's **design** stands as stated, minus the claim that its second axis contrasted a filtered against
an unfiltered gene set *for scGPT*. Those embeddings were also generated unseeded. See
[Step 05](05-multitask-results.md#multi-task-masked-loss-over-all-545-ctrpv2-drugs-26052026).

### The ~50-d PCA baseline and the `(64,32)` PCA trunk

**Established** in the first matrix builds (through 14.06.2026): `add_pca` kept scanpy's default ~50
components, and PCA used a smaller `(64,32)` trunk than scGPT's `(128,64)`.

**Overturned** 14.06.2026 (trunk) and 27.06.2026 (width) — both settings handicapped PCA on axes that
have nothing to do with how genes are encoded, so the comparison was measuring capacity as well as
representation. The earlier scGPT leads on heads-beating (**135 vs 103**, **141 vs 80**) were a capacity
artifact and do not survive matching; with capacity matched, PCA is competitive or better on raw
accuracy. Matching PCA to 512-d *sharpened* the overfitting contrast, because the extra first-layer
capacity lets PCA fit the training set harder while val stays high.

**Replaced by** `add_pca.DEFAULT_N_COMPS = 512` and `DEFAULT_HIDDEN_DIMS = (128,64)` for both
representations; the full matrix was re-run at 512-d (run dirs `runs/20260627_1913xx_*`).

### The 13.07 five-drug numbers

**Established** 13.07.2026 on a 5-drug subset with the harsh gates (coverage ≥ 95 %, `auc_std` ≥ 0.15,
dynamic range ≥ 0.4, ≥ 20 lines killed *and* ≥ 20 spared).

**Overturned** 14.07.2026: the filter was corrected to `learnability = min(#killed, #spared)` with
coverage ≥ 90 % and widened to the top **10**, and every experiment was repeated. The conclusions held;
the margins shrank.

| | 13.07 (5 drugs) | 14.07 (10 drugs) | source |
|---|---|---|---|
| K=10 out-of-fold ρ (PCA / scGPT) | 0.42 / 0.49 | **0.360 / 0.396** | `outputs/legacy/target/target_comparison.csv` |
| K=545 ρ (PCA / scGPT) | 0.378 / 0.430 | **0.316 / 0.328** | `outputs/legacy/target/target_comparison.csv` |
| Ridge line-level (PCA / scGPT) | 0.428 / 0.428 | **0.343 / 0.320** | `outputs/legacy/ablations/ablation_capacity.csv` |
| MLP `(128,64)` (PCA / scGPT) | 0.428 / 0.487 | **0.356 / 0.402** | `outputs/legacy/ablations/ablation_capacity.csv` |
| scGPT − PCA gap | +0.075 | **+0.036 (K=10) / +0.012 (K=545)** | derived |

The consequential change is the last row: **scGPT's edge over PCA shrank to roughly seed-noise size on
the wider set**, so it must not be headlined. Both panels are in any case superseded by the gate defect
above.

**Replaced by** the 14.07 column, itself superseded by the panel void.

**And neither column was ever a generalization number** — a limitation stated when they were produced,
not discovered later. The learnability gates were computed on all 180 lines, val and test included, so
the selection saw held-out labels. Both are **best-case diagnostics**: they answer "does any cross-line
signal exist here?" (yes) and not "how well does this generalize?". Making them reportable needs the
selection run **inside each CV fold on training lines only**, which remains blocking for any headline
number ([TODO](../TODO.md)).

### The first DrEval benchmark — a val-split leak

**Established** 14.07.2026 (`notebooks/analysis/evaluation/dreval_benchmark.ipynb`) with `drevalpy` 1.5.1 under their LCO
protocol, on the 5-drug best-case subset.

**Overturned** the same day: `run_oncomlp` was passing the **test fold** as the validation loader, so
early stopping and best-epoch selection ran on the test fold. Optimistically biased the `OncoMLP` rows
only — the DrEval baselines were never affected. Fixed in **`ee07b00`** by using DrEval's real
`sp['validation']` split (`split_validation=True`).

**Replaced by** the leak-free re-run on the 10-drug set, mean over 5 folds
(`outputs/dreval/dreval_lco_results.csv`): `OncoMLP` scGPT normalized ρ **0.357** / R² **0.114**, PCA
0.340 / 0.086, their `SingleDrugRF` (scgpt) 0.339, `NaiveMeanEffects` 0.000. Still above the naive bar
that half the published field fails, but only faintly above PCA and their RF, and below the paper's best
LCO model (~19 % normalized R²). Use these.

### The same val-split leak, in the code that produced everything else

**Established** as fixed 14.07.2026 — for the DrEval benchmark only. The entry above ends "Fixed in
`ee07b00`", and it was, *in that notebook*.

**Overturned 12.08.2026 (review item 7).** The identical defect sat untouched in
`scripts/training/cv.py::oof_predictions`, which produces every out-of-fold prediction this project
has, and in `train_multitask.cv_evaluate`, which produces every CV metric. `train_model` restores the
lowest-validation-MSE checkpoint; both passed it the held-out fold and then scored that same fold. So
the reported figure is a minimum over up to 25 epochs on the data being reported. Optimistically
biased, in both arms and in the ridge control's absence — the ridge baseline has no early stopping, so
it alone was never selected this way, which flatters the MLP against it.

**Not uniform across the arms, which is why it affects the comparison and not only the level.** Selected
epochs in the last committed fold log (`outputs/legacy/panel_void_8drug/panel_training_folds.csv`, void 8-drug panel —
cited for the mechanism, not as a result): PCA `[1,1,3,1,1]`, scGPT `[10,11,2,21,4]`. The PCA arm's
scored checkpoint was chosen from among near-tied epoch-1 states on the lines it was scored against;
scGPT's came from a real trajectory. The same fact was already on record as the cause of the PCA arm's
`mps` non-reproducibility, without the selection-bias reading being drawn from it.

**Replaced by** `scripts/training/cv.py::inner_holdout` (12.08.2026): 15 % of each fold's training lines,
grouped, become the early-stopping set and the fold is scored clean. Design, the arbitrariness of 0.15
and the separate inner seed: [Step 03](03-model-and-training-design.md#the-early-stopping-set-is-nested-inside-the-training-lines-12082026).
**No number has been re-run under it** — every result on record still carries the bias, and they all move
down at the sweep, from the removed selection and from ~104 training lines per fold where there were 122.

### The evidence that closed model-side tuning

**Established** 13.07.2026 by `notebooks/archive/ablations_and_rescue.ipynb`: four model-side knobs are
flat, `RidgeCV` on 150 line-mean embeddings ties the PCA MLP, and scGPT-plus-a-hidden-layer is the only
configuration that clears ridge (+0.06). It became the project's standing instruction — *don't spend more
time on architecture or hyperparameters, the remaining levers are label-side* — and the reason MIL is
next rather than a bigger model.

**Overturned 12.08.2026 (review item 8).** The notebook's `oof()` calls
`train_model(model, tl, DataLoader(va_ds, ...))`, so **the scored fold is the checkpoint-selection set** —
the same leak recorded above, in the notebook that produced the tuning tables. Three further defects
compound it, each independently sufficient: the target is `auc_z`, [retired](#auc_z-as-the-training-target)
27.07.2026; the five drugs come from the
[learnability gate](#the-learnability-gate-measured-potency-not-rankability), discredited 12.08.2026; and
the notebook is archived under ground 2 — the `auc`/`auc_z` targets file no longer exists, so it cannot
be re-run to check any of this.

**The bias is one-sided, which is what makes it more than bookkeeping.** `RidgeCV` selects its penalty by
generalized cross-validation inside the training fold and has no early stopping, so the ridge column was
never flattered while every MLP column was — and on the PCA arm the two are *equal* (0.428 / 0.428). The
honest ordering could therefore be ridge **above** the PCA MLP, which is a different claim from the one
on record. The +0.06 that carries the project's strongest Q1 statement is likewise a margin between a
flattered number and a clean one.

**Not overturned by this:** the *rescue* result in the same notebook (removing regularization recovers
~70 % of the K=545 collapse) is a within-notebook comparison in which every arm carries the same leak, so
the ranking it establishes survives even though its levels do not. And the general argument — ~150
independent labels cannot support architecture search — is an argument about the data, not a measurement,
and stands on its own.

**Replaced by** nothing yet. What gets re-derived and what does not is [TODO](../TODO.md) item 8C: the
minimal re-run at R4 is trunk `(128,64)` vs a bare linear head, both representations, against ridge on the
same folds and the rebuilt panel — the single comparison the load-bearing claim rests on, not the four-knob
sweep. Until it exists, "model-side tuning is closed" is an instruction without live evidence.

### The two uncentred-target mechanics ran in one training path of three

**Established** 27.07.2026: on a target centred near 0.9 the head biases must start at the per-drug means
and must not be decayed toward 0, both recorded in
[Step 03](03-model-and-training-design.md#the-uncentred-target-is-handled-the-same-way-in-every-training-path)
as mechanics rather than as choices.

**Overturned 12.08.2026 (review item 8) — they were only ever applied in `cv.oof_predictions`.**
`train_multitask.cv_evaluate` and `train_multitask.train_rep`, which together produce the entire 8-run
matrix, initialized no head bias, and `4_training.ipynb` (§B) never set `exclude_output_from_decay`, whose
default was `False`. So the matrix path trained heads that start near 0 and are pulled back toward 0,
while the panel path started them at the drug mean and never decayed them. Two arms of the same project
under different rules, with the docs describing only one of them. `dreval_benchmark.ipynb`, which builds
`OncoMLP` by hand, also has neither — left alone because it is broken twice over and belongs to item 11.

**The flag itself was also wrong.** `exclude_output_from_decay` exempted the *whole output `Linear`,
weight matrix included*, and still decayed the LayerNorm parameters and hidden-layer biases — while
Step 03 described the standard convention ("biases and LayerNorm parameters therefore go in a
`weight_decay=0` group"), which is a different grouping. Neither the code nor the doc matched the other.

**Replaced by** `TrainConfig.no_decay_bias_and_norm` (default on), implemented in
`training_utils._decay_param_groups` — every weight matrix decayed including the output's, every bias and
normalization parameter exempt, which is HuggingFace `transformers`' `Trainer.create_optimizer` grouping
from the BERT reference implementation — and by `init_head_bias=True` on all three paths, via
`OncoMLP.init_head_bias_` over `cv.per_drug_line_mean`. **Takes effect at R4.** The size of either
handicap has never been measured; it is not measured now either, and the change is made on the convention
rather than on evidence.

### The Steps 04–05 numbers as a comparable baseline

**Established** 08.05–27.06.2026 on `mean_pv`, the only target until 13.07.2026.

**Overturned** 13.07.2026 as a basis for comparison — not as results. Absolute MSEs are not comparable
across target scales, and the two scores rank cell lines differently: globally they correlate at
ρ ≈ 0.97, but that is inflated by between-drug potency spread. *Within* a drug, across lines — the only
variation the model must predict — the median Spearman between them is **0.72** (min 0.42). The change
is not cosmetic.

**And unreproducible since 11.08.2026.** `mean_pv` was removed with its reader code when the target
moved to DrEval's reprocessed CTRPv2, so `--score mean_pv` now raises. These numbers cannot be
regenerated by any code in the tree; they stand only as a record of what was measured.

**The Step 04 numbers, moved here 11.08.2026** when that page was reduced to the design that survives
the target change. All are `mean_pv` val MSE on `split_paclitaxel`:

*Split sizes* — 170 cell lines with a paclitaxel label → 119 train / 25 val / 26 test; cells 31,824 /
5,035 / 7,508, 9,146 unassigned, out of 53,513 total with 44,367 labelled. These move with the target
too: the overlap becomes 181 lines and the labelled-cell count changes with it.

| Stage | PCA | scGPT |
|---|---|---|
| Random split (the deliberate mistake) | train 0.0022 / **val 0.0011** — the leak's tell | train 0.0132 / val 0.0137 |
| Cell-line-grouped, unregularized (256-dim MLP) | val 0.0390 | val 0.0437 |
| Grouped + aggressive regularization (hidden 64, dropout 0.5, wd 1e-3) | train 0.0082, val ~0.0380, gap **0.029** | train 0.0260, val ~0.0371, gap **0.013** |
| No-HVG regularized (08.05) | ~0.0375 (ep 10) | ~0.0371 (ep 10) |
| HVG-5000, old model (BatchNorm/ReLU) | 0.0362 (ep 5) | 0.0354 (ep 50) |
| **HVG-5000, upgraded model** | **0.0351 (ep 8)** | **0.0336 (ep 14)** |

Both arms early-stopped under the upgraded setup — PCA at epoch 18 (best 8), scGPT at 24 (best 14).
The `runs/` directories are gitignored and no longer present, so these best-epoch values are all that
survives of those curves.

**What survives the target change is the *shape*, not the numbers:** PCA's val MSE collapsing once the
split became cell-line-grouped is a structural consequence of a label that is constant within a line,
not a property of `mean_pv`. That argument stays in [Step 04](04-single-task-results.md).

---

## Retracted claims

These were written into the docs as findings and later withdrawn. The distinction from *superseded* is
that nothing replaced them — they were simply not true.

### The curve fit preserves signal the dose-average destroys

**Claimed** 13.07.2026, in [Step 03](03-model-and-training-design.md) and
[Step 05](05-multitask-results.md), when the target moved from `mean_pv` to the curve-fit AUC: that the
dose-averaged viability was destroying signal the sigmoid fit preserves.

**Falsified the same day** (`notebooks/archive/target_comparison.ipynb`, re-run with all three targets, 95 %
bootstrap CIs, per-drug dots and a 3-seed check). Trained head-to-head, `mean_pv` and raw `auc` are
statistically identical *everywhere* — K=5: 0.450 vs 0.439 (PCA), 0.481 vs 0.482 (scGPT); K=545: +0.027 /
−0.070 vs +0.016 / −0.087. CIs fully overlap.

**The curve fit buys no measurable accuracy.** Keep it for principled reasons — it is a post-QC sigmoid
fit, it ships per-parameter confidence intervals, and it is the metric family GDSC2 reports, which
cross-database work needs — but never claim it improves prediction. The entire effect attributed to "a
better label" came from **per-drug standardization of the shared loss**, i.e. from fixing the loss, not
the label.

### The learnability filter was validated against the ρ the model achieves

**Claimed** as a rank correlation of **+0.357** between the filter's learnability score and achieved
per-drug ρ, plus "ρ > 0 on 76 % of drugs, median 0.12" across all 545.

**Retracted** 14.07.2026. Both figures came from `scratchpad/learnability_validity.py`, whose output
`learnability_vs_achieved.csv` was never committed and is gone. They are **not reproducible from any
committed artifact** and are therefore dropped, not merely re-labelled.

What survives: the filter selects drugs purely from label statistics (coverage, spread, killed/spared
counts), and it is a legitimate **diagnostic device** — not a measured ranking of learnability.

> ⚠️ **This retraction is not yet fully applied.** `notebooks/README.md` still states that the score
> "correlates only +0.36 with the ρ the model actually reaches" — the retracted figure, live in a tracked
> file. It needs removing.

### `ml210` was rejected on coverage

**Claimed** in an earlier note. **Wrong** — corrected 14.07.2026 against the committed filter output
`notebooks/outputs/legacy/learnability/ctrp_drug_learnability_auc.csv`, where `ml210` has coverage **0.944**,
clears the 0.90 gate, and is `passes_gate=True, selected=True`. It is one of the 10 selected drugs.

### The panel was chosen blind to our labels

**Claimed** briefly 25.07.2026, when the panel moved from our filter to published determinants.

**Corrected the same day.** The candidate list was ranked by `min(kill, spare)` — computed on our own
`auc` values over all 180 lines, val and test included — *before* the literature criterion was applied.
Compounds with published determinants but little spread in our data fell out: `sirolimus`, `neratinib`,
`clofarabine`, `cytarabine hydrochloride`, `gdc-0941`. Correct wording: **literature-anchored,
spread-verified**.

What that correction fixed was the *wording*. The deeper consequence — that the pool itself had been
pre-filtered by a discredited criterion — surfaced only on 27–28.07 and voided the panel
([above](#the-8-drug-literature-panel-and-every-number-computed-on-it)).

### The prediction shrinkage is a defect to fix with lighter regularization

**Claimed** as an open TODO item: predictions are shrunk toward each drug's mean, so loosen the
regularization.

**Withdrawn** 13.07.2026 (`notebooks/archive/ablations_and_rescue.ipynb` §1). `pred_std ≈ ρ × true_std` (scGPT: 0.47
against ρ = 0.48) is exactly what an MSE-optimal predictor **must** do — the conditional mean shrinks
toward the prior in proportion to how little signal exists. It is correct calibration, not timidity, and
loosening dropout *raises* MSE. To report in AUC units, divide by ρ; Spearman is unchanged.

### The scGPT binning invariance was "verified on 200 cells"

**Claimed**, undated, in Step 02's account of what scGPT's embedding path does to its input: binning
"is a rank transform and invariant to any monotone per-cell rescaling — *verified on 200 cells, where
CPM, raw counts, `normalize_total` and `log1p` give identical bins*."

**Retracted 10.08.2026 (review item 4).** No notebook, script or `outputs/` artifact ever backed it;
a search of `notebooks/`, `scripts/` and `docs/` finds nothing. And it could not have been produced in
this project: SCP542 distributes only `CPM_data.txt`, so **the raw counts it names do not exist here**
and there was never anything to compare CPM against.

**The wording was wrong in a second way.** `_digitize` (`scgpt/preprocess.py:239`) resolves values
landing on a repeated bin edge with `np.random.rand`, so binning is **stochastic** wherever a cell has
tied expression values. Without a fixed RNG state the same input does not bin identically to *itself*,
let alone across four input formats. (`gen_embeds.py:243-250` seeds it, which is why the embeddings
are reproducible at all.)

**The conclusion survives; only its status changed.** The invariance follows from the source — bin
edges are quantiles of the cell's own non-zero values, so any strictly monotone transform moves values
and edges together — and Step 02 now states it as
[an argument rather than a measurement](02-preprocessing-and-embeddings.md#what-scgpt-is-fed-and-why-its-scale-does-not-matter).
Nothing downstream depended on the stronger claim, and the input scGPT is actually fed is unchanged.

### "scGPT clears the ridge control" was a first

**Claimed** 27.07.2026 in the first draft of the Step-1 write-up.

**Corrected the same day: it is a replication, not a first.** scGPT MLP over its ridge is **+0.077** on
the 8-drug panel against **+0.082** on the 14.07 10-drug panel (0.402 vs 0.320,
`outputs/legacy/ablations/ablation_capacity.csv`). The replication is the stronger claim of the two, since the
drug identities on the later panel were named by citation rather than by our own labels — though the pool
they were drawn from was not, which is a separate defect.

### PCA prefers the full transcriptome

**Claimed** as a hunch from the earlier matrix runs — and its mirror, that HVG filtering specifically
helps scGPT.

**Not reproduced** 28.06.2026 (then `notebooks/4_training.ipynb` §B4; since 03.08.2026
`notebooks/analysis/qc/verify_variants.ipynb` §9), the gene-set sweep at 1k/2k/3k/5k plus
`all_genes` under identical 5-fold CV: both representations are **flat across the whole axis** (PCA
~203–216 heads beating baseline, scGPT ~184–193), val MSE constant at 0.0105–0.0107. PCA's `all_genes`
value of 204 sits mid-band, *below* `hvg3000`'s 216. There is no sweet spot and no all-genes advantage
for either representation.

> ⚠️ **Narrowed 05.08.2026 — the scGPT half of the last sentence.** The claim being refuted is about
> **PCA**, which read all 22,722 genes, so the refutation itself is untouched. But scGPT never received
> the full transcriptome: at `max_length=1200` the cap binds in every cell of `all_genes`, so it got only
> a random fraction of each cell's expressed genes
> (counts: [Step 02](02-preprocessing-and-embeddings.md#why-hvg-5000-is-the-default-03082026)).
> "No all-genes advantage **for either representation**" therefore holds as written for PCA only. For
> scGPT the flat result supports the narrower statement that roughly twice as many randomly drawn genes
> do no better than the dispersion-selected ones — see
> [Step 05](05-multitask-results.md#gene-set-sweep--heads-beating-vs-gene-count-incl-all_genes-28062026).
> The mirror claim that HVG filtering *helps* scGPT stays refuted, since the four HVG points are all
> below the cap and genuinely comparable.

> ⚠️ **03.08.2026 — the supporting numbers are superseded, the refutation is not.** Those figures came
> from the retired `mean_pv` target; §9 now runs on `auc` and does not read the `mean_pv` cache, so the
> sweep has no live numbers until re-run ([Step 05](05-multitask-results.md#gene-set-sweep--heads-beating-vs-gene-count-incl-all_genes-28062026)).
> The claim stays refuted — nothing has been produced that revives it — but it must not be re-asserted
> from these numbers without the re-run.

---

## Refuted hypotheses

Tested deliberately, and the answer was no. Each one rules something out, which is why they are kept.

### Inverse-density loss weighting improves ranking

**Hypothesis.** Response values are unevenly distributed within each drug, so up-weighting rare values —
the regression analogue of class weighting — should reduce the shrinkage and improve ranking. Implemented
in `scripts/training/density_weighting.py`, fitted per fold on training lines only via
`scripts/training/cv.py`.

**Refuted** 27.07.2026 (`notebooks/4_training.ipynb`): **−0.006 (PCA) / −0.008 (scGPT)** mean
Spearman. Per drug a wash — `selumetinib` +0.09 / +0.06, `tanespimycin` −0.06 / −0.06. The
pre-registered expectation ("MSE worse, Spearman better") also failed; both stayed flat, which is what a
null intervention looks like.

**The mechanism did fire, and that is what makes it a null rather than a bug.** Predicted spread rose
from 0.062 to 0.082 (PCA) and 0.080 (scGPT) — precisely the reduced shrinkage the method is designed to
produce. A weighting that never reached the loss (wrong sign, weights dropped by the mask, an indexing
slip) would leave predictions numerically identical. They are not: the model demonstrably hedges less.
*Did not work* and *was broken* are different claims, and only the first is supported.

**Why it was a null is also understood.** Coherent with
`notebooks/archive/panel_distributions.ipynb`: once the artifacts above `auc` 1.1 are winsorized away, every
drug has |skew| ≤ 0.47 (`topotecan` +2.42 → +0.18), so there was almost no imbalance left for an
imbalance correction to act on.

**What the null buys — this is the reason to keep it.** It removes a standing candidate explanation for
the shrinkage. Predictions span 0.08 against a true spread of 0.171, and one hypothesis was that the
objective is dominated by the crowded middle of each drug's response range. Tested and rejected:
pointing the loss at the sparse extremes does not close the gap. What remains is too little signal
across ~150 independent cell lines — a **label-side** problem, which is the direct argument for MIL and
for more cell lines, and against further objective engineering.

**Decision: do not carry the weighting into Step 2.**

> ⚠️ **Do not read the sign of these deltas.** The PCA unweighted arm is not bit-reproducible on `mps`:
> four identical runs gave 0.313 / 0.315 / 0.317 / 0.320, while every other arm reproduced exactly. The
> cause is that PCA peaks at epoch 1 (best epoch per fold `[1,1,3,1,1]` vs scGPT `[10,11,2,21,4]`), so
> its checkpoint is chosen among near-tied states. The weighting deltas lie inside that band.

### The cell-line effect is largely proliferation

**Hypothesis** (proposed 27.07.2026). A 72 h viability assay confounds drug effect with division rate —
the motivation for GR metrics (Hafner et al., *Nat Methods* 2016; in `references.bib`). Proliferation is
the most visible axis in scRNA-seq, and cell cycle is the top recurrent program in Kinker's analysis of
*this* dataset. If the cell-line effect were largely proliferation, it would explain both why
ridge-on-line-means ties the MLP and why removing the cell-line effect costs ~20 % of the signal.

**Test** (`notebooks/analysis/evaluation/diagnostics.ipynb`). The dataset ships the authors' own published per-cell scores
(`G1/S_score`, `G2/M_score`, plus all 12 heterogeneity programs), so no scoring choice of ours entered.
Averaged per line and correlated against the cell-line effect — the mean over drugs of the per-drug
z-scored AUC, which is DrEval's definition. Outputs:
`outputs/diagnostics/line_effect_vs_programs.csv`, `line_effect_vs_proliferation.png`.

**Refuted.** ρ = **−0.050**, p = 0.50, n = 180, r² = 0.003. The best of the 12 programs is ProtDegra at
−0.113 (p = 0.13, not significant). At n = 180 the 95 % interval is roughly ±0.145, so anything above
|ρ| ≈ 0.2 is excluded — proliferation explains **at most ~4 %** of the effect.

**Both controls pass**, so there was signal to find: the cell-line effect is real and large (std 0.398
z-units, range −2.10 to +1.41), and the cycle scores are not line-centred (between-line std 0.233 vs
within-line 0.593).

**Caveat.** These are transcriptional snapshots of cycling fraction, not measured doubling times. The
definitive test would use DepMap doubling times.

**What the refutation buys.** The cell-line effect is large and explained by **none** of the published
heterogeneity programs of this dataset. That strengthens the case for removing it via a DrEval-aligned
residual target: we would not be discarding a biological program we want to keep.

### The model is over-regularized, or too small

**Hypothesis** (13.07.2026). Predictions are shrunk and the model plateaus, so the network is either
over-regularized or short of capacity.

**Refuted** (`notebooks/archive/ablations_and_rescue.ipynb`), four knobs on the corrected setting, out-of-fold per-drug
Spearman:

| knob | range tested | PCA | scGPT |
|---|---|---|---|
| regularization | none → heavy | 0.42–0.44 | 0.44–**0.49** |
| capacity | 74,629 → 2,565 params | 0.41–0.43 | 0.44–**0.49** |
| batch size | 32 / 128 / 512 | 0.43–0.44 | 0.46–**0.49** |
| sample reweighting | line-balanced, focus-extremes | 0.41–0.43 | 0.48–0.49 |

Every axis flat, defaults at or within noise of the best. With regularization **off**, PCA drives *train*
MSE to ≈ 0.01 — near-perfect memorization of the training lines — and still reaches only 0.42
out-of-fold. A model that can overfit that hard is not being suppressed by its regularizer: it is out of
**signal**, not capacity. Heavy regularization is the only setting that hurts, and it does so via
over-shrinkage (`pred_std` 0.33), not lost ranking.

> **Scope correction (14.07.2026), and it matters.** Those four ablations ran on the **corrected**
> (K=5) setup. They show the knobs do not *improve* the corrected model — they do **not** show the knobs
> could not have *fixed* the K=545 collapse. Tested separately on the broken setting (K=545, raw `auc`,
> scGPT, ρ = −0.063; `outputs/legacy/ablations/rescue_k545.csv`):
>
> | intervention | ρ |
> |---|---|
> | heavy regularization | −0.091 |
> | line-balanced sample reweighting | −0.078 |
> | smaller model (74,629 → 16,645 params) | −0.053 |
> | batch size 32 | +0.027 |
> | **no regularization** | **+0.234** |
> | **per-drug (task) reweighting** | **+0.433** |
>
> **Removing regularization recovers ~70 % of the collapse** — mechanistically consistent with a
> **capacity competition between heads**: at dropout 0.5 the trunk cannot serve both the loud noisy
> drugs and the learnable ones, so the loud ones win. But it is a **symptom fix, not a cause fix**, and
> the interaction proves it: on the corrected setting the *same* regularization is **optimal** (0.488 vs
> 0.456 with none). The model was never over-regularized in absolute terms — it was over-regularized
> *relative to a mis-weighted loss*.

**Decision: stop tuning the model**, on the corrected loss. At ~150 independent labels, architecture
search cannot buy signal.

### Line-balanced reweighting will help

**Hypothesis** (13.07.2026), and a principled one: the entry-pooled loss lets a 500-cell line pull 10×
harder than a 50-cell line, although both carry exactly **one** independent label.

**Refuted** — it changes nothing: scGPT 0.485 against 0.488. In hindsight this is forced, because the
ridge-on-line-means control **is** the fully line-balanced limit, and it ties the PCA MLP.

The underlying imbalance is nonetheless real and is documented as a live defect in
[Step 03](03-model-and-training-design.md) — the per-line loss share spans a factor of 82. What is
refuted is that *reweighting the existing per-cell objective* fixes it. MIL removes it structurally,
because one bag is one line is one example.

---

## Dead ends

Directions that were investigated and abandoned. Kept so they are not re-opened without a reason.

### Migrating the cell-line join to persistent identifiers

**Planned** since 28.07.2026 as review item 2B: replace the `ccl_name_norm` join key (lowercased,
hyphens stripped) with DepMap `ACH-` or Cellosaurus `CVCL_` identifiers, and gene symbols with Ensembl
IDs. The item asserted that the identifiers "exist in the catalog but only in the exploratory
`drug_catalog`, never in the production path", making this a refactor. It also expected the migration
to expose invisible join errors.

**Dropped 10.08.2026 (Selin).** Both halves of the premise failed the audit:

- **The identifiers do not exist.** `v20.meta.per_cell_line.txt` carries `master_ccl_id, ccl_name,
  ccl_availability, ccle_primary_site, ccle_primary_hist, ccle_hist_subtype_1` — no `ACH-`, no
  `CVCL_`. SCP542's `Metadata.txt` carries none either. The only `ACH-` strings in the repository are
  PRISM's own row labels (`ACH-000361::P108::PR500A::REP1M`) appearing as transient cell output in
  `notebooks/analysis/harmonization/drug_catalog.ipynb`; nothing is written to a catalog file. So the
  migration was never a refactor — it required adopting an external resource (DepMap `sample_info.csv`
  or a Cellosaurus release) as a new dataset, with its own version, provenance and licence.
- **The invisible join errors were found without it.** The audit checked the join directly and it
  holds: no normalization collisions on either side, and 189 of 190 name matches independently agree on
  CCLE primary site. The one real defect it did expose, `H292`, was identified and fixed by name — with
  Cellosaurus `CVCL_0455` used as *evidence for one alias*, which is a lookup, not a dependency
  ([Step 01](01-datasets-and-harmonization.md#the-join-dropped-a-screened-cell-line-h292-10082026)).

**What this does not close.** Compounds are a different case: CTRPv2 already ships `broad_cpd_id` (BRD)
in `v20.meta.per_compound.txt` and the production path ignores it, joining on `cpd_name`. That is
harmless today — 545 compound names, zero collisions, and drugs never join across sources — but it is
exactly what the cross-database step needs
([Step 06](06-planned-work.md#a-cross-database-integration)). If PRISM or GDSC are ever wired into
training, the identifier question returns there, with BRD already in hand and only the cell-line side
needing an external mapping.

### Kinker's two named associations do not transfer to this task

**Tried** 25.07.2026, as the most natural source of drug candidates: the SCP542 source paper
([Kinker et al., *Nat Genet* 2020](https://www.nature.com/articles/s41588-020-00726-6)) names two
drug-response associations, so those compounds should be the best-supported panel candidates.

**Abandoned.** The γ-secretase / NOTCH inhibitors (`mk-0752`, `semagacestat`, `l-685458`) and the MDM2
inhibitors (`nutlin-3`, `hli 373`, `sj-172550`, `serdemetan`) have excellent coverage (~0.96–0.98) but
almost no lines below the potency threshold the gate was applying at the time.

**The reason the associations don't transfer is the important part, and it is not about potency.** Kinker
associated the **variability of a heterogeneity program within a cell line** with response. We regress
the **mean AUC level across cell lines**. Those are different quantities, and evidence for one is not
evidence for the other.

> ⚠️ **Read this entry with the gate retraction in mind.** The original note phrased the finding as
> "these compounds kill 0–1 of ~175 lines", which is the discredited criterion
> ([above](#the-learnability-gate-measured-potency-not-rankability)) — `nutlin-3` has spread 0.147,
> comparable to `dasatinib`. So the dead end is specifically that **Kinker's stated associations** do not
> carry over to a cross-line mean-AUC regression. It is **not** that these compounds are unusable, and
> `nutlin-3` in particular is expected to re-enter on the rebuilt pool.

### `kx2-391` carries drug-specific signal

**Believed** 13.07.2026: `kx2-391` was one of the five learnable drugs, and the drug where scGPT most
clearly beat PCA (0.28 vs 0.11) — read at the time as the clearest evidence for the representation.

**Abandoned** 14.07.2026 (`outputs/dreval/dreval_normalized.csv`). Once the cell-line effect is removed,
its correlation collapses from raw **0.283** to **0.006**, against a naive baseline of 0.584 — the entire
apparent signal *was* "this cell line is fragile", with no drug-specific biology. Exactly the artifact
class the DrEval paper describes, found in our own results.

It was excluded from the literature panel on this basis. For contrast, nothing on the 8-drug panel
behaves this way: the weakest are `topotecan` (0.296 → 0.192) and `vincristine` (0.407 → 0.280), while
`dasatinib` (0.546 → 0.558) and `afatinib` (0.325 → 0.368) *gain* — their signal is orthogonal to
fragility.

### scDrugAtlas and ClinTox as data sources

**Tried** 26.03–31.03.2026, while choosing the datasets.

- **scDrugAtlas** (<http://drug.hliulab.tech/scDrugAtlas/>) — abandoned as unusable without provenance.
  Prof. Liu was contacted about the original cell-line IDs in the consolidated downloads, the source
  publication behind each dataset, how the datasets were integrated, whether documentation or a data
  dictionary exists for the consolidated h5ad files, and whether IC50/cell-line annotations exist for
  matching to GDSC. No resolution. Concrete blockers: cell lines cannot be identified reliably, the
  consolidated file format and variable definitions are unclear, and in `breast_cancer_palbociclib` the
  response appears to be encoded as binary. The advisor additionally flagged that the atlas is
  Harmony-processed and warned against direct cross-dataset merging.
- **ClinTox** (<https://tdcommons.ai/single_pred_tasks/tox/#clintox>) — binary toxicity prediction from
  SMILES. Rejected because the labels are binary only, while the project needs a continuous response
  target ([Step 03](03-model-and-training-design.md)).

**Also unresolved:** the GDSC documentation link (`depmap.sanger.ac.uk/documentation/gdsc/`) was dead,
and the DepMap/GDSC team was asked how the "GDSC2 IC50 Data" functional dataset was processed. No
response. GDSC remains downloaded but unused, and is not a modelling priority.

### Considered and never pursued

- **A purpose-built scRNA-seq transformer, or a VAE**, trained from scratch as the representation instead
  of a pretrained model (raised 06.04.2026). Dropped in favour of using **scGPT** as a frozen prior —
  cheaper, and it makes the comparison against PCA a question about representations rather than about
  training budgets.
- **DeepInsight / scDeepInsight** — encoding expression as images for a CNN (raised 06.04.2026). Never
  tried; no result either way.

Both are recorded so they are not re-proposed as though untried. The other 06.04 ideas — bulk or
pseudo-bulk pretraining, fine-tuning on specific cancer types, fine-tuning on clinical data — did become
plan items and live in [Step 06](06-planned-work.md) and [TODO](../TODO.md).

### Retired code paths

Kept because the docs' history refers to them, and because it explains why older notes cite files that
no longer exist at those paths.

| Retired | When | Superseded by |
|---|---|---|
| `train_baseline.py`, `train_scGPT.py` | 26.05.2026 | `train_multitask.py --drugs paclitaxel --use-rep X_pca\|X_scGPT` — K=1 reduces exactly to plain MSE |
| `notebooks/hvg_vs_all_genes_umap.ipynb`, `notebooks/scgpt_umap.ipynb` | 27.06.2026 | consolidated into `notebooks/analysis/qc/verify_variants.ipynb` |
| `notebooks/10_ablations.ipynb` | — | consolidated into `notebooks/archive/ablations_and_rescue.ipynb` |
| `notebooks/01_scDAExploration.ipynb` | 30.07.2026 | renamed to `notebooks/archive/scdrugatlas_exploration.ipynb`. It explores **scDrugAtlas**, not SCP542 — the index's notebook table mislabelled it for months. Archived because the data source itself was [rejected](#scdrugatlas-and-clintox-as-data-sources), not because the notebook is wrong |
| `notebooks/03_analysis.ipynb` | 30.07.2026 | **un-archived** and renamed to `notebooks/archive/ctrp_prism_repurposing.ipynb`. Read-only CTRP→PRISM and clinical-phase mapping; writes nothing, but it is the only notebook that loads `GDSC2_fitted_dose_response_27Oct23.xlsx`, which the "externalize the spread requirement" item needs |
| `notebooks/result_evaluation/target_comparison.ipynb` | 11.08.2026 | **archived to `notebooks/archive/`.** It compared `mean_pv` vs `auc` vs `auc_z`; all three were removed with their reader code on 11.08.2026, so `layout.CTRP_SCORES` now rejects every target it asks for and the notebook cannot run. Its conclusion — retire `auc_z` — is [above](#auc_z-as-the-training-target) |
| `notebooks/data_and_harmonization/replicate_variation.ipynb` | 11.08.2026 | **archived to `notebooks/archive/`.** It measured how far apart two screenings of the same (cell line, drug) fall, and it existed because `_build_drug_table` **averaged** them. It no longer does: DrEval fit one curve across replicates, and `_deduplicate_measurements` raises rather than reconciling. Its measurements are kept, reframed as the evidence for the target switch, in [Step 01](01-datasets-and-harmonization.md#genuine-repeats-are-averaged-and-they-disagree-more-than-the-targets-own-spread-10082026) |
| `notebooks/result_evaluation/ablations_and_rescue.ipynb` | 11.08.2026 | **archived to `notebooks/archive/`.** Its argument — why the 545-head model failed and what rescued it — is a head-to-head between `auc` and `auc_z` across five code cells, and both targets were removed the same day, so it cannot run. Re-wiring it to `auc_cc` / `ln_ic50_cc` would not preserve the argument but replace it with a different question. The σ²-weighting finding it demonstrated is [above](#auc_z-as-the-training-target) |
| `ctrp_to_h5ad._load_score_values`, `._load_ctrp_long`, `._zscore_per_drug` | 11.08.2026 | deleted, not archived — the readers for CTRPv2's own `v20.*` tables and the `auc_z` standardisation. Superseded by `_load_drevalpy_long` + `_deduplicate_measurements` ([Step 01](01-datasets-and-harmonization.md#the-target-moved-to-drevals-reprocessed-ctrpv2-11082026)). In git history; the defect that ended them is [above](#the-auc-target-was-divided-by-the-wrong-quantity) |
| `scratchpad/learnability_validity.py` | never committed | gone; its figures are [retracted](#the-learnability-filter-was-validated-against-the-ρ-the-model-achieves) |

---

## Process failures

Not results — the working-method problems behind several entries above. Kept because each cost real time
and at least three are structurally likely to recur.

**A selection criterion contradicted the evaluation metric for months.** The kill/spare gate filtered on
potency while the target removes potency and the metric reads only ranking. Nothing flagged it; it took
asking *why is `nutlin-3` not in the panel*.
*Check that the selection criterion and the metric measure the same quantity.*

**Two unsourced thresholds set the entire drug panel.** The kill/spare cut-offs of 0.5 and 0.8 had nothing
behind them. Moving them to 0.7/0.8 produced a completely different panel of the same quality — visible
only once the arbitrariness was written down.
*An arbitrary threshold documented as arbitrary is honest; the same threshold stated without comment reads
as principled, and that is the more damaging of the two.*

**A scope change silently made a whole analysis thread irrelevant.** After the 8-drug panel was agreed,
work continued on per-drug variance weighting — estimating assay noise from replicates, deriving
reliability weights, sizing how many drugs are noise-dominated. All correct, all moot: that problem exists
only at K=545, and with eight comparable compounds the variance ratio is 2.5×.
*When the scope changes, re-check which problems the new scope has already dissolved before solving them.*

**Two different weightings were discussed under one word for several rounds.** One party meant weights
*per drug* (columns: some compounds spread wider than others), the other meant *per sample* (rows: some
response values common, some rare — class weighting for regression). The per-drug problem dissolves with
the panel; the per-sample one does not, and it was the one matching the observed failure mode.
*"Reweight the loss" is ambiguous — name the axis first.*

**A rollout plan argued over an option already discarded.** A proposed 1a/1b sequence kept `auc_z` as its
first step after it had been retired — and 1a could not have closed the z-score leak anyway, since μ and σ
are baked into the targets h5ad.
*Before proposing a sequence, check which of its steps the last decision already removed.*

**A broken implementation was invisible in the summary statistics and obvious in the plot.** With α = 1 and
a 10× cap, the weight curve saturated at the cap across whole regions, so the **cap** — an arbitrary safety
limit — set the weights rather than the label density. Plus a real bug: clipping then renormalizing pushed
values back over the cap (observed max **13.1** under a cap of 10). The same figure then revealed what no
amount of theory would have: inverse-density weighting up-weights the sparsest region of the response
range, which here is `auc > 1.1` — *the drug made the cells grow better than the control*. The scheme would
have handed maximum influence to the least trustworthy measurements in the dataset. Hence the
winsorization, which then removed nearly all the skew (`topotecan` +2.42 → +0.18) and reframed the whole
intervention from "repairing a broken distribution" to "deliberately re-emphasising the extremes" — in
hindsight predicting the null result.
*Render the figure and look at it before reporting anything based on it.*

**Two changes landed in one run and the result could not be attributed.** Diagnosing that violation cost
weeks in June 2026, and it is why the 27.07 step moved only the target and the loss.
*Never change the target and the architecture in the same run.*

**Three notebooks and two modules were written in a single sitting.** Their defects were found afterwards
by asking questions, not by reading the code — because a finished artefact arrives with every choice inside
it already made (bandwidth, cap, threshold, aggregation, axis, colour), and those are unreviewable once
buried in three hundred lines that run and produce a plausible figure.
*Work in pieces small enough to be read, questioned and rejected before the next begins.*

**Commits mixed authorship.** On 27–28.07.2026 several of Selin's own edits to `docs/TODO.md` and the
report were swept into commits carrying assistant-authored messages. Nothing was lost, but the history is
wrong about who wrote what, in a repository meant to be citable.
*Stage only the files changed in that piece of work; never `git add -A` or a whole directory.*

**Two claims were overstated and had to be walked back** rather than being hedged when written: that the
panel was chosen blind to our labels
([retraction](#the-panel-was-chosen-blind-to-our-labels)), and that proliferation likely explained the
cell-line effect ([refutation](#the-cell-line-effect-is-largely-proliferation)). Both were corrected within
a day, and the second is kept precisely because it was tested and killed rather than quietly dropped.
*State confidence at the level the evidence supports, the first time.*
