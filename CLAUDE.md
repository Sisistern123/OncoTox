# CLAUDE.md — how to work on this repository

## The hard rule: propose, then wait

**Do nothing before Selin has explicitly agreed to it.** Not code, not data, not docs, not notebooks, not
runs, not commits. Read and inspect freely — that is how a proposal is formed — but every action that
*changes* something or *decides* something waits for an explicit yes.

**No analysis decision is ever made silently.** This is the rule that matters most and the one that has
been broken. It covers every choice that determines what enters the model or how a number is computed:

- which drugs, cell lines, genes or samples are **selected, filtered, ranked or excluded** — and in what
  order those steps are applied
- **thresholds, cut-offs, winsorization, caps, exponents, bin widths, bandwidths**
- how anything is **aggregated** (per cell vs per cell line, pooled vs per fold, mean vs median)
- which **metric, baseline or target** a result is scored against
- how a **figure is computed or displayed**
- what a **loss weights**, and how

If a step feels like plumbing on the way to the thing that was asked for, that does not exempt it. A
selection step *is* an analysis decision. This is exactly how a compromised drug panel was produced on
27.07.2026 and a progress report had to be postponed: the candidate list was ranked on the project's own
response values before the literature criterion was applied, and that ordering was never mentioned.

## What a proposal looks like

Before acting, state in plain language:

1. **What** you intend to do, concretely enough to be vetoed in detail.
2. **What it affects** — which files, which numbers, which conclusions downstream.
3. **Why**, and **what the alternatives are**. If there is a real choice, name it as a choice rather than
   presenting one option as the obvious one.
4. **What would make it wrong** — the assumption it rests on.

Then wait. If Selin says "just do it", that permission covers the thing described, not the next thing
that turns out to be needed along the way.

## When you have changed something

**Explain every change.** Not a list of files touched — what changed, why, and what it means for the
numbers and conclusions that already exist. If a change invalidates something previously reported, say so
in the same message, unprompted.

## When something turns out to be wrong

**Name it as wrong, immediately and in full.** Do not grade a defect down into a caveat, and do not report
the mild version and wait to see whether anyone asks for the rest. If a mistake affects results, say which
results and how badly, in the message where you notice it.

Do not decide unilaterally that a fix should wait — least of all on the basis of a deadline that was never
stated. Whether to fix now or later is Selin's call.

## Working style

- **Verify rather than assert.** Numbers come from code that can be re-run. Anything that exists only as a
  shell command in a chat is not a result — it belongs in a notebook or script before it is reported.
- **Never change the target and the architecture in the same run.** One change at a time, or the outcome
  cannot be attributed. Diagnosing a violation of this cost the project weeks in June 2026.
- **Render figures and look at them** before reporting anything based on them. Two real defects on
  27.07.2026 were invisible in summary statistics and obvious in the plot.
- Prefer correcting the record over adding to it. Superseded conclusions get a dated marker, not a
  silent overwrite — dated log entries are a record of what was believed then.

## Where things live

`docs/TODO.md` — what is next and what blocks it; read the banner at the top first.
`docs/project_progress.md` — index: status table and the reasoned roadmap.
`docs/steps/01`–`08` — the scientific record, one file per stage.
`docs/project_notes.md` — dated decision log, newest first.
`report/` — the written version (LaTeX → `main.pdf`).
`docs/progress_report_*.md` — working record and slide text; **untracked by design**.
`notebooks/` — numbered in pipeline order; outputs under `notebooks/outputs/`.
`scripts/` — the pipeline: `preprocessing/`, `model/`, `training/`, `evaluation/`.
