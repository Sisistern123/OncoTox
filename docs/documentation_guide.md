# Documentation Guide

Two top-level markdown docs with deliberately different roles. Keep them in sync on
numbers, but know which one is authoritative.

## `project_progress.md` — source of truth
Comprehensive, chronological record. Must contain **all** important steps, numbers,
parameters, and results, so that everything is derivable from this file alone. Includes
plan-alignment callouts (✅ on-plan / ⚠️ deviation) against
`~/Desktop/OncoTox/project_plan/project_planning_v2.pdf`.

**When new work is done, update this file** with the hard details: cell/gene counts,
split distributions, hyperparameters, run IDs, MSEs, and any deviations from the plan.

## `project_notes.md` — thought / decision log
Free-form, dated entries where reasoning, advisor updates, ideas, and open questions are
reiterated. An *addition*, not the primary record. Mine it for context, but put
authoritative detail in `project_progress.md`.

## Conventions
- Lead with `project_progress.md` for any factual number; only echo short reasoning into
  `project_notes.md`.
- Keep numbers consistent between the two files.
- **Known cross-doc inconsistency to keep flagging:** SCP542×CTRPv2 cell-line overlap is
  **190** (audit notebook, case-insensitive) vs **180** (pipeline `ctrp_to_h5ad.py`, which
  also strips `-`). Same data, different normalization — pick one per context.
