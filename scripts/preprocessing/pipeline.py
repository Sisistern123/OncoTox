"""The preprocessing steps, one function each -- called in order by the numbered notebooks.

Each function owns the guard and the preconditions that belong to *its* step, so a step cannot be
run out of order or silently clobber an expensive artifact no matter who calls it. What it does
**not** own is the order itself: that lives in the notebooks, which are the pipeline
(`1_data` runs `fetch` and `convert`; `3_representations` runs the remaining four).

Replaces ``run_preprocessing.py``, archived 12.08.2026. That script held the same six steps behind
an ``argparse`` CLI with a ``STEP_ORDER`` list and ``--start-at``. The ordering half became the
notebooks' job when the pipeline was renumbered into five stages, and a second copy of the order in
a CLI is a second thing to keep in step with them. The guards, the preconditions and the scGPT
subprocess bridge are not plumbing and were kept -- they are the reason these are functions rather
than lines in a cell.

There is deliberately no ``__main__`` here: the notebooks are the entry point. To run the pipeline
without a browser, execute them headless (``jupyter nbconvert --execute``).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from scripts.layout import (
    DEFAULT_SCGPT_MODEL_DIR,
    DEFAULT_SCGPT_SCRIPT,
    VARIANT_N_TOP_GENES,
    PipelinePaths,
    guard_output,
)
from scripts.preprocessing import add_pca, create_splits, ctrp_to_h5ad, scp542_conversion
from scripts.sources import fetch_ctrp_response


def fetch(paths: PipelinePaths) -> Path:
    """Download and verify the pinned CTRPv2 response data. Returns the response CSV path.

    Idempotent, and the only step that touches the network: an archive already cached and matching
    its MD5 is neither re-downloaded nor re-extracted. It is first so that a run from an empty data
    directory reproduces the target with no manual step, and so the response data's version is
    fixed by :data:`scripts.layout.ZENODO_RESPONSE_RECORD` rather than by whatever happens to be
    on disk.
    """
    print(f"[fetch] Zenodo record {fetch_ctrp_response.ZENODO_RECORD}")
    fetch_ctrp_response.fetch_ctrp_response(paths.metadata_dir)
    if not paths.ctrp_response_csv.exists():
        raise RuntimeError(
            f"[fetch] Response table missing after fetch:\n  {paths.ctrp_response_csv}"
        )
    return paths.ctrp_response_csv


def convert(
    paths: PipelinePaths,
    *,
    n_top_genes: int | None = None,
    overwrite: bool = False,
) -> Path:
    """SCP542's text matrices -> the raw h5ad, gene-filtered for this variant. Returns ``raw_h5ad``.

    ``n_top_genes`` overrides the variant's own HVG count; leave it ``None`` to follow
    :data:`scripts.layout.VARIANT_N_TOP_GENES`, which is what every recorded run did. **The two
    ``None``\\ s mean different things** and the distinction is load-bearing: ``n_top_genes=None``
    means *take the count from the variant*, whereas a resolved count of ``None`` -- which is what
    the ``all_genes`` variant maps to, and what ``0`` is folded into -- means *apply no HVG filter
    at all*.

    Refuses to overwrite an existing ``raw_h5ad`` unless ``overwrite=True``. It is one of the two
    expensive artifacts, and the cost of a mistaken rebuild is not the minutes: every
    representation and target downstream is derived from this file, so replacing it silently
    invalidates them while leaving them on disk looking current.
    """
    resolved = VARIANT_N_TOP_GENES[paths.variant] if n_top_genes is None else n_top_genes
    hvg = resolved if resolved and resolved > 0 else None
    print(f"[convert] {paths.variant}: " + (f"top-{hvg} HVGs" if hvg else "no HVG filter"))

    guard_output(paths.raw_h5ad, overwrite=overwrite, step="convert")
    # run_preprocessing.py did this once in main() for every step; convert is now the first step
    # that writes into it, so it owns the mkdir.
    paths.processed_dir.mkdir(parents=True, exist_ok=True)
    scp542_conversion.run(
        str(paths.expr_file), str(paths.meta_file), str(paths.raw_h5ad), hvg
    )
    return paths.raw_h5ad


def scgpt(
    paths: PipelinePaths,
    scgpt_python: str | Path | None,
    *,
    model_dir: str | Path = DEFAULT_SCGPT_MODEL_DIR,
    script: str | Path = DEFAULT_SCGPT_SCRIPT,
    overwrite: bool = False,
) -> Path:
    """Embed the raw h5ad with scGPT, in scGPT's own interpreter. Returns ``embed_h5ad``.

    scGPT pins versions this project does not, so it lives in a separate virtualenv and this step
    has to be a subprocess rather than an import -- ``scgpt_python`` is that interpreter.

    **There is no interactive fallback.** ``run_preprocessing._run_scgpt`` printed the command and
    blocked on ``input()`` when no interpreter was given. That was a terminal affordance: run
    headless (``jupyter nbconvert --execute``), it hangs the kernel on an invisible prompt instead
    of failing. Without an interpreter this raises and prints the command to run by hand.

    Note what happens after running it by hand: the embedding now exists, so calling ``scgpt()``
    again **raises** rather than skipping -- ``guard_output`` refuses to overwrite. That is the
    intended behaviour; continue at :func:`targets` instead of re-running this step.
    """
    if not paths.raw_h5ad.exists():
        raise RuntimeError(
            f"[scgpt] Missing the convert output:\n  {paths.raw_h5ad}\n"
            f"Stage 1 (1_data) has not been run for variant {paths.variant!r}."
        )
    guard_output(paths.embed_h5ad, overwrite=overwrite, step="scgpt")

    script, model_dir = Path(script), Path(model_dir)
    if not script.exists():
        raise FileNotFoundError(f"[scgpt] Embedding script not found: {script}")

    manual = (
        f"  <scgpt-python> {script} --input {paths.raw_h5ad} "
        f"--output {paths.embed_h5ad} --model-dir {model_dir}"
    )
    if not scgpt_python:
        raise RuntimeError(
            "[scgpt] needs the scGPT virtualenv's interpreter; there is no fallback.\n"
            "Either pass scgpt_python=<path>, or run this yourself:\n"
            f"{manual}\n"
            "and then continue at targets(). Do not re-run scgpt() afterwards -- it will "
            "refuse, because the embedding it would write now exists."
        )

    cmd = [
        str(scgpt_python), str(script),
        "--input", str(paths.raw_h5ad),
        "--output", str(paths.embed_h5ad),
        "--model-dir", str(model_dir),
    ]
    print("[scgpt] " + " ".join(cmd))
    paths.embed_h5ad.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)

    if not paths.embed_h5ad.exists():
        raise RuntimeError(
            f"[scgpt] Subprocess reported success but wrote no output:\n  {paths.embed_h5ad}"
        )
    return paths.embed_h5ad


def targets(
    paths: PipelinePaths,
    *,
    min_cell_lines: int = ctrp_to_h5ad.DEFAULT_MIN_CELL_LINES,
) -> Path:
    """Map the CTRPv2 response measure onto cells. Returns ``targets_h5ad``.

    Writes ``obsm['Y_ctrp']``, ``obsm['M_ctrp']`` and ``uns['ctrp_drugs']`` for every drug screened
    against at least ``min_cell_lines`` SCP542-overlapping cell lines; ``0`` keeps every drug with
    any overlap at all. The measure is whichever ``paths.score`` names.

    **Deliberately panel-independent (Selin, 12.08.2026).** The drug panel from stage 2 is applied
    at *training* time, not here, so this file does not have to be rebuilt when the panel changes.
    The stronger reason is that :func:`splits` derives eligibility from "this line has at least one
    observed label": restricting the target matrix to the panel would recompute that set over 11
    drugs instead of ~545, and any line thereby dropped would silently re-freeze the split,
    retiring every result scored on the old partition.

    Writes **no** ``viability_<drug>`` / ``train_mask_<drug>`` columns -- that chain was dropped on
    12.08.2026 together with ``ScGPTDrugDataset``, its only consumer (see :func:`splits`).

    Unlike :func:`convert` and :func:`scgpt` this step is not guarded: it is cheap to redo and its
    inputs are immutable, so overwriting the targets h5ad costs minutes rather than a re-embed.
    """
    if not paths.embed_h5ad.exists():
        raise RuntimeError(
            f"[targets] Missing the embedding:\n  {paths.embed_h5ad}\nRun scgpt() first."
        )
    if not paths.ctrp_response_csv.exists():
        raise RuntimeError(
            f"[targets] Missing the response table:\n  {paths.ctrp_response_csv}\n"
            f"Run fetch() in stage 1 first."
        )

    print(f"[targets] score={paths.score}, min_cell_lines={min_cell_lines}")
    ctrp_to_h5ad.run(
        str(paths.embed_h5ad),
        str(paths.targets_h5ad),
        str(paths.ctrp_response_csv),
        min_cell_lines=min_cell_lines,
        extra_single_drug_cols=(),
        score=paths.score,
    )
    return paths.targets_h5ad


def splits(paths: PipelinePaths, *, seed: int = 42, regenerate: bool = False) -> Path:
    """Write the cell-line-grouped 70/15/15 ``split_ctrp`` column. Returns ``targets_h5ad``.

    A cell line is eligible if any of its cells carries at least one observed CTRP label, and the
    assignment is read from the frozen ``splits/split_ctrp.csv`` rather than redrawn -- eligibility
    depends on which drugs survive the filters, so redrawing would move lines between train, val
    and test whenever anything upstream changed, and runs from either side would look comparable
    while being scored on different held-out lines.

    ``regenerate=True`` redraws and overwrites that file. **Every result produced before it was
    regenerated was scored on different held-out lines**, so this is not a repair.

    Writes only ``split_ctrp``. The per-drug ``split_<drug>`` column went on 12.08.2026 with the
    rest of the single-drug chain (:func:`targets`): its consumer ``ScGPTDrugDataset`` had no
    caller, and it cost ``add_pca`` a second 512-component train-only decomposition on every run.
    """
    if not paths.targets_h5ad.exists():
        raise RuntimeError(
            f"[splits] Missing the targets h5ad:\n  {paths.targets_h5ad}\nRun targets() first."
        )
    print(f"[splits] seed={seed}, regenerate={regenerate}")
    create_splits.run_multi(str(paths.targets_h5ad), seed=seed, regenerate=regenerate)
    return paths.targets_h5ad


def pca(
    paths: PipelinePaths,
    *,
    n_comps: int = add_pca.DEFAULT_N_COMPS,
    seed: int = add_pca.DEFAULT_SEED,
    force: bool = False,
) -> Path:
    """Add the PCA baseline representation to the targets h5ad. Returns ``targets_h5ad``.

    Fitted on ``raw_h5ad`` -- the convert output, carrying the full HVG set -- and **not** on the
    targets file's own ``.X``, from which the scGPT step has dropped out-of-vocabulary genes. Both
    representations therefore rest on the same single HVG filter, which is what makes the
    PCA-versus-scGPT comparison a comparison of representations rather than of gene sets.

    Writes the all-cells ``X_pca`` (descriptive: correct for UMAPs, wrong as model input) and, for
    the fixed split, ``X_pca_train_ctrp`` fitted on that split's training cells only. Runs last
    because the train-only fit needs :func:`splits` to have written ``split_ctrp``.
    """
    for required, step in ((paths.targets_h5ad, "targets()"), (paths.raw_h5ad, "convert()")):
        if not required.exists():
            raise RuntimeError(f"[pca] Missing input:\n  {required}\nRun {step} first.")
    print(f"[pca] n_comps={n_comps}, seed={seed}, force={force}")
    add_pca.run(
        str(paths.targets_h5ad),
        force=force,
        counts_h5ad=str(paths.raw_h5ad),
        n_comps=n_comps,
        seed=seed,
    )
    return paths.targets_h5ad
