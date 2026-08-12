"""Derive pipeline file locations from (data_root, variant) — the only place that encodes directory layout."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

# Default locations for this machine (override with --data-root / CLI flags if needed).
DEFAULT_DATA_ROOT = Path("/Users/selin/Desktop/OncoTox/data")
# Vendored next to this file (03.08.2026); it still needs the separate scGPT venv, which
# run_preprocessing.py reaches via --scgpt-python.
# layout.py sits at scripts/ (it is the path contract for training and evaluation too, not only for
# preprocessing), so the embedding script is one directory down rather than a sibling. Moved
# 12.08.2026; `with_name` here silently pointed at a non-existent scripts/gen_embeds.py.
DEFAULT_SCGPT_SCRIPT = Path(__file__).resolve().parent / "preprocessing" / "gen_embeds.py"
DEFAULT_SCGPT_MODEL_DIR = Path("/Users/selin/Desktop/OncoTox/scGPT/scGPT_human")

VARIANTS = ("hvg1000", "hvg2000", "hvg3000", "hvg5000", "all_genes")
DEFAULT_VARIANT = "hvg5000"

#: Zenodo record holding DrEval's reprocessed CTRPv2 response data and the Cellosaurus release it
#: ships with. **Pinned deliberately**: drevalpy's own loader resolves the concept DOI to whatever is
#: latest, so the version it returns changes under you. Bumping this is a target change -- everything
#: downstream has to be re-run, and `docs/steps/01` has to record the new record and retrieval date.
#: Record 21807175 = "Dataset for drevalpy", published 2026-08-05, retrieved 11.08.2026.
ZENODO_RESPONSE_RECORD = 21807175
#: Cache directory carries the record, so an older copy is visible rather than silently overwritten.
FETCH_CACHE_DIRNAME = f"drevalpy_CTRPv2_zenodo_{ZENODO_RESPONSE_RECORD}"

# CTRPv2 response measure used as the training target. Both come from DrEval's re-fit of CTRPv2's raw
# dose-response data with CurveCurator, and both are columns of its published CTRPv2.csv:
#   auc_cc     : `AUC_curvecurator`     -- area under the fitted curve, viability normalised per
#                replicate against the no-drug control. Complete for every curve.
#   ln_ic50_cc : `LN_IC50_curvecurator` -- natural log of the IC50 from the same fit. Missing for
#                ~40 % of curves by construction: DrEval discard an IC50 that falls more than an
#                order of magnitude outside the measured dose range, which is most compounds that
#                never reach half-killing. CTRPv2 itself publishes no IC50 at all.
CTRP_SCORES = ("auc_cc", "ln_ic50_cc")
# The `_cc` suffix is load-bearing. Until 11.08.2026 the default was `auc`, meaning CTRP's published
# `area_under_curve` divided by `conc_pts_fit` -- a defective normalisation (the divisor counts the
# points that survived outlier censoring, not the width of the integral), so every number computed on
# it is void. Reusing the name for CurveCurator's AUC would have made the old and new artifacts
# indistinguishable on disk. Both `auc` and the legacy `mean_pv` were removed with their reader code
# on 11.08.2026 (Selin). See docs/steps/corrections-and-dead-ends.md.
DEFAULT_CTRP_SCORE = "auc_cc"

# hvg1000/2000/3000 added for the HVG-count sweep (find scGPT's filtering sweet spot).
VARIANT_N_TOP_GENES: dict[str, int | None] = {
    "hvg1000": 1000,
    "hvg2000": 2000,
    "hvg3000": 3000,
    "hvg5000": 5000,
    "all_genes": None,
}

H5AD_RAW = "SCP542_CCLE.h5ad"
H5AD_EMBED = "SCP542_CCLE_scGPT_human_embeddings.h5ad"
H5AD_TARGETS = "SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"


def targets_filename(score: str) -> str:
    """Targets h5ad name for a response measure -- one file per measure, never shared.

    The measure is in the filename because it is the only thing distinguishing two otherwise
    identical h5ads, and a target you cannot identify from the path is one that gets mixed into a
    comparison by accident.
    """
    return H5AD_TARGETS.replace(".h5ad", f"_{score}.h5ad")


def resolve_data_root(explicit: Path | str | None = None) -> Path:
    """Resolve the OncoTox data directory (must exist).

    Uses ``explicit`` when passed (e.g. ``--data-root``); otherwise ``DEFAULT_DATA_ROOT``.
    """
    root = Path(explicit).expanduser().resolve() if explicit is not None else DEFAULT_DATA_ROOT
    if not root.is_dir():
        raise SystemExit(f"Data root does not exist or is not a directory: {root}")
    return root


@dataclass(frozen=True)
class PipelinePaths:
    """All inputs/outputs for one preprocessing variant (hvg5000 or all_genes)."""

    data_root: Path
    variant: str
    score: str = DEFAULT_CTRP_SCORE

    def __post_init__(self) -> None:
        if self.variant not in VARIANTS:
            raise ValueError(f"variant must be one of {VARIANTS}, got {self.variant!r}")
        if self.score not in CTRP_SCORES:
            raise ValueError(f"score must be one of {CTRP_SCORES}, got {self.score!r}")

    @property
    def expr_file(self) -> Path:
        return self.data_root / "scRNAseq_SCP542" / "expression" / "CPM_data.txt"

    @property
    def meta_file(self) -> Path:
        return self.data_root / "scRNAseq_SCP542" / "metadata" / "Metadata.txt"

    @property
    def metadata_dir(self) -> Path:
        """Where third-party reference data is cached, one directory per pinned release."""
        return self.data_root / "metadata"

    @property
    def ctrp_dir(self) -> Path:
        """CTRPv2's own 2015 distribution.

        **No longer the source of the training target** (see ``drevalpy_dir``). Kept because it is
        still the only source of CTRP's compound metadata -- ``v20.meta.per_compound.txt`` carries the
        ``broad_cpd_id`` the cross-database step needs -- and several analysis notebooks read it.
        """
        return self.metadata_dir / "CTRPv2.0_2015_ctd2_ExpandedDataset"

    @property
    def drevalpy_dir(self) -> Path:
        """DrEval's reprocessed CTRPv2, pinned to one Zenodo record and fetched by the ``fetch`` step."""
        return self.metadata_dir / FETCH_CACHE_DIRNAME

    @property
    def ctrp_response_csv(self) -> Path:
        """The response table the target is built from: one row per (cell line, drug, experiment)."""
        return self.drevalpy_dir / "CTRPv2" / "CTRPv2.csv"

    @property
    def cellosaurus_file(self) -> Path:
        """Cellosaurus release shipped with the same record -- pinned by it, since Cellosaurus
        publishes no stable per-release download URL of its own."""
        return self.drevalpy_dir / "meta" / "cellosaurus.txt"

    @property
    def processed_dir(self) -> Path:
        return self.data_root / "processed" / "scRNAseq_SCP542" / self.variant

    @property
    def raw_h5ad(self) -> Path:
        return self.processed_dir / H5AD_RAW

    @property
    def embed_h5ad(self) -> Path:
        return self.processed_dir / H5AD_EMBED

    @property
    def targets_h5ad(self) -> Path:
        return self.processed_dir / targets_filename(self.score)

    @classmethod
    def build(
        cls,
        data_root: Path | str | None,
        variant: str = DEFAULT_VARIANT,
        score: str = DEFAULT_CTRP_SCORE,
    ) -> PipelinePaths:
        return cls(resolve_data_root(data_root), variant, score)


def add_data_args(
    parser: argparse.ArgumentParser,
    *,
    variant_default: str = DEFAULT_VARIANT,
    score_default: str = DEFAULT_CTRP_SCORE,
) -> None:
    parser.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help=f"OncoTox data directory (default: {DEFAULT_DATA_ROOT}).",
    )
    parser.add_argument(
        "--variant",
        choices=VARIANTS,
        default=variant_default,
        help="Gene-set variant; outputs go to processed/scRNAseq_SCP542/<variant>/.",
    )
    parser.add_argument(
        "--score",
        choices=CTRP_SCORES,
        default=score_default,
        help=(
            "CTRPv2 response score used as the target. Each score gets its own targets "
            "h5ad, so auc and mean_pv runs can be compared head-to-head."
        ),
    )


def guard_output(path: Path, *, overwrite: bool, step: str) -> None:
    """Refuse to clobber an existing artifact unless ``--overwrite`` is set."""
    if path.exists() and not overwrite:
        raise SystemExit(
            f"[{step}] Output already exists (refusing to overwrite):\n  {path}\n"
            f"Use --overwrite to replace it, or pick another --variant."
        )
