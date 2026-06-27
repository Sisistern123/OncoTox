"""Derive pipeline file locations from (data_root, variant) — the only place that encodes directory layout."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

# Default locations for this machine (override with --data-root / CLI flags if needed).
DEFAULT_DATA_ROOT = Path("/Users/selin/Desktop/OncoTox/data")
DEFAULT_SCGPT_SCRIPT = Path("/Users/selin/PycharmProjects/scGPT/gen_embeds.py")
DEFAULT_SCGPT_MODEL_DIR = Path("/Users/selin/Desktop/OncoTox/scGPT/scGPT_human")

VARIANTS = ("hvg1000", "hvg2000", "hvg3000", "hvg5000", "all_genes")
DEFAULT_VARIANT = "hvg5000"

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

    def __post_init__(self) -> None:
        if self.variant not in VARIANTS:
            raise ValueError(f"variant must be one of {VARIANTS}, got {self.variant!r}")

    @property
    def expr_file(self) -> Path:
        return self.data_root / "scRNAseq_SCP542" / "expression" / "CPM_data.txt"

    @property
    def meta_file(self) -> Path:
        return self.data_root / "scRNAseq_SCP542" / "metadata" / "Metadata.txt"

    @property
    def ctrp_dir(self) -> Path:
        return self.data_root / "metadata" / "CTRPv2.0_2015_ctd2_ExpandedDataset"

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
        return self.processed_dir / H5AD_TARGETS

    @classmethod
    def build(cls, data_root: Path | str | None, variant: str = DEFAULT_VARIANT) -> PipelinePaths:
        return cls(resolve_data_root(data_root), variant)


def add_data_args(
    parser: argparse.ArgumentParser,
    *,
    variant_default: str = DEFAULT_VARIANT,
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


def guard_output(path: Path, *, overwrite: bool, step: str) -> None:
    """Refuse to clobber an existing artifact unless ``--overwrite`` is set."""
    if path.exists() and not overwrite:
        raise SystemExit(
            f"[{step}] Output already exists (refusing to overwrite):\n  {path}\n"
            f"Use --overwrite to replace it, or pick another --variant."
        )
