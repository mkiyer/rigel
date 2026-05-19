"""Inspect loci.feather and locus_stats.feather summaries for two runs."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m")


def describe(label: str, run_dir: Path) -> None:
    loci = pd.read_feather(run_dir / "loci.feather")
    print(f"\n=== {label} :: loci.feather ({len(loci):,} rows) ===")
    print("columns:", list(loci.columns))
    # span = end - start if available
    if "locus_span_bp" in loci.columns:
        print("\nlocus_span_bp percentiles:")
        print(loci["locus_span_bp"].describe(percentiles=[0.5, 0.9, 0.99, 0.999]).to_string())
    if "n_em_fragments" in loci.columns:
        print("\nn_em_fragments percentiles:")
        print(loci["n_em_fragments"].describe(percentiles=[0.5, 0.9, 0.99, 0.999]).to_string())
    print("\ntop 10 by locus_span_bp:")
    cols = ["locus_id", "locus_span_bp", "n_transcripts", "n_nrna_entities", "n_em_fragments", "gdna", "mrna", "nrna"]
    cols = [c for c in cols if c in loci.columns]
    print(loci.nlargest(10, "locus_span_bp")[cols].to_string(index=False))
    print("\ntop 10 by n_em_fragments:")
    print(loci.nlargest(10, "n_em_fragments")[cols].to_string(index=False))


def main() -> None:
    for label, sub in [("with_mm", "with_mm"), ("no_mm", "no_mm")]:
        describe(label, ROOT / sub)


if __name__ == "__main__":
    main()
