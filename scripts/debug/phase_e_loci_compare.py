#!/usr/bin/env python3
"""Compare locus structure between old and new index.

Check: does the new architecture produce different loci sizes?
Are fragments routing differently?
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

BENCH_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune")
OUTPUT_DIR = BENCH_DIR / "phase_e_results"
NEW_INDEX_DIR = OUTPUT_DIR / "rigel_index"
OLD_INDEX_DIR = BENCH_DIR / "rigel_index"

# ss=1.00 condition
COND = "gdna_high_ss_1.00_nrna_default"
BAM = BENCH_DIR / COND / "align_oracle" / "reads_namesort.bam"


def run_and_report(label, index_dir, bam):
    from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig, PipelineConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import run_pipeline

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")

    idx = TranscriptIndex.load(str(index_dir))
    print(f"  Index: {idx.num_transcripts} transcripts")

    cfg = PipelineConfig(
        em=EMConfig(seed=42),
        scan=BamScanConfig(sj_strand_tag="auto", include_multimap=True),
        scoring=FragmentScoringConfig(),
    )
    t0 = time.monotonic()
    pipe = run_pipeline(bam, idx, config=cfg)
    elapsed = time.monotonic() - t0
    print(f"  Run: {elapsed:.1f}s ({pipe.stats.total} fragments)")

    est = pipe.estimator
    stats = pipe.stats

    # Loci
    loci_df = est.get_loci_df()
    sizes = loci_df["n_transcripts"].values
    em_frags = loci_df["n_em_fragments"].values
    gdna_vals = loci_df["gdna"].values

    print(f"\n  Loci summary:")
    print(f"    Total loci: {len(loci_df)}")
    print(f"    Total EM fragments: {em_frags.sum():.0f}")
    print(f"    Total gDNA (from loci): {gdna_vals.sum():.1f}")
    print(f"    gDNA from stats: {stats.n_gdna_total}")
    print(f"    Locus sizes: mean={sizes.mean():.1f}, median={np.median(sizes):.0f}, "
          f"max={sizes.max()}, p99={np.percentile(sizes, 99):.0f}")

    # Size distribution buckets
    for lo, hi in [(1, 1), (2, 5), (6, 10), (11, 50), (51, 100), (101, 500), (501, 1000), (1001, None)]:
        if hi is None:
            mask = sizes >= lo
            label_range = f">={lo}"
        else:
            mask = (sizes >= lo) & (sizes <= hi)
            label_range = f"{lo}-{hi}" if lo != hi else f"{lo}"
        n = mask.sum()
        em_in = em_frags[mask].sum()
        gdna_in = gdna_vals[mask].sum()
        print(f"    Size {label_range:>8s}: {n:6d} loci, {int(em_in):>10,} EM frags, {int(gdna_in):>10,} gDNA")

    # Pool counts
    syn_mask = idx.t_df["is_synthetic_nrna"].values if "is_synthetic_nrna" in idx.t_df.columns else np.zeros(idx.num_transcripts, dtype=bool)
    t_total = est.t_counts.sum(axis=1)
    mrna_total = float(t_total[~syn_mask].sum())
    nrna_total = float(t_total[syn_mask].sum())
    gdna_total = float(stats.n_gdna_total)
    print(f"\n  Pool totals:")
    print(f"    mRNA: {mrna_total:.0f}")
    print(f"    nRNA: {nrna_total:.0f}")
    print(f"    gDNA: {gdna_total:.0f}")

    return loci_df


def main():
    print("Loading old index and running old pipeline...")
    old_loci = run_and_report("BASELINE (pre-Phase-C index)", OLD_INDEX_DIR, BAM)

    print("\n\nLoading new index and running new pipeline...")
    new_loci = run_and_report("PHASE C (new index)", NEW_INDEX_DIR, BAM)


if __name__ == "__main__":
    main()
