#!/usr/bin/env python3
"""Focused investigation of VBEM zero-forcing bug.

Key questions:
1. Which loci contain zero-forced highly-expressed transcripts?
2. How many transcripts/components do those loci have?
3. Which transcripts in those loci received the stolen mass?
4. Is the mega-locus the source of the problem?
"""

import pandas as pd
import numpy as np
from pathlib import Path

BENCHMARK_DIR = Path("/scratch/mkiyer_root/mkiyer0/shared_data/rigel_benchmarks/ccle_vcap_prostate")
CONDITION = "gdna_high_ss_0.90_nrna_none"

def main():
    vbem_dir = BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "vbem"
    map_dir = BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "map"

    vbem_q = pd.read_feather(vbem_dir / "quant.feather")
    map_q = pd.read_feather(map_dir / "quant.feather")
    vbem_loci = pd.read_feather(vbem_dir / "loci.feather")
    map_loci = pd.read_feather(map_dir / "loci.feather")

    # Load truth
    truth = pd.read_csv(BENCHMARK_DIR / "truth_abundances_nrna_none.tsv", sep="\t")

    print(f"VBEM quant columns: {list(vbem_q.columns)}")
    print(f"VBEM loci columns: {list(vbem_loci.columns)}")
    print(f"Truth columns: {list(truth.columns)}")

    # Check if locus_id is in quant
    has_locus = "locus_id" in vbem_q.columns
    print(f"\nlocus_id in quant: {has_locus}")

    # Compute TPM from truth
    truth_tpm = truth.copy()
    total_mrna = truth_tpm["mrna_abundance"].sum()
    truth_tpm["truth_tpm"] = truth_tpm["mrna_abundance"] / total_mrna * 1e6

    # Merge VBEM quant with truth 
    merged = vbem_q.merge(truth_tpm[["transcript_id", "truth_tpm", "mrna_abundance", "n_exons", "spliced_length", "genomic_span"]], 
                          on="transcript_id", how="left")
    merged["truth_tpm"] = merged["truth_tpm"].fillna(0)

    # Find zero-forced highly expressed transcripts
    zero_forced = merged[(merged["truth_tpm"] > 100) & (merged["tpm"] < 0.1)]
    n_zf = len(zero_forced)
    total_truth_tpm_lost = zero_forced["truth_tpm"].sum()
    print(f"\n{'='*80}")
    print(f"ZERO-FORCED TRANSCRIPTS (truth_tpm > 100, vbem < 0.1)")
    print(f"{'='*80}")
    print(f"Count: {n_zf}")
    print(f"Total truth TPM of zero-forced: {total_truth_tpm_lost:.0f}")
    print(f"  (that's {total_truth_tpm_lost/1e6*100:.1f}% of total TPM)")

    if has_locus:
        # Group zero-forced by locus
        zf_per_locus = zero_forced.groupby("locus_id").agg(
            n_zero_forced=("transcript_id", "count"),
            total_lost_tpm=("truth_tpm", "sum"),
        ).sort_values("total_lost_tpm", ascending=False)

        print(f"\nLoci with zero-forced transcripts: {len(zf_per_locus)}")
        print(f"\nTop 20 loci by lost TPM:")
        for locus_id, row in zf_per_locus.head(20).iterrows():
            locus_info = vbem_loci[vbem_loci["locus_id"] == locus_id]
            locus_txs = merged[merged["locus_id"] == locus_id]
            n_txs_in_locus = len(locus_txs)
            n_nonzero_truth = (locus_txs["truth_tpm"] > 0).sum()

            # Check locus details
            locus_cols = ""
            if len(locus_info) > 0:
                li = locus_info.iloc[0]
                for c in locus_info.columns:
                    if c != "locus_id":
                        locus_cols += f", {c}={li[c]}"

            print(f"\n  Locus {locus_id}: "
                  f"zero-forced={row['n_zero_forced']}, "
                  f"lost_tpm={row['total_lost_tpm']:.0f}, "
                  f"total_txs={n_txs_in_locus}, "
                  f"expressed={n_nonzero_truth}"
                  f"{locus_cols}")

            # Show transcripts in this locus (top 10 by truth)
            locus_txs_sorted = locus_txs.nlargest(10, "truth_tpm")
            for _, tx in locus_txs_sorted.iterrows():
                map_tx = map_q[map_q["transcript_id"] == tx["transcript_id"]]
                map_tpm = map_tx["tpm"].values[0] if len(map_tx) > 0 else "?"
                print(f"    {tx['transcript_id']}: truth={tx['truth_tpm']:.1f}, "
                      f"vbem={tx['tpm']:.2f}, map={map_tpm if isinstance(map_tpm, str) else f'{map_tpm:.2f}'}, "
                      f"n_exons={tx.get('n_exons', '?')}")

        # THE BIG QUESTION: Is this all in the mega-locus?
        print(f"\n{'='*80}")
        print("MEGA-LOCUS ANALYSIS")
        print(f"{'='*80}")

        # How many loci have > 100 transcripts?
        locus_sizes = merged.groupby("locus_id").size().sort_values(ascending=False)
        print(f"\nLocus size distribution:")
        print(f"  Total loci: {len(locus_sizes)}")
        print(f"  Max locus size: {locus_sizes.max()}")
        for thresh in [10, 50, 100, 500, 1000, 5000, 10000]:
            n_above = (locus_sizes >= thresh).sum()
            if n_above > 0:
                print(f"  Loci with >= {thresh} transcripts: {n_above}")

        # Zero-forced in mega-loci vs non-mega
        mega_threshold = 100
        mega_locus_ids = locus_sizes[locus_sizes >= mega_threshold].index
        zf_in_mega = zero_forced[zero_forced["locus_id"].isin(mega_locus_ids)]
        zf_not_mega = zero_forced[~zero_forced["locus_id"].isin(mega_locus_ids)]
        print(f"\nZero-forced in mega-loci (>={mega_threshold} txs): {len(zf_in_mega)}, "
              f"lost TPM: {zf_in_mega['truth_tpm'].sum():.0f}")
        print(f"Zero-forced in small loci (<{mega_threshold} txs): {len(zf_not_mega)}, "
              f"lost TPM: {zf_not_mega['truth_tpm'].sum():.0f}")

        # THE BIGGEST LOCUS
        biggest_locus = locus_sizes.index[0]
        biggest_txs = merged[merged["locus_id"] == biggest_locus]
        print(f"\nBiggest locus ({biggest_locus}): {len(biggest_txs)} transcripts")
        biggest_expressed = biggest_txs[biggest_txs["truth_tpm"] > 0]
        biggest_zf = biggest_txs[(biggest_txs["truth_tpm"] > 1) & (biggest_txs["tpm"] < 0.1)]
        print(f"  Expressed (truth > 0): {len(biggest_expressed)}")
        print(f"  Zero-forced (truth > 1, vbem < 0.1): {len(biggest_zf)}")
        print(f"  Total truth TPM in locus: {biggest_txs['truth_tpm'].sum():.0f}")
        print(f"  Total VBEM TPM in locus: {biggest_txs['tpm'].sum():.0f}")

        # Compare to MAP for the biggest locus
        biggest_map = map_q[map_q["locus_id"] == biggest_locus] if "locus_id" in map_q.columns else pd.DataFrame()
        if len(biggest_map) > 0:
            print(f"  Total MAP TPM in locus: {biggest_map['tpm'].sum():.0f}")
            map_zf = biggest_map.merge(truth_tpm[["transcript_id", "truth_tpm"]], on="transcript_id")
            map_zf = map_zf[(map_zf["truth_tpm"] > 1) & (map_zf["tpm"] < 0.1)]
            print(f"  MAP zero-forced (truth > 1, map < 0.1): {len(map_zf)}")

    # Also: check over-estimation
    print(f"\n{'='*80}")
    print("OVER-ESTIMATED TRANSCRIPTS (truth=0, vbem > 10)")
    print(f"{'='*80}")
    overest = merged[(merged["truth_tpm"] == 0) & (merged["tpm"] > 10)]
    print(f"Count: {len(overest)}")
    print(f"Total false TPM: {overest['tpm'].sum():.0f}")
    if has_locus:
        overest_per_locus = overest.groupby("locus_id").agg(
            n_overest=("transcript_id", "count"),
            total_false_tpm=("tpm", "sum"),
        ).sort_values("total_false_tpm", ascending=False)
        print(f"Loci with over-est: {len(overest_per_locus)}")
        overest_in_mega = overest[overest["locus_id"].isin(mega_locus_ids)]
        print(f"Over-est in mega-loci: {len(overest_in_mega)}, "
              f"false TPM: {overest_in_mega['tpm'].sum():.0f}")

    # single-exon vs multi-exon
    if "n_exons" in merged.columns:
        print(f"\n{'='*80}")
        print("SINGLE-EXON vs MULTI-EXON ZERO-FORCING")
        print(f"{'='*80}")
        se_zf = zero_forced[zero_forced["n_exons"] == 1]
        me_zf = zero_forced[zero_forced["n_exons"] > 1]
        print(f"Single-exon zero-forced: {len(se_zf)}, lost TPM: {se_zf['truth_tpm'].sum():.0f}")
        print(f"Multi-exon zero-forced: {len(me_zf)}, lost TPM: {me_zf['truth_tpm'].sum():.0f}")


if __name__ == "__main__":
    main()
