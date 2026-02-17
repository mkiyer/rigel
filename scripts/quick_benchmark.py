#!/usr/bin/env python3
"""Quick benchmark: run hulkrna with a given VBEM prior on all 10 regions."""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.index import HulkIndex
from hulkrna.pipeline import run_pipeline


def load_truth(per_tx_csv):
    truth, salmon, kallisto = {}, {}, {}
    with open(per_tx_csv) as f:
        for row in csv.DictReader(f):
            tid = row["transcript_id"]
            truth[tid] = int(row["truth"])
            salmon[tid] = float(row["salmon"])
            kallisto[tid] = float(row["kallisto"])
    return truth, salmon, kallisto


def mae(truth, pred):
    return np.mean([abs(truth[t] - pred.get(t, 0.0)) for t in truth])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench-dir", type=Path, required=True)
    parser.add_argument("--prior", type=float, default=0.5)
    parser.add_argument("--gdna-threshold", type=float, default=0.0)
    parser.add_argument("--overlap-min-frac", type=float, default=0.99,
                        help="Min fraction of best overlap to keep (0.99=within 1%% of best)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    regions = [
        "EGFR", "FGFR2", "HBB_cluster", "CD44", "BRCA1",
        "HOXA_cluster", "chr12_GAPDH_cluster", "chr19_dense",
        "chr22_CRKL_cluster", "chr1_dense",
    ]

    print(f"\n{'Region':<25} {'hulkrna':>10} {'salmon':>10} {'kallisto':>10} {'ratio_s':>10}")
    print("-" * 65)

    hulk_maes, sal_maes, kal_maes = [], [], []
    for region in regions:
        reg_dir = args.bench_dir / region
        if not reg_dir.exists():
            print(f"{region:<25} MISSING")
            continue

        truth, salmon_c, kallisto_c = load_truth(reg_dir / "per_transcript_counts.csv")
        index = HulkIndex.load(reg_dir / "hulkrna_index")
        bam = reg_dir / "align" / "reads.bam"

        pipe = run_pipeline(
            bam, index, seed=args.seed, sj_strand_tag="auto",
            gdna_threshold=args.gdna_threshold,
            em_pseudocount=args.prior,
            overlap_min_frac=args.overlap_min_frac,
        )
        counts_df = pipe.counter.get_counts_df(index)
        pred = {r.transcript_id: float(r.count) for r in counts_df.itertuples(index=False)}

        h = mae(truth, pred)
        s = mae(truth, salmon_c)
        k = mae(truth, kallisto_c)
        hulk_maes.append(h)
        sal_maes.append(s)
        kal_maes.append(k)
        print(f"{region:<25} {h:>10.1f} {s:>10.1f} {k:>10.1f} {h/s:>10.1f}x")

    print("-" * 65)
    print(f"{'MEAN':<25} {np.mean(hulk_maes):>10.1f} {np.mean(sal_maes):>10.1f} "
          f"{np.mean(kal_maes):>10.1f} {np.mean(hulk_maes)/np.mean(sal_maes):>10.1f}x")


if __name__ == "__main__":
    main()
