#!/usr/bin/env python3
"""Identify transcripts with largest accuracy gap: rigel-minimap2 vs others.

Find transcripts where salmon, kallisto, and rigel-oracle perform well
but rigel-minimap2 performs poorly. Rank by the delta between rigel-minimap2
error and the best-of-others error.
"""

import pandas as pd
import numpy as np
from pathlib import Path

V4 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none")


def main():
    ora = pd.read_csv(V4 / "per_transcript_counts_oracle.csv")
    mm2 = pd.read_csv(V4 / "per_transcript_counts_minimap2.csv")

    # Merge oracle (has rigel_oracle, salmon, kallisto) with minimap2
    m = ora[["transcript_id", "gene_id", "gene_name", "mrna_truth",
             "rigel_oracle", "salmon", "kallisto"]].merge(
        mm2[["transcript_id", "rigel_minimap2"]],
        on="transcript_id",
    )

    # Compute absolute errors for each tool
    m["err_oracle"] = np.abs(m["rigel_oracle"] - m["mrna_truth"])
    m["err_salmon"] = np.abs(m["salmon"] - m["mrna_truth"])
    m["err_kallisto"] = np.abs(m["kallisto"] - m["mrna_truth"])
    m["err_mm2"] = np.abs(m["rigel_minimap2"] - m["mrna_truth"])

    # Best-of-others: min error across oracle, salmon, kallisto
    m["err_best_other"] = m[["err_oracle", "err_salmon", "err_kallisto"]].min(axis=1)

    # Gap = mm2 error minus best-of-others error
    # Positive gap = mm2 is uniquely bad
    m["gap"] = m["err_mm2"] - m["err_best_other"]

    # Also compute signed error (for understanding direction)
    m["signed_mm2"] = m["rigel_minimap2"] - m["mrna_truth"]

    # Filter to transcripts where gap is meaningful
    # (mm2 is significantly worse than at least one other tool)
    candidates = m[m["gap"] > 10].sort_values("gap", ascending=False)

    print("=" * 140)
    print("TOP 50 TRANSCRIPTS: Largest rigel-minimap2 error gap vs best-of-others")
    print("=" * 140)
    print(f"{'rank':>4s}  {'transcript_id':>25s}  {'gene':>12s}  {'truth':>8s}  "
          f"{'oracle':>8s}  {'salmon':>8s}  {'kall':>8s}  {'mm2':>8s}  "
          f"{'err_mm2':>8s}  {'err_best':>8s}  {'gap':>8s}  {'dir':>5s}")
    for i, (_, r) in enumerate(candidates.head(50).iterrows()):
        direction = "over" if r.signed_mm2 > 0 else "under"
        print(f"  {i+1:>2d}  {r.transcript_id:>25s}  {str(r.gene_name):>12s}  {r.mrna_truth:>8.0f}  "
              f"{r.rigel_oracle:>8.0f}  {r.salmon:>8.1f}  {r.kallisto:>8.1f}  {r.rigel_minimap2:>8.0f}  "
              f"{r.err_mm2:>8.0f}  {r.err_best_other:>8.0f}  {r.gap:>+8.0f}  {direction:>5s}")

    # For the top 5, show gene-level context
    print("\n\n" + "=" * 140)
    print("GENE-LEVEL CONTEXT FOR TOP 5 GAP TRANSCRIPTS")
    print("=" * 140)

    top5 = candidates.head(5)
    for rank, (_, r) in enumerate(top5.iterrows(), 1):
        gene_id = r.gene_id
        gene_txs = m[m["gene_id"] == gene_id].sort_values("mrna_truth", ascending=False)
        gene_truth = gene_txs["mrna_truth"].sum()
        gene_mm2 = gene_txs["rigel_minimap2"].sum()
        gene_ora = gene_txs["rigel_oracle"].sum()
        gene_sal = gene_txs["salmon"].sum()
        gene_kal = gene_txs["kallisto"].sum()

        print(f"\n--- #{rank}: {r.transcript_id} ({r.gene_name}) | gap={r.gap:+.0f} ---")
        print(f"  Gene total: truth={gene_truth:.0f}  oracle={gene_ora:.0f}  "
              f"salmon={gene_sal:.0f}  kallisto={gene_kal:.0f}  mm2={gene_mm2:.0f}")
        print(f"  {'transcript_id':>30s}  {'truth':>8s}  {'oracle':>8s}  {'salmon':>8s}  "
              f"{'kall':>8s}  {'mm2':>8s}  {'err_mm2':>8s}  {'gap':>8s}")
        for _, tr in gene_txs.iterrows():
            print(f"  {tr.transcript_id:>30s}  {tr.mrna_truth:>8.0f}  {tr.rigel_oracle:>8.0f}  "
                  f"{tr.salmon:>8.1f}  {tr.kallisto:>8.1f}  {tr.rigel_minimap2:>8.0f}  "
                  f"{tr.err_mm2:>8.0f}  {tr.gap:>+8.0f}")

        # Check: which other genes have transcripts that GAINED counts in mm2?
        # (to find where the counts leaked TO or FROM)
        if r.signed_mm2 < 0:
            # This transcript LOST counts in mm2 — find where they went
            pass  # We'll analyze per-fragment

    # Summary statistics for the gap distribution
    print("\n\n" + "=" * 140)
    print("GAP DISTRIBUTION SUMMARY")
    print("=" * 140)
    print(f"  Total transcripts with gap > 0: {(m['gap'] > 0).sum()}")
    print(f"  Total transcripts with gap > 10: {(m['gap'] > 10).sum()}")
    print(f"  Total transcripts with gap > 100: {(m['gap'] > 100).sum()}")
    print(f"  Total transcripts with gap > 1000: {(m['gap'] > 1000).sum()}")
    print(f"  Sum of all gaps > 0: {m.loc[m['gap'] > 0, 'gap'].sum():.0f}")
    print(f"  Top 5 gaps account for: {top5['gap'].sum():.0f}")
    print(f"  Top 10 gaps account for: {candidates.head(10)['gap'].sum():.0f}")
    print(f"  Top 50 gaps account for: {candidates.head(50)['gap'].sum():.0f}")

    # Direction breakdown
    over = m[(m["gap"] > 100) & (m["signed_mm2"] > 0)]
    under = m[(m["gap"] > 100) & (m["signed_mm2"] < 0)]
    print(f"\n  Among gap > 100:")
    print(f"    Over-estimated (mm2 > truth): {len(over)}, total excess = {over['signed_mm2'].sum():.0f}")
    print(f"    Under-estimated (mm2 < truth): {len(under)}, total shortfall = {under['signed_mm2'].sum():.0f}")

    # Export top 5 for downstream analysis
    top5_ids = top5["transcript_id"].tolist()
    print(f"\n\nTOP 5 TRANSCRIPT IDS FOR EXTRACTION:")
    for tid in top5_ids:
        print(f"  {tid}")


if __name__ == "__main__":
    main()
