#!/usr/bin/env python3
"""Diagnose the strand-channel failure from first principles.

For each VCaP sample, load the saved per-region feather and print:
  1. Fraction of regions at each n-bucket that have k_sense/n in
     {near 0, near 0.5, near 1}.  If SS ≈ 1, RNA should sit at 1,
     gDNA at 0.5, and antisense-biology at 0.  If the data is
     bimodal, the model should trivially separate it.
  2. The three-class empirical counts under *fixed* thresholds
     (k/n >= 0.8 → RNA, k/n <= 0.2 → antisense, else → gDNA).
  3. How many total fragments fall into each bucket (weight by n).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

CACHE = Path("/scratch/mkiyer_root/mkiyer0/shared_data/rigel_scratch/vcap_strand_v5/strand_eval")

N_BINS = [(1, 2), (2, 5), (5, 10), (10, 50), (50, 500), (500, 1_000_000)]


def nominal_fraction(sample: str) -> float:
    x = int(sample[3:].rstrip("m"))
    return x / (20 + x)


def classify(k, n, lo=0.2, hi=0.8):
    """Return (n_rna, n_gdna, n_anti) — counts of REGIONS by k/n."""
    p = k / np.maximum(n, 1.0)
    rna = p >= hi
    anti = p <= lo
    gdna = ~(rna | anti)
    return rna.sum(), gdna.sum(), anti.sum()


def frag_weight(k, n, lo=0.2, hi=0.8):
    """Sum of n over regions in each class — fragment-weighted."""
    p = k / np.maximum(n, 1.0)
    rna = p >= hi
    anti = p <= lo
    gdna = ~(rna | anti)
    return n[rna].sum(), n[gdna].sum(), n[anti].sum()


def analyze_sample(sample: str) -> dict:
    df = pd.read_feather(CACHE / f"regions_{sample}.feather")
    # Recover k_sense from tx_strand convention
    n_u = df["n_u"].to_numpy(dtype=np.float64)
    n_pos = df["n_pos"].to_numpy(dtype=np.float64)
    tx = df["tx_strand"].to_numpy()
    valid = (n_u >= 1) & (tx != 0) & df["eligible"].to_numpy()
    k = np.zeros_like(n_u)
    k[tx == 1] = n_u[tx == 1] - n_pos[tx == 1]
    k[tx == -1] = n_pos[tx == -1]
    k = k[valid]
    n = n_u[valid]

    print(f"\n=== {sample}  (nominal {nominal_fraction(sample):.3f}) "
          f"— N_valid={len(n):,}, Σn={int(n.sum()):,} ===")
    print(f"  mean k/n = {(k/np.maximum(n,1)).mean():.4f}")
    print(f"  median k/n = {np.median(k/np.maximum(n,1)):.4f}")

    # Fragment-weighted fraction of "gDNA-like" reads (p ≈ 0.5):
    # A region with tx_strand ≠ 0 under SS≈1 RNA should have k_sense ≈ n.
    # gDNA should have k_sense ≈ n/2.
    # So total ANTISENSE reads = Σ (n - k) ≈ 0.5·Σ_gDNA n + ε·Σ_RNA n
    # Under SS=0.9999, ε = 1e-4, so:
    #   n_antisense = 0.5 · n_gdna + 1e-4 · n_rna
    #   n_total = n_gdna + n_rna
    # Solve for n_gdna:
    #   n_gdna = (n_antisense - 1e-4·n_total) / (0.5 - 1e-4) ≈ 2·n_antisense
    # This is a MOMENT-OF-METHOD estimator, no BB mixture.
    total_sense = k.sum()
    total_n = n.sum()
    total_anti = total_n - total_sense
    # Method-of-moments: k ~ (1 - pi_g)·n·SS + pi_g·n·0.5
    # => pi_g = (SS - mean(k/n)) / (SS - 0.5)  for n-weighted mean
    # n-weighted mean of k/n = total_sense / total_n
    SS = 0.9997  # from dna40m; varies slightly per sample but ≈ 1
    mean_p = total_sense / total_n
    pi_g_mom = (SS - mean_p) / (SS - 0.5)
    print(f"  total_sense = {total_sense:.0f}  total_anti = {total_anti:.0f}")
    print(f"  fragment-weighted mean p(sense) = {mean_p:.4f}")
    print(f"  method-of-moments π_gDNA = {pi_g_mom:.4f}   "
          f"(nominal {nominal_fraction(sample):.3f})")

    # Three-class classification with hard thresholds
    r, g, a = classify(k, n)
    rf, gf, af = frag_weight(k, n)
    total_r = r + g + a
    total_f = rf + gf + af
    print(f"  regions:  RNA={r:,} ({r/total_r:.1%})  "
          f"gDNA={g:,} ({g/total_r:.1%})  anti={a:,} ({a/total_r:.1%})")
    print(f"  frags:    RNA={int(rf):,} ({rf/total_f:.1%})  "
          f"gDNA={int(gf):,} ({gf/total_f:.1%})  anti={int(af):,} ({af/total_f:.1%})")

    # Stratify by n
    print(f"  Stratified by n (fraction of regions in each bucket):")
    print(f"    n-bin     N_reg      %RNA    %gDNA    %anti   Σn       %RNA_f  %gDNA_f %anti_f")
    for lo, hi in N_BINS:
        mask = (n >= lo) & (n < hi)
        if not mask.any():
            continue
        kb, nb = k[mask], n[mask]
        r, g, a = classify(kb, nb)
        rf, gf, af = frag_weight(kb, nb)
        tot = r + g + a
        totf = rf + gf + af
        print(f"    {lo:>4}-{hi:<5} {tot:>8,}  {r/tot:>6.1%}  {g/tot:>6.1%}  "
              f"{a/tot:>6.1%}  {int(nb.sum()):>8,} {rf/totf:>6.1%}  "
              f"{gf/totf:>6.1%}  {af/totf:>6.1%}")

    return {
        "sample": sample,
        "nominal": nominal_fraction(sample),
        "mean_p": mean_p,
        "pi_g_mom": pi_g_mom,
        "total_sense": int(total_sense),
        "total_anti": int(total_anti),
        "frac_anti_frags": total_anti / total_n,
    }


def main():
    samples = ["dna00m", "dna01m", "dna02m", "dna05m",
               "dna10m", "dna20m", "dna40m", "dna80m"]
    rows = []
    for s in samples:
        rows.append(analyze_sample(s))
    df = pd.DataFrame(rows)
    print("\n\n=== SUMMARY (method-of-moments vs nominal) ===")
    print(df.to_string(index=False,
                       float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
