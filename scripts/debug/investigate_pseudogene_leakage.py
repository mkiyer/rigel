#!/usr/bin/env python3
"""Investigate whether the FL change in v4 affects pseudogene leakage.

The hypothesis: transcript-space FL computation changes the fragment length 
assigned to reads that land on pseudogenes (intronless copies), which changes
EM resolution between real genes and their pseudogenes.

Look at the EIF3C/EIF3CL case in detail:
- v3: EIF3CL got ~3 counts (correct: truth=33)
- v4: EIF3CL got 1751 counts (wrong: truth=33), and EIF3C lost ~1723

Also look at the RPL genes (RPL5, RPL7A, RPL9, etc.) which are known to have
many processed pseudogenes.
"""

import pandas as pd
import numpy as np
from pathlib import Path

V3 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3/gdna_none_ss_0.95_nrna_none")
V4 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none")


def main():
    v3_mm2 = pd.read_csv(V3 / "per_transcript_counts_minimap2.csv")
    v4_mm2 = pd.read_csv(V4 / "per_transcript_counts_minimap2.csv")
    v4_ora = pd.read_csv(V4 / "per_transcript_counts_oracle.csv")

    # Merge all
    m = v3_mm2[["transcript_id", "gene_id", "gene_name", "mrna_truth", "rigel_minimap2"]].merge(
        v4_mm2[["transcript_id", "rigel_minimap2"]],
        on="transcript_id", suffixes=("_v3", "_v4"),
    ).merge(
        v4_ora[["transcript_id", "rigel_oracle", "salmon", "kallisto"]],
        on="transcript_id",
    )

    # Pattern: gene-level leakage between real gene and retrotransposed paralog
    # In v4, counts leak FROM real genes TO pseudogenes
    #
    # Is it possible that the FL computation change made pseudogene single-exon 
    # transcripts score more similarly to multi-exon real transcripts?

    # Check: genes with largest GAIN in v4 that have truth=0 or very low truth
    m["err_v3"] = np.abs(m["rigel_minimap2_v3"] - m["mrna_truth"])
    m["err_v4"] = np.abs(m["rigel_minimap2_v4"] - m["mrna_truth"])
    m["delta"] = m["err_v4"] - m["err_v3"]
    m["gained_v4"] = m["rigel_minimap2_v4"] - m["rigel_minimap2_v3"]

    # Transcripts that gained counts in v4 (false positives)
    false_gainers = m[(m["mrna_truth"] < 10) & (m["gained_v4"] > 100)].sort_values("gained_v4", ascending=False)

    print("=" * 120)
    print("TRANSCRIPTS WITH LOW TRUTH THAT GAINED >100 COUNTS IN V4")
    print("=" * 120)
    print(f"{'transcript_id':>30s}  {'gene':>15s}  {'truth':>8s}  {'oracle':>8s}  {'mm2_v3':>8s}  {'mm2_v4':>8s}  {'gained':>8s}  {'salmon':>8s}  {'kall':>8s}")
    for _, r in false_gainers.head(30).iterrows():
        print(f"  {r.transcript_id:>28s}  {str(r.gene_name):>15s}  {r.mrna_truth:>8.0f}  {r.rigel_oracle:>8.0f}  {r.rigel_minimap2_v3:>8.0f}  {r.rigel_minimap2_v4:>8.0f}  {r.gained_v4:>+8.0f}  {r.salmon:>8.1f}  {r.kallisto:>8.1f}")

    # Total false positive counts
    fp_v3 = m[m["mrna_truth"] == 0]["rigel_minimap2_v3"].sum()
    fp_v4 = m[m["mrna_truth"] == 0]["rigel_minimap2_v4"].sum()
    print(f"\n  Total false positive counts (truth=0):")
    print(f"    v3: {fp_v3:.0f}")
    print(f"    v4: {fp_v4:.0f}")
    print(f"    delta: {fp_v4 - fp_v3:+.0f}")

    # Check: are these pseudogenes? Look at gene_name patterns
    print("\n\n" + "=" * 120)
    print("PSEUDOGENE PATTERN ANALYSIS")
    print("=" * 120)

    # Common pseudogene naming patterns: ends in 'P' + number, 
    # or gene family with 'L' suffix (like EIF3CL)
    import re
    pseudo_pattern = re.compile(r'P\d+$|^LOC\d+|^ENSG\d+')

    for _, r in false_gainers.head(30).iterrows():
        gene = str(r.gene_name)
        is_pseudo = bool(pseudo_pattern.search(gene))
        # Find the "parent" gene by looking at gene families
        parent_candidates = m[(m["gene_name"].str.startswith(gene[:3], na=False)) & 
                            (m["mrna_truth"] > 1000)]["gene_name"].unique()[:5]
        print(f"  {r.transcript_id:>30s}  gene={gene:>15s}  pseudo={is_pseudo}  truth={r.mrna_truth:>6.0f}  v4_mm2={r.rigel_minimap2_v4:>6.0f}  parents={list(parent_candidates)}")

    # Overall: compute total inter-gene leakage for "gene family clusters"
    # where one member has high expression and another has low/zero expression
    print("\n\n" + "=" * 120)
    print("TOTAL LEAKAGE TO LOW-TRUTH TRANSCRIPTS")
    print("=" * 120)
    
    for threshold in [0, 1, 5, 10]:
        low = m["mrna_truth"] <= threshold
        v3_leak = m.loc[low, "rigel_minimap2_v3"].sum()
        v4_leak = m.loc[low, "rigel_minimap2_v4"].sum()
        ora_leak = m.loc[low, "rigel_oracle"].sum()
        sal_leak = m.loc[low, "salmon"].sum()
        kal_leak = m.loc[low, "kallisto"].sum()
        print(f"  truth<={threshold:>2d}: v3={v3_leak:>10.0f}  v4={v4_leak:>10.0f}  delta={v4_leak-v3_leak:>+8.0f}  oracle={ora_leak:>10.0f}  salmon={sal_leak:>10.0f}  kallisto={kal_leak:>10.0f}")

    # How much TOTAL mRNA was assigned correctly?
    print(f"\n  Total assigned (truth>0):")
    active = m["mrna_truth"] > 0
    print(f"    Truth:          {m.loc[active, 'mrna_truth'].sum():>12.0f}")
    print(f"    Oracle:         {m.loc[active, 'rigel_oracle'].sum():>12.0f}")
    print(f"    mm2 v3:         {m.loc[active, 'rigel_minimap2_v3'].sum():>12.0f}")
    print(f"    mm2 v4:         {m.loc[active, 'rigel_minimap2_v4'].sum():>12.0f}")
    print(f"    Salmon:         {m.loc[active, 'salmon'].sum():>12.0f}")
    print(f"    Kallisto:       {m.loc[active, 'kallisto'].sum():>12.0f}")


if __name__ == "__main__":
    main()
