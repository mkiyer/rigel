#!/usr/bin/env python3
"""Follow-up analysis: rigel siphon patterns and salmon 5% bias."""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none"
oracle_df = pd.read_csv(f"{BASE}/per_transcript_counts_oracle.csv")
mm2_df = pd.read_csv(f"{BASE}/per_transcript_counts_minimap2.csv")

expressed = oracle_df[oracle_df["mrna_truth"] > 0].copy()

print("=" * 80)
print("FOLLOW-UP ANALYSIS")
print("=" * 80)

# ── A. Salmon ~5% systematic undercount ──────────────────────
print("\n## A. Salmon Systematic Undercount (~5%)")
print("\n  Salmon loses ~5% of counts vs truth. Since this is an oracle-aligned")
print("  simulation, salmon is doing quasi-mapping from FASTQ and likely losing")
print("  multi-mapping reads or reads that don't quasi-map uniquely.")

# Check if it's uniform across expression bins
total_truth = expressed["mrna_truth"].sum()
total_salmon = expressed["salmon"].sum()
print(f"\n  Total truth:  {total_truth:,.0f}")
print(f"  Total salmon: {total_salmon:,.1f}")
print(f"  Ratio:        {total_salmon / total_truth:.4f}")

# Check mitochondrial transcripts specifically (salmon top errors are all MT-*)
mt_mask = expressed["gene_name"].str.startswith("MT-")
mt = expressed[mt_mask]
print(f"\n  Mitochondrial transcripts: {len(mt)}")
print(f"    Truth: {mt['mrna_truth'].sum():,.0f}")
print(f"    Salmon: {mt['salmon'].sum():,.1f}")
print(f"    Ratio: {mt['salmon'].sum() / mt['mrna_truth'].sum():.4f}")
print(f"    Salmon undercount on MT: {mt['mrna_truth'].sum() - mt['salmon'].sum():,.1f}")

nonmt = expressed[~mt_mask]
print(f"\n  Non-MT transcripts:")
print(f"    Truth: {nonmt['mrna_truth'].sum():,.0f}")
print(f"    Salmon: {nonmt['salmon'].sum():,.1f}")
print(f"    Ratio: {nonmt['salmon'].sum() / nonmt['mrna_truth'].sum():.4f}")

# ── B. Rigel LDHB and MDH1 etc FP analysis ──────────────────
print("\n\n## B. Rigel FP on Unexpressed Transcripts (LDHB, MDH1 siphon)")
unexpressed = oracle_df[oracle_df["mrna_truth"] == 0]
fp = unexpressed[unexpressed["rigel_oracle"] > 100].sort_values("rigel_oracle", ascending=False)
print(f"\n  Unexpressed transcripts with rigel > 100 fragments: {len(fp)}")
for _, row in fp.head(15).iterrows():
    # Check if other isoforms of same gene are expressed
    gene = row["gene_id"]
    same_gene = oracle_df[oracle_df["gene_id"] == gene]
    expressed_isoforms = same_gene[same_gene["mrna_truth"] > 0]
    total_gene_truth = same_gene["mrna_truth"].sum()
    total_gene_rigel = same_gene["rigel_oracle"].sum()
    print(f"\n    {row['transcript_id']} ({row['gene_name']}): rigel={row['rigel_oracle']:.0f}")
    print(f"      Gene {gene}: {len(same_gene)} isoforms, {len(expressed_isoforms)} expressed")
    print(f"      Gene truth: {total_gene_truth:.0f}, Gene rigel: {total_gene_rigel:.0f}")
    if len(expressed_isoforms) > 0:
        for _, iso in expressed_isoforms.iterrows():
            print(f"        {iso['transcript_id']}: truth={iso['mrna_truth']:.0f} "
                  f"rigel={iso['rigel_oracle']:.0f} salmon={iso['salmon']:.1f}")

# ── C. Rigel ~50% error pattern on top failures ─────────────
print("\n\n## C. Rigel ~50% Undercount Pattern")
print("  Several top rigel errors show exactly ~50% undercounting.")
print("  This is suspicious — could indicate EM splitting between mRNA and nRNA.\n")

# Check the top rigel failures for this pattern
rigel_err = expressed["rigel_oracle"] - expressed["mrna_truth"]
rigel_rel_err = rigel_err / expressed["mrna_truth"].clip(lower=1)
half_pattern = expressed[
    (rigel_rel_err < -0.4) & (rigel_rel_err > -0.6) & (expressed["mrna_truth"] > 100)
].copy()
half_pattern["rigel_rel_err"] = rigel_rel_err.loc[half_pattern.index]
half_pattern = half_pattern.sort_values("mrna_truth", ascending=False)
print(f"  Transcripts with rigel ~50% undercount (40-60% underprediction, truth > 100): {len(half_pattern)}")
for _, row in half_pattern.head(20).iterrows():
    print(f"    {row['transcript_id']:25s} ({row['gene_name']:15s}): "
          f"truth={row['mrna_truth']:8.0f}  rigel={row['rigel_oracle']:8.0f}  "
          f"ratio={row['rigel_oracle']/row['mrna_truth']:.3f}")

# ── D. Gene-level accuracy ──────────────────────────────────
print("\n\n## D. Gene-Level Accuracy (summing isoforms)")
gene_truth = expressed.groupby("gene_id").agg({
    "mrna_truth": "sum", "rigel_oracle": "sum", "salmon": "sum",
    "kallisto": "sum", "gene_name": "first"
}).reset_index()
gene_truth = gene_truth[gene_truth["mrna_truth"] > 0]

for tool in ["rigel_oracle", "salmon", "kallisto"]:
    err = gene_truth[tool] - gene_truth["mrna_truth"]
    mae = err.abs().mean()
    print(f"\n  {tool} gene-level:")
    print(f"    MAE:  {mae:.3f}")
    print(f"    Total: {gene_truth[tool].sum():,.1f} (truth: {gene_truth['mrna_truth'].sum():,.0f})")
    r_p, _ = pearsonr(gene_truth["mrna_truth"], gene_truth[tool])
    r_s, _ = spearmanr(gene_truth["mrna_truth"], gene_truth[tool])
    print(f"    Pearson:  {r_p:.6f}")
    print(f"    Spearman: {r_s:.6f}")

# ── E. Quantify nRNA siphon ─────────────────────────────────
print("\n\n## E. Quantifying Rigel's nRNA/gDNA Siphon")
# The pool-level summary showed 3661 nRNA + 1840 gDNA.
# At transcript level, rigel lost ~11,545 fragments from truth.
# But also has 40,886 FP fragments on unexpressed.
# Net: lost 11,545 from truth PLUS 40,886 on unexpressed = some redistribution

total_rigel_expressed = expressed["rigel_oracle"].sum()
total_rigel_all = oracle_df["rigel_oracle"].sum()
total_rigel_unexpressed = oracle_df[oracle_df["mrna_truth"] == 0]["rigel_oracle"].sum()
print(f"  Total rigel (all transcripts): {total_rigel_all:,.0f}")
print(f"  Total rigel (expressed only): {total_rigel_expressed:,.0f}")
print(f"  Total rigel (unexpressed FP): {total_rigel_unexpressed:,.0f}")
print(f"  Truth total: {total_truth:,.0f}")
print(f"  Deviation from truth (all): {total_rigel_all - total_truth:+,.0f}")
print(f"  nRNA+gDNA siphon from benchmark log: nRNA=3661 + gDNA=1840 = 5501")
print(f"  So ~{100*(total_truth - total_rigel_all)/total_truth:.3f}% of fragments attributed away from mRNA")

# ── F. Minimap2 analysis ────────────────────────────────────
print("\n\n## F. Minimap2 Alignment Degradation")
mm2_expressed = mm2_df[mm2_df["mrna_truth"] > 0].copy()
print(f"  Total truth: {mm2_expressed['mrna_truth'].sum():,.0f}")
print(f"  Total rigel_minimap2: {mm2_expressed['rigel_minimap2'].sum():,.1f}")
print(f"  Loss: {mm2_expressed['mrna_truth'].sum() - mm2_expressed['rigel_minimap2'].sum():,.1f}")
print(f"  Loss%: {100*(mm2_expressed['mrna_truth'].sum() - mm2_expressed['rigel_minimap2'].sum())/mm2_expressed['mrna_truth'].sum():.2f}%")

# The benchmark log showed: 17,685,478 buffered (vs 10M input) 
# and 2,815,188 multimapper molecules + nRNA=215,943 + gDNA=3,005
print("\n  From benchmark log:")
print("    Buffered fragments: 17,685,478 (multi-mapping inflates this)")
print("    Multimapper molecules: 2,815,188")
print("    nRNA attributed: 215,943")
print("    gDNA attributed: 3,005")
print("    The nRNA siphon is MUCH worse with minimap2 (215K vs 3.6K with oracle)")
print("    This is because minimap2 generates secondary alignments for multimappers")
print("    which create more ambiguity, allowing nRNA components to absorb more fragments")

# ── G. Rigel error decomposition ────────────────────────────
print("\n\n## G. Rigel Error Decomposition")
# Errors come from: (1) nRNA siphon, (2) gDNA siphon, (3) EM redistribution, (4) FP on unexpressed
print("  Sources of error in rigel (oracle alignment):")
print(f"    1. nRNA siphon:    ~3,661 frags (from pool-level)")
print(f"    2. gDNA siphon:    ~1,840 frags (from pool-level)")
print(f"    3. Isoform misallocation (FP on unexpressed): ~{total_rigel_unexpressed:,.0f} frags")
print(f"    4. Net mRNA tracking loss: {total_truth - total_rigel_all:,.0f} frags")
print(f"  Note: FP on unexpressed = EM redistributing from expressed to unexpressed isoforms")
print(f"  of the same gene. This is a classic EM identifiability problem, not a bug.")

print("\n" + "=" * 80)
print("END OF FOLLOW-UP ANALYSIS")
print("=" * 80)
