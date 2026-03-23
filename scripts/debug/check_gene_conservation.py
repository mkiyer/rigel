#!/usr/bin/env python3
"""Check gene-level conservation for 50% undercount transcripts."""
import pandas as pd

BASE = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none"
df = pd.read_csv(f"{BASE}/per_transcript_counts_oracle.csv")

# Check ENSG00000285976 — 5338 FP fragments, 0 expressed isoforms
g = df[df["gene_id"] == "ENSG00000285976.4"]
print("ENSG00000285976 (ghost gene with 5338 FP):")
for _, r in g.iterrows():
    print(f"  {r.transcript_id}: truth={r.mrna_truth} rigel={r.rigel_oracle} "
          f"salmon={r.salmon:.1f} kallisto={r.kallisto:.1f}")

# Check PIGY — 335 FP, 0 expressed
g = df[df["gene_name"] == "PIGY"]
print(f"\nPIGY (335 FP):")
for _, r in g.iterrows():
    print(f"  {r.transcript_id}: truth={r.mrna_truth} rigel={r.rigel_oracle} "
          f"salmon={r.salmon:.1f} kallisto={r.kallisto:.1f}")

# Summary: how many of the 122 ~50% undercount transcripts have gene-level conservation?
expressed = df[df["mrna_truth"] > 0].copy()
rigel_rel_err = (expressed["rigel_oracle"] - expressed["mrna_truth"]) / expressed["mrna_truth"].clip(lower=1)
half = expressed[(rigel_rel_err < -0.4) & (rigel_rel_err > -0.6) & (expressed["mrna_truth"] > 100)]
gene_ids = half["gene_id"].unique()
print(f"\n122 ~50% undercount transcripts span {len(gene_ids)} genes")

conserved = 0
not_conserved = []
for gid in gene_ids:
    g = df[df["gene_id"] == gid]
    gene_diff = abs(g["rigel_oracle"].sum() - g["mrna_truth"].sum())
    if gene_diff < 5:
        conserved += 1
    else:
        not_conserved.append((gid, g["gene_name"].iloc[0], g["mrna_truth"].sum(), 
                              g["rigel_oracle"].sum(), gene_diff))

print(f"  Gene-level conserved (|gene_rigel - gene_truth| < 5): {conserved}/{len(gene_ids)}")

if not_conserved:
    print(f"\n  Non-conserved genes ({len(not_conserved)}):")
    for gid, gname, truth, rigel, diff in sorted(not_conserved, key=lambda x: -x[4]):
        print(f"    {gid} ({gname}): truth={truth:.0f} rigel={rigel:.0f} diff={diff:.0f}")

# Check: out of all 122 undercount transcripts, what fraction have an
# unexpressed isoform that got the "other half"?
print(f"\n## Checking unexpressed isoform siblings for 50%-undercount transcripts:")
n_has_fp_sibling = 0
for _, row in half.iterrows():
    siblings = df[(df["gene_id"] == row["gene_id"]) & (df["mrna_truth"] == 0)]
    fp_siblings = siblings[siblings["rigel_oracle"] > 10]
    if len(fp_siblings) > 0:
        n_has_fp_sibling += 1
print(f"  {n_has_fp_sibling}/{len(half)} have an unexpressed sibling with rigel > 10")

# What fraction of the ~50% undercount is explained by nRNA siphon vs isoform splitting?
print(f"\n## Cross-check: do these genes appear in multi-isoform loci?")
for _, row in half.head(10).iterrows():
    g = df[df["gene_id"] == row["gene_id"]]
    n_isoforms = len(g)
    n_expressed = (g["mrna_truth"] > 0).sum()
    n_fp = ((g["mrna_truth"] == 0) & (g["rigel_oracle"] > 0)).sum()
    gene_conserved = abs(g["rigel_oracle"].sum() - g["mrna_truth"].sum()) < 5
    print(f"  {row['transcript_id']:25s} ({row['gene_name']:12s}): "
          f"{n_isoforms} isoforms, {n_expressed} expressed, {n_fp} FP, "
          f"gene_conserved={gene_conserved}")
