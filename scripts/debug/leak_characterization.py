"""Characterize the gDNA->RNA leak by transcript/locus TYPE (net_from_gdna per transcript).

Answers: where is the largest leak (absolute), what transcript types leak most (rate), and — the reframe —
which classes have ~zero leak. Breaks down by: mature vs nascent; single- vs multi-exon; AMBIG-locus
(opposite-strand overlap) vs clean single-strand; expression bucket; spliced length.

    python scripts/debug/leak_characterization.py [condition]
"""
import sys
import numpy as np, pandas as pd
pd.set_option('display.width', 220)

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
C = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
tx = pd.read_csv(f"{S}/net_flow_per_transcript.tsv", sep="\t")
tx = tx[tx.condition == C].copy()
tx["is_nrna"] = tx.transcript_id.str.startswith("RIGEL_NRNA_")

# AMBIG-locus flag: does the locus carry MATURE transcripts on BOTH strands? (opposite-strand overlap)
mat = tx[~tx.is_nrna]
strands = mat.groupby("locus_id").strand.agg(lambda s: set(s.dropna()))
ambig_locus = {lid: ({"+", "-"} <= st) for lid, st in strands.items()}
tx["ambig_locus"] = tx.locus_id.map(ambig_locus).fillna(False).astype(bool)

tot = tx.net_from_gdna.sum()
print(f"=== {C} ===")
print(f"TOTAL net_from_gdna (gDNA absorbed by RNA) = {tot:,.0f}   over {len(tx)} transcripts, {tx.locus_id.nunique()} loci\n")

def block(df, name):
    g = df.net_from_gdna
    absleak = g.sum()
    # relative rate: leak / observed mass (how hard is this class pulled per unit of its own signal)
    obs = df.observed.clip(lower=1).sum()
    n = len(df); nz = (g.abs() < 1).sum()
    print(f"  {name:34s} n={n:>5}  leak={absleak:>+11,.0f} ({100*absleak/tot:>5.1f}%)  "
          f"rate=leak/obs={absleak/max(obs,1):>+.4f}  |leak|>0 in {n-nz}/{n}")

print("--- A) MATURE vs NASCENT ---")
block(tx[~tx.is_nrna], "mature (mRNA)")
block(tx[tx.is_nrna], "nascent (nRNA)")

print("\n--- B) MATURE by exon structure ---")
m = tx[~tx.is_nrna]
block(m[m.single_exon == True], "mature single-exon")
block(m[(m.single_exon == False) & (m.n_exons <= 3)], "mature multi-exon (2-3 ex)")
block(m[(m.single_exon == False) & (m.n_exons > 3)], "mature multi-exon (>3 ex)")

print("\n--- C) MATURE by locus ambiguity (opposite-strand overlap) ---")
block(m[m.ambig_locus], "mature @ AMBIG locus (both strands)")
block(m[~m.ambig_locus], "mature @ clean single-strand locus")

print("\n--- D) MATURE by EXPRESSION (mrna_abundance quartile) ---")
m2 = m[m.mrna_abundance.notna()].copy()
m2["qexpr"] = pd.qcut(m2.mrna_abundance.rank(method="first"), 4, labels=["Q1 low", "Q2", "Q3", "Q4 high"])
for q in ["Q1 low", "Q2", "Q3", "Q4 high"]:
    block(m2[m2.qexpr == q], f"mature expr {q}")

print("\n--- E) MATURE by spliced length ---")
block(m[m.spliced_length < 2000], "mature splen <2kb")
block(m[(m.spliced_length >= 2000) & (m.spliced_length < 5000)], "mature splen 2-5kb")
block(m[m.spliced_length >= 5000], "mature splen >=5kb")

print("\n--- F) THE REFRAME: where is there ~NO leak? (|net_from_gdna| < 200, mature) ---")
clean = m[m.net_from_gdna.abs() < 200]; dirty = m[m.net_from_gdna >= 200]
print(f"  clean mature (|leak|<200): n={len(clean)}  mean_expr={clean.mrna_abundance.mean():,.0f}  "
      f"mean_nex={clean.n_exons.mean():.1f}  %single_exon={100*clean.single_exon.mean():.0f}  %AMBIGlocus={100*clean.ambig_locus.mean():.0f}")
print(f"  LEAKY mature (leak>=200): n={len(dirty)}  mean_expr={dirty.mrna_abundance.mean():,.0f}  "
      f"mean_nex={dirty.n_exons.mean():.1f}  %single_exon={100*dirty.single_exon.mean():.0f}  %AMBIGlocus={100*dirty.ambig_locus.mean():.0f}")

print("\n--- top-15 individual leak-absorbing transcripts ---")
top = tx.reindex(tx.net_from_gdna.sort_values(ascending=False).index).head(15)
print(top[["transcript_id", "locus_id", "net_from_gdna", "expected", "observed", "n_exons", "spliced_length", "strand", "single_exon", "is_nrna", "ambig_locus"]].to_string(index=False))
