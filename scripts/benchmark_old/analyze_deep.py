#!/usr/bin/env python3
"""Deep analysis of benchmark results — focus on nRNA siphon and gDNA leakage."""
import pandas as pd
import numpy as np

golden = "scripts/benchmark/golden"

# Load baseline
df = pd.read_csv(f"{golden}/locus_simple_baseline/sweep_results.tsv", sep="\t")

print("=" * 74)
print("BASELINE ANALYSIS: nRNA SIPHON EFFECT")
print("=" * 74)

# Focus on runs where TA1 is expressed
expressed = df[df["TA1"] > 0].copy()
expressed["nrna_siphon_frac"] = np.where(
    expressed["nrna_expected"] > 0,
    expressed["nrna_abs_diff"] / expressed["nrna_expected"],
    0.0
)

# Group by NTA1 level
print("\nnRNA accuracy by NTA1 abundance level:")
print(f"{'NTA1':>6} {'n':>4} {'nrna_rel_err':>12} {'mRNA_err':>10} {'nrna_diff':>10}")
for nta1, grp in expressed.groupby("NTA1"):
    nrna_runs = grp[grp["nrna_expected"] > 0]
    if len(nrna_runs) == 0:
        print(f"{nta1:>6} {len(grp):>4} {'n/a':>12} {grp['total_mrna_rel_err'].mean():>10.3f} {'n/a':>10}")
    else:
        print(f"{nta1:>6} {len(nrna_runs):>4} {nrna_runs['nrna_rel_err'].mean():>12.3f} "
              f"{nrna_runs['total_mrna_rel_err'].mean():>10.3f} {nrna_runs['nrna_abs_diff'].mean():>10.0f}")

# When NTA1 >> TA1, the nRNA component siphons mRNA
print("\n\nNTA1/TA1 ratio effect on mRNA error:")
print(f"{'NTA1':>6} {'TA1':>6} {'ratio':>8} {'mRNA_err':>10} {'nRNA_diff':>10} {'gDNA_diff':>10}")
for _, row in expressed[expressed["gdna_fraction"] == 0.0].sort_values(["NTA1", "TA1"]).iterrows():
    if row["NTA1"] == 0:
        continue
    ratio = row["NTA1"] / max(row["TA1"], 1)
    print(f"{int(row['NTA1']):>6} {int(row['TA1']):>6} {ratio:>8.1f} "
          f"{row['total_mrna_rel_err']:>10.3f} {row['nrna_abs_diff']:>10.0f} {row['gdna_abs_diff']:>10.0f}")

# gDNA leakage: how much gDNA gets misassigned
print("\n\n" + "=" * 74)
print("GDNA LEAKAGE ANALYSIS")
print("=" * 74)
gdna_runs = expressed[expressed["gdna_expected"] > 0].copy()
gdna_runs["gdna_rel_err"] = gdna_runs["gdna_abs_diff"].abs() / gdna_runs["gdna_expected"]

print(f"\n{'ss':>4} {'gDNA_frac':>9} {'n':>4} {'gDNA_relErr':>11} {'mRNA_err':>10}")
for (ss, gf), grp in gdna_runs.groupby(["strand_specificity", "gdna_fraction"]):
    print(f"{ss:>4.1f} {gf:>9.1f} {len(grp):>4} {grp['gdna_rel_err'].mean():>11.3f} "
          f"{grp['total_mrna_rel_err'].mean():>10.3f}")

# EM prior sensitivity
print("\n\n" + "=" * 74)
print("EM PRIOR PSEUDOCOUNT SENSITIVITY")
print("=" * 74)
df_prior = pd.read_csv(f"{golden}/locus_simple_em_prior/sweep_results.tsv", sep="\t")
# Focus on runs with both mRNA and nRNA
both = df_prior[(df_prior["TA1"] > 0) & (df_prior["NTA1"] > 0)]
print(f"\n{'prior':>8} {'n':>4} {'mRNA_err':>10} {'nRNAerr':>10} {'gDNA_err':>10}")
for prior, grp in both.groupby("prior_pseudocount"):
    gdna_sub = grp[grp["gdna_expected"] > 0]
    gdna_re = (gdna_sub["gdna_abs_diff"].abs() / gdna_sub["gdna_expected"]).mean() if len(gdna_sub) > 0 else 0
    print(f"{prior:>8.1f} {len(grp):>4} {grp['total_mrna_rel_err'].mean():>10.3f} "
          f"{grp['nrna_rel_err'].mean():>10.3f} {gdna_re:>10.3f}")

# MAP vs VBEM
print("\n\n" + "=" * 74)
print("MAP vs VBEM COMPARISON")
print("=" * 74)
df_mode = pd.read_csv(f"{golden}/locus_simple_em_mode/sweep_results.tsv", sep="\t")
both = df_mode[(df_mode["TA1"] > 0) & (df_mode["NTA1"] > 0)]
print(f"\n{'mode':>6} {'n':>4} {'mRNA_err':>10} {'nRNAerr':>10}")
for mode, grp in both.groupby("mode"):
    print(f"{mode:>6} {len(grp):>4} {grp['total_mrna_rel_err'].mean():>10.3f} "
          f"{grp['nrna_rel_err'].mean():>10.3f}")

# Strand specificity sensitivity
print("\n\n" + "=" * 74)
print("STRAND SPECIFICITY SENSITIVITY")
print("=" * 74)
df_strand = pd.read_csv(f"{golden}/locus_simple_strand/sweep_results.tsv", sep="\t")
both = df_strand[(df_strand["TA1"] > 0) & (df_strand["NTA1"] > 0)]
print(f"\n{'ss':>4} {'n':>4} {'mRNA_err':>10} {'nRNAerr':>10}")
for ss, grp in both.groupby("strand_specificity"):
    print(f"{ss:>4.1f} {len(grp):>4} {grp['total_mrna_rel_err'].mean():>10.3f} "
          f"{grp['nrna_rel_err'].mean():>10.3f}")

# Scoring penalties
print("\n\n" + "=" * 74)
print("SCORING PENALTY SENSITIVITY")
print("=" * 74)
df_score = pd.read_csv(f"{golden}/locus_simple_scoring/sweep_results.tsv", sep="\t")
both = df_score[(df_score["TA1"] > 0) & (df_score["NTA1"] > 0)]
print(f"\n{'overhang':>10} {'mismatch':>10} {'n':>4} {'mRNA_err':>10} {'nRNAerr':>10}")
for (oh, mm), grp in both.groupby(["overhang_log_penalty", "mismatch_log_penalty"]):
    print(f"{oh:>10.3f} {mm:>10.3f} {len(grp):>4} {grp['total_mrna_rel_err'].mean():>10.3f} "
          f"{grp['nrna_rel_err'].mean():>10.3f}")
