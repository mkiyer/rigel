"""Diagnostics for the Phase 2 Bayesian-prior synthetic-sim run."""
from pathlib import Path

import numpy as np
import pandas as pd

base = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic")

print("=== Per-locus calibration: gdna_prior_count vs realized gdna ===")
print(f"  {'Condition':<38s}  {'Pearson':>8s}  {'Slope':>7s}  {'Sum(a)/Sum(gdna)':>17s}  {'n':>4s}")
for cond in [
    "gdna_low_ss_0.99_nrna_none",
    "gdna_med_ss_0.99_nrna_none",
    "gdna_equal_ss_0.99_nrna_none",
    "gdna_high_ss_0.99_nrna_none",
]:
    loc = pd.read_feather(base / cond / "rigel_out" / "loci.feather")
    valid = loc[loc["n_em_fragments"] > 100]
    rho = float(np.corrcoef(valid["gdna_prior_count"], valid["gdna"])[0, 1])
    slope, intercept = np.polyfit(valid["gdna_prior_count"], valid["gdna"], 1)
    s_alpha = float(valid["gdna_prior_count"].sum())
    s_gdna = float(valid["gdna"].sum())
    print(f"  {cond:<38s}  {rho:>8.4f}  {slope:>7.3f}  {s_alpha/max(s_gdna,1):>17.3f}  {len(valid):>4d}")

print()
print("=== Sparse-locus pathology: gdna_prior_count > n_em_fragments ===")
for cond in ["gdna_low_ss_0.99_nrna_none", "gdna_high_ss_0.99_nrna_none"]:
    loc = pd.read_feather(base / cond / "rigel_out" / "loci.feather")
    loc = loc.copy()
    loc["ratio"] = loc["gdna_prior_count"] / loc["n_em_fragments"].clip(lower=1)
    over = loc[loc["ratio"] > 1.0].sort_values("ratio", ascending=False)
    print(f"\n  {cond}: {len(over)} loci with gdna_prior_count > n_em_fragments")
    if len(over):
        print(over[
            ["locus_id", "n_transcripts", "locus_span_bp", "n_em_fragments",
             "gdna_prior_count", "gdna", "mrna", "ratio"]
        ].head(8).to_string(index=False))

print()
print("=== nRNA false positives (true nRNA = 0 across ALL sims) ===")
print(f"  {'Condition':<38s}  {'nRNA_count':>10s}  {'%_total':>7s}  {'n_loci_with_nrna':>17s}")
for cond in [
    "gdna_none_ss_0.99_nrna_none", "gdna_none_ss_0.50_nrna_none",
    "gdna_low_ss_0.99_nrna_none",  "gdna_low_ss_0.50_nrna_none",
    "gdna_med_ss_0.99_nrna_none",  "gdna_med_ss_0.50_nrna_none",
    "gdna_equal_ss_0.99_nrna_none","gdna_equal_ss_0.50_nrna_none",
    "gdna_high_ss_0.99_nrna_none", "gdna_high_ss_0.50_nrna_none",
]:
    loc = pd.read_feather(base / cond / "rigel_out" / "loci.feather")
    nrna = float(loc["nrna"].sum())
    total = float((loc["mrna"] + loc["nrna"] + loc["gdna"]).sum())
    n_loci_w = int((loc["nrna"] > 1).sum())
    print(f"  {cond:<38s}  {nrna:>10.0f}  {100*nrna/max(total,1):>6.3f}%  {n_loci_w:>17d}")

print()
print("=== gDNA prior calibration efficiency by gDNA level ===")
print("  Σα = expected gDNA pseudocount across ALL loci")
print("  Σgdna = realized gDNA assignments")
print(f"  {'Condition':<38s}  {'Σα':>10s}  {'Σgdna':>10s}  {'Σα/Σgdna':>10s}  {'true_n_gdna':>12s}")
true_n = {"none": 0, "low": 100_000, "med": 500_000, "equal": 1_000_000, "high": 2_000_000}
for cond in [
    "gdna_none_ss_0.99_nrna_none",
    "gdna_low_ss_0.99_nrna_none",
    "gdna_med_ss_0.99_nrna_none",
    "gdna_equal_ss_0.99_nrna_none",
    "gdna_high_ss_0.99_nrna_none",
]:
    loc = pd.read_feather(base / cond / "rigel_out" / "loci.feather")
    s_alpha = float(loc["gdna_prior_count"].sum())
    s_gdna = float(loc["gdna"].sum())
    label = cond.split("_")[1]
    print(f"  {cond:<38s}  {s_alpha:>10.0f}  {s_gdna:>10.0f}  {s_alpha/max(s_gdna,1):>10.3f}  {true_n[label]:>12d}")

print()
print("=== False-positive transcript spotlight (gdna_high_ss_0.99) ===")
import json
truth = pd.read_csv(base / "truth_abundances.tsv", sep="\t")
quant = pd.read_feather(base / "gdna_high_ss_0.99_nrna_none" / "rigel_out" / "quant.feather")
m = truth.merge(quant[["transcript_id", "count", "tpm"]], on="transcript_id", how="left").fillna(0)
fp = m[m["mrna_abundance"] == 0].copy()
fp["gene_id_short"] = fp["gene_id"]
top = fp.nlargest(8, "count")
expressed = m[m["mrna_abundance"] > 0].copy()
print("  These unexpressed transcripts get high false counts in gdna_high condition.")
print("  Hypothesis: they belong to genes that ARE expressed via other isoforms.")
for _, row in top.iterrows():
    gid = row["gene_id"]
    siblings = expressed[expressed["gene_id"] == gid]
    sib_count = float(siblings["count"].sum())
    sib_true = float(siblings["mrna_abundance"].sum())
    n_sib = len(siblings)
    print(f"    {row['transcript_id']:<14s}  gene={gid:<12s}  fp_count={row['count']:>6.0f}  "
          f"n_sibling_isoforms={n_sib:2d}  sibling_true_total={sib_true:>9.1f}  sibling_obs_total={sib_count:>9.1f}")

print()
print("=== nRNA false positives — which loci? (gdna_high_ss_0.99) ===")
loc = pd.read_feather(base / "gdna_high_ss_0.99_nrna_none" / "rigel_out" / "loci.feather")
nrna_loci = loc[loc["nrna"] > 5].nlargest(10, "nrna")
print(nrna_loci[["locus_id","n_transcripts","n_nrna_entities","n_em_fragments","mrna","nrna","gdna","gdna_prior_count"]].to_string(index=False))
