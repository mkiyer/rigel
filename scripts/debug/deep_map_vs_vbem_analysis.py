#!/usr/bin/env python3
"""Deep analysis of MAP vs VBEM benchmark results.

Produces:
1. Overall wall-clock timing comparison  
2. Per-condition accuracy comparison (MAP vs VBEM vs salmon)
3. Mega-locus analysis: transcripts in mega-locus vs rest
4. Pool-level (mRNA/nRNA/gDNA) error analysis
5. Expression-bin stratified comparison
6. Per-transcript scatter for MAP vs VBEM agreement
"""
import pandas as pd
import numpy as np
import json
import os

BDIR = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks/runs"
RDIR = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks/results/map_vs_vbem"

CONDITIONS = ["gdna_none_ss_1.00_nrna_none", "gdna_high_ss_0.90_nrna_none"]
COND_LABELS = {"gdna_none_ss_1.00_nrna_none": "CLEAN", "gdna_high_ss_0.90_nrna_none": "DIRTY"}
TOOLS = ["rigel/map", "rigel/vbem", "rigel/default", "salmon"]

# ─── 1. Wall-clock timing ──────────────────────────────────────────
print("=" * 80)
print("1. WALL-CLOCK TIMING (from run.log)")
print("=" * 80)
with open(os.path.join(RDIR, "run.log")) as f:
    for line in f:
        if "completed in" in line or "Run Summary" in line or "Total runtime" in line:
            print(f"  {line.strip()}")
print()

# ─── 2. Transcript-level metrics side-by-side ──────────────────────
print("=" * 80)
print("2. TRANSCRIPT-LEVEL ACCURACY (MAP vs VBEM vs default vs salmon)")
print("=" * 80)
tm = pd.read_csv(os.path.join(RDIR, "transcript_metrics.csv"))
cols = ["condition", "tool", "pearson_r", "spearman_r", "rmse", "mae", "mape", "ware"]
for cond in CONDITIONS:
    label = COND_LABELS[cond]
    print(f"\n  [{label}] {cond}")
    sub = tm[tm["condition"] == cond][cols].copy()
    sub = sub[sub["tool"].isin(TOOLS)]
    sub = sub.sort_values("tool")
    for _, row in sub.iterrows():
        print(f"    {row['tool']:20s}  R={row['pearson_r']:.4f}  ρ={row['spearman_r']:.4f}  "
              f"RMSE={row['rmse']:.2f}  MAE={row['mae']:.3f}  MAPE={row['mape']:.1f}%  WARE={row['ware']:.4f}")

# ─── 3. Gene-level metrics side-by-side ────────────────────────────
print(f"\n{'=' * 80}")
print("3. GENE-LEVEL ACCURACY (MAP vs VBEM vs default vs salmon)")
print("=" * 80)
gm = pd.read_csv(os.path.join(RDIR, "gene_metrics.csv"))
cols_g = ["condition", "tool", "pearson_r", "spearman_r", "rmse", "mae", "mape", "ware"]
for cond in CONDITIONS:
    label = COND_LABELS[cond]
    print(f"\n  [{label}] {cond}")
    sub = gm[gm["condition"] == cond][cols_g].copy()
    sub = sub[sub["tool"].isin(TOOLS)]
    sub = sub.sort_values("tool")
    for _, row in sub.iterrows():
        print(f"    {row['tool']:20s}  R={row['pearson_r']:.4f}  ρ={row['spearman_r']:.4f}  "
              f"RMSE={row['rmse']:.2f}  MAE={row['mae']:.3f}  MAPE={row['mape']:.1f}%  WARE={row['ware']:.4f}")

# ─── 4. Pool-level summary ─────────────────────────────────────────
print(f"\n{'=' * 80}")
print("4. POOL-LEVEL SUMMARY (mRNA / nRNA / gDNA fragment counts)")
print("=" * 80)
ps = pd.read_csv(os.path.join(RDIR, "pool_summary.csv"))
for cond in CONDITIONS:
    label = COND_LABELS[cond]
    print(f"\n  [{label}] {cond}")
    sub = ps[ps["condition"] == cond].copy()
    sub = sub[sub["tool"].isin(["rigel/map", "rigel/vbem", "rigel/default"])]
    for _, row in sub.iterrows():
        tool = row["tool"]
        mrna_t = row["mrna_frag_truth"]
        mrna_p = row["mrna_pred"]
        mrna_e = row["mrna_rel_error"]
        gdna_t = row["gdna_frag_truth"]
        gdna_p = row["gdna_pred"]
        gdna_e = row["gdna_rel_error"] if pd.notna(row["gdna_rel_error"]) else "N/A"
        nrna_p = row["nrna_pred"]
        print(f"    {tool:20s}  mRNA: {mrna_p:,.0f}/{mrna_t:,.0f} (err={mrna_e:+.4f})  "
              f"gDNA: {gdna_p:,.0f}/{gdna_t:,.0f} (err={gdna_e})  "
              f"nRNA_leak: {nrna_p:,.0f}")

# ─── 5. Mega-locus analysis ────────────────────────────────────────
print(f"\n{'=' * 80}")
print("5. MEGA-LOCUS ANALYSIS")
print("=" * 80)

detail = pd.read_parquet(os.path.join(RDIR, "per_transcript_detail.parquet"))

for cond in CONDITIONS:
    label = COND_LABELS[cond]
    # Load loci to identify mega-locus transcripts
    # Use map run loci (same locus structure for map/vbem since same input)
    loci = pd.read_feather(f"{BDIR}/{cond}/rigel/map/loci.feather")
    
    # Find mega-locus (largest by n_transcripts)
    mega = loci.nlargest(1, "n_transcripts").iloc[0]
    mega_id = mega["locus_id"]
    mega_n_tx = mega["n_transcripts"]
    mega_n_frag = mega["n_em_fragments"]
    
    print(f"\n  [{label}] Mega-locus #{mega_id}: {mega_n_tx} transcripts, {mega_n_frag:,.0f} fragments")
    
    # Load quant to get locus_id per transcript
    for mode in ["map", "vbem"]:
        quant = pd.read_feather(f"{BDIR}/{cond}/rigel/{mode}/quant.feather")
        if "locus_id" in quant.columns:
            mega_tx = set(quant[quant["locus_id"] == mega_id]["transcript_id"].values)
            sub = detail[(detail["condition"] == cond) & (detail["tool"] == f"rigel/{mode}")]
            in_mega = sub[sub["transcript_id"].isin(mega_tx)]
            not_mega = sub[~sub["transcript_id"].isin(mega_tx)]
            
            # Accuracy for mega vs non-mega
            for label2, df in [("MEGA-LOCUS", in_mega), ("OTHER LOCI", not_mega)]:
                expressed = df[df["truth_tpm"] > 0]
                if len(expressed) > 0:
                    corr = np.corrcoef(expressed["truth_tpm"], expressed["predicted"])[0, 1]
                    mae = (expressed["truth_tpm"] - expressed["predicted"]).abs().mean()
                    ware = ((expressed["truth_tpm"] - expressed["predicted"]).abs().sum() /
                            expressed["truth_tpm"].sum())
                    print(f"    {mode:6s} {label2:12s}: n_tx={len(df):6d} n_expressed={len(expressed):6d}  "
                          f"R={corr:.4f}  MAE={mae:.3f}  WARE={ware:.4f}")
        else:
            print(f"    {mode}: no locus_id in quant.feather, skipping mega-locus breakdown")

# ─── 6. MAP vs VBEM per-transcript agreement ───────────────────────
print(f"\n{'=' * 80}")
print("6. MAP vs VBEM PER-TRANSCRIPT AGREEMENT")
print("=" * 80)

for cond in CONDITIONS:
    label = COND_LABELS[cond]
    map_sub = detail[(detail["condition"] == cond) & (detail["tool"] == "rigel/map")][
        ["transcript_id", "predicted", "truth_tpm"]
    ].set_index("transcript_id")
    vbem_sub = detail[(detail["condition"] == cond) & (detail["tool"] == "rigel/vbem")][
        ["transcript_id", "predicted"]
    ].set_index("transcript_id")
    
    merged = map_sub.join(vbem_sub, rsuffix="_vbem", how="inner")
    merged.columns = ["map_pred", "truth", "vbem_pred"]
    merged = merged.fillna(0)
    
    # Correlation between MAP and VBEM predictions
    corr_mv = np.corrcoef(merged["map_pred"], merged["vbem_pred"])[0, 1]
    
    # Mean absolute difference between MAP and VBEM
    mad = (merged["map_pred"] - merged["vbem_pred"]).abs().mean()
    
    # Max absolute difference
    max_diff = (merged["map_pred"] - merged["vbem_pred"]).abs().max()
    max_diff_idx = (merged["map_pred"] - merged["vbem_pred"]).abs().idxmax()
    
    # Fraction where MAP is closer to truth
    expressed = merged[merged["truth"] > 0].copy()
    map_err = (expressed["map_pred"] - expressed["truth"]).abs()
    vbem_err = (expressed["vbem_pred"] - expressed["truth"]).abs()
    map_wins = (map_err < vbem_err).sum()
    vbem_wins = (map_err > vbem_err).sum()
    ties = (map_err == vbem_err).sum()
    
    print(f"\n  [{label}] {cond}")
    print(f"    MAP-VBEM correlation: {corr_mv:.6f}")
    print(f"    Mean |MAP - VBEM|: {mad:.4f}")
    print(f"    Max  |MAP - VBEM|: {max_diff:.2f} (transcript: {max_diff_idx})")
    print(f"    Expressed transcripts: MAP closer to truth: {map_wins} ({100*map_wins/len(expressed):.1f}%)  "
          f"VBEM closer: {vbem_wins} ({100*vbem_wins/len(expressed):.1f}%)  ties: {ties}")

# ─── 7. Stratified MAP vs VBEM ─────────────────────────────────────
print(f"\n{'=' * 80}")
print("7. STRATIFIED MAP vs VBEM (by expression bins)")
print("=" * 80)

sm = pd.read_csv(os.path.join(RDIR, "stratified_metrics.csv"))
for cond in CONDITIONS:
    label = COND_LABELS[cond]
    print(f"\n  [{label}] {cond}")
    for mode in ["rigel/map", "rigel/vbem"]:
        sub = sm[(sm["condition"] == cond) & (sm["tool"] == mode)]
        print(f"    {mode}:")
        for _, row in sub.iterrows():
            bin_name = row.get("expression_bin", "?")
            pr = row.get("pearson_r", float("nan"))
            ware = row.get("ware", float("nan"))
            mae = row.get("mae", float("nan"))
            pr_str = f"{pr:.4f}" if pd.notna(pr) else "N/A"
            ware_str = f"{ware:.4f}" if pd.notna(ware) else "N/A"
            mae_str = f"{mae:.3f}" if pd.notna(mae) else "N/A"
            print(f"      {bin_name:20s}  R={pr_str}  MAE={mae_str}  WARE={ware_str}")

# ─── 8. nRNA leakage breakdown ─────────────────────────────────────
print(f"\n{'=' * 80}")
print("8. nRNA LEAKAGE (truth=0 nRNA, how much leaks into nRNA component)")
print("=" * 80)

for cond in CONDITIONS:
    label = COND_LABELS[cond]
    for mode in ["map", "vbem"]:
        summary = json.load(open(f"{BDIR}/{cond}/rigel/{mode}/summary.json"))
        q = summary["quantification"]
        total = q["mrna_total"] + q["nrna_total"] + q["gdna_total"]
        print(f"  [{label}] {mode:6s}: nRNA = {q['nrna_total']:>10,.0f} "
              f"({100*q['nrna_total']/total:.3f}% of genic)  "
              f"gDNA = {q['gdna_total']:>12,.0f}")

print(f"\n{'=' * 80}")
print("ANALYSIS COMPLETE")
print("=" * 80)
