#!/usr/bin/env python3
"""Detailed diagnostic: trace every fragment, decompose errors, check for regressions.

Questions to answer:
1. Is fragment accounting internally consistent (no leak/double-count)?
2. Where does alignment loss vs model misclassification explain the gap?
3. Did the refactoring change pool-level results vs old runs (rigel/default)?
4. Are annotated nRNA equivalents creating TPM contamination?
5. Is the analysis comparing the right columns?
"""
import pandas as pd
import numpy as np
import json
import os

BDIR = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks"
RUNS = os.path.join(BDIR, "runs")

manifest = json.load(open(os.path.join(BDIR, "manifest.json")))
conditions_meta = {c["name"]: c for c in manifest["conditions"]}

truth_none = pd.read_csv(os.path.join(BDIR, "truth_abundances_nrna_none.tsv"), sep="\t")

CONDITIONS = ["gdna_none_ss_1.00_nrna_none", "gdna_high_ss_0.90_nrna_none"]

print("=" * 100)
print("DIAGNOSTIC 1: Fragment accounting consistency")
print("=" * 100)

for cond in CONDITIONS:
    meta = conditions_meta[cond]
    n_rna_truth = meta["n_rna"]
    n_gdna_truth = meta["n_gdna"]
    total_truth = n_rna_truth + n_gdna_truth

    for mode in ["vbem", "map", "default"]:
        outdir = os.path.join(RUNS, cond, "rigel", mode)
        if not os.path.exists(os.path.join(outdir, "summary.json")):
            continue
        summary = json.load(open(os.path.join(outdir, "summary.json")))
        q = summary["quantification"]
        fstats = summary["fragment_stats"]
        astats = summary["alignment_stats"]

        total_detected = fstats["total"]
        alignment_loss = total_truth - total_detected
        genic = fstats["genic"]
        intergenic = fstats["intergenic"]
        chimeric = fstats["chimeric"]

        # Pool totals from summary
        mrna = q["mrna_total"]
        nrna = q["nrna_total"]
        gdna = q["gdna_total"]
        inter = q["intergenic_total"]
        pool_total = mrna + nrna + gdna + inter
        unaccounted = genic - (mrna + nrna + gdna)

        print(f"\n  [{cond[:20]:>20s}] {mode:>7s}:")
        print(f"    truth={total_truth:>12,}  detected={total_detected:>12,}  "
              f"alignment_loss={alignment_loss:>10,} ({100*alignment_loss/total_truth:.2f}%)")
        print(f"    genic={genic:>12,}  intergenic={intergenic:>8,}  chimeric={chimeric:>6,}")
        print(f"    mrna={mrna:>12,.0f}  nrna={nrna:>10,.0f}  gdna_em={gdna:>10,.0f}  "
              f"inter={inter:>10,}  pool_total={pool_total:>12,.0f}")
        print(f"    genic - (mrna+nrna+gdna) = {unaccounted:,.0f} (should be ≈0)")

print(f"\n{'=' * 100}")
print("DIAGNOSTIC 2: Old (default/star) vs New (vbem/map) comparison")
print("=" * 100)

for cond in CONDITIONS:
    print(f"\n  [{cond}]")
    rows = []
    for mode in ["default", "star", "vbem", "map", "oracle"]:
        outdir = os.path.join(RUNS, cond, "rigel", mode)
        sf = os.path.join(outdir, "summary.json")
        if not os.path.exists(sf):
            continue
        summary = json.load(open(sf))
        q = summary["quantification"]
        fstats = summary["fragment_stats"]
        ver = summary.get("rigel_version", "?")
        ts = summary.get("timestamp", "?")
        rows.append({
            "mode": mode,
            "version": ver,
            "timestamp": ts[:16] if len(ts) > 16 else ts,
            "mrna": q["mrna_total"],
            "nrna": q["nrna_total"],
            "gdna": q["gdna_total"],
            "inter": q["intergenic_total"],
            "n_loci": q["n_loci"],
            "detected": fstats["total"],
        })
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))

print(f"\n{'=' * 100}")
print("DIAGNOSTIC 3: Annotated nRNA equivalent contamination check")
print("=" * 100)

for cond in CONDITIONS[:1]:  # Just clean for brevity
    for mode in ["vbem"]:
        outdir = os.path.join(RUNS, cond, "rigel", mode)
        quant = pd.read_feather(os.path.join(outdir, "quant.feather"))

        # How many is_nrna=True transcripts?
        nrna_equiv = quant[quant["is_nrna"] == True]
        nrna_nonzero = nrna_equiv[nrna_equiv["count"] > 0]

        print(f"\n  [{cond[:20]}] {mode}:")
        print(f"    is_nrna=True in quant.feather: {len(nrna_equiv)}")
        print(f"    is_nrna=True with count > 0: {len(nrna_nonzero)}")
        print(f"    count sum (is_nrna=True): {nrna_equiv['count'].sum():,.1f}")
        print(f"    tpm sum (is_nrna=True): {nrna_equiv['tpm'].sum():,.1f}")
        print(f"    tpm sum (is_nrna=False): {quant.loc[~quant['is_nrna'], 'tpm'].sum():,.1f}")

        # Check truth for these transcripts
        merged = nrna_equiv.merge(
            truth_none[["transcript_id", "mrna_abundance"]],
            on="transcript_id",
            how="left",
        )
        n_with_truth = (merged["mrna_abundance"] > 0).sum()
        n_without_truth = (merged["mrna_abundance"] == 0).sum()
        n_na = merged["mrna_abundance"].isna().sum()
        print(f"    nRNA equiv with truth mrna > 0: {n_with_truth}")
        print(f"    nRNA equiv with truth mrna = 0: {n_without_truth}")
        print(f"    nRNA equiv not in truth: {n_na}")

        # What TPM budget goes to nRNA equiv with mrna_truth = 0?
        spurious = merged[(merged["mrna_abundance"] == 0) & (merged["count"] > 0)]
        print(f"    Spurious (truth=0, pred>0): n={len(spurious)}, "
              f"count={spurious['count'].sum():,.1f}, tpm={spurious['tpm'].sum():,.1f}")

print(f"\n{'=' * 100}")
print("DIAGNOSTIC 4: Error decomposition (alignment loss vs model error)")
print("=" * 100)

for cond in CONDITIONS:
    meta = conditions_meta[cond]
    n_rna_truth = meta["n_rna"]
    n_gdna_truth = meta["n_gdna"]

    for mode in ["vbem"]:
        outdir = os.path.join(RUNS, cond, "rigel", mode)
        summary = json.load(open(os.path.join(outdir, "summary.json")))
        q = summary["quantification"]
        fstats = summary["fragment_stats"]

        total_detected = fstats["total"]
        total_truth = n_rna_truth + n_gdna_truth
        alignment_loss = total_truth - total_detected

        # Pool from summary
        mrna_pred = q["mrna_total"]
        nrna_pred = q["nrna_total"]
        gdna_pred = q["gdna_total"] + q["intergenic_total"]

        # If we had oracle alignment, what would pool look like?
        oracle_dir = os.path.join(RUNS, cond, "rigel", "oracle")
        if os.path.exists(os.path.join(oracle_dir, "summary.json")):
            oracle_sum = json.load(open(os.path.join(oracle_dir, "summary.json")))
            oq = oracle_sum["quantification"]
            of = oracle_sum["fragment_stats"]
            oracle_mrna = oq["mrna_total"]
            oracle_nrna = oq["nrna_total"]
            oracle_gdna = oq["gdna_total"] + oq["intergenic_total"]
            oracle_detected = of["total"]

            print(f"\n  [{cond}] {mode}")
            print(f"    {'':30s}  {'TRUTH':>12s}  {'ORACLE':>12s}  {'STAR':>12s}  {'STAR ERR':>10s}  {'ALIGN LOSS':>10s}  {'MODEL ERR':>10s}")
            print(f"    {'mRNA frags':30s}  {n_rna_truth:>12,}  {oracle_mrna:>12,.0f}  {mrna_pred:>12,.0f}  "
                  f"{mrna_pred - n_rna_truth:>+10,.0f}  {oracle_mrna - n_rna_truth:>+10,.0f}  {mrna_pred - oracle_mrna:>+10,.0f}")
            print(f"    {'gDNA frags':30s}  {n_gdna_truth:>12,}  {oracle_gdna:>12,.0f}  {gdna_pred:>12,.0f}  "
                  f"{gdna_pred - n_gdna_truth:>+10,.0f}  {oracle_gdna - n_gdna_truth:>+10,.0f}  {gdna_pred - oracle_gdna:>+10,.0f}")
            print(f"    {'nRNA leak (truth=0)':30s}  {'0':>12s}  {oracle_nrna:>12,.0f}  {nrna_pred:>12,.0f}  "
                  f"{nrna_pred:>+10,.0f}  {oracle_nrna:>+10,.0f}  {nrna_pred - oracle_nrna:>+10,.0f}")
            print(f"    {'Total detected':30s}  {total_truth:>12,}  {oracle_detected:>12,}  {total_detected:>12,}")
            print(f"    {'Alignment loss':30s}  {'':>12s}  {total_truth - oracle_detected:>+12,}  {alignment_loss:>+12,}")

print(f"\n{'=' * 100}")
print("DIAGNOSTIC 5: TPM column correctness check")
print("=" * 100)

for cond in CONDITIONS[:1]:
    for mode in ["vbem"]:
        outdir = os.path.join(RUNS, cond, "rigel", mode)
        quant = pd.read_feather(os.path.join(outdir, "quant.feather"))

        # Check TPM calculation manually
        eff_len = quant["effective_length"].values
        count = quant["count"].values
        rpk = count / eff_len
        rpk_sum = rpk.sum()
        manual_tpm = rpk / rpk_sum * 1e6

        tpm_diff = np.abs(quant["tpm"].values - manual_tpm)
        print(f"\n  [{cond[:20]}] {mode}:")
        print(f"    TPM sum: {quant['tpm'].sum():.1f} (should be 1M)")
        print(f"    Manual TPM recomputation max error: {tpm_diff.max():.2e}")
        print(f"    tpm_total_rna sum: {quant['tpm_total_rna'].sum():.1f}")

        # What is the denominator difference?
        # tpm uses rpk_sum from annotated transcripts only
        # tpm_total_rna uses rpk_sum from ALL transcripts (annotated + synthetic nRNA)
        # But wait - quant.feather doesn't have synthetics (they're filtered out)
        # So the difference is that tpm_total_rna has count/eff_len from synthetics
        # included in its denominator
        total_rna_tpm = quant["tpm_total_rna"].values
        # If tpm_total_rna < tpm for all transcripts, synthetic nRNAs inflate denominator
        n_lower = (total_rna_tpm < quant["tpm"].values - 0.001).sum()
        print(f"    Transcripts where tpm_total_rna < tpm: {n_lower} (synthetics inflate denominator)")

        # What fraction of TPM budget goes to synthetics?
        tpm_ratio = quant["tpm_total_rna"].sum() / quant["tpm"].sum()
        print(f"    tpm_total_rna / tpm ratio: {tpm_ratio:.6f} "
              f"(1.0 means synthetics have zero RPK)")
        synthetic_tpm_budget = 1.0 - tpm_ratio
        print(f"    Synthetic nRNA TPM budget: {synthetic_tpm_budget*100:.3f}% of total RNA TPM")

print(f"\n{'=' * 100}")
print("DIAGNOSTIC 6: Analysis truth vs predicted comparison audit")
print("=" * 100)

for cond in CONDITIONS[:1]:
    for mode in ["vbem"]:
        outdir = os.path.join(RUNS, cond, "rigel", mode)
        quant = pd.read_feather(os.path.join(outdir, "quant.feather"))

        # Replicate exactly what analysis.py does
        merged = truth_none.merge(
            quant[["transcript_id", "tpm"]].rename(columns={"tpm": "predicted"}),
            on="transcript_id",
            how="left",
        )
        merged["predicted"] = merged["predicted"].fillna(0.0)

        truth_raw = merged["mrna_abundance"].values.astype(np.float64)
        truth_sum = truth_raw.sum()
        truth_arr = truth_raw * (1e6 / truth_sum) if truth_sum > 0 else truth_raw

        pred_arr = merged["predicted"].values.astype(np.float64)

        print(f"\n  [{cond[:20]}] {mode}:")
        print(f"    truth_arr sum: {truth_arr.sum():.1f}")
        print(f"    pred_arr sum: {pred_arr.sum():.1f}")
        print(f"    n_transcripts: {len(merged)}")
        print(f"    n_truth > 0: {(truth_arr > 0).sum()}")
        print(f"    n_pred > 0: {(pred_arr > 0).sum()}")

        # What is the total TPM going to transcripts with truth=0?
        zero_truth_mask = truth_arr == 0
        pred_at_zero = pred_arr[zero_truth_mask]
        print(f"    n_transcripts with truth=0: {zero_truth_mask.sum()}")
        print(f"    TPM assigned to truth=0 transcripts: {pred_at_zero.sum():.1f}")
        print(f"    Median pred for truth=0: {np.median(pred_at_zero):.4f}")

        # Conversely, truth > 0 but pred = 0
        zero_pred_mask = pred_arr == 0
        truth_at_zero_pred = truth_arr[zero_pred_mask & ~zero_truth_mask]
        print(f"    n_expressed with pred=0: {len(truth_at_zero_pred)}")
        print(f"    TPM missed (truth>0, pred=0): {truth_at_zero_pred.sum():.1f}")

print(f"\n{'=' * 100}")
print("DIAGNOSTIC COMPLETE")
print("=" * 100)
