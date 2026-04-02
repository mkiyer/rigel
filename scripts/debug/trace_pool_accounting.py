#!/usr/bin/env python3
"""Trace pool-level fragment accounting to find the source of pool errors.

For each condition/mode, account for every fragment:
- Where do fragments go? (mRNA count, nRNA count, gDNA count, intergenic, unassigned)
- What does the truth say?
- Where is the discrepancy?
"""
import pandas as pd
import numpy as np
import json
import os

BDIR = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks"
RUNS = os.path.join(BDIR, "runs")

# Load manifest for truth
manifest = json.load(open(os.path.join(BDIR, "manifest.json")))
sim = manifest["simulation"]
conditions_meta = {c["name"]: c for c in manifest["conditions"]}

# Load truth abundances
truth_none = pd.read_csv(
    os.path.join(BDIR, "truth_abundances_nrna_none.tsv"), sep="\t"
)

CONDITIONS = ["gdna_none_ss_1.00_nrna_none", "gdna_high_ss_0.90_nrna_none"]

for cond in CONDITIONS:
    meta = conditions_meta[cond]
    n_rna_truth = meta["n_rna"]
    n_gdna_truth = meta["n_gdna"]
    total_frags_truth = n_rna_truth + n_gdna_truth

    print(f"\n{'='*80}")
    print(f"CONDITION: {cond}")
    print(f"{'='*80}")
    print(f"  Truth: n_rna={n_rna_truth:,}, n_gdna={n_gdna_truth:,}, total={total_frags_truth:,}")

    # Truth TPM sums
    truth_mrna_sum = truth_none["mrna_abundance"].sum()
    truth_nrna_sum = truth_none["nrna_abundance"].sum()
    print(f"  Truth TPM: mrna_sum={truth_mrna_sum:.1f}, nrna_sum={truth_nrna_sum:.1f}")

    for mode in ["vbem", "map"]:
        print(f"\n  --- {mode.upper()} ---")

        outdir = os.path.join(RUNS, cond, "rigel", mode)
        summary = json.load(open(os.path.join(outdir, "summary.json")))
        quant_df = pd.read_feather(os.path.join(outdir, "quant.feather"))
        nrna_df = pd.read_feather(os.path.join(outdir, "nrna_quant.feather"))
        loci_df = pd.read_feather(os.path.join(outdir, "loci.feather"))

        # Print quant.feather structure
        print(f"    quant.feather: {quant_df.shape}, columns: {list(quant_df.columns)}")
        print(f"    nrna_quant.feather: {nrna_df.shape}, columns: {list(nrna_df.columns)}")

        # Pool from summary.json
        q = summary["quantification"]
        print(f"\n    summary.json pool totals:")
        print(f"      mrna_total:       {q['mrna_total']:>12,.1f}")
        print(f"      nrna_total:       {q['nrna_total']:>12,.1f}")
        print(f"      gdna_total:       {q['gdna_total']:>12,.1f}")
        print(f"      intergenic_total: {q['intergenic_total']:>12,}")
        total_accounted = q["mrna_total"] + q["nrna_total"] + q["gdna_total"] + q["intergenic_total"]
        print(f"      TOTAL accounted:  {total_accounted:>12,.1f}")
        print(f"      TRUTH total:      {total_frags_truth:>12,}")
        print(f"      MISSING:          {total_frags_truth - total_accounted:>12,.1f}")

        # Alignment stats
        astats = summary["alignment_stats"]
        fstats = summary["fragment_stats"]
        print(f"\n    Alignment stats:")
        print(f"      total_reads:       {astats['total_reads']:>12,}")
        print(f"      mapped_reads:      {astats['mapped_reads']:>12,}")
        print(f"      unique_reads:      {astats['unique_reads']:>12,}")
        print(f"      multimapping:      {astats['multimapping_reads']:>12,}")
        print(f"      proper_pairs:      {astats['proper_pairs']:>12,}")
        print(f"\n    Fragment stats:")
        print(f"      total:             {fstats['total']:>12,}")
        print(f"      genic:             {fstats['genic']:>12,}")
        print(f"      intergenic:        {fstats['intergenic']:>12,}")
        print(f"      chimeric:          {fstats['chimeric']:>12,}")

        # N fragments entering EM
        n_unambig = q["n_unambig_assigned"]
        n_em = q["n_em_assigned"]
        print(f"\n    Fragment routing:")
        print(f"      n_unambig_assigned: {n_unambig:>12,}")
        print(f"      n_em_assigned:      {n_em:>12,}")
        print(f"      total assigned:     {n_unambig + n_em:>12,}")

        # quant.feather count sums
        total_count = quant_df["count"].sum()
        is_nrna_mask = quant_df["is_nrna"] == True
        mrna_count = quant_df.loc[~is_nrna_mask, "count"].sum()
        nrna_equiv_count = quant_df.loc[is_nrna_mask, "count"].sum()
        print(f"\n    quant.feather counts:")
        print(f"      total count:                 {total_count:>12,.1f}")
        print(f"      count (is_nrna=False):       {mrna_count:>12,.1f}")
        print(f"      count (is_nrna=True):        {nrna_equiv_count:>12,.1f}")

        # nrna_quant.feather
        if "count" in nrna_df.columns:
            nrna_count_total = nrna_df["count"].sum()
            print(f"      nrna_quant count total:      {nrna_count_total:>12,.1f}")

        # TPM sums
        tpm_sum = quant_df["tpm"].sum() if "tpm" in quant_df.columns else 0
        tpm_total_rna = (
            quant_df["tpm_total_rna"].sum() if "tpm_total_rna" in quant_df.columns else 0
        )
        print(f"\n    TPM sums:")
        print(f"      tpm sum:           {tpm_sum:>12,.1f}")
        print(f"      tpm_total_rna sum: {tpm_total_rna:>12,.1f}")

        # Loci pool breakdown
        loci_mrna = loci_df["mrna"].sum()
        loci_nrna = loci_df["nrna"].sum()
        loci_gdna = loci_df["gdna"].sum()
        loci_total = loci_df["total"].sum()
        print(f"\n    loci.feather pool sums:")
        print(f"      mrna:  {loci_mrna:>12,.1f}")
        print(f"      nrna:  {loci_nrna:>12,.1f}")
        print(f"      gdna:  {loci_gdna:>12,.1f}")
        print(f"      total: {loci_total:>12,.1f}")

        # Check: does quant count + nrna count + gdna count = total?
        print(f"\n    ACCOUNTING CHECK:")
        print(f"      quant count (all):           {total_count:>12,.1f}")
        print(f"      loci mrna:                   {loci_mrna:>12,.1f}")
        print(f"      loci nrna:                   {loci_nrna:>12,.1f}")
        print(f"      loci gdna:                   {loci_gdna:>12,.1f}")
        print(f"      intergenic:                  {fstats['intergenic']:>12,}")
        print(f"      chimeric:                    {fstats['chimeric']:>12,}")
        reconstruction = loci_mrna + loci_nrna + loci_gdna + fstats["intergenic"]
        print(f"      loci pools + intergenic:     {reconstruction:>12,.1f}")
        print(f"      truth total frags:           {total_frags_truth:>12,}")
        print(f"      gap:                         {total_frags_truth - reconstruction:>12,.1f}")
