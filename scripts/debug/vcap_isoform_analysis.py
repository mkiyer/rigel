#!/usr/bin/env python3
"""Investigate zero-forced transcripts: are they single-isoform genes or isoform competition?

Key question: When VBEM zeros a transcript with truth=2000 TPM, does the mass
go to another isoform of the same gene, or to a completely different gene?
"""

import pandas as pd
import numpy as np
from pathlib import Path

BENCHMARK_DIR = Path("/scratch/mkiyer_root/mkiyer0/shared_data/rigel_benchmarks/ccle_vcap_prostate")
CONDITION = "gdna_high_ss_0.90_nrna_none"

def main():
    vbem_q = pd.read_feather(BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "vbem" / "quant.feather")
    map_q = pd.read_feather(BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "map" / "quant.feather")

    truth = pd.read_csv(BENCHMARK_DIR / "truth_abundances_nrna_none.tsv", sep="\t")
    total_mrna = truth["mrna_abundance"].sum()
    truth["truth_tpm"] = truth["mrna_abundance"] / total_mrna * 1e6

    # Merge for VBEM
    merged = vbem_q.merge(truth[["transcript_id", "truth_tpm", "gene_id", "n_exons"]], 
                          on="transcript_id", suffixes=("", "_truth"), how="left")
    merged["truth_tpm"] = merged["truth_tpm"].fillna(0)

    # Use gene_id from quant (already present)
    # Find zero-forced highly expressed transcripts
    zero_forced = merged[(merged["truth_tpm"] > 10) & (merged["tpm"] < 0.1)].copy()
    
    print(f"Zero-forced transcripts (truth > 10 TPM, vbem < 0.1): {len(zero_forced)}")
    print(f"Total lost TPM: {zero_forced['truth_tpm'].sum():.0f}")
    
    # For each zero-forced transcript, check same-gene behavior
    print(f"\n{'='*80}")
    print("ZERO-FORCED TRANSCRIPTS: SAME-GENE ANALYSIS")
    print(f"{'='*80}")
    
    results = []
    for _, zf_row in zero_forced.iterrows():
        tx_id = zf_row["transcript_id"]
        gene_id = zf_row["gene_id"]
        truth_tpm = zf_row["truth_tpm"]
        
        # All transcripts of same gene in VBEM
        gene_txs = merged[merged["gene_id"] == gene_id]
        gene_map = map_q[map_q["gene_id"] == gene_id]
        
        # How many isoforms?
        n_isoforms = len(gene_txs)
        gene_truth_tpm = gene_txs["truth_tpm"].sum()
        gene_vbem_tpm = gene_txs["tpm"].sum()
        gene_map_tpm = gene_map["tpm"].sum() if len(gene_map) > 0 else 0
        
        # Is the gene total preserved?
        results.append({
            "transcript_id": tx_id,
            "gene_id": gene_id,
            "gene_name": zf_row["gene_name"],
            "truth_tpm": truth_tpm,
            "vbem_tpm": zf_row["tpm"],
            "n_isoforms": n_isoforms,
            "gene_truth_tpm": gene_truth_tpm,
            "gene_vbem_tpm": gene_vbem_tpm,
            "gene_map_tpm": gene_map_tpm,
            "gene_preserved": abs(gene_vbem_tpm - gene_truth_tpm) / max(gene_truth_tpm, 1) < 0.3,
        })
    
    results_df = pd.DataFrame(results)
    
    # Summary
    n_preserved = results_df["gene_preserved"].sum()
    n_not_preserved = (~results_df["gene_preserved"]).sum()
    print(f"\nGene-level assessment:")
    print(f"  Gene total preserved (within 30%): {n_preserved} / {len(results_df)}")
    print(f"  Gene total NOT preserved: {n_not_preserved} / {len(results_df)}")
    
    # Single-isoform genes that are zero-forced (most damning)
    single_iso = results_df[results_df["n_isoforms"] == 1]
    print(f"\nSingle-isoform genes zero-forced: {len(single_iso)}")
    if len(single_iso) > 0:
        print(f"  Total lost TPM: {single_iso['truth_tpm'].sum():.0f}")
        for _, row in single_iso.nlargest(20, "truth_tpm").iterrows():
            print(f"    {row['transcript_id']} ({row['gene_name']}): "
                  f"truth={row['truth_tpm']:.1f}, vbem={row['vbem_tpm']:.2f}, "
                  f"gene_vbem={row['gene_vbem_tpm']:.2f}")
    
    # Multi-isoform genes where total is preserved (isoform competition)
    multi_preserved = results_df[(results_df["n_isoforms"] > 1) & results_df["gene_preserved"]]
    print(f"\nMulti-isoform, gene total preserved: {len(multi_preserved)}")
    if len(multi_preserved) > 0:
        print(f"  → These are isoform-level redistribution, not fundamental failures")
    
    # Multi-isoform genes where total is NOT preserved (mass leaking elsewhere)
    multi_lost = results_df[(results_df["n_isoforms"] > 1) & ~results_df["gene_preserved"]]
    print(f"\nMulti-isoform, gene total NOT preserved: {len(multi_lost)}")
    if len(multi_lost) > 0:
        print(f"  Total lost TPM: {multi_lost['truth_tpm'].sum():.0f}")
        for _, row in multi_lost.nlargest(20, "truth_tpm").iterrows():
            print(f"    {row['transcript_id']} ({row['gene_name']}): "
                  f"truth={row['truth_tpm']:.1f}, "
                  f"gene_truth={row['gene_truth_tpm']:.1f}, "
                  f"gene_vbem={row['gene_vbem_tpm']:.1f}, "
                  f"gene_map={row['gene_map_tpm']:.1f}, "
                  f"n_iso={row['n_isoforms']}")

    # Check nRNA siphon per gene
    print(f"\n{'='*80}")
    print("nRNA SIPHON ANALYSIS FOR ZERO-FORCED GENES")
    print(f"{'='*80}")
    
    vbem_nrna = pd.read_feather(BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "vbem" / "nrna_quant.feather")
    map_nrna = pd.read_feather(BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "map" / "nrna_quant.feather")
    
    print(f"VBEM nRNA columns: {list(vbem_nrna.columns)}")
    print(f"VBEM nRNA rows: {len(vbem_nrna)}")
    print(f"Total VBEM nRNA count: {vbem_nrna['count'].sum():.0f}")
    print(f"Total MAP nRNA count: {map_nrna['count'].sum():.0f}")
    
    # Where is the TPM mass going overall? 
    # Gene-level comparison
    print(f"\n{'='*80}")
    print("GENE-LEVEL COMPARISON (top genes by VBEM vs MAP difference)")
    print(f"{'='*80}")
    
    vbem_gene = pd.read_feather(BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "vbem" / "gene_quant.feather")
    map_gene = pd.read_feather(BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "map" / "gene_quant.feather")
    
    print(f"VBEM gene quant columns: {list(vbem_gene.columns)}")
    
    gene_merged = vbem_gene[["gene_id", "gene_name", "tpm"]].merge(
        map_gene[["gene_id", "tpm"]],
        on="gene_id", suffixes=("_vbem", "_map")
    )
    
    # Add truth
    gene_truth = truth.groupby("gene_id")["truth_tpm"].sum().reset_index()
    gene_merged = gene_merged.merge(gene_truth, on="gene_id", how="left")
    gene_merged["truth_tpm"] = gene_merged["truth_tpm"].fillna(0)
    
    gene_merged["diff"] = gene_merged["tpm_vbem"] - gene_merged["tpm_map"]
    
    # Biggest over-estimates by VBEM vs MAP
    print("\nTop 20 genes where VBEM >> MAP (mass sinks):")
    for _, row in gene_merged.nlargest(20, "diff").iterrows():
        print(f"  {row['gene_name']} ({row['gene_id']}): "
              f"truth={row['truth_tpm']:.1f}, "
              f"vbem={row['tpm_vbem']:.1f}, "
              f"map={row['tpm_map']:.1f}, "
              f"diff={row['diff']:.1f}")
    
    print("\nTop 20 genes where VBEM << MAP (mass lost):")
    for _, row in gene_merged.nsmallest(20, "diff").iterrows():
        print(f"  {row['gene_name']} ({row['gene_id']}): "
              f"truth={row['truth_tpm']:.1f}, "
              f"vbem={row['tpm_vbem']:.1f}, "
              f"map={row['tpm_map']:.1f}, "
              f"diff={row['diff']:.1f}")


if __name__ == "__main__":
    main()
