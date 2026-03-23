#!/usr/bin/env python3
"""Investigate minimap2 v4 regression: which transcript pairs are swapping counts?

This script looks at the top minimap2 regressions and checks whether they are
part of gene families where minimap2 multimapping causes EM confusion.
"""

import pandas as pd
import numpy as np
from pathlib import Path

V3 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3/gdna_none_ss_0.95_nrna_none")
V4 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none")


def main():
    # Load all data
    v3_mm2 = pd.read_csv(V3 / "per_transcript_counts_minimap2.csv")
    v4_mm2 = pd.read_csv(V4 / "per_transcript_counts_minimap2.csv")
    v4_ora = pd.read_csv(V4 / "per_transcript_counts_oracle.csv")

    # Merge
    m = v3_mm2[["transcript_id", "gene_id", "gene_name", "mrna_truth", "rigel_minimap2"]].merge(
        v4_mm2[["transcript_id", "rigel_minimap2"]],
        on="transcript_id", suffixes=("_v3", "_v4"),
    ).merge(
        v4_ora[["transcript_id", "rigel_oracle"]],
        on="transcript_id",
    )

    m["err_v3"] = np.abs(m["rigel_minimap2_v3"] - m["mrna_truth"])
    m["err_v4"] = np.abs(m["rigel_minimap2_v4"] - m["mrna_truth"])
    m["delta"] = m["err_v4"] - m["err_v3"]
    m["oracle_err"] = np.abs(m["rigel_oracle"] - m["mrna_truth"])

    # Focus on regressions (delta > 100)
    regressed = m[m["delta"] > 100].sort_values("delta", ascending=False).head(50)

    print("=" * 120)
    print("TOP 50 MINIMAP2 REGRESSIONS (v4 worse than v3, delta > 100)")
    print("=" * 120)
    print(f"{'transcript_id':>30s}  {'gene':>15s}  {'truth':>8s}  {'ora':>8s}  {'mm2_v3':>8s}  {'mm2_v4':>8s}  {'v3err':>8s}  {'v4err':>8s}  {'delta':>8s}")
    for _, r in regressed.iterrows():
        print(f"  {r.transcript_id:>28s}  {str(r.gene_name):>15s}  {r.mrna_truth:>8.0f}  {r.rigel_oracle:>8.0f}  {r.rigel_minimap2_v3:>8.0f}  {r.rigel_minimap2_v4:>8.0f}  {r.err_v3:>8.0f}  {r.err_v4:>8.0f}  {r.delta:>+8.0f}")

    # Now look at the gene-level picture for the top regressed genes
    print("\n\n" + "=" * 120)
    print("GENE-LEVEL ANALYSIS: For each regressed gene, show ALL its transcripts")
    print("=" * 120)

    # Get genes involved in top regressions
    regressed_genes = regressed["gene_id"].unique()[:15]

    for gene_id in regressed_genes:
        gene_txs = m[m["gene_id"] == gene_id].sort_values("mrna_truth", ascending=False)
        gene_name = gene_txs.iloc[0]["gene_name"]
        total_truth = gene_txs["mrna_truth"].sum()
        total_ora = gene_txs["rigel_oracle"].sum()
        total_v3 = gene_txs["rigel_minimap2_v3"].sum()
        total_v4 = gene_txs["rigel_minimap2_v4"].sum()

        print(f"\n--- Gene: {gene_name} ({gene_id}) ---")
        print(f"  Total: truth={total_truth:.0f}  oracle={total_ora:.0f}  mm2_v3={total_v3:.0f}  mm2_v4={total_v4:.0f}")
        print(f"  {'transcript_id':>30s}  {'truth':>8s}  {'oracle':>8s}  {'mm2_v3':>8s}  {'mm2_v4':>8s}  {'delta':>8s}")
        for _, r in gene_txs.iterrows():
            print(f"  {r.transcript_id:>30s}  {r.mrna_truth:>8.0f}  {r.rigel_oracle:>8.0f}  {r.rigel_minimap2_v3:>8.0f}  {r.rigel_minimap2_v4:>8.0f}  {r.delta:>+8.0f}")

    # Summary: are the gene-level totals stable?
    print("\n\n" + "=" * 120)
    print("GENE-LEVEL TOTALS FOR REGRESSED GENES")
    print("=" * 120)
    gene_totals = m.groupby("gene_id").agg(
        gene_name=("gene_name", "first"),
        n_tx=("transcript_id", "count"),
        truth=("mrna_truth", "sum"),
        oracle=("rigel_oracle", "sum"),
        mm2_v3=("rigel_minimap2_v3", "sum"),
        mm2_v4=("rigel_minimap2_v4", "sum"),
    ).reset_index()
    gene_totals["gene_err_v3"] = np.abs(gene_totals["mm2_v3"] - gene_totals["truth"])
    gene_totals["gene_err_v4"] = np.abs(gene_totals["mm2_v4"] - gene_totals["truth"])
    gene_totals["gene_delta"] = gene_totals["gene_err_v4"] - gene_totals["gene_err_v3"]

    # Overall gene-level MAE
    print(f"\n  Gene-level MAE v3: {gene_totals['gene_err_v3'].mean():.2f}")
    print(f"  Gene-level MAE v4: {gene_totals['gene_err_v4'].mean():.2f}")
    print(f"  Gene-level delta:  {gene_totals['gene_delta'].mean():+.2f}")

    # Top gene-level regressions
    gene_reg = gene_totals.sort_values("gene_delta", ascending=False).head(30)
    print(f"\n  Top 30 gene-level regressions:")
    print(f"  {'gene_id':>25s}  {'gene_name':>15s}  {'n_tx':>5s}  {'truth':>8s}  {'oracle':>8s}  {'mm2_v3':>8s}  {'mm2_v4':>8s}  {'g_delta':>8s}")
    for _, r in gene_reg.iterrows():
        print(f"  {r.gene_id:>25s}  {str(r.gene_name):>15s}  {r.n_tx:>5d}  {r.truth:>8.0f}  {r.oracle:>8.0f}  {r.mm2_v3:>8.0f}  {r.mm2_v4:>8.0f}  {r.gene_delta:>+8.0f}")

    # Check: is the regression mainly intra-gene (isoform redistribution) or inter-gene (leakage)?
    print("\n\n" + "=" * 120)
    print("INTRA-GENE vs INTER-GENE ANALYSIS")
    print("=" * 120)

    # Transcript-level MAE
    tx_mae_v3 = np.mean(m["err_v3"])
    tx_mae_v4 = np.mean(m["err_v4"])

    # Gene-level MAE
    g_mae_v3 = gene_totals["gene_err_v3"].mean()
    g_mae_v4 = gene_totals["gene_err_v4"].mean()

    print(f"  Transcript-level MAE: v3={tx_mae_v3:.2f}  v4={tx_mae_v4:.2f}  delta={tx_mae_v4-tx_mae_v3:+.2f}")
    print(f"  Gene-level MAE:       v3={g_mae_v3:.2f}  v4={g_mae_v4:.2f}  delta={g_mae_v4-g_mae_v3:+.2f}")
    print(f"  Intra-gene component: v3={tx_mae_v3-g_mae_v3:.2f}  v4={tx_mae_v4-g_mae_v4:.2f}")

    # Check if top regressed transcripts are in multi-isoform genes
    print(f"\n  Among top 50 regressions:")
    for _, r in regressed.head(20).iterrows():
        n_isoforms = len(m[m["gene_id"] == r.gene_id])
        print(f"    {r.transcript_id:>30s} ({str(r.gene_name):>12s}): {n_isoforms} isoforms in gene, delta={r.delta:+.0f}")


if __name__ == "__main__":
    main()
