"""Estimate per-locus γ (gdna_prior) under the fixed get_loci_df pathway.

Re-implements compute_locus_priors using only artifacts available
post-run:

  * regions.feather (has mappable_effective_length and region coords)
  * loci.feather + quant.feather (locus → transcript membership)
  * transcripts.feather (genomic coords)
  * summary.json (calibration.lambda_gdna)

It cannot reproduce ``region_n_total`` exactly (which requires
region accumulator output from the BAM scan), so it substitutes the
locus-observed ``n_em_fragments`` as a proxy for n_total within the
locus's overlapping regions.  This is a **lower bound** on Σ n_total
(which also counts unique mappers, intergenic, etc), so the γ
estimates here are **upper bounds** on the true γ.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    run_dir = Path(
        sys.argv[1] if len(sys.argv) > 1
        else "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna80m"
    )
    index_dir = Path(
        sys.argv[2] if len(sys.argv) > 2
        else "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index"
    )

    loci = pd.read_feather(run_dir / "rigel" / "loci.feather")
    summary = json.load(open(run_dir / "rigel" / "summary.json"))
    lam_G = float(summary["calibration"]["lambda_gdna"])
    total_expected_gdna = float(summary["calibration"]["total_expected_gdna"])
    global_frac = float(summary["calibration"]["gdna_fraction"])
    total_n_global = total_expected_gdna / max(global_frac, 1e-12)
    print(f"[calib] lam_G = {lam_G:.4e}")
    print(f"[calib] global E[gdna] = {total_expected_gdna:,.0f}")
    print(f"[calib] global γ  = {global_frac:.4f}")
    print(f"[calib] implied total_n_global = {total_n_global:,.0f}")

    # Reuse the diagnostic TSV we built earlier (has mappable_bp per locus)
    diag_path = run_dir / "rigel" / "locus_mbp_diagnostic.tsv"
    if not diag_path.exists():
        print(f"[err] missing {diag_path}; run locus_mappable_vs_gdna_span.py first")
        sys.exit(1)
    diag = pd.read_csv(diag_path, sep="\t")
    df = loci.merge(diag[["locus_id", "mappable_bp", "span"]], on="locus_id", how="left")

    # Proxy γ = lam_G × mappable_bp / n_em_fragments (upper bound)
    df["E_gdna_locus"] = lam_G * df["mappable_bp"]
    df["gamma_proxy"] = df["E_gdna_locus"] / df["n_em_fragments"].clip(lower=1)
    df["gamma_proxy"] = df["gamma_proxy"].clip(upper=1.0)

    print()
    print("[proxy γ] (upper bound — uses n_em_fragments as n_total proxy)")
    print(df["gamma_proxy"].describe())

    print()
    print("[size-bin]")
    df["size_bin"] = pd.cut(
        df["locus_span_bp"],
        bins=[0, 1e4, 1e5, 1e6, 1e7, 1e10],
        labels=["<10kb", "10-100kb", "100kb-1Mb", "1-10Mb", ">10Mb"],
    )
    g = df.groupby("size_bin", observed=True).agg(
        n_loci=("locus_id", "size"),
        frags=("n_em_fragments", "sum"),
        E_gdna=("E_gdna_locus", "sum"),
        gamma_proxy_mean=("gamma_proxy", "mean"),
        gamma_proxy_median=("gamma_proxy", "median"),
        gdna_rate_wmean=(
            "gdna_rate",
            lambda s: np.average(s, weights=df.loc[s.index, "n_em_fragments"]),
        ),
    )
    g["gamma_bulk"] = g["E_gdna"] / g["frags"]
    print(g.to_string())

    print()
    print("[mega locus]")
    mega = df[df["locus_span_bp"] > 1e7].iloc[0]
    print(f"  locus_id          = {int(mega['locus_id'])}")
    print(f"  locus_span_bp     = {mega['locus_span_bp']:,.0f}")
    print(f"  mappable_bp       = {mega['mappable_bp']:,.0f}")
    print(f"  n_em_fragments    = {mega['n_em_fragments']:,.0f}")
    print(f"  E[gdna_locus]     = λ_G × mappable_bp = {mega['E_gdna_locus']:,.0f}")
    print(f"  γ_proxy (upper)   = {mega['gamma_proxy']:.4f}")
    print(f"  gdna_rate observed= {mega['gdna_rate']:.4f}")
    print(f"  gdna_prior in file= {mega['gdna_prior']:.4f}  <-- DISPLAY BUG (see doc)")


if __name__ == "__main__":
    main()
