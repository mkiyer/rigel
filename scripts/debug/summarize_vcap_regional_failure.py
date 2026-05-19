#!/usr/bin/env python3
"""Summarize VCaP regional-exposure failure diagnostics from generated TSVs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m")
OUT = BASE / "regional_v3_confusion"


def main() -> int:
    locus_false = pd.read_csv(OUT / "gdna_before_correct_after_rna_by_locus.tsv", sep="\t")
    off_loci = pd.read_feather(BASE / "regional_off_v3_with_mm" / "loci.feather")
    auto_loci = pd.read_feather(BASE / "regional_auto_v3_with_mm" / "loci.feather")

    bins = [0.0, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2.5e-1, 5e-1, 1.0, np.inf]
    labels = ["<=1e-4", "1e-4..1e-3", "1e-3..1e-2", "1e-2..5e-2", "5e-2..0.1",
              "0.1..0.25", "0.25..0.5", "0.5..1", ">1"]
    locus_false["ratio_bin"] = pd.cut(
        locus_false["gdna_eff_len_weight_ratio"], bins=bins, labels=labels, include_lowest=True
    )
    by_bin = locus_false.groupby("ratio_bin", observed=True).agg(
        n_loci=("locus_id", "count"),
        new_false=("new_false_gdna_to_rna", "sum"),
        to_mrna=("new_false_to_mrna", "sum"),
        to_nrna=("new_false_to_nrna", "sum"),
        median_span=("locus_span_bp", "median"),
        median_n_tx=("n_transcripts", "median"),
    ).reset_index()
    by_bin["fraction_of_new_false"] = by_bin["new_false"] / by_bin["new_false"].sum()
    by_bin.to_csv(OUT / "gdna_before_correct_after_rna_by_ratio_bin.tsv", sep="\t", index=False)

    merged = auto_loci[[
        "locus_id", "mrna", "nrna", "gdna", "gdna_eff_len", "gdna_eff_len_unweighted",
        "gdna_eff_len_weight_ratio", "gdna_prior_count", "gdna_prior_count_regional",
        "n_em_fragments", "n_transcripts", "locus_span_bp",
    ]].merge(
        off_loci[["locus_id", "mrna", "nrna", "gdna", "gdna_eff_len"]],
        on="locus_id",
        suffixes=("_auto", "_off"),
    )
    merged["delta_gdna"] = merged["gdna_auto"] - merged["gdna_off"]
    merged["delta_nrna"] = merged["nrna_auto"] - merged["nrna_off"]
    merged["delta_mrna"] = merged["mrna_auto"] - merged["mrna_off"]
    merged["gdna_eff_len_log_change"] = np.log(
        np.maximum(merged["gdna_eff_len_auto"], 1.0)
    ) - np.log(np.maximum(merged["gdna_eff_len_off"], 1.0))
    top_lids = locus_false.head(50)["locus_id"].astype(int)
    top_delta = merged[merged["locus_id"].isin(top_lids)].copy()
    top_delta = top_delta.merge(
        locus_false[["locus_id", "new_false_gdna_to_rna", "new_false_to_mrna", "new_false_to_nrna"]],
        on="locus_id",
        how="left",
    ).sort_values("new_false_gdna_to_rna", ascending=False)
    top_delta.to_csv(OUT / "top_false_loci_before_after_delta.tsv", sep="\t", index=False)

    print("new false gDNA->RNA by after-locus gDNA effective-length ratio")
    print(by_bin.to_string(index=False))
    print("\ntop before/after locus deltas")
    print(top_delta.head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())