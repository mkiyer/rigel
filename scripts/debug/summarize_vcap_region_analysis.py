#!/usr/bin/env python3
"""Summarize VCaP region-analysis TSVs for quick inspection."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

OUT_DIR = Path("results/vcap_gdna_false_rna_by_region_kappa_units_fix_2026-05-19")


def pct(numerator: float, denominator: float) -> str:
    if denominator == 0:
        return "n/a"
    return f"{100.0 * numerator / denominator:.2f}%"


def main() -> None:
    mask = pd.read_csv(OUT_DIR / "regional_confusion_by_mask.tsv", sep="\t")
    bins = pd.read_csv(OUT_DIR / "dominant_region_bins.tsv", sep="\t")
    dominant = pd.read_csv(OUT_DIR / "dominant_region_type_confusion.tsv", sep="\t")
    top = pd.read_csv(OUT_DIR / "dominant_region_gdna_false_rna.tsv", sep="\t")

    gdna_rna_mask = mask[(mask.true_source == "gdna") & (mask.predicted_binary == "rna")]
    total_false = int(gdna_rna_mask["count"].sum())

    print("gDNA->RNA by exact mask and detail")
    mask_detail = (
        gdna_rna_mask.groupby(["region_mask", "predicted_detail"], as_index=False)["count"].sum()
        .sort_values("count", ascending=False)
    )
    for row in mask_detail.itertuples(index=False):
        print(f"{row.region_mask}\t{row.predicted_detail}\t{row.count}\t{pct(row.count, total_false)}")

    print("\ngDNA->RNA by dominant type")
    gdna_rna_dom = dominant[(dominant.true_source == "gdna") & (dominant.predicted_binary == "rna")]
    for row in gdna_rna_dom.groupby("dominant_region_type", as_index=False)["count"].sum().itertuples(index=False):
        print(f"{row.dominant_region_type}\t{row.count}\t{pct(row.count, total_false)}")

    print("\nDominant EXON gDNA->RNA by length bin")
    exon_bins = bins[bins.region_type == "EXON"].groupby("length_bin", as_index=False).agg(
        gdna_source_total=("gdna_source_total", "sum"),
        gdna_false_rna=("gdna_false_rna", "sum"),
    )
    order = {"1-100": 0, "101-250": 1, "251-500": 2, "501-1k": 3, "1k-5k": 4, "5k-10k": 5, "10k-50k": 6, ">50k": 7}
    exon_bins = exon_bins.sort_values("length_bin", key=lambda s: s.map(order))
    exon_total_false = int(exon_bins.gdna_false_rna.sum())
    for row in exon_bins.itertuples(index=False):
        print(
            f"{row.length_bin}\t{int(row.gdna_false_rna)}\t{pct(row.gdna_false_rna, exon_total_false)}"
            f"\trate={pct(row.gdna_false_rna, row.gdna_source_total)}"
        )

    print("\nDominant EXON gDNA->RNA by strand")
    exon_strand = bins[bins.region_type == "EXON"].groupby("region_strand", as_index=False).agg(
        gdna_source_total=("gdna_source_total", "sum"),
        gdna_false_rna=("gdna_false_rna", "sum"),
    )
    for row in exon_strand.sort_values("gdna_false_rna", ascending=False).itertuples(index=False):
        print(
            f"{row.region_strand}\t{int(row.gdna_false_rna)}\t{pct(row.gdna_false_rna, exon_total_false)}"
            f"\trate={pct(row.gdna_false_rna, row.gdna_source_total)}"
        )

    print("\nTop-region concentration")
    for n in [20, 100, 200]:
        count = int(top.head(n).gdna_false_rna.sum())
        print(f"top{n}\t{count}\t{pct(count, total_false)}")

    print("\nTop 200 dominant EXON complexity")
    top_exon = top[top.region_type == "EXON"].copy()
    for label, query in [
        ("len>=5kb", top_exon.length_bp >= 5000),
        ("len>=10kb", top_exon.length_bp >= 10000),
        ("ambig_strand", top_exon.region_strand == "AMBIG"),
        ("exon_genes>=2", top_exon.exon_gene_count >= 2),
        ("exon_tx>=10", top_exon.exon_transcript_count >= 10),
    ]:
        count = int(top_exon.loc[query, "gdna_false_rna"].sum())
        print(f"{label}\t{count}\t{pct(count, int(top_exon.gdna_false_rna.sum()))}")


if __name__ == "__main__":
    main()
