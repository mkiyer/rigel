"""Probe Phase 2 Bayesian-prior benchmark issues against oracle BAM tags.

This script is intentionally diagnostic-only. It reads the synthetic benchmark
outputs under ~/Downloads and prints compact tables that separate prior mass,
EM posterior mass, and sampled annotated-BAM assignments.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic")

CONDITIONS_SS99 = [
    "gdna_none_ss_0.99_nrna_none",
    "gdna_low_ss_0.99_nrna_none",
    "gdna_med_ss_0.99_nrna_none",
    "gdna_equal_ss_0.99_nrna_none",
    "gdna_high_ss_0.99_nrna_none",
]

CONDITIONS_ALL = [
    "gdna_none_ss_0.99_nrna_none",
    "gdna_none_ss_0.50_nrna_none",
    "gdna_low_ss_0.99_nrna_none",
    "gdna_low_ss_0.50_nrna_none",
    "gdna_med_ss_0.99_nrna_none",
    "gdna_med_ss_0.50_nrna_none",
    "gdna_equal_ss_0.99_nrna_none",
    "gdna_equal_ss_0.50_nrna_none",
    "gdna_high_ss_0.99_nrna_none",
    "gdna_high_ss_0.50_nrna_none",
]


def _zf(read: pysam.AlignedSegment) -> int:
    try:
        return int(read.get_tag("ZF"))
    except KeyError:
        return 0


def _zl(read: pysam.AlignedSegment) -> int:
    try:
        return int(read.get_tag("ZL"))
    except KeyError:
        return -1


def _zc(read: pysam.AlignedSegment) -> str:
    try:
        return str(read.get_tag("ZC"))
    except KeyError:
        return ""


def summarize_loci(cond: str) -> dict[str, float]:
    loci = pd.read_feather(BASE / cond / "rigel_out" / "loci.feather")
    summary = json.loads((BASE / cond / "rigel_out" / "summary.json").read_text())
    manifest = json.loads((BASE / "manifest.json").read_text())
    condition_meta = {row["name"]: row for row in manifest["conditions"]}[cond]

    genic_fragments = float(summary["fragment_stats"]["genic"])
    n_rna = float(condition_meta["n_rna"])
    true_genic_gdna_proxy = max(genic_fragments - n_rna, 0.0)

    alpha = float(loci["alpha_gdna"].sum())
    gdna = float(loci["gdna"].sum())
    nrna = float(loci["nrna"].sum())
    mrna = float(loci["mrna"].sum())
    n_em = float(loci["n_em_fragments"].sum())
    return {
        "alpha": alpha,
        "gdna": gdna,
        "nrna": nrna,
        "mrna": mrna,
        "n_em": n_em,
        "true_genic_gdna_proxy": true_genic_gdna_proxy,
        "alpha_over_gdna": alpha / gdna if gdna else np.nan,
        "alpha_over_gdna_plus_nrna": alpha / (gdna + nrna) if (gdna + nrna) else np.nan,
        "alpha_over_true_genic": alpha / true_genic_gdna_proxy
        if true_genic_gdna_proxy
        else np.nan,
        "gdna_over_true_genic": gdna / true_genic_gdna_proxy
        if true_genic_gdna_proxy
        else np.nan,
        "mean_prior_rate": alpha / n_em if n_em else np.nan,
    }


def density_summary(cond: str) -> dict[str, float]:
    summary = json.loads((BASE / cond / "rigel_out" / "summary.json").read_text())
    densities = summary["calibration"]["global_densities"]
    rho_ig = float(densities["INTERGENIC"]["rho"])
    rho_in = float(densities["INTRON"]["rho"])
    rho_ex = float(densities["EXON-INTRON"]["rho"])
    return {
        "rho_ig": rho_ig,
        "rho_in": rho_in,
        "rho_ex": rho_ex,
        "rho_ex_over_ig": rho_ex / rho_ig if rho_ig else np.nan,
        "rho_in_over_ig": rho_in / rho_ig if rho_ig else np.nan,
    }


def scan_annotated_bam(cond: str) -> tuple[dict[str, int], pd.DataFrame]:
    bam_path = BASE / cond / "annotated.bam"
    counters: Counter[str] = Counter()
    by_locus: dict[int, Counter[str]] = defaultdict(Counter)
    gdna_name_tokens: Counter[str] = Counter()
    gdna_categories: Counter[str] = Counter()
    gdna_zf_categories: Counter[int] = Counter()
    zl_values: Counter[int] = Counter()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            is_true_gdna = read.query_name.startswith("gdna")
            zf = _zf(read)
            zc = _zc(read)
            zl = _zl(read)
            assigned_gdna = bool(zf & 0x04) or zc == "intergenic"
            assigned_tx = str(read.get_tag("ZT")) if read.has_tag("ZT") else ""
            assigned_nrna = assigned_tx.startswith("RIGEL_NRNA")

            zl_values[zl] += 1
            if is_true_gdna:
                counters["true_gdna"] += 1
                token = read.query_name.split(":", 2)[1] if ":" in read.query_name else ""
                gdna_name_tokens[token] += 1
                gdna_categories[zc] += 1
                gdna_zf_categories[zf] += 1
                if assigned_gdna:
                    counters["true_gdna_assigned_gdna"] += 1
                else:
                    counters["true_gdna_assigned_rna"] += 1
                    if assigned_nrna:
                        counters["true_gdna_assigned_nrna"] += 1
                    else:
                        counters["true_gdna_assigned_mrna"] += 1
                if zl >= 0:
                    counters["true_gdna_with_locus"] += 1
                    by_locus[zl]["true_gdna_with_locus"] += 1
                    if assigned_gdna:
                        by_locus[zl]["true_gdna_assigned_gdna"] += 1
                    else:
                        by_locus[zl]["true_gdna_assigned_rna"] += 1
            else:
                counters["true_rna"] += 1
                if assigned_gdna:
                    counters["true_rna_assigned_gdna"] += 1
                    if zl >= 0:
                        by_locus[zl]["true_rna_assigned_gdna"] += 1
                else:
                    counters["true_rna_assigned_rna"] += 1
                    if assigned_nrna:
                        counters["true_rna_assigned_nrna"] += 1

    print(f"\nOracle/BAM tag diagnostics for {cond}")
    print(f"  gDNA qname token top: {gdna_name_tokens.most_common(8)}")
    print(f"  true-gDNA ZC top:     {gdna_categories.most_common(8)}")
    print(f"  true-gDNA ZF top:     {gdna_zf_categories.most_common(8)}")
    print(f"  ZL top:               {zl_values.most_common(8)}")

    rows = []
    for locus_id, values in by_locus.items():
        row = {"locus_id": locus_id}
        row.update(values)
        rows.append(row)
    by_locus_df = pd.DataFrame(rows).fillna(0)
    if not by_locus_df.empty:
        int_cols = [c for c in by_locus_df.columns if c != "locus_id"]
        by_locus_df[int_cols] = by_locus_df[int_cols].astype(np.int64)
    return dict(counters), by_locus_df


def print_global_tables() -> None:
    print("\n=== Locus EM totals versus prior totals ===")
    print(
        f"{'condition':<34s} {'alpha':>10s} {'true_g':>10s} {'gdna_em':>10s} "
        f"{'nrna_em':>9s} {'a/em':>7s} {'a/true':>8s} {'em/true':>8s} {'rho_ex/ig':>9s}"
    )
    for cond in CONDITIONS_SS99:
        loc = summarize_loci(cond)
        den = density_summary(cond)
        print(
            f"{cond:<34s} {loc['alpha']:>10.0f} {loc['true_genic_gdna_proxy']:>10.0f} "
            f"{loc['gdna']:>10.0f} {loc['nrna']:>9.0f} {loc['alpha_over_gdna']:>7.3f} "
            f"{loc['alpha_over_true_genic']:>8.3f} {loc['gdna_over_true_genic']:>8.3f} "
            f"{den['rho_ex_over_ig']:>9.3f}"
        )


def print_bam_comparison(cond: str) -> None:
    loci = pd.read_feather(BASE / cond / "rigel_out" / "loci.feather")
    counters, by_locus = scan_annotated_bam(cond)

    print("\n=== Oracle assignment counters ===")
    for key in [
        "true_rna",
        "true_gdna",
        "true_rna_assigned_gdna",
        "true_gdna_assigned_gdna",
        "true_gdna_assigned_rna",
        "true_gdna_assigned_nrna",
        "true_gdna_assigned_mrna",
        "true_rna_assigned_nrna",
        "true_gdna_with_locus",
    ]:
        label = "true_gdna_with_zl_tag" if key == "true_gdna_with_locus" else key
        print(f"  {label:<28s} {int(counters.get(key, 0)):>10d}")

    if by_locus.empty:
        return

    merged = loci.merge(by_locus, on="locus_id", how="left").fillna(0)
    for col in [
        "true_gdna_with_locus",
        "true_gdna_assigned_gdna",
        "true_gdna_assigned_rna",
        "true_rna_assigned_gdna",
    ]:
        if col not in merged:
            merged[col] = 0
    merged["em_minus_alpha"] = merged["gdna"] - merged["alpha_gdna"]
    merged["oracle_locus_minus_alpha"] = merged["true_gdna_with_locus"] - merged["alpha_gdna"]
    merged["alpha_rate"] = merged["alpha_gdna"] / merged["n_em_fragments"].clip(lower=1)

    print("\n=== Sums merged by locus ===")
    print("  Note: true_gdna_with_locus uses ZL tags; gDNA-assigned reads generally carry ZL=-1.")
    sum_cols = [
        "alpha_gdna",
        "gdna",
        "nrna",
        "n_em_fragments",
        "true_gdna_with_locus",
        "true_gdna_assigned_gdna",
        "true_gdna_assigned_rna",
        "true_rna_assigned_gdna",
    ]
    for col in sum_cols:
        print(f"  {col:<28s} {float(merged[col].sum()):>12.1f}")

    print("\nTop loci where EM gDNA exceeds alpha:")
    cols = [
        "locus_id",
        "n_transcripts",
        "n_em_fragments",
        "alpha_gdna",
        "gdna",
        "nrna",
        "mrna",
        "true_gdna_with_locus",
        "true_rna_assigned_gdna",
        "em_minus_alpha",
        "alpha_rate",
    ]
    print(merged.nlargest(10, "em_minus_alpha")[cols].to_string(index=False))

    print("\nTop sparse prior loci by alpha/n_em:")
    print(merged.nlargest(10, "alpha_rate")[cols].to_string(index=False))


def print_nrna_scaling() -> None:
    print("\n=== nRNA false positives versus gDNA level ===")
    print(f"{'condition':<34s} {'nrna':>9s} {'gdna':>10s} {'nrna/gdna':>10s} {'nrna/total':>10s}")
    for cond in CONDITIONS_ALL:
        loc = summarize_loci(cond)
        total = loc["mrna"] + loc["nrna"] + loc["gdna"]
        print(
            f"{cond:<34s} {loc['nrna']:>9.0f} {loc['gdna']:>10.0f} "
            f"{(loc['nrna'] / loc['gdna'] if loc['gdna'] else np.nan):>10.4f} "
            f"{(loc['nrna'] / total if total else np.nan):>10.4f}"
        )


def main() -> None:
    print_global_tables()
    print_nrna_scaling()
    print_bam_comparison("gdna_high_ss_0.99_nrna_none")


if __name__ == "__main__":
    main()
