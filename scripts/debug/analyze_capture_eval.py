#!/usr/bin/env python
"""Analyze Rigel outputs for the synthetic hybrid-capture mini-genome runs."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb")
OUTBASE = BASE / "rigel_eval_20260525"
TRUTH = BASE / "truth_abundances_nrna_none.tsv"
GTF = BASE / "reference" / "genes.gtf"
CONDITIONS = [
    "gdna_none_ss_0.99_nrna_none",
    "gdna_none_ss_0.50_nrna_none",
    "gdna_high_ss_0.99_nrna_none",
    "gdna_high_ss_0.50_nrna_none",
]


def parse_gene_spans() -> dict[str, tuple[str, int, int]]:
    spans: dict[str, tuple[str, int, int]] = {}
    with GTF.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = {}
            for item in fields[8].split(";"):
                item = item.strip()
                if not item:
                    continue
                key, value = item.split(" ", 1)
                attrs[key] = value.strip('"')
            gene_id = attrs.get("gene_id")
            if gene_id is None:
                continue
            ref = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            if gene_id not in spans:
                spans[gene_id] = (ref, start, end)
            else:
                old_ref, old_start, old_end = spans[gene_id]
                spans[gene_id] = (old_ref, min(old_start, start), max(old_end, end))
    return spans


def oracle_counts(
    bam_path: Path,
    gene_spans: dict[str, tuple[str, int, int]],
) -> tuple[Counter[str], Counter[str], int, int]:
    seen: set[str] = set()
    transcript_counts: Counter[str] = Counter()
    gdna_gene_counts: Counter[str] = Counter()
    gdna_count = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.query_name in seen:
                continue
            seen.add(read.query_name)
            if read.query_name.startswith("gdna:"):
                gdna_count += 1
                _, ref, coord, *_ = read.query_name.split(":")
                start_text, end_text = coord.split("-", 1)
                frag_start = int(start_text)
                frag_end = int(end_text)
                assigned = False
                for gene_id, (gene_ref, gene_start, gene_end) in gene_spans.items():
                    if ref == gene_ref and frag_start < gene_end and frag_end > gene_start:
                        gdna_gene_counts[gene_id] += 1
                        assigned = True
                if not assigned:
                    gdna_gene_counts["__intergenic__"] += 1
            else:
                transcript_id = read.query_name.split(":", 1)[0]
                transcript_counts[transcript_id] += 1
    return transcript_counts, gdna_gene_counts, gdna_count, len(seen)


def summary_get(summary: dict, path: str, default=np.nan):
    current = summary
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            return default
        current = current[part]
    return current


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def analyze_condition(
    condition: str,
    truth_df: pd.DataFrame,
    gene_spans: dict[str, tuple[str, int, int]],
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    out = OUTBASE / condition
    quant = pd.read_csv(out / "quant.tsv", sep="\t")
    loci = pd.read_csv(out / "loci.tsv", sep="\t")
    nrna = pd.read_csv(out / "nrna_quant.tsv", sep="\t")
    summary = json.loads((out / "summary.json").read_text())
    true_tx_counter, true_gdna_gene_counter, true_gdna, true_total = oracle_counts(
        BASE / condition / "sim_oracle.bam", gene_spans
    )

    truth = truth_df[["transcript_id", "gene_id", "mrna_abundance"]].copy()
    truth["true_count"] = truth["transcript_id"].map(true_tx_counter).fillna(0.0).astype(float)
    truth_ids = set(truth["transcript_id"])
    observed = quant.loc[
        quant["transcript_id"].isin(truth_ids), ["transcript_id", "gene_id", "count"]
    ]
    observed = observed.rename(columns={"count": "rigel_count"})
    merged = truth.merge(observed, on=["transcript_id", "gene_id"], how="left")
    merged["rigel_count"] = merged["rigel_count"].fillna(0.0)
    merged["error"] = merged["rigel_count"] - merged["true_count"]
    merged["abs_error"] = merged["error"].abs()
    merged["condition"] = condition

    expressed = merged["true_count"] > 0
    zero_truth = ~expressed
    transcript_mae_all = float(merged["abs_error"].mean())
    transcript_mae_expressed = float(merged.loc[expressed, "abs_error"].mean())
    transcript_rmse_all = float(np.sqrt(np.mean(np.square(merged["error"].to_numpy()))))
    false_positive_mass = float(merged.loc[zero_truth, "rigel_count"].sum())
    true_mrna = float(merged["true_count"].sum())
    est_mrna = float(merged["rigel_count"].sum())
    est_gdna = float(loci["gdna"].sum())
    est_nrna = float(nrna["count"].sum()) if "count" in nrna else 0.0

    gene_truth = merged.groupby("gene_id", as_index=False)["true_count"].sum()
    gene_est = merged.groupby("gene_id", as_index=False)["rigel_count"].sum()
    gene_merged = gene_truth.merge(gene_est, on="gene_id", how="outer").fillna(0.0)
    gene_merged["error"] = gene_merged["rigel_count"] - gene_merged["true_count"]

    locus_gene = (
        quant.loc[quant["transcript_id"].isin(set(truth["transcript_id"])), ["gene_id", "locus_id"]]
        .drop_duplicates()
        .groupby("locus_id", as_index=False)["gene_id"]
        .first()
    )
    gdna_by_gene = loci[["locus_id", "gdna", "gdna_prior_count_em", "n_em_fragments"]].merge(
        locus_gene, on="locus_id", how="left"
    )
    gdna_by_gene["true_gdna"] = gdna_by_gene["gene_id"].map(true_gdna_gene_counter).fillna(0.0)
    gdna_by_gene["gdna_error"] = gdna_by_gene["gdna"] - gdna_by_gene["true_gdna"]
    gdna_by_gene["condition"] = condition

    row = {
        "condition": condition,
        "true_total": true_total,
        "true_mrna": true_mrna,
        "true_gdna": float(true_gdna),
        "true_gdna_genic": float(
            sum(v for key, v in true_gdna_gene_counter.items() if key != "__intergenic__")
        ),
        "true_gdna_intergenic": float(true_gdna_gene_counter.get("__intergenic__", 0)),
        "est_mrna": est_mrna,
        "est_nrna": est_nrna,
        "est_gdna": est_gdna,
        "est_total": est_mrna + est_nrna + est_gdna,
        "gdna_error": est_gdna - float(true_gdna),
        "gdna_rel_error": (est_gdna - float(true_gdna)) / max(float(true_gdna), 1.0),
        "mrna_error": est_mrna - true_mrna,
        "mrna_rel_error": (est_mrna - true_mrna) / max(true_mrna, 1.0),
        "nrna_false_mass": est_nrna,
        "false_positive_tx_mass": false_positive_mass,
        "tx_mae_all": transcript_mae_all,
        "tx_mae_expressed": transcript_mae_expressed,
        "tx_rmse_all": transcript_rmse_all,
        "tx_corr": safe_corr(merged["true_count"].to_numpy(), merged["rigel_count"].to_numpy()),
        "gene_mae": float(gene_merged["error"].abs().mean()),
        "gene_rmse": float(np.sqrt(np.mean(np.square(gene_merged["error"].to_numpy())))),
        "gene_corr": safe_corr(
            gene_merged["true_count"].to_numpy(), gene_merged["rigel_count"].to_numpy()
        ),
        "n_regions": summary_get(summary, "calibration.density_evidence.n_regions"),
        "rho_ref_source": summary_get(summary, "calibration.density_evidence.rho_ref_source"),
        "rho_ref": summary_get(summary, "calibration.density_evidence.rho_ref"),
        "density_intergenic_fit": summary_get(
            summary, "calibration.density_evidence.priors.INTERGENIC.fit_status", "none"
        ),
        "density_intron_fit": summary_get(
            summary, "calibration.density_evidence.priors.INTRON.fit_status", "none"
        ),
        "density_all_fit": summary_get(
            summary, "calibration.density_evidence.priors.ALL.fit_status", "none"
        ),
        "density_fallback_regions": summary_get(
            summary, "calibration.density_evidence.flags.n_fallback_used", 0
        ),
        "density_high_tail_regions": summary_get(
            summary, "calibration.density_evidence.flags.n_high_tail_tension", 0
        ),
        "intergenic_contained": summary_get(
            summary, "calibration.diagnostics.fl_pool_total.INTERGENIC_CONTAINED", 0.0
        ),
        "intronic_contained": summary_get(
            summary, "calibration.diagnostics.fl_pool_total.INTRONIC_CONTAINED", 0.0
        ),
        "exonic_contained": summary_get(
            summary, "calibration.diagnostics.fl_pool_total.EXONIC_CONTAINED", 0.0
        ),
        "fl_rna_quality": summary_get(summary, "calibration.fl_models.rna_quality"),
        "fl_gdna_quality": summary_get(summary, "calibration.fl_models.gdna_quality"),
        "fl_n_gdna": summary_get(summary, "calibration.fl_models.n_gdna"),
        "strand_p_r1_sense": summary_get(summary, "strand_model.p_r1_sense"),
        "strand_specificity": summary_get(summary, "strand_model.strand_specificity"),
        "strand_near_unstranded_regions": summary_get(
            summary, "calibration.strand_deconv.n_regions_near_unstranded", 0
        ),
        "kappa_d": summary_get(summary, "calibration.strand_deconv.kappa_d"),
        "fused_mean_region_avg": summary_get(
            summary, "calibration.fused_region_gdna.mean_count.mean"
        ),
        "fused_upper_region_avg": summary_get(
            summary, "calibration.fused_region_gdna.upper_count.mean"
        ),
        "fused_density_only_regions": summary_get(
            summary, "calibration.fused_region_gdna.flags.n_density_only", 0
        ),
        "fused_strand_regions": summary_get(
            summary, "calibration.fused_region_gdna.flags.n_strand_used", 0
        ),
        "fused_boundary_fallback_regions": summary_get(
            summary, "calibration.fused_region_gdna.flags.n_boundary_fallback", 0
        ),
        "prior_gdna_sum": summary_get(
            summary, "calibration.prior_table.sum_gdna_prior_count_em", 0.0
        ),
        "prior_gdna_upper_sum": summary_get(
            summary, "calibration.prior_table.sum_gdna_upper_count", 0.0
        ),
        "prior_unallocated": summary_get(
            summary, "calibration.prior_table.unallocated_expected_count", 0.0
        ),
        "prior_enabled_loci": summary_get(
            summary, "calibration.prior_table.n_loci_enable_gdna_true", 0
        ),
        "loci": len(loci),
        "mega_loci_gt_10k_fragments": int((loci["n_em_fragments"] > 10000).sum()),
    }
    return row, merged.sort_values("abs_error", ascending=False), gdna_by_gene


def main() -> None:
    truth_df = pd.read_csv(TRUTH, sep="\t")
    gene_spans = parse_gene_spans()
    rows: list[dict[str, object]] = []
    top_errors: list[pd.DataFrame] = []
    gdna_gene_errors: list[pd.DataFrame] = []
    for condition in CONDITIONS:
        row, errors, gdna_by_gene = analyze_condition(condition, truth_df, gene_spans)
        rows.append(row)
        top_errors.append(errors.head(8))
        gdna_gene_errors.append(gdna_by_gene.sort_values("gdna_error").head(5))

    summary = pd.DataFrame(rows)
    pd.set_option("display.max_columns", 80)
    pd.set_option("display.width", 220)
    print("# TOTALS_AND_ERRORS")
    print(
        summary[
            [
                "condition",
                "true_mrna",
                "true_gdna",
                "true_gdna_genic",
                "true_gdna_intergenic",
                "est_mrna",
                "est_nrna",
                "est_gdna",
                "mrna_rel_error",
                "gdna_rel_error",
                "tx_corr",
                "gene_corr",
                "tx_mae_expressed",
                "gene_mae",
                "false_positive_tx_mass",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.4f}")
    )
    print("\n# CALIBRATION_SUMMARY")
    print(
        summary[
            [
                "condition",
                "n_regions",
                "rho_ref_source",
                "rho_ref",
                "density_intergenic_fit",
                "density_intron_fit",
                "density_all_fit",
                "density_fallback_regions",
                "density_high_tail_regions",
                "intergenic_contained",
                "intronic_contained",
                "exonic_contained",
                "fl_rna_quality",
                "fl_gdna_quality",
                "fl_n_gdna",
                "strand_specificity",
                "strand_near_unstranded_regions",
                "kappa_d",
                "fused_density_only_regions",
                "fused_strand_regions",
                "fused_boundary_fallback_regions",
                "prior_gdna_sum",
                "prior_gdna_upper_sum",
                "prior_unallocated",
                "prior_enabled_loci",
                "mega_loci_gt_10k_fragments",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.4f}" if np.isfinite(x) else "nan")
    )
    print("\n# TOP_TRANSCRIPT_ERRORS")
    top = pd.concat(top_errors, ignore_index=True)
    print(
        top[
            [
                "condition",
                "transcript_id",
                "gene_id",
                "mrna_abundance",
                "true_count",
                "rigel_count",
                "error",
                "abs_error",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.3f}")
    )
    print("\n# WORST_GDNA_LOCUS_UNDERCALLS")
    gdna_top = pd.concat(gdna_gene_errors, ignore_index=True)
    print(
        gdna_top[
            [
                "condition",
                "gene_id",
                "true_gdna",
                "gdna",
                "gdna_error",
                "gdna_prior_count_em",
                "n_em_fragments",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.3f}")
    )


if __name__ == "__main__":
    main()
