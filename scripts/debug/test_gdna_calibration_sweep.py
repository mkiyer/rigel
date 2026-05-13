#!/usr/bin/env python
"""Synthetic mini-genome gDNA calibration sweep.

Creates a minimal 2-gene genome (one multi-exon, one single-exon @ 0 TPM),
simulates oracle BAMs across a range of gDNA fractions, runs the rigel
pipeline, and performs deep analysis of calibration accuracy and EM prior
settings.

Usage:
    conda activate rigel
    python scripts/debug/test_gdna_calibration_sweep.py \
        --outdir /Users/mkiyer/Downloads/rigel_runs/gdna_sweep_v1
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-7s %(name)s — %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("gdna_sweep")

# ── rigel imports ──────────────────────────────────────────────────────────
from rigel.config import EMConfig, PipelineConfig, BamScanConfig, CalibrationConfig
from rigel.pipeline import run_pipeline
from rigel.calibration.density_global import estimate_global_gdna_fragments
from rigel.sim import Scenario, SimConfig
from rigel.sim.reads import GDNAConfig


# ── Sweep configuration ───────────────────────────────────────────────────

GENOME_LENGTH = 30_000
SEED = 42
N_RNA_FRAGMENTS = 10_000
GDNA_FRACTIONS = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]

# Transcript definitions
MULTI_EXON_TRANSCRIPT = {
    "t_id": "tx_multi",
    "exons": [(10000, 10500), (12000, 13000)],
    "abundance": 100.0,
    "nrna_abundance": 0.0,
}
SINGLE_EXON_TRANSCRIPT = {
    "t_id": "tx_single",
    "exons": [(20000, 22000)],
    "abundance": 0.0,       # NOT expressed
    "nrna_abundance": 0.0,
}

SIM_CONFIG = SimConfig(
    frag_mean=250.0,
    frag_std=50.0,
    frag_min=100,
    frag_max=600,
    read_length=150,
    error_rate=0.0,
    strand_specificity=1.0,
    seed=SEED,
)

GDNA_CONFIG = GDNAConfig(
    abundance=1.0,       # placeholder; pool_split controls actual counts
    frag_mean=350.0,
    frag_std=100.0,
    frag_min=100,
    frag_max=1000,
)


# ── Results container ─────────────────────────────────────────────────────

@dataclass
class SweepRow:
    gdna_fraction: float
    n_rna_truth: int
    n_gdna_truth: int
    n_total_truth: int
    gdna_rate_truth: float

    # Calibration outputs
    cal_n_exon_only: int
    cal_n_intron_only: int
    cal_n_intergenic_only: int
    cal_n_exon_intron: int
    cal_n_unannotated: int
    cal_total_observed: int

    cal_rho_intergenic: float
    cal_rho_intron: float
    cal_rho_exon_intron: float
    cal_gdna_fl_mean: float
    cal_rna_fl_mean: float

    cal_n_fl_global: int
    cal_n_fl_rna: int
    cal_n_fl_gdna: int
    cal_rna_fl_quality: str
    cal_gdna_fl_quality: str

    # Per-locus prior outputs
    cal_n_multi_loci: int
    cal_mean_pi_gdna: float

    # EM outputs (locus-level)
    em_mrna_count: float
    em_nrna_count: float
    em_gdna_count: float
    em_gdna_rate: float

    # Global gDNA (calibration density projection)
    gdna_global_count: float
    gdna_global_rate: float

    # Per-transcript EM
    em_tx_multi_count: float
    em_tx_multi_tpm: float
    em_tx_single_count: float
    em_tx_single_tpm: float

    # Locus-level priors
    locus_gdna_prior_count: str     # JSON array

    # Errors
    mrna_count_error: float   # (estimated - truth) / truth
    gdna_count_error: float   # (estimated - truth) / truth  (NaN if truth=0)
    gdna_global_error: float  # (global_projected - truth) / truth
    gdna_rate_error: float    # absolute: em_rate - true_rate
    gdna_global_rate_error: float  # absolute: global_rate - true_rate


# ── Run one scenario ──────────────────────────────────────────────────────

def run_one(gdna_frac: float, outdir: Path) -> SweepRow:
    """Build oracle BAM, run pipeline, extract metrics for one gDNA fraction."""
    tag = f"gdna_{gdna_frac:.2f}"
    work_dir = outdir / tag
    work_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"═══ {tag}: building scenario ═══")

    # Build scenario
    sc = Scenario(
        tag,
        genome_length=GENOME_LENGTH,
        seed=SEED,
        work_dir=work_dir,
        gdna_config=GDNA_CONFIG if gdna_frac > 0 else None,
    )
    sc.add_gene("gene_multi", "+", [MULTI_EXON_TRANSCRIPT])
    sc.add_gene("gene_single", "+", [SINGLE_EXON_TRANSCRIPT])

    # Build oracle BAM with explicit RNA/gDNA split
    result = sc.build_oracle(
        n_rna_fragments=N_RNA_FRAGMENTS,
        gdna_fraction=gdna_frac if gdna_frac > 0 else None,
        sim_config=SIM_CONFIG,
        gdna_config=GDNA_CONFIG if gdna_frac > 0 else None,
    )

    # Ground truth
    gt_mrna = result.ground_truth_from_bam()
    gt_gdna = result.ground_truth_gdna_count_from_bam()
    gt_nrna = result.ground_truth_nrna_count_from_bam()
    n_rna_truth = sum(gt_mrna.values()) + gt_nrna
    n_total = n_rna_truth + gt_gdna
    gdna_rate_truth = gt_gdna / n_total if n_total > 0 else 0.0

    logger.info(
        f"  Ground truth: mRNA={sum(gt_mrna.values())}, nRNA={gt_nrna}, "
        f"gDNA={gt_gdna}, total={n_total}, rate={gdna_rate_truth:.4f}"
    )

    # Run pipeline
    logger.info(f"  Running rigel pipeline...")
    config = PipelineConfig(
        em=EMConfig(seed=SEED, mode="vbem"),
        scan=BamScanConfig(sj_strand_tag="auto"),
        calibration=CalibrationConfig(),
    )
    pr = run_pipeline(result.bam_path, result.index, config=config)

    # ── Extract calibration metrics ───────────────────────────────────
    cal = pr.calibration
    cal_summary = cal.to_summary_dict()
    diag = cal.diagnostics
    gd = cal.global_densities
    fl = cal.fl_models

    # ── Extract EM / quantification metrics ───────────────────────────
    quant_df = pr.estimator.get_counts_df(result.index)
    loci_df = pr.estimator.get_loci_df(result.index)

    # Transcript-level
    tx_multi = quant_df[quant_df["transcript_id"] == "tx_multi"]
    tx_single = quant_df[quant_df["transcript_id"] == "tx_single"]

    em_tx_multi_count = float(tx_multi["count_em"].iloc[0]) if len(tx_multi) else 0.0
    em_tx_multi_tpm = float(tx_multi["tpm"].iloc[0]) if len(tx_multi) else 0.0
    em_tx_single_count = float(tx_single["count_em"].iloc[0]) if len(tx_single) else 0.0
    em_tx_single_tpm = float(tx_single["tpm"].iloc[0]) if len(tx_single) else 0.0

    # Locus-level totals
    em_mrna = float(loci_df["mrna"].sum())
    em_nrna = float(loci_df["nrna"].sum())
    em_gdna = float(loci_df["gdna"].sum())
    em_total = em_mrna + em_nrna + em_gdna
    em_gdna_rate = em_gdna / em_total if em_total > 0 else 0.0

    # Global projected gDNA (from calibration densities over full genome)
    region_df = getattr(result.index, "region_df", None)
    if region_df is not None and len(region_df) > 0:
        gdna_global = estimate_global_gdna_fragments(cal.global_densities, region_df)
    else:
        gdna_global = em_gdna
    gdna_global_total = em_mrna + em_nrna + gdna_global
    gdna_global_rate = gdna_global / gdna_global_total if gdna_global_total > 0 else 0.0

    # Locus-level priors
    gdna_prior_count = (
        loci_df["gdna_prior_count"].tolist()
        if "gdna_prior_count" in loci_df.columns
        else []
    )

    # ── Compute errors ────────────────────────────────────────────────
    mrna_truth = sum(gt_mrna.values())
    mrna_error = (em_mrna - mrna_truth) / mrna_truth if mrna_truth > 0 else float("nan")
    gdna_error = (em_gdna - gt_gdna) / gt_gdna if gt_gdna > 0 else float("nan")
    gdna_global_error = (gdna_global - gt_gdna) / gt_gdna if gt_gdna > 0 else float("nan")
    gdna_rate_err = em_gdna_rate - gdna_rate_truth
    gdna_global_rate_err = gdna_global_rate - gdna_rate_truth

    # ── Save per-scenario detail ──────────────────────────────────────
    quant_df.to_csv(work_dir / "quant.tsv", sep="\t", index=False)
    loci_df.to_csv(work_dir / "loci.tsv", sep="\t", index=False)
    with open(work_dir / "calibration_summary.json", "w") as f:
        json.dump(cal_summary, f, indent=2, default=str)
    with open(work_dir / "ground_truth.json", "w") as f:
        json.dump({
            "mrna_counts": gt_mrna,
            "gdna_count": gt_gdna,
            "nrna_count": gt_nrna,
        }, f, indent=2)

    # Multi-locus prior DF
    if len(cal.multi_locus_prior_df) > 0:
        cal.multi_locus_prior_df.to_csv(
            work_dir / "multi_locus_prior.tsv", sep="\t", index=False
        )
    if len(cal.per_locus_gdna_df) > 0:
        cal.per_locus_gdna_df.to_csv(
            work_dir / "per_locus_gdna.tsv", sep="\t", index=False
        )

    logger.info(
        f"  Pipeline: mRNA_em={em_mrna:.1f} gDNA_em={em_gdna:.1f} "
        f"gDNA_global={gdna_global:.1f} "
        f"rate_global={gdna_global_rate:.4f} rate_true={gdna_rate_truth:.4f}"
    )

    return SweepRow(
        gdna_fraction=gdna_frac,
        n_rna_truth=n_rna_truth,
        n_gdna_truth=gt_gdna,
        n_total_truth=n_total,
        gdna_rate_truth=gdna_rate_truth,

        cal_n_exon_only=diag.n_exon_only,
        cal_n_intron_only=diag.n_intron_only,
        cal_n_intergenic_only=diag.n_intergenic_only,
        cal_n_exon_intron=diag.n_exon_intron,
        cal_n_unannotated=diag.n_unannotated,
        cal_total_observed=diag.total(),

        cal_rho_intergenic=gd.intergenic.rho,
        cal_rho_intron=gd.intron.rho,
        cal_rho_exon_intron=gd.exon_intron.rho,
        cal_gdna_fl_mean=fl.gdna.mean,
        cal_rna_fl_mean=fl.rna.mean,

        cal_n_fl_global=fl.n_global,
        cal_n_fl_rna=fl.n_rna,
        cal_n_fl_gdna=fl.n_gdna,
        cal_rna_fl_quality=fl.rna_quality,
        cal_gdna_fl_quality=fl.gdna_quality,

        cal_n_multi_loci=cal.n_multi_loci,
        cal_mean_pi_gdna=cal_summary.get("mean_pi_gdna", 0.0),

        em_mrna_count=em_mrna,
        em_nrna_count=em_nrna,
        em_gdna_count=em_gdna,
        em_gdna_rate=em_gdna_rate,

        gdna_global_count=gdna_global,
        gdna_global_rate=gdna_global_rate,

        em_tx_multi_count=em_tx_multi_count,
        em_tx_multi_tpm=em_tx_multi_tpm,
        em_tx_single_count=em_tx_single_count,
        em_tx_single_tpm=em_tx_single_tpm,

        locus_gdna_prior_count=json.dumps([float(x) for x in gdna_prior_count]),

        mrna_count_error=mrna_error,
        gdna_count_error=gdna_error,
        gdna_global_error=gdna_global_error,
        gdna_rate_error=gdna_rate_err,
        gdna_global_rate_error=gdna_global_rate_err,
    )


# ── Deep analysis ─────────────────────────────────────────────────────────

def deep_analysis(df: pd.DataFrame, outdir: Path) -> str:
    """Perform deep analysis of calibration and EM results across sweep."""
    lines: list[str] = []
    hr = "─" * 80

    lines.append("")
    lines.append("=" * 80)
    lines.append("  DEEP ANALYSIS: gDNA Calibration Sweep")
    lines.append("=" * 80)

    # ── 1. Ground truth summary ───────────────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("1. GROUND TRUTH SUMMARY")
    lines.append(hr)
    lines.append(f"{'gDNA_frac':>10} {'n_RNA':>8} {'n_gDNA':>8} {'n_total':>8} {'rate_true':>10}")
    for _, r in df.iterrows():
        lines.append(
            f"{r.gdna_fraction:>10.2f} {r.n_rna_truth:>8} {r.n_gdna_truth:>8} "
            f"{r.n_total_truth:>8} {r.gdna_rate_truth:>10.4f}"
        )

    # ── 2. Calibration diagnostics ────────────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("2. CALIBRATION FRAGMENT CLASSIFICATION (diagnostics)")
    lines.append(hr)
    lines.append(
        f"{'gDNA_frac':>10} {'exon':>8} {'intron':>8} {'intergenic':>10} "
        f"{'exon_intr':>10} {'unannot':>8} {'total':>8}"
    )
    for _, r in df.iterrows():
        lines.append(
            f"{r.gdna_fraction:>10.2f} {r.cal_n_exon_only:>8} {r.cal_n_intron_only:>8} "
            f"{r.cal_n_intergenic_only:>10} {r.cal_n_exon_intron:>10} "
            f"{r.cal_n_unannotated:>8} {r.cal_total_observed:>8}"
        )
    # Interpretation
    lines.append("\n  → Exon-only should be ≈ n_RNA (all mRNA is exonic).")
    lines.append("  → Intron + intergenic + exon-intron ≈ gDNA (uniform genome sampling).")
    lines.append("  → Unannotated should be 0 (all reads fall on the single reference).")

    # ── 3. Global gDNA density estimates ──────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("3. GLOBAL gDNA DENSITY ESTIMATES (ρ = frags/bp)")
    lines.append(hr)
    lines.append(
        f"{'gDNA_frac':>10} {'ρ_intergenic':>14} {'ρ_intron':>14} "
        f"{'ρ_exon_intron':>14} {'gDNA_FL_mean':>12} {'RNA_FL_mean':>12}"
    )
    for _, r in df.iterrows():
        lines.append(
            f"{r.gdna_fraction:>10.2f} {r.cal_rho_intergenic:>14.6f} "
            f"{r.cal_rho_intron:>14.6f} {r.cal_rho_exon_intron:>14.6f} "
            f"{r.cal_gdna_fl_mean:>12.1f} {r.cal_rna_fl_mean:>12.1f}"
        )
    lines.append("\n  → Densities should scale proportionally with gDNA fraction.")
    lines.append("  → All three density types should be roughly consistent (same genome).")
    lines.append("  → gDNA FL mean should be ≈350 (simulated); RNA FL mean ≈250.")

    # ── 4. Fragment-length model quality ──────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("4. FRAGMENT-LENGTH MODEL QUALITY")
    lines.append(hr)
    lines.append(
        f"{'gDNA_frac':>10} {'n_global':>10} {'n_rna':>10} {'n_gdna':>10} "
        f"{'rna_qual':>10} {'gdna_qual':>10}"
    )
    for _, r in df.iterrows():
        lines.append(
            f"{r.gdna_fraction:>10.2f} {r.cal_n_fl_global:>10} {r.cal_n_fl_rna:>10} "
            f"{r.cal_n_fl_gdna:>10} {r.cal_rna_fl_quality:>10} {r.cal_gdna_fl_quality:>10}"
        )
    lines.append("\n  → At low gDNA, gdna_quality may be 'fallback' (few gDNA-classifiable frags).")
    lines.append("  → RNA pool should always be 'good' (10k RNA fragments).")

    # ── 5. Per-locus gDNA priors ──────────────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("5. PER-LOCUS gDNA PRIORS")
    lines.append(hr)
    lines.append(
        f"{'gDNA_frac':>10} {'n_multi_loci':>13} {'mean_π_gdna':>12} "
        f"{'gdna_prior_count':>40}"
    )
    for _, r in df.iterrows():
        lines.append(
            f"{r.gdna_fraction:>10.2f} {r.cal_n_multi_loci:>13} "
            f"{r.cal_mean_pi_gdna:>12.6f} "
            f"{r.locus_gdna_prior_count:>40}"
        )
    lines.append("\n  → π_gdna should track true gDNA rate closely.")
    lines.append("  → gdna_prior_count is the expected gDNA pseudocount per locus.")

    # ── 6. EM quantification accuracy ─────────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("6. EM QUANTIFICATION ACCURACY")
    lines.append(hr)
    lines.append(
        f"{'gDNA_frac':>10} {'mRNA_em':>10} {'mRNA_true':>10} {'mRNA_err%':>10} "
        f"{'gDNA_em':>10} {'gDNA_glob':>10} {'gDNA_true':>10} "
        f"{'em_err%':>8} {'glob_err%':>10} "
        f"{'rate_glob':>10} {'rate_true':>10} {'Δrate':>8}"
    )
    for _, r in df.iterrows():
        mrna_err_pct = f"{r.mrna_count_error * 100:+.2f}" if not np.isnan(r.mrna_count_error) else "N/A"
        gdna_err_pct = f"{r.gdna_count_error * 100:+.2f}" if not np.isnan(r.gdna_count_error) else "N/A"
        gdna_glob_err_pct = f"{r.gdna_global_error * 100:+.2f}" if not np.isnan(r.gdna_global_error) else "N/A"
        lines.append(
            f"{r.gdna_fraction:>10.2f} {r.em_mrna_count:>10.1f} {r.n_rna_truth:>10} "
            f"{mrna_err_pct:>10} {r.em_gdna_count:>10.1f} {r.gdna_global_count:>10.1f} "
            f"{r.n_gdna_truth:>10} {gdna_err_pct:>8} {gdna_glob_err_pct:>10} "
            f"{r.gdna_global_rate:>10.4f} {r.gdna_rate_truth:>10.4f} "
            f"{r.gdna_global_rate_error:>+8.4f}"
        )

    # ── 7. Per-transcript breakdown ───────────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("7. PER-TRANSCRIPT EM COUNTS")
    lines.append(hr)
    lines.append(
        f"{'gDNA_frac':>10} {'tx_multi_cnt':>13} {'tx_multi_tpm':>13} "
        f"{'tx_single_cnt':>14} {'tx_single_tpm':>14} {'nRNA_em':>10}"
    )
    for _, r in df.iterrows():
        lines.append(
            f"{r.gdna_fraction:>10.2f} {r.em_tx_multi_count:>13.1f} "
            f"{r.em_tx_multi_tpm:>13.1f} {r.em_tx_single_count:>14.3f} "
            f"{r.em_tx_single_tpm:>14.3f} {r.em_nrna_count:>10.1f}"
        )
    lines.append("\n  → tx_multi should get ≈10000 mRNA counts (the only expressed transcript).")
    lines.append("  → tx_single should get ≈0 (not expressed; any count is gDNA leakage).")
    lines.append("  → nRNA should be ≈0 (not simulated).")

    # ── 8. Consistency checks ─────────────────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("8. CONSISTENCY & ISSUE DETECTION")
    lines.append(hr)

    issues = []
    for _, r in df.iterrows():
        tag = f"gDNA_frac={r.gdna_fraction:.2f}"

        # Check: total observed ≈ total truth
        if r.cal_total_observed > 0:
            obs_err = abs(r.cal_total_observed - r.n_total_truth) / r.n_total_truth
            if obs_err > 0.05:
                issues.append(
                    f"  ⚠ {tag}: observed fragments ({r.cal_total_observed}) differ "
                    f"from truth ({r.n_total_truth}) by {obs_err*100:.1f}%"
                )

        # Check: mRNA error
        if not np.isnan(r.mrna_count_error) and abs(r.mrna_count_error) > 0.05:
            issues.append(
                f"  ⚠ {tag}: mRNA count error = {r.mrna_count_error*100:+.2f}% "
                f"(em={r.em_mrna_count:.1f}, truth={r.n_rna_truth})"
            )

        # Check: gDNA rate error
        if abs(r.gdna_rate_error) > 0.02:
            issues.append(
                f"  ⚠ {tag}: gDNA rate error = {r.gdna_rate_error:+.4f} "
                f"(em={r.em_gdna_rate:.4f}, truth={r.gdna_rate_truth:.4f})"
            )

        # Check: tx_single leakage
        if r.em_tx_single_count > 1.0:
            issues.append(
                f"  ⚠ {tag}: tx_single (unexpressed) received {r.em_tx_single_count:.1f} "
                f"EM counts — possible gDNA leakage into mRNA"
            )

        # Check: nRNA phantom
        if r.em_nrna_count > 1.0:
            issues.append(
                f"  ⚠ {tag}: nRNA received {r.em_nrna_count:.1f} EM counts "
                f"despite zero nRNA simulation"
            )

        # Check: density consistency
        if r.gdna_fraction > 0:
            densities = [r.cal_rho_intergenic, r.cal_rho_intron, r.cal_rho_exon_intron]
            densities_nonzero = [d for d in densities if d > 0]
            if len(densities_nonzero) >= 2:
                ratio = max(densities_nonzero) / min(densities_nonzero)
                if ratio > 3.0:
                    issues.append(
                        f"  ⚠ {tag}: gDNA density estimates inconsistent "
                        f"(max/min ratio = {ratio:.1f}x): "
                        f"intergenic={r.cal_rho_intergenic:.6f}, "
                        f"intron={r.cal_rho_intron:.6f}, "
                        f"exon_intron={r.cal_rho_exon_intron:.6f}"
                    )

        # Check: FL model quality
        if r.gdna_fraction >= 0.25 and r.cal_gdna_fl_quality == "fallback":
            issues.append(
                f"  ⚠ {tag}: gDNA FL model is 'fallback' despite substantial gDNA "
                f"(n_gdna_fl={r.cal_n_fl_gdna})"
            )

    if issues:
        lines.append("\nDETECTED ISSUES:")
        lines.extend(issues)
    else:
        lines.append("\n  ✓ No issues detected across all sweep points.")

    # ── 9. Summary statistics ─────────────────────────────────────────
    lines.append(f"\n{hr}")
    lines.append("9. SUMMARY STATISTICS")
    lines.append(hr)

    has_gdna = df[df.gdna_fraction > 0]
    if len(has_gdna) > 0:
        mrna_errs = has_gdna["mrna_count_error"].dropna()
        gdna_errs = has_gdna["gdna_count_error"].dropna()
        gdna_glob_errs = has_gdna["gdna_global_error"].dropna()
        rate_errs = has_gdna["gdna_rate_error"].abs()
        glob_rate_errs = has_gdna["gdna_global_rate_error"].abs()

        lines.append(f"  Mean |mRNA count error|:        {mrna_errs.abs().mean()*100:.2f}%")
        lines.append(f"  Max  |mRNA count error|:        {mrna_errs.abs().max()*100:.2f}%")
        lines.append(f"  Mean |gDNA count error (EM)|:   {gdna_errs.abs().mean()*100:.2f}%")
        lines.append(f"  Max  |gDNA count error (EM)|:   {gdna_errs.abs().max()*100:.2f}%")
        lines.append(f"  Mean |gDNA count error (glob)|: {gdna_glob_errs.abs().mean()*100:.2f}%")
        lines.append(f"  Max  |gDNA count error (glob)|: {gdna_glob_errs.abs().max()*100:.2f}%")
        lines.append(f"  Mean |gDNA rate error (EM)|:    {rate_errs.mean():.4f}")
        lines.append(f"  Max  |gDNA rate error (EM)|:    {rate_errs.max():.4f}")
        lines.append(f"  Mean |gDNA rate error (glob)|:  {glob_rate_errs.mean():.4f}")
        lines.append(f"  Max  |gDNA rate error (glob)|:  {glob_rate_errs.max():.4f}")

    zero_gdna = df[df.gdna_fraction == 0.0]
    if len(zero_gdna) > 0:
        r0 = zero_gdna.iloc[0]
        lines.append(f"\n  Zero-gDNA baseline:")
        lines.append(f"    mRNA EM count: {r0.em_mrna_count:.1f} (truth: {r0.n_rna_truth})")
        lines.append(f"    gDNA EM count: {r0.em_gdna_count:.1f} (truth: 0)")
        lines.append(f"    Phantom gDNA rate: {r0.em_gdna_rate:.6f}")

    lines.append("")
    lines.append("=" * 80)
    lines.append("  END OF ANALYSIS")
    lines.append("=" * 80)
    lines.append("")

    return "\n".join(lines)


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="gDNA calibration sweep test")
    parser.add_argument(
        "--outdir", type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/gdna_sweep_v1"),
        help="Output directory",
    )
    args = parser.parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    rows: list[SweepRow] = []
    for gf in GDNA_FRACTIONS:
        row = run_one(gf, outdir)
        rows.append(row)

    # Build summary DataFrame
    df = pd.DataFrame([asdict(r) for r in rows])
    summary_path = outdir / "sweep_summary.tsv"
    df.to_csv(summary_path, sep="\t", index=False)
    logger.info(f"Summary saved to {summary_path}")

    # Deep analysis
    report = deep_analysis(df, outdir)
    report_path = outdir / "analysis_report.txt"
    with open(report_path, "w") as f:
        f.write(report)
    logger.info(f"Analysis report saved to {report_path}")

    # Print report
    print(report)

    return 0


if __name__ == "__main__":
    sys.exit(main())
