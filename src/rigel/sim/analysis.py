#!/usr/bin/env python3
"""Run rigel quantification on synthetic simulation conditions and analyze results.

Builds rigel index, runs quant with annotated BAM on all conditions,
then performs comprehensive analysis of:
  1. Calibration accuracy (gDNA FL, pool separation, locus-level gDNA)
  2. Transcript abundance accuracy
  3. Fragment-level misallocation using oracle BAM ground truth
"""

import argparse
import json
import logging
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .manifest import condition_manifest_map, load_manifest
from .truth import parse_origin

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("rigel_analysis")

DEFAULT_SIM_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic")

# gDNA fractions from simulation config
GDNA_FRACS = {
    "gdna_none": 0.0,
    "gdna_low": 0.1,
    "gdna_med": 0.5,
    "gdna_equal": 1.0,
    "gdna_high": 2.0,
}


@dataclass(frozen=True)
class CalibrationAcceptanceThresholds:
    """Post-fix guardrails for the synthetic calibration report."""

    min_rho_ex_over_ig: float = 0.95
    max_nrna_none_count: float = 10.0
    max_gdna_to_rna_leak_rate: float = 0.015


def parse_condition(cond_name: str) -> dict:
    """Parse condition name into components."""
    parts = cond_name.split("_")
    # gdna_XXX_ss_YYY_nrna_ZZZ
    gdna_key = f"{parts[0]}_{parts[1]}"
    ss_val = float(parts[3])
    nrna_label = parts[5] if len(parts) > 5 else "none"
    return {
        "condition": cond_name,
        "gdna_frac": GDNA_FRACS.get(gdna_key, 0.0),
        "strand_specificity": ss_val,
        "gdna_label": gdna_key,
        "nrna_label": nrna_label,
    }


def discover_conditions(sim_base: Path, selected: list[str] | None = None) -> list[str]:
    """Return requested conditions or manifest-discovered conditions."""
    if selected:
        return selected
    manifest = load_manifest(sim_base)
    mapped = condition_manifest_map(manifest)
    if mapped:
        return list(mapped)
    raise ValueError(f"No conditions selected and no manifest conditions found in {sim_base}")


def simulated_fragment_length_means(manifest: dict) -> tuple[float, float]:
    """Return true RNA/gDNA fragment-length means from manifest metadata."""
    sim = manifest.get("simulation", {}) if isinstance(manifest, dict) else {}
    gdna = manifest.get("gdna", {}) if isinstance(manifest, dict) else {}
    return float(sim.get("frag_mean", 250.0)), float(gdna.get("frag_mean", 350.0))


def _format_bp(value: float) -> str:
    """Compact bp value for report prose."""
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def _format_rel_err(value: float) -> str:
    if np.isnan(value):
        return "n/a"
    return f"{value:+.3f}"


# ══════════════════════════════════════════════════════════════════════════
# STEP 1: Build rigel index
# ══════════════════════════════════════════════════════════════════════════

def build_index(sim_base: Path) -> Path:
    """Build rigel index from synthetic reference."""
    ref_dir = sim_base / "reference"
    index_dir = sim_base / "rigel_index"

    if (index_dir / "transcripts.feather").exists():
        logger.info(f"Index already exists at {index_dir}")
        return index_dir

    logger.info("Building rigel index...")
    cmd = [
        "rigel", "index",
        "--fasta", str(ref_dir / "genome.fa"),
        "--gtf", str(ref_dir / "genes.gtf"),
        "-o", str(index_dir),
        "--gtf-parse-mode", "warn-skip",
        "--no-mappability",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Index build failed:\n{result.stderr}")
        sys.exit(1)
    logger.info(f"Index built at {index_dir}")
    return index_dir


# ══════════════════════════════════════════════════════════════════════════
# STEP 2: Run rigel quant on each condition
# ══════════════════════════════════════════════════════════════════════════

def run_quant(sim_base: Path, index_dir: Path, condition: str) -> Path:
    """Run rigel quant on one condition with annotated BAM."""
    cond_dir = sim_base / condition
    out_dir = cond_dir / "rigel_out"
    annotated_bam = cond_dir / "annotated.bam"

    if (out_dir / "quant.feather").exists():
        logger.info(f"  [{condition}] Already quantified, skipping")
        return out_dir

    logger.info(f"  [{condition}] Running rigel quant...")
    cmd = [
        "rigel", "quant",
        "--bam", str(cond_dir / "sim_oracle.bam"),
        "--index", str(index_dir),
        "-o", str(out_dir),
        "--annotated-bam", str(annotated_bam),
        "--sj-strand-tag", "auto",
        "--emit-locus-stats",
        "--tsv",
    ]
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    if result.returncode != 0:
        logger.error(f"  [{condition}] quant FAILED ({elapsed:.1f}s):\n{result.stderr[-2000:]}")
        return None

    logger.info(f"  [{condition}] Done ({elapsed:.1f}s)")
    return out_dir


# ══════════════════════════════════════════════════════════════════════════
# STEP 3: Analysis functions
# ══════════════════════════════════════════════════════════════════════════

def load_truth(sim_base: Path, truth_name: str | None = None) -> pd.DataFrame:
    """Load a ground-truth abundance table."""
    manifest = load_manifest(sim_base)
    if truth_name is None:
        truth_name = manifest.get("truth_abundances", "truth_abundances.tsv")
    truth_path = sim_base / truth_name
    if not truth_path.exists():
        condition_map = condition_manifest_map(manifest)
        truth_names = [
            str(row.get("truth_abundances"))
            for row in condition_map.values()
            if row.get("truth_abundances")
        ]
        if truth_names:
            truth_path = sim_base / truth_names[0]
    return pd.read_csv(truth_path, sep="\t")


def load_condition_truth(
    sim_base: Path,
    condition: str,
    condition_meta: dict[str, dict],
    fallback_truth: pd.DataFrame,
) -> pd.DataFrame:
    """Load condition-specific truth, falling back to the supplied table."""
    truth_name = condition_meta.get(condition, {}).get("truth_abundances")
    if not truth_name:
        return fallback_truth
    truth_path = sim_base / str(truth_name)
    if not truth_path.exists():
        return fallback_truth
    return pd.read_csv(truth_path, sep="\t")


def load_quant(out_dir: Path) -> pd.DataFrame:
    """Load rigel quant output."""
    # Try feather first, fall back to TSV
    feather_path = out_dir / "quant.feather"
    tsv_path = out_dir / "quant.tsv"
    if feather_path.exists():
        return pd.read_feather(feather_path)
    elif tsv_path.exists():
        return pd.read_csv(tsv_path, sep="\t")
    else:
        raise FileNotFoundError(f"No quant output found in {out_dir}")


def load_loci(out_dir: Path) -> pd.DataFrame:
    """Load loci output."""
    feather_path = out_dir / "loci.feather"
    tsv_path = out_dir / "loci.tsv"
    if feather_path.exists():
        return pd.read_feather(feather_path)
    elif tsv_path.exists():
        return pd.read_csv(tsv_path, sep="\t")
    else:
        return pd.DataFrame()


def load_summary(out_dir: Path) -> dict:
    """Load pipeline summary JSON."""
    summary_path = out_dir / "summary.json"
    if summary_path.exists():
        with open(summary_path) as f:
            return json.load(f)
    return {}


# ── Calibration Analysis ──────────────────────────────────────────────────

def analyze_calibration(sim_base: Path, conditions: list[str], truth: pd.DataFrame) -> str:
    """Analyze calibration accuracy across all conditions."""
    lines = []
    hr = "═" * 100
    manifest = load_manifest(sim_base)
    condition_meta = condition_manifest_map(manifest)
    true_rna_fl_mean, true_gdna_fl_mean = simulated_fragment_length_means(manifest)

    lines.append(f"\n{hr}")
    lines.append("  CALIBRATION ACCURACY ANALYSIS")
    lines.append(hr)

    # Table: condition → calibration metrics
    lines.append(f"\n{'Condition':<35} {'gDNA_frac':>9} {'SS':>5} "
                 f"{'ρ_ig':>10} {'ρ_in':>10} {'ρ_ex':>10} "
                 f"{'FL_rna':>7} {'FL_gdna':>8} "
                 f"{'n_loci':>6} {'gdna_rate':>10} {'gdna_true':>10}")
    lines.append("─" * 130)

    cal_rows = []
    for cond in conditions:
        info = parse_condition(cond)
        meta = condition_meta.get(cond, {})
        if meta:
            info["gdna_frac"] = float(meta.get("gdna_rate", info["gdna_frac"]))
            info["strand_specificity"] = float(
                meta.get("strand_specificity", info["strand_specificity"]),
            )
            info["nrna_label"] = str(meta.get("nrna_label", info["nrna_label"]))
        out_dir = sim_base / cond / "rigel_out"
        if not out_dir.exists():
            continue

        summary = load_summary(out_dir)
        if not summary:
            continue

        cal = summary.get("calibration", {})
        density_priors = cal.get("density_evidence", {}).get("priors", {}) or {}
        fl = cal.get("fl_models", {})

        def _prior_mean(name: str) -> float:
            prior = density_priors.get(name) or {}
            try:
                return float(prior.get("mean_density", 0.0))
            except (TypeError, ValueError):
                return 0.0

        rho_ig = _prior_mean("INTERGENIC")
        rho_in = _prior_mean("INTRON")
        rho_ex = 0.0  # v4: EXON-INTRON density family removed
        fl_rna = fl.get("rna_fl_mean", 0)
        fl_gdna = fl.get("gdna_fl_mean", 0)
        n_loci = cal.get("n_multi_loci", 0)

        # Compute true gDNA rate
        gdna_frac = info["gdna_frac"]
        n_rna = int(meta.get("n_rna", 1_000_000))
        n_gdna_true = int(meta.get("n_gdna", round(n_rna * gdna_frac)))
        n_total_true = n_rna + n_gdna_true
        gdna_rate_true = n_gdna_true / n_total_true if n_total_true > 0 else 0

        # Estimated gDNA rate from quantification output
        quant_out = summary.get("quantification", {})
        gdna_rate_est = quant_out.get("gdna_fraction", 0)

        lines.append(
            f"{cond:<35} {gdna_frac:>9.2f} {info['strand_specificity']:>5.2f} "
            f"{rho_ig:>10.6f} {rho_in:>10.6f} {rho_ex:>10.6f} "
            f"{fl_rna:>7.1f} {fl_gdna:>8.1f} "
            f"{n_loci:>6} {gdna_rate_est:>10.4f} {gdna_rate_true:>10.4f}"
        )

        cal_rows.append({
            **info,
            "rho_ig": rho_ig, "rho_in": rho_in, "rho_ex": rho_ex,
            "fl_rna": fl_rna, "fl_gdna": fl_gdna,
            "n_loci": n_loci,
            "gdna_rate_est": gdna_rate_est,
            "gdna_rate_true": gdna_rate_true,
            "n_rna_true": n_rna,
            "n_gdna_true": n_gdna_true,
            "intergenic_total": quant_out.get("intergenic_total", 0),
            "gdna_total": quant_out.get("gdna_total", 0),
            "mrna_total": quant_out.get("mrna_total", 0),
        })

    # ── Density ratio consistency ──
    lines.append(f"\n\n{'─' * 100}")
    lines.append("  DENSITY RATIO CONSISTENCY")
    lines.append(f"{'─' * 100}")
    lines.append("  gDNA density should be uniform across intergenic/intron/exon regions.")
    lines.append("  Ratios close to 1.0 indicate consistent calibration.\n")
    lines.append(f"  {'Condition':<35} {'ρ_in/ρ_ig':>10} {'ρ_ex/ρ_ig':>10} {'ρ_ex/ρ_in':>10}")
    for row in cal_rows:
        if row["rho_ig"] > 0:
            r_in_ig = row["rho_in"] / row["rho_ig"]
            r_ex_ig = row["rho_ex"] / row["rho_ig"]
            r_ex_in = row["rho_ex"] / row["rho_in"] if row["rho_in"] > 0 else float("nan")
            lines.append(f"  {row['condition']:<35} {r_in_ig:>10.3f} {r_ex_ig:>10.3f} {r_ex_in:>10.3f}")

    # ── gDNA FL distribution ──
    lines.append(f"\n\n{'─' * 100}")
    lines.append("  FRAGMENT LENGTH MODELS")
    lines.append(f"{'─' * 100}")
    lines.append(
        "  Simulated FL: "
        f"RNA={_format_bp(true_rna_fl_mean)}bp (mean), "
        f"gDNA={_format_bp(true_gdna_fl_mean)}bp (mean). "
        "Check calibration recovery.\n"
    )
    lines.append(f"  {'Condition':<35} {'FL_rna':>8} {'FL_gdna':>8} {'FL_rna_err':>10} {'FL_gdna_err':>11}")
    for row in cal_rows:
        rna_err = (
            (row["fl_rna"] - true_rna_fl_mean) / true_rna_fl_mean
            if row["fl_rna"] > 0 and true_rna_fl_mean > 0
            else float("nan")
        )
        gdna_err = (
            (row["fl_gdna"] - true_gdna_fl_mean) / true_gdna_fl_mean
            if row["n_gdna_true"] > 0 and row["fl_gdna"] > 0 and true_gdna_fl_mean > 0
            else float("nan")
        )
        lines.append(
            f"  {row['condition']:<35} {row['fl_rna']:>8.1f} {row['fl_gdna']:>8.1f} "
            f"{_format_rel_err(rna_err):>10} {_format_rel_err(gdna_err):>11}"
        )

    return "\n".join(lines)


# ── Transcript Abundance Accuracy ─────────────────────────────────────────

def analyze_abundance(sim_base: Path, conditions: list[str], truth: pd.DataFrame) -> str:
    """Analyze transcript-level abundance accuracy."""
    lines = []
    hr = "═" * 100
    condition_meta = condition_manifest_map(load_manifest(sim_base))

    lines.append(f"\n{hr}")
    lines.append("  TRANSCRIPT ABUNDANCE ACCURACY")
    lines.append(hr)

    # Per-condition summary
    lines.append(f"\n{'Condition':<35} {'Spearman':>8} {'Pearson':>8} "
                 f"{'MARD':>8} {'Med_RE':>8} "
                 f"{'n_FP':>5} {'n_FN':>5} "
                 f"{'gdna_leak':>9}")
    lines.append("─" * 100)

    from scipy.stats import spearmanr, pearsonr

    abundance_details = {}

    for cond in conditions:
        out_dir = sim_base / cond / "rigel_out"
        if not out_dir.exists():
            continue

        quant = load_quant(out_dir)
        cond_truth = load_condition_truth(sim_base, cond, condition_meta, truth)
        # Merge with truth
        merged = cond_truth.merge(
            quant[["transcript_id", "count", "count_em", "tpm"]],
            on="transcript_id", how="left"
        ).fillna(0)

        # Correlation on expressed transcripts only
        expressed = merged[merged["mrna_abundance"] > 0]
        unexpressed = merged[merged["mrna_abundance"] == 0]

        if len(expressed) > 2:
            sp_r, _ = spearmanr(expressed["mrna_abundance"], expressed["tpm"])
            pe_r, _ = pearsonr(
                np.log2(expressed["mrna_abundance"] + 1),
                np.log2(expressed["tpm"] + 1)
            )
        else:
            sp_r = pe_r = float("nan")

        # MARD (mean absolute relative difference) on expressed
        if len(expressed) > 0:
            rel_err = np.abs(expressed["tpm"] - expressed["mrna_abundance"]) / (
                expressed["mrna_abundance"] + 1e-6
            )
            mard = float(rel_err.mean())
            med_re = float(rel_err.median())
        else:
            mard = med_re = float("nan")

        # False positives: unexpressed txs with count > 5
        n_fp = int((unexpressed["count"] > 5).sum())
        # False negatives: expressed txs with count == 0
        n_fn = int((expressed["count"] == 0).sum())

        # gDNA leakage: total counts assigned to unexpressed transcripts
        gdna_leak = float(unexpressed["count"].sum())

        lines.append(
            f"{cond:<35} {sp_r:>8.4f} {pe_r:>8.4f} "
            f"{mard:>8.3f} {med_re:>8.3f} "
            f"{n_fp:>5} {n_fn:>5} "
            f"{gdna_leak:>9.1f}"
        )

        abundance_details[cond] = {
            "merged": merged,
            "spearman": sp_r,
            "pearson": pe_r,
            "mard": mard,
            "n_fp": n_fp,
            "n_fn": n_fn,
            "gdna_leak": gdna_leak,
        }

    # ── Per-gene breakdown for worst cases ──
    lines.append(f"\n\n{'─' * 100}")
    lines.append("  TOP FALSE POSITIVES (unexpressed txs with highest erroneous counts)")
    lines.append(f"{'─' * 100}")

    # Pick the highest-gDNA stranded condition for FP analysis
    worst_cond = "gdna_high_ss_0.99_nrna_none"
    if worst_cond in abundance_details:
        merged = abundance_details[worst_cond]["merged"]
        fp = merged[merged["mrna_abundance"] == 0].nlargest(15, "count")
        lines.append(f"\n  Condition: {worst_cond}")
        lines.append(f"  {'transcript_id':<20} {'gene_id':<12} {'n_exons':>7} "
                     f"{'spliced_len':>11} {'count':>8} {'tpm':>8}")
        for _, row in fp.iterrows():
            lines.append(
                f"  {row['transcript_id']:<20} {row['gene_id']:<12} "
                f"{int(row['n_exons']):>7} {int(row['spliced_length']):>11} "
                f"{row['count']:>8.1f} {row['tpm']:>8.2f}"
            )

    # ── Abundance accuracy by expression level ──
    lines.append(f"\n\n{'─' * 100}")
    lines.append("  ACCURACY BY EXPRESSION LEVEL (gdna_none_ss_0.99 — cleanest condition)")
    lines.append(f"{'─' * 100}")

    clean_cond = "gdna_none_ss_0.99_nrna_none"
    if clean_cond in abundance_details:
        merged = abundance_details[clean_cond]["merged"]
        expressed = merged[merged["mrna_abundance"] > 0].copy()
        expressed["log_abundance"] = np.log10(expressed["mrna_abundance"] + 1)
        bins = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
        expressed["bin"] = pd.cut(expressed["log_abundance"], bins=bins)

        lines.append(f"\n  {'log10(TPM) bin':<18} {'n_tx':>5} {'mean_RE':>8} {'med_RE':>8} "
                     f"{'mean_count':>10} {'mean_true':>10}")
        for b, grp in expressed.groupby("bin", observed=True):
            if len(grp) == 0:
                continue
            re = np.abs(grp["tpm"] - grp["mrna_abundance"]) / (grp["mrna_abundance"] + 1e-6)
            lines.append(
                f"  {str(b):<18} {len(grp):>5} {re.mean():>8.3f} {re.median():>8.3f} "
                f"{grp['count'].mean():>10.1f} {grp['mrna_abundance'].mean():>10.1f}"
            )

    return "\n".join(lines)


# ── Locus-level gDNA Analysis ─────────────────────────────────────────────

def analyze_locus_gdna(sim_base: Path, conditions: list[str]) -> str:
    """Analyze locus-level gDNA estimates."""
    lines = []
    hr = "═" * 100

    lines.append(f"\n{hr}")
    lines.append("  LOCUS-LEVEL gDNA ANALYSIS")
    lines.append(hr)

    lines.append(f"\n{'Condition':<35} {'n_loci':>6} {'tot_mrna':>9} "
                 f"{'tot_nrna':>9} {'tot_gdna':>9} {'tot_unambig':>11} "
                 f"{'gdna_rate':>9} {'max_gdna':>9}")
    lines.append("─" * 110)

    for cond in conditions:
        out_dir = sim_base / cond / "rigel_out"
        if not out_dir.exists():
            continue
        loci = load_loci(out_dir)
        if loci.empty:
            continue

        tot_mrna = loci["mrna"].sum()
        tot_nrna = loci["nrna"].sum() if "nrna" in loci.columns else 0
        tot_gdna = loci["gdna"].sum()
        tot_unambig = loci["count_unambig"].sum() if "count_unambig" in loci.columns else 0
        total = tot_mrna + tot_nrna + tot_gdna
        gdna_rate = tot_gdna / total if total > 0 else 0
        max_gdna = loci["gdna"].max()

        lines.append(
            f"{cond:<35} {len(loci):>6} {tot_mrna:>9.0f} "
            f"{tot_nrna:>9.0f} {tot_gdna:>9.0f} {tot_unambig:>11.0f} "
            f"{gdna_rate:>9.4f} {max_gdna:>9.1f}"
        )

    # Locus details for high-gDNA condition
    lines.append(f"\n\n{'─' * 100}")
    lines.append("  PER-LOCUS DETAIL (gdna_high_ss_0.99_nrna_none — top 20 loci by gDNA)")
    lines.append(f"{'─' * 100}")

    cond = "gdna_high_ss_0.99_nrna_none"
    out_dir = sim_base / cond / "rigel_out"
    if out_dir.exists():
        loci = load_loci(out_dir)
        if not loci.empty:
            top = loci.nlargest(20, "gdna")
            cols = ["locus_id", "n_transcripts", "locus_span_bp",
                    "count_unambig", "mrna", "nrna", "gdna", "gdna_rate", "gdna_prior"]
            avail_cols = [c for c in cols if c in top.columns]
            lines.append(f"\n  {top[avail_cols].to_string(index=False)}")

    return "\n".join(lines)


def collect_fragment_assignment_rows(sim_base: Path, conditions: list[str]) -> list[dict]:
    """Collect fragment-level assignment counts from annotated BAM tags."""
    import pysam

    overview_rows = []
    for cond in conditions:
        cond_dir = sim_base / cond
        annotated_bam = cond_dir / "annotated.bam"

        if not annotated_bam.exists():
            continue

        # Counters
        correct_tx = 0          # RNA → correct transcript
        correct_gene = 0        # RNA → wrong tx but correct gene
        wrong_tx = 0            # RNA → wrong transcript, wrong gene
        rna_as_gdna = 0         # RNA → classified as gDNA/intergenic
        gdna_as_rna = 0         # gDNA → assigned to a transcript
        gdna_correct = 0        # gDNA → classified as gDNA/intergenic
        total_rna = 0
        total_mrna = 0
        total_nrna = 0
        total_gdna = 0
        nrna_as_rna = 0
        nrna_as_gdna = 0
        mrna_as_gdna = 0
        mrna_as_nrna = 0

        with pysam.AlignmentFile(str(annotated_bam), "rb") as bam:
            for read in bam:
                if read.is_read2 or read.is_secondary or read.is_supplementary:
                    continue

                qname = read.query_name
                origin = parse_origin(qname)

                # Get assigned transcript from ZT tag
                try:
                    assigned_tx = read.get_tag("ZT")
                except KeyError:
                    assigned_tx = ""
                try:
                    assigned_gene = read.get_tag("ZG")
                except KeyError:
                    assigned_gene = ""
                try:
                    category = read.get_tag("ZC")
                except KeyError:
                    category = ""

                # Determine if rigel assigned to gDNA pool.
                # Use ZF bit 2 (AF_GDNA_BIT = 0x04). Rigel composes
                # AF_GDNA_EM = 0x05 and AF_GDNA_INTERGENIC = 0x25; both
                # have bit 2 set. ZL is the locus_id, NOT a gDNA flag.
                try:
                    zf = int(read.get_tag("ZF"))
                except KeyError:
                    zf = 0
                is_assigned_gdna = bool(zf & 0x04)
                is_assigned_nrna = bool(zf & 0x08)

                if origin.kind == "gdna":
                    total_gdna += 1
                    if is_assigned_gdna or category == "intergenic":
                        gdna_correct += 1
                    else:
                        gdna_as_rna += 1
                elif origin.kind == "nrna":
                    total_rna += 1
                    total_nrna += 1
                    if is_assigned_gdna or category == "intergenic":
                        rna_as_gdna += 1
                        nrna_as_gdna += 1
                    else:
                        nrna_as_rna += 1
                else:
                    total_rna += 1
                    total_mrna += 1
                    true_tx = origin.transcript_id or ""

                    if is_assigned_gdna or category == "intergenic":
                        rna_as_gdna += 1
                        mrna_as_gdna += 1
                    elif is_assigned_nrna:
                        mrna_as_nrna += 1
                    elif assigned_tx == true_tx:
                        correct_tx += 1
                    elif assigned_gene and true_tx.rsplit(".", 1)[0] == assigned_gene:
                        # Same gene, different isoform
                        correct_gene += 1
                    else:
                        wrong_tx += 1

        total = total_rna + total_gdna
        overview_rows.append({
            "condition": cond,
            "total": total,
            "total_rna": total_rna,
            "total_mrna": total_mrna,
            "total_nrna": total_nrna,
            "total_gdna": total_gdna,
            "correct_tx": correct_tx,
            "correct_gene": correct_gene,
            "wrong_tx": wrong_tx,
            "rna_as_gdna": rna_as_gdna,
            "mrna_as_gdna": mrna_as_gdna,
            "mrna_as_nrna": mrna_as_nrna,
            "nrna_as_rna": nrna_as_rna,
            "nrna_as_gdna": nrna_as_gdna,
            "gdna_as_rna": gdna_as_rna,
            "gdna_correct": gdna_correct,
        })

    return overview_rows


# ── Fragment-level Misallocation Analysis ─────────────────────────────────

def analyze_fragment_assignment(
    sim_base: Path,
    conditions: list[str],
    overview_rows: list[dict] | None = None,
) -> str:
    """Analyze fragment-level assignment using annotated BAM tags vs oracle truth.

    Oracle BAM read name encodes true origin:
      - RNA: {tx_id}:{start}-{end}:{strand}:{frag_num}
      - gDNA: gdna:{region}:{start}-{end}:{strand}:{frag_num}

    Annotated BAM tags:
      - ZT: assigned transcript_id
      - ZC: category (unambig, ambig_same_strand, ambig_opp_strand, multimapper, etc.)
      - ZW: posterior weight
      - ZG: assigned gene
    """
    lines = []
    hr = "═" * 100

    lines.append(f"\n{hr}")
    lines.append("  FRAGMENT-LEVEL ASSIGNMENT ANALYSIS")
    lines.append(hr)

    if overview_rows is None:
        overview_rows = collect_fragment_assignment_rows(sim_base, conditions)

    # Print overview table
    lines.append(f"\n{'Condition':<35} {'Total':>8} "
                 f"{'mRNA_ok':>8} {'mRNA_gene':>9} {'mRNA_bad':>9} "
                 f"{'nRNA→RNA':>9} {'nRNA→gDNA':>10} "
                 f"{'gDNA_ok':>8} {'gDNA→RNA':>9} {'Accuracy':>8}")
    lines.append("─" * 120)

    for row in overview_rows:
        total = row["total"]
        overall_correct = row["correct_tx"] + row["correct_gene"] + row["gdna_correct"]
        accuracy = overall_correct / total if total > 0 else 0

        lines.append(
            f"{row['condition']:<35} {total:>8,} "
            f"{row['correct_tx']:>8,} {row['correct_gene']:>9,} "
            f"{row['wrong_tx']:>9,} {row.get('nrna_as_rna', 0):>9,} "
            f"{row.get('nrna_as_gdna', 0):>10,} "
            f"{row['gdna_correct']:>8,} {row['gdna_as_rna']:>9,} "
            f"{accuracy:>8.4f}"
        )

    # Detailed breakdown for key conditions
    lines.append(f"\n\n{'─' * 100}")
    lines.append("  DETAILED BREAKDOWN")
    lines.append(f"{'─' * 100}")

    for row in overview_rows:
        if row["total_gdna"] == 0 and row["total_rna"] == 0:
            continue
        cond = row["condition"]
        lines.append(f"\n  {cond}:")

        total_mrna = row.get("total_mrna", row["total_rna"])
        if total_mrna > 0:
            rna_precision = row["correct_tx"] / total_mrna
            gene_precision = (row["correct_tx"] + row["correct_gene"]) / total_mrna
            mrna_as_gdna = row.get("mrna_as_gdna", row["rna_as_gdna"])
            rna_misclass = mrna_as_gdna / total_mrna
            lines.append(f"    mRNA ({total_mrna:,} frags):")
            lines.append(f"      Exact transcript: {rna_precision:.4f}")
            lines.append(f"      Correct gene:     {gene_precision:.4f}")
            lines.append(f"      Misclass as gDNA: {rna_misclass:.4f} ({mrna_as_gdna:,} frags)")
            lines.append(f"      Wrong gene:       {row['wrong_tx'] / total_mrna:.4f} ({row['wrong_tx']:,} frags)")

        total_nrna = row.get("total_nrna", 0)
        if total_nrna > 0:
            nrna_to_rna = row.get("nrna_as_rna", 0) / total_nrna
            nrna_to_gdna = row.get("nrna_as_gdna", 0) / total_nrna
            lines.append(f"    nRNA ({total_nrna:,} frags):")
            lines.append(f"      Routed to RNA-compatible tags:  {nrna_to_rna:.4f} ({row.get('nrna_as_rna', 0):,} frags)")
            lines.append(f"      Misclass as gDNA/intergenic:    {nrna_to_gdna:.4f} ({row.get('nrna_as_gdna', 0):,} frags)")

        if row["total_gdna"] > 0:
            gdna_precision = row["gdna_correct"] / row["total_gdna"]
            gdna_leak = row["gdna_as_rna"] / row["total_gdna"]
            lines.append(f"    gDNA ({row['total_gdna']:,} frags):")
            lines.append(f"      Correctly identified: {gdna_precision:.4f} ({row['gdna_correct']:,} frags)")
            lines.append(f"      Leaked to RNA:        {gdna_leak:.4f} ({row['gdna_as_rna']:,} frags)")

    return "\n".join(lines)


# ── Post-fix synthetic acceptance checks ──────────────────────────────────

def analyze_postfix_acceptance(
    sim_base: Path,
    conditions: list[str],
    assignment_rows: list[dict] | None = None,
    thresholds: CalibrationAcceptanceThresholds = CalibrationAcceptanceThresholds(),
) -> str:
    """Summarize the cheap regression checks that gate calibration readiness."""
    manifest = load_manifest(sim_base)
    condition_meta = condition_manifest_map(manifest)
    lines = []
    rows: list[dict] = []
    hr = "═" * 100

    def add_row(
        check: str,
        condition: str,
        value: str,
        threshold: str,
        passed: bool | None,
    ) -> None:
        if passed is None:
            status = "SKIP"
        else:
            status = "PASS" if passed else "FAIL"
        rows.append({
            "check": check,
            "condition": condition,
            "value": value,
            "threshold": threshold,
            "status": status,
        })

    # Boundary-density consistency: the implicit-splice bug mostly showed up as
    # depressed EXON|INTRON density relative to intergenic gDNA density.
    for cond in conditions:
        meta = condition_meta.get(cond, {})
        info = parse_condition(cond)
        n_mrna = int(meta.get("n_mrna", meta.get("n_rna", 1_000_000)))
        n_gdna = int(meta.get("n_gdna", round(n_mrna * info["gdna_frac"])))
        if n_gdna <= 0:
            continue

        summary = load_summary(sim_base / cond / "rigel_out")
        # v4: EXON-INTRON density family removed; the implicit-splice probe
        # used to compare exon-vs-intergenic, but the new density model only
        # fits INTERGENIC/INTRON anchors. Surface a stable "n/a" row.
        _ = summary
        add_row("rho_ex/rho_ig", cond, "n/a (v4)", ">= 0.950", None)

    # nRNA should stay essentially off in synthetic nrna_none conditions after
    # the implicit-splice false positives stop routing true gDNA into nRNA.
    for cond in conditions:
        meta = condition_meta.get(cond, {})
        nrna_label = str(meta.get("nrna_label", "none" if cond.endswith("_nrna_none") else ""))
        if nrna_label not in {"none", "zero"}:
            continue
        loci = load_loci(sim_base / cond / "rigel_out")
        if loci.empty or "nrna" not in loci.columns:
            add_row("nRNA in nrna_none", cond, "n/a", "<= 10", None)
            continue
        total_nrna = float(loci["nrna"].sum())
        add_row(
            "nRNA in nrna_none",
            cond,
            f"{total_nrna:.0f}",
            f"<= {thresholds.max_nrna_none_count:.0f}",
            total_nrna <= thresholds.max_nrna_none_count,
        )

    # Fragment-level gDNA leak requires annotated BAMs. Keep this bounded but do
    # not treat unavailable BAMs as a failed calibration run.
    assignment_rows = assignment_rows or []
    if assignment_rows:
        for row in assignment_rows:
            if row["total_gdna"] <= 0:
                continue
            leak_rate = row["gdna_as_rna"] / row["total_gdna"]
            add_row(
                "gDNA->RNA leak",
                row["condition"],
                f"{leak_rate:.4f}",
                f"<= {thresholds.max_gdna_to_rna_leak_rate:.4f}",
                leak_rate <= thresholds.max_gdna_to_rna_leak_rate,
            )
    else:
        add_row("gDNA->RNA leak", "all", "n/a", "annotated BAM required", None)

    lines.append(f"\n{hr}")
    lines.append("  POST-FIX CALIBRATION ACCEPTANCE CHECKS")
    lines.append(hr)
    lines.append(
        "  These checks pin the implicit-splice fix: boundary density should be "
        "coherent, nrna_none should stay near zero, and true gDNA should not "
        "leak substantially into RNA.\n"
    )
    lines.append(
        f"  {'Check':<20} {'Condition':<35} {'Value':>10} "
        f"{'Threshold':>14} {'Status':>7}"
    )
    lines.append("─" * 96)
    for row in rows:
        lines.append(
            f"  {row['check']:<20} {row['condition']:<35} {row['value']:>10} "
            f"{row['threshold']:>14} {row['status']:>7}"
        )

    failures = [row for row in rows if row["status"] == "FAIL"]
    passes = [row for row in rows if row["status"] == "PASS"]
    lines.append("")
    if failures:
        lines.append(f"  FAIL: {len(failures)} acceptance check(s) failed.")
    elif passes:
        lines.append(f"  PASS: all {len(passes)} evaluated acceptance checks passed.")
    else:
        lines.append("  SKIP: no acceptance checks could be evaluated from available outputs.")

    return "\n".join(lines)


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Run rigel analysis on synthetic simulation")
    parser.add_argument("--sim-base", type=Path, default=DEFAULT_SIM_BASE)
    parser.add_argument("--skip-quant", action="store_true",
                        help="Skip quantification, only run analysis")
    parser.add_argument("--skip-frag-analysis", action="store_true",
                        help="Skip fragment-level analysis (slow)")
    parser.add_argument("--conditions", nargs="*", default=None,
                        help="Optional subset of condition names to evaluate")
    args = parser.parse_args()

    sim_base = args.sim_base
    assert sim_base.exists(), f"Simulation directory not found: {sim_base}"
    conditions = discover_conditions(sim_base, args.conditions)

    print("=" * 100)
    print("  RIGEL SYNTHETIC SIMULATION ANALYSIS")
    print("=" * 100)

    # ── Step 1: Build index ──
    print(f"\n{'─' * 100}")
    print("  STEP 1: Building rigel index")
    print(f"{'─' * 100}")
    index_dir = build_index(sim_base)

    # ── Step 2: Run quant ──
    if not args.skip_quant:
        print(f"\n{'─' * 100}")
        print("  STEP 2: Running rigel quant on all conditions")
        print(f"{'─' * 100}")
        t0 = time.time()
        for cond in conditions:
            run_quant(sim_base, index_dir, cond)
        elapsed = time.time() - t0
        print(f"\n  All conditions quantified in {elapsed:.1f}s")
    else:
        print("\n  [SKIP] Quantification skipped (--skip-quant)")

    # ── Step 3: Analysis ──
    print(f"\n{'─' * 100}")
    print("  STEP 3: Comprehensive Analysis")
    print(f"{'─' * 100}")

    truth = load_truth(sim_base)

    # 3a: Calibration
    cal_report = analyze_calibration(sim_base, conditions, truth)
    print(cal_report)

    # 3b: Abundance accuracy
    abundance_report = analyze_abundance(sim_base, conditions, truth)
    print(abundance_report)

    # 3c: Locus-level gDNA
    locus_report = analyze_locus_gdna(sim_base, conditions)
    print(locus_report)

    # 3d: Fragment-level (optional, slow)
    assignment_rows = []
    if not args.skip_frag_analysis:
        assignment_rows = collect_fragment_assignment_rows(sim_base, conditions)
        frag_report = analyze_fragment_assignment(sim_base, conditions, assignment_rows)
        print(frag_report)
    else:
        frag_report = "\n  [SKIP] Fragment analysis skipped (--skip-frag-analysis)"
        print(frag_report)

    # 3e: Cheap post-fix calibration acceptance checks
    acceptance_report = analyze_postfix_acceptance(sim_base, conditions, assignment_rows)
    print(acceptance_report)

    # ── Save full report ──
    report_path = sim_base / "analysis_report.txt"
    full_report = "\n".join([
        "=" * 100,
        "  RIGEL SYNTHETIC SIMULATION — FULL ANALYSIS REPORT",
        "=" * 100,
        cal_report,
        abundance_report,
        locus_report,
        frag_report,
        acceptance_report,
    ])
    report_path.write_text(full_report)
    print(f"\n\n  Full report saved to: {report_path}")
    print("=" * 100)


if __name__ == "__main__":
    main()
