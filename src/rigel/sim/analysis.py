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
from .read_name import parse_origin

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


def _load_condition_truth_summary(sim_base: Path, meta: dict) -> dict:
    summary_name = meta.get("truth_summary")
    if not summary_name:
        return {}
    summary_path = Path(str(summary_name))
    if not summary_path.is_absolute():
        summary_path = sim_base / summary_path
    if not summary_path.exists():
        return {}
    with open(summary_path) as handle:
        return json.load(handle)


def _truth_fl_mean(summary: dict, kind: str, fallback: float) -> float:
    fragment_lengths = summary.get("fragment_lengths", {}) if summary else {}
    stats = fragment_lengths.get(kind, {}) if isinstance(fragment_lengths, dict) else {}
    mean = stats.get("mean") if isinstance(stats, dict) else None
    return float(mean) if mean is not None else fallback


def _format_bp(value: float) -> str:
    """Compact bp value for report prose."""
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def _format_rel_err(value: float) -> str:
    if np.isnan(value):
        return "n/a"
    return f"{value:+.3f}"


def _as_float(value: object, default: float = float("nan")) -> float:
    try:
        if value is None:
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0 or not np.isfinite(denominator):
        return float("nan")
    return numerator / denominator


def _relative_error(estimated: float, truth: float) -> float:
    return _safe_div(estimated - truth, truth)


def _nested(data: dict | None, *keys: str, default: object = None) -> object:
    current: object = data or {}
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def _stat(summary: object, key: str, default: float = float("nan")) -> float:
    if not isinstance(summary, dict):
        return default
    return _as_float(summary.get(key), default)


def _fmt_float(value: object, digits: int = 3) -> str:
    number = _as_float(value)
    if not np.isfinite(number):
        return "n/a"
    return f"{number:.{digits}f}"


def _fmt_count(value: object) -> str:
    number = _as_float(value)
    if not np.isfinite(number):
        return "n/a"
    return f"{number:,.0f}"


def _fmt_delta(value: object) -> str:
    number = _as_float(value)
    if not np.isfinite(number):
        return "n/a"
    return f"{number:+,.0f}"


def _fmt_pct(value: object, digits: int = 2) -> str:
    number = _as_float(value)
    if not np.isfinite(number):
        return "n/a"
    return f"{100.0 * number:.{digits}f}%"


def _fmt_bool(value: object) -> str:
    if value is True:
        return "yes"
    if value is False:
        return "no"
    return "n/a"


def _fmt_text(value: object, default: str = "none") -> str:
    if value is None:
        return default
    if isinstance(value, float) and np.isnan(value):
        return default
    text = str(value)
    return default if text.lower() in {"", "nan", "none"} else text


def _markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    if not rows:
        return "_No rows available._"
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def _condition_truth_fl_mean(
    truth_summary: dict,
    kind: str,
    fallback: float,
    *,
    has_fragments: bool,
) -> float:
    stats = (truth_summary.get("fragment_lengths", {}) or {}).get(kind, {})
    if isinstance(stats, dict) and stats.get("mean") is not None:
        return float(stats["mean"])
    return fallback if has_fragments else float("nan")


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

    # Table: condition → acyclic-calibrator scalars (one global gDNA density + strand model;
    # the old per-region-type ρ_ig/ρ_in/ρ_ex density family was removed in the acyclic redesign).
    lines.append(f"\n{'Condition':<35} {'gDNA_frac':>9} {'SS':>5} "
                 f"{'dens_glob':>10} {'rna_sense':>9} {'od_gdna':>8} {'od_rna':>8} "
                 f"{'FL_rna':>7} {'FL_gdna':>8} "
                 f"{'n_loci':>6} {'gdna_rate':>10} {'gdna_true':>10}")
    lines.append("─" * 140)

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

        cal = summary.get("calibration", {}) or {}
        density_global = _as_float(cal.get("gdna_density_global"), 0.0)
        rna_sense_frac = _as_float(cal.get("rna_sense_frac"), float("nan"))
        od_gdna = _as_float(cal.get("gdna_strand_overdispersion"), float("nan"))
        od_rna = _as_float(cal.get("rna_strand_overdispersion"), float("nan"))

        fl = summary.get("fragment_length", {}) or {}
        fl_rna = _as_float(_nested(fl, "rna", "summary", "mean", default=0.0), 0.0)
        fl_gdna = _as_float(_nested(fl, "gdna", "summary", "mean", default=0.0), 0.0)
        quant_out = summary.get("quantification", {}) or {}
        n_loci = int(quant_out.get("n_loci", 0) or 0)

        # Compute true gDNA rate
        gdna_frac = info["gdna_frac"]
        n_rna = int(meta.get("n_rna", 1_000_000))
        n_gdna_true = int(meta.get("n_gdna", round(n_rna * gdna_frac)))
        n_total_true = n_rna + n_gdna_true
        gdna_rate_true = n_gdna_true / n_total_true if n_total_true > 0 else 0
        truth_summary = _load_condition_truth_summary(sim_base, meta)
        truth_rna_fl_mean = _truth_fl_mean(truth_summary, "mrna", true_rna_fl_mean)
        truth_gdna_fl_mean = _truth_fl_mean(truth_summary, "gdna", true_gdna_fl_mean)

        gdna_rate_est = quant_out.get("gdna_fraction", 0)

        lines.append(
            f"{cond:<35} {gdna_frac:>9.2f} {info['strand_specificity']:>5.2f} "
            f"{density_global:>10.6f} {rna_sense_frac:>9.3f} {od_gdna:>8.3f} {od_rna:>8.3f} "
            f"{fl_rna:>7.1f} {fl_gdna:>8.1f} "
            f"{n_loci:>6} {gdna_rate_est:>10.4f} {gdna_rate_true:>10.4f}"
        )

        cal_rows.append({
            **info,
            "density_global": density_global,
            "rna_sense_frac": rna_sense_frac,
            "od_gdna": od_gdna, "od_rna": od_rna,
            "fl_rna": fl_rna, "fl_gdna": fl_gdna,
            "n_loci": n_loci,
            "gdna_rate_est": gdna_rate_est,
            "gdna_rate_true": gdna_rate_true,
            "n_rna_true": n_rna,
            "n_gdna_true": n_gdna_true,
            "truth_rna_fl_mean": truth_rna_fl_mean,
            "truth_gdna_fl_mean": truth_gdna_fl_mean,
            "intergenic_total": quant_out.get("intergenic_total", 0),
            "gdna_total": quant_out.get("gdna_total", 0),
            "mrna_total": quant_out.get("mrna_total", 0),
        })

    # ── gDNA FL distribution ──
    lines.append(f"\n\n{'─' * 100}")
    lines.append("  FRAGMENT LENGTH MODELS")
    lines.append(f"{'─' * 100}")
    lines.append(
        "  Default simulated FL: "
        f"RNA={_format_bp(true_rna_fl_mean)}bp (mean), "
        f"gDNA={_format_bp(true_gdna_fl_mean)}bp (mean). "
        "Condition truth_summary files override these means when available.\n"
    )
    lines.append(f"  {'Condition':<35} {'FL_rna':>8} {'FL_gdna':>8} {'FL_rna_err':>10} {'FL_gdna_err':>11}")
    for row in cal_rows:
        rna_err = (
            (row["fl_rna"] - row["truth_rna_fl_mean"]) / row["truth_rna_fl_mean"]
            if row["fl_rna"] > 0 and row["truth_rna_fl_mean"] > 0
            else float("nan")
        )
        gdna_err = (
            (row["fl_gdna"] - row["truth_gdna_fl_mean"]) / row["truth_gdna_fl_mean"]
            if (
                row["n_gdna_true"] > 0
                and row["fl_gdna"] > 0
                and row["truth_gdna_fl_mean"] > 0
            )
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
                    "count_unambig", "mrna", "nrna", "gdna", "gdna_rate", "gdna_prior_count"]
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


# ── Detailed per-condition report ────────────────────────────────────────

def _abundance_diagnostics(
    sim_base: Path,
    condition: str,
    condition_meta: dict[str, dict],
    fallback_truth: pd.DataFrame,
) -> tuple[dict[str, float], list[dict[str, object]]]:
    out_dir = sim_base / condition / "rigel_out"
    if not out_dir.exists():
        return {}, []
    try:
        quant = load_quant(out_dir)
    except FileNotFoundError:
        return {}, []
    cond_truth = load_condition_truth(sim_base, condition, condition_meta, fallback_truth)
    if "transcript_id" not in cond_truth.columns or "transcript_id" not in quant.columns:
        return {}, []

    quant_cols = [col for col in ("transcript_id", "count", "count_em", "tpm") if col in quant]
    if "count" not in quant_cols:
        return {}, []
    merged = cond_truth.merge(quant[quant_cols], on="transcript_id", how="left").fillna(0)
    true = merged["mrna_abundance"].astype(float)
    estimated = merged["count"].astype(float)
    expressed = true > 0
    unexpressed = ~expressed
    error = estimated - true
    abs_error = error.abs()

    metrics = {
        "n_transcripts": float(len(merged)),
        "n_expressed": float(expressed.sum()),
        "n_unexpressed": float(unexpressed.sum()),
        "mRNA_count_true": float(true.sum()),
        "mRNA_count_est": float(estimated.sum()),
        "mRNA_count_abs_error_total": float(abs_error.sum()),
        "mRNA_count_mae": float(abs_error.mean()) if len(abs_error) else float("nan"),
        "false_positive_transcripts_gt5": float((estimated[unexpressed] > 5).sum()),
        "false_positive_count_total": float(estimated[unexpressed].sum()),
        "false_negative_transcripts_lt1": float((estimated[expressed] < 1).sum()),
    }
    if expressed.any():
        rel_error = abs_error[expressed] / (true[expressed] + 1e-9)
        metrics["mRNA_count_mard"] = float(rel_error.mean())
        metrics["mRNA_count_medard"] = float(rel_error.median())
        if true[expressed].nunique() > 1 and estimated[expressed].nunique() > 1:
            metrics["mRNA_count_spearman"] = float(
                true[expressed].corr(estimated[expressed], method="spearman")
            )
        else:
            metrics["mRNA_count_spearman"] = float("nan")
    else:
        metrics["mRNA_count_mard"] = float("nan")
        metrics["mRNA_count_medard"] = float("nan")
        metrics["mRNA_count_spearman"] = float("nan")

    detail_cols = [
        col
        for col in ("transcript_id", "gene_id", "mrna_abundance", "count", "tpm")
        if col in merged
    ]
    merged = merged.assign(abs_error=abs_error, signed_error=error)
    top_errors = (
        merged.nlargest(5, "abs_error")[[*detail_cols, "signed_error", "abs_error"]]
        .to_dict("records")
    )
    return metrics, top_errors


def _pool_rows(row: dict[str, object]) -> list[list[str]]:
    pools = [
        ("mRNA", row["true_mrna"], row["est_mrna"]),
        ("nRNA", row["true_nrna"], row["est_nrna"]),
        ("RNA total", row["true_rna"], row["est_rna"]),
        ("gDNA genic EM", float("nan"), row["est_gdna_genic"]),
        ("gDNA intergenic", float("nan"), row["est_intergenic"]),
        ("gDNA total", row["true_gdna"], row["est_gdna_total"]),
        ("All fragments", row["true_total"], row["est_total"]),
    ]
    rows = []
    for label, truth, estimated in pools:
        truth_value = _as_float(truth)
        estimated_value = _as_float(estimated)
        delta = estimated_value - truth_value if np.isfinite(truth_value) else float("nan")
        rel = _relative_error(estimated_value, truth_value)
        rows.append([
            label,
            _fmt_count(truth_value),
            _fmt_count(estimated_value),
            _fmt_delta(delta),
            _fmt_pct(rel) if np.isfinite(rel) else "n/a",
        ])
    return rows


def _collect_condition_metrics(
    sim_base: Path,
    conditions: list[str],
    assignment_rows: list[dict] | None,
) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    manifest = load_manifest(sim_base)
    condition_meta = condition_manifest_map(manifest)
    fallback_truth = load_truth(sim_base)
    default_rna_fl_mean, default_gdna_fl_mean = simulated_fragment_length_means(manifest)
    assignment_by_condition = {
        row["condition"]: row for row in (assignment_rows or []) if "condition" in row
    }

    rows: list[dict[str, object]] = []
    details: dict[str, dict[str, object]] = {}
    for condition in conditions:
        meta = condition_meta.get(condition, {})
        info = parse_condition(condition)
        gdna_label = str(meta.get("gdna_label", info["gdna_label"]))
        capture_label = str(meta.get("capture_label", "none"))
        strand_specificity = float(meta.get("strand_specificity", info["strand_specificity"]))
        capture_config = meta.get("capture_config", {}) or {}
        capture_probe_panel = meta.get("capture_probe_panel")
        if not capture_probe_panel:
            capture_probe_panel = (
                meta.get("capture_probe_bed")
                or meta.get("capture_probe_tsv")
                or capture_config.get("probes")
                if isinstance(capture_config, dict)
                else None
            )
        capture_probe_source = meta.get("capture_probe_source")
        if capture_probe_source is None and capture_probe_panel:
            capture_probe_source = (
                "generated"
                if meta.get("capture_probe_bed") or meta.get("capture_probe_tsv")
                else "provided"
            )

        out_dir = sim_base / condition / "rigel_out"
        summary = load_summary(out_dir)
        quant = summary.get("quantification", {}) if summary else {}
        cal = summary.get("calibration", {}) if summary else {}
        region = cal.get("region_calibration", {}) if isinstance(cal, dict) else {}
        diagnostics = cal.get("diagnostics", {}) if isinstance(cal, dict) else {}
        # Estimated FL models moved to the top-level "fragment_length" block (acyclic schema).
        fl_models = summary.get("fragment_length", {}) if summary else {}
        background = cal.get("background_model", {}) if isinstance(cal, dict) else {}
        boundary_local = cal.get("boundary_local", {}) if isinstance(cal, dict) else {}
        boundary_sweep = cal.get("boundary_sweep", {}) if isinstance(cal, dict) else {}
        strand_channels = cal.get("strand_channels") or {}
        prior_table = cal.get("prior_table") or summary.get("prior_policy", {}) or {}

        true_mrna = float(meta.get("n_mrna", meta.get("n_rna", 0)))
        true_nrna = float(meta.get("n_nrna", 0))
        true_rna = float(meta.get("n_rna", true_mrna + true_nrna))
        true_gdna = float(meta.get("n_gdna", 0))
        true_total = float(meta.get("n_total", true_rna + true_gdna))

        est_mrna = _as_float(quant.get("mrna_total"))
        est_nrna = _as_float(quant.get("nrna_total"))
        est_gdna_genic = _as_float(quant.get("gdna_total"))
        est_intergenic = _as_float(quant.get("intergenic_total"), 0.0)
        est_rna = est_mrna + est_nrna
        est_gdna_total = est_gdna_genic + est_intergenic
        est_total = est_rna + est_gdna_total

        truth_summary = _load_condition_truth_summary(sim_base, meta)
        truth_rna_fl_mean = _condition_truth_fl_mean(
            truth_summary,
            "mrna",
            default_rna_fl_mean,
            has_fragments=true_mrna > 0,
        )
        truth_gdna_fl_mean = _condition_truth_fl_mean(
            truth_summary,
            "gdna",
            default_gdna_fl_mean,
            has_fragments=true_gdna > 0,
        )

        assignment = assignment_by_condition.get(condition, {})
        total_mrna_assign = float(assignment.get("total_mrna", 0))
        total_gdna_assign = float(assignment.get("total_gdna", 0))
        gdna_as_rna = float(assignment.get("gdna_as_rna", 0))
        mrna_as_gdna = float(assignment.get("mrna_as_gdna", 0))
        wrong_tx = float(assignment.get("wrong_tx", 0))
        correct_tx = float(assignment.get("correct_tx", 0))
        correct_gene = float(assignment.get("correct_gene", 0))

        abundance_metrics, top_errors = _abundance_diagnostics(
            sim_base,
            condition,
            condition_meta,
            fallback_truth,
        )

        row = {
            "condition": condition,
            "gdna_label": gdna_label,
            "strand_specificity": strand_specificity,
            "capture_label": capture_label,
            "capture_enabled": bool(meta.get("capture_enabled", False)),
            "capture_probe_source": capture_probe_source,
            "capture_probe_panel": capture_probe_panel,
            "true_mrna": true_mrna,
            "true_nrna": true_nrna,
            "true_rna": true_rna,
            "true_gdna": true_gdna,
            "true_total": true_total,
            "true_gdna_fraction": _safe_div(true_gdna, true_total),
            "est_mrna": est_mrna,
            "est_nrna": est_nrna,
            "est_rna": est_rna,
            "est_gdna_genic": est_gdna_genic,
            "est_intergenic": est_intergenic,
            "est_gdna_total": est_gdna_total,
            "est_total": est_total,
            "est_gdna_fraction": _safe_div(est_gdna_total, est_total),
            "rna_delta": est_rna - true_rna,
            "gdna_delta": est_gdna_total - true_gdna,
            "gdna_fraction_delta": _safe_div(est_gdna_total, est_total)
            - _safe_div(true_gdna, true_total),
            "truth_rna_fl_mean": truth_rna_fl_mean,
            "truth_gdna_fl_mean": truth_gdna_fl_mean,
            "est_rna_fl_mean": _as_float(_nested(fl_models, "rna", "summary", "mean")),
            "est_gdna_fl_mean": _as_float(_nested(fl_models, "gdna", "summary", "mean")),
            "calibration_n_observed": _as_float(diagnostics.get("n_observed")),
            "calibration_n_excluded_multimap": _as_float(
                diagnostics.get("n_excluded_multimap")
            ),
            "rho_off": _as_float(region.get("rho_off")),
            "kappa_d": _as_float(region.get("kappa_d")),
            "region_n_passes": _as_float(region.get("n_passes")),
            "region_converged": region.get("converged"),
            "region_n_regions": _as_float(region.get("n_regions")),
            "p_unexpressed_p50": _stat(region.get("p_unexpressed"), "p50"),
            "p_unexpressed_p95": _stat(region.get("p_unexpressed"), "p95"),
            "p_expressed_p50": _stat(region.get("p_expressed"), "p50"),
            "p_expressed_p95": _stat(region.get("p_expressed"), "p95"),
            "mu_gdna_p50": _stat(region.get("mu_gdna"), "p50"),
            "mu_gdna_p95": _stat(region.get("mu_gdna"), "p95"),
            "upper_gdna_p95": _stat(region.get("upper_gdna"), "p95"),
            "background_rho_off_mean": _as_float(background.get("rho_off_mean")),
            "background_fit_status": background.get("fit_status"),
            "background_n_seed_regions": _as_float(background.get("n_seed_regions")),
            "boundary_local_regions_with_evidence": _as_float(
                boundary_local.get("n_regions_with_evidence")
            ),
            "boundary_sweep_regions_with_evidence": _as_float(
                boundary_sweep.get("n_regions_with_swept_evidence")
            ),
            "strand_model_est": _as_float(_nested(summary, "strand_model", "strand_specificity")),
            "strand_protocol": _nested(summary, "strand_model", "protocol"),
            "strand_channel_kappa_d": _as_float(strand_channels.get("kappa_d")),
            "fragment_total": _as_float(_nested(summary, "fragment_stats", "total")),
            "fragment_intergenic": _as_float(_nested(summary, "fragment_stats", "intergenic")),
            "fragment_chimeric": _as_float(_nested(summary, "fragment_stats", "chimeric")),
            "n_loci": _as_float(quant.get("n_loci")),
            "gdna_to_rna_rate": _safe_div(gdna_as_rna, total_gdna_assign),
            "mrna_to_gdna_rate": _safe_div(mrna_as_gdna, total_mrna_assign),
            "mrna_exact_rate": _safe_div(correct_tx, total_mrna_assign),
            "mrna_gene_rate": _safe_div(correct_tx + correct_gene, total_mrna_assign),
            "mrna_wrong_gene_rate": _safe_div(wrong_tx, total_mrna_assign),
            "gdna_as_rna": gdna_as_rna,
            "mrna_as_gdna": mrna_as_gdna,
            "prior_global_gdna": _as_float(_nested(prior_table, "global_counts", "gdna")),
            "prior_global_rna": _as_float(_nested(prior_table, "global_counts", "rna")),
            **abundance_metrics,
        }
        rows.append(row)
        details[condition] = {
            "summary": summary,
            "meta": meta,
            "truth_summary": truth_summary,
            "assignment": assignment,
            "top_errors": top_errors,
        }

    return pd.DataFrame(rows), details


def build_condition_report(
    sim_base: Path,
    conditions: list[str],
    assignment_rows: list[dict] | None = None,
) -> tuple[str, pd.DataFrame]:
    metrics, details = _collect_condition_metrics(sim_base, conditions, assignment_rows)
    lines: list[str] = []
    lines.append("# Rigel Synthetic Capture Suite Condition Report")
    lines.append("")
    lines.append(f"Base folder: `{sim_base}`")
    lines.append(f"Conditions evaluated: `{len(conditions)}`")
    lines.append("")
    lines.append(
        "Primary criterion here is pool-level deconvolution: total estimated gDNA "
        "is Rigel genic gDNA plus intergenic fragments, compared with the simulator "
        "post-capture truth for each condition."
    )

    overview_rows = []
    for _, row in metrics.iterrows():
        overview_rows.append([
            row["condition"],
            row["gdna_label"],
            _fmt_float(row["strand_specificity"], 2),
            row["capture_label"],
            _fmt_pct(row["true_gdna_fraction"]),
            _fmt_pct(row["est_gdna_fraction"]),
            _fmt_delta(row["gdna_delta"]),
            _fmt_delta(row["rna_delta"]),
            _fmt_pct(row["gdna_to_rna_rate"]),
            _fmt_pct(row["mrna_to_gdna_rate"]),
            _fmt_bool(row["region_converged"]),
        ])
    lines.append("\n## Pool Deconvolution Overview")
    lines.append(_markdown_table([
        "Condition",
        "gDNA",
        "SS",
        "Capture",
        "True gDNA",
        "Rigel gDNA",
        "gDNA delta",
        "RNA delta",
        "gDNA->RNA",
        "mRNA->gDNA",
        "Cal converged",
    ], overview_rows))

    calibration_rows = []
    for _, row in metrics.iterrows():
        rna_fl_err = _relative_error(row["est_rna_fl_mean"], row["truth_rna_fl_mean"])
        gdna_fl_err = _relative_error(row["est_gdna_fl_mean"], row["truth_gdna_fl_mean"])
        calibration_rows.append([
            row["condition"],
            _fmt_float(row["rho_off"], 6),
            _fmt_float(row["kappa_d"], 3),
            _fmt_float(row["background_rho_off_mean"], 6),
            str(row["background_fit_status"]),
            _fmt_count(row["background_n_seed_regions"]),
            _fmt_float(row["est_rna_fl_mean"], 1),
            _fmt_pct(rna_fl_err),
            _fmt_float(row["est_gdna_fl_mean"], 1),
            _fmt_pct(gdna_fl_err),
            _fmt_float(row["p_unexpressed_p95"], 3),
            _fmt_float(row["p_expressed_p95"], 3),
        ])
    lines.append("\n## Calibration Diagnostics Overview")
    lines.append(_markdown_table([
        "Condition",
        "rho_off",
        "kappa_d",
        "bg rho",
        "bg fit",
        "seed regions",
        "RNA FL",
        "RNA FL err",
        "gDNA FL",
        "gDNA FL err",
        "p_unexpressed p95",
        "p_expressed p95",
    ], calibration_rows))

    for _, row in metrics.iterrows():
        condition = str(row["condition"])
        detail = details.get(condition, {})
        meta = detail.get("meta", {})
        assignment = detail.get("assignment", {})
        summary = detail.get("summary", {})
        cal = summary.get("calibration", {}) if isinstance(summary, dict) else {}
        diagnostics = cal.get("diagnostics", {}) if isinstance(cal, dict) else {}
        coarse = diagnostics.get("mass_by_coarse_class", {}) if isinstance(diagnostics, dict) else {}
        fl_pool = diagnostics.get("fl_pool_total", {}) if isinstance(diagnostics, dict) else {}
        state_mass = _nested(cal, "region_calibration", "state_mass", default={})

        lines.append(f"\n## {condition}")
        lines.append("")
        lines.append(_markdown_table(["Field", "Value"], [
            ["gDNA label", row["gdna_label"]],
            ["strand specificity truth", _fmt_float(row["strand_specificity"], 2)],
            ["strand specificity estimated", _fmt_float(row["strand_model_est"], 3)],
            ["strand protocol", row["strand_protocol"]],
            ["capture label", row["capture_label"]],
            ["capture enabled", _fmt_bool(row["capture_enabled"])],
            ["capture source", _fmt_text(row["capture_probe_source"])],
            ["capture panel", _fmt_text(row["capture_probe_panel"])],
            ["Rigel output", str(sim_base / condition / "rigel_out")],
        ]))

        lines.append("\n### Pool Deconvolution")
        lines.append(_markdown_table(
            ["Pool", "Truth", "Rigel estimate", "Delta", "Relative delta"],
            _pool_rows(row.to_dict()),
        ))
        lines.append("")
        lines.append(
            f"True gDNA fraction: `{_fmt_pct(row['true_gdna_fraction'])}`; "
            f"Rigel gDNA fraction: `{_fmt_pct(row['est_gdna_fraction'])}`. "
            f"Genic EM gDNA: `{_fmt_count(row['est_gdna_genic'])}`; "
            f"intergenic gDNA: `{_fmt_count(row['est_intergenic'])}`."
        )

        lines.append("\n### Calibration")
        lines.append(_markdown_table(["Metric", "Value"], [
            ["calibration observed fragments", _fmt_count(row["calibration_n_observed"])],
            ["excluded multimappers", _fmt_count(row["calibration_n_excluded_multimap"])],
            ["region count", _fmt_count(row["region_n_regions"])],
            ["rho_off", _fmt_float(row["rho_off"], 8)],
            ["kappa_d", _fmt_float(row["kappa_d"], 3)],
            ["passes", _fmt_count(row["region_n_passes"])],
            ["converged", _fmt_bool(row["region_converged"])],
            ["p_unexpressed p50 / p95", f"{_fmt_float(row['p_unexpressed_p50'])} / {_fmt_float(row['p_unexpressed_p95'])}"],
            ["p_expressed p50 / p95", f"{_fmt_float(row['p_expressed_p50'])} / {_fmt_float(row['p_expressed_p95'])}"],
            ["mu_gdna p50 / p95", f"{_fmt_float(row['mu_gdna_p50'])} / {_fmt_float(row['mu_gdna_p95'])}"],
            ["upper_gdna p95", _fmt_float(row["upper_gdna_p95"])],
            ["background fit", row["background_fit_status"]],
            ["background rho mean", _fmt_float(row["background_rho_off_mean"], 8)],
            ["background seed regions", _fmt_count(row["background_n_seed_regions"])],
            ["boundary local evidence regions", _fmt_count(row["boundary_local_regions_with_evidence"])],
            ["boundary sweep evidence regions", _fmt_count(row["boundary_sweep_regions_with_evidence"])],
            ["prior global gDNA / RNA", f"{_fmt_count(row['prior_global_gdna'])} / {_fmt_count(row['prior_global_rna'])}"],
        ]))

        lines.append("\n### Fragment Length Calibration")
        lines.append(_markdown_table(["Pool", "Truth mean", "Rigel mean", "Relative error"], [
            [
                "RNA",
                _fmt_float(row["truth_rna_fl_mean"], 1),
                _fmt_float(row["est_rna_fl_mean"], 1),
                _fmt_pct(_relative_error(row["est_rna_fl_mean"], row["truth_rna_fl_mean"])),
            ],
            [
                "gDNA",
                _fmt_float(row["truth_gdna_fl_mean"], 1),
                _fmt_float(row["est_gdna_fl_mean"], 1),
                _fmt_pct(_relative_error(row["est_gdna_fl_mean"], row["truth_gdna_fl_mean"])),
            ],
        ]))

        lines.append("\n### Calibration Payload Mix")
        lines.append(_markdown_table(["Payload field", "Value"], [
            ["coarse INTERGENIC mass", _fmt_count(_as_float(coarse.get("INTERGENIC")))],
            ["coarse INTRON mass", _fmt_count(_as_float(coarse.get("INTRON")))],
            ["coarse EXON mass", _fmt_count(_as_float(coarse.get("EXON")))],
            ["FL pool global", _fmt_count(_as_float(fl_pool.get("GLOBAL")))],
            ["FL pool RNA", _fmt_count(_as_float(fl_pool.get("RNA")))],
            ["FL pool gDNA", _fmt_count(_as_float(fl_pool.get("GDNA")))],
        ]))

        if isinstance(state_mass, dict) and state_mass:
            state_rows = []
            for name, values in state_mass.items():
                if isinstance(values, dict):
                    state_rows.append([
                        name,
                        _fmt_float(values.get("sum"), 2),
                        _fmt_float(values.get("mean"), 4),
                    ])
            lines.append("\n### Region Latent State Mass")
            lines.append(_markdown_table(["State", "Sum", "Mean"], state_rows))

        lines.append("\n### Fragment Assignment Diagnostics")
        if assignment:
            lines.append(_markdown_table(["Metric", "Value"], [
                ["mRNA exact transcript", _fmt_pct(row["mrna_exact_rate"])],
                ["mRNA correct gene", _fmt_pct(row["mrna_gene_rate"])],
                ["mRNA wrong gene", _fmt_pct(row["mrna_wrong_gene_rate"])],
                ["mRNA routed to gDNA", f"{_fmt_pct(row['mrna_to_gdna_rate'])} ({_fmt_count(row['mrna_as_gdna'])})"],
                ["gDNA routed to RNA", f"{_fmt_pct(row['gdna_to_rna_rate'])} ({_fmt_count(row['gdna_as_rna'])})"],
                ["total mRNA fragments checked", _fmt_count(assignment.get("total_mrna"))],
                ["total gDNA fragments checked", _fmt_count(assignment.get("total_gdna"))],
            ]))
        else:
            lines.append("_Annotated BAM fragment diagnostics were not available._")

        lines.append("\n### Transcript Count Accuracy")
        lines.append(_markdown_table(["Metric", "Value"], [
            ["expressed transcripts", _fmt_count(row.get("n_expressed"))],
            ["unexpressed transcripts", _fmt_count(row.get("n_unexpressed"))],
            ["mRNA count Spearman", _fmt_float(row.get("mRNA_count_spearman"), 4)],
            ["mRNA count MARD", _fmt_pct(row.get("mRNA_count_mard"))],
            ["mRNA count median ARD", _fmt_pct(row.get("mRNA_count_medard"))],
            ["total absolute count error", _fmt_count(row.get("mRNA_count_abs_error_total"))],
            ["false-positive transcripts >5 frags", _fmt_count(row.get("false_positive_transcripts_gt5"))],
            ["false-positive fragment total", _fmt_count(row.get("false_positive_count_total"))],
            ["false-negative expressed tx <1 frag", _fmt_count(row.get("false_negative_transcripts_lt1"))],
        ]))

        top_errors = detail.get("top_errors", [])
        if top_errors:
            top_rows = []
            for item in top_errors:
                top_rows.append([
                    item.get("transcript_id", ""),
                    item.get("gene_id", ""),
                    _fmt_count(item.get("mrna_abundance")),
                    _fmt_count(item.get("count")),
                    _fmt_delta(item.get("signed_error")),
                    _fmt_count(item.get("abs_error")),
                ])
            lines.append("\nTop transcript count errors:")
            lines.append(_markdown_table([
                "Transcript",
                "Gene",
                "Truth",
                "Rigel",
                "Delta",
                "Abs error",
            ], top_rows))

        lines.append("\n### Run Files")
        lines.append(_markdown_table(["Artifact", "Path"], [
            ["summary", str(sim_base / condition / "rigel_out" / "summary.json")],
            ["quant", str(sim_base / condition / "rigel_out" / "quant.feather")],
            ["loci", str(sim_base / condition / "rigel_out" / "loci.feather")],
            ["truth", str(sim_base / str(meta.get("truth_abundances", "")))],
            ["truth summary", str(sim_base / str(meta.get("truth_summary", "")))],
        ]))

    return "\n".join(lines) + "\n", metrics


def write_condition_report(
    sim_base: Path,
    conditions: list[str],
    assignment_rows: list[dict] | None = None,
) -> tuple[Path, Path]:
    report, metrics = build_condition_report(sim_base, conditions, assignment_rows)
    report_path = sim_base / "condition_report.md"
    metrics_path = sim_base / "condition_metrics.tsv"
    report_path.write_text(report)
    metrics.to_csv(metrics_path, sep="\t", index=False)
    return report_path, metrics_path


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

    # 3f: Detailed per-condition Markdown report
    condition_report_path, condition_metrics_path = write_condition_report(
        sim_base,
        conditions,
        assignment_rows,
    )

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
        f"\nDetailed per-condition report: {condition_report_path}",
        f"Per-condition metrics TSV: {condition_metrics_path}",
    ])
    report_path.write_text(full_report)
    print(f"\n\n  Full report saved to: {report_path}")
    print(f"  Detailed condition report saved to: {condition_report_path}")
    print(f"  Condition metrics TSV saved to: {condition_metrics_path}")
    print("=" * 100)


if __name__ == "__main__":
    main()
