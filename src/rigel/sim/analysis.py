#!/usr/bin/env python3
"""Run rigel quantification on synthetic simulation conditions and analyze results.

Builds the rigel index, runs quant (with annotated BAM) on every condition, then reports:
  1. Net fragment-flow deconvolution — the PRIMARY gDNA↔RNA accuracy metric (per-locus net
     leak, per-transcript Δ decomposed into gDNA-source vs RNA-isoform-source; see
     :func:`analyze_net_flow`).
  2. Transcript abundance accuracy (correlation / MARD / false positives).
  3. Calibration scalars + gDNA fragment-length models.
  4. Gross per-fragment confusion (hard labels — secondary) and cheap acceptance checks.
"""

import argparse
import json
import logging
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
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


@dataclass(frozen=True)
class CalibrationAcceptanceThresholds:
    """Guardrails for the cheap synthetic acceptance checks."""

    max_nrna_none_count: float = 10.0
    max_gdna_to_rna_leak_rate: float = 0.015


def parse_condition(cond_name: str) -> dict:
    """Parse a ``gdna_<lbl>_ss_<val>_nrna_<lbl>[...]`` condition name into its parts.

    Structural only — the gDNA *rate* comes from the manifest (``gdna_rate``), never the name.
    """
    parts = cond_name.split("_")
    return {
        "condition": cond_name,
        "strand_specificity": float(parts[3]),
        "gdna_label": f"{parts[0]}_{parts[1]}",
        "nrna_label": parts[5] if len(parts) > 5 else "none",
    }


def _contamination_key(meta: dict) -> tuple[float, float]:
    """Sort key over conditions: (gDNA rate, −strand_specificity).

    ``max`` → dirtiest + least stranded (hardest); ``min`` → cleanest + most stranded.
    """
    return (float(meta.get("gdna_rate", 0.0)), -float(meta.get("strand_specificity", 0.0)))


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


def _nested(data: dict | None, *keys: str, default: object = None) -> object:
    current: object = data or {}
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


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
        "rigel",
        "index",
        "--fasta",
        str(ref_dir / "genome.fa"),
        "--gtf",
        str(ref_dir / "genes.gtf"),
        "-o",
        str(index_dir),
        "--gtf-parse-mode",
        "warn-skip",
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
        "rigel",
        "quant",
        "--bam",
        str(cond_dir / "sim_oracle.bam"),
        "--index",
        str(index_dir),
        "-o",
        str(out_dir),
        "--annotated-bam",
        str(annotated_bam),
        "--sj-strand-tag",
        "auto",
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
    lines.append(
        f"\n{'Condition':<35} {'gDNA_frac':>9} {'SS':>5} "
        f"{'dens_glob':>10} {'rna_sense':>9} {'od_gdna':>8} {'od_rna':>8} "
        f"{'FL_rna':>7} {'FL_gdna':>8} "
        f"{'n_loci':>6} {'gdna_rate':>10} {'gdna_true':>10}"
    )
    lines.append("─" * 140)

    cal_rows = []
    for cond in conditions:
        info = parse_condition(cond)
        meta = condition_meta.get(cond, {})
        gdna_rate = float(meta.get("gdna_rate", 0.0))
        if meta:
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

        # Compute true gDNA rate (manifest counts; fall back to the configured rate).
        n_rna = int(meta.get("n_rna", 1_000_000))
        n_gdna_true = int(meta.get("n_gdna", round(n_rna * gdna_rate)))
        n_total_true = n_rna + n_gdna_true
        gdna_rate_true = n_gdna_true / n_total_true if n_total_true > 0 else 0
        truth_summary = _load_condition_truth_summary(sim_base, meta)
        truth_rna_fl_mean = _truth_fl_mean(truth_summary, "mrna", true_rna_fl_mean)
        truth_gdna_fl_mean = _truth_fl_mean(truth_summary, "gdna", true_gdna_fl_mean)

        gdna_rate_est = quant_out.get("gdna_fraction", 0)

        lines.append(
            f"{cond:<35} {gdna_rate:>9.2f} {info['strand_specificity']:>5.2f} "
            f"{density_global:>10.6f} {rna_sense_frac:>9.3f} {od_gdna:>8.3f} {od_rna:>8.3f} "
            f"{fl_rna:>7.1f} {fl_gdna:>8.1f} "
            f"{n_loci:>6} {gdna_rate_est:>10.4f} {gdna_rate_true:>10.4f}"
        )

        cal_rows.append(
            {
                **info,
                "density_global": density_global,
                "rna_sense_frac": rna_sense_frac,
                "od_gdna": od_gdna,
                "od_rna": od_rna,
                "fl_rna": fl_rna,
                "fl_gdna": fl_gdna,
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
            }
        )

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
    lines.append(
        f"  {'Condition':<35} {'FL_rna':>8} {'FL_gdna':>8} {'FL_rna_err':>10} {'FL_gdna_err':>11}"
    )
    for row in cal_rows:
        rna_err = (
            (row["fl_rna"] - row["truth_rna_fl_mean"]) / row["truth_rna_fl_mean"]
            if row["fl_rna"] > 0 and row["truth_rna_fl_mean"] > 0
            else float("nan")
        )
        gdna_err = (
            (row["fl_gdna"] - row["truth_gdna_fl_mean"]) / row["truth_gdna_fl_mean"]
            if (row["n_gdna_true"] > 0 and row["fl_gdna"] > 0 and row["truth_gdna_fl_mean"] > 0)
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
    lines.append(
        f"\n{'Condition':<35} {'Spearman':>8} {'Pearson':>8} "
        f"{'MARD':>8} {'Med_RE':>8} "
        f"{'n_FP':>5} {'n_FN':>5} "
        f"{'gdna_leak':>9}"
    )
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
            quant[["transcript_id", "count", "count_em", "tpm"]], on="transcript_id", how="left"
        ).fillna(0)

        # Correlation on expressed transcripts only
        expressed = merged[merged["mrna_abundance"] > 0]
        unexpressed = merged[merged["mrna_abundance"] == 0]

        if len(expressed) > 2:
            sp_r, _ = spearmanr(expressed["mrna_abundance"], expressed["tpm"])
            pe_r, _ = pearsonr(
                np.log2(expressed["mrna_abundance"] + 1), np.log2(expressed["tpm"] + 1)
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

    if abundance_details:
        pd.DataFrame(
            [
                {
                    "condition": cond,
                    "spearman": d["spearman"],
                    "pearson": d["pearson"],
                    "mard": d["mard"],
                    "n_false_pos": d["n_fp"],
                    "n_false_neg": d["n_fn"],
                    "gdna_leak_count": d["gdna_leak"],
                }
                for cond, d in abundance_details.items()
            ]
        ).to_csv(sim_base / "abundance_per_condition.tsv", sep="\t", index=False)

    # Dirtiest / cleanest conditions, chosen from the manifest rather than hard-coded.
    evaluated = [c for c in conditions if c in abundance_details]
    worst_cond = (
        max(evaluated, key=lambda c: _contamination_key(condition_meta.get(c, {})))
        if evaluated
        else None
    )
    clean_cond = (
        min(evaluated, key=lambda c: _contamination_key(condition_meta.get(c, {})))
        if evaluated
        else None
    )

    # ── Per-gene breakdown for worst cases ──
    lines.append(f"\n\n{'─' * 100}")
    lines.append("  TOP FALSE POSITIVES (unexpressed txs with highest erroneous counts)")
    lines.append(f"{'─' * 100}")

    if worst_cond is not None:
        merged = abundance_details[worst_cond]["merged"]
        fp = merged[merged["mrna_abundance"] == 0].nlargest(15, "count")
        lines.append(f"\n  Condition: {worst_cond}")
        lines.append(
            f"  {'transcript_id':<20} {'gene_id':<12} {'n_exons':>7} "
            f"{'spliced_len':>11} {'count':>8} {'tpm':>8}"
        )
        for _, row in fp.iterrows():
            lines.append(
                f"  {row['transcript_id']:<20} {row['gene_id']:<12} "
                f"{int(row['n_exons']):>7} {int(row['spliced_length']):>11} "
                f"{row['count']:>8.1f} {row['tpm']:>8.2f}"
            )

    # ── Abundance accuracy by expression level ──
    lines.append(f"\n\n{'─' * 100}")
    lines.append(f"  ACCURACY BY EXPRESSION LEVEL ({clean_cond or 'n/a'} — cleanest condition)")
    lines.append(f"{'─' * 100}")

    if clean_cond is not None:
        merged = abundance_details[clean_cond]["merged"]
        expressed = merged[merged["mrna_abundance"] > 0].copy()
        expressed["log_abundance"] = np.log10(expressed["mrna_abundance"] + 1)
        bins = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
        expressed["bin"] = pd.cut(expressed["log_abundance"], bins=bins)

        lines.append(
            f"\n  {'log10(TPM) bin':<18} {'n_tx':>5} {'mean_RE':>8} {'med_RE':>8} "
            f"{'mean_count':>10} {'mean_true':>10}"
        )
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

    lines.append(
        f"\n{'Condition':<35} {'n_loci':>6} {'tot_mrna':>9} "
        f"{'tot_nrna':>9} {'tot_gdna':>9} {'tot_unambig':>11} "
        f"{'gdna_rate':>9} {'max_gdna':>9}"
    )
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

    # Locus details for the dirtiest condition (highest gDNA, least stranded).
    condition_meta = condition_manifest_map(load_manifest(sim_base))
    evaluated = [c for c in conditions if (sim_base / c / "rigel_out").exists()]
    cond = (
        max(evaluated, key=lambda c: _contamination_key(condition_meta.get(c, {})))
        if evaluated
        else None
    )
    lines.append(f"\n\n{'─' * 100}")
    lines.append(f"  PER-LOCUS DETAIL ({cond or 'n/a'} — top 20 loci by gDNA)")
    lines.append(f"{'─' * 100}")

    if cond is not None:
        loci = load_loci(sim_base / cond / "rigel_out")
        if not loci.empty:
            top = loci.nlargest(20, "gdna")
            cols = [
                "locus_id",
                "n_transcripts",
                "locus_span_bp",
                "count_unambig",
                "mrna",
                "nrna",
                "gdna",
                "gdna_rate",
                "gdna_prior_count",
            ]
            avail_cols = [c for c in cols if c in top.columns]
            lines.append(f"\n  {top[avail_cols].to_string(index=False)}")

    return "\n".join(lines)


# ══════════════════════════════════════════════════════════════════════════
# Net fragment-flow deconvolution (the primary gDNA-vs-RNA accuracy metric)
# ══════════════════════════════════════════════════════════════════════════
#
# Hard per-fragment label recovery is the WRONG target: an unspliced RNA fragment
# and a gDNA fragment from the same locus can be sequence-identical and are
# genuinely unrecoverable. What matters is the *net* deconvolution error per
# component. So we build a per-locus fragment-flow matrix flow[true][assigned]
# from each fragment's TRUE origin (oracle read name) and the tool's HARD MAP
# label (annotated-BAM ZT/ZF), then reduce to NET flow net(a→b)=flow[a][b]-flow[b][a].
# Symmetric (unrecoverable) misassignment cancels to ~0; only systematic bias
# survives — which is exactly the soft/pool-equivalent quantity. The per-fragment
# hard label is just the datum; the net aggregate is the answer.
#
# Components in a locus L = its transcript rows (mRNA + synthetic nRNA spans) plus
# one gdna@L component; intergenic gDNA is gdna@-1. Then for every component
#     observed[c] - expected[c] = Σ_a net(a→c)   (net inflow),
# and each transcript's surplus/deficit decomposes into net flow from other RNA
# isoforms (RNA↔RNA reallocation) vs net flow from gDNA (contamination).


@dataclass
class FlowData:
    """Per-condition fragment-flow matrix over integer-keyed components."""

    condition: str
    flow: dict[tuple[int, int], int]  # (true_cid, assigned_cid) -> fragment count
    comp_name: dict[int, str]  # cid -> transcript_id | "gdna@<locus>" | "unassigned"
    comp_kind: dict[int, str]  # cid -> "rna" | "gdna" | "unassigned"
    comp_locus: dict[int, int]  # cid -> locus_id (gdna@L -> L; unknown -> -1)
    total_gdna_true: int = 0  # true gDNA fragments (oracle), all loci + intergenic
    # 3-pool fragment flow (true_pool, assigned_pool) -> count, pools in {gdna, nrna, mrna}. The
    # net reduction (net(a→b)=flow[a,b]-flow[b,a]) cancels sequence-identical unrecoverable
    # misassignment; per-pool net surplus = assigned-true is the un-conflated leak/FP measure.
    pool_flow: dict[tuple[str, str], int] = field(default_factory=dict)


def collect_fragment_flows(
    sim_base: Path, conditions: list[str]
) -> tuple[dict[str, FlowData], list[dict]]:
    """Single pass over each condition's annotated BAM.

    Returns ``(flows_by_condition, overview_rows)``. ``overview_rows`` reproduces the
    legacy gross-confusion schema (consumed by the condition report + acceptance checks);
    ``flows_by_condition`` carries the sparse per-locus flow matrix for the net analysis.
    """
    import pysam

    flows_by_cond: dict[str, FlowData] = {}
    overview_rows: list[dict] = []

    for cond in conditions:
        cond_dir = sim_base / cond
        annotated_bam = cond_dir / "annotated.bam"
        if not annotated_bam.exists():
            continue

        # Component registry (lazy integer ids).
        key2cid: dict[tuple, int] = {}
        comp_name: dict[int, str] = {}
        comp_kind: dict[int, str] = {}
        comp_locus: dict[int, int] = {}
        # RNA components' home locus is the modal ZL across fragments touching them
        # (BAM-derived, self-consistent with the ZL-keyed gDNA components). This puts
        # nRNA-shadow spans — which carry locus_id=-1 in quant.feather — into their
        # true locus, so mRNA→nRNA flow stays within-locus instead of leaking to cross_locus.
        comp_zl: dict[int, Counter] = defaultdict(Counter)

        def cid_for_transcript(tid: str) -> int:
            key = ("t", tid)
            cid = key2cid.get(key)
            if cid is None:
                cid = len(key2cid)
                key2cid[key] = cid
                comp_name[cid] = tid
                comp_kind[cid] = "rna"
                comp_locus[cid] = -1  # finalized to modal ZL after the pass
            return cid

        def cid_for_gdna(locus: int) -> int:
            key = ("g", locus)
            cid = key2cid.get(key)
            if cid is None:
                cid = len(key2cid)
                key2cid[key] = cid
                comp_name[cid] = f"gdna@{locus}"
                comp_kind[cid] = "gdna"
                comp_locus[cid] = locus
            return cid

        def cid_unassigned() -> int:
            key = ("u",)
            cid = key2cid.get(key)
            if cid is None:
                cid = len(key2cid)
                key2cid[key] = cid
                comp_name[cid] = "unassigned"
                comp_kind[cid] = "unassigned"
                comp_locus[cid] = -1
            return cid

        flow: Counter = Counter()
        pool_flow: Counter = Counter()  # (true_pool, assigned_pool) -> count, pools gdna/nrna/mrna

        # Legacy gross-confusion counters (back-compat for downstream consumers).
        correct_tx = correct_gene = wrong_tx = 0
        rna_as_gdna = gdna_as_rna = gdna_correct = 0
        total_rna = total_mrna = total_nrna = total_gdna = 0
        nrna_as_rna = nrna_as_gdna = mrna_as_gdna = mrna_as_nrna = 0

        with pysam.AlignmentFile(str(annotated_bam), "rb") as bam:
            for read in bam:
                if read.is_read2 or read.is_secondary or read.is_supplementary:
                    continue

                origin = parse_origin(read.query_name)

                def _tag(name: str, default):
                    try:
                        return read.get_tag(name)
                    except KeyError:
                        return default

                zt = str(_tag("ZT", "") or "")
                zg = str(_tag("ZG", "") or "")
                category = str(_tag("ZC", "") or "")
                zf = int(_tag("ZF", 0) or 0)
                zl = int(_tag("ZL", -1))
                # ZF bit 2 (0x04) = gDNA pool; intergenic also carries it but guard on ZC too.
                is_assigned_gdna = bool(zf & 0x04) or category == "intergenic"
                is_assigned_nrna = bool(zf & 0x08)

                # Assigned component (the tool's hard MAP destination).
                if is_assigned_gdna:
                    a_cid = cid_for_gdna(zl)
                elif zt and zt != ".":
                    a_cid = cid_for_transcript(zt)
                else:
                    a_cid = cid_unassigned()

                # True component (oracle) + legacy gross counters.
                if origin.kind == "gdna":
                    total_gdna += 1
                    t_cid = cid_for_gdna(zl)
                    if is_assigned_gdna:
                        gdna_correct += 1
                    else:
                        gdna_as_rna += 1
                else:
                    true_tx = origin.transcript_id or ""
                    t_cid = cid_for_transcript(true_tx)
                    total_rna += 1
                    if origin.kind == "nrna":
                        total_nrna += 1
                        if is_assigned_gdna:
                            rna_as_gdna += 1
                            nrna_as_gdna += 1
                        else:
                            nrna_as_rna += 1
                    else:
                        total_mrna += 1
                        if is_assigned_gdna:
                            rna_as_gdna += 1
                            mrna_as_gdna += 1
                        elif is_assigned_nrna:
                            mrna_as_nrna += 1
                        elif zt == true_tx:
                            correct_tx += 1
                        elif zg and true_tx.rsplit(".", 1)[0] == zg:
                            correct_gene += 1
                        else:
                            wrong_tx += 1

                flow[(t_cid, a_cid)] += 1
                assigned_pool = (
                    "gdna" if is_assigned_gdna else ("nrna" if is_assigned_nrna else "mrna")
                )
                pool_flow[(origin.kind, assigned_pool)] += 1
                # Vote the fragment's genomic locus (ZL) for any RNA endpoint.
                if comp_kind[t_cid] == "rna":
                    comp_zl[t_cid][zl] += 1
                if comp_kind[a_cid] == "rna":
                    comp_zl[a_cid][zl] += 1

        # Finalize RNA components' home locus to their modal ZL.
        for cid, counter in comp_zl.items():
            if counter:
                comp_locus[cid] = counter.most_common(1)[0][0]

        flows_by_cond[cond] = FlowData(
            condition=cond,
            flow=dict(flow),
            comp_name=comp_name,
            comp_kind=comp_kind,
            comp_locus=comp_locus,
            total_gdna_true=total_gdna,
            pool_flow=dict(pool_flow),
        )
        overview_rows.append(
            {
                "condition": cond,
                "total": total_rna + total_gdna,
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
            }
        )

    return flows_by_cond, overview_rows


def collect_fragment_assignment_rows(sim_base: Path, conditions: list[str]) -> list[dict]:
    """Legacy gross-confusion rows (thin wrapper over :func:`collect_fragment_flows`)."""
    return collect_fragment_flows(sim_base, conditions)[1]


# ── Net-flow reduction primitives ─────────────────────────────────────────


def _flow_marginals(flow: dict[tuple[int, int], int]) -> tuple[dict[int, int], dict[int, int]]:
    """Return (expected[c]=row sum from c, observed[c]=col sum into c)."""
    expected: Counter = Counter()
    observed: Counter = Counter()
    for (a, b), n in flow.items():
        expected[a] += n
        observed[b] += n
    return dict(expected), dict(observed)


def _net_flow_rows(fd: FlowData) -> tuple[list[dict], list[dict]]:
    """Reduce one condition's flow matrix to per-transcript and per-locus net-flow rows.

    Per transcript T (home locus L):
      delta = observed - expected
            = net_from_gdna(gdna@L → T) + net_from_rna_isoforms(Σ T'≠T in L) + cross_locus
    """
    flow = fd.flow
    expected, observed = _flow_marginals(flow)

    def net(a: int, b: int) -> int:
        return flow.get((a, b), 0) - flow.get((b, a), 0)

    # Index components by home locus.
    rna_by_locus: dict[int, list[int]] = {}
    gdna_by_locus: dict[int, int] = {}
    for cid, kind in fd.comp_kind.items():
        loc = fd.comp_locus[cid]
        if kind == "rna":
            rna_by_locus.setdefault(loc, []).append(cid)
        elif kind == "gdna":
            gdna_by_locus[loc] = cid

    tx_rows: list[dict] = []
    for loc, rna_cids in rna_by_locus.items():
        gdna_cid = gdna_by_locus.get(loc)
        for cid in rna_cids:
            exp = expected.get(cid, 0)
            obs = observed.get(cid, 0)
            delta = obs - exp
            net_gdna = net(gdna_cid, cid) if gdna_cid is not None else 0
            net_iso = sum(net(other, cid) for other in rna_cids if other != cid)
            tx_rows.append(
                {
                    "condition": fd.condition,
                    "transcript_id": fd.comp_name[cid],
                    "locus_id": loc,
                    "expected": exp,
                    "observed": obs,
                    "delta": delta,
                    "net_from_gdna": net_gdna,
                    "net_from_rna_isoforms": net_iso,
                    "cross_locus": delta - net_gdna - net_iso,
                }
            )

    locus_rows: list[dict] = []
    all_loci = set(rna_by_locus) | set(gdna_by_locus)
    for loc in all_loci:
        rna_cids = rna_by_locus.get(loc, [])
        gdna_cid = gdna_by_locus.get(loc)
        rna_exp = sum(expected.get(c, 0) for c in rna_cids)
        rna_obs = sum(observed.get(c, 0) for c in rna_cids)
        gdna_exp = expected.get(gdna_cid, 0) if gdna_cid is not None else 0
        gdna_obs = observed.get(gdna_cid, 0) if gdna_cid is not None else 0
        locus_rows.append(
            {
                "condition": fd.condition,
                "locus_id": loc,
                "n_transcripts": len(rna_cids),
                "rna_expected": rna_exp,
                "rna_observed": rna_obs,
                "rna_delta": rna_obs - rna_exp,
                "gdna_expected": gdna_exp,
                "gdna_observed": gdna_obs,
                # + => gDNA mass leaked to RNA; - => RNA siphoned into gDNA.
                "net_gdna_to_rna": gdna_exp - gdna_obs,
            }
        )

    return tx_rows, locus_rows


_POOLS_3 = ("gdna", "nrna", "mrna")


def _pool_flow_3way_row(fd: FlowData) -> dict:
    """One per-condition 3-pool (gDNA/nRNA/mRNA) net-flow row.

    Reports each pool's true/assigned totals + **net surplus** (assigned − true; + = the pool is
    net-inflated, − = net-deficit) and the three **net fluxes** between pool pairs. Net cancels the
    sequence-identical, truly-unidentifiable misassignment, so the per-pool net surplus is the
    un-conflated leak / false-positive measure (a gross sum over-counts unrecoverable exchange).
    """
    pf = fd.pool_flow

    def cnt(a: str, b: str) -> int:
        return pf.get((a, b), 0)

    true = {p: sum(cnt(p, b) for b in _POOLS_3) for p in _POOLS_3}
    asg = {p: sum(cnt(a, p) for a in _POOLS_3) for p in _POOLS_3}
    row = {"condition": fd.condition}
    for p in _POOLS_3:
        row[f"{p}_true"] = true[p]
        row[f"{p}_assigned"] = asg[p]
        row[f"{p}_net_surplus"] = asg[p] - true[p]
    row["net_gdna_to_nrna"] = cnt("gdna", "nrna") - cnt("nrna", "gdna")
    row["net_gdna_to_mrna"] = cnt("gdna", "mrna") - cnt("mrna", "gdna")
    row["net_nrna_to_mrna"] = cnt("nrna", "mrna") - cnt("mrna", "nrna")
    return row


def _spearman(x: "pd.Series", y: "pd.Series") -> float:
    """Robust Spearman that returns nan on degenerate input."""
    mask = x.notna() & y.notna()
    if mask.sum() < 5 or x[mask].nunique() < 2 or y[mask].nunique() < 2:
        return float("nan")
    from scipy.stats import spearmanr

    r, _ = spearmanr(x[mask], y[mask])
    return float(r)


def analyze_net_flow(
    sim_base: Path,
    conditions: list[str],
    flows: dict[str, FlowData] | None = None,
) -> str:
    """Primary gDNA-vs-RNA accuracy: net fragment-flow deconvolution.

    Writes ``net_flow_per_transcript.tsv`` and ``net_flow_per_locus.tsv`` (covariate-joined)
    and returns a text report: pool net-leak & direction, per-transcript Δ distribution,
    covariate root-cause ranking, and the identifiability diagnostic.
    """
    if flows is None:
        flows, _ = collect_fragment_flows(sim_base, conditions)

    cmeta = condition_manifest_map(load_manifest(sim_base))
    fallback_truth = load_truth(sim_base)

    tx_frames: list[pd.DataFrame] = []
    locus_frames: list[pd.DataFrame] = []
    for cond in conditions:
        fd = flows.get(cond)
        if fd is None:
            continue
        tx_rows, locus_rows = _net_flow_rows(fd)
        meta = cmeta.get(cond, {})
        tag = {
            "gdna_label": meta.get("gdna_label", parse_condition(cond)["gdna_label"]),
            "ss": meta.get("strand_specificity", parse_condition(cond)["strand_specificity"]),
            "capture": meta.get("capture_label", "off"),
        }

        if tx_rows:
            tdf = pd.DataFrame(tx_rows)
            # Transcript covariates from truth (n_exons, spliced_length, gene_id, strand).
            ctruth = load_condition_truth(sim_base, cond, cmeta, fallback_truth)
            cov_cols = [
                c
                for c in (
                    "transcript_id",
                    "gene_id",
                    "n_exons",
                    "spliced_length",
                    "strand",
                    "mrna_abundance",
                )
                if c in ctruth.columns
            ]
            tdf = tdf.merge(ctruth[cov_cols], on="transcript_id", how="left")
            if "n_exons" in tdf.columns:
                tdf["single_exon"] = tdf["n_exons"].fillna(0) <= 1
            for k, v in tag.items():
                tdf[k] = v
            tx_frames.append(tdf)

        if locus_rows:
            ldf = pd.DataFrame(locus_rows)
            loci_out = load_loci(sim_base / cond / "rigel_out")
            cov_cols = [
                c
                for c in (
                    "locus_id",
                    "locus_span_bp",
                    "n_em_fragments",
                    "gdna_prior_count",
                    "rna_prior_count",
                    "gdna_eff_len_em",
                )
                if c in loci_out.columns
            ]
            if cov_cols and not loci_out.empty:
                ldf = ldf.merge(loci_out[cov_cols], on="locus_id", how="left")
            for k, v in tag.items():
                ldf[k] = v
            locus_frames.append(ldf)

    lines: list[str] = []
    hr = "═" * 100
    lines.append(f"\n{hr}")
    lines.append("  NET FRAGMENT-FLOW DECONVOLUTION  (primary gDNA↔RNA accuracy)")
    lines.append(hr)
    lines.append(
        "  Net flow cancels symmetric (sequence-identical, unrecoverable) misassignment;\n"
        "  only systematic bias remains. + net_gdna_to_rna ⇒ gDNA leaked into RNA;\n"
        "  - ⇒ RNA siphoned into gDNA. Per transcript, Δ = net(gDNA→T) + net(RNA isoforms→T).\n"
    )

    if not tx_frames:
        lines.append("  [no annotated BAMs / loci available — run quant with --annotated-bam]")
        return "\n".join(lines)

    tx_all = pd.concat(tx_frames, ignore_index=True)
    locus_all = pd.concat(locus_frames, ignore_index=True) if locus_frames else pd.DataFrame()

    tx_path = sim_base / "net_flow_per_transcript.tsv"
    locus_path = sim_base / "net_flow_per_locus.tsv"
    tx_all.to_csv(tx_path, sep="\t", index=False)
    if not locus_all.empty:
        locus_all.to_csv(locus_path, sep="\t", index=False)

    # ── Pool net-leak & direction, per condition ──
    lines.append("  POOL NET LEAK (signed; fraction of true gDNA), by condition:")
    lines.append(
        f"  {'gdna':9}{'cap':5}{'ss':6} | {'true_gDNA':>10} {'net→RNA':>9} "
        f"{'leak%':>7} | {'in_locus':>9} {'intergenic':>11}"
    )
    lines.append("  " + "-" * 74)
    pool_rows = []
    for cond in conditions:
        fd = flows.get(cond)
        if fd is None:
            continue
        meta = cmeta.get(cond, {})
        sub = locus_all[locus_all["condition"] == cond] if not locus_all.empty else pd.DataFrame()
        in_locus_net = (
            int(sub[sub["locus_id"] >= 0]["net_gdna_to_rna"].sum()) if not sub.empty else 0
        )
        intergenic_net = (
            int(sub[sub["locus_id"] < 0]["net_gdna_to_rna"].sum()) if not sub.empty else 0
        )
        net_total = in_locus_net + intergenic_net
        true_gdna = fd.total_gdna_true
        leak_frac = net_total / true_gdna if true_gdna else 0.0
        pool_rows.append(
            {
                "gdna_label": meta.get("gdna_label", "?"),
                "gdna_rate": float(meta.get("gdna_rate", 0.0)),
                "capture": meta.get("capture_label", "off"),
                "ss": meta.get("strand_specificity", float("nan")),
                "true_gdna": true_gdna,
                "net_to_rna": net_total,
                "leak_frac": leak_frac,
                "in_locus": in_locus_net,
                "intergenic": intergenic_net,
            }
        )
    pool_rows.sort(key=lambda r: (r["gdna_rate"], r["capture"], -r["ss"]))
    for r in pool_rows:
        lines.append(
            f"  {r['gdna_label']:9}{r['capture']:5}{r['ss']:<6.2f} | {r['true_gdna']:>10,} "
            f"{r['net_to_rna']:>+9,} {r['leak_frac'] * 100:>6.2f}% | "
            f"{r['in_locus']:>+9,} {r['intergenic']:>+11,}"
        )
    cond_path = sim_base / "net_flow_per_condition.tsv"
    pd.DataFrame(pool_rows).to_csv(cond_path, sep="\t", index=False)

    # ── 3-pool net flow: gDNA ↔ nRNA ↔ mRNA (the resurrected nRNA pool needs 3 pools) ──
    lines.append(
        "\n  3-POOL NET FLOW (gDNA/nRNA/mRNA): per-pool net surplus (assigned−true) + net pair fluxes."
    )
    lines.append(
        "  + surplus ⇒ pool net-inflated (false gain); − ⇒ net-deficit. Watch mRNA surplus (mature FP)."
    )
    lines.append(
        f"  {'gdna':9}{'cap':4}{'ss':5}{'nrna':5} | {'gDNA_surp':>10}{'nRNA_surp':>10}{'mRNA_surp':>10}"
        f" | {'g→n':>8}{'g→m':>8}{'n→m':>8}"
    )
    lines.append("  " + "-" * 86)
    pool3_rows = []
    for cond in conditions:
        fd = flows.get(cond)
        if fd is None or not fd.pool_flow:
            continue
        meta = cmeta.get(cond, {})
        row = _pool_flow_3way_row(fd)
        row.update(
            {
                "gdna_label": meta.get("gdna_label", "?"),
                "gdna_rate": float(meta.get("gdna_rate", 0.0)),
                "capture": meta.get("capture_label", "off"),
                "ss": meta.get("strand_specificity", float("nan")),
                "nrna": meta.get("nrna_label", "rnd" if "nrna_rnd" in cond else "none"),
            }
        )
        pool3_rows.append(row)
    pool3_rows.sort(key=lambda r: (r["gdna_rate"], r["capture"], -r["ss"], r["nrna"]))
    for r in pool3_rows:
        lines.append(
            f"  {r['gdna_label']:9}{r['capture']:4}{r['ss']:<5.2f}{r['nrna']:5} | "
            f"{r['gdna_net_surplus']:>+10,}{r['nrna_net_surplus']:>+10,}{r['mrna_net_surplus']:>+10,}"
            f" | {r['net_gdna_to_nrna']:>+8,}{r['net_gdna_to_mrna']:>+8,}{r['net_nrna_to_mrna']:>+8,}"
        )
    if pool3_rows:
        pd.DataFrame(pool3_rows).to_csv(
            sim_base / "net_flow_3pool_per_condition.tsv", sep="\t", index=False
        )

    # ── Per-transcript Δ distribution + source decomposition ──
    lines.append("\n  PER-TRANSCRIPT Δ (observed − expected) and its decomposition:")
    lines.append(
        f"  {'gdna':9}{'cap':5}{'ss':6} | {'n_tx':>5} {'meanΔ':>7} {'sdΔ':>7} "
        f"{'mean|Δ|':>8} | {'fromGDNA':>9} {'fromISO':>8} | {'|Δ|>10':>6}"
    )
    lines.append("  " + "-" * 78)
    for (gl, cap, ss), grp in tx_all.groupby(["gdna_label", "capture", "ss"], dropna=False):
        lines.append(
            f"  {str(gl):9}{str(cap):5}{float(ss):<6.2f} | {len(grp):>5} "
            f"{grp['delta'].mean():>7.2f} {grp['delta'].std():>7.2f} "
            f"{grp['delta'].abs().mean():>8.2f} | "
            f"{grp['net_from_gdna'].mean():>+9.2f} {grp['net_from_rna_isoforms'].mean():>+8.2f} | "
            f"{int((grp['delta'].abs() > 10).sum()):>6}"
        )

    # ── Root cause: covariate ranking against gDNA contamination inflow ──
    lines.append(
        "\n  ROOT CAUSE — Spearman(net_from_gdna, covariate) over transcripts in gDNA>0 conditions:"
    )
    contam = tx_all[tx_all["gdna_label"] != "none"].copy()
    if not contam.empty:
        cov_candidates = {
            "true_RNA_depth(expected)": contam.get("expected"),
            "n_exons": contam.get("n_exons"),
            "spliced_length": contam.get("spliced_length"),
            "single_exon": contam.get("single_exon").astype(float)
            if "single_exon" in contam.columns
            else None,
            "mrna_abundance": contam.get("mrna_abundance"),
        }
        ranked = []
        for name, series in cov_candidates.items():
            if series is None:
                continue
            ranked.append((name, _spearman(contam["net_from_gdna"], series)))
        ranked = [r for r in ranked if not np.isnan(r[1])]
        ranked.sort(key=lambda r: abs(r[1]), reverse=True)
        for name, r in ranked:
            lines.append(f"    {name:<28} ρ = {r:+.3f}")
        if not ranked:
            lines.append("    (insufficient variation to rank covariates)")

    # ── Identifiability diagnostic: gross confusion vs net (expected-unrecoverable vs bias) ──
    lines.append("\n  IDENTIFIABILITY — single-exon (gDNA-identical) vs multi-exon transcripts:")
    if "single_exon" in contam.columns and not contam.empty:
        lines.append(
            f"    {'class':<12} {'n_tx':>6} {'mean|Δ|':>8} {'meanΔ':>8} {'mean net_from_gdna':>18}"
        )
        for label, mask in (
            ("single-exon", contam["single_exon"]),
            ("multi-exon", ~contam["single_exon"]),
        ):
            grp = contam[mask]
            if grp.empty:
                continue
            lines.append(
                f"    {label:<12} {len(grp):>6} {grp['delta'].abs().mean():>8.2f} "
                f"{grp['delta'].mean():>8.2f} {grp['net_from_gdna'].mean():>18.2f}"
            )
        lines.append(
            "    (single-exon transcripts are sequence-identical to gDNA and have no isoforms, so\n"
            "     their Δ ≈ net_from_gdna; compare gDNA inflow across classes to localize where\n"
            "     contamination concentrates. meanΔ near 0 with large spread ⇒ unbiased-but-noisy;\n"
            "     meanΔ ≈ mean net_from_gdna ⇒ systematic gDNA→RNA leak on that class.)"
        )

    lines.append(f"\n  Wrote {cond_path}")
    lines.append(f"  Wrote {tx_path}")
    if not locus_all.empty:
        lines.append(f"  Wrote {locus_path}")
    return "\n".join(lines)


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
    lines.append(
        f"\n{'Condition':<35} {'Total':>8} "
        f"{'mRNA_ok':>8} {'mRNA_gene':>9} {'mRNA_bad':>9} "
        f"{'nRNA→RNA':>9} {'nRNA→gDNA':>10} "
        f"{'gDNA_ok':>8} {'gDNA→RNA':>9} {'Accuracy':>8}"
    )
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
            lines.append(
                f"      Wrong gene:       {row['wrong_tx'] / total_mrna:.4f} ({row['wrong_tx']:,} frags)"
            )

        total_nrna = row.get("total_nrna", 0)
        if total_nrna > 0:
            nrna_to_rna = row.get("nrna_as_rna", 0) / total_nrna
            nrna_to_gdna = row.get("nrna_as_gdna", 0) / total_nrna
            lines.append(f"    nRNA ({total_nrna:,} frags):")
            lines.append(
                f"      Routed to RNA-compatible tags:  {nrna_to_rna:.4f} ({row.get('nrna_as_rna', 0):,} frags)"
            )
            lines.append(
                f"      Misclass as gDNA/intergenic:    {nrna_to_gdna:.4f} ({row.get('nrna_as_gdna', 0):,} frags)"
            )

        if row["total_gdna"] > 0:
            gdna_precision = row["gdna_correct"] / row["total_gdna"]
            gdna_leak = row["gdna_as_rna"] / row["total_gdna"]
            lines.append(f"    gDNA ({row['total_gdna']:,} frags):")
            lines.append(
                f"      Correctly identified: {gdna_precision:.4f} ({row['gdna_correct']:,} frags)"
            )
            lines.append(
                f"      Leaked to RNA:        {gdna_leak:.4f} ({row['gdna_as_rna']:,} frags)"
            )

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
        rows.append(
            {
                "check": check,
                "condition": condition,
                "value": value,
                "threshold": threshold,
                "status": status,
            }
        )

    # nRNA should stay essentially off in synthetic nrna_none conditions.
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
        "  Cheap regression gates: nrna_none should stay near zero, and true gDNA should "
        "not leak substantially into RNA (hard-label rate; see net-flow for the net metric).\n"
    )
    lines.append(f"  {'Check':<20} {'Condition':<35} {'Value':>10} {'Threshold':>14} {'Status':>7}")
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
    parser.add_argument(
        "--skip-quant", action="store_true", help="Skip quantification, only run analysis"
    )
    parser.add_argument(
        "--skip-frag-analysis", action="store_true", help="Skip fragment-level analysis (slow)"
    )
    parser.add_argument(
        "--conditions",
        nargs="*",
        default=None,
        help="Optional subset of condition names to evaluate",
    )
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

    # Single BAM pass: net-flow matrices + legacy gross-confusion rows.
    flows: dict[str, FlowData] = {}
    assignment_rows: list[dict] = []
    if not args.skip_frag_analysis:
        flows, assignment_rows = collect_fragment_flows(sim_base, conditions)

    # 3a: PRIMARY — net fragment-flow deconvolution (gDNA↔RNA).
    if not args.skip_frag_analysis:
        net_flow_report = analyze_net_flow(sim_base, conditions, flows)
    else:
        net_flow_report = "\n  [SKIP] Net-flow analysis skipped (--skip-frag-analysis)"
    print(net_flow_report)

    # 3b: Transcript abundance accuracy (RNA quant view).
    abundance_report = analyze_abundance(sim_base, conditions, truth)
    print(abundance_report)

    # 3c: Calibration scalars.
    cal_report = analyze_calibration(sim_base, conditions, truth)
    print(cal_report)

    # 3d: Locus-level gDNA (tool-internal aggregates).
    locus_report = analyze_locus_gdna(sim_base, conditions)
    print(locus_report)

    # 3e: Secondary — gross (hard-label) confusion table.
    if not args.skip_frag_analysis:
        frag_report = analyze_fragment_assignment(sim_base, conditions, assignment_rows)
        print(frag_report)
    else:
        frag_report = "\n  [SKIP] Fragment analysis skipped (--skip-frag-analysis)"
        print(frag_report)

    # 3f: Cheap acceptance checks.
    acceptance_report = analyze_postfix_acceptance(sim_base, conditions, assignment_rows)
    print(acceptance_report)

    # ── Save full report ──
    report_path = sim_base / "analysis_report.txt"
    full_report = "\n".join(
        [
            "=" * 100,
            "  RIGEL SYNTHETIC SIMULATION — FULL ANALYSIS REPORT",
            "=" * 100,
            net_flow_report,
            abundance_report,
            cal_report,
            locus_report,
            frag_report,
            acceptance_report,
        ]
    )
    report_path.write_text(full_report)
    print(f"\n\n  Full report saved to: {report_path}")
    print("  Per-condition / per-locus / per-transcript TSVs written alongside it.")
    print("=" * 100)


if __name__ == "__main__":
    main()
