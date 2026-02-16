#!/usr/bin/env python3
"""
Isoform stress test: hulkrna vs salmon on a complex multi-isoform gene.

Tests both tools on salmon's home turf:
  - strand_specificity = 1.0 (perfect, salmon's assumption)
  - gDNA = 0 (no gDNA contamination, also salmon's assumption)

Scenario:
  - 20 kb genome
  - 1 gene on + strand with 11 exons (varying sizes 300–1000 bp)
  - 10 isoforms using different subsets of those 11 exons
  - 50,000 total fragments

Two abundance distributions:
  1. Dominant isoform: one isoform at ~95% of counts, others share ~5%
  2. Equal abundance: all 10 isoforms at equal relative abundance

Usage:
    conda run -n hulkrna python scripts/isoform_stress_test.py
"""

import json
import logging
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import (
    Scenario,
    ScenarioResult,
    SimConfig,
    reverse_complement,
    run_benchmark,
)
from hulkrna.types import Strand

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


# =====================================================================
# Constants
# =====================================================================

N_FRAGMENTS = 50_000
SIM_SEED = 42
PIPELINE_SEED = 42
SALMON_THREADS = 1

# =====================================================================
# Exon layout: 11 exons within 20 kb genome
# =====================================================================
#
# Exon sizes vary: roughly 300–1000 bp each.
# Introns are ~300–700 bp between exons.
# All exons fit within [500, 18500] of the 20 kb genome.
#
# Exon  Size   Interval
# E1    500    [500,  1000]
# E2    300    [1300, 1600]
# E3    800    [1900, 2700]
# E4    400    [3100, 3500]
# E5    1000   [3900, 4900]
# E6    300    [5300, 5600]
# E7    700    [5900, 6600]
# E8    500    [7000, 7500]
# E9    300    [7800, 8100]
# E10   600    [8500, 9100]
# E11   400    [9500, 9900]
#
# Total exonic span: 5800 bp across 11 exons, gene span 9400 bp.
#

EXONS = [
    (500, 1000),    # E1  500 bp
    (1300, 1600),   # E2  300 bp
    (1900, 2700),   # E3  800 bp
    (3100, 3500),   # E4  400 bp
    (3900, 4900),   # E5  1000 bp
    (5300, 5600),   # E6  300 bp
    (5900, 6600),   # E7  700 bp
    (7000, 7500),   # E8  500 bp
    (7800, 8100),   # E9  300 bp
    (8500, 9100),   # E10 600 bp
    (9500, 9900),   # E11 400 bp
]

EXON_LENGTHS = [e - s for s, e in EXONS]

# =====================================================================
# 10 isoforms using subsets of the 11 exons
# =====================================================================
#
# Each isoform is defined by its exon indices (0-based).
# Design ensures a mix of:
#   - Full-length (all exons)
#   - Long isoforms (8-10 exons)
#   - Medium isoforms (5-7 exons)
#   - Short isoforms (3-4 exons)
# And importantly: many pairs share most exons, testing EM resolution.
#

ISOFORM_EXON_INDICES = {
    "iso01": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],  # full length (5800 bp)
    "iso02": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],       # skip E11 (5400 bp)
    "iso03": [0, 2, 3, 4, 5, 6, 7, 8, 9, 10],       # skip E2 (5500 bp)
    "iso04": [0, 1, 3, 4, 6, 7, 9, 10],             # skip E3,E6,E9 (4200 bp)
    "iso05": [0, 2, 4, 6, 8, 10],                   # odd exons (3000 bp)
    "iso06": [1, 3, 5, 7, 9],                       # even exons (2100 bp)
    "iso07": [0, 1, 2, 3, 4],                       # first 5 exons (3000 bp)
    "iso08": [6, 7, 8, 9, 10],                      # last 5 exons (2500 bp)
    "iso09": [2, 3, 4, 5, 6],                       # middle 5 exons (3200 bp)
    "iso10": [0, 4, 10],                            # widely separated (1900 bp)
}


def compute_isoform_lengths() -> dict[str, int]:
    """Compute exonic length of each isoform."""
    return {
        iso_id: sum(EXON_LENGTHS[i] for i in exon_idxs)
        for iso_id, exon_idxs in ISOFORM_EXON_INDICES.items()
    }


def describe_isoforms():
    """Print isoform structure summary."""
    lengths = compute_isoform_lengths()
    print("\n=== Isoform Structure ===")
    print(f"{'Isoform':<8} {'Exons':<35} {'Length':>6}  Exon details")
    print("-" * 85)
    for iso_id, exon_idxs in ISOFORM_EXON_INDICES.items():
        exon_names = ", ".join(f"E{i+1}" for i in exon_idxs)
        exon_sizes = "+".join(str(EXON_LENGTHS[i]) for i in exon_idxs)
        print(f"{iso_id:<8} {exon_names:<35} {lengths[iso_id]:>5}bp  ({exon_sizes})")
    print()

    # Show exon sharing between isoforms
    print("=== Exon Uniqueness ===")
    for i, (s, e) in enumerate(EXONS):
        users = [iso for iso, idxs in ISOFORM_EXON_INDICES.items() if i in idxs]
        print(f"  E{i+1:>2} ({EXON_LENGTHS[i]:>4}bp): used by {len(users)}/10 isoforms: {', '.join(users)}")
    print()


# =====================================================================
# Scenario builder
# =====================================================================


def make_isoform_scenario(
    abundances: dict[str, float],
    work_dir: Path,
) -> Scenario:
    """Build a Scenario with 10 isoforms at specified abundances."""
    sc = Scenario(
        "isoform_stress",
        genome_length=20_000,
        seed=SIM_SEED,
        work_dir=work_dir,
    )
    transcripts = []
    for iso_id, exon_idxs in ISOFORM_EXON_INDICES.items():
        exons = [EXONS[i] for i in exon_idxs]
        transcripts.append({
            "t_id": iso_id,
            "exons": exons,
            "abundance": abundances[iso_id],
        })
    sc.add_gene("g1", "+", transcripts)
    return sc


# =====================================================================
# Salmon runner
# =====================================================================


def write_transcript_fasta(result: ScenarioResult, output_path: Path) -> None:
    """Write transcript sequences as FASTA for salmon indexing."""
    with open(output_path, "w") as f:
        for t in result.transcripts:
            exon_seqs = [result.genome[e.start:e.end] for e in t.exons]
            mrna_seq = "".join(exon_seqs)
            if t.strand == Strand.NEG:
                mrna_seq = reverse_complement(mrna_seq)
            f.write(f">{t.t_id}\n{mrna_seq}\n")


def run_salmon(
    result: ScenarioResult,
    work_dir: Path,
) -> dict[str, float]:
    """Run salmon index + quant and return per-transcript counts."""
    work_dir.mkdir(parents=True, exist_ok=True)

    tx_fasta = work_dir / "transcripts.fa"
    write_transcript_fasta(result, tx_fasta)

    idx_dir = work_dir / "salmon_idx"
    cmd_index = [
        "salmon", "index",
        "-t", str(tx_fasta),
        "-i", str(idx_dir),
        "-k", "15",
    ]
    proc = subprocess.run(cmd_index, capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        raise RuntimeError(f"salmon index failed: {proc.stderr[:500]}")

    quant_dir = work_dir / "salmon_quant"
    cmd_quant = [
        "salmon", "quant",
        "-i", str(idx_dir),
        "-l", "A",
        "-1", str(result.fastq_r1),
        "-2", str(result.fastq_r2),
        "-o", str(quant_dir),
        "-p", str(SALMON_THREADS),
        "--validateMappings",
    ]
    proc = subprocess.run(cmd_quant, capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        raise RuntimeError(f"salmon quant failed: {proc.stderr[:500]}")

    quant_sf = quant_dir / "quant.sf"
    df = pd.read_csv(quant_sf, sep="\t")
    return dict(zip(df["Name"], df["NumReads"]))


# =====================================================================
# Core runner
# =====================================================================


def run_scenario(
    label: str,
    abundances: dict[str, float],
    tmp_dir: Path,
    dir_name: str = "",
) -> dict:
    """Run hulkrna and salmon on one abundance configuration."""
    logger.info("=" * 70)
    logger.info(f"Running scenario: {label}")
    logger.info("=" * 70)

    safe_name = dir_name or "".join(c if c.isalnum() or c in "_-" else "_" for c in label)
    sc = make_isoform_scenario(abundances, tmp_dir / safe_name)
    try:
        sim_config = SimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=SIM_SEED,
        )
        result = sc.build(n_fragments=N_FRAGMENTS, sim_config=sim_config)
        ground_truth = result.ground_truth_from_fastq()

        # --- hulkrna ---
        t0 = time.monotonic()
        pipeline_result = run_pipeline(
            result.bam_path, result.index,
            sj_strand_tag="ts",
            seed=PIPELINE_SEED,
        )
        hulk_elapsed = time.monotonic() - t0

        bench = run_benchmark(result, pipeline_result, scenario_name=label)
        hulk_per_t = {}
        for ta in bench.transcripts:
            hulk_per_t[ta.t_id] = float(ta.observed)

        # --- salmon ---
        t0 = time.monotonic()
        salmon_counts = run_salmon(result, tmp_dir / f"{safe_name}_salmon")
        salmon_elapsed = time.monotonic() - t0

        # --- Build result ---
        isoform_lengths = compute_isoform_lengths()
        result_data = {
            "label": label,
            "n_fragments": N_FRAGMENTS,
            "abundances": abundances,
            "hulk_elapsed": round(hulk_elapsed, 3),
            "salmon_elapsed": round(salmon_elapsed, 3),
            "transcripts": [],
        }

        total_hulk_mae = 0.0
        total_salmon_mae = 0.0

        for iso_id in sorted(ISOFORM_EXON_INDICES.keys()):
            truth = ground_truth.get(iso_id, 0)
            hulk = hulk_per_t.get(iso_id, 0.0)
            sal = salmon_counts.get(iso_id, 0.0)
            hulk_err = abs(hulk - truth)
            sal_err = abs(sal - truth)
            total_hulk_mae += hulk_err
            total_salmon_mae += sal_err

            result_data["transcripts"].append({
                "t_id": iso_id,
                "length": isoform_lengths[iso_id],
                "abundance": abundances[iso_id],
                "truth": truth,
                "hulkrna": round(hulk, 1),
                "salmon": round(sal, 1),
                "hulk_err": round(hulk_err, 1),
                "salmon_err": round(sal_err, 1),
            })

        n_iso = len(ISOFORM_EXON_INDICES)
        result_data["hulk_total_mae"] = round(total_hulk_mae, 1)
        result_data["salmon_total_mae"] = round(total_salmon_mae, 1)
        result_data["hulk_mean_mae"] = round(total_hulk_mae / n_iso, 1)
        result_data["salmon_mean_mae"] = round(total_salmon_mae / n_iso, 1)

        return result_data
    finally:
        sc.cleanup()


# =====================================================================
# Report
# =====================================================================


def print_result_table(data: dict):
    """Print formatted comparison table for one scenario."""
    print(f"\n{'=' * 90}")
    print(f"  Scenario: {data['label']}")
    print(f"  Fragments: {data['n_fragments']:,}")
    print(f"  hulkrna elapsed: {data['hulk_elapsed']:.1f}s | salmon elapsed: {data['salmon_elapsed']:.1f}s")
    print(f"{'=' * 90}")

    print(f"\n{'Isoform':<8} {'Len':>5} {'Abund':>7} {'Truth':>7} {'hulkrna':>8} {'salmon':>8} {'H_err':>7} {'S_err':>7} {'Winner':>8}")
    print("-" * 80)

    for t in data["transcripts"]:
        h_err = t["hulk_err"]
        s_err = t["salmon_err"]
        if abs(h_err - s_err) < 0.5:
            winner = "tie"
        elif h_err < s_err:
            winner = "hulkrna"
        else:
            winner = "salmon"
        print(
            f"{t['t_id']:<8} {t['length']:>5} {t['abundance']:>7.1f} "
            f"{t['truth']:>7} {t['hulkrna']:>8.1f} {t['salmon']:>8.1f} "
            f"{h_err:>7.1f} {s_err:>7.1f} {winner:>8}"
        )

    print("-" * 80)
    print(
        f"{'TOTAL':<8} {'':>5} {'':>7} {'':>7} {'':>8} {'':>8} "
        f"{data['hulk_total_mae']:>7.1f} {data['salmon_total_mae']:>7.1f} "
        f"{'hulkrna' if data['hulk_total_mae'] < data['salmon_total_mae'] else 'salmon':>8}"
    )
    print(
        f"{'MEAN':<8} {'':>5} {'':>7} {'':>7} {'':>8} {'':>8} "
        f"{data['hulk_mean_mae']:>7.1f} {data['salmon_mean_mae']:>7.1f}"
    )


def write_report(results: list[dict], output_path: Path):
    """Write markdown report."""
    isoform_lengths = compute_isoform_lengths()

    with open(output_path, "w") as f:
        f.write("# Isoform Stress Test: hulkrna vs salmon\n\n")
        f.write(f"**Fragments**: {N_FRAGMENTS:,}\n")
        f.write(f"**Strand specificity**: 1.0 (perfect)\n")
        f.write(f"**gDNA**: 0 (none)\n")
        f.write(f"**Genome**: 20 kb, 1 gene, 11 exons, 10 isoforms\n\n")

        # Isoform structure table
        f.write("## Isoform Structure\n\n")
        f.write("| Isoform | Exons | Length (bp) |\n")
        f.write("| --- | --- | --- |\n")
        for iso_id, exon_idxs in ISOFORM_EXON_INDICES.items():
            exon_names = ", ".join(f"E{i+1}" for i in exon_idxs)
            f.write(f"| {iso_id} | {exon_names} | {isoform_lengths[iso_id]} |\n")
        f.write("\n")

        # Exon uniqueness
        f.write("## Exon Usage\n\n")
        f.write("| Exon | Size (bp) | Used by (N isoforms) |\n")
        f.write("| --- | --- | --- |\n")
        for i, (s, e) in enumerate(EXONS):
            users = [iso for iso, idxs in ISOFORM_EXON_INDICES.items() if i in idxs]
            f.write(f"| E{i+1} | {EXON_LENGTHS[i]} | {len(users)} |\n")
        f.write("\n")

        # Results for each scenario
        for data in results:
            f.write(f"## {data['label']}\n\n")
            f.write(f"hulkrna: {data['hulk_elapsed']:.1f}s | salmon: {data['salmon_elapsed']:.1f}s\n\n")

            f.write("| Isoform | Length | Abundance | Truth | hulkrna | salmon | hulk err | salmon err |\n")
            f.write("| --- | --- | --- | --- | --- | --- | --- | --- |\n")
            for t in data["transcripts"]:
                f.write(
                    f"| {t['t_id']} | {t['length']} | {t['abundance']:.1f} "
                    f"| {t['truth']} | {t['hulkrna']:.1f} | {t['salmon']:.1f} "
                    f"| {t['hulk_err']:.1f} | {t['salmon_err']:.1f} |\n"
                )

            f.write(f"\n**Total MAE**: hulkrna={data['hulk_total_mae']:.1f}, salmon={data['salmon_total_mae']:.1f}\n")
            f.write(f"**Mean MAE**: hulkrna={data['hulk_mean_mae']:.1f}, salmon={data['salmon_mean_mae']:.1f}\n\n")

        # Summary
        f.write("## Summary\n\n")
        f.write("| Scenario | hulkrna Total MAE | salmon Total MAE | Winner |\n")
        f.write("| --- | --- | --- | --- |\n")
        for data in results:
            winner = "hulkrna" if data["hulk_total_mae"] < data["salmon_total_mae"] else "salmon"
            if abs(data["hulk_total_mae"] - data["salmon_total_mae"]) < 1.0:
                winner = "tie"
            f.write(
                f"| {data['label']} | {data['hulk_total_mae']:.1f} "
                f"| {data['salmon_total_mae']:.1f} | {winner} |\n"
            )


# =====================================================================
# Main
# =====================================================================


def main():
    results = []

    with tempfile.TemporaryDirectory(prefix="iso_stress_") as tmpdir:
        tmp = Path(tmpdir)

        # Scenarios 1-10: Each isoform as dominant (95%, others share 5%)
        isoform_ids = sorted(ISOFORM_EXON_INDICES.keys())
        for dominant_iso in isoform_ids:
            dominant_abundances = {}
            for iso_id in ISOFORM_EXON_INDICES:
                if iso_id == dominant_iso:
                    dominant_abundances[iso_id] = 95.0
                else:
                    dominant_abundances[iso_id] = 5.0 / 9  # ~0.556 each

            label = f"Dominant {dominant_iso} (95%)"
            dir_name = f"dominant_{dominant_iso}"
            result = run_scenario(label, dominant_abundances, tmp, dir_name=dir_name)
            results.append(result)
            print_result_table(result)

        # Scenario 11: Equal abundance (all isoforms equal)
        equal_abundances = {iso_id: 10.0 for iso_id in ISOFORM_EXON_INDICES}

        result_eq = run_scenario("Equal abundance (all 10%)", equal_abundances, tmp)
        results.append(result_eq)
        print_result_table(result_eq)

    # Write reports
    output_dir = Path("docs")
    output_dir.mkdir(exist_ok=True)

    json_path = output_dir / "isoform_stress_results.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"JSON results saved to {json_path}")

    md_path = output_dir / "isoform_stress_test.md"
    write_report(results, md_path)
    logger.info(f"Report saved to {md_path}")


if __name__ == "__main__":
    main()
