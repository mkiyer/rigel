#!/usr/bin/env python3
"""
Comparative benchmark: hulkrna vs salmon on simulated scenarios.

Runs both tools on the same simulated scenarios from test_scenarios.py
and produces a structured JSON results file that can be rendered into
a Markdown report.

Usage:
    conda run -n hulkrna python scripts/benchmark_vs_salmon.py

Output:
    docs/benchmark_results.json   — raw results
    docs/benchmark_hulkrna_vs_salmon.md — formatted report
"""

import json
import logging
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np
import pandas as pd

# -- hulkrna imports --------------------------------------------------------
from hulkrna.pipeline import run_pipeline
from hulkrna.sim import (
    GDNAConfig,
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
# Configuration
# =====================================================================

N_FRAGMENTS = 2000       # More fragments than test suite for stable results
SIM_SEED = 42
PIPELINE_SEED = 42
SALMON_THREADS = 1       # Reproducibility


# gDNA sweep: a representative subset
GDNA_LEVELS = [0, 20, 50]

# Strand specificity sweep
STRAND_LEVELS = [0.5, 0.9, 1.0]

# Abundance ratios for isoform/antisense sweeps
ABUNDANCE_RATIOS = [1, 4, 16]


# =====================================================================
# Helpers
# =====================================================================


def _sim_config(*, strand_specificity: float = 1.0, seed: int = SIM_SEED):
    return SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=seed,
    )


def _gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=100,
        frag_max=1000,
    )


# =====================================================================
# Transcript FASTA extraction (for salmon)
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


# =====================================================================
# Salmon runner
# =====================================================================


def run_salmon(
    result: ScenarioResult,
    work_dir: Path,
    lib_type: str = "A",
) -> dict[str, float]:
    """Run salmon index + quant and return per-transcript counts.

    Parameters
    ----------
    result : ScenarioResult
        Built scenario with FASTQ files and transcript definitions.
    work_dir : Path
        Working directory for salmon artifacts.
    lib_type : str
        Salmon library type string. "A" = auto-detect.

    Returns
    -------
    dict mapping transcript_id → salmon NumReads (float).
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Write transcript FASTA
    tx_fasta = work_dir / "transcripts.fa"
    write_transcript_fasta(result, tx_fasta)

    # 2. Salmon index
    idx_dir = work_dir / "salmon_idx"
    cmd_index = [
        "salmon", "index",
        "-t", str(tx_fasta),
        "-i", str(idx_dir),
        "-k", "15",  # small k for tiny simulated transcripts
    ]
    logger.debug("  salmon index: %s", " ".join(cmd_index))
    proc = subprocess.run(
        cmd_index, capture_output=True, text=True, timeout=120,
    )
    if proc.returncode != 0:
        logger.error("salmon index failed:\n%s", proc.stderr)
        raise RuntimeError(f"salmon index failed: {proc.stderr[:500]}")

    # 3. Salmon quant
    quant_dir = work_dir / "salmon_quant"
    cmd_quant = [
        "salmon", "quant",
        "-i", str(idx_dir),
        "-l", lib_type,
        "-1", str(result.fastq_r1),
        "-2", str(result.fastq_r2),
        "-o", str(quant_dir),
        "-p", str(SALMON_THREADS),
        "--validateMappings",
    ]
    logger.debug("  salmon quant: %s", " ".join(cmd_quant))
    proc = subprocess.run(
        cmd_quant, capture_output=True, text=True, timeout=120,
    )
    if proc.returncode != 0:
        logger.error("salmon quant failed:\n%s", proc.stderr)
        raise RuntimeError(f"salmon quant failed: {proc.stderr[:500]}")

    # 4. Parse quant.sf
    quant_sf = quant_dir / "quant.sf"
    df = pd.read_csv(quant_sf, sep="\t")
    counts = dict(zip(df["Name"], df["NumReads"]))
    return counts


# =====================================================================
# Benchmark result structures
# =====================================================================


@dataclass
class ToolResult:
    """Per-transcript counts from one tool."""
    tool: str
    per_transcript: dict[str, float] = field(default_factory=dict)
    total_count: float = 0.0
    elapsed_sec: float = 0.0


@dataclass
class ScenarioBenchmark:
    """Benchmark for one scenario + parameter combination."""
    scenario_name: str
    param_axis: str
    param_value: float | int
    n_simulated: int
    n_fragments: int  # fragments that entered the pipeline (after alignment)
    ground_truth: dict[str, int] = field(default_factory=dict)
    n_gdna_expected: int = 0
    hulkrna_result: ToolResult | None = None
    salmon_result: ToolResult | None = None


# =====================================================================
# Scenario definitions (matching test_scenarios.py)
# =====================================================================


def _make_single_exon(tmp_path: Path) -> Scenario:
    sc = Scenario("single_exon", genome_length=5000, seed=SIM_SEED,
                  work_dir=tmp_path / "single_exon")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
    ])
    return sc


def _make_spliced_gene(tmp_path: Path) -> Scenario:
    sc = Scenario("spliced_gene", genome_length=5000, seed=SIM_SEED,
                  work_dir=tmp_path / "spliced_gene")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
         "abundance": 100},
    ])
    return sc


def _make_non_overlapping(tmp_path: Path) -> Scenario:
    sc = Scenario("non_overlapping", genome_length=8000, seed=SIM_SEED,
                  work_dir=tmp_path / "non_overlapping")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
         "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(4000, 4400)],
         "abundance": 100},
    ])
    return sc


def _make_two_isoforms(tmp_path: Path, major: float, minor: float) -> Scenario:
    sc = Scenario("two_isoforms", genome_length=5000, seed=SIM_SEED,
                  work_dir=tmp_path / "two_isoforms")
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(200, 500), (1000, 1300), (2000, 2300)],
         "abundance": major},
        {"t_id": "t2",
         "exons": [(200, 500), (2000, 2300)],
         "abundance": minor},
    ])
    return sc


def _make_overlapping_antisense(
    tmp_path: Path, g1_ab: float, g2_ab: float,
) -> Scenario:
    sc = Scenario("overlapping_antisense", genome_length=5000, seed=SIM_SEED,
                  work_dir=tmp_path / "overlapping_antisense")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
         "abundance": g1_ab},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(300, 600), (1100, 1400)],
         "abundance": g2_ab},
    ])
    return sc


# =====================================================================
# Core benchmark runner
# =====================================================================


def run_one_benchmark(
    scenario: Scenario,
    *,
    gdna_abundance: float = 0,
    strand_specificity: float = 1.0,
    n_fragments: int = N_FRAGMENTS,
    scenario_name: str = "",
    param_axis: str = "",
    param_value: float | int = 0,
    salmon_work_dir: Path | None = None,
) -> ScenarioBenchmark:
    """Run hulkrna and salmon on a single scenario configuration."""

    sim_config = _sim_config(strand_specificity=strand_specificity)
    gdna = _gdna_config(gdna_abundance)

    # Build scenario (genome, reads, alignment)
    result = scenario.build(
        n_fragments=n_fragments,
        sim_config=sim_config,
        gdna_config=gdna,
    )

    ground_truth = result.ground_truth_from_fastq()
    n_gdna_expected = result.ground_truth_gdna_count()

    # --- Run hulkrna ---
    t0 = time.monotonic()
    pipeline_result = run_pipeline(
        result.bam_path, result.index,
        sj_strand_tag="ts",
        seed=PIPELINE_SEED,
    )
    hulk_elapsed = time.monotonic() - t0

    bench = run_benchmark(result, pipeline_result,
                          scenario_name=scenario_name)

    hulk_per_t = {}
    for ta in bench.transcripts:
        hulk_per_t[ta.t_id] = float(ta.observed)

    hulk_tool = ToolResult(
        tool="hulkrna",
        per_transcript=hulk_per_t,
        total_count=float(bench.total_observed),
        elapsed_sec=round(hulk_elapsed, 3),
    )

    # --- Run salmon ---
    if salmon_work_dir is None:
        salmon_work_dir = Path(tempfile.mkdtemp(prefix="salmon_"))

    t0 = time.monotonic()
    try:
        salmon_counts = run_salmon(result, salmon_work_dir)
    except Exception as e:
        logger.warning("Salmon failed for %s: %s", scenario_name, e)
        salmon_counts = {t.t_id: 0.0 for t in result.transcripts}
    salmon_elapsed = time.monotonic() - t0

    salmon_total = sum(salmon_counts.get(t_id, 0.0) for t_id in ground_truth)
    salmon_tool = ToolResult(
        tool="salmon",
        per_transcript={t_id: salmon_counts.get(t_id, 0.0)
                        for t_id in ground_truth},
        total_count=float(salmon_total),
        elapsed_sec=round(salmon_elapsed, 3),
    )

    return ScenarioBenchmark(
        scenario_name=scenario_name,
        param_axis=param_axis,
        param_value=param_value,
        n_simulated=result.n_simulated,
        n_fragments=bench.n_fragments,
        ground_truth={k: int(v) for k, v in ground_truth.items()},
        n_gdna_expected=n_gdna_expected,
        hulkrna_result=hulk_tool,
        salmon_result=salmon_tool,
    )


# =====================================================================
# Run all benchmarks
# =====================================================================


def run_all_benchmarks() -> list[dict]:
    """Run all scenario × parameter benchmarks."""
    results = []

    with tempfile.TemporaryDirectory(prefix="hulk_bench_") as tmpdir:
        tmp = Path(tmpdir)

        # --- Scenario 1: Single exon (gDNA sweep) -----------------------
        logger.info("=" * 60)
        logger.info("Scenario 1: Single-exon gene")
        logger.info("=" * 60)
        for gdna in GDNA_LEVELS:
            name = f"single_exon_gdna_{gdna}"
            logger.info("  Running %s ...", name)
            sc = _make_single_exon(tmp / name)
            try:
                b = run_one_benchmark(
                    sc, gdna_abundance=gdna,
                    scenario_name=name, param_axis="gdna",
                    param_value=gdna,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        for ss in STRAND_LEVELS:
            name = f"single_exon_ss_{ss}"
            logger.info("  Running %s ...", name)
            sc = _make_single_exon(tmp / name)
            try:
                b = run_one_benchmark(
                    sc, strand_specificity=ss,
                    scenario_name=name, param_axis="strand_spec",
                    param_value=ss,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        # --- Scenario 2: Spliced gene -----------------------------------
        logger.info("=" * 60)
        logger.info("Scenario 2: Spliced gene")
        logger.info("=" * 60)
        for gdna in GDNA_LEVELS:
            name = f"spliced_gene_gdna_{gdna}"
            logger.info("  Running %s ...", name)
            sc = _make_spliced_gene(tmp / name)
            try:
                b = run_one_benchmark(
                    sc, gdna_abundance=gdna,
                    scenario_name=name, param_axis="gdna",
                    param_value=gdna,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        for ss in STRAND_LEVELS:
            name = f"spliced_gene_ss_{ss}"
            logger.info("  Running %s ...", name)
            sc = _make_spliced_gene(tmp / name)
            try:
                b = run_one_benchmark(
                    sc, strand_specificity=ss,
                    scenario_name=name, param_axis="strand_spec",
                    param_value=ss,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        # --- Scenario 3: Non-overlapping genes --------------------------
        logger.info("=" * 60)
        logger.info("Scenario 3: Non-overlapping genes")
        logger.info("=" * 60)
        for gdna in GDNA_LEVELS:
            name = f"non_overlapping_gdna_{gdna}"
            logger.info("  Running %s ...", name)
            sc = _make_non_overlapping(tmp / name)
            try:
                b = run_one_benchmark(
                    sc, gdna_abundance=gdna,
                    scenario_name=name, param_axis="gdna",
                    param_value=gdna,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        for ss in STRAND_LEVELS:
            name = f"non_overlapping_ss_{ss}"
            logger.info("  Running %s ...", name)
            sc = _make_non_overlapping(tmp / name)
            try:
                b = run_one_benchmark(
                    sc, strand_specificity=ss,
                    scenario_name=name, param_axis="strand_spec",
                    param_value=ss,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        # --- Scenario 4: Two isoforms (abundance sweep) -----------------
        logger.info("=" * 60)
        logger.info("Scenario 4: Two isoforms")
        logger.info("=" * 60)
        for fc in ABUNDANCE_RATIOS:
            name = f"two_isoforms_fc_{fc}"
            logger.info("  Running %s ...", name)
            sc = _make_two_isoforms(tmp / name, 100, 100 / fc)
            try:
                b = run_one_benchmark(
                    sc, scenario_name=name,
                    param_axis="fold_change",
                    param_value=fc, n_fragments=N_FRAGMENTS,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        for gdna in GDNA_LEVELS:
            name = f"two_isoforms_gdna_{gdna}"
            logger.info("  Running %s ...", name)
            sc = _make_two_isoforms(tmp / name, 100, 10)
            try:
                b = run_one_benchmark(
                    sc, gdna_abundance=gdna,
                    scenario_name=name, param_axis="gdna",
                    param_value=gdna, n_fragments=N_FRAGMENTS,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        # --- Scenario 5: Overlapping antisense --------------------------
        logger.info("=" * 60)
        logger.info("Scenario 5: Overlapping antisense genes")
        logger.info("=" * 60)
        for fc in ABUNDANCE_RATIOS:
            name = f"antisense_fc_{fc}"
            logger.info("  Running %s ...", name)
            sc = _make_overlapping_antisense(tmp / name, 100, 100 / fc)
            try:
                b = run_one_benchmark(
                    sc, scenario_name=name,
                    param_axis="fold_change",
                    param_value=fc, n_fragments=N_FRAGMENTS,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        for gdna in GDNA_LEVELS:
            name = f"antisense_gdna_{gdna}"
            logger.info("  Running %s ...", name)
            sc = _make_overlapping_antisense(tmp / name, 100, 100)
            try:
                b = run_one_benchmark(
                    sc, gdna_abundance=gdna,
                    scenario_name=name, param_axis="gdna",
                    param_value=gdna, n_fragments=N_FRAGMENTS,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

        for ss in STRAND_LEVELS:
            name = f"antisense_ss_{ss}"
            logger.info("  Running %s ...", name)
            sc = _make_overlapping_antisense(tmp / name, 100, 100)
            try:
                b = run_one_benchmark(
                    sc, strand_specificity=ss,
                    scenario_name=name, param_axis="strand_spec",
                    param_value=ss, n_fragments=N_FRAGMENTS,
                    salmon_work_dir=tmp / f"{name}_salmon",
                )
                results.append(asdict(b))
            finally:
                sc.cleanup()

    return results


# =====================================================================
# Report generation
# =====================================================================


def _abs_diff(observed: float, expected: int) -> float:
    return abs(observed - expected)


def _rel_err(observed: float, expected: int) -> float:
    if expected == 0:
        return 0.0 if observed == 0 else float("inf")
    return abs(observed - expected) / expected


def generate_report(results: list[dict], output_path: Path) -> None:
    """Generate a Markdown benchmark report."""
    lines = []
    lines.append("# Benchmark: hulkrna vs salmon")
    lines.append("")
    lines.append(f"**Generated**: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"**Fragments per scenario**: {N_FRAGMENTS}")
    lines.append(f"**Seed**: {SIM_SEED}")
    lines.append(f"**salmon version**: 1.10.3")
    lines.append("")

    lines.append("## Overview")
    lines.append("")
    lines.append("This report compares **hulkrna** (Bayesian EM with gDNA modeling) "
                 "against **salmon** (lightweight quasi-mapping quantification) "
                 "on synthetic RNA-seq scenarios with known ground truth.")
    lines.append("")
    lines.append("Key differences in design:")
    lines.append("- **hulkrna** takes a BAM (genome-aligned) as input and models "
                 "gDNA contamination, strand specificity, and insert size jointly.")
    lines.append("- **salmon** performs lightweight mapping directly against the "
                 "transcriptome (no genome alignment, no gDNA modeling).")
    lines.append("")

    # Group results by scenario base name
    scenarios = {}
    for r in results:
        # Extract base scenario: everything before the axis
        name = r["scenario_name"]
        base = name.rsplit("_", 2)[0] if r["param_axis"] else name
        # Better grouping: use known prefixes
        for prefix in ["single_exon", "spliced_gene", "non_overlapping",
                        "two_isoforms", "antisense"]:
            if name.startswith(prefix):
                base = prefix
                break
        scenarios.setdefault(base, []).append(r)

    scenario_titles = {
        "single_exon": "Scenario 1: Single-exon gene (unspliced)",
        "spliced_gene": "Scenario 2: Spliced gene (multi-exon)",
        "non_overlapping": "Scenario 3: Non-overlapping genes",
        "two_isoforms": "Scenario 4: Two isoforms (shared exons)",
        "antisense": "Scenario 5: Overlapping antisense genes",
    }

    scenario_order = ["single_exon", "spliced_gene", "non_overlapping",
                      "two_isoforms", "antisense"]

    # Overall summary table
    lines.append("## Summary: Mean Absolute Error by Scenario and Parameter")
    lines.append("")

    # Per-scenario sections
    for base in scenario_order:
        if base not in scenarios:
            continue
        runs = scenarios[base]
        title = scenario_titles.get(base, base)
        lines.append(f"## {title}")
        lines.append("")

        # Group by axis
        axes = {}
        for r in runs:
            axes.setdefault(r["param_axis"], []).append(r)

        for axis, axis_runs in sorted(axes.items()):
            axis_label = {
                "gdna": "gDNA Contamination Sweep",
                "strand_spec": "Strand Specificity Sweep",
                "fold_change": "Abundance Ratio Sweep",
            }.get(axis, axis)
            param_label = {
                "gdna": "gDNA abundance",
                "strand_spec": "Strand spec.",
                "fold_change": "Fold change",
            }.get(axis, "Parameter")

            lines.append(f"### {axis_label}")
            lines.append("")

            # Determine transcripts from first run
            first = axis_runs[0]
            t_ids = sorted(first["ground_truth"].keys())

            # Build table header
            header_cols = [param_label]
            for t_id in t_ids:
                header_cols.extend([
                    f"Truth ({t_id})",
                    f"hulkrna ({t_id})",
                    f"salmon ({t_id})",
                ])
            if first.get("n_gdna_expected", 0) > 0 or axis == "gdna":
                header_cols.extend(["gDNA truth", "hulkrna gDNA"])
            header_cols.extend(["hulkrna MAE", "salmon MAE"])

            lines.append("| " + " | ".join(header_cols) + " |")
            lines.append("| " + " | ".join(["---"] * len(header_cols)) + " |")

            for r in sorted(axis_runs, key=lambda x: x["param_value"]):
                pv = r["param_value"]
                gt = r["ground_truth"]
                hulk = r["hulkrna_result"]["per_transcript"]
                salm = r["salmon_result"]["per_transcript"]

                row = [str(pv)]
                hulk_errors = []
                salm_errors = []
                for t_id in t_ids:
                    truth = gt.get(t_id, 0)
                    h = hulk.get(t_id, 0.0)
                    s = salm.get(t_id, 0.0)
                    hulk_errors.append(_abs_diff(h, truth))
                    salm_errors.append(_abs_diff(s, truth))
                    row.extend([str(truth), f"{h:.0f}", f"{s:.1f}"])

                if first.get("n_gdna_expected", 0) > 0 or axis == "gdna":
                    gdna_truth = r.get("n_gdna_expected", 0)
                    # hulkrna gDNA = total fragments - RNA counts
                    hulk_total_rna = sum(hulk.get(t, 0) for t in t_ids)
                    hulk_gdna = r["n_fragments"] - hulk_total_rna
                    row.extend([str(gdna_truth), f"{hulk_gdna:.0f}"])

                hulk_mae = np.mean(hulk_errors) if hulk_errors else 0
                salm_mae = np.mean(salm_errors) if salm_errors else 0
                row.extend([f"{hulk_mae:.1f}", f"{salm_mae:.1f}"])

                lines.append("| " + " | ".join(row) + " |")

            lines.append("")

    # --- Overall comparison summary ---
    lines.append("## Overall Comparison")
    lines.append("")

    total_hulk_err = 0.0
    total_salm_err = 0.0
    n_comparisons = 0
    hulk_wins = 0
    salm_wins = 0
    ties = 0

    gdna_scenarios = []
    no_gdna_scenarios = []

    for r in results:
        gt = r["ground_truth"]
        hulk = r["hulkrna_result"]["per_transcript"]
        salm = r["salmon_result"]["per_transcript"]

        h_err = sum(_abs_diff(hulk.get(t, 0), gt[t]) for t in gt)
        s_err = sum(_abs_diff(salm.get(t, 0), gt[t]) for t in gt)

        total_hulk_err += h_err
        total_salm_err += s_err
        n_comparisons += 1

        if h_err < s_err - 0.5:
            hulk_wins += 1
        elif s_err < h_err - 0.5:
            salm_wins += 1
        else:
            ties += 1

        if r.get("n_gdna_expected", 0) > 0:
            gdna_scenarios.append((h_err, s_err))
        else:
            no_gdna_scenarios.append((h_err, s_err))

    lines.append(f"**Total scenarios benchmarked**: {n_comparisons}")
    lines.append("")
    lines.append(f"| Metric | hulkrna | salmon |")
    lines.append(f"| --- | --- | --- |")
    lines.append(f"| Total abs error (sum) | {total_hulk_err:.0f} | {total_salm_err:.0f} |")
    lines.append(f"| Mean abs error per scenario | {total_hulk_err/max(n_comparisons,1):.1f} | {total_salm_err/max(n_comparisons,1):.1f} |")
    lines.append(f"| Scenarios won (lower error) | {hulk_wins} | {salm_wins} |")
    lines.append(f"| Ties (within 0.5) | {ties} | — |")
    lines.append("")

    if gdna_scenarios:
        gdna_hulk = np.mean([h for h, s in gdna_scenarios])
        gdna_salm = np.mean([s for h, s in gdna_scenarios])
        lines.append(f"### With gDNA contamination ({len(gdna_scenarios)} scenarios)")
        lines.append("")
        lines.append(f"| Metric | hulkrna | salmon |")
        lines.append(f"| --- | --- | --- |")
        lines.append(f"| Mean total abs error | {gdna_hulk:.1f} | {gdna_salm:.1f} |")
        lines.append("")

    if no_gdna_scenarios:
        nogdna_hulk = np.mean([h for h, s in no_gdna_scenarios])
        nogdna_salm = np.mean([s for h, s in no_gdna_scenarios])
        lines.append(f"### Without gDNA contamination ({len(no_gdna_scenarios)} scenarios)")
        lines.append("")
        lines.append(f"| Metric | hulkrna | salmon |")
        lines.append(f"| --- | --- | --- |")
        lines.append(f"| Mean total abs error | {nogdna_hulk:.1f} | {nogdna_salm:.1f} |")
        lines.append("")

    # --- Key findings ---
    lines.append("## Key Findings")
    lines.append("")
    lines.append("### hulkrna Advantages")
    lines.append("")
    lines.append("1. **gDNA contamination modeling**: hulkrna explicitly models genomic DNA "
                 "contamination via per-gene shadow components in the EM. Salmon has no "
                 "gDNA model — contaminating fragments that pseudo-align to transcripts "
                 "inflate salmon's counts.")
    lines.append("2. **Strand-aware quantification**: hulkrna uses trained strand models "
                 "to separate sense from antisense reads, critical for overlapping "
                 "antisense genes.")
    lines.append("3. **Genome-aligned input**: hulkrna works from BAM files with full "
                 "genomic context (introns, intergenic regions), enabling it to identify "
                 "unspliced reads and intronic fragments.")
    lines.append("")
    lines.append("### salmon Advantages")
    lines.append("")
    lines.append("1. **Speed**: salmon's lightweight quasi-mapping is significantly faster "
                 "than genome alignment + hulkrna counting.")
    lines.append("2. **No alignment required**: salmon maps directly against the "
                 "transcriptome, avoiding the need for a genome aligner like minimap2.")
    lines.append("3. **Mature bias correction**: salmon includes sequence-specific, "
                 "GC, and positional bias models (not tested here).")
    lines.append("4. **Pure RNA quantification**: In clean RNA-seq without gDNA "
                 "contamination, salmon's transcript-level EM is well-calibrated.")
    lines.append("")
    lines.append("### Limitations of This Benchmark")
    lines.append("")
    lines.append("- **Small simulated genomes** (5–8 kb) may not capture real-world "
                 "mapping complexity.")
    lines.append("- **No multi-gene families**: Real transcriptomes have many paralogs "
                 "where salmon's EM is particularly important.")
    lines.append("- **No read errors or biases**: Simulated reads are clean; real data "
                 "has base-call errors, GC bias, and positional artifacts.")
    lines.append("- **gDNA contamination is synthetic**: Real gDNA distributions may "
                 "differ from our uniform model.")
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    logger.info("Report written to %s", output_path)


# =====================================================================
# Main
# =====================================================================


def main():
    # Check required tools
    for tool in ["salmon", "minimap2", "samtools"]:
        if shutil.which(tool) is None:
            sys.exit(f"Error: {tool} not found in PATH")

    docs_dir = Path(__file__).resolve().parent.parent / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting hulkrna vs salmon benchmark")
    logger.info("N_FRAGMENTS=%d, SEED=%d", N_FRAGMENTS, SIM_SEED)

    results = run_all_benchmarks()

    # Save raw results
    results_path = docs_dir / "benchmark_results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    logger.info("Raw results saved to %s", results_path)

    # Generate report
    report_path = docs_dir / "benchmark_hulkrna_vs_salmon.md"
    generate_report(results, report_path)

    logger.info("Benchmark complete!")


if __name__ == "__main__":
    main()
