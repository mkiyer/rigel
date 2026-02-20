#!/usr/bin/env python3
"""Run all test_scenarios.py scenario classes and produce an accuracy audit.

Outputs a Markdown report with per-transcript accuracy, gDNA accuracy,
nRNA detection, and negative-control false positives for every scenario
and condition combination.

Usage:
    python scripts/scenario_accuracy_audit.py [--outdir /path/to/dir]
"""
from __future__ import annotations

import argparse
import csv
import io
import logging
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Ensure src is on path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

logger = logging.getLogger(__name__)

SIM_SEED = 42
PIPELINE_SEED = 42
N_FRAGMENTS = 500
N_FRAGMENTS_ISO = 1000  # isoform / paralog scenarios use more fragments

GDNA_LEVELS = [0, 5, 20, 50, 100]
STRAND_LEVELS = [0.65, 0.8, 0.9, 0.95, 1.0]
NRNA_FRACTIONS = [0.0, 0.1, 0.3, 0.5, 0.7]


def _sim_config(*, strand_specificity: float = 1.0):
    return SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=strand_specificity,
        seed=SIM_SEED,
    )


def _gdna_config(abundance: float):
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance, frag_mean=350, frag_std=100,
        frag_min=100, frag_max=1000,
    )


@dataclass
class Row:
    scenario: str
    condition: str
    t_id: str
    expected: float
    observed: float
    abs_diff: float
    rel_error: float
    # Aggregate fields
    total_expected: float = 0
    total_observed: float = 0
    gdna_expected: float = 0
    gdna_pipeline: float = 0
    gdna_abs_diff: float = 0
    nrna_expected: float = 0
    nrna_pipeline: float = 0
    nrna_abs_diff: float = 0
    ctrl_observed: float = 0
    n_fragments: int = 0
    alignment_rate: float = 0.0
    n_intergenic: int = 0


def _run_one(scenario, *, n_fragments=N_FRAGMENTS, gdna_abundance=0,
             strand_specificity=1.0, nrna_fraction=0.0,
             include_multimap=False, scenario_name=""):
    sim_config = _sim_config(strand_specificity=strand_specificity)
    gdna = _gdna_config(gdna_abundance)
    result = scenario.build(
        n_fragments=n_fragments, sim_config=sim_config,
        gdna_config=gdna, nrna_fraction=nrna_fraction,
    )
    pr = run_pipeline(
        result.bam_path, result.index,
        sj_strand_tag="ts", seed=PIPELINE_SEED,
        include_multimap=include_multimap,
    )
    bench = run_benchmark(result, pr, scenario_name=scenario_name)
    return bench


def _collect_rows(bench, scenario_name: str, condition: str) -> list[Row]:
    rows = []
    ctrl_obs = 0.0
    for ta in bench.transcripts:
        if ta.t_id == "t_ctrl":
            ctrl_obs = ta.observed
    for ta in bench.transcripts:
        rows.append(Row(
            scenario=scenario_name,
            condition=condition,
            t_id=ta.t_id,
            expected=ta.expected,
            observed=ta.observed,
            abs_diff=ta.abs_diff,
            rel_error=ta.rel_error,
            total_expected=bench.total_expected,
            total_observed=bench.total_observed,
            gdna_expected=bench.n_gdna_expected,
            gdna_pipeline=bench.n_gdna_pipeline,
            gdna_abs_diff=bench.gdna_abs_diff,
            nrna_expected=bench.n_nrna_expected,
            nrna_pipeline=bench.n_nrna_pipeline,
            nrna_abs_diff=bench.nrna_abs_diff,
            ctrl_observed=ctrl_obs,
            n_fragments=bench.n_fragments,
            alignment_rate=bench.alignment_rate,
            n_intergenic=bench.n_intergenic,
        ))
    return rows


# =====================================================================
# Scenario builders
# =====================================================================

def make_single_exon(work_dir: Path):
    sc = Scenario("single_exon", genome_length=8000, seed=SIM_SEED,
                   work_dir=work_dir / "single_exon")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
    ])
    sc.add_gene("g_helper", "+", [
        {"t_id": "t_helper",
         "exons": [(2500, 3000), (3500, 4000)], "abundance": 50},
    ])
    sc.add_gene("g_ctrl", "-", [
        {"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0},
    ])
    return sc


def make_spliced(work_dir: Path):
    sc = Scenario("spliced_gene", genome_length=5000, seed=SIM_SEED,
                   work_dir=work_dir / "spliced_gene")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
         "abundance": 100},
    ])
    sc.add_gene("g_ctrl", "-", [
        {"t_id": "t_ctrl", "exons": [(3500, 3800)], "abundance": 0},
    ])
    return sc


def make_non_overlapping(work_dir: Path):
    sc = Scenario("non_overlapping", genome_length=10000, seed=SIM_SEED,
                   work_dir=work_dir / "non_overlapping")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
         "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(4000, 4400)], "abundance": 100},
    ])
    sc.add_gene("g_ctrl", "+", [
        {"t_id": "t_ctrl", "exons": [(8000, 8300)], "abundance": 0},
    ])
    return sc


def make_two_isoforms(work_dir: Path, major=100, minor=100, suffix=""):
    sc = Scenario("two_isoforms" + suffix, genome_length=6000,
                   seed=SIM_SEED,
                   work_dir=work_dir / ("two_isoforms" + suffix))
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(200, 500), (1000, 1300), (2000, 2300)],
         "abundance": major},
        {"t_id": "t2",
         "exons": [(200, 500), (2000, 2300)],
         "abundance": minor},
    ])
    sc.add_gene("g_ctrl", "-", [
        {"t_id": "t_ctrl", "exons": [(4500, 4800)], "abundance": 0},
    ])
    return sc


def make_antisense(work_dir: Path, g1=100, g2=100, suffix=""):
    sc = Scenario("overlap_anti" + suffix, genome_length=8000,
                   seed=SIM_SEED,
                   work_dir=work_dir / ("overlap_anti" + suffix))
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
         "abundance": g1},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(300, 600), (1100, 1400)],
         "abundance": g2},
    ])
    sc.add_gene("g_ctrl", "+", [
        {"t_id": "t_ctrl", "exons": [(5500, 5800)], "abundance": 0},
    ])
    return sc


def make_paralogs(work_dir: Path, g1=100, g2=100, *, spliced=False, suffix=""):
    sc = Scenario("paralogs" + suffix, genome_length=12000,
                   seed=SIM_SEED,
                   work_dir=work_dir / ("paralogs" + suffix))
    if spliced:
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": g1},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2", "exons": [(5000, 5300), (5700, 6000)],
             "abundance": g2},
        ])
        sc.genome.edit(5000, sc.genome[500:1500])
    else:
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1000)], "abundance": g1},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2", "exons": [(5000, 5500)], "abundance": g2},
        ])
        sc.genome.edit(5000, sc.genome[500:1000])
    sc.add_gene("g_ctrl", "-", [
        {"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0},
    ])
    return sc


def make_dist_paralogs(work_dir: Path, g1=100, g2=100, suffix=""):
    sc = Scenario("dist_paralogs" + suffix, genome_length=12000,
                   seed=SIM_SEED,
                   work_dir=work_dir / ("dist_paralogs" + suffix))
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
         "abundance": g1},
    ])
    sc.add_gene("g2", "+", [
        {"t_id": "t2", "exons": [(5000, 5300), (5700, 5900)],
         "abundance": g2},
    ])
    sc.genome.edit(5000, sc.genome[500:800])
    sc.add_gene("g_ctrl", "-", [
        {"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0},
    ])
    return sc


# =====================================================================
# Main sweep
# =====================================================================

def run_sweep(work_dir: Path) -> list[Row]:
    all_rows: list[Row] = []

    def _sweep(make_fn, name, *, n_frag=N_FRAGMENTS, mm=False,
               gdna_list=None, ss_list=None, nrna_list=None,
               make_kwargs=None):
        gdna_list = gdna_list or [0]
        ss_list = ss_list or [1.0]
        nrna_list = nrna_list or [0.0]
        make_kwargs = make_kwargs or {}

        for gdna in gdna_list:
            for ss in ss_list:
                for nrna in nrna_list:
                    cond = f"gdna={gdna}_ss={ss}_nrna={int(nrna*100)}"
                    sc = make_fn(work_dir, **make_kwargs)
                    try:
                        bench = _run_one(
                            sc, n_fragments=n_frag,
                            gdna_abundance=gdna,
                            strand_specificity=ss,
                            nrna_fraction=nrna,
                            include_multimap=mm,
                            scenario_name=f"{name}_{cond}",
                        )
                        rows = _collect_rows(bench, name, cond)
                        all_rows.extend(rows)
                        # Log progress
                        t1_row = next(
                            (r for r in rows if r.t_id not in ("t_ctrl", "t_helper")),
                            None,
                        )
                        if t1_row:
                            logger.info(
                                "  %s / %s: exp=%.0f obs=%.0f diff=%.0f "
                                "gdna_exp=%.0f gdna_pipe=%.0f ctrl=%.0f",
                                name, cond, t1_row.expected, t1_row.observed,
                                t1_row.abs_diff, t1_row.gdna_expected,
                                t1_row.gdna_pipeline, t1_row.ctrl_observed,
                            )
                    except Exception as e:
                        logger.error("FAILED %s / %s: %s", name, cond, e)
                    finally:
                        sc.cleanup()

    # -- 1. Single exon --
    logger.info("=== Single Exon ===")
    _sweep(make_single_exon, "single_exon",
           gdna_list=GDNA_LEVELS, ss_list=STRAND_LEVELS, nrna_list=[0.0, 0.3, 0.5])

    # -- 2. Spliced gene --
    logger.info("=== Spliced Gene ===")
    _sweep(make_spliced, "spliced",
           gdna_list=GDNA_LEVELS, ss_list=STRAND_LEVELS, nrna_list=[0.0, 0.3, 0.5])

    # -- 3. Non-overlapping genes --
    logger.info("=== Non-Overlapping ===")
    _sweep(make_non_overlapping, "nonoverlap",
           gdna_list=GDNA_LEVELS, ss_list=STRAND_LEVELS, nrna_list=[0.0, 0.3])

    # -- 4. Two isoforms (equal, 10:1) --
    logger.info("=== Two Isoforms (equal) ===")
    _sweep(make_two_isoforms, "iso_equal", n_frag=N_FRAGMENTS_ISO,
           gdna_list=[0, 20, 50], ss_list=[0.9, 1.0], nrna_list=[0.0, 0.3],
           make_kwargs={"major": 100, "minor": 100})

    logger.info("=== Two Isoforms (10:1) ===")
    _sweep(make_two_isoforms, "iso_10_1", n_frag=N_FRAGMENTS_ISO,
           gdna_list=[0, 20, 50], ss_list=[0.9, 1.0], nrna_list=[0.0, 0.3],
           make_kwargs={"major": 100, "minor": 10, "suffix": "_10_1"})

    # -- 5. Overlapping antisense --
    logger.info("=== Overlapping Antisense ===")
    _sweep(make_antisense, "antisense",
           gdna_list=[0, 20, 50], ss_list=[0.8, 0.9, 0.95, 1.0],
           nrna_list=[0.0, 0.3])

    # -- 6. Identical paralogs --
    logger.info("=== Identical Paralogs (unspliced) ===")
    _sweep(make_paralogs, "paralogs_unspliced", mm=True,
           gdna_list=[0, 20, 50], ss_list=[0.9, 1.0], nrna_list=[0.0, 0.3])

    logger.info("=== Identical Paralogs (spliced) ===")
    _sweep(make_paralogs, "paralogs_spliced", mm=True,
           gdna_list=[0, 20, 50], ss_list=[0.9, 1.0], nrna_list=[0.0],
           make_kwargs={"spliced": True})

    # -- 7. Distinguishable paralogs --
    logger.info("=== Distinguishable Paralogs ===")
    _sweep(make_dist_paralogs, "dist_paralogs", mm=True,
           gdna_list=[0, 20, 50], ss_list=[0.9, 1.0], nrna_list=[0.0, 0.3])

    return all_rows


# =====================================================================
# Report writing
# =====================================================================

def write_csv(rows: list[Row], path: Path):
    fields = [
        "scenario", "condition", "t_id", "expected", "observed",
        "abs_diff", "rel_error", "total_expected", "total_observed",
        "gdna_expected", "gdna_pipeline", "gdna_abs_diff",
        "nrna_expected", "nrna_pipeline", "nrna_abs_diff",
        "ctrl_observed", "n_fragments", "alignment_rate", "n_intergenic",
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            d = {k: getattr(r, k) for k in fields}
            # Round floats
            for k in ["observed", "abs_diff", "rel_error", "total_observed",
                       "gdna_pipeline", "gdna_abs_diff", "nrna_pipeline",
                       "nrna_abs_diff", "ctrl_observed", "alignment_rate"]:
                d[k] = round(d[k], 2)
            w.writerow(d)


def write_report(rows: list[Row], path: Path):
    """Write a Markdown accuracy audit report."""
    lines: list[str] = []
    lines.append("# Scenario Accuracy Audit Report\n")

    # Exclude controls/helpers for transcript accuracy
    tx_rows = [r for r in rows if r.t_id not in ("t_ctrl", "t_helper")]
    ctrl_rows = [r for r in rows if r.t_id == "t_ctrl"]

    # ── Overall summary ──
    lines.append("## Overall Summary\n")
    if tx_rows:
        abs_diffs = np.array([r.abs_diff for r in tx_rows])
        rel_errs = np.array([r.rel_error for r in tx_rows
                             if r.expected > 0])
        lines.append(f"- **Total transcript observations**: {len(tx_rows)}")
        lines.append(f"- **Mean abs error**: {abs_diffs.mean():.1f}")
        lines.append(f"- **Median abs error**: {np.median(abs_diffs):.1f}")
        lines.append(f"- **Max abs error**: {abs_diffs.max():.1f}")
        lines.append(f"- **90th pctile abs error**: {np.percentile(abs_diffs, 90):.1f}")
        lines.append(f"- **Mean rel error** (where expected>0): "
                      f"{rel_errs.mean():.3f}")
        lines.append(f"- **Median rel error**: {np.median(rel_errs):.3f}")
        lines.append(f"- **Max rel error**: {rel_errs.max():.3f}")
        lines.append("")

    # ── Negative control false positives ──
    lines.append("## Negative Control False Positives\n")
    if ctrl_rows:
        fp = np.array([r.observed for r in ctrl_rows])
        lines.append(f"- **N conditions**: {len(ctrl_rows)}")
        lines.append(f"- **Median FP**: {np.median(fp):.1f}")
        lines.append(f"- **Mean FP**: {fp.mean():.1f}")
        lines.append(f"- **Max FP**: {fp.max():.1f}")
        lines.append(f"- **Conditions with FP > 5**: "
                      f"{sum(1 for x in fp if x > 5)}")
        lines.append("")
        # Worst false positive cases
        worst_ctrl = sorted(ctrl_rows, key=lambda r: -r.observed)[:20]
        lines.append("### Worst FP cases\n")
        lines.append("| Scenario | Condition | FP count | gDNA_exp | nRNA_exp |")
        lines.append("| --- | --- | ---: | ---: | ---: |")
        for r in worst_ctrl:
            lines.append(
                f"| {r.scenario} | {r.condition} | {r.observed:.0f} "
                f"| {r.gdna_expected:.0f} | {r.nrna_expected:.0f} |"
            )
        lines.append("")

    # ── Per-scenario breakdown ──
    lines.append("## Per-Scenario Accuracy\n")
    scenarios = sorted(set(r.scenario for r in tx_rows))
    for sc_name in scenarios:
        sc_rows = [r for r in tx_rows if r.scenario == sc_name]
        abs_d = np.array([r.abs_diff for r in sc_rows])
        re = np.array([r.rel_error for r in sc_rows if r.expected > 0])
        lines.append(f"### {sc_name}\n")
        lines.append(f"- N observations: {len(sc_rows)}")
        lines.append(f"- Mean abs error: {abs_d.mean():.1f}, "
                      f"Median: {np.median(abs_d):.1f}, "
                      f"Max: {abs_d.max():.1f}")
        if len(re):
            lines.append(f"- Mean rel error: {re.mean():.3f}, "
                          f"Median: {np.median(re):.3f}, "
                          f"Max: {re.max():.3f}")
        lines.append("")

    # ── gDNA accuracy ──
    lines.append("## gDNA Accuracy (where gDNA expected > 0)\n")
    gdna_rows = [r for r in tx_rows if r.gdna_expected > 0]
    # Deduplicate by (scenario, condition) since all transcripts share same gDNA stats
    seen = set()
    gdna_conds = []
    for r in gdna_rows:
        key = (r.scenario, r.condition)
        if key not in seen:
            seen.add(key)
            gdna_conds.append(r)

    if gdna_conds:
        lines.append("| Scenario | Condition | gDNA exp | gDNA pipe | gDNA diff | gDNA rel err |")
        lines.append("| --- | --- | ---: | ---: | ---: | ---: |")
        for r in sorted(gdna_conds, key=lambda x: -x.gdna_abs_diff)[:30]:
            gdna_rel = r.gdna_abs_diff / max(r.gdna_expected, 1)
            lines.append(
                f"| {r.scenario} | {r.condition} "
                f"| {r.gdna_expected:.0f} | {r.gdna_pipeline:.0f} "
                f"| {r.gdna_abs_diff:.0f} | {gdna_rel:.2f} |"
            )
        lines.append("")

    # ── nRNA accuracy ──
    lines.append("## nRNA Accuracy (where nRNA expected > 0)\n")
    nrna_rows = [r for r in tx_rows if r.nrna_expected > 0]
    seen = set()
    nrna_conds = []
    for r in nrna_rows:
        key = (r.scenario, r.condition)
        if key not in seen:
            seen.add(key)
            nrna_conds.append(r)

    if nrna_conds:
        lines.append("| Scenario | Condition | nRNA exp | nRNA pipe | nRNA diff | nRNA rel err |")
        lines.append("| --- | --- | ---: | ---: | ---: | ---: |")
        for r in sorted(nrna_conds, key=lambda x: -x.nrna_abs_diff)[:30]:
            nrna_rel = r.nrna_abs_diff / max(r.nrna_expected, 1)
            lines.append(
                f"| {r.scenario} | {r.condition} "
                f"| {r.nrna_expected:.0f} | {r.nrna_pipeline:.0f} "
                f"| {r.nrna_abs_diff:.0f} | {nrna_rel:.2f} |"
            )
        lines.append("")

    # ── Worst per-transcript errors ──
    lines.append("## Worst Per-Transcript Errors (top 30)\n")
    lines.append("| Scenario | Condition | t_id | Expected | Observed | AbsDiff | RelErr |")
    lines.append("| --- | --- | --- | ---: | ---: | ---: | ---: |")
    for r in sorted(tx_rows, key=lambda x: -x.abs_diff)[:30]:
        lines.append(
            f"| {r.scenario} | {r.condition} | {r.t_id} "
            f"| {r.expected:.0f} | {r.observed:.0f} "
            f"| {r.abs_diff:.0f} | {r.rel_error:.2f} |"
        )
    lines.append("")

    # ── gDNA sweep: RNA count inflation ──
    lines.append("## RNA Count Inflation by gDNA Level\n")
    lines.append("When gDNA is present, do RNA counts inflate (observed > expected)?\n")
    for sc_name in scenarios:
        sc_tx = [r for r in tx_rows if r.scenario == sc_name and
                 r.expected > 10]
        if not sc_tx:
            continue
        # Group by gDNA level
        gdna_groups = {}
        for r in sc_tx:
            # Extract gdna from condition string
            parts = r.condition.split("_")
            gdna_str = parts[0]  # e.g. "gdna=50"
            gdna_groups.setdefault(gdna_str, []).append(r)

        if len(gdna_groups) <= 1:
            continue
        lines.append(f"### {sc_name}\n")
        lines.append("| gDNA level | Mean(obs-exp) | Mean abs_diff | Mean rel_err |")
        lines.append("| --- | ---: | ---: | ---: |")
        for g in sorted(gdna_groups):
            gr = gdna_groups[g]
            bias = np.mean([r.observed - r.expected for r in gr])
            mad = np.mean([r.abs_diff for r in gr])
            mre = np.mean([r.rel_error for r in gr if r.expected > 0])
            lines.append(f"| {g} | {bias:+.1f} | {mad:.1f} | {mre:.3f} |")
        lines.append("")

    # ── Strand specificity impact ──
    lines.append("## Strand Specificity Impact on Accuracy\n")
    for sc_name in scenarios:
        sc_tx = [r for r in tx_rows if r.scenario == sc_name and
                 r.expected > 10]
        if not sc_tx:
            continue
        ss_groups = {}
        for r in sc_tx:
            parts = r.condition.split("_")
            ss_str = parts[1]  # e.g. "ss=0.9"
            ss_groups.setdefault(ss_str, []).append(r)

        if len(ss_groups) <= 1:
            continue
        lines.append(f"### {sc_name}\n")
        lines.append("| SS | Mean(obs-exp) | Mean abs_diff | Mean rel_err | Ctrl FP |")
        lines.append("| --- | ---: | ---: | ---: | ---: |")
        for s in sorted(ss_groups):
            gr = ss_groups[s]
            bias = np.mean([r.observed - r.expected for r in gr])
            mad = np.mean([r.abs_diff for r in gr])
            mre = np.mean([r.rel_error for r in gr if r.expected > 0])
            ctrl_fp = np.mean([r.ctrl_observed for r in gr])
            lines.append(
                f"| {s} | {bias:+.1f} | {mad:.1f} | {mre:.3f} | {ctrl_fp:.1f} |"
            )
        lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("audit_output"))
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        logger.error("minimap2 and/or samtools not in PATH")
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="hulkrna_audit_") as tmpdir:
        work_dir = Path(tmpdir)
        logger.info("Running scenario sweep (work_dir=%s)...", work_dir)
        rows = run_sweep(work_dir)

    csv_path = args.outdir / "scenario_accuracy.csv"
    report_path = args.outdir / "scenario_accuracy_report.md"

    write_csv(rows, csv_path)
    logger.info("Wrote CSV: %s (%d rows)", csv_path, len(rows))

    write_report(rows, report_path)
    logger.info("Wrote report: %s", report_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
