#!/usr/bin/env python3
"""Evaluate the impact of filter_by_overlap on accuracy.

Compares three modes:
  1. filter_by_overlap ON  (default, overlap_min_frac=0.99)
  2. filter_by_overlap OFF (overlap_min_frac=0.0)
  3. filter_by_overlap STRICT (overlap_min_frac=0.5)

For each mode, runs a representative subset of scenarios and
measures per-transcript accuracy, gDNA accuracy, and EM problem size.
"""
from __future__ import annotations

import sys
import os
import time
import logging
from pathlib import Path
from collections import defaultdict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

logging.basicConfig(level=logging.WARNING)

SIM_SEED = 42
PIPELINE_SEED = 42

def _sim_config(*, ss: float = 1.0):
    return SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=ss, seed=SIM_SEED,
    )

def _gdna_config(abundance: float):
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance, frag_mean=350, frag_std=100,
        frag_min=100, frag_max=1000,
    )


# --- Scenario builders ---

def build_spliced(tmp):
    sc = Scenario("spliced", genome_length=5000, seed=SIM_SEED, work_dir=tmp / "spliced")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300), (1800, 2100)], "abundance": 100},
    ])
    sc.add_gene("g_ctrl", "+", [
        {"t_id": "t_ctrl", "exons": [(3500, 3800)], "abundance": 0},
    ])
    return sc

def build_antisense(tmp):
    sc = Scenario("antisense", genome_length=8000, seed=SIM_SEED, work_dir=tmp / "antisense")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300)], "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(300, 600), (1100, 1400)], "abundance": 100},
    ])
    sc.add_gene("g_ctrl", "+", [
        {"t_id": "t_ctrl", "exons": [(5500, 5800)], "abundance": 0},
    ])
    return sc

def build_single_exon(tmp):
    sc = Scenario("single_exon", genome_length=5000, seed=SIM_SEED, work_dir=tmp / "single_exon")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 700)], "abundance": 100},
    ])
    sc.add_gene("g_ctrl", "+", [
        {"t_id": "t_ctrl", "exons": [(3500, 3800)], "abundance": 0},
    ])
    return sc

def build_isoform(tmp):
    sc = Scenario("iso_equal", genome_length=8000, seed=SIM_SEED, work_dir=tmp / "iso_equal")
    sc.add_gene("g1", "+", [
        {"t_id": "t1a", "exons": [(200, 500), (1000, 1300), (1800, 2100)], "abundance": 100},
        {"t_id": "t1b", "exons": [(200, 500), (1200, 1500), (1800, 2100)], "abundance": 100},
    ])
    sc.add_gene("g_ctrl", "+", [
        {"t_id": "t_ctrl", "exons": [(5500, 5800)], "abundance": 0},
    ])
    return sc

def build_nonoverlap(tmp):
    sc = Scenario("nonoverlap", genome_length=8000, seed=SIM_SEED, work_dir=tmp / "nonoverlap")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300), (1800, 2100)], "abundance": 100},
    ])
    sc.add_gene("g2", "+", [
        {"t_id": "t2", "exons": [(3000, 3800)], "abundance": 50},
    ])
    sc.add_gene("g_ctrl", "+", [
        {"t_id": "t_ctrl", "exons": [(5500, 5800)], "abundance": 0},
    ])
    return sc


SCENARIOS = [
    ("spliced", build_spliced),
    ("antisense", build_antisense),
    ("single_exon", build_single_exon),
    ("isoform", build_isoform),
    ("nonoverlap", build_nonoverlap),
]

CONDITIONS = [
    {"gdna": 0, "ss": 1.0, "nrna": 0.0},
    {"gdna": 0, "ss": 0.9, "nrna": 0.0},
    {"gdna": 20, "ss": 1.0, "nrna": 0.0},
    {"gdna": 20, "ss": 0.9, "nrna": 0.3},
    {"gdna": 50, "ss": 0.9, "nrna": 0.0},
    {"gdna": 0, "ss": 1.0, "nrna": 0.5},
    {"gdna": 0, "ss": 0.65, "nrna": 0.0},
]

OVERLAP_MODES = {
    "OFF (0.0)":     0.0,   # No filtering at all
    "LOOSE (0.5)":   0.5,   # Keep candidates within 50% of best
    "DEFAULT (0.99)": 0.99,  # Current: keep within 1% of best
}


def run_one(sc_name, sc_builder, cond, overlap_min_frac, tmp_root):
    """Build scenario, run pipeline, return metrics dict."""
    cond_name = f"gdna={cond['gdna']}_ss={cond['ss']}_nrna={int(cond['nrna']*100)}"
    work = tmp_root / f"{sc_name}_{cond_name}"
    work.mkdir(parents=True, exist_ok=True)

    sc = sc_builder(work)
    try:
        gdna_cfg = _gdna_config(cond["gdna"])
        sim_cfg = _sim_config(ss=cond["ss"])
        result = sc.build(
            n_fragments=500, sim_config=sim_cfg,
            gdna_config=gdna_cfg, nrna_fraction=cond["nrna"],
        )
        pr = run_pipeline(
            result.bam_path, result.index,
            sj_strand_tag="ts", seed=PIPELINE_SEED,
            overlap_min_frac=overlap_min_frac,
        )
        bench = run_benchmark(result, pr, scenario_name=f"{sc_name}_{cond_name}")
        stats = pr.stats

        # Gather per-transcript errors
        t_errors = []
        for t in bench.transcripts:
            if t.t_id == "t_ctrl":
                continue
            t_errors.append({
                "t_id": t.t_id,
                "expected": t.expected,
                "observed": t.observed,
                "abs_diff": t.abs_diff,
                "rel_err": t.rel_error,
            })

        ctrl = next((t for t in bench.transcripts if t.t_id == "t_ctrl"), None)
        ctrl_fp = ctrl.observed if ctrl else 0

        return {
            "scenario": sc_name,
            "condition": cond_name,
            "total_expected": bench.total_expected,
            "total_rna_obs": bench.total_rna_observed,
            "gdna_exp": bench.n_gdna_expected,
            "gdna_pipe": bench.n_gdna_pipeline,
            "ctrl_fp": ctrl_fp,
            "mean_abs": np.mean([e["abs_diff"] for e in t_errors]) if t_errors else 0,
            "mean_rel": np.mean([e["rel_err"] for e in t_errors]) if t_errors else 0,
            "em_unique": stats.em_routed_unique_units,
            "em_iso": stats.em_routed_isoform_ambig_units,
            "em_gene": stats.em_routed_gene_ambig_units,
            "det_unique": stats.deterministic_unique_units,
            "t_errors": t_errors,
        }
    finally:
        sc.cleanup()


def main():
    import tempfile
    tmp_root = Path(tempfile.mkdtemp(prefix="overlap_eval_"))

    # Collect results per mode
    all_results: dict[str, list[dict]] = {}
    for mode_name, frac in OVERLAP_MODES.items():
        print(f"\n{'='*60}")
        print(f"  MODE: {mode_name} (overlap_min_frac={frac})")
        print(f"{'='*60}")
        results = []
        for sc_name, sc_builder in SCENARIOS:
            for cond in CONDITIONS:
                cond_name = f"gdna={cond['gdna']}_ss={cond['ss']}_nrna={int(cond['nrna']*100)}"
                r = run_one(sc_name, sc_builder, cond, frac, tmp_root)
                results.append(r)
                print(f"  {sc_name:15s} {cond_name:30s} "
                      f"abs={r['mean_abs']:6.1f}  rel={r['mean_rel']:.3f}  "
                      f"ctrl_fp={r['ctrl_fp']:.0f}  "
                      f"em_units={r['em_unique']+r['em_iso']+r['em_gene']}")
        all_results[mode_name] = results

    # --- Comparison ---
    print(f"\n\n{'='*80}")
    print("COMPARISON: filter_by_overlap modes")
    print(f"{'='*80}")

    # Per-scenario aggregate
    for sc_name, _ in SCENARIOS:
        print(f"\n--- {sc_name} ---")
        header = f"{'Mode':<20s} {'MeanAbs':>8s} {'MeanRel':>8s} {'CtrlFP':>7s} {'EMunits':>8s} {'DetUniq':>8s}"
        print(header)
        for mode_name in OVERLAP_MODES:
            sc_results = [r for r in all_results[mode_name] if r["scenario"] == sc_name]
            mean_abs = np.mean([r["mean_abs"] for r in sc_results])
            mean_rel = np.mean([r["mean_rel"] for r in sc_results])
            total_ctrl = sum(r["ctrl_fp"] for r in sc_results)
            total_em = sum(r["em_unique"] + r["em_iso"] + r["em_gene"] for r in sc_results)
            total_det = sum(r["det_unique"] for r in sc_results)
            print(f"{mode_name:<20s} {mean_abs:8.1f} {mean_rel:8.3f} {total_ctrl:7.0f} {total_em:8d} {total_det:8d}")

    # Overall aggregate
    print(f"\n--- OVERALL ---")
    header = f"{'Mode':<20s} {'MeanAbs':>8s} {'MeanRel':>8s} {'CtrlFP':>7s} {'EMunits':>8s} {'DetUniq':>8s}"
    print(header)
    for mode_name in OVERLAP_MODES:
        results = all_results[mode_name]
        mean_abs = np.mean([r["mean_abs"] for r in results])
        mean_rel = np.mean([r["mean_rel"] for r in results])
        total_ctrl = sum(r["ctrl_fp"] for r in results)
        total_em = sum(r["em_unique"] + r["em_iso"] + r["em_gene"] for r in results)
        total_det = sum(r["det_unique"] for r in results)
        print(f"{mode_name:<20s} {mean_abs:8.1f} {mean_rel:8.3f} {total_ctrl:7.0f} {total_em:8d} {total_det:8d}")

    # Detailed condition-level diff: OFF vs DEFAULT
    print(f"\n\n--- Per-condition diff: OFF vs DEFAULT ---")
    off_results = all_results["OFF (0.0)"]
    def_results = all_results["DEFAULT (0.99)"]
    print(f"{'Scenario':15s} {'Condition':30s} {'OFF_abs':>8s} {'DEF_abs':>8s} {'OFF_em':>7s} {'DEF_em':>7s} {'delta_abs':>10s}")
    for r_off, r_def in zip(off_results, def_results):
        em_off = r_off["em_unique"] + r_off["em_iso"] + r_off["em_gene"]
        em_def = r_def["em_unique"] + r_def["em_iso"] + r_def["em_gene"]
        delta = r_off["mean_abs"] - r_def["mean_abs"]
        marker = "  <-- OFF better" if delta < -1 else ("  <-- DEF better" if delta > 1 else "")
        print(f"{r_off['scenario']:15s} {r_off['condition']:30s} "
              f"{r_off['mean_abs']:8.1f} {r_def['mean_abs']:8.1f} "
              f"{em_off:7d} {em_def:7d} {delta:+10.1f}{marker}")

    # Cleanup
    import shutil
    shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    main()
