"""Explore g20_n0_s65 failure across fragment counts, gating, density floor.

Three axes:
  1. Fragment count: 2000 vs 10000
  2. gDNA gating: on vs off
  3. _MIN_DENSITY_REGIONS: 20 (current) vs 1 (allow any)
"""

import logging
import numpy as np

from rigel.config import EMConfig, PipelineConfig, BamScanConfig, CalibrationConfig
from rigel.sim import Scenario, GDNAConfig, SimConfig, run_benchmark
from rigel.pipeline import run_pipeline

logging.basicConfig(level=logging.WARNING)

SEED = 42


def make_scenario(work_dir):
    sc = Scenario(
        "diag_g20_n0_s65",
        genome_length=20000,
        seed=SEED,
        work_dir=work_dir,
    )
    sc.add_gene("g1", "+", [{
        "t_id": "t1",
        "exons": [(2000, 4000), (8000, 10000)],
        "abundance": 100,
    }])
    sc.add_gene("g_ctrl", "-", [{
        "t_id": "t_ctrl",
        "exons": [(14000, 16000), (18000, 19000)],
        "abundance": 0,
    }])
    return sc


def run_one(scenario, n_frags, ss=0.65, gdna_abundance=20, nrna_abundance=0, label=""):
    sim_cfg = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=ss, seed=SEED,
    )
    gdna_cfg = GDNAConfig(
        abundance=gdna_abundance, frag_mean=350, frag_std=100,
        frag_min=100, frag_max=1000,
    ) if gdna_abundance > 0 else None

    result = scenario.build_oracle(
        n_fragments=n_frags, sim_config=sim_cfg,
        gdna_config=gdna_cfg, nrna_abundance=nrna_abundance,
    )

    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=config)
    bench = run_benchmark(result, pr, scenario_name=label)

    total_rna_expected = bench.total_expected + bench.n_nrna_expected
    total_rna_observed = bench.total_rna_observed
    err = abs(total_rna_observed - total_rna_expected) / max(total_rna_expected, 1)

    cal = pr.calibration
    gamma = 0.0
    if cal is not None and cal.region_n_total.sum() > 0:
        gamma = cal.region_e_gdna.sum() / cal.region_n_total.sum()

    return {
        "label": label,
        "n_frags": n_frags,
        "gdna_expected": bench.n_gdna_expected,
        "gdna_pipeline": bench.n_gdna_pipeline,
        "nrna_expected": bench.n_nrna_expected,
        "nrna_pipeline": bench.n_nrna_pipeline,
        "mrna_expected": bench.total_expected,
        "mrna_observed": bench.total_observed,
        "total_rna_expected": total_rna_expected,
        "total_rna_observed": total_rna_observed,
        "rna_rel_err": err,
        "gamma": gamma,
        "ss_trained": cal.strand_specificity if cal else 0,
    }


# ── Experiment 1: Fragment count scaling ──
print("=" * 70)
print("EXPERIMENT 1: Fragment Count Scaling (g20_n0_s65)")
print("=" * 70)

for n in [500, 1000, 2000, 5000, 10000]:
    sc = make_scenario(f"/tmp/diag_frag_{n}")
    r = run_one(sc, n, label=f"n={n}")
    print(f"  n={n:5d}: SS_trained={r['ss_trained']:.3f} γ={r['gamma']:.4f} "
          f"gdna_exp={r['gdna_expected']:4d} gdna_pipe={r['gdna_pipeline']:6.0f} "
          f"nrna_pipe={r['nrna_pipeline']:5.0f} "
          f"RNA_err={r['rna_rel_err']:.4f}")
    sc.cleanup()

print()
