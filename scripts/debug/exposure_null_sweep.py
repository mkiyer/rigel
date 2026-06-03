#!/usr/bin/env python
"""Empirical-null sweep for the exposure dispersion (theory: noise-floor doc §4/§6).

Runs the pipeline on GENUINELY-UNIFORM libraries (ω≡1: uniform gDNA, no capture), where
the true exposure dispersion φ=0, and records the dispersion φ̂_null the estimator
*hallucinates*. By construction this is the total noise floor (Poisson + deconvolution +
collapse). We characterize Φ₀ along the two sparsity axes the theory unifies:

  - REGION sweep:   scale K gene-blocks with coverage-per-block held constant (R grows).
  - COVERAGE sweep: fix K, scale the fragment budget (gDNA/region grows).

For each scenario we also record Σμ² = ρ₀²·Σ L_r² and the Poisson prediction
SE_pois = √(2/Σμ²): if φ̂_null ≫ SE_pois the deconvolution term (theory §3) dominates and a
closed-form Poisson floor would under-shrink.

Run with the regularizer OFF (κ=0) to see the RAW floor the prior must subtract.

Usage:  python scripts/debug/exposure_null_sweep.py [--axis region|coverage|both]
"""

from __future__ import annotations

import argparse
import tempfile

import numpy as np

from rigel.calibration.region_arrays import RegionArrays
from rigel.config import BamScanConfig, CalibrationConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark

_T1 = [(2000, 2500), (3000, 3500), (4000, 4500), (5000, 5500),
       (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)]
_TC = [(12000, 12500), (13000, 13500), (14000, 14500), (15000, 15500),
       (16000, 16500), (17000, 17500), (18000, 18500), (19000, 19500)]


def _run(K: int, n_frag: int):
    """Build a uniform-gDNA K-block scenario, run at κ=0, return null metrics."""
    tmp = tempfile.mkdtemp()
    sc = Scenario("null", genome_length=K * 20000, seed=42, work_dir=tmp)
    for i in range(K):
        o = i * 20000
        sc.add_gene(f"g1_{i}", "+", [{"t_id": f"t1_{i}", "exons": [(o + a, o + b) for a, b in _T1], "abundance": 100}])
        sc.add_gene(f"gc_{i}", "-", [{"t_id": f"tc_{i}", "exons": [(o + a, o + b) for a, b in _TC], "abundance": 0}])
    cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                        read_length=100, strand_specificity=0.65, seed=42)
    res = sc.build_oracle(n_fragments=n_frag, sim_config=cfg,
                          gdna_config=GDNAConfig(abundance=20, frag_mean=350, frag_std=100,
                                                 frag_min=100, frag_max=1000),
                          nrna_abundance=70)
    pr = run_pipeline(res.bam_path, res.index,
                      config=PipelineConfig(em=EMConfig(seed=42),
                                            scan=BamScanConfig(sj_strand_tag="auto"),
                                            calibration=CalibrationConfig(exposure_prior_pseudocount=0.0)))
    cal = pr.calibration
    b = run_benchmark(res, pr, scenario_name="null")
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    lengths = (ra.end - ra.start).astype(np.float64)
    m_g = cal.mass_g_contained + cal.mass_g_left + cal.mass_g_right
    sum_mu2 = cal.rho_0**2 * float(np.dot(lengths, lengths))
    se_pois = float(np.sqrt(2.0 / sum_mu2)) if sum_mu2 > 0 else float("inf")
    rna_exp = b.total_expected + b.n_nrna_expected
    return dict(
        R=cal.n_regions, phi_null=cal.exposure_dispersion, rho_0=cal.rho_0,
        gdna_tot=float(m_g.sum()), gdna_per_region=float(m_g.sum()) / cal.n_regions,
        sum_mu2=sum_mu2, se_pois=se_pois,
        rna_rec=(b.total_rna_observed / rna_exp if rna_exp else 0.0),
    )


def _hdr():
    print(f"{'K':>3} {'R':>5} {'nfrag':>7} {'rho_0':>9} {'gDNA/reg':>9} {'phi_null':>9} "
          f"{'Sum_mu2':>10} {'SE_pois':>9} {'phi/SE':>8} {'RNArec':>7}")


def _row(K, nfrag, m):
    ratio = m["phi_null"] / m["se_pois"] if m["se_pois"] > 0 else float("inf")
    print(f"{K:>3} {m['R']:>5} {nfrag:>7} {m['rho_0']:>9.3g} {m['gdna_per_region']:>9.2f} "
          f"{m['phi_null']:>9.3g} {m['sum_mu2']:>10.3g} {m['se_pois']:>9.3g} {ratio:>8.1f} "
          f"{m['rna_rec']:>6.0%}")


def region_sweep():
    print("\n=== EMPIRICAL NULL — REGION sweep (coverage/block constant: nfrag = K*2000) ===")
    _hdr()
    for K in (1, 2, 5, 10, 20, 50):
        _row(K, K * 2000, _run(K, K * 2000))
    print("  Theory: φ_null should fall as R grows; if φ/SE ≫ 1, deconvolution floor dominates Poisson.")


def coverage_sweep():
    print("\n=== EMPIRICAL NULL — COVERAGE sweep (R fixed at K=5; scale fragment budget) ===")
    _hdr()
    for nfrag in (2500, 5000, 10000, 20000, 40000, 80000):
        _row(5, nfrag, _run(5, nfrag))
    print("  Theory: φ_null should fall as coverage grows; compare decay rate to SE_pois ∝ 1/√Σμ².")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--axis", choices=["region", "coverage", "both"], default="both")
    args = ap.parse_args()
    if args.axis in ("region", "both"):
        region_sweep()
    if args.axis in ("coverage", "both"):
        coverage_sweep()


if __name__ == "__main__":
    main()
