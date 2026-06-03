#!/usr/bin/env python
"""Prototype 2 (proper acyclic deconvolution) — demonstrate recovery without collapse.

Adds the three components Prototype 1 proved load-bearing:
  - §3 strand POSTERIOR (regularized) instead of the raw method-of-moments split;
  - §8 total mass per region = contained + left + right (boundary-crossing gDNA);
  - tier-F impute-and-exclude: RNA-dominated regions don't set ρ₀ (they'd dilute it);
    their baseline gDNA is imputed at ρ₀·L.

Strand posterior (per region): a fragment is + (genome) with prob
    p₊(π_g) = ½·π_g + r_rna·(1−π_g),     r_rna = κ_rna (POS gene) | 1−κ_rna (NEG gene)
n₊ ~ Binomial(N, p₊); flat prior on π_g ∈ [0,1]; report posterior mean (and the z=0 default).
Intergenic (TS_NONE) → π_g=1 (gDNA by annotation). AMBIG → imputed at baseline.

ρ₀ from the gDNA-bearing pool (π_g ≥ ½), impute the rest. (The ½ split is a throwaway-proto
stand-in for the production smooth confidence weight; this is a concept check, not shipping.)

Usage:  python scripts/debug/proto_acyclic_deconv2.py
"""

from __future__ import annotations

import tempfile

import numpy as np

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG, TS_NEG, TS_NONE
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark

_T1 = [(2000, 2500), (3000, 3500), (4000, 4500), (5000, 5500),
       (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)]
_TC = [(12000, 12500), (13000, 13500), (14000, 14500), (15000, 15500),
       (16000, 16500), (17000, 17500), (18000, 18500), (19000, 19500)]
_CFG = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                     read_length=100, strand_specificity=0.65, seed=42)
_GD = GDNAConfig(abundance=20, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)


def build(K):
    tmp = tempfile.mkdtemp(prefix="proto2_")
    sc = Scenario("dc", genome_length=K * 20000, seed=42, work_dir=tmp)
    for i in range(K):
        o = i * 20000
        sc.add_gene(f"g1_{i}", "+", [{"t_id": f"t1_{i}", "exons": [(o + a, o + b) for a, b in _T1], "abundance": 100}])
        sc.add_gene(f"gc_{i}", "-", [{"t_id": f"tc_{i}", "exons": [(o + a, o + b) for a, b in _TC], "abundance": 0}])
    return sc.build_oracle(n_fragments=K * 2000, sim_config=_CFG, gdna_config=_GD, nrna_abundance=70)


def strand_posterior_pi_g(pos, neg, ts, kappa):
    """Posterior-mean gDNA fraction per region (flat prior), vectorized over a π_g grid."""
    grid = np.linspace(1e-4, 1 - 1e-4, 200)
    r_rna = np.where(ts == TS_NEG, 1.0 - kappa, kappa)  # RNA '+' rate by gene strand
    out = np.full(pos.shape, 0.5)
    for i in range(pos.size):
        n = pos[i] + neg[i]
        if ts[i] == TS_NONE:
            out[i] = 1.0  # intergenic: gDNA by annotation
            continue
        if ts[i] == TS_AMBIG or n < 0.5:
            out[i] = np.nan  # undecodable → imputed later
            continue
        p_plus = 0.5 * grid + r_rna[i] * (1.0 - grid)
        loglik = pos[i] * np.log(p_plus) + neg[i] * np.log1p(-p_plus)
        w = np.exp(loglik - loglik.max())
        out[i] = float(np.dot(grid, w) / w.sum())
    return out


def run(K):
    res = build(K)
    scan = BamScanConfig(sj_strand_tag="auto")
    _, sm, _, _, payload = scan_and_buffer(str(res.bam_path), res.index, scan)
    sm.finalize()
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(payload, ra)
    kappa = fit_strand_balance(sm).kappa_rna

    # total per-region flux + mass = contained + both boundary sides (§8)
    pos = (sub.contained.n_unspliced_pos + sub.left.n_unspliced_pos + sub.right.n_unspliced_pos).astype(np.float64)
    neg = (sub.contained.n_unspliced_neg + sub.left.n_unspliced_neg + sub.right.n_unspliced_neg).astype(np.float64)
    mass = sub.contained.mass_unspliced + sub.left.mass_unspliced + sub.right.mass_unspliced
    L = sub.L_eff
    ts = sub.ts_class

    pi_g = strand_posterior_pi_g(pos, neg, ts, kappa)
    measurable = ~np.isnan(pi_g)
    g_meas = np.where(measurable, np.nan_to_num(pi_g) * mass, 0.0)

    # ρ₀ from the gDNA-bearing pool (π_g ≥ ½); impute the rest at ρ₀·L (tier F)
    pool = measurable & (pi_g >= 0.5)
    rho0 = float(g_meas[pool].sum()) / float(L[pool].sum()) if L[pool].sum() > 0 else 0.0
    g_final = np.where(pool, g_meas, rho0 * L)  # excluded/RNA-dominated → baseline impute

    pr = run_pipeline(res.bam_path, res.index, config=PipelineConfig(em=EMConfig(seed=42), scan=scan))
    b = run_benchmark(res, pr, scenario_name=f"K{K}")
    return dict(K=K, R=sub.n_regions, kappa=kappa, rho0=rho0, rho0_em=pr.calibration.rho_0,
                g_total=float(g_final.sum()), g_pool=float(g_meas[pool].sum()),
                n_pool=int(pool.sum()), g_exp=b.n_gdna_expected)


def main():
    print(f"{'K':>3} {'R':>5} {'kappa':>6} {'pool':>5} {'rho0_acyc':>10} {'rho0_EM':>10} "
          f"{'G_acyc':>8} {'G_exp':>7} {'recov':>6}")
    for K in (1, 5):
        m = run(K)
        rec = m["g_total"] / m["g_exp"] if m["g_exp"] else 0
        print(f"{m['K']:>3} {m['R']:>5} {m['kappa']:>6.3f} {m['n_pool']:>5} {m['rho0']:>10.4g} "
              f"{m['rho0_em']:>10.4g} {m['g_total']:>8.0f} {m['g_exp']:>7.0f} {rec:>5.0%}")
    print("\n  rho0_acyc (acyclic, strand-posterior + boundary + impute) vs rho0_EM (collapsed).")
    print("  Headline: rho0_acyc bounded & ~uniform-density; G_acyc ≈ G_exp; NO 5e-5 collapse.")


if __name__ == "__main__":
    main()
