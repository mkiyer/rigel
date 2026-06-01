#!/usr/bin/env python
"""Why don't exposure weights shrink to 1.0 under ZERO gDNA?

Builds single_exon (ss=1.0, NO gDNA, with a spliced helper), scans, and re-runs
the calibrate() loop with instrumentation focused on the exposure-shrinkage
question: the gDNA evidence (Σ M_g), the per-region exposure ω and its
uncertainty, the dispersion φ, and — for the gene regions — where any *false*
gDNA comes from (the E-step count vs strand log-Bayes-factors). Prints the
shrinkage arithmetic (1/φ prior pseudocount vs ρ₀·L data scale).
"""

from __future__ import annotations

import argparse
import tempfile

import numpy as np

from rigel.calibration.density import estimate_gdna_density
from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.estep import _llr_count, _llr_strand, _sense_flux, estep_view
from rigel.calibration.exposure import exposure_posterior
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.mstep import (
    fit_rho_d_bb,
    update_exposure_dispersion,
    update_pi_g_prior,
    update_rho_0,
)
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.sweep import sweep_ambig_exposure
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.splice import SpliceType

_TS = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gdna", type=float, default=0.0)
    ap.add_argument("--iters", type=int, default=6)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="zero_gdna_")
    sc = Scenario("single_exon", genome_length=8000, seed=42, work_dir=tmp)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1500)], "abundance": 100}])
    sc.add_gene(
        "g_helper", "+",
        [{"t_id": "t_helper", "exons": [(2500, 3000), (3500, 4000)], "abundance": 50}],
    )
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0}])
    cfg = ReadSimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=1.0, seed=42,
    )
    gdna = None if args.gdna == 0 else GDNAConfig(
        abundance=args.gdna, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000
    )
    result = sc.build_oracle(n_fragments=500, sim_config=cfg, gdna_config=gdna)
    index = result.index
    scan = BamScanConfig(sj_strand_tag="auto")
    _stats, strand_model, fl_scan, buffer, payload = scan_and_buffer(
        str(result.bam_path), index, scan
    )
    try:
        ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
        fl = build_fl_models(
            global_counts=fl_scan.global_model.counts,
            rna_counts=fl_scan.category_models[SpliceType.SPLICED_ANNOT].counts,
            gdna_counts=gdna_fl_mass(payload),
            max_size=fl_scan.max_size,
        )
        cfg_cal = CalibrationConfig()
        strand = fit_strand_balance(CalibrationSubstrate.from_payload(payload, ra), strand_model)
        print("=== strand model (from spliced channel) ===")
        print(f"  strand_specificity = {strand_model.strand_specificity:.4f}  "
              f"n_spliced_obs = {getattr(strand_model, 'n_observations', '?')}")
        print(f"  p_r1_sense = {strand_model.p_r1_sense:.4g}  read1_sense={strand_model.read1_sense}")
        print(f"  fit_strand_balance → kappa_rna={strand.kappa_rna:.4g}  rho_r_bb={strand.rho_r_bb:.4g}  "
              f"(n_obs={strand.n_observations}, n_reads={strand.n_total_reads}, "
              f"fallback={strand.fallback_used}:{strand.fallback_reason})")
        _trace(payload, ra, strand_model, fl.gdna_pmf, cfg_cal, args.iters)
    finally:
        buffer.cleanup()


def _trace(payload, ra, strand_model, gdna_fl_pmf, config, n_iters):
    sub = CalibrationSubstrate.from_payload(payload, ra)
    r = sub.n_regions
    ts_class = sub.ts_class
    rtype = coarse_type_array(ra.signature)
    l_phys = sub.L_eff
    nu = sub.contained.n_unspliced
    ns = sub.contained.n_spliced

    rho_0 = estimate_gdna_density(sub, ra).rho_0
    from rigel.calibration.calibrate import initial_hyperparameters

    disp, rho_d_bb = initial_hyperparameters(sub, config)
    strand = fit_strand_balance(sub, strand_model)
    kappa_rna, rho_r_bb = strand.kappa_rna, strand.rho_r_bb
    region_eff_len = region_eff_length(l_phys, gdna_fl_pmf)
    mu_fl = boundary_eff_length(gdna_fl_pmf)
    boundary_eff = np.full(r, mu_fl)

    # gene regions of interest (t1 exon ~500-1500, t_helper exons ~2500-4000)
    gene = np.where((nu > 0) | (ns > 0))[0]

    omega = np.ones(r)
    pi_prior = np.full(r, 0.5)
    m_d_c = np.zeros(r)
    m_d_l = np.zeros(r)
    m_d_r = np.zeros(r)

    for it in range(1, n_iters + 1):
        shared = dict(omega=omega, rho_0=rho_0, exposure_dispersion=disp, kappa_rna=kappa_rna,
                      rho_r_bb=rho_r_bb, rho_d_bb=rho_d_bb, pi_g_prior=pi_prior)
        cont = estep_view(sub.contained, ts_class, L_eff=region_eff_len, m_d_unspl_prev=m_d_c, **shared)
        left = estep_view(sub.left, ts_class, L_eff=boundary_eff, m_d_unspl_prev=m_d_l, **shared)
        right = estep_view(sub.right, ts_class, L_eff=boundary_eff, m_d_unspl_prev=m_d_r, **shared)
        exp = exposure_posterior(cont.m_g, left.m_g, right.m_g, rho_0=rho_0, L_eff=l_phys,
                                 exposure_dispersion=disp)
        swept = sweep_ambig_exposure(sub, ra, alloc_contained=cont, alloc_left=left, alloc_right=right,
                                     region_eff_len=region_eff_len, mu_fl=mu_fl, rho_0=rho_0,
                                     exposure_dispersion=disp, base_omega=exp.omega,
                                     base_log_omega_var=exp.log_omega_var)
        omega = swept.omega

        # E-step channel diagnostics on the contained view
        mu_g = omega * rho_0 * region_eff_len
        llr_c = _llr_count(nu.astype(np.float64), mu_g, m_d_c, 1.0 / disp)
        k_sense = _sense_flux(sub.contained, ts_class)
        llr_s = _llr_strand(k_sense, nu, ts_class, kappa_rna=kappa_rna, rho_r_bb=rho_r_bb, rho_d_bb=rho_d_bb)

        m_g_tot_total = float(exp.m_g_tot.sum())
        rho_0 = update_rho_0(exp.m_g_tot, omega, l_phys)
        disp = update_exposure_dispersion(exp.omega, exp.log_omega_var, floor=config.exposure_dispersion_floor)
        k_plus_g = np.maximum(cont.k_sense.astype(np.float64) - kappa_rna * cont.m_d_unspl, 0.0)
        rho_d_bb = fit_rho_d_bb(k_plus_g, cont.m_g_unspl)
        pi_prior = update_pi_g_prior(omega, rho_0, region_eff_len, nu)
        m_d_c, m_d_l, m_d_r = cont.m_d_unspl, left.m_d_unspl, right.m_d_unspl

        print(f"\n===== ITER {it} =====  rho_0={rho_0:.5g}  disp(phi)={disp:.4g}  "
              f"total_M_g={m_g_tot_total:.2f}  (gDNA evidence)")
        print("  region table (gene regions):  idx ts rtype start n_u n_s  pi_g  llr_c  llr_s  "
              "M_g(cont)  omega  sd(log w)")
        for i in gene:
            print(f"    {i:3d} {_TS.get(int(ts_class[i]), '?'):>5s} {_RTYPE.get(int(rtype[i]), '?'):>10s}"
                  f" {int(ra.start[i]):6d} {int(nu[i]):4d} {int(ns[i]):3d} {cont.pi_g[i]:.4f}"
                  f" {llr_c[i]:+.2f} {llr_s[i]:+.2f}  {cont.m_g[i]:7.2f}  {exp.omega[i]:7.3f}"
                  f"  {np.sqrt(exp.log_omega_var[i]):.2f}")

    # Shrinkage arithmetic for t1's exon region (the one that should be ω≈1).
    print("\n=== shrinkage arithmetic (final) ===")
    print("  ω = (1/φ + M_g) / (1/φ + ρ₀·L);  shrinks to 1.0 only if 1/φ ≫ ρ₀·L")
    print(f"  φ = {disp:.4g}  →  prior pseudocount 1/φ = {1.0 / disp:.4g}")
    for i in gene:
        if rtype[i] == 2 and ra.start[i] < 2000:  # t1's exon
            print(f"  t1 exon (idx {i}, L={int(l_phys[i])}): ρ₀·L = {rho_0 * l_phys[i]:.4g}  "
                  f"M_g = {(exp.m_g_tot[i]):.3f}  →  ω = {exp.omega[i]:.4f}")


if __name__ == "__main__":
    main()
