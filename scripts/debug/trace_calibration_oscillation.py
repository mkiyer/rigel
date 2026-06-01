#!/usr/bin/env python
"""Trace the calibration outer-loop oscillation in excruciating per-region detail.

Builds a tiny oracle scenario (the pipeline smoke fixture, which oscillates), scans
it, then re-runs the ``calibrate()`` outer loop verbatim with full per-iteration /
per-region instrumentation — the E-step count/strand log-Bayes-factors, ``π_g``, the
deconvolved masses, the exposure, and exactly what ``fit_phi`` / ``fit_rho_d_bb`` are
fitting at each cycle state. Diagnoses the period-2 limit cycle (φ flipping
floor ↔ ~4).

Usage:
    python scripts/debug/trace_calibration_oscillation.py [--iters 10] [--gdna 0.0]
"""

from __future__ import annotations

import argparse

import numpy as np

from rigel.calibration.density import estimate_gdna_density
from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.estep import estep_view
from rigel.calibration.exposure import exposure_posterior
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.mstep import (
    fit_rho_d_bb,
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
from rigel.sim import ReadSimConfig, Scenario
from rigel.splice import SpliceType

_TS = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}


def _nb_fit_dispersion(n_u, mu_g, m_d_unspl, floor, ceil=100.0):
    """The FORMER count-NB dispersion fit (kept only so the nb_* comparison modes
    still reproduce the old oscillating behaviour for the diagnosis)."""
    from scipy.optimize import minimize_scalar
    from scipy.stats import nbinom

    n = np.asarray(n_u, dtype=np.float64)
    mu = np.maximum(np.asarray(mu_g, np.float64) + np.asarray(m_d_unspl, np.float64), 1e-12)
    if float(n.sum()) <= 0.0:
        return floor

    def nll(phi):
        inv = 1.0 / phi
        return -float(np.sum(nbinom.logpmf(n, inv, inv / (inv + mu))))

    res = minimize_scalar(nll, bounds=(floor, ceil), method="bounded")
    return float(np.clip(res.x, floor, ceil))


def build_scenario(tmp_dir, *, gdna_fraction: float, n_fragments: int):
    sc = Scenario("osc", genome_length=5000, seed=42, work_dir=tmp_dir)
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50}])
    cfg = ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=42,
    )
    # gdna_fraction is a Scenario knob if supported; fall back to no-gDNA.
    try:
        return sc, sc.build_oracle(
            n_fragments=n_fragments, sim_config=cfg, gdna_fraction=gdna_fraction
        )
    except TypeError:
        return sc, sc.build_oracle(n_fragments=n_fragments, sim_config=cfg)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--iters", type=int, default=10)
    ap.add_argument("--gdna", type=float, default=0.0, help="gDNA fraction (if scenario supports)")
    ap.add_argument("--nfrag", type=int, default=500)
    ap.add_argument(
        "--phi-mode",
        default="nb_omega",
        choices=["nb_omega", "nb_active", "nb_marginal", "gamma"],
        help="phi M-step variant: nb_omega=current; nb_active=NB on n_u>0 only; "
        "nb_marginal=NB with mean rho_0*L_eff (no omega); gamma=Gamma-variance MLE from omega posteriors",
    )
    ap.add_argument("--quiet", action="store_true", help="trajectory only, no per-region dump")
    args = ap.parse_args()

    import tempfile

    tmp = tempfile.mkdtemp(prefix="osc_trace_")
    sc, result = build_scenario(tmp, gdna_fraction=args.gdna, n_fragments=args.nfrag)
    index = result.index
    scan = BamScanConfig(sj_strand_tag="auto")
    _stats, strand_model, fl_scan, buffer, payload = scan_and_buffer(
        str(result.bam_path), index, scan
    )
    try:
        region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
        fl = build_fl_models(
            global_counts=fl_scan.global_model.counts,
            rna_counts=fl_scan.category_models[SpliceType.SPLICED_ANNOT].counts,
            gdna_counts=gdna_fl_mass(payload),
            max_size=fl_scan.max_size,
        )
        _trace(
            payload,
            region_arrays,
            strand_model,
            fl.gdna_pmf,
            CalibrationConfig(),
            args.iters,
            args.phi_mode,
            args.quiet,
        )
    finally:
        buffer.cleanup()


def _gamma_variance_mstep(m_g_tot, rho_0, l_phys, phi_old):
    """Proper EM M-step for phi as the Gamma exposure-prior variance (s=1/phi).

    Given per-region omega posteriors Gamma(alpha_post, beta_post) with
    alpha_post = 1/phi_old + M_g, beta_post = 1/phi_old + rho_0*L_eff, solve the
    Gamma-shape MLE  log s - psi(s) = mean(E[omega] - E[log omega]) - 1  for s,
    then phi = 1/s. Uses only the exposure posteriors — no count NLL, so no
    empty-region count-misfit feedback.
    """
    from scipy.special import digamma

    inv = 1.0 / phi_old
    alpha_post = inv + np.asarray(m_g_tot, dtype=np.float64)
    beta_post = inv + rho_0 * np.asarray(l_phys, dtype=np.float64)
    e_omega = alpha_post / beta_post
    e_log_omega = digamma(alpha_post) - np.log(beta_post)
    c = float(np.mean(e_omega) - np.mean(e_log_omega)) - 1.0
    if not np.isfinite(c) or c <= 0.0:
        return phi_old
    # log s - psi(s) is positive, strictly decreasing on (0, inf): bisection.
    lo, hi = 1e-8, 1e8
    for _ in range(200):
        mid = (lo * hi) ** 0.5  # geometric bisection (s spans many decades)
        f = np.log(mid) - float(digamma(mid))
        if f > c:
            lo = mid
        else:
            hi = mid
    s = (lo * hi) ** 0.5
    return float(1.0 / s)


def _trace(payload, region_arrays, strand_model, gdna_fl_pmf, config, n_iters, phi_mode, quiet):
    sub = CalibrationSubstrate.from_payload(payload, region_arrays)
    r = sub.n_regions
    ts_class = sub.ts_class
    rtype = coarse_type_array(region_arrays.signature)
    l_phys = sub.L_eff

    rho_0 = estimate_gdna_density(sub, region_arrays).rho_0
    from rigel.calibration.calibrate import initial_hyperparameters

    phi, rho_d_bb = initial_hyperparameters(sub, config)
    strand = fit_strand_balance(sub, strand_model)
    kappa_rna, rho_r_bb = strand.kappa_rna, strand.rho_r_bb

    region_eff_len = region_eff_length(l_phys, gdna_fl_pmf)
    mu_fl = boundary_eff_length(gdna_fl_pmf)
    boundary_eff = np.full(r, mu_fl)

    print(
        f"# regions={r}  seed-rho_0={rho_0:.6g}  init phi={phi:.4g} rho_d_bb={rho_d_bb:.4g} "
        f"kappa_rna={kappa_rna:.4g} rho_r_bb={rho_r_bb:.4g}"
    )
    print(f"# mu_FL(boundary)={mu_fl:.2f}  phi_floor={config.exposure_dispersion_floor:.3g}")
    # Which regions carry contained unspliced/spliced flux?
    nu = sub.contained.n_unspliced
    ns = sub.contained.n_spliced
    active = np.where((nu > 0) | (ns > 0))[0]
    print(
        f"# contained-active regions: {len(active)}  total n_u={int(nu.sum())} n_s={int(ns.sum())}"
    )

    omega = np.ones(r)
    pi_g_prior = np.full(r, 0.5)
    m_d_cont = np.zeros(r)
    m_d_left = np.zeros(r)
    m_d_right = np.zeros(r)
    m_g_tot_prev = None

    for it in range(1, n_iters + 1):
        shared = dict(
            omega=omega,
            rho_0=rho_0,
            exposure_dispersion=phi,
            kappa_rna=kappa_rna,
            rho_r_bb=rho_r_bb,
            rho_d_bb=rho_d_bb,
            pi_g_prior=pi_g_prior,
        )
        cont = estep_view(
            sub.contained, ts_class, L_eff=region_eff_len, m_d_unspl_prev=m_d_cont, **shared
        )
        left = estep_view(sub.left, ts_class, L_eff=boundary_eff, m_d_unspl_prev=m_d_left, **shared)
        right = estep_view(
            sub.right, ts_class, L_eff=boundary_eff, m_d_unspl_prev=m_d_right, **shared
        )
        exposure = exposure_posterior(
            cont.m_g, left.m_g, right.m_g, rho_0=rho_0, L_eff=l_phys, exposure_dispersion=phi
        )
        swept = sweep_ambig_exposure(
            sub,
            region_arrays,
            alloc_contained=cont,
            alloc_left=left,
            alloc_right=right,
            region_eff_len=region_eff_len,
            mu_fl=mu_fl,
            rho_0=rho_0,
            exposure_dispersion=phi,
            base_omega=exposure.omega,
            base_log_omega_var=exposure.log_omega_var,
        )
        omega = swept.omega

        # --- M-step (verbatim from calibrate.py) ---
        rho_0_new = update_rho_0(exposure.m_g_tot, omega, l_phys)
        mu_g_cont = omega * rho_0_new * region_eff_len
        floor = config.exposure_dispersion_floor
        if phi_mode == "nb_omega":
            phi_new = _nb_fit_dispersion(nu, mu_g_cont, cont.m_d_unspl, floor)
        elif phi_mode == "nb_active":
            act = nu > 0
            phi_new = _nb_fit_dispersion(nu[act], mu_g_cont[act], cont.m_d_unspl[act], floor)
        elif phi_mode == "nb_marginal":
            phi_new = _nb_fit_dispersion(nu, rho_0_new * region_eff_len, cont.m_d_unspl, floor)
        else:  # gamma — proper exposure-variance M-step from the omega posteriors
            phi_new = _gamma_variance_mstep(exposure.m_g_tot, rho_0, l_phys, phi)
        k_plus_g = np.maximum(cont.k_sense.astype(np.float64) - kappa_rna * cont.m_d_unspl, 0.0)
        rho_d_bb_new = fit_rho_d_bb(k_plus_g, cont.m_g_unspl)
        pi_g_prior_new = update_pi_g_prior(omega, rho_0_new, region_eff_len, nu)

        m_g_tot = exposure.m_g_tot
        if m_g_tot_prev is None:
            delta = float(m_g_tot.max())
        else:
            delta = float(np.max(np.abs(m_g_tot - m_g_tot_prev) / (m_g_tot_prev + 1.0)))

        if quiet:
            print(
                f"iter {it:2d}  phi_in={phi:10.5g}  ->  phi_out={phi_new:10.5g}   "
                f"rho_0={rho_0_new:.5g}  delta={delta:.4g}"
            )
            rho_0, phi, rho_d_bb = rho_0_new, phi_new, rho_d_bb_new
            pi_g_prior = pi_g_prior_new
            m_d_cont, m_d_left, m_d_right = cont.m_d_unspl, left.m_d_unspl, right.m_d_unspl
            m_g_tot_prev = m_g_tot.copy()
            continue

        print(
            f"\n===== ITER {it}  (E-step uses phi={phi:.4g}, rho_d_bb={rho_d_bb:.4g}) "
            f"delta={delta:.4g} ====="
        )
        print(
            f"  M-step OUT: rho_0={rho_0_new:.5g}  phi={phi_new:.5g}  rho_d_bb={rho_d_bb_new:.4g}"
        )
        # Per-region table — EVERY region that contributes to fit_phi (n_u>0 OR a
        # nonzero modeled gDNA mean μ_g). μ_tot = μ_g + m_d is the NB mean.
        mu_tot_all = mu_g_cont + cont.m_d_unspl
        contrib = np.where((nu > 0) | (mu_g_cont > 0.05))[0]
        print(
            "  region table:  idx rtype       ts     L_eff   n_u    omega    mu_g   m_d_uns   mu_tot   pi_g"
        )
        for i in contrib:
            print(
                f"    {i:3d} {_RTYPE.get(int(rtype[i]), '?'):>10s} {_TS.get(int(ts_class[i]), '?'):>5s}"
                f" {region_eff_len[i]:8.1f} {nu[i]:5d} {omega[i]:7.3f} {mu_g_cont[i]:7.3f}"
                f" {cont.m_d_unspl[i]:8.3f} {mu_tot_all[i]:8.3f}  {cont.pi_g[i]:.4f}"
            )
        # Decompose the fit_phi objective by region at phi=floor vs phi=4.
        _phi_nll_curve(nu, mu_g_cont, cont.m_d_unspl, config.exposure_dispersion_floor)
        _phi_nll_by_region(
            nu, mu_g_cont, cont.m_d_unspl, config.exposure_dispersion_floor, contrib, rtype
        )

        rho_0, phi, rho_d_bb = rho_0_new, phi_new, rho_d_bb_new
        pi_g_prior = pi_g_prior_new
        m_d_cont, m_d_left, m_d_right = cont.m_d_unspl, left.m_d_unspl, right.m_d_unspl
        m_g_tot_prev = m_g_tot.copy()


def _phi_nll_curve(n_u, mu_g, m_d_unspl, phi_floor):
    """Print the NB NLL as a function of phi (the fit_phi objective) on a grid."""
    from scipy.stats import nbinom

    mu = np.maximum(mu_g + m_d_unspl, 1e-12)
    n = n_u.astype(np.float64)

    def nll(phi):
        inv_phi = 1.0 / phi
        p = inv_phi / (inv_phi + mu)
        return -float(np.sum(nbinom.logpmf(n, inv_phi, p)))

    grid = [phi_floor, 1e-4, 1e-2, 0.1, 0.5, 1.0, 2.0, 4.0, 10.0, 50.0, 100.0]
    vals = [(g, nll(g)) for g in grid]
    best = min(vals, key=lambda gv: gv[1])
    cells = "  ".join(f"{g:g}:{v:.2f}" for g, v in vals)
    print(f"  NB-NLL(phi): {cells}   argmin~{best[0]:g}")


def _phi_nll_by_region(n_u, mu_g, m_d_unspl, phi_floor, contrib, rtype):
    """Per-region NB-NLL contribution at phi=floor vs phi=4 — reveals which
    regions drive the fit toward small vs large phi (Δ = NLL@4 − NLL@floor;
    Δ<0 ⇒ region prefers large phi, Δ>0 ⇒ prefers small phi)."""
    from scipy.stats import nbinom

    mu = np.maximum(mu_g + m_d_unspl, 1e-12)
    n = n_u.astype(np.float64)

    def per_region(phi):
        inv_phi = 1.0 / phi
        p = inv_phi / (inv_phi + mu)
        return -nbinom.logpmf(n, inv_phi, p)

    nll_lo = per_region(phi_floor)
    nll_hi = per_region(4.0)
    print("  per-region NLL:  idx rtype       n_u    mu_tot   NLL@floor  NLL@4    Δ(4-floor)")
    for i in contrib:
        print(
            f"    {i:3d} {_RTYPE.get(int(rtype[i]), '?'):>10s} {n_u[i]:5d} {mu[i]:8.3f}"
            f"  {nll_lo[i]:9.3f} {nll_hi[i]:8.3f}  {nll_hi[i] - nll_lo[i]:+9.3f}"
        )
    print(
        f"    TOTAL Δ(4-floor) over ALL regions = {float((nll_hi - nll_lo).sum()):+.3f}  "
        f"(negative ⇒ fit_phi picks large phi)"
    )


if __name__ == "__main__":
    main()
