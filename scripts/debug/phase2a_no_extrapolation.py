"""Phase-2a no-extrapolation check — the NEW IMPUTATION var~mean (fit_imputation_varmean_current).

Replicates the calibrate() setup for the capon (capture-ON + nascent) scenario from
phase1_m2_calibrate_internals.py, then, at the converged gDNA estimate, reports for BOTH the OLD
(fit_imputation_varmean, boundary-crossing axis) and NEW (fit_imputation_varmean_current, region
CURRENT-density axis) imputation models:

  * the model's fit mu-range [exp(x_lo), exp(x_hi)];
  * the EXON (non-count-observable, active) node current-density mu-range;
  * the fraction of exon nodes whose mu lies ABOVE the model's fit-range max (= flat-extrapolation).

The Phase-1 baseline was 71.0% extrapolation for the OLD model. The NEW model should bring the exon
mu INSIDE its fit range (extrapolation -> ~0%).

    python scripts/debug/phase2a_no_extrapolation.py [scenario=capon]
"""
import sys
import tempfile

import numpy as np

# reuse the EXACT scenario + setup from the Phase-1 internals harness
from phase1_m2_calibrate_internals import build_scenario  # noqa: E402

from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import node_gdna_density, count_observable_masks
from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
    overdispersion_for_beta,
)
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_deconv import cleaned_gdna_count, deconv_sides, strand_deconvolve
from rigel.calibration.splice_junction import region_splice_gdna_frac
from rigel.calibration.variance_model import (
    fit_direct_varmean,
    fit_imputation_varmean,
    fit_imputation_varmean_current,
    varmean_points,
)
from rigel.calibration.simplex_sweep import deconv_regions_sweep


def run(kind, pl, ra, sm, gdna_fl_pmf, rna_fl_pmf, cfg, R, gdna_present):
    substrate = CalibrationSubstrate.from_payload(pl, ra)
    region_eff_len = region_eff_length(ra.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)
    rna_fl_mean = boundary_eff_length(rna_fl_pmf)
    region_eff_len_rna = region_eff_length(ra.region_size_bp, rna_fl_pmf)

    rna_sense_frac = float(fit_strand_balance(sm).rna_sense_frac)
    nqv = float(cfg.gdna_deconv_quantile) != 0.5
    nd_raw = node_gdna_density(substrate, ra, region_eff_len, fl_mean, need_count_variance=nqv)
    gst = fit_gdna_strand_from_substrate(
        substrate, ra, nd_raw, boundary_eff_len, rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(cfg.gdna_strand_prior_alpha_beta),
        prior_weight=cfg.gdna_strand_prior_weight)
    god = gst.gdna_strand_overdispersion
    rod = fit_rna_strand_from_substrate(
        substrate, rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(cfg.rna_strand_prior_alpha_beta),
        prior_weight=cfg.rna_strand_prior_weight).rna_strand_overdispersion
    _, ls_, rs_ = strand_deconvolve(
        substrate, ra, rna_sense_frac=rna_sense_frac, gdna_strand_overdispersion=god,
        rna_strand_overdispersion=rod, deconv_quantile=cfg.gdna_deconv_quantile, n_grid=cfg.n_grid)
    i0 = cfg.gdna_strand_info_scale

    def rc(v):
        return v.n_unspliced_pos.astype(np.float64) + v.n_unspliced_neg.astype(np.float64)

    cl = cleaned_gdna_count(ls_, rc(substrate.left), i0)
    cr = cleaned_gdna_count(rs_, rc(substrate.right), i0)
    nd = node_gdna_density(substrate, ra, region_eff_len, fl_mean, need_count_variance=nqv,
                           gdna_counts=(rc(substrate.contained), cl, cr))
    rcf, _ = region_splice_gdna_frac(
        substrate, ra, nd.count_gdna_frac, eff_gdna=fl_mean, eff_rna=rna_fl_mean,
        eff_gdna_region=region_eff_len, eff_rna_region=region_eff_len_rna,
        left_gdna_unspl=cl, right_gdna_unspl=cr)
    left, right = deconv_sides(
        substrate, ra, nd, boundary_eff_len, rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=god, rna_strand_overdispersion=rod,
        deconv_quantile=cfg.gdna_deconv_quantile, n_grid=cfg.n_grid, info_scale=i0)

    c = substrate.contained
    u_tot = c.n_unspliced_pos.astype(np.float64) + c.n_unspliced_neg.astype(np.float64)
    eff_len = np.maximum(np.asarray(region_eff_len, dtype=np.float64), 1e-9)
    mass_u = np.maximum(np.asarray(c.mass_unspliced, dtype=np.float64), 1e-9)
    geom2 = (eff_len / mass_u) ** 2
    obs = np.asarray(nd.region_count_observable, dtype=bool)
    gdna_left = np.asarray(left.gdna_mass, dtype=np.float64)
    gdna_right = np.asarray(right.gdna_mass, dtype=np.float64)
    gdna_c = u_tot.copy()

    # run the loop to convergence (production semantics, including the new builder)
    prev = None
    for _ in range(int(cfg.sweep_max_passes)):
        rho = float(gdna_c[obs].sum() / max(float(eff_len[obs].sum()), 1e-9)) if obs.any() else 0.0
        pts = varmean_points(substrate, ra, region_eff_len, fl_mean,
                             gdna_views=(gdna_c, gdna_left, gdna_right))
        direct = fit_direct_varmean(pts)
        imp_new = fit_imputation_varmean_current(
            substrate, ra, region_eff_len, fl_mean, gdna_views=(gdna_c, gdna_left, gdna_right))
        mu = gdna_c / eff_len
        var_d = np.where(obs, direct.predict(mu), imp_new.predict(mu))
        cf = np.clip(rcf, 0.0, 1.0)
        sig2_frac = np.maximum(np.minimum(var_d * geom2, cf * (1.0 - cf)), 1e-12)
        sig2_glob = np.maximum(direct.predict(np.full_like(mu, rho)) * geom2, 1e-12)
        tau_count = np.minimum(1.0 / sig2_frac, mass_u)
        tau_global = np.clip(1.0 / sig2_glob, 1.0, mass_u)
        regions = deconv_regions_sweep(
            substrate, ra, rna_sense_frac=rna_sense_frac, gdna_strand_overdispersion=god,
            rna_strand_overdispersion=rod, count_gdna_frac=rcf, count_precision=tau_count,
            n_grid=cfg.sweep_n_grid, rho_global=rho, region_eff_len=region_eff_len,
            info_scale=i0, global_tau=tau_global)
        fg = np.asarray(regions.gdna_frac, dtype=np.float64)
        if prev is not None and float(np.mean(np.abs(fg - prev))) < cfg.sweep_convergence_delta:
            gdna_c = np.asarray(regions.gdna_mass, dtype=np.float64)
            break
        prev = fg
        gdna_c = np.asarray(regions.gdna_mass, dtype=np.float64)

    # at the converged estimate, compare OLD vs NEW imputation model fit-range coverage
    mu = gdna_c / eff_len
    sig = np.asarray(ra.signature)
    ref_id = np.asarray(ra.ref_id)
    rco, _ = count_observable_masks(sig, ref_id)
    exon = (~rco) & (u_tot > 0.0)
    ex_mu = mu[exon]

    pts = varmean_points(substrate, ra, region_eff_len, fl_mean,
                         gdna_views=(gdna_c, gdna_left, gdna_right))
    old = fit_imputation_varmean(pts)
    new = fit_imputation_varmean_current(
        substrate, ra, region_eff_len, fl_mean, gdna_views=(gdna_c, gdna_left, gdna_right))

    print(f"\n=== {kind}: NO-EXTRAPOLATION check (R={R}, exon nodes n={int(exon.sum())}) ===")
    print(f"  EXON current-density mu: min={ex_mu.min():.4g} med={np.median(ex_mu):.4g} "
          f"max={ex_mu.max():.4g}")
    for label, m in (("OLD fit_imputation_varmean (boundary-crossing axis)", old),
                     ("NEW fit_imputation_varmean_current (region CURRENT axis)", new)):
        lo, hi = np.exp(m.x_lo), np.exp(m.x_hi)
        above = float(np.mean(ex_mu > hi)) if ex_mu.size else float("nan")
        below = float(np.mean(ex_mu < lo)) if ex_mu.size else float("nan")
        inside = float(np.mean((ex_mu >= lo) & (ex_mu <= hi))) if ex_mu.size else float("nan")
        print(f"  {label}:")
        print(f"    fit mu-range [{lo:.4g}, {hi:.4g}]  (n_fit_points={m.fit_mean.size})")
        print(f"    exon mu ABOVE fit max (extrapolation): {100*above:5.1f}%  "
              f"BELOW: {100*below:4.1f}%  INSIDE: {100*inside:5.1f}%")


def main():
    kind = sys.argv[1] if len(sys.argv) > 1 else "capon"
    with tempfile.TemporaryDirectory() as work:
        sc, pl, ra, sm, gpmf, rpmf, cfg, R, gpres = build_scenario(kind, work)
        try:
            run(kind, pl, ra, sm, gpmf, rpmf, cfg, R, gpres)
        finally:
            try:
                sc.cleanup()
            except Exception:
                pass


if __name__ == "__main__":
    main()
