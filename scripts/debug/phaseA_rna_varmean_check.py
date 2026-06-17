"""Phase-A standalone curve check — the RNA `var~mean` reliability fit (CALIBRATION_PLAN_v5 §7).

Runs `calibrate()` to convergence on the capon (capture-on ss0.99 +nascent) and capoff scenarios
(reusing phase1_m2_calibrate_internals' scenario builders), then calls the NEW
`fit_imputation_rna_varmean` on the converged estimates and reports, for each scenario:

  * #fit points  vs the 18-pt spline-vs-power-law threshold (a real spline >= 18; activation-rate
    found ~52 pairs on these so expect a real spline).
  * whether the fit is monotone over a dense grid.
  * the fit mu-range [exp(x_lo), exp(x_hi)] vs the eligible regions' RNA-density range — the fit must
    SPAN the data (no extrapolation, the 2a no-extrapolation contract).
  * a few (mean, predicted var) samples across the fit range.

It must behave like the gDNA `fit_imputation_varmean_current` (sane, monotone, spans its data).

    python scripts/debug/phaseA_rna_varmean_check.py [scenario=all|capon|capoff]
"""
import importlib.util
import os
import sys
import tempfile

import numpy as np

from rigel.calibration.calibrate import calibrate
from rigel.calibration.effective_length import (
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, TS_NEG, TS_POS
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_deconv import cleaned_gdna_count, strand_deconvolve
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import count_observable_masks
from rigel.calibration.run_fill import same_ref_left_right
from rigel.calibration.rna_density_model import fit_rna_imputation_varmean, rna_strand_densities


def _load_phase1():
    here = os.path.dirname(__file__)
    spec = importlib.util.spec_from_file_location(
        "phase1_m2", os.path.join(here, "phase1_m2_calibrate_internals.py")
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _eligible_rna_density(substrate, ra, region_eff_len_rna, fg):
    """The CURRENT RNA density on s for the eligible single-strand exon regions (the fit's mean axis)."""
    sig = np.asarray(ra.signature).astype(np.int64)
    ts = np.asarray(ra.strand_class)
    ref_id = np.asarray(ra.ref_id)
    r = sig.shape[0]
    rco, bco = count_observable_masks(sig, ref_id)
    ls, rs = same_ref_left_right(ref_id)
    la = np.zeros(r, bool)
    rb = np.zeros(r, bool)
    if r > 1:
        la[1:] = bco[:-1] & ls[1:]
        rb[:-1] = bco[:-1] & rs[:-1]
    left_spl = np.asarray(substrate.left.n_spliced_sense, dtype=np.float64)
    right_spl = np.asarray(substrate.right.n_spliced_sense, dtype=np.float64)
    has_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    single = (ts == TS_POS) | (ts == TS_NEG)
    eligible = single & has_exon & (~rco) & (la & (left_spl > 0.0)) & (rb & (right_spl > 0.0))
    c = substrate.contained
    c_pos = np.asarray(c.n_unspliced_pos, dtype=np.float64)
    c_neg = np.asarray(c.n_unspliced_neg, dtype=np.float64)
    region_unspl_s = np.where(ts == TS_NEG, c_neg, c_pos)
    L = np.maximum(np.asarray(region_eff_len_rna, dtype=np.float64), 1e-12)
    dens = region_unspl_s * (1.0 - np.clip(fg, 0.0, 1.0)) / L
    return eligible, dens[eligible]


def run(kind, work, mod):
    sc, pl, ra, sm, gpmf, rpmf, cfg, R, gpres = mod.build_scenario(kind, work)
    print(f"\n{'='*90}\n=== PHASE-A RNA var~mean CHECK: {kind}  (R={R}, gDNA_present={gpres}) ===\n{'='*90}")

    # Run calibrate to convergence — gives the converged per-region gDNA mass.
    result = calibrate(pl, ra, sm, gpmf, rpmf, cfg)
    substrate = CalibrationSubstrate.from_payload(pl, ra)
    # f_g of the contained UNSPLICED mass = gdna_contained / mass_unspliced (exactly the sweep's f_g).
    mass_unspl = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    fg = np.where(
        mass_unspl > 1e-12,
        np.asarray(result.mass_gdna_contained, dtype=np.float64) / np.maximum(mass_unspl, 1e-12),
        0.0,
    )
    fg = np.clip(fg, 0.0, 1.0)

    region_eff_len_rna = region_eff_length(ra.region_size_bp, rpmf)
    rna_boundary_side_eff_len = boundary_side_eff_length(rpmf, ra.region_size_bp)

    # Recompute the per-side StrandSplit.gdna_frac + cleaned counts exactly as calibrate does.
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    od_g = result.gdna_strand_overdispersion
    od_r = result.rna_strand_overdispersion
    _, left_split, right_split = strand_deconvolve(
        substrate, ra, rna_sense_frac=kappa, gdna_strand_overdispersion=od_g,
        rna_strand_overdispersion=od_r, deconv_quantile=cfg.gdna_deconv_quantile, n_grid=cfg.n_grid,
    )
    i0 = cfg.gdna_strand_info_scale

    def _raw(view):
        return view.n_unspliced_pos.astype(np.float64) + view.n_unspliced_neg.astype(np.float64)

    cleaned_left = cleaned_gdna_count(left_split, _raw(substrate.left), i0)
    cleaned_right = cleaned_gdna_count(right_split, _raw(substrate.right), i0)

    fit = fit_rna_imputation_varmean(
        rna_strand_densities(
            substrate, ra, region_eff_len_rna, rna_boundary_side_eff_len,
            gdna_frac=fg, left_gdna_frac=left_split.gdna_frac, right_gdna_frac=right_split.gdna_frac,
            cleaned_left=cleaned_left, cleaned_right=cleaned_right,
        )
    )

    n_pts = fit.fit_mean.size
    is_spline = n_pts >= 18
    grid = np.logspace(np.log10(max(np.exp(fit.x_lo), 1e-12)),
                       np.log10(max(np.exp(fit.x_hi), 1e-11)), 200)
    pred = fit.predict(grid)
    monotone = bool(np.all(np.diff(pred) >= -1e-9))
    fit_lo, fit_hi = float(np.exp(fit.x_lo)), float(np.exp(fit.x_hi))

    eligible, elig_dens = _eligible_rna_density(substrate, ra, region_eff_len_rna, fg)
    elig_dens = elig_dens[np.isfinite(elig_dens) & (elig_dens > 1e-12)]

    print(f"  fit points = {n_pts}  ({'SPLINE (>=18)' if is_spline else 'POWER-LAW FALLBACK (<18)'})")
    print(f"  monotone over a 200-pt grid? {monotone}")
    print(f"  fit mu-range [exp(x_lo), exp(x_hi)] = [{fit_lo:.4g}, {fit_hi:.4g}]")
    if elig_dens.size:
        print(f"  eligible-region RNA-density range  = [{elig_dens.min():.4g}, {elig_dens.max():.4g}]"
              f"  (n_eligible={int(eligible.sum())})")
        spans = (fit_lo <= elig_dens.min() * (1 + 1e-9)) and (fit_hi >= elig_dens.max() * (1 - 1e-9))
        print(f"  fit SPANS the eligible RNA-density range (no extrapolation)? {spans}")
    else:
        print("  (no finite eligible RNA densities)")
    # samples
    samp_x = np.unique(np.clip(np.array([fit_lo, np.sqrt(fit_lo * fit_hi), fit_hi]), fit_lo, fit_hi))
    samp_y = fit.predict(samp_x)
    print("  (mean, predicted var) samples:")
    for mx, vy in zip(samp_x, samp_y):
        print(f"      mean={mx:10.4g}  var={vy:12.4g}")
    print(f"  [edf={fit.edf:.2f}  lam={fit.lam:.4g}]")
    sc.cleanup()


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    kinds = ["capon", "capoff"] if which == "all" else [which]
    mod = _load_phase1()
    for kind in kinds:
        with tempfile.TemporaryDirectory() as work:
            run(kind, work, mod)


if __name__ == "__main__":
    main()
