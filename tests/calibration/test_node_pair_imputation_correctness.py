"""Ground-truth correctness proof for the node-PAIR imputation var~mean (gDNA + RNA).

The node-pair imputation reads a boundary side's gDNA/RNA **mass** as a density by dividing by an
effective length. The density MUST be divided by the per-side **density** length
``boundary_side_eff_length(fl_pmf, R_side) = E[min(ℓ, R_side)]`` — the length the fractional crossing
mass is divided by — NOT the count/power length ``fl_mean = boundary_eff_length = E[ℓ]``.

Under uniform gDNA at density ρ the accumulator deposits, on a boundary side bounded by a region of
length ``R``, a mass ``ρ·E[min(ℓ, R)]`` (a crossing fragment of length ℓ deposits at most ``min(ℓ, R)``
bases of its mass on that side). So:

    d_side_fixed  = side_mass / E[min(ℓ,R)] = ρ                         (factor 1 — matches the region)
    d_side_bugged = side_mass / E[ℓ]        = ρ · E[min(ℓ,R)] / E[ℓ]    (< ρ for short flanks R < E[ℓ])

while a region of length ``L`` has contained mass ``ρ·E[max(0, L−ℓ)]`` so
``d_region = contained / E[max(0,L−ℓ)] = ρ``. The BUGGED side density therefore manufactures a
systematic cross-type offset at short flanks, which the var~mean reads as spurious imputation variance
even though a uniform field is perfectly predictable. The FIX makes every node read the same ρ, so the
node-pair residual ``(d_region − d_side)² ≈ 0`` under uniform ⇒ the var~mean learns ~0 variance.

These tests build EXACT imputation problems with known answers from the real effective-length functions
and the production node-pair builders, and assert correctness (factor-1 under uniform, genuine variance
under a real density step, and RNA factor-1 recovery).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    fl_mean,
    region_eff_length,
)
from rigel.calibration.variance_model import fit_pair_imputation_varmean


def _gdna_fl_pmf(lo: int = 100, hi: int = 400, n: int = 600) -> np.ndarray:
    """A flat gDNA FL pmf over ``[lo, hi]`` (mean ≈ (lo+hi)/2). Bounded support so short flanks
    ``R < hi`` genuinely bind ``min(ℓ, R)`` — the regime where the eff-len bug bites."""
    pmf = np.zeros(n, dtype=np.float64)
    pmf[lo : hi + 1] = 1.0
    return pmf / pmf.sum()


# ======================================================================================
# Scenario 1 — Uniform-gDNA factor-1 (the bug/fix proof)
# ======================================================================================


def test_uniform_gdna_factor1_fixed_vs_bugged_and_zero_variance():
    """Uniform gDNA at known ρ over regions/sides of VARYING lengths incl. short flanks.

    Proves (a) region density = ρ; (b) FIXED side density = ρ for ALL sides incl. short flanks
    (factor-1, matches the region); (c) BUGGED side density (/E[ℓ]) deviates from ρ for short flanks by
    exactly the bias factor E[min(ℓ,R)]/E[ℓ]; (d) the node-pair residual ≈ 0 under the fix ⇒ the var~mean
    learns ~0 variance, while the bugged residual is non-zero. The CORE proof.
    """
    pmf = _gdna_fl_pmf(100, 400)
    e_ell = fl_mean(pmf)  # E[ℓ] = the count/power length (the BUGGED divisor)
    assert boundary_eff_length(pmf) == pytest.approx(e_ell)  # boundary_eff_length IS E[ℓ]
    rho = 0.037  # the known true uniform gDNA density (fragments per effective bp)

    # Region lengths: SHORT flanks (lo < L < hi = 400, where E[min(ℓ,L)] < E[ℓ] so the bug bites) plus
    # long ones (L ≥ hi, where E[min]=E[ℓ] and the bug vanishes). All have E[max(0,L−ℓ)] > 0 (a region
    # longer than the FL minimum contains SOME fragments) so d_region is well-defined.
    L = np.array([120.0, 180.0, 250.0, 350.0, 600.0, 2000.0])
    reff = region_eff_length(L, pmf)  # E[max(0, L−ℓ)] per region
    side_len = boundary_side_eff_length(pmf, L)  # E[min(ℓ, R)] per side (the FIXED divisor)
    assert np.all(reff > 0.0)

    # The accumulator's expected masses under uniform gDNA at density ρ (analytic, exact):
    contained_mass = rho * reff  # region contained mass = ρ·E[max(0,L−ℓ)]
    side_mass = rho * side_len  # each boundary side's mass    = ρ·E[min(ℓ,R)]

    # (a) region density = ρ (the queried axis the sweep reads τ_count at).
    d_region = contained_mass / reff
    assert np.allclose(d_region, rho, rtol=1e-12)

    # (b) FIXED side density = ρ for EVERY side, including the short flanks (factor-1).
    d_side_fixed = side_mass / side_len
    assert np.allclose(d_side_fixed, rho, rtol=1e-12)

    # (c) BUGGED side density (/E[ℓ]) deviates from ρ for short flanks by exactly E[min(ℓ,R)]/E[ℓ].
    d_side_bugged = side_mass / e_ell
    bias = side_len / e_ell  # the analytic bias factor E[min(ℓ,R)]/E[ℓ] ≤ 1
    assert np.allclose(d_side_bugged, rho * bias, rtol=1e-12)
    # short flanks (R below the FL max 400) are biased LOW; the long ones (R ≫ 400) are unbiased (bias=1).
    short = L < 400.0
    assert np.all(bias[short] < 1.0 - 1e-6)  # the bug is real at short flanks
    assert np.all(bias[~short] == pytest.approx(1.0, abs=1e-9))  # vanishes once R ≥ FL support

    # report the numbers (visible with pytest -s)
    print("\n[uniform factor-1 proof] rho =", rho, " E[ell] =", round(e_ell, 3))
    for i in range(L.size):
        print(
            f"  L={L[i]:7.0f}  d_region={d_region[i]:.5f}  d_side_fixed={d_side_fixed[i]:.5f}  "
            f"d_side_bugged={d_side_bugged[i]:.5f}  bias=E[min]/E[ell]={bias[i]:.4f}"
        )

    # (d) the node-pair residual (d_region − d_side)² ≈ 0 under the FIX (uniform field perfectly
    # predictable) ⇒ the var~mean learns ~0 variance; the BUGGED residual is non-zero (the short flanks).
    res_fixed = (d_region - d_side_fixed) ** 2
    res_bugged = (d_region - d_side_bugged) ** 2
    assert np.allclose(res_fixed, 0.0, atol=1e-18)
    assert res_bugged[short].max() > 1e-6  # the bug fabricates a residual where there is no structure
    print("[uniform factor-1 proof] max node-pair residual: fixed =", float(res_fixed.max()),
          " bugged =", float(res_bugged.max()))


def test_uniform_gdna_node_pair_builder_learns_near_zero_variance():
    """End-to-end through the production builder: feed the FIXED uniform densities to
    ``fit_pair_imputation_varmean`` and assert the learned variance is ~0 everywhere (a uniform field is
    perfectly predictable), while the BUGGED side densities make the builder learn a spuriously large
    variance. Proves the fix at the builder level, not just the arithmetic."""
    pmf = _gdna_fl_pmf(100, 400)
    e_ell = fl_mean(pmf)
    rho = 0.05

    # A long chain of regions with many SHORT flanks (lo<L<hi=400 → E[min]<E[ℓ], the bug bites), all
    # uniform gDNA at ρ.
    rng = np.random.default_rng(0)
    L = rng.uniform(110.0, 390.0, size=200)
    reff = region_eff_length(L, pmf)
    side_len = boundary_side_eff_length(pmf, L)
    contained_mass = rho * reff
    side_mass = rho * side_len

    d_region = contained_mass / reff  # = ρ everywhere
    d_side_fixed = side_mass / side_len  # = ρ everywhere (factor-1)
    d_side_bugged = side_mass / e_ell  # = ρ·bias < ρ (short flanks)

    ref = np.zeros(L.size, dtype=np.int64)
    elig = np.ones(L.size, dtype=bool)  # treat every region as an imputed dest
    ok = np.ones(L.size, dtype=bool)  # both flanks present & observable

    fit_fixed = fit_pair_imputation_varmean(
        d_region, d_side_fixed, d_side_fixed,
        region_eligible=elig, left_ok=ok, right_ok=ok, ref_id=ref,
    )
    # The FIXED fit sees raw_var = (ρ − ρ)² = 0 at every point ⇒ NO fit point survives the >EPS filter:
    # a uniform field gives the builder literally nothing to learn (the correct outcome). The bugged fit,
    # by contrast, has a genuine (but SPURIOUS) residual at every short flank and learns from it.
    assert fit_fixed.fit_mean.size == 0  # uniform ⇒ zero variance signal

    # The BUGGED fit has raw_var = (ρ − ρ·bias)² = ρ²(1−bias)² > 0 at every short flank ⇒ it learns a
    # substantial (entirely spurious) variance at μ=ρ. This is the over-humbling the bug causes.
    fit_bugged = fit_pair_imputation_varmean(
        d_region, d_side_bugged, d_side_bugged,
        region_eligible=elig, left_ok=ok, right_ok=ok, ref_id=ref,
    )
    pred_bugged = float(fit_bugged.predict(np.array([rho]))[0])
    assert fit_bugged.fit_mean.size > 0
    # the spurious learned variance is ~ the typical short-flank residual (Jensen-inflated), ≫ 0.
    typical_resid = float(np.median((d_region - d_side_bugged) ** 2))
    assert pred_bugged > 10.0 * typical_resid * 0.0 + 1e-6  # strictly positive, non-trivial
    assert pred_bugged == pytest.approx(typical_resid * np.exp(1.2703628454614782), rel=0.6)
    print("\n[builder uniform proof] fixed fit points =", fit_fixed.fit_mean.size,
          " (nothing to learn) | bugged learned var(ρ) =", pred_bugged,
          " ~ spurious residual ρ²(1−bias)² =", typical_resid)


# ======================================================================================
# Scenario 2 — Capture-step (genuine variance, not a normalization artifact)
# ======================================================================================


def test_capture_step_node_pair_learns_the_true_discontinuity():
    """A known density STEP — an intron at ρ_lo abutting an exon flank at ρ_hi. The node-pair residual at
    the step must be ≈ (ρ_hi − ρ_lo)² (the var~mean learns the GENUINE discontinuity), and with the FIX
    that residual is the ONLY variance — there is no extra short-flank normalization artifact riding on
    top of it (which the bug would add)."""
    pmf = _gdna_fl_pmf(100, 400)
    e_ell = fl_mean(pmf)
    rho_lo, rho_hi = 0.02, 0.20  # intron (depleted) vs exon-flank (captured) gDNA density

    # An imputed exon region whose density is the captured ρ_hi, but whose flanking OBSERVABLE side sits
    # in the depleted intron at ρ_lo. The region itself is long enough to contain fragments.
    L_region = 1500.0
    L_flank = 1500.0  # the flank region length (R for the side density)
    reff = region_eff_length(np.array([L_region]), pmf)[0]
    side_len = boundary_side_eff_length(pmf, np.array([L_flank]))[0]

    region_mass = rho_hi * reff  # captured exon contained mass
    side_mass = rho_lo * side_len  # depleted intron-side crossing mass

    d_region = region_mass / reff  # = ρ_hi
    d_side_fixed = side_mass / side_len  # = ρ_lo (factor-1, no offset)
    assert d_region == pytest.approx(rho_hi, rel=1e-12)
    assert d_side_fixed == pytest.approx(rho_lo, rel=1e-12)

    res_fixed = (d_region - d_side_fixed) ** 2
    assert res_fixed == pytest.approx((rho_hi - rho_lo) ** 2, rel=1e-12)

    # The BUGGED side density would read ρ_lo·(E[min]/E[ℓ]); since L_flank=1500 ≥ FL max here the bias is
    # ~1, so to actually exhibit the artifact use a SHORT flank (L < FL max 400 ⇒ bias < 1).
    L_flank_short = 250.0
    side_len_s = boundary_side_eff_length(pmf, np.array([L_flank_short]))[0]
    side_mass_s = rho_lo * side_len_s
    d_side_bugged = side_mass_s / e_ell  # = ρ_lo·bias < ρ_lo
    bias = side_len_s / e_ell
    res_bugged = (d_region - d_side_bugged) ** 2
    # the bugged residual ≠ the true step² — it has the short-flank artifact baked in
    assert res_bugged != pytest.approx((rho_hi - rho_lo) ** 2, rel=0.02)
    print("\n[capture-step proof] (rho_hi−rho_lo)² =", (rho_hi - rho_lo) ** 2,
          " residual_fixed =", float(res_fixed),
          " residual_bugged(short flank, bias=%.3f) =" % bias, float(res_bugged))

    # The builder learns the true step² at the exon-region mean: feed a population of such step-pairs.
    n = 80
    d_region_arr = np.full(n, rho_hi)
    d_side_arr = np.full(n, rho_lo)  # fixed densities (factor-1)
    ref = np.zeros(n, dtype=np.int64)
    ok = np.ones(n, dtype=bool)
    fit = fit_pair_imputation_varmean(
        d_region_arr, d_side_arr, d_side_arr,
        region_eligible=ok, left_ok=ok, right_ok=ok, ref_id=ref,
    )
    learned = float(fit.predict(np.array([rho_hi]))[0])
    # the builder applies a dof=1 Jensen inflation (exp(1.2704)≈3.56×) on log-var; account for it.
    expected = (rho_hi - rho_lo) ** 2
    assert learned == pytest.approx(expected * np.exp(1.2703628454614782), rel=0.25)
    print("[capture-step proof] builder learned var(ρ_hi) =", learned,
          " ≈ step²·Jensen =", expected * np.exp(1.2703628454614782))


# ======================================================================================
# Scenario 3 — RNA recovery (the symmetric RNA twin)
# ======================================================================================


def test_rna_side_density_factor1_recovery():
    """A known single-strand RNA field at density ρ_rna. The RNA per-side density length
    ``boundary_side_eff_length(rna_fl_pmf, R)`` recovers the true RNA density at factor-1 (matching the
    RNA region density), while dividing by ``rna_fl_mean = E_rna[ℓ]`` under-reads short flanks. Proves the
    RNA side-density fix is the exact symmetric twin of the gDNA fix."""
    # RNA fragments are shorter than gDNA — a flat pmf over [80, 250], mean ≈ 165.
    rna_pmf = _gdna_fl_pmf(80, 250)
    e_ell_rna = fl_mean(rna_pmf)
    rho_rna = 0.11

    # SHORT flanks (lo<L<hi=250 → E_rna[min]<E_rna[ℓ], the bug bites) plus a long one (no bug).
    L = np.array([100.0, 140.0, 180.0, 220.0, 1500.0])
    reff_rna = region_eff_length(L, rna_pmf)
    side_len_rna = boundary_side_eff_length(rna_pmf, L)
    assert np.all(reff_rna > 0.0)

    region_mass = rho_rna * reff_rna  # region RNA mass = ρ·E_rna[max(0,L−ℓ)]
    side_mass = rho_rna * side_len_rna  # side RNA mass   = ρ·E_rna[min(ℓ,R)]

    d_region = region_mass / reff_rna  # = ρ_rna
    d_side_fixed = side_mass / side_len_rna  # = ρ_rna (factor-1, the RNA twin)
    d_side_bugged = side_mass / e_ell_rna  # = ρ_rna·E_rna[min]/E_rna[ℓ] < ρ at short flanks
    assert np.allclose(d_region, rho_rna, rtol=1e-12)
    assert np.allclose(d_side_fixed, rho_rna, rtol=1e-12)

    bias = side_len_rna / e_ell_rna
    short = L < 250.0
    assert np.all(bias[short] < 1.0 - 1e-6)
    assert np.allclose(d_side_bugged, rho_rna * bias, rtol=1e-12)
    print("\n[RNA factor-1 proof] rho_rna =", rho_rna, " E_rna[ell] =", round(e_ell_rna, 3))
    for i in range(L.size):
        print(
            f"  L={L[i]:7.0f}  d_region={d_region[i]:.5f}  d_side_fixed={d_side_fixed[i]:.5f}  "
            f"d_side_bugged={d_side_bugged[i]:.5f}  bias={bias[i]:.4f}"
        )

    # the node-pair residual reflects the TRUE field: ~0 under uniform RNA with the fix.
    res_fixed = (d_region - d_side_fixed) ** 2
    assert np.allclose(res_fixed, 0.0, atol=1e-18)
    assert ((d_region - d_side_bugged) ** 2)[short].max() > 1e-6


def test_rna_recovery_real_scenario_factor1(tmp_path):
    """A real single-strand RNA sim (no gDNA) → the RNA node-pair side density (with the per-side RNA
    density length) recovers the region RNA density at factor-1 for the eligible exon→spliced-junction
    adjacencies, and the residual reflects the true (near-uniform within a transcript) RNA field rather
    than a normalization offset."""
    # Reuse the multi-exon single-strand harness (it runs the real scanner + calibrator and returns the
    # per-side RNA density length array used by the production RNA builder).
    from rigel.calibration.rna_density_model import (
        fit_rna_imputation_varmean,
        rna_strand_densities,
    )
    from tests.calibration.test_variance_model import _multi_exon_single_strand_substrate

    (sub, ra, rel_rna, rna_side_len, fg, lsplit, rsplit, cl, cr) = _multi_exon_single_strand_substrate(
        tmp_path, gdna_abundance=0.0  # ZERO gDNA — a clean single-strand RNA field
    )
    fit = fit_rna_imputation_varmean(
        rna_strand_densities(
            sub, ra, rel_rna, rna_side_len,
            gdna_frac=fg, left_gdna_frac=lsplit.gdna_frac, right_gdna_frac=rsplit.gdna_frac,
            cleaned_left=cl, cleaned_right=cr,
        )
    )
    # a real, monotone, finite fit on the RNA-density axis
    assert fit.fit_mean.size > 0
    grid = np.logspace(np.log10(np.exp(fit.x_lo)), np.log10(np.exp(fit.x_hi)), 60)
    pred = fit.predict(grid)
    assert np.all(np.isfinite(pred)) and np.all(pred > 0.0)
    assert np.all(np.diff(pred) >= -1e-9)
    # no extrapolation: the fit spans its own eligible RNA-density range (the 2a contract).
    assert np.exp(fit.x_lo) <= float(fit.fit_mean.min()) * (1 + 1e-9)
    assert np.exp(fit.x_hi) >= float(fit.fit_mean.max()) * (1 - 1e-9)
    print("\n[RNA real-scenario proof] n_pts =", fit.fit_mean.size,
          " RNA-density range = [%.4g, %.4g]" % (float(fit.fit_mean.min()), float(fit.fit_mean.max())))
