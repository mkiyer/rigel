"""Frame invariants of the pass-0 message layer, pinned against GROUND TRUTH rather than against goldens.

Every fact here was silently wrong in production at some point, and none of them was caught by the benchmark
— two were held in place by *compensating* errors and "validated" that way. The ground truth is the
accumulator's own deposit rule (``docs/accumulator/00_design.md`` §4.3, reference implementation
``tests/native/_accumulator_reference.py``), which is enumerable in closed form, so these tests do not depend
on any solver behaviour and cannot drift with it.

Derivations: ``docs/calibration/archive/message_layer_derivation.md`` §11 (A1/A2, the spliced effective length) and
the ``c_b`` denominator note in ``bp_solver`` (A3).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.effective_length import (
    boundary_side_eff_length,
    region_eff_length,
    spliced_side_eff_length,
)


def _delta_pmf(ell: int) -> np.ndarray:
    p = np.zeros(ell + 1, dtype=np.float64)
    p[ell] = 1.0
    return p


def _gauss_pmf(mu: float, sd: float, n: int = 1200) -> np.ndarray:
    x = np.arange(float(n))
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    return p / p.sum()


# ---------------------------------------------------------------------------
# A1 / A2 — the ONE-SIDED SPLICED effective length: which divisor, and when.
# ---------------------------------------------------------------------------


def _brute_spliced_deposit(ell: int, R: int, *, continues: bool) -> float:
    """Expected one-sided deposited mass per unit mature density, by direct enumeration of the accumulator
    rule. A spliced fragment credits only its exon flank; ``a`` = its bases on that side.

        a <= R : the near slice lies inside the flank    -> END slice,      n_cross=1 -> a/l
        a >  R : it overruns into the NEXT region        -> INTERIOR slice, n_cross=2 -> R/(2l)

    ``continues=False`` means the transcript ENDS at the flank's far edge, so no fragment can overrun and the
    ``a > R`` term never arises.
    """
    total = 0.0
    for a in range(1, ell):
        if a <= R:
            total += a / ell
        elif continues:
            total += R / (2 * ell)
    return total


@pytest.mark.parametrize("ell", [200, 300])
@pytest.mark.parametrize("R", [25, 50, 100, 150])
def test_spliced_efflen_continuing_exon_is_boundary_side_not_half_triangle(ell, R):
    """When the exon CONTINUES past the flank, the correct divisor is ``E[min(l,R)]/2`` — i.e.
    :func:`boundary_side_eff_length` — and the half-triangle is low by exactly ``l/R``.

    This is the A1 result. The half-triangle was production's divisor for EVERY face, which under-stated the
    spliced density (and so the mature peel) by up to 199x on short flanks."""
    pmf = _delta_pmf(ell)
    brute = _brute_spliced_deposit(ell, R, continues=True)
    two_sided = float(boundary_side_eff_length(pmf, np.array([float(R)]))[0])
    half_tri = float(spliced_side_eff_length(pmf, np.array([float(R)]))[0])

    assert brute == pytest.approx(two_sided, rel=1e-9), "E[min(l,R)]/2 must reproduce the deposit rule"
    # and the half-triangle is wrong here, by the derived factor l/R
    assert two_sided / half_tri == pytest.approx(ell / R, rel=1e-9)


@pytest.mark.parametrize("ell", [200])
@pytest.mark.parametrize("R", [25, 50, 100])
def test_spliced_efflen_terminating_exon_is_the_half_triangle(ell, R):
    """When the transcript ENDS at the flank's far edge no fragment can overrun, and the half-triangle is the
    correct divisor. This is why the selector must be structural rather than global."""
    pmf = _delta_pmf(ell)
    brute = _brute_spliced_deposit(ell, R, continues=False)
    half_tri = float(spliced_side_eff_length(pmf, np.array([float(R)]))[0])
    # brute sums a=1..R while the closed form integrates 0..R — a half-step of discretisation
    assert brute == pytest.approx(half_tri, abs=0.5 + 0.001 * half_tri)


def test_spliced_efflen_divisors_coincide_once_flank_exceeds_fl_support():
    """The selector is SELF-LIMITING: for ``R >> FL support`` the two divisors agree, so the choice only
    matters on short flanks — which is exactly where the geometry matters."""
    pmf = _gauss_pmf(200.0, 50.0)
    for R in (400.0, 1000.0, 5000.0):
        a = float(spliced_side_eff_length(pmf, np.array([R]))[0])
        b = float(boundary_side_eff_length(pmf, np.array([R]))[0])
        assert a == pytest.approx(b, rel=1e-3), f"divisors must coincide at R={R}"


# ---------------------------------------------------------------------------
# A3 — the c_b denominator. c_b = log1p(S_B/D_B) converts an exon's mature-INCLUSIVE composition into a
# boundary's mature-FREE crossing composition. Spliced mass lands on ONE face, so S_B must be divided by THAT
# face's opportunity, not by the sum of both faces'.
# ---------------------------------------------------------------------------


def test_cb_denominator_must_be_per_face_not_summed():
    """With known densities the true dilution is ``r = rho_mu/(rho_g+rho_nu)``. The PER-FACE denominator
    recovers it exactly; the SUMMED-EFF denominator (production's, before A3) under-states it ~2x, which
    under-peels the mature and biases f_g downward.

    Uses an ideal ``D_B`` so the test isolates the denominator question from the separate mixed-FL frame
    error in ``D_B`` (which divides a nascent-containing mass by the gDNA-FL opportunity)."""
    rp = _gauss_pmf(200.0, 50.0)  # the spliced channel is RNA, so only the RNA FL enters
    R_e, R_i = np.array([1000.0]), np.array([2000.0])
    rho_g, rho_nu, rho_mu = 0.5, 0.3, 1.0

    # one-sided spliced opportunity: the exon continues, so E[min(l,R)]/2 on the EXON face (A1/A2)
    esp_e = float(boundary_side_eff_length(rp, R_e)[0])
    esp_i = float(boundary_side_eff_length(rp, R_i)[0])
    spl = rho_mu * esp_e  # deposited on ONE face only

    r_true = rho_mu / (rho_g + rho_nu)
    r_per_face = (spl / esp_e) / (rho_g + rho_nu)
    r_summed = (spl / (esp_e + esp_i)) / (rho_g + rho_nu)

    assert r_per_face == pytest.approx(r_true, rel=1e-12), "per-face denominator must be exact"
    assert r_summed < 0.75 * r_true, "the summed denominator under-states S_B (it did so by ~2x)"
    # and in the log space the correction actually lives in:
    assert np.log1p(r_per_face) == pytest.approx(np.log1p(r_true), rel=1e-12)


def test_region_and_boundary_efflen_are_the_same_physical_quantity():
    """A region's CONTAINED opportunity and a boundary's CROSSING opportunity are two views of one thing —
    the chance a molecule of the local density is observed at that node. That is what makes a DENSITY
    comparable between node kinds, and hence usable as the message currency (§1).

    Sanity: both are positive, monotone in region length, and the crossing length saturates at fl_mean/2
    while the contained length grows without bound."""
    pmf = _gauss_pmf(200.0, 50.0)
    L = np.array([500.0, 1000.0, 5000.0, 20000.0])
    contained = region_eff_length(L, pmf)
    crossing = boundary_side_eff_length(pmf, L)
    assert np.all(np.diff(contained) > 0), "contained opportunity grows with length"
    assert np.all(contained > 0) and np.all(crossing > 0)
    fl_mean = float(np.dot(np.arange(pmf.shape[0], dtype=np.float64), pmf))
    assert crossing[-1] == pytest.approx(fl_mean / 2.0, rel=1e-3), "crossing saturates at fl_mean/2"
    assert contained[-1] > 10.0 * crossing[-1], "contained keeps growing; crossing does not"
