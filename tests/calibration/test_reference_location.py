"""ψ's REFERENCE HAS A MEAN, AND UNTIL NOW NOBODY CHOSE IT — the term that lets it be set, pinned.

``a·log f_g + b·log(1−f_g)`` on the λ grid is exactly ``Beta(a, b)`` in ``f_g``, so ``a`` and ``b`` are
PSEUDO-COUNTS with a STRENGTH ``a + b`` and a MEAN ``a/(a+b)``. Jeffreys fixes both, and its mean of ½
**asserts the library is half gDNA**. :func:`~rigel.calibration.simplex_logodds._location_term` is the
term that moves the mean while leaving the strength — and, critically, the TAILS — alone::

    log p_i(λ) = ½·log f_g + ½·log(1−f_g) − log[ (1−m_i)·f_g + m_i·(1−f_g) ]

⛔⛔ **THE ALTERNATIVE — moving ``a`` and ``b`` — IS WHAT THESE GATES EXIST TO RULE OUT.** The exponents
set the tail slopes and the location TOGETHER, so buying a location costs tail mass: at ``b = 0.03``,
**57 %** of the prior's mass falls outside the shipped ``L = 10`` window and the answer becomes a function
of the grid. Every gate below that mentions ``L`` or a tail is testing that this term does NOT do that.

⭐ Each gate carries its own perturbation: an assertion that cannot fail is not a gate
(`TRAPS: perturb-every-gate`).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.simplex_logodds import (
    _JEFFREYS_REF,
    _gdna_arm,
    _location_term,
    _log1m_fg,
    _log_fg,
    _logodds_grid,
    _solve_regions_logodds_all,
    CompositionPriors,
)

#: the shipped lattice, and two others so nothing below can be a property of one grid
_GRIDS = (60, 256)
_L = 10.0


def _reference(lam, m):
    """ψ's reference alone, as a normalised density over the λ grid."""
    lp = _JEFFREYS_REF * _log_fg(lam) + _JEFFREYS_REF * _log1m_fg(lam)
    if m is not None:
        lp = lp + _location_term(lam, np.array([m]))[0]
    lp = lp - lp.max()
    p = np.exp(lp)
    return p / np.trapezoid(p, lam)


def _median(lam, p):
    c = np.cumsum(p) * (lam[1] - lam[0])
    return float(1.0 / (1.0 + np.exp(-lam[np.searchsorted(c / c[-1], 0.5)])))


# ── the reduction: this is a strict GENERALISATION, not a replacement ───────────────────────────────


@pytest.mark.parametrize("n_grid", _GRIDS)
def test_m_one_half_is_the_shipped_constant_exactly(n_grid):
    """⭐ ``m = ½`` makes the bracket the CONSTANT ½, so the term is constant and cancels in ψ's
    normalisation. The rule agrees with the shipped behaviour precisely where the shipped assumption
    happens to be true, which is the strongest internal consistency check available."""
    lam, _ = _logodds_grid(n_grid, _L)
    term = _location_term(lam, np.array([0.5, 0.5, 0.5]))
    assert term.shape == (3, n_grid)
    assert np.allclose(term, term[0, 0], rtol=0.0, atol=1e-15)
    # ⛔ the perturbation: any other m must NOT be constant, or the term is inert and moves nothing
    assert not np.allclose(
        _location_term(lam, np.array([0.2])), _location_term(lam, np.array([0.2]))[0, 0]
    )


def test_none_writes_no_term_at_all():
    """⛔ ``None`` is the bit-identity contract the whole stage rests on: not "a neutral term" but NO
    term, so every caller reproduces the pre-existing path exactly."""
    lam, _ = _logodds_grid(60, _L)
    assert _location_term(lam, None) == 0.0
    assert isinstance(_location_term(lam, None), float)
    # ⛔ and the arm it is added beside is untouched by it
    assert np.array_equal(_gdna_arm(lam, None), _gdna_arm(lam, None) + _location_term(lam, None))


# ── the median identity, in closed form ────────────────────────────────────────────────────────────


@pytest.mark.parametrize("m", (0.02, 0.2, 0.5, 0.8, 0.98))
def test_the_prior_median_is_m_in_closed_form(m):
    """⭐⭐ ``u = f/(f + r(1−f)) ~ Beta(a,b)`` is symmetric at ``a = b``, so the median sits at ``u = ½``,
    which is ``f = m``. Measured error **4e-6 … 5e-5**. **This is what re-derives
    ``test_relay_mass_rescale``'s hard-coded ``R_own = 0.5`` rather than widening it**: the neutral default
    stops being a number and becomes a formula whose ``m = ½`` case is that number.

    ⚠ **On the UNTRUNCATED prior.** Inside the shipped ``L = 10`` window the median is displaced by up to
    **0.0019** at the extreme locations, because truncation removes more mass from the tilted side —
    that is a property of the window, is pinned by :func:`test_l_invariance_the_answer_does_not_move`,
    and is two orders of magnitude smaller than the effect this term exists to produce.
    """
    lam, _ = _logodds_grid(200001, 60.0)
    assert _median(lam, _reference(lam, m)) == pytest.approx(m, abs=1e-3)
    # ⛔ the perturbation: the SHIPPED reference has median ½ whatever m we wanted
    assert _median(lam, _reference(lam, None)) == pytest.approx(0.5, abs=1e-6)


# ── the tails, which is the whole reason this is a third term and not a different (a, b) ────────────


@pytest.mark.parametrize("m", (0.01, 0.2, 0.5, 0.9, 0.99))
def test_the_tail_slope_stays_at_the_jeffreys_exponent_for_every_m(m):
    """⭐⭐⭐ **ONLY THE LOCATION MOVES.** In the ``+λ`` tail the bracket → ``m·(1−f_g)`` and the term
    contributes ``−log m − log(1−f_g)``; the ``log(1−f_g)`` cancels ONE of the two Jeffreys halves, so
    ψ's reference decays as ``e^(−½|λ|)`` at BOTH ends for every ``m``. ⛔ That is exactly what moving
    ``a``/``b`` cannot do, and it is why ``L``-invariance survives this change."""
    lam, _ = _logodds_grid(40001, 30.0)
    lp = np.log(_reference(lam, m))
    for sl in (slice(0, 200), slice(-200, None)):
        slope = np.polyfit(lam[sl], lp[sl], 1)[0]
        assert abs(abs(slope) - _JEFFREYS_REF) < 1e-3, (m, slope)


@pytest.mark.parametrize("m", (0.01, 0.2, 0.5, 0.9, 0.99))
def test_l_invariance_the_answer_does_not_move(m):
    """⭐⭐ The module's stated acceptance test, tested on the ANSWER rather than on a proxy: the prior's
    median at ``L`` = 10 / 20 / 40 must agree. Measured, this term moves **0.0019** at the extreme
    locations and **exactly 0** at ``m = ½``.

    ⛔⛔ **The perturbation is the design that was REFUSED.** Buying the same location by moving the
    exponents instead moves the answer **0.021** (``b = 0.03``) to **0.041** (``b = 0.25``) across the
    same windows — 11–22x more — because the exponents set the TAILS, and the window then decides how
    much of the prior exists. That is the whole reason this is a third term and not a different ``(a,b)``.
    """
    med = [
        _median(*(lambda g: (g[0], _reference(g[0], m)))(_logodds_grid(200001, L)))
        for L in (10.0, 20.0, 40.0)
    ]
    assert max(med) - min(med) < 5e-3, (m, med)

    # ⛔ the (a, b) route, on the same windows and the same estimator
    def _ab(L, b):
        lam, _ = _logodds_grid(200001, L)
        lp = _JEFFREYS_REF * _log_fg(lam) + b * _log1m_fg(lam)
        lp -= lp.max()
        q = np.exp(lp)
        return _median(lam, q / np.trapezoid(q, lam))

    for b in (0.03, 0.25):
        moved = [_ab(L, b) for L in (10.0, 20.0, 40.0)]
        assert max(moved) - min(moved) > 5.0 * (max(med) - min(med) + 1e-9), (b, moved, med)


@pytest.mark.parametrize("m", (0.01, 0.2, 0.5, 0.9, 0.99))
def test_the_mass_outside_the_window_stays_bounded(m):
    """⭐ Because the tails are untouched, the reference's mass outside ``L = 10`` stays **0.009–0.043**
    across three decades of location — it grows mildly with the tilt and never becomes the answer.
    ⛔ The refused ``(a, b)`` route at ``b = 0.03`` leaves **0.658** outside the same window, so the gain
    there could not be distinguished from prior mass piling at the grid edge."""
    wide, _ = _logodds_grid(200001, 60.0)
    outside = lambda p: float(np.trapezoid(p * (np.abs(wide) > _L), wide))  # noqa: E731
    here = outside(_reference(wide, m))
    lp = _JEFFREYS_REF * _log_fg(wide) + 0.03 * _log1m_fg(wide)
    lp -= lp.max()
    q = np.exp(lp)
    assert here < 10.0 * outside(_reference(wide, 0.5)), (m, here)
    assert outside(q / np.trapezoid(q, wide)) > 10.0 * here


@pytest.mark.parametrize("m", (1e-6, 0.5, 1.0 - 1e-6))
def test_psi_stays_proper_at_every_m_including_the_clamped_ends(m):
    """⛔ Properness is the property the 2026-08-15 ablation LACKED and is why the zero control doubled
    there. The term keeps it for every ``m`` in (0,1), and the clamp is what makes the closed endpoints
    representable rather than a special case."""
    lam, _ = _logodds_grid(40001, 60.0)
    p = _reference(lam, m)
    assert np.all(np.isfinite(p))
    assert float(np.trapezoid(p, lam)) == pytest.approx(1.0, rel=1e-6)
    # ⛔ the perturbation: Haldane (both exponents → 0) is improper and the same integral runs away
    flat = np.ones_like(lam)
    assert float(np.trapezoid(flat, lam)) > 100.0


def test_the_term_is_bounded_and_survives_the_float32_cast():
    """⚠ The AMBIG path sums the arms, casts to float32 and broadcasts across the tilt. A term with
    float64-only dynamic range would work on single-strand objects and silently degrade on AMBIG ones —
    the half of the panel where the answer is already worst. This term's range is
    ``log(max(m,1−m)/min(m,1−m))``, i.e. bounded."""
    lam, _ = _logodds_grid(256, _L)
    m = np.array([1e-5, 0.5, 1 - 1e-5])
    term = _location_term(lam, m)
    assert np.all(np.isfinite(term))
    assert np.abs(term).max() < 25.0
    lost = np.abs(term.astype(np.float32).astype(np.float64) - term).max()
    assert lost < 1e-4, lost


# ── the monotone direction, so a sign error cannot pass ─────────────────────────────────────────────


def test_a_higher_m_moves_prior_mass_toward_the_gdna_vertex():
    """⛔ A sign error here would be invisible in every gate above — the median identity would still hold
    with ``m`` and ``1−m`` swapped only if the reference were symmetric, which it is. This pins the
    DIRECTION against the meaning: ``m`` is the expected gDNA share, so raising it must move mass toward
    ``f_g = 1``."""
    lam, fg = _logodds_grid(20001, _L)
    frac_high = lambda m: float(  # noqa: E731
        np.trapezoid(_reference(lam, m) * (fg > 0.9), lam)
    )
    lo, mid, hi = frac_high(0.1), frac_high(0.5), frac_high(0.9)
    assert lo < mid < hi, (lo, mid, hi)
    # ⛔ and the U-shape SURVIVES: both vertices keep mass at every m — a tilted two-spike prior, not a
    #   bump in the middle, which is what a spike-and-slab population actually wants.
    for m in (0.1, 0.9):
        assert float(np.trapezoid(_reference(lam, m) * (fg < 0.01), lam)) > 0.005
        assert float(np.trapezoid(_reference(lam, m) * (fg > 0.99), lam)) > 0.005


# ── the dataclass plumbing, where the bug class lives ───────────────────────────────────────────────


def test_select_slices_the_location_with_the_arms():
    """⛔ ``select`` restricts every member to one solve block's rows. A location left unsliced would
    broadcast against the wrong rows — the exact bug class this dataclass exists to make
    unrepresentable."""
    p = CompositionPriors(
        gdna=np.arange(12.0).reshape(4, 3), rna=None, location=np.array([0.1, 0.2, 0.3, 0.4])
    )
    q = p.select(np.array([True, False, True, False]))
    assert np.array_equal(q.location, np.array([0.1, 0.3]))
    assert q.location.shape[0] == q.gdna.shape[0]
    assert q.rna is None


def test_regrid_carries_the_arms_and_leaves_the_location_alone():
    """⭐ The location is a per-slot SCALAR, so it has no lattice to be carried between and
    ``regrid-in-the-right-variable`` cannot apply to it. Storing it as a scalar rather than a
    pre-evaluated ``(m, K)`` grid is what makes that true by construction rather than by discipline."""
    loc = np.array([0.25, 0.75])
    p = CompositionPriors(gdna=np.zeros((2, 60)), rna=None, location=loc)
    q = p.regrid(60, 256, _L)
    assert q.gdna.shape == (2, 256)
    assert q.location is loc
    # ⛔ equal grids stay a passthrough on every member
    r = p.regrid(60, 60, _L)
    assert r.gdna is p.gdna and r.location is loc


def test_defaults_are_the_shipped_prior_free_solve():
    """``CompositionPriors()`` is pass-0: no fitted arm, no location, and therefore no new term."""
    p = CompositionPriors()
    assert p.gdna is None and p.rna is None and p.location is None


# ── both ψ paths, end to end: the term is inert at None and LIVE otherwise ──────────────────────────


def _solve(priors, n_grid=256, scale=1.0):
    """One single-strand object and one AMBIG object through the real dispatcher.

    ⭐ **Both ψ call sites in one call.** The single-strand path adds the arms on the fine λ grid; the
    AMBIG path sums them, casts to float32 and broadcasts across the tilt. A term wired into one and not
    the other would pass every scalar gate above and silently do nothing on half the panel.

    ``scale`` multiplies every count, which is how the gates below vary the ONE thing that decides
    whether a one-pseudo-fragment prior is audible: how much real data the object has.
    """
    u_pos = np.array([4.0, 2.5]) * scale
    u_neg = np.array([0.0, 1.5]) * scale
    allow_pos = np.array([True, True])
    allow_neg = np.array([False, True])  # object 0 single-strand, object 1 AMBIG
    return np.asarray(
        _solve_regions_logodds_all(
            u_pos,
            u_neg,
            allow_pos,
            allow_neg,
            np.array([6.0, 5.5]) * scale,
            np.array([0.0, 0.0]),
            kappa=0.9,
            od_g=0.0,
            od_r=0.0,
            n_grid=n_grid,
            priors=priors,
        ).gdna_frac,
        np.float64,
    )


def test_a_none_location_is_bit_identical_through_both_psi_paths():
    """⛔⛔ **THE STAGE-1 CONTRACT.** ``CompositionPriors()`` — the pass-0 default — must reproduce the
    path that existed before this term did, on BOTH ψ solves, exactly.

    ⭐ Proven on real data as well: one ``calibrate`` over ``g50 ss0.99 capture_off`` with the term
    present-but-``None`` against the same run with ``_location_term`` removed entirely is **19/19 output
    arrays bit-identical**, and forcing a nonzero term moves **10/19** with ``max|Δ| = 48.47``. That
    second half is what stops this gate being vacuous, and it is reproduced below in-process.
    """
    for scale in (1.0, 40.0):
        base = _solve(None, scale=scale)
        assert np.array_equal(_solve(CompositionPriors(), scale=scale), base)
        assert np.array_equal(
            _solve(CompositionPriors(gdna=None, rna=None, location=None), scale=scale), base
        )


def test_a_set_location_moves_both_psi_paths():
    """⛔ The perturbation the gate above needs: if a set location changed nothing, the term would be
    unreachable and every number in this file would be about dead code. ⭐ Both objects must move — the
    single-strand one AND the AMBIG one, which is the float32 path — and in the DIRECTION the meaning
    demands, since ``m`` is the expected gDNA share."""
    up = _solve(CompositionPriors(location=np.array([0.98, 0.98])))
    down = _solve(CompositionPriors(location=np.array([0.02, 0.02])))
    assert np.all(up > down), (up, down)
    base = _solve(None)
    for i, kind in enumerate(("single-strand", "AMBIG")):
        assert up[i] > base[i] > down[i], (kind, down[i], base[i], up[i])


def test_evidence_swamps_the_location_where_evidence_exists():
    """⭐⭐⭐ **THE STRENGTH PROPERTY, MEASURED — and it is the reason the location cannot be a back door
    for a wrong assumption.** The reference is worth ``a + b = 1`` pseudo-fragment, so on an object with
    its own composition evidence a maximally-wrong location is a tie-breaker at low depth and inaudible
    at high depth. Measured on the single-strand object at ``kappa = 0.9``, swinging the location from
    0.02 to 0.98 moves ``f_g`` by::

        6 frags 0.7700 · 24 0.0565 · 120 0.0090 · 600 0.0009 · 15,000 0.0000

    ⛔ This is what bounds the term's overlap with the fitted arms, the density factor and the messages:
    wherever any of them speaks, it is negligible by construction."""
    swing = lambda scale: float(  # noqa: E731
        abs(
            _solve(CompositionPriors(location=np.array([0.98, 0.98])), scale=scale)[0]
            - _solve(CompositionPriors(location=np.array([0.02, 0.02])), scale=scale)[0]
        )
    )
    sparse, mid, rich = swing(1.0), swing(20.0), swing(500.0)
    assert sparse > 0.5, sparse
    assert mid < 0.1 * sparse, (sparse, mid)
    assert rich < 0.01 * sparse, (sparse, rich)


def test_where_there_is_no_composition_evidence_the_location_is_the_answer():
    """⛔⛔ **THE OTHER HALF, AND IT IS THE DESIGN'S WHOLE RISK SURFACE — pinned so it cannot be
    forgotten.** An object with no composition evidence of its own has posterior = prior at ANY depth,
    so the location does not merely nudge it, it DETERMINES it. Measured: an AMBIG object, and any object
    once ``kappa = ½`` kills the strand channel, keeps a swing of ~0.95 at 6 fragments and at 3,000.

    ⭐ That is correct Bayes, not a defect — but it means the location's accuracy on evidence-free
    objects is load-bearing in a way it is not anywhere else, and that those objects are exactly the
    population the gDNA landscape is fitted to serve. ⛔ So the reference and the landscape do NOT have
    a bounded overlap there; they are the only two voices, and only a solve arm can price which should
    carry it.
    """
    lo = CompositionPriors(location=np.array([0.02, 0.02]))
    hi = CompositionPriors(location=np.array([0.98, 0.98]))
    for scale in (1.0, 500.0):
        # the AMBIG object: both strands admissible, so strand cannot split gDNA from RNA at any kappa
        assert _solve(hi, scale=scale)[1] - _solve(lo, scale=scale)[1] > 0.5, scale


def test_the_location_is_per_slot_not_shared():
    """⭐ Two objects, two different locations, and each must feel its own — the whole point of a
    per-object reference. A shared scalar would make these two solves equal."""
    a = _solve(CompositionPriors(location=np.array([0.02, 0.98])))
    b = _solve(CompositionPriors(location=np.array([0.98, 0.02])))
    assert not np.allclose(a, b), (a, b)
    assert a[0] < b[0] and a[1] > b[1], (a, b)
