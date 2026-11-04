"""What a certified-RNA count may and may not claim about a slot's composition.

A spliced fragment cannot be gDNA, so ``boundary_spliced`` (and ``sj_count``) are certified RNA, the one
observation in the tool that needs no deconvolution. Turning it into per-slot composition evidence would add a
coefficient ``S`` to ψ's RNA arm, ``(½ + S)·log(1 − f_g)`` in place of the bare Jeffreys reference
``½·log(1 − f_g)``. Half of that derivation holds and half does not; these tests pin where the line falls, as
brute-force arithmetic against scipy's own pmfs. The model is one population of contiguous crossings at one
boundary, split by whether a splice was visible::

    E[S | f_g, q]  =  (q / (1 − q)) · (1 − f_g) · M

    log P(S)  =  S·log(1 − f_g)   −   (q/(1−q))·(1 − f_g)·M   +   const
                 \\___ retained ___/    \\_______ dropped _______/

===  ==========================================================================================
C1   sound: the retained coefficient is the raw count, and no opportunity ratio can reach it
C2   not sound: with ``q`` unknown the profile likelihood in ``f_g`` is exactly flat on [0,1)
C2b  and ``S = 0`` is flat on the closed interval, so a zero count is not vertex evidence either
C3   the raw-count term is one-sided and unbounded: at ``S = 1000`` it answers ~0 whatever the truth
C4   sound: reference + term is exactly Beta(½, ½+S), so ``S = 0`` recovers the reference ψ identically
===  ==========================================================================================

The realised ``q`` is not small enough for the dropped term to be a correction, so keeping only the retained term
would be an assumption about ``q`` rather than evidence, and it is worst where a slot is gDNA-rich. ψ is blind to
this channel.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import expit, log_expit
from scipy.stats import beta as _Beta, poisson as _Poisson

from rigel.calibration.simplex_logodds import _logodds_grid
from rigel.native import transfer_rows as R

#: Every test scores against ``scipy.stats``' own Poisson / Beta rather than a log-pmf written here, since the claim
#: under test is an algebraic decomposition of exactly this pmf.
_C_GRID = (1e-4, 0.03, 1.0, 40.0)
_M_GRID = (5.0, 900.0, 20000.0)


def _residual(S, c, M, lam):
    """``log P(S | c·(1−σ(λ))·M)`` with the retained term ``S·log(1−σ(λ))`` removed."""
    return _Poisson.logpmf(S, c * (1.0 - expit(lam)) * M) - S * log_expit(-lam)


# ── C1, the sound half: no opportunity ratio can reach the retained coefficient ─────────────


@pytest.mark.parametrize("S", [1.0, 7.0, 500.0])
def test_C1_the_retained_coefficient_is_the_RAW_count_whatever_the_two_banks_DIVISORS_are(S):
    """``E[S] = c·(1−f_g)·M`` for some positive ``c``, and every opportunity ratio between the certified bank's divisor
    and the unspliced bank's lives inside ``c``. A multiplicative factor on the mean is an additive constant in log
    space, so it cannot touch the coefficient of ``log(1−f_g)``; neither ``eff_rna`` nor the sj divisor appears in the
    retained term, so there is nothing for two divisors to disagree about.

    Strip ``S·log(1−f_g)`` off scipy's exact Poisson log-pmf, and what remains must be ``−c·M·(1−f_g)`` plus a λ-free
    constant for every ``c``, ``M`` and ``S``: adding the linear term back must leave something flat in λ, to machine
    precision.
    """
    lam = np.linspace(-6.0, 6.0, 4001)
    for c in _C_GRID:
        for M in _M_GRID:
            flat = _residual(S, c, M, lam) + c * M * (1.0 - expit(lam))
            assert np.allclose(flat, flat[0], rtol=0.0, atol=1e-9), (
                c,
                M,
                float(np.abs(flat - flat[0]).max()),
            )


def test_C1_perturbation_an_OPPORTUNITY_SCALED_coefficient_FIRES():
    """``c·S``, an opportunity-scaled count, leaves a λ-dependent remainder, so the test above is not a truism about
    subtraction."""
    lam = np.linspace(-6.0, 6.0, 4001)
    S, c, M = 7.0, 0.03, 900.0
    bad = _Poisson.logpmf(S, c * (1.0 - expit(lam)) * M) - (c * S) * log_expit(-lam)
    flat = bad + c * M * (1.0 - expit(lam))
    assert not np.allclose(flat, flat[0], rtol=0.0, atol=1e-9)


# ── C2, the unsound half: q is not a nuisance that may be dropped ────────────────────────────────


@pytest.mark.parametrize("S", [1.0, 3.0, 40.0, 1000.0])
@pytest.mark.parametrize("M", [4.0, 350.0, 20000.0])
def test_C2_with_q_UNKNOWN_the_profile_likelihood_in_f_g_is_EXACTLY_FLAT(S, M):
    """``q``, the chance a crossing RNA fragment shows a visible splice, is a free parameter on (0,1), so ``q/(1−q)``
    spans (0, ∞). For any ``f_g < 1`` some ``q`` makes ``E[S]`` exactly ``S``, so the profile likelihood ``sup_q log P(S
    | f_g, q)`` is the same at every ``f_g`` in [0,1). The certified count carries one bit about the unspliced split
    ("f_g ≠ 1"); anything more the ``S·log(1−f_g)`` coefficient claims is an assumption about ``q``.

    Brute force over a dense grid of ``ratio = q/(1−q)``, log-spaced over 26 decades so it reaches every mean the model
    can produce, scored with scipy's own Poisson rather than the analytic argument, so an error in the argument itself
    would show.
    """
    fg = np.array([0.0, 0.05, 0.4, 0.5, 0.9, 0.99, 0.99999])
    exact = float(_Poisson.logpmf(S, S))  # the unconstrained Poisson maximum

    def spread(n):
        ratio = np.exp(np.linspace(-30.0, 30.0, n))  # q/(1−q) over q ∈ (0,1), 26 decades
        p = np.array([_Poisson.logpmf(S, ratio * (1.0 - f) * M).max() for f in fg])
        return float(p.max() - p.min()), float(np.abs(p - exact).max())

    # Stated as convergence, not a loosened tolerance: the claim is exact for continuous ``q``, so on a grid the only
    # residual is the quadratic gap at the maximum, which must fall ~4x per halved spacing. A constant residual would
    # mean the profile is tilted.
    coarse, coarse_err = spread(30001)
    fine, fine_err = spread(120001)
    assert fine < 1e-4 and fine_err < 1e-4, (fine, fine_err)
    # 4x the grid density; the ideal quadratic gain is 16x, lower in practice because the maximum does not sit on a grid
    # point. Checked at 4x, far above the 1.0 a tilted profile would give.
    assert fine < 0.25 * coarse, (coarse, fine)
    assert fine_err < 0.25 * coarse_err, (coarse_err, fine_err)
    # the one bit it does carry: f_g = 1 exactly is excluded, since a mean of 0 cannot produce S ≥ 1.
    assert float(_Poisson.logpmf(S, 0.0)) == -np.inf


def test_C2_perturbation_FIXING_q_restores_an_interior_maximum_and_FIRES():
    """The flatness comes from ``q`` being unknown, not from the observation being worthless: with ``q`` pinned the same
    likelihood has a sharp interior maximum at the true split, so the channel is blocked on an opportunity model, not on
    physics."""
    M, q_true, fg_true = 360.0, 0.25, 0.6
    S = (q_true / (1.0 - q_true)) * (1.0 - fg_true) * M  # a consistent world; 48 exactly
    assert S == 48.0
    fg = np.linspace(0.0, 0.999, 4001)
    ll = _Poisson.logpmf(S, (q_true / (1.0 - q_true)) * (1.0 - fg) * M)
    assert float(ll.max() - ll.min()) > 1.0  # no longer flat: the check fires
    assert abs(float(fg[int(np.argmax(ll))]) - fg_true) < 1e-2


@pytest.mark.parametrize("M", [4.0, 350.0, 20000.0])
def test_C2b_a_ZERO_certified_count_is_flat_on_the_CLOSED_interval_so_it_is_NOT_vertex_evidence(M):
    """A zero certified count is not evidence for the gDNA vertex. With ``q`` free, ``S = 0`` is perfectly explained at
    every ``f_g``, ``f_g = 1`` included (take ``q → 0``). The channel's information lives entirely in ``S > 0`` and
    points away from the ``f_g → 1`` vertex, never toward it.
    """
    ratio = np.exp(np.linspace(-30.0, 30.0, 240001))
    profile = np.array(
        [_Poisson.logpmf(0.0, ratio * (1.0 - f) * M).max() for f in (0.0, 0.5, 0.9, 1.0)]
    )
    # ``sup_q log P(0 | q) = 0`` at every f_g, the vertex included, reached in the q → 0 limit, so the grid's smallest
    # ratio is the only gap and it is 1e-13·M.
    assert np.allclose(profile, 0.0, rtol=0.0, atol=1e-8), profile


# ── C3: the raw-count term's damage, as a number ────────────────────────────────────────────────────


@pytest.mark.parametrize("S,ceiling", [(10.0, 0.07), (100.0, 8e-3), (1000.0, 8e-4)])
def test_C3_the_raw_count_term_is_ONE_SIDED_and_UNBOUNDED_in_S(S, ceiling):
    """The term can only lower ``f_g``, never raise it, and does so without limit: the answer is ``Beta(½, ½+S)``'s
    median whatever gDNA the slot holds, and ``M`` never enters. At a gDNA-rich boundary with a large certified count it
    answers near zero against a truth near one.
    """
    med = float(_Beta.ppf(0.5, 0.5, 0.5 + S))
    assert med < ceiling, (S, med)
    # strictly monotone down in S, the one-sidedness itself
    lo = float(_Beta.ppf(0.5, 0.5, 0.5 + S))
    hi = float(_Beta.ppf(0.5, 0.5, 0.5 + 2.0 * S))
    assert hi < lo


# ── C4: the sound half, reusable by any build of this channel ─────────────────────────────


@pytest.mark.parametrize("S", [0.0, 1.0, 10.0, 100.0])
def test_C4_psi_reference_plus_the_term_is_EXACTLY_Beta_half_half_plus_S(S):
    """ψ's two arms on the λ axis, ``½·log f_g + (½+S)·log(1−f_g)``, are the log-density of ``f_g ~ Beta(½, ½+S)``: the
    Jacobian ``df_g/dλ = σ(λ)σ(−λ)`` supplies the missing ``−1`` on each exponent exactly. So the certified count is the
    Bayes update of the reference, and ``S = 0`` returns the reference Beta(½,½) identically, with no new normalisation
    and no ``S = 0`` branch. The coefficient is unlicensed (C2), not the algebra.
    """
    lam, fg = _logodds_grid(1 << 16, 20.0)
    lp = R.JEFFREYS_REF * log_expit(lam) + (R.JEFFREYS_REF + S) * log_expit(-lam)
    w = np.exp(lp - lp.max())
    w /= w.sum()
    grid_median = float(fg[int(np.searchsorted(np.cumsum(w), 0.5))])
    assert grid_median == pytest.approx(float(_Beta.ppf(0.5, 0.5, 0.5 + S)), abs=2e-4)
