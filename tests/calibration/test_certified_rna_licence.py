"""⛔⛔ WHAT A CERTIFIED-RNA COUNT MAY AND MAY NOT CLAIM — the licence, as brute-force arithmetic.

A spliced fragment cannot be gDNA, so ``boundary_spliced`` (and ``sj_count``) are **certified RNA**: the one
observation in the tool that needs no deconvolution. The standing proposal is to turn that into per-slot
composition evidence by adding a coefficient ``S`` to ψ's RNA arm — ``(½ + S)·log(1 − f_g)`` in place of
the bare Jeffreys reference ``½·log(1 − f_g)``.

**Half of that derivation is right and half of it is not**, and these gates pin exactly where the line is.
The model is one population of contiguous crossings at one line, split by whether a splice was VISIBLE::

    E[S | f_g, q]  =  (q / (1 − q)) · (1 − f_g) · M

    log P(S)  =  S·log(1 − f_g)   −   (q/(1−q))·(1 − f_g)·M   +   const
                 \\___ retained ___/    \\_______ dropped _______/

===  ==========================================================================================
TRAPS: a-purity-filter-is-a-length-filter   ⭐ SOUND — the retained coefficient is the RAW COUNT, and no opportunity ratio can reach it
TRAPS: pure-and-length-censored   ⛔ NOT SOUND — with ``q`` unknown the profile likelihood in ``f_g`` is EXACTLY FLAT on [0,1)
C2b  ⛔ and ``S = 0`` is flat on the CLOSED interval, so a zero count is not vertex evidence either
TRAPS: divide-by-a-probability   ⛔ the raw-count term is one-sided and UNBOUNDED — at ``S = 1000`` it answers ~0 whatever the truth
TRAPS: opposite-tilts-must-not-pool   ⭐ SOUND — reference + term is exactly Beta(½, ½+S), so ``S = 0`` recovers today's ψ identically
TRAPS: fractional-mass-is-the-problem   ⛔ ``density_factor_precision`` must NOT price this factor — on a monotone factor it reads the WINDOW
===  ==========================================================================================

⭐ **The measured consequence, which is why this file exists rather than a mechanism**
(`scripts/design/certified_q_census.py`, full 36-condition ladder against origin-split truth): the realised
``q`` has a mass-weighted median of **0.19–0.71** with 60–98 % of the mass below 0.9, so the dropped term
is the same size as the retained one. Scored against truth the raw-count term is **worse than the
uninformative reference on 12 of 36 conditions**, worst **+0.4578** mwae at ``g90 ss0.50 capture_on``. Its
benefit tracks the ANSWER (−0.49 where the truth is all-RNA, +0.45 where it is all-gDNA), not the
evidence — TRAPS: honesty-metrics-reward-ignorance's shape.

⚠ ψ is still blind to the channel and ``test_vertex_reference.py``'s certified-RNA-blindness gate is therefore still GREEN and still
correct. It is the observable a fix must move; nothing here moves it.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import expit, log_expit
from scipy.stats import beta as _Beta, poisson as _Poisson

from rigel.calibration.density_deconv import density_factor_precision
from rigel.calibration.simplex_logodds import _JEFFREYS_REF, _logodds_grid

#: ⛔ Every gate below scores against ``scipy.stats``' OWN Poisson / Beta rather than a log-pmf written
#: here. A gate that re-derives the thing it is checking cannot detect drift in it (TRAPS: a-test-that-redefines), and
#: the whole claim under test is an algebraic decomposition of exactly this pmf.
_C_GRID = (1e-4, 0.03, 1.0, 40.0)
_M_GRID = (5.0, 900.0, 20000.0)


def _residual(S, c, M, lam):
    """``log P(S | c·(1−σ(λ))·M)`` with the RETAINED term ``S·log(1−σ(λ))`` removed."""
    return _Poisson.logpmf(S, c * (1.0 - expit(lam)) * M) - S * log_expit(-lam)


# ── TRAPS: a-purity-filter-is-a-length-filter — the half that IS sound: no opportunity ratio can reach the retained coefficient ─────────────


@pytest.mark.parametrize("S", [1.0, 7.0, 500.0])
def test_C1_the_retained_coefficient_is_the_RAW_count_whatever_the_two_banks_DIVISORS_are(S):
    """⭐ ``E[S] = c·(1−f_g)·M`` for SOME positive ``c``, and every opportunity ratio between the certified
    bank's divisor and the unspliced bank's lives inside ``c``. A multiplicative factor on the MEAN is an
    ADDITIVE CONSTANT in log space, so it cannot touch the coefficient of ``log(1−f_g)``.

    ⛔ **This is why TRAPS: two-divisors-opposite-sign's trap is structurally absent here, not merely avoided.** TRAPS: two-divisors-opposite-sign is two
    divisors built from one pmf responding to it with OPPOSITE SIGN; here neither ``eff_rna`` nor
    ``eff_junction`` appears in the retained term at all, so there is nothing for them to disagree about.
    The junction half's missing "frame statement" is therefore answered: **it does not need one.**

    The gate: strip ``S·log(1−f_g)`` off scipy's exact Poisson log-pmf, and what remains must be
    ``−c·M·(1−f_g)`` plus a λ-free constant — for every ``c``, every ``M``, and every ``S``. Adding the
    linear term back must leave something FLAT in λ, to machine precision.
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
    """⚠ The falsification for TRAPS: a-purity-filter-is-a-length-filter. ``c·S`` — an opportunity-scaled count, which is the other candidate
    §2f offered — leaves a λ-dependent remainder, so the gate above is not a truism about subtraction."""
    lam = np.linspace(-6.0, 6.0, 4001)
    S, c, M = 7.0, 0.03, 900.0
    bad = _Poisson.logpmf(S, c * (1.0 - expit(lam)) * M) - (c * S) * log_expit(-lam)
    flat = bad + c * M * (1.0 - expit(lam))
    assert not np.allclose(flat, flat[0], rtol=0.0, atol=1e-9)


# ── TRAPS: pure-and-length-censored — the half that is NOT sound: q is not a nuisance you may drop ────────────────────────────────


@pytest.mark.parametrize("S", [1.0, 3.0, 40.0, 1000.0])
@pytest.mark.parametrize("M", [4.0, 350.0, 20000.0])
def test_C2_with_q_UNKNOWN_the_profile_likelihood_in_f_g_is_EXACTLY_FLAT(S, M):
    """⛔⛔ **THE GATE THAT FORBIDS THE RAW-COUNT TERM.**

    ``q`` — the chance a crossing RNA fragment shows a *visible* splice — is not a small correction; it
    is a free parameter spanning (0,1), and ``q/(1−q)`` therefore spans (0, ∞). So for **any** ``f_g < 1``
    there is a ``q`` making ``E[S]`` exactly ``S``, and the profile likelihood ``sup_q log P(S | f_g, q)``
    is the SAME number at every ``f_g`` in [0,1).

    ⭐ **The certified count therefore carries exactly ONE BIT about the unspliced split: "f_g ≠ 1".**
    Everything the ``S·log(1−f_g)`` coefficient claims beyond that bit is an assumption about ``q``,
    not evidence — which is exactly why it measured +0.4578 mwae at ``g90 capture_on``.

    ⚠ Brute force over a dense grid of the natural nuisance parameter ``ratio = q/(1−q)`` — log-spaced
    over 26 decades, so it reaches every mean the model can produce — and scored with scipy's own
    Poisson, not the analytic argument, so the gate would catch an error in the argument itself.
    """
    fg = np.array([0.0, 0.05, 0.4, 0.5, 0.9, 0.99, 0.99999])
    exact = float(_Poisson.logpmf(S, S))  # the unconstrained Poisson maximum

    def spread(n):
        ratio = np.exp(np.linspace(-30.0, 30.0, n))  # q/(1−q) over q ∈ (0,1), 26 decades
        p = np.array([_Poisson.logpmf(S, ratio * (1.0 - f) * M).max() for f in fg])
        return float(p.max() - p.min()), float(np.abs(p - exact).max())

    # ⭐ Stated as CONVERGENCE, not as a loosened tolerance: the claim is exact for continuous ``q``, so
    #   on a grid the only residual is the quadratic gap at the maximum — which must fall ~4x per halved
    #   spacing. A constant residual would mean the profile really is tilted (TRAPS: byte-identity-gate).
    coarse, coarse_err = spread(30001)
    fine, fine_err = spread(120001)
    assert fine < 1e-4 and fine_err < 1e-4, (fine, fine_err)
    # ⚠ 4x the grid density; the ideal quadratic gain is 16x and the measured one is 9.6-9.9x (the
    #   maximum does not sit on a grid point). Gated at 4x — far above the 1.0 a genuinely tilted
    #   profile would give, and far below what is observed.
    assert fine < 0.25 * coarse, (coarse, fine)
    assert fine_err < 0.25 * coarse_err, (coarse_err, fine_err)
    # and the ONE BIT it does carry: f_g = 1 exactly is excluded — a mean of 0 cannot produce S ≥ 1.
    assert float(_Poisson.logpmf(S, 0.0)) == -np.inf


def test_C2_perturbation_FIXING_q_restores_an_interior_maximum_and_FIRES():
    """⚠ The falsification for TRAPS: pure-and-length-censored, and the constructive half of it: the flatness is a statement about
    ``q`` being unknown, NOT about the observation being worthless. Pin ``q`` and the same likelihood
    acquires a sharp interior maximum at the true split — so the channel is blocked on an opportunity
    model, not on physics."""
    M, q_true, fg_true = 360.0, 0.25, 0.6
    S = (q_true / (1.0 - q_true)) * (1.0 - fg_true) * M  # a consistent world; 48 exactly
    assert S == 48.0
    fg = np.linspace(0.0, 0.999, 4001)
    ll = _Poisson.logpmf(S, (q_true / (1.0 - q_true)) * (1.0 - fg) * M)
    assert float(ll.max() - ll.min()) > 1.0  # NOT flat any more — the gate fires
    assert abs(float(fg[int(np.argmax(ll))]) - fg_true) < 1e-2


@pytest.mark.parametrize("M", [4.0, 350.0, 20000.0])
def test_C2b_a_ZERO_certified_count_is_flat_on_the_CLOSED_interval_so_it_is_NOT_vertex_evidence(M):
    """⛔ **AND THIS INVERTS THE HOPE THAT OPENED THE INVESTIGATION.** ``test_vertex_reference``'s certified-RNA-blindness gate argues
    that at a silent gene the certified channel is "exactly zero, which is the strongest possible evidence
    for the vertex". It is not: with ``q`` free, ``S = 0`` is perfectly explained at **every** ``f_g``
    *including* ``f_g = 1`` — take ``q → 0``. The information in this channel lives entirely in ``S > 0``
    and it points AWAY from the ``f_g → 1`` vertex, never toward it.

    ⭐ So the vertex population (a measured negative) gets nothing from here either.
    """
    ratio = np.exp(np.linspace(-30.0, 30.0, 240001))
    profile = np.array(
        [_Poisson.logpmf(0.0, ratio * (1.0 - f) * M).max() for f in (0.0, 0.5, 0.9, 1.0)]
    )
    # ``sup_q log P(0 | q) = 0`` at every f_g, the vertex INCLUDED — reached in the q → 0 limit, so the
    # grid's own smallest ratio is the only gap and it is 1e-13·M.
    assert np.allclose(profile, 0.0, rtol=0.0, atol=1e-8), profile


# ── TRAPS: divide-by-a-probability — the raw-count term's damage, as a number ────────────────────────────────────────────────────


@pytest.mark.parametrize("S,ceiling", [(10.0, 0.07), (100.0, 8e-3), (1000.0, 8e-4)])
def test_C3_the_raw_count_term_is_ONE_SIDED_and_UNBOUNDED_in_S(S, ceiling):
    """⛔ The term can only lower ``f_g``, never raise it, and it does so without limit — the answer is
    ``Beta(½, ½+S)``'s median regardless of how much gDNA the slot actually holds. ``M`` never enters.

    ⭐ That is the mechanism behind the panel measurement: 2,565 BOUNDARIES at ``g90 capture_on`` carrying
    **98.3 %** of the certified-bearing mass have a true ``f_g`` of 0.84, and this answers 0.04.
    """
    med = float(_Beta.ppf(0.5, 0.5, 0.5 + S))
    assert med < ceiling, (S, med)
    # strictly monotone DOWN in S — the one-sidedness itself
    lo = float(_Beta.ppf(0.5, 0.5, 0.5 + S))
    hi = float(_Beta.ppf(0.5, 0.5, 0.5 + 2.0 * S))
    assert hi < lo


# ── TRAPS: opposite-tilts-must-not-pool — the half that survives, and is what a future build should reuse ─────────────────────────────


@pytest.mark.parametrize("S", [0.0, 1.0, 10.0, 100.0])
def test_C4_psi_reference_plus_the_term_is_EXACTLY_Beta_half_half_plus_S(S):
    """⭐ ψ's two arms on the λ axis, ``½·log f_g + (½+S)·log(1−f_g)``, are the log-density of
    ``f_g ~ Beta(½, ½+S)`` — the Jacobian ``df_g/dλ = σ(λ)σ(−λ)`` supplies the missing ``−1`` on each
    exponent exactly. So the certified count is the *Bayes update of the reference*, and ``S = 0``
    returns today's Beta(½,½) identically: whatever a future build does with this channel, it needs no
    new normalisation and no ``S = 0`` branch.

    ⛔ It is the COEFFICIENT that is unlicensed (TRAPS: pure-and-length-censored), not the algebra."""
    lam, fg = _logodds_grid(1 << 16, 20.0)
    lp = _JEFFREYS_REF * log_expit(lam) + (_JEFFREYS_REF + S) * log_expit(-lam)
    w = np.exp(lp - lp.max())
    w /= w.sum()
    grid_median = float(fg[int(np.searchsorted(np.cumsum(w), 0.5))])
    assert grid_median == pytest.approx(float(_Beta.ppf(0.5, 0.5, 0.5 + S)), abs=2e-4)


# ── TRAPS: fractional-mass-is-the-problem — the precision the plan reached for reports the WINDOW, not the information ──────────────────


def test_C5_density_factor_precision_reads_the_GRID_WINDOW_on_a_MONOTONE_factor():
    """⛔⛔ The obvious wiring — ``tau_lam += density_factor_precision(cert, lam_grid)``, exactly as
    ``tau_len`` is wired — is out of contract for this factor and the tell is measurable.

    ``density_factor_precision`` reads ``1/Var_λ`` under the NORMALIZED factor. For a peaked factor that
    is the Laplace precision and it is a property of the factor: gated below at exactly 1.0 and 25.0,
    unchanged from ``L = 6`` to ``L = 20``. The certified factor ``S·log σ(−λ)`` is MONOTONE — it has no
    peak and no scale — so its normalized variance is the window's, and the reported "precision" moves by
    orders of magnitude with ``L``.

    ⭐ ``simplex_logodds``' own stated acceptance test is **L-invariance**; a τ that scales with ``L`` is
    disqualified by it. The honest λ-axis information of the term is analytic and window-free —
    ``I = −∂²/∂λ²[S·log σ(−λ)] = S·f_g·(1−f_g)`` — which is what `strand_evidence` returns for its own
    channel; that is the form a future build must use.
    """
    got_cert, got_peak = [], []
    for L in (6.0, 10.0, 20.0):
        lam, _ = _logodds_grid(4096, L)
        got_cert.append(float(density_factor_precision((1e4 * log_expit(-lam))[None, :], lam)[0]))
        got_peak.append(
            float(density_factor_precision((-0.5 * 25.0 * (lam - 1.0) ** 2)[None, :], lam)[0])
        )
    assert np.allclose(got_peak, 25.0, rtol=1e-6), got_peak  # a real factor: L-invariant
    assert max(got_cert) / min(got_cert) > 100.0, got_cert  # the certified one: it IS the window
    # ⭐ and the analytic form it should be replaced by has no grid in it at all: evaluated at a given
    #   composition it is one number, identical whatever window the solver happens to be gridded on.
    i_cert = lambda f: 1e4 * f * (1.0 - f)  # noqa: E731
    assert i_cert(0.5) == pytest.approx(2500.0, rel=0.0, abs=0.0)
    assert i_cert(0.03) == pytest.approx(291.0, rel=1e-12)
