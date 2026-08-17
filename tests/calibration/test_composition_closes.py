"""⭐⭐⭐ ψ's COMPOSITION IS A POINT ON THE SIMPLEX, AND THAT IS NOW STRUCTURAL.

`~rigel.calibration.simplex_logodds._compose` builds ``(f_pos, f_neg)`` as the IMAGE of ψ's two
parameters rather than reading three coordinates off the posterior independently::

    lambda  — the gDNA-vs-RNA LEVEL, read out as f_g (the posterior median over the lambda grid)
    theta   — the RNA-internal TILT, a SHARE with no absolute scale, read out as w_pos

    f_pos = (1 - f_g)*w_pos        f_neg = (1 - f_g)*(1 - w_pos)

so ``f_g + f_pos + f_neg = 1`` identically. **Closure cannot fail**, which is why this file asserts it
rather than measuring it.

⛔⛔ **WHAT IT REPLACED, AND WHY THAT WAS NOT A DECISION.** ``f_g`` was the posterior MEDIAN while
``f_pos``/``f_neg`` were posterior MEANS of the grid quantity ``1 - f_g``. Mixing a quantile with
expectations does not land on the simplex::

    SUM = median(f_g) + (1 - mean(f_g)) = 1 + median(f_g) - mean(f_g)

— the closure error was EXACTLY the posterior's skew — and only **74.7 % of REGIONs and 77.2 % of
BOUNDARYs** closed on a real condition, p5 0.869/0.850. ⭐ The asymmetry was never chosen: the median
for ``f_g`` was argued for in `_posterior_median_fg`; the RNA fractions fell out as expectations of a
grid array nobody revisited.

⛔ **AND THE OBVIOUS ALTERNATIVE IS MEASURED CATASTROPHIC.** Taking MEANS everywhere also closes, by
linearity of expectation, and scores **1.352 / 1.573 / 3.756** on the three in-scope strata and
**1.801** on the zero control (`vertex_ceiling.py --arm psi_mean`, all 16 conditions): the median is
closer to truth at both simplex vertices, where 49-83 % of in-scope error lives.

⚠ **Nor is this the renormalisation `ROADMAP.md` refuses.** Nothing is rescaled to hide a residual —
``f_g`` and the RNA total are exact complements *by parametrisation*, and the tilt is estimated as a
share because a share is what it is.

⭐ Each gate carries its own perturbation (`TRAPS: perturb-every-gate`).
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from rigel.calibration.simplex_logodds import _compose, _solve_regions_logodds_all

_L = 10.0


def _solve(u_pos, u_neg, *, kappa, allow_pos, allow_neg, n_grid=60, n_grid_ss=256, n_tilt=None):
    d = _solve_regions_logodds_all(
        np.asarray(u_pos, np.float64),
        np.asarray(u_neg, np.float64),
        np.asarray(allow_pos, bool),
        np.asarray(allow_neg, bool),
        np.asarray(u_pos, np.float64) + np.asarray(u_neg, np.float64),
        np.zeros_like(np.asarray(u_pos, np.float64)),
        kappa=kappa,
        od_g=0.0,
        od_r=0.0,
        n_grid=n_grid,
        n_grid_ss=n_grid_ss,
        n_tilt=n_tilt,
        L=_L,
    )
    return (
        np.asarray(d.gdna_frac, np.float64),
        np.asarray(d.rna_pos_frac, np.float64),
        np.asarray(d.rna_neg_frac, np.float64),
    )


_DEPTHS = (1.0, 3.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0)
_KAPPAS = (0.5, 0.75, 0.9, 0.99)
_SPLITS = (0.0, 0.01, 0.2, 0.5)


@pytest.mark.parametrize("kappa", _KAPPAS)
@pytest.mark.parametrize(
    "allow_pos,allow_neg", ((True, False), (False, True), (True, True)), ids=("ss+", "ss-", "ambig")
)
def test_the_composition_closes_exactly_on_every_solved_slot(kappa, allow_pos, allow_neg):
    """⭐⭐⭐ **THE PROPERTY THE PARAMETRISATION BUYS.** Not "closes to a tolerance" — closes to float64
    round-off, on both ψ paths, at every depth and strand split, because the read-out is a map onto the
    simplex rather than three independent summaries.

    ⛔ The perturbation is below: a composition assembled the OLD way, from three independent read-outs
    on the same posterior, does NOT close — so this gate is not passing because the quantity is trivial.
    """
    u_pos = np.array([n * (1.0 - f) for n, f in itertools.product(_DEPTHS, _SPLITS)])
    u_neg = np.array([n * f for n, f in itertools.product(_DEPTHS, _SPLITS)])
    ap = np.full(u_pos.shape, allow_pos)
    an = np.full(u_pos.shape, allow_neg)
    f_g, f_pos, f_neg = _solve(u_pos, u_neg, kappa=kappa, allow_pos=ap, allow_neg=an)
    total = f_g + f_pos + f_neg
    assert np.max(np.abs(total - 1.0)) < 1e-12, (kappa, np.max(np.abs(total - 1.0)))
    # ⛔ and the dead strand really is dead, so "closes" is not being bought by leaking RNA onto it
    if not allow_pos:
        assert np.all(f_pos == 0.0)
    if not allow_neg:
        assert np.all(f_neg == 0.0)
    # ⛔ non-vacuous: BOTH components are actually being placed, so closure is not holding because one
    #   side is identically zero. ⚠ Deliberately not a span check on ``f_g`` — at κ = ½ the strand
    #   channel is dead by derivation and every slot sits near the reference, which is correct.
    assert f_g.max() > 0.0 and (f_pos + f_neg).max() > 0.0, (f_g.max(), (f_pos + f_neg).max())


def test_the_old_read_out_did_not_close_and_the_gap_was_the_skew():
    """⛔ **THE PERTURBATION FOR THE GATE ABOVE, AND THE RECORD OF WHAT WAS WRONG.** Reconstruct the
    retired read-out — ``f_g`` from the median, the RNA total from the posterior MEAN — and show it
    misses the simplex by exactly ``median − mean``.

    ⭐ The posterior mean is recoverable from the shipped output without re-solving: the RNA total is
    ``1 − f_g`` now, so feeding the OLD RNA total means replacing it with ``1 − E[f_g]``."""
    from rigel.calibration.simplex_logodds import (
        _local_loglik_logodds,
        _logodds_grid,
        _lse,
        _posterior_median_fg,
    )

    lam, fg = _logodds_grid(60, _L)
    u_pos, u_neg = np.array([27.0, 9.0, 3.0]), np.array([3.0, 1.0, 0.0])
    psi = _local_loglik_logodds(
        u_pos,
        u_neg,
        np.array([True, True, True]),
        np.array([False, False, False]),
        0.99,
        0.0,
        0.0,
        lam,
        fg,
        np.full(3, 0.5),
        np.full(3, 0.5),
        np.zeros(3),
    )
    post = np.exp(psi - _lse(psi, axis=1, keepdims=True))
    median = _posterior_median_fg(post, lam, fg)
    mean = np.sum(post * fg[None, :], axis=1)
    old_sum = median + (1.0 - mean)  # the retired three-read-out composition
    assert np.max(np.abs(old_sum - (1.0 + median - mean))) < 1e-12
    # ⛔ it genuinely misses the simplex on this fixture — otherwise the gate above is vacuous
    assert np.max(np.abs(old_sum - 1.0)) > 0.02, old_sum
    # ⭐ and the shipped read-out on the same counts closes
    f_g, f_pos, f_neg = _solve(
        u_pos, u_neg, kappa=0.99, allow_pos=np.full(3, True), allow_neg=np.full(3, False)
    )
    assert np.max(np.abs(f_g + f_pos + f_neg - 1.0)) < 1e-12


def test_compose_enforces_admissibility_so_rna_cannot_leak_onto_a_forbidden_strand():
    """⛔⛔ **THE ONE WAY A MAP ONTO THE SIMPLEX CAN STILL BE WRONG, AND THIS GATE FOUND IT.** If the
    tilt share were applied blind, a slot with only the + strand admissible and ``w_pos = ½`` would place
    half its RNA on the forbidden − strand, where it is zeroed — and the RNA would simply VANISH, giving
    ``SUM = f_g + (1−f_g)/2``. `_compose` therefore restricts the share to the admissible strands itself
    rather than trusting the caller, so the single-strand path has no tilt to supply at all.

    ⭐ Closure must hold for EVERY ``w_pos``, including ones that are wrong for the slot."""
    f_g = np.array([0.0, 0.05, 0.5, 0.95, 1.0])
    T, F = np.full(5, True), np.full(5, False)
    for w in (0.0, 0.25, 0.5, 0.75, 1.0):
        for ap, an in ((T, T), (T, F), (F, T)):
            p, n = _compose(f_g, np.full(5, w), ap, an)
            assert np.max(np.abs(f_g + p + n - 1.0)) < 1e-15, (w, ap[0], an[0])
            assert np.all(p >= 0.0) and np.all(n >= 0.0)
            assert np.all(p[~ap] == 0.0) and np.all(n[~an] == 0.0)
    # ⭐ a single-strand slot takes the WHOLE RNA total on its admissible strand, whatever w_pos says
    for w in (0.0, 0.5, 1.0):
        p, n = _compose(f_g, np.full(5, w), T, F)
        assert np.allclose(p, 1.0 - f_g) and np.all(n == 0.0), w
        p, n = _compose(f_g, np.full(5, w), F, T)
        assert np.allclose(n, 1.0 - f_g) and np.all(p == 0.0), w
    # ⛔ neither strand admissible ⇒ no RNA to place; the composition is f_g alone, which is a
    #   statement about opportunity and not a closure failure
    p, n = _compose(f_g, np.full(5, 0.5), F, F)
    assert np.all(p == 0.0) and np.all(n == 0.0)


def test_a_share_outside_the_unit_interval_cannot_produce_a_negative_fraction():
    """⛔⛔ **THE ONE WAY A CLOSING COMPOSITION CAN STILL BE NONSENSE, FOUND BY HUNTING FOR IT.** An
    unclamped share yields a NEGATIVE fraction that still sums to 1 — it would pass every closure gate in
    this file. Measured before the clamp: ``w_pos = −1`` gave ``f_pos = −0.8, f_neg = 1.6``, SUM = 1.

    ⚠ Not reachable from the shipped AMBIG caller — ``w_pos = m_pos/(m_pos+m_neg)`` is a ratio of two
    non-negative expectations — which is precisely why the constraint is asserted here rather than assumed
    at the call site."""
    f_g = np.array([0.0, 0.3, 0.9])
    T = np.full(3, True)
    for w in (-1.0, -1e-9, 1.0 + 1e-9, 2.0, 1e9):
        p, n = _compose(f_g, np.full(3, w), T, T)
        assert np.all(p >= 0.0) and np.all(n >= 0.0), (w, p, n)
        assert np.max(np.abs(f_g + p + n - 1.0)) < 1e-15, w
    # ⛔ the perturbation: an IN-range share is passed through untouched, so the clamp is not flattening
    #   the tilt into a constant
    p_lo, _ = _compose(f_g, np.full(3, 0.25), T, T)
    p_hi, _ = _compose(f_g, np.full(3, 0.75), T, T)
    assert np.all(p_hi > p_lo)


def test_the_ambig_tilt_share_stays_a_share():
    """⭐ The AMBIG caller derives ``w_pos`` from two cube expectations of non-negative grid quantities, so
    it is in [0,1] by construction. This asserts that on the SHIPPED grid rather than trusting the
    algebra — an f32 cube reduction is where such a guarantee would quietly fail."""
    for kappa in (0.5, 0.9, 0.99):
        for n in (1.0, 10.0, 1e3, 1e6):
            for frac in (0.0, 0.1, 0.5, 0.9, 1.0):
                f_g, f_pos, f_neg = _solve(
                    np.array([n * (1 - frac)]),
                    np.array([n * frac]),
                    kappa=kappa,
                    allow_pos=np.array([True]),
                    allow_neg=np.array([True]),
                )
                rna = f_pos + f_neg
                assert np.all(rna >= -1e-15) and np.all(f_pos >= -1e-15) and np.all(f_neg >= -1e-15)
                assert np.allclose(rna, 1.0 - f_g, atol=1e-15), (kappa, n, frac)


def test_a_zero_count_slot_publishes_no_data_rather_than_a_composition():
    """⚠ **THE ONE POPULATION THAT DOES NOT CLOSE, AND IT IS DELIBERATE.** A slot with no fragments, or
    with neither strand admissible, is not dispatched to either ψ solve; it publishes ``(0, 0, 0)``,
    which is "no data" and not a composition claim. `region_init` replaces it with the signature-binary
    init, which IS a simplex point.

    ⛔ Recorded so the closure gate above is read with its true scope: closure is asserted over SOLVED
    slots, and this is what the complement looks like."""
    f_g, f_pos, f_neg = _solve(
        np.array([0.0, 30.0]),
        np.array([0.0, 1.0]),
        kappa=0.99,
        allow_pos=np.array([True, False]),
        allow_neg=np.array([False, False]),
    )
    assert np.all(f_g == 0.0) and np.all(f_pos == 0.0) and np.all(f_neg == 0.0)
