"""⛔⛔⛔ A ZERO COUNT OVER A KNOWN OPPORTUNITY IS A MEASUREMENT — TRAPS: a-zero-count-is-a-measurement.

The message currency is a **log**-density, and ``own_precision``'s count term was ``1/n``: the
LARGE-COUNT LIMIT of the Poisson log-rate variance. At ``n = 0`` that limit diverges, so an object with
no counts emitted nothing — even when it is structurally pure gDNA and therefore composition-CERTAIN.
Measured consequence at ``g00`` (zero gDNA by construction): all **1,298** intergenic regions hold exactly
zero counts over **50.7 Mb** of opportunity and **all 1,298 are silent**, pass-0 invents 34–38 % gDNA,
and 100 % of the error is ``relay_only`` with 0 % ``struct_lock``.

⭐ **THE UNIFYING FACT.** Every count channel here is a Poisson rate observed over an opportunity, and
under the Jeffreys prior ψ is already built on (``_JEFFREYS_REF``) the exact posterior is
``Gamma(a + ½, E)``, whose log has variance ``trigamma(a + ½)`` — proper and finite at ``a = 0``, and
asymptotically ``1/a``. So ``1/n`` and the ``n > 0`` branch are one expression's limit.

⚠ **ONLY the variance moves.** The relay fuses in LINEAR density space, so ``ρ = 0`` was always
expressible; the location stays ``f_c·M/E_c`` because the three components are SHARES of one total and
``Σ_c ρ_c·E_c = M`` must hold exactly — Z5.

===  ===========================================================================================
Z1   ``own_precision``'s count term is ``trigamma(n + ½)`` exactly (vs scipy, not re-derived)
Z2   ⭐ ASYMPTOTIC SAFETY — it agrees with the old ``1/n`` to <0.1 % for ``n ≥ 10``, so the change
     is confined to the low-count population. This is what bounds the blast radius
Z3   ⭐⭐ FINITE AND POSITIVE AT ``n = 0`` — the defect, as a number. Precision ``1/4.9348``
Z4   a STRUCTURALLY CERTAIN object with zero counts now EMITS (``struct_lock`` ⇒ ``Var(log f) = 0``)
Z5   ⛔ the LOCATION must NOT move — ``Σ_c ρ_c·E_c = M`` is exact, and the first version of this fix
     broke it by 3/2. The half belongs to the rate's VARIANCE, not to a share of a total
Z6   monotonicity: more counts ⇒ more precision, at every count including 0 → 1
===  ===========================================================================================
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import polygamma

from rigel.calibration.region_init import (
    count_logvar,
    own_composition_logvar,
    own_precision,
)


# ── Z1 — the exact expression ───────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("n", [0.0, 1.0, 3.0, 10.0, 250.0, 112_333.0])
def test_Z1_the_count_term_is_the_exact_Poisson_log_rate_variance(n):
    """``Var(log ρ) = trigamma(a + ½)`` for ``ρ ~ Gamma(a + ½, E)`` — scored against scipy's own
    ``polygamma``, which is a different implementation from ours (TRAPS: a-test-that-redefines)."""
    assert count_logvar(np.array([n]))[0] == pytest.approx(polygamma(1, n + 0.5), rel=1e-12)


def test_Z1_it_does_not_depend_on_the_exposure():
    """⭐ ``Var(log ρ)`` is exposure-FREE — ``E`` only shifts the location, it cannot sharpen the
    claim. A count term that moved with ``E`` would be double-counting the opportunity."""
    v = count_logvar(np.array([0.0, 7.0]))
    assert np.array_equal(v, count_logvar(np.array([0.0, 7.0])))  # pure function of the count


# ── Z2 — the property that bounds the blast radius ──────────────────────────────────────────────


@pytest.mark.parametrize("n", [10.0, 50.0, 1_000.0, 100_000.0])
def test_Z2_it_agrees_with_the_retired_one_over_n_to_better_than_a_tenth_of_a_percent(n):
    """⭐⭐ THE SAFETY PROPERTY. ``trigamma(n+½) → 1/n``, so every object with a real count keeps its
    old answer to <0.1 % and the change lands ONLY where the old form was broken. Without this the
    fix would be a whole-panel perturbation dressed as a bug fix."""
    assert count_logvar(np.array([n]))[0] == pytest.approx(1.0 / n, rel=1e-3)


def test_Z2_and_it_DIVERGES_from_one_over_n_exactly_where_the_defect_was():
    """⚠ The other half of Z2, or Z2 would pass on a no-op: at small counts the two must genuinely
    disagree, which is the whole reason for the change."""
    assert count_logvar(np.array([1.0]))[0] == pytest.approx(0.9348, abs=1e-3)  # vs 1/1
    assert count_logvar(np.array([2.0]))[0] == pytest.approx(0.4903, abs=1e-3)  # vs 1/2


# ── Z3/Z4 — the defect, as numbers ──────────────────────────────────────────────────────────────


def test_Z3_a_zero_count_carries_FINITE_POSITIVE_precision():
    """⭐⭐ THE FIX, IN ONE ASSERTION. ``trigamma(½) = π²/2``, so a zero-count object's log-density
    claim has sd 2.22 nats — loose, but finite, and infinitely more than the nothing it emitted
    before."""
    v = float(count_logvar(np.array([0.0]))[0])
    assert v == pytest.approx(np.pi**2 / 2.0, rel=1e-12)
    assert np.isfinite(v) and v > 0.0
    assert np.sqrt(v) == pytest.approx(2.221, abs=1e-3)


def test_Z4_a_structurally_CERTAIN_object_with_ZERO_counts_now_EMITS():
    """⭐⭐⭐ THE OBJECT THE DISSECTION FOUND: an intergenic region, ``struct_lock`` ⇒ composition
    certain (``Var(log f_g) = 0``), holding zero counts over a large opportunity. It was
    *certain and silent at the same time*."""
    v_fg, _ = own_composition_logvar(np.array([1.0]), np.array([0.0]), np.array([True]))
    assert v_fg[0] == 0.0  # certain, as before
    p = own_precision(np.array([0.0]), v_fg, np.array([True]))
    assert p[0] == pytest.approx(1.0 / (np.pi**2 / 2.0), rel=1e-12)
    assert p[0] > 0.0, "the anchor is still silent — the defect is not fixed"


def test_Z4_no_evidence_at_all_still_emits_NOTHING():
    """⚠ The control that stops Z4 becoming "everything now speaks". An object with no composition
    evidence has ``Var(log f) = ∞``, and ``1/(∞ + trigamma)`` is still exactly 0. The fix restores
    the CERTAIN objects, not the ignorant ones."""
    v_fg, _ = own_composition_logvar(np.array([0.4]), np.array([0.0]), np.array([False]))
    assert not np.isfinite(v_fg[0])
    assert own_precision(np.array([0.0]), v_fg, np.array([True]))[0] == 0.0
    assert own_precision(np.array([5_000.0]), v_fg, np.array([True]))[0] == 0.0


# ── Z5 — ⭐ THE LOCATION MUST NOT MOVE (the design error the mass pin caught) ────────────────────


def test_Z5_the_composition_identity_is_EXACT_so_the_location_may_not_gain_the_half():
    """⛔⛔ **THE GATE THAT KILLED THE FIRST VERSION OF THIS FIX.** The three components are SHARES of one
    observed total, so ``Σ_c ρ_c·E_c = M`` exactly — that identity is what the relay's mass pin enforces
    and what makes three relayed densities a composition rather than three unrelated numbers.

    Giving each component the ``Gamma(a+½, E)`` posterior MEAN breaks it by exactly ``3/2``: correct for
    one rate in isolation, wrong for a share of a total. ``test_relay_mass_pin`` caught it as
    ``R_own = 0.5 + 1/M = 0.5000082``. ⭐ The half belongs to the rate's VARIANCE (`count_logvar`); the
    relay fuses in LINEAR density space, so ``ρ = 0`` was always expressible and only its precision was
    broken."""
    M, E_g, E_r = 121_000.0, 13_284.0, 11_500.0
    f_g, f_pos, f_neg = 0.5, 0.25, 0.25
    rho = [f_g * M / E_g, f_pos * M / E_r, f_neg * M / E_r]
    assert rho[0] * E_g + (rho[1] + rho[2]) * E_r == pytest.approx(M, rel=1e-12)
    # and the rejected form does NOT close, by exactly 3/2 — so this gate is not vacuous
    bad = [(f_g * M + 0.5) / E_g, (f_pos * M + 0.5) / E_r, (f_neg * M + 0.5) / E_r]
    assert bad[0] * E_g + (bad[1] + bad[2]) * E_r == pytest.approx(M + 1.5, rel=1e-12)


# ── Z6 — monotone across the join ───────────────────────────────────────────────────────────────


def test_Z6_precision_is_MONOTONE_in_the_count_across_zero():
    """⭐ There must be no step at ``n = 0 → 1``: the old code jumped from "silent" to "as precise as
    the count allows", and a discontinuity there is what let a whole population fall off the edge."""
    n = np.arange(0.0, 40.0)
    v0 = np.zeros_like(n)
    p = own_precision(n, v0, np.ones(n.shape, bool))
    assert np.all(np.diff(p) > 0.0), p[:6]
    assert p[0] > 0.0


# ── Z7 — ⭐⭐ THE PRODUCTION PATH, not just the arithmetic (the hole P3 exposed) ──────────────────


def _init_with(geometry, chain, statics, belief):
    from rigel.calibration.region_init import build_region_init

    return build_region_init(
        chain,
        statics,
        geometry,
        kappa=0.9,
        od_g=0.2,
        od_r=0.1,
        n_gdna_obs=230.0,
        n_rna_obs=85.0,
        n_grid=60,
        logodds_window=10.0,
        n_tilt=None,
        n_grid_ss=256,
        belief=belief,
    )


def test_Z7_build_region_init_lets_a_ZERO_COUNT_STRUCTURALLY_CERTAIN_slot_EMIT():
    """⛔⛔ **THE GATE THE FIRST PERTURBATION PASS WAS MISSING.** Z1–Z6 all call ``own_precision``
    directly with ``live=True`` handed in, so re-gating ``build_region_init``'s ``live`` back onto the
    COUNT (``rho_g > _EPS``) fired **nothing** — and that is the half of the fix that actually reaches
    the intergenic anchor. TRAPS: perturb-every-gate/TRAPS: name-the-observable-per-site: name the observable for *each place* the change was made.

    ⭐ The fixture is the shipped one with ONE thing varied — an intergenic slot's counts zeroed — so it
    is the ``g00`` intergenic region in miniature: ``struct_lock`` (composition certain), zero counts,
    positive gDNA opportunity. Under the retired count-keyed predicate its ``rho_g`` is 0, so
    ``rho_g > _EPS`` is False and it is silent. That silence, times 1,298, is TRAPS: a-zero-count-is-a-measurement.
    """
    import dataclasses

    import test_region_init as TNI  # the shipped fixture, never a second one

    chain, statics, geometry, belief, _ = TNI._scenario(0.9)
    base = _init_with(geometry, chain, statics, belief)

    count = np.asarray(geometry.unspliced_count, np.float64)
    eff_g = np.asarray(geometry.eff_gdna, np.float64)
    cand = np.flatnonzero(np.asarray(base.struct_lock, bool) & (count.sum(1) > 0) & (eff_g > 0.0))
    assert cand.size, "the fixture has no structurally-certain slot with counts to zero out"
    slot = int(cand[0])
    assert base.prec_g[slot] > 0.0, "the control slot was already silent WITH counts"

    emptied = count.copy()
    emptied[slot] = 0.0
    ni = _init_with(dataclasses.replace(geometry, unspliced_count=emptied), chain, statics, belief)

    assert bool(np.asarray(ni.struct_lock, bool)[slot]), (
        "zeroing the count unlocked it — bad fixture"
    )
    assert ni.rho_g[slot] == 0.0, "the location must stay exactly zero (Z5)"
    assert ni.prec_g[slot] > 0.0, (
        "a certain, zero-count, positive-opportunity slot is still silent — `live` is keyed on the "
        "count again"
    )
    # ⭐ and it speaks QUIETLY: a zero count must be the weakest statement, not a free certainty.
    assert ni.prec_g[slot] < base.prec_g[slot]
