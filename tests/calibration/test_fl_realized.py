"""The TWO ESTIMANDS of the gDNA fragment-length model, and the routing guarantee between them.

`gdna_pmf` is the UNIFORM-FRAME law (what the chemistry makes; what the opportunity/prior mathematics
assumes). `gdna_realized_pmf` is the LIBRARY-CENSUS law (what a sequenced gDNA fragment looks like,
capture selection included; what the EM's per-fragment scorer conditions on). Off capture they coincide;
under capture they split, and feeding either consumer the other one's estimand was measured at
+188,208 misassigned transcripts on one ladder row.

⛔ **The gates here protect the ROUTING, not just the estimator**: the geometry consumers must be
BIT-identical whether or not the realized law is computed, the realized law must fall back to the
uniform one exactly whenever it cannot be estimated, and the on-target correction must vanish
IDENTICALLY when the boundaries carry no enrichment excess.

⚠ **Perturbation record** (`TRAPS: perturb-every-gate`): leaking the realized law into `gdna_pmf`,
breaking the fallback, and dropping the excess's `−1` each fire two gates. Removing the estimator's
early no-density guard fires NOTHING here — with `rho_off = 0` every mass term multiplies to zero and
the mass-zero decline catches it downstream — so that guard's load-bearing work is the honest
diagnostic string and the (near-unreachable) unbracketed-fit branch, and it is kept for those, with
this note standing in place of a gate that cannot be made to fire on this fixture.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.fl import _fl_models_from_histograms


def _models(**kw):
    rng = np.random.default_rng(5)
    g = rng.random(101)
    r = rng.random(101)
    return _fl_models_from_histograms(
        global_counts=g + r, rna_counts=r, gdna_counts=g, max_size=100, **kw
    )


# ── the field exists and is never None: the scorer must be able to read it unconditionally ───────


def test_the_realized_law_is_always_present():
    m = _models()
    assert m.gdna_realized_pmf is not None
    assert m.gdna_realized_pmf.shape == m.gdna_pmf.shape


def test_without_a_realized_estimate_the_two_estimands_coincide_exactly():
    """⭐ The fallback IS the old behaviour: one law, byte-equal, so a consumer that reads the realized
    field on an off-capture or input-starved build gets exactly what shipped before."""
    m = _models()
    np.testing.assert_array_equal(m.gdna_realized_pmf, m.gdna_pmf)


def test_a_supplied_realized_histogram_is_shrunk_like_its_sibling():
    rng = np.random.default_rng(9)
    realized = rng.random(101) * 50
    m = _models(gdna_realized_counts=realized)
    assert not np.array_equal(m.gdna_realized_pmf, m.gdna_pmf)
    assert m.gdna_realized_pmf.min() > 0.0  # EB-shrunk toward the global anchor, like gdna_pmf
    assert m.gdna_realized_pmf.sum() == pytest.approx(1.0)


def test_the_uniform_law_is_bit_identical_with_and_without_the_realized_input():
    """⛔⛔ THE ROUTING GUARANTEE. Computing the realized law must not perturb `gdna_pmf` by one ULP —
    every geometry consumer reads it, and the +188k regression was geometry eating the wrong estimand."""
    rng = np.random.default_rng(9)
    a = _models()
    b = _models(gdna_realized_counts=rng.random(101) * 50)
    np.testing.assert_array_equal(a.gdna_pmf, b.gdna_pmf)
    np.testing.assert_array_equal(a.rna_pmf, b.rna_pmf)
    np.testing.assert_array_equal(a.global_pmf, b.global_pmf)


# ── the realized estimator's own invariants, on constructed payloads ─────────────────────────────


def _payload_fixture(boundary_excess: float):
    """A two-reference toy: intergenic/intron/exon regions with uniform-consistent contained counts,
    and intron|exon + intergenic|exon boundaries whose counts carry ``boundary_excess`` TIMES the
    uniform expectation. At 1.0 the boundaries say "no capture"."""
    from tests.calibration.fl_realized_fixture import build_fixture

    return build_fixture(boundary_excess)


def test_no_enrichment_excess_means_no_on_target_correction():
    """⭐ The closure property: boundaries consistent with the uniform field ⇒ the (eps−1)+ term is
    identically zero and the realized law is the sampled blend alone."""
    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=1.0)
    counts, _uniform_out, diag = _realized_gdna_counts(payload, opp, rl, rt, rna_pmf, uniform)
    assert diag.applied
    assert diag.ontarget_share == pytest.approx(0.0, abs=1e-9)


def test_enriched_boundaries_raise_the_on_target_share():
    """At 50x enrichment the fixture's own arithmetic puts the excess classes near
    ``rho·49·E_contained(200)·2`` of a ~2.6k total — about 0.27. Gate the ORDER, derived, not a guess:
    well clear of the no-excess case's exact 0, and the boundary share alive beside it."""
    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=50.0)
    counts, _uniform_out, diag = _realized_gdna_counts(payload, opp, rl, rt, rna_pmf, uniform)
    assert diag.applied
    assert diag.ontarget_share > 0.2
    assert diag.boundary_share > 0.05


def test_the_on_target_correction_rises_smoothly_with_enrichment():
    """⭐ The capture SPECTRUM, not a switch: the correction is EXACTLY 0 with no excess and rises
    monotonically from there. Its weight is its OWN resolution — how far an exon's enrichment sits from
    1 against its own sampling error — and NOT the strata-split weight `lam`. Gating it on `lam` was a
    real design error caught by this fixture: `lam` asks whether the two strata's LAWS differ, which is
    a different question, and using it here suppressed a correction that was genuinely resolved."""
    from rigel.calibration.fl import _realized_gdna_counts

    shares = [
        _realized_gdna_counts(*_payload_fixture(boundary_excess=x))[2].ontarget_share
        for x in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0)
    ]
    assert shares[0] == 0.0
    assert all(b > a for a, b in zip(shares, shares[1:]))


def test_an_unresolved_enrichment_weight_collapses_with_its_precision():
    """⛔ The excess's own weight, gated at the MECHANISM rather than end to end — and the reason is
    worth recording. Every way of starving the enrichment's precision in the fixture (thinning the
    boundaries, or the whole library at fixed enrichment) also trips an EARLIER guard, so the
    end-to-end share reads 0 whether or not this weight exists: that perturbation cannot isolate it.
    What CAN be isolated is the weight itself, which is where the logic lives.

    ⚠ Traced on the fixture's own numbers: 98 gDNA crossings at a ratio of 10 give w = 0.988, and the
    same ratio on 1e-3 crossings gives w = 0.0008. The ratio alone would fire identically at both.
    """
    from rigel.calibration.fl import _resolution_weight

    eps = 10.0
    strong = _resolution_weight((eps - 1.0) ** 2, eps * eps / 98.0)
    weak = _resolution_weight((eps - 1.0) ** 2, eps * eps / 1e-3)
    assert strong > 0.95
    assert weak < 0.01
    # and it is monotone in the evidence, with no step anywhere along it
    ws = [_resolution_weight((eps - 1.0) ** 2, eps * eps / n) for n in np.logspace(-4, 4, 200)]
    assert all(b >= a for a, b in zip(ws, ws[1:]))
    assert max(abs(b - a) for a, b in zip(ws, ws[1:])) < 0.05


def test_zero_gdna_declines_rather_than_fabricating_a_law():
    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=1.0)
    import dataclasses

    rc = np.zeros_like(np.asarray(payload.region_contained_count))
    pl = np.asarray(payload.pool_lengths, dtype=np.float64).copy()
    pl[:4] = 0.0
    starved = dataclasses.replace(payload, region_contained_count=rc, pool_lengths=pl)
    counts, _uniform_out, diag = _realized_gdna_counts(starved, opp, rl, rt, rna_pmf, uniform)
    assert counts is None and not diag.applied


# ── the CONVERGENCE law: no cliffs, and the two estimands merge when the split is unresolvable ───


def test_resolution_weight_is_a_smooth_signal_to_noise_ratio():
    """⭐ `lam = S/(S+N)`: 0 when the split is pure noise, 1 when noise vanishes, ½ at S = N, and
    MONOTONE in between. No threshold, so there is no value of the data at which behaviour jumps."""
    from rigel.calibration.fl import _resolution_weight

    assert _resolution_weight(0.0, 1.0) == 0.0
    assert _resolution_weight(1.0, 0.0) == 1.0
    assert _resolution_weight(1.0, 1.0) == pytest.approx(0.5)
    prev = -1.0
    for s in np.linspace(0.0, 10.0, 60):
        lam = _resolution_weight(float(s), 1.0)
        assert 0.0 <= lam <= 1.0 and lam >= prev
        prev = lam


def test_the_estimands_converge_when_the_boundary_stratum_is_starved():
    """⛔ THE OWNER'S CASE, capture-OFF side: sparse boundary data must not let the two laws diverge —
    the difference is unmeasurable there, so the honest answer is that they agree."""
    from rigel.calibration.fl import _couple_estimands

    rng = np.random.default_rng(3)
    g_c = np.abs(rng.random(60)) + 0.01
    g_c /= g_c.sum()
    g_b = np.roll(g_c, 6)  # a real shift, but measured on almost nothing
    uni, real, lam = _couple_estimands(g_c, 1e6, g_b, 1e-3)
    assert lam < 1e-3
    # ⭐ The guarantee is that the residual disagreement is BOUNDED BY lam times the split, not that it
    # is bitwise zero — asserting equality would be asserting more than the mathematics gives.
    assert float(np.abs(uni - real).sum()) <= lam * float(np.abs(g_b - g_c).sum()) + 1e-12


def test_the_estimands_converge_when_the_CONTAINED_stratum_is_starved():
    """⛔ THE OWNER'S CASE, infinite-capture side: with no off-target data the chemistry law is not
    estimable, so it must borrow the only law that IS — not hold a stale value."""
    from rigel.calibration.fl import _couple_estimands

    rng = np.random.default_rng(4)
    g_c = np.abs(rng.random(60)) + 0.01
    g_c /= g_c.sum()
    g_b = np.roll(g_c, 6)
    uni, real, lam = _couple_estimands(g_c, 1e-3, g_b, 1e6)
    assert lam < 1e-3
    assert float(np.abs(uni - real).sum()) <= lam * float(np.abs(g_b - g_c).sum()) + 1e-12
    # and the common law it converged on is the well-measured one, not the starved one
    assert float(np.abs(uni - g_b).sum()) < float(np.abs(uni - g_c).sum())


def test_a_well_measured_split_keeps_the_two_estimands_apart():
    from rigel.calibration.fl import _couple_estimands

    rng = np.random.default_rng(5)
    g_c = np.abs(rng.random(60)) + 0.01
    g_c /= g_c.sum()
    g_b = np.roll(g_c, 6)
    uni, real, lam = _couple_estimands(g_c, 1e7, g_b, 1e7)
    assert lam > 0.99
    # lam -> 1 recovers the uncoupled behaviour, to within (1 - lam) times the split
    assert float(np.abs(uni - g_c).sum()) <= (1.0 - lam) * float(np.abs(g_b - g_c).sum()) + 1e-12
    assert float(np.abs(uni - real).sum()) > 0.1


def test_there_is_no_cliff_anywhere_along_the_capture_spectrum():
    """⛔⛔ THE ANTI-CLIFF GATE. Sweep the boundary stratum's mass over eight orders of magnitude and
    assert the uniform law moves CONTINUOUSLY — a binary decline would show as a step."""
    from rigel.calibration.fl import _couple_estimands

    rng = np.random.default_rng(6)
    g_c = np.abs(rng.random(60)) + 0.01
    g_c /= g_c.sum()
    g_b = np.roll(g_c, 6)
    masses = np.logspace(-4, 4, 200)
    unis = np.array([_couple_estimands(g_c, 1e3, g_b, float(m))[0] for m in masses])
    steps = np.abs(np.diff(unis, axis=0)).sum(axis=1)
    assert steps.max() < 0.05, f"a step of {steps.max():.3f} is a cliff, not a fade"
