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
    counts, diag = _realized_gdna_counts(payload, opp, rl, rt, rna_pmf, uniform)
    assert diag.applied
    assert diag.ontarget_share == pytest.approx(0.0, abs=1e-9)


def test_enriched_boundaries_raise_the_on_target_share():
    """At 50x enrichment the fixture's own arithmetic puts the excess classes near
    ``rho·49·E_contained(200)·2`` of a ~2.6k total — about 0.27. Gate the ORDER, derived, not a guess:
    well clear of the no-excess case's exact 0, and the boundary share alive beside it."""
    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=50.0)
    counts, diag = _realized_gdna_counts(payload, opp, rl, rt, rna_pmf, uniform)
    assert diag.applied
    assert diag.ontarget_share > 0.2
    assert diag.boundary_share > 0.05
    _, diag1 = _realized_gdna_counts(
        *_payload_fixture(boundary_excess=1.0)[:1], *_payload_fixture(boundary_excess=1.0)[1:]
    )
    assert diag.ontarget_share > 100 * max(diag1.ontarget_share, 1e-6)


def test_zero_gdna_declines_rather_than_fabricating_a_law():
    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=1.0)
    import dataclasses

    rc = np.zeros_like(np.asarray(payload.region_contained_count))
    pl = np.asarray(payload.pool_lengths, dtype=np.float64).copy()
    pl[:4] = 0.0
    starved = dataclasses.replace(payload, region_contained_count=rc, pool_lengths=pl)
    counts, diag = _realized_gdna_counts(starved, opp, rl, rt, rna_pmf, uniform)
    assert counts is None and not diag.applied
