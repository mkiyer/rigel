"""``apply_grouped_prior_update`` — the identity the per-component prior rests on.

For a GIVEN ``raw_counts`` vector the function must leave the gDNA:RNA split exactly where the two
pseudocounts put it::

    out[gdna]           == raw[gdna] + gdna_prior
    Σ_{i ≠ gdna} out[i] == Σ_{i ≠ gdna} raw[i] + rna_prior

so that however the RNA share is allocated INSIDE the pool, it never moves the number calibration
exists to produce. It is a PER-CALL identity and not an end-to-end one, and conflating the two has
produced a wrong test more than once: the EM iterates around this function, so a larger ``rna_prior``
shifts ``theta``, hence the E-step, hence the next iteration's ``raw_counts``, hence the converged
gDNA total. "The library gDNA fraction cannot move" is false by design. Everything below fixes
``raw_counts`` and calls the function ONCE, which is the only regime the identity claims. The function
is ``static`` in the C++, so the identity was untestable until ``_apply_grouped_prior_update_test``,
the test-only binding these gates reach it through.

THE ALLOCATION (owner, 2026-09-19; ``EQUATIONS.md`` §9b): the RNA pseudocount is shared over the
locus's RNA components in proportion to the evidence each already carries, and NO component is singled
out for zero. Its predecessor withheld the share from SYNTHETIC nascent entities, which made the
prior's factor un-common over the pool and let the prior alone redistribute RNA; the eligibility test
and the per-component flag that fed it are gone. What survives is the property that made the weights
the right ones in the first place, and it is now unconditional: the allocation echoes the EM's own
belief, so it says nothing about the within-RNA split, and ``raw[i] == 0`` is absorbing.

PERTURBATION: every edit is rebuilt and re-run, and each is first scored against a random configuration
battery so that "the tests did not catch it" is separated from "it could not be caught"
(TRAPS: could-the-arm-have-fired). Re-run against the restored allocation, 2026-09-19: re-introducing an
eligibility test fires 8 gates, an equal share per component fires 13, leaking a tenth of the RNA
pseudocount onto the gDNA component fires 5 (the identity below among them), reading a weight without its
guard fires 1, and naming the wrong total in the recipients gate fires 7.

⛔ ONE PERTURBATION IS PROVABLY INERT AND THAT IS A FINDING ABOUT THE CODE, not a hole here: dropping the
no-recipient gate moves nothing, because whenever it would zero `rna_prior` the total it divides by is
itself zero and the branch below never runs. It is defensive redundancy, kept because it states the intent
where the intent is decided — do not "cover" it by weakening a gate, since no input distinguishes the two.
The battery needed falsifying too: its first version generated only well-formed weights, so a perturbation
removing the negative-weight guard moved nothing and read as inert — a verdict about the generator wearing
the costume of a verdict about the code.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel._em_impl import EM_LOG_EPSILON, _apply_grouped_prior_update_test

_EMPTY_F64 = np.zeros(0, dtype=np.float64)


def _update(
    raw,
    *,
    carried=None,
    weight=None,
    gdna_prior=0.0,
    rna_prior=0.0,
    gdna_index=0,
    has_gdna=True,
):
    return _apply_grouped_prior_update_test(
        raw_counts=np.ascontiguousarray(raw, dtype=np.float64),
        carried_state=_EMPTY_F64 if carried is None else np.ascontiguousarray(carried, np.float64),
        rna_prior_weight=(
            _EMPTY_F64 if weight is None else np.ascontiguousarray(weight, np.float64)
        ),
        gdna_prior_fragments=float(gdna_prior),
        rna_prior_fragments=float(rna_prior),
        gdna_index=int(gdna_index),
        has_gdna_candidate=bool(has_gdna),
    )


def _split(out, gdna_index):
    """``(gdna_total, rna_total)`` — the two quantities the identity is about."""
    mask = np.ones(out.shape[0], dtype=bool)
    mask[gdna_index] = False
    return float(out[gdna_index]), float(out[mask].sum())


# ──────────────────────────────────────────────────────────────────────────────
# The identity itself
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "raw",
    [
        [100.0, 30.0, 70.0],
        [0.0, 30.0, 70.0],
        [100.0, 1e-9, 70.0],
        [100.0, 0.0, 70.0],
        [5.0, 1234.5, 0.0, 17.25, 900.0],
    ],
    ids=["plain", "no_gdna_count", "tiny", "one_zero", "five"],
)
def test_the_gDNA_RNA_split_lands_exactly_where_the_two_pseudocounts_put_it(raw):
    """The whole point: gDNA gets ``+gdna_prior`` and the RNA POOL gets ``+rna_prior``, whatever the
    allocation does INSIDE the pool.

    ⭐ UNCONDITIONAL. It used to carry a ``not eligible`` row — an all-synthetic locus, where the
    prior had nowhere to land and the pool got ``+0`` — and that row is gone with the eligibility test
    it described. Every RNA component is a recipient now, so the only way the pool can fail to receive
    the prior is to have no RNA evidence at all, which is a different gate below.
    """
    raw = np.asarray(raw, dtype=np.float64)
    gdna_prior, rna_prior = 12.5, 40.0
    out = _update(raw, gdna_prior=gdna_prior, rna_prior=rna_prior)

    g_out, r_out = _split(out, 0)
    g_raw, r_raw = _split(raw, 0)
    assert g_out == pytest.approx(g_raw + gdna_prior, rel=1e-12, abs=1e-12)
    assert r_out == pytest.approx(r_raw + rna_prior, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize(
    "raw",
    [
        [100.0, 30.0, 70.0],
        [100.0, 1.0, 999.0],
        [5.0, 1234.5, 17.25, 900.0, 0.25],
    ],
    ids=["plain", "skewed", "five"],
)
def test_the_prior_changes_NO_component_s_SHARE_of_the_RNA_pool(raw):
    """The restored rule, stated as the property that makes it a rule rather than a special case:
    the allocation echoes the EM's own current belief, so it moves the gDNA:RNA split and NOTHING
    else. Every RNA component's share of the pool is exactly what it was.

    This is the gate no eligibility test can pass. Hold ONE component out and the prior's factor
    stops being common over the pool, so the prior alone redistributes RNA — which is the defect the
    restoration removed (`EQUATIONS.md` §9b).
    """
    raw = np.asarray(raw, dtype=np.float64)
    out = _update(raw, gdna_prior=7.0, rna_prior=250.0)
    _g_out, r_out = _split(out, 0)
    _g_raw, r_raw = _split(raw, 0)
    np.testing.assert_allclose(out[1:] / r_out, raw[1:] / r_raw, rtol=1e-12, atol=0.0)


def test_a_zero_count_component_CANNOT_be_revived_by_the_prior():
    """``out[i]`` is PROPORTIONAL to ``raw[i]``, so ``raw[i] == 0`` stays 0 — an ABSORBING STATE no
    prior magnitude escapes. That is the structural guard against reviving a shadow entity, it is a
    property of the WEIGHTS rather than of any eligibility test, and it is why admitting every RNA
    component at these weights costs nothing on that axis (``EQUATIONS.md`` §9b.1). The
    coverage-weighted warm start remains the only spark such an entity gets."""
    raw = np.array([100.0, 0.0, 70.0], dtype=np.float64)
    out = _update(raw, gdna_prior=0.0, rna_prior=500.0)
    assert out[1] == 0.0


def test_the_carried_state_path_is_live_and_obeys_the_SAME_split():
    """The zero-evidence branch, reachable under VBEM (the shipped default), which passes ``alpha`` as
    the carried state. Its gate must name the denominator ITS branch divides by: testing one total
    while dividing by another lets a locus keep a live ``rna_prior``, multiply it by ``inv = 0`` and
    drop it — RNA summing to ``rna_count`` while gDNA keeps ``gdna_prior``."""
    raw = np.array([100.0, 0.0, 0.0], dtype=np.float64)  # no RNA evidence at all
    carried = np.array([1.0, 2.0, 6.0], dtype=np.float64)
    out = _update(raw, carried=carried, gdna_prior=3.0, rna_prior=40.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(100.0 + 3.0, rel=1e-12)
    assert r_out == pytest.approx(0.0 + 40.0, rel=1e-12)
    # shared out in proportion to the CARRIED state, not the (zero) counts
    assert out[1] / out[2] == pytest.approx(2.0 / 6.0, rel=1e-12)


def test_the_carried_state_branch_shares_the_pool_over_EVERY_RNA_component():
    """The prior-only pool obeys the same rule as the evidence one: shared over every RNA component,
    here in proportion to the carried alpha, none singled out for zero. A component with zero carried
    alpha still receives nothing — the absorbing state, in the carried coordinate."""
    raw = np.array([100.0, 0.0, 0.0, 0.0], dtype=np.float64)
    carried = np.array([1.0, 3.0, 5.0, 2.0], dtype=np.float64)
    out = _update(raw, carried=carried, gdna_prior=0.0, rna_prior=40.0)

    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(100.0, rel=1e-12)
    assert r_out == pytest.approx(40.0, rel=1e-12), "the pool did not sum to rna_count + rna_prior"
    # shared in proportion to the carried alpha (3 : 5 : 2)
    np.testing.assert_allclose(out[1:] / r_out, carried[1:] / carried[1:].sum(), rtol=1e-12)


def test_a_locus_with_NO_rna_evidence_AT_ALL_drops_the_rna_prior():
    """The configuration the C++ ``prior_recipients`` gate exists for: zero RNA COUNT and zero RNA
    CARRIED alpha, so neither branch below it can place the prior. It must be dropped rather than
    computed against a zero denominator and silently multiplied away.

    ⚠ gDNA KEEPS its own pseudocount here, and that is correct: the two are separate pseudocounts on
    separate components, and the identity this file gates is that each lands where it was aimed."""
    raw = np.array([100.0, 0.0, 0.0], dtype=np.float64)
    carried = np.array([1.0, 0.0, 0.0], dtype=np.float64)  # the gDNA component's alpha only
    out = _update(raw, carried=carried, gdna_prior=3.0, rna_prior=40.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(103.0, rel=1e-12)
    assert r_out == pytest.approx(0.0, abs=1e-12), "a prior with no recipient at all was paid out"


# ──────────────────────────────────────────────────────────────────────────────
# The allocation is UNIFORM over the pool — a property, over many shapes
# ──────────────────────────────────────────────────────────────────────────────


def test_NO_configuration_of_counts_lets_the_prior_move_the_within_RNA_split():
    """The structural replacement for the byte-identity pair that pinned the synthetic mask
    (TRAPS: byte-identity-gate). Those two could only compare a mask against no mask, and the mask is
    gone — so the gate becomes the property they were protecting, asserted over the space instead of
    at one point: for ANY counts, ANY pool size and ANY prior magnitude, the RNA prior scales the
    whole pool by one factor and singles nobody out.

    ⛔ EXACT, not ``approx``. The prior enters as one common multiply, so the ratio between any two
    RNA components survives it to the bit; a rule that reads anything per component would not."""
    rng = np.random.default_rng(0)
    for _ in range(300):
        n = int(rng.integers(2, 9))
        raw = rng.uniform(0.0, 1000.0, size=n)
        raw[0] = rng.uniform(0.0, 1000.0)  # the gDNA component
        rna_prior = float(rng.choice([0.0, 1e-6, 1.0, 37.0, 1e5]))
        out = _update(raw, gdna_prior=11.0, rna_prior=rna_prior)
        _g_out, r_out = _split(out, 0)
        _g_raw, r_raw = _split(raw, 0)
        if r_raw <= 0.0:
            assert r_out == 0.0
            continue
        assert np.all(out[1:][raw[1:] == 0.0] == 0.0), "a zero-evidence component was paid"
        # Every component's SHARE of the pool survives the prior. Compared as shares rather than as
        # the raw ratio `out[i]/raw[i]`: that ratio reintroduces a division this function never
        # performs, so its last bit is the TEST's arithmetic and not the code's.
        np.testing.assert_allclose(out[1:] / r_out, raw[1:] / r_raw, rtol=1e-12, atol=0.0)


# ──────────────────────────────────────────────────────────────────────────────
# Degenerate inputs
# ──────────────────────────────────────────────────────────────────────────────


def test_no_gdna_candidate_withholds_only_the_gDNA_pseudocount():
    """``has_gdna_candidate == False`` withholds the gDNA pseudocount and nothing else.

    The gDNA one must be withheld — there is no component to put it on — but the RNA one lands on the
    RNA components, which exist regardless, so zeroing both withholds the whole RNA prior from every
    locus whose fragments are ALL SPLICED (a gDNA candidate is appended to every unspliced unit).

    That defect is invisible under the shipped evidence-proportional weights, which make the RNA prior
    a COMMON factor over the eligible components that cancels under normalisation — so at such a locus
    it genuinely cannot move ``theta``. An informative per-component weight cancels nothing, which is
    why the gate is worth keeping.
    """
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    out = _update(raw, gdna_prior=12.0, rna_prior=40.0, has_gdna=False)

    assert out[0] == 0.0, "a gDNA pseudocount was paid with no gDNA component to receive it"
    assert float(out[1:].sum()) == pytest.approx(30.0 + 70.0 + 40.0, rel=1e-12), (
        "the RNA prior was discarded at a locus with no gDNA candidate"
    )
    # shared out by evidence, as everywhere else
    assert out[1] / out[2] == pytest.approx(30.0 / 70.0, rel=1e-12)


def test_a_NEGATIVE_or_NONFINITE_raw_count_is_read_as_zero():
    """``nonnegative_finite`` guards every read. A NaN must not propagate into the pool total."""
    raw = np.array([100.0, np.nan, -5.0, 70.0], dtype=np.float64)
    out = _update(raw, gdna_prior=0.0, rna_prior=30.0)
    assert np.all(np.isfinite(out))
    _g, r_out = _split(out, 0)
    assert r_out == pytest.approx(70.0 + 30.0, rel=1e-12)


def test_the_binding_REFUSES_a_mismatched_auxiliary_length():
    """An auxiliary array of the wrong length would read past the end of the components. Empty means
    'this locus has none'; any other wrong length is an error, not a shrug."""
    raw = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    with pytest.raises(Exception):
        _update(raw, carried=[1.0, 2.0])
    with pytest.raises(Exception):
        _update(raw, weight=[1.0, 2.0])


# ──────────────────────────────────────────────────────────────────────────────
# THE WEIGHTED ALLOCATION — the per-transcript lane
# ──────────────────────────────────────────────────────────────────────────────


def test_a_weight_PROPORTIONAL_TO_EVIDENCE_reproduces_the_shipped_rule():
    """The identity that shows the two rules are one rule. The shipped update

        out[i] = raw[i] · (1 + rna_prior/annotated_count)

    is exactly ``raw[i] + rna_prior·w_i/Σw`` at ``w_i = raw[i]``. So handing this branch the raw counts
    as weights must reproduce the shipped answer. Not byte-identical — the two arrive at it by a
    different operation order — but equal to floating tolerance, and that is the point: the
    per-transcript lane GENERALISES the shipped rule rather than replacing it, and the only thing a new
    allocation changes is the weights.
    """
    rng = np.random.default_rng(7)
    raw = rng.uniform(1.0, 1000.0, size=6)
    shipped = _update(raw, gdna_prior=5.0, rna_prior=250.0)
    weighted = _update(raw, weight=raw, gdna_prior=5.0, rna_prior=250.0)
    np.testing.assert_allclose(weighted, shipped, rtol=1e-12, atol=1e-9)


def test_the_weighted_branch_obeys_the_SAME_conservation_identity():
    """Conservation is what makes an allocation safe to change: whatever the weights say, the RNA pool
    still sums to ``rna_count + rna_prior``, so the gDNA:RNA split does not move with the allocation."""
    rng = np.random.default_rng(11)
    for _ in range(200):
        n = int(rng.integers(2, 8))
        raw = rng.uniform(0.0, 500.0, size=n)
        w = rng.uniform(0.0, 10.0, size=n)
        w[rng.random(n) < 0.3] = 0.0  # a weight vector may nominate only some components
        out = _update(raw, weight=w, gdna_prior=3.0, rna_prior=77.0)
        g_out, r_out = _split(out, 0)
        g_raw, r_raw = _split(raw, 0)
        assert g_out == pytest.approx(g_raw + 3.0, rel=1e-11)
        # A weight vector nominating NOBODY falls back to the evidence-proportional rule rather
        # than dropping the prior — the gate below states why — so the pool still receives it
        # whenever it has evidence to receive it with.
        expect = r_raw + (77.0 if (w[1:].sum() > 0.0 or r_raw > EM_LOG_EPSILON) else 0.0)
        assert r_out == pytest.approx(expect, rel=1e-11, abs=1e-9)


def test_a_ZERO_COUNT_component_CAN_be_revived_by_a_weighted_prior():
    """The one place the two allocations are not interchangeable, and it is the consequential one.

    Under the shipped evidence-proportional weights ``out[i] = 0`` is an ABSORBING STATE: the prior is
    proportional to ``raw[i]``, so no prior magnitude whatsoever can revive a component with no
    warm-start evidence. A strictly positive weight has no such state.

    This is a capability, not automatically a good thing — it is precisely the mechanism that could
    revive a shadow entity the data does not support, so it is what an allocation rule must earn.
    """
    raw = np.array([100.0, 0.0, 70.0], dtype=np.float64)
    assert _update(raw, rna_prior=500.0)[1] == 0.0, "the shipped rule should not revive it"
    revived = _update(raw, weight=[0.0, 1.0, 1.0], rna_prior=500.0)
    assert revived[1] == pytest.approx(250.0, rel=1e-12)


def test_a_weight_of_ZERO_is_how_a_caller_withholds_the_prior_from_ONE_component():
    """There is no eligibility test any more, so the weight vector IS the only place a caller can say
    "not this one" — and it says it by naming a weight of zero, per component, on the record.

    ⚠ That is deliberately not the same statement as the rule this replaced. Withholding here is the
    CALLER's claim about one transcript, visible in the array it passed; the rule that went withheld
    from a whole class on the solver's own initiative, for a reason the caller never saw.
    """
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    out = _update(raw, weight=[0.0, 1.0, 0.0], rna_prior=50.0)
    assert out[2] == pytest.approx(70.0, rel=1e-12), "a zero-weight component took prior mass"
    assert out[1] == pytest.approx(30.0 + 50.0, rel=1e-12), "the whole prior did not land"


def test_an_ALL_ZERO_weight_vector_FALLS_BACK_rather_than_dropping_the_prior():
    """A weight vector that nominates nobody has nothing to say, so the shipped rule stands. The
    alternative — paying gDNA its pseudocount while the RNA pool silently loses one — is the defect
    the ``prior_recipients`` gate exists for, and it MOVES the split."""
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    fallback = _update(raw, weight=np.zeros(3), gdna_prior=4.0, rna_prior=50.0)
    shipped = _update(raw, gdna_prior=4.0, rna_prior=50.0)
    np.testing.assert_array_equal(fallback, shipped)


def test_a_NEGATIVE_or_NONFINITE_weight_is_read_as_zero():
    """Weights go through the same ``nonnegative_finite`` guard as counts — twice, so a NaN cannot
    poison the normaliser and then be divided by it."""
    raw = np.array([100.0, 30.0, 70.0, 10.0], dtype=np.float64)
    out = _update(raw, weight=[0.0, np.nan, -3.0, 2.0], rna_prior=50.0)
    assert np.all(np.isfinite(out))
    assert out[1] == pytest.approx(30.0, rel=1e-12)
    assert out[2] == pytest.approx(70.0, rel=1e-12)
    assert out[3] == pytest.approx(10.0 + 50.0, rel=1e-12)


def test_the_weighted_branch_needs_NO_carried_state_to_place_a_prior_only_pool():
    """The carried state answers only one question — who receives the prior when there is no evidence?
    An explicit weight answers it directly, so the weighted branch has no zero-evidence special case,
    and it works where the shipped rule needs VBEM's alpha threaded through."""
    raw = np.zeros(4, dtype=np.float64)
    raw[0] = 100.0  # gDNA only
    out = _update(raw, weight=[0.0, 3.0, 1.0, 0.0], gdna_prior=2.0, rna_prior=40.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(102.0, rel=1e-12)
    assert r_out == pytest.approx(40.0, rel=1e-12)
    assert out[1] == pytest.approx(30.0, rel=1e-12)
    assert out[2] == pytest.approx(10.0, rel=1e-12)
    assert out[3] == 0.0
