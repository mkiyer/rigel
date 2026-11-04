"""``apply_grouped_prior_update``: the identity the per-component prior rests on.

For a given ``raw_counts`` vector the function leaves the gDNA:RNA split exactly where the two pseudocounts put
it::

    out[gdna]           == raw[gdna] + gdna_prior
    Σ_{i ≠ gdna} out[i] == Σ_{i ≠ gdna} raw[i] + rna_prior

so however the RNA share is allocated inside the pool, it never moves the number calibration exists to produce.
The identity holds per call, not end to end: the EM iterates around this function, so a larger ``rna_prior``
shifts ``theta`` and with it the converged gDNA total. Every test here fixes ``raw_counts`` and calls the function
once. The tests reach it through its test binding, ``_em_impl.testing.apply_grouped_prior_update``.

The RNA pseudocount is shared over the locus's RNA components in proportion to the evidence each already carries,
with no component singled out for zero. The allocation echoes the EM's own belief, so it says nothing about the
within-RNA split, and ``raw[i] == 0`` is absorbing.

The C++ no-recipient gate is redundant with the zero-denominator branch below it: whenever it would zero
``rna_prior``, the total it divides by is zero too, so no input can distinguish the two.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel._em_impl.testing import EM_LOG_EPSILON, apply_grouped_prior_update

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
    return apply_grouped_prior_update(
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
def test_the_gdna_rna_split_lands_exactly_where_the_two_pseudocounts_put_it(raw):
    """gDNA gets ``+gdna_prior`` and the RNA pool gets ``+rna_prior``, whatever the allocation does inside the
    pool. Every RNA component is a recipient, so the pool misses the prior only when it has no RNA evidence at
    all, which a separate test covers."""
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
def test_the_prior_changes_no_components_share_of_the_rna_pool(raw):
    """The allocation echoes the EM's current belief, so it moves the gDNA:RNA split and nothing else: every RNA
    component's share of the pool is unchanged. Holding one component out would make the prior's factor uneven
    over the pool, letting the prior alone redistribute RNA."""
    raw = np.asarray(raw, dtype=np.float64)
    out = _update(raw, gdna_prior=7.0, rna_prior=250.0)
    _g_out, r_out = _split(out, 0)
    _g_raw, r_raw = _split(raw, 0)
    np.testing.assert_allclose(out[1:] / r_out, raw[1:] / r_raw, rtol=1e-12, atol=0.0)


def test_a_zero_count_component_cannot_be_revived_by_the_prior():
    """``out[i]`` is proportional to ``raw[i]``, so ``raw[i] == 0`` stays 0 whatever the prior's size. This
    property of the weights is what keeps the prior from reviving a synthetic span; the coverage-weighted warm
    start is the only source of mass such a component gets."""
    raw = np.array([100.0, 0.0, 70.0], dtype=np.float64)
    out = _update(raw, gdna_prior=0.0, rna_prior=500.0)
    assert out[1] == 0.0


def test_the_carried_state_path_is_live_and_obeys_the_same_split():
    """The zero-evidence branch, reached under VBEM (the default), which passes ``alpha`` as the carried state.
    Its guard must test the denominator this branch divides by; testing another total would let ``rna_prior``
    be multiplied by ``inv = 0`` and lost while gDNA keeps ``gdna_prior``."""
    raw = np.array([100.0, 0.0, 0.0], dtype=np.float64)  # no RNA evidence at all
    carried = np.array([1.0, 2.0, 6.0], dtype=np.float64)
    out = _update(raw, carried=carried, gdna_prior=3.0, rna_prior=40.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(100.0 + 3.0, rel=1e-12)
    assert r_out == pytest.approx(0.0 + 40.0, rel=1e-12)
    # shared out in proportion to the carried state, not the (zero) counts
    assert out[1] / out[2] == pytest.approx(2.0 / 6.0, rel=1e-12)


def test_the_carried_state_branch_shares_the_pool_over_every_rna_component():
    """The prior-only pool follows the same rule as the evidence one: shared over every RNA component in
    proportion to the carried alpha. A component with zero carried alpha receives nothing."""
    raw = np.array([100.0, 0.0, 0.0, 0.0], dtype=np.float64)
    carried = np.array([1.0, 3.0, 5.0, 2.0], dtype=np.float64)
    out = _update(raw, carried=carried, gdna_prior=0.0, rna_prior=40.0)

    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(100.0, rel=1e-12)
    assert r_out == pytest.approx(40.0, rel=1e-12), "the pool did not sum to rna_count + rna_prior"
    # shared in proportion to the carried alpha (3 : 5 : 2)
    np.testing.assert_allclose(out[1:] / r_out, carried[1:] / carried[1:].sum(), rtol=1e-12)


def test_a_locus_with_no_rna_evidence_at_all_drops_the_rna_prior():
    """Zero RNA count and zero carried RNA alpha, the case the C++ ``prior_recipients`` guard handles: neither
    branch can place the RNA prior, so it is dropped rather than divided by zero. gDNA keeps its own
    pseudocount, because the two are separate pseudocounts on separate components."""
    raw = np.array([100.0, 0.0, 0.0], dtype=np.float64)
    carried = np.array([1.0, 0.0, 0.0], dtype=np.float64)  # the gDNA component's alpha only
    out = _update(raw, carried=carried, gdna_prior=3.0, rna_prior=40.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(103.0, rel=1e-12)
    assert r_out == pytest.approx(0.0, abs=1e-12), "a prior with no recipient at all was paid out"


# ──────────────────────────────────────────────────────────────────────────────
# The allocation is uniform over the pool, over many shapes
# ──────────────────────────────────────────────────────────────────────────────


def test_no_configuration_of_counts_lets_the_prior_move_the_within_rna_split():
    """For any counts, pool size and prior magnitude, the RNA prior scales the whole pool by one common factor,
    so every component's share of the pool is preserved. The tolerance is at the last bit because the prior
    enters as one common multiply; a rule that read anything per component would miss it."""
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
        # Compared as shares rather than as `out[i]/raw[i]`: that ratio adds a division the function never
        # performs, so its last bit would measure the test's arithmetic, not the code's.
        np.testing.assert_allclose(out[1:] / r_out, raw[1:] / r_raw, rtol=1e-12, atol=0.0)


# ──────────────────────────────────────────────────────────────────────────────
# Degenerate inputs
# ──────────────────────────────────────────────────────────────────────────────


def test_no_gdna_candidate_withholds_only_the_gdna_pseudocount():
    """``has_gdna_candidate == False`` withholds the gDNA pseudocount and nothing else.

    There is no gDNA component to receive it, but the RNA components exist regardless, so zeroing both would
    withhold the RNA prior from every locus whose fragments are all spliced. Evidence-proportional weights hide
    that defect (the common factor cancels under normalisation); an informative per-component weight does not.
    """
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    out = _update(raw, gdna_prior=12.0, rna_prior=40.0, has_gdna=False)

    assert out[0] == 0.0, "a gDNA pseudocount was paid with no gDNA component to receive it"
    assert float(out[1:].sum()) == pytest.approx(30.0 + 70.0 + 40.0, rel=1e-12), (
        "the RNA prior was discarded at a locus with no gDNA candidate"
    )
    # shared out by evidence, as everywhere else
    assert out[1] / out[2] == pytest.approx(30.0 / 70.0, rel=1e-12)


def test_a_negative_or_nonfinite_raw_count_is_read_as_zero():
    """``nonnegative_finite`` guards every read. A NaN must not propagate into the pool total."""
    raw = np.array([100.0, np.nan, -5.0, 70.0], dtype=np.float64)
    out = _update(raw, gdna_prior=0.0, rna_prior=30.0)
    assert np.all(np.isfinite(out))
    _g, r_out = _split(out, 0)
    assert r_out == pytest.approx(70.0 + 30.0, rel=1e-12)


def test_the_binding_refuses_a_mismatched_auxiliary_length():
    """An auxiliary array of the wrong length would read past the end of the components. Empty means "this
    locus has none"; any other wrong length is an error."""
    raw = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    with pytest.raises(Exception):
        _update(raw, carried=[1.0, 2.0])
    with pytest.raises(Exception):
        _update(raw, weight=[1.0, 2.0])


# ──────────────────────────────────────────────────────────────────────────────
# The weighted allocation: the per-transcript lane
# ──────────────────────────────────────────────────────────────────────────────


def test_a_weight_proportional_to_evidence_reproduces_the_shipped_rule():
    """The default update

        out[i] = raw[i] · (1 + rna_prior/annotated_count)

    equals ``raw[i] + rna_prior·w_i/Σw`` at ``w_i = raw[i]``, so passing the raw counts as weights reproduces
    it, to floating tolerance since the operation order differs. The weighted lane generalises the default
    rule; a new allocation changes only the weights.
    """
    rng = np.random.default_rng(7)
    raw = rng.uniform(1.0, 1000.0, size=6)
    shipped = _update(raw, gdna_prior=5.0, rna_prior=250.0)
    weighted = _update(raw, weight=raw, gdna_prior=5.0, rna_prior=250.0)
    np.testing.assert_allclose(weighted, shipped, rtol=1e-12, atol=1e-9)


def test_the_weighted_branch_obeys_the_same_conservation_identity():
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
        # A weight vector nominating nobody falls back to the evidence-proportional rule, so the pool still
        # receives the prior whenever it has evidence.
        expect = r_raw + (77.0 if (w[1:].sum() > 0.0 or r_raw > EM_LOG_EPSILON) else 0.0)
        assert r_out == pytest.approx(expect, rel=1e-11, abs=1e-9)


def test_a_zero_count_component_can_be_revived_by_a_weighted_prior():
    """Where the two allocations differ. Under evidence-proportional weights ``out[i] = 0`` is absorbing; a
    strictly positive weight has no such state. That is also how a weighted rule could revive a synthetic span
    the data do not support, so any weight rule must justify it.
    """
    raw = np.array([100.0, 0.0, 70.0], dtype=np.float64)
    assert _update(raw, rna_prior=500.0)[1] == 0.0, "the shipped rule should not revive it"
    revived = _update(raw, weight=[0.0, 1.0, 1.0], rna_prior=500.0)
    assert revived[1] == pytest.approx(250.0, rel=1e-12)


def test_a_weight_of_zero_is_how_a_caller_withholds_the_prior_from_one_component():
    """The weight vector is the only place a caller can exclude a component from the prior, by giving it a
    weight of zero; the exclusion is visible in the array the caller passed."""
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    out = _update(raw, weight=[0.0, 1.0, 0.0], rna_prior=50.0)
    assert out[2] == pytest.approx(70.0, rel=1e-12), "a zero-weight component took prior mass"
    assert out[1] == pytest.approx(30.0 + 50.0, rel=1e-12), "the whole prior did not land"


def test_an_all_zero_weight_vector_falls_back_rather_than_dropping_the_prior():
    """A weight vector that nominates nobody falls back to the default rule. Dropping the RNA prior while gDNA
    keeps its pseudocount would move the split."""
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    fallback = _update(raw, weight=np.zeros(3), gdna_prior=4.0, rna_prior=50.0)
    shipped = _update(raw, gdna_prior=4.0, rna_prior=50.0)
    np.testing.assert_array_equal(fallback, shipped)


def test_a_negative_or_nonfinite_weight_is_read_as_zero():
    """Weights pass through the same ``nonnegative_finite`` guard as counts, so a NaN cannot reach the
    normaliser."""
    raw = np.array([100.0, 30.0, 70.0, 10.0], dtype=np.float64)
    out = _update(raw, weight=[0.0, np.nan, -3.0, 2.0], rna_prior=50.0)
    assert np.all(np.isfinite(out))
    assert out[1] == pytest.approx(30.0, rel=1e-12)
    assert out[2] == pytest.approx(70.0, rel=1e-12)
    assert out[3] == pytest.approx(10.0 + 50.0, rel=1e-12)


def test_the_weighted_branch_needs_no_carried_state_to_place_a_prior_only_pool():
    """The carried state only decides who receives the prior when there is no evidence. An explicit weight
    decides that directly, so the weighted branch needs no zero-evidence special case."""
    raw = np.zeros(4, dtype=np.float64)
    raw[0] = 100.0  # gDNA only
    out = _update(raw, weight=[0.0, 3.0, 1.0, 0.0], gdna_prior=2.0, rna_prior=40.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(102.0, rel=1e-12)
    assert r_out == pytest.approx(40.0, rel=1e-12)
    assert out[1] == pytest.approx(30.0, rel=1e-12)
    assert out[2] == pytest.approx(10.0, rel=1e-12)
    assert out[3] == 0.0
