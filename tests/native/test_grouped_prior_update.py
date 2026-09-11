"""``apply_grouped_prior_update`` — the identity the per-component prior rests on.

For a GIVEN ``raw_counts`` vector the function must leave the gDNA:RNA split exactly where the two
pseudocounts put it::

    out[gdna]           == raw[gdna] + gdna_prior
    Σ_{i ≠ gdna} out[i] == Σ_{i ≠ gdna} raw[i] + rna_prior

so that withholding the RNA prior from SYNTHETIC components redistributes mass strictly WITHIN the RNA
pool and never moves the number calibration exists to produce. It is a PER-CALL identity and not an
end-to-end one, and conflating the two has produced a wrong test more than once: the EM iterates around
this function, so a larger ``rna_prior`` shifts ``theta``, hence the E-step, hence the next iteration's
``raw_counts``, hence the converged gDNA total. "The library gDNA fraction cannot move" is false by
design. Everything below fixes ``raw_counts`` and calls the function ONCE, which is the only regime the
identity claims. The function is ``static`` in the C++, so the identity was untestable until
``_apply_grouped_prior_update_test``, the test-only binding these gates reach it through.

PERTURBATION: twelve edits to the C++, each rebuilt and re-run and each first scored against a random
configuration battery so that "the tests did not catch it" is separated from "it could not be caught"
(TRAPS: could-the-arm-have-fired). Eleven were caught. The twelfth, dropping the
no-eligible-recipient gate, is provably INERT and that is a finding about the code rather than a hole
here: ``prior_recipients`` cannot change ``out_counts``, because whenever it would zero ``rna_prior``
the eligible total it divides by is itself zero and the ``inv = 0`` guard downstream has already
dropped the prior. It is defensive redundancy, kept because it states the intent where the intent is
decided — do not "cover" it by weakening a gate, since no input distinguishes the two. The battery
needed falsifying too: its first version generated only well-formed weights, so a perturbation removing
the negative-weight guard moved nothing and read as inert — a verdict about the generator wearing the
costume of a verdict about the code.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel._em_impl import _apply_grouped_prior_update_test

_EMPTY_F64 = np.zeros(0, dtype=np.float64)
_EMPTY_U8 = np.zeros(0, dtype=np.uint8)


def _update(
    raw,
    *,
    carried=None,
    is_synthetic=None,
    weight=None,
    gdna_prior=0.0,
    rna_prior=0.0,
    gdna_index=0,
    has_gdna=True,
):
    return _apply_grouped_prior_update_test(
        raw_counts=np.ascontiguousarray(raw, dtype=np.float64),
        carried_state=_EMPTY_F64 if carried is None else np.ascontiguousarray(carried, np.float64),
        is_synthetic=(
            _EMPTY_U8 if is_synthetic is None else np.ascontiguousarray(is_synthetic, np.uint8)
        ),
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
    "raw, is_synthetic, eligible",
    [
        ([100.0, 30.0, 70.0], None, True),
        ([100.0, 30.0, 70.0], [0, 0, 1], True),
        # NO eligible recipient — every RNA component is synthetic. Written here, in the identity's
        # own table, because the identity is CONDITIONAL and a version of this test asserting the
        # unconditional form fails on exactly this row. The pool gets +0, and gDNA still gets its own.
        ([100.0, 30.0, 70.0], [0, 1, 1], False),
        ([0.0, 30.0, 70.0], [0, 0, 1], True),
        ([100.0, 1e-9, 70.0], [0, 1, 0], True),
        ([5.0, 1234.5, 0.0, 17.25, 900.0], [0, 0, 1, 0, 1], True),
    ],
)
def test_the_gDNA_RNA_split_lands_exactly_where_the_two_pseudocounts_put_it(
    raw, is_synthetic, eligible
):
    """The whole point: gDNA gets ``+gdna_prior`` and the RNA POOL gets ``+rna_prior`` — whatever the
    synthetic mask does INSIDE the pool.

    Conditional on there being an eligible recipient. With every RNA component synthetic the prior has
    nowhere to land, so the pool gets ``+0`` — and gDNA must NOT keep its own pseudocount in a way that
    moves the split. That asymmetry is what ``em_solver.cpp``'s ``prior_recipients`` gate exists for.
    """
    raw = np.asarray(raw, dtype=np.float64)
    gdna_prior, rna_prior = 12.5, 40.0
    out = _update(raw, is_synthetic=is_synthetic, gdna_prior=gdna_prior, rna_prior=rna_prior)

    g_out, r_out = _split(out, 0)
    g_raw, r_raw = _split(raw, 0)
    assert g_out == pytest.approx(g_raw + gdna_prior, rel=1e-12, abs=1e-12)
    assert r_out == pytest.approx(r_raw + (rna_prior if eligible else 0.0), rel=1e-12, abs=1e-12)


def test_a_SYNTHETIC_component_receives_none_of_the_rna_prior():
    """A synthetic component is a shadow span the INDEX manufactured; nothing asserts it exists, so
    the null is that it is ABSENT. It keeps its own evidence UNSCALED and takes no pseudocount."""
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    out = _update(raw, is_synthetic=[0, 0, 1], gdna_prior=0.0, rna_prior=50.0)
    assert out[2] == pytest.approx(raw[2], rel=1e-12), "the synthetic component was scaled"
    # ...and the whole prior landed on the one eligible recipient.
    assert out[1] == pytest.approx(raw[1] + 50.0, rel=1e-12)


def test_a_zero_count_component_CANNOT_be_revived_by_the_prior():
    """``out[i]`` is PROPORTIONAL to ``raw[i]``, so ``raw[i] == 0`` stays 0. That is why the
    coverage-weighted warm start is the only spark a synthetic entity ever gets, and why ``theta = 0``
    is a fixed point under ``alpha = 0``."""
    raw = np.array([100.0, 0.0, 70.0], dtype=np.float64)
    out = _update(raw, gdna_prior=0.0, rna_prior=500.0)
    assert out[1] == 0.0


def test_an_ALL_SYNTHETIC_locus_takes_NO_rna_prior_and_the_split_still_holds():
    """The degenerate locus of the derivation: no eligible recipient ⇒ ``rna_prior`` is dropped ⇒ the
    RNA pool must outcompete gDNA unaided. The gDNA side must NOT keep its prior while the RNA side
    silently loses one — that asymmetry is what the gate in ``em_solver.cpp`` exists for."""
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    out = _update(raw, is_synthetic=[0, 1, 1], gdna_prior=9.0, rna_prior=500.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(100.0 + 9.0, rel=1e-12)
    assert r_out == pytest.approx(30.0 + 70.0, rel=1e-12), (
        "an ineligible prior reached the RNA pool"
    )


def test_the_carried_state_path_is_live_and_obeys_the_SAME_split():
    """The zero-evidence branch, reachable under VBEM (the shipped default), which passes ``alpha`` as
    the carried state. Its gate must name the denominator ITS branch divides by: testing
    ``annotated_count && annotated_carried`` while dividing by ``annotated_count`` alone lets a locus
    keep a live ``rna_prior``, multiply it by ``inv = 0`` and drop it — RNA summing to ``rna_count``
    while gDNA keeps ``gdna_prior``."""
    raw = np.array([100.0, 0.0, 0.0], dtype=np.float64)  # no RNA evidence at all
    carried = np.array([1.0, 2.0, 6.0], dtype=np.float64)
    out = _update(raw, carried=carried, gdna_prior=3.0, rna_prior=40.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(100.0 + 3.0, rel=1e-12)
    assert r_out == pytest.approx(0.0 + 40.0, rel=1e-12)
    # shared out in proportion to the CARRIED state, not the (zero) counts
    assert out[1] / out[2] == pytest.approx(2.0 / 6.0, rel=1e-12)


def test_the_carried_state_branch_ALSO_withholds_from_a_synthetic_component():
    """The case a perturbation found was untested. The other carried-state test passes no synthetic
    mask at all, and the all-synthetic one has an empty eligible set — so ``inv = 0`` makes the right
    and wrong answers coincide. Neither can see the branch PAY a synthetic component.

    This is the configuration that distinguishes them: zero RNA evidence, a MIX of annotated and
    synthetic carried alpha. The pool is prior-only, so a synthetic component has earned nothing and
    must receive nothing. Live under VBEM, which is the shipped default and passes ``alpha`` as the
    carried state.
    """
    raw = np.array([100.0, 0.0, 0.0, 0.0], dtype=np.float64)
    carried = np.array([1.0, 3.0, 5.0, 2.0], dtype=np.float64)
    out = _update(raw, carried=carried, is_synthetic=[0, 0, 1, 0], gdna_prior=0.0, rna_prior=40.0)

    assert out[2] == 0.0, "a synthetic component was paid out of the prior-only pool"
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(100.0, rel=1e-12)
    assert r_out == pytest.approx(40.0, rel=1e-12), "the pool did not sum to rna_count + rna_prior"
    # shared between the two ELIGIBLE components in proportion to their carried alpha (3 : 2)
    assert out[1] / out[3] == pytest.approx(3.0 / 2.0, rel=1e-12)


def test_a_locus_with_NO_annotated_carried_alpha_drops_the_prior_from_BOTH_sides():
    """The configuration the C++ gate exists for: zero annotated COUNT AND zero annotated CARRIED
    alpha. The RNA prior must be dropped, and gDNA must not keep its own."""
    raw = np.array([100.0, 0.0, 0.0], dtype=np.float64)
    carried = np.array([1.0, 0.0, 5.0], dtype=np.float64)  # all carried alpha is SYNTHETIC
    out = _update(raw, carried=carried, is_synthetic=[0, 1, 1], gdna_prior=3.0, rna_prior=40.0)
    g_out, r_out = _split(out, 0)
    assert g_out == pytest.approx(103.0, rel=1e-12)
    assert r_out == pytest.approx(0.0, abs=1e-12), "a prior with no eligible recipient was paid out"


# ──────────────────────────────────────────────────────────────────────────────
# Byte-identity: the no-synthetic path must be UNCHANGED by the synthetic branch
# ──────────────────────────────────────────────────────────────────────────────


def test_an_ALL_ANNOTATED_locus_is_BYTE_IDENTICAL_with_and_without_the_mask():
    """The inert-mask control (TRAPS: byte-identity-gate). An all-zero ``is_synthetic`` must reproduce
    the no-mask arm BIT FOR BIT — the C++ is written in the shipped operation order precisely so that
    it does. ``approx`` would not test this; the comparison is exact."""
    rng = np.random.default_rng(0)
    raw = rng.uniform(0.0, 1000.0, size=7)
    a = _update(raw, gdna_prior=11.0, rna_prior=37.0)
    b = _update(raw, is_synthetic=np.zeros(7, np.uint8), gdna_prior=11.0, rna_prior=37.0)
    np.testing.assert_array_equal(a, b)


def test_the_mask_COULD_have_fired():
    """TRAPS: could-the-arm-have-fired. The byte-identity test above is only a control if the mask is
    capable of changing the answer at all — so prove it does, on the same input."""
    rng = np.random.default_rng(0)
    raw = rng.uniform(0.0, 1000.0, size=7)
    inert = _update(raw, is_synthetic=np.zeros(7, np.uint8), gdna_prior=11.0, rna_prior=37.0)
    live = _update(raw, is_synthetic=[0, 0, 1, 0, 0, 1, 0], gdna_prior=11.0, rna_prior=37.0)
    assert not np.array_equal(inert, live), "the synthetic mask cannot move this input"


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
        _update(raw, is_synthetic=[0, 1])
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
        syn = (rng.random(n) < 0.3).astype(np.uint8)
        syn[0] = 0
        out = _update(raw, weight=w, is_synthetic=syn, gdna_prior=3.0, rna_prior=77.0)
        g_out, r_out = _split(out, 0)
        g_raw, r_raw = _split(raw, 0)
        assert g_out == pytest.approx(g_raw + 3.0, rel=1e-11)
        eligible = np.flatnonzero((np.arange(n) != 0) & (syn == 0))
        expect = r_raw + (77.0 if w[eligible].sum() > 0 else 0.0)
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


def test_a_SYNTHETIC_component_is_STILL_ineligible_under_the_weighted_branch():
    """The weighted lane changes HOW the prior is shared out, not WHO may receive it. A synthetic
    component's weight is excluded from the normalisation and it keeps its evidence unscaled, so a
    manufactured shadow span stays absent until the data proves otherwise."""
    raw = np.array([100.0, 30.0, 70.0], dtype=np.float64)
    out = _update(raw, weight=[0.0, 1.0, 9.0], is_synthetic=[0, 0, 1], rna_prior=50.0)
    assert out[2] == pytest.approx(70.0, rel=1e-12), (
        "a synthetic component took weighted prior mass"
    )
    assert out[1] == pytest.approx(30.0 + 50.0, rel=1e-12), "its weight was not excluded"


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
