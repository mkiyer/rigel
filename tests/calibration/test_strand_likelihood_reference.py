"""⭐⭐⭐ THE GATE THE SOURCE CLAIMED TO HAVE AND DID NOT — and it is what makes a "dead" module a REFERENCE.

`simplex.py`'s docstring asserted, in the source's own words, that its three-component strand likelihood
*"collapses to ``strand_loglik`` exactly (the no-regression guard)"*. **No such guard existed anywhere.**
Measured 2026-08-07 by `scripts/design/module_census.py`:

* `strand_likelihood.strand_loglik` had **no production consumer at all**. Its own docstring said "Used by
  the per-region strand module (`strand_deconv`)", and `strand_deconv` does not import it.
* the only importers were two test files, and both exercised it against **its own properties** — never
  against the three-component form it is supposedly the special case of.

⛔ So the two-component function was a **second home for a predicate** whose one home is ψ's
`_mixture_strand_loglik`, with nothing detecting drift between them. That is the shape of TRAPS' "a test that
re-derives a definition cannot detect drift in it", one layer up: the *definition* was duplicated and the
*claim* that the two agree was prose.

⭐ **Deleting it was the other option and this is the better one.** With this file the two-component form
becomes what the source already called it — an executable REFERENCE, gated, in the pattern
`tests/native/_accumulator_reference.py` established. If they ever diverge, this fires instead of a comment
being quietly wrong. And a reference is worth having precisely because ψ's version is 3-component,
vectorised over a lattice and hard to read: the special case is the one a human can check by hand.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.simplex_logodds import _mixture_strand_loglik
from rigel.calibration.strand_likelihood import strand_loglik

#: the candidate gDNA fractions both forms are evaluated on
GRID = np.linspace(0.02, 0.98, 25)


def _three(u_pos, n, f_g, kappa, od_g, od_r, *, tilt=1.0):
    """ψ's three-component form with ALL the RNA on the ``+`` strand — the single-strand special case.

    ⭐ ``tilt = 1`` means ``f_neg = 0``, which is the structural state of a single-strand region: one strand is
    not admissible, so the tilt is not a free nuisance and the mixture reduces to two components.
    ⚠ The variance REFERENCE is passed equal to the live composition, because `strand_loglik` has no
    count-zero-information freeze — it evaluates the variance at the same composition as the mean. Passing a
    different reference would be comparing two different estimators, not two forms of one.
    """
    f_rna = 1.0 - f_g
    f_pos, f_neg = f_rna * tilt, f_rna * (1.0 - tilt)
    return _mixture_strand_loglik(u_pos, n, f_g, f_pos, f_neg, kappa, od_g, od_r, f_g, f_pos, f_neg)


@pytest.mark.parametrize("sense,antisense", [(30.0, 10.0), (126.0, 26.0), (5.0, 4.0), (0.0, 12.0)])
@pytest.mark.parametrize("kappa", [0.5, 0.75, 0.95, 0.99])
@pytest.mark.parametrize("od_g,od_r", [(0.0, 0.0), (0.1, 0.0), (0.0, 0.1), (0.1, 0.1)])
def test_the_three_component_form_collapses_onto_the_two_component_reference(
    sense, antisense, kappa, od_g, od_r
):
    """⭐⭐ **THE GUARD, as the source described it.** With one RNA strand dead the mixture has two
    components, so ψ's form must equal the reference — not approximately, to floating-point identity of the
    same algebra.

    ⚠ Both are log-likelihoods up to an additive constant in the DATA (neither normalises the binomial
    coefficient), so they are compared **after removing a per-call constant offset**: what has to match is
    the SHAPE over the gDNA grid, which is the whole information content. A constant offset cannot move a
    posterior; a shape difference can.
    """
    n = sense + antisense
    ref = strand_loglik(
        GRID,
        sense,
        antisense,
        kappa,
        gdna_strand_overdispersion=od_g,
        rna_strand_overdispersion=od_r,
    )
    got = _three(sense, n, GRID, kappa, od_g, od_r)
    ref, got = np.asarray(ref, np.float64), np.asarray(got, np.float64)
    assert ref.shape == got.shape == GRID.shape
    np.testing.assert_allclose(got - got.mean(), ref - ref.mean(), rtol=0, atol=1e-9)


def test_PERTURBATION_a_wrong_mixture_rate_BREAKS_the_collapse():
    """⛔ TRAPS: perturb-every-gate's second half: the gate above is worth nothing until it is shown to fire.

    The mixture plus-strand rate is ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``. Perturb the gDNA arm's rate from ½ —
    the one number that makes gDNA *unstranded* — and the collapse must fail."""
    sense, antisense, kappa = 30.0, 10.0, 0.9
    n = sense + antisense
    f_rna = 1.0 - GRID

    def bad(u_pos, nn, f_g, f_pos, f_neg, kap, od_g, od_r, *_ref):
        p = 0.55 * f_g + kap * f_pos + (1.0 - kap) * f_neg  # ⛔ 0.55, not ½
        var = np.maximum(nn * p * (1.0 - p), 1e-9)
        return -0.5 * (u_pos - nn * p) ** 2 / var - 0.5 * np.log(var)

    ref = strand_loglik(GRID, sense, antisense, kappa, gdna_strand_overdispersion=0.0)
    got = bad(sense, n, GRID, f_rna, 0.0 * f_rna, kappa, 0.0, 0.0)
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(got - got.mean(), ref - ref.mean(), rtol=0, atol=1e-9)


def test_PERTURBATION_giving_the_dead_strand_mass_BREAKS_the_collapse():
    """The collapse is a claim about the SINGLE-STRAND case specifically. Split the RNA across both strands
    and the two forms must part company — otherwise the test would pass for a three-component input and
    would not be testing the special case at all (TRAPS: could-the-arm-have-fired: check the gate could have failed)."""
    sense, antisense, kappa = 30.0, 10.0, 0.9
    n = sense + antisense
    ref = strand_loglik(GRID, sense, antisense, kappa, gdna_strand_overdispersion=0.0)
    both = _three(sense, n, GRID, kappa, 0.0, 0.0, tilt=0.5)
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(both - both.mean(), ref - ref.mean(), rtol=0, atol=1e-9)


def test_at_kappa_one_half_the_strand_says_nothing_about_composition():
    """⭐ The domain fact both forms must encode: on a genuinely unstranded library the strand channel
    carries **exactly zero** information about the gDNA fraction — the mixture rate is ½ whatever ``f_g``
    is, so the log-likelihood is FLAT over the grid. That is why an unstranded slot has no own composition
    evidence and why the message layer's whole value sits in that stratum."""
    n = 40.0
    got = _three(24.0, n, GRID, 0.5, 0.0, 0.0)
    assert np.ptp(np.asarray(got, np.float64)) < 1e-12
