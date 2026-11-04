"""The two-component strand likelihood is an executable reference, and this module checks the production form against
it.

ψ's strand term is the production form: three components, native (`native/transfer_rows.h`, read through
`_psi_reference.psi_cube`), hard to read by hand. `_psi_reference.strand_loglik` is the two-component special case it
must collapse to, a readable statement of the same predicate kept beside the test oracles and never in `src/`; this
module keeps the two from drifting apart, as `tests/native/_accumulator_reference.py` does for the accumulator.
"""

from __future__ import annotations

import numpy as np
import pytest
from _psi_reference import jeffreys_arms, psi_cube, strand_loglik, strand_loglik_mixture

#: the candidate gDNA fractions both forms are evaluated on
GRID = np.linspace(0.02, 0.98, 25)


def _three(u_pos, n, f_g, kappa):
    """ψ's production strand term with all the RNA on the ``+`` strand, the single-strand special case, read out of the
    solver's own cube: one slot per grid point, each with its variance frozen at that point's live composition
    (`strand_loglik` evaluates the variance at the composition of the mean, so a different reference would compare two
    estimators, not two forms of one), the λ grid the logit of ``f_g``, the arms subtracted, and the diagonal ``psi[j,
    j]`` the term at point ``j``."""
    f_g = np.asarray(f_g, np.float64)
    m = f_g.shape[0]
    lam = np.log(f_g / (1.0 - f_g))
    psi, _fp, _fn, _tau = psi_cube(
        np.full(m, float(u_pos)),
        np.full(m, float(n - u_pos)),
        np.ones(m, bool),
        np.zeros(m, bool),
        f_g,
        1.0 - f_g,
        np.zeros(m),
        kappa=kappa,
        lam=lam,
        ambig=False,
    )
    return psi[np.arange(m), np.arange(m), 0] - jeffreys_arms(lam)


def _three_readable(u_pos, n, f_g, kappa, *, tilt=1.0):
    """The readable three-component term (`_psi_reference`) at an arbitrary tilt, needed because the solver's
    single-strand cube has no tilt to give the dead strand mass."""
    f_rna = 1.0 - f_g
    f_pos, f_neg = f_rna * tilt, f_rna * (1.0 - tilt)
    return strand_loglik_mixture(u_pos, n, f_g, f_pos, f_neg, kappa, f_g, f_pos, f_neg)


@pytest.mark.parametrize("sense,antisense", [(30.0, 10.0), (126.0, 26.0), (5.0, 4.0), (0.0, 12.0)])
@pytest.mark.parametrize("kappa", [0.5, 0.75, 0.95, 0.99])
def test_the_three_component_form_collapses_onto_the_two_component_reference(
    sense, antisense, kappa
):
    """With one RNA strand dead the mixture has two components, so ψ's form must equal the reference to floating-point
    identity of the same algebra. Both are log-likelihoods up to an additive constant in the data (neither normalises
    the binomial coefficient), so they are compared after removing a per-call offset: the shape over the gDNA grid is
    the whole information content, and a constant cannot move a posterior.
    """
    n = sense + antisense
    ref = strand_loglik(GRID, sense, antisense, kappa)
    got = _three(sense, n, GRID, kappa)
    ref, got = np.asarray(ref, np.float64), np.asarray(got, np.float64)
    assert ref.shape == got.shape == GRID.shape
    np.testing.assert_allclose(got - got.mean(), ref - ref.mean(), rtol=0, atol=1e-9)


def test_perturbation_a_wrong_mixture_rate_breaks_the_collapse():
    """The mixture plus-strand rate is ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``. Perturbing the gDNA arm's rate from ½, the
    number that makes gDNA unstranded, must break the collapse."""
    sense, antisense, kappa = 30.0, 10.0, 0.9
    n = sense + antisense
    f_rna = 1.0 - GRID

    def bad(u_pos, nn, f_g, f_pos, f_neg, kap, *_ref):
        p = 0.55 * f_g + kap * f_pos + (1.0 - kap) * f_neg  # 0.55, not ½
        var = np.maximum(nn * p * (1.0 - p), 1e-9)
        return -0.5 * (u_pos - nn * p) ** 2 / var - 0.5 * np.log(var)

    ref = strand_loglik(GRID, sense, antisense, kappa)
    got = bad(sense, n, GRID, f_rna, 0.0 * f_rna, kappa)
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(got - got.mean(), ref - ref.mean(), rtol=0, atol=1e-9)


def test_perturbation_giving_the_dead_strand_mass_breaks_the_collapse():
    """The collapse is specific to the single-strand case: with the RNA split across both strands the two forms must
    part, or the test above would not be testing the special case."""
    sense, antisense, kappa = 30.0, 10.0, 0.9
    n = sense + antisense
    ref = strand_loglik(GRID, sense, antisense, kappa)
    both = _three_readable(sense, n, GRID, kappa, tilt=0.5)
    # the readable term at the single-strand tilt is the production term, so the perturbation is about the same form
    np.testing.assert_allclose(
        _three_readable(sense, n, GRID, kappa) - _three(sense, n, GRID, kappa),
        0.0,
        atol=1e-9,
    )
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(both - both.mean(), ref - ref.mean(), rtol=0, atol=1e-9)


def test_at_kappa_one_half_the_strand_says_nothing_about_composition():
    """On an unstranded library the strand channel carries no information about the gDNA fraction: the mixture rate is ½
    whatever ``f_g`` is, so the log-likelihood is flat over the grid. That is why an unstranded slot has no own
    composition evidence and why the message layer's value lies in that stratum."""
    n = 40.0
    got = _three(24.0, n, GRID, 0.5)
    assert np.ptp(np.asarray(got, np.float64)) < 1e-12
