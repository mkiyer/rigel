"""The two-component strand likelihood is an executable REFERENCE, and this is its gate.

ψ's strand term is the production form: three components, native (`native/transfer_rows.h`, read through
`simplex_logodds.psi_cube`), hard to read by hand. `_psi_reference.strand_loglik` is the two-component special
case it must collapse to — a second, readable statement of one predicate that lives with the gates' oracles and
never in `src/` (one production path), and this file is what stops the two drifting apart
(TRAPS: a-test-that-redefines). The pattern is the one `tests/native/_accumulator_reference.py`
established: keep the readable form beside the gates, and gate the production form against it.
"""

from __future__ import annotations

import numpy as np
import pytest
from _psi_reference import jeffreys_arms, strand_loglik, strand_loglik_mixture

from rigel.calibration.simplex_logodds import psi_cube

#: the candidate gDNA fractions both forms are evaluated on
GRID = np.linspace(0.02, 0.98, 25)


def _three(u_pos, n, f_g, kappa, od_g, od_r):
    """ψ's production strand term with ALL the RNA on the ``+`` strand — the single-strand special case,
    read out of the solver's own cube: one slot per grid point, each with its variance frozen at the live
    composition of that point (`strand_loglik` has no count-zero-information freeze — it evaluates the
    variance at the same composition as the mean, so a different reference would be comparing two
    estimators, not two forms of one), the λ grid the logit of ``f_g``, the arms subtracted, the diagonal
    ``psi[j, j]`` the term at point ``j``."""
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
        od_g=od_g,
        od_r=od_r,
        lam=lam,
        ambig=False,
    )
    return psi[np.arange(m), np.arange(m), 0] - jeffreys_arms(lam)


def _three_readable(u_pos, n, f_g, kappa, od_g, od_r, *, tilt=1.0):
    """The readable three-component term (`_psi_reference`) at an arbitrary tilt — the perturbation's
    handle, since the solver's single-strand cube has no tilt to give the dead strand mass with."""
    f_rna = 1.0 - f_g
    f_pos, f_neg = f_rna * tilt, f_rna * (1.0 - tilt)
    return strand_loglik_mixture(u_pos, n, f_g, f_pos, f_neg, kappa, od_g, od_r, f_g, f_pos, f_neg)


@pytest.mark.parametrize("sense,antisense", [(30.0, 10.0), (126.0, 26.0), (5.0, 4.0), (0.0, 12.0)])
@pytest.mark.parametrize("kappa", [0.5, 0.75, 0.95, 0.99])
@pytest.mark.parametrize("od_g,od_r", [(0.0, 0.0), (0.1, 0.0), (0.0, 0.1), (0.1, 0.1)])
def test_the_three_component_form_collapses_onto_the_two_component_reference(
    sense, antisense, kappa, od_g, od_r
):
    """With one RNA strand dead the mixture has two components, so ψ's form must equal the reference — not
    approximately, to floating-point identity of the same algebra.

    Both are log-likelihoods up to an additive constant in the data (neither normalises the binomial
    coefficient), so they are compared after removing a per-call constant offset: what has to match is the
    shape over the gDNA grid, which is the whole information content. A constant offset cannot move a
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
    """The gate above is worth nothing until it is shown to fire (TRAPS: perturb-every-gate).

    The mixture plus-strand rate is ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``. Perturb the gDNA arm's rate from ½ —
    the one number that makes gDNA *unstranded* — and the collapse must fail."""
    sense, antisense, kappa = 30.0, 10.0, 0.9
    n = sense + antisense
    f_rna = 1.0 - GRID

    def bad(u_pos, nn, f_g, f_pos, f_neg, kap, od_g, od_r, *_ref):
        p = 0.55 * f_g + kap * f_pos + (1.0 - kap) * f_neg  # 0.55, not ½
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
    both = _three_readable(sense, n, GRID, kappa, 0.0, 0.0, tilt=0.5)
    # and the readable term at the single-strand tilt IS the production term, so the perturbation is
    # a statement about the same form
    np.testing.assert_allclose(
        _three_readable(sense, n, GRID, kappa, 0.0, 0.0) - _three(sense, n, GRID, kappa, 0.0, 0.0),
        0.0,
        atol=1e-9,
    )
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(both - both.mean(), ref - ref.mean(), rtol=0, atol=1e-9)


def test_at_kappa_one_half_the_strand_says_nothing_about_composition():
    """The domain fact both forms must encode: on a genuinely unstranded library the strand channel
    carries exactly zero information about the gDNA fraction — the mixture rate is ½ whatever ``f_g``
    is, so the log-likelihood is FLAT over the grid. That is why an unstranded slot has no own composition
    evidence and why the message layer's whole value sits in that stratum."""
    n = 40.0
    got = _three(24.0, n, GRID, 0.5, 0.0, 0.0)
    assert np.ptp(np.asarray(got, np.float64)) < 1e-12
