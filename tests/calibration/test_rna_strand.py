"""Tests for the RNA strand Beta-Binomial overdispersion (docs/em_strand/05).

Two properties:

* **Estimator** — RNA overdispersion drawn at a known level (per splice junction) is recovered by
  :func:`fit_rna_strand_overdispersion` (its component mean is κ, not ½, so the excess-variance
  scale must be κ(1−κ), not ¼).
* **Deconv symmetry** — with the RNA strand modelled as Beta-Binomial too (same default prior as
  gDNA), an *unstranded* region (κ = ½) is **uninformative**, and a balanced region grows more
  gDNA-like as the library becomes more stranded. The earlier gDNA-only overdispersion broke this
  (it pulled balanced unstranded regions toward RNA).
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.gdna_strand import (
    _MAX_OVERDISPERSION,
    RnaStrandModel,
    fit_rna_strand_from_sj_table,
    fit_rna_strand_overdispersion,
    overdispersion_for_beta,
)
from rigel.calibration.strand_likelihood import strand_loglik


def _beta_binom_regions(rng, n_regions, depth, overdispersion, mean):
    """Per-region (sense, total) ~ BetaBinom(depth, mean, overdispersion); symmetric Beta on `mean`."""
    conc = (1.0 - overdispersion) / overdispersion  # a + b
    p = rng.beta(mean * conc, (1.0 - mean) * conc, size=n_regions)
    total = np.full(n_regions, depth, dtype=np.float64)
    sense = rng.binomial(depth, p).astype(np.float64)
    return sense, total


def _decoded_gdna_frac(sense, antisense, kappa, *, gdna_od, rna_od, n_grid=4000):
    """Posterior median gDNA fraction of one region under a FLAT count prior (strand-only deconv).

    Mirrors the strand module (``strand_deconv._deconv_per_region`` strand branch): a weak prior ×
    the strand likelihood, isolating the strand likelihood's effect on the deconvolution.
    """
    grid = np.linspace(1e-6, 1.0 - 1e-6, n_grid)
    ll = strand_loglik(
        grid,
        sense,
        antisense,
        kappa,
        gdna_strand_overdispersion=gdna_od,
        rna_strand_overdispersion=rna_od,
    )
    w = np.exp(ll - ll.max())
    w /= w.sum()
    return float(np.interp(0.5, np.cumsum(w), grid))


# --------------------------------------------------------------------------- estimator


@pytest.mark.parametrize("true_od", [0.02, 0.05, 0.10, 0.20])
@pytest.mark.parametrize("kappa", [0.7, 0.9, 0.99])
def test_recovers_rna_overdispersion(true_od, kappa):
    """RNA overdispersion at a known level is recovered for a stranded library (mean κ ≠ ½).

    Guards the κ(1−κ) excess-variance scale: a ¼ scale would bias the estimate by 4κ(1−κ).
    """
    rng = np.random.default_rng(2024)
    sense, total = _beta_binom_regions(
        rng, n_regions=6000, depth=150, overdispersion=true_od, mean=kappa
    )
    # Weak prior so the 6000 seed sides dominate (tests the estimator, not the shrinkage).
    model = fit_rna_strand_overdispersion(
        sense, total, kappa, prior_overdispersion=1 / 7, prior_weight=5.0
    )
    assert not model.fallback_used
    assert model.n_seed_regions == 6000
    assert model.rna_strand_overdispersion == pytest.approx(true_od, rel=0.20, abs=0.01)


def test_binomial_limit_recovers_near_zero():
    """A clean stranded library (no overdispersion, shared rate exactly κ) → od ≈ 0, not spurious."""
    rng = np.random.default_rng(11)
    kappa = 0.95
    total = np.full(4000, 150.0)
    sense = rng.binomial(150, kappa, size=4000).astype(np.float64)
    model = fit_rna_strand_overdispersion(
        sense, total, kappa, prior_overdispersion=0.0, prior_weight=0.0
    )
    assert model.rna_strand_overdispersion < 0.01


def test_fit_clamped_to_ceiling():
    """Extreme overdispersion is clamped to the Beta(2,2) ceiling (od ≤ 0.2)."""
    rng = np.random.default_rng(3)
    sense, total = _beta_binom_regions(rng, n_regions=4000, depth=150, overdispersion=0.45, mean=0.85)
    model = fit_rna_strand_overdispersion(
        sense, total, 0.85, prior_overdispersion=1 / 7, prior_weight=5.0
    )
    assert model.rna_strand_overdispersion <= _MAX_OVERDISPERSION + 1e-12


def test_no_spliced_data_falls_back_to_prior():
    """No spliced fragments anywhere → fallback to the prior overdispersion (not 0)."""
    total = np.zeros(50)
    sense = np.zeros(50)
    prior = overdispersion_for_beta(3.0)  # 1/7
    model = fit_rna_strand_overdispersion(
        sense, total, 0.9, prior_overdispersion=prior, prior_weight=30.0
    )
    assert model.fallback_used
    assert model.rna_strand_overdispersion == pytest.approx(prior)
    assert model.n_seed_regions == 0


def test_shrinks_to_prior_when_sparse():
    """One thin seed side ⇒ the prior dominates (continuous shrinkage, no hard gate).

    ⚠ The prior weight is now in INFORMATION units (see `gdna_strand._prior_information`), so it must be
    given on that scale; a single 10-fragment side carries only ~4 information units at κ = 0.9."""
    sense = np.array([7.0])
    total = np.array([10.0])
    prior = 1 / 7
    model = fit_rna_strand_overdispersion(
        sense, total, 0.9, prior_overdispersion=prior, prior_weight=909.0
    )
    assert model.rna_strand_overdispersion == pytest.approx(prior, abs=0.03)


# ---------------------------------------------------------------- fit-from-SJ-table wrapper


def _sj_table(sense, antisense):
    """Minimal duck-typed SJStrandTable exposing the two count arrays the fit reads."""
    return SimpleNamespace(
        n_sense=np.asarray(sense, dtype=np.int64),
        n_antisense=np.asarray(antisense, dtype=np.int64),
    )


def test_fit_from_sj_table_uses_every_junction():
    """The wrapper fits one seed per junction; junctions with no fragments drop out."""
    rng = np.random.default_rng(5)
    kappa = 0.9
    sense, total = _beta_binom_regions(rng, 6000, 120, 0.10, kappa)
    model = fit_rna_strand_from_sj_table(
        _sj_table(sense, total - sense),
        rna_sense_frac=kappa,
        prior_overdispersion=1 / 7,
        prior_weight=5.0,
    )
    assert model.n_seed_regions == 6000
    assert model.rna_strand_overdispersion == pytest.approx(0.10, rel=0.25, abs=0.01)


def test_fit_from_sj_table_empty_is_fallback():
    """No junctions carrying fragments → fallback to prior (a library with no spliced reads)."""
    model = fit_rna_strand_from_sj_table(
        _sj_table([0, 0], [0, 0]),
        rna_sense_frac=0.9,
        prior_overdispersion=1 / 7,
        prior_weight=30.0,
    )
    assert model.fallback_used


def test_fit_from_sj_table_matches_the_real_table_type():
    """The production ``SJStrandTable`` satisfies the wrapper's contract (no duck typing)."""
    from rigel.strand_model import SJStrandTable

    assert fit_rna_strand_from_sj_table(SJStrandTable.empty(), rna_sense_frac=0.9).fallback_used


# --------------------------------------------------------------------------- deconv symmetry


def test_unstranded_is_uninformative_with_symmetric_overdispersion():
    """κ = ½, equal gDNA/RNA overdispersion ⇒ a balanced region deconvolves to gdna_frac ≈ ½ (flat)."""
    od = overdispersion_for_beta(3.0)
    frac = _decoded_gdna_frac(50, 50, 0.5, gdna_od=od, rna_od=od)
    assert frac == pytest.approx(0.5, abs=0.02)


def test_asymmetric_overdispersion_biases_unstranded_toward_rna():
    """The OLD asymmetry (gDNA Beta-Binomial, RNA Binomial) spuriously pulls a balanced unstranded
    region toward RNA — the bug this change fixes. Symmetric overdispersion removes the pull."""
    od = overdispersion_for_beta(3.0)
    asym = _decoded_gdna_frac(50, 50, 0.5, gdna_od=od, rna_od=0.0)
    symm = _decoded_gdna_frac(50, 50, 0.5, gdna_od=od, rna_od=od)
    assert asym < 0.4  # materially pulled toward RNA
    assert symm == pytest.approx(0.5, abs=0.02)  # pull removed


def test_graded_information_balanced_region_more_gdna_as_library_stranded():
    """A balanced (50/50) region reads more gDNA-like as κ rises from ½ (unstranded) toward 1.

    At κ = ½ it is uninformative (½); as the library becomes more stranded, a *symmetric* split
    looks increasingly like the symmetric gDNA component, so gdna_frac rises monotonically — the
    'unstranded → weakly → strongly stranded' information gradient.
    """
    od = overdispersion_for_beta(3.0)
    kappas = [0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
    fracs = [_decoded_gdna_frac(50, 50, k, gdna_od=od, rna_od=od) for k in kappas]
    assert fracs[0] == pytest.approx(0.5, abs=0.02)
    assert all(b >= a - 1e-6 for a, b in zip(fracs, fracs[1:]))  # monotone non-decreasing
    assert fracs[-1] > fracs[0] + 0.2  # materially more gDNA at high κ


def test_rna_overdispersion_zero_recovers_binomial_decode():
    """rna_strand_overdispersion = 0 ⇒ strand_loglik is byte-identical to the gDNA-only formula."""
    grid = np.linspace(1e-6, 1 - 1e-6, 200)
    base = strand_loglik(grid, 30, 10, 0.9, gdna_strand_overdispersion=0.1)
    with_rna0 = strand_loglik(
        grid, 30, 10, 0.9, gdna_strand_overdispersion=0.1, rna_strand_overdispersion=0.0
    )
    np.testing.assert_array_equal(base, with_rna0)


def test_rna_model_beta_concentration_roundtrip():
    """RnaStrandModel.beta_concentration inverts overdispersion_for_beta."""
    m = RnaStrandModel(
        rna_strand_overdispersion=overdispersion_for_beta(3.0),
        n_seed_regions=1,
        n_seed_fragments=1,
        fallback_used=False,
    )
    assert m.beta_concentration() == pytest.approx(3.0)
    m0 = RnaStrandModel(
        rna_strand_overdispersion=0.0, n_seed_regions=0, n_seed_fragments=0, fallback_used=True
    )
    assert m0.beta_concentration() == float("inf")
