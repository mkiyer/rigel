"""E-step soft allocation: count + strand log-Bayes-factors, D1, D7."""

from __future__ import annotations

import numpy as np

from rigel.calibration.estep import Allocation, estep_view
from rigel.calibration.signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from rigel.calibration.substrate import SubstrateView


def _view(pos, neg, sense=None, anti=None, mass_u=None, mass_s=None) -> SubstrateView:
    pos = np.asarray(pos, dtype=np.int64)
    neg = np.asarray(neg, dtype=np.int64)
    r = pos.shape[0]
    sense = np.zeros(r, dtype=np.int64) if sense is None else np.asarray(sense, dtype=np.int64)
    anti = np.zeros(r, dtype=np.int64) if anti is None else np.asarray(anti, dtype=np.int64)
    # Contained-style default: mass == count.
    mu = (pos + neg).astype(np.float64) if mass_u is None else np.asarray(mass_u, dtype=np.float64)
    ms = (
        (sense + anti).astype(np.float64)
        if mass_s is None
        else np.asarray(mass_s, dtype=np.float64)
    )
    return SubstrateView(pos, neg, sense, anti, mu, ms)


def _params(r, **overrides):
    base = dict(
        omega=np.ones(r),
        rho_0=1.0,
        L_eff=np.full(r, 100.0),
        exposure_dispersion=0.1,
        kappa_rna=0.95,
        rho_r_bb=0.05,
        rho_d_bb=0.01,
        eps_s=1.0e-3,
        pi_g_prior=np.full(r, 0.5),
        m_d_unspl_prev=np.zeros(r),
    )
    base.update(overrides)
    return base


def test_strand_favoring_rna_drives_pi_g_to_zero():
    # Perfectly stranded library, all-sense region → RNA Beta-Binomial fits, gDNA
    # (symmetric about 0.5) does not → π_g → 0, mass → RNA.
    view = _view([95], [5])
    alloc = estep_view(view, np.array([TS_POS], dtype=np.int8), **_params(1))
    assert isinstance(alloc, Allocation)
    assert alloc.pi_g[0] < 0.01
    assert alloc.m_d[0] > alloc.m_g[0]


def test_neutral_strand_gives_prior_pi_g():
    # Unstranded library (κ_rna = 0.5) with ρ_r == ρ_d → RNA and gDNA BB coincide,
    # LLR_strand = 0; with the count channel silent π_g stays at the prior 0.5.
    view = _view([60], [40])
    alloc = estep_view(
        view, np.array([TS_POS], dtype=np.int8), **_params(1, kappa_rna=0.5, rho_r_bb=0.01)
    )
    np.testing.assert_allclose(alloc.pi_g, [0.5])


def test_ambiguous_region_skips_strand_term():
    # Same asymmetric counts, but AMBIG (both strands) → strand term omitted (D7),
    # count silent → π_g stays at the prior. A decodable POS region with the same
    # counts gets a non-prior π_g.
    view = _view([95], [5])
    ambig = estep_view(view, np.array([TS_AMBIG], dtype=np.int8), **_params(1))
    decodable = estep_view(view, np.array([TS_POS], dtype=np.int8), **_params(1))
    np.testing.assert_allclose(ambig.pi_g, [0.5])
    assert decodable.pi_g[0] < 0.01


def test_none_keeps_strand_term_neutrally():
    # Intergenic (TS_NONE) keeps the strand term; gDNA is unstranded so an
    # arbitrary sense (genome '+') is used. With a balanced count it is neutral.
    view = _view([50], [50])
    alloc = estep_view(
        view, np.array([TS_NONE], dtype=np.int8), **_params(1, kappa_rna=0.5, rho_r_bb=0.01)
    )
    np.testing.assert_allclose(alloc.pi_g, [0.5])


def test_count_channel_activates_with_prior_rna_mass():
    # With a non-zero μ_d the count channel discriminates (κ_rna=0.5 so strand is
    # silent): a count consistent with the gDNA-only mean favours gDNA; a count
    # matching the larger mixture mean favours RNA.
    ts = np.array([TS_POS], dtype=np.int8)
    consistent = estep_view(
        _view([100], [0]),
        ts,
        **_params(1, kappa_rna=0.5, rho_r_bb=0.01, rho_0=1.0, m_d_unspl_prev=np.array([50.0])),
    )
    excess = estep_view(
        _view([100], [0]),
        ts,
        **_params(1, kappa_rna=0.5, rho_r_bb=0.01, rho_0=0.5, m_d_unspl_prev=np.array([50.0])),
    )
    assert consistent.pi_g[0] > 0.5  # n_u ≈ μ_g → gDNA-favoured
    assert excess.pi_g[0] < 0.5  # n_u ≈ μ_g + μ_d → RNA-favoured


def test_count_channel_silent_when_no_prior_rna_mass():
    # μ_d = 0 (the PR 4a single pass): gDNA-only and mixture coincide → LLR=0.
    view = _view([100], [0])
    alloc = estep_view(
        view, np.array([TS_POS], dtype=np.int8), **_params(1, kappa_rna=0.5, rho_r_bb=0.01)
    )
    np.testing.assert_allclose(alloc.pi_g, [0.5])


def test_mass_is_conserved():
    # m_g + m_d == mass_unspliced + mass_spliced, for any allocation.
    view = _view([20, 5, 0], [3, 30, 0], sense=[4, 0, 0], anti=[0, 6, 0])
    ts = np.array([TS_POS, TS_NEG, TS_NONE], dtype=np.int8)
    alloc = estep_view(view, ts, **_params(3))
    total = view.mass_unspliced + view.mass_spliced
    np.testing.assert_allclose(alloc.m_g + alloc.m_d, total)
    assert np.all(alloc.m_g >= 0.0) and np.all(alloc.m_d >= 0.0)


def test_flux_drives_pi_mass_is_allocated():
    # D1: the integer flux (95/5 → strand says RNA) drives π_g, but the fractional
    # boundary mass (10.0, not 100) is what gets allocated.
    view = _view([95], [5], mass_u=[10.0], mass_s=[0.0])
    alloc = estep_view(view, np.array([TS_POS], dtype=np.int8), **_params(1))
    assert alloc.pi_g[0] < 0.01  # π from the 100-read flux
    np.testing.assert_allclose(alloc.m_g_unspl, alloc.pi_g * 10.0)  # mass from 10.0
    assert alloc.m_g[0] < 10.0


def test_spliced_mass_routes_to_rna():
    # Spliced mass is deterministic RNA up to the ε_s gDNA-artifact failsafe.
    view = _view([0], [0], sense=[40], anti=[0])
    alloc = estep_view(view, np.array([TS_POS], dtype=np.int8), **_params(1, eps_s=1.0e-3))
    np.testing.assert_allclose(alloc.m_g, [1.0e-3 * 40])
    np.testing.assert_allclose(alloc.m_d, [(1.0 - 1.0e-3) * 40])
