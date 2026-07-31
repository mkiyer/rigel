"""Unit tests for the aggregate DNA-background reference (`measure_background`)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.background_reference import BackgroundReference, measure_background
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS


def _substrate(sig, pos, neg, eff):
    """Minimal stand-ins for the ``(substrate, region_arrays, node_eff_len)`` inputs.

    ⭐ ``node_contained`` is ONE array with the two GENOME-strand columns side by side, not a pair of
    ``n_unspliced_pos``/``n_unspliced_neg`` fields — and ``measure_background`` sums them, because gDNA
    is strand-symmetric and a background rate is a total.
    """
    ra = SimpleNamespace(signature=np.asarray(sig))
    sub = SimpleNamespace(
        node_contained=SimpleNamespace(
            count=np.stack([np.asarray(pos, float), np.asarray(neg, float)], axis=1)
        )
    )
    return sub, ra, np.asarray(eff, float)


def test_pools_intergenic_only():
    """Only signature-0 (intergenic) regions enter the pool; exon and intron regions are excluded (introns
    carry nascent RNA — the whole point of intergenic-only)."""
    sub, ra, eff = _substrate(
        sig=[0, 0, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_NEG | BIT_INTRON_NEG],
        pos=[10, 20, 500, 300, 400],
        neg=[10, 20, 500, 300, 400],
        eff=[1000, 1000, 1000, 1000, 1000],
    )
    bg = measure_background(sub, ra, eff)
    # only the two intergenic regions: Σg = 20 + 40 = 60 over ΣE = 2000
    assert bg.n_counts == 60.0
    assert bg.eff_total == 2000.0
    assert np.isclose(bg.log_rho_bg, np.log(60.0 / 2000.0))
    assert np.isclose(bg.sigma_bg, 1.0 / np.sqrt(60.0))


def test_zero_counts_is_dormant():
    """A DNA-free (or fully depleted) pool ⇒ log ρ_bg = −inf (the consuming floor is dormant), σ_bg = +inf."""
    sub, ra, eff = _substrate(sig=[0, 0, 0], pos=[0, 0, 0], neg=[0, 0, 0], eff=[1000, 2000, 3000])
    bg = measure_background(sub, ra, eff)
    assert bg.log_rho_bg == -np.inf
    assert bg.sigma_bg == np.inf
    assert bg.n_counts == 0.0
    assert (
        bg.eff_total == 6000.0
    )  # support is still measured (the aggregate exists; the signal is zero)


def test_intron_with_nascent_does_not_leak_into_background():
    """An intron heavy with nascent RNA must NOT raise ρ_bg — the regression guard for the false-positive that
    all-count-observable (intergenic ∪ intron) would cause. Here the only intergenic region is DNA-free."""
    sub, ra, eff = _substrate(
        sig=[0, BIT_INTRON_POS, BIT_INTRON_NEG],
        pos=[0, 9999, 9999],
        neg=[0, 9999, 9999],
        eff=[5000, 1000, 1000],
    )
    bg = measure_background(sub, ra, eff)
    assert bg.log_rho_bg == -np.inf  # intergenic is empty ⇒ background zero, introns ignored
    assert bg.n_counts == 0.0


def test_no_intergenic_regions():
    """No signature-0 region at all ⇒ empty pool ⇒ dormant (no crash)."""
    sub, ra, eff = _substrate(
        sig=[BIT_EXON_POS, BIT_INTRON_POS], pos=[5, 5], neg=[5, 5], eff=[100, 100]
    )
    bg = measure_background(sub, ra, eff)
    assert bg.log_rho_bg == -np.inf
    assert bg.eff_total == 0.0


def test_zero_efflength_regions_excluded():
    """Intergenic regions with no usable eff-length are dropped from the pool (no divide-by-zero)."""
    sub, ra, eff = _substrate(sig=[0, 0], pos=[10, 999], neg=[10, 999], eff=[1000, 0.0])
    bg = measure_background(sub, ra, eff)
    assert bg.n_counts == 20.0  # the E=0 region's 1998 counts are excluded
    assert bg.eff_total == 1000.0


def test_deterministic():
    sub, ra, eff = _substrate(
        sig=[0, 0, BIT_EXON_POS], pos=[10, 20, 5], neg=[11, 21, 5], eff=[1e3, 2e3, 1e3]
    )
    a = measure_background(sub, ra, eff)
    b = measure_background(sub, ra, eff)
    assert a == b
    assert isinstance(a, BackgroundReference)


def test_include_introns_adds_intron_span():
    """``include_introns=True`` pools intergenic + intron (real-data resolution); exons stay excluded."""
    sub, ra, eff = _substrate(
        sig=[0, BIT_INTRON_POS, BIT_INTRON_NEG, BIT_EXON_POS],
        pos=[5, 50, 50, 999],
        neg=[5, 50, 50, 999],
        eff=[1000, 1000, 1000, 1000],
    )
    intergenic = measure_background(sub, ra, eff)  # default: Σg = 10
    assert intergenic.n_counts == 10.0
    with_introns = measure_background(
        sub, ra, eff, include_introns=True
    )  # + two introns: 10 + 100 + 100
    assert with_introns.n_counts == 210.0  # exon (1998) still excluded


def test_robust_trim_drops_nascent_outlier_intron():
    """The MAD fence drops a single nascent-heavy intron (a high per-region density outlier) so it cannot
    inflate the pooled background — the real-data safeguard for sparse nascent."""
    bg = [
        8.0,
        9.0,
        10.0,
        11.0,
        12.0,
        8.0,
        9.0,
        10.0,
        11.0,
        12.0,
    ]  # background ~ density 0.01, varied
    sig = [0] * len(bg) + [BIT_INTRON_POS]  # + one nascent-contaminated intron
    counts = bg + [5000.0]  # the outlier: density 5.0 ≫ 0.01
    eff = [1000.0] * (len(bg) + 1)
    sub, ra, e = _substrate(
        sig=sig, pos=[c / 2 for c in counts], neg=[c / 2 for c in counts], eff=eff
    )
    plain = measure_background(sub, ra, e, include_introns=True)
    trimmed = measure_background(sub, ra, e, include_introns=True, robust_trim_mad=5.0)
    assert (
        trimmed.log_rho_bg < plain.log_rho_bg
    )  # the outlier is removed ⇒ lower, correct background
    assert trimmed.n_counts == float(
        sum(bg)
    )  # only the background regions remain (outlier dropped)
