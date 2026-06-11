"""Phase-1 strand module: the standalone likelihood emitter (`strand_deconvolve`).

Pins the redesign's Phase-1 contract: per node, the gDNA fraction read at the FP-quantile + the
unspliced mass split + the **likelihood** information `I = N·(2κ−1)²` that the count module weights by
(0 / NaN where strand-uninformative). Built alongside the live blend; no count input.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from rigel.calibration.strand_deconv import (
    StrandSplit,
    _grid_posterior_quantile,
    strand_deconvolve,
    strand_posterior_gdna_frac,
)
from rigel.calibration.substrate import CalibrationSubstrate, SubstrateView


def _view(pos, neg, mass=None):
    pos = np.asarray(pos, dtype=np.int64)
    neg = np.asarray(neg, dtype=np.int64)
    n = pos.shape[0]
    mass = (pos + neg).astype(np.float64) if mass is None else np.asarray(mass, dtype=np.float64)
    z = np.zeros(n, dtype=np.int64)
    return SubstrateView(
        n_unspliced_pos=pos, n_unspliced_neg=neg, n_spliced_sense=z, n_spliced_antisense=z,
        mass_unspliced=mass, mass_spliced=np.zeros(n),
    )


def _substrate(strand_class, c_pos, c_neg, c_mass=None, left=None, right=None):
    n = len(strand_class)
    zero = _view(np.zeros(n), np.zeros(n))
    return CalibrationSubstrate(
        n_regions=n, region_len=np.full(n, 1000.0),
        strand_class=np.asarray(strand_class, dtype=np.int8),
        contained=_view(c_pos, c_neg, c_mass),
        left=zero if left is None else left,
        right=zero if right is None else right,
    )


def _ra(strand_class, ref_id=None):
    n = len(strand_class)
    ref_id = np.zeros(n, dtype=np.int64) if ref_id is None else np.asarray(ref_id)
    return SimpleNamespace(strand_class=np.asarray(strand_class), ref_id=ref_id)


def _deconv(sub, ra, *, ss=0.99, q=0.5):
    return strand_deconvolve(
        sub, ra, rna_sense_frac=ss,
        gdna_strand_overdispersion=0.0, rna_strand_overdispersion=0.0, deconv_quantile=q,
    )


def test_confident_balanced_region_reads_gdna():
    # A balanced POS region (sense=antisense) at high SS is gDNA's signature → all gDNA; info = N·(2ss−1)².
    ts = [TS_POS]
    c, _, _ = _deconv(_substrate(ts, [50], [50], [100.0]), _ra(ts), ss=0.99)
    assert c.gdna_frac[0] > 0.9
    assert c.gdna_mass[0] > 90.0
    assert c.info[0] == pytest.approx(100 * (2 * 0.99 - 1) ** 2)
    assert c.gdna_mass[0] + c.rna_mass[0] == pytest.approx(100.0)  # mass conserved where split


def test_rna_skewed_region_reads_rna():
    # Sense fraction == ss (pos=90/100, ss=0.9) is RNA's signature → ~all RNA.
    ts = [TS_POS]
    c, _, _ = _deconv(_substrate(ts, [90], [10], [100.0]), _ra(ts), ss=0.9)
    assert c.gdna_frac[0] < 0.1


def test_ambig_region_not_split():
    ts = [TS_AMBIG]
    c, _, _ = _deconv(_substrate(ts, [40], [40], [80.0]), _ra(ts), ss=0.99)
    assert np.isnan(c.gdna_frac[0])
    assert c.info[0] == 0.0
    assert np.isnan(c.gdna_mass[0]) and np.isnan(c.rna_mass[0])  # count module imputes these


def test_unstranded_gives_no_info():
    ts = [TS_POS, TS_NEG]
    c, _, _ = _deconv(_substrate(ts, [50, 50], [50, 50]), _ra(ts), ss=0.5)
    assert np.all(c.info == 0.0)  # (2·½−1)² = 0
    assert np.all(np.isnan(c.gdna_frac))


def test_info_uses_count_not_mass():
    # info = N·(2ss−1)² with N the unspliced COUNT (pos+neg=50), independent of the mass (500).
    ts = [TS_POS]
    c, _, _ = _deconv(_substrate(ts, [30], [20], [500.0]), _ra(ts), ss=0.8)
    assert c.info[0] == pytest.approx(50 * (2 * 0.8 - 1) ** 2)


def test_quantile_monotone_in_q():
    ts = [TS_POS]
    sub, ra = _substrate(ts, [30], [20], [50.0]), _ra(ts)
    fracs = [_deconv(sub, ra, ss=0.9, q=q)[0].gdna_frac[0] for q in (0.1, 0.5, 0.9)]
    assert fracs[0] <= fracs[1] <= fracs[2]


def test_quantile_half_is_the_posterior_median():
    ts = [TS_POS]
    g_q = _deconv(_substrate(ts, [30], [20], [50.0]), _ra(ts), ss=0.9, q=0.5)[0].gdna_frac[0]
    med, _ = strand_posterior_gdna_frac(
        np.array([30.0]), np.array([20.0]), 0.9,
        gdna_strand_overdispersion=0.0, rna_strand_overdispersion=0.0, n_grid=200, deconv_quantile=0.5,
    )
    assert g_q == pytest.approx(med[0])


def test_grid_quantile_half_matches_np_interp_median():
    rng = np.random.default_rng(0)
    post = rng.random((50, 200))
    post /= post.sum(axis=1, keepdims=True)
    grid = np.linspace(1e-4, 1 - 1e-4, 200)
    mine = _grid_posterior_quantile(post, grid, 0.5)
    ref = np.array([np.interp(0.5, np.cumsum(post[i]), grid) for i in range(post.shape[0])])
    np.testing.assert_allclose(mine, ref, atol=1e-12)


def test_boundary_side_oriented_and_informative():
    # Region 1 (POS) left side neighbours region 0 (NONE intergenic) → strand-observable (wildcard);
    # its crossing count (40) drives info; region 0 (no left neighbour) is silent.
    ts = [TS_NONE, TS_POS]
    left = _view([0, 30], [0, 10], [0.0, 40.0])  # region 1's left side carries the crossing flux
    sub = _substrate(ts, [10, 10], [10, 10], left=left)
    _, lft, _ = _deconv(sub, _ra(ts, ref_id=[0, 0]), ss=0.95)
    assert isinstance(lft, StrandSplit) and lft.info.shape == (2,)
    assert lft.info[0] == 0.0  # reference edge: no left neighbour
    assert lft.info[1] == pytest.approx(40 * (2 * 0.95 - 1) ** 2)
