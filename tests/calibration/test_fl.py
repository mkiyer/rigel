"""gDNA / RNA FL distributions: gDNA pool aggregation + smooth-EB build (PR 4c)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.fl import N_FL_POOLS


def _spike(at: int, total: float, n: int = 1001) -> np.ndarray:
    c = np.zeros(n, dtype=np.float64)
    c[at] = total
    return c


def test_gdna_fl_mass_sums_intergenic_and_intronic_excludes_exonic():
    # Pools 0,1 = intergenic (contained, boundary); 2,3 = intronic; 4,5 = exonic.
    fp = np.zeros((N_FL_POOLS, 5), dtype=np.float64)
    fp[0, 2] = 1.0
    fp[1, 2] = 2.0  # intergenic
    fp[2, 3] = 3.0
    fp[3, 3] = 4.0  # intronic
    fp[4, 2] = 99.0
    fp[5, 3] = 99.0  # exonic — must be EXCLUDED
    g = gdna_fl_mass(SimpleNamespace(fl_pool_mass=fp))
    np.testing.assert_allclose(g, [0.0, 0.0, 3.0, 7.0, 0.0])  # bin2: 1+2; bin3: 3+4


def test_gdna_fl_mass_empty_when_pools_absent():
    assert gdna_fl_mass(SimpleNamespace(fl_pool_mass=None)).size == 0


def test_build_fl_large_pool_is_empirical():
    glob = _spike(100, 1.0e6)
    gdna = _spike(300, 1.0e7)  # huge gDNA pool ⇒ dominated by its own evidence
    fl = build_fl_models(
        global_counts=glob, rna_counts=glob, gdna_counts=gdna, max_size=1000, prior_ess=1000.0
    )
    assert int(np.argmax(fl.gdna_pmf)) == 300
    np.testing.assert_allclose(fl.gdna_pmf.sum(), 1.0)
    np.testing.assert_allclose(fl.global_pmf.sum(), 1.0)


def test_build_fl_empty_pool_collapses_to_global():
    glob = _spike(100, 1.0e6)
    fl = build_fl_models(
        global_counts=glob,
        rna_counts=glob,
        gdna_counts=np.zeros(1001, dtype=np.float64),
        max_size=1000,
    )
    np.testing.assert_allclose(fl.gdna_pmf, fl.global_pmf)  # no gDNA evidence → global
    assert fl.n_gdna == 0.0


def test_build_fl_small_pool_shrinks_toward_global_no_cliff():
    glob = _spike(100, 1.0e6)  # global mode at 100
    gdna = _spike(300, 1.0)  # a single gDNA fragment at 300
    fl = build_fl_models(
        global_counts=glob, rna_counts=glob, gdna_counts=gdna, max_size=1000, prior_ess=1000.0
    )
    # Smooth shrinkage (no threshold): with pool_total=1 ≪ ρ_ess=1000 the global
    # anchor (bin 100) dominates the lone gDNA fragment (bin 300).
    assert fl.gdna_pmf[100] > fl.gdna_pmf[300]
    assert fl.gdna_pmf[300] > 0.0  # but the empirical fragment still nudges its bin


def test_build_fl_just_below_5000_is_not_a_cliff():
    # Decision 5: 4999 vs 5001 differ only marginally (no GOOD/WEAK jump).
    glob = _spike(100, 1.0e6)
    lo = build_fl_models(
        global_counts=glob, rna_counts=glob, gdna_counts=_spike(300, 4999.0), max_size=1000
    )
    hi = build_fl_models(
        global_counts=glob, rna_counts=glob, gdna_counts=_spike(300, 5001.0), max_size=1000
    )
    np.testing.assert_allclose(lo.gdna_pmf[300], hi.gdna_pmf[300], rtol=1e-3)


# ---------------------------------------------------------------------------
# Raw empirical counts + FragmentLengthModel accessors (QC report views)
# ---------------------------------------------------------------------------


def test_build_fl_stores_raw_aligned_counts():
    """Raw (unsmoothed) counts are stored aligned to max_size, totals = n_*."""
    glob = _spike(100, 500.0)
    rna = _spike(200, 300.0)
    gdna = _spike(300, 40.0)
    fl = build_fl_models(
        global_counts=glob, rna_counts=rna, gdna_counts=gdna, max_size=1000, prior_ess=1000.0
    )
    assert fl.global_counts.shape == (1001,)
    assert int(np.argmax(fl.rna_counts)) == 200
    assert int(np.argmax(fl.gdna_counts)) == 300
    # Raw counts are the empirical evidence — NOT the EB-smoothed pmf.
    assert fl.rna_counts[200] == 300.0
    assert fl.gdna_counts[300] == 40.0
    assert fl.n_global == 500.0
    assert fl.n_rna == 300.0
    assert fl.n_gdna == 40.0


def test_build_fl_counts_overflow_folds_into_last_bin():
    """Counts beyond max_size fold into the overflow bin when aligned."""
    over = np.zeros(1201, dtype=np.float64)
    over[1150] = 7.0  # beyond max_size=1000
    fl = build_fl_models(
        global_counts=_spike(100, 10.0),
        rna_counts=_spike(100, 10.0),
        gdna_counts=over,
        max_size=1000,
    )
    assert fl.gdna_counts.shape == (1001,)
    assert fl.gdna_counts[1000] == 7.0
    assert fl.n_gdna == 7.0


def test_accessors_return_empirical_models():
    """rna_model()/gdna_model()/global_model() expose the raw empirical FL."""
    fl = build_fl_models(
        global_counts=_spike(100, 500.0),
        rna_counts=_spike(200, 300.0),
        gdna_counts=_spike(325, 40.0),
        max_size=1000,
        prior_ess=1000.0,
    )
    rna_m = fl.rna_model()
    gdna_m = fl.gdna_model()
    glob_m = fl.global_model()
    # Empirical modes/means track the raw counts (not smoothed toward global).
    assert rna_m.mode == 200
    assert gdna_m.mode == 325
    assert glob_m.mode == 100
    assert rna_m.mean == pytest.approx(200.0)
    # n_observations = pool total.
    assert rna_m.n_observations == 300
    assert gdna_m.n_observations == 40
    # to_dict() works (drives the summary QC report).
    assert gdna_m.to_dict()["summary"]["mode"] == 325
