"""gDNA / RNA FL distributions: the five PURE pools + the smooth-EB build.

The pools are 's, and purity is the whole point: a length model is
fitted only from populations known to be ONE component, so nothing is ever estimated from the fragments
it will later explain.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.fl import (
    _fl_models_from_histograms,
    gdna_fl_mass,
    rna_fl_mass,
    splash_fl_mass,
)
from rigel.scan_payload import (
    N_FRAGMENT_POOLS,
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTERGENIC_EXON,
    POOL_DNA_INTRONIC,
    POOL_DNA_INTRON_EXON,
    POOL_RNA_SPLICED,
)


def _spike(at: int, total: float, n: int = 1001) -> np.ndarray:
    c = np.zeros(n, dtype=np.float64)
    c[at] = total
    return c


def _pools(n_bins: int = 5) -> np.ndarray:
    """A ``pool_lengths`` block with a different value in every pool, so a wrong index cannot pass."""
    pools = np.zeros((N_FRAGMENT_POOLS, n_bins), dtype=np.int64)
    pools[POOL_DNA_INTERGENIC, 2] = 1
    pools[POOL_DNA_INTRONIC, 3] = 3
    pools[POOL_DNA_INTRON_EXON, 2] = 700
    pools[POOL_DNA_INTERGENIC_EXON, 3] = 900
    pools[POOL_RNA_SPLICED, 4] = 11
    return pools


def test_the_pool_indices_ARE_the_specifications_FragmentPool():
    """⛔ The three-way contract: the reference enum, the C++ enum and these constants are one axis.

    A silent disagreement here re-labels every pool — the gDNA model would be fitted from the RNA pool
    and nothing would look wrong. Checked against the executable specification itself, not a written-out
    list, for the same reason the payload schema test does.
    """
    from tests.native._accumulator_reference import FragmentPool

    assert N_FRAGMENT_POOLS == len(FragmentPool)
    assert POOL_DNA_INTERGENIC == FragmentPool.DNA_INTERGENIC
    assert POOL_DNA_INTRONIC == FragmentPool.DNA_INTRONIC
    assert POOL_DNA_INTRON_EXON == FragmentPool.DNA_INTRON_EXON
    assert POOL_DNA_INTERGENIC_EXON == FragmentPool.DNA_INTERGENIC_EXON
    assert POOL_RNA_SPLICED == FragmentPool.RNA_SPLICED


def test_gdna_fl_mass_is_the_two_PURE_CONTAINED_pools_and_NOTHING_else():
    """⭐ The splash pools are deliberately NOT folded in, and this is not cosmetic.

    They are the only ON-TARGET gDNA population, so they sit between
    the pure gDNA and RNA means — 139 and 212 against 88 on LBX0190. The shipped model summed four
    differently-tilted pools and read **146.05** where the pure intergenic pool says **88.0**: biased
    long by ~40 %, by pooling exactly these. Keeping them out is what makes the comparison a QC output
    instead of a guess.
    """
    g = gdna_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(g, [0.0, 0.0, 1.0, 3.0, 0.0])


def test_gdna_fl_mass_excludes_the_RNA_pool():
    """The circularity guard: the gDNA length model must never see a certified-RNA fragment."""
    pools = np.zeros((N_FRAGMENT_POOLS, 5), dtype=np.int64)
    pools[POOL_RNA_SPLICED, 1] = 12345
    assert float(gdna_fl_mass(SimpleNamespace(pool_lengths=pools)).sum()) == 0.0


def test_rna_fl_mass_is_the_ANNOTATED_JUNCTION_pool_alone():
    """gDNA cannot be spliced, so an observed annotated junction certifies RNA — and only that pool
    does. ⚠ ``sj_implicit`` fragments are already excluded by the accumulator, because a splice that was
    never sequenced is a product of the very model this pool is used to fit."""
    r = rna_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(r, [0.0, 0.0, 0.0, 0.0, 11.0])


def test_the_splash_pools_are_reachable_SEPARATELY_for_QC():
    """Named pools, not folded in — so 'is the off-target model mis-centred for the fragments that
    actually leak?' is answerable rather than assumed."""
    s = splash_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(s, [0.0, 0.0, 700.0, 900.0, 0.0])


def test_every_pool_reaches_EXACTLY_ONE_accessor():
    """⛔ Teeth on the partition itself: no pool may be double-counted, and none may be unreachable.

    A pool that no accessor returns is silently discarded evidence; one that two return is
    double-counted. Both are invisible in any single-accessor test.
    """
    payload = SimpleNamespace(pool_lengths=_pools())
    totals = [
        float(gdna_fl_mass(payload).sum()),
        float(rna_fl_mass(payload).sum()),
        float(splash_fl_mass(payload).sum()),
    ]
    assert sum(totals) == float(_pools().sum()), "every pool must land in exactly one accessor"
    for pool in range(N_FRAGMENT_POOLS):
        one = np.zeros((N_FRAGMENT_POOLS, 5), dtype=np.int64)
        one[pool, 1] = 5
        p = SimpleNamespace(pool_lengths=one)
        reached = [
            float(gdna_fl_mass(p).sum()),
            float(rna_fl_mass(p).sum()),
            float(splash_fl_mass(p).sum()),
        ]
        assert sorted(reached) == [0.0, 0.0, 5.0], (
            f"pool {pool} reaches {reached}, expected exactly one"
        )


def test_build_fl_large_pool_is_empirical():
    glob = _spike(100, 1.0e6)
    gdna = _spike(300, 1.0e7)  # huge gDNA pool ⇒ dominated by its own evidence
    fl = _fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=gdna, max_size=1000, prior_ess=1000.0
    )
    assert int(np.argmax(fl.gdna_pmf)) == 300
    np.testing.assert_allclose(fl.gdna_pmf.sum(), 1.0)
    np.testing.assert_allclose(fl.global_pmf.sum(), 1.0)


def test_build_fl_empty_pool_collapses_to_global():
    glob = _spike(100, 1.0e6)
    fl = _fl_models_from_histograms(
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
    fl = _fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=gdna, max_size=1000, prior_ess=1000.0
    )
    # Smooth shrinkage (no threshold): with pool_total=1 ≪ ρ_ess=1000 the global
    # anchor (bin 100) dominates the lone gDNA fragment (bin 300).
    assert fl.gdna_pmf[100] > fl.gdna_pmf[300]
    assert fl.gdna_pmf[300] > 0.0  # but the empirical fragment still nudges its bin


def test_build_fl_just_below_5000_is_not_a_cliff():
    # Decision 5: 4999 vs 5001 differ only marginally (no GOOD/WEAK jump).
    glob = _spike(100, 1.0e6)
    lo = _fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=_spike(300, 4999.0), max_size=1000
    )
    hi = _fl_models_from_histograms(
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
    fl = _fl_models_from_histograms(
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
    fl = _fl_models_from_histograms(
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
    fl = _fl_models_from_histograms(
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
