"""Regression tests for the §1.4 fused locus-prior helpers.

The optimisation hoists ``boundary_crossing_exposure`` and the per-Locus
``region_index.overlap`` + ``contained_exposure_clipped`` work out of
``estimate_locus_gdna`` / ``expected_gdna_count_global`` via two new
optional kwargs (``b_cross`` and ``_scratch``).

These tests pin the contract that:

* passing the precomputed ``b_cross`` reproduces the legacy result
* passing the precomputed ``_scratch`` reproduces the legacy result
* both helpers consuming the *same* ``_scratch`` produce identical
  outputs to two independent calls
"""

from __future__ import annotations

import pytest

from rigel.calibration._exposure import boundary_crossing_exposure
from rigel.calibration.locus_prior import (
    _compute_locus_scratch,
    estimate_locus_gdna,
    expected_gdna_count_global,
)

# Reuse the heavy fixture builders from the per-locus mass test module
# so we exercise the same code paths the integration tests cover.
from test_per_locus_gdna_mass import (  # type: ignore[import-not-found]
    _build_arrays,
    _delta_fl,
    _gdt,
)
from rigel.calibration.density_global import GlobalDensityTable, GlobalGdnaDensity, KappaEstimate  # noqa: E402,F401
from rigel.calibration.regions import RegionType  # noqa: E402,F401
from rigel.calibration._arrays import PayloadArrays  # noqa: E402,F401
from rigel.calibration._region_index_py import RegionIndexPy  # noqa: E402,F401
from rigel.locus import Locus  # noqa: E402


def _setup():
    ra, pa, idx, _ = _build_arrays(
        regions=[
            ("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False),
            ("chr1", 1000, 2000, int(RegionType.EXON), True, True),
            ("chr1", 2000, 3000, int(RegionType.INTRON), False, False),
        ],
        counts_intergenic=[5],
        counts_intron=[3],
        u_left=[0, 2, 0],
        u_right=[0, 1, 0],
    )
    gdna_fl = _delta_fl(200)
    gdt = _gdt(rho_ig=0.0, fl_mean=200)
    locus = Locus(ref="chr1", ref_id=0, start=500, end=2500)
    return ra, pa, idx, gdt, gdna_fl, locus


def _est_equal(a, b):
    assert a.n_gdna == pytest.approx(b.n_gdna)
    assert a.n_gdna_intergenic == pytest.approx(b.n_gdna_intergenic)
    assert a.n_gdna_intron == pytest.approx(b.n_gdna_intron)
    assert a.n_gdna_boundary_observed == pytest.approx(b.n_gdna_boundary_observed)
    assert a.n_gdna_exon_only == pytest.approx(b.n_gdna_exon_only)
    assert a.pi_gdna == pytest.approx(b.pi_gdna)
    assert a.fallback_flags == b.fallback_flags


def _eg_equal(a, b):
    assert a.total == pytest.approx(b.total)
    assert a.intergenic_contained == pytest.approx(b.intergenic_contained)
    assert a.intron_contained == pytest.approx(b.intron_contained)
    assert a.boundary_crossing_expected == pytest.approx(b.boundary_crossing_expected)
    assert a.exon_contained_expected == pytest.approx(b.exon_contained_expected)


def test_b_cross_override_matches_default_estimate():
    ra, pa, idx, gdt, gdna_fl, locus = _setup()
    legacy = estimate_locus_gdna(
        locus=locus, n_obs=10, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl=gdna_fl,
    )
    b = boundary_crossing_exposure(gdna_fl)
    overridden = estimate_locus_gdna(
        locus=locus, n_obs=10, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl=gdna_fl,
        b_cross=b,
    )
    _est_equal(overridden, legacy)


def test_b_cross_override_matches_default_expected():
    ra, pa, idx, gdt, gdna_fl, locus = _setup()
    legacy = expected_gdna_count_global(
        locus=locus, region_index=idx,
        region_arrays=ra, global_densities=gdt, gdna_fl=gdna_fl,
    )
    b = boundary_crossing_exposure(gdna_fl)
    overridden = expected_gdna_count_global(
        locus=locus, region_index=idx,
        region_arrays=ra, global_densities=gdt, gdna_fl=gdna_fl,
        b_cross=b,
    )
    _eg_equal(overridden, legacy)


def test_shared_scratch_matches_legacy_paths():
    ra, pa, idx, gdt, gdna_fl, locus = _setup()
    # Legacy: each helper computes its own scratch internally.
    legacy_est = estimate_locus_gdna(
        locus=locus, n_obs=10, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl=gdna_fl,
    )
    legacy_eg = expected_gdna_count_global(
        locus=locus, region_index=idx,
        region_arrays=ra, global_densities=gdt, gdna_fl=gdna_fl,
    )

    # Optimised: build scratch once, share between both helpers, and
    # hoist b_cross.
    scratch = _compute_locus_scratch(
        locus, idx, ra, gdna_fl,
        intergenic_flank_bp=5_000, ref_length=None,
    )
    assert scratch is not None
    b = boundary_crossing_exposure(gdna_fl)
    fast_est = estimate_locus_gdna(
        locus=locus, n_obs=10, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl=gdna_fl,
        b_cross=b, _scratch=scratch,
    )
    fast_eg = expected_gdna_count_global(
        locus=locus, region_index=idx,
        region_arrays=ra, global_densities=gdt, gdna_fl=gdna_fl,
        b_cross=b, _scratch=scratch,
    )
    _est_equal(fast_est, legacy_est)
    _eg_equal(fast_eg, legacy_eg)


def test_scratch_with_zero_flank_matches_unflanked_path():
    """When ``intergenic_flank_bp == 0`` the scratch's flanked overlap
    collapses to the unflanked one; both helpers must still match.
    """
    ra, pa, idx, gdt, gdna_fl, locus = _setup()
    legacy_est = estimate_locus_gdna(
        locus=locus, n_obs=10, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl=gdna_fl,
        intergenic_flank_bp=0,
    )
    legacy_eg = expected_gdna_count_global(
        locus=locus, region_index=idx,
        region_arrays=ra, global_densities=gdt, gdna_fl=gdna_fl,
    )
    scratch = _compute_locus_scratch(
        locus, idx, ra, gdna_fl,
        intergenic_flank_bp=0, ref_length=None,
    )
    fast_est = estimate_locus_gdna(
        locus=locus, n_obs=10, region_index=idx,
        region_arrays=ra, payload_arrays=pa,
        global_densities=gdt, gdna_fl=gdna_fl,
        intergenic_flank_bp=0, _scratch=scratch,
    )
    fast_eg = expected_gdna_count_global(
        locus=locus, region_index=idx,
        region_arrays=ra, global_densities=gdt, gdna_fl=gdna_fl,
        _scratch=scratch,
    )
    _est_equal(fast_est, legacy_est)
    _eg_equal(fast_eg, legacy_eg)
