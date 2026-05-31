"""Global gDNA-density seed: intergenic + intronic (+ boundaries), no clip."""

from __future__ import annotations

import numpy as np
import pandas as pd

from _synthetic import make_synthetic_payload

from rigel.calibration.density import GdnaDensity, estimate_gdna_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.scan_payload import AccumulatorPayload


def _substrate_from_signatures(contained, signatures) -> tuple[CalibrationSubstrate, RegionArrays]:
    """A boundary-free substrate over regions with the given signatures."""
    contained = np.asarray(contained, dtype=np.uint32)
    r = contained.shape[0]
    zf = np.zeros((r + 1, 4), dtype=np.float32)
    zu = np.zeros((r + 1, 4), dtype=np.uint32)
    payload = AccumulatorPayload(
        boundaries=np.arange(r + 1, dtype=np.int64) * 100,
        ref_pos_offsets=np.array([0, r + 1], dtype=np.int64),
        ref_region_offsets=np.array([0, r], dtype=np.int64),
        ref_boundary_offsets=np.array([0, r + 1], dtype=np.int64),
        region_contained=contained,
        boundary_mass_left=zf,
        boundary_mass_right=zf.copy(),
        boundary_flux_left=zu,
        boundary_flux_right=zu.copy(),
        n_refs=1,
    )
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(r, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * r, dtype="string"),
            "start": np.arange(r, dtype=np.int64) * 100,
            "end": (np.arange(r, dtype=np.int64) + 1) * 100,
            "length": np.full(r, 100, dtype=np.int64),
            "signature": np.asarray(signatures, dtype=np.uint8),
        }
    )
    ra = RegionArrays.from_region_df(region_df, {"chr1": 0})
    return CalibrationSubstrate.from_payload(payload, ra), ra


def test_seed_is_intergenic_only_in_basic_synthetic():
    # make_synthetic_payload: r0 +exon, r1 −exon, r2 intergenic. Only r2 seeds.
    # seed mass = r2 contained (15) + left boundary mass (1.5) + right (0) = 16.5.
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)
    d = estimate_gdna_density(sub, ra)
    assert isinstance(d, GdnaDensity)
    assert d.n_seed_regions == 1
    assert not d.fallback_used
    np.testing.assert_allclose(d.seed_mass, 16.5)  # includes the boundary mass
    np.testing.assert_allclose(d.seed_length, 100.0)
    np.testing.assert_allclose(d.rho_0, 17.5 / 101.0)  # (1 + 16.5) / (1 + 100)


def test_intronic_seeds_exonic_excluded():
    # +exon (excluded), +intron (seed), intergenic (seed).
    sub, ra = _substrate_from_signatures(
        [[100, 0, 0, 0], [40, 0, 0, 0], [30, 0, 0, 0]],
        [BIT_EXON_POS, BIT_INTRON_POS, 0],
    )
    d = estimate_gdna_density(sub, ra)
    assert d.n_seed_regions == 2  # intronic + intergenic, NOT the exon
    np.testing.assert_allclose(d.seed_mass, 70.0)  # 40 + 30 (exon's 100 excluded)
    np.testing.assert_allclose(d.seed_length, 200.0)
    np.testing.assert_allclose(d.rho_0, 71.0 / 201.0)


def test_fallback_when_all_exonic():
    # No non-exonic region anywhere → half-and-half fallback over all regions.
    sub, ra = _substrate_from_signatures(
        [[10, 0, 0, 0], [30, 0, 0, 0]], [BIT_EXON_POS, BIT_EXON_POS]
    )
    d = estimate_gdna_density(sub, ra)
    assert d.fallback_used
    assert d.n_seed_regions == 0
    # α = 1 + 0.5*(10+30) = 21; β = 1 + 200 = 201.
    np.testing.assert_allclose(d.rho_0, 21.0 / 201.0)
    assert d.rho_0 > 0.0


def test_rho_0_strictly_positive_with_no_data():
    sub, ra = _substrate_from_signatures([[0, 0, 0, 0]], [0])
    d = estimate_gdna_density(sub, ra)
    assert d.rho_0 > 0.0  # unit prior floors it: 1 / (1 + 100)
    np.testing.assert_allclose(d.rho_0, 1.0 / 101.0)
