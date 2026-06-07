"""End-to-end: calibrate() recovers a planted gDNA strand overdispersion (docs/em_strand/03 §4).

Builds a payload of many intergenic (pure-gDNA) seed regions whose contained sense/antisense
counts are drawn from BetaBinom(½, od_true), runs the real calibrator, and checks the fitted
gdna_strand_overdispersion. This exercises the full Phase-2 path: substrate → count clue →
seed extraction → fit → result.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from _synthetic import make_gdna_fl_pmf

from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import CalibrationConfig
from rigel.scan_payload import AccumulatorPayload

_STRAND_MODEL = SimpleNamespace(p_r1_sense=0.95, n_observations=200)


def _intergenic_betabinom_payload(n_regions, depth, overdispersion, seed):
    """1-ref payload of `n_regions` intergenic regions; contained gDNA ~ BetaBinom(½, od)."""
    rng = np.random.default_rng(seed)
    a = 0.5 * (1.0 - overdispersion) / overdispersion if overdispersion > 0 else 1e9
    p = rng.beta(a, a, size=n_regions)
    pos = rng.binomial(depth, p)
    neg = depth - pos

    region_contained = np.zeros((n_regions, 4), dtype=np.uint32)
    region_contained[:, 0] = pos
    region_contained[:, 1] = neg

    n_b = n_regions + 1
    payload = AccumulatorPayload(
        boundaries=(np.arange(n_b, dtype=np.int64) * 100),
        ref_pos_offsets=np.array([0, n_b], dtype=np.int64),
        ref_region_offsets=np.array([0, n_regions], dtype=np.int64),
        ref_boundary_offsets=np.array([0, n_b], dtype=np.int64),
        region_contained=region_contained,
        boundary_mass_left=np.zeros((n_b, 4), dtype=np.float32),
        boundary_mass_right=np.zeros((n_b, 4), dtype=np.float32),
        boundary_flux_left=np.zeros((n_b, 4), dtype=np.uint32),
        boundary_flux_right=np.zeros((n_b, 4), dtype=np.uint32),
        n_refs=1,
    )
    starts = np.arange(n_regions, dtype=np.int64) * 100
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(n_regions, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * n_regions, dtype="string"),
            "start": starts,
            "end": starts + 100,
            "length": np.full(n_regions, 100, dtype=np.int64),
            "signature": np.zeros(n_regions, dtype=np.uint8),  # intergenic ⇒ count-observable
        }
    )
    return payload, RegionArrays.from_region_df(region_df, {"chr1": 0})


@pytest.mark.parametrize("od_true", [0.05, 0.10, 0.20])
def test_calibrate_recovers_overdispersion(od_true):
    payload, ra = _intergenic_betabinom_payload(
        n_regions=400, depth=150, overdispersion=od_true, seed=int(od_true * 1000)
    )
    result = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=_STRAND_MODEL,
        gdna_fl_pmf=make_gdna_fl_pmf(),
        config=CalibrationConfig(),
    )
    assert result.gdna_strand_overdispersion == pytest.approx(od_true, rel=0.25, abs=0.02)


def test_calibrate_binomial_gdna_floors_to_zero():
    """Non-overdispersed (50/50) gDNA → the identifiability gate floors od to 0 (Binomial)."""
    payload, ra = _intergenic_betabinom_payload(
        n_regions=400, depth=150, overdispersion=0.0, seed=7
    )
    result = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=_STRAND_MODEL,
        gdna_fl_pmf=make_gdna_fl_pmf(),
        config=CalibrationConfig(),
    )
    assert result.gdna_strand_overdispersion < 0.02
