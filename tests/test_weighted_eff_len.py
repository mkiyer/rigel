"""Tests for weighted_gdna_eff_len_for_loci."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import RegionArrays
from rigel.calibration._exposure import (
    footprint_exposure_weight,
    gdna_eff_len_for_loci,
    transcript_exposure_weights,
    weighted_gdna_eff_len_for_loci,
)
from rigel.calibration._regional_exposure import RegionalGdnaExposure
from rigel.calibration.regions import RegionStrand, RegionType
from rigel.frag_length_model import FragmentLengthModel
from rigel.locus import Locus


def _fl_delta(ell: int = 100, max_size: int = 200) -> FragmentLengthModel:
    fl = FragmentLengthModel(max_size=max_size)
    for _ in range(2000):
        fl.observe(ell)
    fl.finalize()
    return fl


def _three_region_arrays(span: int = 1000) -> RegionArrays:
    df = pd.DataFrame(
        {
            "ref_name": ["chr1", "chr1", "chr1"],
            "start": [0, span, 2 * span],
            "end": [span, 2 * span, 3 * span],
            "type": [int(RegionType.INTERGENIC)] * 3,
            "strand": [int(RegionStrand.NONE)] * 3,
            "boundary_flux_left": [False, False, False],
            "boundary_flux_right": [False, False, False],
        }
    )
    return RegionArrays.from_region_df(df, {"chr1": 0})


def _exposure_with_weights(region_arrays: RegionArrays, weights: np.ndarray) -> RegionalGdnaExposure:
    return RegionalGdnaExposure(
        rho_hat=np.zeros(weights.size),
        log_weight=np.log(weights),
        weight=weights.astype(np.float64),
        mode="regional",
        rho_ref=1e-5,
        n_at_floor=0,
        per_class={},
        ref_offsets=region_arrays.ref_offsets.copy(),
        ref_id=region_arrays.ref_id.copy(),
        start=region_arrays.start.copy(),
        end=region_arrays.end.copy(),
    )


def test_uniform_fast_path_bit_exact():
    region_arrays = _three_region_arrays()
    exp = RegionalGdnaExposure.uniform(region_arrays)
    fl = _fl_delta(ell=80)
    locus = Locus(ref="chr1", ref_id=0, start=500, end=1500)
    ref_lengths = {0: 100_000}
    a = gdna_eff_len_for_loci((locus,), ref_lengths, fl)
    b = weighted_gdna_eff_len_for_loci((locus,), ref_lengths, fl, exp)
    assert a == b


def test_all_unit_weights_equals_unweighted():
    region_arrays = _three_region_arrays()
    exp = _exposure_with_weights(region_arrays, np.array([1.0, 1.0, 1.0]))
    fl = _fl_delta(ell=80)
    locus = Locus(ref="chr1", ref_id=0, start=500, end=1500)
    ref_lengths = {0: 100_000}
    a = gdna_eff_len_for_loci((locus,), ref_lengths, fl)
    b = weighted_gdna_eff_len_for_loci((locus,), ref_lengths, fl, exp)
    assert b == pytest.approx(a, rel=1e-12)


def test_two_state_geometry():
    """Locus spans regions [0,1000) (weight 1.0) and [1000,2000) (weight 0.01)."""
    region_arrays = _three_region_arrays(span=1000)
    weights = np.array([1.0, 0.01, 1.0])
    exp = _exposure_with_weights(region_arrays, weights)
    # Use a tiny FL so midpoint offset is 0 — clean geometry.
    fl = FragmentLengthModel(max_size=10)
    for _ in range(1000):
        fl.observe(1)
    fl.finalize()
    locus = Locus(ref="chr1", ref_id=0, start=0, end=2000)
    ref_lengths = {0: 100_000}
    # Reference: the unweighted result is the FL-marginal count of valid
    # start positions inside [start, end). Weighted should be smaller by
    # exactly the per-region weight ratio.
    unweighted = gdna_eff_len_for_loci((locus,), ref_lengths, fl)
    result = weighted_gdna_eff_len_for_loci((locus,), ref_lengths, fl, exp)
    # For ell=1 (half=0) the midpoint window equals the start window which
    # equals [0, 2000): 1000 bp at weight 1.0 + 1000 bp at weight 0.01.
    # Allow for small leak in the PMF.
    expected_ratio = (1000 * 1.0 + 1000 * 0.01) / 2000.0
    assert result / unweighted == pytest.approx(expected_ratio, rel=0.01)


def test_min_value_floor():
    region_arrays = _three_region_arrays()
    exp = _exposure_with_weights(region_arrays, np.array([0.001, 0.001, 0.001]))
    fl = _fl_delta(ell=80)
    locus = Locus(ref="chr1", ref_id=0, start=500, end=550)  # tiny locus
    ref_lengths = {0: 100_000}
    result = weighted_gdna_eff_len_for_loci(
        (locus,), ref_lengths, fl, exp, min_value=5.0
    )
    assert result >= 5.0


def test_overlapping_intervals_no_double_count():
    """Two overlapping Locus intervals merge to one weighted window."""
    region_arrays = _three_region_arrays(span=1000)
    exp = _exposure_with_weights(region_arrays, np.array([1.0, 1.0, 1.0]))
    fl = _fl_delta(ell=80)
    ref_lengths = {0: 100_000}
    locus_a = Locus(ref="chr1", ref_id=0, start=500, end=1500)
    locus_b = Locus(ref="chr1", ref_id=0, start=1000, end=2000)
    locus_merged = Locus(ref="chr1", ref_id=0, start=500, end=2000)
    a = weighted_gdna_eff_len_for_loci((locus_a, locus_b), ref_lengths, fl, exp)
    b = weighted_gdna_eff_len_for_loci((locus_merged,), ref_lengths, fl, exp)
    assert a == pytest.approx(b, rel=1e-12)


def test_footprint_exposure_weight_merges_overlapping_blocks():
    region_arrays = _three_region_arrays(span=1000)
    exp = _exposure_with_weights(region_arrays, np.array([1.0, 0.5, 1.0]))
    weight = footprint_exposure_weight(
        [(0, 0, 1500), (0, 1000, 2000)],
        exp,
    )
    assert weight == pytest.approx(0.75)


def test_transcript_exposure_weights_use_nrna_span_and_mrna_exons():
    region_arrays = _three_region_arrays(span=1000)
    exp = _exposure_with_weights(region_arrays, np.array([1.0, 0.5, 1.0]))

    t_df = pd.DataFrame(
        {
            "is_nrna": [False, True],
            "start": np.array([0, 0], dtype=np.int64),
            "end": np.array([2000, 2000], dtype=np.int64),
        }
    )

    def build_exon_csr():
        return (
            np.array([0, 2, 3], dtype=np.int32),
            np.array([0, 1000, 0], dtype=np.int32),
            np.array([500, 1500, 2000], dtype=np.int32),
            np.array([0, 500, 0], dtype=np.int32),
        )

    index = SimpleNamespace(
        num_transcripts=2,
        t_df=t_df,
        t_to_ref_arr=np.array([0, 0], dtype=np.int32),
        build_exon_csr=build_exon_csr,
    )

    weights = transcript_exposure_weights(index, exp)
    assert weights[0] == pytest.approx(0.75)
    assert weights[1] == pytest.approx(0.75)
