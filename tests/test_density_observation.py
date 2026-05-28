"""Tests for `rigel.calibration.density_observation` (v4 Phase 2 §6)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration._exposure import (
    fractional_boundary_side_exposure,
    l_eff_contained,
)
from rigel.calibration.density_observation import (
    DensityObservation,
    build_density_observation,
)
from rigel.calibration.region_count_ledger import build_region_count_ledger
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_BOUNDARY_RIGHT,
    COMPARTMENT_CONTAINED,
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)
from rigel.frag_length_model import FragmentLengthModel


def _gdna_fl(mean_fl: int = 100) -> FragmentLengthModel:
    counts = np.zeros(mean_fl * 2 + 1, dtype=np.float64)
    counts[mean_fl] = 1_000_000.0
    return FragmentLengthModel.from_counts(counts, max_size=mean_fl * 2)


def _region_df() -> pd.DataFrame:
    rows = [
        ("chr1", 0, 500, pack_signature()),  # INTERGENIC
        ("chr1", 500, 1000, pack_signature(exon_pos=True)),  # exon-only (not anchor)
        ("chr1", 1000, 2000, pack_signature(intron_pos=True)),  # INTRON
        ("chr1", 2000, 2600, pack_signature(intron_pos=True, intron_neg=True)),  # INTRON
        ("chr1", 2600, 3000, pack_signature(intron_pos=True, exon_pos=True)),  # mixed (not anchor)
    ]
    df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    df.index = df["region_id"].to_numpy()
    return df


def _payload(region_counts: np.ndarray) -> CalibrationScanPayload:
    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    for row, signature in enumerate(_region_df()["signature"].to_numpy()):
        signature_mass[int(signature)] += float(region_counts[row].sum(dtype=np.float64))
    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    return CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        region_contained_unspliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_boundary_left_unspliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_boundary_right_unspliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_contained_spliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_boundary_left_spliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_boundary_right_spliced_support=np.zeros(region_counts.shape[0], dtype=np.uint32),
        region_unspliced_support=np.zeros(region_counts.shape[0], dtype=np.uint64),
        region_spliced_support=np.zeros(region_counts.shape[0], dtype=np.uint64),
        n_observed=int(region_counts.sum(dtype=np.float64)),
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=region_counts.shape[0],
    )


def _build() -> tuple[RegionArrays, PayloadArrays, DensityObservation]:
    region_df = _region_df()
    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0})

    rc = np.zeros((5, N_CHANNELS), dtype=np.float32)
    cup = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    cun = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    blp = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    bln = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    brp = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    brn = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    csp = channel_index(COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_POS)

    # Row 0 (INTERGENIC): contained = 60+40, no boundary, no spliced
    rc[0, cup] = 60.0
    rc[0, cun] = 40.0
    # Row 2 (INTRON): contained = 90+110, left = 7+3, right = 4+6
    rc[2, cup] = 90.0
    rc[2, cun] = 110.0
    rc[2, blp] = 7.0
    rc[2, bln] = 3.0
    rc[2, brp] = 4.0
    rc[2, brn] = 6.0
    # Row 1 (exon-only): irrelevant content, just spliced diagnostic
    rc[1, csp] = 25.0

    payload = _payload(rc)
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    ledger = build_region_count_ledger(payload_arrays)
    obs = build_density_observation(region_arrays, ledger, _gdna_fl())
    return region_arrays, payload_arrays, obs


# ---- tests -----------------------------------------------------------------


def test_anchor_masks_match_signatures() -> None:
    _, _, obs = _build()
    assert obs.anchor_intergenic.tolist() == [True, False, False, False, False]
    assert obs.anchor_intron.tolist() == [False, False, True, True, False]
    assert obs.is_anchor.tolist() == [True, False, True, True, False]


def test_counts_match_ledger_unspliced_totals() -> None:
    _, _, obs = _build()
    # Row 0: 60 + 40 = 100
    assert obs.contained_count[0] == 100.0
    # Row 2: 90 + 110 = 200
    assert obs.contained_count[2] == 200.0
    # Row 2: boundary_left = 7+3 = 10, boundary_right = 4+6 = 10
    assert obs.boundary_left_count[2] == 10.0
    assert obs.boundary_right_count[2] == 10.0
    assert obs.boundary_count[2] == 20.0
    assert obs.observed_compatible_count[2] == 220.0


def test_contained_leff_matches_helper() -> None:
    region_arrays, _, obs = _build()
    spans = (region_arrays.end - region_arrays.start).astype(np.int64)
    expected = l_eff_contained(spans, _gdna_fl())
    np.testing.assert_allclose(obs.contained_leff, expected)


def test_boundary_leff_uses_per_side_helper() -> None:
    region_arrays, _, obs = _build()
    spans = (region_arrays.end - region_arrays.start).astype(np.int64)
    side = fractional_boundary_side_exposure(spans, _gdna_fl())
    np.testing.assert_allclose(obs.boundary_left_leff, side)
    np.testing.assert_allclose(obs.boundary_right_leff, side)
    np.testing.assert_allclose(obs.boundary_leff, 2.0 * side)


def test_short_region_gets_near_zero_contained_opportunity() -> None:
    # Very short region with FL mass far beyond the span -> L_c is ~zero.
    region_df = pd.DataFrame(
        {
            "ref_name": ["chr1"],
            "start": [0],
            "end": [10],
            "signature": np.array([pack_signature()], dtype=np.uint8),
        }
    )
    region_df["region_id"] = [0]
    region_df.index = region_df["region_id"].to_numpy()

    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0})
    rc = np.zeros((1, N_CHANNELS), dtype=np.float32)
    payload = CalibrationScanPayload(
        region_counts=rc,
        channel_mass=rc.sum(axis=0, dtype=np.float64),
        signature_mass=np.zeros(N_SIGNATURES, dtype=np.float64),
        fl_pool_mass=np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64),
        fl_pool_total=np.zeros(N_FL_POOLS, dtype=np.float64),
        region_contained_unspliced_support=np.zeros(1, dtype=np.uint32),
        region_boundary_left_unspliced_support=np.zeros(1, dtype=np.uint32),
        region_boundary_right_unspliced_support=np.zeros(1, dtype=np.uint32),
        region_contained_spliced_support=np.zeros(1, dtype=np.uint32),
        region_boundary_left_spliced_support=np.zeros(1, dtype=np.uint32),
        region_boundary_right_spliced_support=np.zeros(1, dtype=np.uint32),
        region_unspliced_support=np.zeros(1, dtype=np.uint64),
        region_spliced_support=np.zeros(1, dtype=np.uint64),
        n_observed=0,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=1,
    )
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    ledger = build_region_count_ledger(payload_arrays)
    obs = build_density_observation(region_arrays, ledger, _gdna_fl(mean_fl=100))

    # FL has its mass near 100 -> any span much smaller than 100 gives
    # ~zero contained opportunity and a small bounded boundary opportunity.
    assert obs.contained_leff[0] < 1.0e-3
    assert obs.boundary_left_leff[0] >= 0.0
    assert obs.boundary_left_leff[0] <= 5.0  # S / 2 with S = 10


def test_spliced_count_is_diagnostic_only() -> None:
    _, _, obs = _build()
    # Row 1 carries a contained spliced count and no unspliced counts.
    assert obs.spliced_count[1] == 25.0
    # That row is not an anchor and has zero unspliced counts.
    assert obs.contained_count[1] == 0.0
    assert obs.boundary_count[1] == 0.0
    assert not obs.is_anchor[1]


def test_observation_dtypes_are_locked() -> None:
    _, _, obs = _build()
    assert isinstance(obs, DensityObservation)
    for arr in (
        obs.contained_count,
        obs.boundary_left_count,
        obs.boundary_right_count,
        obs.boundary_count,
        obs.observed_compatible_count,
        obs.spliced_count,
    ):
        assert arr.dtype == np.float32
    for arr in (
        obs.contained_leff,
        obs.boundary_left_leff,
        obs.boundary_right_leff,
        obs.boundary_leff,
        obs.region_length,
    ):
        assert arr.dtype == np.float64
    for arr in (obs.anchor_intergenic, obs.anchor_intron, obs.is_anchor):
        assert arr.dtype == np.bool_
