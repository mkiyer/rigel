"""Compatibility coverage for the current FL-aware density surfaces.

The former ``density_global`` module was retired during the fractional
accumulator cutover. These tests keep the old surface name as a regression
anchor while asserting against the live density-observation/model APIs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration._exposure import fractional_boundary_side_exposure, l_eff_contained
from rigel.calibration.density_model import (
    FLAG_FALLBACK_USED,
    FLAG_NON_ANCHOR,
    FLAG_PRIOR_DOMINATED,
    PRIOR_FAMILY_DETERMINISTIC_ZERO,
    fit_density_evidence,
)
from rigel.calibration.density_observation import build_density_observation
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
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)
from rigel.frag_length_model import FragmentLengthModel


def _delta_fl(length: int, *, max_size: int = 1024) -> FragmentLengthModel:
    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[length] = 10_000.0
    return FragmentLengthModel.from_counts(counts, max_size=max_size)


def _region_df() -> pd.DataFrame:
    rows = [
        ("chr1", 0, 1000, pack_signature()),
        ("chr1", 1000, 2000, pack_signature(intron_pos=True)),
        ("chr1", 2000, 2500, pack_signature(exon_pos=True)),
    ]
    df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    df.index = df["region_id"].to_numpy()
    return df


def _payload(region_counts: np.ndarray, region_df: pd.DataFrame) -> CalibrationScanPayload:
    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    for row_idx, signature in enumerate(region_df["signature"].to_numpy(dtype=np.uint8)):
        signature_mass[int(signature)] += float(region_counts[row_idx].sum(dtype=np.float64))
    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    return CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
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


def _build_observation(region_counts: np.ndarray, gdna_fl: FragmentLengthModel):
    region_df = _region_df()
    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0})
    payload_arrays = PayloadArrays.from_payload(_payload(region_counts, region_df), region_arrays)
    ledger = build_region_count_ledger(payload_arrays)
    return region_arrays, build_density_observation(region_arrays, ledger, gdna_fl)


def _channel(compartment: int, strand: int) -> int:
    return channel_index(compartment, SPLICE_UNSPLICED, strand)


def test_hand_counted_observation_uses_fractional_payload_channels() -> None:
    gdna_fl = _delta_fl(200)
    region_counts = np.zeros((3, N_CHANNELS), dtype=np.float32)
    region_counts[0, _channel(COMPARTMENT_CONTAINED, CHANNEL_STRAND_POS)] = 30.0
    region_counts[0, _channel(COMPARTMENT_CONTAINED, CHANNEL_STRAND_NEG)] = 20.0
    region_counts[1, _channel(COMPARTMENT_CONTAINED, CHANNEL_STRAND_POS)] = 12.0
    region_counts[1, _channel(COMPARTMENT_BOUNDARY_LEFT, CHANNEL_STRAND_POS)] = 5.0
    region_counts[1, _channel(COMPARTMENT_BOUNDARY_RIGHT, CHANNEL_STRAND_NEG)] = 7.0
    region_counts[2, _channel(COMPARTMENT_CONTAINED, CHANNEL_STRAND_POS)] = 100.0

    region_arrays, observation = _build_observation(region_counts, gdna_fl)

    assert observation.anchor_intergenic.tolist() == [True, False, False]
    assert observation.anchor_intron.tolist() == [False, True, False]
    assert observation.is_anchor.tolist() == [True, True, False]
    np.testing.assert_allclose(observation.contained_count, [50.0, 12.0, 100.0])
    np.testing.assert_allclose(observation.boundary_left_count, [0.0, 5.0, 0.0])
    np.testing.assert_allclose(observation.boundary_right_count, [0.0, 7.0, 0.0])
    np.testing.assert_allclose(observation.boundary_count, [0.0, 12.0, 0.0])

    spans = (region_arrays.end - region_arrays.start).astype(np.int64)
    np.testing.assert_allclose(observation.contained_leff, l_eff_contained(spans, gdna_fl))
    side = fractional_boundary_side_exposure(spans, gdna_fl)
    np.testing.assert_allclose(observation.boundary_left_leff, side)
    np.testing.assert_allclose(observation.boundary_right_leff, side)


def test_density_evidence_uses_sparse_intergenic_fallback_for_small_anchor_sets() -> None:
    gdna_fl = _delta_fl(200)
    region_counts = np.zeros((3, N_CHANNELS), dtype=np.float32)
    region_counts[0, _channel(COMPARTMENT_CONTAINED, CHANNEL_STRAND_POS)] = 40.0
    region_counts[1, _channel(COMPARTMENT_CONTAINED, CHANNEL_STRAND_POS)] = 10.0
    region_counts[1, _channel(COMPARTMENT_BOUNDARY_LEFT, CHANNEL_STRAND_POS)] = 30.0

    _region_arrays, observation = _build_observation(region_counts, gdna_fl)
    evidence = fit_density_evidence(observation)

    assert evidence.rho_ref > 0.0
    assert evidence.priors["INTERGENIC"].fit_status == "sparse_intergenic_fallback"
    assert evidence.mean_unbounded[1] > evidence.mean_unbounded[0]
    assert evidence.relative_exposure[1] > evidence.relative_exposure[0]
    assert (evidence.flags[2] & FLAG_NON_ANCHOR) != 0


def test_density_evidence_zero_anchors_is_deterministic_zero() -> None:
    gdna_fl = _delta_fl(200)
    region_counts = np.zeros((3, N_CHANNELS), dtype=np.float32)
    region_counts[2, _channel(COMPARTMENT_CONTAINED, CHANNEL_STRAND_POS)] = 25.0

    _region_arrays, observation = _build_observation(region_counts, gdna_fl)
    evidence = fit_density_evidence(observation)

    assert evidence.rho_ref_source == "ZERO"
    np.testing.assert_array_equal(evidence.mean_unbounded, np.zeros(3, dtype=np.float64))
    assert np.all(evidence.prior_family == PRIOR_FAMILY_DETERMINISTIC_ZERO)
    assert np.all((evidence.flags & (FLAG_FALLBACK_USED | FLAG_PRIOR_DOMINATED)) != 0)


def test_density_evidence_validates_confidence() -> None:
    gdna_fl = _delta_fl(200)
    region_counts = np.zeros((3, N_CHANNELS), dtype=np.float32)
    _region_arrays, observation = _build_observation(region_counts, gdna_fl)

    with pytest.raises(ValueError, match="confidence"):
        fit_density_evidence(observation, confidence=1.0)