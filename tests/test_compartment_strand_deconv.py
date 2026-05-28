"""Tests for v6 compartment-aware strand deconvolution."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import PayloadArrays, RegionArrays
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
from rigel.calibration.strand_deconv import (
    build_compartment_strand_counts,
    deconvolve_compartments_by_strand,
)
from rigel.calibration.strand_summary import StrandSummary


def _strand_summary(p_r1_sense: float, n_observations: int) -> StrandSummary:
    n_same = int(round(float(p_r1_sense) * int(n_observations)))
    return StrandSummary(
        p_r1_sense=float(p_r1_sense),
        n_observations=int(n_observations),
        n_same=n_same,
        n_opposite=int(n_observations) - n_same,
    )


def _payload(region_counts: np.ndarray, signatures: np.ndarray) -> CalibrationScanPayload:
    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    for row_idx, signature in enumerate(signatures):
        signature_mass[int(signature)] += float(region_counts[row_idx].sum(dtype=np.float64))
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


def _build_arrays() -> tuple[RegionArrays, PayloadArrays]:
    rows = [
        ("chr1", 0, 100, pack_signature(exon_pos=True)),
        ("chr1", 100, 200, pack_signature(exon_neg=True)),
        ("chr1", 200, 300, pack_signature()),
    ]
    region_df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    region_df["region_id"] = np.arange(len(region_df), dtype=np.int64)
    region_df.index = region_df["region_id"].to_numpy()

    rc = np.zeros((3, N_CHANNELS), dtype=np.float32)
    cup = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    cun = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    csp = channel_index(COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_POS)
    blp = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    bln = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    brp = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    brn = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)

    rc[0, cup] = 80.0
    rc[0, cun] = 20.0
    rc[0, csp] = 99.0
    rc[0, blp] = 3.0
    rc[0, bln] = 7.0
    rc[0, brp] = 5.0
    rc[0, brn] = 1.0

    rc[1, cup] = 30.0
    rc[1, cun] = 70.0
    rc[1, blp] = 2.0
    rc[1, bln] = 8.0

    rc[2, cup] = 11.0
    rc[2, cun] = 13.0

    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0})
    payload_arrays = PayloadArrays.from_payload(
        _payload(rc, region_df["signature"].to_numpy()), region_arrays
    )
    return region_arrays, payload_arrays


def test_build_compartment_strand_counts_folds_unspliced_channels_only() -> None:
    region_arrays, payload_arrays = _build_arrays()

    counts = build_compartment_strand_counts(region_arrays, payload_arrays, p_r1_sense=0.95)

    assert counts.contained_sense.tolist() == [80.0, 70.0, 0.0]
    assert counts.contained_antisense.tolist() == [20.0, 30.0, 0.0]
    assert counts.contained_total.tolist() == [100.0, 100.0, 24.0]
    assert counts.boundary_left_sense.tolist() == [3.0, 8.0, 0.0]
    assert counts.boundary_left_antisense.tolist() == [7.0, 2.0, 0.0]
    assert counts.boundary_right_sense.tolist() == [5.0, 0.0, 0.0]
    assert counts.boundary_right_antisense.tolist() == [1.0, 0.0, 0.0]
    assert counts.eligible.tolist() == [True, True, False]


def test_deconvolve_compartments_is_conservative_for_ineligible_regions() -> None:
    region_arrays, payload_arrays = _build_arrays()
    counts = build_compartment_strand_counts(region_arrays, payload_arrays, p_r1_sense=0.95)

    estimates = deconvolve_compartments_by_strand(
        counts,
        kappa_d=10.0,
        strand_summary=_strand_summary(0.95, 10_000),
    )

    assert estimates.contained_mean[2] == pytest.approx(24.0)
    assert estimates.contained_rna_lower[2] == pytest.approx(0.0)
    assert estimates.boundary_left_mean[2] == pytest.approx(0.0)
    assert estimates.kappa_d == pytest.approx(10.0)
    assert estimates.p_r1_sense == pytest.approx(0.95)


def test_compartment_deconvolution_tracks_channels_independently() -> None:
    region_arrays, payload_arrays = _build_arrays()
    counts = build_compartment_strand_counts(region_arrays, payload_arrays, p_r1_sense=1.0)

    estimates = deconvolve_compartments_by_strand(
        counts,
        kappa_d=1.0e4,
        strand_summary=_strand_summary(1.0, 10_000),
    )

    assert estimates.contained_mean[0] < estimates.boundary_left_mean[0] + 100.0
    assert estimates.boundary_right_mean[0] <= counts.boundary_right_total[0]
    assert estimates.boundary_left_mean[1] <= counts.boundary_left_total[1]
    assert estimates.contained_rna_lower[0] > estimates.boundary_right_rna_lower[0]