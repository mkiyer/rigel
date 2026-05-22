"""Fractional global gDNA density tests."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration.density_global import (
    STRAND_KAPPA_MAX,
    compute_global_densities,
    estimate_strand_balance,
    l_eff_contained,
)
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_CONTAINED,
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)
from rigel.frag_length_model import FragmentLengthModel


def _gdna_fl() -> FragmentLengthModel:
    counts = np.zeros(201, dtype=np.float64)
    counts[100] = 1_000_000.0
    return FragmentLengthModel.from_counts(counts, max_size=200)


def _region_df() -> pd.DataFrame:
    rows = [
        ("chr1", 0, 500, pack_signature()),
        ("chr1", 500, 1000, pack_signature(exon_pos=True)),
        ("chr1", 1000, 2000, pack_signature(intron_pos=True)),
        ("chr1", 2000, 2600, pack_signature(intron_pos=True, intron_neg=True)),
        ("chr1", 2600, 3000, pack_signature(intron_pos=True, exon_pos=True)),
    ]
    df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    df.index = df["region_id"].to_numpy()
    return df


def _payload() -> CalibrationScanPayload:
    region_counts = np.zeros((5, N_CHANNELS), dtype=np.float32)
    pos = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    neg = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)

    region_counts[0, pos] = 60.0
    region_counts[0, neg] = 40.0
    region_counts[1, pos] = 1000.0
    region_counts[2, pos] = 90.0
    region_counts[2, neg] = 110.0
    region_counts[3, neg] = 25.0
    region_counts[4, pos] = 999.0

    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    for row, signature in enumerate(_region_df()["signature"].to_numpy()):
        signature_mass[int(signature)] += float(region_counts[row].sum(dtype=np.float64))

    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    n_observed = int(region_counts.sum(dtype=np.float64))
    return CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        n_observed=n_observed,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=5,
    )


def _arrays() -> tuple[RegionArrays, PayloadArrays]:
    region_arrays = RegionArrays.from_region_df(_region_df(), {"chr1": 0})
    payload_arrays = PayloadArrays.from_payload(_payload(), region_arrays)
    return region_arrays, payload_arrays


def test_contained_density_uses_intergenic_and_intron_only_effective_lengths():
    region_arrays, payload_arrays = _arrays()
    gdna_fl = _gdna_fl()
    table = compute_global_densities(region_arrays, payload_arrays, gdna_fl=gdna_fl)

    spans = region_arrays.end - region_arrays.start
    eff = l_eff_contained(spans, gdna_fl)

    assert table.intergenic.n_fragments == pytest.approx(100.0)
    assert table.intergenic.n_regions_used == 1
    assert table.intergenic.eff_length_bp == pytest.approx(float(eff[0]))
    assert table.intergenic.rho == pytest.approx(100.0 / float(eff[0]))

    assert table.intron.n_fragments == pytest.approx(225.0)
    assert table.intron.n_regions_used == 2
    assert table.intron.eff_length_bp == pytest.approx(float(eff[2] + eff[3]))
    assert table.intron.rho == pytest.approx(225.0 / float(eff[2] + eff[3]))

    assert table.exon_intron.n_fragments == 0.0
    assert table.contained_n_fragments == pytest.approx(325.0)
    assert table.contained_eff_length_bp == pytest.approx(float(eff[0] + eff[2] + eff[3]))
    assert table.contained_rho == pytest.approx(325.0 / float(eff[0] + eff[2] + eff[3]))

    assert table.strand_balance is not None
    assert table.strand_balance.n_fragments == pytest.approx(325.0)
    assert table.strand_balance.n_pos == pytest.approx(150.0)
    assert table.strand_balance.n_neg == pytest.approx(175.0)


def test_strand_balance_estimates_symmetric_beta_binomial_kappa():
    estimate = estimate_strand_balance(
        np.array([60.0, 90.0]),
        np.array([40.0, 110.0]),
        np.array([True, True]),
    )

    expected_residual_sum = 200.0
    expected_binomial_sum = 75.0
    expected_max_sum = 12_500.0
    expected_kappa = (expected_max_sum - expected_residual_sum) / (
        expected_residual_sum - expected_binomial_sum
    )

    assert estimate.fallback_used is False
    assert estimate.kappa == pytest.approx(expected_kappa)
    assert estimate.overdispersion_factor == pytest.approx(
        expected_residual_sum / expected_binomial_sum
    )
    assert estimate.observed_pos_fraction == pytest.approx(0.5)
    assert estimate.alpha == pytest.approx(expected_kappa / 2.0)
    assert estimate.beta == pytest.approx(expected_kappa / 2.0)


def test_strand_balance_keeps_mean_fixed_and_falls_back_when_not_overdispersed():
    estimate = estimate_strand_balance(
        np.array([50.0, 100.0]),
        np.array([50.0, 100.0]),
        np.array([True, True]),
    )

    assert estimate.kappa == pytest.approx(STRAND_KAPPA_MAX)
    assert estimate.fallback_used is True
    assert estimate.fallback_reason == "residual variance <= binomial expectation"
    assert estimate.observed_pos_fraction == pytest.approx(0.5)
    assert estimate.to_summary_dict()["mean_fixed"] == 0.5


def test_global_density_summary_exposes_contained_and_strand_blocks():
    region_arrays, payload_arrays = _arrays()
    table = compute_global_densities(region_arrays, payload_arrays, gdna_fl=_gdna_fl())
    summary = table.to_summary_dict()

    assert summary["contained_global"]["n_fragments"] == pytest.approx(325.0)
    assert summary["intergenic"]["type"] == "INTERGENIC"
    assert summary["intron"]["type"] == "INTRON"
    assert summary["strand_balance"]["model"] == "symmetric_beta_binomial"


def test_payload_arrays_match_calibrate_call_shape():
    index = SimpleNamespace(region_df=_region_df(), ref_name_to_id={"chr1": 0})
    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(_payload(), region_arrays)

    assert region_arrays.signature.shape == (5,)
    assert payload_arrays.contained_unspliced_total.shape == (5,)
