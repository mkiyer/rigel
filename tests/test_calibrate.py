"""Phase 4 calibration wiring tests."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration import calibrate
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_CONTAINED,
    FL_POOL_INTERGENIC_CONTAINED,
    FL_POOL_INTRONIC_CONTAINED,
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)
from rigel.calibration.strand_summary import StrandSummary
from rigel.frag_length_model import FragmentLengthModels
from rigel.splice import SpliceType


def _region_df() -> pd.DataFrame:
    rows = [
        ("chr1", 0, 500, pack_signature()),
        ("chr1", 500, 1000, pack_signature()),
        ("chr1", 1000, 1500, pack_signature(intron_pos=True)),
        ("chr1", 1500, 2000, pack_signature(intron_neg=True)),
        ("chr1", 2000, 2500, pack_signature(exon_pos=True)),
        ("chr1", 2500, 3000, pack_signature(exon_pos=True)),
        ("chr1", 3000, 3500, pack_signature(exon_neg=True)),
    ]
    df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    df.index = df["region_id"].to_numpy()
    return df


def _payload() -> CalibrationScanPayload:
    region_df = _region_df()
    region_counts = np.zeros((len(region_df), N_CHANNELS), dtype=np.float32)
    pos = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    neg = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    assignments = [
        (52.0, 48.0),
        (49.0, 51.0),
        (51.0, 49.0),
        (48.0, 52.0),
        (100.0, 100.0),
        (95.0, 105.0),
        (103.0, 97.0),
    ]
    for row_idx, (pos_count, neg_count) in enumerate(assignments):
        region_counts[row_idx, pos] = pos_count
        region_counts[row_idx, neg] = neg_count

    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    for row_idx, signature in enumerate(region_df["signature"].to_numpy()):
        signature_mass[int(signature)] += float(region_counts[row_idx].sum(dtype=np.float64))

    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    fl_pool_mass[FL_POOL_INTERGENIC_CONTAINED, 100] = 400.0
    fl_pool_mass[FL_POOL_INTRONIC_CONTAINED, 100] = 400.0

    n_observed = int(region_counts.sum(dtype=np.float64))
    return CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        region_unspliced_support=np.zeros(len(region_df), dtype=np.uint64),
        region_spliced_support=np.zeros(len(region_df), dtype=np.uint64),
        n_observed=n_observed,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=len(region_df),
    )


def _scan_trained() -> FragmentLengthModels:
    models = FragmentLengthModels(max_size=200)
    models.global_model.counts[100] = 1000.0
    models.global_model._total_weight = 1000.0
    spliced = models.category_models[SpliceType.SPLICED_ANNOT]
    spliced.counts[100] = 200.0
    spliced._total_weight = 200.0
    return models


def _index():
    return SimpleNamespace(region_df=_region_df(), ref_name_to_id={"chr1": 0})


def _strand_summary() -> StrandSummary:
    return StrandSummary(p_r1_sense=0.95, n_observations=10_000, n_same=9500, n_opposite=500)


def test_calibrate_populates_region_calibration_and_diagnostics():
    result = calibrate(
        index=_index(),
        payload=_payload(),
        scan_trained=_scan_trained(),
        strand_summary=_strand_summary(),
    )

    rc = result.region_calibration
    assert rc.p_states.shape == (7, 2)
    np.testing.assert_allclose(rc.p_states.sum(axis=1), 1.0, rtol=1e-6, atol=1e-6)
    assert rc.mu_gdna.shape == (7,)
    assert rc.upper_gdna.shape == (7,)
    assert rc.region_exposure.omega.shape == (7,)
    assert np.all(rc.upper_gdna >= rc.mu_gdna)
    assert result.strand_channels is not None
    assert result.strand_channels.contained_mean.shape == (7,)
    assert result.strand_channels.p_r1_sense == pytest.approx(0.95)
    assert result.strand_channels.internal_rna_lower_ci == pytest.approx(0.95)
    assert result.background_model.seed_mask.shape == (7,)
    assert result.boundary_local.alpha_excess.shape == (7,)
    assert result.boundary_sweep.mu_sweep.shape == (7,)

    for old_field in ("region_gdna", "region_exposure", "density_evidence", "fused_region_gdna"):
        assert not hasattr(result, old_field)


def test_calibrate_summary_has_region_calibration_blocks():
    result = calibrate(
        index=_index(),
        payload=_payload(),
        scan_trained=_scan_trained(),
        strand_summary=_strand_summary(),
    )

    summary = result.to_summary_dict()

    strand = summary["strand_channels"]
    assert strand["internal_rna_lower_ci"] == pytest.approx(0.95)
    assert strand["p_r1_sense"] == pytest.approx(0.95)
    assert strand["n_regions"] == 7
    assert "contained_mean" in strand

    region_cal = summary["region_calibration"]
    assert region_cal["n_regions"] == 7
    assert region_cal["n_passes"] >= 1
    assert "state_mass" in region_cal
    assert "region_exposure" in region_cal
    assert "A_r" not in region_cal
    assert "p_unexpressed" in region_cal
    assert "p_expressed" in region_cal
    assert "p_captured" not in region_cal
    assert "gamma_r" not in region_cal
    assert summary["background_model"]["n_regions"] == 7
    assert summary["boundary_local"]["n_regions"] == 7
    assert summary["boundary_sweep"]["n_regions"] == 7
    assert summary["prior_table"] is None

    for forbidden in (
        "density_evidence",
        "region_exposure",
        "fused_region_gdna",
        "strand_deconv",
        "regional_exposure",
        "regional_weighting_stats",
        "multi_locus_prior_df",
        "per_locus_gdna_df",
    ):
        assert forbidden not in summary
    json.dumps(summary)
