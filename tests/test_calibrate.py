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


def test_calibrate_populates_region_gdna_and_uniform_exposure():
    result = calibrate(
        index=_index(),
        payload=_payload(),
        scan_trained=_scan_trained(),
        strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=10_000),
        rna_lower_confidence=0.97,
    )

    assert result.region_gdna.n_total.shape == (7,)
    assert result.region_gdna.p_r1_sense == pytest.approx(0.95)
    assert result.region_gdna.rna_lower_confidence == pytest.approx(0.97)
    assert result.region_gdna.kappa_d_n_seed_regions == 4
    assert result.region_gdna.kappa_d_n_exon_self_training >= 1

    assert result.region_exposure.mode == "density"
    assert result.region_exposure.A_r.shape == result.region_gdna.n_total.shape
    assert result.density_evidence is not None
    assert result.fused_region_gdna is not None
    assert result.fused_region_gdna.mean_count.shape == result.region_gdna.n_total.shape
    np.testing.assert_array_equal(
        result.region_exposure.A_r,
        result.density_evidence.relative_exposure.astype(np.float32),
    )


def test_calibrate_summary_has_strand_deconv_and_uniform_exposure_blocks():
    result = calibrate(
        index=_index(),
        payload=_payload(),
        scan_trained=_scan_trained(),
        strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=10_000),
        rna_lower_confidence=0.95,
    )

    summary = result.to_summary_dict()
    assert summary["calibration_config"] == {"rna_lower_confidence": 0.95}

    strand = summary["strand_deconv"]
    assert strand["rna_lower_confidence"] == pytest.approx(0.95)
    assert strand["p_r1_sense"] == pytest.approx(0.95)
    assert strand["kappa_d_n_seed_regions"] == 4
    assert strand["kappa_d_n_exon_self_training"] >= 1
    assert strand["kappa_d_fallback_used"] is False
    assert strand["n_regions"] == 7
    assert strand["n_regions_eligible"] > 0

    exposure = summary["region_exposure"]
    assert exposure["mode"] == "density"
    assert exposure["n_regions"] == 7
    assert "density_evidence" in summary
    assert "priors" in summary["density_evidence"]
    assert summary["fused_region_gdna"] is not None
    assert summary["fused_region_gdna"]["n_regions"] == 7
    assert summary["prior_table"] is None

    for forbidden in (
        "regional_exposure",
        "regional_weighting_stats",
        "multi_locus_prior_df",
        "per_locus_gdna_df",
    ):
        assert forbidden not in summary
    json.dumps(summary)
