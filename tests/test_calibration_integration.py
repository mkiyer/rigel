"""Tests for density/strand calibration fusion."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import nbinom

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration.density_model import DensityEvidence, density_logpmf_grid
from rigel.calibration.density_observation import DensityObservation
from rigel.calibration.integration import (
    FUSED_BOUNDARY_FALLBACK,
    FUSED_DEGENERATE,
    FUSED_DENSITY_ONLY,
    FUSED_STRAND_USED,
    fuse_density_and_strand,
)
from rigel.calibration.region_count_ledger import build_region_count_ledger
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
from rigel.calibration.strand_deconv import build_strand_region_counts
from rigel.calibration.strand_summary import StrandSummary


def _strand_summary(p_r1_sense: float, n_observations: int) -> StrandSummary:
    if n_observations == 0:
        return StrandSummary.uninformative()
    n_same = int(round(float(p_r1_sense) * int(n_observations)))
    assert n_same / int(n_observations) == pytest.approx(float(p_r1_sense), abs=1e-12)
    return StrandSummary(
        p_r1_sense=float(p_r1_sense),
        n_observations=int(n_observations),
        n_same=n_same,
        n_opposite=int(n_observations) - n_same,
    )


def _single_region_inputs(
    signature: int,
    *,
    pos: float,
    neg: float,
    p_r1_sense: float,
) -> tuple[RegionArrays, object, DensityObservation, object]:
    df = pd.DataFrame(
        {
            "region_id": np.array([0], dtype=np.int64),
            "ref_name": pd.array(["chr1"], dtype="string"),
            "start": np.array([0], dtype=np.int64),
            "end": np.array([1000], dtype=np.int64),
            "signature": np.array([signature], dtype=np.uint8),
        }
    )
    region_arrays = RegionArrays.from_region_df(df, {"chr1": 0})
    region_counts = np.zeros((1, N_CHANNELS), dtype=np.float32)
    region_counts[0, channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)] = (
        pos
    )
    region_counts[0, channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)] = (
        neg
    )
    payload = CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=np.zeros(N_SIGNATURES, dtype=np.float64),
        fl_pool_mass=np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64),
        fl_pool_total=np.zeros(N_FL_POOLS, dtype=np.float64),
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
        n_regions=1,
    )
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    ledger = build_region_count_ledger(payload_arrays)
    observed = np.array([pos + neg], dtype=np.float32)
    observation = DensityObservation(
        contained_count=observed,
        boundary_left_count=np.zeros(1, dtype=np.float32),
        boundary_right_count=np.zeros(1, dtype=np.float32),
        boundary_count=np.zeros(1, dtype=np.float32),
        observed_compatible_count=observed,
        contained_leff=np.array([1000.0], dtype=np.float64),
        boundary_left_leff=np.zeros(1, dtype=np.float64),
        boundary_right_leff=np.zeros(1, dtype=np.float64),
        boundary_leff=np.zeros(1, dtype=np.float64),
        anchor_intergenic=np.zeros(1, dtype=bool),
        anchor_intron=np.zeros(1, dtype=bool),
        is_anchor=np.zeros(1, dtype=bool),
        spliced_count=np.zeros(1, dtype=np.float32),
        region_length=np.array([1000.0], dtype=np.float64),
    )
    strand_counts = build_strand_region_counts(
        region_arrays,
        payload_arrays,
        p_r1_sense=p_r1_sense,
    )
    return region_arrays, ledger, observation, strand_counts


def _density_evidence(
    *, alpha: float, beta: float, leff: float, confidence: float = 0.95
) -> DensityEvidence:
    p_nb = beta / (beta + leff)
    mean = alpha / beta * leff
    return DensityEvidence(
        rho_post=np.array([alpha / beta], dtype=np.float64),
        relative_exposure=np.ones(1, dtype=np.float64),
        mean_unbounded=np.array([mean], dtype=np.float64),
        upper_unbounded=np.array([nbinom.ppf(confidence, alpha, p_nb)], dtype=np.float64),
        prior_family=np.zeros(1, dtype=np.uint8),
        fallback_depth=np.zeros(1, dtype=np.uint8),
        flags=np.zeros(1, dtype=np.uint8),
        confidence=confidence,
        priors={},
        rho_ref=alpha / beta,
        rho_ref_source="TEST",
        alpha_post=np.array([alpha], dtype=np.float64),
        beta_post=np.array([beta], dtype=np.float64),
        contained_leff=np.array([leff], dtype=np.float64),
        boundary_count=np.zeros(1, dtype=np.float64),
        variance_unbounded=np.array([alpha * (1.0 - p_nb) / (p_nb * p_nb)], dtype=np.float64),
        tail_probability=np.zeros(1, dtype=np.float64),
        expected_tail_count=np.zeros(1, dtype=np.float64),
    )


def test_neutral_strand_fusion_equals_bounded_density_posterior() -> None:
    region_arrays, ledger, observation, strand_counts = _single_region_inputs(
        pack_signature(),
        pos=10.0,
        neg=0.0,
        p_r1_sense=0.5,
    )
    evidence = _density_evidence(alpha=4.0, beta=100.0, leff=1000.0)
    d_grid = np.arange(11, dtype=np.int64)
    log_density = density_logpmf_grid(evidence, 0, d_grid)
    probs = np.exp(log_density - np.max(log_density))
    probs /= probs.sum()

    fused = fuse_density_and_strand(
        region_arrays=region_arrays,
        ledger=ledger,
        density_observation=observation,
        density_evidence=evidence,
        strand_counts=strand_counts,
        strand_summary=StrandSummary.uninformative(),
        kappa_d=10.0,
        confidence=0.95,
    )

    assert (fused.flags[0] & FUSED_DENSITY_ONLY) != 0
    assert not bool(fused.strand_applicable[0])
    assert fused.strand_weight[0] == 0.0
    assert fused.mean_count[0] == pytest.approx(float(np.sum(d_grid * probs)), abs=1e-6)


def test_strand_evidence_can_overcome_high_density_prior_in_pure_rna_region() -> None:
    region_arrays, ledger, observation, strand_counts = _single_region_inputs(
        pack_signature(exon_pos=True),
        pos=100.0,
        neg=0.0,
        p_r1_sense=1.0,
    )
    evidence = _density_evidence(alpha=2.0, beta=0.025, leff=1.0)

    fused = fuse_density_and_strand(
        region_arrays=region_arrays,
        ledger=ledger,
        density_observation=observation,
        density_evidence=evidence,
        strand_counts=strand_counts,
        strand_summary=_strand_summary(1.0, 1000),
        kappa_d=1.0e6,
        confidence=0.95,
    )

    assert (fused.flags[0] & FUSED_STRAND_USED) != 0
    assert bool(fused.strand_applicable[0])
    assert fused.mean_count[0] < 5.0
    assert fused.upper_count[0] < 10.0


def test_large_count_boundary_mode_sets_boundary_fallback_without_nan() -> None:
    region_arrays, ledger, observation, strand_counts = _single_region_inputs(
        pack_signature(exon_pos=True),
        pos=250.0,
        neg=0.0,
        p_r1_sense=1.0,
    )
    evidence = _density_evidence(alpha=2.0, beta=0.01, leff=1.0)

    fused = fuse_density_and_strand(
        region_arrays=region_arrays,
        ledger=ledger,
        density_observation=observation,
        density_evidence=evidence,
        strand_counts=strand_counts,
        strand_summary=_strand_summary(1.0, 1000),
        kappa_d=1.0e6,
        confidence=0.95,
    )

    assert (fused.flags[0] & FUSED_BOUNDARY_FALLBACK) != 0
    assert np.isfinite(fused.mean_count[0])
    assert np.isfinite(fused.upper_count[0])
    assert 0.0 <= fused.mean_count[0] <= 250.0
    assert 0.0 <= fused.upper_count[0] <= 250.0


def test_large_count_deterministic_zero_density_stays_zero() -> None:
    region_arrays, ledger, observation, strand_counts = _single_region_inputs(
        pack_signature(exon_pos=True),
        pos=0.0,
        neg=2000.0,
        p_r1_sense=0.5,
    )
    evidence = DensityEvidence(
        rho_post=np.zeros(1, dtype=np.float64),
        relative_exposure=np.ones(1, dtype=np.float64),
        mean_unbounded=np.zeros(1, dtype=np.float64),
        upper_unbounded=np.zeros(1, dtype=np.float64),
        prior_family=np.zeros(1, dtype=np.uint8),
        fallback_depth=np.zeros(1, dtype=np.uint8),
        flags=np.zeros(1, dtype=np.uint8),
        confidence=0.95,
        priors={},
        rho_ref=0.0,
        rho_ref_source="ZERO",
        alpha_post=np.zeros(1, dtype=np.float64),
        beta_post=np.ones(1, dtype=np.float64),
        contained_leff=np.array([1000.0], dtype=np.float64),
        boundary_count=np.zeros(1, dtype=np.float64),
        variance_unbounded=np.zeros(1, dtype=np.float64),
        tail_probability=np.zeros(1, dtype=np.float64),
        expected_tail_count=np.zeros(1, dtype=np.float64),
    )

    fused = fuse_density_and_strand(
        region_arrays=region_arrays,
        ledger=ledger,
        density_observation=observation,
        density_evidence=evidence,
        strand_counts=strand_counts,
        strand_summary=StrandSummary.uninformative(),
        kappa_d=2.0,
        confidence=0.95,
    )

    assert fused.mean_count[0] == pytest.approx(0.0)
    assert fused.upper_count[0] == pytest.approx(0.0)
    assert (int(fused.flags[0]) & FUSED_DEGENERATE) != 0
