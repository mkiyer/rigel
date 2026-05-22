"""Tests for ``rigel.calibration._orchestrator.calibrate`` (M8a).

These tests exercise the v6 orchestrator in isolation: synthetic
``CalibrationScanPayload`` + ``FragmentLengthModels`` + a stub index
carrying a hand-built ``region_df``.  No BAM, no real index, no real
``MultiLocus`` graph.

Per-locus priors are filled in by the pipeline (M8b) via
``CalibrationResult.with_priors``; the orchestrator produces a result
with an empty ``PriorTable``.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration import calibrate
from rigel.calibration._result import CalibrationResult
from rigel.calibration.locus_prior import (
    LocusGdnaEstimate,
    MultiLocusPrior,
    PriorTable,
)
from rigel.calibration.regions import RegionType
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
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
from rigel.frag_length_model import FragmentLengthModels
from rigel.locus import Locus
from rigel.splice import SpliceType


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def _make_region_df() -> pd.DataFrame:
    """Three regions on chr1: intergenic / exon / intron tiling [0, 2000)."""
    rows = [
        {
            "ref_name": "chr1",
            "start": 0,
            "end": 500,
            "type": int(RegionType.INTERGENIC),
            "strand": 0,
            "tx_pos_bp": 0,
            "tx_neg_bp": 0,
            "exon_pos_bp": 0,
            "exon_neg_bp": 0,
            "boundary_flux_left": False,
            "boundary_flux_right": False,
            "signature": pack_signature(),
        },
        {
            "ref_name": "chr1",
            "start": 500,
            "end": 1000,
            "type": int(RegionType.EXON),
            "strand": 1,
            "tx_pos_bp": 500,
            "tx_neg_bp": 0,
            "exon_pos_bp": 500,
            "exon_neg_bp": 0,
            "boundary_flux_left": True,
            "boundary_flux_right": True,
            "signature": pack_signature(exon_pos=True),
        },
        {
            "ref_name": "chr1",
            "start": 1000,
            "end": 2000,
            "type": int(RegionType.INTRON),
            "strand": 1,
            "tx_pos_bp": 1000,
            "tx_neg_bp": 0,
            "exon_pos_bp": 0,
            "exon_neg_bp": 0,
            "boundary_flux_left": False,
            "boundary_flux_right": False,
            "signature": pack_signature(intron_pos=True),
        },
    ]
    df = pd.DataFrame(rows)
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    for c in ("start", "end", "tx_pos_bp", "tx_neg_bp", "exon_pos_bp", "exon_neg_bp"):
        df[c] = df[c].astype(np.int64)
    df["type"] = df["type"].astype(np.uint8)
    df["strand"] = df["strand"].astype(np.uint8)
    df["signature"] = df["signature"].astype(np.uint8)
    df.index = df["region_id"].to_numpy()
    return df


def _stub_index(region_df: pd.DataFrame | None = None) -> SimpleNamespace:
    """Tiny TranscriptIndex stand-in carrying only ``region_df``."""
    return SimpleNamespace(region_df=region_df, ref_name_to_id={"chr1": 0})


def _payload(
    *,
    n_intron: int = 200,
    n_intergenic: int = 100,
    n_exon: int = 50,
) -> CalibrationScanPayload:
    """Synthetic payload with hand-set FL histogram + per-region counts."""
    n_regions = 3
    region_counts = np.zeros((n_regions, N_CHANNELS), dtype=np.float32)
    contained_pos = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    region_counts[0, contained_pos] = n_intergenic
    region_counts[1, contained_pos] = n_exon
    region_counts[2, contained_pos] = n_intron

    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    signature_mass[pack_signature()] = n_intergenic
    signature_mass[pack_signature(exon_pos=True)] = n_exon
    signature_mass[pack_signature(intron_pos=True)] = n_intron

    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    _spread_fl_mass(fl_pool_mass, FL_POOL_INTERGENIC_CONTAINED, n_intergenic)
    _spread_fl_mass(fl_pool_mass, FL_POOL_INTRONIC_CONTAINED, n_intron)

    n_obs = int(n_intergenic + n_intron + n_exon)
    return CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        n_observed=n_obs,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=n_regions,
    )


def _spread_fl_mass(
    fl_pool_mass: np.ndarray,
    pool_idx: int,
    total: int,
    *,
    lo: int = 300,
    hi: int = 400,
) -> None:
    if total <= 0:
        return
    width = hi - lo
    base, rem = divmod(int(total), width)
    fl_pool_mass[pool_idx, lo:hi] = base
    if rem:
        fl_pool_mass[pool_idx, lo : lo + rem] += 1.0


def _scan_trained(
    *,
    global_counts: tuple[int, int] = (200, 350),
    spliced_counts: tuple[int, int] = (250, 300),
) -> FragmentLengthModels:
    """Scanner-trained FL models with hand-set count vectors."""
    flm = FragmentLengthModels(max_size=1024)
    lo, hi = global_counts
    flm.global_model.counts[lo:hi] = 100.0
    flm.global_model._total_weight = float(flm.global_model.counts.sum())
    sm = flm.category_models[SpliceType.SPLICED_ANNOT]
    lo, hi = spliced_counts
    sm.counts[lo:hi] = 50.0
    sm._total_weight = float(sm.counts.sum())
    return flm


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestCalibrateHappyPath:
    def test_returns_v6_result(self):
        """Orchestrator returns a ``CalibrationResult`` with all M3–M7 pieces."""
        result = calibrate(
            index=_stub_index(_make_region_df()),
            payload=_payload(),
            scan_trained=_scan_trained(),
        )
        assert isinstance(result, CalibrationResult)
        # FL surface
        assert result.fl_models.global_.mean > 0.0
        assert result.fl_models.rna.mean > 0.0
        assert result.fl_models.gdna.mean > 0.0
        assert result.global_densities is not None
        assert result.global_densities.intergenic.n_fragments == pytest.approx(100.0)
        assert result.global_densities.intron.n_fragments == pytest.approx(200.0)
        assert result.global_densities.strand_balance is not None
        assert result.regional_exposure is not None
        # No priors yet
        assert result.n_multi_loci == 0
        assert len(result.gdna_prior_count) == 0

    def test_diagnostics_populated_from_payload(self):
        """``Diagnostics.from_payload`` is called and the result is exposed."""
        result = calibrate(
            index=_stub_index(_make_region_df()),
            payload=_payload(),
            scan_trained=_scan_trained(),
        )
        d = result.diagnostics
        assert d.mass_by_coarse_class["INTERGENIC"] == 100
        assert d.mass_by_coarse_class["INTRON"] == 200
        assert d.mass_by_coarse_class["EXON"] == 50


class TestCalibrateLegacyIndexRejection:
    def test_raises_on_missing_region_df(self):
        """Stale indexes (no ``region_df``) raise with a rebuild instruction."""
        with pytest.raises(RuntimeError, match="Rebuild the index"):
            calibrate(
                index=_stub_index(region_df=None),
                payload=_payload(),
                scan_trained=_scan_trained(),
            )


class TestCalibrateWithPriorsRoundtrip:
    def test_with_priors_swaps_table(self):
        """Empty-seed result → ``with_priors`` swap → priors visible."""
        result = calibrate(
            index=_stub_index(_make_region_df()),
            payload=_payload(),
            scan_trained=_scan_trained(),
        )

        # Build a tiny PriorTable with one MultiLocus.
        loc = Locus(ref="chr1", ref_id=0, start=500, end=1000)
        est = LocusGdnaEstimate(
            locus=loc,
            n_obs=10,
            n_gdna_intergenic=1.0,
            n_gdna_intron=2.0,
            n_gdna_boundary_observed=0.5,
            n_gdna_exon_only=0.0,
            n_gdna=3.5,
            pi_gdna=0.35,
            rho_loco=(0.0, 0.0, 0.0),
            leff_loco=(0.0, 0.0, 0.0),
            n_eligible_boundaries=0,
            n_boundary_events=0.5,
            nrna_active=False,
            fallback_flags=0,
        )
        mlp = MultiLocusPrior(
            multi_locus_id=0,
            n_obs=10,
            n_gdna=3.5,
            n_rna=6.5,
            pi_gdna=0.35,
            gdna_prior_count=3.5,
            per_locus=(est,),
        )
        new_table = PriorTable(
            multi_locus_priors=(mlp,),
            gdna_prior_count=np.array([3.5], dtype=np.float64),
            gdna_prior_count_em=np.array([3.5], dtype=np.float64),
            gdna_eff_len=np.array([1000.0], dtype=np.float64),
            enable_gdna=np.array([1], dtype=np.uint8),
        )

        updated = result.with_priors(new_table)
        assert updated is not result  # frozen contract
        assert updated.n_multi_loci == 1
        assert updated.gdna_prior_count.tolist() == [3.5]
        assert updated.gdna_prior_count_em.tolist() == [3.5]
        assert len(updated.multi_locus_prior_df) == 1
        assert len(updated.per_locus_gdna_df) == 1


class TestPriorTableEmpty:
    def test_empty_invariants(self):
        """Empty seed has zero-length arrays for every canonical field."""
        t = PriorTable.empty()
        assert t.multi_locus_priors == ()
        assert t.gdna_prior_count.dtype == np.float64
        assert len(t.gdna_prior_count) == 0
        assert t.gdna_prior_count_em.dtype == np.float64
        assert len(t.gdna_prior_count_em) == 0
        assert t.enable_gdna.dtype == np.uint8
        assert len(t.enable_gdna) == 0


class TestCalibrateFLPriorEssPropagates:
    def test_different_ess_yields_different_fl_models(self):
        """``fl_prior_ess`` actually flows to ``build_fl_models``.

        EB shrinkage changes the finalized probability vector used by
        scoring, summary statistics, and effective-length calculations.
        Compare the smoothed distributions directly.
        """
        common = dict(
            index=_stub_index(_make_region_df()),
            payload=_payload(n_intron=50, n_intergenic=20),
            scan_trained=_scan_trained(),
        )
        weak = calibrate(**common, fl_prior_ess=1.0)
        strong = calibrate(**common, fl_prior_ess=10_000.0)
        # ``prior_ess`` is now a maximum: values above the adaptive cap
        # collapse to the same evidence-derived strength, but values below
        # the cap still reduce shrinkage. The smoothed log-prob vectors must
        # differ when the caller supplies a genuinely smaller cap.
        assert not np.allclose(
            weak.fl_models.gdna._log_prob,
            strong.fl_models.gdna._log_prob,
        )


class TestNoFragmentLengthModelsMutation:
    def test_scan_trained_unchanged(self):
        """Orchestrator does NOT mutate the scanner-trained ``FragmentLengthModels``.

        Pre-M8 the v1 path mutated ``frag_length_models.gdna_model``
        with the calibrated copy (pipeline.py:778).  The v6 path must
        leave the scanner-trained instance untouched: the EM scorer
        consumes ``CalibrationResult.fl_models`` directly.
        """
        scan_trained = _scan_trained()
        before_gdna = scan_trained.gdna_model
        before_rna = scan_trained.rna_model
        calibrate(
            index=_stub_index(_make_region_df()),
            payload=_payload(),
            scan_trained=scan_trained,
        )
        # Identity preserved (no slot reassignment).
        assert scan_trained.gdna_model is before_gdna
        assert scan_trained.rna_model is before_rna


class TestCalibratePoolQualityThresholds:
    """``pool_quality_good`` / ``pool_quality_weak`` flow into ``build_fl_models``."""

    def _kwargs(self):
        return dict(
            index=_stub_index(_make_region_df()),
            payload=_payload(n_intron=200, n_intergenic=100),  # gdna n=300
            scan_trained=_scan_trained(spliced_counts=(250, 300)),  # rna n=550
        )

    def test_default_thresholds_classify_weak(self):
        cal = calibrate(**self._kwargs())
        assert cal.fl_models.gdna_quality == "weak"  # 300 in [1, 5000)
        assert cal.fl_models.rna_quality == "weak"  # 550 in [1, 5000)

    def test_lowering_good_threshold_promotes_to_good(self):
        cal = calibrate(**self._kwargs(), pool_quality_good=200)
        assert cal.fl_models.gdna_quality == "good"
        assert cal.fl_models.rna_quality == "good"

    def test_raising_weak_threshold_demotes_to_fallback(self):
        cal = calibrate(**self._kwargs(), pool_quality_weak=10_000)
        assert cal.fl_models.gdna_quality == "fallback"
        assert cal.fl_models.rna_quality == "fallback"
