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
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_EXON,
    MASK_INTERGENIC,
    MASK_INTRON,
    MASK_N_STATES,
    CalibrationScanPayload,
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
            "ref_name": "chr1", "start": 0, "end": 500,
            "type": int(RegionType.INTERGENIC), "strand": 0,
            "tx_pos_bp": 0, "tx_neg_bp": 0,
            "exon_pos_bp": 0, "exon_neg_bp": 0,
            "boundary_flux_left": False, "boundary_flux_right": False,
        },
        {
            "ref_name": "chr1", "start": 500, "end": 1000,
            "type": int(RegionType.EXON), "strand": 1,
            "tx_pos_bp": 500, "tx_neg_bp": 0,
            "exon_pos_bp": 500, "exon_neg_bp": 0,
            "boundary_flux_left": True, "boundary_flux_right": True,
        },
        {
            "ref_name": "chr1", "start": 1000, "end": 2000,
            "type": int(RegionType.INTRON), "strand": 1,
            "tx_pos_bp": 1000, "tx_neg_bp": 0,
            "exon_pos_bp": 0, "exon_neg_bp": 0,
            "boundary_flux_left": False, "boundary_flux_right": False,
        },
    ]
    df = pd.DataFrame(rows)
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    for c in ("start", "end", "tx_pos_bp", "tx_neg_bp", "exon_pos_bp", "exon_neg_bp"):
        df[c] = df[c].astype(np.int64)
    df["type"] = df["type"].astype(np.uint8)
    df["strand"] = df["strand"].astype(np.uint8)
    df.index = df["region_id"].to_numpy()
    return df


def _stub_index(region_df: pd.DataFrame | None = None) -> SimpleNamespace:
    """Tiny TranscriptIndex stand-in carrying only ``region_df``."""
    return SimpleNamespace(region_df=region_df)


def _payload(
    *,
    n_intron: int = 200,
    n_intergenic: int = 100,
    n_exon: int = 50,
) -> CalibrationScanPayload:
    """Synthetic payload with hand-set FL histogram + per-region counts."""
    n_regions = 3
    gc = np.zeros(MASK_N_STATES, dtype=np.int64)
    gc[MASK_INTERGENIC] = n_intergenic
    gc[MASK_INTRON]     = n_intron
    gc[MASK_EXON]       = n_exon
    per_region = np.zeros((n_regions, MASK_N_STATES), dtype=np.int64)
    per_region[0, MASK_INTERGENIC] = n_intergenic
    per_region[1, MASK_EXON]       = n_exon
    per_region[2, MASK_INTRON]     = n_intron
    h = np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64)
    # gDNA pool (intron + intergenic) clusters at FL ≈ 350
    h[MASK_INTRON,     300:400] = max(1, n_intron     // 100)
    h[MASK_INTERGENIC, 300:400] = max(1, n_intergenic // 100)
    # Mirror n_observed into a benign bin so payload validator passes.
    n_obs = int(gc.sum())
    intron_by_orient = np.zeros((n_regions, 3), dtype=np.int64)
    intron_by_orient[:, 2] = per_region[:, MASK_INTRON]
    return CalibrationScanPayload(
        global_counts=gc,
        per_region_counts=per_region,
        fl_hist=h,
        u_left=np.zeros(n_regions, dtype=np.int64),
        u_right=np.zeros(n_regions, dtype=np.int64),
        intron_counts_by_orient=intron_by_orient,
        u_left_by_orient=np.zeros((n_regions, 3), dtype=np.int64),
        u_right_by_orient=np.zeros((n_regions, 3), dtype=np.int64),
        n_observed=n_obs,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_unobserved=0,
        n_unannotated_ref=0,
    )


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
        # Densities
        assert result.global_densities.gdna_fl_mean == pytest.approx(
            float(result.fl_models.gdna.mean)
        )
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
        # All three single-bit diagnostic counters are exposed.
        d = result.diagnostics
        assert d.n_intergenic_only == 100
        assert d.n_intron_only == 200
        # n_exon_only is the Diagnostics field for MASK_EXON.
        assert d.n_exon_only == 50


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
            locus=loc, n_obs=10,
            n_gdna_intergenic=1.0, n_gdna_intron=2.0,
            n_gdna_boundary_observed=0.5, n_gdna_exon_only=0.0,
            n_gdna=3.5, pi_gdna=0.35,
            rho_loco=(0.0, 0.0, 0.0),
            leff_loco=(0.0, 0.0, 0.0),
            n_eligible_boundaries=0,
            n_boundary_events=0.5,
            nrna_active=False,
            fallback_flags=0,
        )
        mlp = MultiLocusPrior(
            multi_locus_id=0, n_obs=10, n_gdna=3.5, n_rna=6.5,
            pi_gdna=0.35,
            gdna_prior_count=3.5,
            per_locus=(est,),
        )
        new_table = PriorTable(
            multi_locus_priors=(mlp,),
            gdna_prior_count=np.array([3.5], dtype=np.float64),
            enable_gdna=np.array([1], dtype=np.uint8),
        )

        updated = result.with_priors(new_table)
        assert updated is not result          # frozen contract
        assert updated.n_multi_loci == 1
        assert updated.gdna_prior_count.tolist() == [3.5]
        assert len(updated.multi_locus_prior_df) == 1
        assert len(updated.per_locus_gdna_df) == 1


class TestPriorTableEmpty:
    def test_empty_invariants(self):
        """Empty seed has zero-length arrays for every canonical field."""
        t = PriorTable.empty()
        assert t.multi_locus_priors == ()
        assert t.gdna_prior_count.dtype == np.float64
        assert len(t.gdna_prior_count) == 0
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
        weak = calibrate(**common, fl_prior_ess=10.0)
        strong = calibrate(**common, fl_prior_ess=10_000.0)
        # Strong shrinkage pulls the gdna posterior mass toward the
        # global distribution; weak shrinkage keeps it close to the raw
        # gdna pool.  The smoothed log-prob vectors must differ.
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
        assert cal.fl_models.gdna_quality == "weak"  # 300 in [200, 5000)
        assert cal.fl_models.rna_quality == "weak"   # 550 in [200, 5000)

    def test_lowering_good_threshold_promotes_to_good(self):
        cal = calibrate(**self._kwargs(), pool_quality_good=200)
        assert cal.fl_models.gdna_quality == "good"
        assert cal.fl_models.rna_quality == "good"

    def test_raising_weak_threshold_demotes_to_fallback(self):
        cal = calibrate(**self._kwargs(), pool_quality_weak=10_000)
        assert cal.fl_models.gdna_quality == "fallback"
        assert cal.fl_models.rna_quality == "fallback"
