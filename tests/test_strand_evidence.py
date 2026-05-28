"""Tests for ``StrandGdnaEvidence`` and ``build_strand_gdna_evidence`` (PR 07 Phase 2)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.fractional_evidence import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from rigel.calibration.region_count_ledger import RegionCountLedger
from rigel.calibration.strand_deconv import StrandRegionCounts
from rigel.calibration.strand_evidence import (
    STRAND_INFO_PSEUDOCOUNT,
    StrandGdnaEvidence,
    build_strand_gdna_evidence,
    compute_global_info_scale,
)
from rigel.calibration.strand_summary import StrandSummary


def _make_ledger(
    *,
    contained_unspliced: np.ndarray,
    left_unspliced: np.ndarray,
    right_unspliced: np.ndarray,
) -> RegionCountLedger:
    """Minimal ledger carrying only the per-compartment unspliced support
    counters that ``build_strand_gdna_evidence`` reads. All other channel
    arrays are zero-filled with matching length.
    """
    n = int(contained_unspliced.size)
    z32 = np.zeros(n, dtype=np.float32)
    z32u = np.zeros(n, dtype=np.uint32)
    z64u = np.zeros(n, dtype=np.uint64)
    return RegionCountLedger(
        contained_unspliced_pos=z32, contained_unspliced_neg=z32,
        boundary_left_unspliced_pos=z32, boundary_left_unspliced_neg=z32,
        boundary_right_unspliced_pos=z32, boundary_right_unspliced_neg=z32,
        contained_spliced_pos=z32, contained_spliced_neg=z32,
        boundary_left_spliced_pos=z32, boundary_left_spliced_neg=z32,
        boundary_right_spliced_pos=z32, boundary_right_spliced_neg=z32,
        contained_unspliced_support=contained_unspliced.astype(np.uint32, copy=False),
        boundary_left_unspliced_support=left_unspliced.astype(np.uint32, copy=False),
        boundary_right_unspliced_support=right_unspliced.astype(np.uint32, copy=False),
        contained_spliced_support=z32u,
        boundary_left_spliced_support=z32u,
        boundary_right_spliced_support=z32u,
        unspliced_support=z64u,
        spliced_support=z64u,
    )


def _make_counts(k_sense, k_anti, *, p_r1_sense=0.95, eligible=None) -> StrandRegionCounts:
    k_sense = np.asarray(k_sense, dtype=np.float32)
    k_anti = np.asarray(k_anti, dtype=np.float32)
    n_total = k_sense + k_anti
    if eligible is None:
        eligible = np.ones(k_sense.size, dtype=bool)
    return StrandRegionCounts(
        k_sense=k_sense,
        k_antisense=k_anti,
        n_total=n_total,
        eligible=np.asarray(eligible, dtype=bool),
        p_r1_sense=float(p_r1_sense),
    )


class TestGlobalInfoScale:
    def test_perfect_strand_high_n_returns_near_one(self):
        s = StrandSummary(p_r1_sense=1.0, n_observations=100_000, n_same=100_000, n_opposite=0)
        g = compute_global_info_scale(s)
        assert 0.9999 < g <= 1.0

    def test_unstranded_returns_zero(self):
        s = StrandSummary(p_r1_sense=0.5, n_observations=100_000, n_same=50_000, n_opposite=50_000)
        g = compute_global_info_scale(s)
        assert g == pytest.approx(0.0, abs=1e-6)

    def test_no_observations_returns_zero(self):
        s = StrandSummary(p_r1_sense=0.95, n_observations=0)
        g = compute_global_info_scale(s)
        assert g == 0.0

    def test_monotone_in_n_observations(self):
        # Holding p_r1_sense fixed, more observations should increase g.
        prev = -1.0
        for n in [10, 100, 1_000, 10_000, 100_000]:
            ns = int(0.7 * n)
            s = StrandSummary(p_r1_sense=0.7, n_observations=n, n_same=ns, n_opposite=n - ns)
            g = compute_global_info_scale(s)
            assert g > prev
            prev = g
        assert prev < 1.0  # finite-sample never saturates


class TestStructuralApplicable:
    """``applicable`` must be strictly structural, derived from ts_class only."""

    def test_ts_pos_neg_are_applicable(self):
        ts = np.array([TS_POS, TS_NEG, TS_POS], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.zeros(3, dtype=np.uint32),
            left_unspliced=np.zeros(3, dtype=np.uint32),
            right_unspliced=np.zeros(3, dtype=np.uint32),
        )
        counts = _make_counts([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
        ev = build_strand_gdna_evidence(
            strand_counts=counts,
            strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50),
            ledger=ledger,
            ts_class=ts,
        )
        assert np.all(ev.applicable)
        assert not np.any(ev.structural_absent)

    def test_ts_none_ambig_are_not_applicable(self):
        ts = np.array([TS_NONE, TS_AMBIG, TS_POS], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.zeros(3, dtype=np.uint32),
            left_unspliced=np.zeros(3, dtype=np.uint32),
            right_unspliced=np.zeros(3, dtype=np.uint32),
        )
        counts = _make_counts([1.0, 2.0, 3.0], [0.5, 0.5, 0.5])
        ev = build_strand_gdna_evidence(
            strand_counts=counts,
            strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50),
            ledger=ledger,
            ts_class=ts,
        )
        np.testing.assert_array_equal(ev.applicable, [False, False, True])
        np.testing.assert_array_equal(ev.structural_absent, [True, True, False])

    def test_applicable_independent_of_counts(self):
        """Zero counts in TS_POS region still applicable; high counts in TS_NONE still not."""
        ts = np.array([TS_POS, TS_NONE], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.array([0, 1000], dtype=np.uint32),
            left_unspliced=np.zeros(2, dtype=np.uint32),
            right_unspliced=np.zeros(2, dtype=np.uint32),
        )
        counts = _make_counts([0.0, 500.0], [0.0, 500.0])
        ev = build_strand_gdna_evidence(
            strand_counts=counts,
            strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50),
            ledger=ledger,
            ts_class=ts,
        )
        np.testing.assert_array_equal(ev.applicable, [True, False])


class TestRegionInfoGain:
    """``g_region`` is continuous, monotone in eff_support, bounded in [0, 1]."""

    def test_zero_eff_support_yields_zero_gain(self):
        ts = np.array([TS_POS], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.array([5], dtype=np.uint32),
            left_unspliced=np.zeros(1, dtype=np.uint32),
            right_unspliced=np.zeros(1, dtype=np.uint32),
        )
        counts = _make_counts([0.0], [0.0])
        ev = build_strand_gdna_evidence(
            strand_counts=counts,
            strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50),
            ledger=ledger,
            ts_class=ts,
        )
        assert ev.region_info_gain[0] == 0.0

    def test_large_eff_support_approaches_one(self):
        ts = np.array([TS_POS], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.array([10000], dtype=np.uint32),
            left_unspliced=np.zeros(1, dtype=np.uint32),
            right_unspliced=np.zeros(1, dtype=np.uint32),
        )
        counts = _make_counts([5000.0], [5000.0])
        ev = build_strand_gdna_evidence(
            strand_counts=counts,
            strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50),
            ledger=ledger,
            ts_class=ts,
        )
        assert ev.region_info_gain[0] > 0.999

    def test_half_gain_at_pseudocount(self):
        """g_region(eff_support = α0) == 0.5."""
        ts = np.array([TS_POS], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.array([10], dtype=np.uint32),
            left_unspliced=np.zeros(1, dtype=np.uint32),
            right_unspliced=np.zeros(1, dtype=np.uint32),
        )
        # eff_support = sense + anti = α0/2 + α0/2 = α0
        half = STRAND_INFO_PSEUDOCOUNT * 0.5
        counts = _make_counts([half], [half])
        ev = build_strand_gdna_evidence(
            strand_counts=counts,
            strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50),
            ledger=ledger,
            ts_class=ts,
        )
        assert ev.region_info_gain[0] == pytest.approx(0.5, rel=1e-12)

    def test_monotone_in_eff_support(self):
        n = 7
        ts = np.full(n, TS_POS, dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.full(n, 10, dtype=np.uint32),
            left_unspliced=np.zeros(n, dtype=np.uint32),
            right_unspliced=np.zeros(n, dtype=np.uint32),
        )
        sense = np.geomspace(0.1, 1000.0, n).astype(np.float32)
        counts = _make_counts(sense, np.zeros(n, dtype=np.float32))
        ev = build_strand_gdna_evidence(
            strand_counts=counts,
            strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50),
            ledger=ledger,
            ts_class=ts,
        )
        diffs = np.diff(ev.region_info_gain)
        assert np.all(diffs > 0)


class TestSchemaInvariants:
    def test_returns_correct_dataclass_with_all_fields(self):
        ts = np.array([TS_POS, TS_NEG, TS_NONE], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.array([3, 4, 0], dtype=np.uint32),
            left_unspliced=np.array([1, 0, 0], dtype=np.uint32),
            right_unspliced=np.array([0, 2, 0], dtype=np.uint32),
        )
        counts = _make_counts([2.0, 3.0, 0.0], [0.5, 1.0, 0.0])
        s = StrandSummary(p_r1_sense=0.9, n_observations=500, n_same=450, n_opposite=50)
        ev = build_strand_gdna_evidence(
            strand_counts=counts, strand_summary=s, ledger=ledger, ts_class=ts,
        )
        assert isinstance(ev, StrandGdnaEvidence)
        assert ev.n_sense.dtype == np.float64
        assert ev.n_anti.dtype == np.float64
        assert ev.n_total.dtype == np.float64
        assert ev.support_total.dtype == np.uint64
        np.testing.assert_array_equal(ev.support_total, [4, 6, 0])
        assert ev.eff_support.dtype == np.float64
        np.testing.assert_allclose(ev.eff_support, [2.5, 4.0, 0.0])
        assert ev.kappa.dtype == np.float64
        assert np.all(ev.kappa == pytest.approx(min(s.p_r1_sense, 1 - s.p_r1_sense)))
        assert ev.applicable.dtype == bool
        assert ev.structural_absent.dtype == bool
        np.testing.assert_array_equal(ev.applicable, ~ev.structural_absent)
        assert ev.p_r1_sense == pytest.approx(0.9)
        assert 0.0 <= ev.global_info_scale <= 1.0
        assert np.all((ev.region_info_gain >= 0.0) & (ev.region_info_gain <= 1.0))

    def test_kappa_uses_folded_minor_branch(self):
        """When p_r1_sense > 0.5, kappa = 1 - p_r1_sense; symmetric below 0.5."""
        ts = np.array([TS_POS], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.array([1], dtype=np.uint32),
            left_unspliced=np.zeros(1, dtype=np.uint32),
            right_unspliced=np.zeros(1, dtype=np.uint32),
        )
        counts = _make_counts([1.0], [0.0])
        for p in [0.2, 0.3, 0.5, 0.7, 0.95]:
            s = StrandSummary(
                p_r1_sense=p, n_observations=1000,
                n_same=int(p * 1000), n_opposite=1000 - int(p * 1000),
            )
            ev = build_strand_gdna_evidence(
                strand_counts=counts, strand_summary=s, ledger=ledger, ts_class=ts,
            )
            assert ev.kappa[0] == pytest.approx(min(p, 1.0 - p))

    def test_ts_class_must_be_in_known_set(self):
        ts = np.array([99], dtype=np.int8)  # not a known TS_* code
        ledger = _make_ledger(
            contained_unspliced=np.array([1], dtype=np.uint32),
            left_unspliced=np.zeros(1, dtype=np.uint32),
            right_unspliced=np.zeros(1, dtype=np.uint32),
        )
        counts = _make_counts([1.0], [0.0])
        with pytest.raises(AssertionError, match="ts_class contains values outside"):
            build_strand_gdna_evidence(
                strand_counts=counts,
                strand_summary=StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50),
                ledger=ledger,
                ts_class=ts,
            )


class TestContinuity:
    """No threshold cliffs: outputs are continuous in counts and p_r1_sense."""

    def test_continuous_in_eff_support(self):
        ts = np.array([TS_POS], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.array([1], dtype=np.uint32),
            left_unspliced=np.zeros(1, dtype=np.uint32),
            right_unspliced=np.zeros(1, dtype=np.uint32),
        )
        s = StrandSummary(p_r1_sense=0.95, n_observations=1000, n_same=950, n_opposite=50)
        sweep = np.linspace(0.0, 20.0, 201)
        gains = np.empty(sweep.size)
        for i, x in enumerate(sweep):
            counts = _make_counts([x], [0.0])
            ev = build_strand_gdna_evidence(
                strand_counts=counts, strand_summary=s, ledger=ledger, ts_class=ts,
            )
            gains[i] = ev.region_info_gain[0]
        # First differences should be small (smooth curve).
        max_jump = float(np.max(np.abs(np.diff(gains))))
        assert max_jump < 0.05

    def test_continuous_in_p_r1_sense(self):
        ts = np.array([TS_POS], dtype=np.int8)
        ledger = _make_ledger(
            contained_unspliced=np.array([10], dtype=np.uint32),
            left_unspliced=np.zeros(1, dtype=np.uint32),
            right_unspliced=np.zeros(1, dtype=np.uint32),
        )
        counts = _make_counts([5.0], [2.0])

        def sweep_g_kappa(npts: int) -> tuple[np.ndarray, np.ndarray]:
            sweep = np.linspace(0.5, 1.0, npts)
            gs = np.empty(sweep.size)
            ks = np.empty(sweep.size)
            for i, p in enumerate(sweep):
                n_same = int(round(p * 10_000))
                s = StrandSummary(
                    p_r1_sense=float(p), n_observations=10_000,
                    n_same=n_same, n_opposite=10_000 - n_same,
                )
                ev = build_strand_gdna_evidence(
                    strand_counts=counts, strand_summary=s, ledger=ledger, ts_class=ts,
                )
                gs[i] = ev.global_info_scale
                ks[i] = ev.kappa[0]
            return sweep, gs, ks

        sweep_a, gs_a, ks_a = sweep_g_kappa(101)
        _, gs_b, _ = sweep_g_kappa(201)
        # global_info_scale rises monotonically from 0 toward 1.
        assert gs_a[0] == pytest.approx(0.0, abs=1e-6)
        assert gs_a[-1] > 0.99
        assert np.all(np.diff(gs_a) >= -1e-12)
        # kappa = 1 - p continuously.
        np.testing.assert_allclose(ks_a, 1.0 - sweep_a, atol=1e-12)
        # True continuity: doubling resolution must shrink the max jump
        # (a smooth function's first-difference magnitude scales with grid
        # spacing; a discontinuity would leave it unchanged).
        jump_a = float(np.max(np.abs(np.diff(gs_a))))
        jump_b = float(np.max(np.abs(np.diff(gs_b))))
        assert jump_b < 0.75 * jump_a
