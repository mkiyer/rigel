"""Tests for the M3 calibration accumulator path (C++ + Python).

Covers:

* ``BamScanner.set_regions`` binding contract
* End-to-end pipeline payload attachment
* ``CalibrationScanPayload.from_scan_dict`` validation (shape, dtype,
  balance assertion)
* Worker-merge equality (1 vs 4 workers)

See ``docs/calibration/m3_implementation_plan.md`` §5.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.scan_payload import CalibrationScanPayload
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.native import BamScanner
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import Scenario, SimConfig

SEED = 42


# ---------------------------------------------------------------------------
# Fixture: small oracle scenario for end-to-end + worker-merge tests
# ---------------------------------------------------------------------------


@pytest.fixture
def calib_scenario(tmp_path):
    """Mini oracle scenario; enough genes/reads to populate calibration."""
    sc = Scenario(
        "calib_acc_test",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "calib_acc",
    )
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50},
        ],
    )
    sim_config = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=1.0, seed=SEED,
    )
    result = sc.build_oracle(n_fragments=500, sim_config=sim_config)
    yield sc, result
    sc.cleanup()


# ---------------------------------------------------------------------------
# 1. set_regions binding contract
# ---------------------------------------------------------------------------


def _make_scanner(index):
    return BamScanner(index.resolver, "XS", True, False)


def _basic_region_arrays(index):
    """Build a minimal valid (ref_ids, starts, ends, type_masks) input."""
    resolver_map = index.resolver.get_ref_to_id()
    # One INTERGENIC region [0, 100) on every known ref.
    ref_ids = np.array(sorted(resolver_map.values()), dtype=np.int32)
    n = ref_ids.size
    starts = np.zeros(n, dtype=np.int64)
    ends = np.full(n, 100, dtype=np.int64)
    type_masks = np.full(n, 0b100, dtype=np.uint8)  # INTERGENIC bit
    return ref_ids, starts, ends, type_masks


class TestSetRegions:
    def test_basic_install(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        ri, s, e, tm = _basic_region_arrays(result.index)
        n_refs = int(ri.max()) + 1
        scanner.set_regions(ri, s, e, tm, n_refs)  # should not raise

    def test_length_mismatch_rejected(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        ri, s, e, tm = _basic_region_arrays(result.index)
        with pytest.raises(Exception):
            # ends array truncated by one
            scanner.set_regions(ri, s, e[:-1], tm, int(ri.max()) + 1)

    def test_double_set_rejected(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        ri, s, e, tm = _basic_region_arrays(result.index)
        n_refs = int(ri.max()) + 1
        scanner.set_regions(ri, s, e, tm, n_refs)
        with pytest.raises(Exception):
            scanner.set_regions(ri, s, e, tm, n_refs)


# ---------------------------------------------------------------------------
# 2. Pipeline end-to-end payload attachment + balance
# ---------------------------------------------------------------------------


class TestPipelinePayload:
    def test_scan_returns_payload(self, calib_scenario):
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", n_scan_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        assert payload is not None
        # global_counts sums to n_observed by validator construction
        assert int(payload.global_counts.sum()) == payload.n_observed
        # at least some EXON-only observations on a clean strand-1 sim
        EXON_ONLY = 0b001
        assert payload.global_counts[EXON_ONLY] > 0
        # FL histogram has shape (8, 1024) and sums to <= n_observed
        assert payload.fl_hist.shape == (8, 1024)
        assert int(payload.fl_hist.sum()) == payload.n_observed
        # per-region rows == n_regions
        n_regions = len(result.index.region_df)
        assert payload.per_region_counts.shape == (n_regions, 8)
        assert payload.u_left.shape == (n_regions,)
        assert payload.u_right.shape == (n_regions,)

    def test_run_pipeline_attaches_payload(self, calib_scenario, tmp_path):
        _, result = calib_scenario
        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(str(result.bam_path), result.index, config)
        assert pr.calibration_payload is not None
        assert pr.calibration_payload.n_observed > 0


# ---------------------------------------------------------------------------
# 3. Payload validation (shape, dtype, balance)
# ---------------------------------------------------------------------------


def _good_payload_dict(n_regions: int = 3, n_observed: int = 10) -> dict:
    """A schema-correct, internally consistent calibration dict."""
    global_counts = np.zeros(8, dtype=np.int64)
    global_counts[1] = n_observed  # all EXON_ONLY
    fl_hist = np.zeros((8, 1024), dtype=np.int64)
    fl_hist[1, 200] = n_observed
    return {
        "global_counts": global_counts,
        "per_region_counts": np.zeros((n_regions, 8), dtype=np.int64),
        "fl_hist": fl_hist,
        "u_left": np.zeros(n_regions, dtype=np.int64),
        "u_right": np.zeros(n_regions, dtype=np.int64),
        "n_observed": n_observed,
        "n_excluded_multimap": 0,
        "n_excluded_chimera": 0,
        "n_excluded_artifact": 0,
        "n_unobserved": 0,
        "n_oor": 0,
    }


class TestPayloadValidation:
    def test_good_dict_round_trips(self):
        d = _good_payload_dict()
        p = CalibrationScanPayload.from_scan_dict(d, n_total=10)
        assert p.n_observed == 10

    def test_balance_violation_raises(self):
        d = _good_payload_dict(n_observed=10)
        # Pretend total was 15 (5 unaccounted)
        with pytest.raises(ValueError, match="balance assertion failed"):
            CalibrationScanPayload.from_scan_dict(d, n_total=15)

    def test_dtype_mismatch_raises(self):
        d = _good_payload_dict()
        d["global_counts"] = d["global_counts"].astype(np.int32)
        with pytest.raises(ValueError, match="dtype"):
            CalibrationScanPayload.from_scan_dict(d)

    def test_shape_mismatch_raises(self):
        d = _good_payload_dict()
        d["fl_hist"] = np.zeros((8, 512), dtype=np.int64)
        with pytest.raises(ValueError, match="shape"):
            CalibrationScanPayload.from_scan_dict(d)

    def test_global_counts_must_sum_to_n_observed(self):
        d = _good_payload_dict(n_observed=10)
        d["n_observed"] = 5  # tampered
        with pytest.raises(ValueError, match="global_counts"):
            CalibrationScanPayload.from_scan_dict(d)

    def test_n_oor_le_n_observed(self):
        d = _good_payload_dict(n_observed=10)
        d["n_oor"] = 99
        with pytest.raises(ValueError, match="n_oor"):
            CalibrationScanPayload.from_scan_dict(d)

    def test_none_dict_raises(self):
        with pytest.raises(ValueError, match="set_regions"):
            CalibrationScanPayload.from_scan_dict(None)


# ---------------------------------------------------------------------------
# 4. Worker-merge equality (1 vs N workers must be byte-identical)
# ---------------------------------------------------------------------------


class TestWorkerMergeEquality:
    def test_one_vs_four_workers_identical(self, calib_scenario):
        _, result = calib_scenario

        def _run(n: int) -> CalibrationScanPayload:
            cfg = BamScanConfig(sj_strand_tag="auto", n_scan_threads=n)
            _, _, _, _, p = scan_and_buffer(
                str(result.bam_path), result.index, cfg,
            )
            return p

        a = _run(1)
        b = _run(4)
        assert a is not None and b is not None
        np.testing.assert_array_equal(a.global_counts, b.global_counts)
        np.testing.assert_array_equal(a.per_region_counts, b.per_region_counts)
        np.testing.assert_array_equal(a.fl_hist, b.fl_hist)
        np.testing.assert_array_equal(a.u_left, b.u_left)
        np.testing.assert_array_equal(a.u_right, b.u_right)
        assert a.n_observed == b.n_observed
        assert a.n_excluded_multimap == b.n_excluded_multimap
        assert a.n_excluded_chimera == b.n_excluded_chimera
        assert a.n_excluded_artifact == b.n_excluded_artifact
        assert a.n_unobserved == b.n_unobserved
        assert a.n_oor == b.n_oor
