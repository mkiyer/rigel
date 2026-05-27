"""Current calibration scan payload and region-install contract tests."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.regions import BoundaryKind, RegionStrand, RegionType, SIGNATURE_SENTINEL
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
    CHANNEL_STRAND_POS,
    COMPARTMENT_CONTAINED,
    FL_POOL_INTERGENIC_CONTAINED,
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.native import BamScanner
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario


SEED = 42


@pytest.fixture
def calib_scenario(tmp_path):
    """Mini oracle scenario with enough alignments to populate calibration payloads."""
    scenario = Scenario(
        "calib_acc_test",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "calib_acc",
    )
    scenario.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    scenario.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50},
        ],
    )
    sim_config = ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=SEED,
    )
    result = scenario.build_oracle(n_fragments=500, sim_config=sim_config)
    yield scenario, result
    scenario.cleanup()


def _make_scanner(index):
    return BamScanner(index.resolver, "XS", True, False)


def _basic_region_arrays(index) -> tuple[np.ndarray, ...]:
    """Build a minimal valid ``BamScanner.set_regions`` argument bundle."""
    resolver_map = index.resolver.get_ref_to_id()
    ref_ids = np.array(sorted(resolver_map.values()), dtype=np.int32)
    n_regions = ref_ids.size
    starts = np.zeros(n_regions, dtype=np.int64)
    ends = np.full(n_regions, 100, dtype=np.int64)
    signatures = np.full(n_regions, pack_signature(), dtype=np.uint8)
    left_signatures = np.full(n_regions, SIGNATURE_SENTINEL, dtype=np.uint8)
    right_signatures = np.full(n_regions, SIGNATURE_SENTINEL, dtype=np.uint8)
    boundary_kind_left = np.full(n_regions, int(BoundaryKind.NONE), dtype=np.uint8)
    boundary_kind_right = np.full(n_regions, int(BoundaryKind.NONE), dtype=np.uint8)
    type_masks = np.full(n_regions, int(RegionType.INTERGENIC), dtype=np.uint8)
    strands = np.full(n_regions, int(RegionStrand.NONE), dtype=np.uint8)
    return (
        ref_ids,
        starts,
        ends,
        signatures,
        left_signatures,
        right_signatures,
        boundary_kind_left,
        boundary_kind_right,
        type_masks,
        strands,
    )


class TestSetRegions:
    def test_basic_install(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        arrays = _basic_region_arrays(result.index)
        n_refs = int(arrays[0].max()) + 1

        scanner.set_regions(*arrays, n_refs)

    def test_legacy_install_overload_rejected(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        ref_ids, starts, ends, *_rest = _basic_region_arrays(result.index)
        n_refs = int(ref_ids.max()) + 1

        with pytest.raises(TypeError):
            scanner.set_regions(ref_ids, starts, ends, np.zeros_like(ref_ids, dtype=np.uint8), n_refs)

    def test_length_mismatch_rejected(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        arrays = list(_basic_region_arrays(result.index))
        n_refs = int(arrays[0].max()) + 1
        arrays[2] = arrays[2][:-1]

        with pytest.raises(Exception):
            scanner.set_regions(*arrays, n_refs)

    def test_double_set_rejected(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        arrays = _basic_region_arrays(result.index)
        n_refs = int(arrays[0].max()) + 1

        scanner.set_regions(*arrays, n_refs)
        with pytest.raises(Exception):
            scanner.set_regions(*arrays, n_refs)


class TestPipelinePayload:
    def test_scan_returns_fractional_payload(self, calib_scenario):
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=1)
        stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path),
            result.index,
            scan_cfg,
        )

        assert payload is not None
        assert payload.n_observed == stats.n_fragments
        assert payload.region_counts.shape == (len(result.index.region_df), N_CHANNELS)
        assert payload.channel_mass.shape == (N_CHANNELS,)
        assert payload.signature_mass.shape == (N_SIGNATURES,)
        assert payload.fl_pool_mass.shape == (N_FL_POOLS, FL_HIST_N_BINS)
        assert payload.fl_pool_total.shape == (N_FL_POOLS,)
        assert payload.region_counts.dtype == np.float32
        assert payload.channel_mass.dtype == np.float64

        region_total = float(payload.region_counts.sum(dtype=np.float64))
        assert region_total == pytest.approx(payload.channel_mass.sum())
        assert region_total == pytest.approx(payload.signature_mass.sum())
        np.testing.assert_allclose(payload.fl_pool_mass.sum(axis=1), payload.fl_pool_total)
        assert payload.fl_pool_total.sum() <= payload.n_observed

    def test_run_pipeline_attaches_payload(self, calib_scenario):
        _, result = calib_scenario
        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto", total_threads=1),
        )
        pipeline_result = run_pipeline(str(result.bam_path), result.index, config)

        assert pipeline_result.calibration_payload is not None
        assert pipeline_result.calibration_payload.n_observed > 0

    def test_worker_merge_preserves_payload_mass(self, calib_scenario):
        _, result = calib_scenario
        payloads = []
        for n_threads in (1, 4):
            scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=n_threads)
            _stats, _sm, _flm, _buf, payload = scan_and_buffer(
                str(result.bam_path),
                result.index,
                scan_cfg,
            )
            payloads.append(payload)

        one, four = payloads
        np.testing.assert_allclose(one.region_counts, four.region_counts, rtol=0.0, atol=1.0e-6)
        np.testing.assert_allclose(one.channel_mass, four.channel_mass, rtol=0.0, atol=1.0e-6)
        np.testing.assert_allclose(one.signature_mass, four.signature_mass, rtol=0.0, atol=1.0e-6)
        np.testing.assert_allclose(one.fl_pool_mass, four.fl_pool_mass, rtol=0.0, atol=1.0e-6)
        assert one.n_observed == four.n_observed


def _good_payload_dict(n_regions: int = 3, n_observed: int = 10) -> dict:
    region_counts = np.zeros((n_regions, N_CHANNELS), dtype=np.float32)
    contained_pos = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
    region_counts[0, contained_pos] = float(n_observed)
    channel_mass = region_counts.sum(axis=0, dtype=np.float64)
    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    signature_mass[pack_signature()] = float(n_observed)
    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    fl_pool_mass[FL_POOL_INTERGENIC_CONTAINED, 200] = float(n_observed)
    return {
        "n_regions": n_regions,
        "region_counts": region_counts,
        "channel_mass": channel_mass,
        "signature_mass": signature_mass,
        "fl_pool_mass": fl_pool_mass,
        "fl_pool_total": fl_pool_mass.sum(axis=1),
        "n_observed": n_observed,
        "n_excluded_multimap": 0,
        "n_excluded_chimera": 0,
        "n_excluded_artifact": 0,
        "n_excluded_strand_ambig": 0,
        "n_unobserved": 0,
        "n_unannotated_ref": 0,
        "n_fl_unavailable": 0,
        "resolver_splicing_anchor_tolerance": 0,
    }


class TestPayloadValidation:
    def test_good_dict_round_trips(self):
        payload = CalibrationScanPayload.from_scan_dict(_good_payload_dict(), n_total=10)

        assert payload.n_observed == 10
        assert payload.region_counts.shape == (3, N_CHANNELS)

    def test_balance_violation_raises(self):
        with pytest.raises(ValueError, match="balance assertion failed"):
            CalibrationScanPayload.from_scan_dict(_good_payload_dict(n_observed=10), n_total=15)

    def test_dtype_mismatch_raises(self):
        payload_dict = _good_payload_dict()
        payload_dict["region_counts"] = payload_dict["region_counts"].astype(np.float64)

        with pytest.raises(ValueError, match="dtype"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_shape_mismatch_raises(self):
        payload_dict = _good_payload_dict()
        payload_dict["fl_pool_mass"] = np.zeros((N_FL_POOLS, 512), dtype=np.float64)

        with pytest.raises(ValueError, match="shape"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_channel_mass_mismatch_raises(self):
        payload_dict = _good_payload_dict()
        payload_dict["channel_mass"][0] += 1.0

        with pytest.raises(ValueError, match="channel_mass"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_fl_pool_total_mismatch_raises(self):
        payload_dict = _good_payload_dict()
        payload_dict["fl_pool_total"][0] += 1.0

        with pytest.raises(ValueError, match="fl_pool_mass"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_none_payload_raises_helpful_error(self):
        with pytest.raises(ValueError, match="scanner returned None"):
            CalibrationScanPayload.from_scan_dict(None)