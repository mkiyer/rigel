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
    region_unspliced_support = np.zeros(n_regions, dtype=np.uint64)
    region_unspliced_support[0] = np.uint64(n_observed)
    region_spliced_support = np.zeros(n_regions, dtype=np.uint64)
    return {
        "n_regions": n_regions,
        "region_counts": region_counts,
        "channel_mass": channel_mass,
        "signature_mass": signature_mass,
        "fl_pool_mass": fl_pool_mass,
        "fl_pool_total": fl_pool_mass.sum(axis=1),
        "region_unspliced_support": region_unspliced_support,
        "region_spliced_support": region_spliced_support,
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

    def test_support_dtype_mismatch_raises(self):
        payload_dict = _good_payload_dict()
        payload_dict["region_unspliced_support"] = payload_dict[
            "region_unspliced_support"
        ].astype(np.int64)

        with pytest.raises(ValueError, match="region_unspliced_support"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_support_shape_mismatch_raises(self):
        payload_dict = _good_payload_dict(n_regions=3)
        payload_dict["region_spliced_support"] = np.zeros(2, dtype=np.uint64)

        with pytest.raises(ValueError, match="region_spliced_support"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_support_exceeds_n_observed_raises(self):
        payload_dict = _good_payload_dict(n_observed=10)
        payload_dict["region_unspliced_support"][0] = np.uint64(11)

        with pytest.raises(ValueError, match="region_unspliced_support contains a value"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_positive_mass_zero_support_raises(self):
        # Region 1 has positive unspliced fractional mass but zero support;
        # the per-region invariant must trip.
        payload_dict = _good_payload_dict(n_observed=10)
        contained_pos = channel_index(
            COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS
        )
        payload_dict["region_counts"][1, contained_pos] = np.float32(1.0)
        payload_dict["channel_mass"] = payload_dict["region_counts"].sum(
            axis=0, dtype=np.float64
        )
        payload_dict["signature_mass"][pack_signature()] += 1.0
        # Bump n_observed/total so the cap and mass-balance checks still pass.
        payload_dict["n_observed"] = 11

        with pytest.raises(
            ValueError,
            match="region_unspliced_support == 0",
        ):
            CalibrationScanPayload.from_scan_dict(payload_dict, n_total=11)


class TestSupportPayload:
    """Per-region physical support contracts on the live scan payload."""

    def test_shapes_and_dtype(self, calib_scenario):
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        n_regions = len(result.index.region_df)
        assert payload.region_unspliced_support.shape == (n_regions,)
        assert payload.region_spliced_support.shape == (n_regions,)
        assert payload.region_unspliced_support.dtype == np.uint64
        assert payload.region_spliced_support.dtype == np.uint64

    def test_support_bounded_by_n_observed(self, calib_scenario):
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        if payload.region_unspliced_support.size:
            assert int(payload.region_unspliced_support.max()) <= payload.n_observed
        if payload.region_spliced_support.size:
            assert int(payload.region_spliced_support.max()) <= payload.n_observed

    def test_positive_mass_implies_positive_support(self, calib_scenario):
        """Region-level invariant: any region with positive unspliced (resp.
        spliced) fractional mass must have positive unspliced (resp.
        spliced) support, and vice versa. This is checked inside
        ``from_scan_dict`` already; here we assert the live scan satisfies
        it without raising."""
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        from rigel.calibration.signature import (
            CHANNEL_STRAND_NEG,
            COMPARTMENT_BOUNDARY_LEFT,
            COMPARTMENT_BOUNDARY_RIGHT,
            SPLICE_SPLICED,
        )

        comps = (
            COMPARTMENT_CONTAINED,
            COMPARTMENT_BOUNDARY_LEFT,
            COMPARTMENT_BOUNDARY_RIGHT,
        )
        strands = (CHANNEL_STRAND_POS, CHANNEL_STRAND_NEG)

        def _splice_mass(splice: int) -> np.ndarray:
            cols = [channel_index(c, splice, s) for c in comps for s in strands]
            return payload.region_counts[:, cols].sum(axis=1, dtype=np.float64)

        u_mass = _splice_mass(SPLICE_UNSPLICED)
        s_mass = _splice_mass(SPLICE_SPLICED)
        assert np.all((payload.region_unspliced_support > 0) | (u_mass == 0.0))
        assert np.all((payload.region_spliced_support > 0) | (s_mass == 0.0))
        assert np.all((payload.region_unspliced_support == 0) | (u_mass > 0.0))
        assert np.all((payload.region_spliced_support == 0) | (s_mass > 0.0))

    def test_worker_merge_preserves_support(self, calib_scenario):
        _, result = calib_scenario
        payloads = []
        for n_threads in (1, 4):
            scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=n_threads)
            _stats, _sm, _flm, _buf, payload = scan_and_buffer(
                str(result.bam_path), result.index, scan_cfg,
            )
            payloads.append(payload)
        one, four = payloads
        np.testing.assert_array_equal(
            one.region_unspliced_support, four.region_unspliced_support
        )
        np.testing.assert_array_equal(
            one.region_spliced_support, four.region_spliced_support
        )

    def test_sorted_view_round_trips(self, calib_scenario):
        """PayloadArrays reorders support vectors by RegionArrays.order;
        the inverse permutation must recover the native payload vector."""
        from rigel.calibration._arrays import PayloadArrays, RegionArrays

        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        region_arrays = RegionArrays.from_region_df(
            result.index.region_df, result.index.ref_name_to_id
        )
        payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
        inv = np.argsort(region_arrays.order)
        np.testing.assert_array_equal(
            payload_arrays.region_unspliced_support_sorted[inv],
            payload.region_unspliced_support,
        )
        np.testing.assert_array_equal(
            payload_arrays.region_spliced_support_sorted[inv],
            payload.region_spliced_support,
        )