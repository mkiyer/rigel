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
    region_unspliced_support = np.zeros(n_regions, dtype=np.uint32)
    region_unspliced_support[0] = np.uint32(n_observed)
    region_spliced_support = np.zeros(n_regions, dtype=np.uint32)
    return {
        "n_regions": n_regions,
        "region_counts": region_counts,
        "channel_mass": channel_mass,
        "signature_mass": signature_mass,
        "fl_pool_mass": fl_pool_mass,
        "fl_pool_total": fl_pool_mass.sum(axis=1),
        "region_contained_unspliced_support": region_unspliced_support,
        "region_boundary_left_unspliced_support": np.zeros(n_regions, dtype=np.uint32),
        "region_boundary_right_unspliced_support": np.zeros(n_regions, dtype=np.uint32),
        "region_contained_spliced_support": region_spliced_support,
        "region_boundary_left_spliced_support": np.zeros(n_regions, dtype=np.uint32),
        "region_boundary_right_spliced_support": np.zeros(n_regions, dtype=np.uint32),
        # Aggregate: one increment per (region, fragment, splice) — same
        # value as contained_*_support for this synthetic payload where
        # no fragment crosses a boundary.
        "region_unspliced_support": region_unspliced_support.astype(np.uint64),
        "region_spliced_support": region_spliced_support.astype(np.uint64),
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
        payload_dict["region_contained_unspliced_support"] = payload_dict[
            "region_contained_unspliced_support"
        ].astype(np.int64)

        with pytest.raises(ValueError, match="region_contained_unspliced_support"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_support_shape_mismatch_raises(self):
        payload_dict = _good_payload_dict(n_regions=3)
        payload_dict["region_contained_spliced_support"] = np.zeros(2, dtype=np.uint32)

        with pytest.raises(ValueError, match="region_contained_spliced_support"):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_support_exceeds_n_observed_raises(self):
        payload_dict = _good_payload_dict(n_observed=10)
        payload_dict["region_contained_unspliced_support"][0] = np.uint32(11)

        with pytest.raises(
            ValueError,
            match="region_contained_unspliced_support contains a value",
        ):
            CalibrationScanPayload.from_scan_dict(payload_dict)

    def test_positive_mass_zero_support_raises(self):
        # Region 1 has positive contained-unspliced fractional mass but zero
        # contained-unspliced support; the per-compartment invariant must trip.
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
            match="region_contained_unspliced_support == 0",
        ):
            CalibrationScanPayload.from_scan_dict(payload_dict, n_total=11)


class TestSupportPayload:
    """Per-region per-compartment physical support contracts on the live scan payload."""

    _SUPPORT_FIELDS = (
        "region_contained_unspliced_support",
        "region_boundary_left_unspliced_support",
        "region_boundary_right_unspliced_support",
        "region_contained_spliced_support",
        "region_boundary_left_spliced_support",
        "region_boundary_right_spliced_support",
    )
    _AGG_FIELDS = ("region_unspliced_support", "region_spliced_support")

    def test_shapes_and_dtype(self, calib_scenario):
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        n_regions = len(result.index.region_df)
        for field in self._SUPPORT_FIELDS:
            arr = getattr(payload, field)
            assert arr.shape == (n_regions,), field
            assert arr.dtype == np.uint32, field
        for field in self._AGG_FIELDS:
            arr = getattr(payload, field)
            assert arr.shape == (n_regions,), field
            assert arr.dtype == np.uint64, field

    def test_support_bounded_by_n_observed(self, calib_scenario):
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        for field in self._SUPPORT_FIELDS:
            arr = getattr(payload, field)
            if arr.size:
                assert int(arr.max()) <= payload.n_observed, field
        for field in self._AGG_FIELDS:
            arr = getattr(payload, field)
            if arr.size:
                assert int(arr.max()) <= payload.n_observed, field

    def test_positive_mass_implies_positive_support(self, calib_scenario):
        """Per-(region, compartment, splice) invariant: positive fractional
        mass in a (compartment, splice) cell iff positive support in the
        matching support cell. Live-scan version of the static check in
        ``CalibrationScanPayload.from_scan_dict``."""
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

        strands = (CHANNEL_STRAND_POS, CHANNEL_STRAND_NEG)

        def _cell_mass(compartment: int, splice: int) -> np.ndarray:
            cols = [channel_index(compartment, splice, s) for s in strands]
            return payload.region_counts[:, cols].sum(axis=1, dtype=np.float64)

        cells = (
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED,
             "region_contained_unspliced_support"),
            (COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED,
             "region_boundary_left_unspliced_support"),
            (COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED,
             "region_boundary_right_unspliced_support"),
            (COMPARTMENT_CONTAINED, SPLICE_SPLICED,
             "region_contained_spliced_support"),
            (COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED,
             "region_boundary_left_spliced_support"),
            (COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED,
             "region_boundary_right_spliced_support"),
        )
        for compartment, splice, field in cells:
            mass = _cell_mass(compartment, splice)
            support = getattr(payload, field)
            assert np.all((support > 0) | (mass == 0.0)), field
            assert np.all((support == 0) | (mass > 0.0)), field

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
        for field in self._SUPPORT_FIELDS:
            np.testing.assert_array_equal(
                getattr(one, field), getattr(four, field), err_msg=field,
            )
        for field in self._AGG_FIELDS:
            np.testing.assert_array_equal(
                getattr(one, field), getattr(four, field), err_msg=field,
            )

    def test_sorted_view_round_trips(self, calib_scenario):
        """PayloadArrays reorders the per-compartment support vectors by
        ``RegionArrays.order``; the inverse permutation must recover the
        native payload vector for each of the six fields."""
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
        for field in self._SUPPORT_FIELDS:
            sorted_arr = getattr(payload_arrays, f"{field}_sorted")
            np.testing.assert_array_equal(
                sorted_arr[inv], getattr(payload, field), err_msg=field,
            )
        for field in self._AGG_FIELDS:
            sorted_arr = getattr(payload_arrays, f"{field}_sorted")
            np.testing.assert_array_equal(
                sorted_arr[inv], getattr(payload, field), err_msg=field,
            )

    def test_compartment_sum_bounds_aggregate(self, calib_scenario):
        """Per-region invariant for spanning-fragment double counting:

            aggregate <= contained + left + right <= 2 * aggregate

        The lower bound: the aggregate counts each unique (region,
        fragment) once; the per-compartment sum counts the same fragment
        in 1 or 2 compartment cells. The upper bound: a fragment hits at
        most two compartments of a region (the only way to score in two
        is to fully span, hitting left AND right). The difference
        (sum - aggregate) is exactly the count of fragments fully
        spanning each region per splice class.
        """
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", total_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        u_sum = (
            payload.region_contained_unspliced_support.astype(np.uint64)
            + payload.region_boundary_left_unspliced_support.astype(np.uint64)
            + payload.region_boundary_right_unspliced_support.astype(np.uint64)
        )
        s_sum = (
            payload.region_contained_spliced_support.astype(np.uint64)
            + payload.region_boundary_left_spliced_support.astype(np.uint64)
            + payload.region_boundary_right_spliced_support.astype(np.uint64)
        )
        assert np.all(u_sum >= payload.region_unspliced_support)
        assert np.all(s_sum >= payload.region_spliced_support)
        assert np.all(u_sum <= 2 * payload.region_unspliced_support.astype(np.uint64))
        assert np.all(s_sum <= 2 * payload.region_spliced_support.astype(np.uint64))