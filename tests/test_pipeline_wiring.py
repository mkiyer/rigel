"""Tests for pipeline wiring helpers."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd

from rigel.calibration._result import CalibrationResult
from rigel.locus import Locus, MultiLocus
import rigel.pipeline as pipeline
from rigel.calibration.regions import BoundaryKind, RegionStrand, RegionType, SIGNATURE_SENTINEL
from rigel.pipeline import _wire_calibration_regions


class _ScannerSpy:
    def __init__(self) -> None:
        self.calls = []

    def set_regions(self, *args):
        self.calls.append(args)


def _region_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ref_name": pd.array(["chrA", "chrB"], dtype="string"),
            "start": np.array([0, 0], dtype=np.int64),
            "end": np.array([100, 200], dtype=np.int64),
            "signature": np.array([0x0, 0x2], dtype=np.uint8),
            "left_signature": np.array([SIGNATURE_SENTINEL, SIGNATURE_SENTINEL], dtype=np.uint8),
            "right_signature": np.array([SIGNATURE_SENTINEL, SIGNATURE_SENTINEL], dtype=np.uint8),
            "boundary_kind_left": np.array(
                [int(BoundaryKind.NONE), int(BoundaryKind.NONE)], dtype=np.uint8
            ),
            "boundary_kind_right": np.array(
                [int(BoundaryKind.NONE), int(BoundaryKind.NONE)], dtype=np.uint8
            ),
            "type": np.array([int(RegionType.INTERGENIC), int(RegionType.EXON)], dtype=np.uint8),
            "strand": np.array([int(RegionStrand.NONE), int(RegionStrand.POS)], dtype=np.uint8),
        }
    )


def test_wire_calibration_regions_uses_index_ref_ids_not_resolver_map():
    scanner = _ScannerSpy()
    index = SimpleNamespace(ref_name_to_id={"chrA": 0, "chrB": 1})

    _wire_calibration_regions(scanner, index, _region_df())

    assert len(scanner.calls) == 1
    (
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
        n_refs,
    ) = scanner.calls[0]
    np.testing.assert_array_equal(ref_ids, np.array([0, 1], dtype=np.int32))
    np.testing.assert_array_equal(starts, np.array([0, 0], dtype=np.int64))
    np.testing.assert_array_equal(ends, np.array([100, 200], dtype=np.int64))
    np.testing.assert_array_equal(signatures, np.array([0x0, 0x2], dtype=np.uint8))
    np.testing.assert_array_equal(
        left_signatures, np.array([SIGNATURE_SENTINEL, SIGNATURE_SENTINEL], dtype=np.uint8)
    )
    np.testing.assert_array_equal(
        right_signatures, np.array([SIGNATURE_SENTINEL, SIGNATURE_SENTINEL], dtype=np.uint8)
    )
    np.testing.assert_array_equal(boundary_kind_left, np.array([0, 0], dtype=np.uint8))
    np.testing.assert_array_equal(boundary_kind_right, np.array([0, 0], dtype=np.uint8))
    np.testing.assert_array_equal(type_masks, np.array([0b100, 0b001], dtype=np.uint8))
    np.testing.assert_array_equal(strands, np.array([0, 1], dtype=np.uint8))
    assert n_refs == 2


def test_wire_calibration_regions_drops_unknown_refs(caplog):
    scanner = _ScannerSpy()
    index = SimpleNamespace(ref_name_to_id={"chrA": 0})

    _wire_calibration_regions(scanner, index, _region_df())

    assert len(scanner.calls) == 1
    (
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
        n_refs,
    ) = scanner.calls[0]
    np.testing.assert_array_equal(ref_ids, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(starts, np.array([0], dtype=np.int64))
    np.testing.assert_array_equal(ends, np.array([100], dtype=np.int64))
    np.testing.assert_array_equal(signatures, np.array([0x0], dtype=np.uint8))
    np.testing.assert_array_equal(left_signatures, np.array([SIGNATURE_SENTINEL], dtype=np.uint8))
    np.testing.assert_array_equal(right_signatures, np.array([SIGNATURE_SENTINEL], dtype=np.uint8))
    np.testing.assert_array_equal(boundary_kind_left, np.array([0], dtype=np.uint8))
    np.testing.assert_array_equal(boundary_kind_right, np.array([0], dtype=np.uint8))
    np.testing.assert_array_equal(type_masks, np.array([0b100], dtype=np.uint8))
    np.testing.assert_array_equal(strands, np.array([0], dtype=np.uint8))
    assert n_refs == 1
    assert "Dropping 1 calibration regions" in caplog.text


def test_quant_from_buffer_wires_scoring_priors_partition_and_em(monkeypatch):
    calls: list[str] = []
    estimator = SimpleNamespace(
        locus_id_per_transcript=np.full(1, -1, dtype=np.int32),
        locus_results=[],
        _gdna_em_total=0.0,
    )
    em_data = SimpleNamespace(n_units=1)
    locus = MultiLocus(
        multi_locus_id=0,
        transcript_indices=np.array([0], dtype=np.int32),
        unit_indices=np.array([0], dtype=np.int32),
        gdna_span=100,
        loci=(Locus(ref="chr1", ref_id=0, start=0, end=100),),
    )
    prior_table = SimpleNamespace(
        gdna_prior_count_em=np.array([1.0], dtype=np.float64),
        alpha_rna_add=np.array([2.0], dtype=np.float64),
        gdna_eff_len=np.array([100.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
        gdna_eff_len_unweighted=np.array([100.0], dtype=np.float64),
        gdna_expected_count=np.array([1.0], dtype=np.float64),
        rna_expected_count=np.array([2.0], dtype=np.float64),
        prior_unspliced_total=np.array([3.0], dtype=np.float64),
        prior_budget_raw=np.array([3.0], dtype=np.float64),
        prior_budget=np.array([3.0], dtype=np.float64),
        prior_gdna_share_raw=np.array([1.0 / 3.0], dtype=np.float64),
        prior_gdna_share_biased=np.array([1.0 / 3.0], dtype=np.float64),
        gdna_prior_density=np.array([0.01], dtype=np.float64),
        gdna_em_exposure_weight=np.array([1.0], dtype=np.float64),
    )
    calibration = CalibrationResult(
        fl_models=SimpleNamespace(rna=SimpleNamespace(), gdna=SimpleNamespace()),
        diagnostics=SimpleNamespace(),
        region_calibration=SimpleNamespace(),
        strand_channels=None,
        background_model=SimpleNamespace(),
        boundary_local=SimpleNamespace(),
        boundary_sweep=SimpleNamespace(),
    )
    index = SimpleNamespace(region_df=pd.DataFrame({"region_id": []}), ref_name_to_id={})

    monkeypatch.setattr(
        pipeline,
        "_setup_geometry_and_estimator",
        lambda *_args, **_kwargs: (object(), estimator),
    )

    def fake_score(*_args, **_kwargs):
        calls.append("score")
        return em_data

    monkeypatch.setattr(pipeline, "_score_fragments", fake_score)

    import rigel.calibration.prior as prior_mod
    import rigel.locus as locus_mod
    import rigel.locus_partition as partition_mod

    def fake_build_multi_loci(*_args, **_kwargs):
        calls.append("build_multi_loci")
        return [locus]

    def fake_assemble_priors(*_args, **_kwargs):
        calls.append("assemble_priors")
        return prior_table

    def fake_partition_and_free(*_args, **_kwargs):
        calls.append("partition")
        return {0: object()}

    def fake_run_em(*_args, **_kwargs):
        calls.append("run_em")

    monkeypatch.setattr(locus_mod, "build_multi_loci", fake_build_multi_loci)
    monkeypatch.setattr(prior_mod, "assemble_priors", fake_assemble_priors)
    monkeypatch.setattr(partition_mod, "partition_and_free", fake_partition_and_free)
    monkeypatch.setattr(pipeline, "_run_locus_em_partitioned", fake_run_em)

    _estimator, out_calibration = pipeline.quant_from_buffer(
        buffer=SimpleNamespace(),
        index=index,
        strand_models=SimpleNamespace(),
        frag_length_models=SimpleNamespace(),
        stats=SimpleNamespace(),
        calibration=calibration,
        calibration_payload=SimpleNamespace(),
    )

    assert calls == ["score", "build_multi_loci", "assemble_priors", "partition", "run_em"]
    assert estimator.locus_id_per_transcript[0] == 0
    assert out_calibration.prior_table is prior_table
    assert out_calibration.n_multi_loci == 1
