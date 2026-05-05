"""Tests for ``regions.feather`` persistence and load-side validation
introduced at calibration M2.

Covers:
  - Happy-path round trip: ``build_index_artifacts`` → write feather →
    ``load_regions`` → ``validate_against_ref_lengths``.
  - INDEX_FORMAT_VERSION bump and version-gate at ``TranscriptIndex.load``.
  - Three invariant violations in ``validate_against_ref_lengths``:
    unknown ref, gap between regions, stale row index.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.regions import (
    RegionStrand,
    RegionType,
    load_regions,
    validate_against_ref_lengths,
)
from rigel.index import (
    INDEX_FORMAT_VERSION,
    REGIONS_FEATHER,
    TranscriptIndex,
    build_index_artifacts,
)
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mk_tx(t_idx: int, ref: str, strand: Strand,
           exons: list[tuple[int, int]]) -> Transcript:
    return Transcript(
        ref=ref, strand=strand,
        exons=[Interval(s, e) for s, e in exons],
        t_id=f"t{t_idx}", g_id=f"g{t_idx}",
        t_index=t_idx, g_index=t_idx,
    )


# ---------------------------------------------------------------------------
# Format version
# ---------------------------------------------------------------------------

def test_index_format_version_bumped_to_3():
    """M2 bumps the on-disk format version from 2 → 3."""
    assert INDEX_FORMAT_VERSION == 3


# ---------------------------------------------------------------------------
# Round-trip
# ---------------------------------------------------------------------------

def test_round_trip_through_feather(tmp_path: Path):
    """build → write → load → validate yields a frame-equal DataFrame."""
    txs = [
        _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)]),
        _mk_tx(1, "chr2", Strand.NEG, [(50, 150)]),
    ]
    ref_lengths = {"chr1": 1000, "chr2": 500}
    _, region_df = build_index_artifacts(txs, ref_lengths)

    path = tmp_path / REGIONS_FEATHER
    region_df.to_feather(path)

    loaded = load_regions(path)
    pd.testing.assert_frame_equal(loaded, region_df)
    validate_against_ref_lengths(loaded, ref_lengths)


def test_load_regions_missing_columns_raises(tmp_path: Path):
    """A feather file lacking required columns must raise ValueError."""
    bad = pd.DataFrame({"region_id": [0], "start": [0], "end": [10]})
    path = tmp_path / "bad.feather"
    bad.to_feather(path)
    with pytest.raises(ValueError, match="missing required columns"):
        load_regions(path)


# ---------------------------------------------------------------------------
# Validation invariants
# ---------------------------------------------------------------------------

def _make_valid_region_df() -> pd.DataFrame:
    """One INTERGENIC region tiling chr1[0,1000)."""
    return pd.DataFrame({
        "region_id": np.array([0], dtype=np.int64),
        "ref_name": pd.array(["chr1"], dtype="string"),
        "start": np.array([0], dtype=np.int64),
        "end": np.array([1000], dtype=np.int64),
        "type": np.array([int(RegionType.INTERGENIC)], dtype=np.uint8),
        "strand": np.array([int(RegionStrand.NONE)], dtype=np.uint8),
        "tx_pos_bp": np.array([0], dtype=np.int64),
        "tx_neg_bp": np.array([0], dtype=np.int64),
        "exon_pos_bp": np.array([0], dtype=np.int64),
        "exon_neg_bp": np.array([0], dtype=np.int64),
        "boundary_flux_left": np.array([False], dtype=np.bool_),
        "boundary_flux_right": np.array([False], dtype=np.bool_),
    })


def test_validate_unknown_ref_raises():
    df = _make_valid_region_df()
    with pytest.raises(ValueError, match="not found in ref_lengths"):
        validate_against_ref_lengths(df, {"chr2": 500})


def test_validate_gap_between_regions_raises():
    """Two regions on chr1 with a gap between them (no contiguity)."""
    df = pd.concat([
        _make_valid_region_df().assign(end=[400]),
        _make_valid_region_df().assign(
            region_id=[1], start=[500], end=[1000],
        ),
    ], ignore_index=True)
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    with pytest.raises(ValueError, match="gap or overlap"):
        validate_against_ref_lengths(df, {"chr1": 1000})


def test_validate_does_not_cover_full_reference_raises():
    """Last region ends short of ref length."""
    df = _make_valid_region_df().assign(end=[800])
    with pytest.raises(ValueError, match="last region ends at 800"):
        validate_against_ref_lengths(df, {"chr1": 1000})


def test_validate_first_region_must_start_at_zero():
    df = _make_valid_region_df().assign(start=[100])
    with pytest.raises(ValueError, match="first region starts at 100"):
        validate_against_ref_lengths(df, {"chr1": 1000})


def test_validate_region_id_mismatch_raises():
    df = _make_valid_region_df()
    df["region_id"] = np.array([42], dtype=np.int64)
    with pytest.raises(ValueError, match="row index does not match"):
        validate_against_ref_lengths(df, {"chr1": 1000})


def test_validate_zero_length_region_raises():
    df = _make_valid_region_df().assign(start=[100], end=[100])
    with pytest.raises(ValueError, match="non-positive length"):
        validate_against_ref_lengths(df, {"chr1": 1000})


# ---------------------------------------------------------------------------
# Version gate at TranscriptIndex.load
# ---------------------------------------------------------------------------

def test_load_rejects_stale_format_version(tmp_path: Path):
    """A manifest with an old format_version must be rejected hard."""
    idx_dir = tmp_path / "idx"
    idx_dir.mkdir()
    (idx_dir / "manifest.json").write_text(json.dumps({
        "format_version": INDEX_FORMAT_VERSION - 1,
        "rigel_version": "test",
    }))
    with pytest.raises(RuntimeError, match="format_version"):
        TranscriptIndex.load(idx_dir)


def test_load_rejects_missing_manifest(tmp_path: Path):
    """An index dir with no manifest.json must be rejected."""
    idx_dir = tmp_path / "idx"
    idx_dir.mkdir()
    with pytest.raises(RuntimeError, match="no manifest.json"):
        TranscriptIndex.load(idx_dir)


# ---------------------------------------------------------------------------
# End-to-end through TranscriptIndex
# ---------------------------------------------------------------------------

def test_transcript_index_attaches_region_df(mini_index):
    """The mini fixture index must expose a validated region_df + ref_lengths."""
    assert mini_index.region_df is not None
    assert len(mini_index.region_df) > 0
    assert hasattr(mini_index, "ref_lengths")
    assert hasattr(mini_index, "ref_name_to_id")
    assert hasattr(mini_index, "ref_names")
    # ref_name_to_id is a contiguous 0..N-1 mapping.
    assert sorted(mini_index.ref_name_to_id.values()) == list(
        range(len(mini_index.ref_names))
    )
    # Validation against ref_lengths is idempotent.
    validate_against_ref_lengths(mini_index.region_df, mini_index.ref_lengths)
