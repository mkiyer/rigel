"""Tests for rigel.index validation behavior."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rigel.index import (
    TranscriptIndex,
    INTERVALS_FEATHER,
    SJ_FEATHER,
    TRANSCRIPTS_FEATHER,
    compute_tss_groups,
)


def _write_minimal_index(tmp_path: Path, t_df: pd.DataFrame) -> Path:
    idx_dir = tmp_path / "idx"
    idx_dir.mkdir(parents=True, exist_ok=True)

    t_df.to_feather(idx_dir / TRANSCRIPTS_FEATHER)

    iv_df = pd.DataFrame(
        {
            "ref": ["chr1"],
            "start": [0],
            "end": [1000],
            "strand": [0],
            "interval_type": [2],
            "t_index": [-1],
        }
    )
    iv_df.to_feather(idx_dir / INTERVALS_FEATHER)

    sj_df = pd.DataFrame(
        columns=[
            "ref",
            "start",
            "end",
            "strand",
            "interval_type",
            "t_index",
        ]
    )
    sj_df.to_feather(idx_dir / SJ_FEATHER)
    return idx_dir


def test_load_raises_on_t_index_mismatch(tmp_path: Path):
    t_df = pd.DataFrame(
        {
            "ref": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [150, 260],
            "strand": [1, 1],
            "length": [50, 60],
            "t_id": ["t1", "t2"],
            "g_id": ["g1", "g1"],
            "t_index": [0, 2],
            "g_index": [0, 0],
            "g_name": ["G1", "G1"],
            "g_type": ["pc", "pc"],
            "is_basic": [True, True],
            "is_mane": [False, False],
            "is_ccds": [False, False],
            "abundance": [None, None],
        }
    )
    idx_dir = _write_minimal_index(tmp_path, t_df)

    with pytest.raises(ValueError, match="t_index"):
        TranscriptIndex.load(idx_dir)


def test_load_raises_when_t_index_missing(tmp_path: Path):
    t_df = pd.DataFrame(
        {
            "ref": ["chr1"],
            "start": [100],
            "end": [150],
            "strand": [1],
            "length": [50],
            "t_id": ["t1"],
            "g_id": ["g1"],
            "g_index": [0],
            "g_name": ["G1"],
            "g_type": ["pc"],
            "is_basic": [True],
            "is_mane": [False],
            "is_ccds": [False],
            "abundance": [None],
        }
    )
    idx_dir = _write_minimal_index(tmp_path, t_df)

    with pytest.raises(ValueError, match="missing 't_index'"):
        TranscriptIndex.load(idx_dir)


# ======================================================================
# TSS grouping tests
# ======================================================================


def _make_tss_df(**overrides):
    """Build a minimal t_df for TSS group tests."""
    defaults = dict(
        ref=["chr1"] * 4,
        start=[100, 150, 500, 510],
        end=[200, 250, 600, 610],
        strand=[1, 1, 1, 1],
        length=[100, 100, 100, 100],
        t_id=["t0", "t1", "t2", "t3"],
        g_id=["g0", "g0", "g1", "g1"],
        t_index=[0, 1, 2, 3],
        g_index=[0, 0, 1, 1],
    )
    defaults.update(overrides)
    df = pd.DataFrame(defaults)
    df.index = df["t_index"]
    return df


class TestComputeTssGroups:
    """Tests for compute_tss_groups()."""

    def test_empty_dataframe(self):
        df = pd.DataFrame(
            columns=["ref", "start", "end", "strand", "t_index"]
        )
        result = compute_tss_groups(df, tss_window=200)
        assert len(result) == 0
        assert result.dtype == np.int32

    def test_single_transcript(self):
        df = _make_tss_df(
            ref=["chr1"], start=[100], end=[200], strand=[1],
            length=[100], t_id=["t0"], g_id=["g0"],
            t_index=[0], g_index=[0],
        )
        result = compute_tss_groups(df, tss_window=200)
        assert len(result) == 1
        assert result[0] == 0

    def test_plus_strand_within_window(self):
        """+ strand: 5' = start. Starts at 100 and 150 → within 200bp window."""
        df = _make_tss_df()
        result = compute_tss_groups(df, tss_window=200)
        # t0 (start=100) and t1 (start=150) within window → same group
        assert result[0] == result[1]
        # t2 (start=500) and t3 (start=510) within window → same group
        assert result[2] == result[3]
        # But the two clusters are different groups
        assert result[0] != result[2]

    def test_plus_strand_outside_window(self):
        """+ strand: starts 100 and 350 with window=200 → different groups."""
        df = _make_tss_df(
            start=[100, 350, 500, 510],
            end=[200, 450, 600, 610],
        )
        result = compute_tss_groups(df, tss_window=200)
        # t0 (start=100) and t1 (start=350): gap=250 > 200 → different
        assert result[0] != result[1]

    def test_minus_strand_uses_end_as_tss(self):
        """− strand: 5' = end coordinate."""
        df = _make_tss_df(
            start=[100, 200, 500, 600],
            end=[300, 310, 800, 810],
            strand=[2, 2, 2, 2],  # NEG strand
        )
        result = compute_tss_groups(df, tss_window=200)
        # t0 end=300, t1 end=310: gap=10 → same group
        assert result[0] == result[1]
        # t2 end=800, t3 end=810: gap=10 → same group
        assert result[2] == result[3]
        # Cluster 1 vs cluster 2: gap=490 > 200 → different
        assert result[0] != result[2]

    def test_different_chroms_are_separate(self):
        """Transcripts on different chromosomes are always separate groups."""
        df = _make_tss_df(
            ref=["chr1", "chr2", "chr1", "chr2"],
            start=[100, 100, 500, 500],
        )
        result = compute_tss_groups(df, tss_window=200)
        # Same position but different chroms → different groups
        assert result[0] != result[1]
        assert result[2] != result[3]

    def test_different_strands_are_separate(self):
        """Transcripts on different strands are separate groups."""
        df = _make_tss_df(
            strand=[1, 2, 1, 2],
            start=[100, 100, 500, 500],
            end=[200, 200, 600, 600],
        )
        result = compute_tss_groups(df, tss_window=200)
        # t0 (+, start=100) and t1 (−, end=200) → different groups
        assert result[0] != result[1]

    def test_exact_matching_window_zero(self):
        """tss_window=0: only exact coordinate matches group together."""
        df = _make_tss_df(
            start=[100, 100, 101, 500],
            end=[200, 200, 201, 600],
        )
        result = compute_tss_groups(df, tss_window=0)
        # t0 and t1 both start=100 → same group
        assert result[0] == result[1]
        # t2 start=101 → different group (window=0 means exact only)
        assert result[0] != result[2]

    def test_single_linkage_chaining(self):
        """Single-linkage chains: 100, 250, 400 with window=200.

        100→250: gap=150 ≤ 200 → same group
        250→400: gap=150 ≤ 200 → same group
        So all three end up in one group despite 100→400 gap of 300.
        """
        df = _make_tss_df(
            ref=["chr1"] * 3,
            start=[100, 250, 400],
            end=[200, 350, 500],
            strand=[1, 1, 1],
            length=[100, 100, 100],
            t_id=["t0", "t1", "t2"],
            g_id=["g0", "g0", "g0"],
            t_index=[0, 1, 2],
            g_index=[0, 0, 0],
        )
        result = compute_tss_groups(df, tss_window=200)
        # All three chained into one group
        assert result[0] == result[1] == result[2]

    def test_group_ids_are_sequential(self):
        """Group IDs are sequential starting from 0."""
        df = _make_tss_df()
        result = compute_tss_groups(df, tss_window=200)
        unique_groups = np.unique(result)
        assert unique_groups[0] == 0
        assert len(unique_groups) == 2
        assert np.array_equal(unique_groups, [0, 1])
