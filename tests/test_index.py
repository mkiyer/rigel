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



