"""
Tests for splice-junction artifact blacklist ingestion.

The blacklist is now sourced from an alignable Zarr store via
``AlignableStore.splice_blacklist()`` — a list of dicts.  These tests
exercise the records-based aggregator
(:func:`rigel.splice_blacklist.load_splice_blacklist_from_records`) with
in-memory dicts so the test suite does not require alignable to be
installed.  The Zarr-opening wrapper
:func:`load_splice_blacklist_from_zarr` is a thin call into alignable
and is covered by integration tests / manual runs.

Aggregation rules:

* Filter rows by ``count >= min_count`` (default 2).
* Group by ``(chrom, intron_start, intron_end)`` and take ``max`` of
  ``max_anchor_left`` and ``max_anchor_right`` across surviving rows.
* Strand is always ``'.'`` in alignable output and is not carried
  through to the Rigel representation.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rigel.splice_blacklist import (
    BLACKLIST_COLUMNS,
    load_splice_blacklist_from_records,
)


def _row(chrom: str, s: int, e: int, rl: int, count: int,
         al: int, ar: int) -> dict:
    return {
        "chrom": chrom, "intron_start": s, "intron_end": e,
        "strand": ".", "read_length": rl, "count": count,
        "max_anchor_left": al, "max_anchor_right": ar,
    }


class TestLoaderFromRecords:
    def test_single_row(self) -> None:
        df = load_splice_blacklist_from_records(
            [_row("chr1", 100, 500, 100, 5, 15, 20)],
            min_count=2,
        )
        assert list(df.columns) == list(BLACKLIST_COLUMNS)
        assert len(df) == 1
        row = df.iloc[0]
        assert row["ref"] == "chr1"
        assert row["start"] == 100
        assert row["end"] == 500
        assert row["max_anchor_left"] == 15
        assert row["max_anchor_right"] == 20

    def test_min_count_drops_singletons(self) -> None:
        df = load_splice_blacklist_from_records(
            [
                _row("chr1", 100, 500, 100, 1, 15, 20),  # below default 2
                _row("chr2", 200, 300, 100, 2, 5, 6),
            ],
            min_count=2,
        )
        assert len(df) == 1
        assert df.iloc[0]["ref"] == "chr2"

    def test_min_count_one_admits_all(self) -> None:
        df = load_splice_blacklist_from_records(
            [_row("chr1", 100, 500, 100, 1, 15, 20)],
            min_count=1,
        )
        assert len(df) == 1

    def test_aggregation_takes_max_across_read_lengths(self) -> None:
        df = load_splice_blacklist_from_records([
            _row("chr1", 100, 500,  50, 5,  8, 12),
            _row("chr1", 100, 500,  75, 3, 12,  9),
            _row("chr1", 100, 500, 100, 7, 15,  6),
            _row("chr1", 100, 500, 125, 4, 10, 18),
        ], min_count=2)
        assert len(df) == 1
        assert df.iloc[0]["max_anchor_left"] == 15
        assert df.iloc[0]["max_anchor_right"] == 18

    def test_aggregation_after_count_filter(self) -> None:
        """Read-length rows below min_count must NOT contribute to max."""
        df = load_splice_blacklist_from_records([
            _row("chr1", 100, 500,  50, 1, 99, 99),  # dropped
            _row("chr1", 100, 500, 100, 5, 10, 12),
        ], min_count=2)
        assert len(df) == 1
        assert df.iloc[0]["max_anchor_left"] == 10
        assert df.iloc[0]["max_anchor_right"] == 12

    def test_multiple_junctions_sorted(self) -> None:
        df = load_splice_blacklist_from_records([
            _row("chr2", 200, 300, 100, 3,  5,  5),
            _row("chr1", 100, 500, 100, 3, 10, 10),
            _row("chr1",  50,  80, 100, 3,  3,  3),
        ], min_count=2)
        assert len(df) == 3
        assert list(df["ref"]) == ["chr1", "chr1", "chr2"]
        assert list(df["start"]) == [50, 100, 200]

    def test_empty_input(self) -> None:
        df = load_splice_blacklist_from_records([], min_count=2)
        assert len(df) == 0
        assert list(df.columns) == list(BLACKLIST_COLUMNS)

    def test_all_filtered(self) -> None:
        df = load_splice_blacklist_from_records(
            [_row("chr1", 1, 2, 100, 1, 5, 5)],
            min_count=2,
        )
        assert len(df) == 0
        assert list(df.columns) == list(BLACKLIST_COLUMNS)

    def test_dtypes_narrow(self) -> None:
        df = load_splice_blacklist_from_records(
            [_row("chr1", 1, 2, 100, 5, 3, 4)],
            min_count=2,
        )
        assert df["start"].dtype == np.int32
        assert df["end"].dtype == np.int32
        assert df["max_anchor_left"].dtype == np.int32
        assert df["max_anchor_right"].dtype == np.int32

    def test_invalid_min_count(self) -> None:
        with pytest.raises(ValueError):
            load_splice_blacklist_from_records([], min_count=0)


# ----------------------------------------------------------------------
# Index round-trip via in-memory blacklist injection.
#
# We exercise the on-disk feather schema by writing a hand-crafted
# blacklist directly to ``splice_blacklist.feather`` (bypassing the
# alignable Zarr opener) and verifying that ``TranscriptIndex.load``
# wires it through to the C++ resolver.
# ----------------------------------------------------------------------


class TestIndexRoundTrip:
    def test_persisted_feather_loads_into_resolver(
        self,
        tmp_path: Path,
        mini_index_inputs,
    ) -> None:
        from rigel.index import TranscriptIndex, SJ_BLACKLIST_FEATHER

        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(
            fasta_file=str(fasta),
            gtf_file=str(gtf),
            output_dir=str(out_dir),
        )

        # Hand-write a blacklist feather with the canonical schema.
        bl = pd.DataFrame({
            "ref": ["chr1"],
            "start": np.asarray([100], dtype=np.int32),
            "end": np.asarray([200], dtype=np.int32),
            "max_anchor_left": np.asarray([10], dtype=np.int32),
            "max_anchor_right": np.asarray([15], dtype=np.int32),
        })
        bl.to_feather(out_dir / SJ_BLACKLIST_FEATHER)

        on_disk = pd.read_feather(out_dir / SJ_BLACKLIST_FEATHER)
        assert list(on_disk.columns) == list(BLACKLIST_COLUMNS)

        idx = TranscriptIndex.load(str(out_dir))
        assert idx is not None

    def test_build_without_zarr_writes_no_blacklist(
        self,
        tmp_path: Path,
        mini_index_inputs,
    ) -> None:
        from rigel.index import TranscriptIndex, SJ_BLACKLIST_FEATHER

        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(
            fasta_file=str(fasta),
            gtf_file=str(gtf),
            output_dir=str(out_dir),
        )
        assert not (out_dir / SJ_BLACKLIST_FEATHER).exists()
        idx = TranscriptIndex.load(str(out_dir))
        assert idx is not None


# ----------------------------------------------------------------------
# Shared fixture for index inputs.  Uses the same mini GTF/FASTA the
# rest of the test suite relies on.
# ----------------------------------------------------------------------


@pytest.fixture
def mini_index_inputs(tmp_path: Path):
    """Write a minimal GTF + FASTA pair suitable for ``TranscriptIndex.build``."""
    import pysam

    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "ann.gtf"

    with open(fasta, "w") as fh:
        fh.write(">chr1\n")
        fh.write("A" * 2000 + "\n")
    pysam.faidx(str(fasta))

    with open(gtf, "w") as fh:
        fh.write(
            'chr1\tsrc\texon\t1\t300\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tsrc\texon\t401\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tsrc\texon\t1\t1000\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n'
        )
    return fasta, gtf
