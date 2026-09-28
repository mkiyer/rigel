"""Ingesting the splice-junction artifact blacklist, and the aggregation rules it applies.

Rows below ``min_count`` are dropped; the survivors are grouped by ``(chrom, intron_start,
intron_end)`` and take the ``max`` of the two anchor lengths, so a junction seen at several read
lengths keeps its longest anchor; strand is not carried through, because the source always reports
it as unknown. The gates run the shipped aggregation (`aggregate_splice_blacklist`, which the store
loader calls) over rows built in memory, so the suite does not require the store to be installed.
The index gates then check that an index uses a blacklist only when its manifest records the store
it came from: a build from a store loads it into the resolver, a build or rebuild without one leaves
none on disk, `load` ignores one dropped into an index built without a store or whose manifest has no
`sources`, and an index whose recorded blacklist was removed loads with detection off.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from _index_builder import rebuild_with_splice_blacklist
from rigel.index import MANIFEST_JSON, SJ_BLACKLIST_FEATHER, SJ_BLACKLIST_TSV, TranscriptIndex
from rigel.splice_blacklist import BLACKLIST_COLUMNS, aggregate_splice_blacklist


def _row(chrom: str, s: int, e: int, rl: int, count: int, al: int, ar: int) -> dict:
    return {
        "chrom": chrom,
        "intron_start": s,
        "intron_end": e,
        "strand": ".",
        "read_length": rl,
        "count": count,
        "max_anchor_left": al,
        "max_anchor_right": ar,
    }


def load_splice_blacklist_from_records(records, *, min_count: int = 2) -> pd.DataFrame:
    """The shipped aggregation over a DataFrame of alignable rows, built here from dicts."""
    return aggregate_splice_blacklist(pd.DataFrame(list(records)), min_count=min_count)


class TestAggregation:
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
        df = load_splice_blacklist_from_records(
            [
                _row("chr1", 100, 500, 50, 5, 8, 12),
                _row("chr1", 100, 500, 75, 3, 12, 9),
                _row("chr1", 100, 500, 100, 7, 15, 6),
                _row("chr1", 100, 500, 125, 4, 10, 18),
            ],
            min_count=2,
        )
        assert len(df) == 1
        assert df.iloc[0]["max_anchor_left"] == 15
        assert df.iloc[0]["max_anchor_right"] == 18

    def test_aggregation_after_count_filter(self) -> None:
        """Read-length rows below min_count must NOT contribute to max."""
        df = load_splice_blacklist_from_records(
            [
                _row("chr1", 100, 500, 50, 1, 99, 99),  # dropped
                _row("chr1", 100, 500, 100, 5, 10, 12),
            ],
            min_count=2,
        )
        assert len(df) == 1
        assert df.iloc[0]["max_anchor_left"] == 10
        assert df.iloc[0]["max_anchor_right"] == 12

    def test_multiple_sj_sorted(self) -> None:
        df = load_splice_blacklist_from_records(
            [
                _row("chr2", 200, 300, 100, 3, 5, 5),
                _row("chr1", 100, 500, 100, 3, 10, 10),
                _row("chr1", 50, 80, 100, 3, 3, 3),
            ],
            min_count=2,
        )
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
# The blacklist in an index: used exactly when the manifest records its store.
# ----------------------------------------------------------------------

_BLACKLIST = pd.DataFrame(
    {
        "ref": ["chr1", "chr1"],
        "start": np.asarray([100, 300], dtype=np.int32),
        "end": np.asarray([200, 400], dtype=np.int32),
        "max_anchor_left": np.asarray([10, 12], dtype=np.int32),
        "max_anchor_right": np.asarray([15, 9], dtype=np.int32),
    }
)


class TestIndexRoundTrip:
    def test_a_blacklist_built_from_a_store_loads_into_the_resolver(
        self, tmp_path: Path, mini_index_inputs
    ) -> None:
        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta, gtf, out_dir)
        rebuild_with_splice_blacklist(out_dir, _BLACKLIST)

        on_disk = pd.read_feather(out_dir / SJ_BLACKLIST_FEATHER)
        assert list(on_disk.columns) == list(BLACKLIST_COLUMNS)
        assert TranscriptIndex.load(out_dir).sj_blacklist_size == len(_BLACKLIST)

    def test_build_without_a_store_writes_no_blacklist(
        self, tmp_path: Path, mini_index_inputs
    ) -> None:
        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta, gtf, out_dir)

        assert not (out_dir / SJ_BLACKLIST_FEATHER).exists()
        assert TranscriptIndex.load(out_dir).sj_blacklist_size == 0

    def test_a_rebuild_without_a_store_applies_no_blacklist(
        self, tmp_path: Path, mini_index_inputs
    ) -> None:
        """Built with a blacklist, then rebuilt in place without a store: the old one is not applied."""
        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta, gtf, out_dir)
        rebuild_with_splice_blacklist(out_dir, _BLACKLIST)
        TranscriptIndex.build(fasta, gtf, out_dir)

        assert TranscriptIndex.load(out_dir).sj_blacklist_size == 0

    def test_a_rebuild_without_a_store_removes_the_old_blacklist(
        self, tmp_path: Path, mini_index_inputs
    ) -> None:
        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta, gtf, out_dir)
        rebuild_with_splice_blacklist(out_dir, _BLACKLIST)
        assert (out_dir / SJ_BLACKLIST_FEATHER).exists()
        assert (out_dir / SJ_BLACKLIST_TSV).exists()

        TranscriptIndex.build(fasta, gtf, out_dir)
        assert not (out_dir / SJ_BLACKLIST_FEATHER).exists()
        assert not (out_dir / SJ_BLACKLIST_TSV).exists()

    def test_load_ignores_a_blacklist_in_an_index_built_without_a_store(
        self, tmp_path: Path, mini_index_inputs
    ) -> None:
        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta, gtf, out_dir)
        _BLACKLIST.to_feather(out_dir / SJ_BLACKLIST_FEATHER)

        assert TranscriptIndex.load(out_dir).sj_blacklist_size == 0

    def test_a_recorded_blacklist_removed_loads_with_detection_off(
        self, tmp_path: Path, mini_index_inputs
    ) -> None:
        """Removing the feather from a store-built index gives its blacklist-free twin."""
        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta, gtf, out_dir)
        rebuild_with_splice_blacklist(out_dir, _BLACKLIST)
        (out_dir / SJ_BLACKLIST_FEATHER).unlink()

        assert TranscriptIndex.load(out_dir).sj_blacklist_size == 0

    def test_a_manifest_without_sources_loads_with_no_blacklist(
        self, tmp_path: Path, mini_index_inputs
    ) -> None:
        """A manifest written before `sources` existed records no store, so its feather is not applied."""
        fasta, gtf = mini_index_inputs
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(fasta, gtf, out_dir)
        rebuild_with_splice_blacklist(out_dir, _BLACKLIST)
        manifest = json.loads((out_dir / MANIFEST_JSON).read_text())
        del manifest["sources"]
        (out_dir / MANIFEST_JSON).write_text(json.dumps(manifest))

        assert TranscriptIndex.load(out_dir).sj_blacklist_size == 0


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
