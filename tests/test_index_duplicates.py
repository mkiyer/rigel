"""Tests for the duplicate-transcript guard in ``read_transcripts``."""
from __future__ import annotations

from pathlib import Path

import pytest

from rigel.index import read_transcripts


def _write_gtf(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


def _gtf_line(feature: str, start: int, end: int, strand: str, gid: str, tid: str,
              ref: str = "chr1") -> str:
    return (
        f'{ref}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t'
        f'gene_id "{gid}"; transcript_id "{tid}";'
    )


def test_duplicate_exon_structure_raises(tmp_path: Path):
    """Two transcripts with identical exon coordinates must raise ValueError."""
    gtf = tmp_path / "dup.gtf"
    _write_gtf(gtf, [
        _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
        _gtf_line("exon", 100, 200, "+", "G1", "T1"),
        _gtf_line("exon", 400, 500, "+", "G1", "T1"),
        # Identical exon set, different transcript_id
        _gtf_line("transcript", 100, 500, "+", "G1", "T2"),
        _gtf_line("exon", 100, 200, "+", "G1", "T2"),
        _gtf_line("exon", 400, 500, "+", "G1", "T2"),
    ])
    with pytest.raises(ValueError, match=r"identical exon coordinates"):
        read_transcripts(gtf)


def test_same_intron_chain_different_utrs_allowed(tmp_path: Path):
    """Transcripts that share an intron chain but differ in UTR boundaries
    are biologically realistic and must not trigger the guard."""
    gtf = tmp_path / "utrs.gtf"
    _write_gtf(gtf, [
        _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
        _gtf_line("exon", 100, 200, "+", "G1", "T1"),
        _gtf_line("exon", 400, 500, "+", "G1", "T1"),
        # Same intron (200,400) but different 5' and 3' UTR
        _gtf_line("transcript", 120, 480, "+", "G1", "T2"),
        _gtf_line("exon", 120, 200, "+", "G1", "T2"),
        _gtf_line("exon", 400, 480, "+", "G1", "T2"),
    ])
    txs = read_transcripts(gtf)
    assert len(txs) == 2
    assert {t.t_id for t in txs} == {"T1", "T2"}


def test_different_strand_not_a_duplicate(tmp_path: Path):
    """Identical coords on opposite strands are distinct transcripts."""
    gtf = tmp_path / "strand.gtf"
    _write_gtf(gtf, [
        _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
        _gtf_line("exon", 100, 200, "+", "G1", "T1"),
        _gtf_line("exon", 400, 500, "+", "G1", "T1"),
        _gtf_line("transcript", 100, 500, "-", "G2", "T2"),
        _gtf_line("exon", 100, 200, "-", "G2", "T2"),
        _gtf_line("exon", 400, 500, "-", "G2", "T2"),
    ])
    txs = read_transcripts(gtf)
    assert len(txs) == 2


def test_duplicate_message_lists_offending_ids(tmp_path: Path):
    gtf = tmp_path / "dup.gtf"
    _write_gtf(gtf, [
        _gtf_line("transcript", 100, 500, "+", "G1", "ALPHA"),
        _gtf_line("exon", 100, 500, "+", "G1", "ALPHA"),
        _gtf_line("transcript", 100, 500, "+", "G1", "BETA"),
        _gtf_line("exon", 100, 500, "+", "G1", "BETA"),
    ])
    with pytest.raises(ValueError) as excinfo:
        read_transcripts(gtf)
    msg = str(excinfo.value)
    assert "ALPHA" in msg
    assert "BETA" in msg


def test_duplicate_error_mentions_collapse_flag(tmp_path: Path):
    """The error should point users at the --collapse-duplicate-transcripts escape hatch."""
    gtf = tmp_path / "dup.gtf"
    _write_gtf(gtf, [
        _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
        _gtf_line("exon", 100, 500, "+", "G1", "T1"),
        _gtf_line("transcript", 100, 500, "+", "G1", "T2"),
        _gtf_line("exon", 100, 500, "+", "G1", "T2"),
    ])
    with pytest.raises(ValueError, match=r"--collapse-duplicate-transcripts"):
        read_transcripts(gtf)


def test_collapse_keeps_lexicographically_smallest_id(tmp_path: Path):
    """With collapse enabled, a duplicate group keeps only the lexicographically
    smallest transcript ID; distinct transcripts are untouched and t_index is
    reassigned contiguously."""
    gtf = tmp_path / "collapse.gtf"
    _write_gtf(gtf, [
        # duplicate pair (identical exons), IDs deliberately out of lexical order
        _gtf_line("transcript", 100, 500, "+", "G1", "ENST00000002"),
        _gtf_line("exon", 100, 200, "+", "G1", "ENST00000002"),
        _gtf_line("exon", 400, 500, "+", "G1", "ENST00000002"),
        _gtf_line("transcript", 100, 500, "+", "G1", "ENST00000001"),
        _gtf_line("exon", 100, 200, "+", "G1", "ENST00000001"),
        _gtf_line("exon", 400, 500, "+", "G1", "ENST00000001"),
        # a distinct transcript (different exons) — must survive
        _gtf_line("transcript", 1000, 1500, "+", "G2", "ENST00000009"),
        _gtf_line("exon", 1000, 1500, "+", "G2", "ENST00000009"),
    ])
    txs = read_transcripts(gtf, collapse_duplicate_transcripts=True)
    assert {t.t_id for t in txs} == {"ENST00000001", "ENST00000009"}
    # dropped transcript is gone; kept indices are contiguous 0..N-1
    assert sorted(t.t_index for t in txs) == [0, 1]


def test_collapse_three_way_group_keeps_one(tmp_path: Path):
    """A 3-way duplicate group collapses to the single lexicographically smallest ID."""
    gtf = tmp_path / "collapse3.gtf"
    lines: list[str] = []
    for tid in ("ENST0000C", "ENST0000A", "ENST0000B"):
        lines += [
            _gtf_line("transcript", 100, 500, "+", "G1", tid),
            _gtf_line("exon", 100, 200, "+", "G1", tid),
            _gtf_line("exon", 400, 500, "+", "G1", tid),
        ]
    _write_gtf(gtf, lines)
    txs = read_transcripts(gtf, collapse_duplicate_transcripts=True)
    assert {t.t_id for t in txs} == {"ENST0000A"}


def test_collapse_is_noop_without_duplicates(tmp_path: Path):
    """collapse=True must not drop transcripts that are merely similar (shared intron,
    different UTRs)."""
    gtf = tmp_path / "nodup.gtf"
    _write_gtf(gtf, [
        _gtf_line("transcript", 100, 500, "+", "G1", "T1"),
        _gtf_line("exon", 100, 200, "+", "G1", "T1"),
        _gtf_line("exon", 400, 500, "+", "G1", "T1"),
        _gtf_line("transcript", 120, 480, "+", "G1", "T2"),
        _gtf_line("exon", 120, 200, "+", "G1", "T2"),
        _gtf_line("exon", 400, 480, "+", "G1", "T2"),
    ])
    txs = read_transcripts(gtf, collapse_duplicate_transcripts=True)
    assert {t.t_id for t in txs} == {"T1", "T2"}
