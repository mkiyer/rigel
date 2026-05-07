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
