"""Tests for GTF parse modes (strict vs warn-skip)."""

from pathlib import Path

import pytest

from hulkrna.gtf import GTFRecord
from hulkrna.transcript import Transcript


def _lines_with_malformed_middle() -> list[str]:
    return [
        'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        'chr1\ttest\texon\tBAD\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        'chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
    ]


def test_gtf_parse_strict_raises_with_line_number():
    with pytest.raises(ValueError, match="line 2"):
        list(GTFRecord.parse(_lines_with_malformed_middle(), parse_mode="strict"))


def test_gtf_parse_warn_skip_continues():
    rows = list(GTFRecord.parse(_lines_with_malformed_middle(), parse_mode="warn-skip"))
    assert len(rows) == 2
    assert all(r.feature == "exon" for r in rows)


def test_gtf_parse_file_strict_and_warn_skip(tmp_path: Path):
    gtf_path = tmp_path / "test.gtf"
    gtf_path.write_text("".join(_lines_with_malformed_middle()))

    with pytest.raises(ValueError, match="line 2"):
        list(GTFRecord.parse_file(gtf_path, parse_mode="strict"))

    rows = list(GTFRecord.parse_file(gtf_path, parse_mode="warn-skip"))
    assert len(rows) == 2


def test_transcript_read_gtf_respects_parse_mode(tmp_path: Path):
    gtf_path = tmp_path / "test_tx.gtf"
    gtf_path.write_text("".join(_lines_with_malformed_middle()))

    with pytest.raises(ValueError, match="line 2"):
        Transcript.read_gtf(str(gtf_path), parse_mode="strict")

    tx = Transcript.read_gtf(str(gtf_path), parse_mode="warn-skip")
    assert len(tx) == 1
    assert tx[0].t_id == "t1"
