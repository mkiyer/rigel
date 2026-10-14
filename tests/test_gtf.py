"""GTF parsing has two modes, and each must behave as its name promises.

Strict raises on a malformed record and names the line number, because an annotation that silently
loses rows produces an index nobody can account for. Warn-and-skip continues past the same record,
and the choice reaches `Transcript.read_gtf` rather than stopping at the parser.
"""

from pathlib import Path

import pytest

from rigel.gtf import GTFRecord
from rigel.transcript import Transcript


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


# ── what a real annotation carries: a last attribute without its `;`, a record with no transcript, a stray exon ──


@pytest.mark.parametrize(
    "attrs",
    ['gene_id "G1"; transcript_id "T1"', "gene_id G1; transcript_id T1"],
    ids=["quoted", "unquoted"],
)
def test_the_last_attribute_needs_no_trailing_semicolon(attrs):
    """GTF asks for a `;` after every attribute and many files omit the last one; the pattern that required it dropped
    the final attribute — here the transcript id — without a word."""
    r = GTFRecord.from_str(f"chr1\ttest\texon\t100\t200\t.\t+\t.\t{attrs}")
    assert r.attrs["gene_id"] == "G1"
    assert r.attrs["transcript_id"] == "T1"


def _write(tmp_path: Path, lines: list[str]) -> str:
    path = tmp_path / "annotation.gtf"
    path.write_text("".join(lines))
    return str(path)


def test_an_exon_without_a_transcript_id_follows_the_parse_mode(tmp_path: Path):
    """An exon with no `transcript_id` raised a bare KeyError from `read_gtf`, outside the parse mode, so even
    warn-and-skip stopped on it. Strict names it; warn-and-skip drops that exon and keeps the rest."""
    path = _write(
        tmp_path,
        [
            'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
            'chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g2";\n',
        ],
    )
    with pytest.raises(ValueError, match="transcript_id"):
        Transcript.read_gtf(path, parse_mode="strict")
    tx = Transcript.read_gtf(path, parse_mode="warn-skip")
    assert [t.t_id for t in tx] == ["t1"]


@pytest.mark.parametrize(
    "stray",
    [
        'chr2\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        'chr1\ttest\texon\t300\t400\t.\t-\t.\tgene_id "g1"; transcript_id "t1";\n',
    ],
    ids=["another-reference", "another-strand"],
)
def test_an_exon_off_its_transcripts_reference_or_strand_follows_the_parse_mode(
    tmp_path: Path, stray
):
    """A transcript's exons share its reference and strand; a stray exon was merged in silently, building a transcript
    that spans two chromosomes or two strands. Strict names it; warn-and-skip drops that exon and keeps the rest."""
    path = _write(
        tmp_path,
        ['chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n', stray],
    )
    with pytest.raises(ValueError, match="t1"):
        Transcript.read_gtf(path, parse_mode="strict")
    tx = Transcript.read_gtf(path, parse_mode="warn-skip")
    assert len(tx) == 1 and tx[0].ref == "chr1" and len(tx[0].exons) == 1
