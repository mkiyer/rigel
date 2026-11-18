"""GTF parsing has two modes, and each must behave as its name promises.

Strict raises on a malformed record and names the line number, because an annotation that silently
loses rows produces an index nobody can account for. Warn-and-skip continues past the same record,
and the choice reaches `Transcript.read_gtf` rather than stopping at the parser.
"""

import logging
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
    """GTF asks for a `;` after every attribute and many files omit the last one, so the final attribute (here the
    transcript id) is read without it rather than dropped in silence."""
    r = GTFRecord.from_str(f"chr1\ttest\texon\t100\t200\t.\t+\t.\t{attrs}")
    assert r.attrs["gene_id"] == "G1"
    assert r.attrs["transcript_id"] == "T1"


def _write(tmp_path: Path, lines: list[str]) -> str:
    path = tmp_path / "annotation.gtf"
    path.write_text("".join(lines))
    return str(path)


def test_an_exon_without_a_transcript_id_follows_the_parse_mode(tmp_path: Path):
    """An exon with no `transcript_id` is judged by the parse mode, not raised as a bare KeyError that stops even
    warn-and-skip. Strict names it; warn-and-skip drops that exon and keeps the rest."""
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
    """A transcript's exons share its reference and strand; a stray exon merged in would build a transcript that spans
    two chromosomes or two strands. Strict names it; warn-and-skip drops that exon and keeps the rest."""
    path = _write(
        tmp_path,
        ['chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n', stray],
    )
    with pytest.raises(ValueError, match="t1"):
        Transcript.read_gtf(path, parse_mode="strict")
    tx = Transcript.read_gtf(path, parse_mode="warn-skip")
    assert len(tx) == 1 and tx[0].ref == "chr1" and len(tx[0].exons) == 1


def _exon(start, end, gene="g1", transcript="t1", ref="chr1", strand="+", extra=""):
    return f'{ref}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}"; transcript_id "{transcript}";{extra}\n'


def test_an_exon_ending_before_it_starts_follows_the_parse_mode(tmp_path: Path):
    """An exon with its end before its start is judged by the parse mode, rather than building an index on which
    `rigel quant` fails far from the cause."""
    path = _write(tmp_path, [_exon(100, 200), _exon(500, 400)])
    with pytest.raises(ValueError, match="line 2"):
        Transcript.read_gtf(path, parse_mode="strict")
    tx = Transcript.read_gtf(path, parse_mode="warn-skip")
    assert [(e.start, e.end) for e in tx[0].exons] == [(99, 200)]


@pytest.mark.parametrize(
    "second",
    [(150, 300), (100, 200)],
    ids=["overlapping", "duplicated"],
)
def test_overlapping_or_duplicated_exons_follow_the_parse_mode(tmp_path: Path, second):
    """Exons of one transcript are disjoint. Overlapping or duplicated ones would lengthen the transcript and write
    negative or empty introns into the splice graph."""
    path = _write(tmp_path, [_exon(100, 200), _exon(*second), _exon(500, 600)])
    with pytest.raises(ValueError, match="overlaps"):
        Transcript.read_gtf(path, parse_mode="strict")
    tx = Transcript.read_gtf(path, parse_mode="warn-skip")
    assert [(e.start, e.end) for e in tx[0].exons] == [(99, 200), (499, 600)]
    assert tx[0].length == 202


def test_abutting_exons_are_one_exon(tmp_path: Path):
    """Exons with no base between them (RefSeq writes them at genome indels) are one block of the transcript. Kept apart,
    they would make a 0-bp intron: the transcript would turn multi-exon and its synthetic span take half its reads."""
    path = _write(tmp_path, [_exon(100, 300), _exon(301, 500)])
    for mode in ("strict", "warn-skip"):
        (t,) = Transcript.read_gtf(path, parse_mode=mode)
        assert (
            [(e.start, e.end) for e in t.exons] == [(99, 500)]
            and t.length == 401
            and t.is_single_exon
        )


def test_a_transcript_id_shared_by_two_genes_follows_the_parse_mode(tmp_path: Path):
    """Two genes' transcripts sharing an ID on one reference and strand would merge into one transcript behind a fake
    intron, and the second gene would vanish from the tables."""
    path = _write(tmp_path, [_exon(100, 200, gene="g1"), _exon(5000, 5200, gene="g2")])
    with pytest.raises(ValueError, match="g2"):
        Transcript.read_gtf(path, parse_mode="strict")
    tx = Transcript.read_gtf(path, parse_mode="warn-skip")
    assert len(tx) == 1 and tx[0].g_id == "g1" and len(tx[0].exons) == 1


def test_an_ensembl_gene_biotype_is_the_gene_type(tmp_path: Path):
    path = _write(tmp_path, [_exon(100, 200, extra=' gene_biotype "lncRNA";')])
    (t,) = Transcript.read_gtf(path, parse_mode="strict")
    assert t.g_type == "lncRNA"


def _fasta(tmp_path: Path, length: int, refs: tuple[str, ...] = ("chr1",)) -> Path:
    import pysam

    fasta = tmp_path / "genome.fa"
    fasta.write_text("".join(f">{ref}\n" + "A" * length + "\n" for ref in refs))
    pysam.faidx(str(fasta))
    return fasta


@pytest.mark.parametrize(
    "stray",
    [
        _exon(100, 300, gene="g2", transcript="t2", ref="chrUn"),
        _exon(900, 1200, gene="g2", transcript="t2"),
    ],
    ids=["contig-missing-from-the-fasta", "past-the-contig-end"],
)
def test_a_transcript_off_the_fasta_follows_the_parse_mode(tmp_path: Path, stray):
    """A transcript on a contig the FASTA lacks, or past a contig's end, follows the parse mode rather than failing the
    build with a bare error. Strict refuses it by name; warn-and-skip drops it and builds the rest."""
    from rigel.index import TranscriptIndex

    fasta = _fasta(tmp_path, 1000)
    gtf = _write(tmp_path, [_exon(100, 300), _exon(500, 700), stray])
    with pytest.raises(ValueError, match="t2"):
        TranscriptIndex.build(fasta, gtf, tmp_path / "strict", gtf_parse_mode="strict")
    TranscriptIndex.build(fasta, gtf, tmp_path / "skip", gtf_parse_mode="warn-skip")
    assert list(TranscriptIndex.load(tmp_path / "skip").t_df["t_id"])[:1] == ["t1"]


@pytest.mark.parametrize("strand", [".", "?", "x"])
def test_a_transcript_with_no_strand_follows_the_parse_mode(tmp_path: Path, caplog, strand):
    """Only a `+` or `-` transcript can be quantified: the scorer reads a stranded fragment against its transcript's
    strand, and calibration's partition holds only stranded transcripts, so a strandless one's exons would be read as
    intergenic and its RNA booked as gDNA. Strict refuses it and names the option that drops it (StringTie writes `.`
    for a transcript it cannot strand); warn-and-skip drops it with a warning, and a strand GTF has no word for, like
    `x`, is skipped the same way rather than stopping the run."""
    path = _write(
        tmp_path,
        [
            _exon(100, 200),
            _exon(300, 400),
            _exon(1000, 1200, gene="g2", transcript="t2", strand=strand),
            _exon(1500, 1600, gene="g2", transcript="t2", strand=strand),
        ],
    )
    with pytest.raises(ValueError, match=r"t2.*--gtf-parse-mode warn-skip"):
        Transcript.read_gtf(path, parse_mode="strict")
    with caplog.at_level(logging.WARNING):
        tx = Transcript.read_gtf(path, parse_mode="warn-skip")
    assert [t.t_id for t in tx] == ["t1"]
    assert f"strand {strand!r}" in caplog.text


def test_rigel_index_refuses_a_transcript_with_no_strand(tmp_path: Path):
    """The issue's case end to end: a two-exon `.` transcript is refused by the build under strict parsing, and is
    absent from the index under warn-and-skip."""
    from rigel.index import TranscriptIndex

    fasta = _fasta(tmp_path, 5000)
    gtf = _write(
        tmp_path,
        [
            _exon(100, 300),
            _exon(500, 700),
            _exon(3001, 3300, gene="g2", transcript="t2", strand="."),
            _exon(3601, 3900, gene="g2", transcript="t2", strand="."),
        ],
    )
    with pytest.raises(ValueError, match="--gtf-parse-mode warn-skip"):
        TranscriptIndex.build(fasta, gtf, tmp_path / "strict", gtf_parse_mode="strict")
    TranscriptIndex.build(fasta, gtf, tmp_path / "skip", gtf_parse_mode="warn-skip")
    t_df = TranscriptIndex.load(tmp_path / "skip").t_df
    assert "t2" not in set(t_df["t_id"])
    assert set(t_df.loc[~t_df["is_synthetic"], "strand"]) == {1}


def test_an_exon_without_a_gene_id_follows_the_parse_mode(tmp_path: Path):
    """An exon with no `gene_id` is judged like one with no `transcript_id`: every such transcript would otherwise
    share one gene, keyed `None`, whatever its reference or strand. Strict names it; warn-and-skip drops it."""
    path = _write(
        tmp_path,
        [
            _exon(100, 200),
            'chr1\ttest\texon\t300\t400\t.\t+\t.\ttranscript_id "t2";\n',
            'chr2\ttest\texon\t300\t400\t.\t-\t.\ttranscript_id "t3";\n',
        ],
    )
    with pytest.raises(ValueError, match="t2.*no gene_id"):
        Transcript.read_gtf(path, parse_mode="strict")
    tx = Transcript.read_gtf(path, parse_mode="warn-skip")
    assert [t.t_id for t in tx] == ["t1"]


def test_rigel_index_holds_no_gene_without_a_gene_id(tmp_path: Path):
    """The issue's case end to end: two transcripts with no `gene_id`, on chr1(+) and chr2(-), are refused by the
    build under strict parsing, and under warn-and-skip neither is indexed and the gene table has no `None` gene."""
    from rigel.index import TranscriptIndex

    fasta = _fasta(tmp_path, 1000, refs=("chr1", "chr2"))
    gtf = _write(
        tmp_path,
        [
            _exon(100, 300),
            'chr1\ttest\texon\t500\t700\t.\t+\t.\ttranscript_id "t2";\n',
            'chr2\ttest\texon\t500\t700\t.\t-\t.\ttranscript_id "t3";\n',
        ],
    )
    with pytest.raises(ValueError, match="gene_id"):
        TranscriptIndex.build(fasta, gtf, tmp_path / "strict", gtf_parse_mode="strict")
    TranscriptIndex.build(fasta, gtf, tmp_path / "skip", gtf_parse_mode="warn-skip")
    index = TranscriptIndex.load(tmp_path / "skip")
    assert set(index.t_df["t_id"]) == {"t1"}
    assert index.g_df["g_id"].notna().all()


def test_a_gene_on_two_references_follows_the_parse_mode(tmp_path: Path, caplog):
    """A gene's row takes one reference and strand and spans its transcripts' extremes, so a `gene_id` with
    transcripts on two references would be one gene on one reference spanning coordinates of both. Strict refuses it
    by name; warn-and-skip drops the gene, since the annotation cannot say which reference is its own, and builds the
    rest."""
    from rigel.index import TranscriptIndex

    fasta = _fasta(tmp_path, 1000, refs=("chr1", "chr2"))
    gtf = _write(
        tmp_path,
        [
            _exon(100, 300),
            _exon(500, 700),
            _exon(100, 300, transcript="t2", ref="chr2"),
            _exon(100, 300, gene="g2", transcript="t3", ref="chr2", strand="-"),
        ],
    )
    with pytest.raises(ValueError, match=r"g1 .*chr1.*chr2"):
        TranscriptIndex.build(fasta, gtf, tmp_path / "strict", gtf_parse_mode="strict")
    with caplog.at_level(logging.WARNING):
        TranscriptIndex.build(fasta, gtf, tmp_path / "skip", gtf_parse_mode="warn-skip")
    index = TranscriptIndex.load(tmp_path / "skip")
    annotated = index.t_df[~index.t_df["is_synthetic"]]
    assert set(annotated["t_id"]) == {"t3"}
    assert "g1" in caplog.text
