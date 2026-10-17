"""Per-fragment truth is carried in the read name, and these are the helpers that write and read it.

Writing a read name and parsing an origin back out of a FASTQ or BAM read name for each of mature,
nascent and gDNA, including ids that contain ":", and refusing a strand outside the alphabet;
counting origins over a whole FASTQ; writing post-capture truth from the origins actually observed
rather than from what was requested, each quantity in one column; the molar table, written exactly
and read back as an abundance input, which a condition's table is refused as; each truth file under
one key in the manifest and the summary; the condition directory naming; and the manifest round trip.
Every truth-scoring instrument in the project reads through these, so a silent parse failure would
make the oracle agree with anything.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

from rigel.sim.manifest import condition_dir_name, gdna_label_for_rate, write_manifest
from rigel.sim.read_name import (
    STRAND_FORWARD,
    STRAND_REVERSE,
    STRANDS,
    Origin,
    format_gdna_name,
    format_rna_name,
    parse_origin,
)
from rigel.sim.truth import count_origins_from_fastq, write_molar_truth, write_post_capture_truth
from rigel.sim.whole_genome import assign_file_abundances
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _truth_transcript(
    transcript_id: str,
    abundance: float,
    nrna_abundance: float = 0.0,
) -> Transcript:
    transcript = Transcript(
        ref="chr1",
        strand=Strand.POS,
        exons=[Interval(0, 500)],
        t_id=transcript_id,
        g_id=f"g_{transcript_id}",
        abundance=abundance,
    )
    transcript.nrna_abundance = nrna_abundance
    transcript.compute_length()
    return transcript


def test_parse_mrna_fastq_origin():
    origin = parse_origin("@tx1:10-250:f:7/1")

    assert origin.kind == "mrna"
    assert origin.transcript_id == "tx1"
    assert origin.start == 10
    assert origin.end == 250
    assert origin.strand == "f"
    assert origin.index == 7


def test_parse_nrna_bam_origin():
    origin = parse_origin("nrna_tx2:4-300:r:11")

    assert origin.kind == "nrna"
    assert origin.transcript_id == "tx2"
    assert origin.start == 4
    assert origin.end == 300
    assert origin.strand == "r"
    assert origin.index == 11


def test_parse_gdna_origin():
    origin = parse_origin("gdna:chrSynthetic:100-390:r:3/2")

    assert origin.kind == "gdna"
    assert origin.transcript_id is None
    assert origin.ref == "chrSynthetic"
    assert origin.start == 100
    assert origin.end == 390
    assert origin.strand == "r"
    assert origin.index == 3


def test_parse_gdna_origin_requires_a_reference():
    with pytest.raises(ValueError, match="without a reference"):
        parse_origin("gdna:100-390:f:3")


def test_parse_origin_keeps_a_colon_inside_a_contig_or_transcript_id():
    gdna = parse_origin("gdna:HLA-A*01:01:01:01:100-390:f:3")
    assert (gdna.kind, gdna.ref, gdna.start, gdna.end, gdna.strand, gdna.index) == (
        "gdna",
        "HLA-A*01:01:01:01",
        100,
        390,
        "f",
        3,
    )

    mrna = parse_origin("tx:1:10-250:r:7")
    assert (mrna.kind, mrna.transcript_id, mrna.start, mrna.end) == ("mrna", "tx:1", 10, 250)

    nrna = parse_origin("nrna_tx:1:4-300:f:11")
    assert (nrna.kind, nrna.transcript_id, nrna.start, nrna.end) == ("nrna", "tx:1", 4, 300)


def test_the_format_functions_write_the_simulator_names():
    """The exact bytes: every oracle BAM and FASTQ already simulated carries names of this form."""
    assert format_rna_name("tx1", 10, 250, STRAND_FORWARD, 7, nascent=False) == "tx1:10-250:f:7"
    assert format_rna_name("tx1", 4, 300, STRAND_REVERSE, 11, nascent=True) == (
        "nrna_tx1:4-300:r:11"
    )
    assert format_gdna_name("chr1", 100, 390, STRAND_REVERSE, 3) == "gdna:chr1:100-390:r:3"


def test_a_formatted_name_parses_back_to_its_fields():
    """The writer and the parser are inverses: every kind, both strands, a transcript id or contig
    that contains ":", and the same name as a FASTQ header."""
    for strand in STRANDS:
        cases = [
            (
                format_rna_name("tx:1", 10, 250, strand, 7, nascent=False),
                Origin("mrna", "tx:1", None, 10, 250, strand, 7),
            ),
            (
                format_rna_name("tx:1", 4, 300, strand, 11, nascent=True),
                Origin("nrna", "tx:1", None, 4, 300, strand, 11),
            ),
            (
                format_gdna_name("HLA-A*01:01", 100, 390, strand, 3),
                Origin("gdna", None, "HLA-A*01:01", 100, 390, strand, 3),
            ),
        ]
        for name, origin in cases:
            assert parse_origin(name) == origin, name
            assert parse_origin(f"@{name}/2") == origin, name


def test_parse_origin_refuses_a_strand_outside_the_alphabet():
    """``+``, ``-`` and ``.`` are the GTF's strand letters, never the simulator's: a name carrying one,
    or any letter outside the alphabet, was not written by the simulator and is refused rather than
    read as truth."""
    for qname in (
        "tx1:10-250:+:7",
        "nrna_tx1:4-300:-:11",
        "gdna:chr1:100-390:.:3",
        "tx1:10-250:F:7",
        "tx1:10-250::7",
    ):
        with pytest.raises(ValueError, match="invalid strand"):
            parse_origin(qname)


def test_count_origins_from_fastq(tmp_path):
    fastq = tmp_path / "reads.fq"
    fastq.write_text(
        "@tx1:0-100:f:0/1\nACGT\n+\nIIII\n"
        "@nrna_tx1:0-100:f:1/1\nACGT\n+\nIIII\n"
        "@gdna:chr1:0-100:r:2/1\nACGT\n+\nIIII\n"
    )

    counts = count_origins_from_fastq(fastq)

    assert counts["mrna"] == 1
    assert counts["nrna"] == 1
    assert counts["gdna"] == 1


def test_write_post_capture_truth_uses_observed_read_origins(tmp_path):
    fastq = tmp_path / "reads.fq"
    fastq.write_text(
        "@TX1:0-100:f:0/1\nACGT\n+\nIIII\n"
        "@nrna_TX1:20-220:f:1/1\nACGT\n+\nIIII\n"
        "@TX2:0-80:f:2/1\nACGT\n+\nIIII\n"
        "@gdna:chr1:5-155:r:3/1\nACGT\n+\nIIII\n"
    )
    transcripts = [
        _truth_transcript("TX1", 10.0, 2.0),
        _truth_transcript("TX2", 50.0),
    ]

    summary = write_post_capture_truth(
        transcripts,
        tmp_path / "truth_abundances.tsv",
        tmp_path / "truth_fragment_lengths.tsv",
        tmp_path / "truth_summary.json",
        fastq_path=fastq,
        condition="cond",
        pre_capture_abundances="molecular.tsv",
    )

    # Each quantity in one column: the molecules before capture, then the fragments realised. The
    # realised count's copies under molar names (`mrna_abundance`, ...) and as `post_capture_*` are gone.
    header = (tmp_path / "truth_abundances.tsv").read_text().splitlines()[0].split("\t")
    assert header == [
        "transcript_id",
        "gene_id",
        "gene_name",
        "ref",
        "strand",
        "pre_capture_mrna_abundance",
        "pre_capture_nrna_abundance",
        "pre_capture_total_rna",
        "observed_mrna_fragments",
        "observed_nrna_fragments",
        "observed_total_rna_fragments",
        "n_exons",
        "spliced_length",
        "genomic_span",
    ]
    abundances = pd.read_csv(tmp_path / "truth_abundances.tsv", sep="\t")
    tx1 = abundances.set_index("transcript_id").loc["TX1"]
    tx2 = abundances.set_index("transcript_id").loc["TX2"]
    assert tx1["pre_capture_mrna_abundance"] == 10.0
    assert tx1["pre_capture_nrna_abundance"] == 2.0
    assert tx1["pre_capture_total_rna"] == 12.0
    assert tx1["observed_mrna_fragments"] == 1
    assert tx1["observed_nrna_fragments"] == 1
    assert tx1["observed_total_rna_fragments"] == 2
    assert (tx2["observed_mrna_fragments"], tx2["pre_capture_mrna_abundance"]) == (1, 50.0)

    fl = pd.read_csv(tmp_path / "truth_fragment_lengths.tsv", sep="\t")
    gdna_150 = fl[(fl["kind"] == "gdna") & (fl["fragment_length"] == 150)]
    rna_200 = fl[(fl["kind"] == "rna") & (fl["fragment_length"] == 200)]
    assert gdna_150["count"].item() == 1
    assert rna_200["count"].item() == 1
    assert summary["origin_counts"] == {"mrna": 2, "nrna": 1, "gdna": 1}

    # one key per path: the three tables' paths are written once each, never repeated under `files`
    written_summary = json.loads((tmp_path / "truth_summary.json").read_text())
    assert list(written_summary) == [
        "condition",
        "truth_kind",
        "gdna_strand_overdispersion_true",
        "pre_capture_abundances",
        "post_capture_abundances",
        "post_capture_fragment_lengths",
        "origin_counts",
        "fragment_lengths",
    ]
    assert written_summary["truth_kind"] == "post_capture_empirical"
    assert written_summary["pre_capture_abundances"] == "molecular.tsv"
    assert written_summary["post_capture_abundances"] == str(tmp_path / "truth_abundances.tsv")
    assert written_summary["post_capture_fragment_lengths"] == str(
        tmp_path / "truth_fragment_lengths.tsv"
    )
    assert written_summary["fragment_lengths"]["mrna"]["n"] == 2


def _molar_transcripts() -> list[Transcript]:
    """A multi-exon transcript starting at base 0, a single-exon nascent entity with no gene name (as
    the index writes a synthetic span), and a minus-strand transcript."""
    spec = [
        ("spliced", Strand.POS, [(0, 100), (200, 300)], 12.5, 0.0, "gene_spliced"),
        ("span", Strand.POS, [(0, 300)], 0.0, 3.14159, ""),
        ("minus", Strand.NEG, [(1000, 1500)], 1234.56789, 0.0, "gene_minus"),
    ]
    transcripts = []
    for t_id, strand, exons, abundance, nrna, g_name in spec:
        transcript = Transcript(
            ref="chr1",
            strand=strand,
            exons=[Interval(start, end) for start, end in exons],
            t_id=t_id,
            g_id=f"g_{t_id}",
            g_name=g_name,
            abundance=abundance,
        )
        transcript.nrna_abundance = nrna
        transcript.compute_length()
        transcripts.append(transcript)
    return transcripts


def test_the_molar_table_is_written_exactly_and_reads_back_as_an_abundance_input(tmp_path):
    """The suite-level molar table: the molecules each transcript was simulated at, to four decimals,
    under the simulator's own input names, so it reads back as ``abundance.file`` and writes the same
    bytes again. Its writer once had no gate, and a fix reached it unguarded: the genomic span of a
    transcript starting at base 0 (``spliced``), recorded as 0 behind a truthiness test of its start.

    Perturbation: that truthiness test in ``truth._transcript_metadata`` writes ``spliced``'s span as
    0, and ``.3f`` in ``truth._molar_row`` moves every abundance column."""
    path = tmp_path / "truth_abundances_nrna_mid.tsv"
    write_molar_truth(_molar_transcripts(), path)
    assert path.read_text() == (
        "transcript_id\tgene_id\tgene_name\tref\tstrand\tmrna_abundance\tnrna_abundance\t"
        "total_rna\tn_exons\tspliced_length\tgenomic_span\n"
        "spliced\tg_spliced\tgene_spliced\tchr1\t+\t12.5000\t0.0000\t12.5000\t2\t200\t300\n"
        "span\tg_span\t\tchr1\t+\t0.0000\t3.1416\t3.1416\t1\t300\t300\n"
        "minus\tg_minus\tgene_minus\tchr1\t-\t1234.5679\t0.0000\t1234.5679\t1\t500\t500\n"
    )

    read_back = _molar_transcripts()
    for transcript in read_back:
        transcript.abundance, transcript.nrna_abundance = None, 0.0
    assert assign_file_abundances(read_back, path) is True, (
        "the table's nascent column was not read"
    )
    assert [(t.abundance, t.nrna_abundance) for t in read_back] == [
        (12.5, 0.0),
        (0.0, 3.1416),
        (1234.5679, 0.0),
    ]
    again = tmp_path / "again.tsv"
    write_molar_truth(read_back, again)
    assert again.read_text() == path.read_text()


def test_a_condition_truth_table_is_refused_as_an_abundance_input(tmp_path):
    """A condition's table holds the fragments that condition realised, never molecules, so it is no
    abundance input. The table this writer writes carries no molar-named column at all; every panel
    simulated before it carries the realised counts under ``mrna_abundance`` / ``nrna_abundance`` too,
    which the ``sim`` format would read as molecules — so both are refused by the realised column.

    Perturbation: dropping the refusal in ``whole_genome._detect_abundance_format`` reads the second
    table's fragment counts as molar abundances."""
    fastq = tmp_path / "reads.fq"
    fastq.write_text("@spliced:0-100:f:0/1\nACGT\n+\nIIII\n")
    current = tmp_path / "cond" / "truth_abundances.tsv"
    write_post_capture_truth(
        _molar_transcripts(),
        current,
        tmp_path / "cond" / "truth_fragment_lengths.tsv",
        tmp_path / "cond" / "truth_summary.json",
        fastq_path=fastq,
    )
    # the header of every condition table on disk before the duplicate columns were dropped
    on_disk_columns = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "ref",
        "strand",
        "mrna_abundance",
        "nrna_abundance",
        "total_rna",
        "pre_capture_mrna_abundance",
        "pre_capture_nrna_abundance",
        "pre_capture_total_rna",
        "post_capture_mrna_fragments",
        "post_capture_nrna_fragments",
        "post_capture_total_rna_fragments",
        "observed_mrna_fragments",
        "observed_nrna_fragments",
        "observed_total_rna_fragments",
        "n_exons",
        "spliced_length",
        "genomic_span",
    ]
    on_disk_row = ["spliced", "g_spliced", "gene_spliced", "chr1", "+", "1.0", "0.0", "1.0"]
    on_disk_row += ["12.5", "0.0", "12.5", "1", "0", "1", "1", "0", "1", "2", "200", "300"]
    on_disk = tmp_path / "on_disk.tsv"
    on_disk.write_text("\t".join(on_disk_columns) + "\n" + "\t".join(on_disk_row) + "\n")
    for table in (current, on_disk):
        with pytest.raises(ValueError, match="condition's truth table"):
            assign_file_abundances(_molar_transcripts(), table)


def test_condition_dir_name_and_gdna_labels():
    assert condition_dir_name("low", 0.9, "rand") == "gdna_low_ss_0.90_nrna_rand"
    assert condition_dir_name("low", 0.9, "rand", "on") == "gdna_low_ss_0.90_nrna_rand_capture_on"
    assert gdna_label_for_rate(0, None, 0) == "r0"
    assert gdna_label_for_rate(0.3, None, 1) == "r0.3"
    assert gdna_label_for_rate(0.3, ["none", "low"], 1) == "low"


def test_manifest_round_trip(tmp_path):
    @dataclass
    class SimSection:
        frag_mean: int = 250

    class Config:
        simulation = SimSection()
        gdna = {"frag_mean": 350}
        nrna = {"modes": ["none"]}
        abundance = {"kind": "lognormal"}

    conditions = [{"name": "gdna_none_ss_1.00_nrna_none", "n_rna": 10}]
    index = {"path": "/idx", "manifest": {"sources": {"gtf": {"path": "/a.gtf", "sha256": "ab"}}}}

    path = write_manifest(tmp_path, Config(), conditions, index=index)
    loaded = json.loads(path.read_text())

    assert path == tmp_path / "manifest.json"
    assert loaded["simulation"]["frag_mean"] == 250
    assert loaded["conditions"] == conditions
    assert loaded["index"] == index


def test_the_condition_grid_names_each_truth_file_under_one_key(tmp_path):
    """A condition's manifest entry and its ``truth_summary.json`` record where its truth tables are, and
    each path is under one key, the same key in both files. The entry once named every path twice
    (``truth_abundances`` beside ``post_capture_abundances``, ...), the manifest a third time at its top
    level, and the summary repeated its paths under ``files``.

    Perturbation: a second key for a path in ``orchestrator.run_condition_grid``'s entry, a top-level
    copy in ``manifest.write_manifest``, or a summary key spelt otherwise."""
    from rigel.sim.capture import CaptureConfig, CaptureScenario
    from rigel.sim.genome import random_dna_array
    from rigel.sim.orchestrator import run_condition_grid
    from rigel.sim.wgs_config import GDNASimConfig, NRNAConfig, SimulationParams

    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "".join(random_dna_array(3000, np.random.default_rng(5))) + "\n")
    pysam.faidx(str(fasta))
    transcripts = _molar_transcripts()
    outdir = tmp_path / "panel"
    entries = run_condition_grid(
        outdir=outdir,
        genome_path=fasta,
        transcripts=transcripts,
        base_abundances=[(t.abundance, t.nrna_abundance) for t in transcripts],
        sim=SimulationParams(
            n_rna_fragments=200,
            frag_mean=80,
            frag_std=10,
            frag_min=50,
            frag_max=120,
            read_length=40,
        ),
        gdna=GDNASimConfig(frag_mean=80, frag_std=10, frag_min=50, frag_max=120),
        nrna=NRNAConfig(),
        genomic_refs=["chr1"],
        gdna_pairs=[("g50", 1.0)],
        gdna_od_pairs=[("od00", 0.0)],
        strand_specificities=[0.99],
        nrna_pairs=[("file", "file", None, 0)],
        capture_scenarios=[CaptureScenario("off", CaptureConfig())],
        include_capture_in_names=False,
        base_seed=3,
    )
    (entry,) = entries
    tables = ("pre_capture_abundances", "post_capture_abundances", "post_capture_fragment_lengths")
    paths = {key: value for key, value in entry.items() if str(value).endswith((".tsv", ".json"))}
    assert sorted(paths) == sorted([*tables, "truth_summary"])
    assert len(set(paths.values())) == len(paths), "two keys name one file"
    assert all((outdir / path).is_file() for path in paths.values())

    summary = json.loads((outdir / entry["truth_summary"]).read_text())
    assert summary["pre_capture_abundances"] == entry["pre_capture_abundances"]
    for key in tables[1:]:
        assert Path(summary[key]) == outdir / entry[key], key
    assert "files" not in summary

    manifest_path = write_manifest(
        outdir, type("Config", (), {})(), entries, index={"path": "/idx", "manifest": {}}
    )
    manifest = json.loads(manifest_path.read_text())
    assert not any(str(value).endswith(".tsv") for value in manifest.values()), (
        "the manifest names a truth table outside its condition's entry"
    )


def test_a_transcript_at_position_zero_has_its_whole_genomic_span():
    """The span was computed behind a truthiness test, so a transcript starting at base 0 recorded a span of 0."""
    from rigel.sim.truth import _transcript_metadata
    from rigel.transcript import Transcript
    from rigel.types import Interval

    t = Transcript(ref="chr1", exons=[Interval(0, 100), Interval(200, 300)], t_id="t1", g_id="g1")
    t.compute_length()
    assert _transcript_metadata(t)["genomic_span"] == 300
