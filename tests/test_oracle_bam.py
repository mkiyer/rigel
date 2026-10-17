"""The simulator's coordinate helpers, and the orientation of the oracle BAM they write.

Projecting a transcript-space interval back to genomic blocks, projecting a pre-mRNA interval, the
CIGAR those blocks imply, and the take-from-left/right helpers that trim a fragment to its two
reads. Then the whole-genome simulator's oracle BAM, whose read names, mate orientation and strand
flags are what every truth-scoring instrument reads a fragment's origin from — an inverted pair there
is a truth source that disagrees with itself, and a name the engine writes other than as
:mod:`rigel.sim.read_name` defines it is one the parser misreads. Its header names the Rigel that
wrote it, and the FASTQ and the BAM give every base the one simulated quality.
"""

import gzip
from pathlib import Path

import pysam

import rigel
from rigel.sim.bam import (
    BASE_QUALITY,
    blocks_to_cigar,
    transcript_to_genomic_blocks,
)
from rigel.sim.genome import MutableGenome
from rigel.sim.read_name import (
    STRAND_FORWARD,
    STRAND_REVERSE,
    Origin,
    format_gdna_name,
    format_rna_name,
    parse_origin,
)
from rigel.sim.wgs_config import GDNASimConfig, SimulationParams
from rigel.sim.wgs_engine import WholeGenomeSimulator
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


# =====================================================================
# Helpers
# =====================================================================


def _make_pos_transcript(exons: list[tuple[int, int]], t_id="t1", g_id="g1") -> Transcript:
    """Create a POS-strand transcript with given exon coords."""
    t = Transcript(
        ref="ref",
        strand=Strand.POS,
        exons=[Interval(s, e) for s, e in exons],
        t_id=t_id,
        g_id=g_id,
        g_name=g_id,
        g_type="protein_coding",
        t_index=0,
        g_index=0,
        abundance=100.0,
    )
    t.compute_length()
    return t


def _make_neg_transcript(exons: list[tuple[int, int]], t_id="t2", g_id="g2") -> Transcript:
    """Create a NEG-strand transcript with given exon coords."""
    t = Transcript(
        ref="ref",
        strand=Strand.NEG,
        exons=[Interval(s, e) for s, e in exons],
        t_id=t_id,
        g_id=g_id,
        g_name=g_id,
        g_type="protein_coding",
        t_index=1,
        g_index=1,
        abundance=100.0,
    )
    t.compute_length()
    return t


def _reference_sequence_for_record(read: pysam.AlignedSegment, fasta: pysam.FastaFile) -> str:
    ref_pos = read.reference_start
    parts = []
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):
            parts.append(fasta.fetch(read.reference_name, ref_pos, ref_pos + length).upper())
            ref_pos += length
        elif op in (2, 3):
            ref_pos += length
    return "".join(parts)


def _name_from_origin(origin: Origin) -> str:
    """The name ``read_name``'s format functions build from a parsed origin's fields."""
    if origin.kind == "gdna":
        return format_gdna_name(origin.ref, origin.start, origin.end, origin.strand, origin.index)
    return format_rna_name(
        origin.transcript_id,
        origin.start,
        origin.end,
        origin.strand,
        origin.index,
        nascent=origin.kind == "nrna",
    )


# =====================================================================
# Coordinate projection tests
# =====================================================================


class TestTranscriptToGenomicBlocks:
    """Test transcript_to_genomic_blocks for POS and NEG strand."""

    def test_single_exon_pos(self):
        """Fragment within a single exon on + strand."""
        t = _make_pos_transcript([(100, 300)])  # 200bp exon
        blocks = transcript_to_genomic_blocks(10, 50, t)
        assert blocks == [(110, 150)]

    def test_single_exon_neg(self):
        """Fragment within a single exon on − strand.

        NEG strand: tx pos 0 = rightmost base of the exon.
        Fragment at tx [10, 50) on a 200bp exon:
        Mirrored: tx [150, 190) → genomic [250, 290)
        """
        t = _make_neg_transcript([(100, 300)])
        blocks = transcript_to_genomic_blocks(10, 50, t)
        # t_len=200, mirrored: start=200-50=150, end=200-10=190
        # exon is [100,300), offset 150→190 → genomic [250, 290)
        assert blocks == [(250, 290)]

    def test_spanning_intron_pos(self):
        """Fragment spanning an intron on + strand."""
        t = _make_pos_transcript([(100, 200), (300, 400)])
        # Transcript: 100bp exon1 + 100bp exon2 = 200bp total
        # Fragment at tx [80, 120) spans exon boundary at pos 100
        blocks = transcript_to_genomic_blocks(80, 120, t)
        # First 20bp in exon1: genomic [180, 200)
        # Next 20bp in exon2: genomic [300, 320)
        assert blocks == [(180, 200), (300, 320)]

    def test_spanning_intron_neg(self):
        """Fragment spanning an intron on − strand."""
        t = _make_neg_transcript([(100, 200), (300, 400)])
        # t_len = 200bp. Fragment at tx [80, 120):
        # Mirrored: [80, 120) → [80, 120) in genomic-ascending tx coords
        # (mirrored: start=200-120=80, end=200-80=120)
        # tx pos 80-99 → exon1 offset 80 → genomic [180, 200)
        # tx pos 100-119 → exon2 offset 0 → genomic [300, 320)
        blocks = transcript_to_genomic_blocks(80, 120, t)
        assert blocks == [(180, 200), (300, 320)]

    def test_within_first_exon_pos(self):
        """Fragment entirely within the first exon."""
        t = _make_pos_transcript([(100, 200), (300, 400)])
        blocks = transcript_to_genomic_blocks(0, 50, t)
        assert blocks == [(100, 150)]

    def test_within_second_exon_pos(self):
        """Fragment entirely within the second exon."""
        t = _make_pos_transcript([(100, 200), (300, 400)])
        # tx position 100 = start of exon2
        blocks = transcript_to_genomic_blocks(100, 150, t)
        assert blocks == [(300, 350)]

    def test_three_exon_spanning_two_introns(self):
        """Fragment covering parts of all three exons."""
        t = _make_pos_transcript([(0, 50), (100, 150), (200, 250)])
        # t_len = 150. Fragment [40, 120):
        # exon1 [0,50): tx 40-49 → genomic [40, 50)
        # exon2 [100,150): tx 50-99 → genomic [100, 150)
        # exon3 [200,250): tx 100-119 → genomic [200, 220)
        blocks = transcript_to_genomic_blocks(40, 120, t)
        assert blocks == [(40, 50), (100, 150), (200, 220)]


class TestBlocksToCigar:
    """Test blocks_to_cigar."""

    def test_single_block(self):
        cigar = blocks_to_cigar([(100, 250)])
        assert cigar == [(0, 150)]  # 150M

    def test_two_blocks_with_intron(self):
        cigar = blocks_to_cigar([(100, 200), (300, 400)])
        # 100M, 100N, 100M
        assert cigar == [(0, 100), (3, 100), (0, 100)]

    def test_three_blocks(self):
        cigar = blocks_to_cigar([(10, 50), (100, 150), (200, 220)])
        assert cigar == [(0, 40), (3, 50), (0, 50), (3, 50), (0, 20)]


class TestWholeGenomeOracleBamOrientation:
    """Regression tests for the vectorized whole-genome oracle BAM writer."""

    def test_bam_seq_matches_reference_for_reverse_flipped_and_gdna_reads(self, tmp_path):
        genome = MutableGenome(length=3000, seed=7, name="ref")
        transcripts = [
            _make_pos_transcript([(200, 500), (700, 1000)], t_id="tx_pos", g_id="g_pos"),
            _make_neg_transcript([(1400, 1700), (1900, 2200)], t_id="tx_neg", g_id="g_neg"),
        ]
        # nascent lives on single-exon ENTITIES over each span (one per strand here), as the index
        # makes them; the entities' reads carry the `nrna_` tag and contiguous genomic blocks
        pos_entity = Transcript(
            ref="ref",
            strand=Strand.POS,
            exons=[Interval(200, 1000)],
            t_id="n_pos",
            g_id="n_pos",
            t_index=2,
            is_single_exon=True,
            is_synthetic=True,
            abundance=0.0,
            nrna_abundance=100.0,
        )
        neg_entity = Transcript(
            ref="ref",
            strand=Strand.NEG,
            exons=[Interval(1400, 2200)],
            t_id="n_neg",
            g_id="n_neg",
            t_index=3,
            is_single_exon=True,
            is_synthetic=True,
            abundance=0.0,
            nrna_abundance=100.0,
        )
        for t in (pos_entity, neg_entity):
            t.length = t.compute_length()
        transcripts += [pos_entity, neg_entity]

        fasta_path = genome.write_fasta(tmp_path)
        sim = WholeGenomeSimulator(
            fasta_path,
            transcripts,
            SimulationParams(
                sim_seed=19,
                frag_mean=120,
                frag_std=1,
                frag_min=120,
                frag_max=120,
                read_length=50,
                error_rate=0.0,
            ),
            GDNASimConfig(frag_mean=120, frag_std=1, frag_min=120, frag_max=120),
            genomic_refs=[genome.name],
            strand_specificity=0.5,
        )
        try:
            _, _, bam_path = sim.simulate_and_write(
                tmp_path / "out",
                n_rna=160,
                n_gdna=80,
                oracle_bam=True,
                n_workers=1,
            )
        finally:
            sim.close()

        assert bam_path is not None
        strand_letter = {
            t.t_id: STRAND_REVERSE if t.strand == Strand.NEG else STRAND_FORWARD
            for t in transcripts
        }
        seen_categories = set()
        with (
            pysam.FastaFile(str(fasta_path)) as fasta,
            pysam.AlignmentFile(str(bam_path), "rb") as bam,
        ):
            for read in bam.fetch(until_eof=True):
                origin = parse_origin(read.query_name)
                seen_categories.add(origin.kind)
                # the engine writes `read_name`'s names, field for field, with the strand letter the
                # format defines: the transcript's strand for RNA, the fragment's genomic strand for
                # gDNA (R1 reverse-aligned on `f`, forward on `r`)
                assert _name_from_origin(origin) == read.query_name
                if origin.kind == "gdna":
                    assert read.is_reverse == (read.is_read1 == (origin.strand == STRAND_FORWARD))
                else:
                    assert origin.strand == strand_letter[origin.transcript_id]
                assert read.query_sequence == _reference_sequence_for_record(read, fasta)

        assert seen_categories == {"mrna", "nrna", "gdna"}


class TestWholeGenomeOracleBamHeader:
    """What the oracle BAM says about itself, and the quality both writers give a base."""

    def test_the_header_names_the_rigel_version_and_every_base_has_the_one_quality(self, tmp_path):
        """The @PG version is the installed Rigel's, not a number of the simulator's own; the FASTQ writer and
        the BAM writer read ``BASE_QUALITY`` from one definition, so they cannot disagree."""
        genome = MutableGenome(length=2000, seed=3, name="ref")
        fasta_path = genome.write_fasta(tmp_path)
        sim = WholeGenomeSimulator(
            fasta_path,
            [_make_pos_transcript([(200, 500), (700, 1000)])],
            SimulationParams(
                sim_seed=5,
                frag_mean=150,
                frag_std=1,
                frag_min=150,
                frag_max=150,
                read_length=50,
                error_rate=0.0,
            ),
            GDNASimConfig(frag_mean=150, frag_std=1, frag_min=150, frag_max=150),
            genomic_refs=[genome.name],
        )
        try:
            r1_path, r2_path, bam_path = sim.simulate_and_write(
                tmp_path / "out", n_rna=20, n_gdna=10, oracle_bam=True, n_workers=1
            )
        finally:
            sim.close()

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            (pg,) = bam.header.to_dict()["PG"]
            bam_qualities = {q for read in bam.fetch(until_eof=True) for q in read.query_qualities}
        assert pg["VN"] == rigel.__version__
        assert bam_qualities == set(pysam.qualitystring_to_array(BASE_QUALITY))
        for path in (r1_path, r2_path):
            lines = gzip.decompress(Path(path).read_bytes()).decode().splitlines()
            assert lines and {c for quality in lines[3::4] for c in quality} == {BASE_QUALITY}
