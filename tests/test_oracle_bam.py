"""Tests for sim BAM coordinate helpers (rigel.sim.bam) + WholeGenomeSimulator oracle-BAM orientation."""

import pysam

from rigel.sim.bam import (
    blocks_to_cigar,
    premrna_to_genomic_interval,
    take_from_left,
    take_from_right,
    transcript_to_genomic_blocks,
)
from rigel.sim.genome import MutableGenome
from rigel.sim.whole_genome import GDNASimConfig, SimulationParams, WholeGenomeSimulator
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


class TestPremrnaToGenomicInterval:
    """Test premrna_to_genomic_interval."""

    def test_pos_strand(self):
        t = _make_pos_transcript([(100, 200), (300, 400)])
        # Pre-mRNA spans genomic [100, 400), len=300
        gstart, gend = premrna_to_genomic_interval(10, 50, t)
        assert gstart == 110
        assert gend == 150

    def test_neg_strand(self):
        t = _make_neg_transcript([(100, 200), (300, 400)])
        # Pre-mRNA spans [100, 400), len=300
        # Mirrored: start=300-50=250, end=300-10=290
        # Genomic: 100+250=350, 100+290=390
        gstart, gend = premrna_to_genomic_interval(10, 50, t)
        assert gstart == 350
        assert gend == 390


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


class TestTakeFromLeftRight:
    """Test take_from_left and take_from_right."""

    def test_take_left_single_block(self):
        result = take_from_left([(100, 300)], 50)
        assert result == [(100, 150)]

    def test_take_left_spanning(self):
        result = take_from_left([(100, 120), (200, 300)], 50)
        assert result == [(100, 120), (200, 230)]

    def test_take_right_single_block(self):
        result = take_from_right([(100, 300)], 50)
        assert result == [(250, 300)]

    def test_take_right_spanning(self):
        result = take_from_right([(100, 200), (300, 320)], 50)
        assert result == [(170, 200), (300, 320)]

    def test_take_left_exact(self):
        """Take exactly the whole block."""
        result = take_from_left([(100, 200)], 100)
        assert result == [(100, 200)]

    def test_take_right_exact(self):
        result = take_from_right([(100, 200)], 100)
        assert result == [(100, 200)]


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
            is_nrna=True,
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
            is_nrna=True,
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
            seed=19,
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
        seen_categories = set()
        with (
            pysam.FastaFile(str(fasta_path)) as fasta,
            pysam.AlignmentFile(str(bam_path), "rb") as bam,
        ):
            for read in bam.fetch(until_eof=True):
                if read.query_name.startswith("gdna:"):
                    seen_categories.add("gdna")
                elif read.query_name.startswith("nrna_"):
                    seen_categories.add("nrna")
                else:
                    seen_categories.add("mrna")
                assert read.query_sequence == _reference_sequence_for_record(read, fasta)

        assert seen_categories == {"mrna", "nrna", "gdna"}
