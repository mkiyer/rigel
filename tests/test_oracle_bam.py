"""Tests for rigel.sim.oracle_bam — Oracle BAM fragment simulator."""

import numpy as np
import pysam
import pytest
from pathlib import Path

from rigel.sim.oracle_bam import (
    OracleBamSimulator,
    _blocks_to_cigar,
    _premrna_to_genomic_interval,
    _take_from_left,
    _take_from_right,
    _transcript_to_genomic_blocks,
)
from rigel.sim.genome import MutableGenome
from rigel.sim.reads import GDNAConfig, ReadSimulator, SimConfig
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


# =====================================================================
# Coordinate projection tests
# =====================================================================


class TestTranscriptToGenomicBlocks:
    """Test _transcript_to_genomic_blocks for POS and NEG strand."""

    def test_single_exon_pos(self):
        """Fragment within a single exon on + strand."""
        t = _make_pos_transcript([(100, 300)])  # 200bp exon
        blocks = _transcript_to_genomic_blocks(10, 50, t)
        assert blocks == [(110, 150)]

    def test_single_exon_neg(self):
        """Fragment within a single exon on − strand.

        NEG strand: tx pos 0 = rightmost base of the exon.
        Fragment at tx [10, 50) on a 200bp exon:
        Mirrored: tx [150, 190) → genomic [250, 290)
        """
        t = _make_neg_transcript([(100, 300)])
        blocks = _transcript_to_genomic_blocks(10, 50, t)
        # t_len=200, mirrored: start=200-50=150, end=200-10=190
        # exon is [100,300), offset 150→190 → genomic [250, 290)
        assert blocks == [(250, 290)]

    def test_spanning_intron_pos(self):
        """Fragment spanning an intron on + strand."""
        t = _make_pos_transcript([(100, 200), (300, 400)])
        # Transcript: 100bp exon1 + 100bp exon2 = 200bp total
        # Fragment at tx [80, 120) spans exon boundary at pos 100
        blocks = _transcript_to_genomic_blocks(80, 120, t)
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
        blocks = _transcript_to_genomic_blocks(80, 120, t)
        assert blocks == [(180, 200), (300, 320)]

    def test_within_first_exon_pos(self):
        """Fragment entirely within the first exon."""
        t = _make_pos_transcript([(100, 200), (300, 400)])
        blocks = _transcript_to_genomic_blocks(0, 50, t)
        assert blocks == [(100, 150)]

    def test_within_second_exon_pos(self):
        """Fragment entirely within the second exon."""
        t = _make_pos_transcript([(100, 200), (300, 400)])
        # tx position 100 = start of exon2
        blocks = _transcript_to_genomic_blocks(100, 150, t)
        assert blocks == [(300, 350)]

    def test_three_exon_spanning_two_introns(self):
        """Fragment covering parts of all three exons."""
        t = _make_pos_transcript([(0, 50), (100, 150), (200, 250)])
        # t_len = 150. Fragment [40, 120):
        # exon1 [0,50): tx 40-49 → genomic [40, 50)
        # exon2 [100,150): tx 50-99 → genomic [100, 150)
        # exon3 [200,250): tx 100-119 → genomic [200, 220)
        blocks = _transcript_to_genomic_blocks(40, 120, t)
        assert blocks == [(40, 50), (100, 150), (200, 220)]


class TestPremrnaToGenomicInterval:
    """Test _premrna_to_genomic_interval."""

    def test_pos_strand(self):
        t = _make_pos_transcript([(100, 200), (300, 400)])
        # Pre-mRNA spans genomic [100, 400), len=300
        gstart, gend = _premrna_to_genomic_interval(10, 50, t)
        assert gstart == 110
        assert gend == 150

    def test_neg_strand(self):
        t = _make_neg_transcript([(100, 200), (300, 400)])
        # Pre-mRNA spans [100, 400), len=300
        # Mirrored: start=300-50=250, end=300-10=290
        # Genomic: 100+250=350, 100+290=390
        gstart, gend = _premrna_to_genomic_interval(10, 50, t)
        assert gstart == 350
        assert gend == 390


class TestBlocksToCigar:
    """Test _blocks_to_cigar."""

    def test_single_block(self):
        cigar = _blocks_to_cigar([(100, 250)])
        assert cigar == [(0, 150)]  # 150M

    def test_two_blocks_with_intron(self):
        cigar = _blocks_to_cigar([(100, 200), (300, 400)])
        # 100M, 100N, 100M
        assert cigar == [(0, 100), (3, 100), (0, 100)]

    def test_three_blocks(self):
        cigar = _blocks_to_cigar([(10, 50), (100, 150), (200, 220)])
        assert cigar == [(0, 40), (3, 50), (0, 50), (3, 50), (0, 20)]


class TestTakeFromLeftRight:
    """Test _take_from_left and _take_from_right."""

    def test_take_left_single_block(self):
        result = _take_from_left([(100, 300)], 50)
        assert result == [(100, 150)]

    def test_take_left_spanning(self):
        result = _take_from_left([(100, 120), (200, 300)], 50)
        assert result == [(100, 120), (200, 230)]

    def test_take_right_single_block(self):
        result = _take_from_right([(100, 300)], 50)
        assert result == [(250, 300)]

    def test_take_right_spanning(self):
        result = _take_from_right([(100, 200), (300, 320)], 50)
        assert result == [(170, 200), (300, 320)]

    def test_take_left_exact(self):
        """Take exactly the whole block."""
        result = _take_from_left([(100, 200)], 100)
        assert result == [(100, 200)]

    def test_take_right_exact(self):
        result = _take_from_right([(100, 200)], 100)
        assert result == [(100, 200)]


# =====================================================================
# OracleBamSimulator integration tests
# =====================================================================


class TestOracleBamSimulatorBasic:
    """Test OracleBamSimulator writes valid BAM files."""

    @pytest.fixture
    def simple_setup(self):
        """Create a simple genome + transcript for testing."""
        genome = MutableGenome(length=2000, seed=42, name="ref")
        # Inject known splice motifs
        genome.edit(200, "GT")  # donor
        genome.edit(298, "AG")  # acceptor
        genome.edit(400, "GT")  # donor
        genome.edit(498, "AG")  # acceptor

        transcripts = [
            _make_pos_transcript(
                [(100, 200), (300, 400), (500, 600)],
                t_id="tx_pos", g_id="gene_pos",
            ),
            _make_neg_transcript(
                [(800, 900), (1000, 1100)],
                t_id="tx_neg", g_id="gene_neg",
            ),
        ]
        return genome, transcripts

    def test_write_bam_creates_file(self, simple_setup, tmp_path):
        genome, transcripts = simple_setup
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        result = sim.write_bam(bam_path, n_fragments=100)
        assert result.exists()

    def test_bam_is_valid(self, simple_setup, tmp_path):
        """BAM should be readable by pysam with proper paired-end records."""
        genome, transcripts = simple_setup
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=100)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            records = list(bam.fetch(until_eof=True))
            assert len(records) > 0
            # Should be even (paired)
            assert len(records) % 2 == 0

    def test_paired_end_flags(self, simple_setup, tmp_path):
        """All records should have proper paired-end flags."""
        genome, transcripts = simple_setup
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=50)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                # Must be paired
                assert read.flag & 0x1
                # Must be proper pair
                assert read.flag & 0x2
                # Must be R1 or R2
                is_r1 = bool(read.flag & 0x40)
                is_r2 = bool(read.flag & 0x80)
                assert is_r1 ^ is_r2  # exactly one
                # NH tag should be 1
                assert read.get_tag("NH") == 1

    def test_name_sorted(self, simple_setup, tmp_path):
        """Default output should be name-sorted."""
        genome, transcripts = simple_setup
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=50)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            names = [r.query_name for r in bam.fetch(until_eof=True)]
            # Name-sorted: each name appears twice (R1+R2) consecutively
            for i in range(0, len(names), 2):
                assert names[i] == names[i + 1]

    def test_spliced_reads_have_xs_tag(self, simple_setup, tmp_path):
        """Spliced mRNA reads should have the XS tag."""
        genome, transcripts = simple_setup
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=200)

        found_spliced = False
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                qname = read.query_name
                # A read from tx_pos (3 exons) might have N in CIGAR
                cigar_ops = [op for op, _ in read.cigartuples] if read.cigartuples else []
                if 3 in cigar_ops:  # BAM_CREF_SKIP
                    found_spliced = True
                    assert read.has_tag("XS")
                    xs = read.get_tag("XS")
                    if qname.startswith("tx_pos"):
                        assert xs == "+"
                    elif qname.startswith("tx_neg"):
                        assert xs == "-"

        assert found_spliced, "Expected at least one spliced read"

    def test_read_name_format(self, simple_setup, tmp_path):
        """Read names should encode ground truth."""
        genome, transcripts = simple_setup
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=50)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                qname = read.query_name
                # Format: {tid}:{start}-{end}:{strand}:{idx}
                parts = qname.split(":")
                assert len(parts) >= 3
                tid = parts[0]
                assert tid in ("tx_pos", "tx_neg")


class TestOracleBamSimulatorWithGDNA:
    """Test OracleBamSimulator with gDNA contamination."""

    def test_gdna_fragments_present(self, tmp_path):
        genome = MutableGenome(length=5000, seed=42, name="ref")
        transcripts = [
            _make_pos_transcript(
                [(100, 300), (500, 700)],
                t_id="tx1", g_id="g1",
            ),
        ]
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        gdna_config = GDNAConfig(abundance=100.0)

        sim = OracleBamSimulator(
            genome, transcripts, config=config,
            gdna_config=gdna_config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=500)

        gdna_count = 0
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.query_name.startswith("gdna:"):
                    gdna_count += 1
        assert gdna_count > 0, "Expected some gDNA reads"

    def test_gdna_reads_unspliced(self, tmp_path):
        """gDNA reads should have simple M CIGAR (no intron skips)."""
        genome = MutableGenome(length=5000, seed=42, name="ref")
        transcripts = [
            _make_pos_transcript(
                [(100, 300)], t_id="tx1", g_id="g1",
            ),
        ]
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, seed=42,
        )
        gdna_config = GDNAConfig(abundance=200.0)

        sim = OracleBamSimulator(
            genome, transcripts, config=config,
            gdna_config=gdna_config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=500)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.query_name.startswith("gdna:"):
                    # Should be simple M only
                    for op, _ in read.cigartuples:
                        assert op == 0, "gDNA reads should only have M (match) CIGAR"


class TestOracleBamSimulatorWithNRNA:
    """Test OracleBamSimulator with nascent RNA."""

    def test_nrna_fragments_present(self, tmp_path):
        genome = MutableGenome(length=5000, seed=42, name="ref")
        transcripts = [
            _make_pos_transcript(
                [(100, 300), (500, 700)],
                t_id="tx1", g_id="g1",
            ),
        ]
        # Set nRNA abundance
        transcripts[0].nrna_abundance = 200.0

        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=500)

        nrna_count = 0
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.query_name.startswith("nrna_"):
                    nrna_count += 1
        assert nrna_count > 0, "Expected some nRNA reads"


class TestOracleBamSimulatorStrandSpecificity:
    """Test that strand specificity affects read orientations."""

    def test_perfect_strandedness(self, tmp_path):
        """With strand_specificity=1.0, all POS-strand tx reads should
        have R1 on − and R2 on +."""
        genome = MutableGenome(length=2000, seed=42, name="ref")
        transcripts = [
            _make_pos_transcript([(100, 400)], t_id="tx1", g_id="g1"),
        ]
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=350,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=100)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.flag & 0x40:  # R1
                    assert read.is_reverse, "R1 should be reverse for POS tx"
                else:  # R2
                    assert not read.is_reverse, "R2 should be forward for POS tx"

    def test_imperfect_strandedness_has_flipped(self, tmp_path):
        """With strand_specificity=0.5, roughly half should be flipped."""
        genome = MutableGenome(length=2000, seed=42, name="ref")
        transcripts = [
            _make_pos_transcript([(100, 400)], t_id="tx1", g_id="g1"),
        ]
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=350,
            read_length=100, strand_specificity=0.5, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=200)

        r1_reverse = 0
        r1_forward = 0
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.flag & 0x40:  # R1
                    if read.is_reverse:
                        r1_reverse += 1
                    else:
                        r1_forward += 1

        total = r1_reverse + r1_forward
        assert total > 0
        # At 0.5 specificity, roughly 50% should be each way
        ratio = r1_reverse / total
        assert 0.2 < ratio < 0.8, f"Expected ~50% reverse for R1, got {ratio:.2f}"


class TestDirectBamCoordSorted:
    """Test coordinate-sorted BAM output."""

    def test_coord_sorted_and_indexed(self, tmp_path):
        genome = MutableGenome(length=2000, seed=42, name="ref")
        transcripts = [
            _make_pos_transcript([(100, 400)], t_id="tx1", g_id="g1"),
        ]
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=50, frag_max=350,
            read_length=100, seed=42,
        )
        sim = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "test.bam"
        sim.write_bam(bam_path, n_fragments=50, name_sorted=False)

        assert bam_path.exists()
        bai = Path(str(bam_path) + ".bai")
        assert bai.exists(), "Index file should be created for coord-sorted BAM"

        # Verify coordinate order
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            positions = [r.reference_start for r in bam.fetch(until_eof=True)]
            assert positions == sorted(positions)


class TestOracleBamTruthCounts:
    """Test that truth counts match between OracleBam and FASTQ simulator."""

    def test_truth_counts_consistency(self, tmp_path):
        """OracleBamSimulator and ReadSimulator with same seed should
        produce the same number of fragments per transcript."""
        genome = MutableGenome(length=3000, seed=42, name="ref")
        transcripts = [
            _make_pos_transcript(
                [(100, 300), (500, 700)],
                t_id="tx1", g_id="g1",
            ),
            _make_neg_transcript(
                [(1000, 1200), (1400, 1600)],
                t_id="tx2", g_id="g2",
            ),
        ]
        n_frags = 500
        config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=50, frag_max=500,
            read_length=100, strand_specificity=1.0, seed=42,
        )

        # Oracle BAM
        sim1 = OracleBamSimulator(
            genome, transcripts, config=config, ref_name="ref",
        )
        bam_path = tmp_path / "oracle.bam"
        sim1.write_bam(bam_path, n_frags)

        # Count from BAM R1 reads
        bam_counts = {}
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if not (read.flag & 0x40):
                    continue
                tid = read.query_name.split(":")[0]
                bam_counts[tid] = bam_counts.get(tid, 0) + 1

        # FASTQ simulator (same seed)
        fastq_sim = ReadSimulator(
            genome, transcripts, config=config,
        )
        fastq_counts = {}
        for r1_name, *_ in fastq_sim.simulate(n_frags):
            qname = r1_name.rstrip("/1")
            tid = qname.split(":")[0]
            fastq_counts[tid] = fastq_counts.get(tid, 0) + 1

        # Same transcript IDs should appear
        assert set(bam_counts.keys()) == set(fastq_counts.keys())
        # Same total fragments
        assert sum(bam_counts.values()) == sum(fastq_counts.values())
        # Same per-transcript counts
        for tid in bam_counts:
            assert bam_counts[tid] == fastq_counts[tid], (
                f"Count mismatch for {tid}: BAM={bam_counts[tid]} vs "
                f"FASTQ={fastq_counts[tid]}"
            )
