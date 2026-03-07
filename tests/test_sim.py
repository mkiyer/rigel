"""Tests for rigel.sim — simulation framework components."""

import numpy as np
import pytest

from rigel.sim.genome import MutableGenome, reverse_complement
from rigel.sim.annotation import GeneBuilder
from rigel.sim.reads import ReadSimulator, SimConfig
from rigel.transcript import Transcript
from rigel.types import Strand, Interval


# =====================================================================
# MutableGenome
# =====================================================================


class TestMutableGenome:
    def test_length(self):
        g = MutableGenome(1000, seed=1, name="chr1")
        assert len(g) == 1000

    def test_deterministic_with_seed(self):
        g1 = MutableGenome(500, seed=42, name="chr1")
        g2 = MutableGenome(500, seed=42, name="chr1")
        assert g1.seq == g2.seq

    def test_different_seeds_differ(self):
        g1 = MutableGenome(500, seed=1, name="chr1")
        g2 = MutableGenome(500, seed=2, name="chr1")
        assert g1.seq != g2.seq

    def test_only_acgt(self):
        g = MutableGenome(10000, seed=7, name="chr1")
        assert set(g.seq) <= {"A", "C", "G", "T"}

    def test_getitem_int(self):
        g = MutableGenome(100, seed=1, name="chr1")
        assert g[0] in "ACGT"

    def test_getitem_slice(self):
        g = MutableGenome(100, seed=1, name="chr1")
        subseq = g[10:20]
        assert len(subseq) == 10
        assert subseq == g.seq[10:20]

    def test_edit(self):
        g = MutableGenome(100, seed=1, name="chr1")
        g.edit(50, "GATTACA")
        assert g[50:57] == "GATTACA"

    def test_edit_boundary(self):
        g = MutableGenome(100, seed=1, name="chr1")
        g.edit(98, "AG")
        assert g[98:100] == "AG"

    def test_edit_out_of_bounds_raises(self):
        g = MutableGenome(100, seed=1, name="chr1")
        with pytest.raises(IndexError):
            g.edit(99, "AG")  # extends to 101

    def test_edit_negative_raises(self):
        g = MutableGenome(100, seed=1, name="chr1")
        with pytest.raises(IndexError):
            g.edit(-1, "A")

    def test_write_fasta(self, tmp_path):
        g = MutableGenome(500, seed=1, name="testchr")
        fasta = g.write_fasta(tmp_path)
        assert fasta.exists()
        assert fasta.name == "testchr.fa"
        # Check .fai was created
        assert (tmp_path / "testchr.fa.fai").exists()

    def test_write_fasta_readable_by_pysam(self, tmp_path):
        import pysam
        g = MutableGenome(500, seed=1, name="testchr")
        fasta = g.write_fasta(tmp_path)
        with pysam.FastaFile(str(fasta)) as fh:
            assert fh.nreferences == 1
            assert fh.references[0] == "testchr"
            assert fh.lengths[0] == 500
            fetched = fh.fetch("testchr", 0, 500)
            assert fetched == g.seq


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ACGT") == "ACGT"

    def test_single(self):
        assert reverse_complement("A") == "T"
        assert reverse_complement("G") == "C"

    def test_palindrome(self):
        assert reverse_complement("AATT") == "AATT"

    def test_asymmetric(self):
        assert reverse_complement("AACG") == "CGTT"

    def test_lowercase(self):
        assert reverse_complement("acgt") == "acgt"


# =====================================================================
# GeneBuilder
# =====================================================================


class TestGeneBuilder:
    def test_single_gene_single_transcript(self, tmp_path):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 300), (500, 700)], "abundance": 100},
        ])
        transcripts = builder.get_transcripts()
        assert len(transcripts) == 1
        assert transcripts[0].t_id == "t1"
        assert transcripts[0].g_id == "g1"
        assert transcripts[0].strand == Strand.POS
        assert len(transcripts[0].exons) == 2

    def test_splice_motif_positive_strand(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
        ])
        # Intron is (300, 500): donor at 300, acceptor at 498
        assert g[300:302] == "GT"
        assert g[498:500] == "AG"

    def test_splice_motif_negative_strand(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "-", [
            {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
        ])
        # Negative strand: CT at donor, AC at acceptor
        assert g[300:302] == "CT"
        assert g[498:500] == "AC"

    def test_multi_exon_splice_motifs(self):
        g = MutableGenome(3000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 200), (400, 500), (700, 800)]},
        ])
        # Intron 1: (200, 400)
        assert g[200:202] == "GT"
        assert g[398:400] == "AG"
        # Intron 2: (500, 700)
        assert g[500:502] == "GT"
        assert g[698:700] == "AG"

    def test_multi_isoform_gene(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 200), (400, 500), (700, 800)]},
            {"t_id": "t2", "exons": [(100, 200), (700, 800)]},
        ])
        transcripts = builder.get_transcripts()
        assert len(transcripts) == 2
        # Both belong to same gene
        assert all(t.g_id == "g1" for t in transcripts)
        assert all(t.g_index == 0 for t in transcripts)
        # Different t_index
        assert transcripts[0].t_index != transcripts[1].t_index

    def test_overlapping_exons_raises(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        with pytest.raises(ValueError, match="overlap"):
            builder.add_gene("g1", "+", [
                {"t_id": "t1", "exons": [(100, 300), (250, 500)]},
            ])

    def test_exon_out_of_bounds_raises(self):
        g = MutableGenome(500, seed=1, name="chr1")
        builder = GeneBuilder(g)
        with pytest.raises(ValueError, match="outside genome"):
            builder.add_gene("g1", "+", [
                {"t_id": "t1", "exons": [(100, 600)]},
            ])

    def test_t_index_assignment(self):
        g = MutableGenome(3000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 200)]},
        ])
        builder.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(1000, 1200)]},
        ])
        transcripts = builder.get_transcripts()
        assert transcripts[0].t_index == 0
        assert transcripts[1].t_index == 1
        assert transcripts[0].g_index != transcripts[1].g_index

    def test_write_gtf(self, tmp_path):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
        ])
        gtf_path = builder.write_gtf(tmp_path)
        assert gtf_path.exists()

        # Read it back with the GTF parser
        from rigel.gtf import GTFRecord
        features = list(GTFRecord.parse_file(gtf_path))
        assert len(features) == 2  # 2 exon lines
        assert all(f.feature == "exon" for f in features)
        assert all(f.attrs["gene_id"] == "g1" for f in features)
        assert all(f.attrs["transcript_id"] == "t1" for f in features)

    def test_write_gtf_roundtrip(self, tmp_path):
        """GTF → Transcript.read_gtf should reconstruct the transcripts."""
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 300), (500, 700)], "abundance": 50.0},
            {"t_id": "t2", "exons": [(100, 300)], "abundance": 25.0},
        ])
        gtf_path = builder.write_gtf(tmp_path)

        # Read back
        loaded = Transcript.read_gtf(str(gtf_path))
        assert len(loaded) == 2

        t1 = next(t for t in loaded if t.t_id == "t1")
        assert len(t1.exons) == 2
        assert t1.exons[0] == Interval(100, 300)
        assert t1.exons[1] == Interval(500, 700)
        assert t1.strand == Strand.POS
        assert t1.g_id == "g1"

    def test_abundance_in_gtf(self, tmp_path):
        """Abundance is written as GTF score field."""
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 300)], "abundance": 42.5},
        ])
        gtf_path = builder.write_gtf(tmp_path)

        from rigel.gtf import GTFRecord
        features = list(GTFRecord.parse_file(gtf_path))
        assert features[0].score == 42.5

    def test_intron_too_short_raises(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        with pytest.raises(ValueError, match="too short"):
            builder.add_gene("g1", "+", [
                {"t_id": "t1", "exons": [(100, 200), (202, 300)]},
            ])

    def test_gtf_to_bed12(self, tmp_path):
        """GTF → BED12 conversion produces valid 12-column BED."""
        from rigel.index import gtf_to_bed12

        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
        ])
        builder.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(1000, 1200), (1400, 1500), (1700, 1800)]},
        ])
        gtf_path = builder.write_gtf(tmp_path)
        bed_path = tmp_path / "annotation.bed"
        result_path = gtf_to_bed12(gtf_path, bed_path)
        assert result_path.exists()

        lines = bed_path.read_text().strip().split("\n")
        assert len(lines) == 2

        # Parse first transcript (2-exon, + strand)
        fields = lines[0].split("\t")
        assert len(fields) == 12
        assert fields[0] == "chr1"     # chrom
        assert fields[1] == "100"      # chromStart
        assert fields[2] == "700"      # chromEnd
        assert fields[3] == "t1"       # name
        assert fields[5] == "+"        # strand
        assert fields[9] == "2"        # blockCount
        assert fields[10] == "200,200" # blockSizes (300-100, 700-500)
        assert fields[11] == "0,400"   # blockStarts (100-100, 500-100)

        # Parse second transcript (3-exon, - strand)
        fields = lines[1].split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "1000"
        assert fields[2] == "1800"
        assert fields[3] == "t2"
        assert fields[5] == "-"
        assert fields[9] == "3"
        assert fields[10] == "200,100,100"  # sizes
        assert fields[11] == "0,400,700"    # starts relative to 1000


# =====================================================================
# ReadSimulator
# =====================================================================


class TestReadSimulator:
    def _make_sim(self, tmp_path=None, seed=42, abundance=100.0, n_exons=2):
        """Build a simple simulator with one transcript."""
        g = MutableGenome(3000, seed=seed, name="chr1")
        builder = GeneBuilder(g)
        if n_exons == 2:
            exons = [(200, 600), (1000, 1400)]
        elif n_exons == 1:
            exons = [(200, 1200)]
        else:
            exons = [(200, 400), (600, 800), (1000, 1200)]
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": exons, "abundance": abundance},
        ])
        transcripts = builder.get_transcripts()
        config = SimConfig(
            frag_mean=250, frag_std=50, frag_min=100,
            frag_max=700, read_length=100, error_rate=0.0, seed=seed,
        )
        sim = ReadSimulator(g, transcripts, config=config)
        return g, transcripts, sim

    def test_correct_read_count(self):
        _, _, sim = self._make_sim()
        reads = list(sim.simulate(100))
        assert len(reads) == 100

    def test_read_name_format(self):
        _, _, sim = self._make_sim()
        reads = list(sim.simulate(10))
        r1_name = reads[0][0]
        r2_name = reads[0][3]
        assert r1_name.endswith("/1")
        assert r2_name.endswith("/2")
        # Ground truth encoded: t_id:start-end:strand:idx/read
        parts = r1_name.rsplit("/", 1)[0]
        fields = parts.split(":")
        assert fields[0] == "t1"  # transcript id
        assert "-" in fields[1]   # frag_start-frag_end
        assert fields[2] in ("f", "r")  # strand char

    def test_read_length(self):
        _, _, sim = self._make_sim()
        reads = list(sim.simulate(20))
        for r1_name, r1_seq, r1_qual, _, r2_seq, r2_qual in reads:
            assert len(r1_seq) == 100
            assert len(r2_seq) == 100
            assert len(r1_qual) == 100
            assert len(r2_qual) == 100

    def test_sequences_are_dna(self):
        _, _, sim = self._make_sim()
        reads = list(sim.simulate(50))
        for _, r1_seq, _, _, r2_seq, _ in reads:
            assert set(r1_seq) <= {"A", "C", "G", "T"}
            assert set(r2_seq) <= {"A", "C", "G", "T"}

    def test_deterministic_with_seed(self):
        _, _, sim1 = self._make_sim(seed=42)
        _, _, sim2 = self._make_sim(seed=42)
        reads1 = list(sim1.simulate(50))
        reads2 = list(sim2.simulate(50))
        for r1, r2 in zip(reads1, reads2):
            assert r1[1] == r2[1]  # R1 sequences match
            assert r1[4] == r2[4]  # R2 sequences match

    def test_write_fastq(self, tmp_path):
        _, _, sim = self._make_sim()
        r1_path, r2_path = sim.write_fastq(tmp_path, 50)
        assert r1_path.exists()
        assert r2_path.exists()

        # Count records
        r1_lines = r1_path.read_text().strip().split("\n")
        r2_lines = r2_path.read_text().strip().split("\n")
        assert len(r1_lines) == 200  # 50 records × 4 lines
        assert len(r2_lines) == 200

    def test_error_injection(self):
        """With a high error rate, not all reads should match exactly."""
        g = MutableGenome(3000, seed=42, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 1200)]},
        ])
        transcripts = builder.get_transcripts()
        config = SimConfig(
            frag_mean=200, frag_std=20, frag_min=100,
            frag_max=500, read_length=100,
            error_rate=0.5,  # Very high error rate
            seed=42,
        )
        sim = ReadSimulator(g, transcripts, config=config)
        reads = list(sim.simulate(100))
        # With 50% error rate, bases will differ from expected
        # Just verify reads are generated and still DNA
        for _, r1_seq, _, _, r2_seq, _ in reads:
            assert set(r1_seq) <= {"A", "C", "G", "T"}

    def test_negative_strand_transcript(self):
        """Negative-strand transcripts should produce valid reads."""
        g = MutableGenome(3000, seed=42, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "-", [
            {"t_id": "t1", "exons": [(200, 600), (1000, 1400)]},
        ])
        transcripts = builder.get_transcripts()
        config = SimConfig(
            frag_mean=200, frag_std=50, frag_min=100,
            frag_max=700, read_length=100, seed=42,
        )
        sim = ReadSimulator(g, transcripts, config=config)
        reads = list(sim.simulate(50))
        assert len(reads) == 50
        # Strand char should be 'r' for negative strand
        r1_name = reads[0][0]
        parts = r1_name.rsplit("/", 1)[0].split(":")
        assert parts[2] == "r"

    def test_abundance_weighting(self):
        """Higher abundance should produce more reads from that transcript."""
        g = MutableGenome(3000, seed=42, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t_high", "exons": [(100, 500)], "abundance": 1000},
            {"t_id": "t_low", "exons": [(1500, 1900)], "abundance": 1},
        ])
        transcripts = builder.get_transcripts()
        config = SimConfig(
            frag_mean=150, frag_std=30, frag_min=50,
            frag_max=350, read_length=100, seed=42,
        )
        sim = ReadSimulator(g, transcripts, config=config)
        reads = list(sim.simulate(500))

        # Count reads per transcript from read names
        counts = {}
        for r1_name, *_ in reads:
            t_id = r1_name.split(":")[0]
            counts[t_id] = counts.get(t_id, 0) + 1

        # t_high should have vastly more reads than t_low
        assert counts.get("t_high", 0) > counts.get("t_low", 0) * 5

    def test_single_exon_no_splicing(self):
        """Single-exon transcripts should still produce valid reads."""
        _, _, sim = self._make_sim(n_exons=1)
        reads = list(sim.simulate(20))
        assert len(reads) == 20

    def test_single_exon_nrna_zeroed(self):
        """Single-exon transcripts should have nrna_abundance zeroed."""
        g = MutableGenome(3000, seed=42, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t_single", "exons": [(200, 1200)], "abundance": 100},
            {"t_id": "t_multi", "exons": [(200, 600), (1000, 1400)], "abundance": 100},
        ])
        transcripts = builder.get_transcripts()
        # Set nRNA abundance for all transcripts
        for t in transcripts:
            t.nrna_abundance = 50.0

        config = SimConfig(
            frag_mean=250, frag_std=50, frag_min=100,
            frag_max=700, read_length=100, seed=42,
        )
        sim = ReadSimulator(g, transcripts, config=config)

        # Single-exon transcript should have nrna_abundance zeroed
        t_single = [t for t in transcripts if t.t_id == "t_single"][0]
        t_multi = [t for t in transcripts if t.t_id == "t_multi"][0]
        assert t_single.nrna_abundance == 0.0
        assert t_multi.nrna_abundance == 50.0

    def test_single_exon_no_nrna_reads(self):
        """No nRNA reads should be generated for single-exon transcripts."""
        g = MutableGenome(3000, seed=42, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t_single", "exons": [(200, 1200)], "abundance": 100},
        ])
        transcripts = builder.get_transcripts()
        for t in transcripts:
            t.nrna_abundance = 100.0

        config = SimConfig(
            frag_mean=250, frag_std=50, frag_min=100,
            frag_max=700, read_length=100, seed=42,
        )
        sim = ReadSimulator(g, transcripts, config=config)
        reads = list(sim.simulate(200))
        # All reads should be mature (no nrna_ prefix)
        for r1_name, *_ in reads:
            t_id = r1_name.split(":")[0]
            assert not t_id.startswith("nrna_"), (
                f"nRNA read generated for single-exon transcript: {r1_name}"
            )


# =====================================================================
# GeneBuilder + TranscriptIndex integration
# =====================================================================


class TestGeneBuilderIndexIntegration:
    def test_build_rigel_index(self, tmp_path):
        """GeneBuilder output should be compatible with TranscriptIndex.build."""
        from rigel.index import TranscriptIndex

        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
        ])
        builder.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(1200, 1400), (1600, 1800)]},
        ])

        fasta = g.write_fasta(tmp_path)
        gtf = builder.write_gtf(tmp_path)

        index_dir = tmp_path / "index"
        TranscriptIndex.build(fasta, gtf, index_dir, write_tsv=False)

        index = TranscriptIndex.load(index_dir)
        assert index.num_transcripts == 2
        assert index.num_genes == 2
