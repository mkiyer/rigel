"""Tests for rigel.sim — simulation framework components."""

import pytest

from rigel.sim.annotation import GeneBuilder
from rigel.sim.genome import MutableGenome, reverse_complement
from rigel.sim.synthetic_genome import generate_genes
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


def test_generate_genes_respects_isoform_bounds_and_antisense_fraction():
    genes = generate_genes(
        500_000,
        10,
        13,
        min_isoforms=1,
        max_isoforms=5,
        target_transcripts=None,
        antisense_overlap_frac=0.0,
    )

    assert len(genes) == 10
    assert all(1 <= len(gene.transcripts) <= 5 for gene in genes)


class TestGeneBuilder:
    def test_single_gene_single_transcript(self, tmp_path):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 300), (500, 700)], "abundance": 100},
            ],
        )
        transcripts = builder.get_transcripts()
        assert len(transcripts) == 1
        assert transcripts[0].t_id == "t1"
        assert transcripts[0].g_id == "g1"
        assert transcripts[0].strand == Strand.POS
        assert len(transcripts[0].exons) == 2

    def test_splice_motif_positive_strand(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
            ],
        )
        # Intron is (300, 500): donor at 300, acceptor at 498
        assert g[300:302] == "GT"
        assert g[498:500] == "AG"

    def test_splice_motif_negative_strand(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "-",
            [
                {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
            ],
        )
        # Negative strand: CT at donor, AC at acceptor
        assert g[300:302] == "CT"
        assert g[498:500] == "AC"

    def test_multi_exon_splice_motifs(self):
        g = MutableGenome(3000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 200), (400, 500), (700, 800)]},
            ],
        )
        # Intron 1: (200, 400)
        assert g[200:202] == "GT"
        assert g[398:400] == "AG"
        # Intron 2: (500, 700)
        assert g[500:502] == "GT"
        assert g[698:700] == "AG"

    def test_multi_isoform_gene(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 200), (400, 500), (700, 800)]},
                {"t_id": "t2", "exons": [(100, 200), (700, 800)]},
            ],
        )
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
            builder.add_gene(
                "g1",
                "+",
                [
                    {"t_id": "t1", "exons": [(100, 300), (250, 500)]},
                ],
            )

    def test_exon_out_of_bounds_raises(self):
        g = MutableGenome(500, seed=1, name="chr1")
        builder = GeneBuilder(g)
        with pytest.raises(ValueError, match="outside genome"):
            builder.add_gene(
                "g1",
                "+",
                [
                    {"t_id": "t1", "exons": [(100, 600)]},
                ],
            )

    def test_t_index_assignment(self):
        g = MutableGenome(3000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 200)]},
            ],
        )
        builder.add_gene(
            "g2",
            "-",
            [
                {"t_id": "t2", "exons": [(1000, 1200)]},
            ],
        )
        transcripts = builder.get_transcripts()
        assert transcripts[0].t_index == 0
        assert transcripts[1].t_index == 1
        assert transcripts[0].g_index != transcripts[1].g_index

    def test_write_gtf(self, tmp_path):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
            ],
        )
        gtf_path = builder.write_gtf(tmp_path)
        assert gtf_path.exists()

        # Read it back with the GTF parser
        from rigel.gtf import GTFRecord

        features = list(GTFRecord.parse_file(gtf_path))
        assert len(features) == 2  # 2 exon boundaries
        assert all(f.feature == "exon" for f in features)
        assert all(f.attrs["gene_id"] == "g1" for f in features)
        assert all(f.attrs["transcript_id"] == "t1" for f in features)

    def test_write_gtf_roundtrip(self, tmp_path):
        """GTF → Transcript.read_gtf should reconstruct the transcripts."""
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 300), (500, 700)], "abundance": 50.0},
                {"t_id": "t2", "exons": [(100, 300)], "abundance": 25.0},
            ],
        )
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
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 300)], "abundance": 42.5},
            ],
        )
        gtf_path = builder.write_gtf(tmp_path)

        from rigel.gtf import GTFRecord

        features = list(GTFRecord.parse_file(gtf_path))
        assert features[0].score == 42.5

    def test_intron_too_short_raises(self):
        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        with pytest.raises(ValueError, match="too short"):
            builder.add_gene(
                "g1",
                "+",
                [
                    {"t_id": "t1", "exons": [(100, 200), (202, 300)]},
                ],
            )

    def test_gtf_to_bed12(self, tmp_path):
        """GTF → BED12 conversion produces valid 12-column BED."""
        from rigel.index import gtf_to_bed12

        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
            ],
        )
        builder.add_gene(
            "g2",
            "-",
            [
                {"t_id": "t2", "exons": [(1000, 1200), (1400, 1500), (1700, 1800)]},
            ],
        )
        gtf_path = builder.write_gtf(tmp_path)
        bed_path = tmp_path / "annotation.bed"
        result_path = gtf_to_bed12(gtf_path, bed_path)
        assert result_path.exists()

        boundaries = bed_path.read_text().strip().split("\n")
        assert len(boundaries) == 2

        # Parse first transcript (2-exon, + strand)
        fields = boundaries[0].split("\t")
        assert len(fields) == 12
        assert fields[0] == "chr1"  # ref
        assert fields[1] == "100"  # refStart
        assert fields[2] == "700"  # refEnd
        assert fields[3] == "t1"  # name
        assert fields[5] == "+"  # strand
        assert fields[9] == "2"  # blockCount
        assert fields[10] == "200,200"  # blockSizes (300-100, 700-500)
        assert fields[11] == "0,400"  # blockStarts (100-100, 500-100)

        # Parse second transcript (3-exon, - strand)
        fields = boundaries[1].split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "1000"
        assert fields[2] == "1800"
        assert fields[3] == "t2"
        assert fields[5] == "-"
        assert fields[9] == "3"
        assert fields[10] == "200,100,100"  # sizes
        assert fields[11] == "0,400,700"  # starts relative to 1000


# =====================================================================
# GeneBuilder + TranscriptIndex integration
# =====================================================================


class TestGeneBuilderIndexIntegration:
    def test_build_rigel_index(self, tmp_path):
        """GeneBuilder output should be compatible with TranscriptIndex.build."""
        from rigel.index import TranscriptIndex

        g = MutableGenome(2000, seed=1, name="chr1")
        builder = GeneBuilder(g)
        builder.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(100, 300), (500, 700)]},
            ],
        )
        builder.add_gene(
            "g2",
            "-",
            [
                {"t_id": "t2", "exons": [(1200, 1400), (1600, 1800)]},
            ],
        )

        fasta = g.write_fasta(tmp_path)
        gtf = builder.write_gtf(tmp_path)

        index_dir = tmp_path / "index"
        TranscriptIndex.build(fasta, gtf, index_dir, write_tsv=False)

        index = TranscriptIndex.load(index_dir)
        assert index.num_transcripts == 4  # 2 annotated + 2 synthetic nRNA
        assert index.num_annotated_genes == 2
        assert index.num_genes == 4  # 2 annotated + 2 synthetic gene rows
