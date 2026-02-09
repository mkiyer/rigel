"""
Shared test fixtures for hulkrna.

Provides small synthetic GTF data, minimal FASTA files, and helper
functions so that tests can run without real genomic data.
"""

import textwrap
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Minimal GTF content (GENCODE-style, 1-based inclusive coordinates)
# ---------------------------------------------------------------------------
# Two genes on chr1, one on + and one on - strand.
# Gene A (g1): two transcripts (t1 has 3 exons, t2 has 2 exons)
# Gene B (g2): one transcript (t3 has 2 exons), negative strand
#
# Coordinates (1-based inclusive as in GTF):
#   t1: exons at 100-200, 300-400, 500-600
#   t2: exons at 100-200, 500-600
#   t3: exons at 1000-1100, 1200-1300 (neg strand)

MINI_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "GeneA"; gene_type "protein_coding";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "GeneA"; gene_type "protein_coding";
    chr1\ttest\texon\t1000\t1100\t.\t-\t.\tgene_id "g2"; transcript_id "t3"; gene_name "GeneB"; gene_type "lncRNA";
    chr1\ttest\texon\t1200\t1300\t.\t-\t.\tgene_id "g2"; transcript_id "t3"; gene_name "GeneB"; gene_type "lncRNA";
""")

# Expected 0-based half-open coordinates after GTF parsing:
#   t1: exons (99,200), (299,400), (499,600)
#   t2: exons (99,200), (499,600)
#   t3: exons (999,1100), (1199,1300)


@pytest.fixture
def mini_gtf_file(tmp_path: Path) -> Path:
    """Write the minimal GTF to a temp file and return its path."""
    gtf_path = tmp_path / "test.gtf"
    gtf_path.write_text(MINI_GTF)
    return gtf_path


@pytest.fixture
def mini_fasta_file(tmp_path: Path) -> Path:
    """Create a minimal FASTA + .fai so pysam can read reference lengths.

    chr1: 2000 bp, chr2: 500 bp (chr2 has no genes — tests intergenic).
    """
    fasta_path = tmp_path / "genome.fa"
    # Write a simple FASTA with 80-char line width
    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        seq1 = "N" * 2000
        for i in range(0, len(seq1), 80):
            f.write(seq1[i : i + 80] + "\n")
        f.write(">chr2\n")
        seq2 = "N" * 500
        for i in range(0, len(seq2), 80):
            f.write(seq2[i : i + 80] + "\n")

    # Create the .fai index
    import pysam
    pysam.faidx(str(fasta_path))

    return fasta_path
