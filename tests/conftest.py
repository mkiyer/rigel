"""The suite-wide pytest fixtures, and the one custom option.

Three fixtures built on the minimal GTF in `_index_builder`: the GTF written to a temp file, a minimal
all-N FASTA with its `.fai` so pysam can read reference lengths, and a session-scoped
`TranscriptIndex` loaded from both. Together they let a structural test run with no real genomic data.
`--update-golden` is registered here because a pytest option has to be declared in a conftest, and
`tests/test_golden_output.py` reads it.
"""

from pathlib import Path

import pytest

from _index_builder import MINI_GTF, build_test_index

# ---------------------------------------------------------------------------
# Global pytest options
# ---------------------------------------------------------------------------


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Regenerate golden output files instead of comparing.",
    )


# ---------------------------------------------------------------------------
# Fixtures over the minimal GTF (GENCODE-style, 1-based inclusive coordinates)
# ---------------------------------------------------------------------------
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


@pytest.fixture(scope="session")
def mini_index(tmp_path_factory):
    """Standard MINI_GTF TranscriptIndex: t0(3-exon), t1(2-exon), t2(neg strand).

    GTF layout (0-based half-open after parse):
      g1 (+): t0 exons (99,200),(299,400),(499,600)
              t1 exons (99,200),(499,600)
      g2 (-): t2 exons (999,1100),(1199,1300)

    Transcript indices: t0=0, t1=1, t2=2
    Gene indices:       g1=0, g2=1
    """
    return build_test_index(tmp_path_factory, MINI_GTF, name="mini_idx")
