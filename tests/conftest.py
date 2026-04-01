"""
Shared test fixtures and helpers for rigel.

Provides small synthetic GTF data, minimal FASTA files, and helper
functions so that tests can run without real genomic data.
"""

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rigel.scored_fragments import Locus, ScoredFragments
from rigel.index import TranscriptIndex
from rigel.splice import SpliceStrandCol

# ---------------------------------------------------------------------------
# Global pytest options
# ---------------------------------------------------------------------------


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden", action="store_true", default=False,
        help="Regenerate golden output files instead of comparing.",
    )

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


# ---------------------------------------------------------------------------
# Shared index builder (used by test_buffer, test_resolution, etc.)
# ---------------------------------------------------------------------------


def build_test_index(tmp_path_factory, gtf_text, genome_size=2000, name="idx"):
    """Build a TranscriptIndex from a GTF string (session/module-scoped helper).

    Parameters
    ----------
    tmp_path_factory : pytest.TempPathFactory
        Pytest factory for session/module-scoped temp dirs.
    gtf_text : str
        GTF content string (GENCODE-style, 1-based inclusive).
    genome_size : int
        Length of chr1 in the synthetic FASTA (default: 2000).
    name : str
        Sub-directory name for this index (keeps separate indexes apart).

    Returns
    -------
    TranscriptIndex
        Loaded index with C++ FragmentResolver ready.
    """
    import pysam

    base = tmp_path_factory.mktemp(name)
    gtf_path = base / "test.gtf"
    gtf_path.write_text(gtf_text)

    fasta_path = base / "genome.fa"
    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        seq = "N" * genome_size
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")
    pysam.faidx(str(fasta_path))

    idx_dir = base / "index"
    TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)
    return TranscriptIndex.load(idx_dir)


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


# ---------------------------------------------------------------------------
# Shared EM helpers (used by test_estimator and test_gdna)
# ---------------------------------------------------------------------------

# Default column for UNSPLICED_SENSE
_UNSPLICED_SENSE = int(SpliceStrandCol.UNSPLICED_SENSE)


def _make_locus_em_data(
    t_indices_per_unit,
    log_liks_per_unit=None,
    count_cols_per_unit=None,
    num_transcripts=None,
    rc=None,
    include_nrna=False,
    include_gdna=False,
    nrna_log_lik=-2.0,
    gdna_log_lik=0.0,
    locus_gamma=0.1,
):
    """Build (ScoredFragments, [Locus], locus_gammas, index) for batch EM tests.

    Returns a tuple (em_data, loci, locus_gammas, index) suitable for
    ``run_batch_locus_em()``.

    Parameters
    ----------
    t_indices_per_unit : list[list[int]]
        GLOBAL mRNA transcript indices per unit.
    rc : AbundanceEstimator or None
        If provided, used for unambig_counts.  The helper also sets
        ``_transcript_spans`` and ``_exonic_lengths`` if they are None.
    include_nrna : bool
        If True, units are marked as unspliced (``is_spliced=False``)
        so the batch C++ adds nRNA shadow candidates.
    include_gdna : bool
        If True, ``locus_gamma > 0`` and ``gdna_log_liks`` are finite
        so the batch C++ adds a gDNA component.
    """
    if num_transcripts is None:
        all_t = [t for unit in t_indices_per_unit for t in unit]
        num_transcripts = (max(all_t) + 1) if all_t else 1

    n_t = num_transcripts
    n_units = len(t_indices_per_unit)

    # Build mRNA-only CSR (no nRNA/gDNA — batch C++ adds those)
    offsets = [0]
    flat_t = []
    flat_lk = []
    flat_cc = []

    for u, t_list in enumerate(t_indices_per_unit):
        for j, t_idx in enumerate(t_list):
            flat_t.append(t_idx)
            flat_lk.append(
                log_liks_per_unit[u][j] if log_liks_per_unit else 0.0
            )
            flat_cc.append(
                count_cols_per_unit[u][j]
                if count_cols_per_unit
                else _UNSPLICED_SENSE
            )
        offsets.append(len(flat_t))

    n_candidates = len(flat_t)

    # Per-unit locus tracking
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_cc = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        if t_list:
            locus_t[u] = t_list[0]
            locus_cc[u] = (
                count_cols_per_unit[u][0]
                if count_cols_per_unit
                else _UNSPLICED_SENSE
            )

    # is_spliced: True (spliced) → no nRNA/gDNA shadows.
    # For include_nrna or include_gdna, set False (unspliced).
    if include_nrna or include_gdna:
        is_spliced = np.zeros(n_units, dtype=bool)
    else:
        is_spliced = np.ones(n_units, dtype=bool)

    # gDNA log-likelihoods per unit
    if include_gdna:
        gdna_log_liks = np.full(n_units, gdna_log_lik, dtype=np.float64)
    else:
        gdna_log_liks = np.full(n_units, -np.inf, dtype=np.float64)

    em_data = ScoredFragments(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.array(flat_lk, dtype=np.float64),
        count_cols=np.array(flat_cc, dtype=np.uint8),
        coverage_weights=np.ones(n_candidates, dtype=np.float64),
        tx_starts=np.zeros(n_candidates, dtype=np.int32),
        tx_ends=np.ones(n_candidates, dtype=np.int32),
        locus_t_indices=locus_t,
        locus_count_cols=locus_cc,
        is_spliced=is_spliced,
        gdna_log_liks=gdna_log_liks,
        frag_ids=np.arange(n_units, dtype=np.int64),
        frag_class=np.zeros(n_units, dtype=np.int8),
        splice_type=np.zeros(n_units, dtype=np.uint8),
        n_units=n_units,
        n_candidates=n_candidates,
        genomic_footprints=np.full(n_units, 200, dtype=np.int32),
    )

    locus = Locus(
        locus_id=0,
        transcript_indices=np.arange(n_t, dtype=np.int32),
        unit_indices=np.arange(n_units, dtype=np.int32),
        gdna_span=10000,
        merged_intervals=[("chr1", 0, 10000)],
    )
    loci = [locus]

    locus_gammas = np.array([locus_gamma if include_gdna else 0.0],
                          dtype=np.float64)

    index = _MockBatchIndex(n_t)

    # Ensure estimator geometry arrays are set if rc is provided
    if rc is not None:
        _ensure_estimator_geometry(rc)

    return em_data, loci, locus_gammas, index


class _MockBatchIndex:
    """Minimal index mock for batch locus EM tests."""

    def __init__(self, num_transcripts):
        self.num_transcripts = num_transcripts
        self.t_df = pd.DataFrame({
            "t_id": [f"t{i}" for i in range(num_transcripts)],
            "ref": ["chr1"] * num_transcripts,
            "start": np.zeros(num_transcripts, dtype=np.int64),
            "end": np.full(num_transcripts, 10000, dtype=np.int64),
            "length": np.full(num_transcripts, 1000, dtype=np.int64),
            "is_nrna": np.zeros(num_transcripts, dtype=bool),
            "is_synthetic": np.zeros(num_transcripts, dtype=bool),
        })


def _ensure_estimator_geometry(rc):
    """Set required geometry arrays on estimator if not already set."""
    n_t = rc.num_transcripts
    if rc._transcript_spans is None:
        rc._transcript_spans = np.full(n_t, 10000.0, dtype=np.float64)
    if rc._exonic_lengths is None:
        rc._exonic_lengths = np.full(n_t, 1000.0, dtype=np.float64)


def _run_and_assign(rc, em_data, loci=None, index=None, locus_gammas=None,
                    *, em_iterations=10):
    """Run batch locus EM via the partitioned path. Returns pool_counts dict.

    Accepts either the new tuple form (em_data, loci, locus_gammas, index)
    separately, or the tuple returned by ``_make_locus_em_data`` as ``em_data``.
    """
    from rigel.partition import partition_and_free

    # Unpack tuple form from _make_locus_em_data
    if isinstance(em_data, tuple):
        em_data, loci, locus_gammas, index = em_data

    _ensure_estimator_geometry(rc)

    # Partition ScoredFragments into per-locus LocusPartition objects
    partitions = partition_and_free(em_data, loci)

    # Build the 12-tuples and transcript index lists expected by C++
    partition_tuples = [
        (
            p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
            p.tx_starts, p.tx_ends, p.count_cols,
            p.is_spliced, p.gdna_log_liks, p.genomic_footprints,
            p.locus_t_indices, p.locus_count_cols,
        )
        for p in [partitions[i] for i in range(len(loci))]
    ]
    locus_t_lists = [l.transcript_indices for l in loci]
    gdna_spans = np.array([l.gdna_span for l in loci], dtype=np.int64)

    total_gdna, _locus_mrna, _locus_gdna = rc.run_batch_locus_em_partitioned(
        partition_tuples, locus_t_lists,
        locus_gammas, gdna_spans, index,
        em_iterations=em_iterations,
    )
    rc._gdna_em_total += total_gdna

    return {
        "mrna": float(rc.em_counts.sum()),
        "nrna": float(rc.nrna_em_count),
        "gdna": float(total_gdna),
    }
