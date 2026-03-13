"""Tests for rigel.region_evidence — fragment-to-region fractional counting."""

from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

from rigel.region_evidence import (
    assemble_fragment,
    count_fragment,
    count_region_evidence,
)
from rigel.types import STRAND_NEG, STRAND_POS

# Reference names tuple (matches BAM header SQ order)
_REF_NAMES = ("chr1", "chr2")


# ---------------------------------------------------------------------------
# Mini-BAM helpers
# ---------------------------------------------------------------------------

# Header matching two-chromosome layout: chr1=2000bp, chr2=500bp
_HEADER = {
    "HD": {"VN": "1.6", "SO": "queryname"},
    "SQ": [
        {"SN": "chr1", "LN": 2000},
        {"SN": "chr2", "LN": 500},
    ],
}


def _make_read(
    qname: str,
    ref_id: int,
    pos: int,
    cigar: list[tuple[int, int]],
    *,
    is_read2: bool = False,
    is_reverse: bool = False,
    nh: int = 1,
    mate_ref_id: int | None = None,
    mate_pos: int | None = None,
    is_paired: bool = True,
) -> pysam.AlignedSegment:
    """Create a single BAM record for testing."""
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.reference_id = ref_id
    a.reference_start = pos
    a.mapping_quality = 60

    flag = 0
    if is_paired:
        flag |= 0x1 | 0x2  # paired, proper pair
    if is_reverse:
        flag |= 0x10
    if is_read2:
        flag |= 0x80
    else:
        flag |= 0x40  # read1

    a.flag = flag
    a.cigar = cigar

    seq_len = sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
    a.query_sequence = "A" * max(seq_len, 1)
    a.query_qualities = pysam.qualitystring_to_array("I" * max(seq_len, 1))

    if mate_ref_id is not None:
        a.next_reference_id = mate_ref_id
    else:
        a.next_reference_id = ref_id
    if mate_pos is not None:
        a.next_reference_start = mate_pos
    else:
        a.next_reference_start = pos + 100

    a.set_tag("NH", nh, "i")
    return a


def _write_name_sorted_bam(
    path: Path,
    reads: list[pysam.AlignedSegment],
    name: str = "test.bam",
) -> Path:
    """Write a name-sorted BAM from a list of AlignedSegments."""
    bam_path = path / name
    with pysam.AlignmentFile(str(bam_path), "wb", header=_HEADER) as out:
        for r in reads:
            out.write(r)
    return bam_path


# ---------------------------------------------------------------------------
# Mini-index region partition reminder (from conftest.py):
#
# chr1 regions (0-based half-open):
#   R0:  [0, 99)      intergenic
#   R1:  [99, 200)     exon_pos (t0+t1 share exon 99-200)
#   R2:  [200, 299)    intron_pos (tx_pos, not exon; t0 intron)
#   R3:  [299, 400)    exon_pos (t0 exon 299-400)
#   R4:  [400, 499)    intron_pos (t0 intron)
#   R5:  [499, 600)    exon_pos (t0+t1 exon 499-600)
#   R6:  [600, 999)    intergenic
#   R7:  [999, 1100)   exon_neg (t2 exon 999-1100)
#   R8:  [1100, 1199)  intron_neg (t2 intron)
#   R9:  [1199, 1300)  exon_neg (t2 exon 1199-1300)
#   R10: [1300, 2000)  intergenic
#
# chr2 regions (only present in two_chr_index, NOT in mini_index):
#   R11: [0, 500)      intergenic
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def two_chr_index(tmp_path_factory):
    """Index with chr1 (2000bp) + chr2 (500bp) for multi-chromosome tests.

    Same genes as mini_index on chr1, plus chr2 as pure intergenic (R11).
    """
    import pysam as _pysam
    from rigel.index import TranscriptIndex
    from conftest import MINI_GTF

    base = tmp_path_factory.mktemp("two_chr_idx")
    gtf_path = base / "test.gtf"
    gtf_path.write_text(MINI_GTF)

    fasta_path = base / "genome.fa"
    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        seq1 = "N" * 2000
        for i in range(0, len(seq1), 80):
            f.write(seq1[i : i + 80] + "\n")
        f.write(">chr2\n")
        seq2 = "N" * 500
        for i in range(0, len(seq2), 80):
            f.write(seq2[i : i + 80] + "\n")
    _pysam.faidx(str(fasta_path))

    idx_dir = base / "index"
    TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)
    return TranscriptIndex.load(idx_dir)


# ===========================================================================
# Section 1: Fragment assembly tests
# ===========================================================================

class TestAssembleFragment:
    """Test assemble_fragment() — from BAM records to merged blocks."""

    def test_single_read_unspliced(self):
        """Single mapped read → one aligned block."""
        r = _make_read("frag1", 0, 150, [(0, 50)])
        blocks, is_spliced = assemble_fragment([r])

        assert len(blocks) == 1
        ref_id, start, end, strand = blocks[0]
        assert ref_id == 0
        assert start == 150
        assert end == 200
        assert strand == STRAND_POS
        assert not is_spliced

    def test_single_read_reverse(self):
        """Single reverse-strand read → NEG strand block."""
        r = _make_read("frag1", 0, 150, [(0, 50)], is_reverse=True)
        blocks, is_spliced = assemble_fragment([r])

        assert len(blocks) == 1
        assert blocks[0][3] == STRAND_NEG
        assert not is_spliced

    def test_paired_end_merge(self):
        """R1+R2 overlapping on same strand → merged into one block."""
        # R1: pos=100, 50M → [100, 150), forward → POS
        r1 = _make_read("frag1", 0, 100, [(0, 50)])
        # R2: pos=130, 50M → [130, 180), reverse → POS (R2 flip: reverse→POS)
        r2 = _make_read("frag1", 0, 130, [(0, 50)],
                         is_read2=True, is_reverse=True)

        blocks, is_spliced = assemble_fragment([r1, r2])

        assert len(blocks) == 1
        ref_id, start, end, strand = blocks[0]
        assert ref_id == 0
        assert start == 100
        assert end == 180
        assert strand == STRAND_POS
        assert not is_spliced

    def test_paired_end_non_overlapping(self):
        """R1+R2 not overlapping → two separate blocks on same strand."""
        # R1: [100, 150), forward → POS
        r1 = _make_read("frag1", 0, 100, [(0, 50)])
        # R2: [300, 350), reverse → POS (R2 flip)
        r2 = _make_read("frag1", 0, 300, [(0, 50)],
                         is_read2=True, is_reverse=True)

        blocks, is_spliced = assemble_fragment([r1, r2])

        assert len(blocks) == 2
        assert blocks[0][1] == 100
        assert blocks[0][2] == 150
        assert blocks[1][1] == 300
        assert blocks[1][2] == 350
        assert not is_spliced

    def test_spliced_fragment(self):
        """CIGAR with N skip → multiple blocks, is_spliced=True."""
        # 50M 200N 50M: blocks [100,150) and [350,400)
        r = _make_read("frag1", 0, 100, [(0, 50), (3, 200), (0, 50)])
        blocks, is_spliced = assemble_fragment([r])

        assert len(blocks) == 2
        assert blocks[0] == (0, 100, 150, STRAND_POS)
        assert blocks[1] == (0, 350, 400, STRAND_POS)
        assert is_spliced

    def test_r2_strand_flip(self):
        """R2 forward → NEG strand (Illumina PE flip)."""
        # R2, forward alignment → strand flipped to NEG
        r2 = _make_read("frag1", 0, 100, [(0, 50)],
                         is_read2=True, is_reverse=False)
        blocks, _ = assemble_fragment([r2])

        assert len(blocks) == 1
        assert blocks[0][3] == STRAND_NEG

    def test_r2_reverse_flip(self):
        """R2 reverse → POS strand (Illumina PE flip)."""
        # R2, reverse alignment → strand flipped to POS
        r2 = _make_read("frag1", 0, 100, [(0, 50)],
                         is_read2=True, is_reverse=True)
        blocks, _ = assemble_fragment([r2])

        assert len(blocks) == 1
        assert blocks[0][3] == STRAND_POS

    def test_unmapped_read_skipped(self):
        """Unmapped record is ignored → empty blocks."""
        r = _make_read("frag1", 0, 0, [(0, 50)])
        r.flag |= 0x4  # unmapped
        blocks, is_spliced = assemble_fragment([r])

        assert len(blocks) == 0
        assert not is_spliced

    def test_cigar_with_deletions(self):
        """CIGAR D op advances reference but is part of aligned block."""
        # 25M 10D 25M → single contiguous block [100, 160)
        r = _make_read("frag1", 0, 100, [(0, 25), (2, 10), (0, 25)])
        blocks, is_spliced = assemble_fragment([r])

        assert len(blocks) == 1
        assert blocks[0][1] == 100
        assert blocks[0][2] == 160
        assert not is_spliced

    def test_cigar_with_insertions(self):
        """CIGAR I op does not advance reference position."""
        # 25M 5I 25M → block [100, 150) (insertion is query-only)
        r = _make_read("frag1", 0, 100, [(0, 25), (1, 5), (0, 25)])
        blocks, is_spliced = assemble_fragment([r])

        assert len(blocks) == 1
        assert blocks[0][1] == 100
        assert blocks[0][2] == 150
        assert not is_spliced

    def test_cigar_soft_clip(self):
        """Soft clip does not count as aligned bases."""
        # 5S 50M 5S → block [100, 150), seq_len=60
        r = _make_read("frag1", 0, 100, [(4, 5), (0, 50), (4, 5)])
        blocks, _ = assemble_fragment([r])

        assert len(blocks) == 1
        assert blocks[0][1] == 100
        assert blocks[0][2] == 150

    def test_empty_after_filter(self):
        """All records filtered (unmapped) → empty blocks."""
        r1 = _make_read("frag1", 0, 0, [(0, 50)])
        r1.flag |= 0x4
        r2 = _make_read("frag1", 0, 0, [(0, 50)], is_read2=True)
        r2.flag |= 0x4
        blocks, is_spliced = assemble_fragment([r1, r2])
        assert len(blocks) == 0


# ===========================================================================
# Section 2: Fractional counting tests
# ===========================================================================

class TestCountFragment:
    """Test count_fragment() — fractional region overlap accumulation."""

    def test_single_block_single_region(self, mini_index):
        """Fragment fully within one region → frac = 1.0."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        # Block fully inside R6 intergenic [600, 999): [700, 800)
        blocks = [(0, 700, 800, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts)

        # R6 should have unspliced_pos = 1.0
        r6 = mini_index.region_df.loc[
            mini_index.region_df["region_id"] == 6
        ].index[0]
        assert counts[6, 0] == pytest.approx(1.0)
        # All others zero
        assert counts.sum() == pytest.approx(1.0)

    def test_single_block_two_regions(self, mini_index):
        """Fragment spanning a region boundary → correct fractions."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        # Block [180, 220): spans R1 [99,200) and R2 [200,299)
        # Overlap with R1: 200-180 = 20bp; overlap with R2: 220-200 = 20bp
        # Total aligned = 40bp; frac_R1 = 20/40 = 0.5, frac_R2 = 20/40 = 0.5
        blocks = [(0, 180, 220, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts)

        assert counts[1, 0] == pytest.approx(0.5)  # R1 unspliced_pos
        assert counts[2, 0] == pytest.approx(0.5)  # R2 unspliced_pos
        assert counts.sum() == pytest.approx(1.0)

    def test_spliced_two_blocks_three_regions(self, mini_index):
        """Spliced fragment with blocks in different regions."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        # Block A: [150, 200) → 50bp in R1 [99,200) exon_pos
        # Block B: [350, 400) → 50bp in R3 [299,400) exon_pos
        # Total aligned = 100bp
        # frac_R1 = 50/100 = 0.5, frac_R3 = 50/100 = 0.5
        blocks = [(0, 150, 200, STRAND_POS),
                  (0, 350, 400, STRAND_POS)]
        count_fragment(blocks, True, mini_index.region_cr, _REF_NAMES, counts)

        # Spliced POS → column 2
        assert counts[1, 2] == pytest.approx(0.5)  # R1 spliced_pos
        assert counts[3, 2] == pytest.approx(0.5)  # R3 spliced_pos
        assert counts.sum() == pytest.approx(1.0)

    def test_fractions_sum_to_one(self, mini_index):
        """Any fragment's fractional counts sum to 1.0."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        # Asymmetric boundary span: [80, 220) → 140bp total
        # R0 [0,99): 99-80 = 19bp; R1 [99,200): 200-99 = 101bp; R2 [200,299): 220-200 = 20bp
        # fracs: 19/140 + 101/140 + 20/140 = 1.0
        blocks = [(0, 80, 220, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts)

        assert counts.sum() == pytest.approx(1.0)
        assert counts[0, 0] == pytest.approx(19.0 / 140.0)
        assert counts[1, 0] == pytest.approx(101.0 / 140.0)
        assert counts[2, 0] == pytest.approx(20.0 / 140.0)

    def test_intron_not_counted(self, mini_index):
        """Intronic gap between spliced blocks contributes zero overlap."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        # Block A: [150, 200) in R1 → 50bp
        # (intron 200-350 is NOT a block)
        # Block B: [350, 400) in R3 → 50bp
        # Total aligned = 100bp (intron excluded from denominator)
        blocks = [(0, 150, 200, STRAND_POS),
                  (0, 350, 400, STRAND_POS)]
        count_fragment(blocks, True, mini_index.region_cr, _REF_NAMES, counts)

        assert counts.sum() == pytest.approx(1.0)
        # R2 (intron region 200-299) gets zero
        assert counts[2, :].sum() == 0.0

    def test_pos_neg_separation(self, mini_index):
        """POS and NEG strand fragments go to separate columns."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        # POS fragment in R6 [600,999)
        blocks_pos = [(0, 700, 800, STRAND_POS)]
        count_fragment(blocks_pos, False, mini_index.region_cr, _REF_NAMES, counts)

        # NEG fragment in R6 [600,999)
        blocks_neg = [(0, 700, 800, STRAND_NEG)]
        count_fragment(blocks_neg, False, mini_index.region_cr, _REF_NAMES, counts)

        assert counts[6, 0] == pytest.approx(1.0)  # unspliced_pos
        assert counts[6, 1] == pytest.approx(1.0)  # unspliced_neg

    def test_spliced_unspliced_separation(self, mini_index):
        """Spliced and unspliced fragments go to separate columns."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        blocks = [(0, 700, 800, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts)  # unspliced
        count_fragment(blocks, True, mini_index.region_cr, _REF_NAMES, counts)   # spliced

        assert counts[6, 0] == pytest.approx(1.0)  # unspliced_pos
        assert counts[6, 2] == pytest.approx(1.0)  # spliced_pos

    def test_neg_strand_spliced(self, mini_index):
        """NEG strand spliced fragment → spliced_neg column."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        blocks = [(0, 700, 800, STRAND_NEG)]
        count_fragment(blocks, True, mini_index.region_cr, _REF_NAMES, counts)

        assert counts[6, 3] == pytest.approx(1.0)  # spliced_neg

    def test_chimeric_multi_ref_excluded(self, mini_index):
        """Multi-reference blocks (chimeric) → not counted."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        blocks = [(0, 700, 800, STRAND_POS),
                  (1, 100, 200, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts)

        assert counts.sum() == 0.0

    def test_chimeric_multi_strand_excluded(self, mini_index):
        """Multi-strand blocks (chimeric) → not counted."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        blocks = [(0, 700, 750, STRAND_POS),
                  (0, 750, 800, STRAND_NEG)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts)

        assert counts.sum() == 0.0

    def test_empty_blocks(self, mini_index):
        """Empty block list → nothing counted."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        count_fragment([], False, mini_index.region_cr, _REF_NAMES, counts)
        assert counts.sum() == 0.0

    def test_intergenic_fragment_counted(self, mini_index):
        """Fragment in intergenic region is still counted."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        # R0 is intergenic [0, 99)
        blocks = [(0, 20, 80, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts)

        assert counts[0, 0] == pytest.approx(1.0)

    def test_chr2_intergenic(self, two_chr_index):
        """Fragment on chr2 (entirely intergenic) → correct region."""
        n = len(two_chr_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        # R11 on chr2 [0, 500)
        blocks = [(1, 100, 200, STRAND_NEG)]
        count_fragment(blocks, False, two_chr_index.region_cr, _REF_NAMES, counts)

        assert counts[11, 1] == pytest.approx(1.0)  # unspliced_neg
        assert counts.sum() == pytest.approx(1.0)


# ===========================================================================
# Section 3: Fragment-length tabulation tests
# ===========================================================================

class TestFragLengthTabulation:
    """Test FL observation collection in count_fragment()."""

    def test_fl_single_region_fragment(self, mini_index):
        """Contained unspliced fragment → FL observation emitted."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)
        fl_rids: list[int] = []
        fl_lens: list[int] = []

        # Fully inside R6 [600, 999): block [700, 900)
        blocks = [(0, 700, 900, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts,
                       fl_rids, fl_lens)

        assert len(fl_rids) == 1
        assert fl_rids[0] == 6
        assert fl_lens[0] == 200  # 900 - 700

    def test_fl_boundary_spanning_excluded(self, mini_index):
        """Fragment spanning region boundary → no FL observation."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)
        fl_rids: list[int] = []
        fl_lens: list[int] = []

        # [180, 220) spans R1/R2 boundary at 200
        blocks = [(0, 180, 220, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts,
                       fl_rids, fl_lens)

        assert len(fl_rids) == 0

    def test_fl_spliced_excluded(self, mini_index):
        """Spliced fragment → no FL observation (even if single region)."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)
        fl_rids: list[int] = []
        fl_lens: list[int] = []

        blocks = [(0, 700, 900, STRAND_POS)]
        count_fragment(blocks, True, mini_index.region_cr, _REF_NAMES, counts,
                       fl_rids, fl_lens)

        assert len(fl_rids) == 0

    def test_fl_correct_frag_len_pe(self, mini_index):
        """PE fragment with two non-overlapping blocks → genomic_footprint."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)
        fl_rids: list[int] = []
        fl_lens: list[int] = []

        # Two blocks in R6 [600, 999): [700,750) and [800,850)
        # genomic_footprint = 850 - 700 = 150
        blocks = [(0, 700, 750, STRAND_POS),
                  (0, 800, 850, STRAND_POS)]
        # Unspliced (no intron between R1/R2 gap)
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts,
                       fl_rids, fl_lens)

        assert len(fl_rids) == 1
        assert fl_lens[0] == 150  # genomic footprint

    def test_fl_none_lists_no_error(self, mini_index):
        """When fl lists are None, no error occurs."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)

        blocks = [(0, 700, 900, STRAND_POS)]
        # Should not raise
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts,
                       None, None)

    def test_fl_intergenic_fragments(self, mini_index):
        """Intergenic contained fragment produces FL observation."""
        n = len(mini_index.region_df)
        counts = np.zeros((n, 4), dtype=np.float64)
        fl_rids: list[int] = []
        fl_lens: list[int] = []

        # R0 intergenic [0, 99): fragment [10, 80)
        blocks = [(0, 10, 80, STRAND_POS)]
        count_fragment(blocks, False, mini_index.region_cr, _REF_NAMES, counts,
                       fl_rids, fl_lens)

        assert len(fl_rids) == 1
        assert fl_rids[0] == 0
        assert fl_lens[0] == 70


# ===========================================================================
# Section 4: Integration tests — count_region_evidence with BAM files
# ===========================================================================

class TestCountRegionEvidence:
    """Integration tests: full BAM → region evidence pipeline."""

    def test_empty_bam(self, mini_index, tmp_path):
        """Empty BAM → all-zero counts, empty FL table."""
        bam_path = _write_name_sorted_bam(tmp_path, [], "empty.bam")
        region_counts, fl_table = count_region_evidence(bam_path, mini_index)

        assert len(region_counts) == len(mini_index.region_df)
        assert region_counts["n_unspliced_pos"].sum() == 0.0
        assert region_counts["n_unspliced_neg"].sum() == 0.0
        assert region_counts["n_spliced_pos"].sum() == 0.0
        assert region_counts["n_spliced_neg"].sum() == 0.0
        assert len(fl_table) == 0

    def test_all_regions_present(self, mini_index, tmp_path):
        """Output has one row per region even with a single fragment."""
        r = _make_read("frag1", 0, 700, [(0, 50)])
        bam_path = _write_name_sorted_bam(tmp_path, [r])
        region_counts, _ = count_region_evidence(bam_path, mini_index)

        assert len(region_counts) == len(mini_index.region_df)
        assert list(region_counts["region_id"]) == list(
            mini_index.region_df["region_id"]
        )

    def test_single_unspliced_pos_fragment(self, mini_index, tmp_path):
        """One unspliced POS fragment fully in R6 intergenic."""
        # R1 forward: [700, 750) → POS strand
        r = _make_read("frag1", 0, 700, [(0, 50)])
        bam_path = _write_name_sorted_bam(tmp_path, [r])
        region_counts, fl_table = count_region_evidence(bam_path, mini_index)

        row6 = region_counts.loc[region_counts["region_id"] == 6].iloc[0]
        assert row6["n_unspliced_pos"] == pytest.approx(1.0)
        assert row6["n_unspliced_neg"] == 0.0
        assert row6["n_spliced_pos"] == 0.0
        assert row6["n_spliced_neg"] == 0.0

        # FL observation
        assert len(fl_table) == 1
        assert fl_table.iloc[0]["region_id"] == 6
        assert fl_table.iloc[0]["frag_len"] == 50

    def test_paired_end_fragment(self, mini_index, tmp_path):
        """PE pair → merged blocks, correct counting."""
        # R1: [700, 750) forward → POS
        r1 = _make_read("frag1", 0, 700, [(0, 50)])
        # R2: [720, 770) reverse → POS (R2 flip: reverse→POS)
        r2 = _make_read("frag1", 0, 720, [(0, 50)],
                         is_read2=True, is_reverse=True)
        bam_path = _write_name_sorted_bam(tmp_path, [r1, r2])
        region_counts, fl_table = count_region_evidence(bam_path, mini_index)

        row6 = region_counts.loc[region_counts["region_id"] == 6].iloc[0]
        assert row6["n_unspliced_pos"] == pytest.approx(1.0)

        # Merged block: [700, 770) → genomic footprint = 70
        assert len(fl_table) == 1
        assert fl_table.iloc[0]["frag_len"] == 70

    def test_spliced_fragment(self, mini_index, tmp_path):
        """Spliced fragment → spliced count column."""
        # 50M 100N 50M at pos=150 → blocks [150,200) and [300,350)
        # [150,200) → R1 [99,200); [300,350) → R3 [299,400)
        r = _make_read("frag1", 0, 150, [(0, 50), (3, 100), (0, 50)])
        bam_path = _write_name_sorted_bam(tmp_path, [r])
        region_counts, fl_table = count_region_evidence(bam_path, mini_index)

        row1 = region_counts.loc[region_counts["region_id"] == 1].iloc[0]
        row3 = region_counts.loc[region_counts["region_id"] == 3].iloc[0]
        assert row1["n_spliced_pos"] == pytest.approx(0.5)
        assert row3["n_spliced_pos"] == pytest.approx(0.5)

        # No FL observation for spliced fragment
        assert len(fl_table) == 0

    def test_multimapper_excluded(self, mini_index, tmp_path):
        """NH > 1 → fragment not counted."""
        r = _make_read("frag1", 0, 700, [(0, 50)], nh=2)
        bam_path = _write_name_sorted_bam(tmp_path, [r])
        region_counts, fl_table = count_region_evidence(bam_path, mini_index)

        assert region_counts["n_unspliced_pos"].sum() == 0.0
        assert len(fl_table) == 0

    def test_neg_strand_fragment(self, mini_index, tmp_path):
        """Reverse-strand R1 → NEG strand, unspliced_neg column."""
        r = _make_read("frag1", 0, 700, [(0, 50)], is_reverse=True)
        bam_path = _write_name_sorted_bam(tmp_path, [r])
        region_counts, _ = count_region_evidence(bam_path, mini_index)

        row6 = region_counts.loc[region_counts["region_id"] == 6].iloc[0]
        assert row6["n_unspliced_neg"] == pytest.approx(1.0)
        assert row6["n_unspliced_pos"] == 0.0

    def test_boundary_spanning_fragment(self, mini_index, tmp_path):
        """Fragment spanning R1/R2 boundary → fractional counts."""
        # [180, 220): 20bp in R1, 20bp in R2
        r = _make_read("frag1", 0, 180, [(0, 40)])
        bam_path = _write_name_sorted_bam(tmp_path, [r])
        region_counts, fl_table = count_region_evidence(bam_path, mini_index)

        row1 = region_counts.loc[region_counts["region_id"] == 1].iloc[0]
        row2 = region_counts.loc[region_counts["region_id"] == 2].iloc[0]
        assert row1["n_unspliced_pos"] == pytest.approx(0.5)
        assert row2["n_unspliced_pos"] == pytest.approx(0.5)

        # Spans two regions → no FL observation
        assert len(fl_table) == 0

    def test_multi_chromosome(self, two_chr_index, tmp_path):
        """Fragments on different chromosomes → correct regions."""
        # chr1 R6 [600,999): [700,750) POS
        r1 = _make_read("frag1", 0, 700, [(0, 50)])
        # chr2 R11 [0,500): [100,150) NEG
        r2 = _make_read("frag2", 1, 100, [(0, 50)], is_reverse=True)
        bam_path = _write_name_sorted_bam(tmp_path, [r1, r2])
        region_counts, fl_table = count_region_evidence(bam_path, two_chr_index)

        row6 = region_counts.loc[region_counts["region_id"] == 6].iloc[0]
        row11 = region_counts.loc[region_counts["region_id"] == 11].iloc[0]
        assert row6["n_unspliced_pos"] == pytest.approx(1.0)
        assert row11["n_unspliced_neg"] == pytest.approx(1.0)
        assert len(fl_table) == 2

    def test_multiple_fragments_accumulate(self, mini_index, tmp_path):
        """Multiple fragments accumulate in region counts."""
        reads = []
        for i in range(5):
            r = _make_read(f"frag{i}", 0, 700, [(0, 50)])
            reads.append(r)
        bam_path = _write_name_sorted_bam(tmp_path, reads)
        region_counts, fl_table = count_region_evidence(bam_path, mini_index)

        row6 = region_counts.loc[region_counts["region_id"] == 6].iloc[0]
        assert row6["n_unspliced_pos"] == pytest.approx(5.0)
        assert len(fl_table) == 5

    def test_total_fractions_sum_to_n_fragments(self, mini_index, tmp_path):
        """Total fractional counts across all regions equals n_fragments."""
        reads = []
        # 3 fragments: one fully in R6, one spanning R1/R2, one in R10
        reads.append(_make_read("f1", 0, 700, [(0, 50)]))
        reads.append(_make_read("f2", 0, 180, [(0, 40)]))
        reads.append(_make_read("f3", 0, 1400, [(0, 50)]))
        bam_path = _write_name_sorted_bam(tmp_path, reads)
        region_counts, _ = count_region_evidence(bam_path, mini_index)

        total = (
            region_counts["n_unspliced_pos"].sum()
            + region_counts["n_unspliced_neg"].sum()
            + region_counts["n_spliced_pos"].sum()
            + region_counts["n_spliced_neg"].sum()
        )
        assert total == pytest.approx(3.0)

    def test_output_dtypes(self, mini_index, tmp_path):
        """Output DataFrames have correct dtypes."""
        r = _make_read("frag1", 0, 700, [(0, 50)])
        bam_path = _write_name_sorted_bam(tmp_path, [r])
        region_counts, fl_table = count_region_evidence(bam_path, mini_index)

        assert region_counts["region_id"].dtype == np.int32  # explicitly cast
        assert region_counts["n_unspliced_pos"].dtype == np.float32
        assert region_counts["n_unspliced_neg"].dtype == np.float32
        assert region_counts["n_spliced_pos"].dtype == np.float32
        assert region_counts["n_spliced_neg"].dtype == np.float32

        assert fl_table["region_id"].dtype == np.int32
        assert fl_table["frag_len"].dtype == np.int32
