"""Tests for hulkrna.fragment — Fragment.from_reads with mock pysam reads."""

from unittest.mock import MagicMock

import pysam
import pytest

from hulkrna.types import GenomicInterval, Strand
from hulkrna.fragment import Fragment


# =====================================================================
# Mock read helpers
# =====================================================================


def _mock_read(
    ref_name="chr1",
    ref_start=100,
    is_reverse=False,
    cigartuples=None,
    xs_tag=None,
    sj_strand_tag="XS",
    extra_tags=None,
):
    """Create a minimal mock pysam.AlignedSegment.

    Parameters
    ----------
    cigartuples : list of (op, length) or None
        CIGAR operations. Default: single 100bp match.
    xs_tag : str or None
        Splice-junction strand tag value. Only returned when queried.
    sj_strand_tag : str
        Name of the BAM tag that holds the SJ strand value
        (default ``"XS"``).
    extra_tags : dict or None
        Additional ``{tag_name: value}`` entries the mock should return.
    """
    read = MagicMock(spec=pysam.AlignedSegment)
    read.reference_name = ref_name
    read.reference_start = ref_start
    read.is_reverse = is_reverse
    if cigartuples is None:
        cigartuples = [(pysam.CMATCH, 100)]
    read.cigartuples = cigartuples

    _tag_name = sj_strand_tag
    _extras = extra_tags or {}

    def _get_tag(tag):
        if tag == _tag_name and xs_tag is not None:
            return xs_tag
        if tag in _extras:
            return _extras[tag]
        raise KeyError(tag)

    def _has_tag(tag):
        if tag == _tag_name:
            return xs_tag is not None
        return tag in _extras

    read.get_tag = _get_tag
    read.has_tag = _has_tag
    return read


# =====================================================================
# Fragment.from_reads
# =====================================================================


class TestFragmentFromReads:
    def test_single_read_unspliced(self):
        r1 = _mock_read(ref_start=100, cigartuples=[(pysam.CMATCH, 100)])
        frag = Fragment.from_reads([r1], [])

        assert len(frag.exons) == 1
        assert frag.exons[0].ref == "chr1"
        assert frag.exons[0].start == 100
        assert frag.exons[0].end == 200
        assert frag.exons[0].strand == Strand.POS
        assert len(frag.introns) == 0

    def test_single_read_spliced(self):
        """Read with 50M 200N 50M → one intron, two exon blocks."""
        r1 = _mock_read(
            ref_start=100,
            cigartuples=[
                (pysam.CMATCH, 50),
                (pysam.CREF_SKIP, 200),
                (pysam.CMATCH, 50),
            ],
            xs_tag="+",
        )
        frag = Fragment.from_reads([r1], [])

        assert len(frag.exons) == 2
        assert frag.exons[0] == GenomicInterval("chr1", 100, 150, Strand.POS)
        assert frag.exons[1] == GenomicInterval("chr1", 350, 400, Strand.POS)
        assert len(frag.introns) == 1
        assert frag.introns[0] == GenomicInterval("chr1", 150, 350, Strand.POS)

    def test_r2_strand_flip(self):
        """R2 strand is flipped: if R1 is forward and R2 is forward,
        after flip R2 becomes reverse → AMBIGUOUS."""
        r1 = _mock_read(ref_start=100, is_reverse=False)
        r2 = _mock_read(ref_start=300, is_reverse=False)  # same direction

        frag = Fragment.from_reads([r1], [r2])

        # R1 POS, R2 POS → after R2 flip: POS + NEG
        # Since they're keyed by (ref, strand), we get exon blocks for both strands
        strands = {e.strand for e in frag.exons}
        assert Strand.POS in strands
        assert Strand.NEG in strands

    def test_proper_fr_pair(self):
        """Proper FR pair: R1 forward, R2 reverse.
        After R2 flip (reverse → POS), both exon blocks have + strand."""
        r1 = _mock_read(ref_start=100, is_reverse=False)
        r2 = _mock_read(ref_start=300, is_reverse=True)

        frag = Fragment.from_reads([r1], [r2])

        # All exon blocks should be POS
        for exon in frag.exons:
            assert exon.strand == Strand.POS

    def test_overlapping_exons_merged(self):
        """Two reads with overlapping exon blocks → merged into one."""
        r1 = _mock_read(ref_start=100, is_reverse=False,
                        cigartuples=[(pysam.CMATCH, 150)])
        r2 = _mock_read(ref_start=200, is_reverse=True,
                        cigartuples=[(pysam.CMATCH, 100)])

        frag = Fragment.from_reads([r1], [r2])

        # After R2 flip (reverse → POS), both on same strand → merge
        assert len(frag.exons) == 1
        assert frag.exons[0].start == 100
        assert frag.exons[0].end == 300

    def test_both_empty(self):
        frag = Fragment.from_reads([], [])
        assert frag.exons == ()
        assert frag.introns == ()

    def test_r1_only(self):
        r1 = _mock_read(ref_start=50, cigartuples=[(pysam.CMATCH, 75)])
        frag = Fragment.from_reads([r1], [])
        assert len(frag.exons) == 1
        assert frag.exons[0].start == 50
        assert frag.exons[0].end == 125

    def test_r2_only(self):
        r2 = _mock_read(ref_start=50, is_reverse=True,
                        cigartuples=[(pysam.CMATCH, 75)])
        frag = Fragment.from_reads([], [r2])
        # R2 strand flipped: reverse → POS
        assert len(frag.exons) == 1
        assert frag.exons[0].strand == Strand.POS

    def test_introns_deduplicated(self):
        """Same splice junction from R1 and R2 should appear once."""
        cigar = [
            (pysam.CMATCH, 50),
            (pysam.CREF_SKIP, 100),
            (pysam.CMATCH, 50),
        ]
        r1 = _mock_read(ref_start=100, is_reverse=False, cigartuples=cigar, xs_tag="+")
        r2 = _mock_read(ref_start=100, is_reverse=True, cigartuples=cigar, xs_tag="+")

        frag = Fragment.from_reads([r1], [r2])

        # Intron set-based → deduplicated to 1
        assert len(frag.introns) == 1
        assert frag.introns[0] == GenomicInterval("chr1", 150, 250, Strand.POS)


# =====================================================================
# Fragment dataclass basics
# =====================================================================


class TestFragmentBasics:
    def test_default_is_empty(self):
        f = Fragment()
        assert f.exons == ()
        assert f.introns == ()
        assert f.insert_size is None


# =====================================================================
# Multi-tag / ts-tag support
# =====================================================================


class TestMultiTagSupport:
    """Tests for parse_read and Fragment.from_reads with multi-tag."""

    def test_parse_read_with_ts_tag(self):
        """Read carrying a 'ts' tag should be parsed when sj_strand_tag='ts'."""
        from hulkrna.bam import parse_read

        r = _mock_read(
            ref_start=100,
            cigartuples=[
                (pysam.CMATCH, 50),
                (pysam.CREF_SKIP, 200),
                (pysam.CMATCH, 50),
            ],
            xs_tag="-",
            sj_strand_tag="ts",
        )
        _, _, _, sjs = parse_read(r, sj_strand_tag="ts")
        assert len(sjs) == 1
        assert sjs[0][2] == Strand.NEG

    def test_parse_read_multi_tag_first_wins(self):
        """When multiple tags are given, the first present tag wins."""
        from hulkrna.bam import parse_read

        # Read has both XS and ts tags via extra_tags
        r = _mock_read(
            ref_start=100,
            cigartuples=[
                (pysam.CMATCH, 50),
                (pysam.CREF_SKIP, 200),
                (pysam.CMATCH, 50),
            ],
            xs_tag="+",
            sj_strand_tag="XS",
            extra_tags={"ts": "-"},
        )
        _, _, _, sjs = parse_read(r, sj_strand_tag=("XS", "ts"))
        # XS is first and present → should use "+"
        assert sjs[0][2] == Strand.POS

    def test_parse_read_multi_tag_fallback(self):
        """When first tag is missing, fall back to second."""
        from hulkrna.bam import parse_read

        # Read has only ts tag
        r = _mock_read(
            ref_start=100,
            cigartuples=[
                (pysam.CMATCH, 50),
                (pysam.CREF_SKIP, 200),
                (pysam.CMATCH, 50),
            ],
            xs_tag="-",
            sj_strand_tag="ts",
        )
        _, _, _, sjs = parse_read(r, sj_strand_tag=("XS", "ts"))
        # XS missing, ts present → should use "-"
        assert sjs[0][2] == Strand.NEG

    def test_parse_read_empty_tuple_gives_none(self):
        """Empty tag tuple → all SJs get Strand.NONE."""
        from hulkrna.bam import parse_read

        r = _mock_read(
            ref_start=100,
            cigartuples=[
                (pysam.CMATCH, 50),
                (pysam.CREF_SKIP, 200),
                (pysam.CMATCH, 50),
            ],
            xs_tag="+",
        )
        _, _, _, sjs = parse_read(r, sj_strand_tag=())
        assert sjs[0][2] == Strand.NONE

    def test_from_reads_multi_tag(self):
        """Fragment.from_reads respects multi-tag tuple."""
        r1 = _mock_read(
            ref_start=100,
            cigartuples=[
                (pysam.CMATCH, 50),
                (pysam.CREF_SKIP, 200),
                (pysam.CMATCH, 50),
            ],
            xs_tag="+",
            sj_strand_tag="ts",
        )
        frag = Fragment.from_reads([r1], [], sj_strand_tag=("XS", "ts"))
        assert len(frag.introns) == 1
        # XS absent, ts present → "+"
        assert frag.introns[0].strand == Strand.POS

    def test_slots(self):
        """Fragment uses __slots__ for memory efficiency."""
        f = Fragment()
        with pytest.raises(AttributeError):
            f.nonexistent_attr = 42


# =====================================================================
# Supplementary record merging
# =====================================================================


class TestSupplementaryMerging:
    """Fragment.from_reads merges primary + supplementary records."""

    def test_two_r1_records_merge_exons(self):
        """Primary R1 on chr1:100-200, supplementary R1 on chr1:5000-5100."""
        r1_pri = _mock_read(ref_start=100, cigartuples=[(pysam.CMATCH, 100)])
        r1_supp = _mock_read(ref_start=5000, cigartuples=[(pysam.CMATCH, 100)])
        frag = Fragment.from_reads([r1_pri, r1_supp], [])
        # Two non-overlapping exon blocks on same strand
        assert len(frag.exons) == 2
        starts = {e.start for e in frag.exons}
        assert 100 in starts
        assert 5000 in starts

    def test_supplementary_on_different_ref(self):
        """Supplementary on different chromosome gives exons on both refs."""
        r1_pri = _mock_read(ref_name="chr1", ref_start=100,
                            cigartuples=[(pysam.CMATCH, 100)])
        r1_supp = _mock_read(ref_name="chr5", ref_start=50000,
                             cigartuples=[(pysam.CMATCH, 100)])
        frag = Fragment.from_reads([r1_pri, r1_supp], [])
        refs = {e.ref for e in frag.exons}
        assert "chr1" in refs
        assert "chr5" in refs

    def test_supplementary_introns_merged(self):
        """Supplementary with splice junction → introns from both records merged."""
        cigar_spliced = [
            (pysam.CMATCH, 50),
            (pysam.CREF_SKIP, 200),
            (pysam.CMATCH, 50),
        ]
        r1_pri = _mock_read(ref_start=100, cigartuples=[(pysam.CMATCH, 100)])
        r1_supp = _mock_read(ref_start=5000, cigartuples=cigar_spliced, xs_tag="+")
        frag = Fragment.from_reads([r1_pri, r1_supp], [])
        # 3 exon blocks total (1 from primary, 2 from spliced supplementary)
        assert len(frag.exons) == 3
        assert len(frag.introns) == 1
