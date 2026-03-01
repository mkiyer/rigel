"""Tests for hulkrna.bam — BAM parser, hit grouping, and read pairing."""

from unittest.mock import MagicMock

import pysam
import pytest

from hulkrna.bam import (
    parse_bam_file,
    parse_read,
    _group_records_by_hit,
)
from hulkrna.types import Strand


# =====================================================================
# Mock read builder
# =====================================================================


def _mock_bam_read(
    query_name="frag1",
    *,
    ref_name="chr1",
    ref_id=0,
    ref_start=100,
    is_reverse=False,
    is_paired=True,
    is_read1=True,
    is_read2=False,
    is_qcfail=False,
    is_unmapped=False,
    is_secondary=False,
    is_supplementary=False,
    is_duplicate=False,
    is_proper_pair=True,
    mate_is_unmapped=False,
    next_ref_id=0,
    next_ref_start=300,
    nh=None,
    hi=None,
    cigartuples=None,
    xs_tag=None,
    sj_strand_tag="XS",
):
    """Build a mock pysam.AlignedSegment with all necessary attributes."""
    read = MagicMock(spec=pysam.AlignedSegment)
    read.query_name = query_name
    read.reference_name = ref_name
    read.reference_id = ref_id
    read.reference_start = ref_start
    read.is_reverse = is_reverse
    read.is_paired = is_paired
    read.is_read1 = is_read1
    read.is_read2 = is_read2 if not is_read1 else False
    if is_read2:
        read.is_read1 = False
        read.is_read2 = True
    read.is_qcfail = is_qcfail
    read.is_unmapped = is_unmapped
    read.is_secondary = is_secondary
    read.is_supplementary = is_supplementary
    read.is_duplicate = is_duplicate
    read.is_proper_pair = is_proper_pair
    read.mate_is_unmapped = mate_is_unmapped
    read.next_reference_id = next_ref_id
    read.next_reference_start = next_ref_start

    if cigartuples is None:
        cigartuples = [(pysam.CMATCH, 100)]
    read.cigartuples = cigartuples

    _tags = {}
    if nh is not None:
        _tags['NH'] = nh
    if hi is not None:
        _tags['HI'] = hi
    if xs_tag is not None:
        _tags[sj_strand_tag] = xs_tag

    def _has_tag(tag):
        return tag in _tags

    def _get_tag(tag):
        if tag in _tags:
            return _tags[tag]
        raise KeyError(tag)

    read.has_tag = _has_tag
    read.get_tag = _get_tag
    return read


def _make_pair(
    query_name="frag1",
    *,
    r1_start=100,
    r2_start=300,
    ref_id=0,
    nh=1,
    hi=None,
    is_secondary=False,
    is_supplementary=False,
    is_proper_pair=True,
):
    """Create a matched R1+R2 pair with reciprocal mate pointers."""
    r1 = _mock_bam_read(
        query_name, ref_start=r1_start, ref_id=ref_id,
        is_read1=True, is_read2=False,
        next_ref_id=ref_id, next_ref_start=r2_start,
        nh=nh, hi=hi,
        is_secondary=is_secondary,
        is_supplementary=is_supplementary,
        is_proper_pair=is_proper_pair,
    )
    r2 = _mock_bam_read(
        query_name, ref_start=r2_start, ref_id=ref_id,
        is_read1=False, is_read2=True, is_reverse=True,
        next_ref_id=ref_id, next_ref_start=r1_start,
        nh=nh, hi=hi,
        is_secondary=is_secondary,
        is_supplementary=is_supplementary,
        is_proper_pair=is_proper_pair,
    )
    return r1, r2


# =====================================================================
# parse_read
# =====================================================================


class TestParseRead:
    """Test CIGAR → exons/introns parsing."""

    def test_simple_match(self):
        r = _mock_bam_read(ref_start=100, cigartuples=[(pysam.CMATCH, 75)])
        ref, strand, exons, sjs = parse_read(r)
        assert ref == "chr1"
        assert strand == Strand.POS
        assert exons == [(100, 175)]
        assert sjs == []

    def test_reverse_strand(self):
        r = _mock_bam_read(ref_start=0, is_reverse=True,
                           cigartuples=[(pysam.CMATCH, 50)])
        _, strand, _, _ = parse_read(r)
        assert strand == Strand.NEG

    def test_spliced_read(self):
        r = _mock_bam_read(
            ref_start=100, xs_tag="+",
            cigartuples=[
                (pysam.CMATCH, 50),
                (pysam.CREF_SKIP, 200),
                (pysam.CMATCH, 50),
            ],
        )
        _, _, exons, sjs = parse_read(r)
        assert exons == [(100, 150), (350, 400)]
        assert len(sjs) == 1
        assert sjs[0] == (150, 350, Strand.POS)

    def test_deletion_in_exon(self):
        """Deletion (D) advances reference within exon."""
        r = _mock_bam_read(
            ref_start=100,
            cigartuples=[
                (pysam.CMATCH, 30),
                (pysam.CDEL, 5),
                (pysam.CMATCH, 30),
            ],
        )
        _, _, exons, sjs = parse_read(r)
        # One contiguous exon block spanning 30 + 5 + 30 = 65 bases
        assert exons == [(100, 165)]
        assert sjs == []


# =====================================================================
# _group_records_by_hit
# =====================================================================


class TestGroupRecordsByHit:

    def test_unique_mapper_one_hit(self):
        """NH=1 primary pair → one hit with one R1 and one R2."""
        r1, r2 = _make_pair(nh=1)
        hits, sec_r1, sec_r2 = _group_records_by_hit([r1, r2])
        assert len(hits) == 1
        assert sec_r1 == []
        assert sec_r2 == []
        r1_reads, r2_reads = hits[0]
        assert r1 in r1_reads
        assert r2 in r2_reads

    def test_primary_with_supplementary_merged(self):
        """Primary R1 + supplementary R1 → same hit."""
        r1_primary = _mock_bam_read(
            is_read1=True, ref_start=100,
        )
        r1_supp = _mock_bam_read(
            is_read1=True, ref_start=5000, is_supplementary=True,
        )
        r2_primary = _mock_bam_read(
            is_read1=False, is_read2=True, ref_start=300,
        )
        hits, sec_r1, sec_r2 = _group_records_by_hit([r1_primary, r1_supp, r2_primary])
        assert len(hits) == 1
        assert sec_r1 == []
        assert sec_r2 == []

    def test_hi_tag_grouping(self):
        """Records with HI tags are grouped by HI value."""
        r1a, r2a = _make_pair(r1_start=100, r2_start=300, hi=1)
        r1b, r2b = _make_pair(r1_start=500, r2_start=700, hi=2)
        hits, sec_r1, sec_r2 = _group_records_by_hit([r1a, r2a, r1b, r2b])
        assert len(hits) == 2
        assert sec_r1 == []
        assert sec_r2 == []
        # Sorted by HI
        assert r1a in hits[0][0]
        assert r1b in hits[1][0]

    def test_hi_tag_with_supplementary(self):
        """HI groups supplementary with its primary."""
        r1_primary = _mock_bam_read(
            is_read1=True, ref_start=100, hi=1,
        )
        r1_supp = _mock_bam_read(
            is_read1=True, ref_start=5000, is_supplementary=True, hi=1,
        )
        r2 = _mock_bam_read(
            is_read1=False, is_read2=True, ref_start=300, hi=1,
        )
        hits, sec_r1, sec_r2 = _group_records_by_hit([r1_primary, r1_supp, r2])
        assert len(hits) == 1
        assert sec_r1 == []
        assert sec_r2 == []
        r1_reads, r2_reads = hits[0]
        assert len(r1_reads) == 2
        assert len(r2_reads) == 1

    def test_fallback_secondary_returned_separately(self):
        """Without HI tags, all R1/R2 locations (including primary)
        go into sec lists for transcript-aware re-pairing."""
        r1_primary = _mock_bam_read(
            is_read1=True, ref_start=100, ref_id=0,
            next_ref_id=0, next_ref_start=300,
        )
        r2_primary = _mock_bam_read(
            is_read1=False, is_read2=True, ref_start=300, ref_id=0,
            next_ref_id=0, next_ref_start=100, is_reverse=True,
        )
        # Secondaries on a DIFFERENT reference (typical minimap2 multimap)
        r1_sec = _mock_bam_read(
            is_read1=True, ref_start=5000, ref_id=1, is_secondary=True,
            next_ref_id=0, next_ref_start=300,  # points to primary mate
        )
        r2_sec = _mock_bam_read(
            is_read1=False, is_read2=True, ref_start=5300, ref_id=1,
            is_secondary=True, is_reverse=True,
            next_ref_id=0, next_ref_start=100,  # points to primary mate
        )
        hits, sec_r1, sec_r2 = _group_records_by_hit(
            [r1_primary, r2_primary, r1_sec, r2_sec],
        )
        # No pre-formed hits — all locations in sec lists for re-pairing
        assert len(hits) == 0
        # Primary and secondary R1/R2 in sec lists
        assert len(sec_r1) == 2
        assert r1_primary in sec_r1[0]  # primary first
        assert r1_sec in sec_r1[1]
        assert len(sec_r2) == 2
        assert r2_primary in sec_r2[0]  # primary first
        assert r2_sec in sec_r2[1]

    def test_fallback_same_ref_secondaries_separated(self):
        """Without HI tags, secondaries on the same reference as primaries
        are returned in sec lists along with primaries for re-pairing."""
        r1_primary = _mock_bam_read(
            is_read1=True, ref_start=100, ref_id=0,
            next_ref_id=0, next_ref_start=300,
        )
        r2_primary = _mock_bam_read(
            is_read1=False, is_read2=True, ref_start=300, ref_id=0,
            next_ref_id=0, next_ref_start=100, is_reverse=True,
        )
        r1_sec = _mock_bam_read(
            is_read1=True, ref_start=5000, ref_id=0, is_secondary=True,
            next_ref_id=0, next_ref_start=300,
        )
        r2_sec = _mock_bam_read(
            is_read1=False, is_read2=True, ref_start=5300, ref_id=0,
            is_secondary=True, is_reverse=True,
            next_ref_id=0, next_ref_start=100,
        )
        hits, sec_r1, sec_r2 = _group_records_by_hit(
            [r1_primary, r2_primary, r1_sec, r2_sec],
        )
        # No pre-formed hits — all in sec lists for re-pairing
        assert len(hits) == 0
        assert len(sec_r1) == 2
        assert r1_primary in sec_r1[0]
        assert r1_sec in sec_r1[1]
        assert len(sec_r2) == 2
        assert r2_primary in sec_r2[0]
        assert r2_sec in sec_r2[1]

    def test_fallback_supplementary_grouped_with_primary(self):
        """Without HI, supplementary goes with primary (hit 0)."""
        r1_primary = _mock_bam_read(is_read1=True, ref_start=100)
        r1_supp = _mock_bam_read(
            is_read1=True, ref_start=5000, is_supplementary=True,
        )
        r2_primary = _mock_bam_read(
            is_read1=False, is_read2=True, ref_start=300,
        )
        hits, sec_r1, sec_r2 = _group_records_by_hit([r1_primary, r1_supp, r2_primary])
        assert len(hits) == 1
        assert sec_r1 == []
        assert sec_r2 == []
        assert len(hits[0][0]) == 2  # primary + supplementary R1
        assert r1_supp in hits[0][0]

    def test_empty_input(self):
        hits, sec_r1, sec_r2 = _group_records_by_hit([])
        assert hits == []
        assert sec_r1 == []
        assert sec_r2 == []


# =====================================================================
# parse_bam_file — core tests
# =====================================================================


class TestParseBamFile:
    """Core tests for parse_bam_file."""

    def _collect(self, bam_records, **kwargs):
        """Parse records and return (list_of_yielded, stats)."""
        stats = {}
        results = list(parse_bam_file(iter(bam_records), stats, **kwargs))
        return results, stats

    def test_unique_mapper_basic(self):
        """NH=1 primary pair → yield one group with one hit."""
        r1, r2 = _make_pair("frag1", nh=1)
        results, stats = self._collect([r1, r2])
        assert len(results) == 1
        nh, hits, sec_r1, sec_r2 = results[0]
        assert nh == 1
        assert len(hits) == 1
        r1_reads, r2_reads = hits[0]
        assert len(r1_reads) == 1
        assert len(r2_reads) == 1

    def test_unique_mapper_stats(self):
        """Stats for a simple unique mapper."""
        r1, r2 = _make_pair("frag1", nh=1)
        _, stats = self._collect([r1, r2])
        assert stats['total'] == 2
        assert stats['n_read_names'] == 1
        assert stats['unique'] == 1
        assert stats['multimapping'] == 0
        assert stats['proper_pair'] == 1

    def test_multimapper_skipped_by_default(self):
        """NH>1 groups are not yielded when include_multimap=False."""
        r1, r2 = _make_pair("frag1", nh=3)
        results, stats = self._collect([r1, r2], include_multimap=False)
        assert len(results) == 0
        assert stats['multimapping'] == 1
        assert stats['n_read_names'] == 1

    def test_multimapper_included(self):
        """NH>1 groups are yielded when include_multimap=True."""
        r1, r2 = _make_pair("frag1", nh=3)
        results, stats = self._collect([r1, r2], include_multimap=True)
        assert len(results) == 1
        nh, hits, *_ = results[0]
        assert nh == 3

    def test_multimapper_with_hi_tags(self):
        """Multimapper with HI tags → multiple hits."""
        r1a, r2a = _make_pair("frag1", r1_start=100, r2_start=300, nh=2, hi=1)
        r1b, r2b = _make_pair("frag1", r1_start=5000, r2_start=5300,
                               nh=2, hi=2, is_secondary=True)
        results, stats = self._collect(
            [r1a, r2a, r1b, r2b], include_multimap=True,
        )
        assert len(results) == 1
        nh, hits, *_ = results[0]
        assert nh == 2
        assert len(hits) == 2

    def test_multimapper_without_hi_tags_fallback(self):
        """Multimapper without HI → all R1/R2 locations in sec lists
        for transcript-aware re-pairing (no pre-formed primary hit)."""
        r1_pri = _mock_bam_read(
            "frag1", is_read1=True, ref_start=100, ref_id=0,
            next_ref_id=0, next_ref_start=300, nh=2,
        )
        r2_pri = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
            ref_start=300, ref_id=0, is_reverse=True,
            next_ref_id=0, next_ref_start=100, nh=2,
        )
        # Secondary pair on a DIFFERENT reference (typical minimap2 multimap)
        r1_sec = _mock_bam_read(
            "frag1", is_read1=True, ref_start=5000, ref_id=1,
            next_ref_id=0, next_ref_start=300, is_secondary=True, nh=2,
        )
        r2_sec = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
            ref_start=5300, ref_id=1, is_reverse=True,
            next_ref_id=0, next_ref_start=100, is_secondary=True, nh=2,
        )
        results, stats = self._collect(
            [r1_pri, r2_pri, r1_sec, r2_sec], include_multimap=True,
        )
        assert len(results) == 1
        nh, hits, sec_r1, sec_r2 = results[0]
        assert nh == 2
        # No pre-formed hits — all locations in sec lists for re-pairing
        assert len(hits) == 0
        # Primary + secondary R1/R2 in sec lists
        assert len(sec_r1) == 2
        assert r1_pri in sec_r1[0]  # primary first
        assert r1_sec in sec_r1[1]
        assert len(sec_r2) == 2
        assert r2_pri in sec_r2[0]  # primary first
        assert r2_sec in sec_r2[1]

    def test_supplementary_included_in_hit(self):
        """Supplementary records are included with primary in same hit."""
        r1_pri = _mock_bam_read(
            "frag1", is_read1=True, ref_start=100, nh=1,
        )
        r1_supp = _mock_bam_read(
            "frag1", is_read1=True, ref_start=5000,
            is_supplementary=True, nh=1,
        )
        r2_pri = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
            ref_start=300, is_reverse=True, nh=1,
        )
        results, stats = self._collect([r1_pri, r1_supp, r2_pri])
        assert len(results) == 1
        nh, hits, *_ = results[0]
        assert len(hits) == 1
        r1_reads, r2_reads = hits[0]
        assert len(r1_reads) == 2  # primary + supplementary
        assert len(r2_reads) == 1
        assert stats['supplementary'] == 1

    def test_supplementary_stats_counted(self):
        """Supplementary records contribute to stats['supplementary']."""
        r1_pri = _mock_bam_read("frag1", is_read1=True, nh=1)
        r1_supp = _mock_bam_read(
            "frag1", is_read1=True, ref_start=5000,
            is_supplementary=True, nh=1,
        )
        r2_pri = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
            is_reverse=True, nh=1,
        )
        _, stats = self._collect([r1_pri, r1_supp, r2_pri])
        assert stats['supplementary'] == 1
        assert stats['secondary'] == 0
        assert stats['total'] == 3

    def test_secondary_stats_counted(self):
        """Secondary records contribute to stats['secondary']."""
        r1_pri, r2_pri = _make_pair("frag1", nh=2)
        r1_sec = _mock_bam_read(
            "frag1", is_read1=True, ref_start=5000,
            is_secondary=True, nh=2,
        )
        _, stats = self._collect(
            [r1_pri, r2_pri, r1_sec], include_multimap=True,
        )
        assert stats['secondary'] == 1

    # -- Filtering tests --

    def test_qc_fail_filtered(self):
        """QC-failed records are excluded."""
        r1_ok, r2_ok = _make_pair("frag1", nh=1)
        r_qc = _mock_bam_read("frag1", is_qcfail=True)
        results, stats = self._collect([r_qc, r1_ok, r2_ok])
        assert stats['qc_fail'] == 1
        assert stats['total'] == 3
        assert len(results) == 1

    def test_unmapped_filtered(self):
        """Unmapped records are excluded."""
        r1_ok, r2_ok = _make_pair("frag1", nh=1)
        r_unmap = _mock_bam_read("frag1", is_unmapped=True)
        results, stats = self._collect([r_unmap, r1_ok, r2_ok])
        assert stats['unmapped'] == 1
        assert len(results) == 1

    def test_all_unmapped_no_yield(self):
        """If all records are unmapped, nothing is yielded."""
        r1 = _mock_bam_read("frag1", is_unmapped=True)
        r2 = _mock_bam_read("frag1", is_unmapped=True)
        results, stats = self._collect([r1, r2])
        assert len(results) == 0
        assert stats['unmapped'] == 2
        assert stats['n_read_names'] == 0

    def test_duplicate_skipped_by_default(self):
        """Duplicates are skipped when skip_duplicates=True (default)."""
        r1 = _mock_bam_read(
            "frag1", is_read1=True, is_duplicate=True, nh=1,
        )
        r2 = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True, is_duplicate=True, nh=1,
        )
        results, stats = self._collect([r1, r2], skip_duplicates=True)
        assert len(results) == 0
        assert stats['duplicate'] == 2

    def test_duplicate_kept(self):
        """Duplicates are kept when skip_duplicates=False."""
        r1, r2 = _make_pair("frag1", nh=1)
        r1.is_duplicate = True
        r2.is_duplicate = True
        results, stats = self._collect([r1, r2], skip_duplicates=False)
        assert len(results) == 1
        assert stats['duplicate'] == 2

    def test_unpaired_raises(self):
        """Unpaired read raises ValueError."""
        r = _mock_bam_read("frag1", is_paired=False)
        with pytest.raises(ValueError, match="paired-end"):
            self._collect([r])

    # -- Multiple read-name groups --

    def test_multiple_query_names(self):
        """Two distinct query names → two yielded groups."""
        r1a, r2a = _make_pair("fragA", nh=1)
        r1b, r2b = _make_pair("fragB", nh=1)
        results, stats = self._collect([r1a, r2a, r1b, r2b])
        assert len(results) == 2
        assert stats['n_read_names'] == 2
        assert stats['unique'] == 2

    def test_mixed_unique_and_multimap(self):
        """Unique + multimapper groups — only unique yielded by default."""
        r1a, r2a = _make_pair("fragA", nh=1)
        r1b, r2b = _make_pair("fragB", nh=3)
        results, stats = self._collect([r1a, r2a, r1b, r2b])
        assert len(results) == 1  # only unique
        assert stats['unique'] == 1
        assert stats['multimapping'] == 1

    # -- Pair-level stats --

    def test_proper_pair_stat(self):
        r1, r2 = _make_pair("frag1", nh=1, is_proper_pair=True)
        _, stats = self._collect([r1, r2])
        assert stats['proper_pair'] == 1
        assert stats['improper_pair'] == 0

    def test_improper_pair_stat(self):
        r1, r2 = _make_pair("frag1", nh=1, is_proper_pair=False)
        _, stats = self._collect([r1, r2])
        assert stats['improper_pair'] == 1
        assert stats['proper_pair'] == 0

    def test_mate_unmapped_stat(self):
        """R1-only hit → mate_unmapped counted."""
        r1 = _mock_bam_read("frag1", is_read1=True, nh=1)
        _, stats = self._collect([r1])
        assert stats['mate_unmapped'] == 1

    # -- NH tag propagation --

    def test_nh_tag_value(self):
        """The yielded nh matches the NH tag, not the number of records."""
        r1 = _mock_bam_read("frag1", is_read1=True, nh=5)
        r2 = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True, nh=5,
        )
        results, _ = self._collect([r1, r2], include_multimap=True)
        nh, hits, *_ = results[0]
        assert nh == 5

    def test_nh_defaults_to_1(self):
        """If NH tag is absent, nh defaults to 1."""
        r1 = _mock_bam_read(
            "frag1", is_read1=True,
        )  # No nh= kwarg → no NH tag
        r2 = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
        )
        results, stats = self._collect([r1, r2])
        nh, hits, *_ = results[0]
        assert nh == 1
        assert stats['unique'] == 1


# =====================================================================
# Integration: supplementary + secondary in same query name
# =====================================================================


class TestComplexReadGroups:
    """Test realistic scenarios with mixed record types."""

    def _collect(self, bam_records, **kwargs):
        stats = {}
        return list(parse_bam_file(iter(bam_records), stats, **kwargs)), stats

    def test_chimeric_unique_mapper_with_supplementary(self):
        """Fusion / chimeric read: primary on chr1, supplementary on chr5."""
        r1_pri = _mock_bam_read(
            "frag1", is_read1=True, ref_name="chr1", ref_start=100, nh=1,
        )
        r1_supp = _mock_bam_read(
            "frag1", is_read1=True, ref_name="chr5", ref_start=50000,
            is_supplementary=True, nh=1,
        )
        r2_pri = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
            ref_name="chr1", ref_start=300, is_reverse=True, nh=1,
        )
        results, stats = self._collect([r1_pri, r1_supp, r2_pri])

        assert len(results) == 1
        nh, hits, *_ = results[0]
        assert nh == 1
        assert len(hits) == 1
        r1_reads, r2_reads = hits[0]
        assert len(r1_reads) == 2
        # The supplementary is in the same hit
        ref_names = {r.reference_name for r in r1_reads}
        assert "chr1" in ref_names
        assert "chr5" in ref_names

    def test_multimapper_with_hi_and_supplementary(self):
        """Multimapper with HI tags, one hit has supplementary."""
        # Hit 1: primary pair on chr1
        r1a = _mock_bam_read(
            "frag1", is_read1=True, ref_start=100, nh=2, hi=1,
        )
        r2a = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
            ref_start=300, is_reverse=True, nh=2, hi=1,
        )
        # Hit 1: supplementary R1 (chimeric piece)
        r1a_supp = _mock_bam_read(
            "frag1", is_read1=True, ref_start=9000,
            is_supplementary=True, nh=2, hi=1,
        )
        # Hit 2: secondary pair on chr1
        r1b = _mock_bam_read(
            "frag1", is_read1=True, ref_start=5000,
            is_secondary=True, nh=2, hi=2,
        )
        r2b = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
            ref_start=5300, is_secondary=True, is_reverse=True,
            nh=2, hi=2,
        )
        results, stats = self._collect(
            [r1a, r2a, r1a_supp, r1b, r2b], include_multimap=True,
        )
        assert len(results) == 1
        nh, hits, *_ = results[0]
        assert nh == 2
        assert len(hits) == 2

        # Hit 1 (HI=1): 2 R1 records (primary + supp), 1 R2
        hit1_r1, hit1_r2 = hits[0]
        assert len(hit1_r1) == 2
        assert len(hit1_r2) == 1

        # Hit 2 (HI=2): 1 R1, 1 R2
        hit2_r1, hit2_r2 = hits[1]
        assert len(hit2_r1) == 1
        assert len(hit2_r2) == 1

    def test_all_qcfail_no_yield(self):
        """If every record in a query_name group is qcfail, skip."""
        r1 = _mock_bam_read("frag1", is_qcfail=True)
        r2 = _mock_bam_read("frag1", is_qcfail=True)
        results, stats = self._collect([r1, r2])
        assert len(results) == 0
        assert stats['qc_fail'] == 2
        assert stats['n_read_names'] == 0

    def test_secondary_only_r1_produces_singleton_hit(self):
        """Secondary R1 without matching R2 → singleton hit."""
        r1_pri = _mock_bam_read(
            "frag1", is_read1=True, ref_start=100, nh=2,
            next_ref_id=0, next_ref_start=300,
        )
        r2_pri = _mock_bam_read(
            "frag1", is_read1=False, is_read2=True,
            ref_start=300, is_reverse=True, nh=2,
            next_ref_id=0, next_ref_start=100,
        )
        r1_sec = _mock_bam_read(
            "frag1", is_read1=True, ref_start=9000,
            is_secondary=True, nh=2,
            next_ref_id=0, next_ref_start=99999,  # no matching R2
            mate_is_unmapped=True,
        )
        results, stats = self._collect(
            [r1_pri, r2_pri, r1_sec], include_multimap=True,
        )
        nh, hits, *_ = results[0]
        assert len(hits) == 2
        # Hit 1: secondary R1 only (mate_is_unmapped)
        sec_hit_r1, sec_hit_r2 = hits[1]
        assert len(sec_hit_r1) == 1
        assert len(sec_hit_r2) == 0
