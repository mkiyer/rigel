"""The splice classification vocabulary shared by the C++ scanner and the Python side.

`SpliceType` and its census field names, the strand column each type maps to, and the column subset
built on them. A category added in one language and forgotten in the other reads as a zero rather
than as an error, so the enum and its derived names are pinned here.
"""

from rigel.splice import (
    SpliceType,
    SpliceStrandCol,
    NUM_SPLICE_STRAND_COLS,
    SPLICED_COLS,
)


# =====================================================================
# SpliceType
# =====================================================================


class TestSpliceType:
    def test_values(self):
        """One spliced type: a spliced fragment always has one annotated transcript holding every junction it observes
        (the scan's annotated-path rule), so there is no unannotated splice to name."""
        assert [(t.name, int(t)) for t in SpliceType] == [
            ("UNSPLICED", 0),
            ("SPLICED_ANNOT", 1),
            ("SPLICED_IMPLICIT", 2),
            ("SPLICE_ARTIFACT", 3),
        ]

    def test_native_constants_match_python_enum(self):
        from rigel._resolve_impl import testing as native_resolve

        assert native_resolve.SPLICE_UNSPLICED == int(SpliceType.UNSPLICED)
        assert not hasattr(native_resolve, "SPLICE_SPLICED_UNANNOT")
        assert native_resolve.SPLICE_SPLICED_ANNOT == int(SpliceType.SPLICED_ANNOT)
        assert native_resolve.SPLICE_IMPLICIT == int(SpliceType.SPLICED_IMPLICIT)
        assert native_resolve.SPLICE_ARTIFACT == int(SpliceType.SPLICE_ARTIFACT)

    def test_native_splice_labels_are_the_member_names(self):
        """The scanner's label for every splice class, which keys the census and stamps the annotated BAM's ZS, is the
        member name lower-cased, so ``census_field`` derives the census key in Python without a second table."""
        from rigel._bam_impl.testing import splice_type_label

        assert [splice_type_label(int(t)) for t in SpliceType] == [
            t.name.lower() for t in SpliceType
        ]
        assert splice_type_label(len(SpliceType)) == "unknown"


# =====================================================================
# SpliceStrandCol — 4 categories × 2 strands = 8 columns
# =====================================================================


class TestSpliceStrandCol:
    def test_num_splice_strand_cols(self):
        assert NUM_SPLICE_STRAND_COLS == 8

    def test_enum_values(self):
        assert [(c.name, int(c)) for c in SpliceStrandCol] == [
            ("UNSPLICED_SENSE", 0),
            ("UNSPLICED_ANTISENSE", 1),
            ("SPLICED_ANNOT_SENSE", 2),
            ("SPLICED_ANNOT_ANTISENSE", 3),
            ("SPLICED_IMPLICIT_SENSE", 4),
            ("SPLICED_IMPLICIT_ANTISENSE", 5),
            ("SPLICE_ARTIFACT_SENSE", 6),
            ("SPLICE_ARTIFACT_ANTISENSE", 7),
        ]

    def test_category_property(self):
        for c in SpliceStrandCol:
            assert c.category == SpliceType(c.value // 2)


# =====================================================================
# Pre-computed column subsets
# =====================================================================


class TestColumnSubsets:
    def test_spliced_cols(self):
        assert SPLICED_COLS == (2, 3)

    def test_spliced_cols_match_enum(self):
        expected = tuple(c.value for c in SpliceStrandCol if c.category == SpliceType.SPLICED_ANNOT)
        assert SPLICED_COLS == expected
