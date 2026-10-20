"""The splice classification vocabulary shared by the C++ scanner and the Python side.

`SpliceType` and its census field names, the strand column each type maps to, and the column subset
built on them. A category added in one language and forgotten in the other reads as a zero rather
than as an error, so the enum and its derived names are pinned here.
"""

from pathlib import Path


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
        (the scan's annotated-path rule, `tests/test_hit_verdict.py`), so there is no unannotated splice to name."""
        assert [(t.name, int(t)) for t in SpliceType] == [
            ("UNSPLICED", 0),
            ("SPLICED_ANNOT", 1),
            ("SPLICED_IMPLICIT", 2),
            ("SPLICE_ARTIFACT", 3),
        ]

    def test_native_constants_match_python_enum(self):
        import rigel._resolve_impl as native_resolve

        assert native_resolve.SPLICE_UNSPLICED == int(SpliceType.UNSPLICED)
        assert not hasattr(native_resolve, "SPLICE_SPLICED_UNANNOT")
        assert native_resolve.SPLICE_SPLICED_ANNOT == int(SpliceType.SPLICED_ANNOT)
        assert native_resolve.SPLICE_IMPLICIT == int(SpliceType.SPLICED_IMPLICIT)
        assert native_resolve.SPLICE_ARTIFACT == int(SpliceType.SPLICE_ARTIFACT)

    def test_native_bam_splice_labels_cover_extended_codes(self):
        """The native BAM writer's private ZS label helper must not fall through.

        ``splice_type_label`` is a static C++ helper, so this source-level test is
        the lightest direct guard that extended splice classes stay labelled in
        annotated BAM output.
        """
        source = Path(__file__).parents[1] / "src/rigel/native/bam_scanner.cpp"
        text = source.read_text()
        assert 'case SPLICE_IMPLICIT:        return "spliced_implicit";' in text
        assert 'case SPLICE_ARTIFACT:        return "splice_artifact";' in text


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
