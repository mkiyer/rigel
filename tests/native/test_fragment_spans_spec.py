"""Spec tests for ``fragment_genomic_spans`` — the contiguous-genomic-span rule.

Phase A (test-first) of the accumulator span redesign
(``docs/calibration/accumulator_fragment_span_redesign.md``). The reference
(`_fragment_spans_reference.py`) is the canonical contract; these table-driven
cases pin every fragment class. The native-parity test is skipped until the C++
``fragment_genomic_spans`` is bound (Phase B).
"""

from __future__ import annotations

import pytest

from ._fragment_spans_reference import fragment_genomic_spans

# Each case: (id, blocks, cut_introns, expected_spans). Intervals are
# (ref_id, start, end), half-open.
CASES = [
    # --- gDNA / unspliced: fill every gap ---
    (
        "unspliced_single_block",
        [(0, 100, 250)],
        [],
        [(0, 100, 250)],
    ),
    (
        "unspliced_paired_mate_gap_filled",  # R1 [100,200), gap, R2 [400,500)
        [(0, 100, 200), (0, 400, 500)],
        [],
        [(0, 100, 500)],  # the unsequenced insert is contiguous DNA → one span
    ),
    (
        "unspliced_overlapping_blocks_merged",
        [(0, 100, 300), (0, 250, 400)],
        [],
        [(0, 100, 400)],
    ),
    # --- explicitly spliced: intron sits exactly at the block boundary ---
    (
        "explicit_splice_two_exons",
        [(0, 100, 200), (0, 300, 400)],
        [(0, 200, 300)],  # CIGAR-N intron == the gap
        [(0, 100, 200), (0, 300, 400)],  # spans are the exons
    ),
    (
        "explicit_splice_three_exons",
        [(0, 100, 200), (0, 300, 400), (0, 600, 700)],
        [(0, 200, 300), (0, 400, 600)],
        [(0, 100, 200), (0, 300, 400), (0, 600, 700)],
    ),
    # --- implicitly spliced: intron lives inside a larger mate gap; flanks fill ---
    (
        "implicit_splice_flanks_filled",  # R1 [100,200), R2 [4900,5000); intron [300,4800)
        [(0, 100, 200), (0, 4900, 5000)],
        [(0, 300, 4800)],
        [(0, 100, 300), (0, 4800, 5000)],  # exonic flanks [200,300) and [4800,4900) filled
    ),
    # --- mixed: a within-exon mate gap (fill) AND a real intron (cut) ---
    (
        "mate_gap_filled_intron_cut",
        # blocks: exonA two reads [100,180)+[220,300) (mate gap 180-220), then exonB [600,700)
        [(0, 100, 180), (0, 220, 300), (0, 600, 700)],
        [(0, 300, 600)],  # only the intron is cut; the 180-220 mate gap fills
        [(0, 100, 300), (0, 600, 700)],
    ),
    # --- multi-reference (chimeric-like): never bridge refs ---
    (
        "multi_ref_separate_spans",
        [(0, 100, 200), (1, 5000, 5100)],
        [],
        [(0, 100, 200), (1, 5000, 5100)],
    ),
    # --- empty ---
    (
        "no_blocks",
        [],
        [],
        [],
    ),
]


@pytest.mark.parametrize("case", CASES, ids=[c[0] for c in CASES])
def test_reference_span_rule(case):
    _id, blocks, cuts, expected = case
    assert fragment_genomic_spans(blocks, cuts) == expected


def test_span_count_bounded_by_blocks():
    """The redesign guarantees #spans ≤ #blocks (the perf/scratch bound)."""
    for _id, blocks, cuts, _expected in CASES:
        spans = fragment_genomic_spans(blocks, cuts)
        assert len(spans) <= max(len(blocks), 0) or not blocks


def test_unspliced_is_single_extent():
    """Any unspliced (no-cut) fragment collapses to exactly one span per ref."""
    blocks = [(0, 100, 200), (0, 350, 420), (0, 900, 1000)]
    assert fragment_genomic_spans(blocks, []) == [(0, 100, 1000)]


@pytest.mark.skip(reason="Phase B: native fragment_genomic_spans not yet bound")
def test_native_parity():
    """Native C++ must match the reference byte-for-byte on every spec case.

    Un-skip when Phase B binds ``rigel.native.fragment_genomic_spans``.
    """
    from rigel.native import fragment_genomic_spans as native_spans  # type: ignore

    for _id, blocks, cuts, _expected in CASES:
        assert list(native_spans(blocks, cuts)) == fragment_genomic_spans(blocks, cuts)
