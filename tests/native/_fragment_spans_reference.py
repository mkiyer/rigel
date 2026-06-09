"""
Pure-Python reference for ``fragment_genomic_spans`` — the contiguous-genomic-span
decomposition of a fragment (accumulator span redesign).

Canonical spec: ``docs/calibration/accumulator_fragment_span_redesign.md`` §3, §5.
The native C++ ``fragment_genomic_spans`` (landing in Phase B) must reproduce this
module's output on every case in ``tests/native/test_fragment_spans_spec.py``.

The rule, in one sentence:

    A fragment occupies, per reference, its full genomic extent
    ``[min block start, max block end]`` MINUS the introns it splices out
    (explicit CIGAR-N introns ∪ implicit introns detected in a mate gap).

Equivalently: fill every inter-block gap (the unsequenced insert is contiguous
DNA / within-exon sequence) EXCEPT where an intron must be cut. This makes the
deposited mass/flux reflect the molecule, not the sequenced read blocks — which
removes the paired-end over-count at its source.

Consequences that fall out of the rule:
- **Unspliced / gDNA** (no cut introns): one span ``[min, max]`` — mate gap filled.
- **Explicitly spliced**: introns sit exactly at block boundaries → spans are the
  blocks (within-read-merged).
- **Implicitly spliced**: the matched intron is cut out of the (larger) mate gap;
  the unsequenced *exonic* flanks on either side of the intron are filled.
- **Multi-reference** (chimeric): a different ``ref_id`` always starts a new span,
  so no span ever bridges references.

This file is test-only — do not import from ``src/rigel/``.
"""

from __future__ import annotations

# A block / intron / span is a plain ``(ref_id, start, end)`` tuple, half-open
# ``[start, end)``. ``cut_introns`` are the intron intervals to keep open
# (explicit ∪ implicit); the caller assembles them.
Interval = tuple[int, int, int]


def fragment_genomic_spans(
    blocks: list[Interval],
    cut_introns: list[Interval],
) -> list[Interval]:
    """Decompose a fragment's aligned blocks into contiguous genomic spans.

    Per reference, take the full extent ``[min start, max end]`` and subtract the
    cut introns falling inside it. Returns spans sorted by ``(ref_id, start)``.
    The number of spans is at most ``len(blocks)`` on a single reference (each
    subtracted intron splits one span into two only if it lies strictly inside).
    """
    if not blocks:
        return []

    # Group blocks by reference; a different ref never shares a span.
    spans: list[Interval] = []
    refs = sorted({b[0] for b in blocks})
    for ref in refs:
        ref_blocks = [b for b in blocks if b[0] == ref]
        lo = min(b[1] for b in ref_blocks)
        hi = max(b[2] for b in ref_blocks)
        if hi <= lo:
            continue

        # Cut introns on this reference that overlap the extent, sorted by start.
        cuts = sorted(
            (max(c[1], lo), min(c[2], hi))
            for c in cut_introns
            if c[0] == ref and min(c[2], hi) > max(c[1], lo)
        )

        # Subtract the (merged) cut intervals from [lo, hi).
        pos = lo
        for c0, c1 in cuts:
            if c0 > pos:
                spans.append((ref, pos, c0))
            pos = max(pos, c1)
        if pos < hi:
            spans.append((ref, pos, hi))

    spans.sort(key=lambda s: (s[0], s[1]))
    return spans
