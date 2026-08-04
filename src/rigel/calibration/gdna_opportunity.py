"""rigel.calibration.gdna_opportunity — de-tilting the FOUR gDNA fragment-length pools.

The gDNA fragment-length model is fitted from pools that are pure gDNA **by construction**, and every
one of them is a length-dependent *selection*. Two are "contained in one node", which longer fragments
achieve less often; two are "crosses exactly one line", which longer fragments achieve more often. So
the four raw histograms tilt in **opposite directions**, and each has to be divided by its own
opportunity before any of them can be combined.

⭐ **The four pools, and why all four are needed.**

| pool | selection rule | dominant when |
|---|---|---|
| ``DNA_INTERGENIC`` | contained in exactly one intergenic node | capture OFF |
| ``DNA_INTRONIC`` | contained in exactly one intronic node | capture OFF |
| ``DNA_INTRON_EXON`` | crosses exactly one line, flanks {intron, exon} | capture ON |
| ``DNA_INTERGENIC_EXON`` | crosses exactly one line, flanks {intergenic, exon} | capture ON |

Off capture the library is spread over the genome and nearly all gDNA lies wholly inside a large
intergenic or intronic node. Under hybrid capture the surviving gDNA sits beside a probe, and a fragment
beside a probe **reaches** the exon boundary — so it stops being contained and becomes crossing. Fitting
from the contained pair alone therefore measures the short half of one population and reads ~15 % short
under capture; the crossing pair holds the long half.

⭐ **The two opportunity functions, both exact.**

*Contained* in one of a set of nodes with lengths ``ell_n``::

    A(w) = SUM_n (ell_n - w + 1)+

*Crossing exactly one line*, for a line whose flanking node lengths are ``a`` (left) and ``b`` (right)::

    A(w) = (w - 1)+  -  (w - 1 - a)+  -  (w - 1 - b)+  +  (w - 1 - a - b)+

A fragment ``[s, s+w)`` crosses the line for ``w - 1`` starts; of those it also crosses the previous
line for ``(w-1-a)+`` and the next for ``(w-1-b)+``, and both for ``(w-1-a-b)+``. ⭐ **The two nearest
lines are the only ones that need excluding** — a fragment is an interval containing the line, so if it
reaches any line beyond ``p-a`` it must also cross ``p-a``. The inclusion-exclusion over the two
neighbours is therefore exact rather than a truncation. ⚠ And the reference ends need no special case:
the partition cuts at ``0`` and at ``L_ref``, so the outermost node's length *is* the distance to the
wall and the same subtraction removes the impossible starts.

⛔ **DIVIDE BY THE PROBABILITY, NEVER BY ``A`` ALONE.** ``count(w)/A(w)`` recovers the distribution
lengths were *drawn* from; every consumer needs the one the library *realizes*, which is the drawn one
weighted by how many placements each length has. So the divisor is ``pi(w) = A(w)/T(w)`` with::

    T(w) = SUM_refs (L_ref - w + 1)+

the total admissible gDNA starts in the reference. ⚠ On whole chromosomes ``T`` is flat to ~1 part in
10^5 and the two forms coincide numerically; on a short reference they do not, and the probability form
is the correct one either way. (This is the same rule the junction pool obeys —
:mod:`rigel.calibration.junction_opportunity`.)

⭐⭐ **AND THE COMBINATION IS DERIVED, NOT CHOSEN.** Sum the counts, sum the opportunities::

    f(w)  ~  [SUM_p count_p(w)] * T(w) / [SUM_p A_p(w)]

That is algebraically the **opportunity-weighted average** of the four de-tilted pools —
``SUM_p A_p f_p / SUM_p A_p`` — and under Poisson counts ``Var(count_p) ∝ A_p``, so weights
proportional to ``A_p`` are exactly inverse-variance. There is no tunable weight anywhere in it.

⛔ **What this is NOT: pooling the four histograms raw.** Summing four differently-tilted counts and
applying one divisor is the defect the shipped model used to have — it read a gDNA mean of 146.05 where
the pure contained pool said 88.0. Summing counts *and* the matching per-pool opportunities is a
different operation with a different answer, and ``tests/calibration/test_gdna_opportunity.py``
separates the two.

⚠ **The residual this does NOT fix, stated so nobody expects it to.** Both opportunity functions assume
gDNA is placed **uniformly** along the genome. Under hybrid capture it is not: placement is proportional
to the capture landscape, which the tool cannot see. Measured on the pilot, the four-pool model lands at
**−0.01 %** against truth off capture and **+7.9 %** under it, against a contained-pair model's −0.53 %
and −14.8 %. The remaining +7.9 % is the non-uniform density, not the opportunity, and closing it needs
a capture-aware placement model rather than a better divisor.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..index import TranscriptIndex

__all__ = [
    "GdnaOpportunity",
    "contained_opportunity",
    "crossing_opportunity",
    "gdna_opportunity_from_index",
    "total_opportunity",
]

#: Coarse node types, as :func:`rigel.calibration.signature.coarse_type_array` emits them.
_TYPE_INTERGENIC, _TYPE_INTRON, _TYPE_EXON = 0, 1, 2


def _ramp(values: np.ndarray, max_width: int) -> np.ndarray:
    """``SUM_i (w - 1 - values[i])+`` for every ``w`` in ``[0, max_width]``, in O(max_width + n).

    ⚠ Values above ``max_width - 2`` can never contribute — a term is live only where ``x <= w - 2`` —
    so they are clipped into a bin the cumulative sums never read.
    """
    clipped = np.clip(np.asarray(values, dtype=np.int64), 0, max_width + 1)
    histogram = np.bincount(clipped, minlength=max_width + 2).astype(np.float64)
    index = np.arange(len(histogram), dtype=np.float64)
    count_below = np.cumsum(histogram)
    sum_below = np.cumsum(histogram * index)

    width = np.arange(max_width + 1, dtype=np.int64)
    out = np.zeros(max_width + 1, dtype=np.float64)
    live = width >= 2
    m = width[live] - 2
    out[live] = (width[live] - 1) * count_below[m] - sum_below[m]
    return out


def contained_opportunity(node_lengths: np.ndarray, max_width: int) -> np.ndarray:
    """``SUM_n (ell_n - w + 1)+`` — the starts at which a length-``w`` fragment fits inside one node.

    ⭐ O(max_width + n), not O(max_width x n): nodes longer than ``max_width`` always contribute, so
    they reduce to ``sum(ell) - (w - 1) * count``, and the rest come from a length histogram's reverse
    cumulative sums. At human scale ``n`` is a million nodes and the naive form is a billion terms.
    """
    lengths = np.asarray(node_lengths, dtype=np.int64)
    width = np.arange(max_width + 1, dtype=np.float64)

    long = lengths > max_width
    out = float(lengths[long].sum()) - (width - 1.0) * float(long.sum())

    short = lengths[~long]
    if short.size:
        histogram = np.bincount(short, minlength=max_width + 2).astype(np.float64)
        index = np.arange(len(histogram), dtype=np.float64)
        # Reverse cumulative sums over lengths >= w.
        count_above = np.cumsum(histogram[::-1])[::-1]
        sum_above = np.cumsum((histogram * index)[::-1])[::-1]
        slot = np.arange(max_width + 1)
        out = out + (sum_above[slot] - (width - 1.0) * count_above[slot])
    return np.maximum(out, 0.0)


def crossing_opportunity(left: np.ndarray, right: np.ndarray, max_width: int) -> np.ndarray:
    """``SUM_lines`` starts at which a length-``w`` fragment crosses THAT line and no other.

    ``left`` and ``right`` are the flanking node lengths, one entry per line. See the module docstring
    for the inclusion-exclusion and why two neighbours suffice.
    """
    left = np.asarray(left, dtype=np.int64)
    right = np.asarray(right, dtype=np.int64)
    width = np.arange(max_width + 1, dtype=np.float64)
    base = float(len(left)) * np.maximum(width - 1.0, 0.0)
    out = base - _ramp(left, max_width) - _ramp(right, max_width) + _ramp(left + right, max_width)
    return np.maximum(out, 0.0)


def total_opportunity(reference_lengths: np.ndarray, max_width: int) -> np.ndarray:
    """``T(w) = SUM_refs (L_ref - w + 1)+`` — every admissible gDNA start in the reference.

    ⚠ Every reference counts, including RNA-only spike-ins: the tool does not know which references
    carry genomic DNA and must not pretend to. Their contribution is negligible by length anyway.
    """
    return contained_opportunity(reference_lengths, max_width)


@dataclass(frozen=True, slots=True)
class GdnaOpportunity:
    """The four pools' opportunities and the total, all ``float64[max_width + 1]``.

    ⚠ Annotation-derived and condition-independent: build it once per index, not once per library.
    """

    intergenic_contained: np.ndarray
    intronic_contained: np.ndarray
    intron_exon_crossing: np.ndarray
    intergenic_exon_crossing: np.ndarray
    total: np.ndarray

    @property
    def pools(self) -> tuple[np.ndarray, ...]:
        """The four in ``rigel.scan_payload`` pool order, so a caller cannot mis-pair them."""
        return (
            self.intergenic_contained,
            self.intronic_contained,
            self.intron_exon_crossing,
            self.intergenic_exon_crossing,
        )

    def combined_probability(self) -> np.ndarray:
        """``pi(w) = [SUM_p A_p(w)] / T(w)`` — the divisor for the four pools' summed counts.

        ⭐ This is the whole model in one line, and the weighting inside it is inverse-variance rather
        than chosen; see the module docstring.
        """
        summed = np.sum(self.pools, axis=0)
        return np.divide(summed, self.total, out=np.zeros_like(summed), where=self.total > 0.0)


def gdna_opportunity_from_index(index: "TranscriptIndex", max_width: int) -> GdnaOpportunity:
    """Build the four pools' opportunities from an index's node partition alone.

    ⚠ Reads the same ``build_node_partition_arrays`` axis the accumulator deposits onto, so the divisor
    and the deposit rule cannot drift apart — a divisor derived from a *different* view of the partition
    is how a pool comes to be divided by an opportunity it does not have.
    """
    from .splice_graph import build_node_partition_arrays

    cuts, ref_cut_offsets, node_types = build_node_partition_arrays(index)
    cuts = np.asarray(cuts, dtype=np.int64)
    offsets = np.asarray(ref_cut_offsets, dtype=np.int64)
    types = np.asarray(node_types, dtype=np.int64)

    node_lengths: list[np.ndarray] = []
    line_left: list[np.ndarray] = []
    line_right: list[np.ndarray] = []
    line_pairs: list[np.ndarray] = []
    reference_lengths: list[int] = []
    node_base = 0
    for r in range(len(offsets) - 1):
        lo, hi = int(offsets[r]), int(offsets[r + 1])
        if hi - lo < 2:
            # A reference with no nodes contributes no cuts, so it cannot host a gDNA fragment either.
            continue
        reference_cuts = cuts[lo:hi]
        lengths = np.diff(reference_cuts)
        n_nodes = len(lengths)
        reference_lengths.append(int(reference_cuts[-1]))
        node_lengths.append(lengths)
        if n_nodes >= 2:
            # Interior lines only: line i sits between node i-1 and node i.
            line_left.append(lengths[:-1])
            line_right.append(lengths[1:])
            flanks = np.stack(
                [
                    types[node_base : node_base + n_nodes - 1],
                    types[node_base + 1 : node_base + n_nodes],
                ]
            )
            # Sorted, so the flank pair is order-insensitive exactly as the accumulator's own
            # `_SPLASH_POOL` key is.
            line_pairs.append(np.sort(flanks, axis=0))
        node_base += n_nodes

    if not node_lengths:
        zero = np.zeros(max_width + 1, dtype=np.float64)
        return GdnaOpportunity(zero, zero.copy(), zero.copy(), zero.copy(), zero.copy())

    lengths = np.concatenate(node_lengths)
    types_by_node = types[: len(lengths)]
    left = np.concatenate(line_left) if line_left else np.zeros(0, np.int64)
    right = np.concatenate(line_right) if line_right else np.zeros(0, np.int64)
    pairs = np.concatenate(line_pairs, axis=1) if line_pairs else np.zeros((2, 0), np.int64)

    intron_exon = (pairs[0] == _TYPE_INTRON) & (pairs[1] == _TYPE_EXON)
    intergenic_exon = (pairs[0] == _TYPE_INTERGENIC) & (pairs[1] == _TYPE_EXON)

    return GdnaOpportunity(
        intergenic_contained=contained_opportunity(
            lengths[types_by_node == _TYPE_INTERGENIC], max_width
        ),
        intronic_contained=contained_opportunity(lengths[types_by_node == _TYPE_INTRON], max_width),
        intron_exon_crossing=crossing_opportunity(left[intron_exon], right[intron_exon], max_width),
        intergenic_exon_crossing=crossing_opportunity(
            left[intergenic_exon], right[intergenic_exon], max_width
        ),
        total=total_opportunity(np.asarray(reference_lengths, dtype=np.int64), max_width),
    )
