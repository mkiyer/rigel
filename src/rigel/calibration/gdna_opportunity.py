"""The opportunities of the four gDNA fragment-length pools (layer 2, opportunity).

The gDNA length law is read from four pools of unspliced fragments (none of them pure gDNA; ``fl`` handles the
contamination), each a length-dependent selection that must be divided by its own opportunity before any are combined:

| pool | selection rule | dominant when |
|---|---|---|
| ``DNA_INTERGENIC`` | contained in exactly one intergenic region | capture OFF |
| ``DNA_INTRONIC`` | contained in exactly one intronic region | capture OFF |
| ``DNA_INTRON_EXON`` | crosses exactly one boundary, flanks {intron, exon} | capture ON |
| ``DNA_INTERGENIC_EXON`` | crosses exactly one boundary, flanks {intergenic, exon} | capture ON |

Off capture nearly all gDNA lies inside a large intergenic or intronic region; under capture the surviving gDNA sits
beside a probe, reaches the exon boundary and becomes crossing, so the contained pair alone would hold only the short
half of the law. The two opportunity functions are exact.

*Contained* in one of a set of regions with lengths ``ell_n``::

    A(w) = SUM_n (ell_n - w + 1)+

*Crossing exactly one boundary*, for a boundary whose flanking region lengths are ``a`` (left) and ``b`` (right)::

    A(w) = (w - 1)+  -  (w - 1 - a)+  -  (w - 1 - b)+  +  (w - 1 - a - b)+

A fragment ``[s, s+w)`` crosses the boundary for ``w - 1`` starts; of those it also crosses the previous boundary for
``(w-1-a)+``, the next for ``(w-1-b)+`` and both for ``(w-1-a-b)+``. Excluding the two nearest boundaries is exact,
since a fragment that reaches any farther boundary also crosses the nearest one. The reference ends need no special
case: the partition puts a bound at ``0`` and at ``L_ref``.

The divisor is the probability ``pi(w) = A(w)/T(w)``, with ``T(w) = SUM_refs (L_ref - w + 1)+`` the total admissible
gDNA starts, as for the sj pool (:mod:`rigel.calibration.sj_opportunity`): ``count/A`` recovers the law lengths were
drawn from, while consumers need the law the library realizes. The raw histograms must never share one divisor: the
contained opportunity falls with length and the crossing one rises, so a shared divisor reads the gDNA mean far too
long. Both functions assume gDNA is placed uniformly along the genome, which hybrid capture breaks, so the divisor
alone does not de-bias a capture-ON fit.

Gate: ``tests/calibration/test_gdna_opportunity.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .region_arrays import (
    boundary_reaches,
    boundary_region_indices,
    observable_pieces,
    observable_stretches,
)
from .signature import RegionType, coarse_type_array

if TYPE_CHECKING:
    from .region_arrays import RegionArrays

__all__ = [
    "GdnaOpportunity",
    "contained_opportunity",
    "crossing_opportunity",
    "observable_gdna_opportunity",
    "total_opportunity",
]


def _ramp(values: np.ndarray, max_width: int) -> np.ndarray:
    """``SUM_i (w - 1 - values[i])+`` for every ``w`` in ``[0, max_width]``, in O(max_width + n).

    Values above ``max_width - 2`` can never contribute, because a term is live only where
    ``x <= w - 2``, so they are clipped into a bin the cumulative sums never read.
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


def contained_opportunity(
    region_lengths: np.ndarray, max_width: int, weights: np.ndarray | None = None
) -> np.ndarray:
    """``SUM_n weight_n (ell_n - w + 1)+``: the starts at which a length-``w`` fragment fits inside one region, each
    region counted at its weight (once, with no weights).

    O(max_width + n): regions longer than ``max_width`` always contribute, so they reduce to
    ``sum(weight·ell) - (w - 1) * sum(weight)``, and the rest come from a weighted length histogram's reverse
    cumulative sums.
    """
    lengths = np.asarray(region_lengths, dtype=np.int64)
    weight = (
        np.ones(lengths.shape, dtype=np.float64)
        if weights is None
        else np.asarray(weights, dtype=np.float64)
    )
    width = np.arange(max_width + 1, dtype=np.float64)

    long = lengths > max_width
    out = float((weight * lengths)[long].sum()) - (width - 1.0) * float(weight[long].sum())

    short = lengths[~long]
    if short.size:
        histogram = np.bincount(short, weights=weight[~long], minlength=max_width + 2)
        index = np.arange(len(histogram), dtype=np.float64)
        # Reverse cumulative sums over lengths >= w.
        count_above = np.cumsum(histogram[::-1])[::-1]
        sum_above = np.cumsum((histogram * index)[::-1])[::-1]
        slot = np.arange(max_width + 1)
        out = out + (sum_above[slot] - (width - 1.0) * count_above[slot])
    return np.maximum(out, 0.0)


def crossing_opportunity(left: np.ndarray, right: np.ndarray, max_width: int) -> np.ndarray:
    """The starts at which a length-``w`` fragment crosses one boundary and no other, summed over boundaries.

    ``left`` and ``right`` are the flanking region lengths, one entry per boundary; the inclusion-exclusion over the
    two neighbours is the module docstring's.
    """
    left = np.asarray(left, dtype=np.int64)
    right = np.asarray(right, dtype=np.int64)
    width = np.arange(max_width + 1, dtype=np.float64)
    base = float(len(left)) * np.maximum(width - 1.0, 0.0)
    out = base - _ramp(left, max_width) - _ramp(right, max_width) + _ramp(left + right, max_width)
    return np.maximum(out, 0.0)


def total_opportunity(reference_lengths: np.ndarray, max_width: int) -> np.ndarray:
    """``T(w) = SUM_refs (L_ref - w + 1)+``: every admissible gDNA start in the reference.

    Every reference counts, RNA-only spike-ins included: the tool does not know which references carry genomic DNA,
    and their contribution is negligible by length.
    """
    return contained_opportunity(reference_lengths, max_width)


@dataclass(frozen=True, slots=True)
class GdnaOpportunity:
    """The four pools' opportunities and the total, all ``float64[max_width + 1]``.

    Annotation-derived and condition-independent: build it once per index, not once per library.
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
        """``pi(w) = [SUM_p A_p(w)] / T(w)``: the divisor for the four pools' summed counts.

        Summing counts and opportunities gives ``f(w) ~ [SUM_p count_p(w)] · T(w) / [SUM_p A_p(w)]``, the
        opportunity-weighted average of the four de-tilted pools. Under Poisson counts ``Var(count_p) ∝ A_p``, so the
        weights are exactly inverse-variance, with nothing tuned.
        """
        summed = np.sum(self.pools, axis=0)
        return np.divide(summed, self.total, out=np.zeros_like(summed), where=self.total > 0.0)


def observable_gdna_opportunity(region_arrays: "RegionArrays", max_width: int) -> GdnaOpportunity:
    """The four pools' opportunities and the total over the observable genome (`rigel.observable`): a start counts only
    where its fragment touches no masked base.

    The contained pools sum each region's observable pieces by the region's type. A crossing pool's flanks are the
    boundary's two regions up to the nearest masked base: a fragment reaching that base is not observed, just as one
    reaching the next boundary is not in the pool, and the inclusion-exclusion of :func:`crossing_opportunity` removes
    both. The total sums every reference's observable stretches. With nothing masked every pool is its regions'
    whole-length opportunity, to the bit.
    """
    types = coarse_type_array(region_arrays.signature)
    piece_region, piece_len = observable_pieces(region_arrays)
    piece_type = types[piece_region]
    lo, hi = boundary_region_indices(region_arrays.ref_id)
    reach_lo, reach_hi = boundary_reaches(region_arrays)
    size = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    left = np.minimum(size[lo], reach_lo).astype(np.int64)
    right = np.minimum(size[hi], reach_hi).astype(np.int64)
    pairs = np.sort(np.stack([types[lo], types[hi]]), axis=0)
    intron_exon = (pairs[0] == RegionType.INTRON) & (pairs[1] == RegionType.EXON)
    intergenic_exon = (pairs[0] == RegionType.INTERGENIC) & (pairs[1] == RegionType.EXON)
    return GdnaOpportunity(
        intergenic_contained=contained_opportunity(
            piece_len[piece_type == RegionType.INTERGENIC], max_width
        ),
        intronic_contained=contained_opportunity(
            piece_len[piece_type == RegionType.INTRON], max_width
        ),
        intron_exon_crossing=crossing_opportunity(left[intron_exon], right[intron_exon], max_width),
        intergenic_exon_crossing=crossing_opportunity(
            left[intergenic_exon], right[intergenic_exon], max_width
        ),
        total=total_opportunity(observable_stretches(region_arrays), max_width),
    )
