"""The annotated-sj opportunity: the length tilt of the spliced pool (layer 2, opportunity).

The spliced pool (fragments whose path used an annotated sj) is a length-dependent selection: a longer fragment
crosses an sj more often, so the pool reads long against the library's RNA lengths. The tilt is computable from the
annotation. For a transcript with exon lengths ``e_1 .. e_K`` in transcript space and ``L = SUM e_i``, the starts at
which a length-``w`` window crosses at least one sj number::

    A_j(w)  =  (L - w + 1)+  -  SUM_i (e_i - w + 1)+

by the complement: a window crosses no sj iff it lies inside one exon, and the exons partition those placements, so
no inclusion-exclusion is needed. Over transcripts, weighted by molar abundance ``theta``::

    T(w) = SUM_t theta_t (L_t - w + 1)+          every placement
    A(w) = SUM_t theta_t A_j(w, t)               the placements that cross something
    pi(w) = A(w) / T(w)                          the probability a placement lands in the pool

and the corrected pool is ``pool(w) / pi(w)``. The divisor is ``pi``, not ``A``: ``A`` alone recovers the law lengths
were drawn from, while consumers need the law the library realizes, the drawn one weighted by its placements ``T``;
and because ``A`` and ``T`` sum over the same transcripts, the ratio tolerates a wrong ``theta``. ``theta`` must be
molar (copies), never a fragment count, which would apply the length weighting twice.
:func:`crossing_probability_from_index` weights every non-synthetic transcript equally, which needs no expression
estimate because the ratio cancels most of the dependence on ``theta``.

Gate: ``tests/calibration/test_sj_opportunity.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..index import TranscriptIndex

__all__ = [
    "crossing_probability",
    "crossing_probability_from_index",
    "detilt_pool",
    "sj_opportunity",
]


def _ramp_sum(lengths: np.ndarray, weights: np.ndarray, max_length: int) -> np.ndarray:
    """``SUM_j weights_j * (lengths_j - w + 1)+`` for every ``w`` in ``[0, max_length]``.

    ``(len - w + 1)+`` is the count of integers in ``[w, len]``, so the whole sum is
    ``SUM_{v >= w} SUM_{j: len_j >= v} weights_j``: the reverse cumulative sum, applied twice, of the weighted
    length histogram, with no loop over ``w``. The histogram must span both the longest object and ``max_length``: a
    shorter span truncates the inner tail sum, and one that stops before ``max_length`` returns a short, misaligned
    array.
    """
    lengths = np.asarray(lengths, dtype=np.int64)
    span = max(int(lengths.max(initial=0)) + 2, max_length + 1)
    hist = np.bincount(lengths, weights=np.asarray(weights, dtype=np.float64), minlength=span)
    tail = np.cumsum(hist[::-1])[::-1]
    out = np.cumsum(tail[::-1])[::-1][: max_length + 1]
    assert out.size == max_length + 1, (out.size, max_length)
    return out.copy()


def sj_opportunity(
    exon_lengths: np.ndarray,
    transcript_offsets: np.ndarray,
    theta: np.ndarray,
    max_length: int,
) -> tuple[np.ndarray, np.ndarray]:
    """``(T, A)``: total and sj-crossing placements, abundance-weighted, over ``w``.

    Parameters
    ----------
    exon_lengths
        int[n_exons] — every exon's length in transcript space, CSR-packed by transcript.
    transcript_offsets
        int[n_transcripts + 1] — the CSR offsets into ``exon_lengths``.
    theta
        float[n_transcripts] — molar abundance. A zero excludes the transcript entirely.
    max_length
        the largest fragment length to return, inclusive.

    Both curves come from the same two length spectra, one over transcripts and one over exons, so they agree on
    what a length is.
    """
    exon_lengths = np.asarray(exon_lengths, dtype=np.int64)
    transcript_offsets = np.asarray(transcript_offsets, dtype=np.int64)
    theta = np.asarray(theta, dtype=np.float64)
    if theta.size != transcript_offsets.size - 1:
        raise ValueError(
            f"theta has {theta.size} entries for {transcript_offsets.size - 1} transcripts"
        )

    n_exons = np.diff(transcript_offsets)
    # a prefix-sum difference, because `add.reduceat` mishandles an empty slice and a transcript with no cached
    # exons is an ordinary state
    prefix = np.concatenate([[0], np.cumsum(exon_lengths)])
    transcript_lengths = prefix[transcript_offsets[1:]] - prefix[transcript_offsets[:-1]]

    total = _ramp_sum(transcript_lengths, theta, max_length)
    inside_one_exon = _ramp_sum(exon_lengths, np.repeat(theta, n_exons), max_length)
    crossing = total - inside_one_exon

    # Bin 0 is not a fragment length, and the complement identity is negative there: a zero-length window sits on an
    # exon boundary, so `SUM_i (e_i + 1)` counts every internal boundary twice and `A_j(0) = 1 - K`. Both curves are
    # defined on `w >= 1`, so bin 0 is zeroed rather than left a negative divisor.
    total[0] = 0.0
    crossing[0] = 0.0
    return total, crossing


def crossing_probability(
    exon_lengths: np.ndarray,
    transcript_offsets: np.ndarray,
    theta: np.ndarray,
    max_length: int,
) -> np.ndarray:
    """``pi(w) = A(w) / T(w)``: the chance a uniformly placed length-``w`` fragment crosses an sj.

    Zero where there are no placements (``w`` beyond every expressed transcript): the pool cannot hold such a fragment.
    """
    total, crossing = sj_opportunity(exon_lengths, transcript_offsets, theta, max_length)
    return np.divide(crossing, total, out=np.zeros_like(total), where=total > 0.0)


def crossing_probability_from_index(index: "TranscriptIndex", max_length: int) -> np.ndarray:
    """``pi(w)`` for an index, weighting every non-synthetic transcript equally.

    The filter is ``~is_synthetic`` alone. Synthetic spans are not sequenced molecules, while a single-exon
    transcript is a real one and must stay in: it adds placements to ``T`` and none to ``A``.
    """
    offsets, exon_starts, exon_ends, _ = index.build_exon_csr()
    theta = (~index.t_df["is_synthetic"].to_numpy()).astype(np.float64)
    return crossing_probability(
        (np.asarray(exon_ends) - np.asarray(exon_starts)),
        np.asarray(offsets),
        theta,
        max_length,
    )


def detilt_pool(counts: np.ndarray, crossing_probability: np.ndarray) -> np.ndarray:
    """Divide a pool histogram by its selection probability, keeping its total.

    De-tilting changes the pool's shape, not how many fragments stand behind it, so the total is preserved for any
    consumer that reads it as an evidence weight. A bin whose probability is zero cannot legitimately hold mass (a
    fragment is in the pool only through a placement the probability counts); it contributes nothing and the surviving
    mass is renormalised over it. With no opportunity anywhere (every transcript single-exon, or an empty annotation)
    ``counts`` comes back unchanged: a correction with no information is inert.
    """
    counts = np.asarray(counts, dtype=np.float64)
    pi = np.asarray(crossing_probability, dtype=np.float64)
    if pi.size < counts.size:
        pi = np.concatenate([pi, np.zeros(counts.size - pi.size)])
    pi = pi[: counts.size]

    out = np.divide(counts, pi, out=np.zeros_like(counts), where=pi > 0.0)
    scale = out.sum()
    if scale <= 0.0:
        return counts.copy()
    return out * (counts.sum() / scale)
