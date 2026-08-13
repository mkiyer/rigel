"""rigel.calibration.junction_opportunity — de-tilting the annotated-junction length pool.

The RNA fragment-length model is fitted from ``RNA_SPLICED``, the pool of deposited fragments whose
path used an **annotated junction**. That is a selection, and it is length-dependent: a longer
fragment covers more of its transcript, so it crosses a junction more often. Measured against the
simulator's own realized mRNA lengths on the chr22 pilot, the pool runs **+6.2 % to +8.1 %** long
while the unconditional anchor runs +0.00 %.

⭐ **The tilt is exactly computable from the annotation.** For a transcript with exon lengths
``e_1 .. e_K`` in transcript space and total ``L = SUM e_i``, the number of start positions at which a
length-``w`` window crosses **at least one** junction is::

    A_j(w)  =  (L - w + 1)+  -  SUM_i (e_i - w + 1)+

⭐ Work with the **complement** and it decomposes with no inclusion-exclusion: a window crosses no
junction iff it lies wholly inside a single exon, and the exons are disjoint, so those placements
partition by which exon contains the window. Attempting the union of "crosses junction ``j``" events
directly is messy and gets worse with every exon; the complement is exact in one boundary.

The library-level quantities are abundance-weighted sums of that, over transcripts::

    T(w) = SUM_t theta_t (L_t - w + 1)+          every placement
    A(w) = SUM_t theta_t A_j(w, t)               the placements that cross something
    pi(w) = A(w) / T(w)                          the probability a placement lands in the pool

and the corrected pool is ``pool(w) / pi(w)``.

⛔ **Divide by ``pi``, never by ``A`` alone, and the difference is not cosmetic.** ``A`` alone recovers
the distribution lengths were *drawn* from; what every consumer needs is the distribution the library
*realizes*, which is the drawn one weighted by how many placements each length has — and that weight
is ``T``. The ratio also makes the correction **safe under a wrong ``theta``**: ``A`` and ``T`` are
sums over the same transcripts, so a reweighting moves both. Measured over a theta sweep including
deliberately pathological regimes (all the abundance on the single steepest-tilted transcript), the
ratio form never does worse than not correcting and the ``A``-only form does.
``tests/calibration/test_junction_opportunity.py`` pins both halves of that.

⚠ **``theta`` is a molar abundance — copies — not an observed fragment count.** ``A_j`` already counts
start positions, so a count applies the length weighting twice.

⭐ **And production uses a UNIFORM theta over the non-synthetic transcripts, which needs no expression
estimate at all.** That is not a shortcut taken for cheapness: swept on the pilot, uniform lands the
corrected pool within 0.2 pp of what the simulator's own molar abundances achieve, because the ratio
cancels most of the dependence. The remaining sensitivity is measured in the test module.
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
    "junction_opportunity",
]


def _ramp_sum(lengths: np.ndarray, weights: np.ndarray, max_length: int) -> np.ndarray:
    """``SUM_j weights_j * (lengths_j - w + 1)+`` for every ``w`` in ``[0, max_length]``.

    ``(len - w + 1)+`` is the count of integers in ``[w, len]``, so the whole sum is
    ``SUM_{v >= w} SUM_{j: len_j >= v} weights_j`` — the reverse cumulative sum, applied twice, of the
    weighted length histogram. Exact integer geometry with no loop over ``w`` and no float rounding
    beyond the weights themselves.

    ⚠ The histogram must span the longest object present **and** ``max_length``, not whichever is
    smaller: too short a span truncates the inner tail sum and every value comes back too small, and
    a spectrum that stops before ``max_length`` returns a short array that silently misaligns against
    the other curve.
    """
    lengths = np.asarray(lengths, dtype=np.int64)
    span = max(int(lengths.max(initial=0)) + 2, max_length + 1)
    hist = np.bincount(lengths, weights=np.asarray(weights, dtype=np.float64), minlength=span)
    tail = np.cumsum(hist[::-1])[::-1]
    out = np.cumsum(tail[::-1])[::-1][: max_length + 1]
    assert out.size == max_length + 1, (out.size, max_length)
    return out.copy()


def junction_opportunity(
    exon_lengths: np.ndarray,
    transcript_offsets: np.ndarray,
    theta: np.ndarray,
    max_length: int,
) -> tuple[np.ndarray, np.ndarray]:
    """``(T, A)``: total and junction-crossing placements, abundance-weighted, over ``w``.

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

    ⭐ Both curves are built from the same two length spectra — one over transcripts, one over exons —
    so they cannot disagree about what a length is.
    """
    exon_lengths = np.asarray(exon_lengths, dtype=np.int64)
    transcript_offsets = np.asarray(transcript_offsets, dtype=np.int64)
    theta = np.asarray(theta, dtype=np.float64)
    if theta.size != transcript_offsets.size - 1:
        raise ValueError(
            f"theta has {theta.size} entries for {transcript_offsets.size - 1} transcripts"
        )

    n_exons = np.diff(transcript_offsets)
    # ⚠ A prefix-sum difference rather than `add.reduceat`, which has no answer for an EMPTY slice —
    # and a transcript with no cached exons is an ordinary state, not an error.
    prefix = np.concatenate([[0], np.cumsum(exon_lengths)])
    transcript_lengths = prefix[transcript_offsets[1:]] - prefix[transcript_offsets[:-1]]

    total = _ramp_sum(transcript_lengths, theta, max_length)
    inside_one_exon = _ramp_sum(exon_lengths, np.repeat(theta, n_exons), max_length)
    crossing = total - inside_one_exon

    # ⛔ Bin 0 is not a fragment length, and the complement identity is NEGATIVE there: a zero-length
    # window sits exactly ON an exon boundary, so `SUM_i (e_i + 1)` counts every internal boundary
    # twice and `A_j(0) = 1 - K`. Both curves are defined on `w >= 1`; leaving a negative divisor in
    # bin 0 would be a landmine for any consumer that does not happen to guard on it.
    total[0] = 0.0
    crossing[0] = 0.0
    return total, crossing


def crossing_probability(
    exon_lengths: np.ndarray,
    transcript_offsets: np.ndarray,
    theta: np.ndarray,
    max_length: int,
) -> np.ndarray:
    """``pi(w) = A(w) / T(w)`` — the chance a uniformly placed length-``w`` fragment crosses a junction.

    Zero where there are no placements at all (``w`` beyond every expressed transcript), which is the
    honest answer: the pool cannot contain such a fragment, so there is nothing to de-tilt.
    """
    total, crossing = junction_opportunity(exon_lengths, transcript_offsets, theta, max_length)
    return np.divide(crossing, total, out=np.zeros_like(total), where=total > 0.0)


def crossing_probability_from_index(index: "TranscriptIndex", max_length: int) -> np.ndarray:
    """``pi(w)`` for an index, weighting every **real** transcript equally.

    ⛔ The transcript filter is ``~is_synthetic``, alone. The manufactured nascent spans are not
    molecules anybody sequenced, and giving them opportunity puts weight on exon structures the
    library does not contain. ⚠ It is **not** ``~is_synthetic & ~is_nrna``: on a real row ``is_nrna``
    means "single-exon, so mature is nascent", and using it as a realness filter silently deletes
    real transcripts.
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
    """Divide a junction-pool histogram by its own opportunity, keeping its evidence weight.

    ⚠ **The total is preserved on purpose.** ``build_fl_models`` shrinks each pool toward the anchor
    with a Dirichlet pseudo-count whose strength is the pool total, and the pool total means "how many
    fragments stand behind this shape". De-tilting changes the shape, not how much evidence there is;
    letting the total move would weaken the shrinkage as an accidental side effect.

    ⭐ A bin where the opportunity is zero cannot legitimately hold mass — a fragment is only in this
    pool because it crossed an annotated junction, which is a placement the opportunity counts. Such a
    bin contributes nothing and the surviving mass is renormalised over it, so no evidence is lost.

    Returns ``counts`` unchanged when the annotation offers no opportunity anywhere (every transcript
    single-exon, or an empty annotation): a correction with no information must be inert.
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
