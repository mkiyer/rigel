"""rigel.calibration._exposure — FL-aware geometric exposure primitives.

Single owner of the *fragment-length-aware geometry* used by the v6
locoregional gDNA prior. Three primitives:

* :func:`contained_exposure_clipped` — per-region full and clipped
  contained-fragment effective lengths, sharing the salmon-style
  eCDF cache in :class:`FragmentLengthModel`. Used by the locus
  proration `ratio_r = eff_clip_r / eff_full_r`.
* :func:`boundary_crossing_exposure` — :math:`B_\\text{cross} =
  \\sum_\\ell h(\\ell)\\,\\max(\\ell - 1, 0)`, the per-side opportunity
  space for a fragment to *strictly* cross a single boundary.
* :func:`boundary_side_in_window` — boolean per-region flags marking
  whether each region's left and right boundaries lie inside a query
  window. Used to localize boundary events to a Locus.

Owning these in one module keeps numerator/denominator geometry consistent
across :mod:`density_global` and :mod:`locus_prior`. See
``docs/calibration/locoregional_gdna_redesign_2026-05-07.md`` §3.1.
"""

from __future__ import annotations

import numpy as np

from ..frag_length_model import FragmentLengthModel


__all__ = [
    "contained_exposure_clipped",
    "boundary_crossing_exposure",
    "boundary_side_in_window",
]


def contained_exposure_clipped(
    starts: np.ndarray,
    ends: np.ndarray,
    clip_lo: int,
    clip_hi: int,
    fl: FragmentLengthModel,
) -> tuple[np.ndarray, np.ndarray]:
    """FL-PMF-weighted contained-fragment effective length, full and clipped.

    For each region ``[start_r, end_r)``, returns:

    * ``eff_full[r]`` — :func:`l_eff_contained` of the *full* region span.
    * ``eff_clip[r]`` — same, after clipping the span to
      ``[clip_lo, clip_hi]``.

    The ratio ``eff_clip / eff_full`` is the FL-aware fraction of a
    region's contained-fragment opportunity space that falls inside the
    locus window. Pair it with the *whole-region* contained count emitted
    by the C++ scanner to get the FL-aware locus-local contained count
    without double-counting fragments that nominally start outside the
    locus.

    Both outputs floor at 0.0 (a span of 0 has zero opportunity).

    Parameters
    ----------
    starts, ends : np.ndarray
        Per-region half-open intervals (int64 or castable). Must have
        identical shape.
    clip_lo, clip_hi : int
        Closed-open clip window applied element-wise.
    fl : FragmentLengthModel
        Finalized fragment-length model. Must have non-zero mean.

    Returns
    -------
    (eff_full, eff_clip) : tuple of np.ndarray
        Two ``float64`` arrays of the same shape as ``starts``.
    """
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    if starts.shape != ends.shape:
        raise ValueError(
            f"contained_exposure_clipped: starts.shape {starts.shape} != "
            f"ends.shape {ends.shape}"
        )
    spans_full = np.maximum(ends - starts, 0)
    clip_starts = np.maximum(starts, np.int64(clip_lo))
    clip_ends = np.minimum(ends, np.int64(clip_hi))
    spans_clip = np.maximum(clip_ends - clip_starts, 0)
    eff_full = fl.compute_all_transcript_eff_lens(spans_full, min_value=0.0)
    eff_clip = fl.compute_all_transcript_eff_lens(spans_clip, min_value=0.0)
    return eff_full, eff_clip


def boundary_crossing_exposure(fl: FragmentLengthModel) -> float:
    """Expected number of fragment start positions that *strictly* cross
    a single boundary, under the gDNA fragment-length distribution.

    Returns :math:`B_\\text{cross} = \\sum_\\ell h(\\ell)\\,\\max(\\ell - 1, 0)`.
    Equals :math:`E[L] - 1` when :math:`L \\ge 1` almost surely, but the
    explicit ``max(ℓ - 1, 0)`` keeps the estimator robust to PMF mass at
    ``ℓ ∈ {0, 1}``. Returns ``0.0`` if the sum underflows (degenerate
    PMF concentrated at ``ℓ ≤ 1``); callers must treat a zero exposure
    as "no boundary information available" and fall back accordingly.
    """
    pmf = fl.pmf
    ell = np.arange(pmf.size, dtype=np.float64)
    val = float((pmf * np.maximum(ell - 1.0, 0.0)).sum())
    return val if val > 0.0 else 0.0


def boundary_side_in_window(
    starts: np.ndarray,
    ends: np.ndarray,
    clip_lo: int,
    clip_hi: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Boolean per-region flags: does each region's left/right boundary
    lie inside the *closed* window ``[clip_lo, clip_hi]``?

    A region's *left* boundary lives at ``start``; its *right* at ``end``.
    Both ends are tested with closed (``>=``, ``<=``) inequalities so a
    region whose end coincides with ``clip_hi`` still contributes a right
    side to the locus. Used to restrict EXON–INTRON boundary-flux events
    to those whose crossing point falls inside the query Locus interval.

    Returns
    -------
    (left_in, right_in) : tuple of np.ndarray
        Two ``bool`` arrays of the same shape as ``starts``.
    """
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    if starts.shape != ends.shape:
        raise ValueError(
            f"boundary_side_in_window: starts.shape {starts.shape} != "
            f"ends.shape {ends.shape}"
        )
    lo = np.int64(clip_lo)
    hi = np.int64(clip_hi)
    left_in = (starts >= lo) & (starts <= hi)
    right_in = (ends >= lo) & (ends <= hi)
    return left_in, right_in
