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
* :func:`gdna_eff_len_for_loci` — FL-marginal overlap effective length
    for the per-``MultiLocus`` gDNA EM component.

Owning these in one module keeps numerator/denominator geometry consistent
across :mod:`density_global` and :mod:`locus_prior`. See
``docs/calibration/locoregional_gdna_redesign_2026-05-07.md`` §3.1.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import numpy as np

from ..frag_length_model import FragmentLengthModel

if TYPE_CHECKING:  # pragma: no cover - import only for type hints
    from ._regional_exposure import RegionalGdnaExposure


__all__ = [
    "contained_exposure_clipped",
    "boundary_crossing_exposure",
    "boundary_side_in_window",
    "gdna_eff_len_for_loci",
    "weighted_gdna_eff_len_for_loci",
]


def _ref_length_for_locus(
    ref_lengths: Mapping[str | int, int] | Sequence[int],
    locus,
) -> int:
    """Return reference length for a Locus-like object."""
    if isinstance(ref_lengths, Mapping):
        if locus.ref in ref_lengths:
            return int(ref_lengths[locus.ref])
        if locus.ref_id in ref_lengths:
            return int(ref_lengths[locus.ref_id])
        raise KeyError(f"Reference length missing for ref={locus.ref!r}, ref_id={locus.ref_id!r}.")

    ref_id = int(locus.ref_id)
    if ref_id < 0 or ref_id >= len(ref_lengths):
        raise KeyError(f"Reference id {ref_id} out of bounds for ref_lengths.")
    return int(ref_lengths[ref_id])


def gdna_eff_len_for_loci(
    loci: tuple | list,
    ref_lengths: Mapping[str | int, int] | Sequence[int],
    fl: FragmentLengthModel,
    *,
    min_value: float = 1.0,
) -> float:
    """FL-marginal gDNA overlap effective length for a ``MultiLocus``.

    Counts the fragment start positions whose length-``ell`` genomic
    interval overlaps any interval in ``loci`` and averages that count over
    the gDNA fragment-length PMF. For one half-open interval ``[a, b)`` away
    from contig ends, the valid start window is ``[a - ell + 1, b)`` and the
    count is ``(b - a) + ell - 1``.

    Multi-interval loci are handled exactly by merging expanded start windows
    per reference for each fragment length, avoiding double-counts when one
    start position overlaps multiple nearby locus intervals.
    """
    if min_value < 0.0:
        raise ValueError(f"gdna_eff_len_for_loci: min_value must be >= 0, got {min_value}.")
    if not loci:
        return float(min_value)

    pmf = fl.pmf
    positive_ell = np.flatnonzero(pmf[1:] > 0.0) + 1
    if positive_ell.size == 0:
        return float(min_value)

    # Fast path: one interval whose expanded start window never clips against
    # contig boundaries for any observed fragment length.
    if len(loci) == 1:
        locus = loci[0]
        span = max(int(locus.end) - int(locus.start), 0)
        ref_len = _ref_length_for_locus(ref_lengths, locus)
        max_ell = int(positive_ell[-1])
        if int(locus.start) - max_ell + 1 >= 0 and int(locus.end) <= ref_len - max_ell + 1:
            probs = pmf[positive_ell]
            positive_mass = float(probs.sum())
            positive_mean = float(np.dot(positive_ell.astype(np.float64), probs))
            eff = float(span) * positive_mass + positive_mean - positive_mass
            return max(eff, float(min_value))

    total = 0.0
    for ell in positive_ell:
        ell_i = int(ell)
        by_ref: dict[int | str, list[tuple[int, int]]] = {}
        for locus in loci:
            ref_len = _ref_length_for_locus(ref_lengths, locus)
            valid_hi = ref_len - ell_i + 1
            if valid_hi <= 0:
                continue
            lo = max(int(locus.start) - ell_i + 1, 0)
            hi = min(int(locus.end), valid_hi)
            if hi <= lo:
                continue
            ref_key = int(locus.ref_id) if hasattr(locus, "ref_id") else locus.ref
            by_ref.setdefault(ref_key, []).append((lo, hi))

        n_valid = 0
        for windows in by_ref.values():
            windows.sort()
            cur_lo, cur_hi = windows[0]
            for lo, hi in windows[1:]:
                if lo <= cur_hi:
                    if hi > cur_hi:
                        cur_hi = hi
                else:
                    n_valid += cur_hi - cur_lo
                    cur_lo, cur_hi = lo, hi
            n_valid += cur_hi - cur_lo

        total += float(pmf[ell_i]) * float(n_valid)

    return max(total, float(min_value))


def weighted_gdna_eff_len_for_loci(
    loci: tuple | list,
    ref_lengths: Mapping[str | int, int] | Sequence[int],
    fl: FragmentLengthModel,
    exposure: "RegionalGdnaExposure",
    *,
    min_value: float = 1.0,
) -> float:
    """Regional-exposure-weighted analogue of :func:`gdna_eff_len_for_loci`.

    Computes :math:`\\sum_\\ell h(\\ell) \\sum_w \\int A(x)\\,dx`, where ``w``
    ranges over the merged per-reference midpoint windows for fragment
    length ``\\ell`` and ``A(x)`` is the per-region exposure weight from
    ``exposure``. The midpoint of a length-``\\ell`` fragment starting at
    position ``s`` is ``s + \\ell // 2``; the merged start window
    ``[lo, hi)`` therefore maps to the midpoint window
    ``[lo + \\ell // 2, hi + \\ell // 2)``.

    For ``exposure.mode == "uniform"`` this delegates to the unweighted
    function and is bit-exact equal.
    """
    if min_value < 0.0:
        raise ValueError(
            f"weighted_gdna_eff_len_for_loci: min_value must be >= 0, got {min_value}."
        )
    if exposure.mode == "uniform":
        return gdna_eff_len_for_loci(loci, ref_lengths, fl, min_value=min_value)
    if not loci:
        return float(min_value)

    pmf = fl.pmf
    positive_ell = np.flatnonzero(pmf[1:] > 0.0) + 1
    if positive_ell.size == 0:
        return float(min_value)

    total = 0.0
    for ell in positive_ell:
        ell_i = int(ell)
        half = ell_i // 2
        by_ref: dict[int, list[tuple[int, int]]] = {}
        for locus in loci:
            ref_len = _ref_length_for_locus(ref_lengths, locus)
            valid_hi = ref_len - ell_i + 1
            if valid_hi <= 0:
                continue
            lo = max(int(locus.start) - ell_i + 1, 0)
            hi = min(int(locus.end), valid_hi)
            if hi <= lo:
                continue
            ref_id = int(locus.ref_id)
            by_ref.setdefault(ref_id, []).append((lo, hi))

        for ref_id, windows in by_ref.items():
            windows.sort()
            cur_lo, cur_hi = windows[0]
            for lo, hi in windows[1:]:
                if lo <= cur_hi:
                    if hi > cur_hi:
                        cur_hi = hi
                else:
                    total += float(pmf[ell_i]) * exposure.weighted_length_on_ref(
                        ref_id, cur_lo + half, cur_hi + half
                    )
                    cur_lo, cur_hi = lo, hi
            total += float(pmf[ell_i]) * exposure.weighted_length_on_ref(
                ref_id, cur_lo + half, cur_hi + half
            )

    return max(total, float(min_value))


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
            f"contained_exposure_clipped: starts.shape {starts.shape} != ends.shape {ends.shape}"
        )
    spans_full = np.maximum(ends - starts, 0)
    clip_starts = np.maximum(starts, np.int64(clip_lo))
    clip_ends = np.minimum(ends, np.int64(clip_hi))
    spans_clip = np.maximum(clip_ends - clip_starts, 0)
    eff_full = fl.compute_all_transcript_eff_lens(spans_full, min_value=0.0)
    eff_clip = fl.compute_all_transcript_eff_lens(spans_clip, min_value=0.0)
    return eff_full, eff_clip


def boundary_crossing_exposure(
    fl: FragmentLengthModel, *, splicing_anchor_tolerance: int = 0
) -> float:
    """Expected number of fragment start positions that *strictly* cross
    a single boundary, under the gDNA fragment-length distribution.

    Returns :math:`B_\\text{cross}(K) = \\sum_\\ell h(\\ell)\\,\\max(\\ell - 2q(K) + 1, 0)`
    where :math:`q(K) = \\max(K, 1)` and ``K = splicing_anchor_tolerance``.

    The :math:`q(K) = \\max(K, 1)` term preserves the strict-crossing
    semantics at ``K = 0``: a fragment must still have at least 1 bp on
    each side of the boundary, yielding :math:`B_\\text{cross}(0) = E[L] - 1`
    (when :math:`L \\ge 1` almost surely). At ``K \\geq 1`` the formula
    extends to "at least K bp on each side," matching the scanner-side
    qualification predicate in
    :class:`rigel.calibration::CalibrationAccumulator`.

    Returns ``0.0`` if the sum underflows (degenerate PMF concentrated
    at :math:`\\ell \\le 2q(K) - 1`); callers must treat a zero exposure
    as "no boundary information available" and fall back accordingly.

    Parameters
    ----------
    fl : FragmentLengthModel
        Finalized FL model. Its ``pmf`` attribute supplies the per-bp
        probability mass.
    splicing_anchor_tolerance : int, optional
        Minimum bp clearance ``K`` on each side of a boundary required
        for a fragment to count as a boundary-crossing event. Default 0
        reproduces the pre-2026.05 strict-crossing semantics.

    Raises
    ------
    ValueError
        If ``splicing_anchor_tolerance < 0``.
    """
    if splicing_anchor_tolerance < 0:
        raise ValueError(
            f"boundary_crossing_exposure: splicing_anchor_tolerance "
            f"({splicing_anchor_tolerance}) must be >= 0"
        )
    pmf = fl.pmf
    ell = np.arange(pmf.size, dtype=np.float64)
    q = float(max(int(splicing_anchor_tolerance), 1))
    val = float((pmf * np.maximum(ell - 2.0 * q + 1.0, 0.0)).sum())
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
            f"boundary_side_in_window: starts.shape {starts.shape} != ends.shape {ends.shape}"
        )
    lo = np.int64(clip_lo)
    hi = np.int64(clip_hi)
    left_in = (starts >= lo) & (starts <= hi)
    right_in = (ends >= lo) & (ends <= hi)
    return left_in, right_in
