"""rigel.calibration._exposure — FL-aware geometric exposure primitives.

Single owner of the *fragment-length-aware geometry* used by the v6
locoregional gDNA prior. Core helpers:

* :func:`contained_exposure_clipped` — per-region full and clipped
  contained-fragment effective lengths, sharing the salmon-style
  eCDF cache in :class:`FragmentLengthModel`. Used by the locus
  proration `ratio_r = eff_clip_r / eff_full_r`.
* :func:`fractional_boundary_side_exposure` — FL-PMF-weighted per-side
  exposure on a vector of region lengths (one side of one boundary per
  region; replaces the legacy ``boundary_crossing_exposure`` scalar).

* :func:`boundary_side_in_window` — boolean per-region flags marking
  whether each region's left and right boundaries lie inside a query
  window. Used to localize boundary events to a Locus.
* :func:`footprint_exposure_weight` — scalar bp-weighted exposure over a
    merged component footprint, used to attenuate the gDNA EM denominator.
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
    "fractional_boundary_side_exposure",
    "boundary_side_in_window",
    "footprint_exposure_weight",
    "transcript_exposure_weights",
    "gdna_eff_len_for_loci",
]


def _merged_blocks(
    blocks: Sequence[tuple[int, int, int]],
) -> list[tuple[int, int, int]]:
    """Return non-overlapping ``(ref_id, start, end)`` blocks."""
    valid_blocks = [
        (int(ref_id), int(start), int(end))
        for ref_id, start, end in blocks
        if int(end) > int(start)
    ]
    if not valid_blocks:
        return []
    valid_blocks.sort()
    merged: list[tuple[int, int, int]] = []
    cur_ref, cur_start, cur_end = valid_blocks[0]
    for ref_id, start, end in valid_blocks[1:]:
        if ref_id == cur_ref and start <= cur_end:
            if end > cur_end:
                cur_end = end
        else:
            merged.append((cur_ref, cur_start, cur_end))
            cur_ref, cur_start, cur_end = ref_id, start, end
    merged.append((cur_ref, cur_start, cur_end))
    return merged


def footprint_exposure_weight(
    blocks: Sequence[tuple[int, int, int]],
    exposure: "RegionalGdnaExposure",
    *,
    min_weight: float = 1.0e-4,
) -> float:
    """Return bp-weighted mean exposure over a component footprint.

    ``blocks`` are half-open genomic intervals as ``(ref_id, start, end)``.
    Overlapping blocks on the same reference are merged before integration so
    locus footprints and exon blocks are not double-counted.
    """
    if min_weight < 0.0:
        raise ValueError(f"footprint_exposure_weight: min_weight must be >= 0, got {min_weight}.")
    merged = _merged_blocks(blocks)
    if not merged:
        return 1.0
    raw_bp = float(sum(end - start for _, start, end in merged))
    if raw_bp <= 0.0:
        return 1.0
    if exposure.mode == "uniform":
        return 1.0
    weighted_bp = sum(
        exposure.weighted_length_on_ref(ref_id, start, end) for ref_id, start, end in merged
    )
    weight = float(weighted_bp) / raw_bp
    return float(np.clip(weight, min_weight, 1.0))


def transcript_exposure_weights(
    index,
    exposure: "RegionalGdnaExposure",
    *,
    min_weight: float = 1.0e-4,
) -> np.ndarray:
    """Return per-transcript EM exposure weights.

    nRNA rows use the unspliced genomic span. Other transcript rows use exon
    blocks, which skips introns for mature mRNA components.
    """
    n_t = int(index.num_transcripts)
    weights = np.ones(n_t, dtype=np.float64)
    if exposure.mode == "uniform" or n_t == 0:
        return weights
    if getattr(index, "t_to_ref_arr", None) is None:
        raise RuntimeError("TranscriptIndex.t_to_ref_arr not populated")

    is_nrna = index.t_df["is_nrna"].to_numpy(dtype=bool)
    starts = index.t_df["start"].to_numpy(dtype=np.int64)
    ends = index.t_df["end"].to_numpy(dtype=np.int64)
    ref_ids = np.asarray(index.t_to_ref_arr, dtype=np.int32)
    exon_offsets, exon_starts, exon_ends, _ = index.build_exon_csr()

    for t_idx in range(n_t):
        ref_id = int(ref_ids[t_idx])
        if is_nrna[t_idx]:
            blocks = [(ref_id, int(starts[t_idx]), int(ends[t_idx]))]
        else:
            lo = int(exon_offsets[t_idx])
            hi = int(exon_offsets[t_idx + 1])
            if hi > lo:
                blocks = [
                    (ref_id, int(exon_starts[pos]), int(exon_ends[pos])) for pos in range(lo, hi)
                ]
            else:
                blocks = [(ref_id, int(starts[t_idx]), int(ends[t_idx]))]
        weights[t_idx] = footprint_exposure_weight(
            blocks,
            exposure,
            min_weight=min_weight,
        )
    return weights


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


def fractional_boundary_side_exposure(
    lengths_bp: np.ndarray,
    gdna_fl: FragmentLengthModel,
) -> np.ndarray:
    """FL-PMF-weighted per-side boundary exposure for a vector of region lengths.

    For each region length :math:`S_r`, returns

    .. math::

        E^{side}_r = \\sum_\\ell h(\\ell)\\,\\min((\\ell - 1)/2, S_r/2),

    the expected mass emitted to one side of one region boundary under
    the fractional accumulator routing rule. Large regions reduce to
    ``(E[FL] - 1) / 2``; regions shorter than the fragment length cap at
    ``S_r / 2`` because fully spanned regions split mass across both sides.

    Parameters
    ----------
    lengths_bp : np.ndarray
        Per-region span lengths (bp), shape ``(R,)``.
    gdna_fl : FragmentLengthModel
        Finalized gDNA fragment-length model.

    Returns
    -------
    np.ndarray
        ``float64`` per-region exposures, shape ``(R,)``. Entries are
        clipped at 0.
    """
    lengths = np.maximum(np.asarray(lengths_bp, dtype=np.int64), 0)
    pmf = np.asarray(gdna_fl.pmf, dtype=np.float64)
    if pmf.size == 0:
        return np.zeros(lengths.shape, dtype=np.float64)

    ell_minus_one = np.maximum(np.arange(pmf.size, dtype=np.float64) - 1.0, 0.0)
    positive_mass = pmf.copy()
    positive_mass[0] = 0.0
    cdf = np.cumsum(positive_mass)
    first_moment = np.cumsum(positive_mass * ell_minus_one)
    total_mass = float(cdf[-1])

    threshold = np.minimum(lengths + 1, pmf.size - 1)
    low_moment = first_moment[threshold]
    tail_mass = total_mass - cdf[threshold]
    return 0.5 * (low_moment + lengths.astype(np.float64, copy=False) * tail_mass)


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
