"""Effective lengths: the divisor that turns a count into a start density (layer 2, opportunity).

An effective length is the expected number of admissible fragment start positions at an object under a
fragment-length pmf ``f`` over molecule length ``w``. It is pure geometry, applied per component because gDNA and
RNA have different length laws and different templates::

    contained   E_f[ (region_len - w + 1)+ ]                              fits wholly inside a region
    crossing    E_f[ max(0, min(w-1, R_lo, R_hi, R_lo + R_hi - w + 1)) ]  spans a 0-bp boundary

``R_lo`` / ``R_hi`` are the molecule's own remaining template either side of the boundary, so reach is per
component: gDNA's template is the chromosome (:data:`UNBOUNDED_REACH`, and the divisor collapses to ``mu - 1``),
RNA's ends with its transcript, which an sj passes as the annotation's exonic reach. At a contiguous
boundary the unspliced RNA reach is :data:`UNBOUNDED_REACH`, a value rather than a second code path. The conserved
share (:func:`conserved_cut_shares`) reads the same placements a third way: how much of each placement's unit the
deposit rule gives an object, so a template's objects share its fl-marginal length with nothing counted twice. The
``observable_*`` forms read all three over the observable genome (a start whose fragment touches a masked base is no
place at all); with nothing masked each equals its plain form to the bit.

An object with no opportunity returns 0, and the caller must treat 0 as "no evidence" rather than floor it: a short
region cannot measure a long component, and a floored division turns no data into a confident wrong answer.

Gate: ``tests/calibration/test_effective_length.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..native import side_shares
from .region_arrays import boundary_reaches, boundary_region_indices, observable_pieces

if TYPE_CHECKING:
    from .region_arrays import RegionArrays


__all__ = [
    "UNBOUNDED_REACH",
    "conserved_cut_shares",
    "contained_eff_length",
    "crossing_eff_length",
    "observable_contained_eff_length",
    "observable_crossing_eff_length",
    "observable_cut_shares",
]

#: The reach of a molecule whose template does not end: gDNA always, and RNA wherever no taper applies. It exceeds
#: any fragment length by orders of magnitude and is finite, so the arithmetic below needs no case for infinity.
UNBOUNDED_REACH = 1.0e12


def _as_pmf(fl_pmf: np.ndarray) -> np.ndarray:
    """The fragment-length pmf as a normalised 1-D float64 array indexed by length, holding no mass at ``w = 0``. A
    zero-length molecule places nowhere, so its mass is dropped, not shared out over the lengths that place: what a
    transcript's effective length reads too (``FragmentLengthModel.compute_all_transcript_eff_lens``)."""
    p = np.asarray(fl_pmf, dtype=np.float64)
    if p.ndim != 1 or p.shape[0] == 0:
        raise ValueError("fl_pmf must be a non-empty 1-D array indexed by fragment length.")
    total = float(p.sum())
    p = p / total if total > 0.0 else p.copy()
    p[0] = 0.0
    return p


def contained_eff_length(region_len_bp: np.ndarray, fl_pmf: np.ndarray) -> np.ndarray:
    """``E_f[(region_len − w + 1)+]`` per region: the number of starts that place the whole molecule inside.

    The ``+1`` counts start positions: ``[s, s+w)`` sits inside ``[a, a+L)`` iff ``s`` is in ``[a, a+L−w]``.
    Computed as ``(L+1)·F(L) − S(L)`` with ``F`` the pmf's CDF and ``S(L) = Σ_{w≤L} w f(w)``; beyond the support the
    full sums apply, giving ``L + 1 − mean``.
    """
    p = _as_pmf(fl_pmf)
    n = p.shape[0]
    lengths = np.arange(n, dtype=np.float64)
    cdf = np.cumsum(p)  # F(w)
    cum_len = np.cumsum(lengths * p)  # S(w)

    region = np.asarray(region_len_bp, dtype=np.float64)
    idx = np.clip(np.floor(region).astype(np.int64), 0, n - 1)
    return np.maximum((region + 1.0) * cdf[idx] - cum_len[idx], 0.0)


def crossing_eff_length(
    fl_pmf: np.ndarray, reach_lo: np.ndarray, reach_hi: np.ndarray
) -> np.ndarray:
    """``E_f[max(0, min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1))]`` per object, for both boundary kinds.

    A molecule of length ``w`` crossing the boundary has ``a`` bases left of it and ``w − a`` right, with
    ``1 ≤ a ≤ w − 1``, and must fit its own template on each side: ``a ≤ R_lo`` and ``w − a ≤ R_hi``. Counting the
    admissible ``a`` gives the four-way ``min``: ``w − 1`` (both sides need a base), ``R_lo`` / ``R_hi`` (one side
    runs out of template) and ``R_lo + R_hi − w + 1`` (the molecule outgrows both remainders together).

    Computed in closed form from the pmf's cumulative sums in ``O(objects)``, with no ``(objects × lengths)`` matrix;
    :data:`UNBOUNDED_REACH` on both sides gives ``mean − 1`` exactly, but a finite reach must keep its taper: where the
    remaining template is shorter than the mean fragment the divisor is far below ``mean − 1``. ``reach_lo`` and
    ``reach_hi`` broadcast against each other; the result has their broadcast shape.
    """
    p = _as_pmf(fl_pmf)
    n = p.shape[0]
    lengths = np.arange(n, dtype=np.float64)
    cdf = np.cumsum(p)  # F(x) = Σ_{w ≤ x} f(w)
    cum_len = np.cumsum(lengths * p)  # S(x) = Σ_{w ≤ x} w f(w)
    lo = np.asarray(reach_lo, dtype=np.float64)
    hi = np.asarray(reach_hi, dtype=np.float64)
    lo, hi = np.broadcast_arrays(lo, hi)
    # The four-way min is piecewise linear in w, with the two reaches ordered a ≤ b:
    #   w ≤ a + 1:            w − 1              (both sides have room)
    #   a + 1 < w ≤ b + 1:    a                  (the short side alone binds)
    #   b + 1 < w ≤ a + b + 1: a + b + 1 − w     (the molecule is longer than both remainders together)
    #   beyond:               0
    # so its expectation is three sums over the pmf, each read off the two cumulative sums at the segment's end.
    a = np.minimum(lo, hi)
    b = np.maximum(lo, hi)

    def at(x):
        """``(F, S)`` at the largest length ``≤ x`` inside the support; ``(0, 0)`` below it."""
        idx = np.floor(np.minimum(x, float(n - 1))).astype(np.int64)
        below = idx < 0
        idx = np.maximum(idx, 0)
        return np.where(below, 0.0, cdf[idx]), np.where(below, 0.0, cum_len[idx])

    f1, s1 = at(a + 1.0)
    f2, s2 = at(b + 1.0)
    f3, s3 = at(a + b + 1.0)
    # the first sum starts at w = 1, where the pmf's mass starts (``_as_pmf``)
    out = (s1 - f1) + a * (f2 - f1) + (a + b + 1.0) * (f3 - f2) - (s3 - s2)
    return np.maximum(out, 0.0).reshape(lo.shape)


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  The conserved frame: where a template's placements put their unit
#
#  The accumulator splits a crossing fragment's unit over the objects it crosses: the boundaries it crosses (and,
#  spliced, the junctions at its blocks' ends) cut its path into slices, and each slice's share (its bases over the
#  fragment's length) is split equally between the objects bounding it (``Accumulator.deposit`` in
#  `tests/native/_accumulator_reference.py`). Summed over a template's placements, the pieces' contained shares
#  (:func:`contained_eff_length`) and the cuts' shares (:func:`conserved_cut_shares`) add up to the template's
#  fl-marginal length exactly. The capture-contracted length prices each share at its object's capture efficiency.
# ═══════════════════════════════════════════════════════════════════════════════════════════════════


def _side_share(
    p: np.ndarray, a: np.ndarray, near: np.ndarray, far: np.ndarray, n_threads: int = 1
) -> np.ndarray:
    """``Σ_w f(w)/w · Σ_x L_a(x)`` — the share a cut takes from the bases of the piece of length ``a`` on
    one side of it, with ``near`` bases of template on that side (``near ≥ a``) and ``far`` on the other.

    ``x`` is the number of the placement's bases on this side, ``max(1, w − far) ≤ x ≤ min(w − 1, near)``;
    ``L_a(x) = x`` while the placement starts inside the piece (the slice is bounded by this cut alone)
    and ``a/2`` once it runs through it (the slice is bounded by this cut and the piece's other cut). With
    no far side the inner sum is ``(w − 1)w/2`` up to ``w = a + 1``, ``a·w/2`` up to ``w = near + 1`` and
    ``a(near + 1)/2`` beyond; the far side removes the ``x < w − far``: ``(w − far − 1)(w − far)/2`` up to
    ``w = far + 1 + a`` and ``a(w − far)/2`` beyond; and no placement is longer than the template
    (``w ≤ near + far``). Each segment is a difference of three cumulative sums over the pmf — ``f``,
    ``w·f`` and ``f/w`` — so the whole is ``O(objects)`` with no per-length loop.
    """
    n = p.shape[0]
    w = np.arange(n, dtype=np.float64)
    live = np.where(w >= 2.0, p, 0.0)  # a crossing places a base on each side
    cum0 = np.concatenate([[0.0], np.cumsum(live)])
    cum1 = np.concatenate([[0.0], np.cumsum(live * w)])
    cumr = np.concatenate([[0.0], np.cumsum(live / np.maximum(w, 1.0))])

    # per cut in native code (``native.side_shares``)
    shape = np.broadcast(a, near, far).shape
    a, near, far = (
        np.ascontiguousarray(np.broadcast_to(x, shape), dtype=np.float64).ravel()
        for x in (a, near, far)
    )
    return side_shares(cum0, cum1, cumr, a, near, far, n_threads=int(n_threads)).reshape(shape)


def conserved_cut_shares(
    fl_pmf: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    reach_lo: np.ndarray,
    reach_hi: np.ndarray,
    *,
    n_threads: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """``(left, right)`` — a template's conserved share at each cut between two of its consecutive pieces,
    split into the part its left piece's bases carry and the part its right piece's carry.

    A length-``w`` placement with ``x`` bases left of the cut gives the cut ``L_a(x)/w`` from its left slice
    and ``L_b(w − x)/w`` from its right, ``a`` and ``b`` the two pieces' lengths (:func:`_side_share`); the
    reaches are the template's bases left and right of the cut (``reach_lo ≥ a``, ``reach_hi ≥ b``;
    :data:`UNBOUNDED_REACH` for gDNA, whose template does not end). A junction is a cut like any other in a
    spliced template's own coordinates. Where neither reach binds, ``left = ½ E_f[min(a, w − 1)]`` and
    ``right = ½ E_f[min(b, w − 1)]``. With the pieces' contained shares, a template's shares total its fl-marginal
    length ``E_f[(L − w + 1)⁺]`` exactly, and each object's share is what the deposit rule gives it. ``n_threads`` is
    the thread budget (0: every core).
    """
    p = _as_pmf(fl_pmf)
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    lo = np.asarray(reach_lo, dtype=np.float64)
    hi = np.asarray(reach_hi, dtype=np.float64)
    a, b, lo, hi = np.broadcast_arrays(a, b, lo, hi)
    if np.any(lo < a) or np.any(hi < b):
        raise ValueError("a cut's reach on each side must include the piece beside it.")
    return _side_share(p, a, lo, hi, n_threads), _side_share(p, b, hi, lo, n_threads)


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  Observable opportunity: the same three readings, where a fragment touching a masked base is no place
# ═══════════════════════════════════════════════════════════════════════════════════════════════════


def _reach(bases: np.ndarray) -> np.ndarray:
    """A side's reach in this module's terms: the bases to the nearest masked base, ``UNBOUNDED_REACH`` where none."""
    return np.where(np.isfinite(bases), bases, UNBOUNDED_REACH)


def observable_contained_eff_length(
    region_arrays: "RegionArrays", fl_pmf: np.ndarray
) -> np.ndarray:
    """``float64[R]`` — each region's contained support over its observable pieces, ``Σ_pieces E_f[(ℓₚ − w + 1)₊]``: a
    fragment contained in the region and touching no masked base lies in one piece. A region with nothing masked is one
    piece, so its support is :func:`contained_eff_length`'s of its length, to the bit."""
    piece_region, piece_len = observable_pieces(region_arrays)
    return np.bincount(
        piece_region,
        weights=contained_eff_length(piece_len.astype(np.float64), fl_pmf),
        minlength=region_arrays.n_regions,
    )


def observable_crossing_eff_length(region_arrays: "RegionArrays", fl_pmf: np.ndarray) -> np.ndarray:
    """``float64[E]`` — gDNA's crossing support at every contiguous boundary, its template the observable stretch the
    boundary lies in: :func:`crossing_eff_length` at the bases to the nearest masked base on each side, and at
    ``UNBOUNDED_REACH`` on a side with none, which with nothing masked is ``mean − 1`` everywhere, to the bit."""
    reach_lo, reach_hi = boundary_reaches(region_arrays)
    return crossing_eff_length(fl_pmf, _reach(reach_lo), _reach(reach_hi))


def observable_cut_shares(region_arrays: "RegionArrays", fl_pmf: np.ndarray) -> np.ndarray:
    """``float64[E]`` — gDNA's conserved share at every contiguous boundary over its observable stretch: the pieces
    beside the cut are the boundary's two flanking regions up to the nearest masked base, and the reaches the bases to
    it (``UNBOUNDED_REACH`` on a side with none). With nothing masked it is :func:`conserved_cut_shares` of the two
    flanking regions at unbounded reach, to the bit."""
    lo, hi = boundary_region_indices(region_arrays.ref_id)
    reach_lo, reach_hi = boundary_reaches(region_arrays)
    size = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    left, right = conserved_cut_shares(
        fl_pmf,
        np.minimum(size[lo], reach_lo),
        np.minimum(size[hi], reach_hi),
        _reach(reach_lo),
        _reach(reach_hi),
    )
    return left + right
