"""The divisor that turns a count into a start density.

Gate: ``tests/calibration/test_effective_length.py``, which enumerates every formula rather than
restating it.

An effective length is the expected number of admissible fragment START POSITIONS at an object. It
is pure geometry against a fragment-length pmf: identical for any species, and applied per
component, because gDNA and RNA have different length distributions and different templates.

Two frames, one family. With ``w`` the molecule length and ``f`` its pmf::

    contained   E_f[ (region_len - w + 1)+ ]                              fits wholly inside a region
    crossing    E_f[ max(0, min(w-1, R_lo, R_hi, R_lo + R_hi - w + 1)) ]  spans a 0-bp boundary

The CONSERVED SHARE (:func:`conserved_cut_shares`) is the third reading of the same placements: not how
many starts reach an object but how much of each placement's unit the deposit rule gives it, so that a
template's objects share its fl-marginal length with nothing counted twice — the frame the
capture-contracted length is priced in.

The crossing formula covers both boundary kinds and both components, with ``R_lo`` / ``R_hi`` the
molecule's own remaining template either side of the boundary. Mean fragment length is its
large-reach limit, not a separate case: gDNA's template is the chromosome, so its reaches are
:data:`UNBOUNDED_REACH` and the divisor collapses to ``mu - 1``. RNA's template ends where its
transcript ends, so its reaches come from the annotation.

Reach is per component, not per boundary, and that is what makes it awkward. An unspliced crossing
is a gDNA/RNA mixture: the RNA part is bounded by its transcript, the gDNA part is not, so a
boundary does not have "a" reach — each component has its own and only the RNA one is finite. An sj
boundary is the easy case, since only a spliced molecule uses it, and it is where the annotation's
exonic reach is actually passed. At a contiguous boundary the unspliced RNA reach is
:data:`UNBOUNDED_REACH` by ruling, which this module expresses as a value rather than as a second
code path.

The ``+1`` in the contained formula is the discrete count of start positions, not a fudge: a
fragment ``[s, s+w)`` sits inside ``[a, a+L)`` iff ``s`` is in ``[a, a+L-w]``, which is
``L - w + 1`` positions. Dropping it makes the divisor exactly 0 when a region is one fragment
long, and a division by zero floored to an epsilon then produces astronomical densities on every
short region of a fine partition.

An object with no opportunity must return 0, and the caller must treat 0 as "no evidence" rather
than flooring it. A short region genuinely cannot measure a long component; that is physics, and a
floored division turns "no data" into a confident wrong answer.
"""

from __future__ import annotations

import numpy as np


__all__ = [
    "UNBOUNDED_REACH",
    "conserved_cut_shares",
    "contained_eff_length",
    "crossing_eff_length",
    "fl_mean",
]

#: The reach of a molecule whose template does not end — gDNA, always; and RNA wherever the taper is
#: deliberately not applied. Large enough to exceed any fragment length by orders of magnitude, finite so
#: that it survives the integer arithmetic below without special-casing infinity.
UNBOUNDED_REACH = 1.0e12


def _as_pmf(fl_pmf: np.ndarray) -> np.ndarray:
    """The fragment-length pmf as a normalised 1-D float64 array indexed by length."""
    p = np.asarray(fl_pmf, dtype=np.float64)
    if p.ndim != 1 or p.shape[0] == 0:
        raise ValueError("fl_pmf must be a non-empty 1-D array indexed by fragment length.")
    total = float(p.sum())
    return p / total if total > 0.0 else p


def fl_mean(fl_pmf: np.ndarray) -> float:
    """``E_f[w]`` — the mean fragment length, the limit the unbounded crossing length derives (the
    gates use it as that oracle)."""
    p = _as_pmf(fl_pmf)
    return float(np.dot(np.arange(p.shape[0], dtype=np.float64), p))


def contained_eff_length(region_len_bp: np.ndarray, fl_pmf: np.ndarray) -> np.ndarray:
    """``E_f[(region_len − w + 1)+]`` per region — the count of starts placing the whole molecule inside.

    Computed as ``(L+1)·F(L) − S(L)`` with ``F`` the pmf's CDF and ``S(L) = Σ_{w≤L} w f(w)``; beyond the
    support the full sums apply, giving ``L + 1 − mean``.

    This is the frame that needs a length model. Containment probability differs several-fold
    between gDNA and RNA at a short region, so an error in either fitted pmf moves the composition
    directly: the length models are load-bearing here, not hygiene.
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
    """``E_f[max(0, min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1))]`` per object — one formula, both boundary kinds.

    A molecule of length ``w`` crossing the boundary lies ``a`` bases to its left and ``w − a`` to its right,
    with ``1 ≤ a ≤ w − 1``; it must fit in what remains of its own template on each side, so
    ``a ≤ R_lo`` and ``w − a ≤ R_hi``. Counting the admissible ``a`` gives the four-way ``min`` — each
    term is one of the four binding constraints, and none is droppable:

    ==============================  ==================================================================
    ``w − 1``                       both sides need at least one base; the unbounded-reach limit
    ``R_lo`` / ``R_hi``             one side alone runs out of template
    ``R_lo + R_hi − w + 1``         the molecule is longer than BOTH remainders together
    ==============================  ==================================================================

    Computed in closed form from the pmf's cumulative sums, in ``O(objects)`` with no
    ``(objects × lengths)`` matrix — that matrix cost 9 GB of transient on the human sj axis
    (`tests/calibration/test_effective_length.py` holds the enumerated and the matrix brute forces).
    Pass :data:`UNBOUNDED_REACH` on both sides for gDNA and the result is ``mean - 1`` exactly. The
    taper is not a refinement: where the remaining template is shorter than the mean fragment, the
    tapered divisor is an order of magnitude below the untapered one, so using the mean blindly
    under-reads the density by the same factor.

    ``reach_lo`` and ``reach_hi`` broadcast against each other; the result has their broadcast shape.
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
    # so its expectation is three sums over the pmf, each read off the two cumulative sums at the
    # segment's end — O(objects), and no (objects × lengths) matrix. Beyond the support the full sums
    # apply, and UNBOUNDED_REACH on both sides gives mean − 1 exactly.
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
    # the first sum runs from one base: a zero-length molecule places nowhere, and the brute force
    # clamps it per length where the cumulative sums would credit it (0 − 1)·f(0)
    out = (s1 - f1 + p[0]) + a * (f2 - f1) + (a + b + 1.0) * (f3 - f2) - (s3 - s2)
    return np.maximum(out, 0.0).reshape(lo.shape)


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  THE CONSERVED FRAME — where a template's placements put their unit
#
#  The accumulator splits a crossing fragment's unit over the objects it crosses: its path is cut into
#  slices by the boundaries it crosses (and, spliced, by the junctions at its blocks' ends), and each
#  slice's share of the fragment — its bases over the fragment's length — is shared equally by the objects
#  bounding it (`tests/native/_accumulator_reference.py`, ``Accumulator.deposit``). Summed over a template's
#  placements that gives each object the template's CONSERVED SHARE of it: a piece's contained share
#  (:func:`contained_eff_length`) and a cut's share (:func:`conserved_cut_shares`), which add up to the
#  template's fl-marginal length exactly, however many cuts one fragment crosses. The capture-contracted
#  length prices each share at its object's capture efficiency (`capture_eff_length`, `priors`).
# ═══════════════════════════════════════════════════════════════════════════════════════════════════


def _side_share(p: np.ndarray, a: np.ndarray, near: np.ndarray, far: np.ndarray) -> np.ndarray:
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

    def over(cum, lo, hi):
        """``Σ`` of the table's summand over the integer lengths ``lo < w ≤ hi``."""
        i0 = np.clip(np.floor(lo).astype(np.int64) + 1, 0, n)
        i1 = np.clip(np.floor(hi).astype(np.int64) + 1, 0, n)
        return np.where(i1 > i0, cum[np.maximum(i1, i0)] - cum[i0], 0.0)

    top = near + far
    starts_inside = np.minimum(a + 1.0, top)
    runs_through = np.minimum(near + 1.0, top)
    kept = (
        0.5 * (over(cum1, -1.0, starts_inside) - over(cum0, -1.0, starts_inside))
        + 0.5 * a * over(cum0, starts_inside, runs_through)
        + 0.5 * a * (near + 1.0) * over(cumr, runs_through, top)
    )
    lo, hi = far + 1.0, np.minimum(far + 1.0 + a, top)
    removed = 0.5 * (
        over(cum1, lo, hi)
        - (2.0 * far + 1.0) * over(cum0, lo, hi)
        + far * (far + 1.0) * over(cumr, lo, hi)
    ) + 0.5 * a * (over(cum0, hi, top) - far * over(cumr, hi, top))
    return np.maximum(kept - removed, 0.0)


def conserved_cut_shares(
    fl_pmf: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    reach_lo: np.ndarray,
    reach_hi: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """``(left, right)`` — a template's conserved share at each cut between two of its consecutive pieces,
    split into the part its left piece's bases carry and the part its right piece's carry.

    A length-``w`` placement with ``x`` bases left of the cut gives the cut ``L_a(x)/w`` from its left slice
    and ``L_b(w − x)/w`` from its right, ``a`` and ``b`` the two pieces' lengths (:func:`_side_share`); the
    reaches are the template's bases left and right of the cut (``reach_lo ≥ a``, ``reach_hi ≥ b``;
    :data:`UNBOUNDED_REACH` for gDNA, whose template does not end). A junction is a cut like any other in a
    spliced template's own coordinates. Where neither reach binds the sums close to
    ``left = ½ E_f[min(a, w − 1)]`` and ``right = ½ E_f[min(b, w − 1)]``. With the pieces' contained shares
    a template's shares total its fl-marginal length ``E_f[(L − w + 1)⁺]`` exactly, and each object's share
    is what the deposit rule gives it (gate: ``tests/calibration/test_effective_length.py``, per object
    against the reference accumulator).
    """
    p = _as_pmf(fl_pmf)
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    lo = np.asarray(reach_lo, dtype=np.float64)
    hi = np.asarray(reach_hi, dtype=np.float64)
    a, b, lo, hi = np.broadcast_arrays(a, b, lo, hi)
    if np.any(lo < a) or np.any(hi < b):
        raise ValueError("a cut's reach on each side must include the piece beside it.")
    return _side_share(p, a, lo, hi), _side_share(p, b, hi, lo)
