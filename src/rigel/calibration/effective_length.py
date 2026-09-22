"""The divisor that turns a count into a start density.

Gate: ``tests/calibration/test_effective_length.py``, which enumerates every formula rather than
restating it.

An effective length is the expected number of admissible fragment START POSITIONS at an object. It
is pure geometry against a fragment-length pmf: identical for any species, and applied per
component, because gDNA and RNA have different length distributions and different templates.

Two frames, one family. With ``w`` the molecule length and ``f`` its pmf::

    contained   E_f[ (region_len - w + 1)+ ]                              fits wholly inside a region
    crossing    E_f[ max(0, min(w-1, R_lo, R_hi, R_lo + R_hi - w + 1)) ]  spans a 0-bp boundary

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

from dataclasses import dataclass

import numpy as np


__all__ = [
    "UNBOUNDED_REACH",
    "BaseTaper",
    "base_taper",
    "contained_eff_length",
    "crossing_base_shares",
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
#  THE PER-BASE FRAME — the same placement counts, read base by base
#
#  A fragment's capture efficiency is the mean per-base efficiency over the bases it covers, so a
#  template's effective length under capture is a sum over its BASES, each weighted by how many starts
#  cover it (`capture_eff_length`, `priors`; EQUATIONS.md §11). Two functionals of the pmf serve that:
#  the end taper τ(x) on a template's own bases, and the share of a crossing fragment's bases that lie
#  in each piece beside a boundary. Both are geometry against the pmf and nothing else, which is why they
#  live here beside the divisors they are built from.
# ═══════════════════════════════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class BaseTaper:
    """``τ(x) = E_f[ (starts whose fragment covers base x) / w ]`` on a template of length ``L``, and its
    interval sums — the weight a base carries in the template's effective length.

    A fragment of length ``w`` starting at ``s`` covers ``x`` iff ``s ≤ x < s + w`` with ``0 ≤ s ≤ L − w``,
    which is ``min(x + 1, w, L − x, L − w + 1)⁺`` starts; each of those starts spreads its unit over its
    ``w`` bases, so the base carries ``1/w`` of it. Summed over the template's bases, ``Σ_x τ(x)`` is the
    fl-marginal length ``Σ_w f(w)(L − w + 1)⁺`` exactly — the per-base frame partitions the start count,
    it does not re-derive it — and away from both ends ``τ = 1``: only the bases within a fragment of an
    end are tapered.

    For ``L ≥ 2·w_max − 1`` the four-way min collapses to ``min(d, w)`` with ``d = min(x + 1, L − x)``, so
    one cumulative table ``CT(d) = Σ_{d' ≤ d} T(d')``, ``T(d) = Σ_w f(w) min(d, w)/w``, serves every such
    template in O(1) per interval; a shorter template is evaluated base by base. A zero-length fragment
    covers no base and carries no weight: a smoothed length model can put mass at ``w = 0`` (a sparse real
    library's did), and it is dropped here.
    """

    pmf: np.ndarray
    cum: np.ndarray
    wmax: int
    #: the pmf's mass at positive lengths — the interior base's weight, 1 for a pmf with no mass at zero
    mass: float

    def _cum_long(self, d: np.ndarray) -> np.ndarray:
        d = np.asarray(d, dtype=np.int64)
        return np.where(
            d <= self.wmax,
            self.cum[np.minimum(d, self.wmax)],
            self.cum[self.wmax] + (d - self.wmax) * self.mass,
        )

    def _cum_short(self, L: int) -> np.ndarray:
        x = np.arange(L, dtype=np.float64)
        d = np.minimum(x + 1.0, L - x)
        w_all = np.arange(self.pmf.shape[0], dtype=np.float64)
        live = (self.pmf > 0.0) & (w_all >= 1.0)
        w = w_all[live]
        starts = np.maximum(
            np.minimum(np.minimum(d[:, None], w[None, :]), L - w[None, :] + 1.0), 0.0
        )
        tau = (self.pmf[live][None, :] * starts / w[None, :]).sum(1)
        return np.concatenate([[0.0], np.cumsum(tau)])

    def interval_sums(self, x0: np.ndarray, x1: np.ndarray, L: np.ndarray | int) -> np.ndarray:
        """``Σ_{x ∈ [x0, x1)} τ(x)`` for arrays of intervals, each on a template of length ``L`` (one
        length, or one per interval): the long templates in one vectorised pass, the short ones grouped
        by length so each table is built once."""
        x0 = np.asarray(x0, dtype=np.int64)
        x1 = np.asarray(x1, dtype=np.int64)
        L = np.broadcast_to(np.asarray(L, dtype=np.int64), x0.shape)
        out = np.empty(x0.shape, dtype=np.float64)
        is_long = L >= 2 * self.wmax - 1
        if is_long.any():
            Ll, a, b = L[is_long], x0[is_long], x1[is_long]
            half = Ll // 2
            total = self._cum_long(half) + self._cum_long(Ll - half)

            def F(x):  # Σ_{x' < x} τ(x'): the left taper up to the middle, the right one mirrored
                return np.where(x <= half, self._cum_long(x), total - self._cum_long(Ll - x))

            out[is_long] = F(b) - F(a)
        short = np.flatnonzero(~is_long)
        short = short[np.argsort(L[short], kind="stable")]
        lengths, first = np.unique(L[short], return_index=True)
        for Ls, s0, s1 in zip(lengths, first, np.r_[first[1:], short.size]):
            sel = short[s0:s1]
            c = self._cum_short(int(Ls))
            out[sel] = c[x1[sel]] - c[x0[sel]]
        return out


def base_taper(fl_pmf: np.ndarray) -> BaseTaper:
    """The :class:`BaseTaper` of a pmf: its table ``CT(d)`` for ``d ≤ w_max``."""
    p = _as_pmf(fl_pmf)
    w_all = np.arange(p.shape[0], dtype=np.float64)
    live = (p > 0.0) & (w_all >= 1.0)
    if not live.any():
        raise ValueError("fl_pmf has no mass at a positive length.")
    wmax = int(np.flatnonzero(live).max())
    w = w_all[live]
    d = np.arange(1, wmax + 1, dtype=np.float64)
    T = (p[live][None, :] * np.minimum(d[:, None], w[None, :]) / w[None, :]).sum(1)
    return BaseTaper(
        pmf=p, cum=np.concatenate([[0.0], np.cumsum(T)]), wmax=wmax, mass=float(p[live].sum())
    )


def crossing_base_shares(
    region_arrays, fl_pmf: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(boundary, region, share)`` — the expected number of base-starts a crossing fragment places in
    each piece within its reach of a contiguous boundary, at UNBOUNDED reach (gDNA's template).

    A molecule of length ``w`` crossing the boundary lies ``a`` bases to its left and ``w − a`` to its
    right, ``1 ≤ a ≤ w − 1``, one start each; spreading each start over its ``w`` bases, the piece at
    cumulative distance ``(c_{j−1}, c_j]`` on a side receives ``G(c_j) − G(c_{j−1})`` base-starts with

        G(c) = Σ_w f(w) Σ_{a=1}^{w−1} min(a, c) / w,

    and the two sides are symmetric. The shares of one boundary sum to its crossing opportunity
    ``E_f[w − 1]`` exactly (:func:`crossing_eff_length` at ``UNBOUNDED_REACH``) as long as the pieces on
    each side reach a fragment's length; a chromosome's end within reach loses the part that has no
    template. A boundary's crossing count is therefore ``Σ_q share_q · ρ_q`` in expectation over the
    pieces ``q`` its fragments cover, which is what lets a piece too short to contain a fragment be read
    from its edges (`capture_efficiency`).
    """
    from .region_arrays import boundary_region_indices

    p = _as_pmf(fl_pmf)
    w = np.arange(p.shape[0], dtype=np.float64)
    live = (p > 0.0) & (w >= 2.0)
    wmax = int(np.flatnonzero(p > 0.0).max())
    c = np.arange(0, wmax + 1, dtype=np.float64)[:, None]
    wl = w[live][None, :]
    # Σ_{a=1}^{w−1} min(a, c): every a for c ≥ w − 1, else the ramp to c and c thereafter
    ramp = c * (c + 1.0) / 2.0 + c * (wl - 1.0 - c)
    G = (p[live][None, :] * np.where(c >= wl - 1.0, wl * (wl - 1.0) / 2.0, ramp) / wl).sum(1)

    starts = np.asarray(region_arrays.start, dtype=np.int64)
    ends = np.asarray(region_arrays.end, dtype=np.int64)
    ref_id = np.asarray(region_arrays.ref_id)
    length = ends - starts
    lo, hi = boundary_region_indices(ref_id)
    out_e: list[np.ndarray] = []
    out_q: list[np.ndarray] = []
    out_a: list[np.ndarray] = []
    # walk outward from each boundary on each side, one step per iteration for every boundary at once
    for first, step in ((lo, -1), (hi, 1)):
        e = np.arange(lo.size, dtype=np.int64)
        r = first.copy()
        cum = np.zeros(lo.size, dtype=np.int64)
        while e.size:
            nxt = np.minimum(cum + length[r], wmax)
            out_e.append(e)
            out_q.append(r)
            out_a.append(G[nxt] - G[cum])
            cum = nxt
            r = r + step
            keep = (cum < wmax) & (r >= 0) & (r < length.size)
            keep &= ref_id[np.clip(r, 0, length.size - 1)] == ref_id[first[e]]
            e, r, cum = e[keep], r[keep], cum[keep]
    if not out_e:
        z = np.zeros(0, dtype=np.int64)
        return z, z, np.zeros(0, dtype=np.float64)
    return np.concatenate(out_e), np.concatenate(out_q), np.concatenate(out_a)
