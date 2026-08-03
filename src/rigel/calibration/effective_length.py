"""rigel.calibration.effective_length — the divisor that turns a count into a start density.


    Gate: ``tests/calibration/test_effective_length.py`` — every formula enumerated, not restated

An **effective length** is the expected number of admissible fragment START POSITIONS at an object. It is
pure geometry against a fragment-length pmf: identical for any species, and applied **per component**,
because gDNA and RNA have different length distributions *and* different templates.

TWO FRAMES, ONE FAMILY. With ``w`` the molecule length and ``f`` its pmf::

    contained   E_f[ (node_len − w + 1)+ ]                                   fits wholly inside a node
    crossing    E_f[ max(0, min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)) ]     spans a 0-bp line

⭐ **The crossing formula covers BOTH edge kinds and BOTH components**, with ``R_lo``/``R_hi`` the
molecule's own remaining template either side of the line. Mean fragment length is its large-reach limit,
not a separate case: gDNA's template is the chromosome, so ``taper_g = 1``, its reaches are
:data:`UNBOUNDED_REACH` and the divisor collapses to ``mu − 1``. RNA's template ends where its transcript
ends, so its reaches come from the annotation.

⚠ **REACH IS PER COMPONENT, NOT PER EDGE, AND THAT IS WHAT MAKES IT AWKWARD.** An *unspliced* crossing is
a gDNA/RNA MIXTURE: the RNA part is bounded by its transcript, the gDNA part is not. So an edge does not
have "a" reach — each component has its own, and only the RNA one is finite. A junction edge is the easy
case, since only a spliced molecule uses it. ⛔ Production has ignored reach entirely up to now; whether
to keep ignoring it is an open decision, and this module is written so that
"ignore it" is expressible as ``UNBOUNDED_REACH`` rather than as a second code path.

⚠ **The `+1` in the contained formula is the discrete count of start positions, not a fudge.** A fragment
``[s, s+w)`` sits inside ``[a, a+L)`` iff ``s ∈ [a, a+L−w]``, which is ``L − w + 1`` positions. Dropping it
makes the divisor exactly 0 when a node is one fragment long — a division by zero that was floored to an
epsilon and produced densities of ~1e9 on 12.4 % of fine-partition nodes.

⚠ **An object with no opportunity must return 0, and the caller must treat 0 as "no evidence" rather than
flooring it.** A short node genuinely cannot measure a long component; that is physics
and a floored division turns "no data" into a confident wrong answer.

⛔ **The three mass-era divisors are DELETED** — ``boundary_side_eff_length`` (``E[min(l,R)]/2``),
``spliced_side_eff_length`` (``E[min²/2l]``) and ``boundary_side_crossing_count_eff_length``. They divided
a per-FACE mass, and a contiguous edge no longer has faces: it is a 0-bp line carrying one set of numbers.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "UNBOUNDED_REACH",
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
    """``E_f[w]`` — the mean fragment length. Reported as QC; the crossing limit derives it."""
    p = _as_pmf(fl_pmf)
    return float(np.dot(np.arange(p.shape[0], dtype=np.float64), p))


def contained_eff_length(node_len_bp: np.ndarray, fl_pmf: np.ndarray) -> np.ndarray:
    """``E_f[(node_len − w + 1)+]`` per node — the count of starts placing the whole molecule inside.

    Computed as ``(L+1)·F(L) − S(L)`` with ``F`` the pmf's CDF and ``S(L) = Σ_{w≤L} w f(w)``; beyond the
    support the full sums apply, giving ``L + 1 − mean``.

    ⚠ This is the frame that still NEEDS a length model. Containment probability differs 6.6× between
    gDNA and RNA at a 150 bp node, so a 10 % error in the fitted pmf costs 0.010–0.026 of composition
    the length models are load-bearing here, not hygiene.
    """
    p = _as_pmf(fl_pmf)
    n = p.shape[0]
    lengths = np.arange(n, dtype=np.float64)
    cdf = np.cumsum(p)  # F(w)
    cum_len = np.cumsum(lengths * p)  # S(w)

    node = np.asarray(node_len_bp, dtype=np.float64)
    idx = np.clip(np.floor(node).astype(np.int64), 0, n - 1)
    return np.maximum((node + 1.0) * cdf[idx] - cum_len[idx], 0.0)


def crossing_eff_length(
    fl_pmf: np.ndarray, reach_lo: np.ndarray, reach_hi: np.ndarray
) -> np.ndarray:
    """``E_f[max(0, min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1))]`` per object — one formula, both edge kinds.

    A molecule of length ``w`` crossing the line lies ``a`` bases to its left and ``w − a`` to its right,
    with ``1 ≤ a ≤ w − 1``; it must fit in what remains of its own template on each side, so
    ``a ≤ R_lo`` and ``w − a ≤ R_hi``. Counting the admissible ``a`` gives the four-way ``min`` — each
    term is one of the four binding constraints, and none is droppable:

    ==============================  ==================================================================
    ``w − 1``                       both sides need at least one base; the unbounded-reach limit
    ``R_lo`` / ``R_hi``             one side alone runs out of template
    ``R_lo + R_hi − w + 1``         the molecule is longer than BOTH remainders together
    ==============================  ==================================================================

    ⭐ Pass :data:`UNBOUNDED_REACH` on both sides for gDNA and the result is ``mean − 1`` exactly. The
    taper is not a refinement: at ``R = 100`` on RNA N(200,50) the divisor is 19.8 against an untapered
    199, so using the mean blindly under-reads the density **tenfold**.

    ``reach_lo`` and ``reach_hi`` broadcast against each other; the result has their broadcast shape.
    """
    p = _as_pmf(fl_pmf)
    lengths = np.arange(p.shape[0], dtype=np.float64)

    lo = np.asarray(reach_lo, dtype=np.float64)
    hi = np.asarray(reach_hi, dtype=np.float64)
    lo, hi = np.broadcast_arrays(lo, hi)
    # objects on the rows, fragment lengths on the columns
    lo_col, hi_col = lo.reshape(-1, 1), hi.reshape(-1, 1)
    placements = np.minimum(
        np.minimum(lengths - 1.0, np.minimum(lo_col, hi_col)),
        lo_col + hi_col - lengths + 1.0,
    )
    np.maximum(placements, 0.0, out=placements)
    return (placements @ p).reshape(lo.shape)
