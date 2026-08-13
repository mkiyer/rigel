"""rigel.calibration.effective_length — the divisor that turns a count into a start density.


    Gate: ``tests/calibration/test_effective_length.py`` — every formula enumerated, not restated

An **effective length** is the expected number of admissible fragment START POSITIONS at an object. It is
pure geometry against a fragment-length pmf: identical for any species, and applied **per component**,
because gDNA and RNA have different length distributions *and* different templates.

TWO FRAMES, ONE FAMILY. With ``w`` the molecule length and ``f`` its pmf::

    contained   E_f[ (region_len − w + 1)+ ]                                   fits wholly inside a region
    crossing    E_f[ max(0, min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)) ]     spans a 0-bp line

⭐ **The crossing formula covers BOTH boundary kinds and BOTH components**, with ``R_lo``/``R_hi`` the
molecule's own remaining template either side of the line. Mean fragment length is its large-reach limit,
not a separate case: gDNA's template is the chromosome, so ``taper_g = 1``, its reaches are
:data:`UNBOUNDED_REACH` and the divisor collapses to ``mu − 1``. RNA's template ends where its transcript
ends, so its reaches come from the annotation.

⚠ **REACH IS PER COMPONENT, NOT PER BOUNDARY, AND THAT IS WHAT MAKES IT AWKWARD.** An *unspliced* crossing is
a gDNA/RNA MIXTURE: the RNA part is bounded by its transcript, the gDNA part is not. So an boundary does not
have "a" reach — each component has its own, and only the RNA one is finite. A junction boundary is the easy
case, since only a spliced molecule uses it. ⛔ Production has ignored reach entirely up to now; whether
to keep ignoring it is an open decision, and this module is written so that
"ignore it" is expressible as ``UNBOUNDED_REACH`` rather than as a second code path.

⚠ **The `+1` in the contained formula is the discrete count of start positions, not a fudge.** A fragment
``[s, s+w)`` sits inside ``[a, a+L)`` iff ``s ∈ [a, a+L−w]``, which is ``L − w + 1`` positions. Dropping it
makes the divisor exactly 0 when a region is one fragment long — a division by zero that was floored to an
epsilon and produced densities of ~1e9 on 12.4 % of fine-partition regions.

⚠ **An object with no opportunity must return 0, and the caller must treat 0 as "no evidence" rather than
flooring it.** A short region genuinely cannot measure a long component; that is physics
and a floored division turns "no data" into a confident wrong answer.

⛔ **The three mass-era divisors are DELETED** — ``boundary_side_eff_length`` (``E[min(l,R)]/2``),
``spliced_side_eff_length`` (``E[min²/2l]``) and ``boundary_side_crossing_count_eff_length``. They divided
a per-FACE mass, and a contiguous boundary no longer has faces: it is a 0-bp line carrying one set of numbers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .region_chain import BOUNDARY, REGION, RegionChain

__all__ = [
    "UNBOUNDED_REACH",
    "LandedMoments",
    "build_slot_moments",
    "contained_eff_length",
    "contained_moments",
    "crossing_eff_length",
    "crossing_moments",
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


def contained_eff_length(region_len_bp: np.ndarray, fl_pmf: np.ndarray) -> np.ndarray:
    """``E_f[(region_len − w + 1)+]`` per region — the count of starts placing the whole molecule inside.

    Computed as ``(L+1)·F(L) − S(L)`` with ``F`` the pmf's CDF and ``S(L) = Σ_{w≤L} w f(w)``; beyond the
    support the full sums apply, giving ``L + 1 − mean``.

    ⚠ This is the frame that still NEEDS a length model. Containment probability differs 6.6× between
    gDNA and RNA at a 150 bp region, so a 10 % error in the fitted pmf costs 0.010–0.026 of composition
    the length models are load-bearing here, not hygiene.
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


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  THE OPPORTUNITY-TILTED LENGTH MOMENTS
#
#  ⭐ A fragment that LANDED at an object is not a draw from the library pmf — it is a draw from the
#  OPPORTUNITY-TILTED one, ``g(w) = f(w)·A(w) / E_f[A]``, because a length the object had more room for
#  is over-represented among the fragments it caught. These are that tilted pmf's moments, in the same
#  two frames :func:`contained_eff_length` and :func:`crossing_eff_length` use, and ``eff`` IS those two
#  functions' output — which is what lets a consumer assert the two against each other rather than
#  maintain two implementations of one quantity.
#
#  ⚠ They lived in a layer-5 composition channel until 2026-08-10 and were mis-filed there: nothing about
#  a tilted moment is a composition claim. The channel was measured and deleted; the geometry is not the
#  channel and belongs here, beside the opportunities it is a functional of.
# ═══════════════════════════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class LandedMoments:
    """The five moments of the opportunity-tilted length distribution ``g_c`` at each object.

    ``m1 = E[u]``, ``m2 = E[w]``, ``q1 = E[u²]``, ``q2 = E[w²]``, ``q12 = E[u·w]`` — everything the
    conditional mean and covariance of ``(Σu, Σw)`` need, and nothing else. ``eff`` is the tilt's own
    normaliser ``E_c[A]``, carried so a consumer can assert it against the divisor the solver used
    (they are the same quantity, and two implementations of one quantity is
    trap 27).

    All arrays are per object, or scalars broadcastable over objects.
    """

    m1: np.ndarray
    m2: np.ndarray
    q1: np.ndarray
    q2: np.ndarray
    q12: np.ndarray
    eff: np.ndarray

def _pmf_cumulants(fl_pmf: np.ndarray):
    """Cumulative sums of ``f(w)·w^k`` for ``k ∈ {−2,−1,0,1,2,3}`` — the whole of the region frame.

    ⚠ ``w = 0`` contributes 0 to every reciprocal sum: a zero-length fragment does not exist, and the
    pmf is 0 there in every real model. Guarding it here rather than trusting the input is what keeps a
    stray ``f(0) > 0`` from producing an infinity three call frames away.
    """
    p = np.asarray(fl_pmf, dtype=np.float64)
    total = float(p.sum())
    p = p / total if total > 0.0 else p
    w = np.arange(p.shape[0], dtype=np.float64)
    inv = np.zeros_like(w)
    np.divide(1.0, w, out=inv, where=w > 0.0)
    return (
        np.cumsum(p),  # F   = Σ f
        np.cumsum(p * inv),  # C1  = Σ f/w
        np.cumsum(p * inv * inv),  # C2  = Σ f/w²
        np.cumsum(p * w),  # S1  = Σ w f
        np.cumsum(p * w * w),  # S2  = Σ w² f
        np.cumsum(p * w * w * w),  # S3  = Σ w³ f
    )

def contained_moments(region_len_bp: np.ndarray, fl_pmf: np.ndarray) -> LandedMoments:
    """Moments of the tilted pmf for the CONTAINED population: ``A(w) = (ell − w + 1)+``, ``u(w) = 1/w``.

    ⭐ **Closed form, O(n_regions).** Each raw moment is ``(ell+1)·<cumsum> − <next cumsum>``, exactly the
    shape `effective_length.contained_eff_length` uses for the denominator — so no
    ``n_regions × max_len`` array is ever materialised (that would be 8 GB at human scale).

        E[A]     = (ell+1)·F  − S1        <- IS contained_eff_length
        E[A·u]   = (ell+1)·C1 − F
        E[A·w]   = (ell+1)·S1 − S2
        E[A·u²]  = (ell+1)·C2 − C1
        E[A·w²]  = (ell+1)·S2 − S3
        E[A·u·w] = (ell+1)·F  − S1        <- ⭐ identical to E[A], because u(w)·w = 1 at a region

    ⭐ That last identity means ``q12 ≡ 1`` at every region, for both components. It is not a shortcut: it
    says the two channels' cross-moment carries no composition information in the contained frame, and
    it falls out of the deposit rule rather than being imposed.
    """
    F, C1, C2, S1, S2, S3 = _pmf_cumulants(fl_pmf)
    n = F.shape[0]
    region = np.asarray(region_len_bp, dtype=np.float64)
    i = np.clip(np.floor(region).astype(np.int64), 0, n - 1)
    a = region + 1.0

    eff = np.maximum(a * F[i] - S1[i], 0.0)
    return _normalise(
        eff,
        a * C1[i] - F[i],
        a * S1[i] - S2[i],
        a * C2[i] - C1[i],
        a * S2[i] - S3[i],
        eff,  # E[A·u·w] == E[A] because u(w)·w == 1
    )

def crossing_moments(fl_pmf: np.ndarray) -> LandedMoments:
    """Moments for the CROSSING population at UNBOUNDED reach: ``A(w) = (w−1)+``, ``u(w) = 1/(w−1)``.

    Every entry is a scalar — under unbounded reach a line's opportunity does not depend on where it is,
    which is the "every boundary has the same expectation" property stated in moments.

        E[A]     = mu − 1                    <- IS crossing_eff_length at UNBOUNDED_REACH
        E[A·u]   = P(w >= 2)
        E[A·w]   = E[w²] − mu
        E[A·u²]  = Σ f(w)/(w−1)
        E[A·w²]  = E[w³] − E[w²]
        E[A·u·w] = mu                        (u(w)·w = w/(w−1), so Σ f(w)(w−1)·w/(w−1) = mu)

    ⚠ **Unbounded reach only, matching `build_region_geometry`'s default.** With the TRAPS: prove-the-substrate taper switched on
    (``boundary_rna_reach``) the opportunity becomes per-boundary and these moments would have to as well. The
    taper was measured as a null (≤ 0.0002), so the default path is the one wired;
    a consumer that turns the taper on must extend this function rather than silently mismatch.
    """
    p = np.asarray(fl_pmf, dtype=np.float64)
    total = float(p.sum())
    p = p / total if total > 0.0 else p
    w = np.arange(p.shape[0], dtype=np.float64)
    ok = w >= 2.0  # a length-0 or length-1 fragment cannot cross a 0-bp line
    inv = np.zeros_like(w)
    np.divide(1.0, w - 1.0, out=inv, where=ok)

    eff = float((p * np.maximum(w - 1.0, 0.0)).sum())
    return _normalise(
        np.asarray(eff),
        np.asarray(float(p[ok].sum())),
        np.asarray(float((p * np.maximum(w - 1.0, 0.0) * w).sum())),
        np.asarray(float((p * inv)[ok].sum())),
        np.asarray(float((p * np.maximum(w - 1.0, 0.0) * w * w).sum())),
        np.asarray(float((p * w)[ok].sum())),
    )

def _normalise(eff, e_u, e_w, e_uu, e_ww, e_uw) -> LandedMoments:
    """Divide the raw ``E[A··]`` moments by ``E[A]`` to get the tilted-pmf moments.

    ⛔ Zero opportunity ⇒ every moment is 0, never a floored division.
    A slot with no opportunity for a component contributes nothing, and the caller's ``det > 0`` gate
    then makes the whole term inert there.
    """
    eff = np.asarray(eff, dtype=np.float64)
    live = eff > 0.0

    def d(x):
        x = np.broadcast_to(np.asarray(x, dtype=np.float64), eff.shape)
        return np.divide(x, eff, out=np.zeros(eff.shape, dtype=np.float64), where=live)

    return LandedMoments(m1=d(e_u), m2=d(e_w), q1=d(e_uu), q2=d(e_ww), q12=d(e_uw), eff=eff)

def build_slot_moments(chain: RegionChain, region_arrays, fl_pmf: np.ndarray) -> LandedMoments:
    """Scatter the two frames' moments onto the chain: contained at REGION slots, crossing at BOUNDARY slots.

    The same slot layout `build_region_geometry` uses for ``eff_gdna``/``eff_rna``, so ``moments.eff`` is
    that array — asserted by ``test_eff_matches_the_solver_divisor``, because two implementations of one
    quantity is how a ½ went unnoticed for months.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_region, is_boundary = kind == REGION, kind == BOUNDARY
    n = int(chain.n_slots)

    region_len = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    region_m = contained_moments(region_len, fl_pmf) if region_len.shape[0] else None
    boundary_m = crossing_moments(fl_pmf)

    fields = {}
    for name in ("m1", "m2", "q1", "q2", "q12", "eff"):
        out = np.zeros(n, dtype=np.float64)
        if region_m is not None:
            out[is_region] = getattr(region_m, name)[obj[is_region]]
        out[is_boundary] = float(getattr(boundary_m, name))
        fields[name] = out
    return LandedMoments(**fields)
