"""The EM's two fragment-length tables: one shared law, unless the data show a difference capture cannot make.

The E-step scores an unspliced fragment's length by the ratio of the two origins' tables, ``f_R(l) / f_G(l)``. Under
hybrid capture a fragment's capture weight belongs to its placement, which both hypotheses share, so it cancels from
that ratio; capture enters only the normalisers (:mod:`rigel.calibration.capture_eff_length`). The two laws the tool
measures are read in different frames, though. gDNA's (``FLModels.gdna_pmf``) is read off target, uncaptured; RNA's
(``FLModels.rna_pmf``) is the spliced census, captured wherever RNA is. On equal chemistry they still differ by what
capture did to RNA's census, and a scorer handed the two as they stand reads that frame difference as chemistry.

So RNA's law is the table of both origins by default, and gDNA's table moves away from it only as far as the data show
a difference capture cannot have made, in two steps.

1. **The capture mixture.** gDNA's law with a captured fraction ``c`` size-biased, ``G_c = (1 - c) g + c g l / mu_g``,
   ``c`` the fraction that best explains RNA's census, between none and ``kappa``, RNA's captured fraction from
   calibration's own efficiencies (:func:`rigel.calibration.capture_eff_length.rna_captured_fraction`). Its
   log-likelihood is concave in ``c``, with slope ``mu_r / mu_g - 1`` at ``c = 0``: capture only lengthens, so where
   RNA's census is no longer than gDNA's law, ``c = 0``. This step shapes gDNA's table at full trust; it does not
   decide how much to trust it.
2. **The decision.** ``D``, the Hellinger distance between the two measured laws less what sampling noise adds,
   against ``B``, the most capture at ``kappa`` can make: ``w = clip(D / B - 1, 0, 1)``. Hellinger distance is a
   metric, so by the triangle inequality the true difference is at least ``D - B``; ``w`` is that guaranteed-real
   part as a fraction of the largest artifact, whole once it is as large. One law while ``D <= B``; the tilted law
   whole from ``D >= 2B``; continuous between, with no threshold. The table lies on the geometric path from RNA's
   law (``w = 0``) to ``G_c`` (``w = 1``), exact at both ends.

``B`` is a bound, not a guess. If a fragment's capture weight grows at most in proportion to its length — its overlap
with the probes does, and saturates — the captured census's likelihood ratio to the uncaptured law is dominated in
convex order by the size-biased one, so the size-biased mixture at ``kappa`` is the worst case for every f-divergence,
Hellinger's included (gate: ``tests/calibration/test_shared_length_law.py``). A response that rewards length faster
than proportionally lies outside that premise and can exceed it.

The decision is a ramp to full trust and never a blend weighted by evidence: a real gap partly shrunk shifts the
balance between the two origins as surely as a false gap kept. Off capture ``kappa = 0``, so ``B = 0``, ``c = 0`` and
any resolved difference is kept whole: each origin's own law. The sampling noise is the plug-in variance of each
law's own estimator, scaled by its empirical-Bayes weight, read at the null's common law.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "LengthDecision",
    "LengthTables",
    "capture_bound",
    "capture_mixture",
    "geometric_path",
    "hellinger_sq",
    "sampling_noise",
    "shared_length_law",
    "size_biased",
    "trust",
]


@dataclass(frozen=True, slots=True)
class LengthDecision:
    """What the shared law decided, and from what.

    ``kappa`` RNA's captured fraction; ``captured`` the capture mixture's ``c`` in ``[0, kappa]``; ``contrast`` ``D``
    and ``bound`` ``B``, Hellinger distances; ``noise`` the squared distance sampling noise alone adds; ``weight``
    ``w``, 0 for one law and 1 for gDNA's own law whole.
    """

    kappa: float
    captured: float
    contrast: float
    bound: float
    noise: float
    weight: float


@dataclass(frozen=True, slots=True)
class LengthTables:
    """The EM's two length tables, and the law RNA's normaliser carries capture from.

    ``rna`` is ``FLModels.rna_pmf`` itself; ``gdna`` the scorer's gDNA table; ``anchor`` gDNA's uncaptured law
    (``FLModels.gdna_pmf``), the law calibration's capture efficiencies are read on. The scorer and RNA's normaliser
    read the same three (:func:`rigel.calibration.capture_eff_length.transcript_capture_eff_lengths`).
    """

    rna: np.ndarray
    gdna: np.ndarray
    anchor: np.ndarray
    decision: LengthDecision


def hellinger_sq(a: np.ndarray, b: np.ndarray) -> float:
    """``Σ_l (√a_l − √b_l)²``: the squared Hellinger distance of two normalised laws, at most 2."""
    return float(((np.sqrt(a) - np.sqrt(b)) ** 2).sum())


def size_biased(g: np.ndarray, c: float) -> np.ndarray:
    """``(1 − c) g + c g l / mu_g``, normalised: ``g`` with a fraction ``c`` of its fragments size-biased."""
    L = np.arange(g.size, dtype=np.float64)
    out = (1.0 - c) * g + c * g * L / float(g @ L)
    return out / out.sum()


def capture_mixture(g: np.ndarray, r: np.ndarray, kappa: float) -> tuple[np.ndarray, float]:
    """``(G_c, c)``: gDNA's law brought into RNA's captured frame as far as capture can move it.

    ``c`` maximises ``Σ_l r(l) log(1 + c (l / mu_g − 1))`` over ``[0, kappa]``. The objective is concave (the mixture is
    normalised for every ``c``), so its slope falls through ``c`` and bisection to machine precision finds the root; at
    ``c = 0`` the slope is ``mu_r / mu_g − 1``, so ``c = 0`` and ``G_c`` is ``g`` itself wherever RNA's census is no
    longer than ``g``. ``kappa < 1`` always (it is a mean of ``1 − 1/c_region``), so the mixture stays positive.
    """
    L = np.arange(g.size, dtype=np.float64)
    mu = float(g @ L)
    if kappa <= 0.0 or mu <= 0.0:
        return g.copy(), 0.0
    both = (r > 0.0) & (g > 0.0)
    rise = (g * L / mu - g)[both]
    weight, base = r[both], g[both]

    def slope(c: float) -> float:
        return float((weight * rise / (base + c * rise)).sum())

    if slope(0.0) <= 0.0:
        return g.copy(), 0.0
    if slope(kappa) >= 0.0:
        return size_biased(g, kappa), float(kappa)
    lo, hi = 0.0, float(kappa)
    while True:
        mid = 0.5 * (lo + hi)
        if not lo < mid < hi:
            break
        if slope(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return size_biased(g, lo), lo


def capture_bound(g: np.ndarray, kappa: float) -> float:
    """``B²``: the squared Hellinger distance capture at captured fraction ``kappa`` can put between ``g`` and its
    captured census, the size-biased mixture being the worst case (module docstring). 0 with no capture."""
    if kappa <= 0.0:
        return 0.0
    return hellinger_sq(g, size_biased(g, kappa))


def sampling_noise(g: np.ndarray, r: np.ndarray, var_g: np.ndarray, var_r: np.ndarray) -> float:
    """The squared Hellinger distance two noisy estimates of ONE law would show: ``Σ_l (Var g_l + Var r_l) / 4 M_l``,
    ``M = (g + r) / 2``, the delta method at the null's common law. The reference is the null's, not each law's own: at
    a length where one law is all but empty and the other is not, its own delta method divides its noise by nearly
    nothing and would call any gap noise."""
    m = 0.5 * (g + r)
    on = m > 0.0
    return float(((var_g + var_r)[on] / (4.0 * m[on])).sum())


def trust(contrast_sq: float, bound_sq: float) -> float:
    """``w = clip(D / B − 1, 0, 1)``: 0 while the difference fits what capture can make, 1 once its guaranteed-real part
    ``D − B`` is as large as ``B``. No difference is ``w = 0``; a resolved difference with no capture is ``w = 1``."""
    if contrast_sq <= 0.0:
        return 0.0
    if bound_sq <= 0.0:
        return 1.0
    return float(np.clip(np.sqrt(contrast_sq / bound_sq) - 1.0, 0.0, 1.0))


def geometric_path(r: np.ndarray, G: np.ndarray, w: float) -> np.ndarray:
    """``normalize(r^(1−w) G^w)``: exactly ``r`` at ``w = 0`` and exactly ``G`` at ``w = 1``, no logarithm at either
    end. Between them, a length only one law holds keeps the smallest positive double, so no length either law admits
    becomes impossible."""
    if w <= 0.0:
        return r.copy()
    if w >= 1.0:
        return G.copy()
    both = (r > 0.0) & (G > 0.0)
    log_t = np.full(r.size, -np.inf)
    log_t[both] = (1.0 - w) * np.log(r[both]) + w * np.log(G[both])
    t = np.zeros(r.size)
    t[both] = np.exp(log_t[both] - log_t[both].max())
    t[((r > 0.0) | (G > 0.0)) & ~both] = np.finfo(np.float64).tiny
    return t / t.sum()


def shared_length_law(
    rna_pmf: np.ndarray,
    gdna_pmf: np.ndarray,
    rna_var: np.ndarray,
    gdna_var: np.ndarray,
    kappa: float,
) -> LengthTables:
    """The EM's two length tables from the two measured laws, their sampling variances and RNA's captured fraction.

    ``rna_pmf`` is the RNA table untouched. gDNA's is ``rna_pmf`` itself where the laws differ by no more than capture
    can make (``w = 0``), the capture mixture ``G_c`` once the difference is twice that (``w = 1``), and the geometric
    path between them. Both laws are normalised already (``FLModels``) and are read as given.
    """
    r = np.asarray(rna_pmf, dtype=np.float64)
    g = np.asarray(gdna_pmf, dtype=np.float64)
    tilted, captured = capture_mixture(g, r, kappa)
    noise = sampling_noise(
        g, r, np.asarray(gdna_var, dtype=np.float64), np.asarray(rna_var, dtype=np.float64)
    )
    contrast_sq = max(hellinger_sq(g, r) - noise, 0.0)
    bound_sq = capture_bound(g, kappa)
    weight = trust(contrast_sq, bound_sq)
    return LengthTables(
        rna=r,
        gdna=geometric_path(r, tilted, weight),
        anchor=g,
        decision=LengthDecision(
            kappa=float(kappa),
            captured=float(captured),
            contrast=float(np.sqrt(contrast_sq)),
            bound=float(np.sqrt(bound_sq)),
            noise=noise,
            weight=weight,
        ),
    )
