"""The message-reliability ``var~mean`` model: a **monotone-increasing P-spline** fit of the biological
dispersion ``σ²_bio(mean)`` via a **Poisson offset** (:meth:`MonotoneVarMean.fit_offset`).

``bp_solver`` fits this each sweep pass on the frozen-snapshot cross-node density residuals — the observed
squared residual decomposes as ``E[R²] = σ²_bio(μ) + V_p`` (``V_p`` = the computed Poisson sampling floor),
so the fit learns only the biological EXCESS — and queries :meth:`MonotoneVarMean.predict` to set each
message's precision ``τ = 1/(σ²_bio + ρ_src/L_src)``.

**Why monotone P-splines.** The var~mean trend must be *monotone* (dispersion non-decreasing in mean),
*smooth* (no isotonic staircases, no LOESS non-monotone fringe), and *flexible at the fringes* (capture
makes the data bimodal: a depleted low-mean floor and an enriched high-mean plateau). A monotone P-spline
is the only candidate that is all three at once. The smoothing penalty (one λ, chosen by GCV) controls
flexibility; the basis size ``k`` is just an upper bound.

**The SCAM reparameterization (Pya & Wood).** Fit ``σ²_bio = Σ β_j B_j(log mean)`` on a B-spline basis;
enforce ``β`` monotone-increasing by the unconstrained reparameterization ``β = cumsum([α₀, exp(α₁), …,
exp(α_{k-1})])`` — ``exp(·) > 0`` guarantees ``β_j ≥ β_{j-1}`` for ANY real ``α``. Fit ``α`` by penalized
robust IRLS least squares ``‖√w·(z − Bβ)‖² + λ‖D₂β‖²``. Pure numpy + scipy; the fit runs once per sweep
pass over ~10³–10⁴ points in well under a second.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.interpolate import BSpline
from scipy.optimize import minimize

__all__ = ["MonotoneVarMean", "MonotoneMean", "BlendedVarMean"]


_EPS = 1.0e-12
_DEFAULT_K = 18  # B-spline basis size (upper bound on flexibility); λ controls the actual smoothness
_DEFAULT_DEGREE = 3  # cubic


def _bspline_design(xq: np.ndarray, knots: np.ndarray, degree: int, k: int) -> np.ndarray:
    """``(n, k)`` B-spline design matrix evaluating each basis function at ``xq``."""
    xq = np.asarray(xq, dtype=np.float64)
    out = np.empty((xq.shape[0], k), dtype=np.float64)
    coef = np.zeros(k, dtype=np.float64)
    for j in range(k):
        coef[j] = 1.0
        out[:, j] = BSpline(knots, coef, degree, extrapolate=True)(xq)
        coef[j] = 0.0
    return out


def _knots(x_lo: float, x_hi: float, k: int, degree: int) -> np.ndarray:
    """Clamped uniform knot vector for ``k`` basis functions of the given degree."""
    n_interior = k - degree - 1
    interior = (
        np.linspace(x_lo, x_hi, n_interior + 2)[1:-1] if n_interior > 0 else np.empty(0)
    )
    return np.concatenate([[x_lo] * (degree + 1), interior, [x_hi] * (degree + 1)])


def _setup_spline(x_sorted: np.ndarray, k: int, degree: int):
    """The shared monotone-spline scaffolding: clamped knots over ``[x_sorted[0], x_sorted[-1]]``, the
    ``(n, k)`` B-spline design at ``x_sorted``, and the 2nd-difference penalty ``P = D₂ᵀD₂``. Used by BOTH the
    variance fit (:meth:`MonotoneVarMean.fit_offset`) and the mean fit (:meth:`MonotoneMean.fit`)."""
    x_lo, x_hi = float(x_sorted[0]), float(x_sorted[-1])
    kn = _knots(x_lo, x_hi, k, degree)
    B = _bspline_design(x_sorted, kn, degree, k)
    D2 = np.diff(np.eye(k), 2, axis=0)
    return kn, B, D2.T @ D2, x_lo, x_hi


def _select_lambda(B: np.ndarray, y: np.ndarray, P: np.ndarray, n: int, n_lambda: int):
    """Pick the smoothing parameter ``λ`` by GCV on the (unweighted, linear) penalized fit ``min‖y−Bβ‖²+λβᵀPβ``.
    Returns ``(lam, edf)``. Shared by the variance and mean fits."""
    lams = np.logspace(-4, 4, n_lambda)
    BtB = B.T @ B
    Bty = B.T @ y
    best = (np.inf, lams[0], float(B.shape[1]))
    for lam in lams:
        A = BtB + lam * P
        try:
            beta = np.linalg.solve(A, Bty)
            edf = float(np.trace(np.linalg.solve(A, BtB)))
        except np.linalg.LinAlgError:
            continue
        rss = float(np.sum((y - B @ beta) ** 2))
        denom = (n - edf) ** 2
        gcv = n * rss / denom if denom > _EPS else np.inf
        if gcv < best[0]:
            best = (gcv, lam, edf)
    return best[1], best[2]


def _slope_evidence_weight(B: np.ndarray, y: np.ndarray, P: np.ndarray, lam: float, edf: float, n: int) -> float:
    """The CONTINUOUS evidence weight ``w ∈ [0, 1]`` that the fitted spline explains ``y`` BEYOND chance —
    the data-derived shrinkage of the learned mean toward the flat (intercept-only) model: ``ê_eff = w·ê +
    (1−w)·flat``. ``w → 1`` under a real ``x→y`` relationship (apply the fit), ``w → 0`` under none (collapse
    to flat). It replaces the binary permutation significance gate with a smooth function of the fit's OWN
    evidence, so one uniform path handles enriched and non-enriched data alike — the fit self-collapsing where
    ``x`` carries no signal.

    **Empirical-Bayes posterior weight.** Treat the flat (intercept-only) model H0 and the GCV-fitted spline
    H1 as two competing hypotheses with EQUAL prior odds; ``w`` is the posterior probability of the sloped
    model — the textbook ``w = ev(H1) / (ev(H1) + ev(H0))``. With the Schwarz (unit-information) prior, the
    log model-evidence is the BIC, so the log evidence-ratio is

        ΔBIC = n·log(RSS_flat / RSS_spline) − (edf − 1)·log n,        w = σ(½·ΔBIC) = ev1 / (ev1 + ev0),

    where ``RSS_spline`` is the residual SS of the penalized fit at the GCV-selected ``λ`` (the fit's own
    out-of-sample prediction-error surrogate, the GCV score's numerator) and ``RSS_flat`` that of the
    intercept-only model. The first term is the deviance the spline buys (the data evidence ``log RSS_flat /
    RSS_spline``, ``→ 0`` exactly when **flat fits as well as sloped** ⇒ ``w → 0``); the second is the Occam
    penalty for the ``edf − 1`` extra dof the spline spends — without it a microscopic chance slope at large
    ``n`` would be promoted (the directive's large-``n`` failure mode). Every quantity comes from the GCV fit
    (``RSS`` from the penalized β, ``edf`` the effective dof, ``n`` the fit-node count); no new constant —
    the ``log n`` penalty is the Schwarz prior's, not a tuned threshold.

    Behaviour: real enrichment drives ``RSS_flat ≫ RSS_spline`` ⇒ ``ΔBIC ≫ 0`` ⇒ ``w → 1`` regardless of
    effect SIZE (a noisy-but-real enrichment, e.g. unstranded capture, is still applied fully — the posterior
    SATURATES, unlike an effect-size weight); uniform gDNA gives ``RSS_flat ≈ RSS_spline`` and the Occam term
    drives ``ΔBIC < 0`` ⇒ ``w → 0``, recovering the byte-flat global where the GCV spline can hallucinate a
    faint slope."""
    BtB = B.T @ B
    A = BtB + lam * P
    try:
        beta = np.linalg.solve(A, B.T @ y)
    except np.linalg.LinAlgError:
        return 0.0
    rss_spline = float(np.sum((y - B @ beta) ** 2))
    rss_flat = float(np.sum((y - float(np.mean(y))) ** 2))
    df1 = max(edf - 1.0, _EPS)  # the spline's extra dof over the flat (intercept-only) model
    # log Bayes factor (Schwarz / unit-information prior): the deviance the spline buys minus its Occam penalty.
    log_bf = 0.5 * n * np.log(max(rss_flat, _EPS) / max(rss_spline, _EPS)) - 0.5 * df1 * np.log(max(float(n), 2.0))
    return float(1.0 / (1.0 + np.exp(-np.clip(log_bf, -700.0, 700.0))))


@dataclass(frozen=True, slots=True)
class MonotoneVarMean:
    """A fitted monotone-increasing ``σ²_bio(mean)`` curve (the Poisson-offset fit, :meth:`fit_offset`).
    Call :meth:`predict` for the variance at a mean; the inputs are stored for the diagnostics dataframe.
    ``coeffs`` are the monotone curve in LINEAR variance space (≥ 0)."""

    knots: np.ndarray
    degree: int
    coeffs: np.ndarray  # the monotone β coefficients (σ²_bio per basis, linear space)
    x_lo: float  # log-mean fit range (predictions are clipped to it — flat extrapolation)
    x_hi: float
    fit_mean: np.ndarray  # input fit points (for diagnostics)
    fit_var: np.ndarray
    edf: float  # effective degrees of freedom (GCV)
    lam: float  # selected smoothing parameter

    @classmethod
    def fit_offset(
        cls,
        mean: np.ndarray,
        raw: np.ndarray,
        offset: np.ndarray,
        weight: np.ndarray | None = None,
        *,
        k: int = _DEFAULT_K,
        degree: int = _DEFAULT_DEGREE,
        robust_iters: int = 2,
        n_lambda: int = 40,
        population_spread: bool = False,
    ) -> "MonotoneVarMean":
        """Fit the **biological dispersion** ``σ²_bio(μ) ≥ 0`` from a Poisson-OFFSET decomposition
        (``imputation_variance_model.md`` §3): the observed squared residual ``raw = (ρ_dst − ρ_src)²`` has
        ``E[raw] = σ²_bio(μ) + V_p`` where ``V_p = offset`` is the COMPUTED Poisson sampling floor
        (``ρ_src/L_src + ρ_dst/L_dst``). We learn only the EXCESS ``σ²_bio``.

        Fit in **LINEAR** variance space (not log-log): the de-offset response ``z = raw − offset`` may be
        negative per point (a sampling-dominated pair), which is an ordinary least-squares residual here —
        no ``log`` of a negative. The monotone curve is ``σ²_bio(μ) = max(spline(log μ), 0)``; the
        non-negativity floor (in :meth:`predict`) is the genuine zero-dispersion region. Heteroskedastic
        IRLS weights ``w = 1/(σ²_bio + V_p)²`` (a squared-Gaussian-residual response has variance ∝ its
        mean²) make the low-μ regime — where the offset dominates — count, plus bisquare robustness. The
        count re-enters at APPLY as ``τ = 1/(σ²_bio(μ) + ρ_src/L_src)`` (the second §0 count→precision
        channel), NOT here.

        ``population_spread`` switches the estimand from the **reliability-weighted MEDIAN** squared
        residual (the default, for the message σ²_g where the Poisson-offset reliability genuinely down-
        weights thin-count edges) to the **unweighted conditional MEAN** ``E[z | log μ]`` — the unbiased
        population CV⁻² that a per-node *hyperprior* needs (a node's predictive density spread is one
        population draw, so every node counts equally). It drops (a) the heteroskedastic ``1/(σ²+V_p)²``
        IRLS reweight — the response is *already* a variance, so down-weighting its heavy χ²₁ upper tail
        discards the very signal being estimated; and (b) the bisquare ``robust_iters`` — same reason; and
        (c) the supplied reliability ``weight`` (the eff-length weight systematically excludes the short /
        variable nodes whose mis-prediction *is* the between-node spread). Used by the enrichment-aware
        global σ²_bio fit and the symmetric RNA-message σ²_bio fit (which sets the cross-regime AMBIG
        message precision); the gDNA message + seed fits keep the default. No new constant — it only
        changes which statistic of the existing ``z = raw − offset`` is fitted.
        """
        mean = np.asarray(mean, dtype=np.float64)
        raw = np.asarray(raw, dtype=np.float64)
        off = np.asarray(offset, dtype=np.float64)
        if population_spread:
            wt = np.ones_like(off)  # the conditional-MEAN estimand counts every node equally (no eff-len weight)
            robust_iters = 0        # the χ²₁ upper tail IS the spread — never rejected as an outlier
        else:
            wt = np.ones_like(off) if weight is None else np.asarray(weight, dtype=np.float64)
        ok = (
            np.isfinite(mean) & (mean > _EPS)
            & np.isfinite(raw) & (raw >= 0.0)
            & np.isfinite(off) & (off >= 0.0)
            & np.isfinite(wt) & (wt > 0.0)
        )
        mean, raw, off, wt = mean[ok], raw[ok], off[ok], wt[ok]
        if mean.size < max(k, 8):
            return cls._constant_offset(mean, raw, off, degree, wt)
        order = np.argsort(mean)
        mean, raw, off, wt = mean[order], raw[order], off[order], wt[order]
        x = np.log(mean)
        z = raw - off  # de-offset excess (LINEAR space; may be < 0 for sampling-dominated pairs)
        if float(x[-1]) - float(x[0]) < _EPS:
            return cls._constant_offset(mean, raw, off, degree, wt)
        kn, B, P, x_lo, x_hi = _setup_spline(x, k, degree)
        # λ by GCV on the (unweighted, linear) penalized fit of the de-offset response z.
        lam, edf = _select_lambda(B, z, P, mean.size, n_lambda)

        # heteroskedastic IRLS (w = 1/total², total = σ²_bio + V_p) + bisquare robustness. Init s=0 ⇒
        # total = V_p (the Poisson floor weights the fit by the predictor reliability).
        total = np.maximum(off, _EPS)
        coeffs = None
        for it in range(robust_iters + 1):
            # population-spread: w = wt (uniform) — the response z is already a variance, so the
            # heteroskedastic 1/total² (which pins the fit to the offset-dominated low-μ tail) is dropped.
            w = wt if population_spread else wt / np.maximum(total, _EPS) ** 2
            if it > 0:
                r = (z - B @ coeffs) * np.sqrt(w)
                s_mad = 1.4826 * float(np.median(np.abs(r - np.median(r)))) + _EPS
                u = r / (4.685 * s_mad)
                w = w * np.where(np.abs(u) < 1.0, (1.0 - u**2) ** 2, 0.0)
            coeffs = _fit_monotone(B, z, w, P, lam, k)
            if it == robust_iters:
                break
            total = np.maximum(B @ coeffs, 0.0) + off
        return cls(
            knots=kn, degree=degree, coeffs=coeffs, x_lo=x_lo, x_hi=x_hi,
            fit_mean=mean, fit_var=np.maximum(z, 0.0), edf=float(edf), lam=float(lam),
        )

    @classmethod
    def _constant_offset(cls, mean, raw, off, degree, weight=None) -> "MonotoneVarMean":
        """Too-few-points / no-log-range fallback for the offset fit: a FLAT ``σ²_bio`` = the
        reliability-weighted mean excess ``max(Σw·(raw−off)/Σw, 0)``, ``w = weight/V_p²``."""
        mean = np.asarray(mean, dtype=np.float64)
        raw = np.asarray(raw, dtype=np.float64)
        off = np.asarray(off, dtype=np.float64)
        sw = np.ones_like(off) if weight is None else np.asarray(weight, dtype=np.float64)
        if mean.size:
            w = sw / np.maximum(off, _EPS) ** 2
            c = max(float(np.sum(w * (raw - off)) / max(float(np.sum(w)), _EPS)), 0.0)
            lx = float(np.log(np.maximum(np.median(mean), _EPS)))
        else:
            c, lx = 0.0, 0.0
        x_lo, x_hi = lx - 1.0, lx + 1.0
        kn = _knots(x_lo, x_hi, max(degree + 1, 2), degree)
        coeffs = np.full(kn.size - degree - 1, c)  # constant curve (≥0)
        return cls(
            knots=kn, degree=degree, coeffs=coeffs, x_lo=x_lo, x_hi=x_hi,
            fit_mean=mean, fit_var=np.maximum(raw - off, 0.0) if mean.size else mean,
            edf=1.0, lam=np.inf,
        )

    def predict(self, mean: np.ndarray) -> np.ndarray:
        """``σ²_bio`` at each ``mean`` — the monotone curve, clipped to the fit range (flat extrapolation).

        The sweep fits the curve on every snapshot edge, so queries interpolate; the curve learns
        genuinely-high dispersion at RNA-rich means directly, so the message precision ``1/(σ²_bio + ...)``
        already humbles the prior where prediction is unreliable. ``coeffs`` are σ²_bio in linear space,
        floored at 0 (the genuine zero-dispersion regions)."""
        m = np.asarray(mean, dtype=np.float64)
        xq = np.clip(np.log(np.maximum(m, _EPS)), self.x_lo, self.x_hi)
        raw = BSpline(self.knots, self.coeffs, self.degree, extrapolate=True)(xq)
        return np.maximum(raw, 0.0)

    def to_dataframe(self, n_curve: int = 100) -> pd.DataFrame:
        """Diagnostics: the input fit points + the fitted curve sampled on a grid (for the companion
        plotting script). Columns: ``kind`` (point/curve), ``mean``, ``var``."""
        pts = pd.DataFrame({"kind": "point", "mean": self.fit_mean, "var": self.fit_var})
        if self.fit_mean.size:
            grid = np.logspace(
                np.log10(self.fit_mean.min()), np.log10(self.fit_mean.max()), n_curve
            )
            crv = pd.DataFrame({"kind": "curve", "mean": grid, "var": self.predict(grid)})
            return pd.concat([pts, crv], ignore_index=True)
        return pts


@dataclass(frozen=True, slots=True)
class BlendedVarMean:
    """A convex blend ``w·hi + (1−w)·lo`` of two fitted ``σ²_bio(mean)`` curves, evaluated at predict time.
    Used for the RNA-message reliability: ``hi`` is the population-CONDITIONAL-MEAN fit (the unbiased
    hyperprior estimand that capture needs), ``lo`` the reliability-weighted-MEDIAN fit (correct off-capture);
    the SAME continuous evidence weight ``w`` that shrinks the enrichment transfer also selects between the two
    estimands — ``w → 1`` under capture (the mean), ``w → 0`` off-capture (the median) — so there is no
    population_spread flag and no capture detector, just the one uniform path."""

    hi: "MonotoneVarMean"
    lo: "MonotoneVarMean"
    w: float

    def predict(self, mean: np.ndarray) -> np.ndarray:
        if self.w >= 1.0 - _EPS:
            return self.hi.predict(mean)
        if self.w <= _EPS:
            return self.lo.predict(mean)
        return self.w * self.hi.predict(mean) + (1.0 - self.w) * self.lo.predict(mean)


@dataclass(frozen=True, slots=True)
class MonotoneMean:
    """A fitted monotone-increasing **conditional MEAN** ``ê(x) = E[y | x]`` — the enrichment transfer
    ``ê(z) = E[ρ_g | z]`` of the enrichment-aware calibration design (``enrichment_aware_calibration.md`` §J).

    Same monotone-P-spline backend as :class:`MonotoneVarMean` (``_setup_spline`` + GCV-λ + the monotone
    reparameterization ``_fit_monotone``), but fitting the MEAN of a positive response, in LOG-LOG space
    (``log y ~ monotone-spline(log x)``) so the multiplicative-noise response and the near-log-linear
    relationship are handled naturally; the spline reduces to ~linear when the data is log-linear.

    A **count-recalibration scale** ``s`` (stored) makes the prediction net-unbiased on the supplied
    ``recal_weight``: ``s = Σ(rw·y) / Σ(rw·ŷ_raw)`` so ``Σ ŷ·rw = Σ y·rw`` on the fit set (with
    ``rw = E_dst`` ⇒ net-unbiased on the gDNA COUNT — the §F fix, no NB needed). :meth:`predict` returns
    ``s·exp(spline(clip(log x)))`` (≥ 0; flat extrapolation outside the fit range)."""

    knots: np.ndarray
    degree: int
    coeffs: np.ndarray  # the monotone β of the LOG-response spline
    x_lo: float  # log-x fit range (predictions clipped to it)
    x_hi: float
    scale: float  # count-recalibration scale (linear space)
    fit_x: np.ndarray  # input fit points (diagnostics)
    fit_y: np.ndarray
    edf: float
    lam: float
    w_enrich: float  # the continuous slope-evidence weight (:func:`_slope_evidence_weight`); 0 for a flat fit

    @classmethod
    def fit(
        cls,
        x: np.ndarray,
        y: np.ndarray,
        weight: np.ndarray | None = None,
        *,
        recal_weight: np.ndarray | None = None,
        k: int = _DEFAULT_K,
        degree: int = _DEFAULT_DEGREE,
        robust_iters: int = 2,
        n_lambda: int = 40,
    ) -> "MonotoneMean":
        """Fit ``ê(x) = E[y | x]`` (both ``x, y > 0``). ``weight`` = the fit reliability weight; ``recal_weight``
        = the net-balance weight (``E_dst``; defaults to ``weight``). Robust (bisquare) IRLS on the log-response.
        Falls back to a flat ``ê`` (the recal-weighted mean of ``y``) on too-few-points / no-x-range."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        w = np.ones_like(x) if weight is None else np.asarray(weight, dtype=np.float64)
        rw = w if recal_weight is None else np.asarray(recal_weight, dtype=np.float64)
        ok = (
            np.isfinite(x) & (x > _EPS) & np.isfinite(y) & (y > _EPS)
            & np.isfinite(w) & (w > 0.0) & np.isfinite(rw) & (rw >= 0.0)
        )
        x, y, w, rw = x[ok], y[ok], w[ok], rw[ok]
        if x.size < max(k, 8):
            return cls._constant(x, y, w, rw, degree)
        order = np.argsort(x)
        x, y, w, rw = x[order], y[order], w[order], rw[order]
        lx, ly = np.log(x), np.log(y)
        if float(lx[-1]) - float(lx[0]) < _EPS:
            return cls._constant(x, y, w, rw, degree)
        kn, B, P, x_lo, x_hi = _setup_spline(lx, k, degree)
        lam, edf = _select_lambda(B, ly, P, x.size, n_lambda)
        # the continuous slope-evidence weight on the SAME GCV fit (the EB sloped-vs-flat posterior) — the
        # data-derived shrinkage that replaces the binary significance gate (the caller blends ê toward flat by this w).
        w_enrich = _slope_evidence_weight(B, ly, P, lam, edf, x.size)

        coeffs = None
        for it in range(robust_iters + 1):
            wi = w.copy()
            if it > 0:
                r = (ly - B @ coeffs) * np.sqrt(w)
                s_mad = 1.4826 * float(np.median(np.abs(r - np.median(r)))) + _EPS
                u = r / (4.685 * s_mad)
                wi = w * np.where(np.abs(u) < 1.0, (1.0 - u**2) ** 2, 0.0)
            coeffs = _fit_monotone(B, ly, wi, P, lam, k)
            if it == robust_iters:
                break
        yhat = np.exp(B @ coeffs)
        scale = float(np.sum(rw * y) / max(float(np.sum(rw * yhat)), _EPS))
        return cls(
            knots=kn, degree=degree, coeffs=coeffs, x_lo=x_lo, x_hi=x_hi, scale=scale,
            fit_x=x, fit_y=y, edf=float(edf), lam=float(lam), w_enrich=float(w_enrich),
        )

    @classmethod
    def _constant(cls, x, y, w, rw, degree) -> "MonotoneMean":
        """Flat ``ê`` fallback: the recal-weight-weighted mean of ``y`` (coeffs ≡ 0 ⇒ exp = 1, scale = the mean).
        Also the explicit OFF-CAPTURE / no-enrichment value the caller forces via the significance gate."""
        if x.size:
            c = float(np.sum(rw * y) / max(float(np.sum(rw)), _EPS))
            lx = float(np.log(np.maximum(np.median(x), _EPS)))
        else:
            c, lx = 0.0, 0.0
        x_lo, x_hi = lx - 1.0, lx + 1.0
        kn = _knots(x_lo, x_hi, max(degree + 1, 2), degree)
        coeffs = np.zeros(kn.size - degree - 1)  # exp(0)=1 ⇒ predict = scale (flat)
        return cls(
            knots=kn, degree=degree, coeffs=coeffs, x_lo=x_lo, x_hi=x_hi, scale=c,
            fit_x=x, fit_y=y, edf=1.0, lam=np.inf, w_enrich=0.0,  # a flat fallback carries no slope evidence
        )

    @classmethod
    def constant(cls, value: float, x_ref: float = 1.0) -> "MonotoneMean":
        """A flat transfer ``ê(x) ≡ value`` — the off-capture / ``z``-unusable fallback (= the genome-wide
        ``ρ_global``), built without data so the caller can substitute it when the significance gate fails."""
        lx = float(np.log(max(x_ref, _EPS)))
        x_lo, x_hi = lx - 1.0, lx + 1.0
        kn = _knots(x_lo, x_hi, max(_DEFAULT_DEGREE + 1, 2), _DEFAULT_DEGREE)
        coeffs = np.zeros(kn.size - _DEFAULT_DEGREE - 1)
        return cls(
            knots=kn, degree=_DEFAULT_DEGREE, coeffs=coeffs, x_lo=x_lo, x_hi=x_hi,
            scale=float(max(value, 0.0)), fit_x=np.zeros(0), fit_y=np.zeros(0), edf=1.0, lam=np.inf,
            w_enrich=0.0,  # the explicit flat ρ_global transfer — no enrichment evidence
        )

    def predict(self, x: np.ndarray) -> np.ndarray:
        """``ê(x) = s·exp(spline(clip(log x)))`` — the net-calibrated monotone mean, ≥ 0, flat-extrapolated."""
        xq = np.clip(np.log(np.maximum(np.asarray(x, dtype=np.float64), _EPS)), self.x_lo, self.x_hi)
        raw = BSpline(self.knots, self.coeffs, self.degree, extrapolate=True)(xq)
        return self.scale * np.maximum(np.exp(raw), 0.0)


def _fit_monotone(B, y, w, P, lam, k) -> np.ndarray:
    """Penalized weighted least squares with the monotone reparameterization β=cumsum([α₀,exp(α₁..)])."""
    sw = np.sqrt(np.maximum(w, 0.0))

    def beta_of(a):
        # clamp the log-increment before exp: a pure overflow guard for the optimizer's transient large
        # steps (exp(700) overflows float64); exp(50)≈5e21 dwarfs any real increment, so it never
        # constrains a genuine fit — it only keeps the objective finite so L-BFGS doesn't stall on inf.
        return np.cumsum(np.concatenate([[a[0]], np.exp(np.minimum(a[1:], 50.0))]))

    def obj(a):
        beta = beta_of(a)
        r = sw * (y - B @ beta)
        db = P @ beta
        return float(r @ r + lam * (beta @ db))

    # init from the monotonized linear solution
    beta0 = np.maximum.accumulate(np.linalg.solve(B.T @ (w[:, None] * B) + lam * P, B.T @ (w * y)))
    inc0 = np.diff(beta0)
    a0 = np.concatenate([[beta0[0]], np.log(np.maximum(inc0, 1e-6))])
    sol = minimize(obj, a0, method="L-BFGS-B", options={"maxiter": 2000, "ftol": 1e-12})
    return beta_of(sol.x)
