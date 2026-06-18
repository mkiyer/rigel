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

__all__ = ["MonotoneVarMean"]


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
        *,
        k: int = _DEFAULT_K,
        degree: int = _DEFAULT_DEGREE,
        robust_iters: int = 2,
        n_lambda: int = 40,
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
        """
        mean = np.asarray(mean, dtype=np.float64)
        raw = np.asarray(raw, dtype=np.float64)
        off = np.asarray(offset, dtype=np.float64)
        ok = (
            np.isfinite(mean) & (mean > _EPS)
            & np.isfinite(raw) & (raw >= 0.0)
            & np.isfinite(off) & (off >= 0.0)
        )
        mean, raw, off = mean[ok], raw[ok], off[ok]
        if mean.size < max(k, 8):
            return cls._constant_offset(mean, raw, off, degree)
        order = np.argsort(mean)
        mean, raw, off = mean[order], raw[order], off[order]
        x = np.log(mean)
        z = raw - off  # de-offset excess (LINEAR space; may be < 0 for sampling-dominated pairs)
        x_lo, x_hi = float(x[0]), float(x[-1])
        if x_hi - x_lo < _EPS:
            return cls._constant_offset(mean, raw, off, degree)
        kn = _knots(x_lo, x_hi, k, degree)
        B = _bspline_design(x, kn, degree, k)
        D2 = np.diff(np.eye(k), 2, axis=0)
        P = D2.T @ D2

        # λ by GCV on the (unweighted, linear) penalized fit of the de-offset response z.
        lams = np.logspace(-4, 4, n_lambda)
        BtB = B.T @ B
        Btz = B.T @ z
        best = (np.inf, lams[0], float(_DEFAULT_K))
        for lam in lams:
            A = BtB + lam * P
            try:
                beta = np.linalg.solve(A, Btz)
                edf = float(np.trace(np.linalg.solve(A, BtB)))
            except np.linalg.LinAlgError:
                continue
            rss = float(np.sum((z - B @ beta) ** 2))
            denom = (mean.size - edf) ** 2
            gcv = mean.size * rss / denom if denom > _EPS else np.inf
            if gcv < best[0]:
                best = (gcv, lam, edf)
        lam, edf = best[1], best[2]

        # heteroskedastic IRLS (w = 1/total², total = σ²_bio + V_p) + bisquare robustness. Init s=0 ⇒
        # total = V_p (the Poisson floor weights the fit by the predictor reliability).
        total = np.maximum(off, _EPS)
        coeffs = None
        for it in range(robust_iters + 1):
            w = 1.0 / np.maximum(total, _EPS) ** 2
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
    def _constant_offset(cls, mean, raw, off, degree) -> "MonotoneVarMean":
        """Too-few-points / no-log-range fallback for the offset fit: a FLAT ``σ²_bio`` = the
        reliability-weighted mean excess ``max(Σw·(raw−off)/Σw, 0)``, ``w = 1/V_p²``."""
        mean = np.asarray(mean, dtype=np.float64)
        raw = np.asarray(raw, dtype=np.float64)
        off = np.asarray(off, dtype=np.float64)
        if mean.size:
            w = 1.0 / np.maximum(off, _EPS) ** 2
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
