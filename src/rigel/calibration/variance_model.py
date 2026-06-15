"""Consolidated var~mean variance model — a **monotone-increasing P-spline (SCAM)** fitter plus the
DIRECT and IMPUTATION point builders. One authoritative variance model for the calibration's three uses:
the per-node count-prior precision (observable=DIRECT, imputed=IMPUTATION), the global gDNA-baseline
variance, and the per-edge RNA odds-coupling. See ``docs/calibration/per_node_deconv_hierarchy_design.md``
§16-17.

**Why monotone P-splines.** The var~mean trend must be *monotone* (variance non-decreasing in mean —
a model where variance falls as mean rises is wrong), *smooth* (no isotonic staircases, no LOESS
non-monotone fringe), and *flexible at the fringes* (capture makes the data bimodal: a depleted low-mean
floor and an enriched high-mean plateau). A monotone P-spline is the only candidate that is all three at
once (validated head-to-head in ``scripts/debug/scam_var_mean.py``). The smoothing penalty (one λ, chosen
by GCV) controls flexibility; the basis size ``k`` is just an upper bound.

**The SCAM reparameterization (Pya & Wood).** Model ``log var = Σ β_j B_j(log mean)`` on a B-spline basis;
enforce ``β`` monotone-increasing by the unconstrained reparameterization ``β = cumsum([α₀, exp(α₁), …,
exp(α_{k-1})])`` — ``exp(·) > 0`` guarantees ``β_j ≥ β_{j-1}`` for ANY real ``α``, turning a constrained
problem into an unconstrained one. Fit ``α`` by penalized (optionally robust/IRLS) least squares
``‖√w·(y − Bβ)‖² + λ‖D₂β‖²``. Pure numpy + scipy (already a dependency); no new packages, and — since the
fit runs once per calibration over ~10³–10⁴ points in well under a second — no C++ port is warranted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.interpolate import BSpline
from scipy.optimize import minimize

__all__ = [
    "MonotoneVarMean",
    "VarMeanPoints",
    "varmean_points",
    "fit_direct_varmean",
    "fit_imputation_varmean",
]

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
    """A fitted monotone-increasing var~mean curve (log-log). Call :meth:`predict` for the variance at a
    mean; the inputs are stored for the diagnostics dataframe."""

    knots: np.ndarray
    degree: int
    coeffs: np.ndarray  # the monotone β coefficients (log-var per basis)
    x_lo: float  # log-mean fit range (predictions are clipped to it — flat extrapolation)
    x_hi: float
    fit_mean: np.ndarray  # input fit points (for diagnostics)
    fit_var: np.ndarray
    edf: float  # effective degrees of freedom (GCV)
    lam: float  # selected smoothing parameter

    @classmethod
    def fit(
        cls,
        mean: np.ndarray,
        var: np.ndarray,
        *,
        k: int = _DEFAULT_K,
        degree: int = _DEFAULT_DEGREE,
        robust_iters: int = 2,
        n_lambda: int = 40,
    ) -> "MonotoneVarMean":
        """Fit ``log var ~ monotone-increasing spline in log mean`` from scattered ``(mean, var)`` points.

        ``k`` is the basis upper bound (the GCV-selected λ controls smoothness); ``robust_iters`` adds
        bisquare IRLS passes (0 ⇒ plain least squares). Points with non-positive mean/var are dropped.
        """
        mean = np.asarray(mean, dtype=np.float64)
        var = np.asarray(var, dtype=np.float64)
        ok = np.isfinite(mean) & np.isfinite(var) & (mean > _EPS) & (var > _EPS)
        mean, var = mean[ok], var[ok]
        if mean.size < max(k, 8):
            # too few points for a stable spline — fall back to a monotone power-law (var = a·mean^b,
            # b≥0 enforced) so the curve is still usable (e.g. tiny toy scenarios).
            return cls._powerlaw_fallback(mean, var, k, degree)

        order = np.argsort(mean)
        mean, var = mean[order], var[order]
        x, y = np.log(mean), np.log(var)
        x_lo, x_hi = float(x[0]), float(x[-1])
        kn = _knots(x_lo, x_hi, k, degree)
        B = _bspline_design(x, kn, degree, k)
        D2 = np.diff(np.eye(k), 2, axis=0)
        P = D2.T @ D2

        # --- λ by GCV on the (linear, unweighted) penalized fit ---
        lams = np.logspace(-4, 4, n_lambda)
        BtB = B.T @ B
        Bty = B.T @ y
        best = (np.inf, lams[0], _DEFAULT_K)
        for lam in lams:
            A = BtB + lam * P
            try:
                beta = np.linalg.solve(A, Bty)
                edf = float(np.trace(np.linalg.solve(A, BtB)))
            except np.linalg.LinAlgError:
                continue
            rss = float(np.sum((y - B @ beta) ** 2))
            denom = (mean.size - edf) ** 2
            gcv = mean.size * rss / denom if denom > _EPS else np.inf
            if gcv < best[0]:
                best = (gcv, lam, edf)
        lam, edf = best[1], best[2]

        # --- monotone fit (reparameterized), with optional bisquare IRLS robustness ---
        w = np.ones(mean.size)
        coeffs = None
        for it in range(robust_iters + 1):
            coeffs = _fit_monotone(B, y, w, P, lam, k)
            if it == robust_iters:
                break
            resid = y - B @ coeffs
            s = 1.4826 * np.median(np.abs(resid - np.median(resid))) + _EPS
            u = resid / (4.685 * s)
            w = np.where(np.abs(u) < 1.0, (1.0 - u**2) ** 2, 0.0)

        return cls(
            knots=kn, degree=degree, coeffs=coeffs, x_lo=x_lo, x_hi=x_hi,
            fit_mean=mean, fit_var=var, edf=float(edf), lam=float(lam),
        )

    @classmethod
    def _powerlaw_fallback(cls, mean, var, k, degree) -> "MonotoneVarMean":
        """Monotone power-law ``var = a·mean^b`` (b clipped ≥ 0) for too-few-points cases."""
        if mean.size >= 2:
            x, y = np.log(mean), np.log(var)
            b, a = np.polyfit(x, y, 1)
            b = max(b, 0.0)
            x_lo, x_hi = float(x.min()), float(x.max())
        else:
            a, b, x_lo, x_hi = 0.0, 0.0, -1.0, 1.0
        # represent the line as a (degenerate) 2-coef spline so predict() is uniform
        kn = _knots(x_lo, x_hi, max(degree + 1, 2), degree)
        kk = kn.size - degree - 1
        xs = np.linspace(x_lo, x_hi, kk)
        coeffs = np.maximum.accumulate(a + b * xs)  # monotone by construction
        return cls(
            knots=kn, degree=degree, coeffs=coeffs, x_lo=x_lo, x_hi=x_hi,
            fit_mean=np.asarray(mean), fit_var=np.asarray(var), edf=2.0, lam=np.inf,
        )

    def predict(self, mean: np.ndarray) -> np.ndarray:
        """Variance at each ``mean`` (the monotone curve, clipped to the fit range = flat extrapolation).

        No separate "confidence/extrapolation guard" is needed: under the all-gDNA-init iterative paradigm
        the model is fit on EVERY node (all means are seen ⇒ predictions interpolate), and the curve
        **learns** genuinely-high variance at RNA-rich means directly (the boundary→region prediction error
        is large there), so the downstream count precision ``1/var`` already yields the count where it is
        unreliable. See ``CALIBRATION_PLAN_v2.md`` §4-5."""
        m = np.asarray(mean, dtype=np.float64)
        xq = np.clip(np.log(np.maximum(m, _EPS)), self.x_lo, self.x_hi)
        return np.exp(BSpline(self.knots, self.coeffs, self.degree, extrapolate=True)(xq))

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
        return np.cumsum(np.concatenate([[a[0]], np.exp(a[1:])]))

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


# ======================================================================================
# Point builders — DIRECT (region-observable) vs IMPUTATION (region imputed from boundaries)
# ======================================================================================


@dataclass(frozen=True, slots=True)
class VarMeanPoints:
    """Per-region var~mean fit points: ``mean`` ½/⅓-measurement mean density, ``raw_var`` the
    inter-measurement disagreement, ``region_observable`` whether the contained region anchored the
    triplet (DIRECT) or was imputed from its boundaries only (IMPUTATION)."""

    mean: np.ndarray
    raw_var: np.ndarray
    region_observable: np.ndarray  # bool — True ⇒ DIRECT point, False ⇒ IMPUTATION point


def varmean_points(
    substrate, region_arrays, region_eff_len, fl_mean, gdna_views=None
) -> VarMeanPoints:
    """Extract the gDNA var~mean fit points from the boundary↔region↔boundary triplets.

    For each region the available density measurements are: the contained own-density (if the region is
    count-observable), and the left/right boundary-crossing densities (if those boundaries are
    count-observable). Where ≥2 measurements exist, the point is ``mean = average``, ``raw_var = the
    disagreement`` (the unbiased k-sample variance of the mean). Tagged ``region_observable`` so the
    caller can split DIRECT (region anchored the estimate) from IMPUTATION (region imputed from boundaries
    only). This is the consolidated successor to ``density_model.density_variance_curve``'s extraction.

    ``gdna_views`` = ``(contained, left, right)`` per-view gDNA fragment counts: the current iteration's
    gDNA estimate (``CALIBRATION_PLAN_v2`` §2/§5). ``None`` ⇒ the **all-gDNA init** (total unspliced =
    every fragment is gDNA); subsequent passes pass the deconvolved gDNA so the curve refines as RNA is
    removed.
    """
    from .density_model import count_observable_masks
    from .run_fill import same_ref_left_right

    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    r = sig.shape[0]
    L = np.maximum(np.asarray(region_eff_len, dtype=np.float64), _EPS)
    inv_fl = 1.0 / float(fl_mean)

    if gdna_views is None:
        def total(view):
            return view.n_unspliced_pos.astype(np.float64) + view.n_unspliced_neg.astype(np.float64)
        c_cnt, l_cnt, r_cnt = total(substrate.contained), total(substrate.left), total(substrate.right)
    else:
        c_cnt, l_cnt, r_cnt = (np.asarray(v, dtype=np.float64) for v in gdna_views)

    rco, bco = count_observable_masks(sig, ref_id)
    ls, rs = same_ref_left_right(ref_id)
    la = np.zeros(r, bool)
    rb = np.zeros(r, bool)
    if r > 1:
        la[1:] = bco[:-1] & ls[1:]  # left side observable (neighbour's boundary)
        rb[:-1] = bco[:-1] & rs[:-1]  # right side observable
    d_left = np.where(la, l_cnt * inv_fl, np.nan)
    d_right = np.where(rb, r_cnt * inv_fl, np.nan)
    c_ok = rco & (L > _EPS)
    contained = np.where(c_ok, c_cnt / L, np.nan)

    obs = np.stack([d_left, d_right, contained])
    msk = np.stack([la, rb, c_ok])
    kcount = msk.sum(0).astype(np.float64)
    mu = np.where(msk, np.nan_to_num(obs), 0.0).sum(0) / np.maximum(kcount, 1.0)
    dev2 = np.where(msk, (np.nan_to_num(obs) - mu) ** 2, 0.0).sum(0)
    raw = np.where(kcount > 1.0, dev2 / np.maximum(kcount - 1.0, _EPS), np.nan) / np.maximum(kcount, 1.0)
    sel = (kcount >= 2.0) & np.isfinite(mu) & (mu > _EPS) & np.isfinite(raw) & (raw > _EPS)
    return VarMeanPoints(mean=mu[sel], raw_var=raw[sel], region_observable=c_ok[sel])


def fit_direct_varmean(points: VarMeanPoints, **kw) -> MonotoneVarMean:
    """DIRECT var~mean — fit over the region-observable points (the region anchored the estimate). Feeds
    the count precision at count-observable nodes and the global-baseline variance."""
    m = points.region_observable
    return MonotoneVarMean.fit(points.mean[m], points.raw_var[m], **kw)


def fit_imputation_varmean(points: VarMeanPoints, **kw) -> MonotoneVarMean:
    """IMPUTATION var~mean — fit over the region-imputed points (boundaries only, no region truth). Feeds
    the count precision at imputed/AMBIG nodes (properly humble — it measures boundary→region disagreement,
    not the anchor's own confidence — which is the §15 over-trust fix)."""
    m = ~points.region_observable
    return MonotoneVarMean.fit(points.mean[m], points.raw_var[m], **kw)
