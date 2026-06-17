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
from scipy.special import digamma

__all__ = [
    "MonotoneVarMean",
    "pair_imputation_points",
    "fit_pair_imputation_varmean",
]


def _jensen_offset(dof: np.ndarray) -> np.ndarray:
    """Per-point Jensen bias correction for ``log(raw_var)`` of a small-dof variance estimate.

    A sample variance ``s²`` of ``dof`` degrees of freedom has ``E[log s²] = log σ² + ψ(dof/2) −
    log(dof/2)`` (the log of a scaled χ²). So ``log s²`` is **downward**-biased by
    ``Δ = log(dof/2) − ψ(dof/2) ≥ 0`` (e.g. dof=1 → +1.2704, dof=2 → +0.5772, →0 as dof→∞). Adding
    ``Δ`` to the fit target ``log s²`` removes the over-confidence (inflates the recovered variance)
    so a k=2 disagreement (dof=1) is not read as a tighter variance than it is.
    """
    d = np.asarray(dof, dtype=np.float64)
    half = np.maximum(d, _EPS) / 2.0
    return np.log(half) - digamma(half)

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
        dof: np.ndarray | None = None,
        k: int = _DEFAULT_K,
        degree: int = _DEFAULT_DEGREE,
        robust_iters: int = 2,
        n_lambda: int = 40,
    ) -> "MonotoneVarMean":
        """Fit ``log var ~ monotone-increasing spline in log mean`` from scattered ``(mean, var)`` points.

        ``k`` is the basis upper bound (the GCV-selected λ controls smoothness); ``robust_iters`` adds
        bisquare IRLS passes (0 ⇒ plain least squares). Points with non-positive mean/var are dropped.

        ``dof`` (optional, per-point degrees of freedom of each ``var`` estimate) applies the **Jensen
        df-offset** to the fit target: ``log var`` of a small-dof variance is downward-biased, so we add
        the positive correction ``Δ = log(dof/2) − ψ(dof/2)`` (see :func:`_jensen_offset`) before fitting,
        inflating the recovered variance to remove the over-confidence. ``None`` ⇒ no offset (exact
        back-compat).
        """
        mean = np.asarray(mean, dtype=np.float64)
        var = np.asarray(var, dtype=np.float64)
        dof_arr = None if dof is None else np.asarray(dof, dtype=np.float64)
        ok = np.isfinite(mean) & np.isfinite(var) & (mean > _EPS) & (var > _EPS)
        if dof_arr is not None:
            ok &= np.isfinite(dof_arr)
        mean, var = mean[ok], var[ok]
        off = None if dof_arr is None else _jensen_offset(dof_arr[ok])
        if mean.size < max(k, 8):
            # too few points for a stable spline — fall back to a monotone power-law (var = a·mean^b,
            # b≥0 enforced) so the curve is still usable (e.g. tiny toy scenarios).
            return cls._powerlaw_fallback(mean, var, k, degree, offset=off)

        order = np.argsort(mean)
        mean, var = mean[order], var[order]
        x, y = np.log(mean), np.log(var)
        if off is not None:
            y = y + off[order]  # Jensen-correct the fit target (does NOT alter stored fit_var)
        x_lo, x_hi = float(x[0]), float(x[-1])
        if x_hi - x_lo < _EPS:
            # All means coincide (no log-range) — a clamped spline knot vector degenerates (zero internal
            # knots ⇒ BSpline error). Fall back to the monotone power-law, which handles a single distinct
            # x gracefully. (Happens when the count-observable own-count densities are all equal, e.g.
            # uniform-gDNA toy scenarios — the DIRECT own-count extractor.)
            return cls._powerlaw_fallback(mean, var, k, degree, offset=off)
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
    def _powerlaw_fallback(cls, mean, var, k, degree, *, offset=None) -> "MonotoneVarMean":
        """Monotone power-law ``var = a·mean^b`` (b clipped ≥ 0) for too-few-points cases.

        ``offset`` (optional, aligned to ``mean``/``var``) is the per-point Jensen df-offset added to
        ``log var`` before the line fit — same correction as the spline path.
        """
        if mean.size >= 2 and float(np.log(mean).max() - np.log(mean).min()) > _EPS:
            x, y = np.log(mean), np.log(var)
            if offset is not None:
                y = y + np.asarray(offset, dtype=np.float64)
            b, a = np.polyfit(x, y, 1)
            b = max(b, 0.0)
            x_lo, x_hi = float(x.min()), float(x.max())
        elif mean.size >= 1:
            # A single distinct mean (no log-range): a flat curve at the mean log-variance. Widen the
            # range by ±1 (in log-mean) so the clamped knot vector / BSpline are well-defined; predict
            # clips flat to it anyway (every query reads the same constant variance).
            lx = float(np.log(np.maximum(mean[0], _EPS)))
            ly = float(np.mean(np.log(var))) + (0.0 if offset is None else float(np.mean(offset)))
            a, b, x_lo, x_hi = ly, 0.0, lx - 1.0, lx + 1.0
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
# Point builder — the node-PAIR IMPUTATION reliability (the Step-2 imputation substrate)
# ======================================================================================


# --------------------------------------------------------------------------------------
# The node-PAIR imputation reliability (CALIBRATION_PLAN v6 §7) — the GENERIC point builder. This module
# is the var~mean MACHINE (densities in, curve out); it is the Step-2 imputation substrate (no production
# consumer in the Step-1 strand+global solve — the per-component density assemblies arrive with Step 2).
# --------------------------------------------------------------------------------------


def pair_imputation_points(
    region_density,
    left_density,
    right_density,
    *,
    region_eligible,
    left_ok,
    right_ok,
) -> tuple[np.ndarray, np.ndarray]:
    """``(means, raw_var)`` node-pair fit points — the generic, component-agnostic builder (v6 §7).

    A *pair* is two adjacent nodes — an imputed destination region and one observable boundary side (the
    single predictor). For each adjacency ``(observable side → eligible region)`` the point is, in **density
    space** (eff-len-normalized, log-log): ``mean = region_density[r]`` (the *queried* axis — the sweep reads
    the reliability at the region's CURRENT density), ``raw_var = (region_density − side_density)²`` (the
    **full** single-predictor residual, a dof=1 estimate). A region with both flanks eligible contributes
    **two** points; one flank → one point (the densification the deleted both-sides triplet missed). The
    reverse (region→boundary) direction is not emitted (its ``μ = d_side`` is queried only once boundaries
    become solved nodes — Phase C).

    Consumed by the gDNA fit (:func:`fit_pair_imputation_varmean`) and, in Step 2, by the per-component
    (per-strand RNA) imputation fits that will pool their points. All arrays are
    per-region (length R); ``*_ok[r]`` means that flank exists, is observable, and carries component mass.
    Standard filters (finite, ``> _EPS``) applied.
    """
    rd = np.asarray(region_density, dtype=np.float64)
    ld = np.asarray(left_density, dtype=np.float64)
    rrd = np.asarray(right_density, dtype=np.float64)
    elig = np.asarray(region_eligible, dtype=bool)
    lok = np.asarray(left_ok, dtype=bool) & elig
    rok = np.asarray(right_ok, dtype=bool) & elig

    means = np.concatenate([rd[lok], rd[rok]])
    raw = np.concatenate([(rd[lok] - ld[lok]) ** 2, (rd[rok] - rrd[rok]) ** 2])
    sel = np.isfinite(means) & (means > _EPS) & np.isfinite(raw) & (raw > _EPS)
    return means[sel], raw[sel]


def fit_pair_imputation_varmean(
    region_density,
    left_density,
    right_density,
    *,
    region_eligible,
    left_ok,
    right_ok,
    ref_id=None,
    **kw,
) -> MonotoneVarMean:
    """The **gDNA** node-pair imputation-reliability var~mean (v6 §7) — :func:`pair_imputation_points` fit.

    Builds the node-pair points from the per-region gDNA densities and fits the monotone curve (Jensen
    dof=1, the single-observation convention). ``ref_id`` is accepted for interface symmetry (the adjacency
    masks already encode same-ref existence); it is unused.
    """
    means, raw = pair_imputation_points(
        region_density,
        left_density,
        right_density,
        region_eligible=region_eligible,
        left_ok=left_ok,
        right_ok=right_ok,
    )
    return MonotoneVarMean.fit(means, raw, dof=np.ones(means.shape[0]), **kw)
