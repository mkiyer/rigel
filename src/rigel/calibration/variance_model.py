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
    "DirectPoints",
    "direct_points",
    "fit_direct_varmean",
    "fit_pair_imputation_varmean",
    "fit_pair_imputation_rna_varmean",
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
# Point builders — DIRECT (own-count Poisson) + the unified node-PAIR IMPUTATION reliability
# ======================================================================================


@dataclass(frozen=True, slots=True)
class DirectPoints:
    """DIRECT var~mean fit points (the count-observable regions only — the region anchored its own
    density). ``mean`` the region's measurement mean (its own contained density, averaged with any
    observable flanking-crossing densities), ``raw_var`` the inter-measurement disagreement (the unbiased
    k-sample variance of that mean), ``kcount`` the number of density measurements behind the point — the
    variance estimate has ``dof = kcount − 1``, which feeds the Jensen df-offset (:func:`fit_direct_varmean`).

    This is the **thin own-count extractor** for the count-observable region set; the IMPUTATION
    (non-observable) extraction is now the node-PAIR builder (:func:`fit_pair_imputation_varmean`). Side
    densities are eff-len-consistent with the contained density (``E[min(ℓ,L_side)]``, not ``E[ℓ]``), so a
    uniform gDNA field yields ~zero region-vs-side disagreement."""

    mean: np.ndarray
    raw_var: np.ndarray
    kcount: np.ndarray  # float — measurements behind the point (dof = kcount − 1)


def direct_points(
    substrate, region_arrays, region_eff_len, boundary_side_eff_len, gdna_views=None
) -> DirectPoints:
    """Extract the DIRECT var~mean fit points — the count-observable regions' own-density estimate.

    For each count-observable region (intergenic / intron — gDNA by signature) the available density
    measurements are its own contained density and any observable flanking-crossing densities; the point is
    ``mean = average``, ``raw_var = the disagreement`` (the unbiased k-sample variance of the mean), tagged
    with ``kcount`` (``dof = kcount − 1`` for the Jensen offset). One point per count-observable region with
    ≥2 measurements (a disagreement must exist).

    Densities are **eff-len-consistent**: the contained count by ``region_eff_len`` (``E[max(0,L−ℓ)]``) and a
    flanking-side crossing count by the per-side **density** length ``boundary_side_eff_len[r] = E[min(ℓ,L_r)]``
    — NOT the count/power length ``fl_mean = E[ℓ]``. Under uniform gDNA at density ρ a region's contained and
    each flanking side both read ρ (factor-1), so a uniform field yields ~zero disagreement; the bugged
    ``count/E[ℓ]`` side density read ``ρ·E[min(ℓ,L)]/E[ℓ] < ρ`` for short flanks, manufacturing a spurious
    region-vs-side disagreement under uniform gDNA.

    ``gdna_views`` = ``(contained, left, right)`` per-view gDNA fragment counts: the current iteration's
    gDNA estimate. ``None`` ⇒ the **all-gDNA init** (total unspliced = every fragment is gDNA).
    """
    from .density_model import count_observable_masks
    from .run_fill import same_ref_left_right

    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    r = sig.shape[0]
    L = np.maximum(np.asarray(region_eff_len, dtype=np.float64), _EPS)
    # Per-side gDNA DENSITY length E[min(ℓ, L_side)] (= boundary_side_eff_len[r] for both of region r's
    # sides, which lie inside region r), the divisor for a flanking crossing MASS/count, NOT fl_mean.
    inv_side = 1.0 / np.maximum(np.asarray(boundary_side_eff_len, dtype=np.float64), _EPS)

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
    d_left = np.where(la, l_cnt * inv_side, np.nan)
    d_right = np.where(rb, r_cnt * inv_side, np.nan)
    c_ok = rco & (L > _EPS)
    contained = np.where(c_ok, c_cnt / L, np.nan)

    obs = np.stack([d_left, d_right, contained])
    msk = np.stack([la, rb, c_ok])
    kcount = msk.sum(0).astype(np.float64)
    mu = np.where(msk, np.nan_to_num(obs), 0.0).sum(0) / np.maximum(kcount, 1.0)
    dev2 = np.where(msk, (np.nan_to_num(obs) - mu) ** 2, 0.0).sum(0)
    raw = np.where(kcount > 1.0, dev2 / np.maximum(kcount - 1.0, _EPS), np.nan) / np.maximum(kcount, 1.0)
    # DIRECT subset: count-observable regions, ≥2 measurements (so a disagreement exists), finite & positive.
    sel = c_ok & (kcount >= 2.0) & np.isfinite(mu) & (mu > _EPS) & np.isfinite(raw) & (raw > _EPS)
    return DirectPoints(mean=mu[sel], raw_var=raw[sel], kcount=kcount[sel])


def fit_direct_varmean(points: DirectPoints, **kw) -> MonotoneVarMean:
    """DIRECT var~mean — fit over the count-observable regions (the region anchored its own density). Feeds
    the count precision at count-observable nodes and the global-baseline variance. Jensen-corrected via the
    per-point ``dof = kcount − 1`` (a 2-measurement point is dof=1, a triplet dof=2)."""
    return MonotoneVarMean.fit(points.mean, points.raw_var, dof=points.kcount - 1.0, **kw)


# --------------------------------------------------------------------------------------
# The unified node-PAIR imputation reliability (CALIBRATION_PLAN_v5 §3, user-locked)
# --------------------------------------------------------------------------------------


def fit_pair_imputation_varmean(
    region_density,
    left_density,
    right_density,
    *,
    region_eligible,
    left_ok,
    right_ok,
    ref_id,
    **kw,
) -> MonotoneVarMean:
    """The unified **node-PAIR** imputation-reliability var~mean (CALIBRATION_PLAN_v5 §3).

    A *pair* is exactly two adjacent nodes — one imputed-region (destination) and one observable boundary
    side (the single predictor / source). For each adjacency ``(observable boundary side → eligible
    region)`` the fit point is, in **density space** (eff-len-normalized, log-log):

    * ``mean    = region_density[r]`` — the *queried* axis (the sweep reads ``τ_count`` at the region's
      CURRENT density only — fit-and-query-on-the-same-axis, the 2a contract);
    * ``raw_var = (region_density[r] − side_density[r])²`` — the **full** single-predictor imputation error
      (NOT the ``¼(d_L−d_R)²`` variance-of-the-mean of the deleted triplet; the genuine residual of
      predicting the region density from one boundary side), a single-observation (dof=1, Jensen) estimate.

    A region with **both** flanks eligible contributes **TWO** points (one per flanking side); with one
    eligible flank, **one** point — the densification the both-sides triplet missed (it required *both*
    sides). The reverse (region→boundary) direction is **NOT emitted**: its ``μ = d_side`` is never queried
    (until boundaries become solved nodes in Phase C) and emitting it would corrupt the GCV/IRLS via
    duplicated residuals (v5 §3).

    Parameters are per-region (length R) arrays:

    * ``region_density[r]`` — the imputed-region (dest) CURRENT density (eff-len-normalized);
    * ``left_density[r]`` / ``right_density[r]`` — the density at region ``r``'s LEFT / RIGHT observable
      boundary side (the single predictor);
    * ``region_eligible[r]`` — the region is an imputed dest (gDNA: ``~region_count_observable``; RNA:
      single-strand-exon). Only eligible regions contribute fit points;
    * ``left_ok[r]`` / ``right_ok[r]`` — that flanking side exists (same-ref), is observable, and carries
      component mass > 0.

    Standard filters (finite, ``mean > EPS``, ``raw_var > EPS``). ``MonotoneVarMean.fit(dof=ones)`` —
    Jensen dof=1, exactly the convention the deleted imputation builders used. ``ref_id`` is accepted for
    interface symmetry / future per-reference logic; the adjacency masks already encode same-ref existence.

    **One builder, symmetric** for gDNA (this signature) and RNA (:func:`fit_pair_imputation_rna_varmean`,
    which assembles the per-strand densities and concatenates the per-strand POINT arrays into one fit).
    """
    rd = np.asarray(region_density, dtype=np.float64)
    ld = np.asarray(left_density, dtype=np.float64)
    rrd = np.asarray(right_density, dtype=np.float64)
    elig = np.asarray(region_eligible, dtype=bool)
    lok = np.asarray(left_ok, dtype=bool) & elig
    rok = np.asarray(right_ok, dtype=bool) & elig

    # One point per (observable side → eligible region) adjacency: mean = region density (the queried
    # axis); raw_var = the FULL single-predictor residual (region density − that side's density)².
    means = np.concatenate([rd[lok], rd[rok]])
    raw = np.concatenate([(rd[lok] - ld[lok]) ** 2, (rd[rok] - rrd[rok]) ** 2])
    sel = np.isfinite(means) & (means > _EPS) & np.isfinite(raw) & (raw > _EPS)
    m = means[sel]
    return MonotoneVarMean.fit(m, raw[sel], dof=np.ones(m.shape[0]), **kw)


def fit_pair_imputation_rna_varmean(
    substrate,
    region_arrays,
    region_eff_len_rna,
    rna_boundary_side_eff_len,
    *,
    gdna_frac,
    left_gdna_frac,
    right_gdna_frac,
    cleaned_left=None,
    cleaned_right=None,
    **kw,
) -> MonotoneVarMean:
    """The **RNA** node-pair imputation-reliability var~mean (CALIBRATION_PLAN_v5 §3) — the symmetric RNA
    twin of :func:`fit_pair_imputation_varmean`, pooling the ``+`` and ``−`` strands into one fit.

    STANDALONE / INERT: this is **not** wired into the live solve (Phase B/C); it is the RNA reliability the
    unified node solver will consume, validated standalone for now.

    Per strand ``s`` (``+`` then ``−``), in **RNA density space** (FL-consistent — ``region_eff_len_rna`` =
    ``E_rna[max(0,L−ℓ)]`` for the region; ``rna_boundary_side_eff_len`` = ``E_rna[min(ℓ,L_side)]`` per side):

    * ``region_density[r] = region_unspliced_s · (1 − f_g[r]) / region_eff_len_rna[r]`` — the region's
      CURRENT RNA density on strand ``s``;
    * ``side_density[r]   = (unspliced_RNA_s + spliced_s) / rna_boundary_side_eff_len[r]`` per side, with
      ``unspliced_RNA_s = side_unspliced_s · (1 − side_gdna_frac)`` (``side_gdna_frac`` from the per-side
      ``StrandSplit.gdna_frac`` with the strand-cleaned-count fallback where ``NaN``, identical to the
      Phase-A RNA builder), and ``spliced_s`` the side's same-strand spliced crossing count. The side MASS is
      divided by the per-side **density** length ``E_rna[min(ℓ,L_side)]`` — NOT the count/power length
      ``rna_fl_mean = E_rna[ℓ]`` (the symmetric RNA twin of the gDNA node-pair eff-len fix): so under a
      uniform single-strand RNA field a region and its flanking sides read the SAME RNA density (factor-1)
      and the residual reflects only the true RNA structure, not a short-flank normalization offset;
    * ``region_eligible[r]`` — single-strand exon on ``s`` (``TS_POS`` for ``+`` / ``TS_NEG`` for ``−``, has
      the exon bit, ``~region_count_observable``);
    * ``*_ok[r]`` — that flank is count-observable AND carries ``> 0`` same-strand spliced mass (so the side
      is RNA-anchored).

    The two strands are pooled by **concatenating the per-strand POINT arrays** (mean, raw_var) returned by
    the core builder — NOT the per-strand node arrays — so no spurious ``+``↔``−`` seam adjacency is ever
    formed (a region eligible on ``+`` is never paired against a ``−`` side). Each strand calls
    :func:`fit_pair_imputation_varmean`'s *point assembly* and the pooled points are fit once.
    """
    from .density_model import count_observable_masks
    from .run_fill import same_ref_left_right
    from .signature import BIT_EXON_NEG, BIT_EXON_POS, TS_NEG, TS_POS

    sig = np.asarray(region_arrays.signature).astype(np.int64)
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    r = sig.shape[0]
    L = np.maximum(np.asarray(region_eff_len_rna, dtype=np.float64), _EPS)
    # Per-side RNA DENSITY length E_rna[min(ℓ,L_side)] — the divisor for a side's RNA MASS (NOT rna_fl_mean).
    inv_side = 1.0 / np.maximum(np.asarray(rna_boundary_side_eff_len, dtype=np.float64), _EPS)

    fg = np.clip(np.asarray(gdna_frac, dtype=np.float64), 0.0, 1.0)
    lgf = np.asarray(left_gdna_frac, dtype=np.float64)
    rgf = np.asarray(right_gdna_frac, dtype=np.float64)

    c = substrate.contained
    c_pos = np.asarray(c.n_unspliced_pos, dtype=np.float64)
    c_neg = np.asarray(c.n_unspliced_neg, dtype=np.float64)
    left = substrate.left
    right = substrate.right
    l_pos = np.asarray(left.n_unspliced_pos, dtype=np.float64)
    l_neg = np.asarray(left.n_unspliced_neg, dtype=np.float64)
    rt_pos = np.asarray(right.n_unspliced_pos, dtype=np.float64)
    rt_neg = np.asarray(right.n_unspliced_neg, dtype=np.float64)
    left_spl_s = np.asarray(left.n_spliced_sense, dtype=np.float64)
    right_spl_s = np.asarray(right.n_spliced_sense, dtype=np.float64)
    left_total = l_pos + l_neg
    right_total = rt_pos + rt_neg
    cl = None if cleaned_left is None else np.asarray(cleaned_left, dtype=np.float64)
    cr = None if cleaned_right is None else np.asarray(cleaned_right, dtype=np.float64)

    def _rna_frac(side_gdna_frac, side_total, cleaned):
        """Per-side RNA fraction (1 − side_gdna_frac), cleaned-count fallback where the gdna_frac is NaN."""
        out = 1.0 - np.where(np.isfinite(side_gdna_frac), side_gdna_frac, 0.0)
        if cleaned is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                fb = np.where(side_total > _EPS, (side_total - cleaned) / side_total, 0.0)
            out = np.where(np.isfinite(side_gdna_frac), out, np.clip(fb, 0.0, 1.0))
        return out

    left_rna_frac = _rna_frac(lgf, left_total, cl)
    right_rna_frac = _rna_frac(rgf, right_total, cr)

    rco, bco = count_observable_masks(sig, ref_id)
    ls, rs = same_ref_left_right(ref_id)
    la = np.zeros(r, bool)  # left side count-observable
    rb = np.zeros(r, bool)  # right side count-observable
    if r > 1:
        la[1:] = bco[:-1] & ls[1:]
        rb[:-1] = bco[:-1] & rs[:-1]
    has_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0

    rna_frac_region = 1.0 - fg
    all_means, all_raw = [], []
    for ts_s, region_unspl_s, l_unspl_s, r_unspl_s in (
        (TS_POS, c_pos, l_pos, rt_pos),
        (TS_NEG, c_neg, l_neg, rt_neg),
    ):
        region_density = region_unspl_s * rna_frac_region / L
        d_left = (l_unspl_s * left_rna_frac + left_spl_s) * inv_side
        d_right = (r_unspl_s * right_rna_frac + right_spl_s) * inv_side
        eligible = (ts == ts_s) & has_exon & (~rco)
        lok = la & (left_spl_s > 0.0) & eligible
        rok = rb & (right_spl_s > 0.0) & eligible
        # Assemble the per-strand POINTS (same form as the core builder); pool by concatenation below.
        means = np.concatenate([region_density[lok], region_density[rok]])
        raw = np.concatenate(
            [(region_density[lok] - d_left[lok]) ** 2, (region_density[rok] - d_right[rok]) ** 2]
        )
        all_means.append(means)
        all_raw.append(raw)

    means = np.concatenate(all_means)
    raw = np.concatenate(all_raw)
    sel = np.isfinite(means) & (means > _EPS) & np.isfinite(raw) & (raw > _EPS)
    m = means[sel]
    return MonotoneVarMean.fit(m, raw[sel], dof=np.ones(m.shape[0]), **kw)
