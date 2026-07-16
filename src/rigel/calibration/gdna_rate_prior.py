"""The pass-0 gDNA-rate prior — a Fixed-Kernel Poisson(-lognormal) Mixture NPMLE in count/rate space.

Replaces the density KDE (`gdna_density_prior`) + the scalar floor + the count-space global
(`bp_solver._floor_estimate` / `_global_logprior`). The gDNA count of a node is ``g ~ Poisson(ρ·E)`` for a
latent rate ``ρ ~ P(ρ)``; we estimate the population ``P(ρ)`` and project it onto each node's ``f_g``.

Design (docs/calibration/npmle_struggles.md §8-9, npmle_roadmap.md):

* **Substrate = ALL nodes at total density** — pass-0 asserts ``f_g=1`` (the gDNA count = the total unspliced
  count) at *zero belief precision*, so the fit is belief-free and covers every probe-targeted region (a later
  REFIT uses the solved belief: ``g_hat = f_g·count`` with width ``τ = √Var(log f_g)``).
* **Fixed-Kernel Poisson Mixture** — ``P(log ρ) = Σ_j w_j·N(log ρ; log ρ_j, h²)``, a mixture of fixed-width
  (bandwidth ``h``) Gaussian kernels on a log-ρ grid; the weights are fit by plain EM (monotone,
  arithmetic-only, deterministic — NO spline / λ / matrix-inversion, so none of the GCV-λ nondeterminism).
  The fixed kernel forbids any peak sharper than ``h`` ⇒ smooth, never a bed-of-nails (the raw NPMLE's
  Kiefer-Wolfowitz atomicity). ``h`` is the familiar KDE bandwidth knob.
* **Count/rate space** — the Poisson is zero-native (``count=0`` ⇒ ``e^{−ρE}`` ⇒ mass at ρ→0), so the sparse,
  zero-inflated data needs no boundary correction; the per-node likelihood width is set by the physics
  (broad for low-count/short-E — the discreteness comb smooths itself away).
* **Extremely weak** — projected onto a node this prior is worth ``n_eff ≈ 0.15`` pseudo-observations (never
  >1), so the strand likelihood and the boundary messages dominate and can peel RNA out of the ``f_g=1`` start
  (measured: `scripts/debug/npmle_variance.py`).

Performance: the per-node likelihood depends only on the ``(count, log-E, τ)`` cell, so ~10⁶ nodes collapse to
~10³–10⁴ weighted cells → the fit is sub-second and needs no C++/numba.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

__all__ = ["GdnaRatePrior"]

_EPS = 1.0e-12


def _lse(a: np.ndarray, axis: int) -> np.ndarray:
    """Numerically-stable log-sum-exp along ``axis`` (drops the axis)."""
    mx = np.max(a, axis=axis, keepdims=True)
    mx = np.where(np.isfinite(mx), mx, 0.0)
    return (mx + np.log(np.sum(np.exp(a - mx), axis=axis, keepdims=True))).squeeze(axis)


def _collapse(g_hat, eff, var_g, *, dlog, dt, g_eps=1.0e-6):
    """Collapse nodes to weighted cells of (near-)identical likelihood, keyed on
    ``(log-ĝ bin | a shared ĝ≈0 bin, log-E bin, τ bin)``. Representative = the in-cell mean; weight = the
    node count. ``dlog``/``dt`` trade memory/compute for binning error — a perf knob, not an accuracy constant.
    Returns ``(gc, Ec, t2c, wc)``."""
    logE = np.log(np.maximum(eff, _EPS))
    tau = np.sqrt(np.maximum(var_g, 0.0))
    gi = np.where(
        g_hat < g_eps,
        np.int64(-(10**9)),
        np.round(np.log(np.maximum(g_hat, _EPS)) / dlog).astype(np.int64),
    )
    cols = np.stack(
        [gi, np.round(logE / dlog).astype(np.int64), np.round(tau / dt).astype(np.int64)], axis=1
    )
    _u, inv, wc = np.unique(cols, axis=0, return_inverse=True, return_counts=True)
    inv = np.asarray(inv).reshape(-1)
    wc = wc.astype(np.float64)
    mean = lambda v: np.bincount(inv, weights=v) / wc  # noqa: E731
    return mean(g_hat), mean(eff), mean(np.maximum(var_g, 0.0)), wc


def _cell_loglik(g_hat, eff, var_g, log_rho, *, n_gh, chunk=4096):
    """Per-cell log-likelihood over the log-ρ grid — the **Poisson-lognormal**: ``g ~ Poisson(ρ·E)`` with the
    belief placing ``log g ~ N(log ĝ, τ²)``; marginalise g by Gauss-Hermite quadrature. ``ĝ=0`` ⇒ ``g_q=0`` ⇒
    ``e^{−ρE}`` (the exact zero anchor). ``τ=0`` (pass-0) ⇒ all quadrature nodes coincide ⇒ the pure Poisson
    (we drop the quadrature entirely then). Returns ``(n_cell, G)``."""
    pure_poisson = float(np.max(var_g)) <= _EPS
    x, wq = (
        (np.zeros(1), np.array([np.sqrt(np.pi)]))
        if pure_poisson
        else np.polynomial.hermite.hermgauss(n_gh)
    )
    lwq = np.log(wq) - 0.5 * np.log(np.pi)  # ∫ e^{−x²}·/√π  (lognormal expectation)
    n = g_hat.shape[0]
    out = np.empty((n, log_rho.shape[0]), dtype=np.float64)
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        tau = np.sqrt(var_g[lo:hi])[:, None]  # (c,1)
        gq = g_hat[lo:hi][:, None] * np.exp(np.sqrt(2.0) * tau * x[None, :])  # (c,Q)
        logE = np.log(eff[lo:hi])[:, None, None]  # (c,1,1)
        gq3 = gq[:, :, None]  # (c,Q,1)
        lam = np.exp(log_rho[None, None, :] + logE)  # (c,1,G) expected counts ρ·E
        logpois = gq3 * (log_rho[None, None, :] + logE) - lam - gammaln(gq3 + 1.0)  # (c,Q,G)
        out[lo:hi] = _lse(lwq[None, :, None] + logpois, axis=1)
    return out


def _kernel_matrix(log_rho, h):
    """``K[i,j] = N(log ρ_i; log ρ_j, h²)`` normalised so each component column integrates to 1 over the grid.
    ``K @ w`` renders the fixed-kernel mixture density on the grid."""
    d = log_rho[:, None] - log_rho[None, :]
    kk = np.exp(-0.5 * (d / h) ** 2)
    return kk / np.maximum(kk.sum(axis=0, keepdims=True), _EPS)


def _em_weights(logL, logwc, n_iter, tol):
    """Mixture EM for the component weights over collapsed cells (monotone, arithmetic-only, deterministic).
    ``logL`` = (n_cell, G) cell-vs-component log-likelihood; ``logwc`` = (n_cell,1) log cell counts."""
    logw = np.full(logL.shape[1], -np.log(logL.shape[1]))
    prev = np.exp(logw)
    for _ in range(n_iter):
        lp = logw[None, :] + logL
        logw = _lse(logwc + (lp - _lse(lp, axis=1)[:, None]), axis=0)
        logw -= _lse(logw, axis=0)
        cur = np.exp(logw)
        if np.max(np.abs(cur - prev)) < tol:
            break
        prev = cur
    return np.exp(logw)


@dataclass(frozen=True, slots=True)
class GdnaRatePrior:
    """The fitted population gDNA-rate prior ``P(log ρ)``, pre-evaluated (smoothed + interior-floored) on a
    log-ρ grid. :meth:`logprior` projects it onto a node's ``f_g`` axis — the ``(n_nodes, K)`` additive term
    the sweep adds to each node's ψ (replacing the KDE + floor + global)."""

    log_rho: np.ndarray  #: (G,) log-rate grid
    logP: (
        np.ndarray
    )  #: (G,) log P(log ρ) — kernel-rendered, uniform-floored inside the support, real tails
    bandwidth: float  #: h (decades), the kernel width
    n_cells: int  #: collapsed-cell count (diagnostic)

    @classmethod
    def fit(
        cls,
        g_hat: np.ndarray,
        eff: np.ndarray,
        var_g: np.ndarray | None = None,
        *,
        bandwidth: float = 0.15,
        n_grid: int = 200,
        em_iters: int = 150,
        em_tol: float = 1.0e-5,
        log_dlog: float = 0.1,
        tau_dt: float = 0.25,
        n_gh: int = 7,
        floor_eps: float = 0.02,
    ) -> "GdnaRatePrior":
        """Fit the Fixed-Kernel Poisson-lognormal Mixture on the per-node gDNA ``(count ĝ, eff-length E,
        belief width τ²=var_g)``. ``var_g=None`` ⇒ ``τ=0`` (pass-0, pure Poisson at ``f_g=1``, ``ĝ=count``).

        The grid spans the observed density support: the top is the MAX density (so the projection never
        extrapolates ``ρ_g=f_g·M/E`` above the grid — avoiding the clamped-tail false-positive), the bottom
        is below the bulk (zero-anchor room; the low side clamps safely to the depleted level)."""
        g_hat = np.asarray(g_hat, dtype=np.float64)
        eff = np.maximum(np.asarray(eff, dtype=np.float64), _EPS)
        # A belief width of ∞ marks an UNSOLVED node (no information). Such nodes carry no gDNA count
        # (``ĝ=0``) so their likelihood is the exact zero anchor regardless of τ — but ``0·exp(√2·∞·x)`` is
        # NaN, so normalise the width to 0 here rather than propagating a NaN through the EM.
        var_g = np.zeros_like(g_hat) if var_g is None else np.asarray(var_g, dtype=np.float64)
        var_g = np.where(np.isfinite(var_g), np.maximum(var_g, 0.0), 0.0)
        ln10 = np.log(10.0)
        h = float(bandwidth) * ln10  # bandwidth in natural-log units

        dens = g_hat / eff
        ld = np.log(dens[dens > 0.0])
        if ld.size == 0:
            ld = np.array([np.log(1.0 / float(np.max(eff)))])
        lo = float(np.percentile(ld, 0.5)) - 3.0 * h  # zero-anchor room (clamp-safe low side)
        hi = (
            float(np.max(ld)) + 2.0 * h
        )  # ≥ max density ⇒ no upward extrapolation in the projection
        log_rho = np.linspace(lo, hi, int(n_grid))

        gc, ec, tc, wc = _collapse(g_hat, eff, var_g, dlog=log_dlog, dt=tau_dt)
        logL = _cell_loglik(gc, ec, tc, log_rho, n_gh=int(n_gh))
        kk = _kernel_matrix(log_rho, h)
        # convolve each cell likelihood with the fixed kernel, EM the weights, render the mixture density.
        ln = np.exp(logL - logL.max(axis=1, keepdims=True))
        logL_sm = np.log(np.maximum(ln @ kk, _EPS))
        w = _em_weights(logL_sm, np.log(wc)[:, None], em_iters, em_tol)
        dens_grid = kk @ w  # (G,) smooth mixture density (real Gaussian tails, width ≥ h)
        dens_grid = dens_grid / max(float(dens_grid.sum()), _EPS)

        # interior uniform floor (fills −∞ valleys between modes) bounded to the fitted support [q0.5, q99.5];
        # OUTSIDE the support the real Gaussian tails are left intact (the false-positive suppression).
        cdf = np.cumsum(w) / max(float(w.sum()), _EPS)
        blo = log_rho[int(np.searchsorted(cdf, 0.005))]
        bhi = log_rho[min(int(np.searchsorted(cdf, 0.995)), n_grid - 1)]
        in_supp = (log_rho >= blo) & (log_rho <= bhi)
        step = float(log_rho[1] - log_rho[0])
        uni = np.where(in_supp, floor_eps / max(float(bhi - blo), _EPS) * step, 0.0)
        logP = np.log(
            np.maximum(np.where(in_supp, (1.0 - floor_eps) * dens_grid, dens_grid) + uni, _EPS)
        )
        return cls(log_rho=log_rho, logP=logP, bandwidth=float(bandwidth), n_cells=int(gc.shape[0]))

    def logprior(self, fg_grid, mass, eff):
        """Project onto the ``f_g`` solve grid → ``(n_nodes, K)`` additive term = ``log P(log ρ_g)`` at
        ``ρ_g = f_g·M/E``. **Bare** — no reference prior, no measure term, no Jacobian.

        This is the ``ψ_λ = strand + logP_g + logP_r`` result (`prior_ramp_and_bp_roadmap.md` §2). The
        latents are rates; conditioning on ``M`` contributes no ``f_g``-dependent Jacobian; and because
        ``logP`` is a density in **log**-rate, its conversion to a linear-rate density (``−log ρ_g``, i.e.
        ``−log f_g`` up to a constant) **cancels the ``log σ'(λ)`` change-of-variable Jacobian exactly**.
        So the two cancel and neither is written — the caller adds no Jacobian either.

        Until ``logP_r`` is fitted the RNA component is implicitly flat in ``log ρ_r``, which is exactly the
        ``logP_r ≡ 0`` reference. That is deliberate and *declared* — unlike the ``−log(1−f_g)`` this
        replaced, which was the RNA-side measure conversion for a prior that was never written, and which
        (with the Beta and the Jacobian) summed to an improper ``+0.5·λ`` ramp whose strength was set by the
        grid half-width.

        The grid top is the max density, so ``ρ_g = f_g·M/E ≤ M/E`` never extrapolates above it; the low
        side clamps to the depleted level (safe)."""
        eff = np.maximum(np.asarray(eff, dtype=np.float64), _EPS)
        mass = np.maximum(np.asarray(mass, dtype=np.float64), _EPS)
        fg = np.minimum(np.maximum(np.asarray(fg_grid, dtype=np.float64), _EPS), 1.0 - _EPS)  # (K,)
        log_me = np.log(mass) - np.log(eff)  # (n,) = log(M/E)
        log_rho_g = np.log(fg)[None, :] + log_me[:, None]  # (n,K)
        return np.interp(
            log_rho_g.ravel(), self.log_rho, self.logP, left=self.logP[0], right=self.logP[-1]
        ).reshape(log_rho_g.shape)
