"""``DensityNPMLE`` — a Fixed-Kernel Poisson(-lognormal) Mixture NPMLE in count/rate space.

**ONE density-estimation engine, TWO distinct calibration uses.** A node's count is ``g ~ Poisson(ρ·E)`` for a
latent rate ``ρ ~ P(ρ)``; :meth:`~DensityNPMLE.fit` estimates the population ``P(log ρ)`` from per-node
``(count, eff-length)`` observations. Which density it is fit on, and which projection is read, defines the use:

1. **Enrichment NPMLE → σ²_transfer** (:meth:`~DensityNPMLE.project`). Fit on every node's **TOTAL** unspliced
   density BEFORE the pass. It models the hybrid-capture ENRICHMENT/DEPLETION landscape, **not** composition
   (a total-density model is composition-vacuous — count-zero-information). Its projection gives each node's
   ``(mu_proj, var_proj)``, from which the sweep builds the belief-free message transfer variance σ²_transfer
  This is the only NPMLE running in the INITIAL solve; it never
   votes composition. (`calibrate.enrichment_prior`.)

2. **gDNA hyperprior → ψ** (:meth:`~DensityNPMLE.logprior`). Fit AFTER the pass-0 solve on the **DECONVOLVED
   gDNA** density (``g_hat = f_g·count`` at belief width ``τ = √Var(log f_g)``) at the solvable nodes, with the
   aggregate background ``ρ_bg`` pinned as a smooth low-density component. Its projection onto each node's
   ``f_g`` axis is the composition (gDNA) arm of ψ for the REFIT solve. ANCHORED, EXTREMELY WEAK
   (``n_eff ≈ 0.15`` pseudo-obs — the strand likelihood + messages dominate). (`calibrate.gdna_prior`.)

Design (the two roles):

* **Fixed-Kernel Poisson Mixture** — ``P(log ρ) = Σ_j w_j·N(log ρ; log ρ_j, h²)``, a mixture of fixed-width
  (bandwidth ``h``) Gaussian kernels on a log-ρ grid; the weights are fit by plain EM (monotone,
  arithmetic-only, deterministic — NO spline / λ / matrix-inversion, so none of the GCV-λ nondeterminism).
  The fixed kernel forbids any peak sharper than ``h`` ⇒ smooth, never a bed-of-nails (the raw NPMLE's
  Kiefer-Wolfowitz atomicity). ``h`` is the familiar KDE bandwidth knob (`config.npmle_bandwidth`).
* **Count/rate space** — the Poisson is zero-native (``count=0`` ⇒ ``e^{−ρE}`` ⇒ mass at ρ→0), so the sparse,
  zero-inflated data needs no boundary correction; the per-node likelihood width is set by the physics
  (broad for low-count/short-E — the discreteness comb smooths itself away). The pass-0 enrichment fit uses
  ``var_g=None`` (τ=0, pure Poisson at ``f_g=1``); the refit hyperprior passes the solved belief width.

Performance: the per-node likelihood depends only on the ``(count, log-E, τ)`` cell, so ~10⁶ nodes collapse to
~10³–10⁴ weighted cells → the fit is sub-second and needs no C++/numba.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from .background_reference import BackgroundReference

__all__ = ["DensityNPMLE"]

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


def _kde_density(log_rho, a_cells, w_cells, h, *, a_floor=-np.inf, w_floor=0.0, h_floor=None):
    """Additive, occupancy-weighted, fixed-bandwidth KDE on the log-rate grid (the Role-B representation —
    ``Σ_c w_c·N(grid; a_c, h) + w_floor·N(grid; a_floor, h_floor)`` → a
    normalized density on ``log_rho`` (G,).

    Each training cell deposits its OWN fixed-width kernel at its (resolution-floored) deconvolved-gDNA point
    estimate ``a_c``, weighted by OCCUPANCY ``w_c`` (node count) — there is NO EM competition, so a minority
    (e.g. capture-enriched) population can never be competed away, and with a COMMON bandwidth occupancy equals
    height (no per-node τ discounting). The **weak floor** is one pseudo-observation (``w_floor``) at the derived
    background level ``a_floor`` — the anchor a sub-resolution / zero-gDNA node lands on, never an ``n_regions``
    tower. ``h_floor`` defaults to ``h`` (used when ``σ_bg`` is infinite, i.e. Σg=0)."""
    inv = 1.0 / (h * np.sqrt(2.0 * np.pi))
    z = (log_rho[:, None] - a_cells[None, :]) / h  # (G, n_cell)
    dens = (w_cells[None, :] * np.exp(-0.5 * z * z)).sum(axis=1) * inv
    if np.isfinite(a_floor) and w_floor > 0.0:
        hf = float(h if h_floor is None else h_floor)
        zf = (log_rho - a_floor) / hf
        dens = dens + w_floor * np.exp(-0.5 * zf * zf) / (hf * np.sqrt(2.0 * np.pi))
    return dens / max(float(dens.sum()), _EPS)


@dataclass(frozen=True, slots=True)
class DensityNPMLE:
    """A fitted population rate density ``P(log ρ)`` — the Fixed-Kernel Poisson-mixture NPMLE, pre-evaluated
    (smoothed + interior-floored) on a log-ρ grid. Read via :meth:`project` (the enrichment → σ²_transfer use)
    or :meth:`logprior` (the gDNA-hyperprior → ψ use); see the module docstring for the two-use design."""

    log_rho: np.ndarray  #: (G,) log-rate grid = the mixture component locations μ_j
    logP: (
        np.ndarray
    )  #: (G,) log P(log ρ) — kernel-rendered, uniform-floored inside the support, real tails
    weights: (
        np.ndarray
    )  #: (G,) mixture component weights w_j (Σ w_j·N(μ_j, h²)) — for :meth:`project`
    bandwidth: float  #: h (decades), the kernel width
    n_cells: int  #: collapsed-cell count (diagnostic)
    #: The aggregate DNA-background reference (`background_reference.measure_background`) — the ONE-SIDED
    #: log-floor `:meth:`logprior`` applies. ``-inf``/``+inf`` (the default) ⇒ dormant, so an unmeasured
    #: background leaves the prior EXACTLY as before (safe default; the A/B gate + the wire-in are separate).
    log_rho_bg: float = (
        -np.inf
    )  #: natural-log background rate; `-inf` ⇒ no floor (DNA-free / fully depleted)
    sigma_bg: float = np.inf  #: Poisson softness of the floor (√Var(log ρ_bg)); `+inf` ⇒ no floor

    @classmethod
    def fit(
        cls,
        g_hat: np.ndarray,
        eff: np.ndarray,
        var_g: np.ndarray | None = None,
        *,
        background: BackgroundReference | None = None,
        bandwidth: float = 0.15,
        n_grid: int = 200,
        em_iters: int = 150,
        em_tol: float = 1.0e-5,
        log_dlog: float = 0.1,
        tau_dt: float = 0.25,
        n_gh: int = 7,
        floor_eps: float = 0.02,
        rho_floor: float | None = None,
        additive: bool = False,
    ) -> "DensityNPMLE":
        """Fit the Fixed-Kernel Poisson-lognormal Mixture on the per-node gDNA ``(count ĝ, eff-length E,
        belief width τ²=var_g)``. ``var_g=None`` ⇒ ``τ=0`` (pass-0, pure Poisson at ``f_g=1``, ``ĝ=count``).

        The grid spans the observed density support: the top is the MAX density (so the projection never
        extrapolates ``ρ_g=f_g·M/E`` above the grid — avoiding the clamped-tail false-positive), the bottom
        is below the bulk (zero-anchor room; the low side clamps safely to the depleted level).

        ``additive=True`` builds the **Role-B representation** instead of the
        EM mixture: an **occupancy-weighted, fixed-bandwidth KDE** on the deconvolved-gDNA point estimates
        (each cell one kernel, weighted by node count — no EM competition, occupancy=height) plus a **weak
        1-pseudo-observation floor** at ``background.log_rho_floor`` (NOT the ``n_regions`` aggregate cell).
        The output object (``log_rho``/``logP``/``weights``/``project``/``logprior``) is interface-identical;
        ``weights`` carries the per-grid density mass so ``project``'s fallback stays defined. ``additive=False``
        (default) is the EM NPMLE, **byte-identical** to before — Role A (σ²_transfer) uses it unchanged."""
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
        # THE AGGREGATE BACKGROUND CELL — the pooled intergenic[/intron] regions injected as ONE Poisson
        # observation at the DERIVED floor density ``ρ_floor`` (``background.log_rho_floor``) over the genome-scale
        # ΣE. This is the hybrid the design calls for (B; gdna_background_floor_derivation
        # md): the pooled rate ``Σg/ΣE`` Fisher-blended with the per-region resolution wall ``1/harmmean(E_zero)``
        # — NOT the old ``e^{−ρΣE}`` which drove ρ to ``1/ΣE`` (~3 logs too low, the confident-FP seed). The cell's
        # huge ΣE makes it a SHARP low mode AT ρ_floor; it carries the population weight (``n_regions``). The pooled
        # regions are excluded from ``g_hat`` by the caller, so nothing is double-counted.
        use_bg = background is not None and background.eff_total > _EPS and background.n_regions > 0
        # Grid bottom = the derived floor (bounded, honest) − 3h of kernel room; the aggregate cell anchors the
        # zero / sub-floor background there. (Old fallback ran to 1/ΣE — the ~3-log collapse.)
        if rho_floor is not None and rho_floor > 0.0:
            lo = float(np.log(rho_floor))
        elif use_bg and np.isfinite(background.log_rho_floor):
            lo = min(lo, float(background.log_rho_floor) - 3.0 * h)
        log_rho = np.linspace(lo, hi, int(n_grid))

        gc, ec, tc, wc = _collapse(g_hat, eff, var_g, dlog=log_dlog, dt=tau_dt)
        if additive:
            # ── THE ROLE-B ADDITIVE KDE — occupancy kernels + a weak floor, NO EM. ──
            # Each cell's point estimate is its deconvolved-gDNA density, floored at the 1-count resolution wall
            # (a zero-ĝ / pure-RNA node reads "gDNA ≤ 1/E", not −∞). Occupancy ``wc`` is the weight; the
            # background is a SEPARATE 1-pseudo-observation floor at ρ_floor — never the ``n_regions`` cell.
            a_cells = np.log(np.maximum(gc, 1.0)) - np.log(ec)
            if use_bg and np.isfinite(background.log_rho_floor):
                a_floor = float(background.log_rho_floor)
                sig = float(background.sigma_bg)
                h_floor = h if not np.isfinite(sig) else max(sig, h)
                w_floor = 1.0  # ONE pseudo-observation (the KDE cap) — weak vs Σ wc = |T|
            else:
                a_floor, h_floor, w_floor = -np.inf, h, 0.0
            # THE GRID MUST SPAN THE FLOORED KERNEL CENTRES (spec §1: hi=max(â)+2h), NOT the raw density ld.
            # When every deconvolved count gc<1 (the common mostly-RNA case: ĝ=f_g·M) the floored centres sit
            # at 1/E — far ABOVE the raw-density grid — so a raw-ld grid would leave every kernel off-grid,
            # underflowing dens_grid to ~0 (a phantom high-gDNA mode + an unclamped-blo IndexError). Rebuild
            # the grid from a_cells (+ the floor) here; the EM path above keeps its raw-ld grid unchanged.
            hi_a = float(np.max(a_cells)) + 2.0 * h
            lo_a = float(np.min(a_cells)) - 3.0 * h
            if np.isfinite(a_floor):
                lo_a = min(lo_a, a_floor - 3.0 * float(h_floor))
            log_rho = np.linspace(lo_a, hi_a, int(n_grid))
            dens_grid = _kde_density(
                log_rho, a_cells, wc, h, a_floor=a_floor, w_floor=w_floor, h_floor=h_floor
            )
            w = dens_grid  # per-grid mass — keeps ``project``'s fallback well-defined (OI-9)
        else:
            if use_bg:
                # the aggregate cell's count places its density at ρ_floor exactly: gc/ec = exp(log_rho_floor).
                # (n0=0 ⇒ log_rho_floor = log(Σg/ΣE) ⇒ gc = Σg = n_counts, byte-identical to the old case.)
                bg_count = (
                    float(np.exp(background.log_rho_floor) * background.eff_total)
                    if np.isfinite(background.log_rho_floor)
                    else float(background.n_counts)
                )
                gc = np.append(gc, bg_count)
                ec = np.append(ec, float(background.eff_total))
                tc = np.append(tc, 0.0)
                wc = np.append(wc, float(background.n_regions))
            logL = _cell_loglik(gc, ec, tc, log_rho, n_gh=int(n_gh))
            kk = _kernel_matrix(log_rho, h)
            # convolve each cell likelihood with the fixed kernel, EM the weights, render the mixture density.
            ln = np.exp(logL - logL.max(axis=1, keepdims=True))
            logL_sm = np.log(np.maximum(ln @ kk, _EPS))
            w_full = _em_weights(logL_sm, np.log(wc)[:, None], em_iters, em_tol)
            dens_grid = (
                kk @ w_full
            )  # (G,) smooth mixture density (the aggregate cell shaped its low mode)
            dens_grid = dens_grid / max(float(dens_grid.sum()), _EPS)
            w = w_full  # the aggregate is a DATA cell, not an extra component ⇒ no column to strip

        # interior uniform floor (fills −∞ valleys between modes) bounded to the fitted support [q0.5, q99.5];
        # OUTSIDE the support the real Gaussian tails are left intact (the false-positive suppression).
        cdf = np.cumsum(w) / max(float(w.sum()), _EPS)
        blo = log_rho[min(int(np.searchsorted(cdf, 0.005)), n_grid - 1)]
        bhi = log_rho[min(int(np.searchsorted(cdf, 0.995)), n_grid - 1)]
        in_supp = (log_rho >= blo) & (log_rho <= bhi)
        step = float(log_rho[1] - log_rho[0])
        uni = np.where(in_supp, floor_eps / max(float(bhi - blo), _EPS) * step, 0.0)
        logP = np.log(
            np.maximum(np.where(in_supp, (1.0 - floor_eps) * dens_grid, dens_grid) + uni, _EPS)
        )
        return cls(
            log_rho=log_rho,
            logP=logP,
            weights=w,
            bandwidth=float(bandwidth),
            n_cells=int(gc.shape[0]),
            log_rho_bg=(-np.inf if background is None else float(background.log_rho_bg)),
            sigma_bg=(np.inf if background is None else float(background.sigma_bg)),
        )

    def logprior(self, fg_grid, mass, eff):
        """Project onto the ``f_g`` solve grid → ``(n_nodes, K)`` additive term = ``log P(log ρ_g)`` at
        ``ρ_g = f_g·M/E``. **Bare** — no reference prior, no measure term, no Jacobian.

        This is the ``ψ_λ = strand + logP_g + logP_r`` result. The
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
        side clamps to the depleted level.

        **The background enters SMOOTHLY, never as a clamp** (B). When this prior was
        fit with a background, the aggregate ``ρ_bg`` is a **pinned Gaussian component of the mixture** (see
        :meth:`fit`) — it fills the ``[0, 1/E]`` low-density vacuum the NPMLE cannot resolve, as ordinary
        (already-normalized) prior mass in ``logP``. There is **no** one-sided floor / half-Gaussian wall: a
        clamp is a cliff, and a Bayesian prior must be smooth with honest strength (so enriched / CNV-amplified
        DNA can overcome it) — the retired floor is in the archive."""
        eff = np.maximum(np.asarray(eff, dtype=np.float64), _EPS)
        mass = np.maximum(np.asarray(mass, dtype=np.float64), _EPS)
        fg = np.minimum(np.maximum(np.asarray(fg_grid, dtype=np.float64), _EPS), 1.0 - _EPS)  # (K,)
        log_me = np.log(mass) - np.log(eff)  # (n,) = log(M/E)
        log_rho_g = np.log(fg)[None, :] + log_me[:, None]  # (n,K)
        return np.interp(
            log_rho_g.ravel(), self.log_rho, self.logP, left=self.logP[0], right=self.logP[-1]
        ).reshape(log_rho_g.shape)

    def project(self, mass, eff, *, chunk: int = 8192):
        """Belief-free projection of each node's total log-density ``d = log(mass/eff)`` onto the mixture
        ``Σ_j w_j·N(μ_j, h²)`` → the latent log-rate's posterior ``(mu_proj, var_proj)`` (design doc
         role 2). The observed density is soft-assigned to the mixture
        components by responsibilities ``r_j ∝ w_j·N(d; μ_j, h²)``; then::

            mu_proj  = Σ_j r_j·μ_j                              (denoised rate — snaps toward a mode)
            var_proj = Σ_j r_j·(μ_j − mu_proj)²  +  h²          (between-component ambiguity + within-mode floor)

        In a mode the responsibilities concentrate ⇒ ``var_proj ≈ h²`` (the count-zero-information max-precision
        floor); in the valley between two modes they split ⇒ the between-component term (the mode gap²) spikes.
        COUNT-FREE — the enrichment landscape only; the count enters the message precision solely as ``1/M_src``.
        The per-edge transfer variance is ``var_proj[dst] + (mu_proj[dst] − mu_proj[src])²`` (F1). Chunked for
        genome scale (the ``(n, G)`` responsibility matrix never materialises whole)."""
        eff = np.maximum(np.asarray(eff, dtype=np.float64), _EPS)
        mass = np.maximum(np.asarray(mass, dtype=np.float64), _EPS)
        d = np.log(mass) - np.log(eff)  # (n,) belief-free total log-density
        h = float(self.bandwidth) * np.log(
            10.0
        )  # bandwidth in natural-log units (μ_j = log_rho is ln)
        logw = np.log(np.maximum(self.weights, _EPS))  # (G,) — the pre-kernel component weights
        mu_out = np.empty_like(d)
        var_out = np.empty_like(d)
        for lo in range(0, d.shape[0], int(chunk)):
            hi = min(lo + int(chunk), d.shape[0])
            z = (d[lo:hi, None] - self.log_rho[None, :]) / h  # (c,G)
            lr = logw[None, :] - 0.5 * z * z
            lr -= lr.max(axis=1, keepdims=True)
            r = np.exp(lr)
            r /= np.maximum(r.sum(axis=1, keepdims=True), _EPS)
            mu = (r * self.log_rho[None, :]).sum(1)
            mu_out[lo:hi] = mu
            var_out[lo:hi] = (r * (self.log_rho[None, :] - mu[:, None]) ** 2).sum(1) + h * h
        return mu_out, var_out
