"""Fusion NPMLE prototype — the belief-width Poisson-rate ``P(ρ_g)``, watched across passes 0→1→2.

The crux the roadmap gates N3 on (``docs/calibration/npmle_roadmap.md`` §4): fit the gDNA-rate distribution
``P(ρ)`` on ALL nodes using the CURRENT belief, with honest precision entering as observation WIDTH so the
prior STARTS WEAK and SHARPENS only as beliefs converge — never pinning a node prematurely, always keeping
the count-0 zero anchor.

Per node the observation is the believed gDNA count ``ĝ = f_g·k`` over the gDNA exposure ``E`` (``eff_global``)
with belief width ``τ² = Var(log f_g)``. The per-node likelihood is the **Poisson-lognormal**: the gDNA count
is ``g ~ Poisson(ρ·E)`` and the belief places ``log g ~ N(log ĝ, τ²)``; marginalise g by Gauss-Hermite
quadrature. ONE likelihood, no arm-seam:

  * ``ĝ=0`` (every count-0 / believed-pure-RNA node) → ``g_q=0`` → ``e^{−ρE}`` — the exact zero anchor.
  * uncertain low-count node (large τ, e.g. the ĝ≈1 boundary comb at τ≈1.5) → the lognormal smears the count
    over a wide g → a BROAD likelihood → cannot form a false sharp spike (this is the honest-width firewall,
    now active in the low-count regime — the fix over the point-estimate Poisson).
  * high ĝ → the Poisson sharpens, the lognormal width dominates → Gaussian-in-log-rate (the seam is automatic).

Nodes with no information (``Var(log f_g)=∞`` — AMBIG at init, empty nodes) are excluded (the honest-width
firewall for free).

PERFORMANCE. Two decouplings make this fast enough to iterate:
  * the slow part — reproducing ``calibrate``'s genome-scale sweep to get the pass 0/1/2 beliefs — is run ONCE
    and cached (compact per-node arrays; the fit then re-runs in milliseconds);
  * the fit itself uses the ``(ĝ-bin, log-E-bin)`` / ``(log-rate-bin, log-σ-bin)`` **cell collapse** — the EM
    likelihood depends only on the cell, so ~10⁶ nodes collapse to ~10³–10⁴ weighted cells (roadmap §1).

We fit at pass 0 (init belief), pass 1 (solved, no KDE), pass 2 (solved, production KDE) and plot the ``P(ρ)``
trajectory. The projection read (observed density → posterior f_g) uses a SMOOTHED + UNIFORM-FLOORED ``P(ρ)``
(the analogue of the KDE mixture-bridge — no −∞ valleys), directly testing the "jagged interp" concern.

    python scripts/debug/npmle_fusion.py --cache C.pkl[,...] --out DIR [--belief-cache DIR]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gammaln

from rigel.calibration.bp_solver import (
    adjacent_disagreement_variance,
    build_node_geometry,
    build_node_statics,
    init_beliefs,
    node_sweep,
)
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.gdna_density_prior import GdnaDensityPrior, build_training_substrate
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
    overdispersion_for_beta,
)
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import _node_region_type
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, nrna_active_strands
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig

sys.path.insert(0, str(Path(__file__).parent))
from calib_cache import load  # noqa: E402

_EPS = 1.0e-12


# ---------------------------------------------------------------------------
# Reproduce calibrate's sweep to get the belief at each pass (b0=init, b1, b2).
# Cached: the sweep is the slow (genome-scale) part; the fit iterates against the cache.
# ---------------------------------------------------------------------------
def _sj_boundary_mask(bsub, ra):
    """Per-boundary single-strand intron↔exon splice-junction mask (copied from pass0_kde_landscape to avoid
    importing its oracle deps). SJ on strand s ⇔ nrna-active on s on both flanks AND exon XOR intron across
    the seam; single-strand ⇔ exactly one of ±SJ."""
    sig = np.asarray(ra.signature).astype(np.int64)
    lr = np.asarray(bsub.left_region, np.int64)
    rr = np.asarray(bsub.right_region, np.int64)
    sig_l = np.where(lr >= 0, sig[np.clip(lr, 0, None)], 0)
    sig_r = np.where(rr >= 0, sig[np.clip(rr, 0, None)], 0)
    nrp_l, nrn_l = nrna_active_strands(sig_l)
    nrp_r, nrn_r = nrna_active_strands(sig_r)
    ex_p_l, ex_p_r = (sig_l & BIT_EXON_POS) != 0, (sig_r & BIT_EXON_POS) != 0
    ex_n_l, ex_n_r = (sig_l & BIT_EXON_NEG) != 0, (sig_r & BIT_EXON_NEG) != 0
    sj_pos = (nrp_l & nrp_r) & (ex_p_l ^ ex_p_r)
    sj_neg = (nrn_l & nrn_r) & (ex_n_l ^ ex_n_r)
    return sj_pos ^ sj_neg
def _run_sweep(inp, cfg: PipelineConfig):
    """Compact per-node arrays for the fit at each pass — the production sequence (`calibrate`):
    b0=init, b1=sweep(b0, no-KDE), b2=sweep(b1, production KDE). Returns {k, eff_global, e_max, fg, varg}."""
    cc = cfg.calibration
    payload, ra = inp["payload"], inp["region_arrays"]
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    bsub = BoundarySubstrate.from_payload(payload)
    gdna_fl_pmf = np.asarray(inp["gdna_fl_pmf"])
    rna_fl_pmf = np.asarray(inp["rna_fl_pmf"])

    region_eff_len = region_eff_length(ra.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)

    kappa = float(fit_strand_balance(inp["strand_model"]).rna_sense_frac)
    node_density_raw = node_gdna_density(substrate, ra, region_eff_len, fl_mean)
    od_g = fit_gdna_strand_from_substrate(
        substrate, ra, node_density_raw, boundary_eff_len, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(cc.gdna_strand_prior_alpha_beta),
        prior_weight=cc.gdna_strand_prior_weight,
    ).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(
        substrate, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(cc.rna_strand_prior_alpha_beta),
        prior_weight=cc.rna_strand_prior_weight,
    ).rna_strand_overdispersion

    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geom = build_node_geometry(chain, substrate, bsub, ra, gdna_fl_pmf, rna_fl_pmf)
    statics = build_node_statics(chain, substrate, bsub, ra)
    b0 = init_beliefs(
        chain, substrate, bsub, ra,
        rna_sense_frac=kappa, gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r,
        n_grid=cc.sweep_n_grid, n_grid_ss=cc.sweep_n_grid_single_strand,
        logodds_window=cc.sweep_logodds_window, statics=statics,
    )
    sig_total = adjacent_disagreement_variance(chain, geom)

    def _sweep(belief_in, prior):
        return node_sweep(
            chain, statics, geom, belief_in, ra, bsub,
            rna_sense_frac=kappa, gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r,
            n_grid=cc.sweep_n_grid, logodds_window=cc.sweep_logodds_window, n_tilt=cc.sweep_n_tilt,
            n_grid_ss=cc.sweep_n_grid_single_strand, gdna_prior=prior, disagreement_sigma2=sig_total,
        )

    b1 = _sweep(b0, None)
    train_sub = build_training_substrate(chain, b1, geom, statics, ra, bsub, min_eff_length=fl_mean)
    prior = (
        GdnaDensityPrior.fit(
            train_sub, bandwidth=cc.gdna_prior_bandwidth,
            mixture_bridge=cc.gdna_prior_mixture_bridge, bridge_trim_pct=cc.calib_kde_bridge_trim_pct,
        )
        if train_sub.n >= cc.calib_kde_min_training_nodes
        else None
    )
    b2 = _sweep(b1, prior) if prior is not None else b1

    is_reg = np.asarray(chain.kind) == REGION
    egl, egr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)
    eff_global = np.maximum(np.where(is_reg, egl, 0.5 * (egl + egr)), _EPS)
    k = np.asarray(statics.u_pos, float) + np.asarray(statics.u_neg, float)
    # region TYPE (0 intergenic / 1 intron / 2 exon; boundaries → -1) for substrate slicing.
    rtype_node, _ = _node_region_type(chain, ra)
    rtype = np.where(is_reg, rtype_node, -1).astype(np.int64)
    # DUAL-STRUCTURAL enriched sample: single-strand SJ boundaries with spliced=0 (RNA-free) — the
    # mature-excluded, on-target gDNA window. Belief-independent (structural). Mapped onto the chain.
    ref = np.asarray(chain.ref_idx, np.int64)
    B = np.asarray(bsub.left_region).shape[0]
    bi = np.clip(ref, 0, B - 1)
    spl_b = np.asarray(bsub.left.mass_spliced, float) + np.asarray(bsub.right.mass_spliced, float)
    sj0_b = _sj_boundary_mask(bsub, ra) & (spl_b <= 1.0e-9)
    dual_enriched = (~is_reg) & sj0_b[bi]
    return dict(
        k=k, eff_global=eff_global, e_max=float(eff_global.max()),
        is_reg=is_reg, rtype=rtype, dual_enriched=dual_enriched,
        fg=np.stack([np.asarray(b.f_g, float) for b in (b0, b1, b2)]),
        varg=np.stack([np.asarray(b.var_gdna, float) for b in (b0, b1, b2)]),
    )


def belief_cache(cache_path, cfg, cache_dir: Path | None):
    """Load the compact per-pass belief arrays, running (and caching) the slow sweep only on a miss."""
    cond = Path(cache_path).stem
    npz = (cache_dir / f"{cond}.beliefs.npz") if cache_dir else None
    if npz is not None and npz.exists():
        d = np.load(npz)
        return {kk: d[kk] for kk in d.files} | {"e_max": float(d["e_max"])}, 0.0
    t0 = time.time()
    bc = _run_sweep(load(cache_path), cfg)
    dt = time.time() - t0
    if npz is not None:
        npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez(npz, **{kk: np.asarray(v) for kk, v in bc.items()})
    return bc, dt


def node_obs(bc, p, substrate="all"):
    """Per-node gDNA observation at pass p: (ĝ, E, τ²=Var(log f_g), fit_mask). fit_mask = finite Var(log f_g)
    ∩ the SUBSTRATE. τ² is the BELIEF width; the Poisson supplies the count noise (no +1/ĝ — no double-count).

    substrate: ``all`` = every solved/locked node; ``region`` = region nodes only (drops short-E boundaries,
    the discreteness comb); ``structural`` = intergenic+intron regions only (belief-independent zero anchor,
    no exons); ``dual`` = the DUAL-STRUCTURAL substrate: structural depleted (intergenic+intron regions) +
    structural enriched (RNA-free single-strand SJ boundaries). Belief-INDEPENDENT: ĝ=k (f_g=1 asserted by
    structure — depleted regions are gDNA-clean, SJ crossings are mature-excluded) and τ=0, so the fit is a
    fixed prior, identical across passes (no GIGO)."""
    var_g = bc["varg"][p]  # Var(log f_g); ∞ ⇒ no information
    k, E = bc["k"], bc["eff_global"]
    if substrate == "total":
        # EVERY node at f_g=1 (total unspliced density as if gDNA), count precision only (τ=0). Shows the
        # TOTAL-density landscape: depleted (no-RNA regions) + enriched (high-RNA regions) — but the enriched
        # mode is RNA-DOMINATED, so this is NOT the gDNA prior (it is what the belief/solve must PEEL).
        return k.copy(), E, np.zeros_like(k), (E > _EPS)
    if substrate == "dual":
        struct = bc["is_reg"] & ((bc["rtype"] == 0) | (bc["rtype"] == 1))
        fit = (struct | bc["dual_enriched"]) & (E > _EPS)
        return k.copy(), E, np.zeros_like(k), fit  # f_g=1 (ĝ=k), τ=0 — structural, belief-independent
    g_hat = bc["fg"][p] * k
    fit = np.isfinite(var_g) & (E > _EPS)
    if substrate == "region":
        fit &= bc["is_reg"]
    elif substrate == "structural":
        fit &= bc["is_reg"] & ((bc["rtype"] == 0) | (bc["rtype"] == 1))
    elif substrate != "all":
        raise ValueError(f"unknown substrate {substrate!r}")
    return g_hat, E, np.maximum(var_g, 0.0), fit


# ---------------------------------------------------------------------------
# The fusion NPMLE fit — Poisson-lognormal likelihood + (ĝ, log-E, τ) cell collapse.
# ---------------------------------------------------------------------------
def _lse(a, axis):
    """Fast log-sum-exp (manual; scipy's has per-call overhead that dominated the EM loop)."""
    mx = np.max(a, axis=axis, keepdims=True)
    mx = np.where(np.isfinite(mx), mx, 0.0)
    return (mx + np.log(np.sum(np.exp(a - mx), axis=axis, keepdims=True))).squeeze(axis)


def _collapse(g_hat, E, var_g, *, dlog=0.1, dt=0.25, g_eps=1e-6):
    """Collapse nodes to weighted cells whose likelihood is (near-)identical: key on ``(log-ĝ bin | zero,
    log-E bin, τ bin)``. ĝ≈0 nodes share one bin (their likelihood is ``e^{−ρE}`` independent of ĝ/τ). The
    bins trade memory/compute for binning error (a perf knob, not an accuracy constant). Representative =
    the in-cell mean. Returns ``(gc, Ec, t2c, wc)``."""
    logE = np.log(np.maximum(E, _EPS))
    tau = np.sqrt(np.maximum(var_g, 0.0))
    gi = np.where(
        g_hat < g_eps, np.int64(-(10**9)), np.round(np.log(np.maximum(g_hat, _EPS)) / dlog).astype(np.int64)
    )
    ei = np.round(logE / dlog).astype(np.int64)
    ti = np.round(tau / dt).astype(np.int64)
    cols = np.stack([gi, ei, ti], axis=1)
    _uniq, inv, wc = np.unique(cols, axis=0, return_inverse=True, return_counts=True)
    inv = np.asarray(inv).reshape(-1)
    wc = wc.astype(np.float64)
    mean = lambda v: np.bincount(inv, weights=v) / wc  # noqa: E731  in-cell mean (accurate representative)
    return mean(g_hat), mean(E), mean(np.maximum(var_g, 0.0)), wc


def _loglik_pln(g_hat, E, var_g, log_rho, *, n_gh=7, chunk=4000):
    """Per-cell log-likelihood over the log-ρ grid = the **Poisson-lognormal**: the gDNA count is
    ``g ~ Poisson(ρ·E)`` and the belief places ``log g ~ N(log ĝ, τ²)`` (τ²=Var(log f_g)); marginalise g by
    Gauss-Hermite quadrature. ĝ=0 ⇒ ``g_q=0`` ⇒ ``e^{−ρE}`` (the exact zero anchor, no special-casing). High ĝ
    ⇒ the Poisson sharpens and the lognormal width dominates ⇒ Gaussian-in-log-rate (the seam is automatic).
    Returns (n_cell, G)."""
    x, wq = np.polynomial.hermite.hermgauss(n_gh)
    lwq = np.log(wq) - 0.5 * np.log(np.pi)  # weight of ∫e^{-x²}·/√π (lognormal expectation)
    log_rho_E = log_rho[None, None, :]  # (1,1,G)
    n = g_hat.shape[0]
    out = np.empty((n, log_rho.shape[0]), dtype=np.float64)
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        tau = np.sqrt(var_g[lo:hi])[:, None]  # (c,1)
        gq = g_hat[lo:hi][:, None] * np.exp(np.sqrt(2.0) * tau * x[None, :])  # (c,Q) quadrature counts
        logE = np.log(E[lo:hi])[:, None, None]  # (c,1,1)
        gq3 = gq[:, :, None]  # (c,Q,1)
        lam = np.exp(log_rho_E + logE)  # (c? no -> 1,1,G broadcast) expected counts ρ·E  -> (c,1,G)
        # logpois(c,q,j) = g_q·(logρ_j+logE_c) − ρ_j E_c − lgamma(g_q+1)
        logpois = gq3 * (log_rho[None, None, :] + logE) - lam - gammaln(gq3 + 1.0)
        out[lo:hi] = _lse(lwq[None, :, None] + logpois, axis=1)  # marginalise g over the quadrature
    return out


def _kernel_matrix(log_rho, h):
    """``K[i,j] = N(log ρ_i; log ρ_j, h²)`` normalised so each component column integrates to 1 over the grid.
    Column j is the fixed-width kernel of mixture component j; ``K @ w`` renders the mixture density."""
    d = log_rho[:, None] - log_rho[None, :]
    K = np.exp(-0.5 * (d / h) ** 2)
    return K / np.maximum(K.sum(axis=0, keepdims=True), _EPS)


def _em_weights(logL, logwc, n_iter, tol):
    """Standard mixture EM for the component weights over collapsed cells (monotone, arithmetic-only,
    deterministic). ``logL`` = (n_cell, G) cell-vs-component log-likelihood; ``logwc`` = (n_cell,1) log counts."""
    G = logL.shape[1]
    logw = np.full(G, -np.log(G))
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


def fusion_fit(g_hat, E, var_g, log_rho, *, method="kernel", kernel_bw=0.15, n_iter=200, tol=1e-6):
    """Estimate ``P(ρ)`` over the log-ρ grid from the collapsed cells. Returns the grid DENSITY (normalised)
    + n_cells. Estimators:

    * ``kernel`` (the Fixed-Kernel Poisson Mixture, DEFAULT — reviewer's synthesis). Parameterise
      ``P(log ρ) = Σ_j w_j·N(log ρ; log ρ_j, h²)`` — a mixture of FIXED-width Gaussian kernels (``h``=
      ``kernel_bw`` decades) on the grid; fit the weights ``w_j`` by standard EM. The cell-vs-component
      likelihood is the Poisson-lognormal convolved with the kernel (``L_c ⊛ N_h``). EM DECONVOLVES (recovers
      the true modes — needed to separate depleted from enriched under capture) but the fixed kernel forbids
      any peak sharper than ``h`` ⇒ smooth, no bed-of-nails, no spline/λ/matrix-inversion — 100 % deterministic
      (EM is monotone, arithmetic-only). ``h`` is the familiar KDE bandwidth knob.
    * ``avg`` — the population average of the per-node likelihood curves (one EM step from uniform; the
      smooth-but-OVERSMOOTHED extreme — no deconvolution, so it cannot separate modes).
    * ``npmle`` — EM to convergence with δ-atoms (h=0); the sharpest deconvolution but ATOMIC/spiky
      (Kiefer-Wolfowitz). Contrast only.

    Returns ``(P_grid, n_cells)``."""
    gc, _Ec, _t2c, wc = _collapse(g_hat, E, var_g)
    logL = _loglik_pln(gc, _Ec, _t2c, log_rho)  # (n_cell, G) Poisson-lognormal
    logwc = np.log(wc)[:, None]
    if method == "avg":
        logw = _lse(logwc + (logL - _lse(logL, axis=1)[:, None]), axis=0)
        return np.exp(logw - _lse(logw, axis=0)), gc.shape[0]
    if method == "npmle":
        w = _em_weights(logL, logwc, n_iter, tol)
        return w, gc.shape[0]
    if method != "kernel":
        raise ValueError(f"unknown method {method!r}")
    # Fixed-kernel Poisson mixture: convolve each cell likelihood with the fixed kernel, EM the weights,
    # render the smooth mixture density. Per-cell max-shift keeps the exp well-scaled (cancels in the E-step).
    K = _kernel_matrix(log_rho, kernel_bw * np.log(10.0))  # (G,G)
    Ln = np.exp(logL - logL.max(axis=1, keepdims=True))  # (n_cell,G) rescaled likelihood in [0,1]
    logL_sm = np.log(np.maximum(Ln @ K, _EPS))  # (n_cell,G) = log(L_c ⊛ N_h) up to a per-cell constant
    w = _em_weights(logL_sm, logwc, n_iter, tol)
    P = K @ w  # render the mixture density on the grid (smooth, width ≥ h)
    return P / max(float(P.sum()), _EPS), gc.shape[0]


def smooth_logprior(w, log_rho, *, h_dex=0.15, eps=0.02):
    """Smoothed ``log P(ρ)`` on the grid with REAL tails — the projection read. Convolve the discrete NPMLE
    weights with a Gaussian of width ``h_dex`` decades (a resolution knob), giving genuine Gaussian decay
    ABOVE and BELOW the fitted mass (the false-positive suppression that a flat floor would destroy — the
    ``logpdf_kernel`` lesson). The uniform ``eps`` floor is applied ONLY inside the fitted support
    ``[q0.5%, q99.5%]`` (the analogue of the KDE mixture-bridge: fills interior −∞ valleys but leaves the
    tails' real decay intact)."""
    h = h_dex * np.log(10.0)
    kern = np.exp(-0.5 * ((log_rho[:, None] - log_rho[None, :]) / h) ** 2)
    kern /= kern.sum(axis=1, keepdims=True)
    dens = kern @ w
    dens = dens / max(float(dens.sum()), _EPS)
    cdf = np.cumsum(w) / max(float(w.sum()), _EPS)
    blo = log_rho[int(np.searchsorted(cdf, 0.005))]
    bhi = log_rho[min(int(np.searchsorted(cdf, 0.995)), log_rho.shape[0] - 1)]
    in_supp = (log_rho >= blo) & (log_rho <= bhi)
    uni = np.where(in_supp, eps / max(float(bhi - blo), _EPS) * float(log_rho[1] - log_rho[0]), 0.0)
    return np.log(np.maximum(np.where(in_supp, (1.0 - eps) * dens, dens) + uni, _EPS))


def projection_curve(log_rho, logP, ld_grid, fg_grid):
    """Observed log-density → posterior f_g under the (smoothed) prior + RNA Jeffreys. Returns (mean, p10, p90)."""
    fg = np.clip(fg_grid, _EPS, 1 - _EPS)
    jeff = -np.log1p(-fg)
    mean, p10, p90 = (np.empty_like(ld_grid) for _ in range(3))
    for i, ld in enumerate(ld_grid):
        psi = np.interp(np.log(fg) + ld, log_rho, logP, left=logP[0], right=logP[-1]) + jeff
        ww = np.exp(psi - psi.max())
        ww /= ww.sum()
        mean[i] = float((ww * fg).sum())
        cdf = np.cumsum(ww)
        p10[i] = fg[int(np.searchsorted(cdf, 0.10))]
        p90[i] = fg[int(np.searchsorted(cdf, 0.90))]
    return mean, p10, p90


def _summary(w, log_rho, e_max):
    """Weighted mean/std of log10 ρ (SPREAD = weak↔sharp) + zero-atom mass (ρ·E_max < 0.1)."""
    l10 = log_rho / np.log(10.0)
    mu = float((w * l10).sum())
    sd = float(np.sqrt(max((w * (l10 - mu) ** 2).sum(), 0.0)))
    pi0 = float(w[log_rho <= np.log(0.1 / e_max)].sum())
    return mu, sd, pi0


def run(cache_path, outdir: Path, cfg: PipelineConfig, cache_dir: Path | None, substrate: str, method: str):
    cond = Path(cache_path).stem
    bc, t_sweep = belief_cache(cache_path, cfg, cache_dir)
    e_max = float(bc["e_max"])

    # Grid range = where the DATA is (percentiles of the observed nonzero densities), NOT tied to the single
    # largest-E region (which drags the floor to a meaningless ~1e-8). A display/perf knob (D3), not accuracy:
    # 1 decade below the 0.5th pct gives the count-0 zero anchor room; 0.5 decade above the 99.5th pct.
    g2, E2, _s2, f2 = node_obs(bc, 2, substrate)
    nz = f2 & (g2 > 0)
    ldp = np.log(g2[nz] / E2[nz]) if nz.any() else np.array([np.log(1.0 / e_max)])
    ln10 = np.log(10.0)
    lo = float(np.percentile(ldp, 0.5)) - 1.0 * ln10
    hi = float(np.percentile(ldp, 99.5)) + 0.5 * ln10
    log_rho = np.linspace(lo, hi, 300)

    fits, stats = [], []
    t1 = time.time()
    for p in range(3):
        gh, E, sg, fit = node_obs(bc, p, substrate)
        w, ncell = fusion_fit(gh[fit], E[fit], sg[fit], log_rho, method=method)
        fits.append(w)
        stats.append((_summary(w, log_rho, e_max), int(fit.sum()), ncell))
    t_fit = time.time() - t1

    l10 = log_rho / np.log(10.0)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    colors = ["tab:orange", "tab:green", "tab:blue"]
    for p in range(3):
        (mu, sd, pi0), nfit, ncell = stats[p]
        axes[0].plot(
            l10, fits[p], color=colors[p], lw=1.8,
            label=f"pass {p}: n={nfit} cells={ncell} σ={sd:.2f} μ={mu:.2f} π₀={pi0:.2f}",
        )
    axes[0].set_title(f"{cond} [{substrate}/{method}]: P(ρ) across passes")
    axes[0].set_xlabel("log10 gDNA rate ρ")
    axes[0].set_ylabel("weight")
    axes[0].legend(fontsize=8)
    # the kernel/avg fit is ALREADY a smooth density → only fill grid gaps + floor (no tail-widening reconv).
    logP2 = smooth_logprior(fits[2], log_rho, h_dex=0.03)
    ld_grid = np.linspace(lo - 0.5, hi + 1.0, 140)
    # f_g on a LOG-ODDS grid (like the real σ(λ) solve grid, L=10) so the projection can reach the small f_g
    # that a high-density node needs to land in the low-ρ prior mass — a linear grid floored at 1e-3 cannot.
    fg_grid = 1.0 / (1.0 + np.exp(-np.linspace(-12.0, 12.0, 400)))
    fmean, fp10, fp90 = projection_curve(log_rho, logP2, ld_grid, fg_grid)
    xd = ld_grid / np.log(10.0)
    axes[1].fill_between(xd, fp10, fp90, color="red", alpha=0.15, label="10–90% posterior")
    axes[1].plot(xd, fmean, "r-", lw=2, label="posterior mean f_g")
    axes[1].axhline(0.5, color="gray", lw=0.5, ls=":")
    axes[1].set_ylim(-0.05, 1.05)
    axes[1].set_title("pass-2 projection: observed density → f_g (smoothed prior)")
    axes[1].set_xlabel("log10 observed density")
    axes[1].set_ylabel("prior-preferred f_g")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    out = outdir / f"fusion_{cond}_{substrate}_{method}.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)

    spreads = " -> ".join(f"{stats[p][0][1]:.2f}" for p in range(3))
    pi0s = " -> ".join(f"{stats[p][0][2]:.2f}" for p in range(3))
    print(
        f"{cond:10s}[{substrate:10s}/{method:5s}] sweep={t_sweep:.1f}s fit={t_fit:.3f}s cells~{stats[2][2]}  "
        f"spread(σ log10ρ) {spreads}  π₀ {pi0s}  -> {out}",
        flush=True,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--belief-cache", default=None, help="dir to cache the slow per-pass sweep beliefs")
    ap.add_argument("--substrate", default="all", help="fit set: all | region | structural (comma-list ok)")
    ap.add_argument("--method", default="avg", help="estimator: avg (KDE-analog) | npmle (comma-list ok)")
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.belief_cache) if args.belief_cache else None
    cfg = PipelineConfig()
    for c in args.cache.split(","):
        for sub in args.substrate.split(","):
            for meth in args.method.split(","):
                run(c, outdir, cfg, cache_dir, sub, meth)


if __name__ == "__main__":
    main()
