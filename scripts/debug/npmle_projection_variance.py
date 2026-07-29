"""NPMLE PROJECTION VARIANCE — observed density -> (mode, VARIANCE) via the mixture, + the RNA bandwidth diag.

Answers three questions:
  (1) Why is the RNA fit wiggly? -> bandwidth. Panel A refits at h=0.15/0.30/0.50 dex.
  (2) The observed-density -> NPMLE projection exists (`_prior_strength`/`projection_curve`). This builds the
      clean RATE-space version returning MODE and VARIANCE.
  (3) The mixture-projection VARIANCE: treat the latent rate as one of the mixture components; the observed
      density gives RESPONSIBILITIES r_j ∝ w_j·N(logρ_obs; μ_j, h²) (optionally × Poisson(k|μ_j,E)). Then
          mode = μ_argmax,   mean = Σ r_j μ_j,   VAR = Σ r_j (μ_j-mean)² + h²   (between-mode + within-mode)
      In a MODE the responsibilities concentrate -> VAR ≈ h² (small). In the VALLEY between two modes they
      split -> the between-mode term (the mode gap²) DOMINATES -> VAR spikes. That spike is the honest
      "this density is ambiguous" signal. COUNT-FREE (mixture only) or COUNT-AWARE (× Poisson).

    OMP_NUM_THREADS=1 python npmle_projection_variance.py [--cond NAME] [--out DIR]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.special import gammaln

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import node_global_geometry
from rigel.calibration.npmle import _cell_loglik, _collapse, _em_weights, _kernel_matrix
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import build_node_geometry
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12
_LN10 = np.log(10.0)


def _oracle_gr(inp, chain):
    def pools(p):
        g = np.asarray(p["gdna_pos"], float) + np.asarray(p["gdna_neg"], float)
        r = sum(np.asarray(p[k], float) for k in ("mat_uns_pos", "mat_uns_neg", "nas_uns_pos", "nas_uns_neg"))
        return g, r

    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION
    gR, rR = pools(inp["region_pools"])
    gB, rB = pools(inp["boundary_pools"])
    ri, bi = np.clip(idx, 0, gR.shape[0] - 1), np.clip(idx, 0, gB.shape[0] - 1)
    return np.where(isr, gR[ri], gB[bi]), np.where(isr, rR[ri], rB[bi])


def fit_weights(g_hat, eff, log_rho, bandwidth):
    """Production NPMLE core -> the mixture WEIGHTS w_j over the grid components (μ_j = log_rho[j], width h)."""
    var_g = np.zeros_like(g_hat)
    gc, ec, tc, wc = _collapse(g_hat, eff, var_g, dlog=0.1, dt=0.25)
    logL = _cell_loglik(gc, ec, tc, log_rho, n_gh=7)
    kk = _kernel_matrix(log_rho, bandwidth * _LN10)
    ln = np.exp(logL - logL.max(axis=1, keepdims=True))
    logL_sm = np.log(np.maximum(ln @ kk, _EPS))
    w = _em_weights(logL_sm, np.log(wc)[:, None], 150, 1e-5)
    return w, kk @ w  # weights, rendered density


def project(log_rho, w, h, log_rho_obs, k=None, eff=None):
    """Project an observed log-density onto the mixture. Returns (mode, mean, var) of the latent log-rate.
    COUNT-FREE: responsibilities r_j ∝ w_j·N(logρ_obs; μ_j, h²). COUNT-AWARE (k,eff given): × Poisson(k|μ_j·E).
    var = between-mode Σr_j(μ_j-mean)²  +  within-mode h²."""
    mu = log_rho  # component locations
    lr = np.log(np.maximum(w, _EPS)) - 0.5 * ((log_rho_obs - mu) / h) ** 2
    if k is not None:
        lam = np.exp(mu) * eff
        lr = lr + k * (mu + np.log(eff)) - lam - gammaln(k + 1.0)
    r = np.exp(lr - lr.max())
    r /= r.sum()
    mean = float((r * mu).sum())
    var_between = float((r * (mu - mean) ** 2).sum())
    return mu[int(np.argmax(lr))], mean, var_between + h * h


def build(inp, index, bandwidth=0.15):
    ra = RegionArrays.from_index(index)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    M, eff = node_global_geometry(chain, geom)
    e = np.maximum(eff, _EPS)
    G, R = _oracle_gr(inp, chain)
    live = e > 1e-9 * 1.001
    return M[live], G[live], R[live], e[live]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--cond", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    inp = _scan_and_truth(suite, a.cond, index, cfg, work, suite / "_selfsolve_cache")
    M, G, R, e = build(inp, index)

    dens_all = np.concatenate([M / e, G / e, R / e])
    dpos = dens_all[dens_all > 0]
    lo = float(np.log(np.percentile(dpos, 0.2))) - 1.0 * _LN10
    hi = float(np.log(dpos.max())) + 0.4 * _LN10
    log_rho = np.linspace(lo, hi, 320)
    x = log_rho / _LN10

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2))

    # --- Panel A: RNA bandwidth ---
    for h, col in ((0.15, "tab:red"), (0.30, "tab:orange"), (0.50, "tab:green")):
        _w, dens = fit_weights(R, e, log_rho, h)
        axes[0].plot(x, dens, color=col, lw=1.8, label=f"RNA  h={h} dex")
    axes[0].set_title(f"(A) RNA fit vs bandwidth h  [{a.cond.split('_gdna300_')[1]}]\n"
                      "h=0.15 (production) is wiggly on the sparse/broad RNA; larger h smooths", fontsize=9)
    axes[0].set_xlabel("log₁₀ ρ"); axes[0].set_ylabel("P(log ρ)"); axes[0].legend(fontsize=8)

    # --- Panel B: total-density NPMLE (bimodal) + projection VARIANCE vs observed density ---
    h = 0.15
    w_tot, dens_tot = fit_weights(M, e, log_rho, h)
    axes[1].fill_between(x, 0, dens_tot / dens_tot.max(), color="0.8", label="total NPMLE (scaled)")
    obs_grid = np.linspace(lo + 0.3, hi - 0.3, 200)
    var_free = np.array([project(log_rho, w_tot, h * _LN10, o)[2] for o in obs_grid])
    var_k30 = np.array([project(log_rho, w_tot, h * _LN10, o, k=30.0, eff=float(np.median(e)))[2] for o in obs_grid])
    var_k1 = np.array([project(log_rho, w_tot, h * _LN10, o, k=1.0, eff=float(np.median(e)))[2] for o in obs_grid])
    axes[1].plot(obs_grid / _LN10, var_free, "k-", lw=2.2, label="VAR count-free (mixture only)")
    axes[1].plot(obs_grid / _LN10, var_k30, "b--", lw=1.6, label="VAR count-aware k=30")
    axes[1].plot(obs_grid / _LN10, var_k1, "r:", lw=1.6, label="VAR count-aware k=1")
    axes[1].set_title("(B) projection VARIANCE of the latent log-rate vs observed density\n"
                      "count-free spikes in the VALLEY between modes (ambiguous); flat-low inside a mode",
                      fontsize=9)
    axes[1].set_xlabel("log₁₀ observed density"); axes[1].set_ylabel("projected Var(log ρ)")
    axes[1].legend(fontsize=8)

    # --- Panel C: source-only vs pair (does the crossing need the destination?) ---
    modes = log_rho[np.r_[False, (dens_tot[1:-1] > dens_tot[:-2]) & (dens_tot[1:-1] > dens_tot[2:]), False]]
    modes = modes[np.interp(modes, log_rho, dens_tot) > 0.05 * dens_tot.max()]
    dep = modes.min() if modes.size else lo + 1
    enr = modes.max() if modes.size else hi - 1
    val = 0.5 * (dep + enr)  # the valley midpoint
    print(f"\n[{a.cond}]  detected modes (log10 ρ): {np.round(modes / _LN10, 2)}  "
          f"dep={dep / _LN10:.2f} enr={enr / _LN10:.2f} valley≈{val / _LN10:.2f}")
    print("\nprojection VARIANCE (count-free) at representative densities:")
    for name, o in (("deep-depleted", dep), ("valley/ambiguous", val), ("deep-enriched", enr)):
        md, mn, vr = project(log_rho, w_tot, h * _LN10, o)
        print(f"   {name:18s} obs log10ρ={o / _LN10:+.2f}  ->  mode={md / _LN10:+.2f}  Var={vr:.3f}")
    print("\nSOURCE-only vs PAIR (mode assignment) for the three edge types:")
    for name, (os_, od_) in (
        ("dep->dep", (dep, dep)), ("enr->enr", (enr, enr)), ("dep<->enr CROSS", (dep, enr))
    ):
        ms = project(log_rho, w_tot, h * _LN10, os_)[0]
        mdd = project(log_rho, w_tot, h * _LN10, od_)[0]
        src_var = project(log_rho, w_tot, h * _LN10, os_)[2]
        pair_gap = (mdd - ms) ** 2
        print(f"   {name:16s} src Var(proj)={src_var:.3f}   mode gap²(pair)={pair_gap:.3f}   "
              f"-> source-only misses the crossing (gap² is 0 for src alone)")
    # Panel C: bar of the three
    labels = ["dep->dep", "enr->enr", "dep<->enr"]
    srcv = [project(log_rho, w_tot, h * _LN10, o)[2] for o in (dep, enr, dep)]
    pairv = [0.0, 0.0, (project(log_rho, w_tot, h * _LN10, enr)[0] - project(log_rho, w_tot, h * _LN10, dep)[0]) ** 2]
    xb = np.arange(3)
    axes[2].bar(xb - 0.2, srcv, 0.4, label="source proj Var (mode ambiguity)", color="tab:blue")
    axes[2].bar(xb + 0.2, pairv, 0.4, label="pair mode-gap² (crossing)", color="tab:red")
    axes[2].set_xticks(xb); axes[2].set_xticklabels(labels, fontsize=8)
    axes[2].set_title("(C) source projection variance vs the pair crossing term\n"
                      "the crossing (dep<->enr) is invisible to the source alone", fontsize=9)
    axes[2].set_ylabel("variance (nats²)"); axes[2].legend(fontsize=8)

    fig.tight_layout()
    out = outdir / "npmle_projection_variance.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
