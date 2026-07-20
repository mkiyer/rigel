"""The starting point: the total-density NPMLE (convolved) vs the ORACLE gDNA and ORACLE RNA NPMLEs
(deconvolved), superimposed on one common log-rate axis.

Three curves per condition, all fit with the PRODUCTION NPMLE (`npmle` internals: Poisson-lognormal
cell likelihood + fixed-kernel EM), on a COMMON log-ρ grid, all on the SAME gDNA effective length so the
decomposition is exactly additive per node:  ρ_total = M/E = (G+R)/E = ρ_g + ρ_r.

  * TOTAL  : g_hat = M   (the belief-free pass-0 substrate — what production fits at pass-0)
  * gDNA   : g_hat = G   (oracle true gDNA unspliced mass per node)
  * RNA    : g_hat = R   (oracle true RNA unspliced mass = mrna+nrna unspliced, per node)

ORACLE = the production accumulator split by TRUE fragment origin (read-name), sum-to-full validated as a hard
gate in `OracleTruth._validate`. Conservation M == G+R per node is asserted here before plotting (the premise).

    OMP_NUM_THREADS=1 python npmle_deconv_plot.py [--out DIR]
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
    """Per chain-node ORACLE gDNA (G) and RNA (R) unspliced mass — matching transfer_variance_diag/_true_fg."""

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
    G = np.where(isr, gR[ri], gB[bi])
    R = np.where(isr, rR[ri], rB[bi])
    return G, R


def _fit_on_grid(g_hat, eff, log_rho, bandwidth=0.15):
    """Production NPMLE core (npmle.DensityNPMLE.fit) on a GIVEN grid, τ=0 (point oracle masses). Returns the
    normalised grid density P(log ρ)."""
    var_g = np.zeros_like(g_hat)
    gc, ec, tc, wc = _collapse(g_hat, eff, var_g, dlog=0.1, dt=0.25)
    logL = _cell_loglik(gc, ec, tc, log_rho, n_gh=7)
    kk = _kernel_matrix(log_rho, bandwidth * _LN10)
    ln = np.exp(logL - logL.max(axis=1, keepdims=True))
    logL_sm = np.log(np.maximum(ln @ kk, _EPS))
    w = _em_weights(logL_sm, np.log(wc)[:, None], 150, 1e-5)
    dens = kk @ w
    return dens / max(float(dens.sum()), _EPS)


def analyse(inp, index):
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    M, eff = node_global_geometry(chain, geom)
    e = np.maximum(eff, _EPS)
    G, R = _oracle_gr(inp, chain)

    # --- CONSERVATION CHECK: the plot premise. M must equal G+R per node. ---
    live = e > 1e-9 * 1.001  # a face that could observe something
    resid = np.abs(M - (G + R))[live]
    rel = resid / np.maximum(M[live], 1.0)
    cons = dict(
        max_abs=float(resid.max()), max_rel=float(rel.max()),
        frac_off=float((rel > 1e-3).mean()), n=int(live.sum()),
        Mtot=float(M[live].sum()), Gtot=float(G[live].sum()), Rtot=float(R[live].sum()),
    )

    # common grid across the three densities' support + zero-anchor room
    dens_all = np.concatenate([
        (M[live] / e[live]), (G[live] / e[live]), (R[live] / e[live])
    ])
    dpos = dens_all[dens_all > 0]
    lo = float(np.log(np.percentile(dpos, 0.2))) - 1.0 * _LN10  # zero-anchor room
    hi = float(np.log(dpos.max())) + 0.4 * _LN10
    log_rho = np.linspace(lo, hi, 320)

    Mg, Gg, Rg = (v[live] for v in (M, G, R))
    eg = e[live]
    P_tot = _fit_on_grid(Mg, eg, log_rho)
    P_gdna = _fit_on_grid(Gg, eg, log_rho)
    P_rna = _fit_on_grid(Rg, eg, log_rho)
    return log_rho, P_tot, P_gdna, P_rna, cons


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = sorted(p.stem for p in cache.glob("*.pkl"))

    n = len(conds)
    ncol = 2
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 3.4 * nrow))
    axes = np.asarray(axes).ravel()
    for ax, c in zip(axes, conds):
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        log_rho, P_tot, P_gdna, P_rna, cons = analyse(inp, index)
        x = log_rho / _LN10  # log10 ρ
        ax.fill_between(x, 0, P_tot, color="0.75", alpha=0.5, label="TOTAL (M/E) — pass-0 substrate")
        ax.plot(x, P_tot, color="0.35", lw=1.4)
        ax.plot(x, P_gdna, color="tab:blue", lw=2.0, label="oracle gDNA (G/E)")
        ax.plot(x, P_rna, color="tab:red", lw=2.0, label="oracle RNA (R/E)")
        cap = "capON" if "capture_on" in c else "capOFF"
        nas = "+nascent" if "nrna_present" in c else "no-nascent"
        ss = "ss0.99" if "0.99" in c else "ss0.50"
        fg = cons["Gtot"] / max(cons["Gtot"] + cons["Rtot"], 1.0)
        ax.set_title(
            f"{cap} · {nas} · {ss}    [true gDNA frac={fg:.2f}]\n"
            f"conservation |M−(G+R)|: max_abs={cons['max_abs']:.2g}, "
            f"max_rel={cons['max_rel']:.1e}, off>1e-3={cons['frac_off']:.1%}",
            fontsize=8.5,
        )
        ax.set_xlabel("log₁₀ density ρ = mass / eff_gdna", fontsize=8)
        ax.set_ylabel("P(log ρ)  [node-weighted]", fontsize=8)
        ax.legend(fontsize=7, loc="upper left")
        ax.tick_params(labelsize=7)
        print(
            f"{c:52s} gDNAfrac={fg:.2f}  conservation max_rel={cons['max_rel']:.2e} "
            f"frac_off={cons['frac_off']:.3%}  (n={cons['n']})",
            flush=True,
        )
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(
        "NPMLE density landscapes — TOTAL (convolved starting point) vs ORACLE gDNA / RNA (deconvolved). "
        "Same eff_gdna axis ⇒ ρ_total = ρ_g + ρ_r per node.",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out = outdir / "npmle_deconv_landscapes.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
