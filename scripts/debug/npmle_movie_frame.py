"""The fit MOVIE — pass-0 total NPMLE, then the SOLVED gDNA/RNA NPMLEs after one full forward-backward solve,
each overlaid on the ORACLE gDNA/RNA NPMLE. "How does the deconvolution look vs oracle, frame by frame?"

Per condition:
  * fit P_total on the belief-free total density (pass-0 — the enrichment profile);
  * run the real solve (DensityNPMLE + node_sweep, projection σ²_transfer ON) → belief f_g per node;
  * fit the SOLVED gDNA NPMLE on ρ_g=f_g·M/E and the SOLVED RNA NPMLE on ρ_r=(1−f_g)·M/E;
  * fit the ORACLE gDNA/RNA NPMLE on G/E and R/E;
  * overlay all on a common eff_gdna axis (so ρ_total = ρ_g+ρ_r additively).

The gap between the SOLVED (dashed) and ORACLE (solid) curves is the deconvolution error — the next frame's
refit fits P_g/P_r on the solved belief, so this shows what that refit would be fed.

    OMP_NUM_THREADS=1 python npmle_movie_frame.py [--conditions a,b] [--out DIR]
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

from rigel.calibration.bp_solver import REGION, node_global_geometry, node_sweep
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.npmle import (
    DensityNPMLE,
    _cell_loglik,
    _collapse,
    _em_weights,
    _kernel_matrix,
)
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
)
from rigel.calibration.node_chain import build_node_chain
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics, init_beliefs
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12
_LN10 = np.log(10.0)


def _fit_grid(g_hat, eff, log_rho, h=0.15):
    var_g = np.zeros_like(g_hat)
    gc, ec, tc, wc = _collapse(g_hat, eff, var_g, dlog=0.1, dt=0.25)
    logL = _cell_loglik(gc, ec, tc, log_rho, n_gh=7)
    kk = _kernel_matrix(log_rho, h * _LN10)
    ln = np.exp(logL - logL.max(axis=1, keepdims=True))
    w = _em_weights(np.log(np.maximum(ln @ kk, _EPS)), np.log(wc)[:, None], 150, 1e-5)
    dens = kk @ w
    return dens / max(float(dens.sum()), _EPS)


def solve_belief(inp, index, cfg):
    """Same setup as sigma2_transfer_ab.solve but returns the belief + geometry + oracle pools."""
    cc = cfg.calibration
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    gfl, rfl = np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"])
    reg_el = region_eff_length(ra.region_size_bp, gfl)
    bnd_el = boundary_side_eff_length(gfl, ra.region_size_bp)
    fl_mean = boundary_eff_length(gfl)
    kappa = float(fit_strand_balance(inp["strand_model"]).rna_sense_frac)
    nd = node_gdna_density(sub, ra, reg_el, fl_mean)
    od_g = fit_gdna_strand_from_substrate(
        sub, ra, nd, bnd_el, rna_sense_frac=kappa
    ).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(
        sub, rna_sense_frac=kappa
    ).rna_strand_overdispersion
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = build_node_geometry(chain, sub, bsub, ra, gfl, rfl)
    st = build_node_statics(chain, sub, bsub, ra)
    b0 = init_beliefs(
        chain, sub, bsub, ra, rna_sense_frac=kappa, gdna_strand_overdispersion=od_g,
        rna_strand_overdispersion=od_r, n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand, logodds_window=cc.sweep_logodds_window, statics=st,
    )
    mass_g, eff_g = node_global_geometry(chain, geom)
    prior = DensityNPMLE.fit(mass_g, eff_g, bandwidth=cc.npmle_bandwidth)
    cap = {}
    belief = node_sweep(
        chain, st, geom, b0, ra, rna_sense_frac=kappa, gdna_strand_overdispersion=od_g,
        rna_strand_overdispersion=od_r, n_grid=cc.sweep_n_grid, logodds_window=cc.sweep_logodds_window,
        n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand, gdna_prior=prior,
        _capture=cap,
    )
    solve_belief.last = dict(prior=prior, cap=cap, mass_g=mass_g, eff_g=eff_g)
    # oracle G/R per chain node
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION
    pr_, pb_ = inp["region_pools"], inp["boundary_pools"]

    def pools(p):
        g = np.asarray(p["gdna_pos"], float) + np.asarray(p["gdna_neg"], float)
        r = sum(np.asarray(p[k], float) for k in ("mat_uns_pos", "mat_uns_neg", "nas_uns_pos", "nas_uns_neg"))
        return g, r

    gR, rR = pools(pr_)
    gB, rB = pools(pb_)
    ri, bi = np.clip(idx, 0, gR.shape[0] - 1), np.clip(idx, 0, gB.shape[0] - 1)
    G = np.where(isr, gR[ri], gB[bi])
    R = np.where(isr, rR[ri], rB[bi])
    fg = np.asarray(belief.f_g, float)
    return mass_g, eff_g, fg, G, R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))

    n = len(conds)
    ncol = 2
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 3.5 * nrow))
    axes = np.asarray(axes).ravel()
    for ax, c in zip(axes, conds):
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        M, eff, fg, G, R = solve_belief(inp, index, cfg)
        e = np.maximum(eff, _EPS)
        live = e > 1e-9 * 1.001
        M, e, fg, G, R = M[live], e[live], fg[live], G[live], R[live]
        gS, rS = fg * M, (1.0 - fg) * M  # SOLVED gDNA / RNA mass
        dens = np.concatenate([M / e, G / e, R / e, gS / e, rS / e])
        dpos = dens[dens > 0]
        lo = float(np.log(np.percentile(dpos, 0.2))) - 1.0 * _LN10
        hi = float(np.log(dpos.max())) + 0.4 * _LN10
        log_rho = np.linspace(lo, hi, 320)
        x = log_rho / _LN10
        P_tot = _fit_grid(M, e, log_rho)
        P_g_or, P_r_or = _fit_grid(G, e, log_rho), _fit_grid(R, e, log_rho)
        P_g_sv, P_r_sv = _fit_grid(gS, e, log_rho), _fit_grid(rS, e, log_rho)
        ax.fill_between(x, 0, P_tot, color="0.82", label="total (pass-0)")
        ax.plot(x, P_g_or, color="tab:blue", lw=2.0, label="gDNA oracle")
        ax.plot(x, P_g_sv, color="tab:blue", lw=1.6, ls="--", label="gDNA SOLVED")
        ax.plot(x, P_r_or, color="tab:red", lw=2.0, label="RNA oracle")
        ax.plot(x, P_r_sv, color="tab:red", lw=1.6, ls="--", label="RNA SOLVED")
        cap = "capON" if "capture_on" in c else "capOFF"
        nas = "+nas" if "nrna_present" in c else ""
        ss = "ss0.99" if "0.99" in c else "ss0.50"
        fgt = G.sum() / max(G.sum() + R.sum(), 1.0)
        ax.set_title(f"{cap} {nas} {ss}  [true gDNAfrac={fgt:.2f}]", fontsize=9)
        ax.set_xlabel("log₁₀ ρ = mass/eff_gdna", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper left")
        ax.tick_params(labelsize=7)
        print(f"{c:52s} plotted (solve done)", flush=True)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(
        "Fit movie frame 1 — SOLVED gDNA/RNA NPMLE (dashed) vs ORACLE (solid), after one full solve "
        "(projection σ²_transfer ON). Gap = deconvolution error the refit would be fed.",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = outdir / "npmle_movie_frame.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
