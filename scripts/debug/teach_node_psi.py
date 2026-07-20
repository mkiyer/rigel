"""TEACHING artifact — walk the per-node ψ solve for ONE flagship node, on its real data, and show the fix.

ψ(f_g) = strand + gDNA-arm + RNA-arm + messages  (simplex_logodds._local_loglik_logodds).
For an unstranded node with depleted/gagged boundaries, strand≈flat and messages≈0, so ψ ≈ gDNA-arm + RNA-arm.

  gDNA arm = logP_g(ρ_g),  ρ_g = f_g·M/E     — the FITTED total-density NPMLE (informative; has modes)
  RNA arm  = ½·log(1−f_g)                     — the UNFITTED Jeffreys REFERENCE (logP_r never written)

The plot shows, over the f_g grid: the gDNA arm (its depleted + enriched modes), the RNA Jeffreys reference,
the current ψ (their sum → argmax = the CRUSH), and ψ with a DEMONSTRATION logP_r fitted on the (oracle-)
deconvolved RNA density in place of the Jeffreys reference → argmax moves to the truth. The demo logP_r uses the
oracle RNA split ONLY to prove the mechanism; ESTIMATING logP_r non-circularly is the jigsaw (refit ordering).

    OMP_NUM_THREADS=1 python teach_node_psi.py [--node 1909] [--out DIR]
"""

from __future__ import annotations

import argparse
import dataclasses
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
from flagship_interrogate import _oracle_per_node  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.npmle import DensityNPMLE  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_LN10 = np.log(10.0)
_EPS = 1e-12


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--condition", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--node", type=int, default=1909)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    inp = _scan_and_truth(suite, a.condition, index, cfg, work, cache)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

    cc0 = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc0, _debug=dbg)
    chain = dbg["chain"]
    cap = dbg["capture"]
    prior = dbg["gdna_prior"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)

    i = a.node
    M, E = mass[i], eff[i]
    rtot = M / max(E, _EPS)
    fo = G[i] / max(G[i] + R[i], _EPS)
    grid = np.asarray(cap["solve_grid"], float)          # f_g grid (K,)
    gdna_arm = np.asarray(cap["global_lp"], float)[i]    # logP_g(ρ_g) on the grid, incl. floor (K,)
    rna_ref = 0.5 * np.log(np.clip(1 - grid, _EPS, 1.0))  # the Jeffreys reference ½log(1−f_g)

    # DEMONSTRATION logP_r: fit the SAME NPMLE machinery on the (oracle-)deconvolved RNA density R/E across all
    # live nodes, then read it at ρ_r = (1−f_g)·M/E. Oracle split used ONLY to prove the mechanism.
    live = (eff > 1e-9 * 1.001) & (mass > _EPS)
    rna_prior = DensityNPMLE.fit(R[live], eff[live], bandwidth=cfg.calibration.npmle_bandwidth)
    log_rho_r = np.log(np.clip(1 - grid, _EPS, 1.0)) + np.log(max(rtot, _EPS))
    rna_arm_fit = np.interp(log_rho_r, rna_prior.log_rho, rna_prior.logP,
                            left=rna_prior.logP[0], right=rna_prior.logP[-1])

    interior = (grid > 0.005) & (grid < 0.995)  # exclude the two edge grid-quantization/clamp artifacts

    def _argmax(psi):
        pj = np.where(interior, psi, -np.inf)
        return grid[int(np.argmax(pj))]

    def _median(psi):
        post = np.exp(psi - psi.max())
        post /= max(post.sum(), _EPS)
        return grid[int(np.clip((np.cumsum(post) < 0.5).sum(), 0, grid.size - 1))]

    psi_cur = gdna_arm + rna_ref            # what production solves (strand flat, msgs ~0 here)
    psi_fix = gdna_arm + rna_arm_fit        # the fix: fitted logP_r in place of the Jeffreys reference

    # ---- printed step-by-step ----
    print(f"=== node {i}  {a.condition} ===")
    print(f"M={M:.0f} E={E:.0f}  ρtot=M/E={rtot:.2f} (log10 {np.log10(rtot):.2f})   TRUTH f_g={fo:.3f}")
    modes = (np.asarray(prior.log_rho) / _LN10)[
        np.r_[False, (prior.logP[1:-1] > prior.logP[:-2]) & (prior.logP[1:-1] >= prior.logP[2:]), False]
    ]
    print(f"gDNA prior modes (log10 ρ): {np.round(modes, 2)}")
    print("  f_g at which ρ_g hits each mode: " +
          "  ".join(f"{10**m/rtot:.4f}" for m in modes if 10**m <= rtot))
    print(f"{'f_g':>7} {'ρ_g':>8} {'ρ_r':>8} {'gDNAarm':>8} {'RNAref':>8} {'ψ_cur':>8} "
          f"{'logP_r':>8} {'ψ_fix':>8}")
    for f in (0.002, 0.01, 0.1, 0.5, 0.9, 0.98):
        j = int(np.argmin(np.abs(grid - f)))
        print(f"{grid[j]:>7.3f} {f*rtot:>8.3f} {(1-f)*rtot:>8.2f} {gdna_arm[j]:>8.2f} {rna_ref[j]:>8.2f} "
              f"{psi_cur[j]:>8.2f} {rna_arm_fit[j]:>8.2f} {psi_fix[j]:>8.2f}")
    print(f"\ninterior ARGMAX f_g:  current={_argmax(psi_cur):.3f}  FIXED(gDNA+logP_r)={_argmax(psi_fix):.3f}"
          f"   truth={fo:.3f}")
    print(f"posterior MEDIAN f_g (the solver read-out): current={_median(psi_cur):.3f}  "
          f"FIXED={_median(psi_fix):.3f}   truth={fo:.3f}")
    # how much does the background mode out-weigh the enriched mode in the gDNA arm alone?
    jb = int(np.argmin(np.abs(grid - 0.002))); je = int(np.argmin(np.abs(grid - 0.9)))
    print(f"gDNA arm: background reading (f_g≈.002)={gdna_arm[jb]:.2f}  vs enriched (f_g≈.9)={gdna_arm[je]:.2f}"
          f"  → background taller by {gdna_arm[jb]-gdna_arm[je]:+.2f} nats")

    # ---- plot ----
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 5.5))
    xr = np.asarray(prior.log_rho) / _LN10
    axL.plot(xr, np.exp(prior.logP - prior.logP.max()), "b-", label="gDNA prior P(ρ) [total-density NPMLE]")
    axL.axvline(np.log10(rtot), color="k", ls="--", label=f"this node ρtot={rtot:.0f}")
    axL.axvline(np.log10(max(G[i]/max(E,_EPS), 1e-9)), color="tab:red", label="TRUTH ρ_g")
    for m in modes:
        axL.axvline(m, color="0.7", lw=0.8)
    axL.set_xlabel("log10 ρ"); axL.set_ylabel("P(ρ) (norm.)")
    axL.set_title("The gDNA prior: depleted (background) mode + enriched mode.\nBoth are gDNA-consistent readings for this node.")
    axL.legend(fontsize=8)

    # Plot in the solver's own coordinate λ = logit(f_g), so BOTH modes (background at λ≈−6, enriched at
    # λ≈+2) are visible and spread out — on a linear f_g axis the decisive background peak is squished at 0.
    lam = np.log(grid / (1.0 - grid))
    axR.plot(lam, gdna_arm - gdna_arm.max(), color="tab:blue",
             label="gDNA arm  logP_g(ρ_g)  [FITTED, informative]")
    axR.plot(lam, rna_ref - rna_ref.max(), color="tab:orange",
             label="RNA arm = ½log(1−f_g)  [Jeffreys ref, UNFITTED]")
    axR.plot(lam, psi_cur - psi_cur.max(), color="k", lw=2,
             label=f"ψ current → median f_g {_median(psi_cur):.3f} (CRUSH)")
    axR.axvline(np.log(fo / (1 - fo)), color="tab:red", label=f"truth f_g={fo:.2f}")
    axR.axvline(0.0, color="0.8", lw=0.8)
    for xt, lab in ((-6.4, "f_g=.002"), (2.2, "f_g=.9")):
        axR.annotate(lab, (xt, 0.1), fontsize=7, ha="center", color="0.4")
    axR.set_xlabel("λ = logit(f_g)   (← more RNA        more DNA →)")
    axR.set_ylabel("ψ (shifted to max 0)")
    axR.set_ylim(-6, 0.6)
    axR.set_title("ψ in the solver's λ coordinate. The gDNA arm's BACKGROUND peak (λ≈−6, f_g≈0)\n"
                  "out-weighs its enriched peak (λ≈+2, f_g≈0.9) by ~2 nats → the solve sits at f_g≈0 (CRUSH).")
    axR.legend(fontsize=8, loc="lower center")
    fig.tight_layout()
    out = outdir / f"teach_node_{i}_psi.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
