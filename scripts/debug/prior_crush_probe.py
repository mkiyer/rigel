"""PRIOR-CRUSH PROBE — why does the pass-0 prior drive a truly-gDNA (fo≈1) node to f_g=0 in the LOCAL solve?

For a handful of the worst flagship nodes (truly ~gDNA, solved to f_g≈0), print the per-node evidence over the
f_g grid: the global NPMLE logprior term, and where the node's density sits vs the prior's modes. A "weak"
prior (n_eff≈0.15) must NOT be able to move f_g from 0.49 (strand-neutral) to 0.00 — if it does, either the
prior is far stronger than advertised for enriched nodes, or the projection is wrong.

    OMP_NUM_THREADS=1 python prior_crush_probe.py [--nodes 1909,2337,...]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402
from refit_loop_study import _setup  # noqa: E402

from rigel.calibration.bp_solver import node_sweep
from rigel.calibration.npmle import DensityNPMLE
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12
_LN10 = np.log(10.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--condition", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--bandwidth", type=float, default=0.15)
    ap.add_argument("--nodes", default="1909,2337,1055,1523,2999")
    ap.add_argument("--no-floor", action="store_true", help="background=None (pre-P2/P3 regime)")
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    cc = cfg.calibration
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    inp = _scan_and_truth(suite, a.condition, index, cfg, work, cache)
    s = _setup(inp, index, cfg)
    mass_g, eff_g = s["mass_g"], s["eff_g"]
    bg = None if a.no_floor else s["bg"]
    prior = DensityNPMLE.fit(mass_g, eff_g, background=bg, bandwidth=a.bandwidth)

    cap = {}
    node_sweep(
        s["chain"], s["st"], s["geom"], s["b0"], s["ra"], rna_sense_frac=s["kappa"],
        gdna_strand_overdispersion=s["od_g"], rna_strand_overdispersion=s["od_r"], n_grid=cc.sweep_n_grid,
        logodds_window=cc.sweep_logodds_window, n_tilt=cc.sweep_n_tilt,
        n_grid_ss=cc.sweep_n_grid_single_strand, gdna_prior=prior, _capture=cap,
    )
    grid = np.asarray(cap["solve_grid"], float)  # (K,) f_g lattice
    glp = np.asarray(cap["global_lp"], float)  # (n,K) global prior term on the grid
    lam_grid = np.log(np.clip(grid, _EPS, 1 - _EPS) / np.clip(1 - grid, _EPS, 1 - _EPS))

    def _logprior_decay(fraction, M, E, slope=0.5):
        """logP_total at ρ = fraction·M/E with a Jeffreys-consistent DECAYING left tail (slope 0.5 in
        natural-log-rate space), instead of the production clamp `left=logP[0]`."""
        log_rho = np.log(np.clip(fraction, _EPS, 1.0)) + (np.log(M) - np.log(max(E, _EPS)))
        x, P = prior.log_rho, prior.logP
        out = np.interp(log_rho, x, P, right=P[-1])
        left = log_rho < x[0]
        out[left] = P[0] + slope * (log_rho[left] - x[0])
        return out

    def _report(nid, name, psi):
        w = np.exp(psi - psi.max())
        w /= max(w.sum(), _EPS)
        mean = float((w * grid).sum())
        return f"{name:22s} argmax={grid[int(np.argmax(psi))]:.3f} mean={mean:.3f}"

    # prior mode locations (log10 ρ) for reference
    P = np.exp(prior.logP - prior.logP.max())
    ismode = np.r_[False, (P[1:-1] > P[:-2]) & (P[1:-1] >= P[2:]), False]
    modes = prior.log_rho[ismode] / _LN10
    print(f"=== {a.condition}  bw={a.bandwidth} ===")
    print(f"prior modes (log10 ρ): {np.round(modes,2)}   grid support log10ρ "
          f"[{prior.log_rho[0]/_LN10:.2f}, {prior.log_rho[-1]/_LN10:.2f}]\n")

    for nid in [int(x) for x in a.nodes.split(",")]:
        M, E = mass_g[nid], eff_g[nid]
        rho_tot = M / max(E, _EPS)
        lp = glp[nid]  # (K,) global prior term over f_g grid
        # normalize to a probability over the grid to read the prior's implied f_g posterior + strength
        w = np.exp(lp - lp.max())
        w = w / max(w.sum(), _EPS)
        lam = np.log(np.clip(grid, _EPS, 1 - _EPS) / np.clip(1 - grid, _EPS, 1 - _EPS))
        m1 = float((w * lam).sum())
        v1 = float((w * (lam - m1) ** 2).sum())
        neff = (1.0 / max(v1, 1e-9)) / 0.25
        # prior term at a few f_g anchors
        def at(fg):
            j = int(np.argmin(np.abs(grid - fg)))
            return lp[j]
        # the RNA Jeffreys REFERENCE arm ½·log(1−f_g) that the solve ADDS alongside the gDNA prior — the only
        # thing bounding the f_g→1 vertex today. Compare its pull against the truth to the gDNA prior's.
        rna_ref = 0.5 * np.log(np.clip(1 - grid, _EPS, 1.0))

        def rat(fg):
            j = int(np.argmin(np.abs(grid - fg)))
            return rna_ref[j]
        print(f"node {nid}: M={M:.0f} E={E:.0f}  log10 ρ_tot={np.log10(max(rho_tot,1e-12)):.2f}  "
              f"(f_g=1 ⇒ ρ_g=ρ_tot; f_g→0 ⇒ ρ_g→0)")
        print(f"   prior-only implied f_g: mode≈σ({m1:.1f})={1/(1+np.exp(-m1)):.3f}  n_eff={neff:.2f}")
        print(f"   gDNA prior term:  f_g=0.01:{at(0.01):+.2f}  0.1:{at(0.1):+.2f}  0.5:{at(0.5):+.2f}  "
              f"0.9:{at(0.9):+.2f}  0.99:{at(0.99):+.2f}   (favor high f_g by {at(0.9)-at(0.1):+.2f})")
        print(f"   RNA ref ½log(1−f_g): 0.1:{rat(0.1):+.2f}  0.5:{rat(0.5):+.2f}  0.9:{rat(0.9):+.2f}  "
              f"0.98:{rat(0.98):+.2f}   (penalize high f_g by {rat(0.9)-rat(0.1):+.2f})")
        # NET ψ (gDNA prior + RNA ref, strand flat): where is the max, and how does truth compare to f_g≈0.1?
        psi = lp + rna_ref
        jmax = int(np.argmax(psi))
        print(f"   NET ψ (prior+RNAref) argmax f_g={grid[jmax]:.3f}  |  ψ(0.98)−ψ(0.1)="
              f"{psi[int(np.argmin(np.abs(grid-0.98)))]-psi[int(np.argmin(np.abs(grid-0.1)))]:+.2f} nats "
              f"(negative ⇒ truth LOSES to f_g≈0.1)")
        # CANDIDATE ψ FORMS (prior-only; strand ≈ flat for these unstranded nodes). Truth in the header.
        gdna_ref = 0.5 * np.log(np.clip(grid, _EPS, 1.0))
        lp_g_decay = _logprior_decay(grid, M, E)  # gDNA arm, total NPMLE, decaying tail
        lp_r_decay = _logprior_decay(1.0 - grid, M, E)  # RNA arm, SAME total NPMLE at ρ_r, decaying tail
        eps = 0.05
        variants = {
            "V0 current(clamp+Jeff)": lp + rna_ref,
            "V1 decayTail+Jeff": lp_g_decay + rna_ref,
            "V2 SYMM total-NPMLE": lp_g_decay + lp_r_decay,
            "V3 Jeffreys guardrail": lp + rna_ref + eps * (gdna_ref + 0.5 * np.log(np.clip(1 - grid, _EPS, 1.0))),
            "V4 restore ½log(f_g)": lp + rna_ref + gdna_ref,
        }
        for nm, psi in variants.items():
            print("   " + _report(nid, nm, psi))


if __name__ == "__main__":
    main()
