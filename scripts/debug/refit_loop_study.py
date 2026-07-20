"""REFIT-LOOP EXPLORATION — study (not score) the calibration refit iteration.

Production (`calibrate.py`) runs `calib_refit_iters` refits: pass-0 fits the gDNA-rate NPMLE on TOTAL density,
solves; each refit re-fits P_g on the SOLVED gDNA counts (g_hat=f_g·M, width var_gdna) and re-solves,
WARM-CONTINUING from the previous belief. Because `node_sweep` uses the one `gdna_prior` for BOTH the logprior
AND the projection σ²_transfer, a refit updates the enrichment landscape the messages are damped by, too.

This harness mirrors that loop EXACTLY and instruments every pass with several NON-redundant lenses (the point
is to STUDY the trajectory, not reduce it to one number):

  * ORACLE mwae        — mass-weighted |f_g_solved − f_g_oracle| (needs the synthetic truth). Overall.
  * DISAGREEMENT vs the ORACLE FLOOR — the ground-truth-FREE optimization objective the owner posed. The truth
    has a NON-ZERO between-node disagreement (real enrichment crossings + expression steps); the refit should
    drive the SOLVED disagreement TOWARD the oracle's, not to zero. Reported per component (gDNA, RNA) as the
    mean squared adjacent log-density gap, alongside the oracle's own value on the SAME edges. Solved < oracle
    ⇒ OVER-smoothing (the degenerate collapse); solved ≫ oracle ⇒ under-solved.
  * CONVERGENCE Δf_g   — mass-weighted mean |f_g^k − f_g^{k-1}|: does the loop settle or oscillate?
  * P_g MODES          — count + top locations of the refit landscape: does it drift or stabilise?

    OMP_NUM_THREADS=1 python refit_loop_study.py [--conditions a,b] [--bandwidths 0.15,0.25] [--passes 5]
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

from rigel.calibration.background_reference import measure_background
from rigel.calibration.bp_solver import REGION, node_global_geometry, node_sweep
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.npmle import DensityNPMLE
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
    overdispersion_for_beta,
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


def _setup(inp, index, cfg):
    """Everything the sweep needs + the per-node oracle gDNA/RNA mass, in chain-node order."""
    cc = cfg.calibration
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    gfl, rfl = np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"])
    reg_el = region_eff_length(ra.region_size_bp, gfl)
    # Production default background: intergenic-only, no robust trim (matches CalibrationConfig defaults).
    bg = measure_background(sub, ra, reg_el)
    bnd_el = boundary_side_eff_length(gfl, ra.region_size_bp)
    fl_mean = boundary_eff_length(gfl)
    kappa = float(fit_strand_balance(inp["strand_model"]).rna_sense_frac)
    nd = node_gdna_density(sub, ra, reg_el, fl_mean)
    od_g = fit_gdna_strand_from_substrate(
        sub, ra, nd, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(cc.gdna_strand_prior_alpha_beta),
        prior_weight=cc.gdna_strand_prior_weight,
    ).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(
        sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(cc.rna_strand_prior_alpha_beta),
        prior_weight=cc.rna_strand_prior_weight,
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
    return dict(
        chain=chain, st=st, geom=geom, b0=b0, ra=ra, bsub=bsub, kappa=kappa, od_g=od_g, od_r=od_r,
        mass_g=mass_g, eff_g=eff_g, G=G, R=R, cc=cc, bg=bg,
    )


def _sweep(s, belief, prior, transfer_variance=True):
    cc = s["cc"]
    return node_sweep(
        s["chain"], s["st"], s["geom"], belief, s["ra"], s["bsub"], rna_sense_frac=s["kappa"],
        gdna_strand_overdispersion=s["od_g"], rna_strand_overdispersion=s["od_r"], n_grid=cc.sweep_n_grid,
        logodds_window=cc.sweep_logodds_window, n_tilt=cc.sweep_n_tilt,
        n_grid_ss=cc.sweep_n_grid_single_strand, gdna_prior=prior, transfer_variance=transfer_variance,
    )


def _disagreement(u, left, edge_ok):
    """Mean squared adjacent log-density gap over valid edges (u = log component density per node)."""
    i = np.where((left >= 0) & edge_ok)[0]
    j = left[i]
    ok = edge_ok[i] & edge_ok[j]
    i, j = i[ok], j[ok]
    if i.size == 0:
        return 0.0, 0
    d = u[i] - u[j]
    return float(np.mean(d * d)), int(i.size)


def _measure(s, belief, prior, prev_fg):
    mass, eff, G, R = s["mass_g"], s["eff_g"], s["G"], s["R"]
    left = np.asarray(s["chain"].left, np.int64)
    fg = np.asarray(belief.f_g, float)
    live = (eff > 1e-9 * 1.001) & (mass > _EPS)
    # oracle f_g on live nodes with any oracle mass
    tot = G + R
    fo = np.where(tot > _EPS, G / np.maximum(tot, _EPS), 0.0)
    wm = np.where(live & (tot > _EPS), mass, 0.0)
    mwae = float((wm * np.abs(fg - fo)).sum() / max(wm.sum(), _EPS))
    # component log-densities (shared eff so solved & oracle are apples-to-apples on each edge)
    e = np.maximum(eff, _EPS)
    ug_s = np.log(np.maximum(fg * mass, _EPS)) - np.log(e)
    ur_s = np.log(np.maximum((1.0 - fg) * mass, _EPS)) - np.log(e)
    ug_o = np.log(np.maximum(G, _EPS)) - np.log(e)
    ur_o = np.log(np.maximum(R, _EPS)) - np.log(e)
    g_edge = live  # gDNA flows genomically everywhere there is mass
    r_edge = live & (R > 1.0)  # RNA smoothness only where the truth has RNA
    dg_s, ng = _disagreement(ug_s, left, g_edge)
    dg_o, _ = _disagreement(ug_o, left, g_edge)
    dr_s, nr = _disagreement(ur_s, left, r_edge)
    dr_o, _ = _disagreement(ur_o, left, r_edge)
    # convergence
    dfg = float((wm * np.abs(fg - prev_fg)).sum() / max(wm.sum(), _EPS)) if prev_fg is not None else np.nan
    # P_g modes
    if prior is None:
        top = np.array([])
        ismode = np.zeros(1, bool)
    else:
        P = np.exp(prior.logP - prior.logP.max())
        ismode = np.r_[False, (P[1:-1] > P[:-2]) & (P[1:-1] >= P[2:]), False]
        modes = prior.log_rho[ismode] / _LN10
        order = np.argsort(P[ismode])[::-1]
        top = modes[order][:3]
    return dict(
        mwae=mwae, dg_s=dg_s, dg_o=dg_o, dr_s=dr_s, dr_o=dr_o, ng=ng, nr=nr, dfg=dfg,
        nmode=int(ismode.sum()), modes=top, fg=fg,
    )


def run_condition_noprior(s):
    """Single pass-0 sweep with gdna_prior=None → PURE symmetric Beta(½,½) Jeffreys reference (½log f_g +
    ½log(1−f_g)) and NO σ²_transfer. The zero-infrastructure baseline: does the NPMLE+σ²_transfer stack beat
    the plain symmetric reference?"""
    belief = _sweep(s, s["b0"], None, transfer_variance=False)
    return [_measure(s, belief, None, None)]


def run_condition(s, bandwidth, passes, cold=False, conf_frac=1.0, transfer_variance=True, no_floor=False):
    """Refit loop; return per-pass measurements. ``cold`` restarts each sweep from the init belief b0
    (isolates the prior-sharpening feedback from the warm-continue accumulation); default WARM-CONTINUES
    from the previous solve, mirroring calibrate.py's closure. ``conf_frac`` < 1 fits the refit prior ONLY
    on the most-confident fraction of live nodes (smallest Var(log f_g)) — the owner's HOLD-OUT idea: keep
    ambiguous, unreliable nodes OUT of the population fit so the prior is not trained on its own guesses.
    ``transfer_variance=False`` reproduces the OLD regime (σ²_transfer=0) to test whether wiring the
    projection broke the previously-beneficial refit. ``no_floor`` drops the DNA-background floor + pinned
    component (background=None) to A/B the PRODUCTION (floor-on) regime against the pre-P2/P3 regime the
    original study measured."""
    bg = None if no_floor else s["bg"]  # production default = floor+pinned ON
    prior = DensityNPMLE.fit(s["mass_g"], s["eff_g"], background=bg, bandwidth=bandwidth)  # pass-0
    belief = _sweep(s, s["b0"], prior, transfer_variance=transfer_variance)
    recs = [_measure(s, belief, prior, None)]
    mass_g, eff_g = s["mass_g"], s["eff_g"]
    live = (eff_g > 1e-9 * 1.001) & (mass_g > _EPS)
    for _ in range(passes):
        fg = np.asarray(belief.f_g, float)
        vg = np.asarray(belief.var_gdna, float)
        g_hat = fg * mass_g
        if conf_frac < 1.0:
            v = np.where(live, vg, np.inf)
            thr = np.quantile(v[np.isfinite(v)], conf_frac)  # keep smallest-variance conf_frac
            keep = live & (vg <= thr)
        else:
            keep = live
        prior = DensityNPMLE.fit(
            g_hat[keep], eff_g[keep], background=bg, var_g=vg[keep], bandwidth=bandwidth
        )
        belief = _sweep(s, s["b0"] if cold else belief, prior, transfer_variance=transfer_variance)
        recs.append(_measure(s, belief, prior, recs[-1]["fg"]))
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    ap.add_argument("--bandwidths", default="0.15")
    ap.add_argument("--passes", type=int, default=5)
    ap.add_argument("--cold", action="store_true", help="restart each sweep from init (no warm-continue)")
    ap.add_argument("--conf-frac", type=float, default=1.0,
                    help="fit refit prior only on most-confident fraction of nodes (hold out ambiguous)")
    ap.add_argument("--no-tv", action="store_true", help="σ²_transfer OFF (old regime) to isolate its refit effect")
    ap.add_argument("--no-floor", action="store_true",
                    help="drop the DNA-background floor+pinned (background=None) — the pre-P2/P3 regime")
    ap.add_argument("--compact", action="store_true", help="one line per condition: full mwae trajectory")
    ap.add_argument("--no-prior", action="store_true",
                    help="gdna_prior=None → pure symmetric Jeffreys, no σ²_transfer (zero-infra baseline)")
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))
    bws = [float(x) for x in a.bandwidths.split(",")]

    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        s = _setup(inp, index, cfg)
        short = c.replace("gdna_gdna300_", "")
        for bw in bws:
            if a.no_prior:
                recs = run_condition_noprior(s)
                print(f"{short[:40]:40s} NO-PRIOR (sym Jeffreys, no σ²_T)  "
                      f"mwae={recs[0]['mwae']:.4f}  disG {recs[0]['dg_s']:.2f}/orc{recs[0]['dg_o']:.2f}",
                      flush=True)
                break
            recs = run_condition(s, bw, a.passes, cold=a.cold, conf_frac=a.conf_frac,
                                 transfer_variance=not a.no_tv, no_floor=a.no_floor)
            tag = ("  COLD" if a.cold else "") + (f"  conf={a.conf_frac}" if a.conf_frac < 1.0 else "") \
                + ("  noTV" if a.no_tv else "") + ("  noFLOOR" if a.no_floor else "")
            if a.compact:
                traj = " ".join(f"{r['mwae']:.4f}" for r in recs)
                dgt = f"disG {recs[0]['dg_s']:.2f}->{recs[-1]['dg_s']:.2f}/orc{recs[0]['dg_o']:.2f}"
                print(f"{short[:34]:34s} bw{bw}{tag:14s} mwae: {traj}   {dgt}", flush=True)
                continue
            print(f"\n=== {short}   bw={bw}{tag} ===")
            print(f"{'pass':>4} {'mwae':>7} {'disG_sol':>9} {'disG_orc':>9} {'disR_sol':>9} "
                  f"{'disR_orc':>9} {'Δf_g':>7} {'nmode':>5}  modes(log10ρ)")
            for k, r in enumerate(recs):
                dfg = "" if np.isnan(r["dfg"]) else f"{r['dfg']:.4f}"
                modes = ",".join(f"{m:.1f}" for m in r["modes"])
                print(f"{k:>4} {r['mwae']:.4f}  {r['dg_s']:>8.3f} {r['dg_o']:>8.3f} {r['dr_s']:>8.3f} "
                      f"{r['dr_o']:>8.3f} {dfg:>7} {r['nmode']:>5}  {modes}")
            print(f"     edges: gDNA={recs[0]['ng']}  RNA={recs[0]['nr']}", flush=True)


if __name__ == "__main__":
    main()
