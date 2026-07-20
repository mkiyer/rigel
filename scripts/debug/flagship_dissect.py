"""FLAGSHIP DISSECTION — localize the unstranded-capture error and find WHY, before any premise change.

The unstranded+capture condition is the DIAGNOSTIC (it exposes problems that strand-specific data hides by
self-solving). This tool decomposes the pass-0 solve into its THREE stages via the `_capture` hook —

    fg_strand : strand likelihood ALONE (no prior, no messages)
    fg_loc    : strand + global NPMLE prior (NO messages)
    f_g       : + messages (the final belief)

— and compares each to the ORACLE f_g = G/(G+R). This says exactly where the error enters (strand? prior?
messages?) and lets us BUCKET it by node class, so we see WHICH subset breaks rather than assuming it is all of
them. Special attention to the message stage (`f_g` vs `fg_loc`): a node whose |error| GROWS when messages
arrive is "responding incorrectly to messages" — the thing to debug.

Also checks the owner's candidate bugs directly:
  * STRAND deposit — in an unstranded library the deposited (u_pos,u_neg) must be ~50/50 for BOTH gDNA and RNA;
    a systematic asymmetry keyed to the oracle strand would be a wrong-strand deposit bug.
  * MESSAGE direction — does the received gDNA message push f_g toward or away from the oracle?

    OMP_NUM_THREADS=1 python flagship_dissect.py [--condition C] [--top 25]
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

from rigel.calibration.bp_solver import REGION, node_sweep
from rigel.calibration.npmle import DensityNPMLE
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12


def _oracle_strands(inp, chain):
    """Per chain-node oracle (Gp, Gn, Rp, Rn) — gDNA and RNA mass split by strand."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION

    def split(p):
        gp = np.asarray(p["gdna_pos"], float)
        gn = np.asarray(p["gdna_neg"], float)
        rp = np.asarray(p["mat_uns_pos"], float) + np.asarray(p["nas_uns_pos"], float)
        rn = np.asarray(p["mat_uns_neg"], float) + np.asarray(p["nas_uns_neg"], float)
        return gp, gn, rp, rn

    gR = split(inp["region_pools"])
    gB = split(inp["boundary_pools"])
    out = []
    for a, b in zip(gR, gB):
        ri = np.clip(idx, 0, a.shape[0] - 1)
        bi = np.clip(idx, 0, b.shape[0] - 1)
        out.append(np.where(isr, a[ri], b[bi]))
    return out  # Gp, Gn, Rp, Rn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--condition", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--bandwidth", type=float, default=0.15)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--no-floor", action="store_true",
                    help="drop the DNA-background floor+pinned (background=None) — the pre-P2/P3 regime")
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    cc = cfg.calibration
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    inp = _scan_and_truth(suite, a.condition, index, cfg, work, cache)
    s = _setup(inp, index, cfg)
    chain = s["chain"]
    mass_g, eff_g, G, R = s["mass_g"], s["eff_g"], s["G"], s["R"]
    Gp, Gn, Rp, Rn = _oracle_strands(inp, chain)

    # pass-0 solve (production config: σ²_transfer ON, floor+pinned ON) with the full capture hook
    cap = {}
    bg = None if a.no_floor else s["bg"]
    prior = DensityNPMLE.fit(mass_g, eff_g, background=bg, bandwidth=a.bandwidth)
    belief = node_sweep(
        chain, s["st"], s["geom"], s["b0"], s["ra"], s["bsub"], rna_sense_frac=s["kappa"],
        gdna_strand_overdispersion=s["od_g"], rna_strand_overdispersion=s["od_r"], n_grid=cc.sweep_n_grid,
        logodds_window=cc.sweep_logodds_window, n_tilt=cc.sweep_n_tilt,
        n_grid_ss=cc.sweep_n_grid_single_strand, gdna_prior=prior, transfer_variance=True, _capture=cap,
    )
    fg = np.asarray(belief.f_g, float)
    fg_loc = np.asarray(cap["fg_loc"], float)
    fg_str = np.asarray(cap["fg_strand"], float)
    free_p = np.asarray(cap["free_pos"], bool)
    free_n = np.asarray(cap["free_neg"], bool)
    prec_g = np.asarray(cap["prec_g"], float)
    mode_g = np.asarray(cap["mode_g"], float)

    tot = G + R
    fo = np.where(tot > _EPS, G / np.maximum(tot, _EPS), np.nan)
    live = (eff_g > 1e-9 * 1.001) & (mass_g > _EPS) & np.isfinite(fo) & (tot > _EPS)
    isr = np.asarray(chain.kind) == REGION

    # structural class from the emission gate (what the solver actually uses)
    cls = np.where(free_p & free_n, "AMBIG", np.where(free_p | free_n, "single", "gonly"))
    cls = np.where(~isr, "bndry", cls)  # boundaries called out separately

    def wmae(mask, x):
        ok = mask & live & np.isfinite(x) & np.isfinite(fo) & np.isfinite(mass_g)
        w = np.where(ok, mass_g, 0.0)
        num = float((w * np.abs(np.where(ok, x, 0.0) - np.where(ok, fo, 0.0))).sum())
        return num / max(w.sum(), _EPS), float(w.sum())

    print(f"=== {a.condition}  bw={a.bandwidth} ===")
    print(f"live nodes={int(live.sum())}  total gDNA true frac={G[live].sum()/max(tot[live].sum(),1):.3f}\n")
    print("STAGE errors (mass-weighted |f_g − oracle|), overall then by class:")
    print(f"{'class':>8} {'massfrac':>8} {'strand':>8} {'+prior':>8} {'+msgs':>8}  {'prior_eff':>9} {'msg_eff':>8}")
    for name in ("ALL", "AMBIG", "single", "gonly", "bndry"):
        m = live if name == "ALL" else (cls == name) & live
        e_str, w = wmae(m, fg_str)
        e_loc, _ = wmae(m, fg_loc)
        e_fin, _ = wmae(m, fg)
        mf = w / max(np.where(live, mass_g, 0.0).sum(), _EPS)
        print(f"{name:>8} {mf:>8.3f} {e_str:>8.3f} {e_loc:>8.3f} {e_fin:>8.3f}  "
              f"{e_loc - e_str:>+9.3f} {e_fin - e_loc:>+8.3f}")

    # signed error at the final stage: are we UNDER- or OVER-calling gDNA?
    signed = np.where(np.isfinite(fg) & np.isfinite(fo), fg - fo, 0.0)
    print("\nSIGNED final error (f_g − oracle), mass-weighted mean by class (+ = over-call gDNA):")
    for name in ("AMBIG", "single", "gonly", "bndry"):
        m = (cls == name) & live & np.isfinite(fg)
        w = np.where(m, mass_g, 0.0)
        sgn = float((w * signed).sum() / max(w.sum(), _EPS))
        # split by oracle truth: is the node truly gDNA (fo>0.8) or RNA (fo<0.2)?
        print(f"  {name:>7}: mean signed={sgn:+.3f}   "
              f"[of these, truly-gDNA(fo>.8) massfrac={float(np.where(m&(fo>0.8),mass_g,0).sum()/max(w.sum(),_EPS)):.2f}, "
              f"truly-RNA(fo<.2)={float(np.where(m&(fo<0.2),mass_g,0).sum()/max(w.sum(),_EPS)):.2f}]")

    # STRAND DEPOSIT check: unstranded ⇒ deposited pos fraction ~0.5 for both gDNA and RNA truth
    gpos_frac = Gp[live].sum() / max((Gp[live].sum() + Gn[live].sum()), _EPS)
    rpos_frac = Rp[live].sum() / max((Rp[live].sum() + Rn[live].sum()), _EPS)
    print(f"\nSTRAND DEPOSIT (oracle, unstranded ⇒ expect ~0.5): gDNA pos-frac={gpos_frac:.3f}  "
          f"RNA pos-frac={rpos_frac:.3f}")

    # MESSAGE response: for nodes that RECEIVED a gDNA message, did it push f_g toward or away from oracle?
    got = live & (prec_g > _EPS) & np.isfinite(fg) & np.isfinite(fg_loc) & np.isfinite(fo)
    moved = np.where(got, np.abs(fg - fo) - np.abs(fg_loc - fo), 0.0)  # >0 = message HURT
    w = np.where(got, mass_g, 0.0)
    hurt = float((w * np.maximum(moved, 0)).sum())
    help_ = float((w * np.maximum(-moved, 0)).sum())
    print(f"\nMESSAGE RESPONSE on {int(got.sum())} nodes with a gDNA msg: "
          f"mass·hurt={hurt:.0f}  mass·help={help_:.0f}  net={'HURT' if hurt>help_ else 'HELP'}")
    for name in ("AMBIG", "single", "gonly", "bndry"):
        m = got & (cls == name)
        wm = np.where(m, mass_g, 0.0)
        h = float((wm * np.maximum(moved, 0)).sum())
        g = float((wm * np.maximum(-moved, 0)).sum())
        print(f"  {name:>7}: n={int(m.sum()):4d}  mass·hurt={h:>9.0f}  mass·help={g:>9.0f}")

    # WORST nodes by mass·|final error| — the localized subset, with the full trace
    err = np.where(live & np.isfinite(fg), mass_g * np.abs(fg - fo), 0.0)
    order = np.argsort(err)[::-1][: a.top]
    print(f"\nTOP {a.top} error nodes (node: cls  mass  oracle→strand→loc→final  msg(mode,prec)  Gp/Gn Rp/Rn):")
    for i in order:
        print(f"  {i:6d}: {cls[i]:>6} m={mass_g[i]:8.0f}  fo={fo[i]:.2f} str={fg_str[i]:.2f} "
              f"loc={fg_loc[i]:.2f} fin={fg[i]:.2f}  msg=({mode_g[i]:+.2f},{prec_g[i]:.2f})  "
              f"G {Gp[i]:.0f}/{Gn[i]:.0f} R {Rp[i]:.0f}/{Rn[i]:.0f}")


if __name__ == "__main__":
    main()
