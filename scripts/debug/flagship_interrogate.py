"""FLAGSHIP INTERROGATE — the single, production-FAITHFUL dissection of the unstranded-capture error.

Supersedes the pass-0-only reconstructions (`flagship_dissect.py`, `prior_crush_probe.py`, the
`refit_loop_study.py` harness): those rebuild the solve from primitives; this one drives the REAL
`rigel.calibration.calibrate.calibrate()` with its `_debug` hook, so what you read here is exactly what ships.
It reports BOTH the pass-0 belief (`calib_refit_iters=0`) AND the shipped belief (production config, refit on),
because production ships the refit and the refit moves f_g by up to ~1.0 on flagship nodes.

For the worst-error nodes it prints the full per-node trace from the production capture — oracle (per strand),
the three solve stages (strand-only → +prior → final), the prior-arm decomposition (does the gDNA arm favour or
oppose high f_g; how hard does the RNA reference oppose it; is the background floor active), the received gDNA
message, AND the chain NEIGHBOURS (their truth, their solved f_g, their gDNA density) — so we can see WHY the
messages do or do not lift a crushed enriched-gDNA node, rather than inferring it.

    OMP_NUM_THREADS=1 python flagship_interrogate.py [--condition C] [--top 15] [--stage shipped|pass0]
"""

from __future__ import annotations

import argparse
import dataclasses
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION
from rigel.calibration.calibrate import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12


def _oracle_per_node(inp, chain):
    """Per chain-node oracle gDNA/RNA mass, per strand, in chain order."""
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


def _solve(inp, ra, refit_iters):
    """Run the REAL production calibrate() with the given refit count; return (_debug, config)."""
    cfg = PipelineConfig()
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=int(refit_iters))
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return dbg


_LN10 = np.log(10.0)


def _deep_dive(i, chain, index, dbg0, cap0, mass, eff, fo, G, R, Gp, Gn, Rp, Rn, cls,
               left, right, rho, grid, glp, rna_ref, mode_g, prec_g, fg_str, fg_loc, fg0, fgS, live):
    """Error-by-error dig on ONE node: (1) single-exon transcript? (2) boundary densities?
    (3) what is the prior doing — is the enriched mode present, is it DNA or RNA, why is/isn't it pulling up?"""
    print(f"\n===== DEEP DIVE node {i} [{cls[i]}] =====")
    rtot = mass[i] / max(eff[i], 1e-12)
    print(f"mass={mass[i]:.0f} eff={eff[i]:.0f}  ρtot={rtot:.2f} (log10 {np.log10(max(rtot,1e-12)):.2f})")
    print(f"oracle fo={fo[i]:.3f}  G={Gp[i]:.0f}/{Gn[i]:.0f}  R(uns)={Rp[i]:.0f}/{Rn[i]:.0f}   "
          f"solve strand={fg_str[i]:.2f} +prior={fg_loc[i]:.2f} pass0={fg0[i]:.2f} shipped={fgS[i]:.2f}")

    # (1) transcript structure — is it a single-exon transcript?
    kind = np.asarray(chain.kind)
    if kind[i] == REGION:
        rid = int(np.asarray(chain.ref_idx)[i])
        rdf = index.nodes_df.iloc[rid]
        ref, rs, re_, sig = rdf["ref_name"], int(rdf["start"]), int(rdf["end"]), int(rdf["signature"])
        print(f"\n(1) STRUCTURE: region_id={rid} {ref}:{rs}-{re_} len={re_-rs} signature={sig:04b}")
        td = index.t_df
        ov = td[(td["ref"] == ref) & (td["start"] < re_) & (td["end"] > rs)]
        if len(ov) == 0:
            print("    no overlapping transcript (intergenic gDNA region)")
        for _, t in ov.iterrows():
            print(f"    transcript {t['t_id']} strand={t['strand']} n_exons={t['n_exons']} "
                  f"is_nrna={t['is_nrna']} span={t['start']}-{t['end']}  "
                  f"{'*** SINGLE-EXON ***' if t['n_exons'] == 1 else ''}")
    else:
        print("\n(1) STRUCTURE: this is a BOUNDARY node")

    # (2) boundary densities (the node's two chain neighbours)
    print("\n(2) BOUNDARY NEIGHBOURS:")
    for side, j in (("L", int(left[i])), ("R", int(right[i]))):
        if j < 0:
            print(f"    {side}: (reference terminal — propagation sink)")
            continue
        print(f"    {side} node {j} [{cls[j]}]: mass={mass[j]:.0f} oracle fo={fo[j]:.2f} "
              f"solved fg={fgS[j]:.2f} log10ρg={np.log10(max(rho[j],1e-12)):.2f}")
    print(f"    received gDNA message: mode={mode_g[i]:+.2f} prec={prec_g[i]:.3f} "
          f"({'WEAK — depleted/gagged' if prec_g[i] < 0.3 else 'active'})")

    # (3) the prior — modes, where the node sits, and the enriched mode's TRUTH
    prior = dbg0["gdna_prior"]
    x = np.asarray(prior.log_rho) / _LN10
    P = np.exp(np.asarray(prior.logP) - np.asarray(prior.logP).max())
    ismode = np.r_[False, (P[1:-1] > P[:-2]) & (P[1:-1] >= P[2:]), False]
    modes = x[ismode]
    print(f"\n(3) PRIOR modes (log10 ρ): {np.round(modes, 2)}")
    lp = glp[i]
    print("    gDNA arm logP_g over f_g (ρ_g=f_g·ρtot):  " +
          "  ".join(f"{f:.2f}:{lp[int(np.argmin(np.abs(grid-f)))]:+.2f}" for f in (0.01, 0.1, 0.5, 0.9, 0.99)))
    print("    RNA ref ½log(1−f_g):                       " +
          "  ".join(f"{f:.2f}:{rna_ref[int(np.argmin(np.abs(grid-f)))]:+.2f}" for f in (0.1, 0.5, 0.9, 0.99)))
    psi = lp + rna_ref
    print(f"    NET ψ (gDNA arm + RNA ref, strand flat) argmax f_g={grid[int(np.argmax(psi))]:.3f}  "
          f"ψ(.9)−ψ(.1)={psi[int(np.argmin(np.abs(grid-0.9)))]-psi[int(np.argmin(np.abs(grid-0.1)))]:+.2f}")
    gpull = lp[int(np.argmin(np.abs(grid-0.9)))] - lp[int(np.argmin(np.abs(grid-0.1)))]
    rpull = rna_ref[int(np.argmin(np.abs(grid-0.9)))] - rna_ref[int(np.argmin(np.abs(grid-0.1)))]
    print(f"    => gDNA-arm pull to f_g=.9: {gpull:+.2f}   RNA-ref pull: {rpull:+.2f}   "
          f"{'PRIOR LOSES to RNA ref' if gpull < -rpull else 'prior wins'}")

    # enriched-mode TRUTH: of the population mass at the node's density, how much is really gDNA vs RNA?
    if modes.size:
        enr_mode = modes[modes > modes.min() + 0.3]
        emode = enr_mode.max() if enr_mode.size else modes.max()
        ltot = np.log10(np.maximum(mass / np.maximum(eff, 1e-12), 1e-12))
        near = live & (np.abs(ltot - emode) < 0.25)
        w = np.where(near, mass, 0.0)
        tg = float(np.where(near & (fo > 0.8), mass, 0).sum() / max(w.sum(), 1e-12))
        tr = float(np.where(near & (fo < 0.2), mass, 0).sum() / max(w.sum(), 1e-12))
        fbar = float((w * fo).sum() / max(w.sum(), 1e-12))
        print(f"    ENRICHED MODE at log10ρ≈{emode:.2f}: population mass there is {tg:.2f} truly-gDNA, "
              f"{tr:.2f} truly-RNA, mean oracle f_g={fbar:.2f} (n={int(near.sum())})")
        print(f"    => the enriched mode is {'DNA' if fbar > 0.6 else ('RNA' if fbar < 0.4 else 'MIXED')} "
              f"in truth; the prior is fit on TOTAL density so it cannot see this split.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--condition", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--stage", choices=["shipped", "pass0"], default="shipped",
                    help="which belief to rank/trace: pass0 (refit=0) or shipped (production refit)")
    ap.add_argument("--node", type=int, default=None,
                    help="deep-dive ONE node: transcript structure (single-exon?), boundary densities, and the "
                         "prior/enriched-mode behaviour")
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    inp = _scan_and_truth(suite, a.condition, index, cfg, work, cache)
    ra = RegionArrays.from_index(index)

    dbg0 = _solve(inp, ra, 0)          # pass-0
    dbgS = _solve(inp, ra, cfg.calibration.calib_refit_iters)  # shipped (refit on)
    chain = dbg0["chain"]
    cap0 = dbg0["capture"]
    capS = dbgS["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    mass = np.asarray(cap0["mass_global"], float)
    eff = np.asarray(cap0["eff_global"], float)
    left = np.asarray(chain.left, np.int64)
    right = np.full(left.shape, -1, np.int64)
    ok = left >= 0
    right[left[ok]] = np.where(ok)[0]

    tot = G + R
    fo = np.where(tot > _EPS, G / np.maximum(tot, _EPS), np.nan)
    live = (eff > 1e-9 * 1.001) & (mass > _EPS) & np.isfinite(fo) & (tot > _EPS)
    isr = np.asarray(chain.kind) == REGION
    fp = np.asarray(cap0["free_pos"], bool)
    fn = np.asarray(cap0["free_neg"], bool)
    cls = np.where(fp & fn, "AMBIG", np.where(fp | fn, "single", "gonly"))
    cls = np.where(~isr, "bndry", cls)

    fg0 = np.asarray(cap0["f_g"], float)
    fgS = np.asarray(capS["f_g"], float)

    def mwae(fg, m=None):
        mask = live if m is None else (m & live)
        w = np.where(mask, mass, 0.0)
        return float((w * np.abs(fg - np.where(mask, fo, 0.0))).sum() / max(w.sum(), _EPS))

    print(f"=== {a.condition} ===   [PRODUCTION calibrate(), refit={cfg.calibration.calib_refit_iters}]")
    print(f"live={int(live.sum())}  true gDNAfrac={G[live].sum()/max(tot[live].sum(),1):.3f}  "
          f"mwae pass0={mwae(fg0):.4f}  shipped={mwae(fgS):.4f}\n")
    print(f"{'class':>7} {'massfrac':>8} {'mwae0':>7} {'mwaeS':>7} {'signedS':>8}  truly-gDNA  truly-RNA")
    for name in ("ALL", "AMBIG", "single", "gonly", "bndry"):
        m = live if name == "ALL" else (cls == name) & live
        w = np.where(m, mass, 0.0)
        mf = w.sum() / max(np.where(live, mass, 0.0).sum(), _EPS)
        sgn = float((w * (fgS - np.where(m, fo, 0.0))).sum() / max(w.sum(), _EPS))
        tg = float(np.where(m & (fo > 0.8), mass, 0).sum() / max(w.sum(), _EPS))
        tr = float(np.where(m & (fo < 0.2), mass, 0).sum() / max(w.sum(), _EPS))
        print(f"{name:>7} {mf:>8.3f} {mwae(fg0,m):>7.3f} {mwae(fgS,m):>7.3f} {sgn:>+8.3f}  "
              f"{tg:>9.2f}  {tr:>9.2f}")

    # arm decomposition on the solve grid (strand-free global prior term already includes floor+pinned)
    grid = np.asarray(cap0["solve_grid"], float)
    glp = np.asarray(cap0["global_lp"], float)  # (n, K)
    rna_ref = 0.5 * np.log(np.clip(1 - grid, _EPS, 1.0))

    def at(vec, x):
        return float(vec[int(np.argmin(np.abs(grid - x)))])

    fg = fgS if a.stage == "shipped" else fg0
    cap = capS if a.stage == "shipped" else cap0
    mode_g = np.asarray(cap["mode_g"], float)
    prec_g = np.asarray(cap["prec_g"], float)
    fg_str = np.asarray(cap0["fg_strand"], float)
    fg_loc = np.asarray(cap0["fg_loc"], float)
    rho = fg * mass / np.maximum(eff, _EPS)  # solved gDNA density

    # POPULATION summary — quantify the two claims beyond the top-N anecdotes.
    enr = live & (mass / np.maximum(eff, _EPS) > 10 ** 0.5)  # "enriched": log10 ρtot > 0.5
    crushed = enr & (fo > 0.8) & (fgS < 0.2)  # truly-gDNA enriched node crushed to ~0
    # does a crushed node have ANY neighbour solved at an enriched gDNA density (log10 ρg > 0.5)?
    lg = np.log10(np.maximum(rho, 1e-12))
    nbr_enr = np.zeros(mass.shape, bool)
    for j_arr in (left, right):
        nbr_lg = np.take(np.r_[lg, -99.0], np.where(j_arr >= 0, j_arr, -1))
        has = (j_arr >= 0) & (nbr_lg > 0.5)
        nbr_enr |= has
    cm = np.where(crushed, mass, 0.0)
    isl = float(np.where(crushed & ~nbr_enr, mass, 0.0).sum() / max(cm.sum(), _EPS))
    print(f"\nPOPULATION: enriched truly-gDNA nodes crushed to <0.2: n={int(crushed.sum())} "
          f"massfrac={cm.sum()/max(np.where(live,mass,0).sum(),_EPS):.3f}  "
          f"of these {isl:.2f} have NO enriched-gDNA neighbour (islands)")
    go = live & (cls == "gonly")
    if go.any():
        lgo = np.log10(np.maximum(mass[go] / np.maximum(eff[go], _EPS), 1e-12))
        print(f"gonly (structural gDNA-only, solve perfect): n={int(go.sum())} massfrac="
              f"{mass[go].sum()/max(np.where(live,mass,0).sum(),_EPS):.3f}  log10ρtot "
              f"[{np.percentile(lgo,10):.2f}, {np.percentile(lgo,50):.2f}, {np.percentile(lgo,90):.2f}]  "
              f"enriched(>0.5): {float((lgo>0.5).mean()):.2f}")

    if a.node is not None:
        _deep_dive(a.node, chain, index, dbg0, cap0, mass, eff, fo, G, R, Gp, Gn, Rp, Rn, cls,
                   left, right, rho, grid, glp, rna_ref, mode_g, prec_g, fg_str, fg_loc, fg0, fgS, live)
        return

    err = np.where(live, mass * np.abs(fg - fo), 0.0)
    order = np.argsort(err)[::-1][: a.top]
    print(f"\nTOP {a.top} error nodes ({a.stage}); trace = oracle→strand→loc→pass0→shipped, prior arms, "
          f"msg, NEIGHBOURS:")
    for i in order:
        gA = at(glp[i], 0.9) - at(glp[i], 0.1)   # gDNA arm: + favours high f_g
        rA = at(rna_ref, 0.9) - at(rna_ref, 0.1)  # RNA arm: - opposes high f_g
        print(f"\n  node {i} [{cls[i]}]  m={mass[i]:.0f} E={eff[i]:.0f} log10ρtot={np.log10(max(mass[i]/max(eff[i],_EPS),1e-12)):.2f}")
        print(f"     oracle fo={fo[i]:.2f}  G={Gp[i]:.0f}/{Gn[i]:.0f}  R={Rp[i]:.0f}/{Rn[i]:.0f}")
        print(f"     solve: strand={fg_str[i]:.2f} +prior(loc)={fg_loc[i]:.2f} pass0={fg0[i]:.2f} "
              f"shipped={fgS[i]:.2f}")
        print(f"     prior arms: gDNA favours high-fg by {gA:+.2f}  |  RNA ref opposes by {rA:+.2f}  "
              f"|  msg=(mode {mode_g[i]:+.2f}, prec {prec_g[i]:.2f})")
        for side, j in (("L", left[i]), ("R", right[i])):
            if j < 0:
                print(f"     {side}-neighbour: (none)")
                continue
            print(f"     {side}-neighbour {j} [{cls[j]}]: oracle fo={fo[j]:.2f} solved fg={fg[j]:.2f} "
                  f"log10ρg={np.log10(max(rho[j],1e-12)):.2f}  m={mass[j]:.0f}")


if __name__ == "__main__":
    main()
