"""CONFIDENT-WRONG NODE TRACE — the confidently-wrong nodes (large |f_g − oracle| AND low var_gdna = the solver
is SURE and WRONG) are what corrupt the gDNA-hyperprior fit. In the unstranded case the strand is silent, so
these are driven ENTIRELY by the message system (new: the (λ,θ) relay + σ²_transfer). This tool traces each one
back through: init → strand-only → local(+reference) → the received gDNA/RNA messages → final, AND shows the
chain NEIGHBOURS that emitted those messages (their truth, their solved f_g, the gDNA density they emit), so we
can tell whether a wrong-and-confident value is (a) an honest propagation from a wrong neighbour, (b) a message
that shouldn't have applied (gating), or (c) a mis-scaled message (a BUG in the new layer).

    OMP_NUM_THREADS=1 python confident_fp_trace.py [--condition C] [--top 15] [--dir over|under|both]
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
from flagship_interrogate import _oracle_per_node, _solve  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--condition", default="gdna_gdna5_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--dir", choices=["over", "under", "both"], default="over")
    ap.add_argument("--nodes", default=None, help="trace these specific node ids (comma-sep) instead of the top")
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    inp = _scan_and_truth(suite, a.condition, index, cfg, work, cache)
    ra = RegionArrays.from_index(index)

    dbg = _solve(inp, ra, 0)
    chain = dbg["chain"]
    cap = dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fg = np.asarray(cap["f_g"], float)
    var_g = np.asarray(cap["var_g"], float)
    fg_str = np.asarray(cap["fg_strand"], float)
    fg_loc = np.asarray(cap["fg_loc"], float)
    mode_g, prec_g = np.asarray(cap["mode_g"], float), np.asarray(cap["prec_g"], float)
    mode_p, prec_p = np.asarray(cap["mode_p"], float), np.asarray(cap["prec_p"], float)
    mode_n, prec_n = np.asarray(cap["mode_n"], float), np.asarray(cap["prec_n"], float)
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    egl, egr = np.asarray(cap["eff_gdna_l"], float), np.asarray(cap["eff_gdna_r"], float)
    isr = np.asarray(chain.kind) == REGION
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    cls = np.where(fp & fn, "AMBIG", np.where(fp | fn, "single", "gonly"))
    cls = np.where(~isr, "bndry", cls)
    left = np.asarray(chain.left, np.int64)
    right = np.full(left.shape, -1, np.int64)
    ok = left >= 0
    right[left[ok]] = np.where(ok)[0]

    tot = G + R
    fo = np.where(tot > _EPS, G / np.maximum(tot, _EPS), np.nan)
    live = (eff > 1e-9 * 1.001) & (mass > _EPS) & np.isfinite(fo) & (tot > _EPS)
    err = fg - fo  # signed
    vmed = float(np.median(var_g[live]))
    confident = live & (var_g < vmed)  # the solver is SURE (below-median belief variance)
    if a.dir == "over":
        wrong = confident & (err > 0.3) & (fo < 0.3)  # over-called: truly RNA, confidently called gDNA
    elif a.dir == "under":
        wrong = confident & (err < -0.3) & (fo > 0.7)  # under-called: truly gDNA, confidently called RNA
    else:
        wrong = confident & (np.abs(err) > 0.3)

    wm = np.where(wrong, mass, 0.0)
    print(f"=== {a.condition}   [Phase-1a]  dir={a.dir} ===")
    print(f"live={int(live.sum())}  median var_gdna={vmed:.3f}  "
          f"CONFIDENTLY-WRONG (var<median & |err|>0.3): n={int(wrong.sum())} massfrac="
          f"{wm.sum()/max(np.where(live,mass,0).sum(),_EPS):.3f}")
    by = {c: int((wrong & (cls == c)).sum()) for c in ("single", "AMBIG", "gonly", "bndry")}
    print(f"  by class: {by}")

    # who drives them? classify each wrong node's DOMINANT gDNA message source neighbour by its OWN truth.
    def _rho_emit(j):  # the gDNA density a neighbour j emits (its solved f_g over its gDNA eff-len)
        e = max(0.5 * (egl[j] + egr[j]), _EPS)
        return fg[j] * mass[j] / e

    order = ([int(x) for x in a.nodes.split(",")] if a.nodes else np.argsort(wm)[::-1][: a.top])
    print("\nTRACE  (strand→loc→FINAL(var); gDNA msg; RNA msg; NEIGHBOURS):")
    for i in order:
        rho_i = mass[i] / max(eff[i], _EPS)
        print(f"\n  node {i} [{cls[i]}]  m={mass[i]:.0f}  log10ρtot={np.log10(max(rho_i,1e-12)):.2f}")
        print(f"     ORACLE fo={fo[i]:.2f}  G={Gp[i]:.0f}/{Gn[i]:.0f}  R={Rp[i]:.0f}/{Rn[i]:.0f}  (unstr ⇒ ~50/50)")
        print(f"     solve: strand={fg_str[i]:.2f}  loc(+ref)={fg_loc[i]:.2f}  FINAL={fg[i]:.2f}  "
              f"var={var_g[i]:.3f} {'CONFIDENT' if var_g[i] < vmed else ''}")
        print(f"     gDNA msg: mode={mode_g[i]:+.2f} prec={prec_g[i]:.2f}   "
              f"RNA msg: +({mode_p[i]:+.2f},{prec_p[i]:.2f}) −({mode_n[i]:+.2f},{prec_n[i]:.2f})")
        for side, j in (("L", int(left[i])), ("R", int(right[i]))):
            if j < 0:
                print(f"     {side}-nbr: (terminal)")
                continue
            print(f"     {side}-nbr {j} [{cls[j]}]: fo={fo[j] if np.isfinite(fo[j]) else -1:.2f} "
                  f"f_g={fg[j]:.2f} m={mass[j]:.0f}  emits gDNA log10ρ={np.log10(max(_rho_emit(j),1e-12)):.2f}")


if __name__ == "__main__":
    main()
