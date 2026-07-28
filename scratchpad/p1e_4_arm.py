"""P1e (5) — THE MECHANISM TEST. On P1d's residual-regression stratum (unstranded x gDNA-rich x capture-OFF),
does the DERIVED attribution  dv_c = bhat2_cons * (s_c/aSa)^2  put the damping on the arm that is WRONG?

P1b measured on that population: the RNA claim is a near-exact MEASUREMENT (36.6453 claimed vs 36.6734 true)
and the gDNA claim is 47x too big. P1d damps the RNA arm (the graft's premise). P1e is supposed to charge the
premium to the arm whose premise actually failed. This checks it, per message, against the oracle.

    OMP_NUM_THREADS=1 python scratchpad/p1e_4_arm.py
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import p1e_lib as L  # noqa: E402

STRATUM = [
    "gdna100_ss0.50_none_capOFF",
    "gdna300_ss0.50_none_capOFF",
    "gdna300_ss0.50_present_capOFF",
]
CONTROL = ["gdna300_ss0.99_present_capON", "gdna100_ss0.50_present_VERYSTRONG"]
_EPS = 1e-9

for name in STRATUM + CONTROL:
    inp, dbg = L.solve(L.CONDS[name])
    t, nf = L.message_table(inp, dbg)
    cap = dbg["capture"]
    us = cap["_uni_static"]
    # oracle per-component densities at the DESTINATION node, in its own frame
    M, E_g, E_r, fo = nf["M"], nf["E_g"], nf["E_r"], nf["fo"]
    og_t = fo * M / np.maximum(E_g, _EPS)
    or_t = (1.0 - fo) * M / np.maximum(E_r, _EPS)
    # the message's claimed densities: recovered from alpha and S
    d = t["dst"]
    S = t["S"]
    cg = t["alpha_g"] * S / np.maximum(E_g[d], _EPS)
    cr = (t["alpha_p"] + t["alpha_n"]) * S / np.maximum(E_r[d], _EPS)
    lg = np.log(np.maximum(cg, _EPS) / np.maximum(og_t[d], _EPS))
    lr = np.log(np.maximum(cr, _EPS) / np.maximum(or_t[d], _EPS))
    sel = (
        (t["nsup"] > 0) & t["lam_emit"] & (t["cls"] == "exon")
        & np.isfinite(lg) & np.isfinite(lr) & (og_t[d] > _EPS) & (or_t[d] > _EPS)
    )
    e = sel & (t["bhat2"] > 0)  # where P1e would actually fire
    print(f"\n{'=' * 122}\n{name}   ({int(sel.sum())} lambda-emitting graft messages into exons; "
          f"{int(e.sum())} with bhat2_cons > 0)\n{'=' * 122}")
    if sel.sum() < 5:
        continue
    print("  WHICH ARM IS WRONG (oracle)      mean log(claimed/true)      exp of it")
    print(f"    gDNA arm  : {np.mean(lg[sel]):>+8.3f}  ->  {np.exp(np.mean(lg[sel])):>8.2f}x   "
          f"(median {np.exp(np.median(lg[sel])):.2f}x)")
    print(f"    RNA  arm  : {np.mean(lr[sel]):>+8.3f}  ->  {np.exp(np.mean(lr[sel])):>8.2f}x   "
          f"(median {np.exp(np.median(lr[sel])):.2f}x)")
    print(f"    mass-weighted:  gDNA {np.exp(np.average(lg[sel], weights=M[d][sel])):.2f}x   "
          f"RNA {np.exp(np.average(lr[sel], weights=M[d][sel])):.2f}x")
    print("\n  WHERE THE DERIVED ATTRIBUTION PUTS THE DAMPING   dv_c = bhat2*(s_c/aSa)^2")
    dg, dr = t["dv_g"], t["dv_r"]
    sh = dg / np.maximum(dg + dr, _EPS)
    print(f"    {'':<30}{'median':>12}{'mean':>12}{'mass-wtd':>12}")
    for tag, v in (("dv_gDNA", dg), ("dv_RNA", dr)):
        print(f"    {tag:<30}{np.median(v[e]):>12.4f}{np.mean(v[e]):>12.4f}"
              f"{np.average(v[e], weights=M[d][e]):>12.4f}")
    print(f"    {'gDNA SHARE of the damping':<30}{np.median(sh[e]):>12.3f}{np.mean(sh[e]):>12.3f}"
          f"{np.average(sh[e], weights=M[d][e]):>12.3f}")
    print(f"    share of messages with dv_g > dv_r : {100 * np.mean(dg[e] > dr[e]):.1f}%")
    print(f"    the inputs:   median s_g = {np.median(t['s_g'][e]):.4g}   "
          f"median s_R (=s_p+s_n) = {np.median((t['s_p'] + t['s_n'])[e]):.4g}   "
          f"median sigma_cm^2 = {np.median(t['s2cm'][e]):.4g}")
    print(f"                  median alpha_g = {np.median(t['alpha_g'][e]):.3f}   "
          f"median w_g = {np.median(t['w_g'][e]):.4g}   median w_R = {np.median(t['w_r'][e]):.4g}")
    # THE TEST: is the damping share aligned with which arm is actually wrong?
    wrong_g = np.abs(lg) > np.abs(lr)
    print(f"\n  ALIGNMENT: the gDNA arm is the more-wrong one on "
          f"{100 * np.mean(wrong_g[e]):.1f}% of these messages; the damping is majority-gDNA on "
          f"{100 * np.mean((dg > dr)[e]):.1f}%")
    print(f"    Spearman( gDNA damping share , |log gDNA err| - |log RNA err| ) = "
          f"{L.spearman(sh[e], (np.abs(lg) - np.abs(lr))[e]):+.3f}")
    # what P1d does by comparison: it damps the RNA arm only (the graft's premise)
    print("    [P1d, for contrast, charges its omega to the RNA arm of EVERY graft edge]")
