"""DISSECTION 3 — ZERO-gDNA: every non-zero gDNA call is a false positive, so trace them.

Owner: *"that's a much easier thing to dissect because you know there's zero DNA. Any nonzero DNA is a
false positive."* And: the data has to GUIDE us to zero — we cannot assume it — so the question is not
"can we force zero" but "where does the FP mass come from, and does the prior amplify or suppress it".

`gdna_d2_floor.py` established the prior is far more PERMISSIVE than the delta-pin it replaces: at
log10 rho = +1 on a zero-gDNA library the delta-pin sits 25 nats below its peak and the landscape only 5 —
and that gap survives removing the floor entirely, so it is the KERNELS, not the smoothing.

Two candidate explanations, and they call for opposite fixes:
  FAITHFUL   pass-0 really does place many nodes up there, the landscape reports them, and the re-solve
             believes its own output. That is CIRCULARITY, and the fix is upstream in pass-0.
  SMEARED    only a handful of nodes are up there, but the kNN rule gives an ISOLATED node the WIDEST
             kernel, converting one observation into a plateau. That is an estimator defect, and the fix
             is here.

Part A settles which. Part B follows the FP mass into node classes, so whichever it is, it is localized.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d3_zero.py
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n5_roughness import region_substrate  # noqa: E402

import rigel.calibration.gdna_landscape as G  # noqa: E402

LN10 = np.log(10.0)
CACHE = Path("/Users/mkiyer/proj/rigel/scratchpad/gdna_dissect_cache.pkl")
NT = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}


def part_a():
    """FAITHFUL or SMEARED: how much TRAINING WEIGHT actually sits above rho = 0, and how wide are those
    nodes' kernels?"""
    print("=== A. Is the prior faithful to pass-0, or smearing a handful of nodes? (ZERO-gDNA) ===\n")
    print(f"{'condition':44s} {'n_train':>8s} {'n>0':>5s} {'wt>0':>8s} {'med h of those':>15s} "
          f"{'med h overall':>14s} {'logP(+1)':>9s}")
    for s in L.load_scenarios("ambig"):
        if s["group"][1] != "none" or s["group"][2] != "0.50":
            continue
        mk = L.masks(s)
        sel = region_substrate(s, mk)
        count, mass, eff = s["g_hat"][sel], s["mass"][sel], s["eff"][sel]
        grid = G._grid(mass, eff)
        centres = np.clip(np.log10(np.maximum(count, 1.0)) - np.log10(eff), grid[0], grid[-1])
        h = G.knn_widths(centres, float(grid[1] - grid[0]))
        w = G._reliability(count, s["var"][sel], mk["struct_zero"][sel])
        hi = centres > 0.0
        ls = G.fit_gdna_landscape(count, mass, eff, s["var"][sel], anchor=mk["struct_zero"][sel])
        lp = float(np.interp(1.0 * LN10, ls.log_rho, ls.logP) - ls.logP.max())
        print(f"{s['cond'][5:]:44s} {sel.sum():8d} {int(hi.sum()):5d} "
              f"{float(w[hi].sum() / max(w.sum(), 1e-9)):8.4f} "
              f"{(float(np.median(h[hi])) if hi.any() else np.nan):15.3f} "
              f"{float(np.median(h)):14.3f} {lp:9.1f}")
    print("\n  If `n>0` is a handful yet their kernels are the WIDEST, the plateau is an artefact of the kNN")
    print("  rule: an isolated node has no near neighbours BY DEFINITION, so the rule hands it the widest")
    print("  kernel in the fit — one observation rendered as decades of population mass.")


def part_b():
    """Where the surviving FP mass actually is, after the re-solve."""
    scen = pickle.loads(CACHE.read_bytes())
    print("\n=== B. The FP mass after the re-solve, by node class (ZERO-gDNA; all gDNA is false) ===\n")
    print(f"{'condition':40s} {'class':11s} {'nodes':>6s} {'FP pass-0':>11s} {'FP refit':>10s} {'ratio':>7s}")
    tot = {}
    for s in scen:
        if s["group"][1] != "none":
            continue
        live = (s["eff"] > 1e-9) & (s["mass"] > 1e-12)
        for t in (0, 1, 2, 3):
            m = live & (s["ntype"] == t)
            if not m.any():
                continue
            a = float((s["f0"][m] * s["mass"][m]).sum())
            b = float((s["f1"][m] * s["mass"][m]).sum())
            acc = tot.setdefault((s["group"][2], NT[t]), [0.0, 0.0, 0])
            acc[0] += a
            acc[1] += b
            acc[2] += int(m.sum())
    print(f"{'':40s} {'':11s} {'':>6s} {'':>11s} {'':>10s}")
    for (ss, cls), (a, b, n) in sorted(tot.items()):
        print(f"{('unstranded' if ss == '0.50' else 'stranded'):40s} {cls:11s} {n:6d} "
              f"{a:11.0f} {b:10.0f} {b / max(a, 1e-9):7.2f}")
    print("\n  ratio > 1 ⇒ the re-solve ADDS false gDNA in that class; < 1 ⇒ it removes it.")


if __name__ == "__main__":
    part_a()
    part_b()
