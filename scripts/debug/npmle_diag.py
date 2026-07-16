"""Fast diagnostic on the cached per-pass beliefs — what is the fit set made of, and what drives P(ρ)?

Loads the compact belief npz (instant — no sweep) and dissects the NPMLE fit set: the count-0 share (the
would-be zero anchor), the discreteness comb (ĝ≈1 over short boundary E), the total believed gDNA, and the
belief widths for the residual low-ĝ nodes (does the honest-width firewall have teeth?).

    python scripts/debug/npmle_diag.py --beliefs DIR/COND.beliefs.npz[,...]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

_EPS = 1e-12


def _q(v, w=None):
    v = np.asarray(v, float)
    if v.size == 0:
        return "∅"
    qs = np.percentile(v, [10, 50, 90])
    return f"p10={qs[0]:.3g} p50={qs[1]:.3g} p90={qs[2]:.3g}"


def run(npz_path):
    cond = Path(npz_path).stem.replace(".beliefs", "")
    d = np.load(npz_path)
    k = d["k"]
    E = d["eff_global"]
    e_max = float(d["e_max"])
    print(f"\n===== {cond}  (n_nodes={k.shape[0]}, E_max={e_max:.0f}) =====")
    short = E < 500.0  # boundary-ish (short exposure) vs region
    for p in range(3):
        fg = d["fg"][p]
        vg = d["varg"][p]
        fit = np.isfinite(vg) & (E > _EPS)
        g = fg * k
        n = int(fit.sum())
        kf = k[fit]
        gf = g[fit]
        vgf = vg[fit]
        shortf = short[fit]
        comb = (gf >= 0.5) & (gf <= 1.5)  # the ĝ≈1 count-1 spike
        print(
            f"[pass {p}] fit={n}  k==0:{np.mean(kf == 0):.3f}  ĝ<0.1:{np.mean(gf < 0.1):.3f}  "
            f"ĝ∈[.5,1.5]:{np.mean(comb):.3f}  Σĝ/Σk={gf.sum() / max(kf.sum(), 1):.4f}  "
            f"Σĝ={gf.sum():.0f}"
        )
        print(
            f"          comb(ĝ≈1): frac_short(E<500)={np.mean(shortf[comb]) if comb.any() else 0:.2f}  "
            f"medE={np.median(E[fit][comb]) if comb.any() else float('nan'):.0f}  "
            f"var_g[{_q(vgf[comb])}]"
        )
        lo = gf < 0.1
        print(
            f"          residual ĝ<0.1: k==0 within={np.mean(kf[lo] == 0):.3f}  "
            f"f_g[{_q(fg[fit][lo])}]  var_g[{_q(vgf[lo])}]"
        )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--beliefs", required=True)
    args = ap.parse_args()
    for p in args.beliefs.split(","):
        run(p)


if __name__ == "__main__":
    main()
