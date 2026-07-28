"""Paired per-condition deltas on the SUBSTRATE, attributing each landed step separately.

    P1d  = pk2    -> p1dp
    P4   = p1dp   -> p4land
    noown= p4land -> noown2   (the destination-own plug-in removal)
    P1e  = noown2 -> p1e1     (probe, RIGEL_P1E=1.0)

    python scratchpad/subacc_delta.py
"""

from __future__ import annotations

import numpy as np

_EPS = 1e-9
ARMS = ["pk2", "p1dp", "p4land", "noown2", "p1e1"]
D = {a: dict(np.load(f"/tmp/subacc_{a}.npz", allow_pickle=True)) for a in ARMS}


def sel(d):
    return (d["live"].astype(bool) & d["isr"].astype(bool)
            & (d["single"].astype(bool) | d["gonly"].astype(bool)) & ~d["intergenic"].astype(bool)
            & np.isfinite(d["fo"]))


def per_cond(a, kind="fg"):
    d = D[a]
    s = sel(d)
    out = {}
    for c in sorted(set(d["cond"].tolist())):
        m = s & (d["cond"] == c)
        w = d["mass"][m]
        if kind == "fg":
            e = np.abs(d["fg"][m] - d["fo"][m])
        else:
            e = np.abs(np.log(np.maximum(d["fg"][m] * d["mass"][m], 1.0))
                       - np.log(np.maximum(d["fo"][m] * d["mass"][m], 1.0)))
        out[c] = (float(np.average(e, weights=w)), float(w.sum()))
    return out


for kind, lab in (("fg", "substrate mwae |f_g-oracle|"), ("rho", "substrate dlog rho_g")):
    print("=" * 104)
    print(f"PAIRED per-condition A/B on the HYPERPRIOR SUBSTRATE — {lab}")
    print("=" * 104)
    print(f"{'step':<28}{'from':>9}{'to':>9}{'delta':>10}{'rel':>9}{'better':>8}{'worse':>7}{'flat':>7}")
    P = {a: per_cond(a, kind) for a in ARMS}
    for step, (x, y) in (("P1d  (graft premise var)", ("pk2", "p1dp")),
                         ("P4   (delete _far)", ("p1dp", "p4land")),
                         ("no-own-plugin", ("p4land", "noown2")),
                         ("P1e  (RIGEL_P1E=1.0)", ("noown2", "p1e1")),
                         ("NET  pk2 -> HEAD", ("pk2", "noown2")),
                         ("NET  pk2 -> HEAD+P1e", ("pk2", "p1e1"))):
        cs = list(P[x])
        wm = np.array([P[x][c][1] for c in cs])
        ex = np.array([P[x][c][0] for c in cs])
        ey = np.array([P[y][c][0] for c in cs])
        ax, ay = float(np.average(ex, weights=wm)), float(np.average(ey, weights=wm))
        b = int(np.sum(ey < ex - 1e-4))
        w_ = int(np.sum(ey > ex + 1e-4))
        print(f"{step:<28}{ax:>9.4f}{ay:>9.4f}{ay - ax:>+10.4f}{(ay - ax) / ax:>+9.1%}"
              f"{b:>8}{w_:>7}{len(cs) - b - w_:>7}")
    print()
