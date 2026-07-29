"""ADVERSARIAL RE-CHECK of the contamination-holdout track's substrate numbers.

Three things the original ho_subz2.py got structurally wrong or left open:

 (A) SUBSTRATE DEFINITION.  `_fit_gdna_hyperprior` on the SHIPPED path (background_floor=True,
     gdna_prior_additive=False) selects
         live & REGION & (single|gonly) & ~intergenic
     -- the `& ~intergenic` fires because `background is not None`.  ho_subz2.py used
         live & REGION & ~ambig
     which KEEPS intergenic.  Intergenic is exact-by-construction (f_g=1, err=0, never moves), so it
     dilutes both the substrate mwae and every Delta measured on it.  Report both.

 (B) z2 CROSS-ARM COMPARISON.  ho_subz2.py held the THRESHOLD fixed but applied it to each arm's OWN
     var, so the confident SET is re-selected per arm.  Report both the re-selected version and the
     HELD-FIXED-NODE-SET version (set frozen from the reference arm).

 (C) WRONG QUANTITY.  The prior is fitted on g_hat = f_g*mass against eff_global, i.e. on the DENSITY
     rho = f_g*mass/eff in log space.  Rank the arms on |log10 rho_hat/rho_oracle| too.

    python scratchpad/v_ho_subcheck.py Vhd0 Vamb0 Vexcl0
"""

from __future__ import annotations

import sys

import numpy as np

import sys as _sys
_sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from z2 import lin_var  # noqa: E402  -- THE single z2 denominator (log Var -> linear Var)

_EPS = 1e-9
ARMS = sys.argv[1:]
D = {a: dict(np.load(f"/tmp/subacc_{a}.npz", allow_pickle=True)) for a in ARMS}


def masks(d):
    isr = d["isr"].astype(bool)
    live = d["live"].astype(bool)
    amb = d["ambig"].astype(bool)
    inter = d["intergenic"].astype(bool)
    ok = np.isfinite(d["fo"]) & (d["mass"] > _EPS)
    return {
        "ok": ok,
        "FIT (shipped: -intergenic)": live & isr & ~amb & ~inter & ok,
        "agent SUBSTRATE (+interg.)": live & isr & ~amb & ok,
        "  x intergenic only": live & isr & ~amb & inter & ok,
        "EXCLUDED (ambig|bnd)": ((~isr) | (isr & amb)) & ok,
    }


def wmean(x, w):
    return float(np.average(x, weights=w)) if w.sum() > 0 else float("nan")


ref = D[ARMS[0]]
mref = masks(ref)
# a node set frozen ONCE, from the reference arm -> every arm scores the SAME nodes
finref = np.isfinite(ref["var"]) & mref["ok"]
THR = float(np.quantile(ref["var"][finref], 0.25))
CONF_FIXED = finref & (ref["var"] <= THR)

print(f"reference arm = {ARMS[0]}   Var(f_g) Q1 threshold = {THR:.6g}   "
      f"|frozen confident set| = {int(CONF_FIXED.sum()):,}")
print()

pops = [k for k in mref if k != "ok"] + ["ok"]
hdr = f"{'arm':<9}{'population':<30}{'nodes':>9}{'mass%':>8}{'mwae':>9}{'|dlog10rho|':>13}"
hdr += f"{'z2 RESEL':>10}{'z2 FIXED':>10}{'CW FIXED':>12}"
print(hdr)
print("-" * len(hdr))
tot_mass = ref["mass"][mref["ok"]].sum()
for a in ARMS:
    d = D[a]
    m = masks(d)
    err = np.abs(d["fg"] - d["fo"])
    # density in log10 -- both sides need a positive gDNA mass and a positive eff length
    ghat = d["fg"] * d["mass"]
    gorc = d["gor"]
    eff = d["eff"]
    dens_ok = (ghat > _EPS) & (gorc > _EPS) & (eff > _EPS)
    dlr = np.where(dens_ok, np.abs(np.log10(np.maximum(ghat, _EPS) / np.maximum(gorc, _EPS))), np.nan)
    fin = np.isfinite(d["var"]) & m["ok"]
    conf_resel = fin & (d["var"] <= THR)

    def z2(sel, conf):
        k = sel & conf & fin
        den = float(np.sum(d["mass"][k] * lin_var(d["var"][k], d["fg"][k])))
        return float(np.sum(d["mass"][k] * err[k] ** 2)) / den if den > 0 else float("nan")

    for lab in pops:
        sel = m[lab] if lab != "ok" else m["ok"]
        name = "ALL" if lab == "ok" else lab
        if not sel.any():
            continue
        dsel = sel & dens_ok
        print(f"{a:<9}{name:<30}{int(sel.sum()):>9,}{100 * d['mass'][sel].sum() / tot_mass:>8.1f}"
              f"{wmean(err[sel], d['mass'][sel]):>9.4f}"
              f"{wmean(dlr[dsel], d['mass'][dsel]):>13.4f}"
              f"{z2(sel, conf_resel):>10.2f}{z2(sel, CONF_FIXED):>10.2f}"
              f"{float(np.sum(err[sel & CONF_FIXED] * d['mass'][sel & CONF_FIXED])):>12,.0f}")
    print("-" * len(hdr))

if len(ARMS) > 1:
    print(f"\nDELTA vs {ARMS[0]} (negative = better).  Node sets FROZEN from {ARMS[0]}.")
    print(f"  {'population':<30}" + "".join(f"{a:>26}" for a in ARMS[1:]))
    print(f"  {'':<30}" + "".join(f"{'d mwae':>13}{'d|dlog10rho|':>13}" for _ in ARMS[1:]))
    d0 = D[ARMS[0]]
    err0 = np.abs(d0["fg"] - d0["fo"])
    g0 = d0["fg"] * d0["mass"]
    ok0 = (g0 > _EPS) & (d0["gor"] > _EPS) & (d0["eff"] > _EPS)
    for lab in pops:
        sel = mref[lab] if lab != "ok" else mref["ok"]
        name = "ALL" if lab == "ok" else lab
        if not sel.any():
            continue
        row = f"  {name:<30}"
        for a in ARMS[1:]:
            d = D[a]
            e = np.abs(d["fg"] - d["fo"])
            g = d["fg"] * d["mass"]
            dk = ok0 & (g > _EPS) & sel
            dl = np.abs(np.log10(np.maximum(g, _EPS) / np.maximum(d["gor"], _EPS)))
            dl0 = np.abs(np.log10(np.maximum(g0, _EPS) / np.maximum(d0["gor"], _EPS)))
            row += (f"{wmean(e[sel], d['mass'][sel]) - wmean(err0[sel], d0['mass'][sel]):>+13.4f}"
                    f"{wmean(dl[dk], d['mass'][dk]) - wmean(dl0[dk], d0['mass'][dk]):>+13.4f}")
        print(row)
