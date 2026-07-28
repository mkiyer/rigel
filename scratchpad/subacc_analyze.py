"""SUBSTRATE ACCURACY — score the arms on the population the gDNA hyperprior is actually fitted from.

Reads the /tmp/subacc_<arm>.npz dumps from `subacc_dump.py` and answers:
  (1) accuracy on SEL = live & REGION & (single|gonly) & ~intergenic  (the EM-path substrate)
  (2) the same on the gDNA DENSITY rho_g = f_g*M/E_g the prior is actually fitted on
  (3) the spread: what fraction of substrate nodes are off by >0.1 / >0.25
  (4) how far the achievable prior is from the oracle-fitted one

    python scratchpad/subacc_analyze.py pk2 p1dp noown2 p1e1
"""

from __future__ import annotations

import sys

import numpy as np

_EPS = 1e-9
ARMS = sys.argv[1:] or ["noown2"]
D = {a: dict(np.load(f"/tmp/subacc_{a}.npz", allow_pickle=True)) for a in ARMS}


def sets(d):
    isr = d["isr"].astype(bool)
    live = d["live"].astype(bool)
    single = d["single"].astype(bool)
    gonly = d["gonly"].astype(bool)
    inter = d["intergenic"].astype(bool)
    sel_em = live & isr & (single | gonly) & ~inter  # what _fit_gdna_hyperprior uses (background is not None)
    sel_raw = live & isr & (single | gonly)  # the same WITHOUT the hybrid intergenic drop
    return sel_em, sel_raw


def mwae(d, m):
    """mass-weighted mean |f_g - oracle| over mask m."""
    ok = m & np.isfinite(d["fo"])
    w = d["mass"][ok]
    if w.sum() <= 0:
        return float("nan")
    return float(np.average(np.abs(d["fg"][ok] - d["fo"][ok]), weights=w))


def dens_err(d, m, floor="unit"):
    """mass-weighted |log(rho_hat/rho_oracle)|, rho = f_g*M/E.  M and E are belief-free, so this is
    |log(f_g/f_oracle)| — floored at the 1-count resolution wall g=1 (the KDE's own floor, npmle._kde_density
    `a_cells = log(max(gc,1)) - log(ec)`), which is what a zero-gDNA node actually contributes."""
    ok = m & np.isfinite(d["fo"])
    w = d["mass"][ok]
    if w.sum() <= 0:
        return float("nan")
    M = d["mass"][ok]
    gh = np.maximum(d["fg"][ok] * M, 1.0 if floor == "unit" else _EPS)
    go = np.maximum(d["fo"][ok] * M, 1.0 if floor == "unit" else _EPS)
    return float(np.average(np.abs(np.log(gh) - np.log(go)), weights=w))


def spread(d, m, thr):
    ok = m & np.isfinite(d["fo"])
    e = np.abs(d["fg"][ok] - d["fo"][ok])
    w = d["mass"][ok]
    return float((e > thr).mean()), float(w[e > thr].sum() / max(w.sum(), _EPS))


d0 = D[ARMS[0]]
sel_em0, sel_raw0 = sets(d0)
print("=" * 118)
print("SUBSTRATE CENSUS (arm-independent: the selection uses free_pos/free_neg/signature, not the belief)")
print("=" * 118)
isr = d0["isr"].astype(bool)
n = isr.shape[0]
for lab, m in (
    ("ALL nodes", np.ones(n, bool)),
    ("  REGION", isr),
    ("  BOUNDARY", ~isr),
    ("live & REGION & (single|gonly)", sel_raw0),
    ("   ... minus intergenic  = SEL(EM path)", sel_em0),
    ("   ... of which single", sel_em0 & d0["single"].astype(bool)),
    ("   ... of which gonly", sel_em0 & d0["gonly"].astype(bool)),
    ("   ... of which exon", sel_em0 & (d0["rt"] == 2)),
    ("   ... of which intron", sel_em0 & (d0["rt"] == 1)),
    ("EXCLUDED (AMBIG or boundary)", (~isr) | d0["ambig"].astype(bool)),
    ("intergenic (→ aggregate cell)", d0["intergenic"].astype(bool)),
):
    w = d0["mass"][m]
    print(f"{lab:<42}{m.sum():>9,} nodes ({m.mean():>6.1%})   mass {w.sum():>14,.0f} "
          f"({w.sum() / d0['mass'].sum():>6.1%})")

print()
print("=" * 118)
print("(1)+(2)+(3)  ARM SCORES.  mwae = mass-wtd |f_g-oracle|;  dlog = mass-wtd |log(rho_hat/rho_oracle)|")
print("=" * 118)
hdr = (f"{'population':<40}" + "".join(f"{a:>13}" for a in ARMS))
for metric, fn in (("mwae", mwae), ("dlog(rho_g)", dens_err)):
    print(f"\n-- {metric} --")
    print(hdr)
    for lab, msk in (
        ("SEL(EM path) = THE SUBSTRATE", "sel_em"),
        ("   substrate × exon", "sel_exon"),
        ("   substrate × intron", "sel_intron"),
        ("   substrate × gonly", "sel_gonly"),
        ("SEL_raw (substrate + intergenic)", "sel_raw"),
        ("EXCLUDED (AMBIG | boundary)", "excl"),
        ("ALL nodes (the suite headline)", "all"),
    ):
        row = f"{lab:<40}"
        for a in ARMS:
            d = D[a]
            se, sr = sets(d)
            mm = {
                "sel_em": se, "sel_raw": sr,
                "sel_exon": se & (d["rt"] == 2), "sel_intron": se & (d["rt"] == 1),
                "sel_gonly": se & d["gonly"].astype(bool),
                "excl": (~d["isr"].astype(bool)) | d["ambig"].astype(bool),
                "all": np.ones(d["mass"].shape[0], bool),
            }[msk]
            row += f"{fn(d, mm):>13.4f}"
        print(row)

print("\n-- (3) SPREAD on the substrate: share of nodes / of substrate MASS with |f_g-oracle| over thr --")
print(f"{'':<40}" + "".join(f"{a:>26}" for a in ARMS))
print(f"{'':<40}" + "".join(f"{'nodes>thr':>13}{'mass>thr':>13}" for a in ARMS))
for thr in (0.1, 0.25, 0.5):
    row = f"{'|err| > ' + str(thr):<40}"
    for a in ARMS:
        d = D[a]
        se, _ = sets(d)
        fn_, fm_ = spread(d, se, thr)
        row += f"{fn_:>13.1%}{fm_:>13.1%}"
    print(row)

print("\n-- per-CONDITION substrate mwae (sorted by HEAD arm's substrate error mass) --")
d = D[ARMS[-1]]
se, _ = sets(d)
conds = sorted(set(d["cond"].tolist()))
rows = []
for c in conds:
    cm = (d["cond"] == c) & se & np.isfinite(d["fo"])
    em = float(np.sum(d["mass"][cm] * np.abs(d["fg"][cm] - d["fo"][cm])))
    rows.append((em, c))
rows.sort(reverse=True)
print(f"{'condition':<48}{'sub mass':>12}" + "".join(f"{a:>10}" for a in ARMS) + f"{'ERRmass':>12}")
for em, c in rows:
    line = f"{c[5:]:<48}"
    d1 = D[ARMS[0]]
    se1, _ = sets(d1)
    line += f"{d1['mass'][(d1['cond'] == c) & se1].sum():>12,.0f}"
    for a in ARMS:
        da = D[a]
        sa, _ = sets(da)
        line += f"{mwae(da, sa & (da['cond'] == c)):>10.4f}"
    print(line + f"{em:>12,.0f}")
