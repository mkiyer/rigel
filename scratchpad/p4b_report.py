"""P4b — price the `_RHO_ITERS = 2` BP violation. REPORT pass.

    python scratchpad/p4b_report.py                 # the full table
    python scratchpad/p4b_report.py --decomp        # + the "what does iteration 2 buy" decomposition

Node sets are HELD FIXED across arms: the trust view is evaluated on the 2-iteration arm's OWN confident
quartile (by stated Var(log f_g)), applied unchanged to the 1-iteration arm, so no re-selection can hide the
effect. Row alignment ((cond, node-id) identical across arms) is asserted, not assumed.
"""

from __future__ import annotations

import argparse
import numpy as np

_EPS = 1e-9
ap = argparse.ArgumentParser()
ap.add_argument("--decomp", action="store_true")
a = ap.parse_args()

D = {}
for arm in ("rho2", "rho1"):
    for r in (0, 1):
        D[(arm, r)] = dict(np.load(f"/tmp/p4b_{arm}_r{r}.npz", allow_pickle=True))

base = D[("rho2", 0)]
for k, v in D.items():
    assert (v["cond"] == base["cond"]).all() and (v["nid"] == base["nid"]).all(), k
    assert np.allclose(v["mass"], base["mass"]) and np.allclose(v["fo"], base["fo"], equal_nan=True), k
print(f"row alignment OK: {len(base['cond']):,} node-rows x 4 arms, identical (cond, node-id, mass, oracle)")

cond, mass, fo = base["cond"], base["mass"], base["fo"]
cls, amb, solv, fit = base["cls"], base["amb"].astype(bool), base["solv"].astype(bool), base["fit"].astype(bool)
conds = sorted(set(cond.tolist()))

STRATA = {
    "ALL 32": lambda c: True,
    "stranded ss_0.99": lambda c: "ss_0.99" in c,
    "unstranded x capON": lambda c: "ss_0.50" in c and "capture_on" in c,
    "capture OFF": lambda c: "capture_off" in c,
    "capture ON": lambda c: "capture_on" in c,
    "verystrong": lambda c: "verystrong" in c,
    "gdna_none": lambda c: c.startswith("gdna_none"),
    "FIT SUBSTRATE (reg x single|gonly)": None,  # a node mask, not a condition mask
}


def mwae(m, arm, r):
    w = mass[m]
    return float(np.average(np.abs(D[(arm, r)]["fg"][m] - fo[m]), weights=w)) if w.sum() > 0 else np.nan


def cond_mask(f):
    return np.array([f(c) for c in cond])


print("\n" + "=" * 104)
print("(1) THE PRICE — mass-weighted mean |f_g − oracle| (mwae), pooled over the stratum's nodes")
print("=" * 104)
hdr = f"{'stratum':<36}{'r0 rho=2':>10}{'r0 rho=1':>10}{'Δ':>9}{'  ':>2}{'r1 rho=2':>10}{'r1 rho=1':>10}{'Δ':>9}"
print(hdr)
print("-" * len(hdr))
for name, f in STRATA.items():
    m = fit if f is None else cond_mask(f)
    if not m.any():
        continue
    a0, b0, a1, b1 = mwae(m, "rho2", 0), mwae(m, "rho1", 0), mwae(m, "rho2", 1), mwae(m, "rho1", 1)
    print(f"{name:<36}{a0:>10.4f}{b0:>10.4f}{b0 - a0:>+9.4f}  {a1:>10.4f}{b1:>10.4f}{b1 - a1:>+9.4f}")

# per-condition better/worse, at both refits
print("\nper-condition (mwae, mass-weighted within condition): rho=1 vs rho=2")
for r in (0, 1):
    bet = wor = flat = 0
    worst = []
    for c in conds:
        m = cond == c
        x, y = mwae(m, "rho2", r), mwae(m, "rho1", r)
        if abs(y - x) < 5e-5:
            flat += 1
        elif y < x:
            bet += 1
        else:
            wor += 1
        worst.append((y - x, c, x, y))
    worst.sort()
    print(f"  refit={r}:  rho=1 better {bet} / worse {wor} / flat {flat}")
    for d, c, x, y in worst[:3]:
        print(f"      best  {c[5:]:<46} {x:.4f} -> {y:.4f}  ({d:+.4f})")
    for d, c, x, y in worst[-3:][::-1]:
        print(f"      worst {c[5:]:<46} {x:.4f} -> {y:.4f}  ({d:+.4f})")

print("\n" + "=" * 104)
print("(1b) THE TRUST VIEW — held-FIXED node set: the rho=2 arm's own confident quartile (same refit)")
print("=" * 104)
for r in (0, 1):
    var2 = D[("rho2", r)]["var"]
    fin = np.isfinite(var2) & solv
    q1 = float(np.quantile(var2[fin], 0.25))
    conf = fin & (var2 <= q1)  # HELD FIXED: selected by the rho=2 arm, applied to BOTH arms
    print(f"\n  refit={r}   fixed threshold Var(log f_g) <= {q1:.5g}   nodes in set: {conf.sum():,}")
    rows = [("ALL (solvable)", solv), ("FIT SUBSTRATE", fit & solv)]
    for c in ("exon", "intron", "boundary", "intergenic"):
        for lab, mm in ((" single", ~amb), (" AMBIG", amb)):
            rows.append((c + lab, (cls == c) & mm & solv))
    h = (f"  {'population':<26}{'CW rho=2':>12}{'CW rho=1':>12}{'Δ%':>8}"
         f"{'z2 rho=2':>10}{'z2 rho=1':>10}{'Δ%':>8}")
    print(h)
    print("  " + "-" * (len(h) - 2))
    for lab, m in rows:
        if not (m & conf).any():
            continue
        out = []
        for arm in ("rho2", "rho1"):
            d = D[(arm, r)]
            e = np.abs(d["fg"] - fo) * mass
            k = m & conf
            cw = float(e[k].sum())
            den = float(np.sum(mass[k] * d["var"][k]))
            num = float(np.sum(mass[k] * (e[k] / np.maximum(mass[k], _EPS)) ** 2))
            out.append((cw, num / den if den > 0 else np.nan))
        (c2, z2), (c1, z1) = out
        print(f"  {lab:<26}{c2:>12,.0f}{c1:>12,.0f}{(c1 / max(c2, 1) - 1) * 100:>+8.1f}"
              f"{z2:>10.2f}{z1:>10.2f}{(z1 / max(z2, _EPS) - 1) * 100:>+8.1f}")

if a.decomp:
    print("\n" + "=" * 104)
    print("(2) WHAT THE SECOND ITERATION BUYS — decomposition")
    print("=" * 104)
    dlr = D[("rho2", 0)]["dlr"]  # |Δ log ρ_face| between the two lazy-ρ iterations
    for r in (0,):
        e2 = np.abs(D[("rho2", r)]["fg"] - fo) * mass
        e1 = np.abs(D[("rho1", r)]["fg"] - fo) * mass
        dfg = np.abs(D[("rho2", r)]["fg"] - D[("rho1", r)]["fg"])
        print(f"\n  refit={r}:  Σ|Δf_g| mass = {float((dfg * mass).sum()):,.0f} reads "
              f"({float((dfg * mass).sum()) / mass.sum():.2%} of node mass); "
              f"nodes moved >0.01: {float((dfg > 0.01).mean()):.1%}")
        print(f"  {'stratum':<36}{'med|Δlogρ|':>12}{'>1% nodes':>11}{'Σmass·|Δf_g|':>14}"
              f"{'ERR rho2':>12}{'ERR rho1':>12}{'Δ ERR':>11}")
        for name, f in list(STRATA.items()):
            m = fit if f is None else cond_mask(f)
            if not m.any():
                continue
            print(f"  {name:<36}{np.median(dlr[m]):>12.4f}{float((dlr[m] > 0.01).mean()):>11.1%}"
                  f"{float((dfg[m] * mass[m]).sum()):>14,.0f}{e2[m].sum():>12,.0f}{e1[m].sum():>12,.0f}"
                  f"{e1[m].sum() - e2[m].sum():>+11,.0f}")
        # is the gain where the reframe changes most?  bin by |Δ log ρ_face|
        print(f"\n  by |Δ log ρ_face| decile (refit={r}) — does the gain track where the reframe moves?")
        q = np.quantile(dlr, np.linspace(0, 1, 11))
        print(f"  {'decile':<10}{'|Δlogρ| hi':>12}{'mass':>14}{'ERR rho2':>12}{'ERR rho1':>12}"
              f"{'Δ ERR':>11}{'Σm|Δf_g|':>12}")
        for i in range(10):
            m = (dlr >= q[i]) & (dlr <= q[i + 1] if i == 9 else dlr < q[i + 1])
            if not m.any():
                continue
            print(f"  {i + 1:<10}{q[i + 1]:>12.4f}{mass[m].sum():>14,.0f}{e2[m].sum():>12,.0f}"
                  f"{e1[m].sum():>12,.0f}{e1[m].sum() - e2[m].sum():>+11,.0f}"
                  f"{float((dfg[m] * mass[m]).sum()):>12,.0f}")
        # MODE vs PRECISION
        print("\n  MODE vs PRECISION: does iteration 2 change the answer, or only the stated variance?")
        v2, v1 = D[("rho2", r)]["var"], D[("rho1", r)]["var"]
        fin = np.isfinite(v2) & np.isfinite(v1) & solv
        print(f"    median stated Var(log f_g):  rho=2 {np.median(v2[fin]):.5g}   rho=1 {np.median(v1[fin]):.5g}"
              f"   median ratio {np.median(v1[fin] / np.maximum(v2[fin], _EPS)):.4f}")
        print(f"    mass-weighted mean |Δf_g| (MODE)  = {float((dfg * mass).sum() / mass.sum()):.5f}")
        print(f"    mass-wtd mean |Δ log Var| (PREC)  = "
              f"{float(np.average(np.abs(np.log(np.maximum(v1[fin], _EPS) / np.maximum(v2[fin], _EPS))), weights=mass[fin])):.5f}")
        # fit substrate vs excluded
        print("\n  FIT SUBSTRATE vs EXCLUDED (the population the hyperprior never sees):")
        for lab, m in (("fit substrate", fit & solv), ("EXCLUDED (AMBIG+bnd)", ~fit & solv)):
            print(f"    {lab:<24} mass {mass[m].sum():>13,.0f}  ERR rho2 {e2[m].sum():>11,.0f}"
                  f"  ERR rho1 {e1[m].sum():>11,.0f}  Δ {e1[m].sum() - e2[m].sum():>+11,.0f}"
                  f"  ({(e1[m].sum() / max(e2[m].sum(), 1) - 1) * 100:+.1f}%)")
