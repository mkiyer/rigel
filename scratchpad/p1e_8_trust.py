"""P1e (4) — the TRUST view A/B from two `pass0_error_table.py` npz dumps.

Reports, per node class: error reads, confidently-wrong reads, and z2 RESTRICTED TO THE CONFIDENT QUARTILE
(`z2|Q1` — the number the ROADMAP and PASS0_FINISH_PLAN quote). Both with each run's OWN quartile threshold
(the shipped convention) and with HEAD's threshold held FIXED, so the comparison is not confounded by the
threshold moving.

    python scratchpad/p1e_8_trust.py /tmp/p1e_state_head.npz /tmp/p1e_state_on.npz
"""

from __future__ import annotations

import sys

import numpy as np

_EPS = 1e-9
A, B = sys.argv[1], sys.argv[2]


def load(p):
    d = dict(np.load(p, allow_pickle=True))
    d["amb"] = d["amb"].astype(bool)
    d["raw"] = np.where(d["mass"] > _EPS, d["err"] / np.maximum(d["mass"], _EPS), 0.0)
    return d


def z2(d, m):
    v, mass, raw = d["var"], d["mass"], d["raw"]
    k = m & np.isfinite(v)
    if not k.any():
        return float("nan")
    den = float(np.sum(mass[k] * v[k]))
    return float(np.sum(mass[k] * raw[k] ** 2)) / den if den > 0 else float("nan")


da, db = load(A), load(B)
qa = np.quantile(da["var"][np.isfinite(da["var"])], 0.25)
qb = np.quantile(db["var"][np.isfinite(db["var"])], 0.25)
print(f"  confident-quartile threshold on Var(log f_g):  HEAD {qa:.5g}   P1e {qb:.5g}")
for thr_name, ta, tb in (("OWN quartile", qa, qb), ("HEAD threshold FIXED", qa, qa)):
    ca = np.isfinite(da["var"]) & (da["var"] <= ta)
    cb = np.isfinite(db["var"]) & (db["var"] <= tb)
    print(f"\n{'=' * 118}\n  {thr_name}\n{'=' * 118}")
    print(f"  {'class':<18}{'ERR head':>13}{'ERR p1e':>13}{'d%':>8}|"
          f"{'CWRONG head':>13}{'CWRONG p1e':>13}{'d%':>8}|{'z2|Q1 head':>12}{'z2|Q1 p1e':>12}")
    rows = [("exon", False), ("exon", True), ("boundary", False), ("boundary", True),
            ("intron", False), ("intron", True), ("intergenic", False)]
    for cls, amb in rows:
        ma = (da["cls"] == cls) & (da["amb"] == amb)
        mb = (db["cls"] == cls) & (db["amb"] == amb)
        if not ma.any():
            continue
        ea, eb = da["err"][ma].sum(), db["err"][mb].sum()
        wa, wb = da["err"][ma & ca].sum(), db["err"][mb & cb].sum()
        print(f"  {cls + (' AMBIG' if amb else ' single'):<18}{ea:>13,.0f}{eb:>13,.0f}"
              f"{100 * (eb - ea) / max(ea, _EPS):>+7.1f}%|{wa:>13,.0f}{wb:>13,.0f}"
              f"{100 * (wb - wa) / max(wa, _EPS):>+7.1f}%|{z2(da, ma & ca):>12.2f}{z2(db, mb & cb):>12.2f}")
    ea, eb = da["err"].sum(), db["err"].sum()
    wa, wb = da["err"][ca].sum(), db["err"][cb].sum()
    print(f"  {'ALL':<18}{ea:>13,.0f}{eb:>13,.0f}{100 * (eb - ea) / ea:>+7.1f}%|"
          f"{wa:>13,.0f}{wb:>13,.0f}{100 * (wb - wa) / wa:>+7.1f}%|"
          f"{z2(da, ca):>12.2f}{z2(db, cb):>12.2f}")
    print(f"  {'mwae (mass-wtd)':<18}{ea / da['mass'].sum():>13.4f}{eb / db['mass'].sum():>13.4f}")

# and by regime, the axis that decides
print(f"\n{'=' * 118}\n  BY REGIME (own-quartile convention)\n{'=' * 118}")
ca = np.isfinite(da["var"]) & (da["var"] <= qa)
cb = np.isfinite(db["var"]) & (db["var"] <= qb)
AX = {
    "stranded x capOFF": lambda c: "ss_0.99" in c and "capture_off" in c,
    "stranded x capON": lambda c: "ss_0.99" in c and "capture_on" in c,
    "unstranded x capOFF": lambda c: "ss_0.50" in c and "capture_off" in c,
    "unstranded x capON": lambda c: "ss_0.50" in c and "capture_on" in c,
    "unstranded x verystrong": lambda c: "ss_0.50" in c and "verystrong" in c,
}
print(f"  {'regime':<26}{'ERR head':>13}{'ERR p1e':>13}{'d%':>8}|"
      f"{'CWRONG head':>13}{'CWRONG p1e':>13}{'d%':>8}|{'z2|Q1 head':>12}{'z2|Q1 p1e':>12}")
for lab, f in AX.items():
    ma = np.array([f(c) for c in da["cond"]])
    mb = np.array([f(c) for c in db["cond"]])
    ea, eb = da["err"][ma].sum(), db["err"][mb].sum()
    wa, wb = da["err"][ma & ca].sum(), db["err"][mb & cb].sum()
    print(f"  {lab:<26}{ea:>13,.0f}{eb:>13,.0f}{100 * (eb - ea) / max(ea, _EPS):>+7.1f}%|"
          f"{wa:>13,.0f}{wb:>13,.0f}{100 * (wb - wa) / max(wa, _EPS):>+7.1f}%|"
          f"{z2(da, ma & ca):>12.2f}{z2(db, mb & cb):>12.2f}")
