"""Independent A/B roll-up of the P1e bench arms (my own rows, /tmp/adv_p1e_bench.tsv)."""

from __future__ import annotations

import csv
import sys

import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else "/tmp/adv_p1e_bench.tsv"
BASE, NEW = (sys.argv[2], sys.argv[3]) if len(sys.argv) > 3 else ("advHEAD_r0", "advP1E_r0")

D: dict = {}
for r in csv.DictReader(open(PATH), delimiter="\t"):
    D.setdefault(r["arm"], {})[r["cond"]] = r


def strata(cond: str):
    out = ["ALL"]
    out.append("capOFF" if "capture_off" in cond else ("capON" if "capture_on" in cond else "vstrong"))
    out.append("ss0.99" if "ss_0.99" in cond else "ss0.50")
    if "ss_0.50" in cond:
        out.append("unstr_x_" + out[1])
    if "ss_0.99" in cond:
        out.append("str_x_" + out[1])
    return out


def mw(arm, sel=None):
    """mass-weighted-ish aggregate: the bench's own convention is a plain mean of per-cond mwae weighted
    by n_nodes*mass? -- use the file's own weight column if present, else plain mean."""
    rows = D[arm]
    vals, wts = [], []
    for c, r in rows.items():
        if sel and sel not in strata(c):
            continue
        vals.append(float(r["mwae_all"]))
        wts.append(float(r.get("mass", 1.0)))
    return np.average(vals, weights=wts), len(vals)


hdr = next(iter(D[BASE].values())).keys()
print("columns:", list(hdr))
print(f"\n{'stratum':<16}{BASE:>14}{NEW:>14}{'delta':>10}{'n':>5}")
for s in ("ALL", "capOFF", "capON", "vstrong", "ss0.99", "ss0.50", "unstr_x_capOFF",
          "unstr_x_capON", "str_x_capOFF", "str_x_capON"):
    try:
        a, n = mw(BASE, s)
        b, _ = mw(NEW, s)
    except ZeroDivisionError:
        continue
    if n == 0:
        continue
    print(f"{s:<16}{a:>14.4f}{b:>14.4f}{b - a:>+10.4f}{n:>5}")

bet = wor = flat = 0
lines = []
for c in sorted(D[BASE]):
    a, b = float(D[BASE][c]["mwae_all"]), float(D[NEW][c]["mwae_all"])
    if abs(b - a) <= 0.002:
        flat += 1
    elif b < a - 0.002:
        bet += 1
    elif b > a + 0.002:
        wor += 1
    lines.append((b - a, c, a, b))
print(f"\nscenario counts {BASE} -> {NEW}: {bet} better / {wor} worse / {flat} flat")
lines.sort()
print("\n  biggest movers")
for d, c, a, b in lines[:6] + lines[-6:]:
    print(f"    {c:<50}{a:>9.4f}{b:>9.4f}{d:>+10.4f}  ({100 * d / max(a, 1e-9):+6.1f}%)")
