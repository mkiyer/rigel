"""Summarize /tmp/pass0_oracle_bench.tsv arms: mass-weighted mwae over strata.

    python scratchpad/ho_summary.py arm1 arm2 ...
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

OUT = Path("/tmp/pass0_oracle_bench.tsv")
rows: dict = {}
with OUT.open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows.setdefault(r["arm"], {})[r["cond"]] = r

arms = sys.argv[1:] or sorted(rows)
STRATA = {
    "ALL 32": lambda c: True,
    "capOFF": lambda c: "capture_off" in c,
    "capON": lambda c: "capture_on" in c,
    "verystrong": lambda c: "verystrong" in c,
    "stranded .99": lambda c: "ss_0.99" in c,
    "unstranded .50": lambda c: "ss_0.50" in c,
    "unstr x capON": lambda c: "ss_0.50" in c and "capture_on" in c,
    "gdna_none": lambda c: "gdna_none" in c,
}
hdr = f"{'arm':<16}" + "".join(f"{k:>16}" for k in STRATA)
print(hdr)
print("-" * len(hdr))
for a in arms:
    d = rows.get(a)
    if not d:
        print(f"{a:<16}  (missing)")
        continue
    line = f"{a:<16}"
    for _k, sel in STRATA.items():
        cs = [c for c in d if sel(c)]
        if not cs:
            line += f"{'-':>16}"
            continue
        e = np.array([float(d[c]["mwae_all"]) for c in cs])
        w = np.array([float(d[c]["mass"]) for c in cs])
        line += f"{np.average(e, weights=w):>16.4f}"
    print(line)

# pairwise better/worse against the first arm
if len(arms) > 1:
    base = rows.get(arms[0], {})
    print(f"\nvs {arms[0]}:  better / worse / flat  (|d|>0.002), and worst regression")
    for a in arms[1:]:
        d = rows.get(a, {})
        cs = sorted(set(base) & set(d))
        if not cs:
            continue
        db = np.array([float(d[c]["mwae_all"]) - float(base[c]["mwae_all"]) for c in cs])
        w = np.array([float(base[c]["mass"]) for c in cs])
        j = int(np.argmax(db))
        print(f"  {a:<16} {int((db < -0.002).sum()):>3} / {int((db > 0.002).sum()):>3} /"
              f" {int((abs(db) <= 0.002).sum()):>3}   Δmw={np.average(db, weights=w):+.4f}"
              f"   worst {cs[j][5:]:<44} {db[j]:+.4f}")
