"""Adversarial verification of P4: aggregate the arms, mass-weighted, exactly as pass0_oracle_bench --report."""
from __future__ import annotations

import csv
import sys

import numpy as np

path = sys.argv[1]
arms = sys.argv[2:]
d = {}
with open(path) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        d[(r["arm"], r["cond"])] = r
conds = sorted({k[1] for k in d})

STRATA = {
    "ALL 32": lambda c: True,
    "stranded ss_0.99": lambda c: "ss_0.99" in c,
    "unstranded ss_0.50": lambda c: "ss_0.50" in c,
    "unstranded x capON": lambda c: "ss_0.50" in c and "capture_on" in c,
    "capture OFF": lambda c: "capture_off" in c,
    "capture ON": lambda c: "capture_on" in c,
    "verystrong": lambda c: "verystrong" in c,
    "gdna_none": lambda c: c.startswith("gdna_none"),
}

base = arms[0]
print(f"{'stratum':<22}{'n':>4}" + "".join(f"{a:>12}" for a in arms))
for lab, f in STRATA.items():
    cs = [c for c in conds if f(c) and all((a, c) in d for a in arms)]
    w = np.array([float(d[(base, c)]["mass"]) for c in cs])
    vals = [np.average([float(d[(a, c)]["mwae_all"]) for c in cs], weights=w) for a in arms]
    print(f"{lab:<22}{len(cs):>4}" + "".join(f"{v:>12.4f}" for v in vals))

print("\nbetter / worse / flat vs", base, " (threshold 0.002)")
for a in arms[1:]:
    cs = [c for c in conds if (a, c) in d and (base, c) in d]
    eb = np.array([float(d[(base, c)]["mwae_all"]) for c in cs])
    ep = np.array([float(d[(a, c)]["mwae_all"]) for c in cs])
    print(f"  {a:<14} {(ep < eb - 0.002).sum()} better, {(ep > eb + 0.002).sum()} worse, "
          f"{(np.abs(ep - eb) <= 0.002).sum()} flat   maxdelta={np.max(np.abs(ep - eb)):.4f}")

print("\nBIT-IDENTITY check (all numeric fields, exact string equality):")
for a in arms[1:]:
    ident = True
    diffs = 0
    for c in conds:
        rb, rp = d.get((base, c)), d.get((a, c))
        if not rb or not rp:
            continue
        for k in rb:
            if k == "arm":
                continue
            if rb[k] != rp[k]:
                ident = False
                diffs += 1
    print(f"  {a:<14} identical={ident}  differing_fields={diffs}")

print("\nper-condition (sorted by |delta| of arm2 vs base):")
if len(arms) > 1:
    cs = [c for c in conds if all((a, c) in d for a in arms)]
    rowsx = sorted(cs, key=lambda c: -abs(float(d[(arms[1], c)]["mwae_all"]) - float(d[(base, c)]["mwae_all"])))
    print(f"{'cond':<46}" + "".join(f"{a:>11}" for a in arms))
    for c in rowsx[:12]:
        print(f"{c[5:]:<46}" + "".join(f"{float(d[(a, c)]['mwae_all']):>11.4f}" for a in arms))
