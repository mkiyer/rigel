"""W4 — the A/B report. Reads the bench TSV and reports every gate the plan names.

Arms, all recorded from HEAD in one session (HANDOFF_16 §0: a stale baseline reads as a broken change):
    base_*  HEAD before W4            — replacing `_gdna_arm`, DensityNPMLE hyperprior
    add_*   + the ADDITIVE `_gdna_arm` (owner: KEEP)
    land_*  + the GdnaLandscape hyperprior   <- what ships

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_w4_report.py [tsv]
"""
from __future__ import annotations

import csv
import sys

import numpy as np

TSV = sys.argv[1] if len(sys.argv) > 1 else "/tmp/pass0_oracle_bench_n4.tsv"
STRATA = [
    ("ALL 32", lambda c: True),
    ("stranded ss_0.99", lambda c: "ss_0.99" in c),
    ("unstranded ss_0.50", lambda c: "ss_0.50" in c),
    ("capture ON", lambda c: "capture_on" in c),
    ("capture OFF", lambda c: "capture_off" in c),
    ("capture VERYSTRONG", lambda c: "verystrong" in c),
    ("unstranded x capON", lambda c: "ss_0.50" in c and "capture_on" in c),
    ("gdna_none (FP guard)", lambda c: c.startswith("gdna_none_")),
]


def load():
    d = {}
    with open(TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[(r["arm"], r["cond"])] = r
    return d, sorted({k[1] for k in d})


def compare(d, conds, a, b, title):
    rows = [(c, float(d[(a, c)]["mwae_all"]), float(d[(b, c)]["mwae_all"]), float(d[(a, c)]["mass"]))
            for c in conds if (a, c) in d and (b, c) in d]
    if not rows:
        return
    ident = sum(1 for _, x, y, _ in rows if x == y)
    print(f"\n=== {title}   [{a} -> {b}] ===")
    print(f"  bit-identical {ident}/{len(rows)}")
    print(f"  {'stratum':22s} {'n':>3s} {'base':>8s} {'new':>8s} {'delta':>9s}  {'b/w/flat':>12s}")
    for name, sel in STRATA:
        r = [x for x in rows if sel(x[0])]
        if not r:
            continue
        arr = np.array([[x[1], x[2]] for x in r])
        w = np.array([x[3] for x in r])
        bt = sum(1 for _, x, y, _ in r if y < x - 0.002)
        ws = sum(1 for _, x, y, _ in r if y > x + 0.002)
        print(f"  {name:22s} {len(r):3d} {np.average(arr[:, 0], weights=w):8.4f} "
              f"{np.average(arr[:, 1], weights=w):8.4f} "
              f"{np.average(arr[:, 1] - arr[:, 0], weights=w):+9.4f}  "
              f"{bt:>3d}/{ws:>3d}/{len(r) - bt - ws:>3d}")
    return rows


def main():
    d, conds = load()
    arms = sorted({k[0] for k in d})
    print(f"arms present: {', '.join(arms)}")
    for a, b, t in (("base_r0", "land_r0", "refit=0 — THE WIRING GATE (prior is None at pass-0)"),
                    ("base_r1", "land_r1", "refit=1 — W4 END TO END (vs HEAD before any of this)"),
                    ("add_r1", "land_r1", "refit=1 — the LANDSCAPE alone (vs the additive arm)")):
        if any((x, conds[0]) in d for x in (a, b)):
            rows = compare(d, conds, a, b, t)
            if rows and b == "land_r1" and a == "base_r1":
                print("\n  per condition:")
                for c, x, y, _ in rows:
                    mk = "  BETTER" if y < x - 0.002 else ("  worse" if y > x + 0.002 else "")
                    print(f"    {c[5:]:50s} {x:.4f} -> {y:.4f}  {y - x:+.4f}{mk}")


if __name__ == "__main__":
    main()
