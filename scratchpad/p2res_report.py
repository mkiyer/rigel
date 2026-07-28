"""The P-2 RESIDUAL A/B report — generic arm-pair comparison over the 32-condition battery.

Same TSV as `gdna_w4_report.py`, but the arm pair is given on the command line and the strata include the
ones the P-2 residual is defined on (`pin_derivation.md` §11): unstranded × capture-OFF × gDNA-bearing.

Run: OMP_NUM_THREADS=1 python scratchpad/p2res_report.py <tsv> <base_arm> <new_arm> [--per-cond]
"""
from __future__ import annotations

import csv
import sys

import numpy as np

STRATA = [
    ("ALL 32", lambda c: True),
    ("stranded ss_0.99", lambda c: "ss_0.99" in c),
    ("unstranded ss_0.50", lambda c: "ss_0.50" in c),
    ("capture ON", lambda c: "capture_on" in c),
    ("capture OFF", lambda c: "capture_off" in c),
    ("capture VERYSTRONG", lambda c: "verystrong" in c),
    ("unstranded x capON", lambda c: "ss_0.50" in c and "capture_on" in c),
    ("** unstr x capOFF x gDNA", lambda c: "ss_0.50" in c and "capture_off" in c
     and not c.startswith("gdna_none_")),
    ("stranded x capOFF", lambda c: "ss_0.99" in c and "capture_off" in c),
    ("gdna_none (FP guard)", lambda c: c.startswith("gdna_none_")),
]


def load(tsv):
    d = {}
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[(r["arm"], r["cond"])] = r
    return d, sorted({k[1] for k in d})


def main():
    tsv = sys.argv[1]
    a, b = sys.argv[2], sys.argv[3]
    per_cond = "--per-cond" in sys.argv
    d, conds = load(tsv)
    rows = [(c, float(d[(a, c)]["mwae_all"]), float(d[(b, c)]["mwae_all"]), float(d[(a, c)]["mass"]),
             float(d[(a, c)]["corr_all"]), float(d[(b, c)]["corr_all"]))
            for c in conds if (a, c) in d and (b, c) in d]
    if not rows:
        print(f"no overlap for {a} / {b}; arms present: {sorted({k[0] for k in d})}")
        return
    ident = sum(1 for r in rows if r[1] == r[2])
    print(f"=== [{a} -> {b}]  bit-identical {ident}/{len(rows)} ===")
    print(f"  {'stratum':24s} {'n':>3s} {'base':>8s} {'new':>8s} {'delta':>9s}  {'b/w/flat':>10s}  {'dcorr':>7s}")
    for name, sel in STRATA:
        r = [x for x in rows if sel(x[0])]
        if not r:
            continue
        arr = np.array([[x[1], x[2]] for x in r])
        w = np.array([x[3] for x in r])
        cc = np.array([[x[4], x[5]] for x in r])
        ok = np.isfinite(cc).all(axis=1)
        dc = np.average(cc[ok, 1] - cc[ok, 0], weights=w[ok]) if ok.any() else float("nan")
        bt = sum(1 for x in r if x[2] < x[1] - 0.002)
        ws = sum(1 for x in r if x[2] > x[1] + 0.002)
        print(f"  {name:24s} {len(r):3d} {np.average(arr[:, 0], weights=w):8.4f} "
              f"{np.average(arr[:, 1], weights=w):8.4f} "
              f"{np.average(arr[:, 1] - arr[:, 0], weights=w):+9.4f}  "
              f"{bt:>3d}/{ws:>3d}/{len(r) - bt - ws:>3d}  {dc:+7.4f}")
    if per_cond:
        print("\n  per condition:")
        for c, x, y, _, _, _ in sorted(rows, key=lambda t: t[2] - t[1]):
            mk = "  BETTER" if y < x - 0.002 else ("  worse" if y > x + 0.002 else "")
            print(f"    {c[5:]:52s} {x:.4f} -> {y:.4f}  {y - x:+.4f}{mk}")


if __name__ == "__main__":
    main()
