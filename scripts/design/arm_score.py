#!/usr/bin/env python
"""Score two ladder arms against each other, on the DELIVERABLE and on pass-0, per stratum.

⛔ `abs_err_all_final` is the deliverable (TRAPS: the-intermediate-is-not-the-deliverable — a −37.2 % pass-0 win was −3.9 % shipped);
`abs_err_all` is pass-0. BOTH are printed on every row, always.
⛔ `g00` is reported SEPARATELY, never inside the ALL row: its truth is exactly 0, so its relative
change is unbounded and would swamp the other 35 conditions in either direction.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import os

#: where `ladder_arm_ab.py --out` wrote the arms. Override with $RIGEL_ARMS.
D = Path(os.environ.get("RIGEL_ARMS", Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_arms"))


def load(name):
    rows = [json.loads(x) for x in (D / f"{name}.jsonl").read_text().splitlines() if x.strip()]
    return {(r["condition"], r["axis"]): r for r in rows}


def stratum(cond):
    return ("stranded" if "ss_0.99" in cond else "unstranded",
            "capture ON" if "capture_on" in cond else "capture OFF")


def main():
    a_name, b_name = sys.argv[1], sys.argv[2]
    A, B = load(a_name), load(b_name)
    keys = sorted(set(A) & set(B))
    assert len(keys) == len(A) == len(B), f"arms differ in shape: {len(A)} vs {len(B)} vs {len(keys)}"

    def agg(field, pred):
        sa = sum(A[k][field] for k in keys if pred(k))
        sb = sum(B[k][field] for k in keys if pred(k))
        return sa, sb

    def line(label, pred, n):
        out = []
        for field in ("abs_err_all_final", "abs_err_all"):
            sa, sb = agg(field, pred)
            pct = (sb - sa) / sa * 100 if sa > 0 else float("nan")
            out.append((sa, sb, pct))
        better = sum(1 for k in keys if pred(k)
                     and B[k]["abs_err_all_final"] < A[k]["abs_err_all_final"])
        (fa, fb, fp), (pa, pb, pp) = out
        print(f"   {label:<26} {fa:>13,.0f} {fb:>13,.0f} {fp:>+8.1f}%   "
              f"{pa:>13,.0f} {pb:>13,.0f} {pp:>+8.1f}%   {better:>3}/{n}")

    def is_g00(k):
        return "_g00_" in k[0]

    print()
    print(f"   {a_name}  vs  {b_name}      36 conditions x 2 axes")
    print(f"   {'':<26} {'':>13} {'DELIVERABLE':>13} {'':>9}   {'':>13} {'PASS-0':>13}")
    print(f"   {'stratum':<26} {'base':>13} {'arm':>13} {'':>9}   {'base':>13} {'arm':>13} "
          f"{'':>9}  better")
    print("   " + "-" * 112)
    n_ex = sum(1 for k in keys if not is_g00(k))
    line("ALL (g00 excluded)", lambda k: not is_g00(k), n_ex)
    def _in(k, s):
        return stratum(k[0]) == s and not is_g00(k)

    for s in (("stranded", "capture ON"), ("stranded", "capture OFF"),
              ("unstranded", "capture ON"), ("unstranded", "capture OFF")):
        def pred(k, s=s):
            return _in(k, s)
        line(" x ".join(s), pred, sum(1 for k in keys if pred(k)))
    print("   " + "-" * 112)
    for ax in ("region", "edge"):
        def pred(k, ax=ax):
            return k[1] == ax and not is_g00(k)
        line(f"axis: {ax}", pred, sum(1 for k in keys if pred(k)))
    print("   " + "-" * 112)
    line("⛔ g00 ZERO CONTROL", is_g00, sum(1 for k in keys if is_g00(k)))

    # the library-level deliverable at the two controls, where truth is a constant
    print()
    print("   ⭐ library f_gdna at the CONTROLS (truth is a constant, so every deviation is signed)")
    print(f"      {'condition':<44} {'truth':>8} {'base':>8} {'arm':>8}")
    for cond in sorted({k[0] for k in keys}):
        if "_g00_" in cond or "_g98_" in cond:
            k = (cond, "region")
            print(f"      {cond:<44} {A[k]['library_f_gdna_truth']:>8.4f} "
                  f"{A[k]['library_f_gdna_final']:>8.4f} {B[k]['library_f_gdna_final']:>8.4f}")

    # the worst regressions and the best wins, on the deliverable
    delta = sorted(
        ((B[k]["abs_err_all_final"] - A[k]["abs_err_all_final"], k) for k in keys if not is_g00(k)),
    )
    print()
    print("   ⭐ the 6 biggest WINS and the 6 biggest REGRESSIONS on the deliverable (fragments)")
    for tag, rows in (("win", delta[:6]), ("regression", delta[-6:][::-1])):
        for d, k in rows:
            base = A[k]["abs_err_all_final"]
            print(f"      {tag:<11} {k[0]:<44} {k[1]:<5} {base:>12,.0f} -> "
                  f"{B[k]['abs_err_all_final']:>12,.0f}  {d / base * 100 if base else 0:>+8.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
