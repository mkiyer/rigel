#!/usr/bin/env python
"""The message-precision sweep, per stratum — is there a PLATEAU?

⛔ Read the STRATA, never the total. the panel total has hidden a sign flip between strata (`TRAPS: never-pool-the-strata`)
between them, so a pooled curve would average the two regimes the sweep exists to separate.
"""

from __future__ import annotations

import json
from pathlib import Path

import os

#: where `ladder_arm_ab.py --out` wrote the arms. Override with $RIGEL_ARMS.
D = Path(os.environ.get("RIGEL_ARMS", Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_arms"))

ARMS = [("msgfree_all", 0.0), ("msgscale_0.001", 0.001), ("msgscale_0.01", 0.01),
        ("msgscale_0.1", 0.1), ("msgscale_0.5", 0.5), ("base", 1.0)]

STRATA = [("stranded x cap ON", "ss_0.99", "capture_on"),
          ("stranded x cap OFF", "ss_0.99", "capture_off"),
          ("unstranded x cap ON", "ss_0.50", "capture_on"),
          ("unstranded x cap OFF", "ss_0.50", "capture_off")]


def load(name):
    f = D / f"{name}.jsonl"
    if not f.is_file():
        return None
    return [json.loads(x) for x in f.read_text().splitlines() if x.strip()]


def main():
    rows = {n: load(n) for n, _ in ARMS}
    have = [(n, k) for n, k in ARMS if rows[n]]

    def total(name, ss, cap, field="abs_err_all_final", g00=False):
        return sum(r[field] for r in rows[name]
                   if ss in r["condition"] and cap in r["condition"]
                   and (("_g00_" in r["condition"]) == g00))

    print()
    print("   ⭐⭐ MESSAGE-PRECISION SWEEP — Σ|err| on the DELIVERABLE, per stratum, g00 excluded")
    print("   (scale multiplies ALL FOUR psi message precisions; 1.0 = base, 0.0 = msgfree_all)")
    print()
    hdr = f"   {'stratum':<22}" + "".join(f"{k:>12g}" for _, k in have)
    print(hdr)
    print("   " + "-" * (len(hdr) - 3))
    for label, ss, cap in STRATA:
        vals = [total(n, ss, cap) for n, _ in have]
        b = total("base", ss, cap)
        print(f"   {label:<22}" + "".join(f"{v / b * 100 - 100:>+11.1f}%" for v in vals))
    print("   " + "-" * (len(hdr) - 3))
    allv = [sum(total(n, ss, cap) for _, ss, cap in STRATA) for n, _ in have]
    ab = sum(total("base", ss, cap) for _, ss, cap in STRATA)
    print(f"   {'ALL (g00 excluded)':<22}" + "".join(f"{v / ab * 100 - 100:>+11.1f}%" for v in allv))
    g00 = [sum(total(n, ss, cap, g00=True) for _, ss, cap in STRATA) for n, _ in have]
    gb = sum(total("base", ss, cap, g00=True) for _, ss, cap in STRATA)
    print(f"   {'⛔ g00 zero control':<22}" + "".join(f"{v / gb * 100 - 100:>+11.0f}%" for v in g00))

    print()
    print("   ⭐ ABSOLUTE Σ|err| (fragments), so the trade is visible in mass not just in percent")
    print(f"   {'stratum':<22}" + "".join(f"{k:>12g}" for _, k in have))
    for label, ss, cap in STRATA:
        print(f"   {label:<22}" + "".join(f"{total(n, ss, cap):>12,.0f}" for n, _ in have))
    print(f"   {'ALL':<22}" + "".join(f"{v:>12,.0f}" for v in allv))

    # where is each stratum's own optimum?
    print()
    print("   ⭐⭐ EACH STRATUM'S OWN BEST SCALE — a PLATEAU needs one scale to be near-best for all four")
    for label, ss, cap in STRATA:
        vals = [(total(n, ss, cap), k) for n, k in have]
        best = min(vals)
        print(f"      {label:<22} best at scale {best[1]:<8g} ({best[0]:>12,.0f})")
    print()
    print("      and the joint optimum, summed over the four strata:")
    for (n, k), v in zip(have, allv):
        mark = " <-- best" if v == min(allv) else ""
        print(f"        scale {k:<8g} {v:>12,.0f}{mark}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
