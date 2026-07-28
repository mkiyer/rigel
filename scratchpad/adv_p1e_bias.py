"""ADVERSARIAL re-run of P1e §4.5's BIAS-vs-DISPERSION decomposition, with the strata the original omitted.

`p1e_9_bias.py` drops `intergenic` destinations (line 34) and reports only all-live / lambda-emitting /
one-component. The report's sentence -- "the bias share is 0.0-8.3 % in every condition AND STRATUM" -- is
therefore untested on (a) the intergenic class, whose delta is a near-CONSTANT -0.57/-0.70, and (b) the
GRAFT x one-component stratum, which is where the term's effect lives (section 2b).

    OMP_NUM_THREADS=1 python scratchpad/adv_p1e_bias.py
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import p1e_lib as L  # noqa: E402

ORDER = [
    "gdna300_ss0.99_present_capOFF",
    "gdna300_ss0.50_present_capOFF",
    "gdna300_ss0.50_none_capOFF",
    "gdna100_ss0.50_none_capOFF",
    "gdna300_ss0.99_present_capON",
    "gdna100_ss0.50_present_VERYSTRONG",
]
print(f"{'condition':<32}{'stratum':<26}{'n':>6}{'E[d]':>9}{'sd(d)':>9}"
      f"{'E[d^2]':>9}{'BIAS share':>12}{'mean b^2':>10}{'b2*M share':>12}")
print("-" * 125)
for name in ORDER:
    inp, dbg = L.solve(L.CONDS[name])
    t, _ = L.message_table(inp, dbg)
    live = t["nsup"] > 0
    tot = float((t["bhat2"] * t["M"])[live].sum())
    strata = [
        ("all live (incl intergenic)", live),
        ("intergenic dst ONLY", live & (t["cls"] == "intergenic")),
        ("one-component", live & ~t["lam_emit"]),
        ("one-comp x non-intergenic", live & ~t["lam_emit"] & (t["cls"] != "intergenic")),
        ("GRAFT x one-component", live & ~t["lam_emit"] & t["graft"]),
        ("GRAFT x lambda-emitting", live & t["lam_emit"] & t["graft"]),
    ]
    for tag, m in strata:
        if m.sum() < 20:
            continue
        d = t["delta"][m]
        e2 = float(np.mean(d * d))
        b = float(np.mean(d)) ** 2
        share = float((t["bhat2"] * t["M"])[m].sum()) / max(tot, 1e-9)
        print(f"{name:<32}{tag:<26}{int(m.sum()):>6}{np.mean(d):>+9.4f}{np.std(d):>9.4f}"
              f"{e2:>9.4f}{100 * b / max(e2, 1e-12):>11.1f}%{np.mean(t['bhat2'][m]):>10.4f}"
              f"{100 * share:>11.1f}%")
    print()
