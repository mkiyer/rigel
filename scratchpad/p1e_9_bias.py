"""P1e (4b) — WHY capture-ON regresses: is the conservation violation a BIAS or a VARIANCE?

b^2_cons = max(0, delta^2 - declared) treats the whole violation as random dispersion. If E[delta] != 0 the
violation is a systematic MIS-SPECIFICATION, and pricing a bias as a variance only widens messages — it can
never move the mode toward truth. Decompose E[delta^2] = (E delta)^2 + Var(delta) per regime.

    OMP_NUM_THREADS=1 python scratchpad/p1e_9_bias.py
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
    "none_ss0.50_present_capOFF",
]
print(f"{'condition':<36}{'stratum':<18}{'n':>6}{'E[delta]':>10}{'sd(delta)':>11}"
      f"{'E[d^2]':>10}{'(Ed)^2':>9}{'BIAS share':>12}{'mean b^2':>10}")
print("-" * 122)
for name in ORDER:
    inp, dbg = L.solve(L.CONDS[name])
    t, _ = L.message_table(inp, dbg)
    live = (t["nsup"] > 0) & (t["cls"] != "intergenic")
    for tag, m in (("all live", live),
                   ("lambda-emitting", live & t["lam_emit"]),
                   ("one-component", live & ~t["lam_emit"])):
        if m.sum() < 20:
            continue
        d = t["delta"][m]
        e2 = float(np.mean(d * d))
        b = float(np.mean(d)) ** 2
        print(f"{name:<36}{tag:<18}{int(m.sum()):>6}{np.mean(d):>+10.4f}{np.std(d):>11.4f}"
              f"{e2:>10.4f}{b:>9.4f}{100 * b / max(e2, 1e-12):>11.1f}%"
              f"{np.mean(t['bhat2'][m]):>10.4f}")
