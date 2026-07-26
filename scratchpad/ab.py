"""Compact multi-arm A/B reader for /tmp/pass0_oracle_bench.tsv.

    python scratchpad/ab.py NEW_ARM [BASE_ARM]        # default base g5

Mass-weighted aggregate (the same statistic pass0_oracle_bench.py --report prints) plus the strata the
owner asks to see called out explicitly, and the per-condition table.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict

import numpy as np

D: dict = defaultdict(dict)
for r in csv.DictReader(open("/tmp/pass0_oracle_bench.tsv"), delimiter="\t"):
    D[r["arm"]][r["cond"]] = r

new = sys.argv[1] if len(sys.argv) > 1 else "m11a"
base = sys.argv[2] if len(sys.argv) > 2 else "g5"
extra = sys.argv[3:]

STRATA = [
    ("ALL 32", lambda c: True),
    ("stranded ss_0.99", lambda c: "ss_0.99" in c),
    ("unstranded ss_0.50", lambda c: "ss_0.50" in c),
    ("unstranded x capON", lambda c: "ss_0.50" in c and "capture_on" in c),
    ("capture OFF", lambda c: "capture_off" in c),
    ("capture ON", lambda c: "capture_on" in c),
    ("verystrong", lambda c: "verystrong" in c),
    ("gdna_none", lambda c: c.startswith("gdna_none")),
    ("low gDNA (1/5/none)", lambda c: any(c.startswith(f"gdna_{t}") for t in ("gdna1", "gdna5", "none"))),
    ("nrna_none", lambda c: "nrna_none" in c),
]

for suf in ("_r0", "_r1"):
    b, p = base + suf, new + suf
    if b not in D or p not in D:
        continue
    conds = sorted(set(D[b]) & set(D[p]))
    print(f"\n{'=' * 104}\n  refit{suf[-1]}   base={b}   new={p}\n{'=' * 104}")
    hdr = f"  {'stratum':<24}{'n':>4}{'base':>9}{'new':>9}{'delta':>9}   better/worse/flat"
    print(hdr)
    for name, sel in STRATA:
        cs = [c for c in conds if sel(c)]
        if not cs:
            continue
        eb = np.array([float(D[b][c]["mwae_all"]) for c in cs])
        ep = np.array([float(D[p][c]["mwae_all"]) for c in cs])
        w = np.array([float(D[b][c]["mass"]) for c in cs])
        nb = int((ep < eb - 0.002).sum())
        nw = int((ep > eb + 0.002).sum())
        print(f"  {name:<24}{len(cs):>4}{np.average(eb, weights=w):>9.4f}"
              f"{np.average(ep, weights=w):>9.4f}{np.average(ep - eb, weights=w):>+9.4f}"
              f"   {nb} / {nw} / {len(cs) - nb - nw}")
    print(f"\n  {'condition':<44}{'base':>9}{'new':>9}{'delta':>9}"
          + "".join(f"{a:>9}" for a in extra))
    for c in conds:
        eb, ep = float(D[b][c]["mwae_all"]), float(D[p][c]["mwae_all"])
        mk = "  BETTER" if ep < eb - 0.002 else ("  worse" if ep > eb + 0.002 else "")
        ex = "".join(
            f"{float(D[a + suf][c]['mwae_all']):>9.4f}" if (a + suf) in D and c in D[a + suf] else f"{'-':>9}"
            for a in extra
        )
        print(f"  {c[5:]:<44}{eb:>9.4f}{ep:>9.4f}{ep - eb:>+9.4f}{ex}{mk}")
