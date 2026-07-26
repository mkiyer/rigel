"""RC3 — offline A/B of each candidate structural fix, on the bit-faithful offline solver (`rc3_solve`).

    OMP_NUM_THREADS=1 python scratchpad/rc3_ab.py [--conds a,b] [--arms peel_zero_gate,...]
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from rc3_replay import SUITE, build, load, oracle  # noqa: E402
from rc3_solve import SWITCHES, solve  # noqa: E402

_EPS = 1e-9

ap = argparse.ArgumentParser()
ap.add_argument("--conds", default=None)
ap.add_argument("--arms", default=None)
a = ap.parse_args()

conds = (
    a.conds.split(",")
    if a.conds
    else sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
)
arms = ["BASE"] + (a.arms.split(",") if a.arms else list(SWITCHES))

rows: dict[str, list] = {k: [] for k in ("mass", "oracle", "self", "cls", "full", "cond")}
res: dict[str, list] = {k: [] for k in arms}

for ci, cond in enumerate(conds):
    ctx = load(cond)
    S = build(ctx)
    fo = oracle(ctx)
    ok = np.isfinite(fo) & (S["mass"] > _EPS)
    idx = np.flatnonzero(ok)
    rows["mass"].append(S["mass"][idx])
    rows["oracle"].append(fo[idx])
    rows["self"].append(S["fg_loc"][idx])
    rows["cls"].append(S["cls"][idx])
    rows["full"].append((S["tau_own"][idx] > _EPS))
    rows["cond"].append(np.full(idx.size, cond))
    for arm in arms:
        kw = {} if arm == "BASE" else dict.fromkeys(arm.split("+"), True)
        res[arm].append(solve(ctx, S, **kw)[idx])
    print(f"  [{ci + 1:>2}/{len(conds)}] {cond}", flush=True)

D = {k: np.concatenate(v) for k, v in rows.items()}
Rr = {k: np.concatenate(v) for k, v in res.items()}
mass, fo, self_fg = D["mass"], D["oracle"], D["self"]
serr = np.abs(self_fg - fo)
base_err = np.abs(Rr["BASE"] - fo)
full = D["full"]
hurt = full & (base_err > serr + 0.02)

strata = [
    ("ALL", np.ones_like(full, bool)),
    ("full-rank (τ_own>0)", full),
    ("HURT full-rank", hurt),
    ("τ_own=0", ~full),
    ("exon", D["cls"] == 2),
    ("intron", D["cls"] == 1),
    ("boundary", D["cls"] == -1),
]
print(f"\n{'arm':<20}" + "".join(f"{lab:>22}" for lab, _ in strata))
for arm in arms:
    e = np.abs(Rr[arm] - fo)
    cells = []
    for lab, m in strata:
        v = np.average(e[m], weights=mass[m])
        b = np.average(base_err[m], weights=mass[m])
        cells.append(f"{v:.4f}({v - b:+.4f})".rjust(22))
    print(f"{arm:<20}" + "".join(cells))

print("\nper-condition (ALL, mass-weighted):")
cl = sorted(set(D["cond"].tolist()))
hdr = f"{'condition':<48}" + "".join(f"{arm[:14]:>16}" for arm in arms)
print(hdr)
for c in cl:
    m = D["cond"] == c
    b = np.average(base_err[m], weights=mass[m])
    cells = []
    for arm in arms:
        v = np.average(np.abs(Rr[arm] - fo)[m], weights=mass[m])
        cells.append((f"{v:.4f}" if arm == "BASE" else f"{v:.4f}({v - b:+.4f})").rjust(16))
    print(f"{c:<48}" + "".join(cells))
