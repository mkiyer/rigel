"""RC1 — the fix arms over the FULL suite, split by STRATUM so the cost is visible.

    OMP_NUM_THREADS=1 python scratchpad/rc1_strata.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import rc1_replay as R  # noqa: E402

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
ARMS = {
    "base": dict(),
    "A graft@dst": dict(graft_mode="dst", relay_graft="dst"),
    "B post-pin": dict(graft_mode="post_pin", relay_graft="dst"),
    "E B+meas-local": dict(graft_mode="post_pin", relay_graft="dst", meas_local=True),
}
STRATA = ["ALL", "full-rank (tau>0)", "  RC1 target", "  full-rank exon", "  full-rank other",
          "tau_own=0", "exon", "  exon tau>0", "  exon tau=0", "intron", "boundary",
          "  boundary tau>0", "  boundary tau=0", "intergenic"]


def main():
    conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
    acc = {(a, s): [0.0, 0.0] for a in ARMS for s in STRATA}
    per_cond = {}
    for cond in conds:
        ctx = R.Ctx(cond)
        fo, mass = ctx.fo, ctx.mass
        ok = np.isfinite(fo) & (mass > _EPS)
        res = {k: ctx.run(**v)["fg_final"] for k, v in ARMS.items()}
        e0 = np.abs(res["base"] - fo)
        serr = np.abs(ctx.fg_loc - fo)
        full = ctx.tau_own > 0
        T = ok & full & (e0 > serr + 0.02) & (ctx.cls == "exon")
        masks = {
            "ALL": ok, "full-rank (tau>0)": ok & full, "  RC1 target": T,
            "  full-rank exon": ok & full & (ctx.cls == "exon"),
            "  full-rank other": ok & full & (ctx.cls != "exon"),
            "tau_own=0": ok & ~full, "exon": ok & (ctx.cls == "exon"),
            "  exon tau>0": ok & full & (ctx.cls == "exon"),
            "  exon tau=0": ok & ~full & (ctx.cls == "exon"),
            "intron": ok & (ctx.cls == "intron"), "boundary": ok & (ctx.cls == "boundary"),
            "  boundary tau>0": ok & full & (ctx.cls == "boundary"),
            "  boundary tau=0": ok & ~full & (ctx.cls == "boundary"),
            "intergenic": ok & (ctx.cls == "intergenic"),
        }
        per_cond[cond] = {}
        for a in ARMS:
            e = np.abs(res[a] - fo)
            per_cond[cond][a] = float(np.average(e[ok], weights=mass[ok]))
            for s, m in masks.items():
                if m.any():
                    acc[(a, s)][0] += float((e[m] * mass[m]).sum())
                    acc[(a, s)][1] += float(mass[m].sum())
        print(f"  {cond}", flush=True)

    print(f"\n{'stratum':<22}" + "".join(f"{a:>19}" for a in ARMS))
    for s in STRATA:
        row = f"{s:<22}"
        b = acc[("base", s)][0] / max(acc[("base", s)][1], 1)
        for a in ARMS:
            v = acc[(a, s)][0] / max(acc[(a, s)][1], 1)
            row += f"{v:>11.4f}{'' if a == 'base' else f' ({v - b:+.4f})':>8}"
        print(row)

    print(f"\nper-condition wins/losses vs base (delta mwae, - = better)\n{'condition':<48}"
          + "".join(f"{a:>16}" for a in ARMS if a != "base"))
    tally = {a: [0, 0, 0] for a in ARMS if a != "base"}
    for cond in conds:
        line = f"{cond:<48}"
        for a in ARMS:
            if a == "base":
                continue
            d = per_cond[cond][a] - per_cond[cond]["base"]
            line += f"{d:>+16.4f}"
            tally[a][0 if d < -1e-4 else (1 if abs(d) <= 1e-4 else 2)] += 1
        print(line)
    print(f"\n{'better / flat / worse':<48}"
          + "".join(f"{f'{t[0]}/{t[1]}/{t[2]}':>16}" for a, t in tally.items()))


if __name__ == "__main__":
    main()
