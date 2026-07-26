"""RC1 — the fix arms over the FULL 32-scenario suite (pass-0 vs oracle, refit=0), offline replay.

    OMP_NUM_THREADS=1 python scratchpad/rc1_suite.py
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
    "BOTH@dst": dict(graft_mode="dst", relay_graft="dst"),
    "BOTHpostpin": dict(graft_mode="post_pin", relay_graft="dst"),
}


def main():
    conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
    tot = {k: [0.0, 0.0] for k in ARMS}
    tgt = {k: [0.0, 0.0] for k in ARMS}
    print(f"{'condition':<48}" + "".join(f"{k:>13}" for k in ARMS) + "   n_tgt")
    for cond in conds:
        ctx = R.Ctx(cond)
        fo, mass = ctx.fo, ctx.mass
        ok = np.isfinite(fo) & (mass > _EPS)
        res = {k: ctx.run(**v)["fg_final"] for k, v in ARMS.items()}
        err0 = np.abs(res["base"] - fo)
        serr = np.abs(ctx.fg_loc - fo)
        T = ok & (ctx.tau_own > 0) & (err0 > serr + 0.02) & (ctx.cls == "exon")
        line = f"{cond:<48}"
        for k in ARMS:
            e = np.abs(res[k] - fo)
            tot[k][0] += float((e[ok] * mass[ok]).sum())
            tot[k][1] += float(mass[ok].sum())
            tgt[k][0] += float((e[T] * mass[T]).sum())
            tgt[k][1] += float(mass[T].sum())
            line += f"{float(np.average(e[ok], weights=mass[ok])):>13.4f}"
        print(line + f"   {int(T.sum()):>5}", flush=True)
    print("\n" + "-" * 92)
    print(f"{'SUITE mwae (all live nodes)':<48}"
          + "".join(f"{tot[k][0] / tot[k][1]:>13.4f}" for k in ARMS))
    print(f"{'  on the RC1 target set':<48}"
          + "".join(f"{tgt[k][0] / max(tgt[k][1], 1):>13.4f}" for k in ARMS))
    print(f"{'  target err-MASS':<48}"
          + "".join(f"{tgt[k][0]:>13,.0f}" for k in ARMS))


if __name__ == "__main__":
    main()
