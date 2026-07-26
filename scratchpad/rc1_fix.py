"""RC1 — counterfactual ARMS on the graft, scored per condition (pass-0 vs oracle, refit=0).

    OMP_NUM_THREADS=1 python scratchpad/rc1_fix.py [--all] [COND ...]
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import rc1_replay as R  # noqa: E402

_EPS = 1e-9
DEFAULT = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna100_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
]
ARMS = {
    "base": dict(),
    "cmb@dst": dict(graft_mode="dst"),
    "cmb postpin": dict(graft_mode="post_pin"),
    "BOTH@dst": dict(graft_mode="dst", relay_graft="dst"),
    "BOTH postpin": dict(graft_mode="post_pin", relay_graft="dst"),
    "BOTH off": dict(graft_mode="off", relay_graft="off"),
}


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("-")]
    conds = args or DEFAULT
    tot = {k: [0.0, 0.0] for k in ARMS}  # [sum |err|*mass, sum mass]
    tgt_tot = {k: [0.0, 0.0] for k in ARMS}
    print(f"{'condition':<46}" + "".join(f"{k:>16}" for k in ARMS))
    for cond in conds:
        ctx = R.Ctx(cond)
        fo, mass = ctx.fo, ctx.mass
        ok = np.isfinite(fo) & (mass > _EPS)
        res = {k: ctx.run(**v) for k, v in ARMS.items()}
        b = res["base"]["fg_final"]
        err0 = np.abs(b - fo)
        serr = np.abs(ctx.fg_loc - fo)
        tgt = ok & (ctx.tau_own > 0) & (err0 > serr + 0.02) & (ctx.cls == "exon")
        line = f"{cond:<46}"
        for k in ARMS:
            e = np.abs(res[k]["fg_final"] - fo)
            a = float(np.average(e[ok], weights=mass[ok]))
            t = float(np.average(e[tgt], weights=mass[tgt])) if tgt.any() else float("nan")
            tot[k][0] += float((e[ok] * mass[ok]).sum())
            tot[k][1] += float(mass[ok].sum())
            tgt_tot[k][0] += float((e[tgt] * mass[tgt]).sum())
            tgt_tot[k][1] += float(mass[tgt].sum())
            line += f"{a:>8.4f}/{t:<7.4f}"
        print(line + f"   [target n={int(tgt.sum())}]")
    print("\n" + "-" * 120)
    print(f"{'POOLED mwae (all / target)':<46}"
          + "".join(f"{tot[k][0] / tot[k][1]:>8.4f}/{tgt_tot[k][0] / max(tgt_tot[k][1], 1):<7.4f}"
                    for k in ARMS))


if __name__ == "__main__":
    main()
