"""RC1 — is the gDNA-measurement bias (37.2% of the target error mass) the SAME defect?

`_pin_v` forces the message's components to sum to the destination's own mass, so on a two-live-component
message f_g^msg + f_R^msg == 1 EXACTLY: an over-claimed RNA leg IS the under-claimed gDNA leg.  Check the
identity, then re-bin (mo_g - oracle f_g) under base vs the fix.

    OMP_NUM_THREADS=1 python scratchpad/rc1_gchan.py [COND ...]
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import rc1_replay as R  # noqa: E402

_EPS = 1e-9
CONDS = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
]
BINS = [0.0, 0.5, 0.9, 0.99, 1.01]


def main():
    conds = sys.argv[1:] or CONDS
    acc = {}
    for cond in conds:
        ctx = R.Ctx(cond)
        base = ctx.run(trace_node=True)
        fix = ctx.run(graft_mode="post_pin", relay_graft="dst")
        fo, mass = ctx.fo, ctx.mass
        ok = np.isfinite(fo) & (mass > _EPS)
        e0 = np.abs(base["fg_final"] - fo)
        serr = np.abs(ctx.fg_loc - fo)
        T = ok & (ctx.tau_own > 0) & (e0 > serr + 0.02) & (ctx.cls == "exon")

        # (1) the pin identity, per delivered message, on graft edges into hurt exons
        tr = base["trace"][1]
        bad = []
        for i in np.flatnonzero(T):
            for t, kg, kp, kn, pg_, pp_, pn_ in (
                    (tr[0], "ag", "ap", "an", "apg", "app", "apn"),
                    (tr[1], "bg", "bp", "bn", "bpg", "bpp", "bpn")):
                if not (t["r"][i] > 0) or not t["graft"][i]:
                    continue
                if not (base[pg_][i] > 0 and (base[pp_][i] + base[pn_][i]) > 0):
                    continue  # a one-legged message: the pin substitutes the node's own density
                s = (base[kg][i] * ctx.E_g[i] + (base[kp][i] + base[kn][i]) * ctx.E_r[i]) / ctx.M[i]
                bad.append(abs(s - 1.0))
        print(f"### {cond}")
        print(f"  pin identity  max|f_g^msg + f_R^msg - 1| over {len(bad)} two-legged graft "
              f"messages = {max(bad) if bad else float('nan'):.3g}")

        # (2) the delivered gDNA measurement mode vs oracle, binned
        for lab, arm in (("base", base), ("fix(post-pin)", fix)):
            mg = np.exp(arm["mo_g"])
            for lo, hi in zip(BINS[:-1], BINS[1:]):
                m = T & (fo >= lo) & (fo < hi)
                if not m.any():
                    continue
                k = (lab, f"[{lo:.2f},{hi:.2f})")
                acc.setdefault(k, [0.0, 0.0, 0])
                acc[k][0] += float((mg[m] - fo[m]).sum())
                acc[k][1] += float((np.exp(arm["mo_p"])[m] + np.exp(arm["mo_n"])[m]
                                    - (1.0 - fo[m])).sum())
                acc[k][2] += int(m.sum())
    print(f"\n{'oracle f_g bin':<18}{'arm':<16}{'n':>6}{'mean(mo_g - oracle f_g)':>26}"
          f"{'mean(mo_R - oracle f_R)':>26}")
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        b = f"[{lo:.2f},{hi:.2f})"
        for lab in ("base", "fix(post-pin)"):
            v = acc.get((lab, b))
            if not v:
                continue
            print(f"{b:<18}{lab:<16}{v[2]:>6}{v[0] / v[2]:>26.4f}{v[1] / v[2]:>26.4f}")


if __name__ == "__main__":
    main()
