"""ADJ — split RC3's `mismatch_deflate(... & known)` finding BY DIRECTION.

`contradicted = live_msg XOR live_own` conflates two opposite situations:
  (i)  msg = 0, own > 0  — the message ASSERTS AN ABSENCE (the structural-absence bug class);
  (ii) msg > 0, own = 0  — the DESTINATION has no density of its own (unsolved default) and the message is
       the only information there is.
RC3's proposed fix (drop `& known`) kills BOTH at v_own = inf. Count them, weighted by destination error
mass, so the asymmetric form can be adjudicated.
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import rc1_replay as R  # noqa: E402

_EPS = 1e-9
CONDS = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
]


def main():
    tot = {}
    for cond in CONDS:
        ctx = R.Ctx(cond)
        fwd, bwd = ctx.build_relays()  # shipped relay
        _, rho_lf, rho_rf = ctx.rho_faces(ctx.fg_init)
        err = np.abs(ctx.fg_init * 0.0)  # placeholder; use final solved error below
        out = ctx.run()
        err = np.abs(out["fg_final"] - ctx.fo)
        err = np.where(np.isfinite(err), err, 0.0)
        em = err * ctx.mass
        rows = []
        for arrs, src, valid, df, sf, dfv, sfv in (
            (fwd, ctx.sl, ctx.vl, 0, 1, rho_lf, rho_rf),
            (bwd, ctx.sr, ctx.vr, 1, 0, rho_rf, rho_lf),
        ):
            tg, tp, tn, tpg, tpp, tpn, tmg, tmp, tmn, ttau = ctx.transport(
                src, valid, df, sf, arrs, dfv, sfv, no_dl=True
            )
            for nm, t, o, pr in (("g", tg, ctx.og, tpg), ("p", tp, ctx.op, tpp),
                                 ("n", tn, ctx.on, tpn)):
                lt, lo = t > 1e-12, o > 1e-12
                live = valid & (pr > 0.0)
                absent = live & lt.__invert__() & lo  # (i) message asserts ABSENCE
                orphan = live & lt & lo.__invert__()  # (ii) destination has no own density
                unknown = ~np.isfinite(ctx.v_own_g if nm == "g" else ctx.v_own_r)
                rows.append((nm, int((absent & unknown).sum()), float(em[absent & unknown].sum()),
                             int((orphan & unknown).sum()), float(em[orphan & unknown].sum()),
                             int((absent & ~unknown).sum()), int((orphan & ~unknown).sum())))
        agg = {}
        for nm, a_n, a_m, o_n, o_m, ak, ok_ in rows:
            d = agg.setdefault(nm, [0, 0.0, 0, 0.0, 0, 0])
            d[0] += a_n; d[1] += a_m; d[2] += o_n; d[3] += o_m; d[4] += ak; d[5] += ok_
        print(f"\n{cond}")
        print(f"  {'chan':<5}{'(i) msg=0,own>0 @v_own=inf':>32}{'(ii) msg>0,own=0 @v_own=inf':>32}"
              f"{'  (i)@finite':>13}{'(ii)@finite':>12}")
        for nm, d in agg.items():
            print(f"  {nm:<5}{d[0]:>12d} msgs  em={d[1]:>10,.0f}{d[2]:>14d} msgs  em={d[3]:>10,.0f}"
                  f"{d[4]:>13d}{d[5]:>12d}")
            t = tot.setdefault(nm, [0, 0.0, 0, 0.0, 0, 0])
            for k in range(6):
                t[k] += d[k]
    print("\n=== TOTAL over conditions ===")
    for nm, d in tot.items():
        print(f"  {nm:<5}(i) absence-claim @v_own=inf: {d[0]:>7d} msgs, dst err-mass {d[1]:>12,.0f}"
              f"   | (ii) orphan-dst @v_own=inf: {d[2]:>7d} msgs, dst err-mass {d[3]:>12,.0f}"
              f"   | @finite v_own (already killed): (i) {d[4]:>6d}  (ii) {d[5]:>6d}")


if __name__ == "__main__":
    main()
