"""RC1 — why does the DL cliff term not kill the over-claim?  b2 = max(0, G^2 - v_msg - v_own).

    OMP_NUM_THREADS=1 python scratchpad/rc1_dl.py [COND ...]
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


def main():
    conds = sys.argv[1:] or CONDS
    tot = [0, 0, 0.0, 0.0]
    for cond in conds:
        ctx = R.Ctx(cond)
        base = ctx.run(trace_node=True)
        nodl = ctx.run(no_dl=True)
        fo, mass = ctx.fo, ctx.mass
        ok = np.isfinite(fo) & (mass > _EPS)
        e0 = np.abs(base["fg_final"] - fo)
        serr = np.abs(ctx.fg_loc - fo)
        T = ok & (ctx.tau_own > 0) & (e0 > serr + 0.02) & (ctx.cls == "exon")
        # replay the DL arithmetic on the RNA leg of every graft message into a hurt exon
        tr = base["trace"][1]
        n_all = n_inert = 0
        gaps = []
        for i in np.flatnonzero(T):
            for t, kp, kn, pp_, pn_ in ((tr[0], "ap", "an", "app", "apn"),
                                        (tr[1], "bp", "bn", "bpp", "bpn")):
                if not (t["r"][i] > 0) or not t["graft"][i]:
                    continue
                rho_msg = base[kp][i] + base[kn][i]
                rho_own = ctx.op[i] + ctx.on[i]
                if rho_msg <= _EPS or rho_own <= _EPS:
                    continue
                G = np.log(rho_msg / rho_own)
                p = base[pp_][i] + base[pn_][i]
                v_msg = 1.0 / p if p > 0 else np.inf
                v_own = ctx.v_own_r[i]
                b2 = max(0.0, G * G - v_msg - v_own)
                n_all += 1
                n_inert += int(b2 <= 0.0)
                gaps.append((G, np.sqrt(2.0 * v_own), np.exp(G), np.exp(np.sqrt(2.0 * v_own))))
        g = np.array(gaps)
        print(f"### {cond}")
        print(f"  graft RNA messages into hurt exons: {n_all}   DL INERT (b2 = 0): {n_inert} "
              f"({n_inert / max(n_all, 1):.1%})")
        print(f"  |G| (log over-claim)  median {np.median(np.abs(g[:, 0])):.3f} "
              f"(x{np.median(np.exp(np.abs(g[:, 0]))):.2f})   "
              f"sqrt(2)*sigma_own median {np.median(g[:, 1]):.3f} (x{np.median(g[:, 3]):.2f})")
        print(f"  mwae(target) base {np.average(e0[T], weights=mass[T]):.4f}  "
              f"DL off {np.average(np.abs(nodl['fg_final'] - fo)[T], weights=mass[T]):.4f}")
        tot[0] += n_all
        tot[1] += n_inert
    print(f"\nPOOLED: {tot[1]}/{tot[0]} = {tot[1] / max(tot[0], 1):.1%} of the graft RNA messages into "
          "hurt full-rank exons pass DL UNDAMPED")


if __name__ == "__main__":
    main()
