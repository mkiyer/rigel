"""RC1 — END-TO-END arithmetic of the RNA measurement message for one node.

    OMP_NUM_THREADS=1 python scratchpad/rc1_trace.py COND NODE [NODE ...]
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import rc1_replay as R  # noqa: E402

_EPS = 1e-9


def trace(ctx, i, base):
    fo, cls, M, E_g, E_r = ctx.fo, ctx.cls, ctx.M, ctx.E_g, ctx.E_r
    tr = base["trace"][1]  # iteration 2 (the one that produced the shipped answer)
    print("=" * 118)
    print(f"NODE {i}  [{cls[i]}]  cond={ctx.cond}")
    print(f"  M(unspliced mass)={M[i]:,.1f}   E_g={E_g[i]:,.1f}  E_r={E_r[i]:,.1f}   "
          f"free=(+{int(ctx.fp_a[i])},-{int(ctx.fn_a[i])})")
    print(f"  ORACLE   f_g={fo[i]:.4f}   gDNA={ctx.G[i]:,.0f}  RNA+={ctx.Rp[i]:,.1f}  RNA-={ctx.Rn[i]:,.1f}"
          f"   -> oracle f_R(total)={1 - fo[i]:.4f}")
    print(f"  SELF-SOLVE f_g={ctx.fg_loc[i]:.4f}  (own rho: g={ctx.og[i]:.4g} +={ctx.op[i]:.4g} "
          f"-={ctx.on[i]:.4g}; tau_own={ctx.tau_own[i]:.4g})")
    print(f"  SOLVED   f_g={base['fg_final'][i]:.4f}")
    print(f"  delivered to psi:  mo_g={np.exp(base['mo_g'][i]):.4f}/p={base['cm_g'][i]:.4g}   "
          f"mo_p={np.exp(base['mo_p'][i]):.4f}/p={base['cm_p'][i]:.4g}   "
          f"mo_n={np.exp(base['mo_n'][i]):.4f}/p={base['cm_n'][i]:.4g}   "
          f"lam f_g_eq={1 / (1 + np.exp(-base['lam_msg'][i])):.4f}/p={base['c_tau'][i]:.4g}")

    for side, t, pre in (("LEFT (fwd)", tr[0], ("ap", "app", "amp", "an", "apn", "amn", "ag", "apg")),
                         ("RIGHT (bwd)", tr[1], ("bp", "bpp", "bmp", "bn", "bpn", "bmn", "bg", "bpg"))):
        src = int(t["src"][i])
        if not (t["r"][i] > 0):
            print(f"\n  --- {side}: no valid neighbour")
            continue
        sf, df = t["sf"], t["df"]
        SPs, SNs, ESPs = ctx.SP[sf][src], ctx.SN[sf][src], ctx.ESP[sf][src]
        rho_mu_p, rho_mu_n = ctx.spl_p_f[sf][src], ctx.spl_n_f[sf][src]
        print(f"\n  --- {side}  src node {src} [{cls[src]}] oracle f_g={fo[src]:.4f} "
              f"mass={M[src]:,.1f} src-face={sf} dst-face={df}")
        print(f"      MEASURED SPLICED at src face {sf}:  SP={SPs:,.2f}  SN={SNs:,.2f}   "
              f"E_spl={ESPs:,.2f}  ->  rho_mu+ = {rho_mu_p:.5g}  rho_mu- = {rho_mu_n:.5g}")
        print(f"      relayed src density:   rho_g={t['rg_s'][i]:.5g}(p={t['pg_s'][i]:.4g})  "
              f"rho_+={t['rp_s'][i]:.5g}(p={t['pp_s'][i]:.4g})  rho_-={t['rn_s'][i]:.5g}(p={t['pn_s'][i]:.4g})")
        print(f"      REFRAME r = rho_tot(dst face {df})/rho_tot(src face {sf}) = "
              f"{t['dst_face'][i]:.5g} / {t['src_face'][i]:.5g} = {t['r'][i]:.5g}"
              f"   graft={bool(t['graft'][i])} peel={bool(t['peel'][i])} s2t={t['s2t'][i]:.4g}")
        for lab, rel, mu in (("+", t["rp_s"][i], t["gp"][i]), ("-", t["rn_s"][i], t["gn"][i])):
            tot = (rel + mu) * t["r"][i]
            if tot <= 0:
                continue
            print(f"      GRAFT {lab}: (rho_relayed {rel:.5g} + rho_mu {mu:.5g}) x r {t['r'][i]:.5g} "
                  f"= {tot:.5g}   [count share {mu / max(rel + mu, _EPS):.1%}, "
                  f"r contributes x{t['r'][i]:.4g}]")
        print(f"      PIN k = M/(sum rho.E) = {t['pin'][i]:.5g}   pre-pin (g,+,-) = "
              f"({t['pre_pin'][0][i]:.5g}, {t['pre_pin'][1][i]:.5g}, {t['pre_pin'][2][i]:.5g})")
        g_, p_, n_ = (base[pre[6]][i], base[pre[0]][i], base[pre[3]][i])
        print(f"      DELIVERED densities  g={g_:.5g} +={p_:.5g} -={n_:.5g}  ->  as FRACTIONS "
              f"f_g={g_ * E_g[i] / M[i]:.4f} f_+={p_ * E_r[i] / M[i]:.4f} f_-={n_ * E_r[i] / M[i]:.4f}")
        print(f"      precisions (post-DL): mode-fuse pg={base[pre[7]][i]:.4g} pp={base[pre[1]][i]:.4g} "
              f"pn={base[pre[4]][i]:.4g} | measurement mp={base[pre[2]][i]:.4g} mn={base[pre[5]][i]:.4g}")
        print(f"      spliced-count precision seed _spc={t['spc'][i]:.4g}  _snc={t['snc'][i]:.4g} "
              f"(= SP/(1+SP*s2t), SP={t['sp'][i]:.4g} SN={t['sn'][i]:.4g})")
        # the claim decomposition, in FRACTION units
        for lab, mu, rel in (("+", t["gp"][i], t["rp_s"][i]), ("-", t["gn"][i], t["rn_s"][i])):
            if mu <= 0 and rel <= 0:
                continue
            f_raw = mu * E_r[i] / M[i]
            f_r = mu * t["r"][i] * E_r[i] / M[i]
            f_pin = mu * t["r"][i] * t["pin"][i] * E_r[i] / M[i]
            print(f"      CLAIM{lab} from the COUNT alone:  rho_mu*E_r/M = {f_raw:.4f}"
                  f"  -> xr = {f_r:.4f}  -> xpin = {f_pin:.4f}"
                  f"   (oracle f_R = {1 - fo[i]:.4f})")


def main():
    cond = sys.argv[1]
    nodes = [int(x) for x in sys.argv[2:]]
    ctx = R.Ctx(cond)
    d, dd, base = R.fidelity(ctx)
    print(f"# replay fidelity max|dfg| = {d:.3g}  (0 == exact)")
    base = ctx.run(trace_node=True)
    for i in nodes:
        trace(ctx, i, base)


if __name__ == "__main__":
    main()
