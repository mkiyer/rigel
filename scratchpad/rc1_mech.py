"""RC1 — the MECHANISM, quantified: where does the RNA measurement over-claim come from?

For every full-rank exon node in a condition, decompose the delivered RNA claim into
  (i)  the raw measured mature density -> fraction   rho_mu * E_r^dst / M_dst
  (ii) the reframe r
  (iii) the pin k
and test the r-cancellation identity.

    OMP_NUM_THREADS=1 python scratchpad/rc1_mech.py [COND ...]
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
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
]


def analyse(cond):
    ctx = R.Ctx(cond)
    base = ctx.run(trace_node=True)
    r1 = ctx.run(force_r1=True)
    ng = ctx.run(no_graft=True)
    fo = ctx.fo
    fg = base["fg_final"]
    self_fg = ctx.fg_loc
    mass, M, E_g, E_r = ctx.mass, ctx.M, ctx.E_g, ctx.E_r
    err = np.abs(fg - fo)
    serr = np.abs(self_fg - fo)
    ok = np.isfinite(fo) & (mass > _EPS)
    tgt = ok & (ctx.tau_own > 0) & (err > serr + 0.02) & (ctx.cls == "exon")
    print("=" * 110)
    print(f"### {cond}")
    print(f"    live nodes {int(ok.sum())}  mwae={np.average(err[ok], weights=mass[ok]):.4f}   "
          f"target(full-rank exon, msgs hurt) n={int(tgt.sum())} "
          f"err-mass={float((err * mass)[tgt].sum()):,.0f}")

    def mw(x, m):
        return float(np.average(x[m], weights=mass[m])) if m.any() else float("nan")

    print(f"\n    ARM                     mwae(all)   mwae(target)   mean f_g(target)  oracle={mw(fo, tgt):.4f}")
    for lab, arm in (("base", base), ("r forced to 1", r1), ("graft removed", ng)):
        f = arm["fg_final"]
        print(f"    {lab:<22}{mw(np.abs(f - fo), ok):>10.4f}{mw(np.abs(f - fo), tgt):>14.4f}"
              f"{mw(f, tgt):>18.4f}")
    print(f"    {'self-solve':<22}{mw(serr, ok):>10.4f}{mw(serr, tgt):>14.4f}{mw(self_fg, tgt):>18.4f}")

    # ── per-side decomposition on the target nodes ──
    tr = base["trace"][1]
    rows = []
    for i in np.flatnonzero(tgt):
        for side, t, k_ap, k_ag in (("L", tr[0], "ap", "ag"), ("R", tr[1], "bp", "bg")):
            if not (t["r"][i] > 0) or not t["graft"][i]:
                continue
            src = int(t["src"][i])
            sf = t["sf"]
            mu = t["gp"][i] + t["gn"][i]
            rel = t["rp_s"][i] + t["rn_s"][i]
            if mu <= 0 and rel <= 0:
                continue
            rows.append(dict(
                node=i, side=side, src=src, mass=mass[i], oracle_fR=1.0 - fo[i],
                rho_mu=mu, rho_rel=rel, rho_g_src=t["rg_s"][i], r=t["r"][i], pin=t["pin"][i],
                f_count=mu * E_r[i] / M[i],
                f_src_ratio=(rel + mu) * E_r[i] / max((rel + mu) * E_r[i] + t["rg_s"][i] * E_g[i], _EPS),
                f_deliv=(base[k_ap][i] + base["an" if side == "L" else "bn"][i]) * E_r[i] / M[i],
                sp=t["sp"][i] + t["sn"][i],
                E_spl=ctx.ESP[sf][src], E_r=E_r[i], E_g=E_g[i],
                rho_tot_src=t["src_face"][i], rho_tot_dst=t["dst_face"][i],
            ))
    if not rows:
        return
    A = {k: np.array([r[k] for r in rows], float) for k in rows[0] if k not in ("side",)}
    print(f"\n    GRAFT edges into target exons: n={len(rows)}")
    print(f"      raw count claim   f = rho_mu*E_r/M      median {np.median(A['f_count']):.4f}"
          f"   (oracle f_R median {np.median(A['oracle_fR']):.4f})")
    print(f"      source-ratio claim (r-INVARIANT)        median {np.median(A['f_src_ratio']):.4f}")
    print(f"      DELIVERED fraction                      median {np.median(A['f_deliv']):.4f}")
    print(f"      reframe r                               median {np.median(A['r']):.4g}"
          f"  p10 {np.percentile(A['r'], 10):.4g}  p90 {np.percentile(A['r'], 90):.4g}")
    print(f"      pin k                                   median {np.median(A['pin']):.4g}")
    print(f"      rho_mu / rho_g(src)                     median {np.median(A['rho_mu'] / np.maximum(A['rho_g_src'], _EPS)):.4g}")
    print(f"      rho_g(src) / rho_g_own(dst)             median "
          f"{np.median(A['rho_g_src'] / np.maximum(np.array([ctx.og[int(r['node'])] for r in rows]), _EPS)):.4g}")
    print(f"      E_r(dst)/E_spl(src)                     median {np.median(A['E_r'] / A['E_spl']):.4g}")
    # the identity check: delivered == source ratio (r cancels)?
    d = np.abs(A["f_deliv"] - A["f_src_ratio"])
    print(f"      |delivered - source-ratio|  max {d.max():.3g}  median {np.median(d):.3g}"
          "   <- 0 proves r CANCELS in the pin")
    return ctx, base


def main():
    conds = sys.argv[1:] or CONDS
    for c in conds:
        analyse(c)


if __name__ == "__main__":
    main()
