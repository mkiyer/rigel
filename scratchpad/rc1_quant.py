"""RC1 — quantify the amplification A = r*k and the CHANNEL-differential that causes it.

At every graft edge (boundary src -> exon dst) into a hurt full-rank exon, compare
  * the gDNA channel:   rho_g(src, unspliced crossing)  vs  rho_g(dst, contained)   [the capture cliff]
  * the mature channel: rho_mu(src, spliced junction)   vs  rho_R(dst, contained)   [no cliff]
If the first ratio is far from 1 and the second near 1, the source's own component RATIO — which is what the
message transports — is invalid, and the delivered claim is inflated by exactly that difference.

    OMP_NUM_THREADS=1 python scratchpad/rc1_quant.py [COND ...]
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


def q(x, w=None):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not x.size:
        return "n/a"
    return f"p25 {np.percentile(x, 25):>9.3g}  med {np.median(x):>9.3g}  p75 {np.percentile(x, 75):>9.3g}"


def main():
    conds = sys.argv[1:] or CONDS
    for cond in conds:
        ctx = R.Ctx(cond)
        base = ctx.run(trace_node=True)
        fo, mass, M, E_g, E_r = ctx.fo, ctx.mass, ctx.M, ctx.E_g, ctx.E_r
        err = np.abs(base["fg_final"] - fo)
        serr = np.abs(ctx.fg_loc - fo)
        ok = np.isfinite(fo) & (mass > _EPS)
        tgt = ok & (ctx.tau_own > 0) & (err > serr + 0.02) & (ctx.cls == "exon")
        tr = base["trace"][1]
        rows = []
        for i in np.flatnonzero(tgt):
            for t in (tr[0], tr[1]):
                if not (t["r"][i] > 0) or not t["graft"][i]:
                    continue
                s = int(t["src"][i])
                mu = t["gp"][i] + t["gn"][i]
                if mu <= _EPS:
                    continue
                # oracle contained densities at the destination exon
                rho_g_dst_oracle = ctx.G[i] / E_g[i]
                rho_R_dst_oracle = (ctx.Rp[i] + ctx.Rn[i]) / E_r[i]
                # the source boundary's OWN (message-free) unspliced gDNA density + its oracle gDNA density
                rho_g_src_own = ctx.og[s]
                rho_g_src_oracle = ctx.G[s] / max(ctx.E_g[s], _EPS)
                rows.append((
                    mass[i],
                    t["r"][i] * t["pin"][i],                          # amplification A
                    mu * E_r[i] / M[i],                               # count-only claim
                    mu * E_r[i] / M[i] * t["r"][i] * t["pin"][i],      # graft's delivered claim
                    1.0 - fo[i],                                      # oracle f_R
                    rho_g_src_oracle / max(rho_g_dst_oracle, _EPS),    # gDNA channel ratio (ORACLE)
                    mu / max(rho_R_dst_oracle, _EPS),                 # mature channel ratio (ORACLE)
                    rho_g_src_own / max(ctx.og[i], _EPS),              # gDNA channel ratio (as SOLVED)
                    t["r"][i], t["pin"][i],
                    ctx.ESP[t["sf"]][s], E_r[i],
                ))
        if not rows:
            continue
        A = np.array(rows, float)
        w = A[:, 0]
        print("=" * 108)
        print(f"### {cond}    graft edges into hurt full-rank exons: n={len(rows)}")
        print(f"  amplification A = r*k                        {q(A[:, 1])}")
        print(f"  count-only claim  rho_mu*E_r/M               {q(A[:, 2])}")
        print(f"  graft delivered claim  (= count * A)         {q(A[:, 3])}")
        print(f"  ORACLE f_R at the exon                       {q(A[:, 4])}")
        print(f"  count-only / oracle f_R  (eff-len fidelity)  {q(A[:, 2] / np.maximum(A[:, 4], _EPS))}")
        print(f"  delivered / oracle f_R   (the OVER-CLAIM)    {q(A[:, 3] / np.maximum(A[:, 4], _EPS))}")
        print("  --- the CHANNEL differential (oracle truth, no solver in the loop) ---")
        print(f"  gDNA   rho_g(seam)/rho_g(exon)               {q(A[:, 5])}")
        print(f"  mature rho_mu(seam)/rho_R(exon)              {q(A[:, 6])}")
        print(f"  ratio of the two (= the false ratio factor)  {q(A[:, 6] / np.maximum(A[:, 5], _EPS))}")
        print(f"  gDNA seam/exon as SOLVED (own densities)     {q(A[:, 7])}")
        print(f"  reframe r                                    {q(A[:, 8])}")
        print(f"  pin k                                        {q(A[:, 9])}")
        print(f"  E_r(dst)/E_spl(src)                          {q(A[:, 11] / A[:, 10])}")
        print(f"  mass-weighted mean over-claim factor         "
              f"{np.average(A[:, 3] / np.maximum(A[:, 4], _EPS), weights=w):.1f}x   "
              f"(count-only would be {np.average(A[:, 2] / np.maximum(A[:, 4], _EPS), weights=w):.2f}x)")


if __name__ == "__main__":
    main()
