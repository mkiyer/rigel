"""(e) numerator-only vs both-sums; (f) the alpha operating window; (g) the scaling law;
(h) the ratio-of-sums (Jensen) term."""
from __future__ import annotations

import numpy as np
import od_tb_core as T
from scipy.special import betainc
from scipy.optimize import brentq
from od_tb_agg import D, estimator_terms, seed_view, size_hist, pooled_bias
from od_tb_bias import ALPHA_BONF, min_n_rejectable
from od_tb_contam import keep_flags


def part_e():
    print("=" * 104)
    print("(e) WHERE THE REJECTED SEED IS DROPPED FROM.  a_null = 3, real size distributions.")
    print("    both  = drop from numerator AND denominator (correct);  num   = drop from numerator only.")
    print(f"    {'sample':10s} {'alpha':>9s} {'od':>6s} {'q_sc':>9s} {'bias both':>10s} "
          f"{'bias num-only':>14s} {'ratio':>7s}")
    for name, d in D.items():
        s, t, w, k = seed_view(d)
        _, sc, nc = estimator_terms(s, t, w, k)
        u, c = size_hist(nc)
        for alpha in (ALPHA_BONF, 0.01, 0.05):
            for od in (0.1, 0.2):
                r = pooled_bias(u, c, 3.0, alpha, od)
                b1 = 1.0 - r["od_hat"] / od
                b2 = 1.0 - r["od_hat_numonly"] / od
                print(f"    {name:10s} {alpha:9.3g} {od:6.3f} {r['q_scale_weighted']:9.2e} "
                      f"{b1:10.3e} {b2:14.3e} {b2/max(b1,1e-300):7.3f}")
    print()


def part_f():
    print("=" * 104)
    print("(f) THE OPERATING WINDOW: truncation bias vs contamination removed, a_null = 3.")
    print("    delta = od_raw - od_screened (what the screen removes, real data).")
    print("    bias  = exact truncation bias in od units, evaluated at od = min(od_screened, ceiling).")
    print(f"    {'sample':10s} {'alpha':>10s} {'#rej':>6s} {'delta od':>9s} {'bias od':>10s} "
          f"{'bias/delta':>11s}")
    for name, d in D.items():
        s, t, w, k = seed_view(d)
        ex, sc, nc = estimator_terms(s, t, w, k)
        u, c = size_hist(nc)
        od_raw = ex.sum() / sc.sum()
        for alpha in (ALPHA_BONF, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 5e-2):
            keep, _ = keep_flags(nc, s, 3.0, alpha)
            od_scr = ex[keep].sum() / sc[keep].sum()
            delta = od_raw - od_scr
            ev = float(np.clip(od_scr, 1e-4, 0.2))
            bias = (1.0 - pooled_bias(u, c, 3.0, alpha, ev)["od_hat"] / ev) * ev
            print(f"    {name:10s} {alpha:10.3g} {int((~keep).sum()):6d} {delta:+9.4f} "
                  f"{bias:10.2e} {bias/max(abs(delta),1e-12):11.2e}")
        print()


def _x_boundary(a_null, alpha):
    """Continuum acceptance boundary in sense-FRACTION: 2*F_{Beta(a,a)}(x) = alpha."""
    return brentq(lambda x: 2 * betainc(a_null, a_null, x) - alpha, 1e-14, 0.5)


def part_g():
    print("=" * 104)
    print("(g) THE SCALING LAW (large-n continuum).  x solves 2*F_{a_null}(x)=alpha;")
    print("    q(od) = 2*F_{a_true}(x)  and  bias ~ q(1-od)/od.  Note q ~ alpha^(a_true/a_null).")
    for a_null in (2.0, 3.0):
        print(f"\n    a_null = {a_null:.0f}")
        print(f"      {'alpha':>10s} {'x (sense frac)':>15s} " +
              "".join(f"{'q@od='+format(od,'.3g'):>14s}" for od in (0.2, 0.143, 0.1, 0.0345)))
        for alpha in (ALPHA_BONF, 1e-4, 1e-3, 1e-2, 5e-2):
            x = _x_boundary(a_null, alpha)
            row = f"      {alpha:10.3g} {x:15.6f} "
            for od in (0.2, 0.143, 0.1, 0.0345):
                a_t = T.a_of_od(od)
                row += f"{2*betainc(a_t, a_t, x):14.3e}"
            print(row)
    print()


def part_h():
    print("=" * 104)
    print("(h) RATIO-OF-SUMS (Jensen) TERM:  sd(D)/E[D] for D = sum scale*1{keep}, a_null=3.")
    print("    The E[N]/E[D] limit used above is exact to O((sd(D)/E[D])^2).")
    for name, d in D.items():
        s, t, w, k = seed_view(d)
        _, sc, nc = estimator_terms(s, t, w, k)
        n = np.rint(nc).astype(np.int64)
        for alpha, od in ((ALPHA_BONF, 0.2), (0.05, 0.2)):
            nmin = min_n_rejectable(3.0, alpha)
            ED = 0.0
            VD = 0.0
            for ni in np.unique(n[n >= 2]):
                cnt = int((n == ni).sum())
                scale = ni * (ni - 1) / 4.0
                pk = 1.0 if (nmin < 0 or ni < nmin) else T.seed_moments(int(ni), od, 3.0, alpha)[0]
                ED += cnt * scale * pk
                VD += cnt * scale**2 * pk * (1 - pk)
            print(f"    {name:10s} alpha={alpha:9.3g} od={od}: sd(D)/E[D] = "
                  f"{np.sqrt(VD)/ED:.3e}")
    print()


if __name__ == "__main__":
    part_e()
    part_g()
    part_h()
    part_f()
