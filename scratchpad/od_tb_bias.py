"""Per-seed truncation bias r_n(od) = E[excess|kept]/E[excess], exact."""
from __future__ import annotations

import numpy as np
import od_tb_core as T

ALPHA_BONF = 1.0 / 160366.0


def min_n_rejectable(a_null: float, alpha: float, nmax: int = 200000) -> int:
    """Smallest n whose most extreme outcome (all one strand) has two-sided p < alpha."""
    lo, hi = 1, nmax
    if T.min_pvalue(nmax, a_null) >= alpha:
        return -1
    while lo < hi:
        mid = (lo + hi) // 2
        if T.min_pvalue(mid, a_null) < alpha:
            hi = mid
        else:
            lo = mid + 1
    return lo


def table(a_null, alphas, ns, ods):
    print(f"a_null = {a_null}  (ceiling od = {T.od_of_a(a_null):.4f})")
    for alpha in alphas:
        nmin = min_n_rejectable(a_null, alpha)
        print(f"\n  alpha = {alpha:.3g}   smallest rejectable seed n = {nmin}")
        hdr = "    n     " + "".join(f"  od={od:<6.3g}       " for od in ods)
        print(hdr)
        print("           " + "".join("  p_keep   r        " for _ in ods))
        for n in ns:
            row = f"  {n:6d}   "
            for od in ods:
                pk, e_all, e_kept, r = T.seed_moments(n, od, a_null, alpha)
                if od == 0:
                    # report absolute bias in od units instead of the 0/0 ratio
                    row += f"  {pk:.5f}  {e_kept/(n*(n-1)/4.):+8.1e}"
                else:
                    row += f"  {pk:.5f}  {r:8.5f} "
            print(row)


if __name__ == "__main__":
    ns = [2, 5, 10, 30, 100, 264, 300, 1000, 1523]
    ods = [0.0, 0.01, 0.0345, 0.1, 0.143, 0.2]
    for a_null in (2.0, 3.0):
        print("=" * 118)
        table(a_null, [ALPHA_BONF, 1e-4, 1e-3, 0.01, 0.05], ns, ods)
        print()
