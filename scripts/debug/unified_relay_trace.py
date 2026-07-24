"""UNIFIED-SOLVER RELAY TRACE — where along the chain does the message amplitude leave the simplex?

The fused message's implied composition must sum to 1 (`unified_solver_design.md` §2). It does not
(`UNIFIED_SOLVER_HANDOFF.md`). This probe walks the forward relay hop by hop and reports the **drift**

    D(i) = [ rho_g(i)*E_g(i) + (rho_p(i)+rho_n(i))*E_r(i) ] / M(i)      (must be ~1)

per node, plus the per-hop ratio D(i)/D(src) — which isolates WHICH edge type injects the uncancelled
factor. A telescoping relay keeps D constant along a path; anything else names the defect.

    RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/unified_relay_trace.py [cond]
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from unified_message_audit import run  # noqa: E402

_EPS = 1e-12
_NAMES = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cond = args[0] if args else "gdna_gdna300_ss_0.50_nrna_present_capture_on"
    d = run(cond)
    cap, cls, chain = d["cap"], d["cls"], d["chain"]
    st = cap["_uni_static"]
    M, E_g, E_r = st["M"], st["E_g"], st["E_r"]
    left = np.asarray(chain.left)
    order = [int(x) for x in np.asarray(chain.order)]

    def drift(g, p, n):
        return (g * E_g + (p + n) * E_r) / np.maximum(M, _EPS)

    D_own = drift(st["og"], st["op"], st["on"])
    D_fwd = drift(st["fwd_g"], st["fwd_p"], st["fwd_n"])
    D_g_own = st["og"] * E_g / np.maximum(M, _EPS)
    D_g_fwd = st["fwd_g"] * E_g / np.maximum(M, _EPS)
    ok = M > 1e-9

    print(f"# cond={cond}")
    print("\n=== drift D (implied SUM f_c; must be ~1) ===")
    print(f"{'class':<11}{'n':>6}{'D_own':>9}{'D_fwd':>9}{'Dg_own':>9}{'Dg_fwd':>9}")
    for c in (0, 1, 2, 3):
        m = ok & (cls == c)
        if not m.sum():
            continue
        w = M[m]
        print(f"{_NAMES[c]:<11}{int(m.sum()):>6}"
              f"{np.average(D_own[m], weights=w):>9.3f}{np.average(D_fwd[m], weights=w):>9.3f}"
              f"{np.average(D_g_own[m], weights=w):>9.3f}{np.average(D_g_fwd[m], weights=w):>9.3f}")

    # ---- PER-HOP gDNA amplification along the forward relay: does D_g telescope? ----
    # For a node that ADOPTS (own precision 0) the transported density, D_g(i)/D_g(src) isolates the
    # uncancelled factor of THIS edge. Bucket by (src class -> dst class).
    pg_own, fwd_pg = st["pg_own"], st["fwd_pg"]
    rows: dict[tuple[int, int], list[float]] = {}
    adopt: dict[tuple[int, int], list[float]] = {}
    for i in order:
        s = int(left[i])
        if s < 0 or M[i] <= 1e-9 or M[s] <= 1e-9:
            continue
        if D_g_fwd[s] <= _EPS or D_g_fwd[i] <= _EPS:
            continue
        k = (int(cls[s]), int(cls[i]))
        rows.setdefault(k, []).append(D_g_fwd[i] / D_g_fwd[s])
        if pg_own[i] <= 0.0:  # pure adopt — no own term to explain the change
            adopt.setdefault(k, []).append(D_g_fwd[i] / D_g_fwd[s])

    print("\n=== per-hop gDNA drift amplification  D_g(dst)/D_g(src), forward relay ===")
    print(f"{'edge':<26}{'n':>6}{'geomean':>10}{'median':>9}{'p10':>9}{'p90':>9}{'n_adopt':>9}{'geo_adopt':>11}")
    for k in sorted(rows, key=lambda z: -len(rows[z])):
        v = np.asarray(rows[k])
        a = np.asarray(adopt.get(k, []))
        lab = f"{_NAMES[k[0]]} -> {_NAMES[k[1]]}"
        ga = np.exp(np.mean(np.log(np.maximum(a, 1e-12)))) if a.size else np.nan
        print(f"{lab:<26}{v.size:>6}{np.exp(np.mean(np.log(np.maximum(v, 1e-12)))):>10.3f}"
              f"{np.median(v):>9.3f}{np.percentile(v, 10):>9.3f}{np.percentile(v, 90):>9.3f}"
              f"{a.size:>9}{ga:>11.3f}")

    # ---- Walk the single longest adopt-chain from an anchor, printing each hop ----
    anchor = None
    best = 0
    seen_len = np.zeros(M.size, int)
    for i in order:
        s = int(left[i])
        seen_len[i] = (seen_len[s] + 1) if (s >= 0 and pg_own[i] <= 0.0) else 0
        if seen_len[i] > best:
            best, anchor = seen_len[i], i
    if anchor is not None:
        path = []
        j = anchor
        while j >= 0 and len(path) < 24:
            path.append(j)
            j = int(left[j])
        path.reverse()
        print(f"\n=== longest zero-own-precision forward path (len {best}), hop by hop ===")
        print(f"{'node':>6}{'class':<12}{'M':>10}{'E_g':>9}{'E_r':>9}{'rho_tot':>11}"
              f"{'rho_L':>11}{'rho_R':>11}{'relay_rg':>11}{'D_g':>8}{'r_hop':>9}")
        rho0, rl0, rr0 = st["rho_node0"], st["rho_l0"], st["rho_r0"]
        prev = None
        for j in path:
            rh = (rl0[j] / rr0[prev]) if (prev is not None and rr0[prev] > _EPS) else np.nan
            print(f"{j:>6}{_NAMES[int(cls[j])]:<12}{M[j]:>10.1f}{E_g[j]:>9.1f}{E_r[j]:>9.1f}"
                  f"{rho0[j]:>11.4f}{rl0[j]:>11.4f}{rr0[j]:>11.4f}"
                  f"{st['fwd_g'][j]:>11.4f}{D_g_fwd[j]:>8.2f}{rh:>9.3f}")
            prev = j

    # ---- how much of the chain is zero-own-precision? ----
    print(f"\n# nodes with pg_own == 0: {int((pg_own[ok] <= 0).sum())}/{int(ok.sum())}"
          f"   fwd_pg == 0: {int((fwd_pg[ok] <= 0).sum())}/{int(ok.sum())}")
    print(f"# mean adopt-run length (consecutive zero-own-precision): {seen_len[ok].mean():.2f}"
          f"  max {seen_len.max()}")


if __name__ == "__main__":
    main()
