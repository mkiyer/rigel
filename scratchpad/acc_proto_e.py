"""EXPERIMENT E — the owner's geometry, worked end to end.

    A 500bp <-> B 50bp <-> C 10bp <-> D 1bp <-> E 1000bp

Three estimators of the SAME true density at each edge:

  v3        rho = flux / E[w]                       integer count, region-free divisor
  today     rho = (m_right[r] + m_left[r+1]) / (0.5*(sidelen_r + sidelen_r1))
                                                    the shipped pooled-seam mass estimator
  split     rho = (sum_f w_f / k_f) / E[w]          "partition the length across the k edges
                                                     the fragment spans", geometry-blind divisor

Truth is uniform rho, so every edge must read the same number. Anything that does not is
carrying the local node geometry into its answer.
"""

from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(20260730)

LENS = np.array([500, 50, 10, 1, 1000], dtype=np.int64)
NAMES = ["A", "B", "C", "D", "E"]
CUTS = np.concatenate([[0], np.cumsum(LENS)])          # 0, 500, 550, 560, 561, 1561
GENOME = int(CUTS[-1])
EDGES = [(i, i + 1) for i in range(len(LENS) - 1)]     # A|B, B|C, C|D, D|E


def fl(mu, sd, n=1500):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return p / p.sum()


def moments(pmf):
    w = np.arange(pmf.size, dtype=np.float64)
    return np.sum(pmf * w), np.sum(pmf * w * w)


def side_len(pmf, L):
    """boundary_side_eff_length = E[min(w, L)] / 2 — today's per-side DENSITY divisor."""
    w = np.arange(pmf.size, dtype=np.float64)
    return 0.5 * np.sum(pmf * np.minimum(w, L))


def region_eff(pmf, L):
    """region_eff_length = E[max(0, L - w + 1)]."""
    w = np.arange(pmf.size, dtype=np.float64)
    return np.sum(pmf * np.maximum(0.0, L - w + 1.0))


def node_of(pos):
    return np.searchsorted(CUTS, pos, side="right") - 1


def simulate(rho, pmf, n_edges=len(EDGES)):
    """Deposit under all three rules at once. Returns per-edge arrays."""
    flux = np.zeros(n_edges)
    flux_len = np.zeros(n_edges)
    split_len = np.zeros(n_edges)
    mass_r = np.zeros(len(LENS))     # today: region r's RIGHT-side boundary mass
    mass_l = np.zeros(len(LENS))     # today: region r's LEFT-side boundary mass
    contained = np.zeros(len(LENS))
    contained_len = np.zeros(len(LENS))

    n = RNG.poisson(rho * GENOME)
    s = RNG.integers(0, GENOME, n)
    w = RNG.choice(pmf.size, size=n, p=pmf)
    keep = s + w <= GENOME
    s, w = s[keep], w[keep]

    for si, wi in zip(s, w):
        a, z = node_of(si), node_of(si + wi - 1)
        if a == z:
            contained[a] += 1
            contained_len[a] += wi
            continue
        crossed = list(range(a, z))          # edge i is between node i and i+1
        k = len(crossed)
        for e in crossed:
            flux[e] += 1
            flux_len[e] += wi
            split_len[e] += wi / k
        # --- today's rule: per-region slices, share = (slice/L)/n_cross ---
        slices = []
        for r in range(a, z + 1):
            lo, hi = max(si, CUTS[r]), min(si + wi, CUTS[r + 1])
            slices.append((r, hi - lo))
        m = len(slices)
        for i, (r, ln) in enumerate(slices):
            cl, cr = i > 0, i < m - 1
            ncross = int(cl) + int(cr)
            share = (ln / wi) / ncross
            if cr:
                mass_r[r] += share       # deposited on r's RIGHT boundary, r's side
            if cl:
                mass_l[r] += share       # deposited on r's LEFT boundary, r's side
    return flux, flux_len, split_len, mass_l, mass_r, contained, contained_len


def main():
    pmf = fl(200, 50)
    E_w, E_w2 = moments(pmf)
    rho = 0.05
    print("geometry: " + "  ".join(f"{n}={L}" for n, L in zip(NAMES, LENS))
          + f"   genome {GENOME} bp")
    print(f"FL: mean {E_w:.1f}, E[w^2]/E[w] = {E_w2/E_w:.1f}   |   TRUE rho = {rho}\n")

    flux, flux_len, split_len, mass_l, mass_r, cont, cont_len = simulate(rho, pmf)

    sl = np.array([side_len(pmf, L) for L in LENS])
    print("PER-EDGE DENSITY (truth = 0.0500 everywhere)")
    print(f"  {'edge':6} {'flux':>8} {'flux_len':>10} | {'v3':>8} {'today':>8} {'split':>8}")
    for e, (r, r1) in enumerate(EDGES):
        seam_m = mass_r[r] + mass_l[r1]
        seam_S = 0.5 * (sl[r] + sl[r1])
        v3 = flux[e] / E_w
        today = seam_m / seam_S
        split = split_len[e] / E_w
        print(f"  {NAMES[r]}|{NAMES[r1]:3} {flux[e]:8.0f} {flux_len[e]:10.0f} | "
              f"{v3:8.4f} {today:8.4f} {split:8.4f}")

    def spread(v):
        v = np.asarray(v)
        return v.max() / max(v.min(), 1e-12)

    v3s = np.array([flux[e] / E_w for e in range(len(EDGES))])
    tds = np.array([(mass_r[r] + mass_l[r1]) / (0.5 * (sl[r] + sl[r1]))
                    for r, r1 in EDGES])
    sps = np.array([split_len[e] / E_w for e in range(len(EDGES))])
    print(f"\n  max/min across edges (1.00 = geometry-free):  "
          f"v3 {spread(v3s):.2f}   today {spread(tds):.2f}   split {spread(sps):.2f}")

    print("\nPER-NODE CONTAINED (truth = 0.0500)")
    for i, nm in enumerate(NAMES):
        Er = region_eff(pmf, LENS[i])
        est = cont[i] / Er if Er > 1e-9 else float("nan")
        print(f"  node {nm} L={LENS[i]:5d}  contained {cont[i]:7.0f}  E_cont {Er:8.2f}  "
              f"rho {est:8.4f}" + ("   <- no contained opportunity" if Er < 1 else ""))

    # ---- two components: does (count, sum-length) recover f_g at each edge? ----
    print("\nTWO COMPONENTS at each edge — gDNA FL 60, RNA FL 200, true f_g = 0.286")
    pg, pr = fl(60, 15), fl(200, 50)
    rg, rr = 0.02, 0.05
    Eg1, Eg2 = moments(pg)
    Er1, Er2 = moments(pr)
    M = np.array([[Eg1, Er1], [Eg2, Er2]])
    fg_true = rg / (rg + rr)
    acc = {e: [] for e in range(len(EDGES))}
    for _ in range(200):
        f1, l1, *_ = simulate(rg, pg)
        f2, l2, *_ = simulate(rr, pr)
        for e in range(len(EDGES)):
            obs = np.array([f1[e] + f2[e], l1[e] + l2[e]])
            sol = np.maximum(np.linalg.solve(M, obs), 0.0)
            if sol.sum() > 0:
                acc[e].append(sol[0] / sol.sum())
    for e, (r, r1) in enumerate(EDGES):
        v = np.asarray(acc[e])
        print(f"  {NAMES[r]}|{NAMES[r1]:3}  f_g {v.mean():.3f} +/- {v.std():.3f}   "
              f"bias {v.mean()-fg_true:+.4f}")


if __name__ == "__main__":
    main()
