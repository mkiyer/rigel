"""Is the shipped POOLED-SEAM density m_s/S_s equal to rho, or to 2*rho?

Drives the current deposit rule (docs/accumulator/00_design.md 4.3) under a uniform
genomic density and compares the pooled-seam density against the truth, using the
REPO's own effective-length functions.

  m_s = mass_gdna_right[r] + mass_gdna_left[r+1]          (_pooled_seam_arrays)
  S_s = 0.5 * (gdna_boundary_len[r] + gdna_boundary_len[r+1])
      where gdna_boundary_len = boundary_side_eff_length = E[min(l,L)]/2   (calibrate.py:226)
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.effective_length import boundary_side_eff_length

RNG = np.random.default_rng(7)


def fl(mu, sd, n=1200):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return p / p.sum()


def run(region_len, n_reg, pmf, rho=0.05):
    cuts = np.arange(n_reg + 1, dtype=np.int64) * region_len
    genome = int(cuts[-1])
    mass_r = np.zeros(n_reg)   # region r's RIGHT-side boundary mass
    mass_l = np.zeros(n_reg)   # region r's LEFT-side boundary mass
    flux = np.zeros(n_reg)     # integer crossings at the boundary right of region r

    n = RNG.poisson(rho * genome)
    s = RNG.integers(0, genome, n)
    w = RNG.choice(pmf.size, size=n, p=pmf)
    keep = s + w <= genome
    s, w = s[keep], w[keep]

    for si, wi in zip(s, w):
        a = int(np.searchsorted(cuts, si, "right") - 1)
        z = int(np.searchsorted(cuts, si + wi - 1, "right") - 1)
        if a == z:
            continue
        slices = []
        for r in range(a, z + 1):
            lo, hi = max(si, cuts[r]), min(si + wi, cuts[r + 1])
            slices.append((r, hi - lo))
        m = len(slices)
        for i, (r, ln) in enumerate(slices):
            cl, cr = i > 0, i < m - 1
            share = (ln / wi) / (int(cl) + int(cr))
            if cr:
                mass_r[r] += share
                flux[r] += 1
            if cl:
                mass_l[r] += share

    side = boundary_side_eff_length(pmf, np.full(n_reg, float(region_len)))  # E[min(l,L)]/2
    w_arr = np.arange(pmf.size, dtype=np.float64)
    fl_mean = float(np.sum(pmf * w_arr))

    interior = np.arange(2, n_reg - 3)          # avoid genome edges
    m_s = mass_r[interior] + mass_l[interior + 1]
    S_s_code = 0.5 * (side[interior] + side[interior + 1])   # what _pooled_seam_arrays computes
    S_s_sum = side[interior] + side[interior + 1]            # the sum, not the average
    fx = flux[interior]

    print(f"  L={region_len:5d}  FL={fl_mean:5.1f}  seams={interior.size:4d}   true rho={rho}")
    print(f"      m_s / S_s   (SHIPPED, averaged)  = {np.mean(m_s / S_s_code):.4f}"
          f"   -> ratio to truth {np.mean(m_s / S_s_code) / rho:.3f}")
    print(f"      m_s / (S_l+S_r) (summed)         = {np.mean(m_s / S_s_sum):.4f}"
          f"   -> ratio to truth {np.mean(m_s / S_s_sum) / rho:.3f}")
    print(f"      flux / fl_mean  (v3)             = {np.mean(fx / fl_mean):.4f}"
          f"   -> ratio to truth {np.mean(fx / fl_mean) / rho:.3f}")


if __name__ == "__main__":
    print("Uniform genomic density; does the pooled-seam density recover rho?\n")
    for L in (2000, 500, 200, 100, 50):
        run(L, 400, fl(200, 50))
        print()
