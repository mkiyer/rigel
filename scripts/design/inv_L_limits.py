"""Where does Sum 1/L STOP being exact? Two stress tests the earlier (infinite, uniform) sim could not see.

TEST 1 -- TRANSCRIPT END.  A molecule must FIT on its template. Near the end of a transcript the number of
admissible placements stops being L and becomes (distance to the end), the SAME for every L. The 1/L weight
then no longer cancels anything.
TEST 2 -- DENSITY STEP.  Sum 1/L at p averages rho over the molecules covering p, which reach back up to L
bases. Two components with different lengths average over different windows, so at a step the composition
is distorted.
"""

import numpy as np

rng = np.random.default_rng(11)
L = np.arange(1, 1201)


def fl(mu, sd):
    w = np.exp(-0.5 * ((L - mu) / sd) ** 2)
    return w / w.sum()


def placements(w, r_lo, r_hi):
    """admissible start positions for a length-w molecule COVERING a point with r_lo bases of template
    to its left (inclusive of the point) and r_hi to its right."""
    w = np.asarray(w, float)
    r_lo = np.broadcast_to(np.asarray(r_lo, float), w.shape)
    r_hi = np.broadcast_to(np.asarray(r_hi, float), w.shape)
    return np.clip(np.minimum.reduce([w, r_lo, r_hi, r_lo + r_hi - w]), 0, None)


print("TEST 1 -- distance from a transcript end (analytic, RNA N(200,50))")
f = fl(200, 50)
EinvL = float((f / L).sum())
print(
    f"  E[1/L] = {EinvL:.6f}   deep-interior Sum1/L / rho = {(f * placements(L, 10**6, 10**6) / L).sum():.4f}"
)
print(f"  {'dist to end':>12} {'Sum1/L / rho':>13} {'count/E[L] / rho':>18}")
EL = float((f * L).sum())
for d in (10, 20, 50, 100, 150, 200, 300, 500, 1000, 5000):
    pl = placements(L, 10**6, d)
    print(f"  {d:>12} {(f * pl / L).sum():>13.4f} {(f * pl).sum() / EL:>18.4f}")

print("\nTEST 1b -- SIMULATED, transcript of length 1099 (the human median), rho = 0.05")
T, RHO = 1099, 0.05
tot = np.zeros(T)
cnt = np.zeros(T)
for _ in range(4000):
    n = rng.poisson(RHO * T)
    lens = np.clip(rng.normal(200, 50, n).round().astype(int), 30, 600)
    ok = lens <= T
    lens = lens[ok]
    starts = rng.integers(0, T - lens + 1)
    for s, w in zip(starts, lens):
        tot[s : s + w] += 1.0 / w
        cnt[s : s + w] += 1.0
tot /= 4000
cnt /= 4000
for d in (10, 20, 50, 100, 200, 400, 549):
    p = T - d
    print(
        f"  {d:>5} bp from the 3' end: Sum1/L={tot[p]:.5f} ({tot[p] / RHO:.3f}x rho)   "
        f"count={cnt[p]:.3f} ({cnt[p] / (RHO * EL):.3f}x)"
    )

print(
    "\nTEST 2 -- DENSITY STEP. rho doubles at x=0. gDNA N(110,30), RNA N(200,50), true f_g = 0.30"
)
fg_dist, fr_dist = fl(110, 30), fl(200, 50)


def step_read(fdist, x, lo, hi):
    """E[Sum 1/L] at offset x from a step where rho = lo for pos<0 and hi for pos>=0."""
    out = 0.0
    for w, pw in zip(L, fdist):
        if pw <= 0:
            continue
        s = np.arange(x - w + 1, x + 1)  # starts covering x
        # a molecule spans [s, s+w); assign it the density at its START (where it was born)
        rho = np.where(s < 0, lo, hi)
        out += pw * rho.mean()
    return out


print(f"  {'offset':>8} {'rho_g read':>11} {'rho_r read':>11} {'f_g read':>9} {'error':>8}")
for x in (-300, -150, -60, -20, 0, 20, 60, 150, 300):
    rg = step_read(fg_dist, x, 0.30 * 1.0, 0.30 * 2.0)
    rr = step_read(fr_dist, x, 0.70 * 1.0, 0.70 * 2.0)
    fg = rg / (rg + rr)
    print(f"  {x:>8} {rg:>11.4f} {rr:>11.4f} {fg:>9.4f} {fg - 0.30:>+8.4f}")
