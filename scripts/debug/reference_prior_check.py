"""Numerically verify the per-node reference/measure decomposition on the (lambda, tau) grid.

CLAIM under test:
    log p(lam, tau) = sum_c logP_c(log rho_c) + R
where logP_c is each component's LOG-RATE prior density, and R is a pure measure residual:
    R = 0                       for n=2 (single-strand: g, r)
    R = -log((1-tau^2)/4)       for n=3 (AMBIG: g, +, -)
and the Jeffreys reference for an unfitted component is logP_c = +0.5*log(f_c) + const.

Consequences to check:
  A) n=2, both Jeffreys  -> psi_lam = 0.5*log f_g + 0.5*log(1-f_g)   == Beta(1/2,1/2) in f_g
  B) n=3, all Jeffreys   -> psi = 0.5*log f_g + log(1-f_g) - 0.5*log((1-tau^2)/4)
                              == Dirichlet(1/2,1/2,1/2) on (f_g,f_+,f_-)
  C) n=2, gDNA FITTED    -> extra term is +0.5*log(1-f_g)     (NO +0.5*log f_g)
  D) n=3, gDNA FITTED    -> extra term is +log(1-f_g) - 0.5*log((1-tau^2)/4)
"""

import numpy as np

rng = np.random.default_rng(0)


def sigma(x):
    return 1.0 / (1.0 + np.exp(-x))


# ---------------------------------------------------------------- A) n=2 Beta(1/2,1/2)
# Sample f_g ~ Beta(1/2,1/2), map to lam, and check the empirical density in lam
# equals exp(0.5*log f_g + 0.5*log(1-f_g)) up to normalization.
N = 4_000_000
fg = rng.beta(0.5, 0.5, N)
lam = np.log(fg) - np.log1p(-fg)
hist, edges = np.histogram(lam, bins=400, range=(-10, 10), density=True)
ctr = 0.5 * (edges[:-1] + edges[1:])
f_c = sigma(ctr)
psi_A = 0.5 * np.log(f_c) + 0.5 * np.log(1 - f_c)
pred_A = np.exp(psi_A)
pred_A /= np.trapezoid(pred_A, ctr)
# renormalize hist over the same window
hist_n = hist / np.trapezoid(hist, ctr)
errA = np.max(np.abs(hist_n - pred_A)) / np.max(pred_A)
print(f"A) n=2 both-Jeffreys  == Beta(1/2,1/2) in lam : max rel err = {errA:.4f}")

# ---------------------------------------------------------------- B) n=3 Dirichlet(1/2,1/2,1/2)
d = rng.dirichlet([0.5, 0.5, 0.5], N)
fg3, fp3, fn3 = d[:, 0], d[:, 1], d[:, 2]
lam3 = np.log(fg3) - np.log1p(-fg3)
tau3 = (fp3 - fn3) / np.maximum(fp3 + fn3, 1e-300)
# 2-D histogram vs the predicted psi
H, xe, ye = np.histogram2d(
    lam3, tau3, bins=[120, 120], range=[[-6, 6], [-0.995, 0.995]], density=True
)
xc = 0.5 * (xe[:-1] + xe[1:])
yc = 0.5 * (ye[:-1] + ye[1:])
FG = sigma(xc)[:, None]
TAU = yc[None, :]
psi_B = 0.5 * np.log(FG) + np.log(1 - FG) - 0.5 * np.log((1 - TAU**2) / 4.0)
predB = np.exp(psi_B)
dx = xc[1] - xc[0]
dy = yc[1] - yc[0]
predB /= predB.sum() * dx * dy
Hn = H / (H.sum() * dx * dy)
errB = np.max(np.abs(Hn - predB)) / np.max(predB)
print(f"B) n=3 all-Jeffreys   == Dirichlet(1/2,1/2,1/2) in (lam,tau) : max rel err = {errB:.4f}")

# ---------------------------------------------------------------- naive mnemonic comparison (n=3)
psi_naive = 0.5 * np.log(FG) + np.log(1 - FG) + 0.5 * np.log((1 - TAU**2) / 4.0)
predN = np.exp(psi_naive)
predN /= predN.sum() * dx * dy
errN = np.max(np.abs(Hn - predN)) / np.max(predB)
print(f"   naive '+0.5*log f_c per component, no measure residual' : max rel err = {errN:.4f}")
print("   -> the naive mnemonic is WRONG for n=3 (tau-term sign) iff errN >> errB")

# ---------------------------------------------------------------- C/D) symbolic-numeric check of the residual
# Verify: sum_c (-log rho_c) + log|J_chart| == -R, i.e. the measure bookkeeping.
# n=2: chart (f_g); |J| = f_g(1-f_g).  sum_c log f_c = log f_g + log(1-f_g). residual = 0.
# n=3: chart (f_g,f_+); |J| = f_g(1-f_g)^2/2. sum_c log f_c = log f_g + 2log(1-f_g) + log((1-tau^2)/4).
for name, fgv, tauv in [("mid", 0.3, 0.4), ("tilted", 0.7, -0.85), ("balanced", 0.5, 0.0)]:
    u = 1 - fgv
    fp = u * (1 + tauv) / 2
    fn = u * (1 - tauv) / 2
    # n=2
    J2 = fgv * (1 - fgv)
    res2 = np.log(J2) - (np.log(fgv) + np.log(1 - fgv))
    # n=3
    J3 = fgv * (1 - fgv) ** 2 / 2
    res3 = np.log(J3) - (np.log(fgv) + np.log(fp) + np.log(fn))
    want3 = -np.log((1 - tauv**2) / 4.0)
    print(
        f"C/D) {name:9s} residual n=2 = {res2:+.6f} (want 0)   "
        f"n=3 = {res3:+.6f} vs -log((1-t^2)/4) = {want3:+.6f}  "
        f"[match={abs(res3 - want3) < 1e-9}]"
    )

# ---------------------------------------------------------------- what the FITTED-gDNA case adds
# n=2, gDNA fitted: extra = +0.5*log(1-f_g). n=3, gDNA fitted: extra = +log(1-f_g) -0.5*log((1-t^2)/4)
print()
print("With logP_g FITTED (NPMLE) and RNA unfitted (Jeffreys), the term to ADD to psi is:")
print("  1-DOF single-strand : +0.5*log(1-f_g)                        [no +0.5*log f_g!]")
print("  2-DOF AMBIG         : +log(1-f_g) - 0.5*log((1-tau^2)/4)")
print("Both push f_g DOWN (bound the f_g->1 vertex). Neither contains +0.5*log f_g.")
print()
# Beta(1/2,1) sanity: the 1-DOF fitted-gDNA reference ALONE (flat logP_g) -> p(f) ~ f^-1 (1-f)^-1/2 ... check
lamg = np.linspace(-12, 12, 20001)
fgg = sigma(lamg)
p = np.exp(0.5 * np.log(1 - fgg))
p /= np.trapezoid(p, lamg)
cdf = np.cumsum(p) * (lamg[1] - lamg[0])
med = fgg[np.searchsorted(cdf, 0.5)]
print(f"1-DOF, flat logP_g + (+0.5*log(1-f_g)) alone on L=12: median f_g = {med:.4f}")
p0 = np.ones_like(lamg)
p0 /= np.trapezoid(p0, lamg)
cdf0 = np.cumsum(p0) * (lamg[1] - lamg[0])
print(
    f"   vs the SHIPPED Haldane (nothing written)         : median f_g = {fgg[np.searchsorted(cdf0, 0.5)]:.4f}"
)
print(
    "   (Haldane's median is 0.5 only by grid symmetry; it is improper and L-driven once anything tilts.)"
)
