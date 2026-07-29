"""IS 'REACH' WORTH IT?  Simulate the true physics and let the data answer.

The owner's challenge, in three parts:
  (1) near a TES fewer RNA fragments fit, so an observed fragment there is MORE likely gDNA — shouldn't
      gDNA get the advantage rather than RNA?
  (2) doesn't Sum(1/L) already compensate, since only short fragments fit and 1/L is larger for them?
  (3) is the extra machinery worth it at all?

TRUE PHYSICS simulated here, no estimator assumptions:
  * a mature molecule covers the WHOLE transcript; fragmentation then cuts it. A mature fragment is a
    uniformly-placed interval that must lie ENTIRELY inside the transcript.
  * gDNA is a uniformly-placed interval on the genome, unconstrained by the transcript.
  * rho_M and rho_g are both UNIFORM per bp. That is the null the estimator must recover.
"""
from __future__ import annotations
import numpy as np

RNG = np.random.default_rng(20260729)
T_LEN   = 3000      # transcript, single exon [0, T_LEN) inside a longer genome
GENOME  = 9000
FL_MU, FL_SD = 200.0, 50.0
RHO_M, RHO_G = 0.030, 0.010     # molecules per bp -- gDNA is 25% of molecules
N_REP   = 400

W = np.arange(1200)
PMF = np.exp(-0.5*((W-FL_MU)/FL_SD)**2); PMF[:20] = 0.0; PMF /= PMF.sum()
FLMEAN = float((W*PMF).sum())

def draw(n):
    return RNG.choice(W.size, size=n, p=PMF)

def trapezoid_E(Rd, Ra):
    """E_f[min(w-1, Rd, Ra, Rd+Ra-w+1)_+] -- the reach-aware crossing opportunity."""
    v = np.minimum.reduce([W-1.0, np.full(W.size, float(Rd)), np.full(W.size, float(Ra)),
                           Rd+Ra-W+1.0])
    return float((PMF*np.maximum(v, 0.0)).sum())

def simulate():
    """Return per-crossing-point (flux_M, recip_M, flux_g, recip_g) at every position in the transcript."""
    fM = np.zeros(T_LEN); rM = np.zeros(T_LEN); fG = np.zeros(T_LEN); rG = np.zeros(T_LEN)
    for _ in range(N_REP):
        # MATURE: fragmentation events at density rho_M per bp of transcript; an event at s of length w
        # yields a fragment only if it FITS inside the molecule. Drawing starts uniformly on [0, T) and
        # REJECTING the non-fitting ones is what makes the per-position start density exactly rho_M --
        # placing a fixed molecule count into the shrinking valid range instead inflates it by
        # T/(T-w+1) ~ 1.07 and puts a spurious +7% on every column.
        n = RNG.poisson(RHO_M * T_LEN)
        w = draw(n) + 1
        s = (RNG.random(n) * T_LEN).astype(np.int64)
        keep = (s >= 0) & (s + w <= T_LEN)
        for a, b in zip(s[keep], (s+w)[keep]):
            fM[a+1:b] += 1.0; rM[a+1:b] += 1.0/(b-a)      # crosses every interior point
        # gDNA: starts uniform on the genome, no transcript constraint
        n = RNG.poisson(RHO_G * GENOME)
        w = draw(n) + 1
        s = (RNG.random(n) * (GENOME - w + 1)).astype(np.int64)
        for a, b in zip(s, s+w):
            lo, hi = max(a+1, 0), min(b, T_LEN)
            if hi > lo:
                fG[lo:hi] += 1.0; rG[lo:hi] += 1.0/(b-a)
    return fM/N_REP, rM/N_REP, fG/N_REP, rG/N_REP

fM, rM, fG, rG = simulate()
print(f"transcript {T_LEN} bp, FL N({FL_MU:.0f},{FL_SD:.0f}) fl_mean={FLMEAN:.1f}, "
      f"rho_M={RHO_M} rho_g={RHO_G}  ({N_REP} reps)\n")

print("Distance from the TES, at a crossing point INSIDE the transcript:")
print(f"{'d(TES)':>7} {'true f_g':>9} | {'rhoM/fl_mean':>13} {'rhoM/E_J':>10} | "
      f"{'Sum1/L':>9} {'Sum1/L(cor)':>12} | {'reach Rd,Ra':>14}")
print("-"*88)
for d in (1500, 800, 400, 250, 200, 150, 100, 50, 20):
    p = T_LEN - d                      # the crossing point
    Rd, Ra = p, T_LEN - p              # exonic bases either side, within the transcript
    tot_f = fM[p] + fG[p]
    true_fg = fG[p] / tot_f if tot_f > 0 else np.nan
    naive  = fM[p] / FLMEAN            # count / fl_mean  -- NO reach
    EJ = trapezoid_E(Rd, Ra)
    aware  = fM[p] / EJ if EJ > 0 else np.nan   # count / E_J -- WITH reach
    s1L    = rM[p]                      # Sum(1/L): the owner's candidate self-correction
    corr   = float((PMF*np.maximum(np.minimum.reduce(
              [W-1.0, np.full(W.size,float(Rd)), np.full(W.size,float(Ra)), Rd+Ra-W+1.0]),0.0)/np.maximum(W,1)).sum())
    print(f"{d:>7} {true_fg:>9.3f} | {naive/RHO_M:>13.3f} {aware/RHO_M:>10.3f} | "
          f"{s1L/RHO_M:>9.3f} {s1L/max(corr,1e-9)/RHO_M:>12.3f} | {Rd:>6},{Ra:<7}")
print("\ncolumns are rho_hat / rho_TRUE for the MATURE component -- 1.000 is correct.")
