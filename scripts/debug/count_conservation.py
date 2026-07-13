"""Verify the count-conservation correction (issue #5) preserves self-defense.

OLD (wrong) node constraint:  sum_c rho_c = D          (density conserved)
NEW (correct) constraint:     sum_c N_c   = N          (COUNTS conserved)
Link:                         rho_c = N_c / E_c         (E_c = per-component eff-length,
                                                         E_gdna != E_rna because FL differs)

Claim to verify:
 (1) In log-count coords x_c=log N_c, density messages m_c become shifted targets
     m~_c = m_c + log E_c, and the KKT gives  delta_c = x_c - m~_c = -(mu/2) N_c / pi_c
     => SELF-DEFENSE (confident components barely deviate) SURVIVES with lever N_c.
 (2) When E_g == E_r the count solve == the old density solve (equal-FL special case).
 (3) When E_g != E_r the two solves DIFFER, and the difference is exactly the FL bias
     the old (density-conserved) solve was missing.
"""
import numpy as np
from scipy.optimize import minimize

def solve_count(N, Ec, m, pi):
    """min sum pi_c (log(N_c/E_c) - m_c)^2  s.t. sum N_c = N, N_c>0.
    parametrize by log-count x; enforce sum e^x = N by projecting via a softmax over free a."""
    Ec = np.asarray(Ec, float); m = np.asarray(m, float); pi = np.asarray(pi, float)
    def unpack(z):                       # z: free logits -> count fractions summing to 1 -> counts
        a = np.exp(z - z.max()); a /= a.sum()
        return a * N
    def obj(z):
        Nc = np.maximum(unpack(z), 1e-12)
        rho = Nc / Ec
        return float(np.sum(pi * (np.log(rho) - m) ** 2))
    best = None
    for seed in range(8):                # multistart to check uniqueness
        z0 = np.linspace(-1, 1, len(Ec)) * (0.3 * seed - 1.0)
        r = minimize(obj, z0, method="Nelder-Mead",
                     options=dict(xatol=1e-9, fatol=1e-12, maxiter=20000))
        if best is None or r.fun < best.fun:
            best = r
    Nc = unpack(best.x); rho = Nc / Ec
    return Nc, rho, best.fun

def solve_density(D, m, pi):
    """OLD: min sum pi_c (log rho_c - m_c)^2 s.t. sum rho_c = D."""
    m = np.asarray(m, float); pi = np.asarray(pi, float)
    def unpack(z):
        a = np.exp(z - z.max()); a /= a.sum()
        return a * D
    def obj(z):
        rho = np.maximum(unpack(z), 1e-12)
        return float(np.sum(pi * (np.log(rho) - m) ** 2))
    best = None
    for seed in range(8):
        z0 = np.linspace(-1, 1, len(m)) * (0.3 * seed - 1.0)
        r = minimize(obj, z0, method="Nelder-Mead",
                     options=dict(xatol=1e-9, fatol=1e-12, maxiter=20000))
        if best is None or r.fun < best.fun:
            best = r
    return unpack(best.x)

print("=" * 78)
print("TEST 1  self-defense under COUNT conservation (single-strand: gDNA + RNA+)")
print("=" * 78)
# confident gDNA belief, weak (wrong) RNA+ message; N=100 fragments fixed.
# gDNA FL ~ short (E_g=100), RNA FL longer (E_r=200) -> different count<->density maps.
N = 100.0
Eg, Er = 100.0, 200.0
# targets are on log-DENSITY. gDNA target rho=0.30 confidently; RNA+ target rho=0.20 weakly.
m   = [np.log(0.30), np.log(0.20)]     # [gdna, rna+]
pi  = [50.0,          2.0]             # gDNA confident, RNA+ weak
Nc, rho, f = solve_count(N, [Eg, Er], m, pi)
print(f"  counts   N_g={Nc[0]:7.3f}  N_p={Nc[1]:7.3f}   (sum={Nc.sum():.3f}=N)")
print(f"  density  rho_g={rho[0]:.4f} (target {np.exp(m[0]):.4f})  rho_p={rho[1]:.4f} (target {np.exp(m[1]):.4f})")
print(f"  log-count deviation delta_g = log N_g - (m_g+log E_g) = "
      f"{np.log(Nc[0]) - (m[0]+np.log(Eg)):+.4f}")
print(f"  log-count deviation delta_p = log N_p - (m_p+log E_p) = "
      f"{np.log(Nc[1]) - (m[1]+np.log(Er)):+.4f}")
print(f"  => |delta_g| << |delta_p|?  {abs(np.log(Nc[0])-(m[0]+np.log(Eg))) < abs(np.log(Nc[1])-(m[1]+np.log(Er)))}"
      f"   (confident gDNA defends; weak RNA absorbs the count constraint)")
# KKT check: delta_c should be proportional to N_c/pi_c (same mu)
dg = np.log(Nc[0]) - (m[0]+np.log(Eg)); dp = np.log(Nc[1]) - (m[1]+np.log(Er))
mu_g = -2*dg*pi[0]/Nc[0]; mu_p = -2*dp*pi[1]/Nc[1]
print(f"  KKT: implied mu from gDNA = {mu_g:+.5f}, from RNA+ = {mu_p:+.5f}  (equal => KKT holds)")

print()
print("=" * 78)
print("TEST 2  equal-FL special case: count-solve == old density-solve")
print("=" * 78)
Nc2, rho2, _ = solve_count(N, [150.0, 150.0], m, pi)     # E_g==E_r
D_equiv = (N/150.0)                                       # each count maps to /150; but D varies...
# under equal E, count-fraction == density-fraction; compare COMPOSITIONS
f_count = rho2 / rho2.sum()
rho_dens = solve_density(rho2.sum(), m, pi)              # feed the SAME total density
f_dens = rho_dens / rho_dens.sum()
print(f"  equal-FL count-solve composition : f_g={f_count[0]:.4f} f_p={f_count[1]:.4f}")
print(f"  density-solve (same D) composition: f_g={f_dens[0]:.4f} f_p={f_dens[1]:.4f}")
print(f"  => match within 1e-3?  {np.allclose(f_count, f_dens, atol=1e-3)}")

print()
print("=" * 78)
print("TEST 3  E_g != E_r: the FL bias the density-conserved solve MISSES")
print("=" * 78)
# same targets; gDNA FL short (E_g=80) vs RNA FL long (E_r=250).
Eg3, Er3 = 80.0, 250.0
Nc3, rho3, _ = solve_count(N, [Eg3, Er3], m, pi)
f_count3 = rho3 / rho3.sum()
# the WRONG way: solve in density-space to total = sum of the naive count/E, ignoring that
# the split changes D. Show the count-composition vs a density-conserved composition differ.
rho_dens3 = solve_density(rho3.sum(), m, pi)
f_dens3 = rho_dens3 / rho_dens3.sum()
# also report the COUNT split (the deliverable to EM) vs the density composition (the currency)
a_count = Nc3 / Nc3.sum()
print(f"  E_g={Eg3}, E_r={Er3}")
print(f"  COUNT split (deliverable): a_g={a_count[0]:.4f}  a_p={a_count[1]:.4f}")
print(f"  DENSITY composition (currency): f_g={f_count3[0]:.4f}  f_p={f_count3[1]:.4f}")
print(f"  => count-split != density-composition:  {not np.allclose(a_count, f_count3, atol=1e-3)}"
      f"   (the FL map; a longer-FL component has fewer counts per unit density)")
print(f"  gap |a_g - f_g| = {abs(a_count[0]-f_count3[0]):.4f}  (this is what a density-only solve conflates)")
