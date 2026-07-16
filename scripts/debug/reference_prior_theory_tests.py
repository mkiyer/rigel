"""THEORY-SPACE tests of the psi reference. No solver, no messages, no benchmark — pure measure theory.

Each test is a property the reference must satisfy REGARDLESS of whether the rest of calibration works.

  T1  transform consistency : does the (lambda,tau) reference reproduce the intended composition prior?
  T2  properness            : is int exp(psi_ref) finite? (equivalently: is the answer L-invariant?)
  T3  reparam invariance    : same posterior whether we grid lambda or grid f_g?
  T4  CLASS CONSISTENCY     : does a single-strand node and an AMBIG node carry the SAME f_g reference?
  T5  reduction             : does the 2-DOF reference reduce to the 1-DOF one as a strand dies?
"""

import numpy as np
from scipy.stats import beta as B
from scipy.stats import dirichlet as D

rng = np.random.default_rng(7)
N = 3_000_000


def sig(x):
    return 1.0 / (1.0 + np.exp(-x))


print("=" * 78)
print("T4/T5  CLASS CONSISTENCY — the f_g marginal of the 2-DOF reference vs the 1-DOF reference")
print("=" * 78)
print("Dirichlet(a1,a2,a3) marginal of x1 is Beta(a1, a2+a3)  [standard aggregation property]")
print()
print(
    f"{'2-DOF reference':<26}{'f_g marginal':<18}{'1-DOF ref is Beta(1/2,1/2)':<28}{'consistent?'}"
)
print("-" * 78)
for name, a in [
    ("Dirichlet(1/2,1/2,1/2)", (0.5, 0.5, 0.5)),
    ("Dirichlet(1/2,1/4,1/4)", (0.5, 0.25, 0.25)),
]:
    marg = f"Beta({a[0]}, {a[1] + a[2]})"
    ok = "YES" if abs((a[1] + a[2]) - 0.5) < 1e-12 else "NO  <-- class-dependent f_g prior"
    print(f"{name:<26}{marg:<18}{'Beta(0.5, 0.5)':<28}{ok}")
print()

# empirical confirmation of the aggregation property
for name, a in [("Dir(1/2,1/2,1/2)", (0.5, 0.5, 0.5)), ("Dir(1/2,1/4,1/4)", (0.5, 0.25, 0.25))]:
    x = D.rvs(np.array(a), size=N, random_state=rng)[:, 0]
    q = np.percentile(x, [10, 50, 90])
    q_half = B.ppf([0.1, 0.5, 0.9], 0.5, 0.5)
    print(f"  {name:<18} f_g quantiles(10/50/90) = {q[0]:.4f} {q[1]:.4f} {q[2]:.4f}")
print(
    f"  {'Beta(1/2,1/2)':<18} f_g quantiles(10/50/90) = {q_half[0]:.4f} {q_half[1]:.4f} {q_half[2]:.4f}"
)
print("  -> Dir(1/2,1/4,1/4) must match Beta(1/2,1/2); Dir(1/2,1/2,1/2) must NOT.")
print()

print("=" * 78)
print("T1  TRANSFORM CONSISTENCY — the (lambda,tau) form of each reference")
print("=" * 78)
print(
    "General identity (derived):  log p(lam,tau) = SUM_c logP_c  -  SUM_c log f_c  +  log|J_chart|"
)
print(
    "  n=2: |J| = f_g(1-f_g)          SUM log f_c = log f_g + log(1-f_g)        -> residual R = 0"
)
print("  n=3: |J| = f_g(1-f_g)^2 / 2    SUM log f_c = log f_g + 2log(1-f_g) + log((1-t^2)/4)")
print(
    "                                                        -> residual R = -log((1-t^2)/4) - log2"
)
print()


def lam_tau_form(a_g, a_p, a_n, fg, tau):
    """Closed form of log Dirichlet(a_g,a_p,a_n) as a density in (lambda, tau), up to a constant."""
    u = 1.0 - fg
    fp = u * (1 + tau) / 2
    fn = u * (1 - tau) / 2
    return (
        (a_g - 1) * np.log(fg)
        + (a_p - 1) * np.log(fp)
        + (a_n - 1) * np.log(fn)
        + np.log(fg)
        + 2 * np.log(u)  # log|J| (the -log2 is constant)
    )


lam = np.linspace(-8, 8, 601)
tau = np.linspace(-0.97, 0.97, 601)
FG = sig(lam)[:, None]
T = tau[None, :]
T2 = (1 - T**2) / 4.0

cases = [
    ("Dir(1/2,1/2,1/2)", (0.5, 0.5, 0.5), 0.5 * np.log(FG) + np.log(1 - FG) - 0.5 * np.log(T2)),
    (
        "Dir(1/2,1/4,1/4)",
        (0.5, 0.25, 0.25),
        0.5 * np.log(FG) + 0.5 * np.log(1 - FG) - 0.75 * np.log(T2),
    ),
]
for name, a, claimed in cases:
    exact = lam_tau_form(a[0], a[1], a[2], FG, T)
    d = exact - claimed
    print(
        f"  {name}: claimed closed form matches the exact transform up to a constant? "
        f"spread of (exact-claimed) = {float(np.max(d) - np.min(d)):.2e}"
    )
print()

print("=" * 78)
print("T2  PROPERNESS / L-INVARIANCE of the pure reference (no likelihood, no messages, no prior)")
print("=" * 78)
print("A PROPER psi gives an L-invariant posterior median. An IMPROPER one walks with L.\n")


def median_fg(psi_fn, L, K=20001):
    lm = np.linspace(-L, L, K)
    p = np.exp(psi_fn(lm) - np.max(psi_fn(lm)))
    c = np.cumsum(p) / p.sum()
    return sig(lm)[np.searchsorted(c, 0.5)]


variants = {
    "shipped (Haldane both)      ": lambda lm: np.zeros_like(lm),
    "1-DOF Jeffreys  +.5(lg+la)  ": lambda lm: 0.5 * np.log(sig(lm)) + 0.5 * np.log(1 - sig(lm)),
    "gDNA FITTED-flat + .5*la    ": lambda lm: 0.5 * np.log(1 - sig(lm)),
}
Ls = [6, 8, 10, 14, 20, 30]
print(f"{'variant':<30}" + "".join(f"{'L=' + str(Lv):>8}" for Lv in Ls) + f"{'spread':>9}")
print("-" * 78)
for nm, fn in variants.items():
    v = [median_fg(fn, Lv) for Lv in Ls]
    print(f"{nm:<30}" + "".join(f"{x:>8.4f}" for x in v) + f"{max(v) - min(v):>9.4f}")
print()
print(
    "  NOTE the 3rd row: with a FLAT logP_g (which is what np.interp(left=logP[0]) makes the NPMLE"
)
print(
    "  BELOW its grid), adding the RNA reference is improper at the OTHER vertex: f_g -> 0 as L grows."
)
print(
    "  The measure is only proper when BOTH arms are bounded. Fitted-logP_g must have a REAL left tail."
)
print()

print("=" * 78)
print("T3  REPARAMETERIZATION INVARIANCE — grid in lambda vs grid in f_g")
print("=" * 78)
print("The same prior, integrated on two different coordinates, must give the same median.\n")


def med_lambda(a_g, a_r, L=20, K=200001):
    lm = np.linspace(-L, L, K)
    f = sig(lm)
    lp = (a_g - 1) * np.log(f) + (a_r - 1) * np.log(1 - f) + np.log(f) + np.log(1 - f)
    p = np.exp(lp - lp.max())
    return f[np.searchsorted(np.cumsum(p) / p.sum(), 0.5)]


def med_f(a_g, a_r, K=200001):
    f = np.linspace(1e-9, 1 - 1e-9, K)
    lp = (a_g - 1) * np.log(f) + (a_r - 1) * np.log(1 - f)
    p = np.exp(lp - lp.max())
    return f[np.searchsorted(np.cumsum(p) / p.sum(), 0.5)]


print(f"{'prior':<20}{'median via lambda-grid':>24}{'median via f-grid':>20}{'Beta median':>14}")
print("-" * 78)
for a_g, a_r in [(0.5, 0.5), (1.0, 1.0), (0.5, 1.0)]:
    print(
        f"{'Beta(' + str(a_g) + ',' + str(a_r) + ')':<20}{med_lambda(a_g, a_r):>24.5f}"
        f"{med_f(a_g, a_r):>20.5f}{B.median(a_g, a_r):>14.5f}"
    )
print(
    "\n  agreement => the lambda-grid + the written reference correctly represent the intended prior."
)
