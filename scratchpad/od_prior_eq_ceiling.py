"""Adversarial evaluation: does prior od0 == ceiling od hold up?

Pure math / no repo data needed beyond the estimator's formulas.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import betabinom


def od_of(a):
    return 1.0 / (2.0 * a + 1.0)


def a_of(od):
    return 0.5 * (1.0 - od) / od


print("=" * 78)
print("0. od <-> a")
for a in [2, 3, 4, 4.5, 5, 7, 14, 29]:
    print(f"   a={a:<5} od={od_of(a):.4f}")
print(f"   od=0.10 -> a={a_of(0.10):.3f};  od=C/2=0.10 (C=0.2) -> a={a_of(0.1):.3f}")

print()
print("=" * 78)
print("1. JOB 2 (rejection null): two-sided tail p of the vcap top seed N=1523 sense=5")
print("   Bonferroni thresholds: 0.05/160366 = %.2e (LBX0190), 0.05/863246 = %.2e (vcap)"
      % (0.05 / 160366, 0.05 / 863246))
N, k = 1523, 5
for a in [2, 3, 4, 4.5, 5, 6, 8, 14]:
    rv = betabinom(N, a, a)
    p = 2.0 * rv.cdf(k)
    verd_l = "REJECT" if p < 0.05 / 160366 else "keep  "
    verd_v = "REJECT" if p < 0.05 / 863246 else "keep  "
    print(f"   a={a:<5} od={od_of(a):.4f}  p={p:.3e}   LBX0190:{verd_l}  vcap:{verd_v}")

print()
print("   owner's examples (must be KEPT):")
for (n, s) in [(10, 9), (100, 90), (10, 10), (20, 18)]:
    row = []
    for a in [2, 3, 4, 4.5]:
        rv = betabinom(n, a, a)
        p = 2.0 * min(rv.cdf(s - 1) if s > 0 else 0.0, 1.0 - rv.cdf(s - 1))
        p = 2.0 * min(rv.cdf(s), 1.0 - rv.cdf(s - 1))
        row.append(f"a={a}: p={min(p,1.0):.3f}")
    print(f"   {s}/{n}: " + "  ".join(row))

print()
print("=" * 78)
print("2. n_eff = n / [1 + (n-1) od]  -- how much a seed is worth")
ods = [od_of(2), od_of(3), od_of(4.5), od_of(14), 0.02, 0.0]
print("   n      " + "  ".join(f"od={o:.4f}" for o in ods))
for n in [10, 100, 1000, 1523, 10000]:
    row = [n / (1 + (n - 1) * o) if o > 0 else float(n) for o in ods]
    print(f"   {n:<6} " + "  ".join(f"{v:9.1f}" for v in row))

print()
print("=" * 78)
print("3. STRAND FISHER INFORMATION for f_g -- 'destroying our strongest evidence'")
print("   large-N limit: I(g) = (kappa-1/2)^2 / (od * Q(g)),  Q = g^2/4 + (1-g)^2 k(1-k)")
print("   -> se(g_hat) = 1/sqrt(I), INDEPENDENT OF N (od caps the information).")


def I_mean_large_n(kappa, od, g):
    Q = g * g * 0.25 + (1 - g) ** 2 * kappa * (1 - kappa)
    return (kappa - 0.5) ** 2 / (od * Q)


for kappa in [0.99, 0.95, 0.90, 0.75]:
    print(f"   kappa={kappa}")
    for g in [0.25, 0.5, 0.75]:
        row = []
        for a, lbl in [(2, "a=2 (ceil)"), (4.5, "a=4.5 (C/2)"), (14, "a=14 (now)")]:
            I = I_mean_large_n(kappa, od_of(a), g)
            row.append(f"{lbl}: I={I:8.1f} se={1/np.sqrt(I):.3f}")
        print(f"      g={g}:  " + " | ".join(row))

print()
print("   information RATIO a=14 -> a=2 (fraction of strand info retained at the ceiling):")
for kappa in [0.99, 0.95, 0.9, 0.75]:
    r = od_of(14) / od_of(2)
    print(f"      kappa={kappa}: ratio = od(14)/od(2) = {r:.4f}  (independent of kappa & g)")
    break
print("      -> adopting the ceiling as the no-information default DISCARDS "
      f"{100*(1-od_of(14)/od_of(2)):.1f}% of the strand Fisher information.")
print(f"      -> a=14 -> a=4.5 (C/2) discards {100*(1-od_of(14)/od_of(4.5)):.1f}%.")

print()
print("=" * 78)
print("4. THE SPURIOUS VARIANCE CHANNEL: at kappa=1/2 an od>0 makes a BALANCED node INFORMATIVE")
print("   strand_loglik(g) with kappa=1/2, s = N/2 exactly:")
print("   var(g) = N/4 + od*N^2/4*(g^2+(1-g)^2); loglik = -0.5*(s-N/2)^2/var - 0.5*log var")


def dll(N, s, od, kappa=0.5):
    g = np.array([0.0, 0.5, 1.0])
    p = 0.5 * g + kappa * (1 - g)
    mean = N * p
    var = (N * p * (1 - p) + (N * g) ** 2 * 0.25 * od
           + (N * (1 - g)) ** 2 * kappa * (1 - kappa) * od)
    ll = -0.5 * (s - mean) ** 2 / var - 0.5 * np.log(var)
    return ll


print("   BALANCED node (s = N/2). Delta = loglik(g=0.5) - loglik(g=0), in nats:")
print("   N      od=0.2(ceil)  od=0.1(C/2)  od=0.0345(now)   od=0")
for N in [2, 4, 10, 20, 50, 100, 1000]:
    s = N / 2
    row = []
    for od in [0.2, 0.1, od_of(14), 0.0]:
        L = dll(N, s, od)
        row.append(L[1] - L[0])
    print(f"   {N:<6} " + "  ".join(f"{v:11.4f}" for v in row))
print("   (asymptote = 0.5*log2 = %.4f nats; ANY od>0 breaks the 'unstranded node is "
      "uninformative' invariant)" % (0.5 * np.log(2)))

print()
print("   DISPLACED balanced-mean node (s = 0.9N), Delta = loglik(g=0.5)-loglik(g=0):")
for N in [10, 20, 50, 100, 1000]:
    s = 0.9 * N
    row = []
    for od in [0.2, 0.1, od_of(14), 0.0]:
        L = dll(N, s, od)
        row.append(L[1] - L[0])
    print(f"   {N:<6} " + "  ".join(f"{v:11.4f}" for v in row))

print()
print("   the artifact SATURATES once od*N >~ 1, i.e. N >~ 1/od:")
for a in [2, 4.5, 14]:
    print(f"      a={a}: od={od_of(a):.4f} -> saturates above N ~ {1/od_of(a):.0f} fragments")

print()
print("=" * 78)
print("5. UPWARD BIAS of shrink-to-ceiling.  od_final = (I*od_mom + W*C)/(I+W), clipped <= C")
print("   bias = W/(I+W) * (C - od_mom) >= 0 ALWAYS.  Fraction of the ceiling reported when")
print("   the truth is od_mom = 0 (the synthetic suite's ground truth):")
print("   W/I     od_final/C")
for r in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 10.0]:
    print(f"   {r:<7} {r/(1+r):.4f}")

print()
print("=" * 78)
print("6. P2 UNITS FIX: what happens to the prior's weight share if W stays 30?")
rows = [("LBX0190", 160366, 699894), ("LBX0588", 297293, 4914517),
        ("MO_3021", 624913, 10745654), ("vcap", 863246, 100998511)]
print("   sample     nodes      pairs        W/(n+W) now     W/(pairs+W) after P2")
for nm, n, P in rows:
    print(f"   {nm:<10} {n:<10} {P:<12} {30/(n+30)*100:10.4f}%   {30/(P+30)*100:14.6f}%")
print("   -> the units fix ALONE makes the prior 4-100x WEAKER, not stronger. W must be re-derived.")

print()
print("   MoM sampling se: for seeds of n_c=2, per-seed ratio is exactly +/-1 (mean 0, var 1)")
print("   under od=0 => se(od_mom) = 1/sqrt(P_eff). Implied resolution on real samples:")
for nm, n, P in rows:
    print(f"   {nm:<10} P={P:<12} se(od_mom) ~ {1/np.sqrt(P):.5f}")
print("   -> at 0.7M pairs the data resolves od to +/-0.0012. LBX0190 is NOT underpowered for")
print("      od as a VARIANCE; if it saturates, that is BIAS (contamination), which no prior fixes.")

print()
print("   W implied by a proper prior's own precision (W = 1/sigma0^2, pair units):")
for lbl, sig in [("uniform on [0,C], sd=C/sqrt(12)", 0.2 / np.sqrt(12)),
                 ("half of C", 0.1), ("C itself", 0.2)]:
    print(f"      {lbl:<34} sigma0={sig:.4f} -> W={1/sig**2:8.1f} pairs")
print("   -> ANY defensible W is O(10^2-10^3) pairs, i.e. NEGLIGIBLE against 7e5-1e8 real pairs.")
print("      The prior only ever bites on a library with < ~1e3 informative pairs.")

print()
print("=" * 78)
print("7. COHERENCE: can a proper prior on [0,C] have mean C?")
print("   E_pi[od] <= C with equality iff pi = delta_C (degenerate).")
print("   Linear shrinkage toward m with FINITE weight W is the posterior mean of a conjugate")
print("   prior with mean m. m = C therefore requires pi = delta_C => W = infinity, not W = 30.")
print("   With finite W, m=C is the posterior mean of a prior with support ABOVE C.")
print("   The interior ignorance target under a flat prior on the ICC is E = C/2 = %.3f (a = %.1f)."
      % (0.2 / 2, a_of(0.1)))
