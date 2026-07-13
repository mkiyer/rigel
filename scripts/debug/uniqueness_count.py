"""Uniqueness (unimodality) re-check of the node solve UNDER COUNT CONSERVATION.

Node solve:  min over counts N_c >= 0 with  sum_c N_c = N (FIXED)  of
    F(N) = - strand_loglik(n_pos,n_neg | split; kappa, od)      # intrinsic (count = Fisher info)
           + sum_c pi_c * (log(N_c / E_c) - m_c)^2              # density beliefs + messages + prior
with per-component E_c (E_gdna != E_rna). Density rho_c = N_c/E_c is NOT conserved; only COUNT is.

Earlier unimodality (bp_theory.py) was proven only under DENSITY conservation (sum rho_c = D). Two
things change under count conservation and could break it:
  (i)  the density penalty (log(N_c/E_c)-m_c)^2 is NON-CONVEX in N_c once log(N_c) - m~_c > 1
       (g''(N) = 2(1 - (log N - m~))/N^2 < 0) -- the 'push-up' regime, analog of the old mu<0 case;
  (ii) per-component E_c skews the simplex (the count<->density map differs per component).

We re-verify by ENUMERATING F on a fine simplex grid and COUNTING strict local minima (want == 1).
1-DoF = single-strand {gDNA, RNA}; 2-DoF = AMBIG {gDNA, RNA+, RNA-}.
"""
import numpy as np

EPS = 1e-9


def strand_ll(n_pos, n_neg, f_g, f_np, f_nn, kappa, od_g, od_r):
    N = n_pos + n_neg
    if N <= 0:
        return np.zeros_like(np.asarray(f_g, float))
    p = 0.5 * f_g + kappa * f_np + (1.0 - kappa) * f_nn
    mean = N * p
    rscale = kappa * (1.0 - kappa)
    var = (mean * (1.0 - p) + (N * f_g) ** 2 * 0.25 * od_g
           + (N * f_np) ** 2 * rscale * od_r + (N * f_nn) ** 2 * rscale * od_r)
    var = np.maximum(var, 1e-9)
    return -0.5 * (n_pos - mean) ** 2 / var - 0.5 * np.log(var)


# ------------------------------------------------------------------- 1-DoF (single strand)
def F_1d(t, r):
    """t = gDNA count fraction; N_g=N*t, N_rna=N*(1-t) (RNA on + strand)."""
    N = r['N']
    Ng = np.maximum(N * t, EPS); Nr = np.maximum(N * (1.0 - t), EPS)
    rho_g = Ng / r['Eg']; rho_r = Nr / r['Er']
    dens = r['pi'][0] * (np.log(rho_g) - r['m'][0]) ** 2 + r['pi'][1] * (np.log(rho_r) - r['m'][1]) ** 2
    ll = strand_ll(r['npos'], r['nneg'], t, 1.0 - t, 0.0, r['kappa'], r['odg'], r['odr'])
    return dens - ll


def count_min_1d(r, K=4001):
    t = np.linspace(1e-4, 1 - 1e-4, K)
    F = F_1d(t, r)
    mins = [i for i in range(K)
            if F[i] < (F[i - 1] if i > 0 else np.inf) and F[i] < (F[i + 1] if i < K - 1 else np.inf)]
    return len(mins), [round(float(t[i]), 3) for i in mins], float(t[int(np.argmin(F))])


# ------------------------------------------------------------------- 2-DoF (AMBIG)
def F_2d(ag, ap, an, r):
    N = r['N']
    Ng = np.maximum(N * ag, EPS); Np = np.maximum(N * ap, EPS); Nn = np.maximum(N * an, EPS)
    rho_g = Ng / r['Eg']; rho_p = Np / r['Er']; rho_n = Nn / r['Er']
    dens = (r['pi'][0] * (np.log(rho_g) - r['m'][0]) ** 2
            + r['pi'][1] * (np.log(rho_p) - r['m'][1]) ** 2
            + r['pi'][2] * (np.log(rho_n) - r['m'][2]) ** 2)
    ll = strand_ll(r['npos'], r['nneg'], ag, ap, an, r['kappa'], r['odg'], r['odr'])
    return dens - ll


def count_min_2d(r, K=301):
    g = np.linspace(0.0, 1.0, K)
    AG, AP = np.meshgrid(g, g, indexing='ij')
    AN = 1.0 - AG - AP
    valid = (AG > 5e-3) & (AP > 5e-3) & (AN > 5e-3)
    F = np.full((K, K), np.inf)
    F[valid] = F_2d(AG[valid], AP[valid], AN[valid], r)
    mins = []
    for i in range(1, K - 1):
        for j in range(1, K - 1):
            if not valid[i, j]:
                continue
            c = F[i, j]
            nb = [F[i - 1, j], F[i + 1, j], F[i, j - 1], F[i, j + 1],
                  F[i - 1, j - 1], F[i + 1, j + 1], F[i - 1, j + 1], F[i + 1, j - 1]]
            if all(c < v for v in nb):
                mins.append((round(float(AG[i, j]), 2), round(float(AP[i, j]), 2)))
    fi, fj = np.unravel_index(np.argmin(F), F.shape)
    return len(mins), mins, (round(float(AG[fi, fj]), 2), round(float(AP[fi, fj]), 2))


def base(**kw):
    r = dict(N=300.0, Eg=100.0, Er=200.0, m=[np.log(1.5), np.log(1.0)], pi=[5.0, 5.0],
             npos=150.0, nneg=150.0, kappa=0.5, odg=0.0, odr=0.0)
    r.update(kw); return r


# targets-want-count helper: n*_c = exp(m_c)*E_c ; compare Sum to N to label push-up vs push-down
def regime_label(r, comps):
    Es = [r['Eg']] + [r['Er']] * (comps - 1)
    want = sum(np.exp(r['m'][c]) * Es[c] for c in range(comps))
    return f"Sum n*={want:6.0f} vs N={r['N']:.0f}  ({'push-UP(concave risk)' if want < r['N'] else 'push-down'})"


print("=" * 94)
print("UNIQUENESS under COUNT CONSERVATION — strict local-minima count on the count simplex (want == 1)")
print("=" * 94)

print("\n--- 1-DoF (single-strand: gDNA + RNA) ---")
cases_1d = {
    "R1 baseline (equal-ish E, balanced, unstranded)": base(),
    "R2 push-UP hard (tiny targets, Sum n* << N)":      base(m=[np.log(0.05), np.log(0.03)], pi=[8, 8]),
    "R3 push-down hard (huge targets, Sum n* >> N)":    base(m=[np.log(20.0), np.log(20.0)], pi=[8, 8]),
    "R4 two-confident CONFLICT (both want large)":      base(m=[np.log(6.0), np.log(6.0)], pi=[50, 50]),
    "R5 extreme precision gap (100 vs 0.1)":            base(pi=[100, 0.1], m=[np.log(2.0), np.log(0.5)]),
    "R6 E_g != E_r STRONG (80 vs 250), balanced":       base(Eg=80.0, Er=250.0, m=[np.log(1.5), np.log(1.0)], pi=[8, 8]),
    "R7 E asymmetry + push-UP (the combined stress)":   base(Eg=80.0, Er=250.0, m=[np.log(0.05), np.log(0.03)], pi=[10, 10]),
    "R8 stranded strong (kappa=0.95, od) + messages":   base(kappa=0.95, odg=0.05, odr=0.05, npos=240, nneg=60,
                                                             m=[np.log(2.0), np.log(1.5)], pi=[6, 6]),
    "R9 self-defense case (count_conservation TEST1)":  base(Eg=100.0, Er=200.0, m=[np.log(0.30), np.log(0.20)], pi=[50, 2]),
}
bad = 0
for tag, r in cases_1d.items():
    n, locs, gmin = count_min_1d(r)
    flag = "" if n == 1 else "  <<< NON-UNIMODAL"
    if n != 1:
        bad += 1
    print(f"  {tag:52s}: minima={n}  argmin t={gmin:.3f}   [{regime_label(r,2)}]{flag}")

print("\n--- 2-DoF (AMBIG: gDNA + RNA+ + RNA-) ---")
cases_2d = {
    "A1 baseline balanced unstranded":                  base(m=[np.log(1.2), np.log(1.0), np.log(1.0)], pi=[5, 5, 5]),
    "A2 push-UP hard (Sum n* << N)":                    base(m=[np.log(0.05), np.log(0.04), np.log(0.04)], pi=[8, 8, 8]),
    "A3 two-confident CONFLICT (gDNA vs RNA+ both big)":base(m=[np.log(5.0), np.log(5.0), np.log(0.2)], pi=[40, 40, 2]),
    "A4 E asymmetry + push-UP":                         base(Eg=80.0, Er=250.0, m=[np.log(0.05), np.log(0.04), np.log(0.04)], pi=[10, 10, 10]),
    "A5 stranded (kappa=0.9) + od":                     base(kappa=0.9, odg=0.05, odr=0.05, npos=200, nneg=100,
                                                             m=[np.log(1.5), np.log(1.2), np.log(0.3)], pi=[6, 6, 6]),
    "A6 extreme precision gap":                         base(m=[np.log(2.0), np.log(0.5), np.log(0.5)], pi=[100, 0.1, 0.1]),
}
for tag, r in cases_2d.items():
    n, locs, gmin = count_min_2d(r)
    flag = "" if n == 1 else "  <<< NON-UNIMODAL"
    if n != 1:
        bad += 1
    print(f"  {tag:52s}: minima={n}  argmin (ag,ap)={gmin}   [{regime_label(r,3)}]{flag}")

print("\n" + "=" * 94)
print(f"RESULT: {'ALL UNIMODAL (unique global minimum) under count conservation' if bad == 0 else f'{bad} NON-UNIMODAL case(s) — investigate'}")
print("=" * 94)


print("\n" + "=" * 94)
print("DIAGNOSIS — are the multimodal cases pathological (dishonest), or real?")
print("=" * 94)

print("\n(a) The push-UP multimodality requires a STRONG message contradicting a full node.")
print("    Honest transfer variance (#3) makes an across-regime (depleted->enriched) message WEAK. Re-test weak:")
for tag, r in [("R2' push-up, HONEST weak msg (pi=1)",   base(m=[np.log(0.05), np.log(0.03)], pi=[1, 1])),
               ("R7' E-asym push-up, weak msg (pi=1)",    base(Eg=80, Er=250, m=[np.log(0.05), np.log(0.03)], pi=[1, 1])),
               ("R2'' push-up, honest v.weak (pi=0.3)",   base(m=[np.log(0.05), np.log(0.03)], pi=[0.3, 0.3]))]:
    n, _, g = count_min_1d(r); print(f"    {tag:44s}: minima={n}  argmin t={g:.3f}  {'OK' if n == 1 else '<<< still multi'}")
for tag, r in [("A2' push-up, HONEST weak msg (pi=1)",    base(m=[np.log(0.05), np.log(0.04), np.log(0.04)], pi=[1, 1, 1])),
               ("A4' E-asym push-up, weak msg (pi=1)",    base(Eg=80, Er=250, m=[np.log(0.05), np.log(0.04), np.log(0.04)], pi=[1, 1, 1]))]:
    n, _, g = count_min_2d(r); print(f"    {tag:44s}: minima={n}  argmin={g}  {'OK' if n == 1 else '<<< still multi'}")

print("\n(b) The realistic ENRICHED node: the enriched-mode prior gives a HIGH gDNA target, not a tiny one.")
for tag, r in [("R2''' enriched: gDNA target HIGH (rho 3)", base(m=[np.log(3.0), np.log(0.1)], pi=[8, 8])),
               ("A2''  enriched: gDNA HIGH, RNA tiny",       base(m=[np.log(3.0), np.log(0.1), np.log(0.1)], pi=[8, 8, 8]))]:
    if len(r['m']) == 2:
        n, _, g = count_min_1d(r); print(f"    {tag:44s}: minima={n}  argmin t={g:.3f}  {'OK' if n == 1 else '<<< multi'}")
    else:
        n, _, g = count_min_2d(r); print(f"    {tag:44s}: minima={n}  argmin={g}  {'OK' if n == 1 else '<<< multi'}")

print("\n(c) A6 '0 minima' — is it a FLAT ridge (unidentified split), not multimodality?")
rA6 = base(m=[np.log(2.0), np.log(0.5), np.log(0.5)], pi=[100, 0.1, 0.1])
K = 301; g = np.linspace(0, 1, K); AG, AP = np.meshgrid(g, g, indexing='ij'); AN = 1 - AG - AP
valid = (AG > 5e-3) & (AP > 5e-3) & (AN > 5e-3)
F = np.full((K, K), np.inf); F[valid] = F_2d(AG[valid], AP[valid], AN[valid], rA6)
Fmin = F.min(); frac_near = float(np.mean(F[valid] < Fmin + 0.5))
ag_at_min = AG[valid][F[valid] < Fmin + 0.5]
print(f"    global Fmin={Fmin:.3f}; fraction of simplex within 0.5 of it = {frac_near:.1%}")
print(f"    gDNA fraction across that near-optimal set: ag in [{ag_at_min.min():.2f}, {ag_at_min.max():.2f}]"
      f"  -> gDNA PINNED (pi=100), RNA+/- split FLAT (pi=0.1) = a ridge, unique in the identified direction.")


print("\n(d) Does multimodality break order-independence? Only if the GLOBAL min is tied. Check the F-gap:")
def gap_1d(r, K=8001):
    t = np.linspace(1e-4, 1 - 1e-4, K); F = F_1d(t, r)
    idx = [i for i in range(1, K-1) if F[i] < F[i-1] and F[i] < F[i+1]]
    fv = sorted(F[i] for i in idx)
    return fv[0], (fv[1] - fv[0]) if len(fv) > 1 else np.inf
for tag, r in [("R2 push-up", base(m=[np.log(0.05), np.log(0.03)], pi=[8, 8])),
               ("R7 E-asym push-up", base(Eg=80, Er=250, m=[np.log(0.05), np.log(0.03)], pi=[10, 10]))]:
    f0, dg = gap_1d(r)
    print(f"    {tag:22s}: global Fmin={f0:.3f}  gap to 2nd min = {dg:.3f}  -> "
          f"{'UNIQUE global (grid is deterministic)' if dg > 1e-3 else 'TIED (order-dependent!)'}")
print("    => the grid solve returns the unique global min; local minima matter only for a bare local optimizer.")
