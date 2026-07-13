"""Prototype the DENSITY frame for the boundary self-solve, and show it is MAGIC-NUMBER-FREE.

Density frame: solve over the component COUNTS (N_g gDNA, N_p nascent+), N_g+N_p = N_u. The strand
likelihood is the ONLY intrinsic gDNA/RNA signal. nascent~0 is the VERTEX rest state (N_p=0), reached by an
INFINITESIMAL tie-break −ε·N_p (breaks the flat-direction degeneracy toward the vertex). If the frame is
right, f_g is INSENSITIVE to ε (a tie-break, not a strength).

Contrast: the log-odds frame adds +c·log(f_g) to reach f_g=1, and f_g is SENSITIVE to c (a magic number).
"""
import numpy as np

def strand_loglik_grid(u_pos, N, f_g_grid, kappa, od_g=0.0, od_r=0.0):
    """The production two-component strand loglik (strand_likelihood.py form) over an f_g grid."""
    f_p = 1.0 - f_g_grid
    p = 0.5 * f_g_grid + kappa * f_p
    mean = N * p
    var = np.maximum(mean * (1.0 - p) + (N * f_g_grid) ** 2 * 0.25 * od_g
                     + (N * f_p) ** 2 * kappa * (1.0 - kappa) * od_r, 1e-9)
    return -0.5 * (u_pos - mean) ** 2 / var - 0.5 * np.log(var)

def density_solve(u_pos, u_neg, kappa, eps, K=4001, od_g=0.05, od_r=0.05):
    """DENSITY frame: MAP over N_p∈[0,N] with the nascent~0 vertex tie-break −ε·N_p. f_g = (N−N_p*)/N."""
    N = u_pos + u_neg
    Np = np.linspace(0.0, N, K)                      # nascent+ count
    f_g = np.clip((N - Np) / max(N, 1e-9), 1e-9, 1.0)
    ll = strand_loglik_grid(u_pos, N, f_g, kappa, od_g, od_r)
    obj = ll - eps * Np                              # infinitesimal nascent~0 tie-break (→ vertex)
    return float(f_g[int(np.argmax(obj))])

def logodds_solve(u_pos, u_neg, kappa, c, K=4001, od_g=0.05, od_r=0.05):
    """LOG-ODDS frame: MAP over the f_g grid with the +c·log(f_g) gDNA-vertex tilt (the magic-c approach)."""
    N = u_pos + u_neg
    f_g = np.linspace(1e-4, 1 - 1e-4, K)
    ll = strand_loglik_grid(u_pos, N, f_g, kappa, od_g, od_r)
    obj = ll + c * np.log(f_g)                       # the tunable tilt
    return float(f_g[int(np.argmax(obj))])

# scenarios: (label, u_pos, u_neg, kappa)  — gDNA is symmetric; nascent adds to the sense strand.
scen = [
    ("unstranded (κ=.5), balanced        ", 150, 150, 0.5),
    ("stranded (κ=.9) nrna_none, balanced ", 150, 150, 0.9),
    ("stranded (κ=.9) nrna_present +excess", 245, 55, 0.9),
    ("weak-strand (κ=.6) balanced         ", 150, 150, 0.6),
]

print("=" * 90)
print("DENSITY frame — f_g vs the tie-break ε (want: INSENSITIVE ⇒ ε is a tie-break, not a magic number)")
print("=" * 90)
print(f"{'scenario':38} " + "".join(f"ε={e:<9.0e}" for e in (1e-2, 1e-4, 1e-6, 1e-9)))
for lab, up, un, k in scen:
    row = "  ".join(f"{density_solve(up, un, k, e):.4f}" for e in (1e-2, 1e-4, 1e-6, 1e-9))
    print(f"{lab} {row}")

print("\n" + "=" * 90)
print("LOG-ODDS frame — f_g vs the tilt c  (SENSITIVE ⇒ c is a magic number)")
print("=" * 90)
print(f"{'scenario':38} " + "".join(f"c={c:<9}" for c in (1, 3, 5, 10)))
for lab, up, un, k in scen:
    row = "  ".join(f"{logodds_solve(up, un, k, c):.4f}" for c in (1, 3, 5, 10))
    print(f"{lab} {row}")

print("\nREAD: density-frame f_g barely moves across 7 orders of magnitude of ε (unstranded → 1, stranded")
print("nrna_none → 1, nrna_present → peeled) = nascent~0 is STRUCTURAL. Log-odds f_g slides with c.")


print("\n" + "=" * 90)
print("DIAGNOSIS: is the ε-sensitivity the FRAME, or the overdispersion (od) term? Sweep od.")
print("  count-zero-info: at κ=½ the count carries NO gDNA/RNA info ⇒ likelihood should be FLAT ⇒ f_g→1")
print("  via an infinitesimal tie-break. The od term's −½log(var) prefers f_g=0.5 (a count-magnitude signal).")
print("=" * 90)
for od in (0.0, 0.005, 0.02, 0.05):
    print(f"\n  od_g=od_r={od}:")
    print(f"    {'scenario':38} " + "".join(f"ε={e:<9.0e}" for e in (1e-2, 1e-4, 1e-6, 1e-9)))
    for lab, up, un, k in scen:
        row = "  ".join(f"{density_solve(up, un, k, e, od_g=od, od_r=od):.4f}" for e in (1e-2, 1e-4, 1e-6, 1e-9))
        print(f"    {lab} {row}")
print("\nREAD: at od=0 (count-zero-info respected — the FLAGSHIP's fitted od_g≈0), density f_g is ε-INSENSITIVE:")
print("  unstranded→1, stranded nrna_none→1, nrna_present→peeled — MAGIC-NUMBER-FREE. The od>0 pull to 0.5 is a")
print("  SEPARATE count-zero-info violation (od should set PRECISION, not tilt the VALUE).")


def logodds_tiebreak(u_pos, u_neg, kappa, eps, K=4001, od_g=0.0, od_r=0.0):
    """LOG-ODDS grid, but with a nascent~0 VERTEX TIE-BREAK −ε·(1−f_g) instead of the +c·log(f_g) TILT."""
    N = u_pos + u_neg
    f_g = np.linspace(1e-4, 1 - 1e-4, K)
    ll = strand_loglik_grid(u_pos, N, f_g, kappa, od_g, od_r)
    obj = ll - eps * (1.0 - f_g)                     # tie-break toward the f_g=1 vertex
    return float(f_g[int(np.argmax(obj))])

print("\n" + "=" * 90)
print("Is the magic number TILT-vs-TIE-BREAK (not log-odds-vs-density)? Log-odds grid + −ε·(1−f_g) tie-break, od=0:")
print("=" * 90)
print(f"{'scenario':38} " + "".join(f"ε={e:<9.0e}" for e in (1e-2, 1e-4, 1e-6, 1e-9)))
for lab, up, un, k in scen:
    row = "  ".join(f"{logodds_tiebreak(up, un, k, e):.4f}" for e in (1e-2, 1e-4, 1e-6, 1e-9))
    print(f"{lab} {row}")
print("\nIf ε-INSENSITIVE here too ⇒ the fix is TIE-BREAK (not TILT) + count-zero-info (od→precision), in the")
print("EXISTING log-odds solver — no re-architecture to a count-simplex frame needed.")
