"""Pin the nascent~0 value mechanism: can a gDNA-vertex prior term in the grid ψ drive an UNSTRANDED
node's f_g -> 1 (nascent~0), while a STRANDED node still peels via the strand likelihood?
Uses the real _solve_nodes_logodds_all with a custom global_logprior = c*log(f_g)  (favors f_g=1)."""
import numpy as np
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all, _logodds_grid

L, K = 10.0, 256
_, fg = _logodds_grid(K, L)            # the f_g grid the global_logprior is evaluated on
log_fg = np.log(np.clip(fg, 1e-9, 1.0))

def solve(u_pos, u_neg, kappa, c, allow=(True, False)):
    m = 1
    up = np.array([float(u_pos)]); un = np.array([float(u_neg)])
    ap = np.array([allow[0]]); an = np.array([allow[1]])
    so = np.array([allow[0] ^ allow[1]])
    mu = np.array([float(u_pos + u_neg)]); ms = np.zeros(1)
    glp = (c * log_fg)[None, :] if c else None    # gDNA-vertex prior: +c·log f_g (favors f_g→1)
    dc = _solve_nodes_logodds_all(up, un, ap, an, so, mu, ms, kappa=kappa, od_g=0.0, od_r=0.0,
                                  n_grid=K, L=L, n_grid_ss=K, global_logprior=glp)
    return float(dc.gdna_frac[0])

print("nascent~0 prior = +c·log(f_g)  (a gDNA-vertex tilt; c = pseudo-fragments of gDNA)")
print(f"{'c':>4} | unstranded(κ=.5, 150/150)  weak-strand(κ=.6, 165/135)  stranded(κ=.9, 250/50)")
for c in (0, 1, 2, 3, 5, 10):
    f_unstr = solve(150, 150, 0.5, c)
    f_weak  = solve(165, 135, 0.6, c)
    f_str   = solve(250, 50, 0.9, c)
    print(f"{c:>4} |   f_g={f_unstr:.3f}              f_g={f_weak:.3f}            f_g={f_str:.3f}")
print("\nWANT: unstranded -> ~1 (nascent~0);  stranded -> peels DOWN (real nascent on the + excess).")
print("(strand-biased 250/50 at κ=.9 ⇒ real RNA present ⇒ f_g should drop below 1.)")
