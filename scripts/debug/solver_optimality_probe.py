"""Probe the per-node simplex solver in ISOLATION (no messages, no global) vs an EXACT-likelihood reference.

The production strand term (`_mixture_strand_loglik`) is a GAUSSIAN approximation to a Binomial/Beta-Binomial
mixture. Question: on its own evidence, does the solver recover f_g≈0 for a pure-RNA node (zero-gDNA), or does
the Gaussian approx / lattice / median / overdispersion inject PHANTOM gDNA?

Reference = the EXACT convolved Beta-Binomial likelihood of u_pos = (gDNA tagged +) + (+RNA tagged +), with
N_g = N·f_g gDNA fragments at rate ½ (od_g) and N_r = N·(1−f_g) +RNA at rate κ (od_r), over a fine f_g grid
(single-strand: f_neg=0). Same Jeffreys prior + median readout, so only the LIKELIHOOD MODEL differs.

    python scripts/debug/solver_optimality_probe.py
"""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp
from scipy.stats import betabinom, binom

from rigel.calibration.simplex_sweep import _fg_median, _local_loglik, _simplex_lattice

_Z = np.array([0.0])
_T = np.array([True])
_F = np.array([False])


def solver_fg(u_pos, N, kappa, od_g, od_r, *, jeffreys=True, n_grid=60, readout="median"):
    """Production per-node solve: Gaussian strand term + (Jeffreys) over the n_grid lattice, no msgs/global."""
    lat = _simplex_lattice(n_grid)
    up = np.array([float(u_pos)])
    un = np.array([float(N - u_pos)])
    psi = _local_loglik(up, un, _Z, _Z, _T, _F, kappa, od_g, od_r, lat,
                        strand_obs=_T if jeffreys else None)
    if readout == "mean":
        post = np.exp(psi - logsumexp(psi, axis=1, keepdims=True))
        return float((post @ lat[2])[0])
    return float(_fg_median(psi, lat[2])[0])


def _bb_pmf(n, rate, od):
    if n <= 0:
        return np.array([1.0])
    if od <= 1e-9:
        return binom.pmf(np.arange(n + 1), n, rate)
    s = 1.0 / od - 1.0  # a+b ; intra-class corr = 1/(a+b+1) = od
    return betabinom.pmf(np.arange(n + 1), n, rate * s, (1.0 - rate) * s)


def exact_fg(u_pos, N, kappa, od_g, od_r, *, jeffreys=True, fine=800, readout="median"):
    """Exact convolved Beta-Binomial posterior over f_g (single-strand, f_neg=0)."""
    fgs = np.linspace(0.0, 1.0, fine + 1)
    k = int(round(u_pos))
    ll = np.empty(fgs.size)
    for i, fg in enumerate(fgs):
        n_g = int(round(N * fg))
        conv = np.convolve(_bb_pmf(n_g, 0.5, od_g), _bb_pmf(N - n_g, kappa, od_r))
        ll[i] = np.log(max(conv[k] if k < conv.size else 0.0, 1e-300))
    if jeffreys:
        ll = ll + (-0.5) * (np.log(np.clip(fgs, 1e-3, 1)) + np.log(np.clip(1 - fgs, 1e-3, 1)))
    post = np.exp(ll - logsumexp(ll))
    if readout == "mean":
        return float(post @ fgs)
    return float(fgs[np.searchsorted(np.cumsum(post), 0.5)])


def solver_fg_2d(u_pos, N, kappa, od_g, od_r, *, n_grid=60):
    """Production 2-D AMBIG solve: both strands allowed, NO Jeffreys (strand_obs=False), median f_g."""
    lat = _simplex_lattice(n_grid)
    up = np.array([float(u_pos)])
    un = np.array([float(N - u_pos)])
    psi = _local_loglik(up, un, _Z, _Z, _T, _T, kappa, od_g, od_r, lat, strand_obs=None)
    return float(_fg_median(psi, lat[2])[0])


def exact_fg_2d(u_pos, N, kappa, od_g, od_r, *, fine=60):
    """Exact 3-component (gDNA½ / +RNAκ / −RNA(1−κ)) convolved Beta-Binomial posterior, f_g median.
    Flat prior over the 2-simplex (AMBIG gets no Jeffreys), to match the production AMBIG solve."""
    k = int(round(u_pos))
    fg_vals, weights = [], []
    for i in range(fine + 1):
        for j in range(fine + 1 - i):
            fpos, fneg = i / fine, j / fine
            n_pos = int(round(N * fpos))
            n_neg = int(round(N * fneg))
            n_g = N - n_pos - n_neg
            if n_g < 0:
                continue
            conv = np.convolve(np.convolve(_bb_pmf(n_g, 0.5, od_g), _bb_pmf(n_pos, kappa, od_r)),
                               _bb_pmf(n_neg, 1.0 - kappa, od_r))
            fg_vals.append(1.0 - fpos - fneg)
            weights.append(conv[k] if k < conv.size else 0.0)
    fg_vals = np.array(fg_vals)
    w = np.array(weights)
    w = w / max(w.sum(), 1e-300)
    order = np.argsort(fg_vals)
    cw = np.cumsum(w[order])
    return float(fg_vals[order][np.searchsorted(cw, 0.5)])


def main():
    KAP = 0.0094  # the test's fitted κ
    OG, OR = 0.20, 0.101  # the test's fitted overdispersions (od_g = fallback MAX)

    print("=== ZERO-gDNA per-node solve (f_g_true=0, single-strand +RNA): does f_g stay ~0? ===")
    print(f"  κ={KAP}, od_g={OG}, od_r={OR}   (u_pos = round(N·κ); truth f_g=0)")
    print(f"  {'N':>6} {'u_pos':>6} | {'solver+J':>9} {'solver-J':>9} {'solverMEAN':>10} | "
          f"{'EXACT+J':>8} {'EXACTmean':>9} | {'binom+J':>8}")
    for N in (20, 50, 100, 191, 627, 2000, 8000):
        up = round(N * KAP)
        s_j = solver_fg(up, N, KAP, OG, OR, jeffreys=True)
        s_nj = solver_fg(up, N, KAP, OG, OR, jeffreys=False)
        s_mean = solver_fg(up, N, KAP, OG, OR, jeffreys=True, readout="mean")
        e_j = exact_fg(up, N, KAP, OG, OR, jeffreys=True)
        e_mean = exact_fg(up, N, KAP, OG, OR, jeffreys=True, readout="mean")
        b_j = exact_fg(up, N, KAP, 0.0, 0.0, jeffreys=True)  # binomial (no overdispersion)
        print(f"  {N:>6} {up:>6} | {s_j:>9.3f} {s_nj:>9.3f} {s_mean:>10.3f} | "
              f"{e_j:>8.3f} {e_mean:>9.3f} | {b_j:>8.3f}")

    print("\n=== OVERDISPERSION sweep at the boundary (N=191, u_pos=2, κ=0.0094, f_g_true=0) ===")
    print(f"  {'od_g':>6} | {'solver+J':>9} {'EXACT+J':>8} {'binom+J':>8}")
    for og in (0.0, 0.02, 0.05, 0.10, 0.20):
        s = solver_fg(2, 191, KAP, og, OR, jeffreys=True)
        e = exact_fg(2, 191, KAP, og, OR, jeffreys=True)
        b = exact_fg(2, 191, KAP, 0.0, 0.0, jeffreys=True)
        print(f"  {og:>6.2f} | {s:>9.3f} {e:>8.3f} {b:>8.3f}")

    print("\n=== RECOVERY: does the solver read gDNA when it IS present? (N=500, κ=0.0094) ===")
    print(f"  {'f_g_true':>9} {'u_pos':>6} | {'solver+J':>9} {'EXACT+J':>8}")
    for fg in (0.0, 0.1, 0.3, 0.5, 0.8, 1.0):
        p = 0.5 * fg + KAP * (1 - fg)
        up = round(500 * p)
        s = solver_fg(up, 500, KAP, OG, OR, jeffreys=True)
        e = exact_fg(up, 500, KAP, OG, OR, jeffreys=True)
        print(f"  {fg:>9.2f} {up:>6} | {s:>9.3f} {e:>8.3f}")

    print("\n=== 2-D AMBIG node (balanced +/− RNA, f_g_true=0): the 2-D lattice vs exact ===")
    print(f"  κ={KAP}, od_g={OG}, od_r={OR}; balanced u_pos≈u_neg → f_g_true=0 (pure RNA, half + half −)")
    print(f"  {'N':>6} {'u_pos':>6} | {'solver2D':>9} {'EXACT2D':>8}")
    for N in (200, 600, 1220):
        # balanced RNA, f_g=0: + and − each N/2; observed u_pos = N_+·(read-+ rate) + N_-·(read-+ rate)
        p = KAP * 0.5 + (1.0 - KAP) * 0.5  # = 0.5 (balanced RNA looks like p≈0.5)
        up = round(N * p)
        s = solver_fg_2d(up, N, KAP, OG, OR)
        e = exact_fg_2d(up, N, KAP, OG, OR, fine=60)
        print(f"  {N:>6} {up:>6} | {s:>9.3f} {e:>8.3f}")

    print("\n=== UNSTRANDED (κ=0.5): per-node solve is uninformative — what does each give? ===")
    print(f"  {'N':>6} {'u_pos':>6} {'f_g_true':>9} | {'solver+J':>9} {'EXACT+J':>8}")
    for fg in (0.0, 1.0):
        for N in (100, 1000):
            up = round(N * 0.5)
            s = solver_fg(up, N, 0.5, OG, OR, jeffreys=True)
            e = exact_fg(up, N, 0.5, OG, OR, jeffreys=True)
            print(f"  {N:>6} {up:>6} {fg:>9.2f} | {s:>9.3f} {e:>8.3f}")


if __name__ == "__main__":
    main()
