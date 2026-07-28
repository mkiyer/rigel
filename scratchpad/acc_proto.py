"""Accumulator redesign — numerical prototype.

Settles three questions that cannot be settled by argument:

  A. IDENTIFIABILITY. Can the cell counts n(A,Z) alone separate two components with
     DIFFERENT fragment-length distributions (the owner's cfRNA case: gDNA 60 bp,
     RNA 200 bp), with no strand information and no prior? How does the power decay
     as the two FLs converge?

  B. TERMINUS BIAS. Does the naive anchor rule (E_c(A) = |A| for every component)
     over-call gDNA near a transcript 3' end, and does the admissible-start
     correction remove it?

  C. HEAD-TO-HEAD. Cell estimator vs the CURRENT fractional-mass rule, on identical
     simulated data, as a function of region length.

Everything is a direct simulation of the generative model + an exact analytic design
matrix; no calibration code is involved.
"""

from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(20260728)


# ---------------------------------------------------------------------------
# fragment-length pmfs
# ---------------------------------------------------------------------------
def fl_pmf(mu: float, sd: float, n: int = 1200) -> np.ndarray:
    """Discretised positive Gaussian FL pmf, index = length in bp."""
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return p / p.sum()


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
class Partition:
    """A 1-D region partition on [0, end)."""

    def __init__(self, cuts: list[int]):
        self.cuts = np.asarray(sorted(set(cuts)), dtype=np.int64)
        self.starts = self.cuts[:-1]
        self.ends = self.cuts[1:]
        self.n = self.starts.size
        self.lens = self.ends - self.starts

    def region_of(self, pos):
        return np.searchsorted(self.cuts, pos, side="right") - 1


def E_gdna_cells(part: Partition, A: int, pmf: np.ndarray) -> np.ndarray:
    """E_g(A, Z) for every Z: expected count per unit density, exactly.

    E(A,Z) = sum_w f(w) * #{s in A : s + w - 1 in Z}.  Direct enumeration over the
    start positions of A (small in the prototype) — this is the ground truth the
    closed form must match.
    """
    out = np.zeros(part.n, dtype=np.float64)
    s = np.arange(part.starts[A], part.ends[A], dtype=np.int64)
    genome_end = int(part.cuts[-1])
    for w in np.nonzero(pmf)[0]:
        last = s + w - 1
        ok = last < genome_end
        if not ok.any():
            continue
        z = part.region_of(last[ok])
        out += pmf[w] * np.bincount(z, minlength=part.n)
    return out


class Transcript:
    """A spliced transcript: exons in genomic order. Mature fragments live in
    transcript coordinates and must FIT (u + w <= L_tx)."""

    def __init__(self, exons: list[tuple[int, int]]):
        self.exons = exons
        self.g2t = np.concatenate([np.arange(a, b) for a, b in exons])  # tcoord -> gpos
        self.L = self.g2t.size


def E_mature_cells(part: Partition, tx: Transcript, A: int, pmf: np.ndarray) -> np.ndarray:
    """E_r(A, Z) for mature RNA: enumerate admissible (u, w) with u + w <= L_tx."""
    out = np.zeros(part.n, dtype=np.float64)
    reg_of_t = part.region_of(tx.g2t)  # region of every transcript coordinate
    for w in np.nonzero(pmf)[0]:
        if w > tx.L:
            break
        u = np.arange(0, tx.L - w + 1, dtype=np.int64)
        a = reg_of_t[u]
        m = a == A
        if not m.any():
            continue
        z = reg_of_t[u[m] + w - 1]
        out += pmf[w] * np.bincount(z, minlength=part.n)
    return out


# ---------------------------------------------------------------------------
# simulation
# ---------------------------------------------------------------------------
def sim_gdna(part: Partition, rho: float, pmf: np.ndarray):
    """gDNA: start positions Poisson(rho) per bp over the whole partition."""
    genome_end = int(part.cuts[-1])
    n = RNG.poisson(rho * genome_end)
    s = RNG.integers(0, genome_end, n)
    w = RNG.choice(pmf.size, size=n, p=pmf)
    last = s + w - 1
    keep = last < genome_end
    return part.region_of(s[keep]), part.region_of(last[keep])


def sim_mature(part: Partition, tx: Transcript, rho: float, pmf: np.ndarray):
    """Mature RNA: rate(u, w) prop f(w) * 1[u + w <= L_tx]  (rejection sampling)."""
    n = RNG.poisson(rho * tx.L)
    u = RNG.integers(0, tx.L, n)
    w = RNG.choice(pmf.size, size=n, p=pmf)
    keep = u + w <= tx.L
    u, w = u[keep], w[keep]
    reg_of_t = part.region_of(tx.g2t)
    return reg_of_t[u], reg_of_t[u + w - 1]


def cell_counts(part: Partition, a_arr, z_arr) -> np.ndarray:
    """n[A, Z] integer cell counts."""
    m = np.zeros((part.n, part.n), dtype=np.int64)
    np.add.at(m, (a_arr, z_arr), 1)
    return m


# ---------------------------------------------------------------------------
# estimator: nonneg Poisson MLE over the cell design matrix
# ---------------------------------------------------------------------------
def poisson_nnls(n_obs: np.ndarray, M: np.ndarray, iters: int = 4000) -> np.ndarray:
    """min over rho >= 0 of the Poisson deviance for  n_z ~ Poisson(sum_c M[z,c] rho_c).
    Multiplicative EM update (Richardson-Lucy); converges monotonically."""
    col = M.sum(axis=0)
    rho = np.full(M.shape[1], max(n_obs.sum(), 1.0) / max(col.sum(), 1e-12))
    for _ in range(iters):
        lam = M @ rho
        lam = np.maximum(lam, 1e-12)
        rho = rho * (M.T @ (n_obs / lam)) / np.maximum(col, 1e-12)
    return rho


def fisher(M: np.ndarray, rho: np.ndarray) -> np.ndarray:
    lam = np.maximum(M @ rho, 1e-12)
    return (M / lam[:, None]).T @ M


# ===========================================================================
# EXPERIMENT A — identifiability of (rho_g, rho_r) from cells alone
# ===========================================================================
def experiment_A(rna_mu: float, reg_len: int = 150, reps: int = 40, verbose=True):
    """A chain of equal regions, all exonic (a single-exon 'transcript' spanning the
    middle), gDNA everywhere. No strand information at all. Estimate (rho_g, rho_r)
    at the CENTRE region from its own cell counts."""
    n_reg = 41
    cuts = [i * reg_len for i in range(n_reg + 1)]
    part = Partition(cuts)
    A = n_reg // 2
    tx = Transcript([(0, int(part.cuts[-1]))])  # RNA present everywhere, unspliced

    f_g = fl_pmf(60, 15)
    f_r = fl_pmf(rna_mu, rna_mu * 0.25)

    Eg = E_gdna_cells(part, A, f_g)
    Er = E_mature_cells(part, tx, A, f_r)
    live = (Eg + Er) > 1e-9
    M = np.column_stack([Eg[live], Er[live]])

    rho_g_true, rho_r_true = 0.02, 0.05
    est = []
    for _ in range(reps):
        ag, zg = sim_gdna(part, rho_g_true, f_g)
        ar, zr = sim_mature(part, tx, rho_r_true, f_r)
        cg = cell_counts(part, ag, zg)[A]
        cr = cell_counts(part, ar, zr)[A]
        est.append(poisson_nnls((cg + cr)[live], M))
    est = np.asarray(est)

    fg_true = rho_g_true / (rho_g_true + rho_r_true)
    tot = est.sum(axis=1)
    ok = tot > 0
    fg_hat = est[ok, 0] / tot[ok]
    fim = fisher(M, np.array([rho_g_true, rho_r_true]))
    cond = np.linalg.cond(fim)
    corr = fim[0, 1] / np.sqrt(fim[0, 0] * fim[1, 1])

    if verbose:
        print(f"  region {reg_len} bp | gDNA FL 60, RNA FL {rna_mu:.0f} "
              f"| cells used {int(live.sum())}")
        print(f"      design corr(g,r) = {corr:+.4f}   Fisher cond = {cond:8.1f}")
        print(f"      f_g  true {fg_true:.3f}   est {fg_hat.mean():.3f} "
              f"+/- {fg_hat.std():.3f}   bias {fg_hat.mean()-fg_true:+.4f}")
    return dict(mu=rna_mu, cond=cond, corr=corr, bias=fg_hat.mean() - fg_true,
                sd=fg_hat.std(), reg_len=reg_len)


# ===========================================================================
# EXPERIMENT B — terminus bias of the naive anchor rule
# ===========================================================================
def experiment_B(reg_len: int = 100, reps: int = 60):
    """A transcript ending at a TES. Tile the region partition across the last part of
    it. Compare f_g from (i) naive anchor E_c(A) = |A|, (ii) admissible-start E_c(A),
    (iii) the full cell estimator."""
    tx_start, tx_end = 2000, 6000
    cuts = list(range(0, 8000 + reg_len, reg_len))
    part = Partition(cuts)
    tx = Transcript([(tx_start, tx_end)])

    f_g = fl_pmf(60, 15)
    f_r = fl_pmf(200, 50)
    rho_g_true, rho_r_true = 0.02, 0.05
    fg_true = rho_g_true / (rho_g_true + rho_r_true)

    # admissible-start effective lengths, per region, per component
    genome_end = int(part.cuts[-1])
    cdf_r = np.cumsum(f_r)
    E_naive_g = part.lens.astype(float)
    E_naive_r = np.zeros(part.n)
    E_adm_g = np.zeros(part.n)
    E_adm_r = np.zeros(part.n)
    for A in range(part.n):
        s = np.arange(part.starts[A], part.ends[A])
        # gDNA: fits unless it runs off the reference
        E_adm_g[A] = cdf_r.size and np.sum(
            np.clip(np.searchsorted(np.arange(f_g.size), genome_end - s), 0, f_g.size - 1)
            * 0.0 + np.cumsum(f_g)[np.clip(genome_end - s, 0, f_g.size - 1)]
        )
        # mature RNA: needs (tx_end - s) bases remaining
        rem = tx_end - s
        inside = (s >= tx_start) & (s < tx_end)
        E_naive_r[A] = inside.sum()
        E_adm_r[A] = np.sum(cdf_r[np.clip(rem, 0, cdf_r.size - 1)] * inside)

    rows = []
    for A in range(part.n):
        if not (part.starts[A] >= tx_start and part.ends[A] <= tx_end):
            continue
        Eg = E_gdna_cells(part, A, f_g)
        Er = E_mature_cells(part, tx, A, f_r)
        live = (Eg + Er) > 1e-9
        M = np.column_stack([Eg[live], Er[live]])
        naive, adm, cell = [], [], []
        for _ in range(reps):
            ag, zg = sim_gdna(part, rho_g_true, f_g)
            ar, zr = sim_mature(part, tx, rho_r_true, f_r)
            cg = cell_counts(part, ag, zg)[A]
            cr = cell_counts(part, ar, zr)[A]
            ng, nr = cg.sum(), cr.sum()
            n_tot = ng + nr
            # (i) naive: E_g = E_r = |A|  -> the "cancellation" claim
            r_g = n_tot / E_naive_g[A]
            #   with equal E the composition is just the count split, which we cannot
            #   observe; the honest naive estimator solves the 2x2 with E_g = E_r = |A|
            #   using the SAME cell contrast, i.e. it is the cell estimator with the
            #   WRONG (uncorrected) design.  Build that design explicitly:
            Mn = np.column_stack([Eg[live], Er[live] * (E_naive_r[A] / max(E_adm_r[A], 1e-9))])
            rn = poisson_nnls((cg + cr)[live], Mn)
            naive.append(rn[0] / max(rn.sum(), 1e-12))
            # (ii) admissible-start node marginal only (1 equation, 2 unknowns ->
            #      report what you get if you assume the FL-free cancellation)
            adm.append(np.nan if E_adm_r[A] <= 0 else
                       (n_tot - rho_r_true * E_adm_r[A]) / max(n_tot, 1e-9))
            # (iii) full cell estimator with the CORRECT design
            rc = poisson_nnls((cg + cr)[live], M)
            cell.append(rc[0] / max(rc.sum(), 1e-12))
            _ = r_g
        rows.append((part.starts[A], tx_end - part.ends[A],
                     float(np.mean(naive)), float(np.mean(cell)),
                     E_adm_r[A] / max(E_naive_r[A], 1e-9)))
    print(f"  transcript [{tx_start},{tx_end}), TES at {tx_end}; regions {reg_len} bp; "
          f"true f_g = {fg_true:.3f}")
    print(f"  {'dist to TES':>12} {'E_r adm/|A|':>12} {'f_g NAIVE':>11} {'f_g CELL':>10}")
    for st, d, nv, cl, ratio in rows[-8:]:
        print(f"  {d:12d} {ratio:12.3f} {nv:11.3f} {cl:10.3f}")
    return rows, fg_true


if __name__ == "__main__":
    print("=" * 78)
    print("EXPERIMENT A — can cell counts alone separate two components by FL?")
    print("=" * 78)
    for mu in (200.0, 150.0, 100.0, 80.0, 65.0, 60.0):
        experiment_A(mu)
    print()
    print("  region-length sweep at gDNA 60 / RNA 200:")
    for rl in (60, 100, 150, 300, 600, 1500):
        experiment_A(200.0, reg_len=rl, reps=25)

    print()
    print("=" * 78)
    print("EXPERIMENT B — terminus bias of the naive anchor rule")
    print("=" * 78)
    experiment_B()
