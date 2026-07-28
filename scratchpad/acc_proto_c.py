"""EXPERIMENT C — anchor marginal vs coverage marginal at a SHARP density step.

The node-vs-edge allocation question reduces to which summary of the path counts a
node should receive:

  ANCHOR    n(A) = # paths whose first node is A,  divided by the admissible-start
            effective length.  Exact under a Poisson-start process; uses each
            fragment once; few fragments at a short node.

  COVERAGE  V(A) = sum over ALL covering fragments of o_A/L, divided by |A|.
            Uses every fragment that touches the node -> far lower variance,
            but it estimates rho SMOOTHED over one fragment length, and the two
            components smooth over DIFFERENT lengths (gDNA 60 vs RNA 200), so the
            COMPOSITION is biased wherever rho has a step.

A sharp step is exactly an exon boundary under hybrid capture. Measure both.
"""

from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(20260728)


def fl_pmf(mu, sd, n=1200):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return p / p.sum()


def run(reg_len=100, n_reg=61, reps=30, rna_mu=200.0, gdna_mu=60.0, block=(25, 35)):
    """rho_G uniform everywhere; rho_R nonzero only on regions [block) -> a sharp step."""
    cuts = np.arange(n_reg + 1, dtype=np.int64) * reg_len
    starts, ends = cuts[:-1], cuts[1:]
    genome_end = int(cuts[-1])
    f_g, f_r = fl_pmf(gdna_mu, gdna_mu * 0.25), fl_pmf(rna_mu, rna_mu * 0.25)

    rho_g, rho_r = 0.02, 0.05
    lo, hi = block
    rna_lo, rna_hi = int(cuts[lo]), int(cuts[hi])          # the RNA support, contiguous

    def region_of(pos):
        return np.searchsorted(cuts, pos, side="right") - 1

    # ---- admissible-start effective length per region, per component -------
    cdf_g, cdf_r = np.cumsum(f_g), np.cumsum(f_r)
    E_anchor_g = np.zeros(n_reg)
    E_anchor_r = np.zeros(n_reg)
    for a in range(n_reg):
        s = np.arange(starts[a], ends[a])
        E_anchor_g[a] = cdf_g[np.clip(genome_end - s, 0, cdf_g.size - 1)].sum()
        inside = (s >= rna_lo) & (s < rna_hi)
        E_anchor_r[a] = (cdf_r[np.clip(rna_hi - s, 0, cdf_r.size - 1)] * inside).sum()

    acc_anchor, acc_cov = [], []
    for _ in range(reps):
        anchor_n = np.zeros(n_reg)
        cov_v = np.zeros(n_reg)
        for rho, pmf, sup_lo, sup_hi in ((rho_g, f_g, 0, genome_end),
                                         (rho_r, f_r, rna_lo, rna_hi)):
            n = RNG.poisson(rho * (sup_hi - sup_lo))
            s = RNG.integers(sup_lo, sup_hi, n)
            w = RNG.choice(pmf.size, size=n, p=pmf)
            keep = s + w <= sup_hi
            s, w = s[keep], w[keep]
            # anchor: +1 at the region containing the first base
            np.add.at(anchor_n, region_of(s), 1.0)
            # coverage: + o_n / L to every region the fragment touches
            for si, wi in zip(s, w):
                a0, a1 = region_of(si), region_of(si + wi - 1)
                for a in range(a0, a1 + 1):
                    o = min(ends[a], si + wi) - max(starts[a], si)
                    cov_v[a] += o / wi
        acc_anchor.append(anchor_n)
        acc_cov.append(cov_v)

    A = np.mean(acc_anchor, axis=0)
    C = np.mean(acc_cov, axis=0)
    Asd = np.std(acc_anchor, axis=0)
    Csd = np.std(acc_cov, axis=0)

    truth_fg = np.where((np.arange(n_reg) >= lo) & (np.arange(n_reg) < hi),
                        rho_g / (rho_g + rho_r), 1.0)

    # ANCHOR composition: solve rho from the two admissible-start eff lengths.
    # (2 unknowns, 1 equation per node -> use the ORACLE split for a bias-only read:
    #  what f_g does each summary imply if the component counts were separable?)
    #  Here we report the density each summary recovers for the TOTAL, and the
    #  implied f_g using each component's own effective length.
    fg_anchor = np.zeros(n_reg)
    fg_cov = np.zeros(n_reg)
    for a in range(n_reg):
        # anchor: rho_hat_c = n_c / E_anchor_c  -> but n_c is unobserved; use the
        # KNOWN generative split to isolate ATTRIBUTION bias from estimation noise.
        ng = rho_g * E_anchor_g[a]
        nr = rho_r * E_anchor_r[a]
        rg = ng / max(E_anchor_g[a], 1e-9)
        rr = nr / max(E_anchor_r[a], 1e-9)
        fg_anchor[a] = rg / max(rg + rr, 1e-12)
        # coverage: the expected coverage statistic, smoothed by each component's FL
        def cov_expect(rho, pmf, sup_lo, sup_hi):
            tot = 0.0
            xs = np.arange(starts[a], ends[a])
            for w in np.nonzero(pmf)[0]:
                s_lo = np.maximum(xs - w + 1, sup_lo)
                s_hi = np.minimum(xs, sup_hi - w)
                tot += pmf[w] * np.maximum(0, s_hi - s_lo + 1).sum() / w
            return rho * tot
        cg = cov_expect(rho_g, f_g, 0, genome_end)
        cr = cov_expect(rho_r, f_r, rna_lo, rna_hi)
        fg_cov[a] = (cg / reg_len) / max(cg / reg_len + cr / reg_len, 1e-12)

    print(f"regions {reg_len} bp | gDNA FL {gdna_mu:.0f} rho {rho_g} | "
          f"RNA FL {rna_mu:.0f} rho {rho_r} on regions [{lo},{hi})")
    print(f"{'region':>7} {'true f_g':>9} {'ANCHOR f_g':>11} {'COVERAGE f_g':>13} "
          f"{'anchor CV':>10} {'cov CV':>8}")
    for a in range(lo - 3, hi + 3):
        cvA = Asd[a] / max(A[a], 1e-9)
        cvC = Csd[a] / max(C[a], 1e-9)
        mark = "  <- step" if a in (lo - 1, lo, hi - 1, hi) else ""
        print(f"{a:7d} {truth_fg[a]:9.3f} {fg_anchor[a]:11.3f} {fg_cov[a]:13.3f} "
              f"{cvA:10.3f} {cvC:8.3f}{mark}")

    interior = np.arange(lo + 3, hi - 3)
    print()
    print(f"  INTERIOR of the RNA block (regions {lo+3}..{hi-4}):")
    print(f"    anchor   f_g bias {np.mean(fg_anchor[interior]-truth_fg[interior]):+.4f}"
          f"   mean CV {np.mean(Asd[interior]/np.maximum(A[interior],1e-9)):.3f}")
    print(f"    coverage f_g bias {np.mean(fg_cov[interior]-truth_fg[interior]):+.4f}"
          f"   mean CV {np.mean(Csd[interior]/np.maximum(C[interior],1e-9)):.3f}")
    edge = np.array([lo - 1, lo, hi - 1, hi])
    print(f"  AT THE STEP (regions {edge.tolist()}):")
    print(f"    anchor   f_g bias {np.mean(np.abs(fg_anchor[edge]-truth_fg[edge])):+.4f} (abs)")
    print(f"    coverage f_g bias {np.mean(np.abs(fg_cov[edge]-truth_fg[edge])):+.4f} (abs)")


if __name__ == "__main__":
    print("=" * 84)
    run(reg_len=100)
    print()
    print("=" * 84)
    run(reg_len=25, n_reg=241, block=(100, 140), reps=15)
