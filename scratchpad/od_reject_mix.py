"""(1) Correct within-stratum pooled od.  (2) Prototype: TWO-COMPONENT responsibility weighting.

The contamination is not an arbitrary displacement — it is RNA, whose sense rate is the ALREADY-FITTED
kappa (or 1-kappa, depending on the contaminating gene's strand).  So the seed population is a mixture
of two components with KNOWN means:

    G : BetaBinom(n, 1/2,  od_g)      genomic DNA          (mean fixed by biology)
    R : BetaBinom(n, kappa or 1-kappa, od_r)   transcription  (mean fitted from spliced, pure RNA)

Fit (pi, od_g) by EM.  The seed's gDNA responsibility r_s replaces both the hard reject AND the
annotation-derived gdna_weight (which is measured to be 1.000 for 100% of seeds -> zero information).

Key property: the two components separate exactly when kappa is far from 1/2 -- i.e. the screen has
power exactly when the strand likelihood itself has power.  On unstranded data the components coincide,
pi is unidentifiable, responsibilities go flat, and this reduces to the pooled estimator.
"""

from __future__ import annotations

import pickle
import sys

import numpy as np
from scipy.special import logsumexp
from scipy.stats import betabinom

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")

REAL = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_seeds.pkl"
SYNTH = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_synth.pkl"
CEIL = 0.2
_EPS = 1e-12


def ab(mu, rho):
    c = (1.0 - rho) / rho
    return mu * c, (1.0 - mu) * c


def strata(tag, sense, total):
    ok = total > 0
    s, t = sense[ok], total[ok]
    allp = np.sum(t * (t - 1) / 2.0)
    cells = []
    for lo, hi in ((2, 2), (3, 4), (5, 8), (9, 16), (17, 32), (33, 64), (65, 128),
                   (129, 256), (257, 512), (513, 10**9)):
        m = (t >= lo) & (t <= hi)
        if not m.any():
            continue
        n, k = t[m], s[m]
        num = np.sum((k - n / 2.0) ** 2 - n / 4.0)
        den = np.sum(n * (n - 1.0) / 4.0)
        cells.append((f"{lo}-{hi}", int(m.sum()), np.sum(n * (n - 1) / 2) / allp, num / den))
    print(f"\n  {tag}")
    print("     " + " ".join(f"{c[0]:>9}" for c in cells))
    print("     " + " ".join(f"{c[1]:>9d}" for c in cells) + "   #seeds")
    print("     " + " ".join(f"{c[2]:>8.1%} " for c in cells) + "  %pairs")
    print("     " + " ".join(f"{c[3]:>9.4f}" for c in cells) + "   od (pooled within stratum)")


def em_mixture(sense, total, kappa, od_r=0.05, iters=60, od0=0.05, pi0=0.9):
    n = total
    k = sense
    ok = n > 0
    n, k = n[ok], k[ok]
    od_g, pi = od0, pi0
    kap = min(max(kappa, 1e-4), 1 - 1e-4)
    for _ in range(iters):
        a1, b1 = ab(0.5, max(od_g, 1e-4))
        lg = betabinom.logpmf(k, n, a1, b1)
        a2, b2 = ab(kap, od_r)
        a3, b3 = ab(1 - kap, od_r)
        lr = logsumexp(
            np.stack([betabinom.logpmf(k, n, a2, b2), betabinom.logpmf(k, n, a3, b3)]),
            axis=0,
            b=0.5,
        )
        num = np.log(max(pi, 1e-12)) + lg
        den = np.logaddexp(num, np.log(max(1 - pi, 1e-12)) + lr)
        r = np.exp(num - den)
        pi = float(np.mean(r))
        w = r * np.maximum(n * (n - 1.0), 0.0) * 0.25
        ex = r * ((k - n / 2.0) ** 2 - n * 0.25)
        new = float(np.sum(ex) / max(np.sum(w), _EPS))
        new = min(max(new, 1e-4), CEIL)
        if abs(new - od_g) < 1e-7:
            od_g = new
            break
        od_g = new
    return od_g, pi, r, n, k


def main():
    real = pickle.load(open(REAL, "rb"))
    syn = pickle.load(open(SYNTH, "rb"))
    print("=" * 100)
    print("A. CORRECT within-stratum pooled od  (od <= 1 always; >1 earlier was a binning artefact)")
    print("=" * 100)
    for nm in ("LBX0190", "LBX0588", "MO_3021", "vcap"):
        s, t, w, kap, _ = real[nm]["gdna"]
        strata(f"REAL {nm}  (kappa={kap:.4f})", s, t)
    for nm, v in syn.items():
        s, t, w, kap, _ = v["gdna"]
        strata(f"SYNTH {nm[9:]}  (kappa={kap:.4f})", s, t)

    print()
    print("=" * 100)
    print("B. TWO-COMPONENT EM: od_g and the estimated contaminated fraction (1-pi)")
    print("=" * 100)
    print(f"  {'library':44s} {'kappa':>7} {'pooled':>9} {'od_g(EM)':>9} {'1-pi':>8} "
          f"{'%pairs in R':>12}")
    for nm in ("LBX0190", "LBX0588", "MO_3021", "vcap"):
        s, t, w, kap, _ = real[nm]["gdna"]
        pooled = float(np.sum((s - t / 2) ** 2 - t / 4) / np.sum(t * (t - 1) / 4))
        od_g, pi, r, n, k = em_mixture(s, t, kap)
        pr = 1 - np.sum(r * n * (n - 1) / 2) / np.sum(n * (n - 1) / 2)
        print(f"  REAL  {nm:38s} {kap:7.4f} {pooled:9.4f} {od_g:9.4f} {1 - pi:8.3%} {pr:12.2%}")
    for nm, v in syn.items():
        s, t, w, kap, _ = v["gdna"]
        pooled = float(np.sum((s - t / 2) ** 2 - t / 4) / np.sum(t * (t - 1) / 4))
        od_g, pi, r, n, k = em_mixture(s, t, kap)
        pr = 1 - np.sum(r * n * (n - 1) / 2) / np.sum(n * (n - 1) / 2)
        print(f"  SYNTH {nm[9:]:38s} {kap:7.4f} {pooled:9.4f} {od_g:9.4f} {1 - pi:8.3%} "
              f"{pr:12.2%}   TRUTH=0")

    print()
    print("  CONTROL — an UNSTRANDED library (kappa=1/2) must make the mixture unidentifiable:")
    s, t, w, kap, _ = syn["gdna_gdna100_ss_0.50_nrna_none_capture_off"]["gdna"]
    od_g, pi, r, n, k = em_mixture(s, t, 0.5)
    print(f"     forced kappa=0.5 -> od_g={od_g:.4f}, 1-pi={1 - pi:.3%} "
          f"(responsibilities flat: sd(r)={np.std(r):.4f})")


if __name__ == "__main__":
    main()
