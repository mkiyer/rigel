"""WHAT IS THE RNA STRAND OVERDISPERSION ON REAL DATA?  Measured at junctions deep enough to see it.

OWNER'S FRAMING (2026-07-28), which is the correct one:

  During strand-specific prep (dUTP) the reverse-strand synthesis is degraded so fragments can be traced
  to their strand. It is imperfect. The ALIGNER only reports Forward/Reverse — that is not a strand call.
  But IFF a fragment is SPLICED, the genomic motif (GT/AG) gives the TRUE strand independently of prep.
  So over spliced fragments, comparing motif strand (truth) to aligner orientation MEASURES LIBRARY-PREP
  EFFICIENCY. Perfect prep = 100 % agreement; inefficient prep = less.

  `kappa` is that mean. What we need is the DISPERSION around it, because that is what sets how sharply
  the strand channel can separate gDNA from RNA in the UNSPLICED population at each node.

  ⚠ And the owner's power point is the crux: at kappa = 0.99 you need a junction with LOTS of reads to
  see even one disagreeing read. So the estimate can only come from DEEP junctions, and pooling shallow
  ones in can only add noise.

This script therefore reports, per library:
  * the aggregate agreement rate of the seeds vs the `kappa` the fit is handed (are they the same
    observable at all?);
  * the per-junction agreement distribution RESTRICTED BY DEPTH;
  * the implied overdispersion from the deep junctions ALONE, which is the honest estimate.

Run: OMP_NUM_THREADS=1 python scratchpad/od_rna_depth.py
"""
from __future__ import annotations

import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from od_outliers import CF, capture  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

_EPS = 1e-12
DEPTHS = (1, 10, 100, 1000)


def od_at(sense, total, kappa, depth):
    """Pooled MoM overdispersion using ONLY seeds with total >= depth, at mean kappa."""
    m = total >= depth
    if m.sum() == 0:
        return np.nan, 0, np.nan
    s, t = sense[m], total[m]
    ex = (s - t * kappa) ** 2 - t * kappa * (1 - kappa)
    sc = np.maximum(t * (t - 1.0), 0.0) * kappa * (1 - kappa)
    good = sc > 0
    od = ex[good].sum() / max(sc[good].sum(), _EPS) if good.any() else np.nan
    return od, int(m.sum()), float(s.sum() / max(t.sum(), _EPS))


def main():
    cfg = PipelineConfig()
    print("=== R1. IS THE SEED POPULATION THE SAME OBSERVABLE AS kappa? ===\n")
    print(f"{'sample':10s} {'kappa (StrandModel,':>20s} {'seeds aggregate':>17s} {'ratio':>8s} "
          f"{'n seeds':>10s} {'total frags':>13s}")
    print(f"{'':10s} {'read-1 sense MLE)':>20s} {'agreement rate':>17s}")
    keep = {}
    for name in sorted(p.stem for p in CF.glob("*.pkl")):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        g = capture(lambda d=d: calibrate(
            d["payload"], d["region_arrays"], d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        s, t, _, k, _ = g["rna"]
        ok = t > 0
        s, t = s[ok], t[ok]
        keep[name] = (s, t, k)
        agg = s.sum() / max(t.sum(), _EPS)
        print(f"{name:10s} {k:20.5f} {agg:17.5f} {agg / max(k, _EPS):7.1f}x {t.size:10,d} "
              f"{t.sum():13,.0f}")

    print("\n\n=== R2. THE AGREEMENT DISTRIBUTION BY JUNCTION DEPTH ===")
    print("    At kappa near 0 or 1 a shallow junction CANNOT show disagreement -- it is all-or-nothing")
    print("    by arithmetic. Only deep junctions carry information about the dispersion.\n")
    for name, (s, t, k) in keep.items():
        print(f"  {name}   (kappa = {k:.5f})")
        print(f"    {'depth':>10s} {'n junc':>9s} {'frags':>12s} {'agg rate':>9s} "
              f"{'p10':>7s} {'p50':>7s} {'p90':>7s} {'sd':>7s}")
        for dep in DEPTHS:
            m = t >= dep
            if m.sum() == 0:
                continue
            sf = s[m] / t[m]
            print(f"    {'>= ' + str(dep):>10s} {int(m.sum()):9,d} {t[m].sum():12,.0f} "
                  f"{s[m].sum() / max(t[m].sum(), _EPS):9.4f} "
                  f"{np.quantile(sf, .1):7.3f} {np.quantile(sf, .5):7.3f} "
                  f"{np.quantile(sf, .9):7.3f} {sf.std():7.3f}")
        print()

    print("\n=== R3. THE HONEST OVERDISPERSION -- deep junctions only, and at the RIGHT mean ===")
    print("    'at kappa' uses the fitted kappa; 'at own mean' re-centres on the depth-restricted")
    print("    aggregate, which isolates DISPERSION from any MEAN misfit.\n")
    print(f"{'sample':10s} {'depth':>8s} {'n':>8s} {'od at kappa':>13s} {'od at own mean':>16s} "
          f"{'implied Beta(a,a)':>18s}")
    for name, (s, t, k) in keep.items():
        for dep in DEPTHS:
            m = t >= dep
            if m.sum() < 5:
                continue
            own = s[m].sum() / max(t[m].sum(), _EPS)
            od_k, n, _ = od_at(s, t, k, dep)
            od_o, _, _ = od_at(s, t, own, dep)
            a = (1 / od_o - 1) / 2 if (od_o and od_o > 0) else np.nan
            print(f"{name:10s} {dep:8d} {n:8,d} {od_k:13.4f} {od_o:16.4f} {a:18.2f}")
        print()
    print("  od at own mean is the quantity the owner asked for: given the mean, what is the variance?")
    print("  ⚠ od is an intraclass CORRELATION and is bounded by 1. Anything above 1 is not a")
    print("    dispersion at all -- it is a mean misfit being reported as one.")


if __name__ == "__main__":
    main()
