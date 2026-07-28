"""WHAT SHAPE IS THE CONTAMINATION?  The measurement that decides the robust-estimator design.

`od_screen.py` showed a seed-consistency screen fixes 2 of 4 real samples and leaves LBX0190 / MO_3021
pinned at the ceiling even under the MEDIAN of per-seed moments (a 50 % breakdown point).  That means one
of two very different things, and they call for different estimators:

  (A) TAIL contamination — most seeds are clean, a minority are displaced. A robust SCALE (IQR/MAD of z)
      still measures the true dispersion, and the design is "fit the scale from the centre".
  (B) PERVASIVE displacement — most seeds are individually inconsistent with Binomial(1/2). No robust
      estimator can help, because there is no clean majority; the SEED SET itself is wrong.

The discriminator is the distribution of the seed's own standardized deviation

        z = (sense − N·mean) / sqrt(N·mean·(1−mean))

Under the null (pure gDNA, no overdispersion) z ~ N(0,1): median 0, IQR 1.349. Overdispersion inflates
the IQR symmetrically. Contamination by a stranded transcript DISPLACES the mean — a different moment.
That asymmetry is the whole basis of a principled screen, so measure it before designing one.

Run: OMP_NUM_THREADS=1 python scratchpad/od_shape.py
"""
from __future__ import annotations

import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from od_outliers import CF, SUITE, capture  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

_EPS = 1e-12
_IQR_N01 = 1.3489795


def zed(sense, total, weight, kappa):
    node_mean = 0.5 * weight + kappa * (1.0 - weight)
    sd = np.sqrt(np.maximum(total * node_mean * (1.0 - node_mean), _EPS))
    ok = total > 0.0
    return ((sense - total * node_mean) / sd)[ok], total[ok], weight[ok], (sense / np.maximum(total, _EPS))[ok]


def main():
    cfg = PipelineConfig()
    runs = {}
    for name in sorted(p.stem for p in CF.glob("*.pkl")):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        runs[f"REAL  {name}"] = capture(lambda d=d: calibrate(
            d["payload"], d["region_arrays"], d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
    from selfsolve_diag import _scan_and_truth

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.index import TranscriptIndex
    ix = TranscriptIndex.load(str(SUITE / "rigel_index"))
    ra = RegionArrays.from_region_df(ix.region_df, ix.ref_name_to_id)
    for cond in ("gdna_gdna100_ss_0.99_nrna_present_capture_on",
                 "gdna_gdna100_ss_0.50_nrna_none_capture_off"):
        inp = _scan_and_truth(SUITE, cond, ix, cfg, Path("/tmp/rigel_selfsolve"),
                              SUITE / "_selfsolve_cache")
        runs[f"SYNTH {cond[5:26]}"] = capture(lambda inp=inp: calibrate(
            inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
            np.asarray(inp["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))

    print("=== S1. THE z DISTRIBUTION of the gDNA seeds. Null: median 0, IQR 1.349 ===\n")
    print(f"{'sample':30s} {'n':>7s} {'med z':>8s} {'IQR z':>8s} {'IQR/null':>9s} | "
          f"{'|z|>2':>7s} {'|z|>3':>7s} {'|z|>5':>7s} | {'med N':>7s}")
    keep = {}
    for tag, g in runs.items():
        s, t, w, k, _ = g["gdna"]
        z, tt, ww, sf = zed(s, t, w, k)
        keep[tag] = (z, tt, ww, sf)
        iqr = float(np.subtract(*np.percentile(z, [75, 25])))
        print(f"{tag:30s} {z.size:7d} {np.median(z):8.3f} {iqr:8.3f} {iqr / _IQR_N01:9.2f} | "
              f"{np.mean(np.abs(z) > 2):6.1%} {np.mean(np.abs(z) > 3):6.1%} {np.mean(np.abs(z) > 5):6.1%} | "
              f"{np.median(tt):7.0f}")

    print("\n   IQR/null = 1 means the CENTRE of the seed population is exactly Binomial(1/2) — no")
    print("   dispersion at all. >1 means genuine spread. |z|>3 well above ~0.3 % is contamination.")

    print("\n\n=== S2. IS IT DISPLACEMENT OR SPREAD?  the seeds' own sense fraction ===\n")
    print(f"{'sample':30s} {'med sf':>8s} {'IQR sf':>8s} | {'sf<0.1':>7s} {'sf>0.9':>7s} "
          f"{'|sf-.5|>0.3':>12s} | {'their N':>8s} {'their wt':>9s}")
    for tag, (z, tt, ww, sf) in keep.items():
        ext = np.abs(sf - 0.5) > 0.3
        print(f"{tag:30s} {np.median(sf):8.3f} {float(np.subtract(*np.percentile(sf, [75, 25]))):8.3f} | "
              f"{np.mean(sf < 0.1):6.1%} {np.mean(sf > 0.9):6.1%} {np.mean(ext):11.1%} | "
              f"{np.median(tt[ext]) if ext.any() else np.nan:8.0f} "
              f"{np.median(ww[ext]) if ext.any() else np.nan:9.3f}")
    print("\n   gDNA is symmetric BY CONSTRUCTION, so a seed far from sf = 0.5 is displaced, not")
    print("   dispersed. 'their N' is how many fragments those seeds bring into the pooled fit.")

    print("\n\n=== S3. WHERE THE FRAGMENTS ARE (the pooled fit is fragment-weighted, not seed-weighted) ===\n")
    print(f"{'sample':30s} {'seeds |sf-.5|>0.3':>18s} {'their FRAGMENT share':>21s}")
    for tag, (z, tt, ww, sf) in keep.items():
        ext = np.abs(sf - 0.5) > 0.3
        print(f"{tag:30s} {np.mean(ext):17.1%} {tt[ext].sum() / max(tt.sum(), _EPS):20.1%}")
    print("\n   A small SEED share with a large FRAGMENT share is case (A): trimmable.")
    print("   A large seed share is case (B): the seed SET is wrong and no estimator can rescue it.")


if __name__ == "__main__":
    main()
