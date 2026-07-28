"""Split the seeds by CLASS (contained region vs boundary side) and stratify each by size.

`fit_gdna_strand_from_substrate` concatenates region seeds then boundary-side seeds, so the class
boundary is just the length of the region array.  Recompute both from the cached payloads.
"""

from __future__ import annotations

import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")

from rigel.calibration import gdna_strand as GS  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")


def capture_split(run):
    grabbed: dict = {}
    orig_r, orig_b = GS._region_seeds, GS.boundary_side_seeds

    def r_spy(*a, **k):
        out = orig_r(*a, **k)
        grabbed["region"] = tuple(np.asarray(x, float).copy() for x in out)
        return out

    def b_spy(*a, **k):
        out = orig_b(*a, **k)
        grabbed["boundary"] = tuple(np.asarray(x, float).copy() for x in out)
        return out

    GS._region_seeds, GS.boundary_side_seeds = r_spy, b_spy
    try:
        run()
    finally:
        GS._region_seeds, GS.boundary_side_seeds = orig_r, orig_b
    return grabbed


def strat(tag, sense, total):
    ok = total > 0
    s, t = sense[ok], total[ok]
    print(f"\n  {tag}   seeds={s.size}")
    print(f"     {'n':>10} {'#seeds':>9} {'%pairs':>8} {'od_hat':>9}")
    allp = np.sum(t * (t - 1) / 2.0)
    for lo, hi in ((2, 2), (3, 3), (4, 4), (5, 6), (7, 10), (11, 20), (21, 50), (51, 200),
                   (201, 10**9)):
        m = (t >= lo) & (t <= hi)
        if not m.any():
            continue
        n, k = t[m], s[m]
        nb = np.mean(n)
        od = (4.0 * np.mean((k - n / 2.0) ** 2) / nb - 1.0) / max(nb - 1.0, 1e-12)
        print(
            f"     {f'{lo}-{hi}':>10} {int(m.sum()):9d} "
            f"{np.sum(n * (n - 1) / 2.0) / max(allp, 1e-12):8.2%} {od:9.4f}"
        )


def main():
    cfg = PipelineConfig()
    for name in ("LBX0190", "vcap"):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        g = capture_split(
            lambda d=d: calibrate(
                d["payload"],
                d["region_arrays"],
                d["strand_model"],
                np.asarray(d["gdna_fl_pmf"]),
                np.asarray(d["rna_fl_pmf"]),
                dataclasses.replace(cfg.calibration, calib_refit_iters=0),
            )
        )
        print("=" * 80)
        print(name)
        print("=" * 80)
        for cls in ("region", "boundary"):
            s, t, w = g[cls]
            print(f"  [{cls}] n={t.size}, w>=0.999: {np.mean(w > 0.999):.2%}, "
                  f"integral N: {np.mean(np.abs(t - np.rint(t)) < 1e-6):.2%}, "
                  f"sum N={t.sum():,.0f}")
            strat(cls, s, t)


if __name__ == "__main__":
    main()
