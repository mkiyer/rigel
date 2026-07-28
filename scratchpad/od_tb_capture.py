"""Capture the production gDNA-seed arrays from the 4 cached cfRNA libraries, once."""
from __future__ import annotations

import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from od_outliers import CF, capture  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

OUT = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/5a99d9c9-dd1f-4d8a-b7a3-a9d62be385dc/scratchpad/seeds.pkl")


def main():
    cfg = PipelineConfig()
    out = {}
    for name in sorted(p.stem for p in CF.glob("*.pkl")):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        g = capture(lambda d=d: calibrate(
            d["payload"], d["region_arrays"], d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        s, t, w, k, _ = g["gdna"]
        rs, rt, _, rk, _ = g["rna"]
        out[name] = dict(sense=s, total=t, weight=w, kappa=k,
                         rna_sense=rs, rna_total=rt, rna_kappa=rk)
        print(f"{name}: {t.size} seeds, sum total {t.sum():,.0f}, kappa {k:.4f}", flush=True)
    pickle.dump(out, open(OUT, "wb"))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
