"""Dump the EXACT production gDNA/RNA strand-overdispersion seed arrays for the 4 real cfRNA libs.

Written for the rejection-criterion analysis (Bonferroni vs Tarone vs BH vs robust weighting).
Run: OMP_NUM_THREADS=1 python scratchpad/od_reject_dump.py
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

OUT = Path("/Users/mkiyer/proj/rigel/scratchpad/od_reject_seeds.pkl")


def main():
    cfg = PipelineConfig()
    out: dict = {}
    for name in sorted(p.stem for p in CF.glob("*.pkl")):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        g = capture(
            lambda d=d: calibrate(
                d["payload"],
                d["region_arrays"],
                d["strand_model"],
                np.asarray(d["gdna_fl_pmf"]),
                np.asarray(d["rna_fl_pmf"]),
                dataclasses.replace(cfg.calibration, calib_refit_iters=0),
            )
        )
        out[name] = g
        print("captured", name, flush=True)
    with open(OUT, "wb") as fh:
        pickle.dump(out, fh)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
