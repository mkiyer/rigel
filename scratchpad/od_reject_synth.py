"""Dump the synthetic-suite gDNA strand seeds (TRUE od = 0 by construction) for the frame analysis."""

from __future__ import annotations

import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from od_outliers import SUITE, capture  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

OUT = Path("/Users/mkiyer/proj/rigel/scratchpad/od_reject_synth.pkl")
SYNTH = (
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_off",
)


def main():
    from selfsolve_diag import _scan_and_truth

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.index import TranscriptIndex

    cfg = PipelineConfig()
    ix = TranscriptIndex.load(str(SUITE / "rigel_index"))
    ra = RegionArrays.from_region_df(ix.region_df, ix.ref_name_to_id)
    out: dict = {}
    for cond in SYNTH:
        inp = _scan_and_truth(
            SUITE, cond, ix, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
        )
        g = capture(
            lambda inp=inp: calibrate(
                inp["payload"],
                ra,
                inp["strand_model"],
                np.asarray(inp["gdna_fl_pmf"]),
                np.asarray(inp["rna_fl_pmf"]),
                dataclasses.replace(cfg.calibration, calib_refit_iters=0),
            )
        )
        out[cond] = g
        print("captured", cond, flush=True)
    with open(OUT, "wb") as fh:
        pickle.dump(out, fh)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
