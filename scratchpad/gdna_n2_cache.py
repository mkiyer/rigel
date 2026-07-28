"""N2 — build the PAIRED pass-0 / re-solved belief cache for the iterative-AMBIG test.

HANDOFF_17 §0 proposes resolving the circularity-vs-enriched-mass conflict with the architecture that
already exists:  pass-0 -> prior #1 (no AMBIG) -> RE-SOLVE -> prior #2 (including the re-solved AMBIG),
which is non-circular by construction because the AMBIG estimates prior #2 trains on were produced by a
prior that never saw them.

To measure that at all we need the re-solved beliefs, which no cache carries. This adds them, per node,
alongside pass-0's — same conditions, same substrate, same scan cache (`_selfsolve_cache`, warm), so the
only difference between the two belief sets is `calib_refit_iters`.

⚠ **Prior #1 here is the SHIPPED `DensityNPMLE`** — the object this workstream is retiring. So this
measures whether the ARCHITECTURE delivers usable AMBIG estimates, not whether the final prior will. If it
does deliver, the measurement must be repeated once W4 wires the landscape in.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n2_cache.py [--suite <dir>] [--out <pkl>]
"""
from __future__ import annotations

import argparse
import dataclasses
import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
ap.add_argument("--out", default="/Users/mkiyer/proj/rigel/scratchpad/gdna_refit_cache.pkl")
args = ap.parse_args()

suite = Path(args.suite)
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
cfg = PipelineConfig()
work = Path("/tmp/rigel_selfsolve")
cache = suite / "_selfsolve_cache"
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))

out = {}
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    rec = {}
    for iters in (0, 1):
        dbg = {}
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=iters), _debug=dbg)
        cap = dbg["capture"]
        rec[iters] = dict(
            f0=np.asarray(dbg["belief"].f_g, np.float64),
            var=np.asarray(dbg["belief"].var_gdna, np.float64),
            mass=np.asarray(cap["mass_global"], np.float64),
            eff=np.asarray(cap["eff_global"], np.float64),
        )
        if iters == 0:
            rec["is_region"] = np.asarray(dbg["chain"].kind) == REGION
            rec["fp"] = np.asarray(cap["free_pos"], bool)
            rec["fn"] = np.asarray(cap["free_neg"], bool)
    out[cond] = rec
    d0, d1 = rec[0], rec[1]
    amb = rec["fp"] & rec["fn"]
    moved = float(np.mean(np.abs(d1["f0"][amb] - d0["f0"][amb]))) if amb.any() else float("nan")
    print(f"  {cond:52s} AMBIG {int(amb.sum()):5d}  mean |Δf_g| pass0→refit {moved:.4f}")

Path(args.out).write_bytes(pickle.dumps(out))
print(f"\nwrote {len(out)} conditions -> {args.out}")
