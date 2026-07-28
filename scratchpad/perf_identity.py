"""FAST bit-identity gate: hash every solved per-node array, not the mwae summary.

`verify_clean.sh` compares two aggregate statistics per condition (~11 min per refit). This compares the
FULL solved state — f_g / f_pos / f_neg and all three variances, node by node — over a chosen subset of
conditions, so it catches a perturbation the aggregates round away and it runs in seconds per condition.
Use it as the inner loop; `verify_clean.sh` stays the landing gate.

    OMP_NUM_THREADS=1 python scratchpad/perf_identity.py --save base        # record
    OMP_NUM_THREADS=1 python scratchpad/perf_identity.py --check base       # compare
    ... --conds 8   (default: 8 conditions spanning stranded/unstranded x capture, both refits)
"""
from __future__ import annotations

import argparse
import dataclasses
import importlib
import os
import pickle
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
STORE = Path("/tmp/rigel_perf_identity")
STORE.mkdir(exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--save")
ap.add_argument("--check")
ap.add_argument("--conds", type=int, default=8)
a = ap.parse_args()
name = a.save or a.check
assert name, "pass --save NAME or --check NAME"

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
allc = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
# span the axes the owner tracks: stranded/unstranded x capture off/on x gDNA depth
conds = [allc[i] for i in np.linspace(0, len(allc) - 1, a.conds).round().astype(int)]

out: dict = {}
t0 = time.perf_counter()
for cond in conds:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    for refit in (0, 1):
        dbg: dict = {}
        cc = dataclasses.replace(cfg.calibration, calib_refit_iters=refit)
        res = calmod.calibrate(
            inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
            np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
        )
        cap = dbg["capture"]
        rec = {k: np.asarray(cap[k], np.float64).copy() for k in
               ("f_g", "f_pos", "f_neg", "var_g", "vg_loc", "vp_loc", "mode_g", "mode_p", "mode_n",
                "prec_g", "prec_p", "prec_n", "mass_global")}
        for k in ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg"):
            rec["bel_" + k] = np.asarray(getattr(dbg["belief"], k), np.float64).copy()
        for f in dataclasses.fields(res):
            v = getattr(res, f.name)
            if isinstance(v, np.ndarray) and v.dtype.kind == "f":
                rec["res_" + f.name] = v.astype(np.float64).copy()
            elif isinstance(v, float):
                rec["res_" + f.name] = np.array([v])
        out[(cond, refit)] = rec
print(f"{len(conds)} conditions x 2 refits in {time.perf_counter() - t0:.1f} s")

p = STORE / f"{name}.pkl"
if a.save:
    pickle.dump(out, open(p, "wb"))
    print(f"saved -> {p}  ({sum(len(v) for v in out.values())} arrays)")
else:
    ref = pickle.load(open(p, "rb"))
    bad = tot = 0
    for key, rec in out.items():
        r = ref.get(key)
        if r is None:
            print(f"⛔ missing baseline for {key}")
            bad += 1
            continue
        for k, v in rec.items():
            tot += 1
            b = r.get(k)
            if b is None or b.shape != v.shape or not np.array_equal(
                v.view(np.int64), b.view(np.int64)
            ):
                nd = -1 if b is None or b.shape != v.shape else int(
                    (v.view(np.int64) != b.view(np.int64)).sum()
                )
                print(f"⛔ DIFFERS {key} {k}: {nd} of {v.size} elements")
                bad += 1
    print(f"{tot - bad}/{tot} arrays BIT-IDENTICAL" + ("" if not bad else f"   ⛔ {bad} DIFFER"))
    raise SystemExit(1 if bad else 0)
