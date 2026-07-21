"""Precompute the per-node SUBSTRATE for the whole scenario battery (pass-0 once) and pickle it, so landscape
exploration (workflow agents) is pure numpy with NO calibrate re-runs. Includes BOTH region and BOUNDARY nodes.
Run once:  OMP_NUM_THREADS=1 python scripts/debug/gdna_cache_build.py"""
from __future__ import annotations

import dataclasses
import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

import argparse  # noqa: E402

_SCR = "/private/tmp/claude-503/-Users-mkiyer-proj-rigel/4f7a248b-0c78-4b40-9030-462373aefb19/scratchpad"
_ap = argparse.ArgumentParser()
_ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_ap.add_argument("--out", default=f"{_SCR}/gdna_substrate_cache.pkl")
_args = _ap.parse_args()
OUT = Path(_args.out)
suite = Path(_args.suite)
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
cfg = PipelineConfig()
work = Path("/tmp/rigel_selfsolve")
cache = suite / "_selfsolve_cache"
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))
rtype_all = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)

scen = []
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]),
              dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain = dbg["chain"]
    cap = dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    isr = np.asarray(chain.kind) == REGION
    idx = np.asarray(chain.ref_idx, np.int64)
    ntype = np.where(isr, rtype_all[np.clip(idx, 0, len(rtype_all) - 1)], 3)  # 0/1/2 region-type, 3 boundary
    dbkt = next((k for k in ("none", "gdna1000", "gdna300", "gdna100", "gdna5", "gdna1") if k in cond), "?")
    f0 = np.asarray(dbg["belief"].f_g, float)
    scen.append(dict(
        cond=cond,
        group=("ON" if "capture_on" in cond else "OFF", dbkt,
               "0.99" if "0.99" in cond else "0.50",
               "nrna" if ("present" in cond or "rnd" in cond) else "none"),
        is_region=isr,
        ntype=ntype.astype(np.int8),
        fp=np.asarray(cap["free_pos"], bool),
        fn=np.asarray(cap["free_neg"], bool),
        mass=np.asarray(cap["mass_global"], np.float64),
        eff=np.asarray(cap["eff_global"], np.float64),
        f0=f0,
        g_hat=(f0 * np.asarray(cap["mass_global"], np.float64)),
        var=np.asarray(dbg["belief"].var_gdna, np.float64),
        G=(Gp + Gn).astype(np.float64),
        R=(Rp + Rn).astype(np.float64),
    ))
    print(f"  cached {cond}: {isr.sum()} region + {(~isr).sum()} boundary nodes")

OUT.write_bytes(pickle.dumps(scen))
print(f"\nwrote {len(scen)} scenarios -> {OUT}  ({OUT.stat().st_size / 1e6:.1f} MB)")
