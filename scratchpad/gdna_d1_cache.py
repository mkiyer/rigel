"""DISSECTION — the paired pass-0 / re-solved cache AT THE W4 HEAD, plus the fitted landscape itself.

`gdna_refit_cache.pkl` was built when prior #1 was the retiring delta-pin. Everything below needs the
CURRENT object, and needs it per node, because the questions are about MOVEMENT:

    which direction does the re-solve move each node, and is that direction right for what the node IS?

Stores, per condition: the pass-0 and re-solved beliefs, the node's mass/eff, the structural masks, and the
fitted `GdnaLandscape` (so the prior can be interrogated along each node's own ray without refitting).

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d1_cache.py [--out <pkl>]
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
from flagship_interrogate import _oracle_per_node  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.calibrate import _fit_gdna_hyperprior  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
ap.add_argument("--out", default="/Users/mkiyer/proj/rigel/scratchpad/gdna_dissect_cache.pkl")
args = ap.parse_args()

suite = Path(args.suite)
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
cfg = PipelineConfig()
work, cache = Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache"
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))
rtype_all = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)


def _group(cond):
    p = cond.split("_")
    return ({"on": "ON", "off": "OFF", "verystrong": "VSTRONG"}.get(p[7], p[7].upper()),
            p[1], p[3], "none" if p[5] == "none" else "nrna")


out = []
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    beliefs = {}
    for iters in (0, 1):
        dbg = {}
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=iters), _debug=dbg)
        beliefs[iters] = dbg
    d0, d1 = beliefs[0], beliefs[1]
    chain, cap = d0["chain"], d0["capture"]
    mg = np.asarray(cap["mass_global"], np.float64)
    eg = np.asarray(cap["eff_global"], np.float64)
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    isr = np.asarray(chain.kind) == REGION
    idx = np.asarray(chain.ref_idx, np.int64)
    # refit the landscape from the PASS-0 belief — exactly the object the re-solve consumed
    land = _fit_gdna_hyperprior(chain, d0["belief"], d0["statics"], ra, mg, eg,
                                strength=cfg.calibration.gdna_prior_strength)
    out.append(dict(
        cond=cond, group=_group(cond), is_region=isr,
        ntype=np.where(isr, rtype_all[np.clip(idx, 0, len(rtype_all) - 1)], 3).astype(np.int8),
        fp=np.asarray(cap["free_pos"], bool), fn=np.asarray(cap["free_neg"], bool),
        mass=mg, eff=eg,
        f0=np.asarray(d0["belief"].f_g, np.float64), var0=np.asarray(d0["belief"].var_gdna, np.float64),
        f1=np.asarray(d1["capture"]["f_g"], np.float64),
        G=(Gp + Gn).astype(np.float64), R=(Rp + Rn).astype(np.float64),
        landscape=land,
    ))
    fp_mass = float(np.sum(out[-1]["f1"] * mg)) if out[-1]["group"][1] == "none" else np.nan
    print(f"  {cond:52s} mean|df_g|={np.mean(np.abs(out[-1]['f1'] - out[-1]['f0'])):.4f} "
          f"FP mass={fp_mass:.0f}", flush=True)

Path(args.out).write_bytes(pickle.dumps(out))
print(f"\nwrote {len(out)} conditions -> {args.out}")
