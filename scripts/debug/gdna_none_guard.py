"""The hard gate: does a change manufacture PHANTOM gDNA in zero-gDNA libraries? Sum the gDNA-attributed mass
Σ f_g·mass_global over the gdna_none scenarios (oracle gDNA ≈ 0 ⇒ any f_g·mass is phantom), split region vs
boundary. Run twice (baseline vs RIGEL_BND_SPLICED_DENSITY=1) and compare — the fix must NOT increase it."""
from __future__ import annotations
import dataclasses
import importlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.bp_solver import REGION
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_none_"))
print(f"gdna_none scenarios ({len(conds)}): phantom gDNA mass Σ f_g·mass  [flag={'ON' if __import__('os').environ.get('RIGEL_BND_SPLICED_DENSITY')=='1' else 'off'}]\n")
tot_r = tot_b = 0.0
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    fg = np.asarray(cap["f_g"]); mass = np.asarray(cap["mass_global"]); kind = np.asarray(chain.kind)
    isr = kind == REGION
    phr = float(np.sum(fg[isr] * mass[isr])); phb = float(np.sum(fg[~isr] * mass[~isr]))
    tot_r += phr; tot_b += phb
    print(f"  {cond:>44} | region {phr:>10.0f} | boundary {phb:>9.0f}")
print(f"\n  {'TOTAL':>44} | region {tot_r:>10.0f} | boundary {tot_b:>9.0f} | all {tot_r+tot_b:>10.0f}")
