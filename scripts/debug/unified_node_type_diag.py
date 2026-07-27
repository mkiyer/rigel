"""Per-node-TYPE f_g diagnostic for the UNIFIED solver — localize the systematic composition bias.

Runs ONE scenario and prints, per node class (intergenic / intron / exon / boundary), the mass-weighted
mean f_g and |f_g − oracle| — so a base vs unified vs oracle comparison shows WHERE the unified solver
mis-calls the composition. Run the two arms in separate processes (the RIGEL_UNIFIED flag is read at import):

    OMP_NUM_THREADS=1 python scripts/debug/unified_node_type_diag.py base
    RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/unified_node_type_diag.py uni
    RIGEL_UNIFIED=1 RIGEL_UNIFIED_ROUTE=0 OMP_NUM_THREADS=1 python scripts/debug/unified_node_type_diag.py uni-nr

Optional 2nd arg = condition name (default a gDNA-present capture-on unstranded scenario, where the bias is
clearest). See docs/calibration/archive/UNIFIED_SOLVER_HANDOFF.md.
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
lbl = sys.argv[1] if len(sys.argv) > 1 else "arm"
cond = sys.argv[2] if len(sys.argv) > 2 else "gdna_gdna300_ss_0.50_nrna_present_capture_on"

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
calmod.calibrate(inp["payload"], ra, inp["strand_model"],
                 np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
chain = dbg["chain"]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > 1e-9, G / np.maximum(G + R, 1e-9), np.nan)
fg = np.asarray(dbg["capture"]["f_g"])
mass = np.asarray(dbg["capture"]["mass_global"])
rt, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind)
ok = np.isfinite(fo) & (mass > 1e-9)

print(f"[{lbl}] {cond}")
for name, sel in [("intergenic", (kind == 0) & (rt == 0)), ("intron", (kind == 0) & (rt == 1)),
                  ("exon", (kind == 0) & (rt == 2)), ("boundary", kind != 0)]:
    m = ok & sel
    if m.sum():
        err = np.average(np.abs(fg[m] - fo[m]), weights=mass[m])
        mfg = np.average(fg[m], weights=mass[m])
        mo = np.average(fo[m], weights=mass[m])
        print(f"  {name:<11} n={int(m.sum()):5d}  |err|={err:.4f}  meanfg={mfg:.3f}  oracle={mo:.3f}")
