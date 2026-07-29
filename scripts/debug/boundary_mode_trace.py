"""Confirm the boundary-message saturation MECHANISM: for a given boundary node, dump the per-edge density-mode
arithmetic across scenarios. implied f_g = clip(max(ρ_g,1/E_g^dst)·E_g^dst / md, 0, 1). If the imputed gDNA mass
ρ_g·E_g^dst routinely EXCEEDS the dst observed crossing total md, the message clamps to f_g=1 — a saturated,
source-decoupled all-gDNA message. Uses the inert `_edge_modes` capture (per-edge rho_g/egd/md/sm)."""
from __future__ import annotations
import dataclasses
import importlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
from flagship_interrogate import _oracle_per_node
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
NODE = int(sys.argv[1]) if len(sys.argv) > 1 else 3236
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and "0.50" in d.name)

print(f"node {NODE} — per-scenario incoming-edge density-mode arithmetic\n")
print(f"{'scenario':>40} | {'bnd_or':>6} | {'src':>5} | {'rho_g':>8} {'egd':>7} {'imp=rho*egd':>11} "
      f"{'md':>8} | {'imp/md':>7} {'implied':>7}")
for cond in sorted(conds):
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    edges = [e for e in cap.get("_edge_modes", []) if e["dst"] == NODE]
    tag = cond.replace("gdna_", "").replace("_ss_0.50_nrna", "").replace("_capture", "/cap")
    for e in edges:
        imp = max(e["rho_g"], 1.0 / max(e["egd"], _EPS)) * e["egd"]
        ratio = imp / max(e["md"], _EPS)
        implied = min(max(ratio, 0.0), 1.0)
        print(f"{tag:>40} | {fo[NODE]:>6.2f} | {e['src']:>5} | {e['rho_g']:>8.3f} {e['egd']:>7.1f} "
              f"{imp:>11.1f} {e['md']:>8.1f} | {ratio:>7.2f} {implied:>7.2f}")
