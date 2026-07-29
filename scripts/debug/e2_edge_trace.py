"""Trace WHERE the gDNA message rides the shift vs the density mode, and hence where E2's reframe can fire.

E2 reframes the gDNA message ONLY on density-mode edges (``not use_shift_g``); on shift edges the message
keeps the source frame. So if the critical enrichment cliff (seam -> exon, the factory hop INTO the enriched
exon) is classified shift-mode, E2 structurally cannot correct it. This reconstructs ``use_shift_g`` offline
from the calibrate capture (which records use_shift, kind_s/d, exon_s/d) and tabulates edges by type.

Usage: OMP_NUM_THREADS=1 python scripts/debug/e2_edge_trace.py [condition]
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
REGION = 0
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.50_nrna_present_capture_verystrong"

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
calmod.calibrate(inp["payload"], ra, inp["strand_model"],
                 np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
cap = dbg["capture"]
edges = cap.get("_edge_modes", [])

def _type(kind, exon):
    if kind != REGION:
        return "BND"
    return "exon" if exon else "reg"

tab = Counter()
for e in edges:
    ex_s, ex_d = e["exon_s"], e["exon_d"]
    bnd_i = e["kind_d"] != REGION
    bnd_src = e["kind_s"] != REGION
    use_shift = e["use_shift"]
    use_shift_g = use_shift or (bnd_i and ex_s) or (ex_d and bnd_src)
    st = _type(e["kind_s"], ex_s)
    dt = _type(e["kind_d"], ex_d)
    mode = "SHIFT (E2 skips msg)" if use_shift_g else "density (E2 reframes)"
    tab[(f"{st:>4} -> {dt:<4}", mode)] += 1

print(f"\ncondition: {cond}")
print(f"total gDNA edges: {len(edges)}\n")
print(f"  {'edge type':<14} {'mode':<24} {'count':>7}")
for (etype, mode), n in sorted(tab.items()):
    print(f"  {etype:<14} {mode:<24} {n:>7}")

shift = sum(n for (e, m), n in tab.items() if "SHIFT" in m)
dens = sum(n for (e, m), n in tab.items() if "density" in m)
print(f"\n  SHIFT (E2 message reframe SKIPPED): {shift}  ({100*shift/max(len(edges),1):.1f}%)")
print(f"  density (E2 message reframe FIRES): {dens}  ({100*dens/max(len(edges),1):.1f}%)")
print("\n  --> the edges INTO an exon (where the enrichment cliff matters most):")
for (etype, mode), n in sorted(tab.items()):
    if "-> exon" in etype:
        print(f"      {etype:<14} {mode:<24} {n:>7}")
