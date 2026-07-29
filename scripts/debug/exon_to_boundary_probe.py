"""Direct observation of what the exon -> boundary edge emits (docs/CARRY_FORWARD.md §10.1).

HYPOTHESIS (from the derivation): on an exon -> splice-junction-boundary edge the code suppresses the exon's
RNA density to `SPs/_esp - absorb_p`, and a REGION has `spliced_* == 0`, so that term is <= 0. The shift
normalizer `_den = Mg + max(rho_pos,0)*erd + max(rho_neg,0)*erd` therefore collapses to `Mg`, and

        mode_g = log(Mg/_den) = log(1) = 0        i.e.  "boundary, you are PURE gDNA"

is emitted at the exon's full counting precision (plus +c_b, which pushes it ABOVE 1). If so, every
splice-junction boundary receives a high-precision all-gDNA assertion from its exon flank -- which would
explain census finding M3 (boundaries: corr ~ 0 at precision 2-3x the reference scale, r(prec,err) > 0).

This script does NOT change the solver. It reads the per-edge `_edge_modes` records the inert `_capture`
hook already emits, and reports, per edge family, the distribution of exp(mode_g) and prec_g.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/exon_to_boundary_probe.py
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
REF_PREC = 1.0 / 3.289868  # 1/Var_Beta(1/2,1/2)(log f) -- the derived "stronger than the node's own prior" scale

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
conds = sorted(
    d.name
    for d in SUITE.iterdir()
    if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and "0.50" in d.name
)

fam: dict[tuple, list] = defaultdict(list)

for cond in conds:
    capst = "ON " if "capture_on" in cond else ("VSTR" if "verystrong" in cond else "OFF")
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    cap = dbg["capture"]
    for e in cap.get("_edge_modes", []):
        if e["kind_d"] == 0:          # destination is a REGION
            if e["kind_s"] == 0:
                continue
            src = "bnd"
            dst = "exon" if e["exon_d"] else "reg"
            name = f"{src} -> {dst}"
        else:                          # destination is a BOUNDARY
            name = ("exon" if e["exon_s"] else "nonexon") + " -> bnd"
        fam[(capst, name)].append((e["mode_g"], e["prec_g"], e["c_mat"], e["md"]))

print("\nexp(mode_g) is the gDNA FRACTION the message asserts for the destination.")
print(f"prec_g is its precision; the derived reference scale is 1/Var_Beta(1/2,1/2)(log f) = {REF_PREC:.3f}")
print("-- a message above that scale is stronger than the destination's own prior.\n")
hdr = (
    f"{'capture':>7} {'edge family':>16} | {'edges':>7} | {'median':>7} {'>1.0':>6} {'>0.95':>7} |"
    f" {'med prec':>8} {'pr>ref':>7} | {'med md':>7} {'md<1':>6}"
)
print(hdr)
print("-" * len(hdr))
for capst in ("OFF", "ON ", "VSTR"):
    for k in sorted(fam):
        if k[0] != capst:
            continue
        v = np.array([r[:4] for r in fam[k]], dtype=float)
        m, p, _cb, md = v[:, 0], v[:, 1], v[:, 2], v[:, 3]
        # ONLY live messages. `amg`/`apg` are zero-initialised, so an edge whose gDNA channel did not emit
        # reads back as mode_g = 0 => exp = 1.0 -- an artifact, not an assertion.
        live = p > 1e-12
        if not live.any():
            continue
        fr = np.exp(m[live])
        pl, mdl = p[live], md[live]
        print(
            f"{capst:>7} {k[1]:>16} | {live.sum():>7} |"
            f" {np.median(fr):>7.3f} {np.mean(fr > 1.0):>6.1%} {np.mean(fr > 0.95):>7.1%} |"
            f" {np.median(pl):>8.3f}"
            f" {np.mean(pl > REF_PREC):>7.1%} | {np.median(mdl):>7.3f} {np.mean(mdl < 1.0):>6.1%}"
        )
    print()
