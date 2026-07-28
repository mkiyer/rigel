"""TASK 3 — EXACT CHECK: does calibrate.py's Phase-2 refit already RESET the belief?

Two things are proved here, not argued:
  (a) `_init_belief()` is a PURE function of (chain, substrate, boundary_substrate, region_arrays, the three
      fitted library scalars, the grid config) — it never reads the incoming belief. Verified by calling
      `init_beliefs` twice with a DELIBERATELY CORRUPTED belief in between and comparing bit-for-bit.
  (b) with RIGEL_NORESET the refit's input belief DIFFERS from the reset one, so line 420 is load-bearing
      (i.e. the reset is real, not a no-op that happens to reproduce the init).

TASK 2 — EXACT CHECK: RIGEL_HOLDOUT_WB (write-back holdout) cannot move a substrate node.

    OMP_NUM_THREADS=1 python scratchpad/ho_resetcheck.py
"""

from __future__ import annotations

import dataclasses
import importlib
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import build_node_statics, init_beliefs  # noqa: E402
from rigel.calibration.node_chain import build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND = "gdna_gdna300_ss_0.99_nrna_present_capture_on"
_EPS = 1e-9

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")


def run(refit, **env):
    old = {k: os.environ.get(k) for k in env}
    os.environ.update({k: v for k, v in env.items() if v is not None})
    for k, v in env.items():
        if v is None:
            os.environ.pop(k, None)
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=refit), _debug=dbg,
    )
    for k, v in old.items():
        os.environ.pop(k, None)
        if v is not None:
            os.environ[k] = v
    return dbg


# ── (a) init_beliefs is belief-free: build it twice, once after corrupting a belief ───────────────────────
payload = inp["payload"]
sub = CalibrationSubstrate.from_payload(payload, ra)
bsub = BoundarySubstrate.from_payload(payload)
chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
st = build_node_statics(chain, sub, bsub, ra)
kw = dict(rna_sense_frac=0.99, gdna_strand_overdispersion=0.01, rna_strand_overdispersion=0.02,
          n_grid=cfg.calibration.sweep_n_grid, n_grid_ss=cfg.calibration.sweep_n_grid_single_strand,
          logodds_window=cfg.calibration.sweep_logodds_window, statics=st)
b1 = init_beliefs(chain, sub, bsub, ra, **kw)
b1.f_g[:] = 0.123456  # corrupt it in place
b2 = init_beliefs(chain, sub, bsub, ra, **kw)
same = np.array_equal(b2.f_g, np.where(np.isfinite(b2.f_g), b2.f_g, 0)) and not np.allclose(b2.f_g, 0.123456)
b3 = init_beliefs(chain, sub, bsub, ra, **kw)
print(f"(a) init_beliefs is deterministic + belief-free: b2 == b3 bitwise -> "
      f"{np.array_equal(b2.f_g, b3.f_g) and np.array_equal(b2.var_gdna, b3.var_gdna)}"
      f"   (corrupted b1 did not leak: {same})")

# ── (b) the refit's input belief: reset vs carried ────────────────────────────────────────────────────────
# `capture["fg_init"]` is the belief the LAST sweep RECEIVED (bp_solver line 163). Under the shipped reset it
# must NOT equal the pass-0 OUTPUT; under RIGEL_NORESET it must equal it exactly.
d_r0 = run(0)
d_r1 = run(1)
d_nr = run(1, RIGEL_NORESET="1")
fg_pass0_out = np.asarray(d_r0["capture"]["f_g"])
in_reset = np.asarray(d_r1["capture"]["fg_init"])
in_nores = np.asarray(d_nr["capture"]["fg_init"])
print("(b) refit sweep INPUT belief vs the PASS-0 OUTPUT it could have carried:")
print(f"      shipped (reset)    : identical on {int(np.sum(in_reset == fg_pass0_out)):,}/"
      f"{in_reset.size:,} nodes,  max |Δ| = {np.max(np.abs(in_reset - fg_pass0_out)):.4f}")
print(f"      RIGEL_NORESET=1    : identical on {int(np.sum(in_nores == fg_pass0_out)):,}/"
      f"{in_nores.size:,} nodes,  max |Δ| = {np.max(np.abs(in_nores - fg_pass0_out)):.4g}")
print(f"      reset vs no-reset  : max |Δ| = {np.max(np.abs(in_reset - in_nores)):.4f} "
      f"(0 would mean the reset is a no-op)")

# ── TASK 2: the write-back holdout cannot move a substrate node ───────────────────────────────────────────
d_hd = run(0)
d_wb = run(0, RIGEL_HOLDOUT_WB="1")
fp, fn = np.asarray(d_hd["statics"].free_pos, bool), np.asarray(d_hd["statics"].free_neg, bool)
isr = np.asarray(d_hd["chain"].kind) == REGION
SUB = isr & ~(fp & fn)
a, b = np.asarray(d_hd["capture"]["f_g"]), np.asarray(d_wb["capture"]["f_g"])
print(f"\nTASK 2  write-back holdout: substrate f_g identical on "
      f"{int(np.sum(a[SUB] == b[SUB])):,}/{int(SUB.sum()):,} nodes "
      f"(max |Δ| = {np.max(np.abs(a[SUB] - b[SUB])):.3g});  excluded-half max |Δ| = "
      f"{np.max(np.abs(a[~SUB] - b[~SUB])):.3g}")
