"""Verify P4's SECOND-violation claim: in the forward relay, does the message INTO node j carry j's own
self-solve?  The peel at node i uses `_far[0][.][i] = op[right[i]]`; the relay then writes node i's belief,
which node right[i] consumes on the next hop IF right[i] is visited after i in `order_list`.

Also quantifies HOW MANY relay hops are actually affected (peel edges only), not just the ordering.
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
CONDS = sys.argv[1:] or ["gdna_gdna300_ss_0.99_nrna_present_capture_off"]

for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    us = cap["_uni_static"]
    order = np.asarray(cap["order"]) if "order" in cap else np.asarray(chain.order)
    li, ri = np.asarray(us["left"]), np.asarray(us["right"])
    is_bnd, is_exon = np.asarray(us["is_bnd"], bool), np.asarray(us["is_exon"], bool)
    n = li.size
    # position in the FORWARD visiting order
    pos = np.full(n, -1, np.int64)
    pos[order] = np.arange(order.size)
    okL = li >= 0
    okR = ri >= 0
    lt = bool((pos[li[okL]] < pos[np.arange(n)[okL]]).all())
    gt = bool((pos[ri[okR]] > pos[np.arange(n)[okR]]).all())
    # peel edges in the FORWARD relay: i is a boundary, its LEFT neighbour is an exon
    peelF = is_bnd & okL & is_exon[np.clip(li, 0, n - 1)]
    # of those, how many have a right neighbour that is later in the order AND itself receives from i
    selfmsg = peelF & okR & (pos[np.clip(ri, 0, n - 1)] > pos[np.arange(n)]) & (li[np.clip(ri, 0, n - 1)] == np.arange(n))
    print(f"\n{cond}   n={n}")
    print(f"  order is a permutation of 0..n-1 and monotone: {bool((order == np.sort(order)).all())}")
    print(f"  pos[left[i]]  < pos[i] for all i with a left  neighbour : {lt}")
    print(f"  pos[right[i]] > pos[i] for all i with a right neighbour : {gt}")
    print(f"  forward-relay PEEL hops (bnd with exon on the left)     : {int(peelF.sum()):,}")
    print(f"  ... of which right[i] is visited AFTER i and its source is i "
          f"(=> that node's OWN op[] rode into its own incoming message): {int(selfmsg.sum()):,} "
          f"({selfmsg.sum() / max(peelF.sum(), 1):.1%})")
