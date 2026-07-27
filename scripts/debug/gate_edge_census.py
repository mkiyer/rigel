"""Static edge census for the mature-crossing gate — Phase 1 instrument (no production change).

Counts, per directed boundary↔region edge, the quantities the gate's correctness rests on — all STATIC
(functions of signature + geometry only, never of the running belief), so the census is EXACT, not sampled:

  * ``emit``            — edges where the ``+``/``−`` RNA message fires today (``emit_p``/``emit_n``).
  * ``gated``           — edges the asymmetric gate ``send_s = mrna_active_s[dst] or not mrna_active_s[src]``
                          would silence (``n_nasc → 0``).
  * ``gated_live_nmat`` — gated edges that ALSO carry a live spliced measurement ``n_mat = SPs[src] > 0``.
                          **Must be 0** — that is why the gate touches only ``n_nasc`` and leaves the ``n_mat``
                          measurement path intact (``docs/calibration/archive/mature_crossing_gate.md`` §3d). If this
                          is > 0 the gate as written would silence a real measurement — STOP and re-plan.
  * ``absorb`` / ``absorb_gated`` / ``absorb_survives`` — edges where the dst-face mature absorption
                          ``− SPd[dst]/ESPd[dst]`` fires; how many the gate subsumes; and how many still fire
                          (exon→exon-junction). **absorb_survives > 0 means the absorption term must stay**
                          (§3f) — it is not made dead by the gate.

Mirrors ``bp_solver.node_sweep``'s edge geometry exactly (``_adjacent`` scan direction, the per-face
``MS``/``SP``/``SN`` selection). Reuses ``ablate_replay.build`` so the chain/statics/geometry are byte-identical
to what the solver sees.

    OMP_NUM_THREADS=1 python scripts/debug/gate_edge_census.py [--suite DIR] [--conditions a,b]

Expected on ``ambig_dense_10mb`` (summed over all 7 selfsolve conditions):
    emit=37485  gated=5543  gated_live_nmat=0  absorb=5172  absorb_gated=4851  absorb_survives=321
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from ablate_replay import build

from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1.0e-9


def _scan_census(seq, nbr, sf, df, masks):
    """Count the edge classes for ONE scan direction. A directed edge (src → dst) presents the src's facing
    face ``sf`` to the dst's ``df``; the forward scan uses the LEFT neighbour (src presents its right face),
    the backward scan the RIGHT neighbour — matching ``node_sweep``'s two ``_scan`` passes exactly."""
    fp, fn, mrp, mrn, MS, SP, SN = masks
    MSs, MSd = MS[sf], MS[df]
    SPs, SNs = SP[sf], SN[sf]
    SPd, SNd = SP[df], SN[df]
    t = dict(emit=0, gated=0, gated_live_nmat=0, absorb=0, absorb_gated=0, absorb_survives=0)
    for i in seq:
        lsrc = nbr[i]
        if lsrc < 0:
            continue
        sm = MSs[lsrc]
        for free, mr, SPsrc, SPdst in (
            (fp, mrp, SPs, SPd),
            (fn, mrn, SNs, SNd),
        ):
            emit = bool(free[lsrc] and free[i] and (sm > _EPS or SPsrc[lsrc] > _EPS))
            if not emit:
                continue
            t["emit"] += 1
            send = bool(mr[i] or not mr[lsrc])  # asymmetric gate
            if not send:
                t["gated"] += 1
                if SPsrc[lsrc] > _EPS:
                    t["gated_live_nmat"] += 1
            if SPdst[i] > _EPS:  # dst-face mature absorption fires
                t["absorb"] += 1
                if not send:
                    t["absorb_gated"] += 1
                else:
                    t["absorb_survives"] += 1
    return t


def _census(chain, st, geom):
    fp = np.asarray(st.free_pos, bool)
    fn = np.asarray(st.free_neg, bool)
    mrp = np.asarray(st.mrna_active_pos, bool)
    mrn = np.asarray(st.mrna_active_neg, bool)
    MS = (np.asarray(geom.mass_left, float), np.asarray(geom.mass_right, float))
    SP = (np.asarray(geom.spliced_pos_left, float), np.asarray(geom.spliced_pos_right, float))
    SN = (np.asarray(geom.spliced_neg_left, float), np.asarray(geom.spliced_neg_right, float))
    masks = (fp, fn, mrp, mrn, MS, SP, SN)
    order = [int(x) for x in np.asarray(chain.order)]
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    fwd = _scan_census(order, left, 1, 0, masks)  # forward: src=left, presents right face
    bwd = _scan_census(order[::-1], right, 0, 1, masks)  # backward: src=right, presents left face
    return {k: fwd[k] + bwd[k] for k in fwd}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))

    from selfsolve_diag import _scan_and_truth  # noqa: PLC0415  (shared cache loader)

    rows = []
    total = dict(emit=0, gated=0, gated_live_nmat=0, absorb=0, absorb_gated=0, absorb_survives=0)
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, Path("/tmp/rigel_selfsolve"), cache)
        ra, sub, bsub, chain, st, geom, kappa, belief, cap, cc = build(inp, index, cfg)
        t = _census(chain, st, geom)
        rows.append(dict(condition=c, **t))
        for k in total:
            total[k] += t[k]
    rows.append(dict(condition="TOTAL", **total))
    df = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(df.to_string(index=False))
    assert total["gated_live_nmat"] == 0, (
        f"gated_live_nmat = {total['gated_live_nmat']} (expected 0) — the gate would silence a live "
        "measurement; §3d is violated, STOP."
    )
    print("\nOK: gated_live_nmat == 0 (the gate touches only n_nasc; the n_mat measurement path is intact).")


if __name__ == "__main__":
    main()
