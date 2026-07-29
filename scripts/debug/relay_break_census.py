"""B1 DESIGN — how badly does a zero-mass node break the forward-backward relay?

Two distinct defects live at a zero-mass node, and they must be separated:

  D1  EMISSION.  `_scan` builds the message density as `rho_c = f_c * sm / E_c` with `sm` the SOURCE's own
      facing mass, and gates on `emit_g = sm > _EPS`. A node with no mass therefore has no density to report
      and sends NOTHING -- even though the density it *should* report is its NEIGHBOUR's (gDNA is genomically
      continuous; nascent is continuous within a transcript). The chain dies there.

  D2  FALSE CERTAINTY.  `solvable = (fp|fn) & mass>0`, `locked = ~solvable`, and a locked node is reseeded to
      `lam = logit(f_g_init) = +L` (all-gDNA) at `lvar = 0` (CERTAIN). That is correct for a STRUCTURAL
      all-gDNA node (intergenic / a no-RNA-crossing seam: no admissible RNA strand, so its unspliced mass is
      necessarily gDNA). It is WRONG for a node that merely has no data: it admits RNA strands, the init
      itself says `var_gdna = inf` ("no information"), and the sweep overrides that to 0 ("certain").
      So the two are conflated by one predicate.

This script sizes both, and -- the design question -- measures how much of the chain the breaks ISOLATE:
segments of consecutive live nodes separated by breaks, and whether each segment still reaches a structural
gDNA ANCHOR (an intergenic node with mass). A segment with no anchor has no absolute gDNA reference at all.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/relay_break_census.py
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

from rigel.calibration.bp_solver import REGION, node_global_geometry  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
INTERGENIC = 0

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

agg: dict[str, float] = defaultdict(float)
seg_sizes: list[int] = []
n_scen = 0

for cond in conds:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    chain, st, geom = dbg["chain"], dbg["statics"], dbg["geometry"]
    n_scen += 1
    kind = np.asarray(chain.kind)
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    mu = np.asarray(st.mass_unspliced, float)
    fp = np.asarray(st.free_pos, bool)
    fn = np.asarray(st.free_neg, bool)
    rtype, _ = _node_region_type(chain, ra)
    mass_g, eff_g = node_global_geometry(chain, geom)

    has_rna = fp | fn
    live = mu > _EPS                      # has an observation
    struct_g = ~has_rna                   # no admissible RNA strand -> structurally all-gDNA
    # D2: nodes the code marks "locked/certain" that are NOT structurally gDNA
    false_certain = ~live & has_rna
    agg["nodes"] += kind.shape[0]
    agg["live"] += live.sum()
    agg["dead (mass==0)"] += (~live).sum()
    agg["  dead & STRUCTURAL gdna (correctly certain)"] += (~live & struct_g).sum()
    agg["  dead & has RNA strands (FALSELY certain)"] += false_certain.sum()
    agg["    ... of which REGION"] += (false_certain & (kind == REGION)).sum()
    agg["    ... of which BOUNDARY"] += (false_certain & (kind != REGION)).sum()
    # D1: dead interior nodes with a LIVE neighbour -> the relay would have carried something
    interior = (left >= 0) & (right >= 0)
    nbr_live = np.zeros_like(live)
    for nb in (left, right):
        ok = nb >= 0
        nbr_live |= np.where(ok, live[np.clip(nb, 0, None)], False)
    agg["dead interior w/ live neighbour (relay BREAKS)"] += (~live & interior & nbr_live).sum()
    agg["  ... eff_g < 1 (cannot contain a fragment)"] += (
        ~live & interior & nbr_live & (eff_g < 1.0)
    ).sum()

    # ---- SEGMENTATION: split each reference's chain at dead nodes ----
    anchor = live & struct_g & (rtype == INTERGENIC)   # a real, measured pure-gDNA reference
    order = np.asarray(chain.order)
    seg_id = np.full(kind.shape[0], -1, np.int64)
    cur = -1
    prev_ref = None
    for n in order:
        r = int(chain.ref_idx[n]) if False else None  # ref grouping is via left==-1 below
        if left[n] < 0:      # start of a reference path
            cur += 1
            new_seg = True
        else:
            new_seg = False
        if not live[n]:
            cur += 1         # a dead node terminates the current segment
            continue
        if new_seg or seg_id[left[n]] < 0:
            cur += 1
        seg_id[n] = cur
    segs = defaultdict(list)
    for n in np.where(seg_id >= 0)[0]:
        segs[int(seg_id[n])].append(n)
    n_anch = n_noanch = 0
    m_anch = m_noanch = 0.0
    for _s, nodes in segs.items():
        nn = np.array(nodes)
        seg_sizes.append(len(nn))
        if anchor[nn].any():
            n_anch += 1
            m_anch += float(mass_g[nn].sum())
        else:
            n_noanch += 1
            m_noanch += float(mass_g[nn].sum())
    agg["segments (anchored)"] += n_anch
    agg["segments (NO anchor)"] += n_noanch
    agg["mass in anchored segments"] += m_anch
    agg["mass in ANCHOR-LESS segments"] += m_noanch

print(f"\nAveraged over {n_scen} scenarios of the ambig_dense_10mb battery.\n")
print("--- node inventory ---")
N = agg["nodes"] / n_scen
for k in (
    "nodes", "live", "dead (mass==0)",
    "  dead & STRUCTURAL gdna (correctly certain)",
    "  dead & has RNA strands (FALSELY certain)",
    "    ... of which REGION", "    ... of which BOUNDARY",
    "dead interior w/ live neighbour (relay BREAKS)",
    "  ... eff_g < 1 (cannot contain a fragment)",
):
    print(f"   {k:<48} {agg[k]/n_scen:>9,.0f}  ({100*agg[k]/agg['nodes']:>5.1f}% of nodes)")

print("\n--- chain segmentation (segments = maximal runs of LIVE nodes, split at every dead node) ---")
sa, sn = agg["segments (anchored)"] / n_scen, agg["segments (NO anchor)"] / n_scen
ma, mn = agg["mass in anchored segments"], agg["mass in ANCHOR-LESS segments"]
print(f"   segments reaching a structural gDNA anchor : {sa:>8,.0f}   mass {100*ma/max(ma+mn,1):>5.1f}%")
print(f"   segments with NO anchor                    : {sn:>8,.0f}   mass {100*mn/max(ma+mn,1):>5.1f}%")
ss = np.array(seg_sizes)
print(f"   segment size: median {np.median(ss):.0f}  mean {ss.mean():.1f}  p90 {np.percentile(ss,90):.0f}  max {ss.max()}")
print(f"   singleton segments (1 node, fully isolated): {100*np.mean(ss==1):.1f}%")
