"""EXON-EXON BOUNDARIES THAT ARE ALSO SPLICE JUNCTIONS — the neglected topology (owner, 2026-07-23).

A boundary's "is a splice junction" and "mature crosses contiguously" are INDEPENDENT properties, and the
mature-absorption code conflates them. Owner's example:

    TA+  (1000,2000) (10000,11000)
    TB+  (1000,2000) ( 9000,11000)

The signature changes at 10000 (TA's exon starts), so there is an exon<->exon boundary there. It is ALSO a
splice ACCEPTOR: TA splices 2000 -> 10000. So at that one boundary, simultaneously:

    * spliced mass arrives                       (TA's junction)
    * mature crosses UNSPLICED and contiguously  (TB's exon spans 9000-11000)

The c_b mature-dilution term assumes "the boundary's unspliced crossing is mature-FREE". That is true for an
exon<->intron junction and FALSE here -- so c_b over-peels.

Worse, c_b is built from the SUM of both faces:

    _spl_mass = SP[0]+SN[0]+SP[1]+SN[1] ;  _S_B = _spl_mass / (ESP[0]+ESP[1])

while the accumulator already routes spliced mass ONE-SIDED to its exon flank (node_geometry, by
junction_strand x exon bit). In the example the acceptor's exon flank is the RIGHT, so spliced_right > 0 and
spliced_left == 0 -- meaning:

    right exon -> boundary   SHOULD peel  (TA's mature terminates here)
    left  exon -> boundary   should NOT   (TB's mature continues through)

The per-face quantity `absorb_* = SPd[i]/_espd` is already directional and correct; `mature_dilution` is not.

This script measures how much of the chain is in that neglected class.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/exon_exon_junction_census.py
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
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

tot: dict[str, float] = defaultdict(float)
n_scen = 0

for cond in conds:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    chain, geom, st = dbg["chain"], dbg["geometry"], dbg["statics"]
    n_scen += 1
    kind = np.asarray(chain.kind)
    is_bnd = kind != REGION
    mass, _eff = node_global_geometry(chain, geom)

    spl_l = np.asarray(geom.spliced_pos_left) + np.asarray(geom.spliced_neg_left)
    spl_r = np.asarray(geom.spliced_pos_right) + np.asarray(geom.spliced_neg_right)
    has_spl = (spl_l + spl_r) > _EPS
    # "mature crosses contiguously" == exon on BOTH flanks, per strand (node_geometry._boundary_strand_stats)
    mat_x = np.asarray(st.mrna_active_pos, bool) | np.asarray(st.mrna_active_neg, bool)

    classes = {
        "seam            (no spliced, mature does NOT cross)": is_bnd & ~has_spl & ~mat_x,
        "classic junction(spliced,    mature does NOT cross)": is_bnd & has_spl & ~mat_x,
        "exon-exon       (no spliced, mature DOES cross)    ": is_bnd & ~has_spl & mat_x,
        "*** exon-exon + SPLICE JUNCTION (BOTH) ***         ": is_bnd & has_spl & mat_x,
    }
    for k, m in classes.items():
        tot[k + "|n"] += int(m.sum())
        tot[k + "|mass"] += float(mass[m].sum())
    tot["ALL boundaries|n"] += int(is_bnd.sum())
    tot["ALL boundaries|mass"] += float(mass[is_bnd].sum())

    # of the BOTH class: is the spliced mass one-sided (as the accumulator intends) or on both faces?
    both = is_bnd & has_spl & mat_x
    tot["both:onesided"] += int((both & ((spl_l > _EPS) ^ (spl_r > _EPS))).sum())
    tot["both:twosided"] += int((both & (spl_l > _EPS) & (spl_r > _EPS)).sum())
    # how badly does the SUMMED c_b over-state the directional peel? compare summed vs per-face
    sel = both & ((spl_l > _EPS) ^ (spl_r > _EPS))
    if sel.any():
        summed = spl_l[sel] + spl_r[sel]
        facing0 = spl_l[sel]  # the face with NO spliced, for the flank whose mature continues
        tot["both:peel_ratio_num"] += float(np.sum(summed))
        tot["both:peel_ratio_den"] += float(np.sum(np.maximum(facing0, 0.0)))

print(f"\nAveraged over {n_scen} scenarios of the ambig_dense_10mb battery.\n")
hdr = f"{'boundary class':<52} | {'count':>9} | {'% of bnd':>8} | {'mass':>12} | {'% mass':>7}"
print(hdr)
print("-" * len(hdr))
N = tot["ALL boundaries|n"]
M = tot["ALL boundaries|mass"]
for k in (
    "seam            (no spliced, mature does NOT cross)",
    "classic junction(spliced,    mature does NOT cross)",
    "exon-exon       (no spliced, mature DOES cross)    ",
    "*** exon-exon + SPLICE JUNCTION (BOTH) ***         ",
):
    n, m = tot[k + "|n"], tot[k + "|mass"]
    print(f"{k:<52} | {n/n_scen:>9,.0f} | {100*n/max(N,1):>7.1f}% | {m/n_scen:>12,.0f} | {100*m/max(M,1):>6.1f}%")
print(f"{'ALL boundaries':<52} | {N/n_scen:>9,.0f} | {100.0:>7.1f}% | {M/n_scen:>12,.0f} | {100.0:>6.1f}%")

print("\nWithin the BOTH class — is the spliced deposit one-sided, as node_geometry routes it?")
print(f"      one-sided (spliced on exactly ONE face): {tot['both:onesided']/n_scen:>9,.0f}")
print(f"      two-sided (spliced on BOTH faces)      : {tot['both:twosided']/n_scen:>9,.0f}")
print("\n  => where it is one-sided, the SUMMED c_b peels the FULL junction flux from BOTH directions,")
print("     including the flank whose mature demonstrably continues (its own face carries ZERO spliced).")
