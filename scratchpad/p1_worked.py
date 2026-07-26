"""P1 — ONE node, in plain units. Counts, lengths, densities. No solver jargon.

Dumps the single largest confidently-wrong exon on unstranded x capture-OFF and everything that reaches it,
in fragments and base pairs, so the conundrum can be read off directly.

    OMP_NUM_THREADS=1 python scratchpad/p1_worked.py
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

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.50_nrna_present_capture_off"
NODES = [int(x) for x in sys.argv[2:]] or [2651, 2441, 1405]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
calmod.calibrate(
    inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
    np.asarray(inp["rna_fl_pmf"]),
    dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
)
chain, cap, st, geo = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
us, uni = cap["_uni_static"], cap["_uni"][-1]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
T = Gp + Gn + Rp + Rn
fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
kind = np.asarray(chain.kind)
rt, _ = _node_region_type(chain, ra)
n = kind.shape[0]
CLS = {0: "intergenic", 1: "intron", 2: "exon"}
cls = np.array([CLS.get(int(rt[j]), "?") if kind[j] == REGION else "boundary" for j in range(n)])
li, ri = us["left"], us["right"]
M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
n_u = (np.asarray(us["n_unspl_l"], float), np.asarray(us["n_unspl_r"], float))
SPf = (us["SP_l"], us["SP_r"])
SNf = (us["SN_l"], us["SN_r"])
NSP = (np.asarray(us["spl_n_pos_l"]) + np.asarray(us["spl_n_neg_l"]),
       np.asarray(us["spl_n_pos_r"]) + np.asarray(us["spl_n_neg_r"]))
ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
ms = np.asarray(st.mass_spliced, float)
mu = np.asarray(st.mass_unspliced, float)

for node in NODES:
    print(f"\n{'=' * 100}")
    print(f"NODE {node}   ({cls[node]}, {COND[5:]})")
    print(f"{'=' * 100}")
    nn = float(n_u[0][node] + (n_u[1][node] if cls[node] == "boundary" else 0.0))
    print("  WHAT THIS NODE ITSELF SEES")
    print(f"    unspliced fragments (need deconvolution) : {nn:>12,.0f}  over {E_g[node]:>9,.0f} bp "
          f"(gDNA eff-len) / {E_r[node]:>9,.0f} bp (RNA eff-len)")
    print(f"    spliced fragments INSIDE it (pure RNA)   : {ms[node]:>12,.0f}")
    print(f"    unspliced mass                           : {mu[node]:>12,.1f}")
    print(f"    TRUTH: {100 * (1 - fo[node]):.1f}% of its unspliced mass is RNA "
          f"({100 * fo[node]:.1f}% gDNA)")
    print(f"    SOLVER: {100 * (1 - cap['f_g'][node]):.1f}% RNA, and it states a std-dev of "
          f"{np.sqrt(cap['var_g'][node]):.4f} on that -- so it is "
          f"{abs(cap['f_g'][node] - fo[node]) / max(np.sqrt(cap['var_g'][node]), 1e-9):.1f} "
          f"std-devs away from the truth")
    print(f"    (with no messages at all it would say {100 * (1 - cap['fg_loc'][node]):.1f}% RNA "
          f"-- the uninformed default)")

    print("\n  WHAT ARRIVES FROM ITS TWO NEIGHBOURING BOUNDARIES")
    print(f"    {'':<6}{'nbr':>6}{'spliced frags':>15}{'over eff-len':>14}{'= density':>12}"
          f"{'   vs this exon: RNA density if 100% RNA':>42}")
    for side, nb, face in (("left", li[node], 1), ("right", ri[node], 0)):
        if nb < 0:
            continue
        cnt = float(NSP[face][nb])
        mass_s = float(SPf[face][nb] + SNf[face][nb])
        e = float(ESP[face][nb])
        dens = mass_s / max(e, _EPS)
        print(f"    {side:<6}{nb:>6}{cnt:>15,.0f}{e:>14,.0f}{dens:>12.4f}"
              f"{mu[node] / max(E_r[node], _EPS):>42.4f}")
    print(f"\n    the exon spans {E_r[node]:,.0f} bp of RNA effective length; each junction is measured over "
          f"~{np.mean([ESP[1][li[node]], ESP[0][ri[node]]]):,.0f} bp")
    print(f"    ratio: the spliced measurement has to speak for an interval "
          f"{E_r[node] / max(np.mean([ESP[1][li[node]], ESP[0][ri[node]]]), _EPS):.0f}x longer than the one "
          f"it was measured over")

    # what the two channels actually CLAIM, in the node's own units
    og_t = fo[node] * mu[node] / max(E_g[node], _EPS)          # oracle gDNA density
    or_t = (1 - fo[node]) * mu[node] / max(E_r[node], _EPS)    # oracle RNA density
    cg, cr = uni["cg"][node], uni["cp"][node] + uni["cn"][node]
    print("\n  WHAT THE TWO CHANNELS CLAIM, in fragments-per-base (the node's own frame)")
    print(f"    {'':<26}{'claimed':>12}{'TRUTH':>12}{'ratio':>10}")
    print(f"    {'gDNA density':<26}{cg:>12.4f}{og_t:>12.4f}{cg / max(og_t, _EPS):>9.1f}x")
    print(f"    {'RNA density':<26}{cr:>12.4f}{or_t:>12.4f}{cr / max(or_t, _EPS):>9.1f}x")
    print(f"    those imply f_gDNA = {cg * E_g[node] / max(mu[node], _EPS):>6.3f} against a truth of "
          f"{fo[node]:.3f}")
    print(f"    and they account for {(cg * E_g[node] + cr * E_r[node]) / max(mu[node], _EPS):>5.2f}x "
          f"the fragments this node actually observed (must be 1.00 -- the conservation check)")

    print("\n  THE CONFIDENCE EACH CHANNEL BRINGS TO THIS NODE (higher = more sure)")
    print(f"    spliced-RNA measurement : {uni['cm_p'][node] + uni['cm_n'][node]:>10.1f}")
    print(f"    gDNA-level measurement  : {uni['cm_g'][node]:>10.1f}")
    print(f"    composition (ratio) msg : {uni['c_tau'][node]:>10.1f}")
    print(f"    the node's own evidence : {us['tau_own'][node]:>10.1f}   "
          f"(0 = unstranded, it has none of its own)")
