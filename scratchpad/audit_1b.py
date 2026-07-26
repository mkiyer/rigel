"""AUDIT #1b — hop-by-hop sigma^2_transfer damping of the DISTANT (node-1914) spliced+ count as it
telescopes down the backward relay into boundary 1910, then UNDAMPED (graft, s2t=0) into exon 1909.

Re-runs the exact backward mp relay but prints, per hop along the 1914->1910 path, the graft/plain edge
kind, the applied s2t, and the surviving fraction of the origin spliced count.
"""
import os, sys, dataclasses, importlib
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from pathlib import Path
from selfsolve_diag import _scan_and_truth
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND = "gdna_gdna300_ss_0.99_nrna_present_capture_on"
_EPS = 1.0e-9
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
os.environ.pop("RIGEL_S2T_OFF", None)
dbg = {}
cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
calmod.calibrate(inp["payload"], ra, inp["strand_model"],
                 np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
st = dbg["capture"]["_uni_static"]
is_exon = np.asarray(st["is_exon"], bool); is_bnd = np.asarray(st["is_bnd"], bool)
left = np.asarray(st["left"], np.int64); right = np.asarray(st["right"], np.int64)
fp = np.asarray(st["fp"], bool)
logvar_tot = np.asarray(st["logvar_tot"], np.float64)
SP_l = np.asarray(st["SP_l"], np.float64); SP_r = np.asarray(st["SP_r"], np.float64)
mass = np.asarray(dbg["capture"]["mass_global"], np.float64)

# Walk the chain to the RIGHT of 1910 to find the path to origin 1914.
path = [1910]
cur = 1910
for _ in range(20):
    nb = right[cur]
    if nb < 0:
        break
    path.append(nb)
    if nb == 1914:
        break
    cur = nb
print("chain path 1910 -> ... -> 1914 (rightward):")
for node in path:
    print(f"   node {node:>5} is_exon={int(is_exon[node])} is_bnd={int(is_bnd[node])} "
          f"mass={mass[node]:>10,.0f}  logvar_tot={logvar_tot[node]:.4g}  SP_l={SP_l[node]:.3f} SP_r={SP_r[node]:.3f}")

# Now replay the BACKWARD relay (nbr=right, sf=0) but only track the contribution that originates at 1914
# (grafted into an exon from boundary 1914) and follow it hop-by-hop back to 1910.
# In the backward pass the sequence is order[::-1]; node i pulls from right[i]. So the origin spliced at 1914
# is deposited when the relay processes the EXON immediately to the LEFT of 1914 (whose right nbr is 1914),
# then propagates leftward hop by hop to 1910.
print("\nhop-by-hop damping of the origin-1914 spliced+ contribution as it flows leftward to 1910:")
print(f"{'dst':>6} {'src':>6} {'edge':>6} {'s2t':>10} {'in':>10} {'out(damped)':>12} {'+graft':>10} {'mp[dst]':>10}")

def _damp(p, s2t):
    return p / (1.0 + p * s2t) if p > 0.0 else 0.0

# reconstruct mp[dst] along the path from 1914 down to 1910, isolating the 1914 origin.
# find the exon left-of-1914 that grafts SP_l[1914].
seq_nodes = list(reversed(path))  # 1914 ... 1910
mp = 0.0
carried = 0.0  # the running 1914-origin contribution
for k in range(len(seq_nodes) - 1):
    src = seq_nodes[k]
    dst = seq_nodes[k + 1]  # dst = left neighbour along path (dst.right == src)
    assert right[dst] == src, (dst, src, right[dst])
    gr = bool(is_exon[dst] and is_bnd[src])
    s2t = 0.0 if gr else (logvar_tot[dst] + logvar_tot[src])
    mp_in = mp
    mp_out = _damp(mp_in, s2t)
    graft = 0.0
    if gr and SP_l[src] > _EPS:  # sf=0 -> SP_l at the source boundary
        graft = SP_l[src] / (1.0 + SP_l[src] * s2t)
    # isolate carried 1914-origin part (proportional attribution through the damping)
    scale = (mp_out / mp_in) if mp_in > _EPS else 0.0
    carried = carried * scale
    if gr and src == 1914:
        carried += graft
    mp_dst = 0.0
    if fp[dst]:
        mp_dst = mp_out + graft
    mp = mp_dst
    edge = "GRAFT" if gr else "plain"
    print(f"{dst:>6} {src:>6} {edge:>6} {s2t:>10.4f} {mp_in:>10.4f} {mp_out:>12.4f} {graft:>10.4f} {mp_dst:>10.4f}  carried(1914)={carried:.4f}")

print(f"\n=> origin-1914 contribution delivered to mp[1910] = {carried:.4f}")
print(f"   raw spliced+ at 1914 (SP_l) = {SP_l[1914]:.3f}; survived {carried/max(SP_l[1914],_EPS):.1%} of it to 1910,")
print(f"   then the GRAFT 1910->1909 applies s2t=0 (exon<-boundary), so it enters 1909 UNDAMPED.")
