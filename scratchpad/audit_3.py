"""audit_3 — ENRICHMENT LEAK into PRECISION.

Trace cm_p=26.45 at node 1909 to its ORIGINATING discrete spliced counts, and test the invariant
precision_claimed <= originating_count. Decompose the measurement stream cm_p = amp + bmp into
(a) propagated measurement _dv(mp[src]) and (b) the graft _spc = SP/(1+SP*s2t), and measure the
reframe r at each grafting edge.
"""
import os, sys, dataclasses, importlib
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from pathlib import Path
from flagship_interrogate import _oracle_per_node
from selfsolve_diag import _scan_and_truth
from rigel.calibration.bp_solver import REGION
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.node_geometry import _node_region_type
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND = "gdna_gdna300_ss_0.99_nrna_present_capture_on"
index = TranscriptIndex.load(str(SUITE / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
os.environ.pop("RIGEL_S2T_OFF", None)
dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                 np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
chain = dbg["chain"]; cap = dbg["capture"]; uni = cap["_uni"][-1]; st = cap["_uni_static"]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > 1e-9, G / np.maximum(G + R, 1e-9), np.nan)

left = st["left"]; right = st["right"]; is_bnd = st["is_bnd"]; is_exon = st["is_exon"]
M = st["M"]; E_g = st["E_g"]; E_r = st["E_r"]
SP_l = st["SP_l"]; SP_r = st["SP_r"]; SN_l = st["SN_l"]; SN_r = st["SN_r"]
logvar = st["logvar_tot"]
fwd_mp = st["fwd_mp"]; bwd_mp = st["bwd_mp"]
fwd_mn = st["fwd_mn"]; bwd_mn = st["bwd_mn"]
rho_l0 = st["rho_l0"]; rho_r0 = st["rho_r0"]
fp = st["fp"]; fn = st["fn"]
EPS = 1e-9


def transfer(dst, src, graft):
    return 0.0 if graft else (logvar[dst] + logvar[src])


def dv(psrc, s2t):
    return 1.0 / (1.0 / max(psrc, EPS) + s2t) if psrc > 0.0 else 0.0


def decompose(i):
    """Reconstruct cm_p at node i = amp(left transport) + bmp(right transport)."""
    print(f"\n=== node {i}  is_exon={is_exon[i]} is_bnd={is_bnd[i]} M={M[i]:,.0f} "
          f"oracle_fg={fo[i]:.3f} logvar_tot={logvar[i]:.4g}")
    print(f"    cm_p(reported)={uni['cm_p'][i]:.4g}  mo_p->f_pos={np.exp(uni['mo_p'][i]):.3f}")
    # LEFT message: src=left[i], dst face 0, src face 1 ; uses fwd stream ; dst_face_v=rho_l0, src_face_v=rho_r0
    # RIGHT message: src=right[i], dst face 1, src face 0 ; uses bwd stream ; dst_face_v=rho_r0, src_face_v=rho_l0
    total = 0.0
    for tag, src, mp_stream, dstf, srcf, SP_src_face in (
        ("L", int(left[i]), fwd_mp, rho_l0[i], rho_r0, SP_r),
        ("R", int(right[i]), bwd_mp, rho_r0[i], rho_l0, SP_l),
    ):
        if src < 0:
            continue
        graft = bool(is_exon[i] and is_bnd[src])
        s2t = transfer(i, src, graft)
        # reframe r for this edge
        rho_src = srcf[src]
        r = (dstf / rho_src) if (rho_src > EPS and dstf > EPS) else 1.0
        prop = dv(mp_stream[src], s2t)          # propagated measurement _dv(mp[src])
        sp = SP_src_face[src] if graft else 0.0  # graft spliced MASS (used as the count)
        s2t_spl = s2t if np.isfinite(s2t) else 0.0
        spc = sp / (1.0 + sp * s2t_spl) if sp > EPS else 0.0
        tmp = prop + spc  # matches _transport tmp (before fp_a gate)
        gate = fp[i]
        tmp_g = tmp if gate else 0.0
        total += tmp_g
        print(f"  {tag}: src={src} is_bnd={is_bnd[src]} graft={graft} s2t={s2t:.4g}  r={r:,.3g}")
        print(f"       mp_stream[src]={mp_stream[src]:.4g} -> _dv(prop)={prop:.4g}")
        print(f"       SP_srcface={sp:.4g} -> _spc(graft)={spc:.4g}  (s2t={s2t_spl:.3g})")
        print(f"       tmp={tmp:.4g} (fp_gate={gate}) -> contrib={tmp_g:.4g}")
    print(f"    reconstructed cm_p = {total:.4g}   (reported {uni['cm_p'][i]:.4g})")
    return total


for i in (1909,):
    decompose(i)

# --- follow the measurement stream back: what discrete spliced counts feed 1908 and 1910? ---
print("\n\n########## PROVENANCE OF THE MEASUREMENT STREAM AT THE NEIGHBOURS ##########")
for i in (int(left[1909]), int(right[1909])):
    print(f"\nnode {i} is_bnd={is_bnd[i]} is_exon={is_exon[i]} M={M[i]:,.0f} oracle_fg={fo[i]:.3f}")
    print(f"   fwd_mp={fwd_mp[i]:.4g}  bwd_mp={bwd_mp[i]:.4g}")
    print(f"   SP_l={SP_l[i]:.4g} SP_r={SP_r[i]:.4g}  (spliced pos MASS per face)")
    print(f"   logvar_tot={logvar[i]:.4g}")

# --- GLOBAL invariant check: for every node, is cm_p > total spliced-pos count in its 2-hop nbhd? ---
print("\n\n########## HONEST-PRECISION INVARIANT SCAN (cm_p vs local spliced count) ##########")
cm_p = uni["cm_p"]
n = len(cm_p)
SP_tot = SP_l + SP_r  # per-node spliced pos MASS (~count)
worst = []
for i in range(n):
    if not fp[i] or cm_p[i] <= 0:
        continue
    # local originating spliced count = this node + its immediate neighbours
    nb = [i]
    if left[i] >= 0: nb.append(int(left[i]))
    if right[i] >= 0: nb.append(int(right[i]))
    local_sp = sum(SP_tot[j] for j in nb)
    if cm_p[i] > local_sp + 1e-6:
        worst.append((cm_p[i] - local_sp, i, cm_p[i], local_sp))
worst.sort(reverse=True)
print(f"nodes where cm_p > (self+nbr) spliced-pos count: {len(worst)}")
for excess, i, c, s in worst[:15]:
    print(f"  node {i}: cm_p={c:.3g} local_spliced_pos={s:.3g} excess={excess:.3g} "
          f"M={M[i]:,.0f} oracle_fg={fo[i]:.3f} is_exon={is_exon[i]}")
