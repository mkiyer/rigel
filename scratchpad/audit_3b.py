"""audit_3b — trace the measurement-stream ACCUMULATION that produces cm_p=4310 at node 2933
with ZERO local spliced count. Walk the relay hop-by-hop, recording the per-hop s2t damping and
reframe r, and compare the accumulated cm_p to the TOTAL spliced-pos count anywhere in the chain.
"""
import os, sys, dataclasses, importlib
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from pathlib import Path
from flagship_interrogate import _oracle_per_node
from selfsolve_diag import _scan_and_truth
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.region_arrays import RegionArrays
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
M = st["M"]; logvar = st["logvar_tot"]; order = st["order"]
SP_l = st["SP_l"]; SP_r = st["SP_r"]
fwd_mp = st["fwd_mp"]; bwd_mp = st["bwd_mp"]
fp = st["fp"]
EPS = 1e-9
SP_tot = SP_l + SP_r

print(f"TOTAL spliced-pos MASS across ALL nodes = {SP_tot.sum():,.1f}")
print(f"MAX single-node spliced-pos MASS       = {SP_tot.max():,.1f} (node {SP_tot.argmax()})")

# The forward relay recurrence for the measurement stream mp (RNA+):
#   mp[i] = mp_own[i](=0) + _dv(mp[src], s2t) + graft_spc
# where src = left[i] in fwd order, graft only if is_exon[i] & is_bnd[src], s2t=0 on graft else logvar[i]+logvar[src].
# Reconstruct the fwd relay to see accumulation, then walk 2933's contributing chain.


def dv(psrc, s2t):
    return 1.0 / (1.0 / max(psrc, EPS) + s2t) if psrc > 0.0 else 0.0


def relay(seq, nbr, src_is_left):
    mp = np.zeros(len(M))
    contrib_src = [None] * len(M)  # (src, propagated_in, graft_added, s2t, r_proxy)
    for i in seq:
        s = nbr[i]
        if s < 0:
            continue
        graft = bool(is_exon[i] and is_bnd[s])
        s2t = 0.0 if graft else (logvar[i] + logvar[s])
        prop = dv(mp[s], s2t)
        # graft spliced face: fwd uses src face 1 (SP_r), bwd uses src face 0 (SP_l)
        spface = (SP_r if src_is_left else SP_l)[s] if graft else 0.0
        spc = spface / (1.0 + spface * (s2t if np.isfinite(s2t) else 0.0)) if spface > EPS else 0.0
        tmp = prop + spc
        if not fp[i]:
            tmp = 0.0
        mp[i] = 0.0 + tmp
        contrib_src[i] = (int(s), prop, spc, s2t, mp[s])
    return mp, contrib_src


order_list = [int(x) for x in order]
mp_f, csrc_f = relay(order_list, left, True)
mp_b, csrc_b = relay(order_list[::-1], right, False)
print(f"\nreconstructed fwd_mp[2933]={mp_f[2933]:.4g} (captured {fwd_mp[2933]:.4g})")
print(f"reconstructed bwd_mp[2933]={mp_b[2933]:.4g} (captured {bwd_mp[2933]:.4g})")


def walk_back(node, mp, csrc, nbr_arr, tag):
    print(f"\n--- walk back {tag} from node {node} (mp={mp[node]:.4g}) ---")
    hops = 0
    cur = node
    total_graft = 0.0
    while cur is not None and hops < 40:
        c = csrc[cur]
        if c is None:
            print(f"  node {cur}: seed (no src)")
            break
        s, prop, spc, s2t, mp_s = c
        total_graft += spc
        print(f"  node {cur}: is_exon={is_exon[cur]} M={M[cur]:,.0f} <- src {s}: "
              f"prop_in={prop:.4g} graft_spc={spc:.4g} s2t={s2t:.3g} (mp[src]={mp_s:.4g})")
        if spc > 0 or prop == 0:
            pass
        cur = s if mp_s > EPS or spc > 0 else None
        hops += 1
    print(f"  >>> total grafted spliced-pos along this chain = {total_graft:.4g}")


# node 2933 receives L (fwd) and R (bwd). Which dominates?
print(f"\nnode 2933: cm_p={uni['cm_p'][2933]:.4g}  = amp(L,fwd) + bmp(R,bwd)")
# amp uses fwd stream at left[2933]; bmp uses bwd stream at right[2933]
lsrc, rsrc = int(left[2933]), int(right[2933])
print(f"  left={lsrc} (fwd_mp={mp_f[lsrc]:.4g})   right={rsrc} (bwd_mp={mp_b[rsrc]:.4g})")
walk_back(lsrc, mp_f, csrc_f, left, "FWD via left nbr")
walk_back(rsrc, mp_b, csrc_b, right, "BWD via right nbr")
