"""audit_4 — PROPAGATION SEMANTICS of the spliced RNA MEASUREMENT stream (mp/mn).

Claim under test: a spliced measurement at junction J measures J's RNA. Does the measurement stream (mp)
PROPAGATE J's precision to distant nodes and ADD it as if it measured them — inflating cm_p at nodes far
from J? Is the measurement precision DECAYING with distance (honest) or ACCUMULATING (dishonest)?

Method: reconstruct the forward/backward measurement relay for node 1909 exactly, walk the relay path hop
by hop, and show (a) where each spliced count is injected, (b) the s2t damping applied at each hop,
(c) whether the precision decays. Compare precision_CLAIMED (cm_p) vs originating discrete spliced count.
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

def run(s2t_off):
    os.environ.pop("RIGEL_S2T_OFF", None)
    if s2t_off: os.environ["RIGEL_S2T_OFF"] = "1"
    dbg = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return dbg

def analyze(dbg, tag):
    chain = dbg["chain"]; cap = dbg["capture"]; uni = cap["_uni"][-1]
    st = cap["_uni_static"]
    left = np.asarray(chain.left); right = np.asarray(chain.right)
    order = np.asarray(st["order"])
    fwd_mp = np.asarray(st["fwd_mp"]); bwd_mp = np.asarray(st["bwd_mp"])
    SP_l = np.asarray(st["SP_l"]); SP_r = np.asarray(st["SP_r"])
    logvar = np.asarray(st["logvar_tot"])
    mass = np.asarray(cap["mass_global"])
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > 1e-9, G / np.maximum(G + R, 1e-9), np.nan)
    rt, _ = _node_region_type(chain, ra)
    CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
    kind = np.asarray(chain.kind)
    def cl(i):
        return CLS[int(rt[i])] if kind[i] == REGION else "boundary"
    amp = np.asarray(uni["amp"]); bmp = np.asarray(uni["bmp"])
    cm_p = np.asarray(uni["cm_p"])
    i = 1909
    print(f"\n{'='*90}\n[{tag}]  node {i} [{cl(i)}] mass={mass[i]:,.0f} oracle_fg={fo[i]:.3f}")
    print(f"  cm_p (CLAIMED RNA+ measurement precision) = {cm_p[i]:.4f}")
    print(f"    = amp[left-stream]={amp[i]:.4f} + bmp[right-stream]={bmp[i]:.4f}")
    print(f"  fwd_mp[{i}]={fwd_mp[i]:.4f}  bwd_mp[{i}]={bwd_mp[i]:.4f}  (relay accumulators, pre-transport)")

    # Walk the FORWARD relay path INTO node i: the forward pass proceeds along `order`; at each node the
    # source is its LEFT neighbour. Reconstruct the chain of left-predecessors back from i.
    def walk(nbr, mp_arr, SP_face_src, facename):
        # trace back along nbr[] from i, showing mp accumulation + spliced injections
        print(f"  --- {facename} relay path into {i} (source = {facename} neighbour each hop) ---")
        path = []
        j = i
        seen = set()
        while j >= 0 and j not in seen:
            seen.add(j); path.append(j); j = int(nbr[j])
        # path[0]=i, path[1]=its source, ... walk from far end toward i to show accumulation
        for k in range(len(path) - 1, -1, -1):
            node = path[k]
            src = int(nbr[node]) if k < len(path) else -1
            src = path[k + 1] if k + 1 < len(path) else -1
            spc_src = SP_face_src[src] if src >= 0 else 0.0
            lv = logvar[node] + (logvar[src] if src >= 0 else 0.0)
            print(f"    node {node:5d} [{cl(node):10s}] mass={mass[node]:>10,.0f}  mp={mp_arr[node]:8.4f}"
                  f"  logvar_tot={logvar[node]:8.4g}  src_spliced_inject(SP@src)={spc_src:8.2f}  s2t(node,src)={lv:8.4g}")
        return path
    fpath = walk(left, fwd_mp, SP_r, "LEFT(fwd)")   # forward source faces src's RIGHT face -> SP_r
    bpath = walk(right, bwd_mp, SP_l, "RIGHT(bwd)")  # backward source faces src's LEFT face -> SP_l

    # Originating spliced counts anywhere on the two paths — the honest ceiling for cm_p at i is the
    # spliced count at i's OWN flanking junction, NOT distant ones.
    all_path = set(fpath) | set(bpath)
    inj = [(j, SP_l[j] + SP_r[j]) for j in all_path if (SP_l[j] + SP_r[j]) > 1e-9]
    inj.sort(key=lambda x: -x[1])
    print(f"  --- spliced-count injection sites on the relay paths (node, total SP) ---")
    for j, s in inj[:12]:
        d_fwd = fpath.index(j) if j in fpath else None
        d_bwd = bpath.index(j) if j in bpath else None
        d = d_fwd if d_fwd is not None else d_bwd
        print(f"    node {j:5d} [{cl(j):10s}] SP={s:8.2f}  hops_from_1909={d}")
    return dict(cm_p=cm_p[i], amp=amp[i], bmp=bmp[i], inj=inj, fpath=fpath, bpath=bpath,
                fwd_mp=fwd_mp, bwd_mp=bwd_mp, mass=mass, SP_l=SP_l, SP_r=SP_r)

print("RUNNING σ²_transfer ON ...")
dbg_on = run(False)
res_on = analyze(dbg_on, "S2T ON")
print("\n\nRUNNING σ²_transfer OFF ...")
dbg_off = run(True)
res_off = analyze(dbg_off, "S2T OFF")

print(f"\n\n{'#'*90}\nSUMMARY — cm_p at node 1909:")
print(f"  S2T ON : cm_p={res_on['cm_p']:.3f}  (amp={res_on['amp']:.3f} + bmp={res_on['bmp']:.3f})")
print(f"  S2T OFF: cm_p={res_off['cm_p']:.3f}  (amp={res_off['amp']:.3f} + bmp={res_off['bmp']:.3f})")
print(f"\n  node 1909 own flanking spliced (SP_l+SP_r)[1909] = {res_on['SP_l'][1909]+res_on['SP_r'][1909]:.3f}")
print(f"  Left neighbour (1908) SP = {res_on['SP_l'][1908]+res_on['SP_r'][1908]:.3f}")
print(f"  Right neighbour (1910) SP = {res_on['SP_l'][1910]+res_on['SP_r'][1910]:.3f}")

# ── GENOME-WIDE test of the invariant: precision_claimed (cm_p) vs the node's OWN originating spliced count.
cap = dbg_on["capture"]; uni = cap["_uni"][-1]; st = cap["_uni_static"]
chain = dbg_on["chain"]
cm_p = np.asarray(uni["cm_p"]); cm_n = np.asarray(uni["cm_n"])
mo_p = np.asarray(uni["mo_p"])
SP = np.asarray(st["SP_l"]) + np.asarray(st["SP_r"])
SN = np.asarray(st["SN_l"]) + np.asarray(st["SN_r"])
own_rna_ct = SP + SN
logvar = np.asarray(st["logvar_tot"])
mass = np.asarray(cap["mass_global"])
print(f"\n  logvar_tot stats (S2T ON): min={logvar.min():.4g} median={np.median(logvar):.4g} "
      f"max={logvar.max():.4g}  frac==0: {(logvar==0).mean():.3f}")
print(f"  logvar_tot[1909]={logvar[1909]:.4g}  logvar_tot[1910]={logvar[1910]:.4g}  logvar_tot[1963]={logvar[1963]:.4g}")

# Nodes claiming RNA measurement precision with ~ZERO own spliced count = precision entirely IMPORTED.
imported = (cm_p + cm_n > 1.0) & (own_rna_ct < 1e-6)
print(f"\n  INVARIANT TEST — nodes with cm_p+cm_n > 1 but OWN spliced count == 0 (imported measurement):")
print(f"    count = {imported.sum()} nodes")
idx = np.argsort(-(cm_p + cm_n) * imported)[:10]
mo_p_frac = np.exp(mo_p)
for j in idx:
    j = int(j)
    if not imported[j]: continue
    print(f"    node {j:5d} mass={mass[j]:>10,.0f}  own_SP+SN={own_rna_ct[j]:.3f}  "
          f"cm_p={cm_p[j]:.3f} cm_n={cm_n[j]:.3f}  mo_p(f_pos mode)={mo_p_frac[j]:.3f}")

# ── The HONEST CEILING: an exon legitimately receives ONLY its ADJACENT boundary's spliced (1-hop graft).
# Anything above that is precision imported from junctions >1 hop away and ADDED additively.
left = np.asarray(chain.left); right = np.asarray(chain.right)
SPtot = SP + SN  # total spliced count per node (nonzero on boundaries = junctions)
adj_ceiling = np.zeros_like(cm_p)
for i in range(cm_p.shape[0]):
    c = 0.0
    for nb in (left[i], right[i]):
        if nb >= 0:
            c = max(c, SPtot[nb])   # the 1-hop adjacent junction that legitimately grafts into i
    adj_ceiling[i] = c
claimed = cm_p + cm_n
# nodes where the RNA measurement precision EXCEEDS the adjacent-junction count (imported from farther)
excess = (claimed > 1.0) & (claimed > 2.0 * np.maximum(adj_ceiling, 1e-9))
print(f"\n  HONEST-CEILING TEST — cm_p+cm_n vs the ADJACENT-junction spliced count (the only count that")
print(f"  legitimately measures this node's incoming mature RNA). Excess = imported from >1 hop away:")
print(f"    nodes with claimed > 2× adjacent-junction ceiling: {excess.sum()}")
idx2 = np.argsort(-(claimed - adj_ceiling) * excess)[:12]
for j in idx2:
    j = int(j)
    if not excess[j]: continue
    print(f"    node {j:5d} mass={mass[j]:>9,.0f}  claimed(cm_p+cm_n)={claimed[j]:9.1f}  "
          f"adjacent-junction ceiling={adj_ceiling[j]:8.1f}  EXCESS={claimed[j]-adj_ceiling[j]:9.1f}  "
          f"mode f_pos={np.exp(mo_p[j]):.3f}")
print(f"\n  → node 1909: claimed={claimed[1909]:.2f}  adjacent-junction ceiling={adj_ceiling[1909]:.2f}  "
      f"(its own boundary 1910 SP={SPtot[1910]:.2f}, plus imported from junction 1913 SP={SPtot[1912]:.2f}/1913 area)")
