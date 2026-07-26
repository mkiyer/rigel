"""audit_2 — INVARIANT SCAN across ALL nodes of the anchor scenario.

For every node compute per-message precision_CLAIMED (cm_g, cm_p, cm_n, c_tau, th_prec + the fused psi
precisions) and compare to the honest Poisson bound = the discrete fragment COUNT available at the message's
ORIGIN. Rank the worst violations (claimed >> count). Is the violation SYSTEMATIC (tiny-mass sources ->
high precision)? Quantify the distribution.

Origins by stream (from the relay code, bp_solver._unified_solve):
  * cm_g  (anchor gDNA measurement): seed mg_own = struct_lock * prec_g, prec_g <= n_unspl. Transported by
           _damp (<=1) + additive fuse. Origin count = struct_lock nodes' unspliced COUNT n_unspl.
  * cm_p/cm_n (spliced RNA measurement): seed 0; injected ONLY at a boundary->exon GRAFT as _spc = SP/(1+SP*s2t)
           with s2t==0 on the graft (matched-set exemption) => _spc == SP (the spliced MASS, UNDAMPED).
           Origin count = the spliced flux at the grafting boundary (integer spl_n_pos / spl_n_neg).
  * c_tau (composition tau_lam): seed tau_own = i_strand(single-strand) + intron NB curvature. Bounded by the
           strand/unspliced COUNT at the source. Transported damped + additive fuse.

The relay accumulates ADDITIVELY along the chain (fwd via left[], bwd via right[]) with per-hop damping only,
so cm at a node is a SUM of damped origin counts. The invariant is violated if a node's claimed precision
exceeds ANY single originating count (accumulation of foreign junctions), and grossly so if it exceeds the
node's OWN discrete fragment count.
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
dbg = {}
cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                 np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)

chain = dbg["chain"]; cap = dbg["capture"]; uni = cap["_uni"][-1]; st = cap["_uni_static"]
left = np.asarray(chain.left); right = np.asarray(chain.right)
kind = np.asarray(chain.kind)
N = left.shape[0]

# claimed precisions handed to psi
cm_g = np.asarray(uni["cm_g"]); cm_p = np.asarray(uni["cm_p"]); cm_n = np.asarray(uni["cm_n"])
c_tau = np.asarray(uni["c_tau"])
fp = np.asarray(st["fp"], bool); fn = np.asarray(st["fn"], bool)
is_amb = fp & fn
th_prec = np.where(is_amb, cm_p + cm_n, 0.0)
solvable = np.asarray(cap["solvable"], bool)

# seeds / geometry (origin counts)
SP_l = np.asarray(st["SP_l"]); SP_r = np.asarray(st["SP_r"])          # spliced+ MASS per face (boundary)
SN_l = np.asarray(st["SN_l"]); SN_r = np.asarray(st["SN_r"])
spnl = np.asarray(st["spl_n_pos_l"]); spnr = np.asarray(st["spl_n_pos_r"])  # spliced+ integer COUNT
snnl = np.asarray(st["spl_n_neg_l"]); snnr = np.asarray(st["spl_n_neg_r"])
nul = np.asarray(st["n_unspl_l"]); nur = np.asarray(st["n_unspl_r"])        # unspliced integer COUNT
tau_own = np.asarray(st["tau_own"]); mg_own = np.asarray(st["mg_own"])
struct = np.asarray(st["struct_lock"], bool)
mass = np.asarray(cap["mass_global"])
n_self = nul + (nur * (kind != REGION))   # region: contained count = n_unspl_left; boundary: both faces
n_self = np.where(kind == REGION, nul, nul + nur)
SP_node = SP_l + SP_r; SN_node = SN_l + SN_r
spn_node = spnl + spnr; snn_node = snnl + snnr

# oracle + class
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > 1e-9, G / np.maximum(G + R, 1e-9), np.nan)
rt, _ = _node_region_type(chain, ra)
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
def cl(i): return CLS[int(rt[i])] if kind[i] == REGION else "boundary"

# ---------------------------------------------------------------------------------------------------------
# Per-node ORIGIN reconstruction: walk the fwd (left[]) and bwd (right[]) relay chains, collect the graft
# spliced counts (SP/SN at boundary sources feeding an exon) and the tau/mg seeds that reach node i.
# Returns per node: sum + max single of the originating SPLICED counts (for cm_p/cm_n), the originating
# unspliced-count sum/max of struct_lock nodes (for cm_g), and tau seed sum/max.
# ---------------------------------------------------------------------------------------------------------
def relay_origins(nbr):
    """For each node, list of source nodes on its relay chain (excluding itself), nearest first."""
    chains = [None] * N
    for i in range(N):
        path = []
        j = int(nbr[i]); seen = {i}
        while j >= 0 and j not in seen:
            seen.add(j); path.append(j); j = int(nbr[j])
        chains[i] = path
    return chains

fwd_chain = relay_origins(left)   # sources reached going left
bwd_chain = relay_origins(right)  # sources reached going right

# Per the graft rule: forward source faces its RIGHT face (SP_r); backward source faces its LEFT face (SP_l).
def sp_origin(chains, sp_face_src, sn_face_src):
    sump = np.zeros(N); maxp = np.zeros(N); sumn = np.zeros(N); maxn = np.zeros(N)
    for i in range(N):
        for s in chains[i]:
            vp = sp_face_src[s]; vn = sn_face_src[s]
            if vp > 1e-9:
                sump[i] += vp; maxp[i] = max(maxp[i], vp)
            if vn > 1e-9:
                sumn[i] += vn; maxn[i] = max(maxn[i], vn)
    return sump, maxp, sumn, maxn

f_sump, f_maxp, f_sumn, f_maxn = sp_origin(fwd_chain, SP_r, SN_r)
b_sump, b_maxp, b_sumn, b_maxn = sp_origin(bwd_chain, SP_l, SN_l)
orig_sump = f_sump + b_sump                 # total spliced+ MASS reachable (both directions)
orig_maxp = np.maximum(f_maxp, b_maxp)      # largest SINGLE spliced+ count reachable
orig_sumn = f_sumn + b_sumn
orig_maxn = np.maximum(f_maxn, b_maxn)

# integer spliced counts reachable (the honest Poisson ceiling) — same walk over the count arrays
fc = relay_origins(left); bc = relay_origins(right)
def spn_origin(chains, cnt_face):
    s = np.zeros(N); m = np.zeros(N)
    for i in range(N):
        for src in chains[i]:
            v = cnt_face[src]
            if v > 1e-9:
                s[i] += v; m[i] = max(m[i], v)
    return s, m
fp_cnt_s, fp_cnt_m = spn_origin(fc, spnr); bp_cnt_s, bp_cnt_m = spn_origin(bc, spnl)
orig_cnt_sump = fp_cnt_s + bp_cnt_s; orig_cnt_maxp = np.maximum(fp_cnt_m, bp_cnt_m)

# gDNA measurement origin: struct_lock nodes' unspliced count reachable
mg_seed = np.where(struct, n_self, 0.0)   # honest ceiling of mg_own is n_self at struct nodes
def scalar_origin(chains, seed):
    s = np.zeros(N); m = np.zeros(N)
    for i in range(N):
        for src in chains[i]:
            v = seed[src]
            if v > 1e-9:
                s[i] += v; m[i] = max(m[i], v)
    return s, m
fg_s, fg_m = scalar_origin(fc, mg_seed); bg_s, bg_m = scalar_origin(bc, mg_seed)
orig_mg_sum = fg_s + bg_s; orig_mg_max = np.maximum(fg_m, bg_m)
# tau origin
ft_s, ft_m = scalar_origin(fc, tau_own); bt_s, bt_m = scalar_origin(bc, tau_own)
orig_tau_sum = ft_s + bt_s + tau_own; orig_tau_max = np.maximum(np.maximum(ft_m, bt_m), tau_own)

# =========================================================================================================
print("="*100)
print("SEED HONESTY: is the spliced-measurement precision SEED (SP mass) <= the integer spliced COUNT?")
bm = SP_node > 1e-9
print(f"  boundaries with SP>0: {bm.sum()};  fraction with SP_mass <= spn_count+1e-6: "
      f"{(SP_node[bm] <= spn_node[bm]+1e-6).mean():.3f}")
print(f"  max ratio SP_mass/spn_count over those: {np.nanmax(SP_node[bm]/np.maximum(spn_node[bm],1e-9)):.3f}")
print("  => if ratio<=~1 the seed itself is count-honest; inflation (if any) is in ACCUMULATION/reframe.\n")

# =========================================================================================================
print("="*100)
print("ANCHOR node 1909 — full message provenance")
i = 1909
print(f"  [{cl(i)}] mass={mass[i]:,.0f} oracle_fg={fo[i]:.3f} n_self(unspl count)={n_self[i]:.1f} "
      f"own spliced+ count={spn_node[i]:.1f}")
print(f"  CLAIMED:  cm_p={cm_p[i]:.3f}  cm_n={cm_n[i]:.3f}  cm_g={cm_g[i]:.3f}  c_tau={c_tau[i]:.3f}  "
      f"th_prec={th_prec[i]:.3f}")
print(f"  RNA+ measurement origin: sum(SP mass reachable)={orig_sump[i]:.2f}  "
      f"max single SP={orig_maxp[i]:.2f}  sum(integer spliced+ count)={orig_cnt_sump[i]:.2f}  "
      f"max single count={orig_cnt_maxp[i]:.2f}")
print(f"  neighbours: L={int(left[i])}[{cl(int(left[i]))}] SP={SP_node[int(left[i])]:.2f} "
      f"mass={mass[int(left[i])]:,.0f} | R={int(right[i])}[{cl(int(right[i]))}] "
      f"SP={SP_node[int(right[i])]:.2f} mass={mass[int(right[i])]:,.0f}")
print(f"  INVARIANT: cm_p / (max single originating spliced count) = "
      f"{cm_p[i]/max(orig_cnt_maxp[i],1e-9):.2f}x   "
      f"cm_p / (node's OWN spliced+ count {spn_node[i]:.1f}) = {cm_p[i]/max(spn_node[i],1e-9):.2f}x")

# =========================================================================================================
# DISTRIBUTION across all solvable nodes: how often does cm_p exceed (a) the largest single originating
# spliced count, and (b) the node's own spliced count?
print("\n" + "="*100)
print("DISTRIBUTION — RNA+ measurement (cm_p) vs honest origin, over SOLVABLE nodes with cm_p>0.5")
m = solvable & (cm_p > 0.5)
r_max = cm_p[m] / np.maximum(orig_cnt_maxp[m], 1e-9)   # claimed / largest single originating count
r_own = cm_p[m] / np.maximum(spn_node[m], 1e-9)        # claimed / node's OWN spliced count
print(f"  nodes: {m.sum()}")
for lab, arr in (("cm_p / max-single-origin-count", r_max), ("cm_p / node-own-spliced-count", r_own)):
    q = np.nanpercentile(arr, [50, 75, 90, 99, 100])
    print(f"    {lab:34s} p50={q[0]:7.2f} p75={q[1]:7.2f} p90={q[2]:7.2f} p99={q[3]:8.2f} max={q[4]:9.2f}")
print(f"  frac with cm_p > max single originating count : {(r_max > 1.001).mean():.3f}")
print(f"  frac with cm_p > node's own spliced count     : {(r_own > 1.001).mean():.3f}")

# systematic tiny-mass -> high precision? correlate node mass vs (cm_p / own spliced count)
print("\n  systematic check: is a HIGH cm_p concentrated on nodes whose OWN spliced count is ~0 (foreign)?")
foreign = m & (spn_node < 1.0) & (cm_p > 2.0)
print(f"    nodes with own spliced+ count<1 but cm_p>2 (RNA+ precision with ~no local spliced evidence): "
      f"{foreign.sum()}  (of {m.sum()} = {foreign.sum()/max(m.sum(),1):.3f})")
if foreign.sum():
    idx = np.where(foreign)[0]
    idx = idx[np.argsort(-cm_p[idx])][:12]
    print("    top such nodes (cm_p, own_spliced, max-origin-count, mass, class, oracle_fg, solved_fg):")
    fgout = np.asarray(uni["fg_out"])
    for j in idx:
        print(f"      node {j:5d} cm_p={cm_p[j]:7.2f} own_spl={spn_node[j]:4.1f} "
              f"max_orig_cnt={orig_cnt_maxp[j]:6.1f} mass={mass[j]:>9,.0f} [{cl(j):10s}] "
              f"oracle_fg={fo[j]:.3f} solved_fg={fgout[j]:.3f}")

# =========================================================================================================
print("\n" + "="*100)
print("WORST cm_p ABSOLUTE (top 15 solvable) — claimed RNA+ precision, provenance")
fgout = np.asarray(uni["fg_out"])
order = np.argsort(-np.where(solvable, cm_p, -1))[:15]
print(f"  {'node':>6} {'cm_p':>8} {'maxOrigCnt':>10} {'sumOrigCnt':>10} {'ownSpl':>7} {'mass':>10} "
      f"{'class':>10} {'oracle':>7} {'solved':>7}")
for j in order:
    print(f"  {j:6d} {cm_p[j]:8.2f} {orig_cnt_maxp[j]:10.1f} {orig_cnt_sump[j]:10.1f} {spn_node[j]:7.1f} "
          f"{mass[j]:10,.0f} {cl(j):>10} {fo[j]:7.3f} {fgout[j]:7.3f}")

# =========================================================================================================
print("\n" + "="*100)
print("OTHER STREAMS — cm_g (gDNA meas) and c_tau (composition) vs origin counts, over solvable nodes")
for name, claim, omax, osum, ceil in (
    ("cm_g", cm_g, orig_mg_max, orig_mg_sum, orig_mg_max),
    ("c_tau", c_tau, orig_tau_max, orig_tau_sum, orig_tau_max),
):
    mm = solvable & (claim > 0.5)
    r = claim[mm] / np.maximum(ceil[mm], 1e-9)
    q = np.nanpercentile(r, [50, 90, 99, 100]) if mm.sum() else [0,0,0,0]
    print(f"  {name}: nodes>0.5={mm.sum():5d}  claimed/max-single-origin  p50={q[0]:.2f} p90={q[1]:.2f} "
          f"p99={q[2]:.2f} max={q[3]:.2f}  frac(claim>maxOrigin)={ (r>1.001).mean() if mm.sum() else 0:.3f}")

print("\nDONE.")
