"""AUDIT #1 — PROVENANCE of cm_p (the RNA+ measurement precision) at node 1909.

Reconstructs the `mp` measurement-precision stream of `_relay`/`_transport` EXACTLY from the captured
static arrays, tracking per-contribution provenance: which node's discrete spliced+ count SP fed it,
and how sigma^2_transfer damped it hop-by-hop. Then decomposes cm_p[1909] into
  (a) LOCAL: the spliced+ counts at 1909's own two junctions, vs
  (b) PROPAGATED: spliced counts telescoped in from other junctions along the relay.
Compares precision_CLAIMED (cm_p) against the discrete spliced+ count available at 1909's OWN junctions
(Poisson bound: precision of an RNA measurement <= the spliced count that produced it).
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
TARGET = 1909
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
chain = dbg["chain"]
cap = dbg["capture"]
uni = cap["_uni"][-1]
st = cap["_uni_static"]

Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > 1e-9, G / np.maximum(G + R, 1e-9), np.nan)

# ---- static arrays (exactly what the relay/transport consume) ----
is_exon = np.asarray(st["is_exon"], bool)
is_bnd = np.asarray(st["is_bnd"], bool)
left = np.asarray(st["left"], np.int64)
right = np.asarray(st["right"], np.int64)
fp = np.asarray(st["fp"], bool)
fn = np.asarray(st["fn"], bool)
order = np.asarray(st["order"], np.int64)
logvar_tot = np.asarray(st["logvar_tot"], np.float64)
SP_l = np.asarray(st["SP_l"], np.float64)   # spliced+ at LEFT face  (SP[0])
SP_r = np.asarray(st["SP_r"], np.float64)   # spliced+ at RIGHT face (SP[1])
mass = np.asarray(cap["mass_global"], np.float64)
cm_p = np.asarray(uni["cm_p"], np.float64)
fwd_mp_solver = np.asarray(st["fwd_mp"], np.float64)
bwd_mp_solver = np.asarray(st["bwd_mp"], np.float64)
n = left.shape[0]

print(f"=== node {TARGET}: cm_p (solver) = {cm_p[TARGET]:.4f}   mass={mass[TARGET]:,.0f} "
      f"oracle_fg={fo[TARGET]:.3f}  is_exon={is_exon[TARGET]}")
print(f"    left nbr={left[TARGET]} (exon={is_exon[left[TARGET]]},bnd={is_bnd[left[TARGET]]},"
      f"mass={mass[left[TARGET]]:,.0f})   right nbr={right[TARGET]} "
      f"(exon={is_exon[right[TARGET]]},bnd={is_bnd[right[TARGET]]},mass={mass[right[TARGET]]:,.0f})")

# ---------------------------------------------------------------------------
# Reconstruct the forward/backward mp streams EXACTLY, carrying per-origin provenance.
# mp[i] contributions are a dict {origin_boundary_node: delivered_precision}. A graft into exon i from
# boundary s deposits SP_face[s] (tagged origin=s). Propagation damps the WHOLE mp by _damp(., s2t);
# we attribute that damping proportionally across existing contributions (exact when s2t=0, i.e. on grafts).
# ---------------------------------------------------------------------------
def _damp(p, s2t):
    return p / (1.0 + p * s2t) if p > 0.0 else 0.0

def run_relay(seq, nbr, sf, df):
    # sf: source face for the graft spliced (1=right for fwd, 0=left for bwd)
    SP_src_face = SP_r if sf == 1 else SP_l
    mp = np.zeros(n)
    prov = [dict() for _ in range(n)]  # per node: {origin_node: amount}
    for i in seq:
        s = nbr[i]
        if s < 0:
            continue
        gr = bool(is_exon[i] and is_bnd[s])
        s2t = 0.0 if gr else (logvar_tot[i] + logvar_tot[s])
        # propagate source contributions, damped
        src_tot = mp[s]
        damped = _damp(src_tot, s2t)
        scale = (damped / src_tot) if src_tot > _EPS else 0.0
        contrib = {o: a * scale for o, a in prov[s].items()}
        tmp = damped
        if gr:
            spc = SP_src_face[s] / (1.0 + SP_src_face[s] * s2t) if SP_src_face[s] > _EPS else 0.0
            tmp += spc
            if spc > 0.0:
                contrib[s] = contrib.get(s, 0.0) + spc  # origin = the grafting boundary s
        if not fp[i]:
            tmp = 0.0
            contrib = {}
        mp[i] = 0.0 + tmp  # mp_own = 0
        prov[i] = contrib
    return mp, prov

fwd_mp, fwd_prov = run_relay(list(order), left, sf=1, df=0)
bwd_mp, bwd_prov = run_relay(list(order[::-1]), right, sf=0, df=1)

# Validate reconstruction vs solver-captured arrays
print(f"\n[validate] reconstructed fwd_mp[{TARGET}]={fwd_mp[TARGET]:.6f} vs solver={fwd_mp_solver[TARGET]:.6f}  "
      f"max|diff| over chain={np.max(np.abs(fwd_mp - fwd_mp_solver)):.3e}")
print(f"[validate] reconstructed bwd_mp[{TARGET}]={bwd_mp[TARGET]:.6f} vs solver={bwd_mp_solver[TARGET]:.6f}  "
      f"max|diff| over chain={np.max(np.abs(bwd_mp - bwd_mp_solver)):.3e}")

# ---------------------------------------------------------------------------
# Reconstruct cm_p[1909] from the combine _transport (df/sf swapped for the two messages).
#   left  message: src=left[i], sf=1  -> amp = _dv(fwd_mp[src], s2t) + graft _spc(SP_r[src])
#   right message: src=right[i], sf=0 -> bmp = _dv(bwd_mp[src], s2t) + graft _spc(SP_l[src])
# On a graft (exon i <- boundary src) s2t=0, so _dv is UNDAMPED and _spc is the RAW spliced count.
# ---------------------------------------------------------------------------
def transport_mp(i, src, mp_src_arr, prov_src, SP_face):
    if src < 0:
        return 0.0, {}, 0.0, 0.0
    gr = bool(is_exon[i] and is_bnd[src])
    s2t = 0.0 if gr else (logvar_tot[i] + logvar_tot[src])
    # _dv(mp[src], s2t)  (propagated measurement precision from the neighbour)
    dv = (1.0 / (1.0 / mp_src_arr[src] + s2t)) if mp_src_arr[src] > 0.0 else 0.0
    scale = (dv / mp_src_arr[src]) if mp_src_arr[src] > _EPS else 0.0
    prop = {o: a * scale for o, a in prov_src[src].items()}
    # graft spliced at this hop (local junction between i and src)
    spc = 0.0
    if gr and SP_face[src] > _EPS:
        s2t_spl = s2t if np.isfinite(s2t) else 0.0
        spc = SP_face[src] / (1.0 + SP_face[src] * s2t_spl)
    if not fp[i]:
        return 0.0, {}, 0.0, 0.0
    return dv + spc, prop, spc, SP_face[src]

amp, amp_prov, amp_local_spc, amp_local_SP = transport_mp(TARGET, left[TARGET], fwd_mp, fwd_prov, SP_r)
bmp, bmp_prov, bmp_local_spc, bmp_local_SP = transport_mp(TARGET, right[TARGET], bwd_mp, bwd_prov, SP_l)

print(f"\n=== cm_p[{TARGET}] reconstruction ===")
print(f"  amp (from LEFT  nbr {left[TARGET]}) = {amp:.4f}   [ _dv(propagated)={amp-amp_local_spc:.4f} + local graft spliced+ SP_r[{left[TARGET]}]={amp_local_spc:.4f} (raw SP={amp_local_SP:.4f}) ]")
print(f"  bmp (from RIGHT nbr {right[TARGET]}) = {bmp:.4f}   [ _dv(propagated)={bmp-bmp_local_spc:.4f} + local graft spliced+ SP_l[{right[TARGET]}]={bmp_local_spc:.4f} (raw SP={bmp_local_SP:.4f}) ]")
print(f"  cm_p (reconstructed) = {amp + bmp:.4f}   (solver = {cm_p[TARGET]:.4f})")

local = amp_local_spc + bmp_local_spc
propagated = (amp - amp_local_spc) + (bmp - bmp_local_spc)
print(f"\n  LOCAL (1909's OWN two junctions' spliced+ counts) = {local:.4f}")
print(f"  PROPAGATED (telescoped from other nodes)          = {propagated:.4f}")
print(f"  => propagated fraction of cm_p = {propagated / max(amp+bmp, _EPS):.1%}")

# ---------------------------------------------------------------------------
# Provenance: where does the PROPAGATED precision come from? List origin junctions + their raw spliced+.
# ---------------------------------------------------------------------------
merged = {}
for o, a in amp_prov.items():
    merged[o] = merged.get(o, 0.0) + a
for o, a in bmp_prov.items():
    merged[o] = merged.get(o, 0.0) + a
# also add the two local grafts as origins for a full census
if amp_local_spc > 0:
    merged[left[TARGET]] = merged.get(left[TARGET], 0.0) + amp_local_spc
if bmp_local_spc > 0:
    merged[right[TARGET]] = merged.get(right[TARGET], 0.0) + bmp_local_spc

items = sorted(merged.items(), key=lambda kv: -kv[1])
print(f"\n=== FULL ORIGIN CENSUS of cm_p[{TARGET}] (origin boundary node -> delivered precision) ===")
print(f"{'origin':>7} {'delivered':>10} {'raw SP_l':>9} {'raw SP_r':>9} {'mass':>10} {'is_bnd':>6} {'|dist chain-hops from 1909'}")
def hop_dist(node):
    # distance in chain hops walking left/right adjacency from TARGET
    for d, direction in ((0, None),):
        pass
    # BFS
    from collections import deque
    seen = {TARGET}
    dq = deque([(TARGET, 0)])
    while dq:
        cur, dd = dq.popleft()
        if cur == node:
            return dd
        for nb in (left[cur], right[cur]):
            if nb >= 0 and nb not in seen:
                seen.add(nb)
                dq.append((nb, dd + 1))
    return -1
for o, a in items:
    print(f"{o:>7} {a:>10.4f} {SP_l[o]:>9.3f} {SP_r[o]:>9.3f} {mass[o]:>10,.0f} {str(bool(is_bnd[o])):>6}  {hop_dist(o)}")

total_local_SP = amp_local_SP + bmp_local_SP
print(f"\n=== HONEST-PRECISION INVARIANT CHECK ===")
print(f"  precision CLAIMED at node {TARGET}: cm_p = {cm_p[TARGET]:.4f}")
print(f"  discrete spliced+ counts at {TARGET}'s OWN two junctions: SP_r[{left[TARGET]}]+SP_l[{right[TARGET]}] "
      f"= {amp_local_SP:.4f} + {bmp_local_SP:.4f} = {total_local_SP:.4f}")
print(f"  Poisson bound (precision of an RNA+ measurement <= originating spliced count): "
      f"cm_p {'>' if cm_p[TARGET] > total_local_SP else '<='} own-junction count")
if cm_p[TARGET] > total_local_SP + 1e-6:
    print(f"  ** VIOLATION: claimed {cm_p[TARGET]:.2f} > {total_local_SP:.2f} available at own junctions; "
          f"excess {cm_p[TARGET]-total_local_SP:.2f} is propagated from distant nodes. **")
