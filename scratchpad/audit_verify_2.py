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
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
os.environ.pop("RIGEL_S2T_OFF", None)
dbg = {}
cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
calmod.calibrate(inp["payload"], ra, inp["strand_model"],
                 np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)

chain = dbg["chain"]; cap = dbg["capture"]
st = cap["_uni_static"]; uni = cap["_uni"][-1]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > 1e-9, G / np.maximum(G + R, 1e-9), np.nan)

order = st["order"]
left = st["left"].astype(np.int64); right = st["right"].astype(np.int64)
is_bnd = st["is_bnd"].astype(bool); is_exon = st["is_exon"].astype(bool)
fp = st["fp"].astype(bool); fn = st["fn"].astype(bool)
logvar = st["logvar_tot"]
SP_l, SP_r = st["SP_l"], st["SP_r"]
SN_l, SN_r = st["SN_l"], st["SN_r"]
fwd_mp, bwd_mp = st["fwd_mp"], st["bwd_mp"]
fwd_mg, bwd_mg = st["fwd_mg"], st["bwd_mg"]
fwd_tau, bwd_tau = st["fwd_tau"], st["bwd_tau"]
tau_own = st["tau_own"]; mg_own = st["mg_own"]; struct_lock = st["struct_lock"].astype(bool)
n_unspl_l, n_unspl_r = st["n_unspl_l"], st["n_unspl_r"]
# per-node discrete spliced fragment counts (for the count-honesty bound)
snp_l, snp_r = st["spl_n_pos_l"], st["spl_n_pos_r"]
snn_l, snn_r = st["spl_n_neg_l"], st["spl_n_neg_r"]

cm_p = uni["cm_p"]; cm_g = uni["cm_g"]; cm_n = uni["cm_n"]; c_tau = uni["c_tau"]
amp, bmp = uni["amp"], uni["bmp"]
mo_p = uni["mo_p"]; fg_out = uni["fg_out"]
N = cm_p.shape[0]
EPS = 1e-9


def _damp(p, s2t):
    return 1.0 / (1.0 / p + s2t) if p > 0.0 else 0.0


print("=" * 90)
print("PART 0: SEED HONESTY — is a graft's injected precision (SP mass) <= its discrete spliced COUNT?")
print("=" * 90)
# SP is fractional accumulator MASS; the discrete count is spl_n_pos.  precision seed = SP <= count?
for face, spm, spc in (("L+", SP_l, snp_l), ("R+", SP_r, snp_r), ("L-", SN_l, snn_l), ("R-", SN_r, snn_r)):
    m = spm > EPS
    ratio = spm[m] / np.maximum(spc[m], EPS)
    viol = np.mean(ratio > 1.0 + 1e-6) if m.sum() else 0.0
    print(f"  face {face}: n_seed={m.sum():5d}  max(SP_mass/count)={ratio.max() if m.sum() else 0:.4f}  frac(mass>count)={viol:.4f}")
print("  -> if all <=1: the ORIGINATING seed count is NOT inflated; any cm_p excess is TRANSPORT, not seed.")

print()
print("=" * 90)
print("PART 1: is an exon's OWN CONTAINED spliced count structurally ZERO? (audit#1 used it as the ceiling)")
print("=" * 90)
exon = is_exon & (np.asarray(chain.kind) == REGION)
own_spl_exon = (snp_l + snp_r + snn_l + snn_r)[exon]
print(f"  exons: {exon.sum()}   exons with own contained spliced count>0: {(own_spl_exon > EPS).sum()}")
print(f"  mean own contained spliced count on exons = {own_spl_exon.mean():.4f}  (max {own_spl_exon.max():.2f})")
print("  -> spliced fragments live on BOUNDARY nodes; 'exon own spliced count' as a ceiling is structurally wrong.")

print()
print("=" * 90)
print("PART 2: RECONSTRUCT the measurement relay with PROVENANCE; validate vs captured fwd_mp/bwd_mp")
print("=" * 90)


def relay_mp(seq, nbr, sf_sp, free):
    """Replicate the pos measurement stream + track provenance {origin_junction: precision}."""
    mp = np.zeros(N)
    prov = [dict() for _ in range(N)]
    for i in seq:
        s = int(nbr[i])
        if s < 0:
            continue
        gr = is_exon[i] and is_bnd[s]
        s2t = 0.0 if gr else (logvar[i] + logvar[s])
        # damp incoming, propagate provenance
        if mp[s] > 0.0:
            fac = _damp(mp[s], s2t) / mp[s]
            tmp = mp[s] * fac
            pr = {o: v * fac for o, v in prov[s].items()}
        else:
            tmp = 0.0
            pr = {}
        if gr:
            spc = sf_sp[s] / (1.0 + sf_sp[s] * s2t) if sf_sp[s] > EPS else 0.0
            tmp += spc
            pr[s] = pr.get(s, 0.0) + spc  # origin = this adjacent junction s
        if not free[i]:
            tmp = 0.0
            pr = {}
        mp[i] = tmp
        prov[i] = pr
    return mp, prov


rec_fwd, prov_fwd = relay_mp([int(x) for x in order], left, SP_r, fp)     # fwd: nbr=left, sf=1 -> SP_r
rec_bwd, prov_bwd = relay_mp([int(x) for x in order][::-1], right, SP_l, fp)  # bwd: nbr=right, sf=0 -> SP_l
print(f"  fwd_mp reconstruction max|Δ| = {np.max(np.abs(rec_fwd - fwd_mp)):.3e}")
print(f"  bwd_mp reconstruction max|Δ| = {np.max(np.abs(rec_bwd - bwd_mp)):.3e}")

# reconstruct amp/bmp (combine) and cm_p
graft_L = is_exon & is_bnd[np.clip(left, 0, N - 1)] & (left >= 0)
graft_R = is_exon & is_bnd[np.clip(right, 0, N - 1)] & (right >= 0)
s2t_L = np.where(graft_L, 0.0, logvar + logvar[np.clip(left, 0, N - 1)])
s2t_R = np.where(graft_R, 0.0, logvar + logvar[np.clip(right, 0, N - 1)])


def _dv_vec(psrc, srcidx, s2t, valid):
    out = np.zeros(N)
    for i in range(N):
        s = int(srcidx[i])
        if not valid[i] or s < 0 or psrc[s] <= 0.0:
            continue
        out[i] = 1.0 / (1.0 / psrc[s] + s2t[i])
    return out


rec_amp = _dv_vec(fwd_mp, left, s2t_L, left >= 0) + np.where(graft_L & fp, SP_r[np.clip(left, 0, N - 1)], 0.0)
rec_bmp = _dv_vec(bwd_mp, right, s2t_R, right >= 0) + np.where(graft_R & fp, SP_l[np.clip(right, 0, N - 1)], 0.0)
rec_amp = np.where(fp, rec_amp, 0.0); rec_bmp = np.where(fp, rec_bmp, 0.0)
rec_cm_p = rec_amp + rec_bmp
print(f"  cm_p reconstruction max|Δ| = {np.max(np.abs(rec_cm_p - cm_p)):.3e}  (validates the decomposition)")

# honest 1-hop (adjacent junction) contribution vs foreign accumulation
honest_1hop = (np.where(graft_L & fp, SP_r[np.clip(left, 0, N - 1)], 0.0)
               + np.where(graft_R & fp, SP_l[np.clip(right, 0, N - 1)], 0.0))
foreign = cm_p - honest_1hop
# adjacent-junction discrete-count ceiling: the spliced count at the two flanking boundaries facing the exon
adj_count = (np.where(graft_L, snp_r[np.clip(left, 0, N - 1)], 0.0)
             + np.where(graft_R, snp_l[np.clip(right, 0, N - 1)], 0.0))

print()
print("=" * 90)
print("PART 3: honest 1-hop graft vs foreign accumulation (which precisions are count-backed?)")
print("=" * 90)
solv = cm_p > 0.5
print(f"  nodes with cm_p>0.5: {solv.sum()}")
print(f"  1-hop honest precision  <= adjacent-junction COUNT? frac(honest>adj_count+1e-6) = "
      f"{np.mean(honest_1hop[solv] > adj_count[solv] + 1e-6):.4f}   (expect 0 -> honest seed is count-bounded)")
frac_foreign = foreign[solv] / np.maximum(cm_p[solv], EPS)
print(f"  foreign fraction of cm_p: mean={frac_foreign.mean():.3f}  median={np.median(frac_foreign):.3f}  "
      f"p90={np.quantile(frac_foreign, .9):.3f}")
print(f"  nodes where cm_p is MAJORITY foreign (>50%): {np.mean(frac_foreign > 0.5):.4f}")
# the audit#1 ceiling vs the correct ceiling
own_contained = (snp_l + snp_r)  # exon own contained pos spliced (~0)
print(f"  audit#1 ceiling (own contained pos spliced): frac(cm_p>own_contained) among solv = "
      f"{np.mean(cm_p[solv] > own_contained[solv]):.4f}   <-- 68.9% claim reproduced against STRUCTURALLY-ZERO ceiling")
print(f"  correct ceiling (adjacent junction count):   frac(cm_p>adj_count)      among solv = "
      f"{np.mean(cm_p[solv] > adj_count[solv] + 1e-6):.4f}")
print(f"  strict single-origin (audit's own):          frac(honest 1-hop term ALONE > adj_count) = "
      f"{np.mean(honest_1hop[solv] > adj_count[solv] + 1e-6):.4f}")

print()
print("=" * 90)
print("PART 4: are the OTHER streams (c_tau composition, cm_g anchor gDNA) count-honest?")
print("=" * 90)
# composition tau: enters only via strand Fisher (I_strand) or intron NB curvature -> per-node own tau_own,
# accumulated along the relay.  Bound: the max single reachable tau_own on the path.
n_strand = (n_unspl_l + n_unspl_r)  # discrete unspliced count backing the strand Fisher
print(f"  c_tau: max = {c_tau.max():.3f};  max own tau_own = {tau_own.max():.3f}")
print(f"  own tau_own <= n_strand (Fisher <= count)? frac(tau_own>n_strand+1e-6) = "
      f"{np.mean(tau_own > n_strand + 1e-6):.4f}")
ct = c_tau > 0.5
print(f"  c_tau vs max reachable own tau across chain: (accumulation present) c_tau>max(tau_own)? "
      f"frac = {np.mean(c_tau[ct] > tau_own.max() + 1e-6):.4f}")
# cm_g: struct_lock only own count (mg_own), accumulated
print(f"  cm_g: max = {cm_g.max():.3f};  mg_own only on struct_lock: max mg_own = {mg_own.max():.3f}")
print(f"  mg_own <= own unspliced count? frac(mg_own > n_strand+1e-6) = "
      f"{np.mean(mg_own > n_strand + 1e-6):.4f}")

print()
print("=" * 90)
print("PART 5: ANCHOR NODE 1909 — full provenance of cm_p")
print("=" * 90)
i = 1909
print(f"  node {i}: exon={is_exon[i]} oracle_fg={fo[i]:.3f} solved_fg={fg_out[i]:.3f} mode_f_pos={np.exp(mo_p[i]):.3f}")
print(f"  cm_p={cm_p[i]:.3f}  = amp {amp[i]:.3f} + bmp {bmp[i]:.3f}")
print(f"  left nbr {int(left[i])} (bnd={is_bnd[int(left[i])]}) SP_r={SP_r[int(left[i])]:.3f} spl_count_r={snp_r[int(left[i])]:.3f}")
print(f"  right nbr {int(right[i])} (bnd={is_bnd[int(right[i])]}) SP_l={SP_l[int(right[i])]:.3f} spl_count_l={snp_l[int(right[i])]:.3f}")
print(f"  honest 1-hop precision = {honest_1hop[i]:.3f}  (adjacent junction count ceiling = {adj_count[i]:.3f})")
print(f"  foreign accumulation   = {foreign[i]:.3f}  ({100*foreign[i]/max(cm_p[i],EPS):.1f}% of cm_p)")
# top provenance origins
allprov = {}
for o, v in prov_fwd[i].items():
    allprov[o] = allprov.get(o, 0.0) + _damp(fwd_mp[int(left[i])], s2t_L[i]) / max(fwd_mp[int(left[i])], EPS) * v if fwd_mp[int(left[i])] > 0 else allprov.get(o, 0.0)
# simpler: report the relay provenance of the two source arrays directly (pre-combine damp)
print("  provenance of fwd_mp[left] (upstream junctions carried into the left boundary), top 6 by contribution:")
topf = sorted(prov_fwd[int(left[i])].items(), key=lambda kv: -kv[1])[:6]
for o, v in topf:
    hops = "adjacent" if o in (int(left[i]), int(right[i])) else "FOREIGN"
    print(f"     origin junction {o}: raw contrib {v:.3f}  [{hops}]  SP_r_origin={SP_r[o]:.3f}")
print("  provenance of bwd_mp[right], top 6:")
topb = sorted(prov_bwd[int(right[i])].items(), key=lambda kv: -kv[1])[:6]
for o, v in topb:
    hops = "adjacent" if o in (int(left[i]), int(right[i])) else "FOREIGN"
    print(f"     origin junction {o}: raw contrib {v:.3f}  [{hops}]  SP_l_origin={SP_l[o]:.3f}")

print()
print("=" * 90)
print("PART 6: TOP cm_p offenders on gDNA exons (oracle_fg>0.9, own contained spliced ~0)")
print("=" * 90)
mask = exon & (fo > 0.9) & (cm_p > 2.0)
idx = np.argsort(-cm_p * mask)[:8]
print(f"  {'node':>6} {'oracle':>7} {'solved':>7} {'cm_p':>10} {'honest1hop':>11} {'foreign':>10} {'adj_cnt':>8} {'own_spl':>8}")
for j in idx:
    if not mask[j]:
        continue
    print(f"  {j:>6} {fo[j]:>7.3f} {fg_out[j]:>7.3f} {cm_p[j]:>10.2f} {honest_1hop[j]:>11.2f} "
          f"{foreign[j]:>10.2f} {adj_count[j]:>8.2f} {own_contained[j]:>8.2f}")
