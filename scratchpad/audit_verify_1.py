"""VERIFIER #1 — independent reproduction + adversarial test of the cm_p phantom-precision root cause.

Invariant under test: precision is bounded by DISCRETE COUNTS at the message ORIGIN (Poisson: prec<=n).
Enrichment (reframe r) may scale the MODE but must NEVER scale the PRECISION.
Falsifiable claims of the trace agents:
  (T1) at node 1909 cm_p=26.45, own spliced count=0, delivered mode f_pos=0.718 -> f_g collapses.
  (T2) cm_p provenance: adjacent boundary 1910 grafts SP undamped (s2t=0), + additive accumulation of
       far junctions along the relay (mp[i]=mp_own+tmp ; cm_p=amp+bmp).
  (T3) genome-wide: many nodes carry cm_p >> own spliced count / adjacent-junction ceiling.
My additions:
  (V1) exact provenance of bmp[1909] hop-by-hop: which junctions contribute, how much, single-origin bound.
  (V2) MODE-vs-PRECISION ablation: replay psi at 1909 with (mode,prec) vs (mode, prec_capped_to_own_count)
       vs (correct-ish mode, high prec). Which one collapses f_g? -> is the bug the MODE or the PRECISION?
  (V3) adversarial SHOULD/ SHOULD-NOT: high-mass gDNA exon next to tiny high-RNA boundary (trigger) vs a
       matched-composition node (no trigger).
  (V4) deeper: does the graft's undamped s2t=0 exemption, or the additive accumulation, dominate?
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
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND = "gdna_gdna300_ss_0.99_nrna_present_capture_on"
index = TranscriptIndex.load(str(SUITE / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")

cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
def run(s2t_off):
    if s2t_off: os.environ["RIGEL_S2T_OFF"] = "1"
    else: os.environ.pop("RIGEL_S2T_OFF", None)
    dbg = {}
    calmod.calibrate(inp["payload"], ra, inp["strand_model"],
                     np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return dbg

dbg = run(False)
chain = dbg["chain"]; cap = dbg["capture"]; uni = cap["_uni"][-1]; st = cap["_uni_static"]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > 1e-9, G / np.maximum(G + R, 1e-9), np.nan)
rt, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind)
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
def cls(i): return CLS[int(rt[i])] if kind[i] == REGION else "boundary"

order = st["order"]; left = st["left"]; right = st["right"]
SP_l, SP_r = st["SP_l"], st["SP_r"]; SN_l, SN_r = st["SN_l"], st["SN_r"]
logvar = st["logvar_tot"]; is_bnd = st["is_bnd"]; is_exon = st["is_exon"]
fp = st["fp"]; fn = st["fn"]
cm_p = uni["cm_p"]; cm_n = uni["cm_n"]; cm_g = uni["cm_g"]; c_tau = uni["c_tau"]
amp = uni["amp"]; bmp = uni["bmp"]; amn = uni["amn"]; bmn = uni["bmn"]
mo_p = uni["mo_p"]; mo_g = uni["mo_g"]; mo_n = uni["mo_n"]; lam_msg = uni["lam_msg"]
fg_out = uni["fg_out"]; mass = st["M"]
fwd_mp = st["fwd_mp"]; bwd_mp = st["bwd_mp"]
# own spliced count at a node = spliced mass on its own faces (a node observes ITS OWN junctions).
own_spl_p = SP_l + SP_r; own_spl_n = SN_l + SN_r

N = fg_out.shape[0]
print("=" * 90)
print("PART 1 — reproduce anchor node 1909")
print("=" * 90)
i = 1909
print(f"node {i} [{cls(i)}] mass={mass[i]:,.0f} oracle_fg={fo[i]:.3f} solved_fg={fg_out[i]:.3f}")
print(f"  own spliced count (SP_l+SP_r)={own_spl_p[i]:.3f}  free=({int(fp[i])},{int(fn[i])})")
print(f"  RNA+ MEASUREMENT: cm_p={cm_p[i]:.3f}  mode f_pos=exp(mo_p)={np.exp(mo_p[i]):.4f}")
print(f"     decompose cm_p = amp(left)={amp[i]:.4f} + bmp(right)={bmp[i]:.4f}")
print(f"  gDNA meas cm_g={cm_g[i]:.3f}  comp c_tau={c_tau[i]:.3f}  lam_msg fg_eq={1/(1+np.exp(-lam_msg[i])):.3f}")
for tag, j in (("L", int(left[i])), ("R", int(right[i]))):
    if j >= 0:
        print(f"  nbr {tag} node {j} [{cls(j)}] mass={mass[j]:,.1f} oracle={fo[j]:.3f} "
              f"own_spl_p={own_spl_p[j]:.3f} logvar_tot={logvar[j]:.4g}")

print("\n" + "=" * 90)
print("PART 2 — INVARIANT: precision cm_p vs discrete spliced COUNT available at origin")
print("=" * 90)
# origin-adjacency ceiling: the spliced count that actually crosses INTO node i is on its own faces + its
# 1-hop graft neighbours' faces. Build the max single reachable spliced count along each node's relay path.
solv = np.asarray(cap["solvable_mask"], bool)
# adjacent-junction ceiling: own face spliced + the spliced on the immediate boundary neighbours' facing side.
adj_ceiling = own_spl_p.copy()
for i2 in range(N):
    l, r = int(left[i2]), int(right[i2])
    if l >= 0 and is_bnd[l] and is_exon[i2]:
        adj_ceiling[i2] = max(adj_ceiling[i2], SP_r[l])  # left nbr faces us on its right
    if r >= 0 and is_bnd[r] and is_exon[i2]:
        adj_ceiling[i2] = max(adj_ceiling[i2], SP_l[r])
m = solv & (cm_p > 0.5)
print(f"solvable nodes with cm_p>0.5: {int(m.sum())}")
frac_own = float((cm_p[m] > np.maximum(own_spl_p[m], 1e-9)).mean())
frac_adj = float((cm_p[m] > 2.0 * np.maximum(adj_ceiling[m], 1e-9)).mean())
print(f"  frac( cm_p > OWN spliced count )        = {frac_own:.3f}")
print(f"  frac( cm_p > 2x adjacent-junction cap ) = {frac_adj:.3f}")
zero_own = solv & (cm_p + cm_n > 1.0) & (own_spl_p + own_spl_n < 1e-6)
print(f"  nodes claiming cm_p+cm_n>1 with ZERO own spliced count = {int(zero_own.sum())}")
# max single reachable spliced count along the whole chain (global bound the seeds cannot exceed at origin)
max_single = max(float(SP_l.max()), float(SP_r.max()))
print(f"  global max single spliced count anywhere = {max_single:.1f}; "
      f"frac(cm_p > that) = {float((cm_p[m] > max_single).mean()):.3f}")
# top offenders
idx = np.argsort(-cm_p * solv)[:8]
print("  top cm_p offenders (node, cls, cm_p, own_spl, oracle_fg, solved_fg):")
for j in idx:
    print(f"    {j} {cls(j):10s} cm_p={cm_p[j]:9.1f} own_spl={own_spl_p[j]:6.2f} "
          f"oracle={fo[j]:.3f} solved={fg_out[j]:.3f} mass={mass[j]:,.0f}")

print("\n" + "=" * 90)
print("PART 3 — bmp[1909] hop-by-hop provenance (walk the backward relay from 1909)")
print("=" * 90)
# Replicate the _relay backward accumulation restricted to what flows into 1909 through its RIGHT neighbour.
# bmp[1909] = _dv(bwd_mp)[right[1909]] transported. Trace the raw bwd_mp accumulation chain.
i = 1909; r = int(right[i])
print(f"node 1909 right neighbour = {r} [{cls(r)}]  bwd_mp[{r}]={bwd_mp[r]:.4f}")
# walk right-ward following 'right' pointers, printing spliced injected per hop
cur = r; hops = 0; contribs = []
seen = set()
while cur >= 0 and cur not in seen and hops < 4000:
    seen.add(cur); rr = int(right[cur])
    inj = 0.0
    if rr >= 0 and is_bnd[rr] and is_exon[cur]:  # graft into cur from rr (its right nbr on right pass)
        pass
    # spliced injected AT cur when cur is a graft target: source is the neighbour with spliced facing
    spc = own_spl_p[cur]
    if spc > 0.5:
        contribs.append((cur, spc, cls(cur)))
    cur = rr; hops += 1
contribs.sort(key=lambda x: -x[1])
tot = sum(c[1] for c in contribs)
print(f"  reachable spliced-bearing nodes along the R-relay: {len(contribs)}, sum SP={tot:,.1f}, "
      f"largest single={contribs[0][1]:.1f} at node {contribs[0][0]} [{contribs[0][2]}]" if contribs else "  none")
print(f"  raw bwd_mp accumulator at r={r}: {bwd_mp[r]:.2f}  (vs adjacent single graft {adj_ceiling[i]:.2f})")

print("\n" + "=" * 90)
print("PART 4 — MODE vs PRECISION ablation at node 1909 (the deeper question)")
print("=" * 90)
# Replay the psi solve at node 1909 under channel variants to see whether the MODE or the PRECISION collapses f_g.
stt = dbg["statics"]
# solver scalar params: kappa/od from calibrate debug, grids from the config we passed in.
P = dict(_kappa=float(dbg["rna_sense_frac"]),
         _od_g=float(dbg["calibration_priors"].gdna_strand_overdispersion),
         _od_r=float(dbg["calibration_priors"].rna_strand_overdispersion),
         _n_grid=int(cc.sweep_n_grid), _L=float(cc.sweep_logodds_window),
         _n_tilt=cc.sweep_n_tilt, _n_grid_ss=cc.sweep_n_grid_single_strand)
cap.update(P)
def solve_one(idx, gm, gp, rmp, rpp, rmn, rpn, lam, ctau):
    # single-node psi solve using the same backend; pass through the node's own strand data.
    sl = _solve_nodes_logodds_all(
        np.asarray([stt.u_pos[idx]]), np.asarray([stt.u_neg[idx]]),
        np.asarray([fp[idx]]), np.asarray([fn[idx]]),
        np.asarray([stt.mass_unspliced[idx]]), np.asarray([stt.mass_spliced[idx]]),
        kappa=cap["_kappa"], od_g=cap["_od_g"], od_r=cap["_od_r"],
        n_grid=cap["_n_grid"], L=cap["_L"], n_tilt=cap["_n_tilt"], n_grid_ss=cap["_n_grid_ss"],
        global_logprior=None if cap["global_lp"] is None else cap["global_lp"][idx:idx+1],
        gdna_imp_mode=None if gm is None else np.asarray([gm]),
        gdna_imp_prec=None if gp is None else np.asarray([gp]),
        rna_imp_mode=None if rmp is None else (np.asarray([rmp]), np.asarray([rmn])),
        rna_imp_prec=None if rpp is None else (np.asarray([rpp]), np.asarray([rpn])),
        lam_imp_mode=None if lam is None else np.asarray([lam]),
        lam_imp_prec=None if ctau is None else np.asarray([ctau]),
        lam_logprior=None if cap["intron_prior"] is None else None,
        fg_ref=np.asarray([cap["fg_init"][idx]]),
        fpos_ref=np.asarray([cap["fpos_init"][idx]]),
        fneg_ref=np.asarray([cap["fneg_init"][idx]]),
    )
    return float(sl.gdna_frac[0])

have_params = all(k in cap for k in ("_kappa", "_od_g", "_od_r", "_n_grid", "_L", "_n_tilt", "_n_grid_ss"))
if not have_params:
    print("  (solver scalar params not published on capture; skipping single-node replay — see PART 4b)")
else:
    i = 1909
    base = solve_one(i, mo_g[i], cm_g[i], mo_p[i], cm_p[i], mo_n[i], cm_n[i], lam_msg[i], c_tau[i])
    no_rna = solve_one(i, mo_g[i], cm_g[i], mo_p[i], 0.0, mo_n[i], 0.0, lam_msg[i], c_tau[i])
    cap_own = solve_one(i, mo_g[i], cm_g[i], mo_p[i], min(cm_p[i], own_spl_p[i]),
                        mo_n[i], min(cm_n[i], own_spl_n[i]), lam_msg[i], c_tau[i])
    print(f"  node 1909 oracle_fg={fo[i]:.3f}")
    print(f"    (a) full messages           f_g={base:.3f}  [reproduces solver {fg_out[i]:.3f}]")
    print(f"    (b) RNA meas precision -> 0  f_g={no_rna:.3f}")
    print(f"    (c) RNA prec capped to own spliced count(={own_spl_p[i]:.2f}) f_g={cap_own:.3f}")

print("\n" + "=" * 90)
print("PART 4b — genome-wide: does capping RNA-meas precision to own spliced count recover oracle?")
print("=" * 90)
if have_params:
    # batch replay: full vs capped, on solvable nodes, compare mean |f_g - oracle| weighted by mass.
    def batch(cap_prec):
        cpp = np.minimum(cm_p, own_spl_p) if cap_prec else cm_p
        cpn = np.minimum(cm_n, own_spl_n) if cap_prec else cm_n
        sl = _solve_nodes_logodds_all(
            stt.u_pos, stt.u_neg, fp, fn, stt.mass_unspliced, stt.mass_spliced,
            kappa=cap["_kappa"], od_g=cap["_od_g"], od_r=cap["_od_r"],
            n_grid=cap["_n_grid"], L=cap["_L"], n_tilt=cap["_n_tilt"], n_grid_ss=cap["_n_grid_ss"],
            global_logprior=cap["global_lp"],
            gdna_imp_mode=mo_g, gdna_imp_prec=cm_g,
            rna_imp_mode=(mo_p, mo_n), rna_imp_prec=(cpp, cpn),
            lam_imp_mode=lam_msg, lam_imp_prec=c_tau,
            theta_imp_mode=None, theta_imp_prec=None,
            lam_logprior=cap["intron_prior"],
            fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
        )
        return np.asarray(sl.gdna_frac)
    fg_full = batch(False); fg_capped = batch(True)
    ok = solv & np.isfinite(fo)
    w = mass * ok
    def mwae(x): return float(np.sum(w * np.abs(x - np.where(ok, fo, 0.0))) / np.sum(w))
    print(f"  mass-weighted |f_g - oracle| over solvable nodes:")
    print(f"    full (cm_p as-is)                 = {mwae(fg_full):.4f}")
    print(f"    RNA-meas prec capped to own count = {mwae(fg_capped):.4f}")
    moved = ok & (np.abs(fg_full - fg_capped) > 0.05)
    print(f"    nodes moved >0.05 by the cap: {int(moved.sum())} (mass {float(mass[moved].sum()):,.0f})")

print("\n" + "=" * 90)
print("PART 5 — adversarial: SHOULD-trigger vs SHOULD-NOT among real nodes")
print("=" * 90)
# SHOULD trigger: high-mass gDNA exon (oracle_fg high), adjacent to a low-mass boundary with spliced>0.
trig = np.zeros(N, bool)
for i2 in range(N):
    if not (is_exon[i2] and solv[i2] and mass[i2] > 5000 and fo[i2] > 0.9): continue
    for j in (int(left[i2]), int(right[i2])):
        if j >= 0 and is_bnd[j] and mass[j] < 500 and own_spl_p[j] > 1.0:
            trig[i2] = True
print(f"SHOULD-trigger set (hi-mass hi-gDNA exon next to tiny spliced boundary): {int(trig.sum())} nodes")
if trig.sum():
    med_err = np.abs(fg_out[trig] - fo[trig])
    print(f"  median |solved-oracle| = {np.median(med_err):.3f}; median cm_p={np.median(cm_p[trig]):.2f}; "
          f"median own_spl={np.median(own_spl_p[trig]):.3f}")
# SHOULD NOT: matched-composition solvable exon with its own spliced evidence (own_spl_p>2) and moderate mass.
nott = solv & is_exon & (own_spl_p > 2.0) & (mass > 1000) & (fo < 0.5) & np.isfinite(fo)
print(f"SHOULD-NOT set (RNA-rich exon with its OWN spliced evidence): {int(nott.sum())} nodes")
if nott.sum():
    print(f"  median |solved-oracle| = {np.median(np.abs(fg_out[nott]-fo[nott])):.3f}; "
          f"median cm_p={np.median(cm_p[nott]):.2f}; median own_spl={np.median(own_spl_p[nott]):.2f}")
    print(f"  here cm_p is JUSTIFIED by own count: frac(cm_p<=2*own_spl)="
          f"{float((cm_p[nott]<=2*own_spl_p[nott]).mean()):.2f}")

print("\n" + "=" * 90)
print("PART 6 — S2T on/off: does sigma2_transfer damp the accumulation? (graft is exempt either way)")
print("=" * 90)
dbg2 = run(True)
uni2 = dbg2["capture"]["_uni"][-1]
print(f"  cm_p[1909]:  S2T ON = {cm_p[1909]:.2f}   S2T OFF = {uni2['cm_p'][1909]:.2f}")
print(f"  solved f_g[1909]: S2T ON = {fg_out[1909]:.3f}  S2T OFF = {uni2['fg_out'][1909]:.3f}")
print("  -> if OFF cm_p RISES, s2t partially damps accumulation but the GRAFT (s2t=0) is undamped regardless.")
