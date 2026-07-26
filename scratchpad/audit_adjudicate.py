import os, sys, dataclasses, importlib
import numpy as np
sys.path.insert(0,"/Users/mkiyer/proj/rigel/scripts/debug")
from pathlib import Path
from flagship_interrogate import _oracle_per_node
from selfsolve_diag import _scan_and_truth
from rigel.calibration.bp_solver import REGION
calmod=importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.node_geometry import _node_region_type
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

SUITE=Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND="gdna_gdna300_ss_0.99_nrna_present_capture_on"
index=TranscriptIndex.load(str(SUITE/"rigel_index")); cfg=PipelineConfig()
ra=RegionArrays.from_region_df(index.region_df,index.ref_name_to_id)
inp=_scan_and_truth(SUITE,COND,index,cfg,Path("/tmp/rigel_selfsolve"),SUITE/"_selfsolve_cache")
os.environ.pop("RIGEL_S2T_OFF",None)
dbg={}; cc=dataclasses.replace(cfg.calibration,calib_refit_iters=0)
calmod.calibrate(inp["payload"],ra,inp["strand_model"],np.asarray(inp["gdna_fl_pmf"]),np.asarray(inp["rna_fl_pmf"]),cc,_debug=dbg)
chain=dbg["chain"]; cap=dbg["capture"]; st=dbg["statics"]; uni=cap["_uni"][-1]; us=cap["_uni_static"]
Gp,Gn,Rp,Rn=_oracle_per_node(inp,chain); G,R=Gp+Gn,Rp+Rn
fo=np.where(G+R>1e-9,G/np.maximum(G+R,1e-9),np.nan)
rt,_=_node_region_type(chain,ra)
left,right=np.asarray(chain.left),np.asarray(chain.right)
mass=np.asarray(cap["mass_global"]); solv=np.asarray(cap["solvable_mask"],bool)

# discrete spliced fragment COUNTS (n_pos/n_neg per face) vs fractional MASS (SP/SN per face)
npos_l=us["spl_n_pos_l"]; npos_r=us["spl_n_pos_r"]; nneg_l=us["spl_n_neg_l"]; nneg_r=us["spl_n_neg_r"]
SP_l=us["SP_l"]; SP_r=us["SP_r"]; SN_l=us["SN_l"]; SN_r=us["SN_r"]
tau_own=us["tau_own"]; struct=us["struct_lock"]
cm_p=uni["cm_p"]; cm_n=uni["cm_n"]; c_tau=uni["c_tau"]; amp=uni["amp"]; bmp=uni["bmp"]
fwd_mp=us["fwd_mp"]; bwd_mp=us["bwd_mp"]; fwd_tau=us["fwd_tau"]; bwd_tau=us["bwd_tau"]
n_unspl=us["n_unspl_l"]+us["n_unspl_r"]

print("="*90)
print("PART A — node 1909: cm_p vs the ADJACENT junction's DISCRETE spliced count (resolve verifier #2)")
print("="*90)
i=1909; L=int(left[i]); Rn_=int(right[i])
print(f"node {i}: oracle_fg={fo[i]:.3f} mass={mass[i]:,.0f} solved={uni['fg_out'][i]:.3f} fg_loc={cap['fg_loc'][i]:.3f}")
print(f"  cm_p={cm_p[i]:.3f}  (amp={amp[i]:.3f} + bmp={bmp[i]:.3f})   cm_n={cm_n[i]:.3f}")
print(f"  own spliced faces: SP_l={SP_l[i]:.3f} SP_r={SP_r[i]:.3f}  n_pos_l={npos_l[i]:.1f} n_pos_r={npos_r[i]:.1f}")
for tag,j in (("L",L),("R",Rn_)):
    if j<0: continue
    print(f"  nbr {tag} node {j}: SP_l={SP_l[j]:.3f}(n={npos_l[j]:.1f}) SP_r={SP_r[j]:.3f}(n={npos_r[j]:.1f}) "
          f"SN_l={SN_l[j]:.3f}(n={nneg_l[j]:.1f}) SN_r={SN_r[j]:.3f}(n={nneg_r[j]:.1f}) mass={mass[j]:,.0f}")
# the graft face: boundary source s grafts into exon i via SP[sf][s]; fwd uses sf=1 (src RIGHT face), bwd uses sf=0 (src LEFT face)
# fwd relay: dst i receives from left nbr s=left[i] on i's LEFT (df=0), src RIGHT face sf=1
# bwd relay: dst i receives from right nbr s=right[i] on i's RIGHT (df=1), src LEFT face sf=0
print(f"\n  GRAFT provenance of cm_p (fwd from L via src RIGHT face, bwd from R via src LEFT face):")
if L>=0: print(f"    from L={L}: SP_r[L]={SP_r[L]:.3f} discrete n_pos_r[L]={npos_r[L]:.1f}")
if Rn_>=0: print(f"    from R={Rn_}: SP_l[R]={SP_l[Rn_]:.3f} discrete n_pos_l[R]={npos_l[Rn_]:.1f}")
# max discrete count among the two adjacent graft faces
adj_disc=[]
if L>=0: adj_disc.append(npos_r[L])
if Rn_>=0: adj_disc.append(npos_l[Rn_])
adj_max=max(adj_disc) if adj_disc else 0.0
print(f"  => adjacent-junction max DISCRETE spliced-pos count = {adj_max:.1f}")
print(f"     cm_p={cm_p[i]:.3f}   cm_p <= adjacent discrete count? {cm_p[i]<=adj_max}")
print(f"     bmp={bmp[i]:.3f} = graft(R) + relayed upstream; graft(R)=damped SP_l[R]={SP_l[Rn_]:.3f}")

print("\n"+"="*90)
print("PART B — GENOME-WIDE: is cm_p count-honest vs the ADJACENT junction discrete count?")
print("="*90)
# adjacent graft discrete count ceiling per node: from left nbr's RIGHT face + right nbr's LEFT face
adj_ceiling=np.zeros_like(cm_p)
adj_mass_ceiling=np.zeros_like(cm_p)
for i in range(len(cm_p)):
    c=0.0; cm=0.0
    l=int(left[i]); r=int(right[i])
    if l>=0: c+=npos_r[l]+nneg_r[l]; cm+=SP_r[l]+SN_r[l]
    if r>=0: c+=npos_l[r]+nneg_l[r]; cm+=SP_l[r]+SN_l[r]
    adj_ceiling[i]=c; adj_mass_ceiling[i]=cm
cmpn=cm_p+cm_n
sel=solv&(cmpn>0.5)
print(f"solvable nodes with cm_p+cm_n>0.5: {sel.sum()}")
print(f"  frac(cm_p+cm_n > adjacent DISCRETE count):     {np.mean(cmpn[sel]>adj_ceiling[sel]):.3f}")
print(f"  frac(cm_p+cm_n > adjacent FRACTIONAL mass):    {np.mean(cmpn[sel]>adj_mass_ceiling[sel]+1e-9):.3f}")
print(f"  frac(cm_p+cm_n > 2x adjacent DISCRETE count):  {np.mean(cmpn[sel]>2*adj_ceiling[sel]):.3f}")
# global max single discrete spliced count anywhere
gmax_disc=max(npos_l.max(),npos_r.max(),nneg_l.max(),nneg_r.max())
print(f"  global max single DISCRETE spliced count = {gmax_disc:.1f}")
print(f"  frac(cm_p+cm_n > global max single discrete count): {np.mean(cmpn[sel]>gmax_disc):.4f}")
# foreign fraction: cm from beyond the 1-hop adjacent graft
# honest 1-hop = adjacent graft damped spliced; approximate honest as min(cmpn, adj_mass_ceiling)
honest=np.minimum(cmpn,adj_mass_ceiling)
foreign=np.maximum(cmpn-adj_mass_ceiling,0.0)
ff=np.where(cmpn>1e-9,foreign/np.maximum(cmpn,1e-9),0.0)
print(f"  mean foreign-fraction of cm_p+cm_n (beyond adjacent graft mass): {np.mean(ff[sel]):.3f}")
print(f"  frac(nodes with foreign-fraction>0.5): {np.mean(ff[sel]>0.5):.3f}")
# does the graft ALONE (adjacent, count-honest) exceed the discrete count?
print(f"  frac(adjacent graft MASS > adjacent DISCRETE count) [seed honesty]: {np.mean(adj_mass_ceiling[sel]>adj_ceiling[sel]+1e-9):.4f}")

print("\n"+"="*90)
print("PART C — c_tau (composition) stream: additive accumulation vs own? (resolve verifier#1 vs #2)")
print("="*90)
selt=solv&(c_tau>0.5)
print(f"solvable nodes with c_tau>0.5: {selt.sum()}")
print(f"  frac(c_tau > tau_own):        {np.mean(c_tau[selt]>tau_own[selt]):.3f}")
print(f"  frac(c_tau > 2x tau_own):     {np.mean(c_tau[selt]>2*tau_own[selt]):.3f}")
print(f"  frac(c_tau>1 with tau_own==0): {np.mean((c_tau>1.0)&(tau_own<=1e-9)&solv)/max(np.mean(solv),1e-9):.4f}  count={np.sum((c_tau>1.0)&(tau_own<=1e-9)&solv)}")
# is tau_own itself count-bounded? tau_lam = i_strand(~N_eff<=n) + tau_fac. compare to n_unspl count
print(f"  frac(tau_own > n_unspl count): {np.mean(tau_own[solv&(tau_own>0.5)]>n_unspl[solv&(tau_own>0.5)]):.4f}")
# max single reachable tau along relay = max over path. Bound below by max(tau_own over all) as loose proxy;
# tighter: c_tau vs global max tau_own
gmax_tau=tau_own.max()
print(f"  global max single tau_own = {gmax_tau:.4g};  max c_tau = {c_tau.max():.4g}")
print(f"  frac(c_tau > global max single tau_own): {np.mean(c_tau[selt]>gmax_tau):.4f}")
# foreign fraction of c_tau (beyond own)
tf=np.where(c_tau>1e-9,np.maximum(c_tau-tau_own,0.0)/np.maximum(c_tau,1e-9),0.0)
print(f"  mean foreign-fraction of c_tau (beyond own tau): {np.mean(tf[selt]):.3f}")
print(f"  frac(c_tau foreign-fraction>0.5): {np.mean(tf[selt]>0.5):.3f}")
i=1909
print(f"\n  node 1909: c_tau={c_tau[i]:.4g} tau_own={tau_own[i]:.4g} lam_msg f_g_eq={1/(1+np.exp(-uni['lam_msg'][i])):.3f}")

print("\n"+"="*90)
print("PART D — FUNCTIONAL ablation at node 1909 (faithful replay): which stream drives the collapse?")
print("="*90)
fp=np.asarray(st.free_pos,bool); fn=np.asarray(st.free_neg,bool)
is_amb=fp&fn
mo_g=uni["mo_g"]; mo_p=uni["mo_p"]; mo_n=uni["mo_n"]; cm_g=uni["cm_g"]; lam_msg=uni["lam_msg"]
cp=uni["cp"]; cn=uni["cn"]
tau_tilt=np.clip(np.where(cp+cn>1e-9,(cp-cn)/np.maximum(cp+cn,1e-9),0.0),-1,1)
th_msg=np.arcsin(tau_tilt); th_prec=np.where(is_amb,cm_p+cm_n,0.0)
base=dict(kappa=float(dbg["rna_sense_frac"]),
          od_g=float(dbg["calibration_priors"].gdna_strand_overdispersion),
          od_r=float(dbg["calibration_priors"].rna_strand_overdispersion),
          n_grid=cc.sweep_n_grid,L=cc.sweep_logodds_window,n_tilt=cc.sweep_n_tilt,
          n_grid_ss=cc.sweep_n_grid_single_strand,
          fg_ref=cap["fg_init"],fpos_ref=cap["fpos_init"],fneg_ref=cap["fneg_init"],
          global_logprior=cap["global_lp"],lam_logprior=cap["intron_prior"])
pos=(st.u_pos,st.u_neg,st.free_pos,st.free_neg,st.mass_unspliced,st.mass_spliced)
def solve(rna=True,gdna=True,lam=True):
    kw=dict(base)
    kw.update(gdna_imp_mode=mo_g if gdna else None, gdna_imp_prec=cm_g if gdna else None,
              rna_imp_mode=(mo_p,mo_n) if rna else None, rna_imp_prec=(cm_p,cm_n) if rna else None,
              lam_imp_mode=lam_msg if lam else None, lam_imp_prec=c_tau if lam else None,
              theta_imp_mode=th_msg,theta_imp_prec=th_prec)
    return np.asarray(_solve_nodes_logodds_all(*pos,**kw).gdna_frac)
full=solve(); fid=abs(full[1909]-uni["fg_out"][1909])
print(f"  replay fidelity at 1909: |replay-shipped|={fid:.2e}")
print(f"  oracle f_g=0.985  own fg_loc={cap['fg_loc'][1909]:.3f}")
print(f"  ALL channels:        f_g={full[1909]:.3f}")
print(f"  drop RNA-meas:       f_g={solve(rna=False)[1909]:.3f}")
print(f"  drop lambda-comp:    f_g={solve(lam=False)[1909]:.3f}")
print(f"  drop gDNA-meas:      f_g={solve(gdna=False)[1909]:.3f}")
print(f"  drop RNA+lambda:     f_g={solve(rna=False,lam=False)[1909]:.3f}")
print(f"  drop ALL msgs:       f_g={solve(rna=False,gdna=False,lam=False)[1909]:.3f}")

print("\n"+"="*90)
print("PART E — GENOME-WIDE error under stream caps (mass-weighted, solvable, error-mass)")
print("="*90)
ok=np.isfinite(fo)&(mass>1e-9)&solv
w=mass[ok]
def mwae(fgv): return np.average(np.abs(fgv[ok]-fo[ok]),weights=w)
# cap cm_p/cm_n to adjacent discrete count; cap c_tau to tau_own
cm_p_cap=np.minimum(cm_p,adj_ceiling); cm_n_cap=np.minimum(cm_n,np.maximum(adj_ceiling-cm_p,0.0))
def solve_full(cmp_,cmn_,ctau_):
    kw=dict(base)
    kw.update(gdna_imp_mode=mo_g,gdna_imp_prec=cm_g,rna_imp_mode=(mo_p,mo_n),rna_imp_prec=(cmp_,cmn_),
              lam_imp_mode=lam_msg,lam_imp_prec=ctau_,theta_imp_mode=th_msg,
              theta_imp_prec=np.where(is_amb,cmp_+cmn_,0.0))
    return np.asarray(_solve_nodes_logodds_all(*pos,**kw).gdna_frac)
print(f"  baseline (shipped):                 mwae={mwae(solve_full(cm_p,cm_n,c_tau)):.4f}")
print(f"  cap RNA-meas to adjacent count:     mwae={mwae(solve_full(cm_p_cap,cm_n_cap,c_tau)):.4f}")
print(f"  cap c_tau to own:                   mwae={mwae(solve_full(cm_p,cm_n,np.minimum(c_tau,tau_own))):.4f}")
print(f"  cap BOTH:                           mwae={mwae(solve_full(cm_p_cap,cm_n_cap,np.minimum(c_tau,tau_own))):.4f}")
