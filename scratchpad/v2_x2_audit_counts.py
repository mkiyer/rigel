"""Re-harvest x2's edges but ALSO carry the true COUNTS (spliced count arrays; unspliced count/mass ratio)."""
from __future__ import annotations
import dataclasses, importlib, sys
from pathlib import Path
import numpy as np
from scipy.special import polygamma
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
from rigel.calibration.bp_solver import REGION
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS, _TINY = 1e-9, 1e-12
CONDS = ["gdna_gdna300_ss_0.99_nrna_present_capture_off",
         "gdna_gdna300_ss_0.50_nrna_present_capture_off",
         "gdna_gdna100_ss_0.50_nrna_present_capture_off",
         "gdna_gdna300_ss_0.50_nrna_none_capture_off"]
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig(); ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
psi1 = lambda n: polygamma(1, np.maximum(np.asarray(n,float), _TINY))

def harvest(cond):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE/"_selfsolve_cache")
    dbg={}
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]),
                     dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us = cap["_uni_static"]
    kind=np.asarray(chain.kind); idx=np.asarray(chain.ref_idx,np.int64); isr=kind==REGION
    def pool(k):
        a=np.asarray(inp["region_pools"][k],float); b=np.asarray(inp["boundary_pools"][k],float)
        return np.where(isr, a[np.clip(idx,0,a.shape[0]-1)], b[np.clip(idx,0,b.shape[0]-1)])
    G=pool("gdna_pos")+pool("gdna_neg")
    Ru=pool("mat_uns_pos")+pool("nas_uns_pos")+pool("mat_uns_neg")+pool("nas_uns_neg")
    Rs=pool("mat_spl")+pool("nas_spl")
    E_g=np.asarray(us["E_g"],float); E_r=np.asarray(us["E_r"],float)
    li,ri=us["left"],us["right"]; is_bnd,is_exon=us["is_bnd"],us["is_exon"]
    SPf,SNf=(us["SP_l"],us["SP_r"]),(us["SN_l"],us["SN_r"])
    SPn,SNn=(us["spl_n_pos_l"],us["spl_n_pos_r"]),(us["spl_n_neg_l"],us["spl_n_neg_r"])
    ESP=(np.asarray(geo.eff_spl_left,float),np.asarray(geo.eff_spl_right,float))
    # node unspliced MASS and COUNT (mass_left/right + n_unspl_left/right, region=contained, bnd=both faces)
    ml,mr=np.asarray(geo.mass_left,float),np.asarray(geo.mass_right,float)
    nl,nr=np.asarray(geo.n_unspl_left,float),np.asarray(geo.n_unspl_right,float)
    M_un=np.where(isr, ml, ml+mr); N_un=np.where(isr, nl, nl+nr)
    cscale=np.where(M_un>_EPS, N_un/np.maximum(M_un,_EPS), 1.0)   # count per unit mass at this node
    rho_g=np.where(E_g>_EPS,G/np.maximum(E_g,_EPS),np.nan)
    rho_R=np.where(E_r>_EPS,(Ru+Rs)/np.maximum(E_r,_EPS),np.nan)
    rho_nu=np.where(E_r>_EPS,Ru/np.maximum(E_r,_EPS),np.nan)
    rec={k:[] for k in ("phi","n_spl_mass","n_spl_cnt","n_unspl_R_bnd","c_bnd","n_R_exon","c_exon",
                        "n_g_bnd","n_g_exon","w_nu","w_mu","rho_mu","rho_R","share")}
    for face,nbr in ((1,li),(0,ri)):
        for i in np.flatnonzero(is_exon):
            b=nbr[i]
            if b<0 or not is_bnd[b]: continue
            sp,sn=float(SPf[face][b]),float(SNf[face][b]); espl=float(ESP[face][b])
            mu=(sp+sn)/max(espl,_EPS)
            if not (mu>_EPS) or not np.isfinite(rho_R[i]) or rho_R[i]<=_EPS: continue
            if not (rho_g[b]>_EPS and rho_g[i]>_EPS): continue
            step=rho_g[b]/rho_g[i]; phi=(rho_nu[b]+mu)/(rho_R[i]*step)
            if not (np.isfinite(phi) and phi>0): continue
            tot=rho_nu[b]+mu
            rec["phi"].append(phi); rec["n_spl_mass"].append(sp+sn)
            rec["n_spl_cnt"].append(float(SPn[face][b])+float(SNn[face][b]))
            rec["n_unspl_R_bnd"].append(float(Ru[b])); rec["c_bnd"].append(float(cscale[b]))
            rec["n_R_exon"].append(float(Ru[i]+Rs[i])); rec["c_exon"].append(float(cscale[i]))
            rec["n_g_bnd"].append(float(G[b])); rec["n_g_exon"].append(float(G[i]))
            rec["w_nu"].append(rho_nu[b]/tot); rec["w_mu"].append(mu/tot)
            rec["rho_mu"].append(mu); rec["rho_R"].append(float(rho_R[i]))
            rec["share"].append((sp+sn)/max(float(Ru[i]+Rs[i]),_EPS))
    return {k:np.asarray(v,float) for k,v in rec.items()}

D={c:harvest(c) for c in CONDS}
print("\n=== mass vs COUNT at the junction face (project: 'v_mu uses the spliced COUNT, never the mass') ===")
print(f"{'condition':<48}{'med mass':>10}{'med count':>11}{'med cnt/mass':>14}{'mean cnt/mass':>15}")
for c in CONDS:
    d=D[c]; r=d["n_spl_cnt"]/np.maximum(d["n_spl_mass"],_TINY)
    print(f"{c[5:]:<48}{np.median(d['n_spl_mass']):>10.2f}{np.median(d['n_spl_cnt']):>11.1f}"
          f"{np.median(r):>14.2f}{np.mean(r):>15.2f}")
print("\n=== unspliced count/mass scale at boundary and exon nodes ===")
for c in CONDS:
    d=D[c]; print(f"{c[5:]:<48} bnd med={np.median(d['c_bnd']):.2f}  exon med={np.median(d['c_exon']):.2f}")

def budget(d, use_counts):
    ns = d["n_spl_cnt"] if use_counts else d["n_spl_mass"]
    nu = d["n_unspl_R_bnd"]*(d["c_bnd"] if use_counts else 1.0)
    nR = d["n_R_exon"]*(d["c_exon"] if use_counts else 1.0)
    ng_b = d["n_g_bnd"]*(d["c_bnd"] if use_counts else 1.0)
    ng_e = d["n_g_exon"]*(d["c_exon"] if use_counts else 1.0)
    vnum = np.where(d["w_nu"]>0, d["w_nu"]**2/np.maximum(nu,_TINY),0.0)+d["w_mu"]**2/np.maximum(ns,_TINY)
    return vnum+psi1(nR)+psi1(ng_b)+psi1(ng_e)

print("\n=== RESIDUAL: x2's MASS-as-count vs the true COUNTS (EXACT form) ===")
print(f"{'condition':<48}{'Var':>9}{'E[vp] mass':>12}{'RESID mass':>12}{'E[vp] cnt':>11}{'RESID cnt':>11}")
allv=[];allm=[];allc=[]
for c in CONDS:
    d=D[c]; lp=np.log(d["phi"]); V=float(np.var(lp))
    bm=float(np.mean(budget(d,False))); bc=float(np.mean(budget(d,True)))
    allv.append(lp); allm.append(budget(d,False)); allc.append(budget(d,True))
    print(f"{c[5:]:<48}{V:>9.4f}{bm:>12.4f}{V-bm:>12.4f}{bc:>11.4f}{V-bc:>11.4f}")
LP=np.concatenate(allv); BM=np.concatenate(allm); BC=np.concatenate(allc)
print(f"{'POOLED':<48}{np.var(LP):>9.4f}{BM.mean():>12.4f}{np.var(LP)-BM.mean():>12.4f}{BC.mean():>11.4f}{np.var(LP)-BC.mean():>11.4f}")

print("\n=== SHAPE by junction spliced COUNT (not mass), pooled cap-OFF ===")
NS=np.concatenate([D[c]["n_spl_cnt"] for c in CONDS])
for a,b in ((-np.inf,30),(30,100),(100,300),(300,1000),(1000,np.inf)):
    m=(NS>=a)&(NS<b)
    if m.sum()<3: continue
    print(f"  n_spl_cnt {a:g}-{b:g}: n={m.sum():5d} med={np.median(NS[m]):8.1f} Var={np.var(LP[m]):8.4f} "
          f"E[vp_cnt]={BC[m].mean():8.4f} RESID={np.var(LP[m])-BC[m].mean():8.4f}")

print("\n=== SHAPE by the ACTUAL SHARE n_spl/n_R(exon) (the quantity the headline names) ===")
SH=np.concatenate([D[c]["share"] for c in CONDS])
for a,b in ((-np.inf,0.05),(0.05,0.15),(0.15,0.3),(0.3,0.6),(0.6,np.inf)):
    m=(SH>=a)&(SH<b)
    if m.sum()<3: continue
    print(f"  share {a:g}-{b:g}: n={m.sum():5d} med={np.median(SH[m]):6.3f} Var={np.var(LP[m]):8.4f} "
          f"E[vp]={BC[m].mean():8.4f} RESID={np.var(LP[m])-BC[m].mean():8.4f}")
np.savez("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/df34908c-f6e7-45d0-827a-31f2a818e9d9/scratchpad/a2.npz",
         **{f"{c}__{k}":v for c in CONDS for k,v in D[c].items()})
