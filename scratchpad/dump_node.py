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
chain=dbg["chain"]; cap=dbg["capture"]; uni=cap["_uni"][-1]
Gp,Gn,Rp,Rn=_oracle_per_node(inp,chain); G,R=Gp+Gn,Rp+Rn
fo=np.where(G+R>1e-9,G/np.maximum(G+R,1e-9),np.nan)
rt,_=_node_region_type(chain,ra); CLS={0:"intergenic",1:"intron",2:"exon",-1:"boundary"}
fp=np.asarray(dbg["statics"].free_pos,bool); fn=np.asarray(dbg["statics"].free_neg,bool)
tau0=np.asarray(cap["_tau0_lam"]); fgloc=np.asarray(cap["fg_loc"])
left,right=np.asarray(chain.left),np.asarray(chain.right)
mass=np.asarray(cap["mass_global"])
def show(i):
    cl=CLS[int(rt[i])] if np.asarray(chain.kind)[i]==REGION else "boundary"
    print(f"\nnode {i} [{cl}] mass={mass[i]:,.0f} oracle_fg={fo[i]:.3f} gDNA={G[i]:,.0f} RNA={R[i]:,.0f}")
    print(f"   own(_ni): fg_loc={fgloc[i]:.3f} tau0_lam={tau0[i]:.4g} free=({int(fp[i])},{int(fn[i])})")
    print(f"   SOLVED f_g={uni['fg_out'][i]:.3f}")
    print(f"   λ-MSG   mode(f_g_eq)={1/(1+np.exp(-uni['lam_msg'][i])):.3f} prec(c_tau)={uni['c_tau'][i]:.4g}")
    print(f"   gDNA-meas mode={np.exp(uni['mo_g'][i]):.3f} prec(cm_g)={uni['cm_g'][i]:.4g}   full cpg={uni['cpg'][i]:.4g}")
    print(f"   RNA-meas +mode={np.exp(uni['mo_p'][i]):.3f} prec(cm_p)={uni['cm_p'][i]:.4g}  -mode={np.exp(uni['mo_n'][i]):.3f} prec(cm_n)={uni['cm_n'][i]:.4g}")
    for tag,j in (("L",int(left[i])),("R",int(right[i]))):
        if j>=0:
            clj=CLS[int(rt[j])] if np.asarray(chain.kind)[j]==REGION else "boundary"
            print(f"   nbr {tag} node {j} [{clj}] oracle={fo[j]:.3f} solved={uni['fg_out'][j]:.3f} mass={mass[j]:,.0f}")
for i in (1909,1908,1910,1087):
    show(i)
