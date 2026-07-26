import os, sys, dataclasses
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from pathlib import Path
from flagship_interrogate import _oracle_per_node
from selfsolve_diag import _scan_and_truth
from rigel.calibration.bp_solver import REGION
import importlib; calmod=importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
SUITE=Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND="gdna_gdna300_ss_0.99_nrna_present_capture_on"
index=TranscriptIndex.load(str(SUITE/"rigel_index")); cfg=PipelineConfig()
ra=RegionArrays.from_region_df(index.region_df,index.ref_name_to_id)
inp=_scan_and_truth(SUITE,COND,index,cfg,Path("/tmp/rigel_selfsolve"),SUITE/"_selfsolve_cache")
def run(s2t_off):
    os.environ.pop("RIGEL_S2T_OFF",None)
    if s2t_off: os.environ["RIGEL_S2T_OFF"]="1"
    dbg={}
    cc=dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"],ra,inp["strand_model"],np.asarray(inp["gdna_fl_pmf"]),np.asarray(inp["rna_fl_pmf"]),cc,_debug=dbg)
    return dbg
dbg_off=run(True); dbg_on=run(False)
chain=dbg_on["chain"]
Gp,Gn,Rp,Rn=_oracle_per_node(inp,chain); G,R=Gp+Gn,Rp+Rn
fo=np.where(G+R>1e-9,G/np.maximum(G+R,1e-9),np.nan)
fg_off=np.asarray(dbg_off["capture"]["f_g"]); fg_on=np.asarray(dbg_on["capture"]["f_g"])
mass=np.asarray(dbg_on["capture"]["mass_global"])
is_reg=np.asarray(chain.kind)==REGION
ok=np.isfinite(fo)&(mass>1e-9)
# regression = |on-oracle| - |off-oracle|, error-mass-weighted
reg=(np.abs(fg_on-fo)-np.abs(fg_off-fo))
regmass=np.where(ok,reg*mass,0.0)
print(f"COND={COND}")
for tag,m in (("ALL",ok),("region",ok&is_reg),("boundary",ok&~is_reg)):
    w=mass[m]
    e_off=np.average(np.abs(fg_off[m]-fo[m]),weights=w); e_on=np.average(np.abs(fg_on[m]-fo[m]),weights=w)
    print(f"  {tag:<9} mwae σ²t-OFF={e_off:.4f}  σ²t-ON={e_on:.4f}  Δ={e_on-e_off:+.4f}")
order=np.argsort(-regmass)[:10]
print("\nWORST 10 REGRESSED NODES (σ²t ON vs OFF):")
left,right=np.asarray(chain.left),np.asarray(chain.right)
for i in order:
    i=int(i)
    if regmass[i]<=0: continue
    cls="reg" if is_reg[i] else "BND"
    print(f"  node {i} [{cls}] mass={mass[i]:,.0f} oracle={fo[i]:.3f}  f_g OFF={fg_off[i]:.3f} ON={fg_on[i]:.3f}"
          f"  Δregress={regmass[i]:,.0f}  nbrs L={left[i]} R={right[i]}")
