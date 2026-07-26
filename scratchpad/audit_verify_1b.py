"""VERIFIER #1 — PART 7: is cm_p the WHOLE story, or do the composition-lambda and gDNA-meas streams
carry the SAME additive-accumulation / undamped-reframe defect? Multi-channel ablation at node 1909 +
a genuine SHOULD-NOT case (a node whose message precision IS backed by its own discrete spliced count)."""
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
os.environ.pop("RIGEL_S2T_OFF", None)
cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
dbg = {}
calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                 np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
chain = dbg["chain"]; cap = dbg["capture"]; uni = cap["_uni"][-1]; st = cap["_uni_static"]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
fo = np.where(G + R > 1e-9, G / np.maximum(G + R, 1e-9), np.nan)
stt = dbg["statics"]; fp = st["fp"]; fn = st["fn"]
mo_g, mo_p, mo_n = uni["mo_g"], uni["mo_p"], uni["mo_n"]
cm_g, cm_p, cm_n, c_tau, lam_msg = uni["cm_g"], uni["cm_p"], uni["cm_n"], uni["c_tau"], uni["lam_msg"]
own_spl_p = st["SP_l"] + st["SP_r"]; own_spl_n = st["SN_l"] + st["SN_r"]
mass = st["M"]; tau_own = st["tau_own"]

kw = dict(kappa=float(dbg["rna_sense_frac"]),
          od_g=float(dbg["calibration_priors"].gdna_strand_overdispersion),
          od_r=float(dbg["calibration_priors"].rna_strand_overdispersion),
          n_grid=int(cc.sweep_n_grid), L=float(cc.sweep_logodds_window),
          n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand)

def solve_one(i, use_g=True, use_rna=True, use_lam=True):
    return float(_solve_nodes_logodds_all(
        np.asarray([stt.u_pos[i]]), np.asarray([stt.u_neg[i]]),
        np.asarray([fp[i]]), np.asarray([fn[i]]),
        np.asarray([stt.mass_unspliced[i]]), np.asarray([stt.mass_spliced[i]]),
        global_logprior=None,
        gdna_imp_mode=np.asarray([mo_g[i]]) if use_g else None,
        gdna_imp_prec=np.asarray([cm_g[i]]) if use_g else None,
        rna_imp_mode=(np.asarray([mo_p[i]]), np.asarray([mo_n[i]])) if use_rna else None,
        rna_imp_prec=(np.asarray([cm_p[i]]), np.asarray([cm_n[i]])) if use_rna else None,
        lam_imp_mode=np.asarray([lam_msg[i]]) if use_lam else None,
        lam_imp_prec=np.asarray([c_tau[i]]) if use_lam else None,
        fg_ref=np.asarray([cap["fg_init"][i]]),
        fpos_ref=np.asarray([cap["fpos_init"][i]]),
        fneg_ref=np.asarray([cap["fneg_init"][i]]), **kw).gdna_frac[0])

i = 1909
print(f"node 1909 oracle_fg={fo[i]:.3f}  own message-free fg_loc={cap['fg_loc'][i]:.3f}")
print(f"  own tau0={tau_own[i]:.3f}  -> ACCUMULATED c_tau={c_tau[i]:.3f}  (lam mode fg_eq={1/(1+np.exp(-lam_msg[i])):.3f})")
print(f"  cm_g={cm_g[i]:.2f} (mode fg={np.exp(mo_g[i]):.3f})   cm_p={cm_p[i]:.2f} (mode fpos={np.exp(mo_p[i]):.3f})")
print("  channel ablation (which messages collapse f_g?):")
print(f"    all three channels ON     f_g={solve_one(i,1,1,1):.3f}")
print(f"    drop RNA-meas             f_g={solve_one(i,1,0,1):.3f}")
print(f"    drop LAMBDA-comp          f_g={solve_one(i,1,1,0):.3f}")
print(f"    drop gDNA-meas            f_g={solve_one(i,0,1,1):.3f}")
print(f"    drop RNA + LAMBDA         f_g={solve_one(i,1,0,0):.3f}")
print(f"    drop ALL messages (own)   f_g={solve_one(i,0,0,0):.3f}")

# c_tau accumulation genome-wide: does the COMPOSITION precision exceed a node's OWN tau0?
solv = np.asarray(cap["solvable_mask"], bool)
m = solv & (c_tau > 1e-6)
print(f"\ncomposition c_tau vs own tau0 (accumulation check), solvable nodes with c_tau>0: {int(m.sum())}")
print(f"  frac( c_tau > own tau0 )          = {float((c_tau[m] > tau_own[m] + 1e-9).mean()):.3f}")
print(f"  frac( c_tau > 2x own tau0 )       = {float((c_tau[m] > 2*tau_own[m] + 1e-9).mean()):.3f}")
zero_tau = solv & (c_tau > 1.0) & (tau_own < 1e-6)
print(f"  nodes with c_tau>1 but ZERO own tau0 (pure imported composition precision) = {int(zero_tau.sum())}")

# SHOULD-NOT: nodes whose RNA-meas precision IS backed by their OWN spliced count. Expect solved~oracle,
# and cm_p not wildly above own count.
fg_out = uni["fg_out"]
good = solv & (own_spl_p > 10.0) & np.isfinite(fo)
print(f"\nSHOULD-NOT (own spliced count>10, precision self-backed): {int(good.sum())} nodes")
if good.sum():
    print(f"  median |solved-oracle| = {np.median(np.abs(fg_out[good]-fo[good])):.3f}")
    print(f"  median cm_p/own_spl_p  = {np.median(cm_p[good]/np.maximum(own_spl_p[good],1e-9)):.2f}")
    print(f"  (a self-backed node: precision ~ its own count, mode correct -> solve tracks oracle)")
