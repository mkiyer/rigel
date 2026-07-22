"""DERIVE the boundary RNA imputation: for the exon→boundary hop, is the residual RNA density
    rho_r = rna_src − mat_abs  =  (exon RNA density) − (boundary junction-mature density)
collapsing to ~0 and driving f_g→1, and is the cause a FRAME/SCALE mismatch between the two densities (ERs, the
RNA contained eff-len, vs ESPd, the spliced half-triangle)? Traces a boundary node's incoming edges across the 20
unstranded scenarios: the two RNA terms, their ratio, the residual, the eff-lengths, vs the boundary ORACLE f_g."""
from __future__ import annotations
import dataclasses
import importlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
from flagship_interrogate import _oracle_per_node
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.bp_solver import REGION
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
NODE = int(sys.argv[1]) if len(sys.argv) > 1 else 3236
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and "0.50" in d.name)

def _onode(inp, chain, key):
    """per chain-node oracle POOL (mass) for a raw pool key, e.g. nas_uns_pos."""
    kind = np.asarray(chain.kind); idx = np.asarray(chain.ref_idx, np.int64); isr = kind == REGION
    rp = np.asarray(inp["region_pools"][key], float); bp = np.asarray(inp["boundary_pools"][key], float)
    ri = np.clip(idx, 0, rp.shape[0] - 1); bi = np.clip(idx, 0, bp.shape[0] - 1)
    return np.where(isr, rp[ri], bp[bi])


print(f"node {NODE} — exon→boundary RNA imputation vs ORACLE nascent/mature at the boundary\n")
print(f"{'scenario':>30} | {'bnd_or':>6} | {'or_nas':>6} {'or_mat':>6} | {'rna_src':>7} {'mat_abs':>7} | "
      f"{'rho_r':>7} {'rho_g':>7} | {'implF_g':>7}")
agg = {"ratio": [], "collapsed": 0, "n": 0}
for cond in sorted(conds):
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    or_nas = _onode(inp, chain, "nas_uns_pos") + _onode(inp, chain, "nas_uns_neg")
    or_mat = _onode(inp, chain, "mat_uns_pos") + _onode(inp, chain, "mat_uns_neg")
    edges = [e for e in cap.get("_edge_modes", []) if e["dst"] == NODE and e["mat_abs"] > _EPS]  # exon-facing hop
    tag = cond.replace("gdna_", "").replace("_ss_0.50_nrna_present", "").replace("_ss_0.50_nrna_none", "").replace("_capture", "/cap")
    for e in edges:
        rho_r = e["rho_pos"] + e["rho_neg"]
        Mg = e["rho_g"] * e["egd"]; Mr = max(rho_r, 0.0) * e["erd"]
        implf = Mg / max(Mg + Mr, _EPS)  # SHIFT-implied f_g (imputed total, the theory-correct rule)
        agg["ratio"].append(e["mat_abs"] / max(e["rna_src"], _EPS)); agg["n"] += 1
        agg["collapsed"] += int(rho_r <= e["rna_src"] * 0.1)
        print(f"{tag:>30} | {fo[NODE]:>6.2f} | {or_nas[NODE]:>6.1f} {or_mat[NODE]:>6.1f} | "
              f"{e['rna_src']:>7.2f} {e['mat_abs']:>7.2f} | {rho_r:>7.2f} {e['rho_g']:>7.3f} | {implf:>7.2f}")
if agg["n"]:
    r = np.array(agg["ratio"])
    print(f"\nmat_abs / rna_src  ratio: median {np.median(r):.2f}  mean {r.mean():.2f}  "
          f"[{np.percentile(r,10):.2f}–{np.percentile(r,90):.2f}]   collapsed(res<10%): {agg['collapsed']}/{agg['n']}")
    print("If the ratio ≫ 1 systematically, the mature term OVER-scales the RNA — a frame/scale bug, not real mature.")
