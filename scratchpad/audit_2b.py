"""audit_2b — decompose cm_p[1909]: immediate graft neighbour vs ACCUMULATED foreign junctions, and confirm
the graft sigma2_transfer exemption keeps it undamped."""
import os, sys, dataclasses, importlib
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from pathlib import Path
from selfsolve_diag import _scan_and_truth
from rigel.calibration.bp_solver import REGION
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND = "gdna_gdna300_ss_0.99_nrna_present_capture_on"
index = TranscriptIndex.load(str(SUITE / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")

def run(s2t_off):
    os.environ.pop("RIGEL_S2T_OFF", None)
    if s2t_off: os.environ["RIGEL_S2T_OFF"] = "1"
    dbg = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return dbg

for s2t_off in (False, True):
    dbg = run(s2t_off)
    cap = dbg["capture"]; uni = cap["_uni"][-1]; st = cap["_uni_static"]
    left = np.asarray(dbg["chain"].left); right = np.asarray(dbg["chain"].right)
    cm_p = np.asarray(uni["cm_p"]); amp = np.asarray(uni["amp"]); bmp = np.asarray(uni["bmp"])
    SP_l = np.asarray(st["SP_l"]); SP_r = np.asarray(st["SP_r"])
    fwd_mp = np.asarray(st["fwd_mp"]); bwd_mp = np.asarray(st["bwd_mp"])
    i = 1909
    tag = "S2T OFF" if s2t_off else "S2T ON "
    print(f"[{tag}] node 1909: cm_p={cm_p[i]:.3f} = amp(left)={amp[i]:.3f} + bmp(right)={bmp[i]:.3f}")
    print(f"         fwd_mp[1909]={fwd_mp[i]:.3f}  bwd_mp[1909]={bwd_mp[i]:.3f}  (relay accumulators pre-transport)")
    # immediate graft: node 1909 exon receives graft from a boundary neighbour. fwd source=left[1909], its
    # emitting face is RIGHT (SP_r); bwd source=right[1909], face LEFT (SP_l).
    L = int(left[i]); R = int(right[i])
    print(f"         left nbr {L}: SP_r(emit face)={SP_r[L]:.3f} | right nbr {R}: SP_l(emit face)={SP_l[R]:.3f}")
    # count how many distinct upstream junctions inject into the bwd accumulator reaching 1909
    if not s2t_off:
        j = R; hops = 0; injectors = []
        seen = {i}
        while j >= 0 and j not in seen and hops < 4000:
            seen.add(j)
            face = SP_l[j]  # bwd source emits its LEFT face
            if face > 0.5:
                injectors.append((j, face, hops))
            j = int(right[j]); hops += 1
        print(f"         bwd relay into 1909 spans {hops} hops; junctions with SP_l>0.5 injecting: {len(injectors)}")
        tot = sum(x[1] for x in injectors)
        print(f"         sum of those injected spliced masses = {tot:.1f} (vs immediate neighbour {SP_l[R]:.1f});"
              f" nearest 5: {[(jj,round(f,1),h) for jj,f,h in injectors[:5]]}")
    print()
