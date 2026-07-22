"""Owner's question: at an intron|exon boundary, is σ²_transfer lower on the EXON→boundary edge (dense, cliff-
crossed) than the INTRON→boundary edge (clean but sparse/depleted)? And does the intron edge degrade with better
capture? σ²_transfer(src→dst) = var_proj[dst] + (mu_proj[dst] − mu_proj[src])². Split by capture regime; only
intron|exon boundaries (one intron flank, one exon flank)."""
from __future__ import annotations
import dataclasses
import importlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.bp_solver import REGION
from rigel.calibration.node_geometry import _node_region_type
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and "0.50" in d.name)


def capreg(name):
    if "verystrong" in name:
        return "verystrong"
    if "capture_on" in name:
        return "on"
    return "off"


rows = {"off": {"exon": [], "intron": [], "exon_m": [], "intron_m": []},
        "on": {"exon": [], "intron": [], "exon_m": [], "intron_m": []},
        "verystrong": {"exon": [], "intron": [], "exon_m": [], "intron_m": []}}
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    node_type, _ = _node_region_type(chain, ra)
    mu, vp = np.asarray(cap["_mu_proj"]), np.asarray(cap["_var_proj"])
    massf = (np.asarray(cap["mass_l"]), np.asarray(cap["mass_r"]))  # source facing mass, per side
    kind = np.asarray(chain.kind); left, right = np.asarray(chain.left), np.asarray(chain.right)
    reg = capreg(cond)
    for i in np.where(kind != REGION)[0]:
        L, R = int(left[i]), int(right[i])
        if L < 0 or R < 0:
            continue
        lt, rt = int(node_type[L]), int(node_type[R])
        # identify the intron flank and exon flank (need exactly one of each)
        flanks = {}
        if lt == 1:
            flanks["intron"] = L
        elif lt == 2:
            flanks["exon"] = L
        if rt == 1:
            flanks["intron"] = R
        elif rt == 2:
            flanks["exon"] = R
        if "intron" not in flanks or "exon" not in flanks:
            continue
        for tag, f in flanks.items():
            s2t = vp[i] + (mu[i] - mu[f]) ** 2
            rows[reg][tag].append(s2t)
            rows[reg][tag + "_m"].append(cap["mass_global"][f])  # flank mass (sparsity of the intron source)

print("σ²_transfer at intron|exon boundaries, by capture regime "
      "(lower = the message is trusted more):\n")
print(f"{'capture':>12} | {'exon→bnd σ²_T':>26} | {'intron→bnd σ²_T':>26} | {'intron flank mass':>18}")
print(f"{'':>12} | {'median  [p25–p75]':>26} | {'median  [p25–p75]':>26} | {'median':>18}")
for reg in ("off", "on", "verystrong"):
    ex = np.array(rows[reg]["exon"]); it = np.array(rows[reg]["intron"]); im = np.array(rows[reg]["intron_m"])
    if len(ex) < 5:
        continue
    def q(a):
        return f"{np.median(a):>6.2f}  [{np.percentile(a,25):.2f}–{np.percentile(a,75):.2f}]"
    print(f"{reg:>12} | {q(ex):>26} | {q(it):>26} | {np.median(im):>18.1f}  (n={len(ex)})")
print("\nInterpretation: owner expects exon→bnd σ²_T < intron→bnd σ²_T (exon more trusted), and the intron flank"
      "\nmass to SHRINK as capture strengthens (off→on→verystrong) — the depleted-intron sparsity he flagged.")
