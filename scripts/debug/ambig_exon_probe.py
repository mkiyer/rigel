"""Decompose the AMBIG-exon under-call, node by node, vs oracle.

For the top AMBIG-exon regions (production under-calls gDNA vs oracle), dump the per-node solve
decomposition: strand-only f_g, local (strand+global prior) f_g, final f_g, and the proximity
density-prior internals (x_obs, x*, σ_node, ρ_floor). Tests the circular-trap hypothesis: a low
pass-1 f_g → low x_obs → proximity maps to the depleted mode → stuck low.

    OMP_NUM_THREADS=1 python scripts/debug/ambig_exon_probe.py [reg1,reg2,...]
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402
import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
targets = [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else
                            ["242", "1080", "672", "224", "1083", "231"])]

index, blob = build_or_load_cache(COND, False)
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

CAP = {}
PASSES = []
CHAIN = {}
_orig = bp.node_sweep
def _wrapped(*a, **k):
    d = {}
    k["_capture"] = d
    CHAIN["chain"] = a[0] if a else k.get("chain")
    r = _orig(*a, **k)
    PASSES.append(d)
    CAP.update(d)
    return r
bp.node_sweep = _wrapped
import rigel.calibration.calibrate  # noqa: F401  (ensure submodule is in sys.modules)
_calmod = sys.modules["rigel.calibration.calibrate"]  # the MODULE (attr `calibrate` is the func)
_calmod.node_sweep = _wrapped

cal = calibrate(payload=blob["payload_full"], region_arrays=ra,
                strand_model=blob["strand_full"], gdna_fl_pmf=blob["gdna_pmf"],
                rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())

chain = CHAIN["chain"]
kind = np.asarray(chain.kind); rref = np.asarray(chain.ref_idx)
# region index -> node id
reg_to_node = {}
for nid in range(len(kind)):
    if kind[nid] == REGION:
        reg_to_node[int(rref[nid])] = nid

# oracle per region
sub_g = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
sub_r = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
g_or = np.asarray(sub_g.contained.mass_unspliced, float)
r_or = np.asarray(sub_r.contained.mass_unspliced, float)
fg_or = g_or / (g_or + r_or + 1e-9)

def g(key):
    v = CAP.get(key)
    return None if v is None else np.asarray(v)

# ── stage-by-stage gDNA mass loss across ALL AMBIG-exon region nodes ──
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature
sig = index.region_df["signature"].to_numpy()
scls = np.array([coarse_strand_from_signature(int(s)) for s in sig])
tcls = np.array([coarse_type_from_signature(int(s)) for s in sig])
ambig_exon = [r for r in range(len(sig)) if scls[r] == 3 and tcls[r] == 2 and r in reg_to_node]
nodes = np.array([reg_to_node[r] for r in ambig_exon])
regs = np.array(ambig_exon)
M = np.asarray(CAP["mass_global"])[nodes]        # unspliced mass per node
g_or_n = g_or[regs]                               # oracle gDNA frac
P0, P1 = PASSES[0], PASSES[1]
stage = {
    "oracle":      g_or_n,
    "strand":      np.asarray(P0["fg_strand"])[nodes],
    "P0_loc(+floor)": np.asarray(P0["fg_loc"])[nodes],
    "P0_final(+msg)": np.asarray(P0["f_g"])[nodes],
    "P1_loc(+KDE)":   np.asarray(P1["fg_loc"])[nodes],
    "P1_final":       np.asarray(P1["f_g"])[nodes],
}
print(f"=== {len(regs)} AMBIG-exon nodes: mass-weighted gDNA (mass sum={M.sum():,.0f}) ===")
print(f"  oracle gDNA mass = {(M*g_or_n).sum():,.0f}")
for name, fgv in stage.items():
    gm = (M * fgv).sum()
    print(f"  {name:>16}: mean_fg={ (M*fgv).sum()/M.sum():.3f}   gDNA_mass={gm:>12,.0f}   vs_oracle={gm-(M*g_or_n).sum():>+12,.0f}")
print()

print(f"n_passes captured={len(PASSES)}   floor: rho_global={g('rho_global')}")
print("=== per-pass fg_loc / f_g for each target node (pass1=floor-only, pass2=+KDE) ===")
for reg in targets:
    nid = reg_to_node.get(reg)
    if nid is None: continue
    cells = []
    for pi, d in enumerate(PASSES):
        fl = np.asarray(d.get("fg_loc"))[nid]
        ff = np.asarray(d.get("f_g"))[nid]
        fs = np.asarray(d.get("fg_strand"))[nid]
        cells.append(f"P{pi}[strand={fs:.3f} loc={fl:.3f} final={ff:.3f}]")
    print(f"  reg {reg:>5} (or {fg_or[reg]:.3f}): " + "  ".join(cells))
print()
hdr = ["reg", "node", "fg_or", "fg_strand", "fg_loc", "fg_final",
       "M", "E", "M/E", "rho_obs", "x_obs", "x*", "sig_node", "var*", "rho_floor", "in_floor"]
print("  ".join(f"{h:>8}" for h in hdr))
for reg in targets:
    nid = reg_to_node.get(reg)
    if nid is None:
        print(f"{reg}: no node"); continue
    def val(key, i=nid):
        v = g(key)
        if v is None: return float("nan")
        v = np.asarray(v)
        return float(v[i]) if v.ndim else float(v)
    M = val("mass_global"); E = val("eff_global")
    rf = g("rho_floor"); fm = g("floor_mask")
    rho_floor = (float(np.asarray(rf)[nid]) if (rf is not None and np.asarray(rf).ndim) else
                 (float(rf) if rf is not None else float("nan")))
    in_floor = (bool(np.asarray(fm)[nid]) if fm is not None else False)
    # precision breakdown + message anatomy: neighbors, their belief, msg mode/prec
    vgl = g("vg_loc"); pg = g("prec_g"); mg = g("mode_g")
    ch = CHAIN["chain"]; lft = int(np.asarray(ch.left)[nid]); rgt = int(np.asarray(ch.right)[nid])
    fgf = g("f_g"); rho_g_all = None
    def nb(x):
        return "-" if x < 0 else (f"n{x}[kind={int(np.asarray(ch.kind)[x])} "
                f"fg={float(np.asarray(fgf)[x]):.3f} M/E={float(np.asarray(g('mass_global'))[x])/max(float(np.asarray(g('eff_global'))[x]),1e-9):.2f}]")
    print(f"    reg {reg}: vg_loc={float(np.asarray(vgl)[nid]):.4g}  msg_mode_g={float(np.asarray(mg)[nid]) if mg is not None else float('nan'):.3f}"
          f"  msg_prec_g={float(np.asarray(pg)[nid]):.4g}  var*={val('var_star'):.4g}")
    print(f"        LEFT={nb(lft)}  RIGHT={nb(rgt)}")
    row = [reg, nid, fg_or[reg], val("fg_strand"), val("fg_loc"), val("f_g"),
           M, E, (M / E if E else float("nan")), val("rho_obs"), val("x_obs"),
           val("x_star"), val("sigma_node"), val("var_star"), rho_floor, in_floor]
    def fmt(x):
        if isinstance(x, bool): return f"{str(x):>8}"
        if isinstance(x, (int, np.integer)): return f"{int(x):>8}"
        return f"{x:>8.3f}"
    print("  ".join(fmt(x) for x in row))
