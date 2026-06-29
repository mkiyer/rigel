"""Node-level diagnosis of the log-density (logodds) regression vs the lattice.

For a condition, runs BOTH solvers with the _capture hook + the oracle (split BAMs), and decomposes
per REGION node: oracle f_g, the strand-only f_g, the local f_g (strand⊗global), the final f_g (after
messages) — for lattice AND logodds. Aggregates by node class (single-strand exon / AMBIG exon) so we
see WHERE the under-call enters (which class, which stage) and the logodds global TARGET + precision.

    OMP_NUM_THREADS=1 python scripts/debug/logodds_node_diagnosis.py [cond]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
DEFAULT = "gdna_gdna300_ss_0.50_nrna_none_capture_on"


def run_capture(blob, ra, solver):
    cm = importlib.import_module("rigel.calibration.calibrate")
    orig = cm.node_sweep
    cap = {}
    cm.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    try:
        calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                  gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"],
                  config=CalibrationConfig(calibration_solver=solver))
    finally:
        cm.node_sweep = orig
    return cap


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    kind = np.asarray(ch.kind)
    refidx = np.asarray(ch.ref_idx, np.int64)
    is_reg = kind == REGION
    R = len(index.region_df)
    reg_node = {int(refidx[i]): int(i) for i in np.where(is_reg)[0]}

    sub_g = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    sub_r = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    g_or = np.asarray(sub_g.contained.mass_unspliced, float)
    r_or = np.asarray(sub_r.contained.mass_unspliced, float)
    tot = g_or + r_or
    fg_true = np.where(tot > _EPS, g_or / np.maximum(tot, _EPS), np.nan)

    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    lat = run_capture(blob, ra, "lattice")
    lo = run_capture(blob, ra, "logodds")

    # region-id arrays for capture (indexed by chain node)
    def at(cap, key):
        return np.asarray(cap[key], float)

    print(f"\n================= {cond} =================")
    print(f"κ={'?'}  rho_global lat={float(lat['rho_global']):.4f} lo={float(lo['rho_global']):.4f}  "
          f"enrich_w lat={float(lat['enrich_w']):.2f} lo={float(lo['enrich_w']):.2f}")

    # Per-class aggregate over EXON region nodes with real mass.
    classes = [("SS-exon", lambda r: ((scl[r] == 1) | (scl[r] == 2)) & (tcl[r] == 2)),
               ("AMBIG-exon", lambda r: (scl[r] == 3) & (tcl[r] == 2))]
    cols = "class n mass true LAT:str LAT:loc LAT:fin  LO:str LO:loc LO:fin".split()
    print("\n  MASS-WEIGHTED mean f_g by class+stage (strand→local→final):")
    print("  " + "  ".join(f"{c:>8}" for c in cols))
    for name, sel in classes:
        rids = np.array([r for r in range(R) if sel(r) and tot[r] > 1.0])
        if rids.size == 0:
            continue
        nodes = np.array([reg_node[int(r)] for r in rids])
        w = tot[rids]
        sw = max(w.sum(), _EPS)
        def mw(cap, key):
            return float(np.sum(at(cap, key)[nodes] * w) / sw)
        row = [name, f"{rids.size}", f"{w.sum():,.0f}",
               f"{np.sum(fg_true[rids]*w)/sw:.3f}",
               f"{mw(lat,'fg_strand'):.3f}", f"{mw(lat,'fg_loc'):.3f}", f"{mw(lat,'f_g'):.3f}",
               f"{mw(lo,'fg_strand'):.3f}", f"{mw(lo,'fg_loc'):.3f}", f"{mw(lo,'f_g'):.3f}"]
        print("  " + "  ".join(f"{v:>8}" for v in row))

    # logodds GLOBAL target (implied f_g) + N on exon nodes — is the target right? how strong?
    ehat = lo["ehat"]; z = at(lo, "z_enrich"); eff = at(lo, "eff_global"); mass = at(lo, "mass_global")
    rg = float(lo["rho_global"]); w_en = float(lo["enrich_w"]); app = np.asarray(lo["enrich_apply"], bool)
    print("\n  logodds GLOBAL implied-f_g target vs oracle (exon nodes, enrich-applied):")
    for name, sel in classes:
        rids = np.array([r for r in range(R) if sel(r) and tot[r] > 1.0])
        nodes = np.array([reg_node[int(r)] for r in rids]) if rids.size else np.array([], int)
        am = app[nodes] if nodes.size else np.array([], bool)
        if not am.any():
            print(f"    {name}: (none enrich-applied)")
            continue
        zn = z[nodes][am]
        rho_e = np.maximum(ehat.predict(zn), _EPS)
        rho_hat = np.exp(w_en*np.log(rho_e) + (1-w_en)*np.log(max(rg,_EPS)))
        target_frac = np.clip(rho_hat*eff[nodes][am]/np.maximum(mass[nodes][am],_EPS), _EPS, 1.0)
        ot = fg_true[rids][am]
        print(f"    {name}: n={am.sum()}  mean target_f_g={np.nanmean(target_frac):.3f}  "
              f"mean oracle_f_g={np.nanmean(ot):.3f}  (target≈oracle? {'YES' if abs(np.nanmean(target_frac)-np.nanmean(ot))<0.15 else 'NO — global target is wrong'})")

    # ── MESSAGE TRACE on the worst-crashed SS-exon nodes (local→final divergence) ──
    vg_loc = at(lo, "vg_loc"); mode_g = at(lo, "mode_g"); prec_g = at(lo, "prec_g")
    fg_loc = at(lo, "fg_loc"); f_fin = at(lo, "f_g")
    ss_ex = np.array([r for r in range(R) if ((scl[r] == 1) | (scl[r] == 2)) & (tcl[r] == 2) and tot[r] > 100])
    nodes = np.array([reg_node[int(r)] for r in ss_ex])
    crash = (fg_loc[nodes] - f_fin[nodes]) * tot[ss_ex]  # downward message mass
    worst = np.argsort(-crash)[:6]
    print("\n  MESSAGE TRACE — worst SS-exon local→final crashes (logodds):")
    print("   rid  mass  oracle  local  pg_loc  msg_mode(log f) msg_impl_f  prec_g   final")
    for j in worst:
        r = int(ss_ex[j]); n = nodes[j]
        pgl = 1.0 / max(vg_loc[n], _EPS)
        impl = float(np.exp(min(mode_g[n], 0.0)))  # implied message f_g (clip log≤0 for display)
        print(f"   {r:>4} {tot[r]:>6,.0f}  {fg_true[r]:.3f}  {fg_loc[n]:.3f}  {pgl:>6.2f}  "
              f"{mode_g[n]:>9.2f}      {impl:.4f}    {prec_g[n]:>6.2f}   {f_fin[n]:.3f}")


if __name__ == "__main__":
    main()
