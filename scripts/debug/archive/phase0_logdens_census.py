"""Phase 0 measurement for the log-density 1-D solver redesign.

For each cached condition: the NODE-CLASS census (single-strand vs AMBIG vs G1, over ALL chain nodes
— the cost model for the 1-D grid + the AMBIG tilt marginalization), the gDNA-density ρ_g DYNAMIC
RANGE (oracle + solved, log10 — justifies the log grid and sets the span c), and the f_g EXTREMES
distribution (how much true mass sits at f_g→0 depleted / f_g→1 enriched, where the uniform linear
lattice is coarsest). Pure measurement; no solver change.

    OMP_NUM_THREADS=1 python scripts/debug/phase0_logdens_census.py [cond1 cond2 ...]
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
from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
CONDS = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",   # flagship (stranded, capture)
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",  # stranded, off-capture
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",   # unstranded, capture (AMBIG-hard)
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",  # unstranded, off-capture
    "gdna_gdna300_ss_0.50_nrna_rnd_capture_on",    # unstranded + nascent (E4 phantom)
    "gdna_none_ss_0.50_nrna_none_capture_off",     # zero gDNA (floor / pseudocount)
]


def _pct(a, ps=(0, 1, 5, 50, 95, 99, 100)):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {p: float("nan") for p in ps}
    return {p: float(np.percentile(a, p)) for p in ps}


def census_one(cond: str):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]

    calmod = importlib.import_module("rigel.calibration.calibrate")
    orig = calmod.node_sweep
    cap = {}
    calmod.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    try:
        cal = calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
                        gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"],
                        config=CalibrationConfig())
    finally:
        calmod.node_sweep = orig

    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    kind = np.asarray(ch.kind)
    is_reg = kind == REGION
    n_nodes = kind.size

    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    solv = np.asarray(cap["solvable"], bool)
    f_g = np.asarray(cap["f_g"], float)
    eff_g = np.asarray(cap["eff_global"], float)
    mass_g = np.asarray(cap["mass_global"], float)

    single = fp ^ fn
    ambig = fp & fn
    g1 = ~fp & ~fn

    def share(mask):
        return int(mask.sum()), 100.0 * mask.sum() / max(n_nodes, 1)

    print(f"\n{'='*78}\n{cond}")
    print(f"  κ={cal.rna_sense_frac:.3f}  ρ_global={float(getattr(cal,'gdna_density_global',float('nan'))):.4f}")
    print(f"  chain nodes m = {n_nodes:,}   (regions {int(is_reg.sum()):,}  boundaries {int((~is_reg).sum()):,})")
    print("  NODE-CLASS census (all nodes / solvable-with-mass):")
    for nm, mask in (("single-strand", single), ("AMBIG", ambig), ("G1 (no-RNA)", g1)):
        c, p = share(mask)
        cs, ps = share(mask & solv)
        print(f"    {nm:<16} {c:>7,} ({p:5.1f}%)   solvable {cs:>7,} ({ps:5.1f}%)")
    n_ambig_solv = int((ambig & solv).sum())

    # ρ_g dynamic range — SOLVED per-node gDNA density (f_g·M/E), over solvable nodes with gDNA.
    rho_solved = f_g * mass_g / np.maximum(eff_g, _EPS)
    rs = rho_solved[solv & (rho_solved > _EPS)]
    print("  ρ_g SOLVED density (solvable, >0):")
    pr = _pct(rs)
    print(f"    n={rs.size:,}  min={pr[0]:.4g}  p1={pr[1]:.4g}  p50={pr[50]:.4g}  p99={pr[99]:.4g}  max={pr[100]:.4g}")
    if rs.size:
        lr = np.log10(rs)
        print(f"    log10 range = {lr.max()-lr.min():.2f} decades  (p99/p1 = {pr[99]/max(pr[1],_EPS):.4g}×)")

    # ORACLE ρ_g + f_g for REGIONS (split BAMs) — the TRUE range the grid must span.
    sub_g = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    sub_r = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    g_or = np.asarray(sub_g.contained.mass_unspliced, float)
    r_or = np.asarray(sub_r.contained.mass_unspliced, float)
    tot = g_or + r_or
    reg_eff_g = region_eff_length(np.asarray(ra.region_size_bp, float), blob["gdna_pmf"])
    has = tot > _EPS
    fg_true = np.where(has, g_or / np.maximum(tot, _EPS), np.nan)
    rho_true = np.where((g_or > _EPS), g_or / np.maximum(reg_eff_g, _EPS), np.nan)
    rt = rho_true[np.isfinite(rho_true) & (rho_true > _EPS)]
    print("  ρ_g ORACLE region density (gDNA mass / eff_len, >0):")
    prt = _pct(rt)
    print(f"    n={rt.size:,}  p1={prt[1]:.4g}  p50={prt[50]:.4g}  p99={prt[99]:.4g}  max={prt[100]:.4g}")
    if rt.size:
        lrt = np.log10(rt)
        print(f"    log10 range = {lrt.max()-lrt.min():.2f} decades")

    # f_g EXTREMES — mass-weighted, where the linear lattice (spacing 1/60≈0.0167) is coarsest.
    fgt = fg_true[has]
    w = tot[has]
    bins = [(0, 0.01), (0.01, 0.0167), (0.0167, 0.05), (0.05, 0.95),
            (0.95, 0.9833), (0.9833, 0.99), (0.99, 1.0001)]
    print("  TRUE f_g distribution (mass-weighted) — lattice spacing = 1/60 ≈ 0.0167:")
    for lo, hi in bins:
        m = (fgt >= lo) & (fgt < hi)
        wm = w[m].sum()
        print(f"    f_g∈[{lo:.4f},{hi:.4f}): {int(m.sum()):>6,} regions  {100*wm/max(w.sum(),_EPS):5.1f}% of mass")
    # how many regions have true f_g within ONE lattice cell of a vertex (under-resolved)?
    near0 = int(((fgt < 0.0167) & (fgt > _EPS)).sum())
    near1 = int((fgt > 1 - 0.0167).sum())
    print(f"    within 1 lattice cell of a vertex: f_g<0.0167 → {near0:,} regions ; f_g>0.9833 → {near1:,}")
    return {"cond": cond, "m": n_nodes, "n_ambig_solv": n_ambig_solv}


def main():
    conds = sys.argv[1:] or CONDS
    summ = []
    for c in conds:
        try:
            summ.append(census_one(c))
        except Exception as e:  # noqa: BLE001
            print(f"\n!! {c}: {type(e).__name__}: {e}")
    print(f"\n{'='*78}\nMEMORY MODEL (K=60 g-grid, K_t=40 tilt-grid):")
    for s in summ:
        mk = s["m"] * 60 * 8 / 1e9
        amk = s["n_ambig_solv"] * 60 * 40 * 8 / 1e9
        print(f"  {s['cond']:<46} (m,K)={mk:.4f}GB  AMBIG(m,K,K_t)={amk:.4f}GB")
    print("  VCaP projection @ m=1.5M: (m,K)=0.72GB ; AMBIG share from the census above.")


if __name__ == "__main__":
    main()
