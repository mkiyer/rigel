"""In-process A/B of the mixture-bridge Fix-1 across all conditions (calibration deconvolution level).

Per condition: scan ONCE, then run the real calibrate() twice on the SAME payload (bridge off vs on,
toggled via RIGEL_MIX_BRIDGE in-process → deterministic, no cross-process scanner noise). Reports the
per-region gDNA deconvolution error vs the sim oracle: Σ|err| and net by class (exon/intron/ig).

    OMP_NUM_THREADS=1 python calib_bridge_ab.py <suite_dir> [eps] [cond_glob]
"""
import os, sys, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import numpy as np, pysam

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.calibration.signature import coarse_type_from_signature
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

suite = Path(sys.argv[1])
eps = sys.argv[2] if len(sys.argv) > 2 else "0.01"
cond_glob = sys.argv[3] if len(sys.argv) > 3 else "gdna_*"
idx = TranscriptIndex.load(suite / "rigel_index")
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
starts = np.asarray(ra.start, np.int64); ends = np.asarray(ra.end, np.int64)
roff = np.asarray(ra.ref_offsets, np.int64)
n = ra.n_regions
sig = idx.region_df["signature"].to_numpy()
cls = np.array(["ig", "intron", "exon"])[[coarse_type_from_signature(int(s)) for s in sig]]
cfg = CalibrationConfig()

def oracle_true_g(cond):
    tg = np.zeros(n)
    with pysam.AlignmentFile(str(cond / "sim_oracle.bam"), "rb") as f:
        dref = f.references[0]; seen = set()
        for r in f:
            q = r.query_name
            if q in seen: continue
            seen.add(q); o = parse_origin(q)
            if o.kind != "gdna" or o.start is None: continue
            rid = idx.ref_name_to_id.get(str(o.ref if o.ref is not None else dref))
            if rid is None: continue
            lo0, hi0 = int(roff[rid]), int(roff[rid + 1]); a, b = int(o.start), int(o.end)
            rr = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
            if rr < hi0 and starts[rr] <= a and b <= ends[rr]:
                tg[rr] += 1.0
    return tg

from dataclasses import replace as _replace

def calib_obs(pl, fl, bridge):
    c = _replace(cfg, gdna_prior_mixture_bridge=float(bridge))  # config-param bridge (was RIGEL_MIX_BRIDGE env)
    return np.asarray(calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, c).mass_gdna_contained, float)

conds = sorted(Path(p) for p in glob.glob(str(suite / cond_glob)) if (Path(p) / "sim_oracle.bam").exists())
print(f"in-process bridge A/B  eps={eps}  ({len(conds)} conditions)\n")
hdr = f"{'condition':<44} {'Σ|err|_off':>11} {'Σ|err|_on':>11} {'Δ%':>7}   {'exon_off':>9} {'exon_on':>9}"
print(hdr); print("-" * len(hdr))
tot_off = tot_on = 0.0
rows = []
for cond in conds:
    tg = np.minimum(oracle_true_g(cond), 1e18)
    st, sm, flm, buf, pl = scan_and_buffer(str(cond / "sim_oracle.bam"), idx, BamScanConfig())
    sub = CalibrationSubstrate.from_payload(pl, ra)
    raw = np.asarray(sub.contained.mass_unspliced, float)
    tgc = np.minimum(tg, raw)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    off = calib_obs(pl, fl, 0.0)
    on = calib_obs(pl, fl, eps)
    e_off = np.abs(off - tgc); e_on = np.abs(on - tgc)
    s_off, s_on = e_off.sum(), e_on.sum()
    exon = cls == "exon"
    nx_off = float((off - tgc)[exon].sum()); nx_on = float((on - tgc)[exon].sum())
    d = 100 * (s_on - s_off) / max(s_off, 1)
    print(f"{cond.name:<44} {s_off:>11,.0f} {s_on:>11,.0f} {d:>+6.1f}%   {nx_off:>+9,.0f} {nx_on:>+9,.0f}")
    tot_off += s_off; tot_on += s_on
    rows.append((cond.name, s_off, s_on, d))
print("-" * len(hdr))
print(f"{'TOTAL':<44} {tot_off:>11,.0f} {tot_on:>11,.0f} {100*(tot_on-tot_off)/max(tot_off,1):>+6.1f}%")
worse = [r for r in rows if r[3] > 2.0]
if worse:
    print("\nREGRESSED conditions (Σ|err| worse by >2%):")
    for nm, so, sn, d in sorted(worse, key=lambda r: -r[3]):
        print(f"  {nm:<44} {so:,.0f} -> {sn:,.0f} ({d:+.1f}%)")
