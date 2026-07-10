"""Backward LAYER trace of the residual on gdna300 ss99 capON nrna_none (quick). Separates the sources:
  A baseline                     — production masses + production capture eff-length + EM
  B oracle_masses                — cal masses -> ORACLE (isolates calibration-MASS error: A-B)
  C oracle_masses + efflen_OFF   — also disable capture eff-length contraction (isolates eff-length: B-C)
  D oracle_masses + nascent_out  — also price nascent shadows out (x100 eff-len) (isolates nascent-shadow sink)
True nascent = 0 ⇒ siphon = nrna_em_count. gDNA leak = true_g − assigned_g.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from dataclasses import replace as dc
import dataclasses
from pathlib import Path
import numpy as np
import pysam
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from oracle import OracleTruth
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType
import rigel.calibration.capture_eff_length as cel
from rigel.sim.read_name import parse_origin

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
bam = f"{S}/{COND}/sim_oracle.bam"
index = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
wd = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_oracle_split"
_REAL_GREF, _REAL_TCEL = cel._global_reference_density, cel.transcript_capture_eff_lengths
is_nrna = index.t_df["is_nrna"].to_numpy(bool)

# true observed gDNA fragment count
true_g = 0
with pysam.AlignmentFile(bam, "rb") as f:
    seen = set()
    for rd in f:
        if rd.query_name in seen:
            continue
        seen.add(rd.query_name)
        if parse_origin(rd.query_name).kind == "gdna":
            true_g += 1

orc = OracleTruth.from_bam(bam, index, cfg, wd, COND)
override = orc.override_masses(ra)


def run(arm):
    cel._global_reference_density, cel.transcript_capture_eff_lengths = _REAL_GREF, _REAL_TCEL
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, sc)
    sm.finalize()
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
    if arm != "A":
        cal = dataclasses.replace(cal, **override)
    if arm == "C":
        cel._global_reference_density = lambda *a, **k: None
    if arm == "D":
        def _tcel(calibration, region_arrays, idx, fl_eff, _r=_REAL_TCEL, _m=is_nrna):
            e = _r(calibration, region_arrays, idx, fl_eff).copy()
            e[_m] *= 100.0
            return e
        cel.transcript_capture_eff_lengths = _tcel
    try:
        est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                                   em_config=cfg.em, scoring=cfg.scoring)
    finally:
        cel._global_reference_density, cel.transcript_capture_eff_lengths = _REAL_GREF, _REAL_TCEL
    g = float(est.gdna_em_count) + float(stats.n_intergenic)
    return g, float(est.nrna_em_count), float(est.get_counts_df(index)["count"].sum())


print(f"=== {COND}  true_gdna={true_g:,} ===")
print(f"{'arm':32} {'gdna_assigned':>13} {'gdna_leak':>10} {'siphon(nrna)':>13} {'mature':>10}")
for arm, label in [("A", "A baseline (production)"),
                   ("B", "B +oracle masses"),
                   ("C", "C +oracle masses +efflen OFF"),
                   ("D", "D +oracle masses +nascent priced out")]:
    g, n, m = run(arm)
    print(f"{label:32} {g:>13,.0f} {true_g-g:>10,.0f} {n:>13,.1f} {m:>10,.0f}", flush=True)
