"""Warm-start A/B: does removing the ambiguous-fragment warm-start seed kill the nascent siphon?
Runs the flagship (gdna300 ss99 capON nrna_none) with ORACLE calibration (isolate the EM) and the current
RIGEL_EM_WARMSTART mode (read once by the C++). true nascent=0 ⇒ nrna_em_count = pure siphon.
Run this in SEPARATE processes with RIGEL_EM_WARMSTART in {(unset)=default, unambig, uniform}.
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
from rigel.sim.read_name import parse_origin

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
bam = f"{S}/{COND}/sim_oracle.bam"
index = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
wd = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_oracle_split"
mode = os.environ.get("RIGEL_EM_WARMSTART", "default")

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
scan = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, scan)
sm.finalize()
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
cal = dataclasses.replace(cal, **override)
est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                           em_config=cfg.em, scoring=cfg.scoring)
g = float(est.gdna_em_count) + float(stats.n_intergenic)
n = float(est.nrna_em_count)
m = float(est.get_counts_df(index)["count"].sum())
print(f"WARMSTART={mode:8}  true_gdna={true_g:,}  gdna_assigned={g:,.0f}  gdna_leak={true_g-g:+,.0f}  "
      f"siphon(nrna)={n:,.1f}  mature={m:,.0f}", flush=True)
