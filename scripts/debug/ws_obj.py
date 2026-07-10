"""Compare the EM marginal data log-likelihood at the two fixed points (flat=73k vs unambig=23k).
If flat's total loglik > unambig's, the siphon is the MLE (model fix). If unambig's is higher, flat is a
local-optimum trap (optimizer/init fix). Runs the oracle-calibrated flagship with emit_locus_stats and sums
final_data_loglik over loci, for the current RIGEL_EM_WARMSTART mode. Run in separate processes."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from dataclasses import replace as dc
import dataclasses
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from oracle import OracleTruth
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
bam = f"{S}/{COND}/sim_oracle.bam"
index = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
wd = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_oracle_split"
mode = os.environ.get("RIGEL_EM_WARMSTART", "default")

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
                           em_config=cfg.em, scoring=cfg.scoring, emit_locus_stats=True)
ll = sum(float(d.get("final_data_loglik", 0.0)) for d in (est.locus_stats or []))
siphon = float(est.nrna_em_count)
print(f"WARMSTART={mode:8}  total_data_loglik={ll:,.1f}  siphon(nrna)={siphon:,.0f}  "
      f"n_loci_with_stats={len(est.locus_stats or [])}", flush=True)
