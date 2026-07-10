"""Levered oracle: which calibration OUTPUT is wrong? Override oracle truth for only a SUBSET of the
per-node mass arrays and see which subset collapses the capON gDNA leak + siphon.

  gdna    : mass_gdna_{contained,left,right}      (the gDNA side of the split)
  rna     : mass_rna_{contained,left,right}       (the unspliced RNA side)
  spliced : mass_rna_spliced                      (the mature MEASUREMENT hint)
  all     : everything (== full oracle_calib)

Whichever lever moves the 3-pool to truth is the culprit calibration output.
"""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
from dataclasses import replace as dc
import dataclasses
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType
from _metrics import oracle_node_masses

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
bam = f"{S}/{COND}/sim_oracle.bam"
index = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
man = json.load(open(f"{S}/manifest.json"))
truth = next(c for c in man["conditions"] if c["name"] == COND)
T = {"gdna": truth["n_gdna_observed"], "nrna": truth["n_nrna_observed"], "mrna": truth["n_mrna_observed"]}

GDNA = ["mass_gdna_contained", "mass_gdna_left", "mass_gdna_right"]
RNA = ["mass_rna_contained", "mass_rna_left", "mass_rna_right"]
SPL = ["mass_rna_spliced"]
LEVERS = {"baseline": [], "oracle_gdna": GDNA, "oracle_rna": RNA, "oracle_spliced": SPL,
          "oracle_gdna+spl": GDNA + SPL, "oracle_all": GDNA + RNA + SPL}


def run(keys):
    scan = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, scan)
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
    if keys:
        om = oracle_node_masses(bam, ra, index)
        cal = dataclasses.replace(cal, **{k: om[k] for k in keys})
    est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                               em_config=cfg.em, scoring=cfg.scoring)
    return {"gdna": float(est.gdna_em_count) + float(stats.n_intergenic),
            "nrna": float(est.nrna_em_count),
            "mrna": float(est.get_counts_df(index)["count"].sum())}


print(f"=== {COND} ===  (surplus = assigned - true;  siphon = nascent surplus)")
print(f"{'TRUTH':16} gdna={T['gdna']:>10,.0f}  nrna={T['nrna']:>9,.0f}  mrna={T['mrna']:>10,.0f}")
for name, keys in LEVERS.items():
    a = run(keys)
    print(f"{name:16} gdna={a['gdna']-T['gdna']:>+10,.0f}  siphon={a['nrna']-T['nrna']:>+9,.0f}  "
          f"mrna={a['mrna']-T['mrna']:>+10,.0f}")
