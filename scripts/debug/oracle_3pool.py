"""Decisive test: does PERFECT calibration (oracle per-node masses) fix the capON gDNA->RNA leak?

Runs scan+calibrate+EM under baseline vs oracle_calib (cal masses replaced by ground truth from read
origins) and prints the 3-pool EM totals (gDNA / nascent / mature) vs the manifest OBSERVED truth. If
oracle_calib collapses the -136k gDNA leak -> calibration IS the source. If the leak persists -> the
residual is DOWNSTREAM (prior->EM translation or EM assignment), given an accurate calibration.
"""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
from dataclasses import replace as dc
import dataclasses
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))  # for _metrics when run from scratchpad
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


def run(oracle):
    scan = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, scan)
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
    if oracle:
        cal = dataclasses.replace(cal, **oracle_node_masses(bam, ra, index))
    est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                               em_config=cfg.em, scoring=cfg.scoring)
    gdna = float(est.gdna_em_count) + float(stats.n_intergenic)
    nrna = float(est.nrna_em_count)
    mrna = float(est.get_counts_df(index)["count"].sum())
    return {"gdna": gdna, "nrna": nrna, "mrna": mrna}


print(f"=== {COND} ===")
print(f"{'':12} {'gdna':>28} {'nascent(siphon)':>28} {'mature':>28}")
print(f"{'TRUTH':12} {T['gdna']:>28,.0f} {T['nrna']:>28,.0f} {T['mrna']:>28,.0f}")
for name, orc in [("baseline", False), ("oracle_calib", True)]:
    a = run(orc)
    print(f"{name:12} " + " ".join(
        f"{a[p]:>12,.0f} (surplus {a[p]-T[p]:>+10,.0f})" for p in ("gdna", "nrna", "mrna")))
