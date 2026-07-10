"""DECISIVE re-attribution on the CORRECT basis: does perfect calibration (the new, validated
accumulator-consistent override) fix the capON gDNA->RNA leak? Baseline vs oracle-override -> EM -> 3-pool.

If oracle-override => leak ~0: the leak IS calibration (small errors snowball). If it STILL leaks: the leak
is downstream (EM/eff-length), given an accurate calibration. Uses scripts/debug/oracle.py (validated),
NOT the retired oracle_node_masses.
"""
import os, sys, json
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
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
bam = f"{S}/{COND}/sim_oracle.bam"
index = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
wd = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_oracle_split"
man = json.load(open(f"{S}/manifest.json"))
truth = next(c for c in man["conditions"] if c["name"] == COND)
T = {"gdna": truth["n_gdna_observed"], "nrna": truth["n_nrna_observed"], "mrna": truth["n_mrna_observed"]}

orc = OracleTruth.from_bam(bam, index, cfg, wd, COND)
print("oracle validated.")
override = orc.override_masses(ra)


def run(use_oracle):
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, sc)
    sm.finalize()
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
    if use_oracle:
        cal = dataclasses.replace(cal, **override)
    est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                               em_config=cfg.em, scoring=cfg.scoring)
    return {"gdna": float(est.gdna_em_count) + float(stats.n_intergenic),
            "nrna": float(est.nrna_em_count),
            "mrna": float(est.get_counts_df(index)["count"].sum())}


print(f"\n=== {COND} — re-attribution on the CORRECT (accumulator) basis ===")
print(f"{'':16} {'gdna surplus':>16} {'siphon(nrna)':>16} {'mature surplus':>16}")
for name, orcl in [("baseline", False), ("oracle_override", True)]:
    a = run(orcl)
    print(f"{name:16} {a['gdna']-T['gdna']:>+16,.0f} {a['nrna']-T['nrna']:>+16,.0f} {a['mrna']-T['mrna']:>+16,.0f}")
