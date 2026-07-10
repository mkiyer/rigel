"""Oracle EM-lever test for the nascent siphon — quantifies which term drives it.

On the nrna_NONE flagship sibling (true nascent abundance = 0, so ANY nascent count is pure siphon), it
re-runs scan+calibrate+EM under each lever and reports the total nascent count (the siphon), total mature
count, and total gDNA count. All levers are Python-level (NO recompile):

  baseline   : production.
  strand     : RNA strand model p_sense=p_antisense=0.5 ⇒ RNA strand term = gDNA's LOG_HALF (removes the
               +0.68-nat RNA-over-gDNA strand tilt). Tests the strand-asymmetry hypothesis.
  fl         : gDNA FL pmf := RNA FL pmf ⇒ FL gives ZERO gDNA-vs-RNA discrimination. Tests how much the FL
               term is holding the leak back (siphon should GROW).
  strand+fl  : both neutralized ⇒ only prior + eff-length + coverage remain.

A lever that COLLAPSES the nascent siphon is the driver. Deterministic (OMP_NUM_THREADS=1). One re-scan per
lever (the buffer is consumed by a quant pass).

Usage: OMP_NUM_THREADS=1 python scripts/debug/oracle_lever_test.py [condition]
"""
import os
import sys
from dataclasses import replace as dc

import numpy as np

from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType

import rigel.calibration.capture_eff_length as cel_mod
from _metrics import oracle_node_masses, rna_component_breakdown

S = os.environ.get("AUDIT_SUITE", "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb")
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
# efflen_off: disable capture contraction everywhere (nascent eff → full fl, 3x longer → far less
#             competitive). no_nascent: inflate nascent eff x100 (deweight nascent out of the EM) — where
#             does its mass go, gDNA (correct) or mature?
LEVERS = (os.environ["LEVERS"].split(",") if os.environ.get("LEVERS")
          else ["baseline", "oracle_calib", "strand", "fl", "strand+fl", "efflen_off", "no_nascent"])
_REAL_TCEL = cel_mod.transcript_capture_eff_lengths
_REAL_GREF = cel_mod._global_reference_density


# oracle_node_masses + rna_component_breakdown live in _metrics.py (the single source of truth for the
# siphon metric + oracle calibration — see docs/calibration/siphon_measurement.md).


def run(index, bam, cfg, lever):
    scan = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, scan)
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
    if lever == "oracle_calib":  # replace fitted per-node masses with the TRUE masses (perfect calibration)
        import dataclasses
        cal = dataclasses.replace(cal, **oracle_node_masses(bam, ra, index))
    if "strand" in lever:  # neutralize the RNA strand tilt (RNA strand term := gDNA's LOG_HALF)
        sm.exonic_spliced._cached_p_sense = 0.5
        sm.exonic_spliced._cached_p_antisense = 0.5
    if "fl" in lever:  # neutralize FL discrimination (gDNA FL := RNA FL)
        fl = dc(fl, gdna_pmf=fl.rna_pmf.copy())
    # eff-length levers: patch the (lazily-imported) module functions so quant + priors pick them up
    cel_mod._global_reference_density = _REAL_GREF
    cel_mod.transcript_capture_eff_lengths = _REAL_TCEL
    if lever == "efflen_off":  # disable capture contraction everywhere (nascent eff → full fl)
        cel_mod._global_reference_density = lambda *a, **k: None
    if lever == "no_nascent":  # inflate nascent eff x100 ⇒ nascent priced out of the EM competition
        is_n = index.t_df["is_nrna"].to_numpy(bool)
        def _tcel(calibration, region_arrays, idx, fl_eff, _r=_REAL_TCEL, _m=is_n):
            e = _r(calibration, region_arrays, idx, fl_eff).copy()
            e[_m] *= 100.0
            return e
        cel_mod.transcript_capture_eff_lengths = _tcel
    try:
        est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                                   em_config=cfg.em, scoring=cfg.scoring)
    finally:
        cel_mod._global_reference_density = _REAL_GREF
        cel_mod.transcript_capture_eff_lengths = _REAL_TCEL
    # SIPHON = EM mass on the SYNTHETIC shadows (see _metrics / docs/calibration/siphon_measurement.md).
    siphon, single_exon, mature = rna_component_breakdown(est, index)
    syn = index.t_df["is_synthetic"].to_numpy(bool)
    nn = index.t_df["is_nrna"].to_numpy(bool)
    return siphon, siphon + single_exon, mature, int(syn.sum()), int((nn & ~syn).sum())


def main():
    index = TranscriptIndex.load(f"{S}/rigel_index")
    bam = f"{S}/{COND}/sim_oracle.bam"
    cfg = PipelineConfig()
    print(f"=== oracle EM-lever test  cond={COND}  (true nascent=0 ⇒ SIPHON = em_counts on is_synthetic shadows) ===")
    print(f"{'lever':14} {'SIPHON(synthetic)':>18} {'Δ vs base':>12} {'is_nrna(all)':>14} {'mature':>12}   (syn/nonsyn rows)")
    base_s = None
    for lever in LEVERS:
        siphon, nrna_all, mature, n_syn, n_non = run(index, bam, cfg, lever)
        if base_s is None:
            base_s = siphon
        print(f"{lever:14} {siphon:>18,.0f} {siphon - base_s:>+12,.0f} {nrna_all:>14,.0f} {mature:>12,.0f}   ({n_syn}/{n_non})", flush=True)


if __name__ == "__main__":
    main()
