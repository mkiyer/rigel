"""Diagnostic: in a ZERO-gDNA library, WHERE does calibration manufacture false gDNA?

gdna_none => every gDNA call is a false positive. We calibrate and attribute the total gDNA mass
the calibrator assigns across node classes (observable region / exon region / boundary sides) and
SS regimes, to localize the FP source: is it the observable-anchor strand clean, the exon count
prior (local imputation), or the boundary sides? This tells us whether the robust-clean fix belongs
at the anchors, the imputation, or the sides — and whether any global target is even involved.
"""

import os
import numpy as np
from dataclasses import replace as _dc_replace

from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration import calibrate
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.effective_length import region_eff_length, boundary_side_eff_length, boundary_eff_length
from rigel.calibration.signature import TS_POS, TS_NEG, TS_NONE, TS_AMBIG
from rigel.splice import SpliceType

SUITE = os.environ.get("DIAG_SUITE", "/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb")
index = TranscriptIndex.load(f"{SUITE}/rigel_index")


def analyze(cond):
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    cfg = PipelineConfig()
    scan = _dc_replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, scan)
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size,
    )
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, config=cfg.calibration)
    nd = node_gdna_density(ra and None or None, ra, None, None, rna_sense_frac=0.5) if False else None

    # recompute node observability for attribution
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    region_eff_len = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    fl_mean = boundary_eff_length(fl.gdna_pmf)
    from rigel.calibration.strand_balance import fit_strand_balance
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    nd = node_gdna_density(substrate, ra, region_eff_len, fl_mean, rna_sense_frac=kappa,
                           strand_clean_prior_weight=cfg.calibration.strand_clean_prior_weight)
    obs = np.asarray(nd.region_count_observable)
    ts = np.asarray(ra.strand_class)

    g_contained = np.asarray(cal.mass_gdna_contained)
    g_left = np.asarray(cal.mass_gdna_left)
    g_right = np.asarray(cal.mass_gdna_right)
    r_contained = np.asarray(cal.mass_rna_contained)

    total_fp_gdna = g_contained.sum() + g_left.sum() + g_right.sum()
    print(f"\n=== {cond}  (kappa={kappa:.4f}) ===")
    print(f"  TOTAL false gDNA mass = {total_fp_gdna:,.0f}")
    print(f"   region-contained: {g_contained.sum():,.0f}  "
          f"(observable {g_contained[obs].sum():,.0f} / exon+ambig {g_contained[~obs].sum():,.0f})")
    print(f"   boundary-left:    {g_left.sum():,.0f}")
    print(f"   boundary-right:   {g_right.sum():,.0f}")
    # by strand class of the region (contained)
    for name, mask in [("POS", ts == TS_POS), ("NEG", ts == TS_NEG),
                       ("NONE", ts == TS_NONE), ("AMBIG", ts == TS_AMBIG)]:
        if g_contained[mask].sum() > 1:
            print(f"     contained {name:5s}: gdna={g_contained[mask].sum():,.0f}  "
                  f"(observable {g_contained[mask & obs].sum():,.0f})")


if __name__ == "__main__":
    analyze("gdna_none_ss_0.50_nrna_none_capture_on")
    analyze("gdna_none_ss_0.99_nrna_none_capture_on")
