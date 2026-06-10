"""Diagnostic: how LOCAL is the clean gDNA density under capture vs no-capture?

Settles the architecture question (CALIBRATION_TODO Issue #1 redesign): is a pooled/global gDNA
density a viable shrinkage target, or does hybrid capture make density so spatially variable that
only LOCAL imputation can work? Calibrates a capture_on and a capture_off condition and reports the
spatial distribution of the per-region clean gDNA density over COUNT-OBSERVABLE regions (gDNA-pure
by signature), where the density is read directly with no imputation.

We also report, per observable region, the strand-clean's effect: raw count/eff_len vs cleaned
density, to see how much the strand clean is moving things on nrna_none data (where observable mass
is PURE gDNA, so the clean should be a near no-op in expectation — any large moves are noise).
"""

import os
import json
import numpy as np
from dataclasses import replace as _dc_replace

from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import node_gdna_density, count_observable_masks
from rigel.calibration.effective_length import region_eff_length, boundary_side_eff_length, boundary_eff_length
from rigel.calibration.strand_balance import fit_strand_balance
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
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    region_eff_len = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    fl_mean = boundary_eff_length(fl.gdna_pmf)
    balance = fit_strand_balance(sm)
    kappa = float(balance.rna_sense_frac)
    nd = node_gdna_density(substrate, ra, region_eff_len, fl_mean, rna_sense_frac=kappa,
                           strand_clean_prior_weight=cfg.calibration.strand_clean_prior_weight)

    obs = np.asarray(nd.region_count_observable)
    dens = np.asarray(nd.density)
    eff = np.asarray(region_eff_len)
    contained_mass = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    # raw (uncleaned) density at observable regions = total contained mass / eff_len
    with np.errstate(divide="ignore", invalid="ignore"):
        raw_dens = np.where(eff > 0, contained_mass / np.maximum(eff, 1e-9), np.nan)

    # observable regions with real coverage
    seed = obs & (eff > 50) & (contained_mass > 0)
    d = dens[seed]
    rd = raw_dens[seed]
    d_pos = d[d > 0]
    print(f"\n=== {cond}  (kappa={kappa:.4f}) ===")
    print(f"  n_regions={len(obs)}  n_observable={int(obs.sum())}  n_seed(cov>0,eff>50)={int(seed.sum())}")
    if len(d_pos) == 0:
        print("  no positive-density observable seeds")
        return
    print(f"  CLEAN gDNA density over observable seeds (frags/bp):")
    print(f"    mean={d.mean():.4g}  median={np.median(d):.4g}  std={d.std():.4g}  CV={d.std()/max(d.mean(),1e-12):.3f}")
    qs = np.percentile(d_pos, [1, 5, 25, 50, 75, 95, 99])
    print(f"    pctiles[1,5,25,50,75,95,99]={np.array2string(qs, precision=4)}")
    print(f"    max/median ratio={d_pos.max()/np.median(d_pos):.1f}  p99/p1={qs[-1]/max(qs[0],1e-12):.1f}")
    # clean vs raw: how much is the strand clean moving pure-gDNA observable density?
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = d / np.where(rd > 0, rd, np.nan)
    ratio = ratio[np.isfinite(ratio)]
    print(f"  strand-clean effect on observable density (clean/raw; 1.0=no-op on pure gDNA):")
    print(f"    mean={ratio.mean():.4f}  median={np.median(ratio):.4f}  std={ratio.std():.4f}  "
          f"min={ratio.min():.3f}  max={ratio.max():.3f}")
    return d_pos


if __name__ == "__main__":
    on = analyze("gdna_gdna400_ss_0.99_nrna_none_capture_on")
    off = analyze("gdna_gdna400_ss_0.99_nrna_none_capture_off")
    # also the unstranded capture case — where the clean is at kappa~0.5 (its failure mode)
    on50 = analyze("gdna_gdna400_ss_0.50_nrna_none_capture_on")
