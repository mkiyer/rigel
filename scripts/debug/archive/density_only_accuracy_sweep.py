#!/usr/bin/env python3
"""Density-only calibration accuracy sweep across the capture-on benchmark conditions.

For each condition, scan + calibrate, then compare the GENE0037 exon-region deconv gDNA fraction
(calibration output) to the oracle gDNA fraction (from read names). This characterizes the
density-only scheme (Phases 1+3, no DNA-fraction estimator) across gDNA level × strand specificity.

Run:  python scripts/debug/density_only_accuracy_sweep.py
"""
from __future__ import annotations

import numpy as np
import pysam

from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import coarse_type_array
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _check_region_payload_alignment, scan_and_buffer
from rigel.splice import SpliceType
from rigel.sim.read_name import parse_origin

SUITE = "/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb"
REF = "chr_syn"
LOCUS = (1659774, 1701056)  # GENE0037
CONDS = [
    f"gdna_{g}_ss_{s}_nrna_none_capture_on"
    for s in ("0.99", "0.50")
    for g in ("none", "gdna100", "gdna400", "gdna1000")
]


def main() -> None:
    idx = TranscriptIndex.load(SUITE + "/rigel_index")
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    g = idx.ref_name_to_id[REF]
    rid = np.asarray(ra.ref_id)
    S0, E0 = np.asarray(ra.start), np.asarray(ra.end)
    ct = coarse_type_array(np.asarray(ra.signature))
    win = np.where((rid == g) & (E0 > LOCUS[0]) & (S0 < LOCUS[1]) & (ct == 2))[0]

    print(f"\n{'condition':>42} {'kappa':>6} {'cal_gf':>7} {'oracle_gf':>9} {'abs_err':>7}")
    for cond in CONDS:
        bam = f"{SUITE}/{cond}/sim_oracle.bam"
        try:
            _s, sm, fla, _b, pl = scan_and_buffer(bam, idx, BamScanConfig(sj_strand_tag="auto"))
            sm.finalize()
            _check_region_payload_alignment(ra, pl)
            flm = build_fl_models(
                global_counts=fla.global_model.counts,
                rna_counts=fla.category_models[SpliceType.SPLICED_ANNOT].counts,
                gdna_counts=gdna_fl_mass(pl), max_size=fla.max_size,
            )
            cal = calibrate(pl, ra, sm, flm.gdna_pmf, flm.rna_pmf, CalibrationConfig())
        except Exception as exc:  # noqa: BLE001
            print(f"{cond:>42}  calib error: {exc}")
            continue
        # calibration deconv gDNA fraction (contained) per region
        denom = cal.mass_gdna_contained + cal.mass_rna_contained
        cal_gf = np.where(denom > 0, cal.mass_gdna_contained / np.maximum(denom, 1e-9), 0.0)
        # oracle gDNA fraction over GENE0037 exons (contained template span)
        ogd = np.zeros(ra.n_regions)
        otot = np.zeros(ra.n_regions)
        with pysam.AlignmentFile(bam, "rb") as b:
            for r in b:
                if r.is_read2 or r.is_secondary or r.is_supplementary or r.reference_name != REF:
                    continue
                s = min(r.reference_start, r.next_reference_start)
                e = s + abs(r.template_length)
                j = np.searchsorted(E0, s, side="right")
                if not (0 <= j < ra.n_regions and rid[j] == g and S0[j] <= s and e <= E0[j]):
                    continue
                if parse_origin(r.query_name).kind == "gdna":
                    ogd[j] += 1
                otot[j] += 1
        m = win[otot[win] > 0]
        cal_m = float(np.mean(cal_gf[m]))
        orc_m = float(np.mean(ogd[m] / np.maximum(otot[m], 1)))
        print(f"{cond:>42} {cal.rna_sense_frac:>6.3f} {cal_m:>7.3f} {orc_m:>9.3f} {abs(cal_m-orc_m):>7.3f}")


if __name__ == "__main__":
    main()
