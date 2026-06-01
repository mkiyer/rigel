#!/usr/bin/env python
"""Dump the calibration-v6 state for a BAM + index (PR 5 diagnostic).

Runs the BAM scan and the calibrator, then prints the converged library
hyperparameters, the per-iteration mass-change trajectory, and a per-region-type
(intergenic / intron / exon) summary of the deconvoluted gDNA fraction and
exposure. Use ``-v`` to see the per-iteration M-step trajectory (rho_0, phi,
rho_d_bb, eps_s, delta) logged by the outer loop.

Usage:
    python scripts/debug/dump_calibration_state.py --index INDEX_DIR --bam SAMPLE.bam [-v]

Built for debugging convergence / parameter behaviour (e.g. the armis2
regression hunt; docs/acc_caljointmodel/00_implementation_plan.md §6).
"""

from __future__ import annotations

import argparse
import logging

import numpy as np

from rigel.calibration import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.splice import SpliceType


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", required=True, help="rigel index directory")
    ap.add_argument("--bam", required=True, help="name-sorted BAM with NH tags")
    ap.add_argument("-v", "--verbose", action="store_true", help="log per-iteration M-step state")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    index = TranscriptIndex.load(args.index)
    scan = BamScanConfig(sj_strand_tag="auto")
    _stats, strand_models, fl_models_scan, buffer, payload = scan_and_buffer(args.bam, index, scan)
    try:
        region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
        fl = build_fl_models(
            global_counts=fl_models_scan.global_model.counts,
            rna_counts=fl_models_scan.category_models[SpliceType.SPLICED_ANNOT].counts,
            gdna_counts=gdna_fl_mass(payload),
            max_size=fl_models_scan.max_size,
        )
        result = calibrate(
            payload=payload,
            region_arrays=region_arrays,
            strand_model=strand_models,
            gdna_fl_pmf=fl.gdna_pmf,
            config=CalibrationConfig(),
        )
    finally:
        buffer.cleanup()

    print("\n=== calibration-v6 state ===")
    print(f"regions          : {result.n_regions}")
    print(f"iterations       : {result.n_iterations}  converged={result.converged}")
    print(f"mass_change      : {np.array2string(result.mass_change_history, precision=4)}")
    print("hyperparameters  :")
    print(f"  rho_0   = {result.rho_0:.6g}     exposure_dispersion = {result.exposure_dispersion:.6g}")
    print(f"  rho_d_bb= {result.rho_d_bb:.4g}  eps_s = {result.eps_s:.4g}")
    print(f"  kappa_rna={result.kappa_rna:.4g}  rho_r_bb={result.rho_r_bb:.4g}  (fixed, PR 3)")
    print(
        f"FL means         : gdna={_mean(fl.gdna_pmf):.1f}  rna={_mean(fl.rna_pmf):.1f}  "
        f"global={_mean(fl.global_pmf):.1f}  (n_gdna_pool={fl.n_gdna:.0f})"
    )

    # Per-region-type breakdown of the deconvoluted gDNA fraction + exposure.
    rtype = coarse_type_array(region_arrays.signature)
    m_g = result.mass_g_contained + result.mass_g_left + result.mass_g_right
    m_d = result.mass_d_contained + result.mass_d_left + result.mass_d_right
    total = m_g + m_d
    print("per-region-type  : (gDNA fraction of deconvoluted mass; mean ω)")
    for code, name in ((0, "intergenic"), (1, "intron"), (2, "exon")):
        mask = rtype == code
        if not np.any(mask):
            continue
        tot = float(total[mask].sum())
        gfrac = float(m_g[mask].sum()) / tot if tot > 0 else float("nan")
        print(
            f"  {name:11s}: n={int(mask.sum()):6d}  gDNA_frac={gfrac:.3f}  "
            f"omega_mean={float(result.omega[mask].mean()):.3f}"
        )


def _mean(pmf: np.ndarray) -> float:
    return float(np.dot(np.arange(pmf.size, dtype=np.float64), pmf))


if __name__ == "__main__":
    main()
