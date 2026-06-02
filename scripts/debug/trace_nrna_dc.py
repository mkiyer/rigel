#!/usr/bin/env python
"""Full-chain trace for the nrna_dc g20 failures (count/EM investigation).

Builds nrna_double_count (t1 8-exon +, abundance 100; t_ctrl 8-exon −, abundance 0,
physically separated) with (gdna, nrna, ss), runs the full pipeline, and dumps the
whole chain — calibration hyperparameters, per-region gDNA/RNA deconvolution,
per-locus prior (alpha_gdna_add / alpha_rna_add) + EM split, and final
per-transcript / nRNA / gDNA counts vs expected — to localize where reads are
mis-routed (calibration deconvolution vs per-locus prior vs EM).

Usage:  python scripts/debug/trace_nrna_dc.py --gdna 20 --nrna 0 --ss 0.9
"""

from __future__ import annotations

import argparse
import tempfile

import numpy as np

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark

_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}
_TS = {0: "NONE", 1: "POS", -1: "NEG", 2: "AMBIG"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gdna", type=float, default=20.0)
    ap.add_argument("--nrna", type=float, default=0.0)
    ap.add_argument("--ss", type=float, default=0.9)
    ap.add_argument("--nfrag", type=int, default=2000)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="nrna_dc_")
    sc = Scenario("nrna_dc", genome_length=20000, seed=42, work_dir=tmp)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(2000, 2500), (3000, 3500),
        (4000, 4500), (5000, 5500), (6000, 6500), (7000, 7500), (8000, 8500),
        (9500, 10000)], "abundance": 100}])
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(12000, 12500),
        (13000, 13500), (14000, 14500), (15000, 15500), (16000, 16500),
        (17000, 17500), (18000, 18500), (19000, 19500)], "abundance": 0}])
    cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                        read_length=100, strand_specificity=args.ss, seed=42)
    gdna = None if args.gdna == 0 else GDNAConfig(
        abundance=args.gdna, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
    result = sc.build_oracle(n_fragments=args.nfrag, sim_config=cfg,
                             gdna_config=gdna, nrna_abundance=args.nrna)
    index = result.index
    pr = run_pipeline(result.bam_path, index,
                      config=PipelineConfig(em=EMConfig(seed=42),
                                            scan=BamScanConfig(sj_strand_tag="auto")))
    cal = pr.calibration
    bench = run_benchmark(result, pr,
                          scenario_name=f"g{int(args.gdna)}_n{int(args.nrna)}_s{int(args.ss * 100)}")

    print(f"\n===== {bench.scenario_name}  (gdna={args.gdna} nrna={args.nrna} ss={args.ss}) =====")
    print(f"calibration: iters={cal.n_iterations} converged={cal.converged}  rho_0={cal.rho_0:.4g}  "
          f"exp_disp={cal.exposure_dispersion:.4g}  kappa_rna={cal.kappa_rna:.4g}  "
          f"rho_r_bb={cal.rho_r_bb:.4g}  rho_d_bb={cal.rho_d_bb:.4g}")

    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    rtype = coarse_type_array(ra.signature)
    m_g = cal.mass_g_contained + cal.mass_g_left + cal.mass_g_right
    m_d = cal.mass_d_contained + cal.mass_d_left + cal.mass_d_right
    tot = m_g + m_d
    print("\n--- per-region deconvolution (regions with mass) ---")
    print("  idx ts    rtype      start   end    m_g     m_d  gfrac  omega")
    for i in range(ra.n_regions):
        if tot[i] < 0.5:
            continue
        gfrac = m_g[i] / tot[i] if tot[i] > 1e-9 else float("nan")
        gene = "t1" if ra.start[i] < 11000 else ("t_ctrl" if ra.start[i] < 19500 else "")
        print(f"  {i:3d} {_TS.get(int(ra.ts_class[i]), '?'):>5s} {_RTYPE.get(int(rtype[i]), '?'):>9s}"
              f" {int(ra.start[i]):6d}{int(ra.end[i]):6d} {m_g[i]:6.1f}{m_d[i]:6.1f} {gfrac:5.2f}"
              f" {cal.omega[i]:.3f}  {gene}")

    print("\n--- per-locus prior + EM (get_loci_df) ---")
    loci = pr.estimator.get_loci_df(index)
    cols = [c for c in ["locus_id", "n_transcripts", "locus_span_bp", "n_em_fragments",
            "alpha_gdna_add", "alpha_rna_add", "rna_total", "gdna", "gdna_eff_len_em"]
            if c in loci.columns]
    print(loci[cols].to_string(index=False))

    print("\n--- counts (expected vs observed) ---")
    for t in bench.transcripts:
        print(f"  {t.t_id:8s} expected={t.expected:6.0f}  observed={t.observed:7.1f}")
    print(f"  gDNA:      expected={bench.n_gdna_expected:6.0f}  observed={bench.n_gdna_pipeline:7.1f}")
    print(f"  nRNA:      expected={bench.n_nrna_expected:6.0f}  observed={bench.n_nrna_pipeline:7.1f}")
    print(f"  total RNA: expected={bench.total_expected + bench.n_nrna_expected:6.0f}  "
          f"observed={bench.total_rna_observed:7.1f}")


if __name__ == "__main__":
    main()
