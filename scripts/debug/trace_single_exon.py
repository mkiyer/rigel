#!/usr/bin/env python
"""Trace the gDNA over-call on the single-exon scenario (calibration → prior → EM).

Builds the `single_exon` baseline (ss=1.0, no gDNA, with a spliced helper gene),
runs the full pipeline once, and dumps the whole chain for the single-exon gene
t1 (which loses ~298/339 reads to gDNA): the calibration hyperparameters, the
per-region deconvolved gDNA fraction over the gene regions, the per-locus prior
(alpha_gdna_add / alpha_rna_add), and the final per-transcript counts.

Usage:  python scripts/debug/trace_single_exon.py [--gdna 0] [--ss 1.0]
"""

from __future__ import annotations

import argparse
import tempfile

import numpy as np

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}
_TS = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gdna", type=float, default=0.0)
    ap.add_argument("--ss", type=float, default=1.0)
    ap.add_argument("--nfrag", type=int, default=500)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="single_exon_")
    sc = Scenario("single_exon", genome_length=8000, seed=42, work_dir=tmp)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1500)], "abundance": 100}])
    sc.add_gene(
        "g_helper",
        "+",
        [{"t_id": "t_helper", "exons": [(2500, 3000), (3500, 4000)], "abundance": 50}],
    )
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0}])
    cfg = ReadSimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=args.ss, seed=42,
    )
    gdna = None if args.gdna == 0 else GDNAConfig(
        abundance=args.gdna, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000
    )
    result = sc.build_oracle(n_fragments=args.nfrag, sim_config=cfg, gdna_config=gdna)
    index = result.index
    pr = run_pipeline(
        result.bam_path, index,
        config=PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="auto")),
    )
    cal = pr.calibration

    print("\n=== calibration hyperparameters ===")
    print(f"  iters={cal.n_iterations} converged={cal.converged}")
    print(f"  exposure_dispersion = {cal.exposure_dispersion:.5g}   rho_0 = {cal.rho_0:.5g}")
    print(f"  kappa_rna = {cal.kappa_rna:.5g}   rho_r_bb = {cal.rho_r_bb:.5g}   "
          f"rho_d_bb = {cal.rho_d_bb:.4g}")

    # Per-region deconvolved gDNA fraction over the annotated (gene) regions.
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pr.calibration_payload, ra)
    rtype = coarse_type_array(ra.signature)
    m_g = cal.mass_g_contained + cal.mass_g_left + cal.mass_g_right
    m_d = cal.mass_d_contained + cal.mass_d_left + cal.mass_d_right
    tot = m_g + m_d
    nu = sub.contained.n_unspliced
    ns = sub.contained.n_spliced
    print("\n=== gene-region deconvolution (regions with any contained flux) ===")
    print("  idx  ts     rtype       start    end    n_u  n_s    m_g     m_d   gDNA_frac  omega")
    for i in range(ra.n_regions):
        if nu[i] == 0 and ns[i] == 0 and tot[i] < 1e-9:
            continue
        gfrac = m_g[i] / tot[i] if tot[i] > 1e-12 else float("nan")
        print(f"  {i:3d}  {_TS.get(int(ra.ts_class[i]), '?'):>5s}  {_RTYPE.get(int(rtype[i]), '?'):>10s}"
              f"  {int(ra.start[i]):6d} {int(ra.end[i]):6d} {int(nu[i]):4d} {int(ns[i]):4d}"
              f"  {m_g[i]:6.1f} {m_d[i]:6.1f}   {gfrac:6.3f}   {cal.omega[i]:.3f}")

    print("\n=== per-locus prior + EM result (get_loci_df) ===")
    loci = pr.estimator.get_loci_df(index)
    cols = [c for c in
            ["locus_id", "n_transcripts", "n_genes", "locus_span_bp", "n_em_fragments",
             "alpha_gdna_add", "alpha_rna_add", "rna_total", "gdna", "gdna_eff_len_em"]
            if c in loci.columns]
    print(loci[cols].to_string(index=False))

    print("\n=== final transcript counts (expected vs observed) ===")
    counts = pr.estimator.get_counts_df(index)
    for _, row in counts.iterrows():
        if row["transcript_id"] in ("t1", "t_helper", "t_ctrl"):
            print(f"  {row['transcript_id']:10s}  count={row['count']:.1f}")
    print(f"  gdna_em_count = {pr.estimator.gdna_em_count:.1f}   "
          f"n_intergenic = {pr.stats.n_intergenic}")


if __name__ == "__main__":
    main()
