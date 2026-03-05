#!/usr/bin/env python
"""Diagnostic: trace nrna_frac prior and EM to understand nRNA under-estimation."""
import numpy as np
import tempfile
import pathlib
import logging

logging.basicConfig(level=logging.WARNING)

from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
from hulkrna.pipeline import run_pipeline
from hulkrna.sim import Scenario, SimConfig, run_benchmark

with tempfile.TemporaryDirectory() as td:
    td = pathlib.Path(td)
    sc = Scenario("diag", genome_length=20000, seed=42, work_dir=td / "diag")
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(2000, 4000), (8000, 10000)], "abundance": 100},
    ])
    sc.add_gene("g_ctrl", "-", [
        {"t_id": "t_ctrl", "exons": [(14000, 16000), (18000, 19000)], "abundance": 0},
    ])
    sim = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=1.0, seed=42,
    )
    result = sc.build_oracle(n_fragments=10000, sim_config=sim, nrna_abundance=100)
    cfg = PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="auto"))
    pr = run_pipeline(result.bam_path, result.index, config=cfg)
    bench = run_benchmark(result, pr, scenario_name="diag")
    est = pr.estimator

    print("=== Ground truth ===")
    print(f"mRNA gt: {bench.total_expected}")
    print(f"nRNA gt: {bench.n_nrna_expected}")
    print(f"gDNA gt: {bench.n_gdna_expected}")
    print(f"mRNA pipe: {bench.total_observed:.1f}")
    print(f"nRNA pipe: {bench.n_nrna_pipeline:.1f}")
    print(f"gDNA pipe: {bench.n_gdna_pipeline:.1f}")
    print(f"nRNA rel err: {(bench.n_nrna_pipeline - bench.n_nrna_expected)/bench.n_nrna_expected:.4f}")

    print()
    print("=== nrna_frac prior (per transcript) ===")
    for i in range(est.num_transcripts):
        a = est.nrna_frac_alpha[i]
        b = est.nrna_frac_beta[i]
        mean = a / (a + b) if (a + b) > 0 else 0
        kappa = a + b
        print(f"  t{i}: alpha={a:.2f}  beta={b:.2f}  mean={mean:.4f}  kappa={kappa:.1f}")

    # What should nrna_frac be?
    # For t1: exons (2000,4000) + (8000,10000) = 4000 bp exonic, 4000 bp intron, 8000 bp span
    L_exonic = 4000
    L_genomic = 8000
    A_m, A_n = 100, 100
    true_nrna_frac = (A_n * L_genomic) / (A_m * L_exonic + A_n * L_genomic)
    density_nrna_frac = A_n / (A_m + A_n)
    print()
    print(f"True nrna_frac (fragment ratio): {true_nrna_frac:.4f}")
    print(f"Density nrna_frac (abundance ratio): {density_nrna_frac:.4f}")
    print(f"Gap: {true_nrna_frac - density_nrna_frac:.4f}")

    print()
    print("=== Pre-EM accumulators ===")
    for i in range(est.num_transcripts):
        print(f"  t{i}: exon_sense={est.transcript_exonic_sense[i]:.1f}"
              f"  exon_anti={est.transcript_exonic_antisense[i]:.1f}"
              f"  intron_sense={est.transcript_intronic_sense[i]:.1f}"
              f"  intron_anti={est.transcript_intronic_antisense[i]:.1f}")
    if est._exonic_lengths is not None:
        print(f"  L_exonic: {est._exonic_lengths}")
    if est._transcript_spans is not None:
        print(f"  L_span: {est._transcript_spans}")

    # Manual nrna_frac calculation
    print()
    print("=== Manual nrna_frac calculation ===")
    t_idx = 0
    D_exon = est.transcript_exonic_sense[t_idx] / est._exonic_lengths[t_idx]
    L_intron = est._transcript_spans[t_idx] - est._exonic_lengths[t_idx]
    D_intron = est.transcript_intronic_sense[t_idx] / L_intron if L_intron > 0 else 0
    D_nrna = max(0, D_intron)
    D_mrna = max(0, D_exon - D_intron)
    nrna_frac_density = D_nrna / (D_mrna + D_nrna) if (D_mrna + D_nrna) > 0 else 0
    print(f"  D_exon = {D_exon:.4f}")
    print(f"  D_intron = {D_intron:.4f}")
    print(f"  D_mrna = {D_mrna:.4f}")
    print(f"  D_nrna = {D_nrna:.4f}")
    print(f"  nrna_frac (density) = {nrna_frac_density:.4f}")

    # Fragment-ratio corrected
    L_ex = float(est._exonic_lengths[t_idx])
    L_sp = float(est._transcript_spans[t_idx])
    nrna_frac_frag = (D_nrna * L_sp) / (D_mrna * L_ex + D_nrna * L_sp) if (D_mrna * L_ex + D_nrna * L_sp) > 0 else 0
    print(f"  nrna_frac (fragment-corrected) = {nrna_frac_frag:.4f}")
    print(f"  true nrna_frac = {true_nrna_frac:.4f}")

    print()
    print("=== EM convergence ===")
    print(f"  nrna_em_counts: {est.nrna_em_counts}")
    total_em = float(est.em_counts[t_idx].sum() if hasattr(est.em_counts, '__getitem__') else 0)
    total_nrna = float(est.nrna_em_counts[t_idx])
    print(f"  mRNA EM: {total_em:.1f}")
    print(f"  nRNA EM: {total_nrna:.1f}")
    if total_em + total_nrna > 0:
        actual_nrna_frac = total_nrna / (total_em + total_nrna)
        print(f"  actual nrna_frac from EM: {actual_nrna_frac:.4f}")
