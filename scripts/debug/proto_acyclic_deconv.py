#!/usr/bin/env python
"""Prototype 1 (headline): strand-only acyclic deconvolution — show the collapse is GONE.

Computes the per-region gDNA mass Ĝ_r from a *local strand split alone* (no global ρ₀
input), derives ρ₀ = ΣĜ/ΣL, and compares to (a) the truth and (b) the production EM's
collapsed ρ₀. If the design doc is right, Ĝ ≈ the true gDNA and ρ₀ is healthy — with no
feedback loop able to collapse it.

Deconvolution (frame-agnostic): gDNA is unstranded (50/50), so all sense/antisense
*imbalance* is RNA-driven. Per region (unspliced flux n₊, n₋; N=n₊+n₋):

    R = |n₊ − n₋| / |2κ_rna − 1|        (RNA count from the strand asymmetry)
    G = clip(N − R, 0, N),  π_g = G/N    (gDNA fraction)        Ĝ_mass = π_g · mass_unspliced

TS_NONE (intergenic) → G=N (gDNA by annotation). TS_AMBIG → undecodable (skip; would impute).

Usage:  python scripts/debug/proto_acyclic_deconv.py
"""

from __future__ import annotations

import tempfile

import numpy as np

from rigel.calibration.density import estimate_gdna_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG, TS_NONE, coarse_type_array
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark

_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}
_TS = {0: "NONE", 1: "POS", -1: "NEG", 2: "AMBIG"}


def build_nrna_dc(tmp):
    sc = Scenario("nrna_dc", genome_length=20000, seed=42, work_dir=tmp)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(2000, 2500), (3000, 3500),
        (4000, 4500), (5000, 5500), (6000, 6500), (7000, 7500), (8000, 8500),
        (9500, 10000)], "abundance": 100}])
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(12000, 12500),
        (13000, 13500), (14000, 14500), (15000, 15500), (16000, 16500),
        (17000, 17500), (18000, 18500), (19000, 19500)], "abundance": 0}])
    cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                        read_length=100, strand_specificity=0.65, seed=42)
    return sc.build_oracle(n_fragments=2000, sim_config=cfg,
                           gdna_config=GDNAConfig(abundance=20, frag_mean=350, frag_std=100,
                                                  frag_min=100, frag_max=1000),
                           nrna_abundance=70)


def main():
    tmp = tempfile.mkdtemp(prefix="proto_")
    res = build_nrna_dc(tmp)

    # --- raw substrate (per-region sense/antisense flux) from one scan ---
    scan = BamScanConfig(sj_strand_tag="auto")
    _, strand_models, _, _, payload = scan_and_buffer(str(res.bam_path), res.index, scan)
    strand_models.finalize()
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(payload, ra)
    kappa = fit_strand_balance(strand_models).kappa_rna
    power = abs(2.0 * kappa - 1.0)

    c = sub.contained
    pos = c.n_unspliced_pos.astype(np.float64)
    neg = c.n_unspliced_neg.astype(np.float64)
    N = pos + neg
    ts = sub.ts_class
    L = sub.L_eff

    # strand-only deconvolution (NO global ρ₀ input)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.abs(pos - neg) / power
    pi_g = np.where(N > 0, np.clip((N - R) / np.maximum(N, 1e-12), 0.0, 1.0), 0.0)
    pi_g = np.where(ts == TS_NONE, 1.0, pi_g)          # intergenic = gDNA by annotation
    decodable = ts != TS_AMBIG
    g_mass = pi_g * c.mass_unspliced
    g_mass = np.where(decodable, g_mass, 0.0)

    rho0_strand = float(g_mass.sum()) / float(L[decodable].sum())
    rho0_seed = estimate_gdna_density(sub, ra).rho_0

    # --- production EM (for the collapsed ρ₀ + the expected gDNA) ---
    pr = run_pipeline(res.bam_path, res.index,
                      config=PipelineConfig(em=EMConfig(seed=42), scan=scan))
    b = run_benchmark(res, pr, scenario_name="nrna_dc")
    rtype = coarse_type_array(ra.signature)

    print(f"\n=== Prototype 1: strand-only acyclic deconvolution (nrna_dc g20 n70 ss65) ===")
    print(f"kappa_rna = {kappa:.4f}   strand power |2k-1| = {power:.4f}")
    print(f"\n  idx ts    rtype       start  N(flux)  pos  neg   R     pi_g  G_mass  gene")
    for i in range(sub.n_regions):
        if N[i] < 0.5 and c.mass_unspliced[i] < 0.5:
            continue
        gene = "t1" if ra.start[i] < 11000 else ("t_ctrl" if ra.start[i] < 19500 else "")
        print(f"  {i:3d} {_TS.get(int(ts[i]),'?'):>5s} {_RTYPE.get(int(rtype[i]),'?'):>9s}"
              f" {int(ra.start[i]):6d} {N[i]:7.0f} {pos[i]:4.0f} {neg[i]:4.0f} {R[i]:5.1f}"
              f"  {pi_g[i]:.2f}  {g_mass[i]:6.1f}  {gene}")

    print(f"\n--- aggregate ---")
    print(f"  Σ Ĝ_mass (strand-deconvolved gDNA) = {g_mass.sum():7.1f}   expected gDNA = {b.n_gdna_expected:.0f}")
    print(f"  ρ₀ (strand-only, acyclic)  = {rho0_strand:.5g}")
    print(f"  ρ₀ (density seed)          = {rho0_seed:.5g}")
    print(f"  ρ₀ (production EM)         = {pr.calibration.rho_0:.5g}   <- the COLLAPSE (true ≈ 0.03)")
    print(f"\n  HEADLINE: strand-only ρ₀ should be healthy (~0.03), NOT collapsed to ~5e-5,")
    print(f"           and Σ Ĝ should be ≈ {b.n_gdna_expected:.0f} — with no global-ρ₀ feedback at all.")


if __name__ == "__main__":
    main()
