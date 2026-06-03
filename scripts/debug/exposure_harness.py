#!/usr/bin/env python
"""Exposure-prior validation harness — the three regimes (count/EM Mechanism A).

Establishes baselines for the adaptive exposure-prior work
(docs/futureprs/exposure_prior_robustness_theory.md):

  1. SPARSE + phi-downsampling sweep   — the reviewer's acceptance diagnostic:
     recovered exposure-dispersion phi vs total fragments, plus the per-locus
     gDNA-eff-len min and RNA recovery. Pre-fix: phi ~flat, RNA collapses when
     sparse. Post-fix: phi should inflect toward 0, ω→1, RNA graceful.
  2. DENSE-UNIFORM (scale-up)          — many regions, uniform gDNA. Robust by the
     law of large numbers; the fix must not regress it.
  3. DENSE-CAPTURE (synthetic probes)  — gDNA focally enriched at probe targets.
     The critical guard: the prior must RELAX so real enrichment (ω>1 at targets)
     is recovered, not shrunk to ω=1.

Usage:  python scripts/debug/exposure_harness.py [--regime sparse|dense|capture|all]
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import numpy as np

from rigel.calibration.region_arrays import RegionArrays
from rigel.config import BamScanConfig, CalibrationConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import CaptureConfig, GDNAConfig, ReadSimConfig, Scenario, run_benchmark

_KAPPA: float | None = None  # exposure_prior_pseudocount override (set by --kappa)

_T1_EXONS = [(2000, 2500), (3000, 3500), (4000, 4500), (5000, 5500),
             (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)]
_TC_EXONS = [(12000, 12500), (13000, 13500), (14000, 14500), (15000, 15500),
             (16000, 16500), (17000, 17500), (18000, 18500), (19000, 19500)]


def _pipeline(res):
    cal = CalibrationConfig() if _KAPPA is None else CalibrationConfig(exposure_prior_pseudocount=_KAPPA)
    return run_pipeline(
        res.bam_path, res.index,
        config=PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="auto"),
                              calibration=cal),
    )


def _metrics(pr, res, label):
    b = run_benchmark(res, pr, scenario_name=label)
    loci = pr.estimator.get_loci_df(res.index)
    eff = loci["gdna_eff_len_em"].to_numpy() if "gdna_eff_len_em" in loci.columns else np.array([np.nan])
    rna_exp = b.total_expected + b.n_nrna_expected
    return dict(
        phi=pr.calibration.exposure_dispersion, conv=pr.calibration.converged,
        eff_min=float(eff.min()), rna_obs=b.total_rna_observed, rna_exp=rna_exp,
        gdna=b.n_gdna_pipeline, gdna_exp=b.n_gdna_expected, nrna=b.n_nrna_pipeline,
    )


# ---------------------------------------------------------------- regime 1: sparse sweep
def sparse_sweep():
    print("\n=== REGIME 1: SPARSE — phi-vs-downsampling sweep (nrna_dc g20 n70 ss65) ===")
    print(f"{'n_frag':>7} {'phi':>8} {'conv':>5} {'eff_min':>8} {'RNAobs/exp':>14} {'gDNA(exp596)':>13} {'nRNA':>6}")
    for nfrag in (500, 1000, 2000, 4000, 8000, 16000):
        tmp = tempfile.mkdtemp()
        sc = Scenario("sp", genome_length=20000, seed=42, work_dir=tmp)
        sc.add_gene("g1", "+", [{"t_id": "t1", "exons": _T1_EXONS, "abundance": 100}])
        sc.add_gene("gc", "-", [{"t_id": "tc", "exons": _TC_EXONS, "abundance": 0}])
        cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                            read_length=100, strand_specificity=0.65, seed=42)
        res = sc.build_oracle(n_fragments=nfrag, sim_config=cfg,
                              gdna_config=GDNAConfig(abundance=20, frag_mean=350, frag_std=100,
                                                     frag_min=100, frag_max=1000),
                              nrna_abundance=70)
        m = _metrics(_pipeline(res), res, "sp")
        rec = m["rna_obs"] / m["rna_exp"] if m["rna_exp"] else 0
        print(f"{nfrag:>7} {m['phi']:>8.3g} {str(m['conv']):>5} {m['eff_min']:>8.0f} "
              f"{m['rna_obs']:>6.0f}/{m['rna_exp']:<6.0f}({rec:.0%}) {m['gdna']:>13.0f} {m['nrna']:>6.0f}")
    print("  Target post-fix: phi inflects toward 0 as n_frag drops; RNA recovery stays graceful (no collapse).")


# ---------------------------------------------------------------- regime 2: dense-uniform
def dense_uniform():
    print("\n=== REGIME 2: DENSE-UNIFORM — scale-up (K gene-pairs, uniform gDNA, g20 n70 ss65) ===")
    print(f"{'K':>3} {'regions':>8} {'phi':>8} {'eff_min':>8} {'RNAobs/exp':>16} {'gDNA':>8}")
    for K in (1, 5, 20):
        tmp = tempfile.mkdtemp()
        sc = Scenario("du", genome_length=K * 20000, seed=42, work_dir=tmp)
        for i in range(K):
            o = i * 20000
            sc.add_gene(f"g1_{i}", "+", [{"t_id": f"t1_{i}", "exons": [(o + a, o + b) for a, b in _T1_EXONS], "abundance": 100}])
            sc.add_gene(f"gc_{i}", "-", [{"t_id": f"tc_{i}", "exons": [(o + a, o + b) for a, b in _TC_EXONS], "abundance": 0}])
        cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                            read_length=100, strand_specificity=0.65, seed=42)
        res = sc.build_oracle(n_fragments=K * 2000, sim_config=cfg,
                              gdna_config=GDNAConfig(abundance=20, frag_mean=350, frag_std=100,
                                                     frag_min=100, frag_max=1000),
                              nrna_abundance=70)
        pr = _pipeline(res)
        m = _metrics(pr, res, "du")
        rec = m["rna_obs"] / m["rna_exp"] if m["rna_exp"] else 0
        print(f"{K:>3} {pr.calibration.n_regions:>8} {m['phi']:>8.3g} {m['eff_min']:>8.0f} "
              f"{m['rna_obs']:>7.0f}/{m['rna_exp']:<7.0f}({rec:.0%}) {m['gdna']:>8.0f}")
    print("  Robust baseline (law of large numbers); the fix must not regress it.")


# ---------------------------------------------------------------- regime 3: dense-capture
def dense_capture():
    print("\n=== REGIME 3: DENSE-CAPTURE — synthetic probes, gDNA focally enriched at targets ===")
    tmp = tempfile.mkdtemp()
    sc = Scenario("cap", genome_length=120000, seed=42, work_dir=tmp)
    # 6 expressed genes; CAPTURE the first 3 (probes over their full transcripts).
    captured = {0, 1, 2}
    for i in range(6):
        o = i * 20000
        sc.add_gene(f"g{i}", "+", [{"t_id": f"t{i}", "exons": [(o + a, o + b) for a, b in _T1_EXONS], "abundance": 50}])
    probe_path = Path(tmp) / "probes.tsv"
    with open(probe_path, "w") as fh:
        fh.write("transcript_id\tstart\tend\n")
        for i in sorted(captured):
            fh.write(f"t{i}\t0\t4000\n")  # whole 8-exon (4000 bp) transcript captured
    cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                        read_length=100, strand_specificity=0.9, seed=42)
    res = sc.build_oracle(n_fragments=12000, sim_config=cfg,
                          gdna_config=GDNAConfig(abundance=40, frag_mean=350, frag_std=100,
                                                 frag_min=100, frag_max=1000),
                          nrna_abundance=0,
                          capture_config=CaptureConfig(probes=str(probe_path), binding_per_base=30.0))
    pr = _pipeline(res)
    m = _metrics(pr, res, "cap")
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    omega = pr.calibration.omega
    cap_starts = {i * 20000 for i in captured}
    # classify each region's owning gene-block start
    cap_mask = np.array([(int(s) // 20000) in captured for s in ra.start])
    print(f"  phi={m['phi']:.3g}  converged={m['conv']}  gDNA obs/exp={m['gdna']:.0f}/{m['gdna_exp']:.0f}")
    print(f"  omega at CAPTURED-gene regions:   mean={omega[cap_mask].mean():.3f}  "
          f"median={np.median(omega[cap_mask]):.3f}  max={omega[cap_mask].max():.3f}")
    print(f"  omega at UNCAPTURED-gene regions: mean={omega[~cap_mask].mean():.3f}  "
          f"median={np.median(omega[~cap_mask]):.3f}  max={omega[~cap_mask].max():.3f}")
    print("  Critical guard: captured-region omega should be >1 (real enrichment recovered),")
    print("  and the fix must NOT shrink it to 1 (the prior must relax when evidence is dense).")


def main() -> None:
    global _KAPPA
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--regime", choices=["sparse", "dense", "capture", "all"], default="all")
    ap.add_argument("--kappa", type=float, default=None,
                    help="override exposure_prior_pseudocount (default: config default)")
    args = ap.parse_args()
    _KAPPA = args.kappa
    print(f"[exposure_prior_pseudocount = {'default' if _KAPPA is None else _KAPPA}]")
    if args.regime in ("sparse", "all"):
        sparse_sweep()
    if args.regime in ("dense", "all"):
        dense_uniform()
    if args.regime in ("capture", "all"):
        dense_capture()


if __name__ == "__main__":
    main()
