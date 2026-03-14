#!/usr/bin/env python3
"""diagnostic_init_params.py — Instrument rigel initialization parameters.

Uses the same Scenario / nRNA-group infrastructure as synthetic_sim_sweep.py
to run a set of simulation scenarios and capture ALL internal initialization
parameters so we can compare them to ground truth values.

Captured parameters:
  - strand_specificity (estimated vs true)
  - nrna_init per nRNA component
  - gdna_init per locus (from EB shrinkage)
  - nrna_frac_alpha, nrna_frac_beta per nRNA
  - nrna_frac (derived η = α/(α+β))
  - locus_results (mrna, nrna, gdna allocated per locus)
  - intronic_sense/antisense per nRNA
  - exonic_sense/antisense per nRNA
  - unspliced_sense/antisense per nRNA
  - gdna_contamination_rate

Usage:
    conda run -n rigel python scripts/diagnostic_init_params.py
"""

import json
import logging
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

# Ensure src/ is importable when running from the scripts/ directory.
_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "src"))

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark
from synthetic_sim_sweep import parse_nrna_coords

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────
# Scenario definition (same layout as nrna_sweep_config.yaml)
# ──────────────────────────────────────────────────────────────────

GENOME_LENGTH = 50000
SEED = 42

TRANSCRIPTS = {
    "TA1": {"strand": "+", "exons": [[1000, 2000], [5000, 5500], [7000, 7500], [9000, 10000]]},
    "TA2": {"strand": "+", "exons": [[1000, 2000], [9000, 10000]]},
    "TA3": {"strand": "+", "exons": [[1000, 2000], [5000, 5500], [9000, 10000]]},
    "TA4": {"strand": "+", "exons": [[4500, 5500], [9000, 10000]]},
    "TB":  {"strand": "-", "exons": [[20000, 24000]]},
    "TC":  {"strand": "-", "exons": [[30000, 31000], [35000, 36000]]},
}
NRNAS_CONFIG = {"NTA": ["+", 1000, 10000]}

# Pre-compute nRNA groups (same logic as in the sweep script)
NRNA_GROUPS = parse_nrna_coords(NRNAS_CONFIG, TRANSCRIPTS)


def build_scenario(params):
    """Build a Scenario from parameter dict using sweep-style nRNA groups."""
    label = params.get("label", "diag")
    seed = params.get("seed", SEED)

    with tempfile.TemporaryDirectory(prefix=f"rigel_diag_{label}_") as tmpdir:
        sc = Scenario(
            name=label,
            genome_length=GENOME_LENGTH,
            seed=seed,
            work_dir=Path(tmpdir) / "work",
        )

        # Add transcripts grouped by nRNA group (gene_id = group label)
        for grp_label, grp in NRNA_GROUPS.items():
            nrna_ab = float(params.get(grp_label, 0))
            t_list = []
            for t_id in grp.t_ids:
                t_def = TRANSCRIPTS[t_id]
                ab = float(params.get(t_id, 0))
                nrna = nrna_ab if t_id == grp.carrier else 0.0
                t_list.append({
                    "t_id": t_id,
                    "exons": [tuple(e) for e in t_def["exons"]],
                    "abundance": ab,
                    "nrna_abundance": nrna,
                })
            sc.add_gene(grp_label, grp.strand, t_list)

        yield sc, tmpdir

# ──────────────────────────────────────────────────────────────────
# Diagnostic capture
# ──────────────────────────────────────────────────────────────────

def capture_diagnostics(params, label=""):
    """Run one scenario and capture all internal initialization parameters."""
    ss = params.get("strand_specificity", 1.0)
    n_frags = params.get("n_fragments", 100000)
    gdna_ab = float(params.get("gdna", 0))
    seed = params.get("seed", SEED)

    # Build scenario using sweep-style nRNA groups
    for sc, tmpdir in build_scenario(params):
        gdna_cfg = None
        if gdna_ab > 0:
            gdna_cfg = GDNAConfig(
                abundance=gdna_ab,
                frag_mean=350, frag_std=100, frag_min=100, frag_max=1000,
            )
        sim_cfg = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=ss, seed=seed,
        )
        result = sc.build_oracle(
            n_fragments=n_frags,
            sim_config=sim_cfg,
            gdna_config=gdna_cfg,
            nrna_abundance=0.0,
        )

        # Run pipeline
        pipe_cfg = PipelineConfig(
            em=EMConfig(
                seed=seed,
            ),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)

        # Run benchmark for ground truth
        bench = run_benchmark(result, pr, scenario_name=label)

        est = pr.estimator

        # ── Extract internal parameters ──
        diag = {"label": label}
        diag.update(params)

        # Ground truth
        diag["gt_gdna_fragments"] = bench.n_gdna_expected
        diag["gt_nrna_fragments"] = bench.n_nrna_expected
        diag["gt_mrna_fragments"] = bench.total_expected
        diag["gt_strand_specificity"] = ss

        # Pipeline results
        diag["pipe_gdna"] = round(bench.n_gdna_pipeline, 2)
        diag["pipe_nrna"] = round(bench.n_nrna_pipeline, 2)
        diag["pipe_mrna"] = round(bench.total_observed, 2)
        diag["mrna_rel_err"] = (round(abs(bench.total_observed - bench.total_expected) /
                                      bench.total_expected, 4)
                                if bench.total_expected > 0 else 0)
        diag["n_fragments_actual"] = bench.n_fragments
        diag["n_intergenic"] = bench.n_intergenic

        # Strand model
        diag["est_strand_specificity"] = round(pr.strand_models.strand_specificity, 6)

        # Per-nRNA init
        for i in range(est.num_nrna):
            diag[f"nrna_init_{i}"] = round(float(est.nrna_init[i]), 2)

        # Per-nRNA frac alpha/beta → η
        for i in range(est.num_nrna):
            alpha_i = float(est.nrna_frac_alpha[i])
            beta_i = float(est.nrna_frac_beta[i])
            eta_i = alpha_i / (alpha_i + beta_i) if (alpha_i + beta_i) > 0 else 0
            diag[f"nrna_frac_alpha_{i}"] = round(alpha_i, 4)
            diag[f"nrna_frac_beta_{i}"] = round(beta_i, 4)
            diag[f"nrna_frac_eta_{i}"] = round(eta_i, 4)

        # Per-nRNA intronic evidence
        for i in range(est.num_nrna):
            diag[f"intronic_sense_{i}"] = round(float(est.transcript_intronic_sense[i]), 2)
            diag[f"intronic_anti_{i}"] = round(float(est.transcript_intronic_antisense[i]), 2)

        # Per-nRNA exonic evidence (transcript-level; aggregate to nRNA)
        n_t = est.num_transcripts
        n_nrna = est.num_nrna
        t_to_nrna = est._t_to_nrna
        if t_to_nrna is not None and n_nrna > 0:
            nrna_exon_s = np.zeros(n_nrna, dtype=np.float64)
            nrna_exon_a = np.zeros(n_nrna, dtype=np.float64)
            np.add.at(nrna_exon_s, t_to_nrna, est.transcript_exonic_sense)
            np.add.at(nrna_exon_a, t_to_nrna, est.transcript_exonic_antisense)
            for i in range(n_nrna):
                diag[f"exonic_sense_{i}"] = round(float(nrna_exon_s[i]), 2)
                diag[f"exonic_anti_{i}"] = round(float(nrna_exon_a[i]), 2)

        # Per-nRNA unspliced evidence
        for i in range(min(n_nrna, 6)):
            diag[f"unspliced_sense_{i}"] = round(float(est.transcript_unspliced_sense[i]), 2)
            diag[f"unspliced_anti_{i}"] = round(float(est.transcript_unspliced_antisense[i]), 2)

        # Locus results (gdna_init, mrna, nrna, gdna per locus)
        for i, lr in enumerate(est.locus_results):
            diag[f"locus_{i}_gdna_init"] = round(lr.get("gdna_init", 0), 2)
            diag[f"locus_{i}_mrna"] = round(lr.get("mrna", 0), 2)
            diag[f"locus_{i}_nrna"] = round(lr.get("nrna", 0), 2)
            diag[f"locus_{i}_gdna"] = round(lr.get("gdna", 0), 2)
            diag[f"locus_{i}_n_frags"] = lr.get("n_em_fragments", 0)

        # gDNA contamination rate
        diag["gdna_contamination_rate"] = round(est.gdna_contamination_rate, 6)

        # Compute what ground truth nrna_frac should be
        gt_nrna = bench.n_nrna_expected
        gt_mrna = bench.total_expected
        if (gt_nrna + gt_mrna) > 0:
            diag["gt_nrna_frac"] = round(gt_nrna / (gt_nrna + gt_mrna), 4)
        else:
            diag["gt_nrna_frac"] = 0.0

        return diag


# ──────────────────────────────────────────────────────────────────
# Test scenarios
# ──────────────────────────────────────────────────────────────────

SCENARIOS = [
    # Baseline: no contamination
    {"label": "baseline_clean",
     "TA1": 100, "TA4": 3, "NTA": 0, "gdna": 0,
     "strand_specificity": 1.0, "n_fragments": 100000,

    # Moderate nRNA, no gDNA
    {"label": "nrna_only_100",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 0,
     "strand_specificity": 1.0, "n_fragments": 100000,

    # Moderate gDNA, no nRNA
    {"label": "gdna_only_100_k2",
     "TA1": 100, "TA4": 3, "NTA": 0, "gdna": 100,
     "strand_specificity": 1.0, "n_fragments": 100000,

    {"label": "gdna_only_100_k6",
     "TA1": 100, "TA4": 3, "NTA": 0, "gdna": 100,
     "strand_specificity": 1.0, "n_fragments": 100000,

    # Both moderate
    {"label": "both_100_k2",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 100,
     "strand_specificity": 1.0, "n_fragments": 100000,

    {"label": "both_100_k6",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 100,
     "strand_specificity": 1.0, "n_fragments": 100000,

    # Extreme gDNA
    {"label": "gdna_1000_k2",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 1000,
     "strand_specificity": 1.0, "n_fragments": 100000,

    {"label": "gdna_1000_k6",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 1000,
     "strand_specificity": 1.0, "n_fragments": 100000,

    # Very extreme gDNA
    {"label": "gdna_5000_k2",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 5000,
     "strand_specificity": 1.0, "n_fragments": 100000,

    {"label": "gdna_5000_k6",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 5000,
     "strand_specificity": 1.0, "n_fragments": 100000,

    # Extreme nRNA
    {"label": "nrna_1000_k2",
     "TA1": 100, "TA4": 3, "NTA": 1000, "gdna": 100,
     "strand_specificity": 1.0, "n_fragments": 100000,

    {"label": "nrna_1000_k6",
     "TA1": 100, "TA4": 3, "NTA": 1000, "gdna": 100,
     "strand_specificity": 1.0, "n_fragments": 100000,

    # Asymmetric gDNA (SS=0.9) — the penalty should detect imbalance
    {"label": "gdna_1000_ss09_k2",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 1000,
     "strand_specificity": 0.9, "n_fragments": 100000,

    {"label": "gdna_1000_ss09_k6",
     "TA1": 100, "TA4": 3, "NTA": 100, "gdna": 1000,
     "strand_specificity": 0.9, "n_fragments": 100000,

    # Both extreme
    {"label": "both_5000_k6",
     "TA1": 100, "TA4": 3, "NTA": 5000, "gdna": 5000,
     "strand_specificity": 1.0, "n_fragments": 200000,
]


def main():
    out_dir = "sweep_results/diagnostics"
    os.makedirs(out_dir, exist_ok=True)

    results = []
    for i, params in enumerate(SCENARIOS):
        label = params.get("label", f"run_{i}")
        print(f"\n{'='*70}")
        print(f"[{i+1}/{len(SCENARIOS)}] {label}")
        print(f"  TA1={params.get('TA1',0)}, TA4={params.get('TA4',0)}, "
              f"NTA={params.get('NTA',0)}, gDNA={params.get('gdna',0)}, "
              f"SS={params.get('strand_specificity',1.0)}")

        diag = capture_diagnostics(params, label=label)
        results.append(diag)

        # Print key metrics
        print(f"  Results: mRNA_err={diag['mrna_rel_err']:.1%}, "
              f"pipe_gdna={diag['pipe_gdna']:.0f} (gt={diag['gt_gdna_fragments']}), "
              f"pipe_nrna={diag['pipe_nrna']:.0f} (gt={diag['gt_nrna_fragments']})")
        if 'nrna_frac_eta_0' in diag:
            print(f"  nRNA frac: eta={diag['nrna_frac_eta_0']:.4f} "
                  f"(gt={diag['gt_nrna_frac']:.4f}), "
                  f"alpha={diag['nrna_frac_alpha_0']:.4f}, "
                  f"beta={diag['nrna_frac_beta_0']:.4f}")
        if 'nrna_init_0' in diag:
            print(f"  nRNA init: {diag['nrna_init_0']:.1f} "
                  f"(gt nRNA frags={diag['gt_nrna_fragments']})")
        if 'locus_0_gdna_init' in diag:
            print(f"  Locus 0: gdna_init={diag['locus_0_gdna_init']:.1f}, "
                  f"mrna={diag['locus_0_mrna']:.1f}, "
                  f"nrna={diag['locus_0_nrna']:.1f}, "
                  f"gdna={diag['locus_0_gdna']:.1f}")
        if 'intronic_sense_0' in diag:
            print(f"  Intronic: sense={diag['intronic_sense_0']:.0f}, "
                  f"anti={diag['intronic_anti_0']:.0f}")
        print(f"  Est SS: {diag['est_strand_specificity']:.6f}")

    # Save JSON
    json_path = os.path.join(out_dir, "diagnostics.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n\nSaved {len(results)} diagnostic records to {json_path}")

    # Print summary table
    print(f"\n{'='*120}")
    print("DIAGNOSTIC SUMMARY")
    print(f"{'='*120}")
    print(f"{'label':>25s} | {'mRNA%':>7s} | {'pipe_gDNA':>10s} | {'gt_gDNA':>8s} | "
          f"{'pipe_nRNA':>10s} | {'gt_nRNA':>8s} | {'eta':>6s} | {'gt_eta':>7s} | "
          f"{'nrna_init':>10s} | {'gdna_init':>10s} | {'est_SS':>8s}")
    print("-" * 120)
    for d in results:
        eta = d.get('nrna_frac_eta_0', 0)
        gt_eta = d.get('gt_nrna_frac', 0)
        ninit = d.get('nrna_init_0', 0)
        ginit = d.get('locus_0_gdna_init', 0)
        print(f"{d['label']:>25s} | {d['mrna_rel_err']:>6.1%} | {d['pipe_gdna']:>10.0f} | "
              f"{d['gt_gdna_fragments']:>8d} | {d['pipe_nrna']:>10.1f} | "
              f"{d['gt_nrna_fragments']:>8d} | {eta:>6.4f} | {gt_eta:>7.4f} | "
              f"{ninit:>10.1f} | {ginit:>10.1f} | {d['est_strand_specificity']:>8.6f}")

    # Deep-dive: nrna_frac error analysis
    print(f"\n{'='*120}")
    print("nRNA FRAC PRIOR vs GROUND TRUTH")
    print(f"{'='*120}")
    print(f"{'label':>25s} | {'alpha':>8s} | {'beta':>8s} | {'eta_prior':>10s} | "
          f"{'gt_eta':>8s} | {'eta_err':>8s} | {'int_s':>7s} | {'int_a':>7s} | "
          f"{'ex_s':>7s} | {'ex_a':>7s} | {'nrna_init':>10s}")
    print("-" * 120)
    for d in results:
        alpha = d.get('nrna_frac_alpha_0', 0)
        beta = d.get('nrna_frac_beta_0', 0)
        eta = d.get('nrna_frac_eta_0', 0)
        gt_eta = d.get('gt_nrna_frac', 0)
        eta_err = eta - gt_eta
        int_s = d.get('intronic_sense_0', 0)
        int_a = d.get('intronic_anti_0', 0)
        ex_s = d.get('exonic_sense_0', 0)
        ex_a = d.get('exonic_anti_0', 0)
        ninit = d.get('nrna_init_0', 0)
        print(f"{d['label']:>25s} | {alpha:>8.4f} | {beta:>8.4f} | {eta:>10.4f} | "
              f"{gt_eta:>8.4f} | {eta_err:>+8.4f} | {int_s:>7.0f} | {int_a:>7.0f} | "
              f"{ex_s:>7.0f} | {ex_a:>7.0f} | {ninit:>10.1f}")


if __name__ == "__main__":
    main()
