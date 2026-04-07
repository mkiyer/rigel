"""Diagnose the prior imbalance in nRNA test scenarios.

Runs the nRNA double-counting scenario and inspects calibration
output + locus priors to understand the alpha_gdna / alpha_rna
imbalance.
"""
import numpy as np
from pathlib import Path
import tempfile
import logging

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import Scenario, SimConfig
from rigel.locus import compute_locus_priors, build_loci

logging.basicConfig(level=logging.WARNING)

SEED = 42
N_FRAGMENTS = 2000

# Test conditions that fail + passing controls
CONDITIONS = [
    (30, 0.65),
    (30, 0.90),
    (70, 0.65),
    (70, 0.90),
    (30, 1.00),
    (70, 1.00),
]


def run_diagnostic(nrna_abundance, strand_specificity):
    """Run scenario and inspect calibration + priors."""
    print(f"\n{'='*70}")
    print(f"nRNA={nrna_abundance}, SS={strand_specificity}")
    print(f"{'='*70}")

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        sc = Scenario(
            "diag", genome_length=20000, seed=SEED,
            work_dir=tmp_path / "diag",
        )
        sc.add_gene("g1", "+", [{
            "t_id": "t1",
            "exons": [(2000, 4000), (8000, 10000)],
            "abundance": 100,
        }])
        sc.add_gene("g_ctrl", "-", [{
            "t_id": "t_ctrl",
            "exons": [(14000, 16000), (18000, 19000)],
            "abundance": 0,
        }])

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=strand_specificity, seed=SEED,
        )
        result = sc.build_oracle(
            n_fragments=N_FRAGMENTS, sim_config=sim_config,
            gdna_config=None, nrna_abundance=nrna_abundance,
        )

        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )

        # Run full pipeline to get results + capture calibration
        pr = run_pipeline(result.bam_path, result.index, config=config)

        # Extract calibration from pipeline result
        cal = pr.calibration
        print(f"\nCalibration:")
        print(f"  lambda_gdna:       {cal.lambda_gdna:.6f}")
        print(f"  strand_specificity: {cal.strand_specificity:.4f}")
        print(f"  total E[gDNA]:     {cal.region_e_gdna.sum():.1f}")
        print(f"  total N_total:     {cal.region_n_total.sum():.1f}")
        global_gamma = cal.region_e_gdna.sum() / max(cal.region_n_total.sum(), 1.0)
        print(f"  global gamma:      {global_gamma:.4f}")

        # Per-region breakdown (only non-zero)
        print(f"\n  Active regions (N > 0):")
        for i in range(len(cal.region_e_gdna)):
            e = cal.region_e_gdna[i]
            n = cal.region_n_total[i]
            if n > 0:
                gamma = e / n
                print(f"    region {i}: E[gDNA]={e:.1f}, N={n:.0f}, gamma={gamma:.4f}")

        # Compute locus priors: fraction-based C=1.0
        print(f"\n  Prior comparison (per total):")
        total_e = cal.region_e_gdna.sum()
        total_n = cal.region_n_total.sum()
        gamma_global = total_e / max(total_n, 1.0)

        # Current system (fraction-based, C=1.0)
        ag_cur = gamma_global * 1.0
        ar_cur = (1 - gamma_global) * 1.0
        print(f"    FRACTION C=1: alpha_gdna={ag_cur:.4f}, alpha_rna={ar_cur:.4f}, ratio={ag_cur/max(ar_cur,1e-30):.4f}")

        # Pipeline results
        est = pr.estimator
        print(f"\n  Pipeline output:")
        gdna_total = float(est.gdna_locus_counts.sum())
        total_mrna = float(est.t_counts.sum())
        print(f"    gDNA (EM):  {gdna_total:.1f}")
        print(f"    total mRNA: {total_mrna:.1f}")
        if gdna_total > 10:
            print(f"    !! {gdna_total:.0f} frags classified as gDNA in a ZERO gDNA scenario!")

        sc.cleanup()


if __name__ == "__main__":
    for nrna, ss in CONDITIONS:
        run_diagnostic(nrna, ss)
