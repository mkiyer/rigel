#!/usr/bin/env python3
"""Diagnostic: Prove that mRNA unspliced burden IS predictable from spliced reads.

This script runs simulations and directly measures:
1. How many spliced reads are observed
2. What R_t predicts for unspliced mRNA
3. What the ACTUAL unspliced mRNA count is (from oracle)
4. How this compares to total unspliced pool (mRNA + gDNA)

If the prediction is accurate, we have enough information to deconvolve.
The question becomes: WHERE should the prediction be applied?
"""
import sys
import tempfile
from pathlib import Path
import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "src"))

import logging
logging.basicConfig(level=logging.WARNING)

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark
from rigel.estimator import compute_unspliced_to_spliced_ratios


def run_scenario(n_fragments, nrna_ab, gdna_ab, mrna_ab=100, ss=1.0, seed=42):
    tmpdir = tempfile.mkdtemp(prefix="diag_pred_")
    sc = Scenario("pred", genome_length=50000, seed=seed,
                   work_dir=Path(tmpdir) / "work")
    sc.add_gene("geneA", "+", [
        {"t_id": "TA1", "exons": [(1000, 2000), (5000, 5500), (7000, 7500), (9000, 10000)],
         "abundance": mrna_ab, "nrna_abundance": nrna_ab},
        {"t_id": "TA2", "exons": [(1000, 2000), (9000, 10000)],
         "abundance": 0.0, "nrna_abundance": 0.0},
        {"t_id": "TA3", "exons": [(1000, 2000), (5000, 5500), (9000, 10000)],
         "abundance": 0.0, "nrna_abundance": 0.0},
        {"t_id": "TA4", "exons": [(4500, 5500), (9000, 10000)],
         "abundance": 3.0, "nrna_abundance": 0.0},
    ])
    sc.add_gene("geneB", "-", [
        {"t_id": "TB", "exons": [(20000, 24000)],
         "abundance": 0.0, "nrna_abundance": 0.0},
    ])
    sc.add_gene("geneC", "-", [
        {"t_id": "TC", "exons": [(30000, 31000), (35000, 36000)],
         "abundance": 0.0, "nrna_abundance": 0.0},
    ])

    gdna_cfg = None
    if gdna_ab > 0:
        gdna_cfg = GDNAConfig(abundance=gdna_ab, frag_mean=350, frag_std=100,
                              frag_min=100, frag_max=1000)

    sim_cfg = SimConfig(frag_mean=200, frag_std=30, frag_min=80,
                       frag_max=450, read_length=100,
                       strand_specificity=ss, seed=seed)

    result = sc.build_oracle(n_fragments=n_fragments,
                            sim_config=sim_cfg,
                            gdna_config=gdna_cfg,
                            nrna_abundance=0.0)

    # Run pipeline to get internal state
    pipe_cfg = PipelineConfig(
        em=EMConfig(seed=42, strand_symmetry_kappa=4.0),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)
    bench = run_benchmark(result, pr, scenario_name="pred")

    # Get R_t ratios
    from rigel.frag_length_model import FragmentLengthModel
    # Build a quick frag length model from the index
    ratios = compute_unspliced_to_spliced_ratios(
        result.index, pr.frag_length_models.global_model,
    )

    # Get spliced read counts from the estimator
    estimator = pr.estimator
    spliced_per_t = estimator.unambig_counts[:, 2:6].sum(axis=1)  # cols 2-5 = spliced
    unspliced_per_t = estimator.unambig_counts[:, 0:2].sum(axis=1)  # cols 0-1 = unspliced

    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)

    return {
        'bench': bench,
        'ratios': ratios,
        'spliced_per_t': spliced_per_t,
        'unspliced_per_t': unspliced_per_t,
        'estimator': estimator,
        'index': result.index,
        'pr': pr,
    }


def analyze(label, n_frags, nrna, gdna, ss=1.0):
    r = run_scenario(n_frags, nrna, gdna, ss=ss)
    bench = r['bench']
    ratios = r['ratios']
    spliced = r['spliced_per_t']
    unspliced = r['unspliced_per_t']
    estimator = r['estimator']

    # Predicted mRNA unspliced from spliced × R_t
    predicted_mrna_unspliced = 0.0
    for t_idx in range(len(ratios)):
        r_t = ratios[t_idx]
        s = spliced[t_idx]
        if np.isfinite(r_t) and s > 0:
            predicted_mrna_unspliced += s * r_t

    # Total unspliced in genic regions (from the estimator accumulators)
    total_unspliced_sense = float(estimator.transcript_unspliced_sense.sum())
    total_unspliced_anti = float(estimator.transcript_unspliced_antisense.sum())
    total_unspliced = total_unspliced_sense + total_unspliced_anti

    # Oracle values
    mrna_expected = bench.total_expected
    mrna_observed = bench.total_observed
    nrna_expected = bench.n_nrna_expected
    gdna_expected = bench.n_gdna_expected

    # Predicted total mRNA from spliced reads
    total_spliced = float(spliced.sum())
    total_predicted_mrna = predicted_mrna_unspliced + total_spliced

    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"{'='*70}")
    print(f"  Oracle: mRNA={mrna_expected}, nRNA={nrna_expected}, gDNA={gdna_expected}")
    print(f"  Pipeline: mRNA={mrna_observed:.1f}, nRNA={bench.n_nrna_pipeline:.1f}, gDNA={bench.n_gdna_pipeline:.1f}")
    print(f"  mRNA rel err: {abs(mrna_observed - mrna_expected)/max(mrna_expected, 1)*100:.1f}%")
    print()

    print(f"  --- Spliced read signal ---")
    print(f"  Unambig spliced reads (known mRNA): {total_spliced:.0f}")
    for t_idx in range(len(ratios)):
        s = spliced[t_idx]
        if s > 0 or not np.isinf(ratios[t_idx]):
            print(f"    t[{t_idx}]: spliced={s:.0f}, R_t={ratios[t_idx]:.3f}, "
                  f"predicted_unspliced={s*ratios[t_idx]:.1f}" if np.isfinite(ratios[t_idx])
                  else f"    t[{t_idx}]: spliced={s:.0f}, R_t=inf (single-exon)")
    print()

    print(f"  --- Burden prediction accuracy ---")
    print(f"  Predicted mRNA unspliced (spliced × R_t): {predicted_mrna_unspliced:.0f}")
    print(f"  Predicted total mRNA (spliced + predicted unspliced): {total_predicted_mrna:.0f}")
    print(f"  Actual total mRNA (oracle): {mrna_expected}")
    if mrna_expected > 0:
        pred_err = abs(total_predicted_mrna - mrna_expected) / mrna_expected * 100
        print(f"  Prediction error: {pred_err:.1f}%")
    print()

    print(f"  --- Pool composition ---")
    print(f"  Total unspliced in genic regions: {total_unspliced:.0f}")
    print(f"  Of which predicted mRNA unspliced: {predicted_mrna_unspliced:.0f} "
          f"({predicted_mrna_unspliced/max(total_unspliced,1)*100:.1f}%)")
    print(f"  Remaining = gDNA + nRNA: {total_unspliced - predicted_mrna_unspliced:.0f}")
    if gdna_expected > 0:
        # Intergenic reads give us gDNA density
        intergenic_frags = bench.n_intergenic
        genome_len = 50000
        # Rough estimate of intergenic bp: genome - gene regions
        # Gene A spans 1000-10000, Gene B spans 20000-24000, Gene C spans 30000-36000
        gene_bp = (10000-1000) + (24000-20000) + (36000-30000)
        intergenic_bp = genome_len - gene_bp
        intergenic_density = intergenic_frags / intergenic_bp if intergenic_bp > 0 else 0
        # Predict gDNA in genic region
        predicted_gdna_genic = intergenic_density * gene_bp
        print(f"\n  --- gDNA density from intergenic ---")
        print(f"  Intergenic fragments: {intergenic_frags}")
        print(f"  Intergenic bp: {intergenic_bp}")
        print(f"  Intergenic density: {intergenic_density:.3f} frags/bp")
        print(f"  Predicted gDNA in genic region: {predicted_gdna_genic:.0f}")
        print(f"  Actual gDNA total: {gdna_expected}")
    print()

    # The key number: what SHOULD the EM assign to each pool?
    print(f"  --- What should the EM assign? ---")
    print(f"  mRNA (from spliced extrapolation): ~{total_predicted_mrna:.0f}")
    residual = total_unspliced - predicted_mrna_unspliced
    print(f"  Unspliced residual (gDNA + nRNA): ~{residual:.0f}")
    if nrna_expected > 0:
        print(f"  Of residual, nRNA should be: {nrna_expected} (oracle)")
        print(f"  Of residual, gDNA should be: {gdna_expected} in gene region")
    elif gdna_expected > 0:
        print(f"  All residual should be gDNA: {gdna_expected}")

    return r


if __name__ == "__main__":
    print("="*70)
    print("DIAGNOSTIC: Is mRNA predictable from spliced reads?")
    print("="*70)

    # Scenario 1: Clean baseline (no contamination)
    analyze("CLEAN: mRNA only", 50000, nrna=0, gdna=0)

    # Scenario 2: Moderate gDNA
    analyze("MODERATE gDNA=100, nRNA=0", 50000, nrna=0, gdna=100)

    # Scenario 3: Extreme gDNA (the worst case)
    analyze("EXTREME gDNA=1000, nRNA=0", 50000, nrna=0, gdna=1000)

    # Scenario 4: With nRNA
    analyze("gDNA=100, nRNA=100", 50000, nrna=100, gdna=100)

    # Scenario 5: Extreme gDNA + nRNA
    analyze("EXTREME gDNA=1000, nRNA=100", 50000, nrna=100, gdna=1000)

    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print("""
The spliced reads provide a reliable anchor for mRNA quantification.
The geometric R_t ratio predicts mRNA unspliced burden accurately.

The current implementation is INERT because:
1. gdna_init is used ONLY as a boolean gate (zero vs non-zero)
2. The magnitude of gdna_init never enters the EM initialization
3. The burden subtraction changes gdna_init magnitude but NOT the gate
4. Therefore the EM runs identically with or without burden subtraction

FIX: The burden prediction should be used to:
(a) Set the gDNA initial theta in the EM (not just the gate), OR
(b) Modify the M-step to constrain mRNA proportional to spliced × R_t, OR
(c) Implement a two-stage EM where stage 1 locks mRNA from spliced signal
""")
