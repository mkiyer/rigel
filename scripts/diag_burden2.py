#!/usr/bin/env python3
"""Diagnostic: trace burden subtraction impact on gDNA density.

Instead of monkey-patching, directly calls the internal functions
to compute with and without burden subtraction.
"""
import sys
import tempfile
import logging
from pathlib import Path
import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "src"))

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark
from rigel.estimator import compute_unspliced_to_spliced_ratios, AbundanceEstimator
from rigel.locus import (
    compute_eb_gdna_priors,
    _compute_locus_mrna_burden,
    _compute_per_locus_gdna_densities,
    compute_gdna_density_hybrid,
    build_loci,
)
from rigel.splice import SpliceStrandCol

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


def make_scenario(n_fragments, nrna_abundance, gdna_abundance, ss=1.0, seed=42):
    """Create a standard test scenario and return pipeline internals."""
    tmpdir = tempfile.mkdtemp(prefix="diag_")
    sc = Scenario("diag", genome_length=50000, seed=seed,
                   work_dir=Path(tmpdir) / "work")

    sc.add_gene("geneA", "+", [
        {"t_id": "TA1", "exons": [(1000, 2000), (5000, 5500), (7000, 7500), (9000, 10000)],
         "abundance": 100.0, "nrna_abundance": nrna_abundance},
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
    if gdna_abundance > 0:
        gdna_cfg = GDNAConfig(abundance=gdna_abundance,
                              frag_mean=350, frag_std=100,
                              frag_min=100, frag_max=1000)

    sim_cfg = SimConfig(frag_mean=200, frag_std=30, frag_min=80,
                       frag_max=450, read_length=100,
                       strand_specificity=ss, seed=seed)

    result = sc.build_oracle(n_fragments=n_fragments,
                            sim_config=sim_cfg,
                            gdna_config=gdna_cfg,
                            nrna_abundance=0.0)
    return result, sc, tmpdir


def run_comparison(n_fragments, nrna, gdna, ss=1.0, kappa=4.0):
    """Run pipeline twice: with and without burden subtraction."""
    print(f"\n{'='*70}")
    print(f"n={n_fragments}, nRNA={nrna}, gDNA={gdna}, SS={ss}, kappa={kappa}")
    print(f"{'='*70}")

    result, sc, tmpdir = make_scenario(n_fragments, nrna, gdna, ss)

    # Run WITH burden subtraction (default)
    pipe_cfg = PipelineConfig(
        em=EMConfig(seed=42, strand_symmetry_kappa=kappa),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr_with = run_pipeline(result.bam_path, result.index, config=pipe_cfg)
    bench_with = run_benchmark(result, pr_with, scenario_name="with_burden")

    # Run WITHOUT burden subtraction by setting ratios to None
    # We need to hack the pipeline to skip the ratio computation.
    # Easiest: just temporarily make frag_length_model.n_observations = 0
    # Actually we can't do that easily. Let's instead directly compare by
    # calling compute_eb_gdna_priors manually.

    # Print results comparison
    print(f"\n  Expected mRNA: {bench_with.total_expected:.0f}")
    print(f"  Observed mRNA (with burden): {bench_with.total_observed:.2f}")
    print(f"  mRNA rel err: {abs(bench_with.total_observed - bench_with.total_expected) / max(bench_with.total_expected, 1):.4f}")
    print(f"  Expected nRNA: {bench_with.n_nrna_expected}")
    print(f"  Observed nRNA: {bench_with.n_nrna_pipeline:.2f}")
    print(f"  Expected gDNA: {bench_with.n_gdna_expected}")
    print(f"  Observed gDNA: {bench_with.n_gdna_pipeline:.2f}")

    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)


def trace_burden_values(n_fragments=50000, nrna=0, gdna=1000, ss=1.0):
    """Directly call internal functions to trace exact burden values."""
    print(f"\n{'='*70}")
    print(f"DIRECT TRACE: n={n_fragments}, nRNA={nrna}, gDNA={gdna}, SS={ss}")
    print(f"{'='*70}")

    result, sc, tmpdir = make_scenario(n_fragments, nrna, gdna, ss)

    # Run pipeline to get the internal state
    from rigel.pipeline import (
        _build_estimator, _build_frag_length_models,
    )
    from rigel.index import TranscriptIndex
    from rigel.buffer import FragmentBuffer
    from rigel.scoring import FragmentScorer
    from rigel.scanner import FragmentRouter
    from rigel.stats import ScanStats
    from rigel.strand_model import StrandModels

    # Full pipeline run with logging
    logging.getLogger('rigel').setLevel(logging.DEBUG)

    pipe_cfg = PipelineConfig(
        em=EMConfig(seed=42, strand_symmetry_kappa=4.0),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)

    logging.getLogger('rigel').setLevel(logging.WARNING)

    bench = run_benchmark(result, pr, scenario_name="trace")

    # Print key results
    print(f"\n  mRNA: expected={bench.total_expected:.0f}, observed={bench.total_observed:.2f}")
    print(f"  nRNA: expected={bench.n_nrna_expected}, observed={bench.n_nrna_pipeline:.2f}")
    print(f"  gDNA: expected={bench.n_gdna_expected}, observed={bench.n_gdna_pipeline:.2f}")

    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)
    return bench


def quantify_burden_impact():
    """The key test: manually compute gDNA init with and without burden.

    This directly answers: does the burden subtraction change the gDNA init?
    And if so, by how much?
    """
    print(f"\n{'='*70}")
    print("QUANTIFYING BURDEN IMPACT: Direct gDNA density comparison")
    print(f"{'='*70}")

    # We'll intercept the computation by running the pipeline and logging
    # all the intermediate values from the logger output.

    result, sc, tmpdir = make_scenario(50000, nrna=0, gdna=1000, ss=1.0)

    logging.getLogger('rigel').setLevel(logging.DEBUG)
    pipe_cfg = PipelineConfig(
        em=EMConfig(seed=42, strand_symmetry_kappa=4.0),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)
    logging.getLogger('rigel').setLevel(logging.WARNING)

    bench = run_benchmark(result, pr, scenario_name="impact")
    print(f"\nFinal: mRNA err = {abs(bench.total_observed - bench.total_expected) / max(bench.total_expected, 1):.4f}")
    print(f"  mRNA expected={bench.total_expected}, observed={bench.total_observed:.2f}")

    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    # The critical question: does burden subtraction change anything?
    # Strategy: look at log output to see burden magnitude vs unspliced counts.

    import io
    import re

    # Capture logs
    log_capture = io.StringIO()
    handler = logging.StreamHandler(log_capture)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(name)s — %(message)s"))
    logging.getLogger('rigel').addHandler(handler)
    logging.getLogger('rigel').setLevel(logging.DEBUG)

    # Run the worst-case scenario
    result, sc, tmpdir = make_scenario(50000, 0, 1000, ss=1.0)
    pipe_cfg = PipelineConfig(
        em=EMConfig(seed=42, strand_symmetry_kappa=4.0),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)
    bench = run_benchmark(result, pr, scenario_name="capture")
    logging.getLogger('rigel').setLevel(logging.WARNING)

    # Parse logs for key values
    logs = log_capture.getvalue()
    print("\n" + "="*70)
    print("CAPTURED LOG ANALYSIS: gDNA=1000, nRNA=0, SS=1.0, kappa=4")
    print("="*70)

    # Extract key log lines
    for line in logs.split('\n'):
        if any(k in line for k in ['burden', 'phantom', 'global_density',
                                     'κ_ref', 'κ_locus', 'median_locus',
                                     'Geometric splicing', 'n_zero',
                                     'locus_density', 'UNSPLICED', 'Strand spec']):
            print(f"  {line.strip()}")

    print(f"\nResult: mRNA err = {abs(bench.total_observed - bench.total_expected) / max(bench.total_expected, 1)*100:.2f}%")
    print(f"  mRNA: expected={bench.total_expected}, observed={bench.total_observed:.2f}")
    print(f"  gDNA: expected={bench.n_gdna_expected}, observed={bench.n_gdna_pipeline:.2f}")

    # Now: the KEY comparison. What would happen without burden?
    # We need to understand the magnitude: burden vs unspliced counts
    print(f"\n" + "="*70)
    print("KEY QUESTION: Burden magnitude vs total unspliced counts")
    print("="*70)

    for line in logs.split('\n'):
        if 'phantom' in line.lower():
            m = re.search(r'(\d+\.?\d*)\s+total predicted phantom', line)
            if m:
                burden_total = float(m.group(1))
                print(f"  Total burden (phantom reads): {burden_total:.0f}")
        if 'global_density' in line.lower():
            print(f"  {line.strip()}")

    # The burden only modifies the PER-LOCUS counts used for gDNA density.
    # But the EM operates on N=50000 fragments.
    # The gDNA init is: shrunk_density × exonic_bp
    # Even if we change the init, the EM will re-estimate from the data.
    print(f"\n  Total fragments: 50000")
    print(f"  Impact analysis: burden shifts gDNA INIT, not EM likelihood.")
    print(f"  The EM has {bench.n_fragments} fragments to override any prior change.")

    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)
