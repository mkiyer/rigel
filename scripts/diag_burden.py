#!/usr/bin/env python3
"""Diagnostic: trace the geometric splicing burden subtraction end-to-end.

Runs a single simulation with known parameters and instruments
compute_eb_gdna_priors to see exactly what the burden values are.
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
from rigel.estimator import compute_unspliced_to_spliced_ratios

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)-5s %(name)s — %(message)s")
logger = logging.getLogger(__name__)


def make_scenario(n_fragments, nrna_abundance, gdna_abundance, ss=1.0, seed=42):
    """Create a standard test scenario."""
    with tempfile.TemporaryDirectory(prefix="diag_") as tmpdir:
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

        yield result, sc


def run_diagnostic(n_fragments, nrna, gdna, ss=1.0, kappa=4.0):
    """Run pipeline and extract diagnostic info about burden subtraction."""
    print(f"\n{'='*70}")
    print(f"DIAGNOSTIC: n={n_fragments}, nRNA={nrna}, gDNA={gdna}, SS={ss}, kappa={kappa}")
    print(f"{'='*70}")

    for result, sc in make_scenario(n_fragments, nrna, gdna, ss):
        pipe_cfg = PipelineConfig(
            em=EMConfig(seed=42, strand_symmetry_kappa=kappa),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )

        # Monkey-patch compute_eb_gdna_priors to capture intermediate values
        from rigel import locus as locus_mod
        original_fn = locus_mod.compute_eb_gdna_priors
        captured = {}

        def patched_compute_eb_gdna_priors(*args, **kwargs):
            # Capture the ratios argument
            ratios = kwargs.get('unspliced_to_spliced_ratios')
            captured['ratios'] = ratios

            # Instrument _compute_locus_mrna_burden
            original_burden = locus_mod._compute_locus_mrna_burden
            burden_results = []

            def instrumented_burden(locus, em_data, estimator, ratios, strand_spec):
                bs, ba = original_burden(locus, em_data, estimator, ratios, strand_spec)
                # Also capture the inputs
                t_arr = locus.transcript_indices
                n_t = estimator.num_transcripts
                spliced_counts = {}
                for t_idx in t_arr:
                    c_spliced = estimator.unambig_counts[int(t_idx), 2:6].sum()
                    r_t = ratios[int(t_idx)]
                    spliced_counts[int(t_idx)] = {
                        'spliced_unambig': c_spliced,
                        'r_t': r_t,
                        'burden_contribution': c_spliced * r_t if np.isfinite(r_t) and c_spliced > 0 else 0.0,
                    }
                burden_results.append({
                    'locus_id': locus.locus_id,
                    'transcript_indices': list(t_arr),
                    'burden_sense': bs,
                    'burden_anti': ba,
                    'per_transcript': spliced_counts,
                })
                return bs, ba

            locus_mod._compute_locus_mrna_burden = instrumented_burden
            result = original_fn(*args, **kwargs)
            locus_mod._compute_locus_mrna_burden = original_burden
            captured['burden_results'] = burden_results
            captured['gdna_inits'] = result
            return result

        locus_mod.compute_eb_gdna_priors = patched_compute_eb_gdna_priors

        # Also capture the unspliced sense/anti arrays
        original_per_locus = locus_mod._compute_per_locus_gdna_densities
        per_locus_inputs = {}

        def patched_per_locus(loci, t_sense, t_anti, t_exonic, t_refs,
                             strand_spec, intergenic_density, ref_shrunk,
                             global_density, locus_burden_sense=None,
                             locus_burden_anti=None):
            per_locus_inputs['t_sense'] = t_sense.copy()
            per_locus_inputs['t_anti'] = t_anti.copy()
            per_locus_inputs['t_exonic'] = t_exonic.copy()
            per_locus_inputs['burden_sense'] = locus_burden_sense
            per_locus_inputs['burden_anti'] = locus_burden_anti
            return original_per_locus(loci, t_sense, t_anti, t_exonic, t_refs,
                                      strand_spec, intergenic_density, ref_shrunk,
                                      global_density, locus_burden_sense,
                                      locus_burden_anti)

        locus_mod._compute_per_locus_gdna_densities = patched_per_locus

        pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)

        # Restore
        locus_mod.compute_eb_gdna_priors = original_fn
        locus_mod._compute_per_locus_gdna_densities = original_per_locus

        bench = run_benchmark(result, pr, scenario_name="diag")

        # Print diagnostics
        print(f"\n--- Benchmark Summary ---")
        print(bench.summary())

        print(f"\n--- Geometric Splicing Ratios (R_t) ---")
        ratios = captured.get('ratios')
        if ratios is not None:
            for i, r in enumerate(ratios):
                print(f"  t[{i}]: R_t = {r:.4f}" if np.isfinite(r) else f"  t[{i}]: R_t = inf (single-exon)")
        else:
            print("  RATIOS WERE None — burden subtraction was DISABLED!")

        print(f"\n--- Per-Locus Burden Subtraction Details ---")
        burden_results = captured.get('burden_results', [])
        if not burden_results:
            print("  NO BURDEN RESULTS — _compute_locus_mrna_burden was never called!")
        for b in burden_results:
            print(f"\n  Locus {b['locus_id']} (transcripts: {b['transcript_indices']}):")
            print(f"    burden_sense = {b['burden_sense']:.2f}")
            print(f"    burden_anti  = {b['burden_anti']:.2f}")
            for t_idx, info in b['per_transcript'].items():
                print(f"    t[{t_idx}]: spliced_unambig={info['spliced_unambig']:.0f}, "
                      f"R_t={info['r_t']:.4f}" + (f", burden={info['burden_contribution']:.2f}"
                                                    if info['burden_contribution'] > 0 else ""))

        print(f"\n--- Unspliced Sense/Anti Arrays Going Into gDNA Density ---")
        t_sense = per_locus_inputs.get('t_sense')
        t_anti = per_locus_inputs.get('t_anti')
        if t_sense is not None:
            for i in range(len(t_sense)):
                print(f"  t[{i}]: unspliced_sense={t_sense[i]:.1f}, unspliced_anti={t_anti[i]:.1f}")
            print(f"  TOTAL: sense={t_sense.sum():.1f}, anti={t_anti.sum():.1f}")

        burden_s = per_locus_inputs.get('burden_sense')
        burden_a = per_locus_inputs.get('burden_anti')
        if burden_s is not None:
            print(f"\n--- Burden Values Passed to _compute_per_locus_gdna_densities ---")
            for i, (bs, ba) in enumerate(zip(burden_s, burden_a)):
                print(f"  locus[{i}]: burden_sense={bs:.2f}, burden_anti={ba:.2f}")
        else:
            print("\n  BURDEN was None at _compute_per_locus_gdna_densities!")

        print(f"\n--- gDNA Inits ---")
        gdna_inits = captured.get('gdna_inits', [])
        for i, gi in enumerate(gdna_inits):
            print(f"  locus[{i}]: gdna_init = {gi:.4f}")

        return bench, captured, per_locus_inputs


if __name__ == "__main__":
    # Test case 1: Moderate gDNA, no nRNA — the gDNA siphon scenario
    print("\n" + "#"*70)
    print("# TEST 1: gDNA=100, nRNA=0, SS=1.0, kappa=4")
    print("#"*70)
    run_diagnostic(50000, nrna=0, gdna=100, ss=1.0, kappa=4.0)

    # Test case 2: Heavy gDNA
    print("\n" + "#"*70)
    print("# TEST 2: gDNA=1000, nRNA=0, SS=1.0, kappa=4")
    print("#"*70)
    run_diagnostic(50000, nrna=0, gdna=1000, ss=1.0, kappa=4.0)

    # Test case 3: With nRNA
    print("\n" + "#"*70)
    print("# TEST 3: gDNA=100, nRNA=100, SS=0.95, kappa=4")
    print("#"*70)
    run_diagnostic(50000, nrna=100, gdna=100, ss=0.95, kappa=4.0)
