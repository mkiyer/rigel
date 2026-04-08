"""Trace per-fragment log-likelihoods to understand gDNA dominance.

Runs the nRNA scenario and dumps:
- Per-unit candidates (component indices + log-likelihoods)
- gDNA vs nRNA log-likelihood comparison
- Unit breakdown: gDNA-only vs RNA-only vs both
"""
import numpy as np
from pathlib import Path
import tempfile
import logging

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import Scenario, SimConfig
from rigel.locus import build_loci, compute_locus_priors
from rigel.scoring import FragmentScorer
from rigel.scan import FragmentRouter
from rigel.estimator import AbundanceEstimator
from rigel.frag_length_model import FragmentLengthModels
from rigel.calibration import calibrate_gdna

logging.basicConfig(level=logging.WARNING)

SEED = 42
N_FRAGMENTS = 500


def run_trace(nrna_abundance, strand_specificity):
    """Trace EM internals for a single scenario."""
    print(f"\n{'='*70}")
    print(f"nRNA={nrna_abundance}, SS={strand_specificity}")
    print(f"{'='*70}")

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        sc = Scenario(
            "trace", genome_length=20000, seed=SEED,
            work_dir=tmp_path / "trace",
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

        # Run scan phase
        (stats, strand_models, fl_models, buffer,
         region_counts, fl_table) = scan_and_buffer(
            str(result.bam_path), result.index, config.scan,
        )

        # Finalize models (required before scoring!)
        strand_models.finalize()
        fl_models.build_scoring_models()
        fl_models.finalize()

        # Calibrate
        calibration = calibrate_gdna(
            region_counts, fl_table,
            result.index.region_df,
            strand_models.strand_specificity,
            density_percentile=config.calibration.density_percentile,
        )

        print(f"\nCalibration: E[gDNA]={calibration.region_e_gdna.sum():.1f}, "
              f"N_total={calibration.region_n_total.sum():.1f}")
        print(f"SS={calibration.strand_specificity:.4f}")

        # Build scored fragments via the scoring pipeline
        estimator = AbundanceEstimator(
            result.index.num_transcripts, em_config=config.em,
        )
        scorer = FragmentScorer.from_models(
            strand_models, fl_models, result.index, estimator,
        )
        router = FragmentRouter(
            scorer, estimator, stats, result.index, strand_models,
        )
        em_data = router.scan(buffer, log_every=999999)
        buffer.release()

        # Build loci and compute priors
        loci = build_loci(em_data, result.index)
        alpha_gdna, alpha_rna = compute_locus_priors(
            loci, result.index, calibration,
            kappa_min=config.calibration.gdna_prior_kappa_min,
            kappa_max=config.calibration.gdna_prior_kappa_max,
        )

        # Access the raw scored fragment arrays
        offsets = em_data.offsets
        t_indices = em_data.t_indices
        log_liks = em_data.log_liks
        gdna_log_liks = em_data.gdna_log_liks
        is_spliced = em_data.is_spliced
        tx_starts = em_data.tx_starts
        tx_ends = em_data.tx_ends
        genomic_footprints = em_data.genomic_footprints
        n_total_units = len(offsets) - 1

        # Dump FL model info
        print(f"\n  FL model info:")
        print(f"    has RNA FL:  {fl_models.rna_model.total_weight > 0}")
        print(f"    has gDNA FL: {fl_models.gdna_model.total_weight > 0}")
        if fl_models.rna_model.total_weight > 0:
            lut = fl_models.rna_model._log_prob
            print(f"    RNA FL LUT size: {len(lut)}")
            print(f"    RNA FL mean: {fl_models.rna_model.mean:.1f}")
            # Print a few values around 200bp
            for flen in [0, 100, 150, 200, 250, 300, 400]:
                if flen < len(lut):
                    print(f"    RNA FL({flen}) = {lut[flen]:.4f}")
        else:
            print(f"    RNA FL model is EMPTY (no training data)")
        if fl_models.gdna_model.total_weight > 0:
            lut = fl_models.gdna_model._log_prob
            print(f"    gDNA FL LUT size: {len(lut)}")
            for flen in [0, 100, 150, 200, 250, 300, 400]:
                if flen < len(lut):
                    print(f"    gDNA FL({flen}) = {lut[flen]:.4f}")
        else:
            print(f"    gDNA FL model is EMPTY (no training data)")

        # Dump strand model
        print(f"\n  Strand model: SS={strand_models.strand_specificity:.4f}")
        print(f"    log_p_sense={scorer.log_p_sense:.4f}")
        print(f"    log_p_antisense={scorer.log_p_antisense:.4f}")
        print(f"    r1_antisense={scorer.r1_antisense}")

        for li, locus in enumerate(loci):
            if len(locus.transcript_indices) == 0:
                continue

            # Get units for this locus
            unit_ids = locus.unit_indices
            n_units = len(unit_ids)

            if n_units < 10:
                continue

            print(f"\n--- Locus {li}: {n_units} units, "
                  f"{len(locus.transcript_indices)} transcripts ---")
            print(f"  alpha_gdna={alpha_gdna[li]:.6f}, alpha_rna={alpha_rna[li]:.6f}")
            gamma = alpha_gdna[li] / (alpha_gdna[li] + alpha_rna[li] + 1e-30)
            print(f"  gamma={gamma:.6f}")

            # Transcript info
            t_df = result.index.t_df
            for ti in locus.transcript_indices:
                row = t_df.iloc[ti]
                print(f"  transcript {ti}: {row.get('transcript_id', '?')}, "
                      f"length={row.get('length', '?')}")

            # Count breakdown
            n_spliced_count = sum(1 for ui in unit_ids if is_spliced[ui])
            n_unspliced_count = n_units - n_spliced_count
            print(f"  spliced={n_spliced_count}, unspliced={n_unspliced_count}")

            # Log-likelihood statistics per candidate type
            rna_lls = []
            gdna_lls_both = []
            n_gdna_only = 0
            n_rna_only = 0
            n_both = 0
            n_gdna_wins = 0
            n_rna_wins = 0

            for ui in unit_ids:
                start = offsets[ui]
                end = offsets[ui + 1]
                n_cands = end - start
                spliced = bool(is_spliced[ui])
                gdna_ll = float(gdna_log_liks[ui])
                has_gdna = not spliced and np.isfinite(gdna_ll)

                unit_rna_lls = [float(log_liks[j]) for j in range(start, end)]

                if has_gdna and n_cands > 0:
                    n_both += 1
                    best_rna = max(unit_rna_lls)
                    rna_lls.append(best_rna)
                    gdna_lls_both.append(gdna_ll)
                    if gdna_ll > best_rna:
                        n_gdna_wins += 1
                    else:
                        n_rna_wins += 1
                elif has_gdna and n_cands == 0:
                    n_gdna_only += 1
                    gdna_lls_both.append(gdna_ll)
                elif n_cands > 0:
                    n_rna_only += 1

            print(f"\n  Unit breakdown:")
            print(f"    RNA-only units:       {n_rna_only}")
            print(f"    gDNA-only units:      {n_gdna_only}")
            print(f"    Both RNA+gDNA units:  {n_both}")

            if n_both > 0:
                rna_arr = np.array(rna_lls)
                gdna_arr = np.array(gdna_lls_both[:n_both])
                diff = gdna_arr - rna_arr  # positive = gDNA wins
                print(f"\n  Per-fragment log-lik comparison (gDNA vs best RNA) [{n_both} units]:")
                print(f"    gDNA wins: {n_gdna_wins} ({100*n_gdna_wins/n_both:.1f}%)")
                print(f"    RNA wins:  {n_rna_wins} ({100*n_rna_wins/n_both:.1f}%)")
                print(f"    Mean gDNA-RNA diff: {diff.mean():.3f}")
                print(f"    Median diff:        {np.median(diff):.3f}")
                print(f"    Min diff:           {diff.min():.3f}")
                print(f"    Max diff:           {diff.max():.3f}")
                print(f"    Std diff:           {diff.std():.3f}")
                print(f"\n    Mean RNA ll:  {rna_arr.mean():.3f}")
                print(f"    Mean gDNA ll: {gdna_arr.mean():.3f}")

            if n_gdna_only > 0:
                print(f"\n    *** {n_gdna_only} units with gDNA as ONLY candidate! ***")
                print(f"    These fragments have NO RNA candidate to compete with gDNA.")

            # Show first few example units with full detail
            print(f"\n  First 10 units (detail):")
            for idx, ui in enumerate(unit_ids[:10]):
                start = offsets[ui]
                end = offsets[ui + 1]
                n_cands = end - start
                spliced = bool(is_spliced[ui])
                gdna_ll = float(gdna_log_liks[ui])
                has_gdna = not spliced and np.isfinite(gdna_ll)
                gfp = int(genomic_footprints[ui])

                cands = []
                for j in range(start, end):
                    t_idx = int(t_indices[j])
                    ll = float(log_liks[j])
                    ts = int(tx_starts[j])
                    te = int(tx_ends[j])
                    flen = te - ts
                    cands.append(f"t{t_idx}(ll={ll:.2f}, flen={flen}, "
                                 f"ts={ts}, te={te})")

                gdna_str = f"gDNA={gdna_ll:.2f}" if has_gdna else "no-gDNA"
                print(f"    unit {idx}: [{'; '.join(cands)}] {gdna_str} "
                      f"gfp={gfp} {'spliced' if spliced else 'unspliced'}")

        sc.cleanup()


if __name__ == "__main__":
    # Run with same params as test_nrna_double_counting
    run_trace(70, 0.65)

    # Also run full pipeline for the failing case
    print("\n\n" + "=" * 70)
    print("FULL PIPELINE RUN: nRNA=70, SS=0.65, gdna=0")
    print("=" * 70)

    import tempfile
    from pathlib import Path
    from rigel.config import EMConfig, PipelineConfig, BamScanConfig
    from rigel.pipeline import run_pipeline
    from rigel.sim import Scenario, SimConfig, run_benchmark

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        sc = Scenario(
            "fullpipe", genome_length=20000, seed=SEED,
            work_dir=tmp_path / "fullpipe",
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
            read_length=100, strand_specificity=0.65, seed=SEED,
        )
        result = sc.build_oracle(
            n_fragments=N_FRAGMENTS, sim_config=sim_config,
            gdna_config=None, nrna_abundance=70,
        )
        print(f"Simulated fragments: {result.n_simulated}")

        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )

        pr = run_pipeline(result.bam_path, result.index, config=config)

        # Dump estimator results
        est = pr.estimator
        print(f"\nEstimator results:")
        print(f"  em_counts shape: {est.em_counts.shape}")
        em_total = est.em_counts.sum(axis=1)
        for i in range(len(em_total)):
            if em_total[i] > 0.1:
                row = result.index.t_df.iloc[i]
                tid = row.get("transcript_id", f"t{i}")
                is_nrna = row.get("is_nrna", False)
                is_syn = row.get("is_synthetic", False)
                length = row.get("length", 0)
                print(f"  [{i}] {tid}: em={em_total[i]:.1f}, "
                      f"len={length}, nrna={is_nrna}, syn={is_syn}")

        unambig = est.unambig_counts.sum(axis=1)
        for i in range(len(unambig)):
            if unambig[i] > 0.1:
                row = result.index.t_df.iloc[i]
                tid = row.get("transcript_id", f"t{i}")
                print(f"  [{i}] {tid}: unambig={unambig[i]:.1f}")

        print(f"\n  gdna_em_count: {est.gdna_em_count:.1f}")
        print(f"  nrna_em_count: {est.nrna_em_count:.1f}")
        print(f"  intergenic: {pr.stats.n_intergenic}")

        # Run benchmark to see the same metrics as the test
        bench = run_benchmark(result, pr, scenario_name="g0_n70_s65")
        print(f"\nBenchmark:")
        print(f"  total_expected: {bench.total_expected}")
        print(f"  total_observed: {bench.total_observed:.1f}")
        print(f"  n_nrna_expected: {bench.n_nrna_expected}")
        print(f"  n_nrna_pipeline: {bench.n_nrna_pipeline:.1f}")
        print(f"  n_synthetic_observed: {bench.n_synthetic_observed:.1f}")
        print(f"  total_rna_observed: {bench.total_rna_observed:.1f}")
        print(f"  n_gdna_pipeline: {bench.n_gdna_pipeline:.1f}")
        print(f"  n_gdna_expected: {bench.n_gdna_expected}")
        total_rna_exp = bench.total_expected + bench.n_nrna_expected
        print(f"  total_rna_expected: {total_rna_exp}")
        if total_rna_exp > 0:
            print(f"  rna_rel_err: "
                  f"{abs(bench.total_rna_observed - total_rna_exp) / total_rna_exp:.4f}")

        sc.cleanup()
