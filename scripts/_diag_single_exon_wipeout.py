#!/usr/bin/env python
"""Trace the exact logic path that causes single-exon gene wipeout at low SS.

For a single-exon gene with ss=0.65 and gdna=0, traces:
1. How fragments are classified as sense/antisense
2. The exact gDNA rate computed from the strand formula
3. EB shrinkage computation
4. EM init & competition outcome
"""
import sys
import tempfile
import logging
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.sim.scenario import Scenario
from hulkrna.sim.reads import SimConfig
from hulkrna.pipeline import run_pipeline, _compute_gdna_rate_from_strand
from hulkrna.sim.benchmark import run_benchmark

logging.basicConfig(level=logging.WARNING)
SIM_SEED = 42


def trace_single_exon(ss_value, n_frags=500):
    """Full trace for a single-exon gene at a given SS."""
    print(f"\n{'='*72}")
    print(f"SINGLE-EXON GENE, ss={ss_value}, gdna=0, nrna=0, N={n_frags}")
    print(f"{'='*72}\n")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = Scenario("se_trace", genome_length=5000, seed=SIM_SEED,
                       work_dir=tmp / "se_trace")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 700)], "abundance": 100},
        ])

        sim_config = SimConfig(
            read_length=100, frag_mean=200, frag_std=30,
            strand_specificity=ss_value, seed=SIM_SEED,
            frag_min=80, frag_max=450,
        )

        result = sc.build(n_fragments=n_frags, sim_config=sim_config)
        index = result.index

        pr = run_pipeline(
            str(result.bam_path), index,
            sj_strand_tag="ts", seed=SIM_SEED,
        )
        bench = run_benchmark(result, pr, scenario_name=f"se_ss{ss_value}")
        counter = pr.counter

        # --- Step 1: Strand counts ---
        sense_g1 = float(counter.gene_sense_all[0])
        anti_g1 = float(counter.gene_antisense_all[0])
        total_g1 = sense_g1 + anti_g1
        print("Step 1: Fragment strand classification (unique+isoform unspliced)")
        print(f"  sense    = {sense_g1:.0f}")
        print(f"  antisense= {anti_g1:.0f}")
        print(f"  total    = {total_g1:.0f}")
        if total_g1 > 0:
            print(f"  anti_frac= {anti_g1/total_g1:.3f}")
        print(f"  Expected (pure RNA): sense={total_g1*ss_value:.0f}, "
              f"anti={total_g1*(1-ss_value):.0f}")
        print()

        # --- Step 2: gDNA rate formula ---
        print("Step 2: gDNA rate formula")
        denom_ss = 2.0 * ss_value - 1.0
        print(f"  denom_ss = 2*{ss_value} - 1 = {denom_ss:.2f}")
        print(f"  Guard: denom_ss <= 0.2? {denom_ss <= 0.2}")

        if denom_ss > 0.2 and total_g1 > 0:
            g_numerator = 2.0 * (anti_g1 * ss_value - sense_g1 * (1.0 - ss_value))
            g_estimated = max(g_numerator / denom_ss, 0.0)
            raw_rate = min(g_estimated / total_g1, 1.0) if total_g1 > 0 else 0

            print(f"  Numerator = 2*({anti_g1:.0f}*{ss_value} - {sense_g1:.0f}*{1-ss_value:.2f})")
            print(f"            = 2*({anti_g1*ss_value:.1f} - {sense_g1*(1-ss_value):.1f})")
            print(f"            = {g_numerator:.1f}")
            print(f"  G_estimated = max({g_numerator:.1f}/{denom_ss:.2f}, 0) = {g_estimated:.1f}")
            print(f"  rate = {g_estimated:.1f} / {total_g1:.0f} = {raw_rate:.4f}")

            # Cross-check with actual function
            actual_rate = _compute_gdna_rate_from_strand(sense_g1, anti_g1, ss_value)
            print(f"  _compute_gdna_rate_from_strand() = {actual_rate:.4f}")
        elif denom_ss <= 0.2:
            raw_rate = 0.0
            print(f"  Guard triggered → rate = 0.0 (SS too low to disambiguate)")
        else:
            raw_rate = 0.0
            print(f"  No fragments → rate = 0.0")
        print()

        # --- Step 3: nRNA init ---
        print("Step 3: nRNA init")
        print(f"  nrna_init[t1] = {float(counter.nrna_init[0]):.1f}")
        print(f"  (Single-exon: intronic span = 0 → nrna_init = 0)")
        print()

        # --- Step 4: Pipeline output ---
        print("Step 4: Pipeline output")
        for ta in bench.transcripts:
            print(f"  {ta.t_id}: expected={ta.expected:.0f}, "
                  f"observed(mRNA)={ta.observed:.1f}, "
                  f"abs_diff={ta.abs_diff:.1f}, "
                  f"rel_err={ta.rel_error:.1%}")
        print(f"  gDNA: expected=0, pipeline={bench.n_gdna_pipeline:.1f}")
        print(f"  nRNA: expected=0, pipeline={bench.n_nrna_pipeline:.1f}")
        print()

        # --- Step 5: Count breakdown ---
        print("Step 5: Count category breakdown (t1)")
        unique = counter.unique_counts[0]
        em = counter.em_counts[0]
        total_counts = unique + em
        col_labels = [
            "unique_sense_spliced_annot", "unique_sense_spliced_unannot",
            "unique_sense_unspliced",
            "unique_anti_spliced_annot", "unique_anti_spliced_unannot",
            "unique_anti_unspliced",
            "multi_spliced_annot", "multi_spliced_unannot",
        ]
        for i in range(min(8, len(unique))):
            if unique[i] > 0 or em[i] > 0:
                label = col_labels[i] if i < len(col_labels) else f"col{i}"
                print(f"  [{i}] {label}: unique={unique[i]:.1f}, "
                      f"em={em[i]:.1f}, total={total_counts[i]:.1f}")
        print(f"  gdna_total = {counter.gdna_total:.1f}")
        print()

        sc.cleanup()
        return bench


def theoretical_gdna_rate():
    """Show what the gDNA formula predicts for pure RNA at each SS."""
    print("\n" + "="*72)
    print("THEORETICAL: gDNA rate formula on PURE RNA (no actual gDNA)")
    print("="*72)
    print("\nWhen ALL fragments are RNA, sense=N*SS, anti=N*(1-SS).")
    print("The formula should return 0, but what does it return?\n")
    print(f"  {'SS':>5}  {'denom_ss':>8}  {'guard?':>6}  "
          f"{'G_est':>7}  {'rate':>7}  {'comment'}")
    print(f"  {'-----':>5}  {'--------':>8}  {'------':>6}  "
          f"{'-------':>7}  {'-------':>7}  {'-------'}")

    N = 500
    for ss in [1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55]:
        sense = N * ss
        anti = N * (1 - ss)
        denom = 2*ss - 1
        guarded = denom <= 0.2

        if not guarded:
            g_num = 2 * (anti * ss - sense * (1-ss))
            g = max(g_num / denom, 0.0)
            rate = min(g / N, 1.0)
        else:
            g = 0
            rate = 0

        actual = _compute_gdna_rate_from_strand(sense, anti, ss)
        comment = ""
        if guarded:
            comment = "GUARDED (denom<=0.2)"
        elif abs(rate) < 1e-6:
            comment = "OK: formula returns 0"
        else:
            comment = f"FALSE POSITIVE (see below)"

        print(f"  {ss:5.2f}  {denom:8.2f}  {'YES' if guarded else 'no':>6}  "
              f"{g:7.1f}  {actual:7.4f}  {comment}")


def theoretical_explanation():
    """Explain why the formula returns 0 for pure RNA."""
    print("\n" + "="*72)
    print("WHY THE FORMULA RETURNS 0 FOR PURE RNA (theoretical)")
    print("="*72)
    print("""
For pure RNA: sense = N*SS, anti = N*(1-SS)

G = 2*(anti*SS - sense*(1-SS)) / (2*SS-1)
  = 2*(N*(1-SS)*SS - N*SS*(1-SS)) / (2*SS-1)
  = 2*(N*SS*(1-SS) - N*SS*(1-SS)) / (2*SS-1)
  = 0 / (2*SS-1)
  = 0

So theoretically the formula returns EXACTLY 0 for pure RNA.
The issue is SAMPLING VARIANCE: with finite reads, the proportion
of antisense reads won't be exactly (1-SS).

Example at SS=0.65, N=500:
  Expected: sense=325, anti=175
  Binomial std: σ ≈ √(N*p*(1-p)) ≈ √(500*0.65*0.35) ≈ 10.7
  A +1σ fluctuation on antisense: anti=186, sense=314
  G = 2*(186*0.65 - 314*0.35) / 0.30
    = 2*(120.9 - 109.9) / 0.30
    = 22.0 / 0.30
    = 73.3 fragments falsely attributed to gDNA!

The formula amplifies sampling noise by 1/denom_ss:
  - At SS=0.90: denom_ss=0.80, amplification ~1.25x
  - At SS=0.65: denom_ss=0.30, amplification ~3.33x
  - At SS=0.60: denom_ss=0.20, GUARDED (returns 0)

The std dev of the estimated gDNA count scales as:
  σ(G) ∝ 2 * √(N*SS*(1-SS)) / (2SS-1)

At SS=0.65, N=500: σ(G) ≈ 2*√(113.75)/0.30 ≈ 71.1
At SS=0.90, N=500: σ(G) ≈ 2*√(45)/0.80 ≈ 16.8

So at SS=0.65, the noise in the gDNA *estimate* is enormous relative
to a true value of 0. And because we clamp G≥0, the expected value 
of max(G,0) is NOT zero — it's biased positive (half-normal).
    E[max(G,0)] ≈ σ(G) * √(2/π) ≈ 71.1 * 0.798 ≈ 56.7
""")


def trace_em_internals(ss_value=0.9, n_frags=500):
    """Inspect the actual EM data to find the asymmetry."""
    print(f"\n{'='*72}")
    print(f"EM INTERNALS TRACE: ss={ss_value}, N={n_frags}")
    print(f"{'='*72}\n")

    from hulkrna.pipeline import (
        scan_and_buffer,
        count_from_buffer,
        _scan_and_build_em_data,
        _build_loci,
        _build_locus_em_data,
        _compute_eb_gdna_priors,
        _compute_nrna_init,
        _GDNA_SPLICE_PENALTIES,
    )
    from hulkrna.categories import SpliceType
    import pysam

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = Scenario("se_em", genome_length=5000, seed=SIM_SEED,
                       work_dir=tmp / "se_em")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 700)], "abundance": 100},
        ])
        sim_config = SimConfig(
            read_length=100, frag_mean=200, frag_std=30,
            strand_specificity=ss_value, seed=SIM_SEED,
            frag_min=80, frag_max=450,
        )
        result = sc.build(n_fragments=n_frags, sim_config=sim_config)
        index = result.index

        # Run scan phase
        bamfh = pysam.AlignmentFile(str(result.bam_path), "rb")
        stats, strand_models, insert_models, buffer = scan_and_buffer(
            bamfh.fetch(until_eof=True), index,
            sj_strand_tag="ts",
        )
        bamfh.close()

        print(f"Strand models:")
        print(f"  exonic_spliced: n_obs={strand_models.exonic_spliced.n_observations}, "
              f"SS={strand_models.exonic_spliced.strand_specificity:.3f}")
        print(f"  intergenic: n_obs={strand_models.intergenic.n_observations}, "
              f"SS={strand_models.intergenic.strand_specificity:.3f}")
        print(f"Insert models:")
        print(f"  global: n_obs={insert_models.global_model.n_observations}, "
              f"mean={insert_models.global_model.mean:.1f}")
        print()

        # Compute geometry
        exonic_lengths = index.t_df["length"].values.astype(np.float64)
        mean_frag = (
            insert_models.global_model.mean
            if insert_models.global_model.n_observations > 0 else 200.0
        )
        effective_lengths = np.maximum(exonic_lengths - mean_frag + 1.0, 1.0)
        transcript_spans = (
            index.t_df["end"].values - index.t_df["start"].values
        ).astype(np.float64)
        intronic_spans = index.intron_span_per_gene()

        from hulkrna.counter import ReadCounter
        counter = ReadCounter(
            index.num_transcripts, index.num_genes, seed=SIM_SEED,
            alpha=0.5, effective_lengths=effective_lengths,
            t_to_g=index.t_to_g_arr,
            gene_spans=(index.g_df["end"].values - index.g_df["start"].values).astype(np.float64),
            mean_frag=mean_frag,
            intronic_spans=intronic_spans,
            transcript_spans=transcript_spans,
        )

        # Build EM data
        gdna_splice_penalties = dict(_GDNA_SPLICE_PENALTIES)

        em_data = _scan_and_build_em_data(
            buffer, index, strand_models, insert_models, counter, stats,
            log_every=1_000_000,
            gdna_splice_penalties=gdna_splice_penalties,
        )

        print(f"EM data: {em_data.n_units} units, {em_data.n_candidates} candidates")
        print(f"  Spliced: {em_data.is_spliced.sum()}")
        print(f"  Unspliced: {(~em_data.is_spliced).sum()}")
        print()

        # Check gdna log-likelihoods
        gdna_lls = em_data.gdna_log_liks
        valid_gdna = gdna_lls[gdna_lls > -np.inf]
        print(f"gDNA log-likelihoods:")
        print(f"  Valid: {len(valid_gdna)}/{len(gdna_lls)}")
        if len(valid_gdna) > 0:
            print(f"  Mean: {valid_gdna.mean():.4f}")
            print(f"  Min: {valid_gdna.min():.4f}, Max: {valid_gdna.max():.4f}")
        print()

        # Check RNA log-likelihoods
        all_lls = em_data.log_liks
        print(f"RNA log-likelihoods:")
        print(f"  Total candidates: {len(all_lls)}")
        print(f"  Mean: {all_lls.mean():.4f}")
        print(f"  Min: {all_lls.min():.4f}, Max: {all_lls.max():.4f}")
        # Histogram
        import collections
        bins = collections.Counter()
        for ll in all_lls:
            bins[round(ll, 1)] += 1
        print(f"  Histogram (rounded):")
        for k in sorted(bins.keys()):
            print(f"    {k:.1f}: {bins[k]}")
        print()

        # Compare RNA vs gDNA likelihoods per unit
        offsets = em_data.offsets
        diffs = []
        for u in range(em_data.n_units):
            s, e = int(offsets[u]), int(offsets[u+1])
            rna_lls = all_lls[s:e]
            gdna_ll = float(gdna_lls[u])
            if gdna_ll > -np.inf and len(rna_lls) > 0:
                best_rna = float(rna_lls.max())
                diffs.append(best_rna - gdna_ll)

        diffs = np.array(diffs)
        print(f"Per-unit: best_RNA_ll - gDNA_ll:")
        print(f"  N={len(diffs)}")
        if len(diffs) > 0:
            print(f"  Mean: {diffs.mean():.4f}")
            print(f"  Min: {diffs.min():.4f}, Max: {diffs.max():.4f}")
            print(f"  >0 (RNA favored): {(diffs > 0).sum()}")
            print(f"  =0 (tied): {(diffs == 0).sum()}")
            print(f"  <0 (gDNA favored): {(diffs < 0).sum()}")
        print()

        # Build loci
        loci = _build_loci(em_data, index)
        print(f"Loci: {len(loci)}")
        for locus in loci:
            if 0 in locus.gene_indices:
                print(f"  g1 locus: {len(locus.unit_indices)} units, "
                      f"footprint={locus.gdna_footprint_bp}bp")
                print(f"  mRNA eff_len = {effective_lengths[0]:.1f}")
                gdna_eff = max(float(locus.gdna_footprint_bp) - mean_frag + 1, 1)
                print(f"  gDNA eff_len = {gdna_eff:.1f}")
                print(f"  nRNA eff_len = {max(transcript_spans[0] - mean_frag + 1, 1):.1f}")
                print()

                # Compute gDNA init
                eb = _compute_eb_gdna_priors(
                    loci, em_data, counter, index, strand_models,
                )
                gdna_init = eb[0]
                print(f"  gdna_init = {gdna_init:.4f}")

                # Build locus EM
                locus_em = _build_locus_em_data(
                    locus, em_data, counter, index,
                    mean_frag, gdna_init,
                )
                print(f"  EM components: {locus_em.n_components}")
                print(f"  unique_totals = {locus_em.unique_totals}")
                print(f"  prior = {locus_em.prior}")
                print(f"  eff_len = {locus_em.effective_lengths}")
                print()

                # Run EM
                theta, alpha = counter.run_locus_em(locus_em)
                print(f"  Converged theta = {theta}")
                print(f"  Converged alpha = {alpha}")
                n_t = locus_em.n_transcripts
                print(f"  mRNA total = {alpha[:n_t].sum():.1f}")
                print(f"  nRNA total = {alpha[n_t:2*n_t].sum():.1f}")
                print(f"  gDNA total = {alpha[2*n_t]:.1f}")
                break

        buffer.cleanup()
        sc.cleanup()


if __name__ == "__main__":
    # Detailed EM internals trace
    trace_em_internals(0.9)

    # Theoretical analysis
    theoretical_gdna_rate()
    theoretical_explanation()

    # Concrete traces
    trace_single_exon(0.9)
    trace_single_exon(0.65)
