#!/usr/bin/env python
"""Concrete demonstration of filter_by_overlap bug in antisense scenario.

Shows exactly HOW and WHY filter_by_overlap removes the correct gene from
the candidate set for unspliced (nRNA) fragments in the antisense overlap
zone.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.sim.scenario import Scenario
from hulkrna.sim.reads import SimConfig
from hulkrna.resolution import filter_by_overlap
from hulkrna.pipeline import run_pipeline
from hulkrna.sim.benchmark import run_benchmark

SIM_SEED = 42


def print_geometry():
    """Print the antisense overlap geometry."""
    print("=" * 72)
    print("ANTISENSE GENE GEOMETRY")
    print("=" * 72)
    print()
    print("g1 (+strand): t1 with exons at (200,500) and (1000,1300)")
    print("              intron: 500 -> 1000")
    print()
    print("g2 (-strand): t2 with exons at (300,600) and (1100,1400)")
    print("              intron: 600 -> 1100")
    print()
    print("Overlap zones (genomic coordinate):")
    print("  200-300 : g1 exon1 only")
    print("  300-500 : g1 exon1 + g2 exon1  [EXON-EXON overlap]")
    print("  500-600 : g1 INTRON + g2 exon1 [asymmetric!]")
    print("  600-1000: g1 intron + g2 intron [INTRON-INTRON]")
    print("  1000-1100: g1 EXON2 + g2 INTRON [asymmetric! KEY ZONE]")
    print("  1100-1300: g1 exon2 + g2 exon2  [EXON-EXON overlap]")
    print("  1300-1400: g2 exon2 only")
    print()
    print("The zone 1000-1100 is EXON for g1 but INTRON for g2.")
    print("An nRNA fragment from g2 spanning this boundary will have")
    print("MORE exon_bp for g1 than for g2 (its true source gene).")
    print()


def demonstrate_concrete_fragments():
    """Show the bug with hand-traced fragments."""
    print("=" * 72)
    print("CONCRETE EXAMPLES: Fragments in the asymmetric overlap zone")
    print("=" * 72)

    cases = [
        ("nRNA from g2 spanning intron->exon2 boundary",
         800, 1200, "g2 pre-mRNA, mostly intronic"),
        ("nRNA from g2 centered on asymmetric zone",
         950, 1150, "g2 pre-mRNA, straddles 1000-1100"),
        ("nRNA from g1 spanning exon2->right boundary",
         1050, 1250, "g1 pre-mRNA, mostly exonic for g1"),
    ]

    for title, frag_start, frag_end, desc in cases:
        frag_length = frag_end - frag_start
        print()
        print(f"--- {title} ---")
        print(f"Fragment: {frag_start}-{frag_end} ({frag_length}bp) [{desc}]")
        print()

        # t1 (g1+): exons (200,500), (1000,1300); intron (500,1000)
        t1_exon_bp = (max(0, min(500, frag_end) - max(200, frag_start))
                      + max(0, min(1300, frag_end) - max(1000, frag_start)))
        t1_intron_bp = max(0, min(1000, frag_end) - max(500, frag_start))

        # t2 (g2-): exons (300,600), (1100,1400); intron (600,1100)
        t2_exon_bp = (max(0, min(600, frag_end) - max(300, frag_start))
                      + max(0, min(1400, frag_end) - max(1100, frag_start)))
        t2_intron_bp = max(0, min(1100, frag_end) - max(600, frag_start))

        t1_genic = t1_exon_bp + t1_intron_bp
        t2_genic = t2_exon_bp + t2_intron_bp

        print(f"  t1 (g1+): exon_bp={t1_exon_bp:3d}, intron_bp={t1_intron_bp:3d}  "
              f"exon_frac={t1_exon_bp/frag_length:.3f}  "
              f"genic_frac={t1_genic/frag_length:.3f}")
        print(f"  t2 (g2-): exon_bp={t2_exon_bp:3d}, intron_bp={t2_intron_bp:3d}  "
              f"exon_frac={t2_exon_bp/frag_length:.3f}  "
              f"genic_frac={t2_genic/frag_length:.3f}")
        print()

        overlap_bp = {0: (t1_exon_bp, t1_intron_bp),
                      1: (t2_exon_bp, t2_intron_bp)}
        t_inds = frozenset({0, 1})

        best = max(t1_exon_bp, t2_exon_bp) / frag_length
        threshold = best * 0.99
        t1_frac = t1_exon_bp / frag_length
        t2_frac = t2_exon_bp / frag_length

        winner = "t1" if t1_exon_bp >= t2_exon_bp else "t2"
        loser = "t2" if winner == "t1" else "t1"
        loser_frac = t2_frac if winner == "t1" else t1_frac

        result = filter_by_overlap(t_inds, overlap_bp, frag_length,
                                   min_frac_of_best=0.99)

        print(f"  filter_by_overlap(min_frac_of_best=0.99):")
        print(f"    best exon_frac = {best:.3f} ({winner})")
        print(f"    threshold      = {threshold:.3f}")
        print(f"    {loser} exon_frac = {loser_frac:.3f} "
              f"{'< threshold -> REMOVED' if loser_frac < threshold else '>= threshold -> KEPT'}")
        print(f"    Result: {set(result)}")

        removed = t_inds - result
        if removed:
            # Which was removed?
            rem_idx = list(removed)[0]
            rem_name = "t2 (g2)" if rem_idx == 1 else "t1 (g1)"
            rem_genic = t2_genic if rem_idx == 1 else t1_genic
            print()
            print(f"    BUG: {rem_name} removed despite having "
                  f"genic_frac={rem_genic/frag_length:.3f}.")
            print(f"         The filter only sees exon_bp, ignoring intron_bp.")
            print(f"         For nRNA fragments, intron overlap IS the signal.")

    print()
    print("=" * 72)
    print("ROOT CAUSE")
    print("=" * 72)
    print()
    print("filter_by_overlap() uses EXON overlap only to decide which")
    print("transcripts to keep.  For UNSPLICED fragments (nRNA/gDNA),")
    print("most of the overlap is INTRONIC.  In an antisense overlap,")
    print("gene A's exon spatially overlaps gene B's intron.  So an")
    print("nRNA fragment from gene B in the overlap zone will have:")
    print()
    print("  - HIGH exon overlap with gene A (false signal)")
    print("  - LOW exon overlap with gene B  (true source!)")
    print()
    print("The filter removes gene B (the correct gene), leaving only")
    print("gene A.  The EM then has no choice but to assign the fragment")
    print("to gene A, inflating its count and zeroing gene B's nRNA.")
    print()


def demonstrate_pipeline_impact():
    """Run simulation with nrna=50% and compare filter ON vs OFF."""
    print("=" * 72)
    print("PIPELINE IMPACT: antisense scenario, nrna=50%, ss=1.0, gdna=0")
    print("=" * 72)
    print()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = Scenario("antisense_diag", genome_length=8000, seed=SIM_SEED,
                       work_dir=tmp / "antisense_diag")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(300, 600), (1100, 1400)],
             "abundance": 100},
        ])
        sc.add_gene("g_ctrl", "+", [
            {"t_id": "t_ctrl", "exons": [(5500, 5800)], "abundance": 0},
        ])

        sim_config = SimConfig(
            read_length=100,
            frag_mean=300,
            frag_std=30,
            strand_specificity=1.0,
            seed=SIM_SEED,
        )

        result = sc.build(
            n_fragments=1000,
            sim_config=sim_config,
            nrna_fraction=0.50,
        )
        index = result.index
        bam_path = result.bam_path

        print(f"Simulated fragments: {result.n_simulated}")
        print(f"  nRNA ground truth: {result.ground_truth_nrna_count()}")
        print()

        # Run pipeline with filter OFF
        pr_off = run_pipeline(
            str(bam_path), index,
            sj_strand_tag="ts", seed=SIM_SEED,
            overlap_min_frac=0.0,
        )
        bench_off = run_benchmark(result, pr_off, scenario_name="filter_OFF")

        # Run pipeline with filter ON (default)
        pr_on = run_pipeline(
            str(bam_path), index,
            sj_strand_tag="ts", seed=SIM_SEED,
            overlap_min_frac=0.99,
        )
        bench_on = run_benchmark(result, pr_on, scenario_name="filter_ON")

        print(f"{'Transcript':<12} {'Expected':>10} {'Filter OFF':>12} "
              f"{'Filter ON':>12} {'Delta':>8}")
        print("-" * 58)
        for ta_off, ta_on in zip(bench_off.transcripts, bench_on.transcripts):
            delta = ta_on.observed - ta_off.observed
            marker = " <-- LOSS" if delta < -5 else ""
            print(f"{ta_off.t_id:<12} {ta_off.expected:>10.0f} "
                  f"{ta_off.observed:>12.1f} {ta_on.observed:>12.1f} "
                  f"{delta:>+8.1f}{marker}")

        print(f"{'gDNA':<12} {'0':>10} {bench_off.n_gdna_pipeline:>12.1f} "
              f"{bench_on.n_gdna_pipeline:>12.1f} "
              f"{bench_on.n_gdna_pipeline - bench_off.n_gdna_pipeline:>+8.1f}")
        print(f"{'nRNA (total)':<12} "
              f"{result.ground_truth_nrna_count():>10} "
              f"{bench_off.n_nrna_pipeline:>12.1f} "
              f"{bench_on.n_nrna_pipeline:>12.1f} "
              f"{bench_on.n_nrna_pipeline - bench_off.n_nrna_pipeline:>+8.1f}")
        print()

        total_off = sum(ta.observed for ta in bench_off.transcripts)
        total_on = sum(ta.observed for ta in bench_on.transcripts)
        total_exp = sum(ta.expected for ta in bench_off.transcripts
                        if ta.expected > 0)
        print(f"Total RNA observed (filter OFF): {total_off:.1f}")
        print(f"Total RNA observed (filter ON):  {total_on:.1f}")
        print(f"Total RNA expected:              {total_exp:.0f}")
        print(f"RNA loss from filter:            {total_off - total_on:.1f} fragments")
        print()

        sc.cleanup()


if __name__ == "__main__":
    print_geometry()
    demonstrate_concrete_fragments()
    demonstrate_pipeline_impact()
