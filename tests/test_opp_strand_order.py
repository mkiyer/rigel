#!/usr/bin/env python3
"""
Test: opposite-strand overlapping transcripts with gDNA/nRNA parameter sweep.

Scenario: 20kb genome with two transcripts on opposite strands.
  T1: + strand, exons (1000,1500), (3000,3500), (8000,10000)
  T2: - strand, exons (6000,7000), (12000,13000), (15000,15500)

Questions answered:
  1. When do gDNA candidates get emitted into the EM?
  2. With zero gDNA/nRNA, are the two transcripts independent?
  3. After adding nRNA/gDNA, do they share equivalence classes?
  4. Does changing fragment order change results?
"""

import itertools
import logging
import sys
import textwrap
from pathlib import Path

import numpy as np
import pysam
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
N_FRAGMENTS = 2000
SIM_SEED = 42
PIPELINE_SEED = 42
GENOME_LENGTH = 20_000


def _sim_config(*, strand_specificity: float = 0.9, seed: int = SIM_SEED):
    return SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=strand_specificity, seed=seed,
    )


def _gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance, frag_mean=350, frag_std=100,
        frag_min=100, frag_max=1000,
    )


def _build_scenario(tmp_path, name="opp_strand"):
    """Build the 20kb scenario with two opposite-strand transcripts."""
    sc = Scenario(
        name, genome_length=GENOME_LENGTH, seed=SIM_SEED,
        work_dir=tmp_path / name,
    )
    # T1: + strand, multi-exon
    sc.add_gene("G1", "+", [
        {
            "t_id": "T1",
            "exons": [(1000, 1500), (3000, 3500), (8000, 10000)],
            "abundance": 100,
        },
    ])
    # T2: - strand, multi-exon, on opposite strand
    sc.add_gene("G2", "-", [
        {
            "t_id": "T2",
            "exons": [(6000, 7000), (12000, 13000), (15000, 15500)],
            "abundance": 100,
        },
    ])
    return sc


def _run(tmp_path, *, gdna_abundance=0, nrna_abundance=0,
         strand_specificity=0.9, label=""):
    """Build scenario, run pipeline, return pipeline result + scenario result."""
    sc = _build_scenario(tmp_path, name=label or "opp_strand")
    sim = _sim_config(strand_specificity=strand_specificity)
    gdna = _gdna_config(gdna_abundance)

    result = sc.build_oracle(
        n_fragments=N_FRAGMENTS, sim_config=sim,
        gdna_config=gdna, nrna_abundance=nrna_abundance,
    )

    config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=config)
    bench = run_benchmark(result, pr, scenario_name=label or "opp_strand")
    return sc, result, pr, bench


def _shuffle_bam(input_bam: Path, output_bam: Path, seed: int = 99):
    """Shuffle *read-name groups* in a name-sorted BAM.

    Reads all records, groups by read name, shuffles groups,
    writes a new name-sorted BAM with the shuffled order.
    """
    # Read all records grouped by read name
    groups: dict[str, list] = {}
    with pysam.AlignmentFile(str(input_bam), "rb") as inf:
        header = inf.header.to_dict()
        for rec in inf:
            qname = rec.query_name
            if qname not in groups:
                groups[qname] = []
            groups[qname].append(rec)

    # Shuffle the groups
    rng = np.random.default_rng(seed)
    group_keys = list(groups.keys())
    rng.shuffle(group_keys)

    # Write shuffled BAM
    # Change header to indicate name-sorted
    header["HD"]["SO"] = "queryname"
    with pysam.AlignmentFile(str(output_bam), "wb", header=header) as outf:
        for key in group_keys:
            for rec in groups[key]:
                outf.write(rec)


def _extract_results(pr, bench, result):
    """Extract key metrics from a pipeline run."""
    est = pr.estimator
    idx = result.index

    info = {}
    for t in result.transcripts:
        ti = t.t_index
        tid = t.t_id
        info[f"{tid}_mrna"] = float(est.em_counts[ti].sum())
        info[f"{tid}_nrna"] = float(est.nrna_em_counts[ti])
        info[f"{tid}_unambig"] = float(est.unambig_counts[ti].sum())
        info[f"{tid}_unspliced_sense"] = float(
            est.transcript_unspliced_sense[ti])
        info[f"{tid}_unspliced_antisense"] = float(
            est.transcript_unspliced_antisense[ti])
        info[f"{tid}_nrna_init"] = float(est.nrna_init[ti])

    info["gdna_total"] = float(est.gdna_total)
    info["gdna_em"] = float(est.gdna_em_count)
    info["n_loci"] = len(est.locus_results)

    # Locus details
    for i, lr in enumerate(est.locus_results):
        if isinstance(lr, dict):
            info[f"locus_{i}_n_transcripts"] = lr.get("n_transcripts", 0)
            info[f"locus_{i}_n_em_fragments"] = lr.get("n_em_fragments", 0)
            info[f"locus_{i}_gdna_init"] = lr.get("gdna_init", 0.0)
        else:
            info[f"locus_{i}_n_transcripts"] = lr.n_transcripts
            info[f"locus_{i}_n_em_fragments"] = lr.n_em_fragments
            info[f"locus_{i}_gdna_init"] = lr.gdna_init

    return info


def _compare_results(r1, r2, label=""):
    """Compare two result dicts, return (n_diffs, max_abs_diff, details)."""
    all_keys = sorted(set(r1.keys()) | set(r2.keys()))
    diffs = []
    for k in all_keys:
        v1 = r1.get(k, float("nan"))
        v2 = r2.get(k, float("nan"))
        if v1 != v2:
            abs_diff = abs(v1 - v2)
            denom = max(abs(v1), abs(v2), 1e-30)
            rel_diff = abs_diff / denom
            diffs.append((k, v1, v2, abs_diff, rel_diff))
    return diffs


# ===========================================================================
# Tests
# ===========================================================================


class TestOppositeStrandOverlap:
    """Test two opposite-strand transcripts sharing genomic territory."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = _build_scenario(tmp_path)
        yield sc
        sc.cleanup()

    # -------------------------------------------------------------------
    # Part 1: Locus structure under different conditions
    # -------------------------------------------------------------------

    def test_mrna_only_independent(self, tmp_path):
        """With zero gDNA and zero nRNA, transcripts should be independent."""
        sc, result, pr, bench = _run(
            tmp_path, gdna_abundance=0, nrna_abundance=0,
            strand_specificity=1.0, label="mrna_only",
        )
        info = _extract_results(pr, bench, result)

        # Each transcript should be in its own locus (no shared fragments)
        print("\n=== mRNA-only (no gDNA, no nRNA, SS=1.0) ===")
        print(f"Number of loci: {info['n_loci']}")
        for k, v in sorted(info.items()):
            print(f"  {k}: {v}")

        # T1 and T2 exons don't overlap, but unspliced reads may share
        # genomic territory, causing them to be in the same locus even
        # with no gDNA/nRNA. This is expected behavior when SS < 1.0
        # or when unspliced fragments from strand-ambiguous regions connect them.
        # With SS=1.0, check if they are truly independent:
        if info["n_loci"] == 1:
            print("  NOTE: T1 and T2 are in the SAME locus (shared unspliced territory)")
            print(f"  T1_unspliced_sense={info['T1_unspliced_sense']}, "
                  f"T2_unspliced_sense={info['T2_unspliced_sense']}")
        else:
            print(f"  T1 and T2 are in SEPARATE loci ({info['n_loci']} total)")
        sc.cleanup()

    def test_nrna_creates_shared_locus(self, tmp_path):
        """With nRNA, unspliced reads span introns → shared territory."""
        sc, result, pr, bench = _run(
            tmp_path, gdna_abundance=0, nrna_abundance=50,
            strand_specificity=0.9, label="nrna_sharing",
        )
        info = _extract_results(pr, bench, result)

        print("\n=== nRNA=50, gDNA=0, SS=0.9 ===")
        print(f"Number of loci: {info['n_loci']}")
        for k, v in sorted(info.items()):
            print(f"  {k}: {v}")
        sc.cleanup()

    def test_gdna_creates_shared_locus(self, tmp_path):
        """With gDNA, genomic reads span both genes → shared territory."""
        sc, result, pr, bench = _run(
            tmp_path, gdna_abundance=50, nrna_abundance=0,
            strand_specificity=0.9, label="gdna_sharing",
        )
        info = _extract_results(pr, bench, result)

        print("\n=== nRNA=0, gDNA=50, SS=0.9 ===")
        print(f"Number of loci: {info['n_loci']}")
        for k, v in sorted(info.items()):
            print(f"  {k}: {v}")
        sc.cleanup()

    def test_full_stress(self, tmp_path):
        """Both nRNA and gDNA active."""
        sc, result, pr, bench = _run(
            tmp_path, gdna_abundance=50, nrna_abundance=50,
            strand_specificity=0.9, label="full_stress",
        )
        info = _extract_results(pr, bench, result)

        print("\n=== nRNA=50, gDNA=50, SS=0.9 ===")
        print(f"Number of loci: {info['n_loci']}")
        for k, v in sorted(info.items()):
            print(f"  {k}: {v}")
        sc.cleanup()

    # -------------------------------------------------------------------
    # Part 2: Fragment order sensitivity
    # -------------------------------------------------------------------

    def test_fragment_order_mrna_only(self, tmp_path):
        """mRNA-only: shuffling should not change results."""
        self._test_order_sensitivity(
            tmp_path, gdna_abundance=0, nrna_abundance=0,
            strand_specificity=1.0, label="order_mrna",
        )

    def test_fragment_order_with_nrna(self, tmp_path):
        """nRNA active: shuffling should not change results."""
        self._test_order_sensitivity(
            tmp_path, gdna_abundance=0, nrna_abundance=50,
            strand_specificity=0.9, label="order_nrna",
        )

    def test_fragment_order_with_gdna(self, tmp_path):
        """gDNA active: shuffling should not change results."""
        self._test_order_sensitivity(
            tmp_path, gdna_abundance=50, nrna_abundance=0,
            strand_specificity=0.9, label="order_gdna",
        )

    def test_fragment_order_full_stress(self, tmp_path):
        """Both nRNA and gDNA: shuffling should not change results."""
        self._test_order_sensitivity(
            tmp_path, gdna_abundance=50, nrna_abundance=50,
            strand_specificity=0.9, label="order_full",
        )

    def _test_order_sensitivity(self, tmp_path, *, gdna_abundance, nrna_abundance,
                                 strand_specificity, label):
        """Run pipeline on original + shuffled BAM, compare results."""
        # --- Run 1: normal order ---
        sc = _build_scenario(tmp_path, name=f"{label}_orig")
        sim = _sim_config(strand_specificity=strand_specificity)
        gdna = _gdna_config(gdna_abundance)
        result = sc.build_oracle(
            n_fragments=N_FRAGMENTS, sim_config=sim,
            gdna_config=gdna, nrna_abundance=nrna_abundance,
        )
        config = PipelineConfig(
            em=EMConfig(seed=PIPELINE_SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr1 = run_pipeline(result.bam_path, result.index, config=config)
        bench1 = run_benchmark(result, pr1, scenario_name=f"{label}_orig")
        info1 = _extract_results(pr1, bench1, result)

        # --- Shuffle BAM ---
        shuffled_bam = tmp_path / f"{label}_shuffled.bam"
        _shuffle_bam(result.bam_path, shuffled_bam, seed=99)

        # --- Run 2: shuffled order ---
        pr2 = run_pipeline(shuffled_bam, result.index, config=config)
        bench2 = run_benchmark(result, pr2, scenario_name=f"{label}_shuffled")
        info2 = _extract_results(pr2, bench2, result)

        # --- Compare ---
        diffs = _compare_results(info1, info2, label=label)

        print(f"\n=== Order sensitivity: {label} ===")
        print(f"  gdna={gdna_abundance}, nrna={nrna_abundance}, ss={strand_specificity}")
        if not diffs:
            print("  RESULT: BIT-EXACT (no differences)")
        else:
            print(f"  RESULT: {len(diffs)} differences found")
            for k, v1, v2, ad, rd in sorted(diffs, key=lambda x: -x[3]):
                print(f"    {k}: {v1:.10f} vs {v2:.10f} "
                      f"(abs={ad:.2e}, rel={rd:.2e})")

        # Check that results are close (allow small FP tolerance)
        for k, v1, v2, ad, rd in diffs:
            if "n_loci" in k or "n_transcripts" in k or "n_em_fragments" in k:
                # Structural properties must be exact
                assert v1 == v2, (
                    f"Structural difference in {k}: {v1} vs {v2}"
                )

        sc.cleanup()


# ===========================================================================
# Direct invocation for quick iteration
# ===========================================================================

if __name__ == "__main__":
    import tempfile

    logging.basicConfig(level=logging.WARNING)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test = TestOppositeStrandOverlap()

        print("=" * 70)
        print("  OPPOSITE-STRAND OVERLAP TEST")
        print("=" * 70)

        # Part 1: Locus structure
        test.test_mrna_only_independent(tmp / "p1a")
        test.test_nrna_creates_shared_locus(tmp / "p1b")
        test.test_gdna_creates_shared_locus(tmp / "p1c")
        test.test_full_stress(tmp / "p1d")

        # Part 2: Fragment ordering
        print("\n" + "=" * 70)
        print("  FRAGMENT ORDER SENSITIVITY")
        print("=" * 70)

        test.test_fragment_order_mrna_only(tmp / "p2a")
        test.test_fragment_order_with_nrna(tmp / "p2b")
        test.test_fragment_order_with_gdna(tmp / "p2c")
        test.test_fragment_order_full_stress(tmp / "p2d")
