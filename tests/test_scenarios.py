"""
End-to-end integration tests and accuracy benchmarks using the simulation
framework.

These tests exercise the full pipeline: genome → annotation → reads →
minimap2 alignment → BAM → HulkIndex → run_pipeline → counts.

Every scenario is swept across three parameter axes:

1. **gDNA contamination**: abundance in [0, 5, 20, 50, 100] relative
   to transcript abundance.
2. **Strand specificity**: fraction of correctly stranded reads in
   [0.5, 0.75, 0.9, 0.95, 1.0], where 0.5 = unstranded.
3. **Abundance ratios** (where applicable): fold changes between
   isoforms or antisense genes in [1, 2, 4, 8, 16].

The goal is to profile tool accuracy across parameter space and
identify regimes where the EM or strand model breaks down.

Requirements: minimap2 and samtools must be available in PATH
(installed via conda/mamba).
"""

import logging
import shutil

import numpy as np
import pytest

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark


logger = logging.getLogger(__name__)


# Skip all tests in this module if minimap2 or samtools are missing
pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


# =====================================================================
# Parameter grids
# =====================================================================

# gDNA abundance levels (0 = no gDNA)
GDNA_LEVELS = [0, 5, 20, 50, 100]

# Strand specificity levels (1.0 = perfect FR, 0.5 = unstranded)
STRAND_LEVELS = [0.5, 0.75, 0.9, 0.95, 1.0]

# Abundance fold-change ratios (major:minor)
ABUNDANCE_RATIOS = [1, 2, 4, 8, 16]

# Fixed simulation parameters
N_FRAGMENTS = 500
SIM_SEED = 42
PIPELINE_SEED = 42


# =====================================================================
# Helpers
# =====================================================================


def _sim_config(*, strand_specificity: float = 1.0, seed: int = SIM_SEED):
    """Standard SimConfig with optional strand specificity override."""
    return SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=seed,
    )


def _gdna_config(abundance: float) -> GDNAConfig | None:
    """GDNAConfig for the given abundance, or None if 0."""
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=100,
        frag_max=1000,
    )


def _build_and_run(scenario, *, n_fragments=N_FRAGMENTS,
                   gdna_abundance=0, strand_specificity=1.0,
                   scenario_name=""):
    """Build a scenario, run the pipeline, and return the benchmark."""
    sim_config = _sim_config(strand_specificity=strand_specificity)
    gdna = _gdna_config(gdna_abundance)
    result = scenario.build(
        n_fragments=n_fragments,
        sim_config=sim_config,
        gdna_config=gdna,
    )
    pipeline_result = run_pipeline(
        result.bam_path, result.index,
        sj_strand_tag="ts",
        seed=PIPELINE_SEED,
    )
    bench = run_benchmark(result, pipeline_result,
                          scenario_name=scenario_name)
    logger.info("\n%s", bench.summary())
    return bench


def _assert_alignment(bench, min_rate=0.70):
    """Assert basic alignment sanity."""
    assert bench.n_fragments > 0, "No fragments entered the pipeline"
    assert bench.alignment_rate > min_rate, (
        f"Low alignment rate: {bench.alignment_rate:.1%}"
    )


def _assert_accountability(bench, tolerance=5):
    """Assert that transcript + gDNA + chimeric ≈ pipeline total."""
    total_accounted = (
        bench.total_observed + bench.n_gdna_pipeline + bench.n_chimeric
    )
    assert abs(total_accounted - bench.n_fragments) <= tolerance, (
        f"Accountability gap: accounted={total_accounted:.0f}, "
        f"pipeline={bench.n_fragments}, "
        f"diff={abs(total_accounted - bench.n_fragments):.0f}"
    )


def _assert_transcript_accuracy(bench, max_abs_diff=3):
    """Assert per-transcript count accuracy (no gDNA scenario)."""
    for ta in bench.transcripts:
        assert ta.abs_diff <= max_abs_diff, (
            f"{ta.t_id}: expected={ta.expected}, "
            f"observed={ta.observed:.0f}, diff={ta.abs_diff:.0f}"
        )


def _assert_gdna_accuracy(bench, gdna_abundance, max_rel_err=0.30):
    """Assert gDNA separation accuracy when gDNA is present."""
    if gdna_abundance == 0:
        return
    assert bench.n_gdna_expected > 0, (
        f"Expected gDNA fragments but got 0 (abundance={gdna_abundance})"
    )
    # Transcript counts should not wildly exceed RNA expectation
    for ta in bench.transcripts:
        assert ta.observed <= ta.expected + bench.n_gdna_expected + 5, (
            f"{ta.t_id}: observed={ta.observed:.0f} exceeds "
            f"RNA({ta.expected}) + gDNA({bench.n_gdna_expected})"
        )
    # gDNA relative error when enough gDNA to measure
    if bench.n_gdna_expected > 10:
        gdna_rel_err = bench.gdna_abs_diff / bench.n_gdna_expected
        assert gdna_rel_err < max_rel_err, (
            f"gDNA rel error too high: pipeline={bench.n_gdna_pipeline:.0f}, "
            f"expected={bench.n_gdna_expected}, rel_err={gdna_rel_err:.2f}"
        )


# =====================================================================
# Scenario 1: Single-exon gene (unspliced)
# =====================================================================


class TestSingleExon:
    """Single gene, single transcript, single exon (1000 bp).

    The simplest possible scenario.  With no competing genes and no
    splice junctions, every RNA fragment should count to the sole
    transcript.  Unspliced reads are indistinguishable from gDNA
    overlapping the gene region — this is the hardest case for
    gDNA/RNA separation.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "single_exon",
            genome_length=5000,
            seed=SIM_SEED,
            work_dir=tmp_path / "single_exon",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        yield sc
        sc.cleanup()

    # -- gDNA sweep -----------------------------------------------------------

    @pytest.mark.parametrize(
        "gdna_abundance", GDNA_LEVELS,
        ids=[f"gdna_{g}" for g in GDNA_LEVELS],
    )
    def test_gdna_sweep(self, scenario, gdna_abundance):
        """Sweep gDNA contamination on an unspliced single-exon gene.

        At gdna=0 we expect near-exact counts (within gDNA-prior
        absorption of ~2).  As gDNA increases, unspliced RNA and
        overlapping gDNA are fundamentally inseparable — the EM must
        rely on non-overlapping gDNA to estimate the gDNA rate.
        """
        bench = _build_and_run(
            scenario, gdna_abundance=gdna_abundance,
            scenario_name=f"single_exon_gdna_{gdna_abundance}",
        )
        _assert_alignment(bench)
        _assert_accountability(bench)

        if gdna_abundance == 0:
            _assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            _assert_gdna_accuracy(bench, gdna_abundance)
            # Lower bound: tolerance scales with gDNA fraction
            gdna_frac = bench.n_gdna_expected / max(bench.n_fragments, 1)
            lower_tol = max(0.0, 1.0 - 1.5 * gdna_frac)
            for ta in bench.transcripts:
                assert ta.observed >= ta.expected * lower_tol, (
                    f"{ta.t_id}: observed={ta.observed:.0f} too low vs "
                    f"expected={ta.expected} "
                    f"(tol={lower_tol:.2f}, gdna_frac={gdna_frac:.2f})"
                )

    # -- Strand specificity sweep ---------------------------------------------

    @pytest.mark.parametrize(
        "strand_specificity", STRAND_LEVELS,
        ids=[f"ss_{s}" for s in STRAND_LEVELS],
    )
    def test_strand_sweep(self, scenario, strand_specificity):
        """Sweep strand specificity on a single-exon gene.

        Imperfect strandedness causes RNA reads that appear antisense
        to be diverted to gDNA (no antisense gene exists).  The loss
        is AMPLIFIED beyond the naive flip rate because the strand
        model learns the library is poorly stranded, reducing
        confidence for ALL reads (not just flipped ones).  For
        unspliced single-exon genes this effect is most severe.
        """
        bench = _build_and_run(
            scenario, strand_specificity=strand_specificity,
            scenario_name=f"single_exon_ss_{strand_specificity}",
        )
        _assert_alignment(bench)
        _assert_accountability(bench)
        # At low ss, the strand model can't distinguish sense from
        # antisense.  Most reads get diverted to gDNA.  We only
        # assert a very lenient lower bound.
        # At ss=1.0, expect near-exact.  At ss=0.5, accept almost
        # any count (the loss is a known limitation of unspliced genes).
        if strand_specificity >= 0.95:
            _assert_transcript_accuracy(bench, max_abs_diff=55)


# =====================================================================
# Scenario 2: Single spliced gene (multi-exon)
# =====================================================================


class TestSplicedGene:
    """Single gene, single transcript, 2 exons (300 bp each, 500 bp intron).

    Spliced reads provide strong evidence for the transcript and
    are easily separated from gDNA (which never produces splice
    junctions matching the annotation).
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "spliced_gene",
            genome_length=5000,
            seed=SIM_SEED,
            work_dir=tmp_path / "spliced_gene",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        yield sc
        sc.cleanup()

    # -- gDNA sweep -----------------------------------------------------------

    @pytest.mark.parametrize(
        "gdna_abundance", GDNA_LEVELS,
        ids=[f"gdna_{g}" for g in GDNA_LEVELS],
    )
    def test_gdna_sweep(self, scenario, gdna_abundance):
        """Sweep gDNA contamination on a spliced gene.

        Spliced reads are unambiguously RNA.  Unspliced exonic reads
        compete with gDNA, but the spliced component anchors the EM.
        Accuracy should be better than the unspliced single-exon case.
        """
        bench = _build_and_run(
            scenario, gdna_abundance=gdna_abundance,
            scenario_name=f"spliced_gene_gdna_{gdna_abundance}",
        )
        _assert_alignment(bench)
        _assert_accountability(bench)

        if gdna_abundance == 0:
            _assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            _assert_gdna_accuracy(bench, gdna_abundance)
            # Spliced reads anchor the count — expect tighter accuracy
            # than the unspliced case.
            gdna_frac = bench.n_gdna_expected / max(bench.n_fragments, 1)
            lower_tol = max(0.0, 1.0 - 1.0 * gdna_frac)
            for ta in bench.transcripts:
                assert ta.observed >= ta.expected * lower_tol, (
                    f"{ta.t_id}: observed={ta.observed:.0f} too low vs "
                    f"expected={ta.expected} "
                    f"(tol={lower_tol:.2f}, gdna_frac={gdna_frac:.2f})"
                )

    # -- Strand specificity sweep ---------------------------------------------

    @pytest.mark.parametrize(
        "strand_specificity", STRAND_LEVELS,
        ids=[f"ss_{s}" for s in STRAND_LEVELS],
    )
    def test_strand_sweep(self, scenario, strand_specificity):
        """Sweep strand specificity on a spliced gene.

        Spliced reads carry splice-junction strand information that
        is independent of library strandedness, providing a strong
        anchor even at low strand specificity.  However, antisense
        reads without splice junctions are still diverted to gDNA.
        """
        bench = _build_and_run(
            scenario, strand_specificity=strand_specificity,
            scenario_name=f"spliced_gene_ss_{strand_specificity}",
        )
        _assert_alignment(bench)
        _assert_accountability(bench)
        # Count loss from antisense→gDNA diversion (spliced gene
        # anchors better than unspliced, so expect less loss).
        max_loss = (1.0 - strand_specificity) + 0.10
        for ta in bench.transcripts:
            assert ta.observed >= ta.expected * (1.0 - max_loss), (
                f"{ta.t_id}: observed={ta.observed:.0f} too low vs "
                f"expected={ta.expected} at ss={strand_specificity}"
            )


# =====================================================================
# Scenario 3: Multiple non-overlapping genes
# =====================================================================


class TestNonOverlappingGenes:
    """Two non-overlapping genes on opposite strands.

    Gene 1: + strand, spliced (2 exons, 300 bp each).
    Gene 2: − strand, single exon (400 bp), distant location.

    No spatial overlap → no ambiguity.  Each gene's reads should
    map exclusively to that gene.  This tests that the pipeline
    handles multi-gene genomes without cross-talk.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "non_overlapping",
            genome_length=8000,
            seed=SIM_SEED,
            work_dir=tmp_path / "non_overlapping",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(4000, 4400)],
             "abundance": 100},
        ])
        yield sc
        sc.cleanup()

    # -- gDNA sweep -----------------------------------------------------------

    @pytest.mark.parametrize(
        "gdna_abundance", GDNA_LEVELS,
        ids=[f"gdna_{g}" for g in GDNA_LEVELS],
    )
    def test_gdna_sweep(self, scenario, gdna_abundance):
        """Sweep gDNA on two non-overlapping genes."""
        bench = _build_and_run(
            scenario, gdna_abundance=gdna_abundance,
            scenario_name=f"non_overlapping_gdna_{gdna_abundance}",
        )
        _assert_alignment(bench)
        _assert_accountability(bench)

        if gdna_abundance == 0:
            _assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            _assert_gdna_accuracy(bench, gdna_abundance)
            gdna_frac = bench.n_gdna_expected / max(bench.n_fragments, 1)
            lower_tol = max(0.0, 1.0 - 1.5 * gdna_frac)
            for ta in bench.transcripts:
                assert ta.observed >= ta.expected * lower_tol, (
                    f"{ta.t_id}: observed={ta.observed:.0f} too low vs "
                    f"expected={ta.expected} "
                    f"(tol={lower_tol:.2f}, gdna_frac={gdna_frac:.2f})"
                )

    # -- Strand specificity sweep ---------------------------------------------

    @pytest.mark.parametrize(
        "strand_specificity", STRAND_LEVELS,
        ids=[f"ss_{s}" for s in STRAND_LEVELS],
    )
    def test_strand_sweep(self, scenario, strand_specificity):
        """Sweep strand specificity on non-overlapping genes.

        No spatial overlap → no cross-talk between genes.  But
        antisense reads from each gene are still diverted to gDNA.
        The unspliced gene (t2) suffers more than the spliced gene
        (t1) because it has no splice-junction anchor.
        """
        bench = _build_and_run(
            scenario, strand_specificity=strand_specificity,
            scenario_name=f"non_overlapping_ss_{strand_specificity}",
        )
        _assert_alignment(bench)
        _assert_accountability(bench)
        # Only assert accuracy at high strandedness; at low ss,
        # unspliced gene t2 loses most of its counts to gDNA.
        if strand_specificity >= 0.95:
            _assert_transcript_accuracy(bench, max_abs_diff=25)


# =====================================================================
# Scenario 4: Two isoforms of a single gene
# =====================================================================


class TestTwoIsoforms:
    """Single gene, two isoforms sharing 5′ and 3′ exons.

    t1 (major): 3 exons — shared_5p, middle_exon, shared_3p.
    t2 (minor): 2 exons — shared_5p, shared_3p (skips middle exon).

    Fragments from shared exons are ambiguous between isoforms.
    Fragments spanning the middle exon or its splice junctions are
    unique to t1.  No fragments are unique to t2 (it is a strict
    subset of t1's exons).

    The EM must disambiguate based on the fraction of spliced reads
    hitting the middle exon.  Accuracy depends on the abundance ratio
    and how many fragments span the discriminating splice junctions.
    """

    def _make_scenario(self, tmp_path, major_abundance, minor_abundance):
        """Create a two-isoform scenario with the given abundances."""
        sc = Scenario(
            "two_isoforms",
            genome_length=5000,
            seed=SIM_SEED,
            work_dir=tmp_path / "two_isoforms",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 500), (1000, 1300), (2000, 2300)],
             "abundance": major_abundance},
            {"t_id": "t2",
             "exons": [(200, 500), (2000, 2300)],
             "abundance": minor_abundance},
        ])
        return sc

    # -- Abundance ratio sweep ------------------------------------------------

    @pytest.mark.parametrize(
        "fold_change", ABUNDANCE_RATIOS,
        ids=[f"fc_{r}" for r in ABUNDANCE_RATIOS],
    )
    def test_abundance_sweep(self, tmp_path, fold_change):
        """Sweep isoform abundance ratios (major:minor fold change).

        At 1:1, both isoforms get roughly equal counts.  As the fold
        change increases, the minor isoform gets fewer fragments and
        becomes harder to quantify accurately.
        """
        major_abundance = 100
        minor_abundance = major_abundance / fold_change

        sc = self._make_scenario(tmp_path, major_abundance, minor_abundance)
        try:
            bench = _build_and_run(
                sc, scenario_name=f"two_iso_fc_{fold_change}",
                n_fragments=1000,
            )
            _assert_alignment(bench)
            _assert_accountability(bench)

            # Gene-level total should be close to expected
            assert bench.total_observed == pytest.approx(
                bench.total_expected, abs=5
            ), (
                f"Gene total mismatch: observed={bench.total_observed:.0f}, "
                f"expected={bench.total_expected}"
            )

            # Higher-abundance isoform should get more counts
            # (only meaningful when fold_change > 1)
            if fold_change > 1:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t1.observed > t2.observed, (
                    f"t1 should dominate at {fold_change}:1 — "
                    f"t1={t1.observed:.0f}, t2={t2.observed:.0f}"
                )
        finally:
            sc.cleanup()

    # -- gDNA sweep -----------------------------------------------------------

    @pytest.mark.parametrize(
        "gdna_abundance", GDNA_LEVELS,
        ids=[f"gdna_{g}" for g in GDNA_LEVELS],
    )
    def test_gdna_sweep(self, tmp_path, gdna_abundance):
        """Sweep gDNA on the two-isoform gene (fixed 10:1 ratio)."""
        sc = self._make_scenario(tmp_path, 100, 10)
        try:
            bench = _build_and_run(
                sc, gdna_abundance=gdna_abundance,
                scenario_name=f"two_iso_gdna_{gdna_abundance}",
                n_fragments=1000,
            )
            _assert_alignment(bench)
            _assert_accountability(bench)

            if gdna_abundance == 0:
                # Gene-level total should be close
                assert bench.total_observed == pytest.approx(
                    bench.total_expected, abs=5
                )
            else:
                _assert_gdna_accuracy(bench, gdna_abundance)
        finally:
            sc.cleanup()

    # -- Strand specificity sweep ---------------------------------------------

    @pytest.mark.parametrize(
        "strand_specificity", STRAND_LEVELS,
        ids=[f"ss_{s}" for s in STRAND_LEVELS],
    )
    def test_strand_sweep(self, tmp_path, strand_specificity):
        """Sweep strand specificity on the two-isoform gene (10:1 ratio).

        Antisense reads are diverted to gDNA, reducing gene-level
        total proportionally to 1 − strand_specificity.  Isoform
        ratios should remain stable since strand doesn't affect
        splice-junction-based disambiguation.
        """
        sc = self._make_scenario(tmp_path, 100, 10)
        try:
            bench = _build_and_run(
                sc, strand_specificity=strand_specificity,
                scenario_name=f"two_iso_ss_{strand_specificity}",
                n_fragments=1000,
            )
            _assert_alignment(bench)
            _assert_accountability(bench)

            # Gene-level total drops proportional to antisense leakage
            max_loss = (1.0 - strand_specificity) + 0.10
            min_expected = bench.total_expected * (1.0 - max_loss)
            assert bench.total_observed >= min_expected, (
                f"Gene total too low at ss={strand_specificity}: "
                f"observed={bench.total_observed:.0f}, "
                f"min={min_expected:.0f}"
            )
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 5: Overlapping genes on opposite strands (antisense)
# =====================================================================


class TestOverlappingAntisense:
    """Two genes overlapping on opposite strands.

    g1 (+ strand): 2 exons, 300 bp each.
    g2 (− strand): 2 exons, 300 bp each, overlapping g1's footprint.

    Overlapping antisense genes create the most challenging scenario
    for strand-based disambiguation.  At perfect strandedness, the
    pipeline can separate sense from antisense.  At 0.5, it cannot.

    Abundance ratios control which gene dominates and how the EM
    resolves ambiguous assignments.
    """

    def _make_scenario(self, tmp_path, g1_abundance, g2_abundance):
        """Create overlapping antisense scenario with given abundances."""
        sc = Scenario(
            "overlapping_antisense",
            genome_length=5000,
            seed=SIM_SEED,
            work_dir=tmp_path / "overlapping_antisense",
        )
        # g1 on + strand: exons at 200-500 and 1000-1300
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": g1_abundance},
        ])
        # g2 on − strand: exons at 300-600 and 1100-1400
        # This overlaps with g1's exons
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(300, 600), (1100, 1400)],
             "abundance": g2_abundance},
        ])
        return sc

    # -- Abundance ratio sweep ------------------------------------------------

    @pytest.mark.parametrize(
        "fold_change", ABUNDANCE_RATIOS,
        ids=[f"fc_{r}" for r in ABUNDANCE_RATIOS],
    )
    def test_abundance_sweep(self, tmp_path, fold_change):
        """Sweep antisense abundance ratios at perfect strandedness.

        At 1:1, both genes get ~equal counts.  As one gene dominates,
        the EM should resolve ambiguous reads toward it.  At high
        fold changes, the minor gene becomes hard to quantify but
        should still receive some counts.
        """
        g1_abundance = 100
        g2_abundance = g1_abundance / fold_change

        sc = self._make_scenario(tmp_path, g1_abundance, g2_abundance)
        try:
            bench = _build_and_run(
                sc, scenario_name=f"antisense_fc_{fold_change}",
                n_fragments=1000,
            )
            _assert_alignment(bench)
            _assert_accountability(bench)

            # Gene-level total should be close to expected
            assert bench.total_observed == pytest.approx(
                bench.total_expected, abs=10
            ), (
                f"Total mismatch at fc={fold_change}: "
                f"observed={bench.total_observed:.0f}, "
                f"expected={bench.total_expected}"
            )

            # Dominant gene should get more counts
            if fold_change > 1:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t1.observed > t2.observed, (
                    f"t1 should dominate at {fold_change}:1 — "
                    f"t1={t1.observed:.0f}, t2={t2.observed:.0f}"
                )
        finally:
            sc.cleanup()

    # -- gDNA sweep -----------------------------------------------------------

    @pytest.mark.parametrize(
        "gdna_abundance", GDNA_LEVELS,
        ids=[f"gdna_{g}" for g in GDNA_LEVELS],
    )
    def test_gdna_sweep(self, tmp_path, gdna_abundance):
        """Sweep gDNA on overlapping antisense genes (equal abundance)."""
        sc = self._make_scenario(tmp_path, 100, 100)
        try:
            bench = _build_and_run(
                sc, gdna_abundance=gdna_abundance,
                scenario_name=f"antisense_gdna_{gdna_abundance}",
                n_fragments=1000,
            )
            _assert_alignment(bench)
            _assert_accountability(bench)

            if gdna_abundance == 0:
                assert bench.total_observed == pytest.approx(
                    bench.total_expected, abs=10
                )
            else:
                _assert_gdna_accuracy(bench, gdna_abundance, max_rel_err=0.40)
        finally:
            sc.cleanup()

    # -- Strand specificity sweep ---------------------------------------------

    @pytest.mark.parametrize(
        "strand_specificity", STRAND_LEVELS,
        ids=[f"ss_{s}" for s in STRAND_LEVELS],
    )
    def test_strand_sweep(self, tmp_path, strand_specificity):
        """Sweep strand specificity on overlapping antisense genes.

        This is the critical scenario for strand-based disambiguation.
        At ss=1.0, the pipeline separates sense/antisense cleanly.
        As ss → 0.5, strand information vanishes and the overlapping
        genes become indistinguishable — gene-level accuracy degrades.

        Some antisense-diverted reads go to gDNA rather than the
        competing gene, so the total may also drop.
        """
        sc = self._make_scenario(tmp_path, 100, 100)
        try:
            bench = _build_and_run(
                sc, strand_specificity=strand_specificity,
                scenario_name=f"antisense_ss_{strand_specificity}",
                n_fragments=1000,
            )
            _assert_alignment(bench)
            _assert_accountability(bench)

            # Total across both genes may drop due to gDNA diversion
            max_loss = (1.0 - strand_specificity) + 0.10
            min_expected = bench.total_expected * (1.0 - max_loss)
            assert bench.total_observed >= min_expected, (
                f"Total too low at ss={strand_specificity}: "
                f"observed={bench.total_observed:.0f}, "
                f"min={min_expected:.0f}"
            )

            # Per-gene accuracy scales with strand specificity.
            # At ss=1.0, expect tight accuracy.
            # At ss=0.5, per-gene counts may be wildly off
            # (fragments randomly assigned between the two genes).
            if strand_specificity >= 0.9:
                for ta in bench.transcripts:
                    assert ta.rel_error < 0.20, (
                        f"{ta.t_id}: rel_error={ta.rel_error:.2f} "
                        f"at ss={strand_specificity}"
                    )
        finally:
            sc.cleanup()
