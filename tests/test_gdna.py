"""Tests for gDNA contamination modeling in locus-level EM.

Verifies that:
- _score_gdna_candidate computes correct log-likelihoods
- gDNA shadow component competes with mRNA in locus EM
- gDNA assignments propagate to correct counters
- StrandModels container works correctly
- PipelineStats gDNA fields are correct
"""

import numpy as np
import pytest

from rigel.splice import SpliceType
from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator
from scoring_helpers import score_gdna_standalone as _score_gdna_candidate
from rigel.frag_length_model import FragmentLengthModel, FragmentLengthModels
from rigel.strand_model import StrandModels
from rigel.types import Strand

from conftest import _UNSPLICED_SENSE, _make_locus_em_data, _run_and_assign


# =====================================================================
# Helpers
# =====================================================================


def _make_frag_length_models(n_obs=200, size=250):
    """FragmentLengthModels with a known distribution for gDNA scoring tests.

    Builds the models the same way the production pipeline does:
    observations → build_scoring_models → inject gDNA → finalize.
    Since these tests exercise gDNA scoring (not fragment length training),
    we use a simple single-peak distribution.
    """
    im = FragmentLengthModels()
    for _ in range(n_obs):
        im.observe(size, SpliceType.SPLICED_ANNOT)
        im.observe(size, SpliceType.UNSPLICED)
        im.observe(size, SpliceType.SPLICED_UNANNOT)
        im.observe(size, None)  # intergenic
    im.build_scoring_models()
    # Inject gDNA model (simulating calibration providing it)
    gdna_counts = np.zeros(im.max_size + 1, dtype=np.float64)
    gdna_counts[size] = float(n_obs)
    im.gdna_model = FragmentLengthModel.from_counts(gdna_counts, max_size=im.max_size)
    im.finalize()
    return im


def _make_strand_models_default():
    """StrandModels with untrained intergenic model (P_strand=0.5)."""
    return StrandModels()


# =====================================================================
# _score_gdna_candidate
# =====================================================================


class TestScoreGDNA:
    """Tests for pipeline._score_gdna_candidate."""

    def test_unspliced_no_penalty(self):
        """UNSPLICED gets splice_penalty=1.0 → log(1)=0."""
        im = _make_frag_length_models()
        cat = int(SpliceType.UNSPLICED)
        score = _score_gdna_candidate(
            int(Strand.POS),
            cat,
            250,
            im,
        )
        # Uses gDNA insert model (injected from calibration)
        gdna_model = im.gdna_model
        assert score < 0
        assert score == pytest.approx(np.log(0.5) + gdna_model.log_likelihood(250))

    def test_spliced_unannot_penalty(self):
        """SPLICED_UNANNOT gets heavy penalty (default 0.01)."""
        im = _make_frag_length_models()
        score_unannot = _score_gdna_candidate(
            int(Strand.POS),
            int(SpliceType.SPLICED_UNANNOT),
            250,
            im,
        )
        score_unspliced = _score_gdna_candidate(
            int(Strand.POS),
            int(SpliceType.UNSPLICED),
            250,
            im,
        )
        assert score_unannot < score_unspliced
        penalty_diff = score_unspliced - score_unannot
        assert penalty_diff == pytest.approx(-np.log(0.01), rel=0.01)

    def test_custom_penalties(self):
        """Custom penalty dict overrides defaults."""
        im = _make_frag_length_models()
        cat = int(SpliceType.UNSPLICED)
        custom = {SpliceType.UNSPLICED: 0.5}
        score = _score_gdna_candidate(
            int(Strand.POS),
            cat,
            250,
            im,
            gdna_splice_penalties=custom,
        )
        gdna_model = im.gdna_model
        expected = np.log(0.5) + gdna_model.log_likelihood(250) + np.log(0.5)
        assert score == pytest.approx(expected)

    def test_zero_frag_length(self):
        """frag_length=0 → P_insert=1.0 (log=0)."""
        im = _make_frag_length_models()
        score = _score_gdna_candidate(
            int(Strand.POS),
            int(SpliceType.UNSPLICED),
            0,
            im,
        )
        assert score == pytest.approx(np.log(0.5))


# =====================================================================
# gDNA component in locus EM
# =====================================================================


class TestLocusGDNATheta:
    """Verify gDNA competes with mRNA in locus EM."""

    def test_gdna_absorbs_with_high_init(self):
        """With large gdna_init and equal likelihoods, gDNA takes share."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode="fractional"))
        bundle = _make_locus_em_data(
            [[0]] * 10,
            num_transcripts=3,
            include_gdna=True,
            alpha_gdna=2.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=1)
        # gDNA should absorb some fragments
        assert pool_counts["gdna"] > 0

    def test_high_gdna_log_lik_biases_assignment(self):
        """gDNA with much higher log-likelihood absorbs most counts."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        # Give gDNA a much higher log-lik than mRNA (mRNA defaults to 0.0)
        bundle = _make_locus_em_data(
            [[0]] * 100,
            num_transcripts=2,
            include_gdna=True,
            alpha_gdna=1.0,
            gdna_log_lik=5.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)

        gdna_count = pool_counts["gdna"]
        mrna_count = rc.em_counts.sum()
        assert gdna_count > mrna_count, (
            f"gDNA count ({gdna_count:.4f}) should dominate over mRNA ({mrna_count:.4f})"
        )

    def test_no_gdna_init_minimal_absorption(self):
        """With zero gdna_init and weak gDNA likelihood, gDNA count is small."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 1000.0
        bundle = _make_locus_em_data(
            [[0]] * 100,
            num_transcripts=2,
            rc=rc,
            include_gdna=True,
            alpha_gdna=0.0,
            gdna_log_lik=-5.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)

        assert pool_counts["gdna"] < 1.0, f"gDNA count too high: {pool_counts['gdna']:.4f}"


# =====================================================================
# gDNA assignment routing in locus EM
# =====================================================================


class TestLocusGDNAAssignment:
    """Verify gDNA assignments go to pool_counts, not em_counts."""

    def test_gdna_assignment_not_in_em_counts(self):
        """Fragments assigned to gDNA are counted in pool_counts return."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0]] * 100,
            log_liks_per_unit=[[-10.0]] * 100,
            num_transcripts=2,
            include_gdna=True,
            alpha_gdna=5.0,
            gdna_log_lik=0.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert gdna_count > 0
        total = rc.em_counts.sum() + gdna_count
        assert total == pytest.approx(100.0, abs=1.0)

    def test_total_counts_preserved(self):
        """em_counts + gdna == n_units."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=3,
            include_nrna=True,
            include_gdna=True,
            alpha_gdna=2.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        total = rc.em_counts.sum() + gdna_count
        assert total == pytest.approx(200.0, abs=1.0)

    def test_strong_transcript_signal_beats_gdna(self):
        """When transcript likelihood >> gDNA, most go to transcript."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 500.0
        bundle = _make_locus_em_data(
            [[0]] * 500,
            log_liks_per_unit=[[0.0]] * 500,
            num_transcripts=2,
            rc=rc,
            include_gdna=True,
            alpha_gdna=1.0,
            gdna_log_lik=-20.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert rc.em_counts[0].sum() > 490
        assert gdna_count < 10

    def test_high_gdna_log_lik_absorbs_fragments(self):
        """With gDNA log-lik >> mRNA log-lik, gDNA takes larger share."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0]] * 500,
            log_liks_per_unit=[[-10.0]] * 500,
            num_transcripts=2,
            include_gdna=True,
            alpha_gdna=1.0,
            gdna_log_lik=0.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert gdna_count > 100, f"Expected gDNA to absorb many fragments, got {gdna_count}"


# =====================================================================
# gDNA locus attribution
# =====================================================================


class TestGDNALocusAttribution:
    """When a fragment goes to gDNA, locus_counts tracks it."""

    def test_locus_counts_populated(self):
        """gDNA assignments populate gdna_locus_counts."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        cc = _UNSPLICED_SENSE
        bundle = _make_locus_em_data(
            [[0]] * 100,
            log_liks_per_unit=[[-10.0]] * 100,
            count_cols_per_unit=[[cc]] * 100,
            num_transcripts=2,
            include_gdna=True,
            alpha_gdna=5.0,
            gdna_log_lik=0.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert gdna_count > 0
        assert rc.gdna_locus_counts[0].sum() > 0

    def test_locus_counts_zero_when_no_gdna(self):
        """No gDNA assignments → locus counts stay zero."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc.unambig_counts[0, 0] = 1000.0
        bundle = _make_locus_em_data(
            [[0]] * 50,
            log_liks_per_unit=[[0.0]] * 50,
            num_transcripts=2,
            rc=rc,
            include_gdna=True,
            alpha_gdna=0.0,
            gdna_log_lik=-50.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        # With fractional assignment, a vanishingly small posterior may
        # leak to gDNA even when it's extremely unlikely (log_lik=-50).
        assert gdna_count < 1e-10
        assert rc.gdna_locus_counts.sum() < 1e-10


# =====================================================================
# gDNA properties on AbundanceEstimator
# =====================================================================


class TestGDNAProperties:
    def test_gdna_total_from_em_total(self):
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc._gdna_em_total = 80.0
        assert rc.gdna_em_count == 80.0

    def test_gdna_contamination_rate(self):
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc._gdna_em_total = 10.0
        rc.unambig_counts[0, 0] = 80.0
        rc.em_counts[0, 0] = 10.0
        assert rc.gdna_contamination_rate == pytest.approx(0.1)

    def test_contamination_rate_zero_when_empty(self):
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        assert rc.gdna_contamination_rate == 0.0


# =====================================================================
# Stats integration
# =====================================================================


class TestStatsGDNA:
    """Verify PipelineStats gDNA field changes."""

    def test_intergenic_split(self):
        from rigel.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 100
        stats.n_intergenic_spliced = 10
        assert stats.n_intergenic == 110

    def test_gdna_unambig_property(self):
        from rigel.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 80
        stats.n_intergenic_spliced = 20
        assert stats.n_gdna_unambig == 100

    def test_gdna_total_property(self):
        from rigel.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 80
        stats.n_intergenic_spliced = 20
        stats.n_gdna_em = 50
        assert stats.n_gdna_total == 150

    def test_to_dict_includes_gdna(self):
        from rigel.stats import PipelineStats

        stats = PipelineStats()
        stats.n_intergenic_unspliced = 10
        stats.n_intergenic_spliced = 5
        stats.n_gdna_em = 3

        d = stats.to_dict()
        assert d["n_intergenic"] == 15
        assert d["n_gdna_unambig"] == 15
        assert d["n_gdna_total"] == 18


# =====================================================================
# Locus gDNA behavior — integration tests
# =====================================================================


class TestLocusGDNABehavior:
    """Integration tests for gDNA behavior in locus EM."""

    def test_two_transcripts_share_one_gdna_shadow(self):
        """Two transcripts in one locus share a single gDNA component."""
        rc = AbundanceEstimator(3, em_config=EMConfig(seed=42, assignment_mode="fractional"))
        bundle = _make_locus_em_data(
            [[0, 1]] * 200,
            num_transcripts=3,
            include_gdna=True,
            alpha_gdna=3.0,
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        assert gdna_count > 0  # gDNA has non-zero count
        total = rc.em_counts.sum() + gdna_count
        assert total == pytest.approx(200.0, abs=1.0)

    def test_gdna_log_lik_determines_absorption(self):
        """Higher gDNA log-likelihood → more fragments absorbed by gDNA."""
        # Low gDNA log-lik + strong RNA prior → no gDNA absorption
        rc_low = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        rc_low.unambig_counts[0, _UNSPLICED_SENSE] = 500.0
        bundle_low = _make_locus_em_data(
            [[0]] * 200,
            num_transcripts=2,
            rc=rc_low,
            include_gdna=True,
            alpha_gdna=1.0,
            gdna_log_lik=-5.0,
        )
        pc_low = _run_and_assign(rc_low, bundle_low, em_iterations=10)
        gc_low = pc_low["gdna"]

        # High gDNA log-lik → gDNA absorbs fragments
        rc_high = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        bundle_high = _make_locus_em_data(
            [[0]] * 200,
            log_liks_per_unit=[[-10.0]] * 200,
            num_transcripts=2,
            include_gdna=True,
            alpha_gdna=1.0,
            gdna_log_lik=0.0,
        )
        pc_high = _run_and_assign(rc_high, bundle_high, em_iterations=10)
        gc_high = pc_high["gdna"]

        assert gc_high > gc_low
        assert gc_high > 100  # most fragments go to gDNA
        assert gc_low < 0.1  # RNA dominates (tiny fractional leak ok)

    def test_gdna_em_count_property(self):
        """AbundanceEstimator.gdna_em_count accumulates across loci."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        assert rc.gdna_em_count == 0.0

        rc._gdna_em_total = 42.0
        assert rc.gdna_em_count == 42.0


# =====================================================================
# VBEM alpha prior interaction with gDNA (Phase 5 tests)
# =====================================================================


class TestVBEMGDNAAlpha:
    """Verify VBEM naturally suppresses/enables gDNA via alpha priors."""

    def test_tiny_alpha_gdna_suppresses_gdna(self):
        """Very small α_gDNA → gDNA component gets negligible weight."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42, assignment_mode="fractional"))
        bundle = _make_locus_em_data(
            [[0]] * 200,
            num_transcripts=2,
            include_gdna=True,
            alpha_gdna=0.001,  # tiny gDNA prior
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        # With tiny prior and no gDNA signal, gDNA should be negligible
        assert gdna_count < 5.0

    def test_zero_alpha_gdna_disables_gdna(self):
        """α_gDNA = 0 → gDNA component disabled (gate removes it)."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42))
        bundle = _make_locus_em_data(
            [[0]] * 200,
            num_transcripts=2,
            include_gdna=True,
            alpha_gdna=0.0,  # zero → disabled
            gdna_log_lik=5.0,  # high gDNA likelihood (should be irrelevant)
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        # With alpha_gdna=0, gate disables gDNA entirely
        assert gdna_count < 1e-10

    def test_large_alpha_gdna_enables_absorption(self):
        """Large α_gDNA + equal likelihoods → gDNA absorbs fragments."""
        rc = AbundanceEstimator(2, em_config=EMConfig(seed=42, assignment_mode="fractional"))
        bundle = _make_locus_em_data(
            [[0]] * 200,
            num_transcripts=2,
            include_gdna=True,
            alpha_gdna=50.0,  # strong gDNA prior
        )
        pool_counts = _run_and_assign(rc, bundle, em_iterations=10)
        gdna_count = pool_counts["gdna"]

        # Strong gDNA prior should pull fragments toward gDNA
        assert gdna_count > 2.0


# =====================================================================
# StrandModels container
# =====================================================================


class TestStrandModelsContainer:
    """Tests for the simplified StrandModels container."""

    def test_default_construction(self):
        from rigel.strand_model import StrandModels

        sm = StrandModels()
        assert sm.exonic_spliced.n_observations == 0
        assert sm.exonic.n_observations == 0
        assert sm.intergenic.n_observations == 0

    def test_delegation_to_exonic_spliced(self):
        from rigel.strand_model import StrandModels
        from rigel.types import Strand

        sm = StrandModels()
        for _ in range(100):
            sm.exonic_spliced.observe(Strand.POS, Strand.POS)
        for _ in range(2):
            sm.exonic_spliced.observe(Strand.POS, Strand.NEG)

        assert sm.p_r1_sense == sm.exonic_spliced.p_r1_sense
        assert sm.strand_specificity == sm.exonic_spliced.strand_specificity
        assert sm.read1_sense == sm.exonic_spliced.read1_sense
        assert sm.n_observations == sm.exonic_spliced.n_observations

    def test_to_dict_structure(self):
        from rigel.strand_model import StrandModels

        sm = StrandModels()
        d = sm.to_dict()
        assert "exonic_spliced" in d
        assert "diagnostics" in d
        assert "exonic" in d["diagnostics"]
        assert "intergenic" in d["diagnostics"]
        assert "observations" in d["exonic_spliced"]

    def test_write_json(self, tmp_path):
        from rigel.strand_model import StrandModels
        from rigel.types import Strand
        import json

        sm = StrandModels()
        for _ in range(50):
            sm.exonic_spliced.observe(Strand.POS, Strand.POS)

        out = tmp_path / "strand_models.json"
        sm.write_json(out)
        with open(out) as f:
            data = json.load(f)
        assert "strand_models" in data
        assert "exonic_spliced" in data["strand_models"]
        assert "diagnostics" in data["strand_models"]

    def test_independent_models(self):
        """Each region's model is independent."""
        from rigel.strand_model import StrandModels
        from rigel.types import Strand

        sm = StrandModels()
        sm.exonic_spliced.observe(Strand.POS, Strand.POS)
        sm.exonic.observe(Strand.POS, Strand.NEG)
        sm.intergenic.observe(Strand.NEG, Strand.NEG)

        assert sm.exonic_spliced.n_observations == 1
        assert sm.exonic.n_observations == 1
        assert sm.intergenic.n_observations == 1

    def test_mle_on_finalize(self):
        """finalize() uses MLE for exonic_spliced."""
        from rigel.strand_model import StrandModels
        from rigel.types import Strand

        sm = StrandModels()
        for _ in range(10):
            sm.exonic_spliced.observe(Strand.POS, Strand.POS)
        sm.finalize()
        # MLE: p_r1_sense = 10/10 = 1.0
        assert sm.p_r1_sense == pytest.approx(1.0)
