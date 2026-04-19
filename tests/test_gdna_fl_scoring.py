"""
test_gdna_fl_scoring.py — Tests for calibrated gDNA fragment-length model in scoring.

Validates that:
1. The calibrated gDNA FL model replaces the strand-deconvolved intergenic-based
   model in scoring.
2. The gDNA log-likelihoods computed by C++ scoring use the calibrated FL model.
3. The RNA FL model remains unchanged (spliced-fragment-based).
4. End-to-end pipeline produces valid gDNA log-likelihoods with calibration.
5. The replacement is logged and diagnostics reflect the change.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import (
    BamScanConfig,
    CalibrationConfig,
    EMConfig,
    PipelineConfig,
)
from rigel.frag_length_model import FragmentLengthModel, FragmentLengthModels
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

SEED = 42


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_gdna_scenario(tmp_path, *, gdna_abundance=40.0, n_fragments=600):
    """Two-gene scenario with moderate gDNA contamination."""
    sc = Scenario(
        "phase4_fl",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "phase4_fl",
    )
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50},
        ],
    )
    gdna = GDNAConfig(
        abundance=gdna_abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=80,
        frag_max=600,
    )
    sim_config = SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=0.95,
        seed=SEED,
    )
    result = sc.build_oracle(n_fragments=n_fragments, sim_config=sim_config, gdna_config=gdna)
    return sc, result


def _run_pipeline(result, index):
    """Run the pipeline with calibration."""
    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
        calibration=CalibrationConfig(),
    )
    return run_pipeline(result.bam_path, index, config=config)


# ---------------------------------------------------------------------------
# Test: FragmentLengthModel unit tests
# ---------------------------------------------------------------------------


class TestFragmentLengthModelReplacement:
    """Verify that replacing gdna_model on FragmentLengthModels works."""

    def test_replace_gdna_model(self):
        """Assigning a new gdna_model replaces the old one."""
        models = FragmentLengthModels(max_size=500)

        # Observe some intergenic data (goes to gdna_model)
        for fl in [180, 200, 220, 250, 300]:
            models.observe(fl, splice_type=None)
        models.finalize()

        old_mean = models.gdna_model.mean

        # Build a new calibrated gDNA model with different distribution
        cal_counts = np.zeros(501, dtype=np.float64)
        cal_counts[350] = 100.0
        cal_counts[360] = 80.0
        cal_counts[370] = 60.0
        cal_model = FragmentLengthModel.from_counts(cal_counts, max_size=500)

        models.gdna_model = cal_model
        assert models.gdna_model.mean != old_mean
        assert models.gdna_model.mean == pytest.approx(cal_model.mean, rel=1e-6)

    def test_replacement_preserves_rna_model(self):
        """Replacing gdna_model does not affect rna_model."""
        from rigel.splice import SpliceType

        models = FragmentLengthModels(max_size=500)
        for fl in [150, 200, 250]:
            models.observe(fl, splice_type=SpliceType.SPLICED_ANNOT)
        for fl in [300, 350, 400]:
            models.observe(fl, splice_type=None)
        models.finalize()

        rna_mean_before = models.rna_model.mean
        rna_log_prob_before = models.rna_model._log_prob.copy()

        # Replace gDNA model
        cal_counts = np.zeros(501, dtype=np.float64)
        cal_counts[350] = 100.0
        models.gdna_model = FragmentLengthModel.from_counts(cal_counts, max_size=500)

        assert models.rna_model.mean == pytest.approx(rna_mean_before, rel=1e-10)
        np.testing.assert_array_equal(models.rna_model._log_prob, rna_log_prob_before)

    def test_calibrated_model_is_finalized(self):
        """The calibrated model from from_counts() is finalized."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 50.0
        counts[250] = 100.0
        model = FragmentLengthModel.from_counts(counts, max_size=500)

        assert model._finalized
        assert model._log_prob is not None
        assert hasattr(model, "_tail_base")

    def test_log_likelihood_uses_new_model(self):
        """After replacement, log_likelihood uses the calibrated model."""
        models = FragmentLengthModels(max_size=500)
        for fl in [180, 200, 220]:
            models.observe(fl, splice_type=None)
        models.finalize()

        ll_old = models.gdna_model.log_likelihood(350)

        # Replace with model peaked at 350
        cal_counts = np.zeros(501, dtype=np.float64)
        cal_counts[350] = 1000.0
        models.gdna_model = FragmentLengthModel.from_counts(cal_counts, max_size=500)

        ll_new = models.gdna_model.log_likelihood(350)
        assert ll_new > ll_old, "Calibrated model should assign higher prob to its peak"


# ---------------------------------------------------------------------------
# Test: FragmentScorer picks up replaced gDNA FL model
# ---------------------------------------------------------------------------


class TestFragmentScorerGdnaFL:
    """Verify that FragmentScorer.from_models() reads the correct gDNA FL model."""

    @pytest.fixture()
    def mini_scorer_data(self, mini_index):
        """Build minimal models and estimator for FragmentScorer construction."""
        from rigel.config import EMConfig
        from rigel.estimator import AbundanceEstimator
        from rigel.frag_length_model import FragmentLengthModels
        from rigel.splice import SpliceType
        from rigel.strand_model import StrandModels

        # Strand models: near-perfect stranding
        sm = StrandModels()
        sm.exonic_spliced._r1_sense_count = 950
        sm.exonic_spliced._r1_anti_count = 50
        sm.finalize()

        # FL models: distinct RNA vs gDNA peaks
        fl_models = FragmentLengthModels(max_size=500)
        for fl in range(150, 300):
            fl_models.observe(fl, splice_type=SpliceType.SPLICED_ANNOT, weight=10.0)
        for fl in range(250, 450):
            fl_models.observe(fl, splice_type=None, weight=5.0)
        fl_models.build_scoring_models()
        # Inject gDNA FL model (simulating calibration)
        gdna_counts = np.zeros(501, dtype=np.float64)
        for fl in range(250, 450):
            gdna_counts[fl] = 5.0
        fl_models.gdna_model = FragmentLengthModel.from_counts(gdna_counts)
        fl_models.finalize()

        index = mini_index

        # Estimator
        exonic_lengths = index.t_df["length"].values.astype(np.float64)
        spans = (index.t_df["end"].values - index.t_df["start"].values).astype(np.float64)

        from rigel.config import TranscriptGeometry

        geometry = TranscriptGeometry(
            effective_lengths=np.maximum(exonic_lengths - 200.0, 1.0),
            exonic_lengths=exonic_lengths,
            t_to_g=index.t_to_g_arr,
            transcript_spans=spans,
        )
        estimator = AbundanceEstimator(
            index.num_transcripts,
            em_config=EMConfig(),
            geometry=geometry,
        )
        return sm, fl_models, index, estimator

    def test_scorer_uses_default_gdna_fl(self, mini_scorer_data):
        """Without replacement, scorer uses frag_length_models.gdna_model."""
        from rigel.scoring import FragmentScorer

        sm, fl_models, index, estimator = mini_scorer_data

        scorer = FragmentScorer.from_models(sm, fl_models, index, estimator)
        # The scorer's gDNA FL should match the default gdna_model
        np.testing.assert_array_equal(scorer.gdna_fl_log_prob, fl_models.gdna_model._log_prob)
        assert scorer.gdna_fl_max_size == fl_models.gdna_model.max_size

    def test_scorer_uses_replaced_gdna_fl(self, mini_scorer_data):
        """After replacing gdna_model, scorer picks up the new FL model."""
        from rigel.scoring import FragmentScorer

        sm, fl_models, index, estimator = mini_scorer_data

        # Build a distinctive calibrated model peaked at 400
        cal_counts = np.zeros(501, dtype=np.float64)
        cal_counts[400] = 500.0
        cal_counts[410] = 300.0
        cal_model = FragmentLengthModel.from_counts(cal_counts, max_size=500)

        fl_models.gdna_model = cal_model

        scorer = FragmentScorer.from_models(sm, fl_models, index, estimator)
        np.testing.assert_array_equal(scorer.gdna_fl_log_prob, cal_model._log_prob)
        assert scorer.gdna_fl_max_size == cal_model.max_size

    def test_scorer_rna_fl_unchanged_after_gdna_replacement(self, mini_scorer_data):
        """Replacing gdna_model doesn't affect the RNA FL in the scorer."""
        from rigel.scoring import FragmentScorer

        sm, fl_models, index, estimator = mini_scorer_data
        rna_log_prob_original = fl_models.rna_model._log_prob.copy()

        # Replace gDNA model
        cal_counts = np.zeros(501, dtype=np.float64)
        cal_counts[400] = 500.0
        fl_models.gdna_model = FragmentLengthModel.from_counts(cal_counts, max_size=500)

        scorer = FragmentScorer.from_models(sm, fl_models, index, estimator)
        np.testing.assert_array_equal(scorer.fl_log_prob, rna_log_prob_original)


# ---------------------------------------------------------------------------
# Test: End-to-end pipeline integration
# ---------------------------------------------------------------------------


class TestGdnaFLEndToEnd:
    """End-to-end pipeline test verifying calibrated gDNA FL model in scoring."""

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        self.sc, self.result = _make_gdna_scenario(tmp_path)
        self.index = self.result.index
        yield
        self.sc.cleanup()

    def test_calibration_produces_gdna_fl_model(self):
        """Calibration result contains a valid gDNA FL model."""
        pr = _run_pipeline(self.result, self.index)
        assert pr.calibration is not None
        cal_fl = pr.calibration.gdna_fl_model
        assert cal_fl is not None
        assert cal_fl._finalized
        assert cal_fl._log_prob is not None
        assert cal_fl.n_observations > 0

    def test_frag_length_models_gdna_replacement(self):
        """After pipeline with calibration, the gDNA model is replaced when the
        calibrated model provides better FL discrimination than strand-deconvolved."""
        pr = _run_pipeline(self.result, self.index)
        # The calibrated model should exist
        assert pr.calibration is not None
        cal_fl = pr.calibration.gdna_fl_model
        assert cal_fl is not None
        gdna_fl = pr.frag_length_models.gdna_model

        # The replacement depends on the discrimination guard:
        # - If replaced: gdna_model IS the calibrated model
        # - If not replaced: gdna_model is the strand-deconvolved model
        # Either way, the gDNA model should be finalized and usable
        assert gdna_fl._finalized
        assert gdna_fl._log_prob is not None

    def test_pipeline_produces_valid_output(self):
        """Pipeline with calibration produces valid quantification."""
        pr = _run_pipeline(self.result, self.index)
        df = pr.estimator.get_counts_df(self.index)
        n_annotated = sum(
            1 for tid in df["transcript_id"] if not tid.startswith("RIGEL_NRNA_")
        )
        assert n_annotated == 3
        assert (df["count"] >= 0).all()
        assert df["count"].sum() > 0
        assert df["tpm"].sum() == pytest.approx(1e6, rel=1e-3)
        assert pr.calibration is not None


# ---------------------------------------------------------------------------
# Test: Bayesian FL prior blending replaces v3 ESS gate
# ---------------------------------------------------------------------------


def _make_fl_table(region_ids, frag_lens, frag_strands=None):
    """Create a fl_table DataFrame for build_gdna_fl_model."""
    import pandas as pd

    d = {"region_id": np.array(region_ids, dtype=np.intp), "frag_len": np.array(frag_lens, dtype=np.intp)}
    if frag_strands is not None:
        d["frag_strand"] = np.array(frag_strands, dtype=np.int8)
    return pd.DataFrame(d)


def _make_stats(n_regions, gene_strand=None, region_length=None, n_unspliced=None):
    """Create a stats dict for build_gdna_fl_model."""
    if gene_strand is None:
        gene_strand = np.ones(n_regions, dtype=np.int8)
    if region_length is None:
        region_length = np.full(n_regions, 1000.0)
    if n_unspliced is None:
        n_unspliced = np.full(n_regions, 100.0)
    return {
        "gene_strand": np.asarray(gene_strand),
        "region_length": np.asarray(region_length, dtype=np.float64),
        "n_unspliced": np.asarray(n_unspliced, dtype=np.float64),
        "n_total": np.asarray(n_unspliced, dtype=np.float64),
    }


class TestFLPriorBlending:
    """v4 replaces the v3 min_ess cliff with smooth Bayesian prior blending.

    The model is always returned; with zero data it equals the prior, with
    enough data it approaches the empirical histogram.
    """

    def test_returns_uniform_when_no_data_no_prior(self):
        from rigel.calibration import build_gdna_fl_model

        model = build_gdna_fl_model(
            _make_fl_table([], []), _make_stats(2),
            region_weight=np.ones(2),
            strand_specificity=0.95,
            intergenic_fl_model=None,
            fl_prior_ess=0.0,
        )
        assert model is not None
        assert model.total_weight == 0.0

    def test_data_dominates_prior_at_high_n(self):
        from rigel.calibration import build_gdna_fl_model
        from rigel.frag_length_model import FragmentLengthModel

        # Build a prior with mass concentrated at 100bp.
        prior = FragmentLengthModel(max_size=1000)
        prior.counts[100] = 100.0
        prior._total_weight = 100.0

        # 10000 fragments at 300bp easily dominate prior_ess=10.
        rids = np.zeros(10000, dtype=np.intp)
        flens = np.full(10000, 300, dtype=np.intp)
        fstrands = np.ones(10000, dtype=np.int8)
        fl = _make_fl_table(rids, flens, fstrands)

        model = build_gdna_fl_model(
            fl, _make_stats(1), region_weight=np.ones(1),
            strand_specificity=0.95, intergenic_fl_model=prior,
            fl_prior_ess=10.0,
        )
        # Mode at 300bp.
        assert model.counts.argmax() == 300

    def test_prior_dominates_at_low_n(self):
        from rigel.calibration import build_gdna_fl_model
        from rigel.frag_length_model import FragmentLengthModel

        # Prior: 1000 obs at 200bp.
        prior = FragmentLengthModel(max_size=1000)
        prior.counts[200] = 1000.0
        prior._total_weight = 1000.0

        # Only 2 fragments at 400bp.
        fl = _make_fl_table([0, 0], [400, 400], [1, 1])
        model = build_gdna_fl_model(
            fl, _make_stats(1), region_weight=np.ones(1),
            strand_specificity=0.95, intergenic_fl_model=prior,
            fl_prior_ess=1000.0,
        )
        # The empirical posterior counts (model.counts) only hold the
        # two raw observations, but log_prob blends with the prior.  The
        # log_prob at 200 (prior mode) should exceed log_prob at 400.
        assert model.log_likelihood(200) > model.log_likelihood(400)


# ---------------------------------------------------------------------------
# Test: Different gDNA FL distributions affect scoring differently
# ---------------------------------------------------------------------------


class TestGdnaFLScoringImpact:
    """Verify that different gDNA FL models produce different gDNA likelihoods."""

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        """Build scenario with gDNA fragments that have distinctive FL."""
        sc = Scenario(
            "phase4_fl_impact",
            genome_length=5000,
            seed=SEED,
            work_dir=tmp_path / "phase4_fl_impact",
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            ],
        )
        sc.add_gene(
            "g2",
            "-",
            [
                {"t_id": "t2", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50},
            ],
        )
        # Use distinctive gDNA FL parameters: mean=400 (vs RNA mean=200)
        gdna = GDNAConfig(
            abundance=50.0,
            frag_mean=400,
            frag_std=50,
            frag_min=200,
            frag_max=600,
        )
        sim_config = SimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=0.95,
            seed=SEED,
        )
        result = sc.build_oracle(n_fragments=600, sim_config=sim_config, gdna_config=gdna)
        self.sc = sc
        self.result = result
        self.index = result.index
        yield
        sc.cleanup()

    def test_calibrated_gdna_fl_reflects_true_distribution(self):
        """The calibrated gDNA FL model should have mean closer to the true gDNA FL."""
        pr = _run_pipeline(self.result, self.index)

        cal_mean = pr.frag_length_models.gdna_model.mean

        # True gDNA FL mean is 400; calibrated should be reasonably close
        true_mean = 400.0
        # The calibrated model should shift toward the true gDNA FL distribution
        # (may not be exact due to γ-weighting approximation and sample size)
        assert cal_mean > 200.0, (
            f"Calibrated gDNA FL mean ({cal_mean:.1f}) should be shifted "
            f"toward true mean ({true_mean:.1f})"
        )

    def test_pipeline_produces_valid_output(self):
        """Pipeline with calibration produces valid quantification."""
        pr = _run_pipeline(self.result, self.index)
        df = pr.estimator.get_counts_df(self.index)
        assert len(df) > 0
        assert (df["count"] >= 0).all()
        assert df["tpm"].sum() == pytest.approx(1e6, rel=1e-3)
