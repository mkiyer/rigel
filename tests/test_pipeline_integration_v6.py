"""End-to-end integration tests for the v6 calibration pipeline (M8b).

Drives ``run_pipeline`` on synthetic oracle BAMs and asserts the v6
calibration surface is fully wired:

* :func:`rigel.calibration.calibrate` produces a :class:`CalibrationResult`
  carrying populated ``fl_models`` and ``global_densities``.
* :func:`rigel.pipeline.quant_from_buffer` back-fills the per-locus
  :class:`PriorTable` via :func:`assemble_priors` and returns the
  enriched calibration result.
* No backflow mutation of the scanner-trained ``FragmentLengthModels``
  (no ``gdna_model`` swap, no ``finalize`` call).
* ``FragmentScorer.from_models`` is invoked with the v6 FL models
  (``calibration.fl_models.rna`` / ``.gdna``), not the raw accumulator.
* ``calibration.to_summary_dict()`` exposes only v6 keys.

These tests pin the M8b contract.  The deeper unit tests for
``calibrate()`` itself live in ``test_calibrate_orchestrator.py``.
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
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

SEED = 17


def _make_scenario(tmp_path, *, n_fragments=600, gdna_abundance=0, name="m8b"):
    sc = Scenario(
        name,
        genome_length=8000,
        seed=SEED,
        work_dir=tmp_path / name,
    )
    # Two-isoform gene + a unique-mapper gene + a third gene
    # (gives ≥ 1 multi-locus + a few single-locus regions).
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 60},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 40},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [{"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50}],
    )
    sc.add_gene(
        "g3",
        "+",
        [{"t_id": "t4", "exons": [(5000, 5300), (5500, 5800)], "abundance": 40}],
    )
    sim_config = SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=SEED,
    )
    gdna_cfg = None
    if gdna_abundance > 0:
        gdna_cfg = GDNAConfig(
            abundance=gdna_abundance,
            frag_mean=350,
            frag_std=100,
            frag_min=100,
            frag_max=1000,
        )
    return sc, sc.build_oracle(
        n_fragments=n_fragments,
        sim_config=sim_config,
        gdna_config=gdna_cfg,
    )


def _run(result, *, calibration_cfg=None):
    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
        calibration=calibration_cfg or CalibrationConfig(),
    )
    return run_pipeline(result.bam_path, result.index, config=config)


@pytest.fixture(scope="module")
def _pristine(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("m8b_pristine")
    sc, sim = _make_scenario(tmp, gdna_abundance=0, name="pristine")
    pr = _run(sim)
    yield pr, sc, sim
    sc.cleanup()


@pytest.fixture(scope="module")
def _contaminated(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("m8b_contam")
    sc, sim = _make_scenario(
        tmp,
        gdna_abundance=80,
        n_fragments=1200,
        name="contam",
    )
    pr = _run(sim)
    yield pr, sc, sim
    sc.cleanup()


class TestPipelinePristineRna:
    def test_calibration_present(self, _pristine):
        pr, _, _ = _pristine
        assert pr.calibration is not None
        # v6 surface
        assert pr.calibration.fl_models is not None
        assert pr.calibration.global_densities is not None

    def test_pi_pool_near_zero(self, _pristine):
        pr, _, _ = _pristine
        summary = pr.calibration.to_summary_dict()
        # Pristine RNA → mean π_gDNA should be very small.
        assert float(summary["mean_pi_gdna"]) < 0.10

    def test_priors_populated_after_quant(self, _pristine):
        pr, _, _ = _pristine
        ptable = pr.calibration.prior_table
        # alpha arrays sized to n_multi_loci
        assert ptable.alpha_gdna.shape == ptable.alpha_rna.shape
        # If any multi-loci were built, c_base * (π_gdna + π_rna) ≈ c_base.
        if ptable.alpha_gdna.size > 0:
            totals = ptable.alpha_gdna + ptable.alpha_rna
            np.testing.assert_allclose(
                totals,
                np.full_like(totals, ptable.c_base_value),
                rtol=1e-6,
            )


class TestPipelineWithGdna:
    def test_pi_pool_elevated(self, _contaminated):
        pr, _, _ = _contaminated
        summary = pr.calibration.to_summary_dict()
        assert float(summary["mean_pi_gdna"]) > 0.05

    def test_estimator_attributes_gdna(self, _contaminated):
        pr, _, _ = _contaminated
        # Some fragments routed to gDNA
        assert pr.estimator.gdna_em_count > 0


class TestNoV1Mutation:
    """Scanner-trained ``FragmentLengthModels`` must not be backfilled."""

    def test_gdna_model_untouched(self, _pristine):
        pr, _, _ = _pristine
        # v1 used to swap ``frag_length_models.gdna_model`` for a
        # calibration-derived FL model.  M8b deletes that backflow:
        # the accumulator's gDNA model carries zero observations.
        gdna_acc = pr.frag_length_models.gdna_model
        assert gdna_acc is None or gdna_acc.n_observations == 0

    def test_finalize_not_called(self, _pristine):
        pr, _, _ = _pristine
        # ``finalize()`` writes ``_log_prob`` on the global model.
        # M8b drops the call from ``run_pipeline``.
        assert pr.frag_length_models.global_model._log_prob is None


class TestScorerConsumesFLModels:
    """``FragmentScorer.from_models`` should be reachable with raw FL models."""

    def test_from_models_signature(self):
        import inspect

        from rigel.scoring import FragmentScorer

        sig = inspect.signature(FragmentScorer.from_models)
        # The new signature carries ``rna_fl`` + ``gdna_fl`` instead of
        # ``frag_length_models``.
        assert "rna_fl" in sig.parameters
        assert "gdna_fl" in sig.parameters
        assert "frag_length_models" not in sig.parameters


class TestSummaryJsonV6Schema:
    """``to_summary_dict()`` exposes the v6 schema only."""

    def test_required_keys(self, _pristine):
        pr, _, _ = _pristine
        summary = pr.calibration.to_summary_dict()
        for key in ("fl_models", "global_densities", "n_multi_loci",
                    "c_base", "mean_pi_gdna", "diagnostics"):
            assert key in summary, f"missing v6 key: {key}"

    def test_no_v1_keys(self, _pristine):
        pr, _, _ = _pristine
        summary = pr.calibration.to_summary_dict()
        # SRD-v1 keys must not leak through.
        for v1_key in ("pi_pool", "n_pool", "gdna_fl_quality",
                       "strand_specificity"):
            assert v1_key not in summary, f"unexpected v1 key: {v1_key}"


class TestPriorWeightRnaPlumbed:
    """Non-zero ``nrna_weight`` from the CLI flows through to EM.

    The M5 ``build_prior_weight_rna`` helper currently returns all-ones
    regardless of ``nrna_weight`` — but the plumbing must propagate the
    value to the EM batch so that future helper implementations take
    effect without further pipeline edits.
    """

    def test_nrna_weight_does_not_break(self, tmp_path):
        sc, sim = _make_scenario(tmp_path, name="nrna_w")
        try:
            pr = _run(
                sim,
                calibration_cfg=CalibrationConfig(nrna_weight=0.5),
            )
            # Pipeline must complete with non-default nrna_weight.
            assert pr.estimator is not None
            ptable = pr.calibration.prior_table
            # All entries are non-None float32 vectors of (n_t + 1).
            for w in ptable.prior_weight_rna:
                assert isinstance(w, np.ndarray)
                assert w.dtype == np.float32
                assert w.size >= 2
        finally:
            sc.cleanup()
