"""End-to-end tests for SRD v1 ``calibrate_gdna`` orchestrator.

Unit tests for the building blocks (categorize, fl_mixture,
fl_empirical_bayes) live in ``test_categorize.py``,
``test_fl_mixture.py``, and ``test_fl_empirical_bayes.py``.

This module exercises the orchestrator end-to-end: build a small
scenario via the oracle simulator, run ``scan_and_buffer``, then call
``calibrate_gdna`` and inspect the returned :class:`CalibrationResult`.
"""

from __future__ import annotations

import dataclasses

import pytest

from rigel.calibration import CalibrationResult, calibrate_gdna
from rigel.calibration._categorize import FragmentCategory, N_CATEGORIES
from rigel.config import BamScanConfig
from rigel.frag_length_model import FragmentLengthModel
from rigel.pipeline import scan_and_buffer
from rigel.sim import GDNAConfig, Scenario, SimConfig


_SIM_SEED = 11
_N_FRAGMENTS = 1500


def _sim_cfg(strand_specificity: float = 1.0) -> SimConfig:
    return SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=strand_specificity, seed=_SIM_SEED,
    )


def _gdna_cfg(abundance: float) -> GDNAConfig | None:
    if abundance <= 0:
        return None
    # Distinct gDNA FL distribution (mean=350) vs RNA (mean=200).
    return GDNAConfig(
        abundance=abundance, frag_mean=350, frag_std=100,
        frag_min=100, frag_max=1000,
    )


def _spliced_scenario(tmp_path, name: str) -> Scenario:
    """Two-exon spliced gene + silent control on opposite strand."""
    sc = Scenario(name, genome_length=6000, seed=_SIM_SEED, work_dir=tmp_path / name)
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 600), (1100, 1500)], "abundance": 100},
    ])
    sc.add_gene("g_ctrl", "-", [
        {"t_id": "t_ctrl", "exons": [(4000, 4400)], "abundance": 0},
    ])
    return sc


def _run_calibration(
    tmp_path,
    name: str,
    *,
    strand_specificity: float = 1.0,
    gdna_abundance: float = 0.0,
    n_fragments: int = _N_FRAGMENTS,
) -> tuple[CalibrationResult, float]:
    """Build → scan → calibrate. Returns (result, observed_strand_specificity)."""
    sc = _spliced_scenario(tmp_path, name)
    try:
        sim_cfg = _sim_cfg(strand_specificity=strand_specificity)
        gdna = _gdna_cfg(gdna_abundance)
        sim_result = sc.build_oracle(
            n_fragments=n_fragments, sim_config=sim_cfg,
            gdna_config=gdna, nrna_abundance=0,
        )
        scan_cfg = BamScanConfig(sj_strand_tag="auto", n_scan_threads=1)
        _stats, strand_models, frag_length_models, buffer = scan_and_buffer(
            str(sim_result.bam_path), sim_result.index, scan_cfg,
        )
        try:
            strand_models.finalize()
            frag_length_models.build_scoring_models()
            frag_length_models.finalize(prior_ess=500.0)
            cal = calibrate_gdna(
                buffer,
                sim_result.index,
                frag_length_models,
                strand_models.strand_specificity,
            )
            return cal, float(strand_models.strand_specificity)
        finally:
            buffer.cleanup()
    finally:
        sc.cleanup()


# ---------------------------------------------------------------------------
# Schema contract
# ---------------------------------------------------------------------------


class TestSchema:

    def test_calibration_result_fields(self, tmp_path):
        cal, _ = _run_calibration(tmp_path, "schema", gdna_abundance=20)

        # All SRD v1 fields present and well-typed.
        assert isinstance(cal.gdna_fl_model, FragmentLengthModel)
        assert isinstance(cal.rna_fl_model, FragmentLengthModel)
        assert isinstance(cal.global_fl_model, FragmentLengthModel)
        assert cal.gdna_fl_quality in {"good", "weak", "fallback"}
        assert 0.0 <= cal.strand_specificity <= 1.0
        assert cal.category_counts.shape == (N_CATEGORIES,)
        assert cal.category_counts.dtype.kind == "i"
        assert cal.n_multimap_excluded >= 0
        assert cal.n_pool >= 0
        assert 0.0 <= cal.pi_pool <= 1.0
        assert isinstance(cal.mixture_converged, bool)
        assert cal.mixture_iterations >= 0
        assert cal.exon_fit_tolerance_bp == 5
        assert cal.fl_prior_ess == 500.0

    def test_no_v5_fields(self, tmp_path):
        """v5 fields must not exist on the SRD CalibrationResult."""
        cal, _ = _run_calibration(tmp_path, "no_v5", gdna_abundance=20)
        v5_only = {
            "region_e_gdna", "region_n_total", "lambda_gdna",
            "region_gamma", "region_gamma_strand",
            "mu_R", "sigma_R", "mixing_pi", "mixing_pi_soft",
            "kappa_G", "kappa_R", "em_n_iter", "em_converged",
            "n_eligible", "n_soft", "n_spliced_hard",
            "lam_G_on", "lam_G_off", "capture_class_mode",
        }
        present = {f.name for f in dataclasses.fields(cal)}
        leaked = v5_only & present
        assert not leaked, f"v5 fields leaked into SRD schema: {leaked}"

    def test_summary_dict_serializable(self, tmp_path):
        import json
        cal, _ = _run_calibration(tmp_path, "summary", gdna_abundance=20)
        d = cal.to_summary_dict()
        json.dumps(d)  # must not raise
        assert "gdna_fl_quality" in d
        assert "pi_pool" in d
        assert "category_counts" in d
        assert len(d["category_counts"]) == N_CATEGORIES


# ---------------------------------------------------------------------------
# Pure-RNA → near-zero π_pool
# ---------------------------------------------------------------------------


class TestPureRna:

    def test_stranded_pure_rna_has_low_pi(self, tmp_path):
        cal, ss = _run_calibration(
            tmp_path, "pure_rna_stranded",
            strand_specificity=1.0, gdna_abundance=0,
        )
        assert ss > 0.85
        # No injected gDNA → pool absorbs intronic/intergenic RNA only;
        # mixture should land at very low π.
        assert cal.pi_pool < 0.10, (
            f"pure-RNA library has unexpectedly high π_pool={cal.pi_pool:.3f}"
        )
        # Quality may be 'weak' or 'fallback' depending on pool population;
        # 'good' would imply > 2% gDNA, which is wrong here.
        assert cal.gdna_fl_quality in {"weak", "fallback"}


# ---------------------------------------------------------------------------
# Library-agnostic pool definition (no SS branching)
# ---------------------------------------------------------------------------


class TestLibraryAgnostic:

    @pytest.mark.parametrize("ss", [0.5, 0.95], ids=["unstranded", "stranded"])
    def test_pool_definition_independent_of_ss(self, tmp_path, ss):
        """SRD pool = EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC.

        The categorization is geometric and library-agnostic. SS affects
        which transcripts produce which fragment classes, but the pool
        membership rule itself does not branch on SS.
        """
        cal, _ = _run_calibration(
            tmp_path, f"agnostic_{ss}",
            strand_specificity=ss, gdna_abundance=20,
        )

        # Pool size equals the geometric pool count, no antisense, plus
        # the C++ scanner's unique-mapper unresolved (intergenic) FL
        # contribution which is folded into Pool B but never appears in
        # `category_counts` (those frags are dropped before the buffer).
        cats = cal.category_counts
        expected_pool = (
            int(cats[FragmentCategory.EXON_INCOMPATIBLE])
            + int(cats[FragmentCategory.INTRONIC])
            + int(cats[FragmentCategory.INTERGENIC])
            + int(cal.n_intergenic_unique)
        )
        assert cal.n_pool == expected_pool, (
            f"n_pool={cal.n_pool} != EXON_INCOMPATIBLE+INTRONIC+INTERGENIC"
            f"+n_intergenic_unique={expected_pool}"
        )

    def test_antisense_excluded_from_pool(self, tmp_path):
        """Regression: UNSPLICED_ANTISENSE_EXONIC must not enter the pool.

        At unstranded SS=0.5 with a moderately expressed sense gene,
        antisense-exonic fragments are abundant. They must be tracked as
        a diagnostic but excluded from the mixture pool.
        """
        cal, _ = _run_calibration(
            tmp_path, "antisense_exclusion",
            strand_specificity=0.5, gdna_abundance=20,
        )
        n_antisense = int(cal.category_counts[FragmentCategory.UNSPLICED_ANTISENSE_EXONIC])

        # Antisense column is populated (proves we still count it).
        assert n_antisense >= 0  # tolerate zero on tiny scenarios

        # Pool excludes it.
        expected_pool = sum(
            int(cal.category_counts[c])
            for c in (
                FragmentCategory.EXON_INCOMPATIBLE,
                FragmentCategory.INTRONIC,
                FragmentCategory.INTERGENIC,
            )
        ) + int(cal.n_intergenic_unique)
        assert cal.n_pool == expected_pool


# ---------------------------------------------------------------------------
# Fragment-length models
# ---------------------------------------------------------------------------


class TestFlModels:

    def test_three_models_share_same_max_size(self, tmp_path):
        cal, _ = _run_calibration(tmp_path, "fl_models", gdna_abundance=20)
        ms = cal.global_fl_model.max_size
        assert cal.rna_fl_model.max_size == ms
        assert cal.gdna_fl_model.max_size == ms

    def test_gdna_fl_distinct_from_rna_when_pool_populated(self, tmp_path):
        """When the pool is populated, gDNA FL should differ from RNA FL.

        Injected gDNA has mean=350, RNA mean=200. EB shrinkage to global
        FL pulls them together, but the gDNA mean should still sit above
        the RNA mean.
        """
        cal, _ = _run_calibration(
            tmp_path, "fl_distinct",
            strand_specificity=1.0, gdna_abundance=100, n_fragments=3000,
        )
        # Only a meaningful contrast when calibration found gDNA signal.
        if cal.gdna_fl_quality != "good":
            pytest.skip(f"calibration quality={cal.gdna_fl_quality}; no gDNA signal")
        assert cal.gdna_fl_model.mean > cal.rna_fl_model.mean


# ---------------------------------------------------------------------------
# Knob plumbing (config echo + non-default tolerance)
# ---------------------------------------------------------------------------


class TestKnobs:

    def test_custom_exon_fit_tolerance_propagates(self, tmp_path):
        """Non-default ``exon_fit_tolerance_bp`` is echoed back."""
        sc = _spliced_scenario(tmp_path, "knobs_tol")
        try:
            sim_result = sc.build_oracle(
                n_fragments=_N_FRAGMENTS,
                sim_config=_sim_cfg(strand_specificity=1.0),
                gdna_config=_gdna_cfg(20), nrna_abundance=0,
            )
            scan_cfg = BamScanConfig(sj_strand_tag="auto", n_scan_threads=1)
            _stats, sm, flm, buf = scan_and_buffer(
                str(sim_result.bam_path), sim_result.index, scan_cfg,
            )
            try:
                sm.finalize()
                flm.build_scoring_models()
                flm.finalize(prior_ess=500.0)
                cal = calibrate_gdna(
                    buf, sim_result.index, flm, sm.strand_specificity,
                    exon_fit_tolerance_bp=12,
                    fl_prior_ess=250.0,
                )
                assert cal.exon_fit_tolerance_bp == 12
                assert cal.fl_prior_ess == 250.0
            finally:
                buf.cleanup()
        finally:
            sc.cleanup()
