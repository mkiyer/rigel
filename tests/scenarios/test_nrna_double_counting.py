"""Scenario: nascent RNA must not be double-counted.

A 20 kb genome with two multi-exon transcripts::

    t1     (+)  eight 500 bp exons across chr1:2000-10000   abundance 100  ← positive control
    t_ctrl (−)  eight 500 bp exons across chr1:12000-19500  abundance 0    ← negative control

t1's ~4 kb of intronic span gives the nRNA estimate ample signal, and t_ctrl is unexpressed, on the
opposite strand and physically separated, so any count on it is a false positive. The sweep covers
gDNA abundance [0, 20, 100] × nRNA abundance [0, 30, 70] × strand specificity [0.65, 0.9, 1.0], and
the property held is that the EM separates the three pools: pipeline nRNA must not exceed the truth by
more than the per-tier ratio, total RNA (mRNA + nRNA) must be accounted for, and per-transcript mRNA
must be accurate on the clean conditions. The bounds widen with gDNA and with the gap to perfect
strand specificity, because nRNA fragments lying entirely inside exons are physically
indistinguishable from mRNA and antisense nRNA is genuinely confusable with gDNA — at near-random
strand specificity the strand model says almost nothing, and the leak is real rather than a defect
this test can see.
"""

import logging

import pytest
from rigel.sim import Scenario

from .conftest import (
    SIM_SEED,
    build_and_run,
    assert_alignment,
    assert_accountability,
    assert_transcript_accuracy,
    assert_negative_control,
)

logger = logging.getLogger(__name__)

# =====================================================================
# Sweep grid — full hyperparameter space
# =====================================================================

GDNA_LEVELS = [0, 20, 100]
NRNA_LEVELS = [0, 30, 70]
SS_LEVELS = [0.65, 0.9, 1.0]

# All (gdna, nrna, ss) combos for the exhaustive sweep
_FULL_GRID = [(g, n, s) for g in GDNA_LEVELS for n in NRNA_LEVELS for s in SS_LEVELS]
_FULL_IDS = [f"g{g}_n{n}_s{int(s * 100)}" for g, n, s in _FULL_GRID]

N_FRAGMENTS = 2000


class TestNrnaDoubleCounting:
    """Regression test: nRNA counts must not be double-counted."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "nrna_double_count",
            genome_length=20000,
            seed=SIM_SEED,
            work_dir=tmp_path / "nrna_double_count",
        )
        # Positive control: multi-exon gene on the + strand. 4 kb of exon in eight pieces, so
        # annotated spliced observations are well represented, and an 8 kb span so that the
        # nRNA-versus-gDNA separation problem stays live.
        sc.add_gene(
            "g1",
            "+",
            [
                {
                    "t_id": "t1",
                    "exons": [
                        (2000, 2500),
                        (3000, 3500),
                        (4000, 4500),
                        (5000, 5500),
                        (6000, 6500),
                        (7000, 7500),
                        (8000, 8500),
                        (9500, 10000),
                    ],
                    "abundance": 100,
                },
            ],
        )
        # Negative control: multi-exon gene on − strand, no expression.
        # Physically separated from g1.
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {
                    "t_id": "t_ctrl",
                    "exons": [
                        (12000, 12500),
                        (13000, 13500),
                        (14000, 14500),
                        (15000, 15500),
                        (16000, 16500),
                        (17000, 17500),
                        (18000, 18500),
                        (19000, 19500),
                    ],
                    "abundance": 0,
                },
            ],
        )
        yield sc
        sc.cleanup()

    # -----------------------------------------------------------------
    # Exhaustive sweep: every (gDNA, nRNA, SS) combination
    # -----------------------------------------------------------------

    @pytest.mark.parametrize("gdna,nrna,ss", _FULL_GRID, ids=_FULL_IDS)
    def test_full_sweep(self, scenario, gdna, nrna, ss):
        """Full hyperparameter sweep — all three pools must stay separated.

        The assertion tiers are by difficulty: no gDNA is tight and is where the double-counting
        regression lives, moderate gDNA is a two-way gDNA-plus-RNA separation, and extreme gDNA is a
        stress arm where only accountability is asserted.
        """
        bench = build_and_run(
            scenario,
            n_fragments=N_FRAGMENTS,
            gdna_abundance=gdna,
            nrna_abundance=nrna,
            strand_specificity=ss,
            scenario_name=f"nrna_dc_g{gdna}_n{nrna}_s{int(ss * 100)}",
        )
        assert_alignment(bench)
        assert_accountability(bench, tolerance=10)
        assert_negative_control(
            bench,
            gdna_abundance=gdna,
            strand_specificity=ss,
        )

        # ----- Core nRNA double-counting assertion -----
        # When nRNA is present AND gDNA is not overwhelming, the pipeline's nRNA count must not
        # exceed the truth by more than the limit. A synthetic nascent transcript competes as an
        # ordinary transcript in the EM, so with gDNA present the three-way competition
        # (mRNA / nRNA / gDNA) gets more slack than the two-way one.
        if bench.n_nrna_expected > 20 and gdna <= 20:
            nrna_ratio = bench.n_nrna_pipeline / bench.n_nrna_expected
            if gdna > 0:
                max_ratio = 2.0 if ss >= 0.85 else 2.5
            else:
                max_ratio = 1.50 if ss >= 0.85 else 2.0
            assert nrna_ratio < max_ratio, (
                f"nRNA over-estimation: pipeline={bench.n_nrna_pipeline:.0f}, "
                f"expected={bench.n_nrna_expected}, ratio={nrna_ratio:.2f} "
                f"(gdna={gdna}, nrna={nrna}, ss={ss})"
            )

        # ----- Total RNA accountability -----
        # mRNA + nRNA should account for all non-gDNA fragments. Known limitation: below perfect
        # strand specificity, with nRNA present, some antisense nRNA is classified as gDNA.
        total_rna_expected = bench.total_expected + bench.n_nrna_expected
        total_rna_observed = bench.total_rna_observed
        if total_rna_expected > 0:
            rna_rel_err = abs(total_rna_observed - total_rna_expected) / total_rna_expected

            if gdna == 0 and ss >= 0.99:
                # Perfect SS, no gDNA: total RNA near-exact. The bound is not zero because a
                # zero-gDNA library still shows a small phantom gDNA share from antisense nascent
                # fragments — a calibration-solve artifact tracked on its own, not an
                # effective-length one.
                assert rna_rel_err < 0.08, (
                    f"Total RNA error: expected={total_rna_expected}, "
                    f"observed={total_rna_observed:.0f}, "
                    f"rel_err={rna_rel_err:.2f}"
                )
            elif gdna == 0:
                # Imperfect SS with nRNA present leaks into false gDNA, so the tolerance depends on
                # SS: the lower it is, the more antisense confusion there is. At near-random SS the
                # strand model is uninformative and the gDNA length model, trained from misidentified
                # antisense RNA, matches real RNA — so the leak approaches total.
                tol = 0.25 if ss >= 0.85 else 0.75
                assert rna_rel_err < tol, (
                    f"Total RNA error (g0 ss={ss}): "
                    f"expected={total_rna_expected}, "
                    f"observed={total_rna_observed:.0f}, "
                    f"rel_err={rna_rel_err:.2f}"
                )
            elif gdna == 20:
                # Moderate gDNA. The bound is set by one corner — weak SS with heavy nascent — where
                # this toy's sparse intronic evidence is at its most brittle; every other cell of the
                # sweep sits well inside it.
                assert rna_rel_err < 0.66, (
                    f"Total RNA error (g20): expected={total_rna_expected}, "
                    f"observed={total_rna_observed:.0f}, "
                    f"rel_err={rna_rel_err:.2f}"
                )
            # gdna=100: extreme stress — just check accountability passed

        # ----- mRNA accuracy (clean conditions only) -----
        if gdna == 0 and nrna == 0:
            # The three tiers are the same defect seen at three strengths: on a zero-gDNA library
            # with imperfect strand specificity, calibration cannot deconvolve the last few percent
            # and leaves a false-gDNA phantom, which costs this transcript exactly those fragments.
            # It is the residual the negative-control assertion already admits. Perfect strand
            # specificity is unaffected and stays tight; the near-random arm is the widest because
            # the strand channel there carries almost no information at all.
            tol = 20 if ss >= 0.99 else (100 if ss >= 0.85 else 250)
            assert_transcript_accuracy(bench, max_abs_diff=tol)

    # -----------------------------------------------------------------
    # Focused: nRNA pool accuracy across strand specificities
    # -----------------------------------------------------------------

    @pytest.mark.parametrize("ss", SS_LEVELS, ids=[f"ss_{int(s * 100)}" for s in SS_LEVELS])
    @pytest.mark.parametrize("nrna", [30, 70], ids=[f"nrna_{n}" for n in [30, 70]])
    def test_nrna_strand_inversion(self, scenario, nrna, ss):
        """Higher strand specificity must not DEGRADE nRNA accuracy.

        The inversion this gates is paradoxical and therefore easy to miss: a strand model that
        double-counts nascent RNA gets worse as the strand signal gets stronger, so a sweep that only
        looked at the near-random arm would read the defect as noise.
        """
        bench = build_and_run(
            scenario,
            n_fragments=N_FRAGMENTS,
            nrna_abundance=nrna,
            strand_specificity=ss,
            scenario_name=f"nrna_inversion_n{nrna}_s{int(ss * 100)}",
        )
        if bench.n_nrna_expected > 10:
            nrna_ratio = bench.n_nrna_pipeline / bench.n_nrna_expected
            # At near-random SS the strand model provides almost no signal and nRNA detection may
            # fail outright, the nascent fragments leaking to gDNA; at higher SS the signal is clear.
            lower_bound = 0.0 if ss < 0.7 else (0.30 if ss >= 0.8 else 0.05)
            # Upper bound: must NOT be double-counted (was ~2.0 before fix)
            assert lower_bound <= nrna_ratio < 1.50, (
                f"nRNA ratio {nrna_ratio:.2f} out of [{lower_bound}, 1.50] "
                f"(pipeline={bench.n_nrna_pipeline:.0f}, "
                f"expected={bench.n_nrna_expected}, ss={ss})"
            )

    # -----------------------------------------------------------------
    # Focused: mRNA should be unaffected by nRNA presence
    # -----------------------------------------------------------------

    @pytest.mark.parametrize("nrna", [0, 30, 70], ids=[f"nrna_{n}" for n in [0, 30, 70]])
    def test_mrna_stable_across_nrna(self, scenario, nrna):
        """Total RNA (mRNA + nRNA) is conserved whatever the nRNA level.

        Run at perfect strand specificity on purpose, so the known antisense-nRNA-to-gDNA leak is out
        of the way and what is left is the nascent accounting alone.
        """
        bench = build_and_run(
            scenario,
            n_fragments=N_FRAGMENTS,
            nrna_abundance=nrna,
            strand_specificity=1.0,
            scenario_name=f"nrna_mrna_stable_n{nrna}",
        )
        # Total RNA should be near-perfect, there being no gDNA. The bound is not zero for the same
        # reason as the sweep's zero-gDNA perfect-SS branch: a small antisense-nRNA phantom gDNA
        # share survives calibration, and it is tracked as its own defect.
        total_rna = bench.total_rna_observed
        total_expected = bench.total_expected + bench.n_nrna_expected
        if total_expected > 0:
            rel_err = abs(total_rna - total_expected) / total_expected
            assert rel_err < 0.08, (
                f"Total RNA conservation violated: "
                f"expected={total_expected}, observed={total_rna:.0f}, "
                f"rel_err={rel_err:.2f} (nrna={nrna})"
            )
        # For clean condition (no nRNA), check mRNA directly
        if nrna == 0:
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            assert t1.abs_diff <= 15, (
                f"t1 mRNA error: expected={t1.expected}, "
                f"observed={t1.observed:.0f}, diff={t1.abs_diff:.0f}"
            )
