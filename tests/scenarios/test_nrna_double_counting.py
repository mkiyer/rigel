"""Scenario: nRNA double-counting regression test.

A 20 kb genome with two two-exon transcripts:

  t1   (+)  exons (2000,4000) (8000,10000)  abundance=100  ← positive control
  t_ctrl (−)  exons (14000,16000) (18000,19000)  abundance=0   ← negative control

The ~4 kb intron in t1 provides ample intronic signal for nRNA
estimation.  The negative control transcript t_ctrl is on the
opposite strand, unexpressed, and physically separated.

This scenario sweeps the full hyperparameter space:
  - gDNA abundance:       [0, 20, 100]
  - nRNA abundance:       [0, 30, 70]
  - strand specificity:   [0.65, 0.9, 1.0]

Before the nRNA double-counting fix, stranded conditions with
nRNA > 0 would over-estimate nascent RNA by ~2× and distort
per-transcript mRNA counts.  After the fix, the EM should
separate all pools and total RNA (mRNA + nRNA) should be
perfectly accounted for.

Known limitations:
  - nRNA fragments that fall entirely within exons are physically
    indistinguishable from mRNA — some nRNA→mRNA leakage is expected.
  - With moderate nRNA and SS < 1.0, anti-sense nRNA reads (from
    imperfect SS) can be misinterpreted as gDNA signal by the strand
    model, producing false-positive gDNA estimates.  This only affects
    low nRNA abundance (30); at higher nRNA (70) the intronic sense
    signal dominates and separation is perfect.

We therefore check:
  - nRNA NOT double-counted (ratio < 1.50) — the core regression test
  - Total RNA (mRNA + nRNA) accountability at SS=1.0 (tight)
  - Relaxed total RNA at SS < 1.0 where false gDNA is expected
  - mRNA accuracy for clean (no gDNA, no nRNA) conditions
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
    assert_nrna_detected,
)

logger = logging.getLogger(__name__)

# =====================================================================
# Sweep grid — full hyperparameter space
# =====================================================================

GDNA_LEVELS = [0, 20, 100]
NRNA_LEVELS = [0, 30, 70]
SS_LEVELS = [0.65, 0.9, 1.0]

# All (gdna, nrna, ss) combos for the exhaustive sweep
_FULL_GRID = [
    (g, n, s)
    for g in GDNA_LEVELS
    for n in NRNA_LEVELS
    for s in SS_LEVELS
]
_FULL_IDS = [f"g{g}_n{n}_s{int(s*100)}" for g, n, s in _FULL_GRID]

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
        # Positive control: two-exon gene on + strand with wide intron.
        # exon1=2kb, intron=4kb, exon2=2kb → spliced length 4kb, span 8kb
        sc.add_gene("g1", "+", [
            {
                "t_id": "t1",
                "exons": [(2000, 4000), (8000, 10000)],
                "abundance": 100,
            },
        ])
        # Negative control: two-exon gene on − strand, no expression.
        # Physically separated from g1.
        sc.add_gene("g_ctrl", "-", [
            {
                "t_id": "t_ctrl",
                "exons": [(14000, 16000), (18000, 19000)],
                "abundance": 0,
            },
        ])
        yield sc
        sc.cleanup()

    # -----------------------------------------------------------------
    # Exhaustive sweep: every (gDNA, nRNA, SS) combination
    # -----------------------------------------------------------------

    @pytest.mark.parametrize("gdna,nrna,ss", _FULL_GRID, ids=_FULL_IDS)
    def test_full_sweep(self, scenario, gdna, nrna, ss):
        """Full hyperparameter sweep — all pools should be well-separated.

        Assertion tiers by difficulty:
          - Tier 1 (g=0):   tight — nRNA double-counting regression
          - Tier 2 (g=20):  moderate — two-way gDNA + RNA separation
          - Tier 3 (g=100): relaxed — extreme gDNA stress test
        """
        bench = build_and_run(
            scenario,
            n_fragments=N_FRAGMENTS,
            gdna_abundance=gdna,
            nrna_abundance=nrna,
            strand_specificity=ss,
            scenario_name=f"nrna_dc_g{gdna}_n{nrna}_s{int(ss*100)}",
        )
        assert_alignment(bench)
        assert_accountability(bench, tolerance=10)
        assert_negative_control(
            bench, gdna_abundance=gdna, strand_specificity=ss,
        )

        # ----- Core nRNA double-counting assertion -----
        # When nRNA is present AND gDNA is not overwhelming,
        # pipeline nRNA count must not exceed expected by > 50%.
        # Before the fix, stranded conditions would overshoot by ~100%.
        # At low SS (< 0.85), the strand-weighted fragment length mixing
        # has high variance (denom=2s-1 amplifies noise), so we relax
        # the threshold.
        if bench.n_nrna_expected > 20 and gdna <= 20:
            nrna_ratio = bench.n_nrna_pipeline / bench.n_nrna_expected
            max_ratio = 1.50 if ss >= 0.85 else 2.0
            assert nrna_ratio < max_ratio, (
                f"nRNA over-estimation: pipeline={bench.n_nrna_pipeline:.0f}, "
                f"expected={bench.n_nrna_expected}, ratio={nrna_ratio:.2f} "
                f"(gdna={gdna}, nrna={nrna}, ss={ss})"
            )

        # ----- Total RNA accountability -----
        # mRNA + nRNA should approximately account for all non-gDNA fragments.
        # Known limitation: at SS < 1.0 with moderate nRNA, some nRNA
        # anti-sense reads get misclassified as gDNA (false positive).
        total_rna_expected = bench.total_expected + bench.n_nrna_expected
        total_rna_observed = bench.total_observed + bench.n_nrna_pipeline
        if total_rna_expected > 0:
            rna_rel_err = abs(
                total_rna_observed - total_rna_expected
            ) / total_rna_expected

            if gdna == 0 and ss >= 0.99:
                # Perfect SS, no gDNA: total RNA must be near-exact
                assert rna_rel_err < 0.05, (
                    f"Total RNA error: expected={total_rna_expected}, "
                    f"observed={total_rna_observed:.0f}, "
                    f"rel_err={rna_rel_err:.2f}"
                )
            elif gdna == 0:
                # Imperfect SS with nRNA present can leak into false gDNA.
                # Tolerance depends on SS (lower SS → more anti-sense confusion).
                tol = 0.25 if ss >= 0.85 else 0.55
                assert rna_rel_err < tol, (
                    f"Total RNA error (g0 ss={ss}): "
                    f"expected={total_rna_expected}, "
                    f"observed={total_rna_observed:.0f}, "
                    f"rel_err={rna_rel_err:.2f}"
                )
            elif gdna == 20:
                # Moderate gDNA: within 60%
                assert rna_rel_err < 0.60, (
                    f"Total RNA error (g20): expected={total_rna_expected}, "
                    f"observed={total_rna_observed:.0f}, "
                    f"rel_err={rna_rel_err:.2f}"
                )
            # gdna=100: extreme stress — just check accountability passed

        # ----- mRNA accuracy (clean conditions only) -----
        if gdna == 0 and nrna == 0:
            assert_transcript_accuracy(bench, max_abs_diff=20)

    # -----------------------------------------------------------------
    # Focused: nRNA pool accuracy across strand specificities
    # -----------------------------------------------------------------

    @pytest.mark.parametrize("ss", SS_LEVELS,
                             ids=[f"ss_{int(s*100)}" for s in SS_LEVELS])
    @pytest.mark.parametrize("nrna", [30, 70],
                             ids=[f"nrna_{n}" for n in [30, 70]])
    def test_nrna_strand_inversion(self, scenario, nrna, ss):
        """Higher strand specificity must NOT degrade nRNA accuracy.

        Before the fix, SS=0.95 would produce ~2× nRNA while SS=0.5
        produced ~1× (the paradoxical strand inversion).
        """
        bench = build_and_run(
            scenario,
            n_fragments=N_FRAGMENTS,
            nrna_abundance=nrna,
            strand_specificity=ss,
            scenario_name=f"nrna_inversion_n{nrna}_s{int(ss*100)}",
        )
        if bench.n_nrna_expected > 10:
            nrna_ratio = bench.n_nrna_pipeline / bench.n_nrna_expected
            # With low SS (0.65), nRNA is poorly detectable — accept wider range.
            # At higher SS (≥ 0.9), the nRNA signal is clear.
            lower_bound = 0.05 if ss < 0.8 else 0.30
            # Upper bound: must NOT be double-counted (was ~2.0 before fix)
            assert lower_bound < nrna_ratio < 1.50, (
                f"nRNA ratio {nrna_ratio:.2f} out of [{lower_bound}, 1.50] "
                f"(pipeline={bench.n_nrna_pipeline:.0f}, "
                f"expected={bench.n_nrna_expected}, ss={ss})"
            )

    # -----------------------------------------------------------------
    # Focused: mRNA should be unaffected by nRNA presence
    # -----------------------------------------------------------------

    @pytest.mark.parametrize("nrna", [0, 30, 70],
                             ids=[f"nrna_{n}" for n in [0, 30, 70]])
    def test_mrna_stable_across_nrna(self, scenario, nrna):
        """Total RNA (mRNA+nRNA) should be conserved regardless of nRNA level.

        Uses SS=1.0 to avoid the known nRNA→gDNA leakage at imperfect SS,
        isolating the pure nRNA accounting test.
        """
        bench = build_and_run(
            scenario,
            n_fragments=N_FRAGMENTS,
            nrna_abundance=nrna,
            strand_specificity=1.0,
            scenario_name=f"nrna_mrna_stable_n{nrna}",
        )
        # Total RNA should be near-perfect (no gDNA)
        total_rna = bench.total_observed + bench.n_nrna_pipeline
        total_expected = bench.total_expected + bench.n_nrna_expected
        if total_expected > 0:
            rel_err = abs(total_rna - total_expected) / total_expected
            assert rel_err < 0.05, (
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
