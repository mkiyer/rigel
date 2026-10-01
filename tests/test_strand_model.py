"""`rigel.strand_model` — the 2x2 sense/antisense counts, the per-sj table underneath them, and the
posterior they imply.

The counts and their accumulation; the per-sj strand table, which is what the overdispersion fit
reads and which the 2x2 is only the marginal of; the posterior over a fragment's strand; the derived
properties; and the container that holds one model per population. A model that kept only the 2x2
would satisfy every marginal check and silently disable the dispersion estimate.
"""

import numpy as np
import pytest

from rigel.types import Strand
from rigel.strand_model import SJStrandTable, StrandModel, StrandModels


def _labels(pairs):
    """Expand ``[(align, sj, count), ...]`` into the C++ scanner's two parallel label arrays."""
    align, sj = [], []
    for a, s, n in pairs:
        align += [int(a)] * n
        sj += [int(s)] * n
    return np.asarray(align, dtype=np.int8), np.asarray(sj, dtype=np.int8)


def _model(pairs) -> StrandModel:
    return StrandModel.from_labels(*_labels(pairs))


def _table(rows) -> SJStrandTable:
    """Build a table from ``[(motif_strand, n_sense, n_antisense), ...]``."""
    motif = [int(m) for m, _, _ in rows]
    return SJStrandTable(
        ref_id=np.zeros(len(rows), dtype=np.int32),
        start=np.arange(len(rows), dtype=np.int64) * 1000,
        end=np.arange(len(rows), dtype=np.int64) * 1000 + 100,
        motif_strand=np.asarray(motif, dtype=np.int8),
        n_sense=np.asarray([s for _, s, _ in rows], dtype=np.int64),
        n_antisense=np.asarray([a for _, _, a in rows], dtype=np.int64),
    )


class TestStrandModelCounts:
    """The 2×2 contingency table and derived counts."""

    def test_default_counts_are_zero(self):
        sm = StrandModel()
        assert sm.pos_pos == 0
        assert sm.pos_neg == 0
        assert sm.neg_pos == 0
        assert sm.neg_neg == 0
        assert sm.n_observations == 0
        assert sm.sj_table is None

    def test_from_labels_fills_each_cell(self):
        sm = _model(
            [
                (Strand.POS, Strand.POS, 3),
                (Strand.POS, Strand.NEG, 5),
                (Strand.NEG, Strand.POS, 7),
                (Strand.NEG, Strand.NEG, 11),
            ]
        )
        assert (sm.pos_pos, sm.pos_neg, sm.neg_pos, sm.neg_neg) == (3, 5, 7, 11)
        assert sm.n_same == 3 + 11
        assert sm.n_opposite == 5 + 7
        assert sm.n_observations == 26

    def test_is_immutable(self):
        sm = _model([(Strand.POS, Strand.POS, 1)])
        with pytest.raises(Exception):
            sm.pos_pos = 99


class TestSJStrandTable:
    """The per-sj refinement and the marginal identity that licenses it."""

    def test_empty(self):
        t = SJStrandTable.empty()
        assert t.n_sj == 0
        assert t.n_observations == 0
        assert t.contingency() == (0, 0, 0, 0)

    def test_marginal_is_the_2x2(self):
        """The correctness argument: the 2x2 is exactly the table's marginal."""
        t = _table([(Strand.POS, 3, 7), (Strand.POS, 10, 2), (Strand.NEG, 5, 4)])
        pos_pos, pos_neg, neg_pos, neg_neg = t.contingency()
        assert pos_pos == 3 + 10  # sense on motif-POS sj
        assert neg_pos == 7 + 2  # antisense on motif-POS sj
        assert neg_neg == 5  # sense on motif-NEG sj
        assert pos_neg == 4  # antisense on motif-NEG sj

    def test_from_sj_table_agrees_with_from_labels(self):
        """Both constructors, one population: the same fragments give the same 2×2 and κ."""
        rows = [(Strand.POS, 3, 7), (Strand.POS, 10, 2), (Strand.NEG, 5, 4)]
        from_table = StrandModel.from_sj_table(_table(rows))
        # The same fragments as raw labels: motif POS + sense ⇒ (align POS, sj POS), etc.
        from_labels = _model(
            [
                (Strand.POS, Strand.POS, 13),  # sense on motif POS
                (Strand.NEG, Strand.POS, 9),  # antisense on motif POS
                (Strand.NEG, Strand.NEG, 5),  # sense on motif NEG
                (Strand.POS, Strand.NEG, 4),  # antisense on motif NEG
            ]
        )
        assert (from_table.pos_pos, from_table.pos_neg, from_table.neg_pos, from_table.neg_neg) == (
            from_labels.pos_pos,
            from_labels.pos_neg,
            from_labels.neg_pos,
            from_labels.neg_neg,
        )
        assert from_table.p_r1_sense == pytest.approx(from_labels.p_r1_sense)

    def test_depth_and_totals(self):
        t = _table([(Strand.POS, 3, 7), (Strand.NEG, 100, 0)])
        np.testing.assert_array_equal(t.depth, [10, 100])
        assert t.n_observations == 110
        assert t.n_sj == 2
        m = StrandModel.from_sj_table(t)
        assert m.n_same / m.n_observations == pytest.approx(
            103 / 110
        )  # the marginal; κ is TestGenuineKappa's

    def test_to_dict_reports_deep_sj(self):
        t = _table([(Strand.POS, 1, 1), (Strand.POS, 60, 60), (Strand.NEG, 900, 200)])
        d = t.to_dict()
        assert d["n_sj"] == 3
        assert d["n_observations"] == 2 + 120 + 1100
        assert d["n_sj_depth_ge_100"] == 2
        assert d["n_sj_depth_ge_1000"] == 1
        assert d["depth_max"] == 1100

    def test_model_carries_table_through(self):
        t = _table([(Strand.POS, 4, 1)])
        sm = StrandModel.from_sj_table(t)
        assert sm.sj_table is t
        assert sm.n_observations == 5


class TestStrandModelPosterior:
    """MLE probability computation."""

    def test_no_observations(self):
        assert StrandModel().p_r1_sense == 0.5

    def test_strong_fr_library(self):
        """A strongly R1-sense library (most same-direction)."""
        sm = _model([(Strand.POS, Strand.POS, 95), (Strand.POS, Strand.NEG, 5)])
        assert sm.p_r1_sense == pytest.approx(0.95)
        assert sm.strand_specificity > 0.9
        assert sm.read1_sense is True

    def test_strong_rf_library(self):
        """A strongly R1-antisense library (most opposite-direction)."""
        sm = _model([(Strand.POS, Strand.POS, 5), (Strand.POS, Strand.NEG, 95)])
        assert sm.p_r1_sense < 0.1
        assert sm.p_r1_antisense > 0.9
        assert sm.strand_specificity > 0.9
        assert sm.read1_sense is False


class TestStrandModelProperties:
    def test_posterior_variance_with_data(self):
        sm = _model([(Strand.POS, Strand.POS, 80), (Strand.POS, Strand.NEG, 20)])
        # p = 80/100 = 0.8, variance = 0.8 * 0.2 / 100 = 0.0016
        assert sm.posterior_variance() == pytest.approx(0.0016)

    def test_posterior_variance_no_observations(self):
        assert StrandModel().posterior_variance() == 0.25

    def test_posterior_95ci(self):
        sm = _model([(Strand.POS, Strand.POS, 95), (Strand.POS, Strand.NEG, 5)])
        lo, hi = sm.posterior_95ci()
        assert lo < hi
        assert lo > 0.88  # should be heavily skewed towards 1.0
        assert hi <= 1.0


class TestStrandModelsContainer:
    """The StrandModels container and its construction from a scanner dict."""

    def _scan_dict(self, rows, exonic=()):
        t = _table(rows)
        align, sj = _labels(exonic) if exonic else (np.empty(0, np.int8), np.empty(0, np.int8))
        return {
            "sj_ref_id": t.ref_id,
            "sj_start": t.start,
            "sj_end": t.end,
            "sj_motif_strand": t.motif_strand,
            "sj_n_sense": t.n_sense,
            "sj_n_antisense": t.n_antisense,
            "exonic_obs": align,
            "exonic_truth": sj,
        }

    def test_delegation_to_exonic_spliced(self):
        models = StrandModels.from_scan(self._scan_dict([(Strand.POS, 50, 0)]))
        assert models.strand_specificity == models.exonic_spliced.strand_specificity
        assert models.p_r1_sense == models.exonic_spliced.p_r1_sense
        assert models.read1_sense == models.exonic_spliced.read1_sense
        assert models.n_observations == models.exonic_spliced.n_observations

    def test_mle_from_scan(self):
        models = StrandModels.from_scan(self._scan_dict([(Strand.POS, 10, 0)]))
        assert models.p_r1_sense == pytest.approx(1.0)

    def test_spliced_2x2_is_the_table_marginal(self):
        """The container's spliced model has ONE source of truth — the sj table."""
        rows = [(Strand.POS, 30, 3), (Strand.NEG, 12, 5)]
        models = StrandModels.from_scan(self._scan_dict(rows))
        assert models.exonic_spliced.contingency_matches_table()
        assert models.sj_table.n_sj == 2

    def test_zero_observations_warns(self, caplog):
        with caplog.at_level("WARNING"):
            models = StrandModels.from_scan(self._scan_dict([]))
        assert "No spliced strand observations" in caplog.text
        assert models.p_r1_sense == 0.5

    def test_low_observations_warns(self, caplog):
        with caplog.at_level("WARNING"):
            StrandModels.from_scan(self._scan_dict([(Strand.POS, 5, 0)]))
        assert "Only 5 spliced strand observations" in caplog.text

    def test_diagnostic_model_built_from_labels_and_has_no_table(self):
        models = StrandModels.from_scan(
            self._scan_dict([(Strand.POS, 1, 0)], exonic=[(Strand.POS, Strand.POS, 10)])
        )
        assert models.exonic.n_observations == 10
        assert models.exonic.sj_table is None

    def test_default_container_has_an_empty_table(self):
        assert StrandModels().sj_table.n_sj == 0


class TestGenuineKappa:
    """κ from the genuine junctions only (owner, 2026-09-30).

    In a stranded library a junction is genuine RNA (wrong-strand rate κ), a splice artifact —
    misaligned gDNA, wrong-strand rate ½ — or reversed (1 − κ). κ is the genuine class's rate, the
    posterior mean under the shipped Beta(1, 1) prior with the three-class likelihood, so it reduces to
    ``(n_same + 1)/(n_obs + 2)`` wherever no junction needs another class. The strand-live gate reads
    the pooled counts, and an unstranded library keeps its pooled rate.
    """

    @staticmethod
    def _contaminated():
        genuine = [(Strand.POS, 1, 99)] * 200  # 200 deep junctions at a wrong-strand rate of 0.01
        artifacts = (
            [(Strand.POS, 1, 1)] * 300 + [(Strand.POS, 2, 1)] * 50 + [(Strand.POS, 1, 2)] * 50
        )
        reversed_ = [(Strand.POS, 20, 0)] * 10  # every read on the motif-opposite side
        return _table(genuine + artifacts + reversed_)

    def test_kappa_excludes_splice_artifacts_and_reversed_junctions(self):
        from rigel.calibration.strand_balance import fit_strand_balance

        table = self._contaminated()
        model = StrandModel.from_sj_table(table)
        pooled = model.n_same / model.n_observations
        kappa = fit_strand_balance(StrandModels(exonic_spliced=model)).rna_sense_frac
        assert pooled > 0.03  # the artifacts and reversed junctions inflate the pooled rate
        assert kappa == pytest.approx(0.01, abs=0.0015)
        assert model.p_r1_sense == pytest.approx(kappa, rel=0.01)  # the EM reads the same κ
        assert model.contingency_matches_table()  # the 2x2 stays the observed record

    def test_without_contamination_kappa_is_the_shipped_posterior_mean(self):
        from rigel.calibration.strand_balance import fit_strand_balance

        model = StrandModel.from_sj_table(
            _table([(Strand.POS, 3, 297)] * 40 + [(Strand.NEG, 2, 198)] * 30)
        )
        expect = (model.n_same + 1.0) / (model.n_observations + 2.0)
        kappa = fit_strand_balance(StrandModels(exonic_spliced=model)).rna_sense_frac
        assert kappa == pytest.approx(expect, rel=1e-6)

    def test_a_deep_clean_library_resolves_its_narrow_posterior(self):
        """A million spliced reads, every junction genuine: κ is finite and is the shipped posterior mean."""
        model = StrandModel.from_sj_table(_table([(Strand.POS, 10, 990)] * 1000))
        expect = (model.n_same + 1.0) / (model.n_observations + 2.0)
        assert np.isfinite(model.p_r1_sense)
        assert (model.genuine_n_same + 1.0) / (model.n_observations + 2.0) == pytest.approx(
            expect, rel=1e-6
        )

    @staticmethod
    def _independent_maximum(table):
        """The three-class mixture's maximum by a route the fit does not share: a grid over κ and the class shares, then
        a Nelder-Mead polish from its best points. Returns ``(the fit's log-likelihood, the search's, the tolerance)``;
        the tolerance adds what the fit's κ floor can cost, ``floor × reads``."""
        from scipy.optimize import minimize
        from scipy.special import logsumexp

        from rigel.config import CONSTANTS
        from rigel.junction_fit import _genuine_mixture

        kn, cnt = np.unique(
            np.stack([table.n_sense, table.depth], axis=1).astype(float), axis=0, return_counts=True
        )
        k, n, c = kn[:, 0], kn[:, 1], cnt.astype(float)
        if c @ k > 0.5 * (c @ n):
            k = n - k
        fit_loglik = _genuine_mixture(k, n, c)[4]

        def classes(kappa):
            return np.stack(
                [
                    k * np.log(kappa) + (n - k) * np.log1p(-kappa),
                    n * np.log(0.5),
                    k * np.log1p(-kappa) + (n - k) * np.log(kappa),
                ]
            )

        def loglik(kappa, w):
            return float(
                c @ logsumexp(classes(kappa) + np.log(np.maximum(w, 1e-300))[:, None], axis=0)
            )

        g = np.linspace(0.0, 1.0, 51)
        simplex = np.array([(a, b, 1.0 - a - b) for a in g for b in g if a + b <= 1.0 + 1e-12])
        log_w = np.log(np.maximum(simplex, 1e-300))[:, :, None]
        grid = sorted(
            (
                (float(v), kap, w)
                for kap in np.geomspace(1e-12, 0.499, 80)
                for v, w in zip(logsumexp(log_w + classes(kap)[None], axis=1) @ c, simplex)
            ),
            key=lambda r: r[0],
        )
        best = grid[-1][0]
        for _, kap, w in grid[-3:]:
            x0 = np.r_[np.log(kap / (1.0 - kap)), np.log(np.maximum(w, 1e-12))]

            def neg(x):
                kp = 1.0 / (1.0 + np.exp(-x[0]))
                return -loglik(min(kp, 0.5), np.exp(x[1:] - logsumexp(x[1:])))

            fit = minimize(
                neg,
                x0,
                method="Nelder-Mead",
                options={"xatol": 1e-12, "fatol": 1e-13, "maxiter": 3000},
            )
            best = max(best, -fit.fun)
        return fit_loglik, best, 1e-8 + CONSTANTS.junction_fit.kappa_floor * float(c @ n)

    def test_the_fit_reaches_the_mixture_maximum(self):
        """The fit is the three-class mixture's maximum likelihood, by an independent search."""
        fit, best, tol = self._independent_maximum(self._contaminated())
        assert fit >= best - tol

    def test_a_maximum_at_the_kappa_boundary_is_found(self):
        """A near-perfectly stranded library whose only wrong-strand reads are five single-read junctions (a real 1 %
        VCaP draw's shape): the maximum calls them reversed, at κ → 0. EM stalled at a saddle 1.1 units short of it."""
        fit, best, tol = self._independent_maximum(
            _table(
                [(Strand.POS, 0, 1)] * 20000
                + [(Strand.POS, 0, 3)] * 10000
                + [(Strand.POS, 1, 0)] * 5
            )
        )
        assert fit >= best - tol

    def test_a_weakly_stranded_library_is_fit_exactly(self):
        """κ near ½ merges the classes and EM crawled for seconds; the profile fit is exact and bounded there."""
        rows = (
            [(Strand.POS, 9, 11)] * 150
            + [(Strand.POS, 4, 6)] * 150
            + [(Strand.POS, 2, 3)] * 150
            + [(Strand.POS, 1, 1)] * 60
        )
        fit, best, tol = self._independent_maximum(_table(rows))
        assert fit >= best - tol

    def test_an_unstranded_library_keeps_its_pooled_rate(self):
        """The strand-live gate reads the pooled counts. On this unstranded table the three-class fit, applied anyway,
        would read its one-sided junctions as genuine and reversed RNA and put κ near 0."""
        model = StrandModel.from_sj_table(
            _table(
                [(Strand.POS, 3, 0)] * 40
                + [(Strand.POS, 0, 3)] * 40
                + [(Strand.POS, 1, 1)] * 100
                + [(Strand.POS, 1, 0)] * 200
                + [(Strand.POS, 0, 1)] * 200
            )
        )
        assert model.p_r1_sense == model.n_same / model.n_observations == 0.5
