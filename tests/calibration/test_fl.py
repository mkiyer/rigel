"""gDNA / RNA fragment-length laws: the pools, RNA's own-count law, gDNA's smooth-EB build and two-pool contrast, and
each law's sampling variance.

The pools are defined by structure and none is pure: only the RNA pool is certified, since gDNA cannot splice, while RNA
inside introns and unannotated transcription reach the gDNA pools too. The partition tests hold the pool accessors
exhaustive and disjoint. RNA's law is its own counts', read at their junction-crossing chance: a small pool is not
pulled toward the majority's law, a length the pool could reach and did not is evidence, an unreachable one is not, and
an empty pool is the global anchor exactly. gDNA's build shrinks a sparse pool smoothly toward the global anchor with no
threshold, and the accessor tests cover the raw empirical views the QC report reads. ``gdna_pmf`` is the uniform-frame
law the chemistry makes and the opportunity and prior mathematics assume: the contrast's alone. Each component pmf
carries its per-length sampling variance (RNA's from its posterior, gDNA's its raw law's plug-in at its empirical-Bayes
weight squared), which the EM's two length tables read to tell a difference from noise.
"""

from __future__ import annotations

from functools import partial

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.fl import (
    fl_models_from_histograms,
    gdna_fl_mass,
    rna_fl_mass,
)
from rigel.scan_payload import (
    N_FRAGMENT_POOLS,
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTERGENIC_EXON,
    POOL_DNA_INTRONIC,
    POOL_DNA_INTRON_EXON,
    POOL_RNA_SPLICED,
)


def _spike(at: int, total: float, n: int = 1001) -> np.ndarray:
    c = np.zeros(n, dtype=np.float64)
    c[at] = total
    return c


def _pools(n_bins: int = 5) -> np.ndarray:
    """A ``pool_lengths`` block with a different value in every pool, so a wrong index cannot pass."""
    pools = np.zeros((N_FRAGMENT_POOLS, n_bins), dtype=np.int64)
    pools[POOL_DNA_INTERGENIC, 2] = 1
    pools[POOL_DNA_INTRONIC, 3] = 3
    pools[POOL_DNA_INTRON_EXON, 2] = 700
    pools[POOL_DNA_INTERGENIC_EXON, 3] = 900
    pools[POOL_RNA_SPLICED, 4] = 11
    return pools


def test_the_pool_indices_are_the_specifications_fragment_pool():
    """The reference enum, the C++ enum and these constants are one axis. A disagreement would re-label every pool (the
    gDNA model fitted from the RNA pool) with nothing looking wrong, so this checks against the executable specification
    itself.
    """
    from native._accumulator_reference import FragmentPool

    assert N_FRAGMENT_POOLS == len(FragmentPool)
    assert POOL_DNA_INTERGENIC == FragmentPool.DNA_INTERGENIC
    assert POOL_DNA_INTRONIC == FragmentPool.DNA_INTRONIC
    assert POOL_DNA_INTRON_EXON == FragmentPool.DNA_INTRON_EXON
    assert POOL_DNA_INTERGENIC_EXON == FragmentPool.DNA_INTERGENIC_EXON
    assert POOL_RNA_SPLICED == FragmentPool.RNA_SPLICED


def test_gdna_fl_mass_is_all_four_gdna_pools():
    """All four, because the contained pair alone measures the short half of one population. Under hybrid capture the
    surviving off-target gDNA sits beside a probe, and a fragment beside a probe reaches the exon boundary, so it
    becomes crossing rather than contained; the contained pair alone reads short under capture, while all four, each
    divided by its own opportunity, stay close to truth in both regimes. The histogram means something only with the
    matching divisor: the four pools tilt in opposite directions, so the raw sum is biased long, and nothing consumes it
    alone.
    """
    g = gdna_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(g, [0.0, 0.0, 701.0, 903.0, 0.0])


def test_gdna_fl_mass_excludes_the_rna_pool():
    """The gDNA length model must never see a certified-RNA fragment."""
    pools = np.zeros((N_FRAGMENT_POOLS, 5), dtype=np.int64)
    pools[POOL_RNA_SPLICED, 1] = 12345
    assert float(gdna_fl_mass(SimpleNamespace(pool_lengths=pools)).sum()) == 0.0


def test_rna_fl_mass_is_the_annotated_sj_pool_alone():
    """gDNA cannot be spliced, so a fragment on an annotated sj certifies RNA, and only that pool does. The accumulator
    admits a splice whether sequenced or implied by the one hypothesis that survived arbitration: pool membership is by
    determinacy, not provenance."""
    r = rna_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(r, [0.0, 0.0, 0.0, 0.0, 11.0])


def test_the_two_component_accessors_partition_every_pool():
    """``gdna_fl_mass`` and ``rna_fl_mass`` are the two component accessors and must be exhaustive and disjoint: a pool
    no accessor returns is discarded evidence, and one two return is double-counted; neither shows in a single-accessor
    test.
    """
    payload = SimpleNamespace(pool_lengths=_pools())
    total = float(gdna_fl_mass(payload).sum()) + float(rna_fl_mass(payload).sum())
    assert total == float(_pools().sum()), "the two component accessors must be exhaustive"

    for pool in range(N_FRAGMENT_POOLS):
        one = np.zeros((N_FRAGMENT_POOLS, 5), dtype=np.int64)
        one[pool, 1] = 5
        p = SimpleNamespace(pool_lengths=one)
        reached = [float(gdna_fl_mass(p).sum()), float(rna_fl_mass(p).sum())]
        assert sorted(reached) == [0.0, 5.0], (
            f"pool {pool} reaches {reached}, expected exactly one component"
        )


def test_build_fl_large_pool_is_empirical():
    glob = _spike(100, 1.0e6)
    gdna = _spike(300, 1.0e7)  # huge gDNA pool, dominated by its own evidence
    fl = fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=gdna, max_size=1000
    )
    assert int(np.argmax(fl.gdna_pmf)) == 300
    np.testing.assert_allclose(fl.gdna_pmf.sum(), 1.0)
    np.testing.assert_allclose(fl.global_pmf.sum(), 1.0)


def test_build_fl_empty_pool_collapses_to_global():
    glob = _spike(100, 1.0e6)
    fl = fl_models_from_histograms(
        global_counts=glob,
        rna_counts=glob,
        gdna_counts=np.zeros(1001, dtype=np.float64),
        max_size=1000,
    )
    np.testing.assert_allclose(fl.gdna_pmf, fl.global_pmf)  # no gDNA evidence → global
    assert fl.n_gdna == 0.0


def test_build_fl_small_pool_shrinks_toward_global_no_cliff():
    glob = _spike(100, 1.0e6)  # global mode at 100
    gdna = _spike(300, 1.0)  # a single gDNA fragment at 300
    fl = fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=gdna, max_size=1000
    )
    # Smooth shrinkage (no threshold): with pool_total=1 ≪ the prior's effective size the global
    # anchor (bin 100) dominates the lone gDNA fragment (bin 300).
    assert fl.gdna_pmf[100] > fl.gdna_pmf[300]
    assert fl.gdna_pmf[300] > 0.0  # but the empirical fragment still nudges its bin


# ---------------------------------------------------------------------------
# Raw empirical counts + FragmentLengthModel accessors (QC report views)
# ---------------------------------------------------------------------------


def test_build_fl_stores_raw_aligned_counts():
    """Raw (unsmoothed) counts are stored aligned to max_size, totals = n_*."""
    rna = _spike(200, 300.0)
    gdna = _spike(300, 40.0)
    glob = _spike(100, 160.0) + rna + gdna  # every pool count is a deposited fragment
    fl = fl_models_from_histograms(
        global_counts=glob, rna_counts=rna, gdna_counts=gdna, max_size=1000
    )
    assert fl.global_counts.shape == (1001,)
    assert int(np.argmax(fl.rna_counts)) == 200
    assert int(np.argmax(fl.gdna_counts)) == 300
    # Raw counts are the empirical evidence, not the EB-smoothed pmf.
    assert fl.rna_counts[200] == 300.0
    assert fl.gdna_counts[300] == 40.0
    assert float(fl.global_counts.sum()) == 500.0
    assert fl.n_rna == 300.0
    assert fl.n_gdna == 40.0


def test_build_fl_counts_overflow_folds_into_last_bin():
    """Counts beyond max_size fold into the overflow bin when aligned."""
    over = np.zeros(1201, dtype=np.float64)
    over[1150] = 7.0  # beyond max_size=1000
    fl = fl_models_from_histograms(
        global_counts=_spike(100, 10.0),
        rna_counts=_spike(100, 10.0),
        gdna_counts=over,
        max_size=1000,
    )
    assert fl.gdna_counts.shape == (1001,)
    assert fl.gdna_counts[1000] == 7.0
    assert fl.n_gdna == 7.0


def test_accessors_return_empirical_models():
    """rna_model()/gdna_model()/global_model() expose the raw empirical FL."""
    fl = fl_models_from_histograms(
        global_counts=_spike(100, 500.0) + _spike(200, 300.0) + _spike(325, 40.0),
        rna_counts=_spike(200, 300.0),
        gdna_counts=_spike(325, 40.0),
        max_size=1000,
    )
    rna_m = fl.rna_model()
    gdna_m = fl.gdna_model()
    glob_m = fl.global_model()
    # Empirical modes/means track the raw counts (not smoothed toward global).
    assert rna_m.mode == 200
    assert gdna_m.mode == 325
    assert glob_m.mode == 100
    assert rna_m.mean == pytest.approx(200.0)
    # n_observations = pool total.
    assert rna_m.n_observations == 300
    assert gdna_m.n_observations == 40
    # the summary QC report reads them as they are
    from rigel.cli import length_summary

    assert length_summary(gdna_m)["mode"] == 325


# ── the two estimands of the gDNA length law, and the routing guarantee between them ─────────


def _models(**kw):
    rng = np.random.default_rng(5)
    g = rng.random(101)
    r = rng.random(101)
    return fl_models_from_histograms(
        global_counts=g + r, rna_counts=r, gdna_counts=g, max_size=100, **kw
    )


def _eb_weight(m) -> float:
    """The empirical-Bayes weight the build gave gDNA's own law, read off what the build did to its pmf, not restated:
    the shrinkage is ``gdna_pmf = w·law + (1 − w)·anchor`` at every length, asserted here, so ``w`` is one number the
    variance must carry too."""
    law = m.gdna_counts / m.gdna_counts.sum()
    gap = law - m.global_pmf
    at = int(np.argmax(np.abs(gap)))
    w = float((m.gdna_pmf[at] - m.global_pmf[at]) / gap[at])
    np.testing.assert_allclose(
        m.gdna_pmf, w * law + (1.0 - w) * m.global_pmf, rtol=1e-9, atol=1e-15
    )
    assert 0.0 < w < 1.0, w
    return w


# ── the field exists and is never None: the scorer must be able to read it unconditionally ───────


def test_gdnas_pmf_carries_its_raw_laws_variance_at_its_eb_weight_squared():
    """``Var(pmf_l) = (n / (n + prior_ess))² · law_l² / count_l``: a length's relative variance is one over its
    count, and the EB step shrinks the variance with the law."""
    m = _models()
    law = m.gdna_counts / m.gdna_counts.sum()
    weight = _eb_weight(m)
    np.testing.assert_allclose(
        m.gdna_pmf_var, weight * weight * law * law / m.gdna_counts, rtol=1e-12
    )


def test_a_gdna_length_with_no_count_carries_no_variance():
    counts = np.zeros(101)
    counts[[20, 40, 60]] = (5.0, 50.0, 500.0)
    m = fl_models_from_histograms(
        global_counts=counts + 1.0, rna_counts=counts, gdna_counts=counts, max_size=100
    )
    assert np.count_nonzero(m.gdna_pmf_var) == 3
    np.testing.assert_array_equal(m.gdna_pmf_var == 0.0, counts == 0.0)


def test_more_evidence_carries_proportionally_less_variance():
    """gDNA's raw law's variance is one over its counts, so ten times the counts is a tenth of it, read through each
    model's own EB weight."""
    from rigel.config import CONSTANTS

    rng = np.random.default_rng(3)
    counts = rng.random(101) * 100.0
    one = fl_models_from_histograms(
        global_counts=counts, rna_counts=counts, gdna_counts=counts, max_size=100
    )
    ten = fl_models_from_histograms(
        global_counts=counts, rna_counts=10.0 * counts, gdna_counts=10.0 * counts, max_size=100
    )

    def raw_var(m):
        # gDNA's law is its anchor's here, so no shrinkage shows: the weight is read from the code's own constant
        weight = m.n_gdna / (m.n_gdna + CONSTANTS.fragment_length.pool_prior_ess)
        return m.gdna_pmf_var / (weight * weight)

    np.testing.assert_allclose(raw_var(ten), raw_var(one) / 10.0, rtol=1e-12)
    # RNA's is its posterior's: more evidence, less variance at every length its law holds
    on = one.rna_pmf > 1e-4
    assert np.all(ten.rna_pmf_var[on] < one.rna_pmf_var[on])


def test_the_contrasts_own_variance_is_carried_at_its_eb_weight():
    rng = np.random.default_rng(4)
    supplied = rng.random(101) * 1e-6
    m = _models(gdna_var=supplied)
    weight = _eb_weight(m)
    np.testing.assert_allclose(m.gdna_pmf_var, weight * weight * supplied, rtol=1e-12)


def test_de_tilting_adds_no_evidence_to_the_rna_law():
    """The builder reads the spliced pool's RAW counts at their sj opportunity: de-tilting rescales a length and adds
    no evidence, so the law and its variance are the own-count law of the raw counts, never of the de-tilted ones."""
    from rigel.calibration.fl import own_count_law

    rng = np.random.default_rng(8)
    pools = np.zeros((N_FRAGMENT_POOLS, 101))
    pools[POOL_RNA_SPLICED] = np.round(rng.random(101) * 200.0)
    pools[POOL_DNA_INTERGENIC] = np.round(rng.random(101) * 50.0)
    pools[POOL_RNA_SPLICED, :5] = 0.0
    payload = SimpleNamespace(pool_lengths=pools, deposited_lengths=pools.sum(0), max_length=100)
    pi = np.linspace(0.0, 1.0, 101)
    m = _rna_law_models(payload, pi)
    law, var, _ = own_count_law(pools[POOL_RNA_SPLICED], pi, m.global_pmf)
    np.testing.assert_array_equal(m.rna_pmf, law)
    np.testing.assert_array_equal(m.rna_pmf_var, var)


def _rna_law_models(payload, pi):
    """The models' kernel on a payload's anchor and spliced pool at the sj opportunity ``pi``, as `build_fl_models`
    hands them to it (RNA's law reads nothing of gDNA's)."""
    pools = np.asarray(payload.pool_lengths, dtype=np.float64)
    return fl_models_from_histograms(
        global_counts=payload.deposited_lengths,
        rna_counts=pools[POOL_RNA_SPLICED],
        rna_selection=pi,
        gdna_counts=pools[POOL_DNA_INTERGENIC],
        max_size=int(payload.max_length),
    )


# ── RNA's law from its own counts, read as counts ─────────────────────────────────────────────────────────────

_L = np.arange(501, dtype=np.float64)


def _gauss(mu: float, sd: float, lo: int = 50) -> np.ndarray:
    p = np.exp(-0.5 * ((_L - mu) / sd) ** 2)
    p[:lo] = 0.0
    return p / p.sum()


def _rna_payload(n_rna: int, seed: int = 20261026, *, anchor_from: int = 50):
    """A 98 %-gDNA library: the all-fragment law is mostly gDNA's N(250, 150), RNA's own is N(100, 50), and the
    spliced pool is RNA's law drawn through a junction-crossing chance that rises with length and is zero below 30 bp."""
    rng = np.random.default_rng(seed)
    truth = _gauss(100.0, 50.0)
    anchor = 0.98 * _gauss(250.0, 150.0, anchor_from) + 0.02 * truth
    pi = np.clip((_L - 30.0) / 400.0, 0.0, 1.0)
    pools = np.zeros((N_FRAGMENT_POOLS, 501))
    if n_rna:
        pools[POOL_RNA_SPLICED] = rng.multinomial(n_rna, (truth * pi) / (truth * pi).sum())
    pools[POOL_DNA_INTERGENIC] = np.round(anchor * 1e6)
    payload = SimpleNamespace(
        pool_lengths=pools, deposited_lengths=np.round(anchor * 5e6), max_length=500
    )
    return payload, pi, truth


def _h2(a: np.ndarray, b: np.ndarray) -> float:
    return float(((np.sqrt(a / a.sum()) - np.sqrt(b / b.sum())) ** 2).sum())


def test_a_small_rna_pool_reads_its_own_law_not_the_all_fragment_one():
    """2,000 spliced fragments in a 98 %-gDNA library: RNA's law is its own counts', not gDNA's shape at a weight no
    measurement supplied."""
    payload, pi, truth = _rna_payload(2_000)
    m = _rna_law_models(payload, pi)
    assert abs(m.rna_pmf @ _L - truth @ _L) < 3.0
    assert _h2(m.rna_pmf, truth) < 0.02


def test_an_empty_rna_pool_is_the_all_fragment_law_exactly():
    payload, pi, _ = _rna_payload(0)
    m = _rna_law_models(payload, pi)
    np.testing.assert_array_equal(m.rna_pmf, m.global_pmf)


def test_a_deep_rna_pool_is_its_own_law():
    payload, pi, truth = _rna_payload(2_000_000)
    m = _rna_law_models(payload, pi)
    assert _h2(m.rna_pmf, truth) < 1e-3


def test_a_length_the_pool_could_reach_and_did_not_is_evidence():
    """Past RNA's support the pool reached every length (its crossing chance is high there) and drew nothing, so the law
    is low there. The all-fragment law's long tail is gDNA's and says nothing about RNA."""
    payload, pi, truth = _rna_payload(2_000)
    m = _rna_law_models(payload, pi)
    assert m.rna_pmf[300:].sum() < 10.0 * truth[300:].sum() + 1e-3


def test_where_the_pool_cannot_reach_the_law_keeps_the_anchors_shape():
    """A length no spliced fragment can have (crossing chance zero) carries no evidence: the law there is the
    all-fragment law's shape, held at one ratio."""
    payload, pi, _ = _rna_payload(2_000, anchor_from=10)
    m = _rna_law_models(payload, pi)
    dark = (pi == 0.0) & (m.global_pmf > 0.0)
    assert dark.sum() >= 10
    ratio = m.rna_pmf[dark] / m.global_pmf[dark]
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-9)


def test_a_one_length_law_is_certain():
    """Normalised, a law of one length holds probability one whatever its count: no variance. The deviation's level,
    which only the pool's total fixes, is not a property of the law and must not reach its variance."""
    from rigel.calibration.fl import own_count_law

    law, var, _ = own_count_law(np.array([100.0]), np.array([1.0]), np.array([1.0]))
    np.testing.assert_array_equal(law, [1.0])
    np.testing.assert_array_equal(var, [0.0])


def test_a_two_length_laws_two_variances_are_equal():
    """``f_1 + f_2 = 1``, so ``Var f_1 = Var f_2``: the normalised law's variances, not the deviation's."""
    from rigel.calibration.fl import own_count_law

    _law, var, _ = own_count_law(np.array([30.0, 70.0]), None, np.array([0.5, 0.5]))
    assert var[0] > 0.0
    assert var[0] == pytest.approx(var[1], rel=1e-12)


def test_the_rna_laws_variance_is_its_normalised_delta_method():
    """Brute force: the posterior's Hessian rebuilt densely at the returned mode, ``H = diag(mu) + K / step`` with ``K``
    the path's first-difference precision and ``mu`` the fitted Poisson means (at the mode they sum to the pool's
    total, since the walk's penalty has no level), its inverse pushed through the normalisation's Jacobian
    ``J = diag(f) - f f'``: the reported variance is ``diag(J H^-1 J')`` at every length, over a window holding an
    unreachable length and an empty one, and zero outside the window."""
    from rigel.calibration.fl import own_count_law

    raw = np.array([0, 0, 4, 9, 0, 15, 22, 18, 11, 6, 2, 1, 0, 0], dtype=np.float64)
    reach = np.array([1, 1, 0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1, 1, 1, 1], dtype=np.float64)
    reach[2] = 0.0
    raw[2] = 0.0
    anchor = np.array([0, 0.02, 0.05, 0.08, 0.1, 0.12, 0.14, 0.13, 0.11, 0.09, 0.07, 0.05, 0.04, 0])
    law, var, step = own_count_law(raw, reach, anchor)
    lo, hi = 1, 13
    f = law[lo:hi]
    live = (reach[lo:hi] > 0.0) & (anchor[lo:hi] > 0.0)
    exposure = np.where(live, reach[lo:hi] * f, 0.0)
    mu = raw.sum() * exposure / exposure.sum()
    n = hi - lo
    walk = 2.0 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    walk[0, 0] = walk[-1, -1] = 1.0
    covariance = np.linalg.inv(np.diag(mu) + walk / step)
    jacobian = np.diag(f) - np.outer(f, f)
    np.testing.assert_allclose(var[lo:hi], np.diag(jacobian @ covariance @ jacobian.T), rtol=1e-6)
    assert np.all(var[:lo] == 0.0) and np.all(var[hi:] == 0.0)


def test_the_variance_solve_holds_a_block_of_lengths_not_the_support_squared(monkeypatch):
    """The variance's right-hand side is dense, one column per length: solved whole, a 10,000-length support (a
    configured ``max_frag_length``) would hold 800 MB. Each solve holds at most ``law_variance_block`` columns, and
    the block moves no number: the lengths' solves are independent, so one length at a time, a few, and the default
    give the same bits."""
    import scipy.linalg

    import rigel.calibration.fl as fl
    from _constants import replaced
    from rigel.config import CONSTANTS

    widths = []
    solve = scipy.linalg.solve_banded

    def recording(l_and_u, ab, b, *args, **kwargs):
        widths.append(np.shape(b)[1] if np.ndim(b) == 2 else 1)
        return solve(l_and_u, ab, b, *args, **kwargs)

    monkeypatch.setattr(scipy.linalg, "solve_banded", recording)
    anchor = _gauss(250.0, 80.0)
    raw = np.random.default_rng(3).poisson(2_000 * anchor).astype(np.float64)
    _law, var, _ = fl.own_count_law(raw, None, anchor)
    block = CONSTANTS.fragment_length.law_variance_block
    assert np.count_nonzero(anchor > 0.0) > block
    assert widths and max(widths) <= block
    for other in (1, 7):
        monkeypatch.setattr(
            fl, "CONSTANTS", replaced(CONSTANTS, "fragment_length.law_variance_block", other)
        )
        np.testing.assert_array_equal(fl.own_count_law(raw, None, anchor)[1], var)


def test_the_rna_laws_variance_predicts_what_replicate_libraries_show():
    """Repeated libraries at one law: the Hellinger noise the variance predicts, ``Σ Var f / 4 f``, against the replicate
    mean of ``H²(f̂, f)`` — what the length decision reads it for. At an exact null (RNA's law the anchor's, so the
    deviation is flat) the deviation's level variance would be most of the prediction; away from it the posterior
    stays on the safe side of what replicates show."""
    from rigel.calibration.fl import own_count_law

    pi = np.clip((_L - 30.0) / 400.0, 0.0, 1.0)
    rng = np.random.default_rng(20261027)
    for truth, anchor, bounds in (
        (_gauss(200.0, 60.0), _gauss(200.0, 60.0), (0.5, 2.0)),
        (_gauss(100.0, 50.0), 0.98 * _gauss(250.0, 150.0) + 0.02 * _gauss(100.0, 50.0), (0.8, 2.5)),
    ):
        draw = truth * pi / (truth * pi).sum()
        for n in (300, 10_000):
            seen, predicted = [], []
            for _ in range(40):
                law, var, _ = own_count_law(rng.multinomial(n, draw).astype(np.float64), pi, anchor)
                on = law > 0.0
                seen.append(float(((np.sqrt(law) - np.sqrt(truth)) ** 2).sum()))
                predicted.append(float((var[on] / (4.0 * law[on])).sum()))
            ratio = np.mean(predicted) / np.mean(seen)
            assert bounds[0] < ratio < bounds[1], (n, ratio)


def test_the_crossing_chance_is_an_exposure_up_to_a_constant():
    payload, pi, _ = _rna_payload(2_000)
    a = _rna_law_models(payload, pi)
    b = _rna_law_models(payload, 0.5 * pi)
    np.testing.assert_allclose(a.rna_pmf, b.rna_pmf, rtol=1e-6, atol=1e-12)


def test_the_smoothing_is_read_off_the_data():
    """The random walk's step is the evidence maximum, inside its search range on an ordinary pool: neither the flat
    deviation (the anchor's shape) nor no smoothing at all."""
    from rigel.calibration.fl import own_count_law
    from rigel.config import CONSTANTS

    payload, pi, _ = _rna_payload(2_000)
    anchor = np.asarray(payload.deposited_lengths, dtype=np.float64)
    _pmf, _var, step = own_count_law(
        np.asarray(payload.pool_lengths)[POOL_RNA_SPLICED], pi, anchor / anchor.sum()
    )
    fl = CONSTANTS.fragment_length
    assert fl.law_step_min * 10.0 < step < fl.law_step_max / 10.0


def test_a_pool_count_the_anchor_does_not_hold_is_refused():
    """The anchor is every deposited fragment and the pool is drawn from it: a pool count at a length the anchor has
    none at means the two came from different payloads, and the law refuses rather than drop it."""
    from rigel.calibration.fl import own_count_law

    anchor = _spike(100, 1.0)
    with pytest.raises(ValueError, match="different payloads"):
        own_count_law(_spike(200, 5.0), None, anchor)


def test_the_contrasts_variance_falls_as_its_evidence_grows():
    """The two-pool contrast's per-length variance runs through its own coefficients: ten times the contained
    counts, at the same purities and the same shapes, is a tenth of the variance."""
    import dataclasses

    from _fl_fixture import build_fixture

    from rigel.calibration.fl import deconvolved_gdna_counts

    payload, opp, regions = build_fixture()
    counts, var = deconvolved_gdna_counts(payload, opp, regions)
    assert counts is not None
    assert np.all(np.isfinite(var)) and np.all(var >= 0.0) and var.sum() > 0.0
    scaled = dataclasses.replace(
        payload,
        pool_lengths=np.asarray(payload.pool_lengths, dtype=np.float64) * 10.0,
        region_contained_count=np.asarray(payload.region_contained_count, dtype=np.float64) * 10.0,
    )
    counts10, var10 = deconvolved_gdna_counts(scaled, opp, regions)
    assert counts10 is not None
    np.testing.assert_allclose(
        counts10 / counts10.sum(), counts / counts.sum(), rtol=1e-9, atol=1e-15
    )
    np.testing.assert_allclose(var10, var / 10.0, rtol=1e-9, atol=1e-30)


def _clipping_contrast_payload(scale: float = 1.0):
    """The fixture with RNA in both contained pools at different shares and with different laws, so the two-pool
    inversion's premise (one contaminant shape) fails and it runs negative, and clips, at the intron's RNA lengths:
    intergenic region 0 carries RNA, region 4 is the clean object the one-sided rate reads, and the intron carries more
    and longer RNA."""
    import dataclasses

    from _fl_fixture import MAX_SIZE, RHO, build_fixture

    payload, opp, regions = build_fixture()
    L = np.arange(MAX_SIZE + 1, dtype=np.float64)
    g_law = np.exp(-0.5 * ((L - 45.0) / 10.0) ** 2)
    r_ig = np.exp(-0.5 * ((L - 55.0) / 10.0) ** 2)
    r_in = np.exp(-0.5 * ((L - 80.0) / 10.0) ** 2)
    g_law[0] = r_ig[0] = r_in[0] = 0.0
    g_law, r_ig, r_in = g_law / g_law.sum(), r_ig / r_ig.sum(), r_in / r_in.sum()
    counts = np.array(payload.region_contained_count, dtype=np.float64, copy=True) * 400.0
    gdna_ig = float(counts[0].sum() + counts[4].sum())
    gdna_in = float(counts[2].sum())
    rna_ig, rna_in = 0.25 * gdna_ig, 1.0 * gdna_in
    counts[0] += rna_ig / 2.0
    counts[2] += rna_in / 2.0
    pools = np.array(payload.pool_lengths, dtype=np.float64, copy=True)
    pools[0] = g_law * gdna_ig + r_ig * rna_ig
    pools[1] = g_law * gdna_in + r_in * rna_in
    expected = dataclasses.replace(
        payload, pool_lengths=pools * scale, region_contained_count=counts * scale
    )
    return expected, opp, regions, RHO


def test_the_contrasts_variance_is_its_own_estimators_through_the_clip():
    """Redrawn many times from the pools' own Poisson counts and pushed through the contrast, the law's per-length
    variance matches the plug-in the contrast reports: the inverse's coefficient where the inversion is open, the
    mixture's alone where the clip holds it at zero. Charging the inverse's coefficient through the clip would overstate
    a clipped length's variance many times over.
    """
    import dataclasses

    from rigel.calibration.fl import contained_purities, deconvolved_gdna_counts

    payload, opp, regions, _rho = _clipping_contrast_payload()
    counts, var = deconvolved_gdna_counts(payload, opp, regions)
    _rate, (a0, a1), _n = contained_purities(payload, regions)
    assert counts is not None and a0 > a1
    law = counts / counts.sum()
    pools = np.asarray(payload.pool_lengths, dtype=np.float64)
    inverse = (1.0 - a1) * pools[0] / pools[0].sum() - (1.0 - a0) * pools[1] / pools[1].sum()
    open_ = (inverse > 0.0) & (law > 1e-3)
    clipped = (inverse < -1e-3) & (pools[1] > 0.0)  # well inside the clip, away from the kink
    assert clipped.sum() >= 3 and open_.sum() >= 5, (
        "the fixture must clip somewhere and stay open elsewhere"
    )

    rng = np.random.default_rng(17)
    draws = []
    for _ in range(400):
        drawn = pools.copy()
        drawn[:2] = rng.poisson(pools[:2]).astype(np.float64)
        c, _v = deconvolved_gdna_counts(
            dataclasses.replace(payload, pool_lengths=drawn), opp, regions
        )
        assert c is not None
        draws.append(c / c.sum())
    mc = np.var(np.array(draws), axis=0, ddof=1)
    ratio = mc[open_] / var[open_]
    assert np.all((ratio > 0.6) & (ratio < 1.6)), (
        f"Monte Carlo over plug-in where open: {np.round(ratio, 2)}"
    )
    # through the clip only the mixture's small share is left, and the plug-in says so: both are at the law's own
    # scale there, many orders below what the inverse's coefficient would charge
    assert np.all(var[clipped] <= 4.0 * mc[clipped] + 1e-15) and np.all(
        mc[clipped] <= 4.0 * var[clipped] + 1e-15
    )


# ── the RNA law the census reads: a law, whatever its scale, and the one FLModels carries ───────────


def _normal_law(L, mu) -> np.ndarray:
    """A normal length law of sd 30 centred on *mu*, over the lengths *L*, with no mass at 0."""
    p = np.exp(-0.5 * ((L - mu) / 30.0) ** 2)
    p[0] = 0.0
    return p / p.sum()


def _tie_fixture(rna_in_introns: float, rna_in_intergenic: float = 0.0):
    """Four intergenic and four intronic 20-kb regions of pure gDNA at the uniform expectation, and four probed exons:
    the contained pools carry the uniform-frame law ``g`` (mean 150), the two crossing pools a capture-selected law
    (mean 200) at ten times their mass. Every off-target region holds exactly its expected count, so the one-sided rate
    lands just above it and both purities clip at 1, the exact tie of a near-pure captured library. ``rna_in_introns``
    adds RNA (mean 120) to the introns per unit of opportunity; enough of it pulls the intronic purity below 1 and
    resolves the tie. The same amount in both pools (``rna_in_intergenic``, its own RNA law, mean 100) keeps the
    purities tied below 1 with the two pools' shapes apart. Returns ``(payload, opportunity, regions, L)``."""
    from _fl_fixture import regions_of

    from rigel.calibration.effective_length import contained_eff_length

    n = 400
    L = np.arange(n + 1, dtype=np.float64)
    law = partial(_normal_law, L)
    g, g_captured, r, r_ig = law(150.0), law(200.0), law(120.0), law(100.0)
    lengths = np.full(12, 20_000.0)
    lengths[8:] = 300.0
    types = np.array([0] * 4 + [1] * 4 + [2] * 4, dtype=np.uint8)  # intergenic / intron / exon
    e = contained_eff_length(lengths, g)
    counts = np.zeros((12, 2))
    counts[:8, 0] = 0.05 * e[:8]
    # the RNA sits in two regions of a pool at twice the pool rate, so the one-sided rate still reads
    # the gDNA level off the other two and the pool's purity falls below 1
    counts[4:6, 0] += 2.0 * rna_in_introns * e[4:6]
    counts[0:2, 0] += 2.0 * rna_in_intergenic * e[0:2]
    counts[8:, 0] = 1e4
    a_ig = np.clip(lengths[:4].sum() - L + 1.0, 0.0, None)
    a_in = np.clip(lengths[4:8].sum() - L + 1.0, 0.0, None)
    a_x = 8.0 * np.clip(L - 1.0, 0.0, None)
    total = a_ig + a_in + 2.0 * a_x
    rna = 2.0 * rna_in_introns * float(e[4:6].sum())
    rna_ig = 2.0 * rna_in_intergenic * float(e[0:2].sum())
    pools = np.zeros((N_FRAGMENT_POOLS, n + 1))
    pools[POOL_DNA_INTERGENIC] = g * a_ig / (g * a_ig).sum() * (counts[:4].sum() - rna_ig)
    if rna_ig > 0.0:
        pools[POOL_DNA_INTERGENIC] += r_ig * a_ig / (r_ig * a_ig).sum() * rna_ig
    pools[POOL_DNA_INTRONIC] = g * a_in / (g * a_in).sum() * (counts[4:8].sum() - rna)
    if rna > 0.0:
        pools[POOL_DNA_INTRONIC] += r * a_in / (r * a_in).sum() * rna
    pools[POOL_DNA_INTRON_EXON] = pools[POOL_DNA_INTERGENIC_EXON] = (
        g_captured * a_x / (g_captured * a_x).sum() * 5e5
    )
    pools[POOL_RNA_SPLICED] = r * 1e5
    payload = SimpleNamespace(
        pool_lengths=pools,
        region_contained_count=counts,
        max_length=n,
        deposited_lengths=pools.sum(axis=0),
    )
    opportunity = SimpleNamespace(
        pools=(a_ig, a_in, a_x, a_x),
        total=total,
        combined_probability=lambda: (a_ig + a_in + 2.0 * a_x) / total,
    )
    return payload, opportunity, regions_of(lengths, types), L


def _mean_length(p, L) -> float:
    p = np.asarray(p, dtype=np.float64)
    return float((L[: p.size] * p).sum() / p.sum())


def test_a_purity_tie_returns_the_contained_pairs_own_mixture():
    """At an exact purity tie (``a_0 = a_1``) the contrast's resolution weight is 0, so its answer is the contained
    pair's own de-tilted mixture, the limit its fade approaches and what the boundary pair returns at its own tie. It
    does not decline: declining would hand the uniform-frame law to the capture-selected four-pool census, which reads
    long.
    """
    from rigel.calibration.fl import contained_purities, deconvolved_gdna_counts

    payload, opp, regions, L = _tie_fixture(0.0)
    with np.errstate(divide="raise", invalid="raise"):  # a tie never divides by its zero separation
        counts, _var = deconvolved_gdna_counts(payload, opp, regions)
    _rate, (a0, a1), _n = contained_purities(payload, regions)
    assert a0 == a1 == 1.0
    assert counts is not None
    assert _mean_length(counts, L) == pytest.approx(150.0, abs=0.5)


def test_the_uniform_law_is_continuous_through_a_purity_tie():
    """No cliff at the purity axis, end to end through ``build_fl_models``: the uniform-frame law just past the tie
    (intronic purity below 1) and at the tie must agree, and both must read the uniform-frame law, not the
    capture-selected census 50 bp longer."""
    from rigel.calibration.fl import contained_purities, build_fl_models

    means, seps = [], []
    for rna in (0.0005, 0.001):
        payload, opp, regions, L = _tie_fixture(rna)
        _rate, (a0, a1), _n = contained_purities(payload, regions)
        seps.append(a0 - a1)
        m = build_fl_models(
            payload, sj_opportunity=np.ones(L.size), gdna_opportunity=opp, regions=regions
        )
        means.append(_mean_length(m.gdna_pmf, L))
    assert seps[0] == 0.0 < seps[1]  # the fixture straddles the tie
    assert abs(means[0] - means[1]) < 0.5, (
        f"{means[0]:.1f} at the tie against {means[1]:.1f} just past it"
    )
    assert means[0] == pytest.approx(150.0, abs=1.0)


def test_a_tie_below_full_purity_is_the_pairs_mixture_too():
    """Two pools of equal purity carry no contrast whatever that purity is, so a tie below 1 also answers with the
    pair's de-tilted mixture, contaminant included. The ``errstate`` guards hold that a tie never divides by its zero
    separation: at ``sep = 0`` the inversion's ``0/0`` would only fall through to the mixture by accident.
    """
    from rigel.calibration.fl import contained_purities, deconvolved_gdna_counts
    from rigel.calibration.sj_opportunity import detilt_pool

    payload, opp, regions, L = _tie_fixture(0.005, 0.005)
    with np.errstate(divide="raise", invalid="raise"):
        counts, _var = deconvolved_gdna_counts(payload, opp, regions)
    _rate, (a0, a1), _n = contained_purities(payload, regions)
    assert a0 - a1 == 0.0 and a0 < 1.0
    assert counts is not None and np.isfinite(counts).all()
    pools = np.asarray(payload.pool_lengths)
    f = [
        detilt_pool(pools[p], np.asarray(opp.pools[i]) / opp.total)
        for i, p in enumerate((POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC))
    ]
    n = [float(payload.region_contained_count[m].sum()) for m in (slice(0, 4), slice(4, 8))]
    mixture = n[0] * f[0] / f[0].sum() + n[1] * f[1] / f[1].sum()
    np.testing.assert_allclose(counts / counts.sum(), mixture / mixture.sum(), atol=1e-12)


# ── the convergence law: no cliffs, and the two estimands merge when the split is unresolvable ───


def test_resolution_weight_is_a_smooth_signal_to_noise_ratio():
    """`lam = S/(S+N)`: 0 when the split is pure noise, 1 when noise vanishes, ½ at S = N, and monotone in between. No
    threshold, so behaviour never jumps."""
    from rigel.calibration.fl import resolution_weight

    assert resolution_weight(0.0, 1.0) == 0.0
    assert resolution_weight(1.0, 0.0) == 1.0
    assert resolution_weight(1.0, 1.0) == pytest.approx(0.5)
    prev = -1.0
    for s in np.linspace(0.0, 10.0, 60):
        lam = resolution_weight(float(s), 1.0)
        assert 0.0 <= lam <= 1.0 and lam >= prev
        prev = lam


# ── the length law over observable bases ───────────────────────────────────────────────────────────

#: The masked genome's unknown runs: the first intergenic region ([0, 20 000)) N end to end, and the middle
#: half of the second ([25 000, 35 000)).
_UNKNOWN = [(0, 20_000), (25_000, 35_000)]


#: The observable pieces of the masked genome's eight off-target regions, written out by hand from ``_UNKNOWN``:
#: the all-N region has none, the half-N one two 5-kb pieces, the rest are whole.
_PIECES = [[], [5_000, 5_000]] + [[20_000]] * 6


def _masked_fixture(rna_intergenic: float, rna_intronic: float = 0.05, mask=_UNKNOWN):
    """`_tie_fixture`'s genome (four intergenic and four intronic 20-kb regions and four probed exons) with the unknown
    sequence a real assembly has. No fragment is drawn on N, so every count is per observable start: gDNA (mean 150) at
    ``RHO``, RNA of one law (mean 120) at *rna_intronic* in two introns and at *rna_intergenic* in two intergenic
    regions, so the contrast's premise (one contaminant law in both pools) holds exactly. The opportunity is the frame's
    own (`observable_gdna_opportunity`), so each contained pool is its exact expectation. Returns ``(payload,
    opportunity, regions, L)``."""
    from _fl_fixture import regions_of

    from rigel.calibration.effective_length import contained_eff_length
    from rigel.calibration.gdna_opportunity import observable_gdna_opportunity

    rho = 0.05
    n = 400
    L = np.arange(n + 1, dtype=np.float64)
    law = partial(_normal_law, L)
    g, g_captured, r = law(150.0), law(200.0), law(120.0)
    lengths = np.full(12, 20_000)
    lengths[8:] = 300
    types = np.array([0] * 4 + [1] * 4 + [2] * 4)
    regions = regions_of(lengths, types, mask)

    def support(pmf):
        return np.array([contained_eff_length(np.array(ps, float), pmf).sum() for ps in _PIECES])

    e_g, e_r = support(g), support(r)
    counts = np.zeros((12, 2))
    counts[:8, 0] = rho * e_g
    counts[2:4, 0] += rna_intergenic * e_r[2:4]
    counts[4:6, 0] += rna_intronic * e_r[4:6]
    counts[8:, 0] = 1e4
    opportunity = observable_gdna_opportunity(regions, n)
    whole = np.clip(20_000.0 - L + 1.0, 0.0, None)  # the RNA regions are unmasked 20-kb regions
    pools = np.zeros((N_FRAGMENT_POOLS, n + 1))
    pools[POOL_DNA_INTERGENIC] = rho * g * opportunity.intergenic_contained
    pools[POOL_DNA_INTERGENIC] += 2.0 * rna_intergenic * r * whole
    pools[POOL_DNA_INTRONIC] = rho * g * opportunity.intronic_contained
    pools[POOL_DNA_INTRONIC] += 2.0 * rna_intronic * r * whole
    pools[POOL_DNA_INTRON_EXON] = pools[POOL_DNA_INTERGENIC_EXON] = g_captured * 1e3
    pools[POOL_RNA_SPLICED] = r * 1e5
    payload = SimpleNamespace(
        pool_lengths=pools,
        region_contained_count=counts,
        max_length=n,
        deposited_lengths=pools.sum(axis=0),
    )
    # the pools' own arithmetic: every contained count is in exactly one pool
    assert pools[POOL_DNA_INTERGENIC].sum() == pytest.approx(counts[:4, 0].sum(), rel=1e-12)
    assert pools[POOL_DNA_INTRONIC].sum() == pytest.approx(counts[4:8, 0].sum(), rel=1e-12)
    return payload, opportunity, regions, L


def test_the_one_sided_rate_and_the_purities_read_observable_support():
    """A region of N carries no fragment. Read over its whole length it would be the cleanest region of all (a zero
    count against 20 kb of exposure) and pull the one-sided rate to a hundredth of the truth; read over its pieces it is
    no object at all, and a half-N region is the density it holds. So the rate is the estimator on each region's
    observable support (written out by hand here), and each pool's gDNA share reads the same support; the half-N region
    read whole would overstate the intergenic share by a fifth. On exact-expectation counts the estimator itself reads a
    few percent above ``RHO``, because its dispersion term expects scatter that exact counts lack; that is independent
    of the genome and not what this test is about.
    """
    from rigel.calibration.fl import (
        contained_purities,
        deconvolved_gdna_counts,
    )
    from rigel.calibration.effective_length import contained_eff_length
    from rigel.calibration.gdna_density import one_sided_rate

    payload, opportunity, regions, _L = _masked_fixture(rna_intergenic=0.005)
    counts_law, _var = deconvolved_gdna_counts(payload, opportunity, regions)
    got_rate, (a_ig, a_in), _n = contained_purities(payload, regions)
    pooled = np.asarray(payload.pool_lengths, float)[[POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC]].sum(
        axis=0
    )
    pooled /= pooled.sum()
    support = np.array([contained_eff_length(np.array(ps, float), pooled).sum() for ps in _PIECES])
    counts = payload.region_contained_count[:8].sum(axis=1)
    live = support > 0.0  # the all-N region is no object
    rate = one_sided_rate(counts[live], support[live]).rate
    assert counts_law is not None
    assert got_rate == pytest.approx(rate, rel=1e-12)
    for share, pool in ((a_ig, slice(0, 4)), (a_in, slice(4, 8))):
        assert share == pytest.approx(
            min(rate * support[pool].sum() / counts[pool].sum(), 1.0), rel=1e-12
        )
    assert 0.9 < a_ig < 1.0  # RNA in the intergenic pool: the share does not clip


def test_over_observable_bases_the_contrast_recovers_the_uniform_law():
    """End to end through the contrast: with the intergenic pool pure, its gDNA share clips at 1 and the contrast is
    that pool's de-tilted histogram, the uniform-frame law (mean 150 bp), moved toward the pools' mixture only by the
    resolution weight, which resolves the shares' separation almost whole. With the rate read over whole lengths, both
    shares would fall near zero, their separation drown in its error, the resolution weight go to zero, and the law
    become the pools' mixture, pulled short by the introns' RNA (mean 120).
    """
    from rigel.calibration.fl import contained_purities, deconvolved_gdna_counts

    payload, opportunity, regions, L = _masked_fixture(rna_intergenic=0.0)
    counts, _var = deconvolved_gdna_counts(payload, opportunity, regions)
    _rate, (a_ig, _a_in), _n = contained_purities(payload, regions)
    assert counts is not None
    assert a_ig == 1.0
    assert _mean_length(counts, L) == pytest.approx(150.0, abs=0.1)


def test_the_length_law_reads_the_index_mask(tmp_path, monkeypatch):
    """`pipeline.library_fl_models` hands `build_fl_models` the index's frame, carrying its unobservable mask, and
    the gDNA opportunity over that observable genome: an index whose genome holds an N run gets both, and the
    pools' total ``T(w)`` counts no start on N."""
    import rigel.calibration.fl as fl
    from _index_builder import MINI_GTF, write_fasta

    from rigel.calibration.gdna_opportunity import observable_gdna_opportunity, total_opportunity
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.index import TranscriptIndex
    from rigel.pipeline import library_fl_models

    fasta = write_fasta(tmp_path / "g.fa", {"chr1": 2000}, {"chr1": [(1500, 1800)]})
    gtf = tmp_path / "a.gtf"
    gtf.write_text(MINI_GTF)
    TranscriptIndex.build(fasta, gtf, tmp_path / "idx", write_tsv=False)
    index = TranscriptIndex.load(tmp_path / "idx")
    seen = {}
    monkeypatch.setattr(fl, "build_fl_models", lambda payload, **kw: seen.update(kw))
    library_fl_models(SimpleNamespace(max_length=200), index)

    assert seen["regions"].unobservable.tolist() == [[0, 1500, 1800]]
    want = observable_gdna_opportunity(RegionArrays.from_index(index), 200)
    got = seen["gdna_opportunity"]
    for got_pool, want_pool in zip(got.pools + (got.total,), want.pools + (want.total,)):
        np.testing.assert_array_equal(got_pool, want_pool)
    assert np.all(
        got.total < total_opportunity(np.array([2000]), 200)
    )  # the whole reference's starts


def test_a_frame_of_another_partition_is_refused():
    """The contrast reads each region's count against that region's support, so the frame must be the partition the
    payload was deposited into: one with another region count is refused, never truncated."""
    from _fl_fixture import regions_of

    from rigel.calibration.fl import deconvolved_gdna_counts

    payload, opportunity, _regions, _L = _tie_fixture(0.0)  # twelve regions
    other = regions_of([20_000] * 11, [0] * 11)
    with pytest.raises(ValueError, match="different partitions"):
        deconvolved_gdna_counts(payload, opportunity, other)
