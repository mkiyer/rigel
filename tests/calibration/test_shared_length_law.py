"""``shared_length_law``: one length table for both origins unless the data show a difference capture cannot make.

The gates, in the order the design argues them. On equal chemistry the scorer's two tables are one array, however
strong the capture, because a captured census of gDNA's law sits inside the bound; a real gap is kept whole, never
partly shrunk; capture only lengthens, so an RNA census no longer than gDNA's law moves nothing; the mixture stays
within RNA's captured fraction and more capture never widens what is left; the bound is the worst case over every
capture response that rewards length at most proportionally, and only over those; the weight and the mixture are
continuous in every input; with no capture each origin keeps its own law; an RNA law no spliced fragment measured is
one law; and the geometric path is exact at both ends. On laws BUILT from nested pools and handed over as the pipeline
hands them, equal chemistry under capture is one law at every RNA pool size, a gap the spliced fragments measure is
kept whole, and RNA with no spliced fragment is one law.

PERTURBATION (`TRAPS: perturb-every-gate`): the no-measurement branch deleted, that branch keeping gDNA's own law, and
the pipeline reading every RNA law as measured each fire at least one gate.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from rigel.calibration.shared_length_law import (
    capture_bound,
    capture_mixture,
    geometric_path,
    hellinger_sq,
    sampling_noise,
    shared_length_law,
    size_biased,
    trust,
)

L = np.arange(1001, dtype=np.float64)
ZERO = np.zeros(L.size)


def _law(mean: float, sd: float) -> np.ndarray:
    p = stats.gamma.pdf(L, a=(mean / sd) ** 2, scale=sd**2 / mean)
    return p / p.sum()


def _captured(g: np.ndarray, kappa: float, m: np.ndarray) -> np.ndarray:
    """``g``'s census under a capture response ``m(l)`` with captured fraction ``kappa``: ``(1 − kappa) g + kappa g m /
    <m>``, the response's weight ``1 + beta m(l)`` written at the fraction it captures."""
    out = (1.0 - kappa) * g + kappa * g * m / float(g @ m)
    return out / out.sum()


G = _law(225.0, 82.0)


def _tables(r, g, kappa, var_r=ZERO, var_g=ZERO):
    return shared_length_law(r, g, var_r, var_g, kappa, rna_measured=True)


# ── one law where capture explains the difference ──────────────────────────────────────────────────────


@pytest.mark.parametrize("kappa", [0.2, 0.6, 0.98])
def test_equal_chemistry_under_capture_is_one_table(kappa):
    """RNA's census is gDNA's law captured at every fraction up to ``kappa``: the gDNA table IS the RNA table."""
    for fraction in (0.0, 0.5 * kappa, kappa):
        r = size_biased(G, fraction)
        t = _tables(r, G, kappa)
        assert t.decision.weight == 0.0
        np.testing.assert_array_equal(t.gdna, r)
        assert t.rna is not t.gdna


def test_a_saturating_capture_response_is_one_table_too():
    """Overlap with a probe saturates at the probe's length: a response inside the bound's premise."""
    for probe in (60.0, 150.0, 400.0):
        r = _captured(G, 0.9, np.minimum(L, probe))
        assert _tables(r, G, 0.9).decision.weight == 0.0


def test_sampling_noise_alone_is_one_table_under_capture():
    """Two independent draws of one law, each with its plug-in variance: the noise term absorbs the difference."""
    rng = np.random.default_rng(11)
    for n_g, n_r in ((2e4, 1e5), (2e5, 1e6)):
        for _ in range(20):
            cg = rng.multinomial(int(n_g), G).astype(np.float64)
            cr = rng.multinomial(int(n_r), size_biased(G, 0.5)).astype(np.float64)
            g, r = cg / cg.sum(), cr / cr.sum()
            var_g = np.divide(g * g, cg, out=np.zeros_like(g), where=cg > 0)
            var_r = np.divide(r * r, cr, out=np.zeros_like(r), where=cr > 0)
            assert _tables(r, g, 0.9, var_r, var_g).decision.weight == 0.0


# ── a real gap is kept whole ───────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "g_mean, r_mean", [(130.0, 260.0), (260.0, 130.0), (200.0, 300.0), (300.0, 200.0)]
)
def test_a_real_gap_is_never_partly_shrunk(g_mean, r_mean):
    """A gap far beyond what capture can make takes the tilted law whole: ``w = 1`` exactly, the table ``G_c``."""
    g, r = _law(g_mean, 70.0), _law(r_mean, 80.0)
    t = _tables(r, g, 0.9)
    assert t.decision.weight == 1.0
    tilted, _ = capture_mixture(g, r, 0.9)
    np.testing.assert_array_equal(t.gdna, tilted)
    np.testing.assert_array_equal(t.rna, r)


def test_the_weight_is_the_guaranteed_real_part_of_the_difference():
    """``w = D/B − 1``: 0 up to the bound, 1 from twice the bound, linear in the distance between."""
    assert trust(0.0, 0.04) == 0.0
    assert trust(0.04, 0.04) == 0.0
    assert trust(0.09, 0.04) == pytest.approx(0.5)
    assert trust(0.16, 0.04) == 1.0
    assert trust(1.0, 0.04) == 1.0


# ── capture only lengthens ─────────────────────────────────────────────────────────────────────────────


def test_an_rna_census_no_longer_than_gdna_moves_nothing():
    """The mixture's slope at ``c = 0`` is ``mu_r / mu_g − 1``: an RNA census no longer than gDNA's law gives ``c = 0``
    and ``G_c`` is gDNA's own law, bit for bit."""
    for r_mean in (150.0, 200.0, 224.0):
        tilted, c = capture_mixture(G, _law(r_mean, 82.0), 0.9)
        assert c == 0.0
        np.testing.assert_array_equal(tilted, G)


def test_the_mixture_never_exceeds_kappa_and_reaches_it_when_rna_is_far_longer():
    for kappa in (0.1, 0.5, 0.9):
        _, c = capture_mixture(G, _law(400.0, 90.0), kappa)
        assert c == kappa
        _, c = capture_mixture(G, size_biased(G, 0.3 * kappa), kappa)
        assert 0.0 < c <= kappa


def test_the_mixture_is_the_maximum_likelihood_fraction():
    """A brute-force grid over ``[0, kappa]`` never beats the returned ``c``."""
    r = _captured(G, 0.5, np.minimum(L, 180.0))
    kappa = 0.9
    _, c = capture_mixture(G, r, kappa)
    on = r > 0

    def loglik(x):
        return float((r[on] * np.log(size_biased(G, x)[on])).sum())

    best = max(loglik(x) for x in np.linspace(0.0, kappa, 2001))
    assert loglik(c) >= best - 1e-12


def test_more_capture_never_leaves_a_wider_contrast():
    """The family is nested in ``kappa``: a larger captured fraction can only bring gDNA's law nearer RNA's census."""
    r = _law(300.0, 90.0)

    def kl(a, b):
        on = (a > 0) & (b > 0)
        return float((a[on] * np.log(a[on] / b[on])).sum())

    previous = np.inf
    for kappa in np.linspace(0.0, 0.95, 39):
        tilted, _ = capture_mixture(G, r, kappa)
        here = kl(r, tilted)
        assert here <= previous + 1e-12
        previous = here


# ── the bound is a bound, and only under its premise ─────────────────────────────────────────────────────


def _proportional_responses(rng, n):
    """Random capture responses that reward length at most proportionally: ``m(l) = l s(l)``, ``s`` non-increasing and
    ``m`` non-decreasing — a probe overlap that grows with length and saturates. Built directly: each step keeps
    ``s(l + 1) / s(l)`` in ``[l / (l + 1), 1]``, the whole class, with half the draws concave."""
    out = []
    lo = np.divide(L[1:-1], L[2:])
    for k in range(n):
        if k % 2:
            steps = np.sort(rng.uniform(0.0, 1.0, L.size - 1))[::-1] ** rng.uniform(0.2, 5.0)
            out.append(np.concatenate(([0.0], np.cumsum(steps))))
            continue
        u = rng.uniform(0.0, 1.0, lo.size) ** rng.uniform(0.2, 5.0)
        s = np.concatenate(([1.0, 1.0], np.cumprod(lo + u * (1.0 - lo))))
        out.append(L * s)
    return out


def test_the_bound_is_the_worst_case_over_proportional_capture_responses():
    """For every response whose capture weight grows at most in proportion to length, the captured census is no
    further from the law than the size-biased one at the same captured fraction — in Hellinger distance, KL and
    chi-square alike (convex order)."""
    rng = np.random.default_rng(5)
    on = G > 0

    def kl(a, b):
        return float((a[on] * np.log(a[on] / b[on])).sum())

    def chi2(a, b):
        return float(((a[on] - b[on]) ** 2 / b[on]).sum())

    for m in _proportional_responses(rng, 300):
        for kappa in (0.2, 0.6, 0.95):
            census, linear = _captured(G, kappa, m), size_biased(G, kappa)
            assert hellinger_sq(G, census) <= capture_bound(G, kappa) * (1.0 + 1e-9)
            assert kl(census, G) <= kl(linear, G) * (1.0 + 1e-9)
            assert chi2(census, G) <= chi2(linear, G) * (1.0 + 1e-9)


def test_a_response_that_rewards_length_faster_than_proportionally_can_exceed_the_bound():
    """The premise is necessary: capture that takes only fragments above a length breaks the bound."""
    step = np.where(L > 300.0, L, 0.0)
    assert hellinger_sq(G, _captured(G, 0.6, step)) > 5.0 * capture_bound(G, 0.6)


# ── continuous in every input ──────────────────────────────────────────────────────────────────────────


def test_the_weight_is_continuous_in_the_contrast_and_the_bound():
    bound = 0.04
    contrasts = np.linspace(0.0, 0.25, 20001)
    w = np.array([trust(x, bound) for x in contrasts])
    assert np.max(np.abs(np.diff(w))) < 1e-3
    bounds = np.linspace(1e-4, 0.3, 20001)
    w = np.array([trust(0.1, b) for b in bounds])
    assert np.max(np.abs(np.diff(w))) < 1e-2


def test_the_mixture_is_continuous_in_kappa_and_through_equal_means():
    r = _law(300.0, 90.0)
    cs = np.array([capture_mixture(G, r, k)[1] for k in np.linspace(0.0, 0.95, 2001)])
    assert np.max(np.abs(np.diff(cs))) < 1e-3
    # through equal means: c falls to 0 continuously as RNA's census shortens to gDNA's law
    cs = np.array(
        [capture_mixture(G, size_biased(G, x), 0.9)[1] for x in np.linspace(0.05, 0.0, 501)]
    )
    assert cs[-1] == 0.0 and np.max(np.abs(np.diff(cs))) < 1e-3


# ── no capture, and the ends of the path ─────────────────────────────────────────────────────────────────


def test_with_no_capture_each_origin_keeps_its_own_law():
    """``kappa = 0``: no mixture, no bound, and any resolved difference whole — gDNA's own law, bit for bit."""
    g, r = _law(200.0, 70.0), _law(215.0, 75.0)
    t = _tables(r, g, 0.0)
    assert t.decision.captured == 0.0 and t.decision.bound == 0.0 and t.decision.weight == 1.0
    np.testing.assert_array_equal(t.gdna, g)
    np.testing.assert_array_equal(t.anchor, g)


def test_noise_that_covers_the_difference_is_one_table_even_without_capture():
    g, r = _law(200.0, 70.0), _law(201.0, 70.0)
    var = np.full(L.size, 1e-6) * g
    t = _tables(r, g, 0.0, var_r=var, var_g=var)
    assert hellinger_sq(g, r) < sampling_noise(g, r, var, var)
    assert t.decision.weight == 0.0
    np.testing.assert_array_equal(t.gdna, r)


def test_a_noisy_near_empty_tail_of_one_law_does_not_hide_a_real_gap():
    """Where gDNA's law is all but empty (its EB floor) and RNA's is not, a plug-in variance carried there must not read
    as a reason for one law: the noise is the null's, read at the two laws' common mass (cfRNA-B's case: gDNA's law
    near 100 bp, RNA's near 260, and gDNA's contrast variance spread over RNA's lengths)."""
    g, r = _law(110.0, 30.0), _law(260.0, 60.0)
    g = (g + 1e-12) / (g + 1e-12).sum()
    # at RNA's own lengths gDNA's law sits at its floor, while the contrast's intronic pool, which carries RNA there,
    # gives it a plug-in variance in proportion to RNA's mass, inflated by the inversion's coefficient
    var_g = np.where(L > 200.0, 1e-3 * r, 1e-9 * g)
    var_r = 1e-9 * r
    t = _tables(r, g, 0.7, var_r=var_r, var_g=var_g)
    assert t.decision.weight == 1.0
    assert t.decision.noise < 0.5 * hellinger_sq(g, r)


def test_the_geometric_path_is_exact_at_both_ends_and_a_law_between():
    r, g = _law(200.0, 70.0), _law(260.0, 80.0)
    np.testing.assert_array_equal(geometric_path(r, g, 0.0), r)
    np.testing.assert_array_equal(geometric_path(r, g, 1.0), g)
    mid = geometric_path(r, g, 0.5)
    assert mid.sum() == pytest.approx(1.0)
    assert np.all(mid[(r > 0) | (g > 0)] > 0.0)
    expected = np.sqrt(r * g)
    np.testing.assert_allclose(mid, expected / expected.sum(), rtol=1e-12, atol=1e-300)


def test_an_rna_law_no_spliced_fragment_measured_shows_no_difference():
    """An empty spliced pool leaves RNA's law at the anchor, a fallback no count measured, and its zero variance is no
    precision: the decision reads no RNA measurement as no difference and keeps one law, however far the anchor sits
    from gDNA's law, capture or none."""
    g, anchor = _law(200.0, 70.0), _law(260.0, 80.0)
    for kappa in (0.0, 0.7):
        t = shared_length_law(anchor, g, ZERO, ZERO, kappa, rna_measured=False)
        d = t.decision
        assert d.weight == 0.0 and d.contrast == 0.0 and d.captured == 0.0 and d.noise == 0.0
        assert d.kappa == kappa and d.bound == pytest.approx(np.sqrt(capture_bound(g, kappa)))
        np.testing.assert_array_equal(t.gdna, anchor)
        np.testing.assert_array_equal(t.rna, anchor)
        np.testing.assert_array_equal(t.anchor, g)


# ── the decision on built laws: pools nested in the anchor as production's are ─────────────────────────────

_BL = np.arange(601, dtype=np.float64)
_PI = np.clip((_BL - 30.0) / 400.0, 0.0, 1.0)


def _built_law(mean: float, sd: float) -> np.ndarray:
    p = stats.gamma.pdf(_BL, a=(mean / sd) ** 2, scale=sd**2 / mean)
    p[:20] = 0.0
    return p / p.sum()


def _built_tables(monkeypatch, g_law, r_law, n_rna, kappa, seed, *, spliced=True):
    """``build_fl_models`` on a payload whose gDNA pool is every gDNA fragment, whose spliced pool is RNA's fragments
    thinned by their junction-crossing chance (none with ``spliced`` False: RNA all unspliced), and whose anchor is
    every fragment of both; then the decision exactly as the pipeline hands it over (``pipeline.length_tables``), at
    a captured fraction ``kappa``."""
    from types import SimpleNamespace

    import rigel.calibration.capture_eff_length as cel
    from rigel.calibration.fl import build_fl_models
    from rigel.pipeline import length_tables
    from rigel.scan_payload import N_FRAGMENT_POOLS, POOL_DNA_INTERGENIC, POOL_RNA_SPLICED

    rng = np.random.default_rng(seed)
    g_all = rng.multinomial(400_000, g_law).astype(np.float64)
    r_all = rng.multinomial(n_rna, r_law).astype(np.float64)
    pools = np.zeros((N_FRAGMENT_POOLS, _BL.size))
    pools[POOL_DNA_INTERGENIC] = g_all
    if spliced:
        pools[POOL_RNA_SPLICED] = rng.binomial(r_all.astype(np.int64), _PI).astype(np.float64)
    payload = SimpleNamespace(
        pool_lengths=pools, deposited_lengths=g_all + r_all, max_length=_BL.size - 1
    )
    monkeypatch.setattr(cel, "rna_captured_fraction", lambda calibration, region_arrays: kappa)
    return length_tables(build_fl_models(payload, sj_opportunity=_PI), None, None)


def test_equal_chemistry_under_capture_is_one_law_at_every_rna_pool_size(monkeypatch):
    """RNA's census gDNA's law captured at half the captured fraction, from an empty spliced pool to a deep one: the
    built laws' own variances keep it one law in every library."""
    g = _built_law(225.0, 82.0)
    for n_rna in (0, 1_000, 10_000, 100_000):
        for seed in range(6):
            t = _built_tables(monkeypatch, g, size_biased(g, 0.3), n_rna, 0.6, seed)
            assert t.decision.weight == 0.0, (n_rna, seed)


@pytest.mark.parametrize("kappa", [0.0, 0.6])
def test_a_real_gap_measured_by_spliced_fragments_is_kept_whole(monkeypatch, kappa):
    g = _built_law(225.0, 82.0)
    for r_law in (_built_law(130.0, 50.0), _built_law(300.0, 90.0)):
        for n_rna in (1_000, 100_000):
            for seed in range(3):
                assert (
                    _built_tables(monkeypatch, g, r_law, n_rna, kappa, seed).decision.weight == 1.0
                )


@pytest.mark.parametrize("kappa", [0.0, 0.6])
def test_rna_with_no_spliced_fragment_shows_no_length_difference(monkeypatch, kappa):
    """RNA in the library but none of it spliced: its law is the anchor, which carries unspliced RNA of another length,
    and nothing measured RNA's own law — so one law, never the anchor read as RNA's census against gDNA's law."""
    g = _built_law(225.0, 82.0)
    t = _built_tables(monkeypatch, g, _built_law(130.0, 50.0), 200_000, kappa, 0, spliced=False)
    assert t.decision.weight == 0.0
    np.testing.assert_array_equal(t.gdna, t.rna)


def test_the_decision_reports_what_it_read():
    g, r = _law(130.0, 60.0), _law(260.0, 80.0)
    t = _tables(r, g, 0.7)
    d = t.decision
    assert d.kappa == 0.7 and 0.0 <= d.captured <= 0.7
    assert d.contrast == pytest.approx(np.sqrt(hellinger_sq(g, r)))
    assert d.bound == pytest.approx(np.sqrt(capture_bound(g, 0.7)))
    assert d.noise == 0.0 and d.weight == 1.0
