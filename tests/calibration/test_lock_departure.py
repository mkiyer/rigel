"""The locks' departure: what a locked object's two strands can share.

An object the annotation admits no RNA strand at is booked as gDNA. gDNA lands on either genome strand evenly, a
transcribed molecule on its own strand's column, so a stranded library's locked objects say which of them hold
something else, and how much of each one's count its two strands can share. An unstranded library's say nothing.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import binom

from rigel.calibration.lock_departure import LockFit, fit_lock_departure

KAPPA = 0.99


def _gdna(rng, mean, size):
    """All-gDNA objects: a Poisson count at ``mean``, its strands balanced."""
    g = rng.poisson(mean, size)
    u = rng.binomial(g, 0.5)
    return u, g - u


def _with_rna(rng, u, v, rna, on_plus, kappa=KAPPA):
    """The same objects with ``rna`` transcribed fragments each, on one strand at ``kappa``."""
    s = rng.binomial(rna, kappa, u.shape[0])
    return (u + s, v + rna - s) if on_plus else (u + rna - s, v + s)


# --- the three properties the form has by construction ----------------------------------------------------------------


def test_a_dead_strand_channel_is_the_lock():
    """At κ = ½ transcribed fragments are balanced too: the strands say nothing, and every object keeps its count,
    however one-sided."""
    fit = fit_lock_departure([40, 3, 0, 17], [0, 3, 0, 2], 0.5)
    assert fit.departure == 0.0 and fit.evidence == 0.0 and fit.n_objects == 3
    assert np.all(fit.share == 1.0) and np.all(fit.variance == 0.0)
    assert np.all(fit.rna_pos == 0.0) and np.all(fit.rna_neg == 0.0)


@pytest.mark.parametrize("seed", range(8))
def test_an_all_gdna_library_is_the_lock(seed):
    """Three thousand objects holding only gDNA, sixty fragments each on average. Binomial noise is not departure: the
    fit keeps their gDNA to a part in a thousand, where the bare shared part, twice the smaller strand, gives up a
    tenth of it."""
    rng = np.random.default_rng(seed)
    u, v = _gdna(rng, 60.0, 3000)
    fit = fit_lock_departure(u, v, KAPPA)
    n = (u + v).astype(float)
    assert (fit.share * n).sum() > 0.999 * n.sum(), (fit.departure, (fit.share * n).sum() / n.sum())
    assert (2.0 * np.minimum(u, v)).sum() < 0.92 * n.sum()


def test_most_all_gdna_libraries_are_the_lock_exactly():
    """Under the lock the likelihood's slope at ε = 0 has mean zero and a long right tail, so its maximum is at 0 in
    most all-gDNA libraries, and the fit is then the lock to the bit."""
    exact = 0
    for seed in range(40):
        u, v = _gdna(np.random.default_rng(100 + seed), 30.0, 2000)
        fit = fit_lock_departure(u, v, KAPPA)
        exact += fit.departure == 0.0 and bool(np.all(fit.share == 1.0))
    assert exact >= 30, exact


def test_transcription_on_both_strands_alike_stays_gdna():
    """A locked object holding as much transcription from one strand as from the other is balanced. Its strands cannot
    tell it from gDNA, and it keeps its count."""
    rng = np.random.default_rng(7)
    u, v = _gdna(rng, 20.0, 2000)
    u, v = _with_rna(rng, u, v, 30, True)
    u, v = _with_rna(rng, u, v, 30, False)
    fit = fit_lock_departure(u, v, KAPPA)
    n = (u + v).astype(float)
    assert (fit.share * n).sum() > 0.995 * n.sum(), fit.departure


# --- gDNA until proven otherwise ----------------------------------------------------------------------------------------


def test_a_lock_is_gdna_until_its_own_strands_say_otherwise():
    """A library where one lock in three is one-sided, so its rate of departed locks is high. A lock whose own strands
    do not favour the departed reading is still all gDNA, exactly: one fragment, which says nothing either way; a
    balanced pair; three and three; thirty and thirty-one. The library's rate tempers how far a one-sided lock is
    believed. It never moves a lock by itself."""
    u = np.concatenate([[1, 0, 1, 2, 3, 10, 30, 301], np.full(200, 40), np.full(400, 20)]).astype(
        float
    )
    v = np.concatenate([[0, 1, 1, 1, 3, 10, 31, 300], np.zeros(200), np.full(400, 20)]).astype(
        float
    )
    fit = fit_lock_departure(u, v, KAPPA)
    assert fit.departure > 0.2, fit.departure
    assert np.all(fit.share[:8] == 1.0), fit.share[:8]
    assert np.all(fit.variance[:8] == 0.0)
    assert np.all(fit.rna_pos[:8] == 0.0) and np.all(fit.rna_neg[:8] == 0.0)
    assert np.all(fit.share[8:208] < 0.2)


@pytest.mark.parametrize(("gdna", "kept_at_least"), [(13.0, 0.98), (60.0, 0.995)])
def test_all_gdna_locks_among_departed_ones_keep_their_gdna(gdna, kept_at_least):
    """Two thousand all-gDNA locks, holding a dozen fragments each or sixty, beside a hundred that also hold eighty
    transcribed ones. An all-gDNA lock is one-sided only by chance, and the deeper it is the less chance gives: the
    dozen-fragment locks keep 98 % of their gDNA and the sixty-fragment ones 99.5 %."""
    rng = np.random.default_rng(17)
    u0, v0 = _gdna(rng, gdna, 2000)
    u1, v1 = _with_rna(rng, *_gdna(rng, gdna, 100), 80, on_plus=True)
    fit = fit_lock_departure(np.concatenate([u0, u1]), np.concatenate([v0, v1]), KAPPA)
    n = (u0 + v0).astype(float)
    kept = fit.share[:2000] * n
    assert kept.sum() > kept_at_least * n.sum(), (fit.departure, kept.sum() / n.sum())


def test_where_every_lock_is_shallow_each_is_read_on_its_own_strands():
    """Two thousand all-gDNA locks of a fragment or two beside a hundred departed ones. Shallow locks cannot pin the
    library's rate, which reads high, so each lock stands on its own strands: one fragment, or a balanced two or three,
    is all gDNA exactly; two or three on one strand are believed departed in part. The all-gDNA locks keep 19 fragments
    in 20."""
    rng = np.random.default_rng(17)
    u0, v0 = _gdna(rng, 1.3, 2000)
    u1, v1 = _with_rna(rng, *_gdna(rng, 1.3, 100), 80, on_plus=True)
    fit = fit_lock_departure(np.concatenate([u0, u1]), np.concatenate([v0, v1]), KAPPA)
    n, share = (u0 + v0).astype(float), fit.share[:2000]
    assert fit.departure > 0.5, fit.departure
    assert np.all(share[n == 1.0] == 1.0)
    balanced = (n >= 2.0) & (np.abs(u0 - v0) <= 1)
    assert balanced.sum() > 100 and np.all(share[balanced] == 1.0)
    one_sided = (n >= 2.0) & (np.minimum(u0, v0) == 0)
    assert one_sided.sum() > 100 and np.all(share[one_sided] < 1.0)
    pair = one_sided & (n == 2.0)
    assert pair.sum() > 50 and np.all(share[pair] > 0.8), share[pair].min()
    assert (share * n).sum() > 0.95 * n.sum(), (share * n).sum() / n.sum()


# --- a departed lock --------------------------------------------------------------------------------------------------


def _mixed(seed, *, gdna=60.0, rna=80, n_null=2000, n_departed=100):
    rng = np.random.default_rng(seed)
    u0, v0 = _gdna(rng, gdna, n_null)
    g = _gdna(rng, gdna, n_departed)
    u1, v1 = _with_rna(rng, *g, rna, on_plus=True)
    return np.concatenate([u0, u1]), np.concatenate([v0, v1]), n_null, (g[0] + g[1]).astype(float)


def test_a_departed_lock_books_what_its_strands_share():
    """One object in twenty holds eighty transcribed fragments on its sixty of gDNA. The departed objects book their
    gDNA, within a tenth, and not their count; the rest keep theirs to a part in a hundred."""
    u, v, n_null, gdna = _mixed(1)
    fit = fit_lock_departure(u, v, KAPPA)
    n = (u + v).astype(float)
    booked = fit.share * n
    assert 0.9 * gdna.sum() < booked[n_null:].sum() < 1.1 * gdna.sum(), (
        booked[n_null:].sum(),
        gdna.sum(),
        n[n_null:].sum(),
    )
    assert booked[:n_null].sum() > 0.99 * n[:n_null].sum(), booked[:n_null].sum() / n[:n_null].sum()
    assert fit.departure > 0.0 and fit.evidence > 100.0


def test_the_remainder_is_on_the_strand_that_holds_it():
    """The simplex closes, and the transcribed part is booked to the strand whose column holds the excess: the plus
    strand's under a sense protocol, the minus strand's when read 1 is antisense."""
    u, v, n_null, _ = _mixed(2)
    for kappa, pos_heavy in ((KAPPA, True), (1.0 - KAPPA, False)):
        fit = fit_lock_departure(u, v, kappa)
        np.testing.assert_allclose(fit.share + fit.rna_pos + fit.rna_neg, 1.0, rtol=0, atol=1e-12)
        pos, neg = fit.rna_pos[n_null:].sum(), fit.rna_neg[n_null:].sum()
        assert (pos > 20.0 * neg) if pos_heavy else (neg > 20.0 * pos), (kappa, pos, neg)


def test_a_library_with_no_gdna_reads_every_lock_as_departed():
    """Every locked object one-sided: the library's rate is 1, and each books a few fragments of its forty, two or
    three on average: twice its smaller strand and the width of that."""
    rng = np.random.default_rng(3)
    s = rng.binomial(40, KAPPA, 200)
    fit = fit_lock_departure(s, 40 - s, KAPPA)
    assert fit.departure > 0.99
    booked = fit.share * 40.0
    assert booked.mean() < 3.0 and booked.max() < 10.0, (booked.mean(), booked.max())


def test_a_one_sided_lock_is_believed_more_as_its_count_grows():
    """In a library where one lock in six is one-sided, a lock holding fragments on one strand only is believed
    departed by degrees: one fragment proves nothing, two a little, and each further one more. A spectrum, with no
    count at which the reading switches."""
    kept = []
    for n in (1, 2, 3, 5, 8, 13, 21):
        u = np.concatenate([[n], np.full(300, 20), np.full(60, 40)]).astype(float)
        v = np.concatenate([[0], np.full(300, 20), np.zeros(60)]).astype(float)
        kept.append(float(fit_lock_departure(u, v, KAPPA).share[0]))
    assert kept[0] == 1.0, kept
    assert all(a > b for a, b in zip(kept, kept[1:])), kept
    assert kept[1] > 0.9 and kept[-1] < 0.2, kept


def test_a_lock_that_shares_about_one_fragment_has_no_location():
    """A one-sided lock's shared part is a bound near zero, a fragment or two give or take as much. Its log share is
    wider than one nat², the landscape's location floor, so it does not train the landscape; a lock that is all gDNA
    is certain."""
    u = np.concatenate([np.full(300, 30), np.full(60, 40)])
    v = np.concatenate([np.full(300, 30), np.zeros(60)])
    fit = fit_lock_departure(u, v, KAPPA)
    assert np.all(fit.variance[300:] > 1.0), fit.variance[300]
    assert np.all(fit.variance[:300] < 0.01), fit.variance[0]


def test_the_two_strands_are_read_alike():
    """Swapping an object's two columns swaps the strand its remainder is booked to and changes nothing else, from one
    fragment to two hundred on one strand: the share and its spread are the same to rounding, and finite."""
    u = np.array([1, 0, 5, 40, 200, 200, 3, 60, 31, 0], float)
    v = np.array([0, 2, 0, 1, 0, 2, 3, 20, 29, 200], float)
    one, other = fit_lock_departure(u, v, KAPPA), fit_lock_departure(v, u, KAPPA)
    assert one.departure == pytest.approx(other.departure, rel=1e-9)
    np.testing.assert_allclose(one.share, other.share, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(one.rna_pos, other.rna_neg, rtol=1e-9, atol=1e-12)
    assert np.all(np.isfinite(one.variance)) and np.all(np.isfinite(other.variance))
    np.testing.assert_allclose(one.variance, other.variance, rtol=1e-6, atol=1e-9)


def test_an_object_with_no_fragment_keeps_the_lock():
    fit = fit_lock_departure([0, 40, 0], [0, 0, 0], KAPPA)
    assert fit.n_objects == 1
    assert fit.share[0] == 1.0 and fit.share[2] == 1.0 and fit.variance[0] == 0.0


def test_counts_are_whole_fragments():
    with pytest.raises(ValueError, match="whole fragments"):
        fit_lock_departure([1.5], [1.0], KAPPA)
    with pytest.raises(ValueError, match="aligned"):
        fit_lock_departure([1.0, 2.0], [1.0], KAPPA)


# --- the closed forms against brute force -----------------------------------------------------------------------------


def _brute(u, n, kappa, eps):
    """One object's null probability, share and Var(log share) by direct integration over the departed share."""
    a = binom.pmf(u, n, 0.5)
    parts = []
    for q in (kappa, 1.0 - kappa):

        def lik(f, q=q):
            return binom.pmf(u, n, 0.5 * f + (1.0 - f) * q)

        b = quad(lik, 0.0, 1.0, limit=200)[0]
        m1 = quad(lambda f: f * lik(f), 0.0, 1.0, limit=200)[0]
        l1 = quad(lambda f: np.log(f) * lik(f), 0.0, 1.0, limit=200)[0]
        l2 = quad(lambda f: np.log(f) ** 2 * lik(f), 0.0, 1.0, limit=200)[0]
        parts.append((b, m1, l1, l2))
    b = 0.5 * (parts[0][0] + parts[1][0])
    w = min(1.0, a / ((1.0 - eps) * a + eps * b))
    total = parts[0][0] + parts[1][0]
    mean = (parts[0][1] + parts[1][1]) / total
    first = (1.0 - w) * (parts[0][2] + parts[1][2]) / total
    second = (1.0 - w) * (parts[0][3] + parts[1][3]) / total
    return w, w + (1.0 - w) * mean, second - first * first


@pytest.mark.parametrize("kappa", [0.99, 0.9, 0.7, 0.02])
def test_the_closed_forms_match_direct_integration(kappa):
    """Each object's share and the spread of its log share, against quadrature over the departed share at the fitted
    rate: the share to 1e-9, the variance to the midpoint quantiles' discretization."""
    u = np.array([0, 1, 3, 5, 10, 12, 20, 38, 40, 7, 30, 30], float)
    v = np.array([1, 0, 3, 0, 10, 2, 22, 2, 0, 33, 30, 28], float)
    fit: LockFit = fit_lock_departure(u, v, kappa)
    assert 0.0 < fit.departure < 1.0
    for j in range(u.shape[0]):
        w, share, var = _brute(u[j], u[j] + v[j], kappa, fit.departure)
        assert fit.share[j] == pytest.approx(share, rel=1e-9, abs=1e-12), (j, u[j], v[j])
        assert fit.variance[j] == pytest.approx(var, rel=0.05, abs=2e-4), (j, u[j], v[j], w)
