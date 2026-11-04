"""The sj pool's opportunity function, and the de-tilt it feeds.

``RNA_SPLICED`` is selected on "the path used an annotated sj", and that population is longer than the library, so the
pool the RNA fragment-length model is fitted from is tilted long; dividing it by its own opportunity ``A_j(w) = (L−w+1)+
− Σ_i (e_i−w+1)+`` removes the tilt. The correction divides by the theta-weighted crossing probability ``A(w)/T(w)``
rather than by ``A`` alone, because ``A`` recovers the distribution lengths were drawn from, while the pool must become
the one the library realizes: the drawn distribution weighted by each length's placements, which is ``T``. Under ``A``
alone a badly wrong theta makes the correction worse than none. Every oracle here enumerates and none calls the module
under test.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from rigel.calibration.sj_opportunity import (
    crossing_probability,
    crossing_probability_from_index,
    detilt_pool,
    sj_opportunity,
)

# ---------------------------------------------------------------------------
# oracles: enumerating, and sharing no code with the module under test
# ---------------------------------------------------------------------------


def _oracle_a_j(exon_lengths, w: int) -> int:
    """Count start positions whose ``[s, s+w)`` window strictly contains a sj coordinate."""
    total = int(sum(exon_lengths))
    region_bounds = list(itertools.accumulate(exon_lengths))[:-1]
    return sum(1 for s in range(0, total - w + 1) if any(s < int(c) < s + w for c in region_bounds))


def _oracle_total(exon_lengths, w: int) -> int:
    """Count every start position, the placements a length-``w`` window has at all."""
    return max(0, int(sum(exon_lengths)) - w + 1)


def _csr(transcripts):
    """``(exon_lengths, offsets)`` for a list of per-transcript exon-length lists."""
    offsets = np.zeros(len(transcripts) + 1, dtype=np.int64)
    for i, exons in enumerate(transcripts):
        offsets[i + 1] = offsets[i] + len(exons)
    lengths = np.array([e for exons in transcripts for e in exons], dtype=np.int64)
    return lengths, offsets


def _aggregate(transcripts, theta, max_length):
    """``(T, A)`` from :mod:`sj_opportunity`, for a list-of-lists transcriptome."""
    lengths, offsets = _csr(transcripts)
    return sj_opportunity(lengths, offsets, np.asarray(theta, dtype=np.float64), max_length)


# ---------------------------------------------------------------------------
# the formula itself
# ---------------------------------------------------------------------------


def test_the_opportunity_matches_an_enumerating_oracle_on_every_small_configuration():
    """Exhaustive over 1-4 exons x exon lengths 1-7 x every ``w`` up to ``L + 2``. The oracle walks every start position
    and tests it against every sj coordinate; the formula works with the complement (a window crossing nothing lies
    wholly inside one exon, and exons are disjoint). They share no code, so agreement is evidence.
    """
    checked = 0
    for k in range(1, 5):
        for exons in itertools.product(range(1, 8), repeat=k):
            _, a = _aggregate([list(exons)], [1.0], int(sum(exons)) + 2)
            for w in range(1, int(sum(exons)) + 3):
                assert a[w] == _oracle_a_j(exons, w), (exons, w, a[w])
                checked += 1
    assert checked > 48_000, checked


@pytest.mark.parametrize(
    "exons",
    [[1, 2000], [2000, 1], [50] * 20, [1, 1, 1], [500, 500], [713], [1], [200, 1, 200]],
)
def test_the_opportunity_matches_the_oracle_at_realistic_scales(exons):
    """The exhaustive sweep is all tiny. A 1 bp exon beside a 2 kb one is the shape that breaks an off-by-one, and 20 x
    50 bp is the many-short-exons transcript the tilt is steepest on."""
    total = int(sum(exons))
    _, a = _aggregate([exons], [1.0], total + 2)
    for w in {1, 2, 49, 50, 51, 200, 499, 500, 501, total, total + 1, total + 2} - {0}:
        if w <= total + 2:
            assert a[w] == _oracle_a_j(exons, w), (exons, w)


def test_the_five_properties_the_derivation_claims_are_checked_not_stated():
    """P1-P5 of the derivation, each a property the correction depends on."""
    # P1 — a single-exon transcript can never populate the sj pool
    _, a = _aggregate([[500]], [1.0], 600)
    assert not a.any()

    # P2 — a 1 bp fragment crosses no 0-bp boundary
    _, a = _aggregate([[10, 10, 10]], [1.0], 40)
    assert a[1] == 0

    # P3: A_j rises with w up to the longest exon. This is the tilt being corrected.
    _, a = _aggregate([[100, 100]], [1.0], 200)
    rising = a[1:101]
    assert np.all(np.diff(rising) >= 0) and rising[-1] > rising[0]

    # P4: beyond the longest exon every placement crosses, so the tilt saturates
    total_a, a = _aggregate([[100, 100]], [1.0], 200)
    for w in range(101, 201):
        assert a[w] == total_a[w] == max(0, 200 - w + 1)

    # P5: two transcripts of equal length but different structure differ, so the aggregate depends on which transcripts
    # are expressed, which is why theta appears.
    _, few = _aggregate([[100, 100]], [1.0], 200)
    _, many = _aggregate([[20] * 10], [1.0], 200)
    assert not np.array_equal(few, many)
    assert many[30] > few[30]


def test_the_vectorised_aggregate_is_the_scalar_formula_summed():
    """The kernel builds ``T`` and ``A`` from two length spectra with a double reverse-cumsum; this re-derives them one
    transcript at a time by direct enumeration on a mixed transcriptome.
    """
    transcripts = [[7], [3, 11], [50] * 4, [1, 400], [120, 5, 9, 60]]
    theta = [3.0, 0.0, 17.5, 1.0, 250.0]
    t, a = _aggregate(transcripts, theta, 500)
    for w in (1, 2, 5, 60, 121, 200, 401, 500):
        want_t = sum(th * _oracle_total(e, w) for th, e in zip(theta, transcripts))
        want_a = sum(th * _oracle_a_j(e, w) for th, e in zip(theta, transcripts))
        assert t[w] == pytest.approx(want_t)
        assert a[w] == pytest.approx(want_a)


# ---------------------------------------------------------------------------
# the crossing probability
# ---------------------------------------------------------------------------


def test_the_crossing_probability_is_a_probability_and_saturates_at_one():
    pi = crossing_probability(*_csr([[10, 10]]), np.array([1.0]), 30)
    assert np.all((pi >= 0.0) & (pi <= 1.0))
    # Bin 0 included: the complement identity is negative there (a zero-length window sits on a boundary, so every
    # internal one is counted twice and `A_j(0) = 1 - K`). It is not a fragment length, but a negative divisor must not
    # be left in the array.
    assert pi[0] == 0.0
    assert pi[1] == 0.0
    # every placement of a window longer than the longest exon must cross
    for w in range(11, 21):
        assert pi[w] == pytest.approx(1.0)
    # past the transcript there are no placements at all, so the probability is undefined -> 0
    assert pi[21] == 0.0


def test_a_single_exon_annotation_makes_the_correction_inert_not_a_division_by_zero():
    """No sj anywhere means no opportunity anywhere, so the correction has nothing to say and must neither divide by
    zero nor delete the pool."""
    pi = crossing_probability(*_csr([[500], [900]]), np.array([1.0, 1.0]), 600)
    assert not pi.any()
    counts = np.zeros(601)
    counts[100:400] = 7.0
    np.testing.assert_array_equal(detilt_pool(counts, pi), counts)


# ---------------------------------------------------------------------------
# Does the correction recover the library distribution?
# ---------------------------------------------------------------------------


def _enumerate_library(transcripts, theta, pmf):
    """Every placement, enumerated: ``(library, pool)`` histograms over fragment length. A fragment of length ``w`` is
    drawn with probability ``pmf[w]`` and placed uniformly over the start positions available to it on a transcript
    chosen with probability proportional to ``theta``. ``library`` counts every placement; ``pool`` only those crossing
    a sj. This is the generative model the opportunity function claims to invert.
    """
    n = len(pmf)
    library = np.zeros(n)
    pool = np.zeros(n)
    for th, exons in zip(theta, transcripts):
        total = sum(exons)
        region_bounds = list(itertools.accumulate(exons))[:-1]
        for w in range(1, n):
            if pmf[w] == 0.0 or w > total:
                continue
            for s in range(0, total - w + 1):
                library[w] += th * pmf[w]
                if any(s < int(c) < s + w for c in region_bounds):
                    pool[w] += th * pmf[w]
    return library, pool


def _mean_length(h) -> float:
    """The mean of a length histogram."""
    return float(np.dot(np.arange(h.size), h) / h.sum())


def test_the_correction_recovers_the_library_distribution_from_an_enumerated_pool():
    """Enumerate every placement to get the library histogram and its sj-crossing subset, then de-tilt the subset: it
    must come back to the library histogram, shape-exactly since the de-tilt preserves the total, with nothing tuned and
    no tolerance beyond float.
    """
    transcripts = [[30, 40], [12] * 6, [70], [25, 8, 61], [100, 3]]
    theta = [5.0, 2.0, 9.0, 1.0, 4.0]
    pmf = np.zeros(91)
    pmf[10:81] = np.exp(-0.5 * ((np.arange(10, 81) - 45) / 12.0) ** 2)
    pmf /= pmf.sum()

    library, pool = _enumerate_library(transcripts, theta, pmf)
    assert pool.sum() > 0 and pool.sum() < library.sum()

    pi = crossing_probability(*_csr(transcripts), np.array(theta), len(pmf) - 1)
    corrected = detilt_pool(pool, pi)

    want = library * (pool.sum() / library.sum())
    np.testing.assert_allclose(corrected, want, rtol=1e-12, atol=1e-12)

    # and the raw pool is not already the answer, or the test would pass with no correction at all
    assert not np.allclose(pool, want, rtol=1e-3)


def test_the_correction_moves_the_pool_mean_onto_the_librarys():
    """The same statement in the one number the calibrator consumes."""
    transcripts = [[30, 40], [12] * 6, [70], [25, 8, 61], [100, 3]]
    theta = [5.0, 2.0, 9.0, 1.0, 4.0]
    pmf = np.zeros(91)
    pmf[10:81] = 1.0 / 71.0
    library, pool = _enumerate_library(transcripts, theta, pmf)
    pi = crossing_probability(*_csr(transcripts), np.array(theta), len(pmf) - 1)

    assert _mean_length(pool) > _mean_length(library)  # the tilt exists, and it is long
    assert _mean_length(detilt_pool(pool, pi)) == pytest.approx(_mean_length(library), rel=1e-12)


def test_the_ratio_form_is_what_makes_a_wrong_theta_safe():
    """With a deliberately wrong theta (all the weight on the transcript with the shortest exons, which has the steepest
    tilt), dividing by ``A`` alone overshoots past the library mean, worse than no correction. The ratio form ``A / T``
    must not, because a reweighting moves ``A`` and ``T`` together.
    """
    transcripts = [[30, 40], [6] * 12, [70], [25, 8, 61], [100, 3]]
    theta = [5.0, 2.0, 9.0, 1.0, 4.0]
    pmf = np.zeros(91)
    pmf[10:81] = 1.0 / 71.0
    library, pool = _enumerate_library(transcripts, theta, pmf)

    wrong = np.array([0.0, 1.0, 0.0, 0.0, 0.0])  # every copy is the many-short-exon transcript
    lengths, offsets = _csr(transcripts)
    pi_wrong = crossing_probability(lengths, offsets, wrong, len(pmf) - 1)
    _, a_wrong = sj_opportunity(lengths, offsets, wrong, len(pmf) - 1)

    ratio_err = abs(_mean_length(detilt_pool(pool, pi_wrong)) - _mean_length(library))
    a_only_err = abs(
        _mean_length(detilt_pool(pool, a_wrong / max(a_wrong.max(), 1.0))) - _mean_length(library)
    )
    raw_err = abs(_mean_length(pool) - _mean_length(library))

    assert ratio_err < raw_err, (ratio_err, raw_err)
    assert a_only_err > raw_err, (a_only_err, raw_err)


# ---------------------------------------------------------------------------
# the de-tilt's own contract
# ---------------------------------------------------------------------------


def test_the_detilt_preserves_the_evidence_count():
    """The de-tilt changes the pool's shape, never how much evidence it represents: ``build_fl_models`` shrinks each
    pool toward the anchor with a Dirichlet pseudo-count whose strength is the pool total, so inflating the total would
    weaken the shrinkage as a side effect.
    """
    pi = crossing_probability(*_csr([[40, 60], [10] * 10]), np.array([3.0, 1.0]), 120)
    counts = np.zeros(121)
    counts[20:110] = np.arange(90, dtype=np.float64)
    out = detilt_pool(counts, pi)
    assert out.sum() == pytest.approx(counts.sum())
    assert not np.allclose(out, counts)


def test_the_detilt_cannot_invent_mass_where_the_pool_has_none():
    pi = crossing_probability(*_csr([[40, 60]]), np.array([1.0]), 120)
    counts = np.zeros(121)
    counts[50] = 100.0
    out = detilt_pool(counts, pi)
    assert out[50] == pytest.approx(100.0)
    assert not out[np.arange(121) != 50].any()


# ---------------------------------------------------------------------------
# the index adapter
# ---------------------------------------------------------------------------


def test_the_index_adapter_uses_the_real_transcripts_and_only_those(mini_index):
    """``~is_synthetic`` alone: the synthetic nascent spans are not sequenced molecules, and weighting them would put
    opportunity on exon structures the library does not contain. ``is_single_exon`` is not a realness filter: a real
    single-exon transcript is a real transcript.
    """
    max_length = 400
    pi = crossing_probability_from_index(mini_index, max_length)

    real = ~mini_index.t_df["is_synthetic"].to_numpy()
    assert real.sum() == 3, "MINI_GTF has three real transcripts"
    transcripts = [
        [int(b - a) for a, b in mini_index.get_exon_intervals(int(i))] for i in np.nonzero(real)[0]
    ]
    want = crossing_probability(*_csr(transcripts), np.ones(len(transcripts)), max_length)
    np.testing.assert_allclose(pi, want, rtol=0, atol=0)

    # and it differs from weighting everything, or the filter would be doing nothing
    assert int((~real).sum()) > 0, "MINI_GTF grew no synthetic rows — this test proves nothing"
    every = []
    for i in range(len(mini_index.t_df)):
        intervals = mini_index.get_exon_intervals(int(i))
        if intervals is not None and len(intervals):
            every.append([int(b - a) for a, b in intervals])
    assert not np.allclose(pi, crossing_probability(*_csr(every), np.ones(len(every)), max_length))


def test_the_index_adapter_agrees_with_the_enumerating_oracle_on_a_real_index(mini_index):
    """The adapter, the kernel and a hand enumeration on the same three transcripts, so a wrong exon-length unpacking
    cannot hide behind the kernel."""
    pi = crossing_probability_from_index(mini_index, 400)
    real = np.nonzero(~mini_index.t_df["is_synthetic"].to_numpy())[0]
    transcripts = [[int(b - a) for a, b in mini_index.get_exon_intervals(int(i))] for i in real]
    for w in (1, 2, 100, 101, 102, 202, 303, 304):
        want_a = sum(_oracle_a_j(e, w) for e in transcripts)
        want_t = sum(_oracle_total(e, w) for e in transcripts)
        assert pi[w] == pytest.approx(want_a / want_t if want_t else 0.0)


# ---------------------------------------------------------------------------
# the wiring
# ---------------------------------------------------------------------------


def test_the_sj_divisor_reaches_the_rna_pool_alone_and_keeps_its_evidence():
    """The models' kernel reads the spliced pool at its sj opportunity: the divisor moves RNA's law, leaves gDNA's
    untouched, and adds no evidence weight."""
    from rigel.calibration.fl import fl_models_from_histograms
    from rigel.scan_payload import N_FRAGMENT_POOLS, POOL_DNA_INTERGENIC, POOL_RNA_SPLICED

    pools = np.zeros((N_FRAGMENT_POOLS, 51), dtype=np.int64)
    pools[POOL_RNA_SPLICED, 20] = 300
    pools[POOL_RNA_SPLICED, 40] = 100
    pools[POOL_DNA_INTERGENIC, 30] = 500
    anchor = pools.sum(axis=0)  # every pool count is a deposited fragment
    anchor[25] += 900

    def models(selection):
        return fl_models_from_histograms(
            global_counts=anchor,
            rna_counts=pools[POOL_RNA_SPLICED],
            rna_selection=selection,
            gdna_counts=pools[POOL_DNA_INTERGENIC],
            max_size=50,
        )

    pi = np.zeros(51)
    pi[1:] = np.linspace(0.1, 0.9, 50)
    plain, tilted = models(None), models(pi)
    assert not np.allclose(tilted.rna_pmf, plain.rna_pmf), "the divisor reached nothing"
    np.testing.assert_array_equal(
        tilted.gdna_pmf, plain.gdna_pmf
    )  # the sj divisor must not touch the gDNA pool
    assert tilted.n_rna == pytest.approx(plain.n_rna), "the EB evidence weight moved"
