"""The junction pool's opportunity function, and the de-tilt it feeds.

``RNA_SPLICED`` is selected on *"the path used an annotated junction"*, and that population is
genuinely longer than the library — seeing a splice is roughly length-independent while having an
unsequenced mate gap is a pure length threshold. So the pool the RNA fragment-length model is FITTED
FROM is tilted long, and dividing it by its own opportunity is what removes the tilt.

The quantity is

    A_j(w)  =  (L - w + 1)+  -  SUM_i (e_i - w + 1)+

for a transcript of exon lengths ``e_1 .. e_K`` and total ``L``, and the correction divides the pool
by the theta-weighted CROSSING PROBABILITY ``pi(w) = A(w) / T(w)`` rather than by ``A(w)`` alone.

⛔ **Why the ratio and not ``A`` alone, which is what the derivation's own worked examples used.**
``A(w)`` recovers the distribution fragment lengths were DRAWN from; the pool has to be turned into
the distribution the library REALIZES, which is the drawn one weighted by how many placements each
length has. That is ``T(w)``. Empirically the difference is not cosmetic: dividing by ``A`` alone,
a badly wrong theta makes the correction WORSE than no correction, and dividing by ``pi`` it never
does — see :func:`test_the_ratio_form_is_what_makes_a_WRONG_theta_SAFE`.

⚠ Every oracle here ENUMERATES. None of them calls the module under test — a validator that shares
the implementation's helper validates nothing.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from rigel.calibration.junction_opportunity import (
    crossing_probability,
    crossing_probability_from_index,
    detilt_pool,
    junction_opportunity,
)

# ---------------------------------------------------------------------------
# oracles — enumerating, and sharing no code with the module under test
# ---------------------------------------------------------------------------


def _oracle_a_j(exon_lengths, w: int) -> int:
    """Count start positions whose ``[s, s+w)`` window strictly contains a junction coordinate."""
    total = int(sum(exon_lengths))
    region_bounds = list(itertools.accumulate(exon_lengths))[:-1]
    return sum(1 for s in range(0, total - w + 1) if any(s < int(c) < s + w for c in region_bounds))


def _oracle_total(exon_lengths, w: int) -> int:
    """Count every start position — the placements a length-``w`` window has at all."""
    return max(0, int(sum(exon_lengths)) - w + 1)


def _csr(transcripts):
    """``(exon_lengths, offsets)`` for a list of per-transcript exon-length lists."""
    offsets = np.zeros(len(transcripts) + 1, dtype=np.int64)
    for i, exons in enumerate(transcripts):
        offsets[i + 1] = offsets[i] + len(exons)
    lengths = np.array([e for exons in transcripts for e in exons], dtype=np.int64)
    return lengths, offsets


def _aggregate(transcripts, theta, max_length):
    """``(T, A)`` from :mod:`junction_opportunity`, for a list-of-lists transcriptome."""
    lengths, offsets = _csr(transcripts)
    return junction_opportunity(lengths, offsets, np.asarray(theta, dtype=np.float64), max_length)


# ---------------------------------------------------------------------------
# the formula itself
# ---------------------------------------------------------------------------


def test_A_j_matches_an_ENUMERATING_oracle_over_every_small_configuration():
    """⭐ The derivation, re-proven in the suite rather than cited from a document.

    Exhaustive over 1-4 exons x exon lengths 1-7 x every ``w`` up to ``L + 2``. The oracle walks every
    start position and tests it against every junction coordinate; the formula works with the
    complement (a window crossing nothing lies wholly inside ONE exon, and exons are disjoint). They
    share no code, so agreement is evidence.
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
def test_A_j_matches_the_oracle_at_REALISTIC_scales_too(exons):
    """⚠ The exhaustive sweep is all tiny. A 1 bp exon beside a 2 kb one is the shape that breaks
    an off-by-one, and 20 x 50 bp is the many-short-exons transcript the tilt is steepest on."""
    total = int(sum(exons))
    _, a = _aggregate([exons], [1.0], total + 2)
    for w in {1, 2, 49, 50, 51, 200, 499, 500, 501, total, total + 1, total + 2} - {0}:
        if w <= total + 2:
            assert a[w] == _oracle_a_j(exons, w), (exons, w)


def test_the_five_properties_the_derivation_CLAIMS_are_CHECKED_not_stated():
    """P1-P5 of the derivation, each one a thing the correction silently depends on."""
    # P1 — a single-exon transcript can never populate the junction pool
    _, a = _aggregate([[500]], [1.0], 600)
    assert not a.any()

    # P2 — a 1 bp fragment crosses no 0-bp boundary
    _, a = _aggregate([[10, 10, 10]], [1.0], 40)
    assert a[1] == 0

    # P3 — A_j RISES with w up to the longest exon. ⭐ This IS the tilt being corrected.
    _, a = _aggregate([[100, 100]], [1.0], 200)
    rising = a[1:101]
    assert np.all(np.diff(rising) >= 0) and rising[-1] > rising[0]

    # P4 — beyond the longest exon every placement crosses, so the tilt SATURATES
    total_a, a = _aggregate([[100, 100]], [1.0], 200)
    for w in range(101, 201):
        assert a[w] == total_a[w] == max(0, 200 - w + 1)

    # P5 — ⛔ two transcripts of EQUAL length but different structure differ, so the aggregate
    #      depends on WHICH transcripts are expressed. This is the whole reason theta appears.
    _, few = _aggregate([[100, 100]], [1.0], 200)
    _, many = _aggregate([[20] * 10], [1.0], 200)
    assert not np.array_equal(few, many)
    assert many[30] > few[30]


def test_the_vectorised_aggregate_is_the_SCALAR_formula_SUMMED():
    """Two implementations of one quantity, checked against each other on a mixed transcriptome.

    The kernel builds ``T`` and ``A`` from two length spectra with a double reverse-cumsum; this
    re-derives them one transcript at a time by direct enumeration. Trap 27: two implementations of
    one quantity that nobody ever diffed.
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


def test_the_crossing_probability_is_a_PROBABILITY_and_saturates_at_one():
    pi = crossing_probability(*_csr([[10, 10]]), np.array([1.0]), 30)
    assert np.all((pi >= 0.0) & (pi <= 1.0))
    # ⛔ Bin 0 included: the complement identity is NEGATIVE there (a zero-length window sits ON a
    # boundary, so every internal one is counted twice and `A_j(0) = 1 - K`). It is not a fragment
    # length, but a negative divisor left in the array is a landmine for the next consumer.
    assert pi[0] == 0.0
    assert pi[1] == 0.0
    # every placement of a window longer than the longest exon must cross
    for w in range(11, 21):
        assert pi[w] == pytest.approx(1.0)
    # past the transcript there are no placements at all, so the probability is undefined -> 0
    assert pi[21] == 0.0


def test_a_SINGLE_EXON_annotation_makes_the_correction_INERT_not_a_division_by_zero():
    """⛔ No junction anywhere means no opportunity anywhere. The correction must have NOTHING to
    say, and in particular must not divide by zero or delete the pool."""
    pi = crossing_probability(*_csr([[500], [900]]), np.array([1.0, 1.0]), 600)
    assert not pi.any()
    counts = np.zeros(601)
    counts[100:400] = 7.0
    np.testing.assert_array_equal(detilt_pool(counts, pi), counts)


# ---------------------------------------------------------------------------
# ⭐ THE FALSIFICATION: does the correction actually recover the library distribution?
# ---------------------------------------------------------------------------


def _enumerate_library(transcripts, theta, pmf):
    """Every placement, enumerated: ``(library, pool)`` histograms over fragment length.

    A fragment of length ``w`` is drawn with probability ``pmf[w]`` and placed uniformly over the
    start positions available to it on a transcript chosen with probability proportional to
    ``theta``. ``library`` counts every placement; ``pool`` counts only those crossing a junction.
    This is the generative model the opportunity function claims to invert, written out by hand.
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


def test_the_correction_RECOVERS_the_library_distribution_from_an_ENUMERATED_pool():
    """⭐⭐ The gate TRAPS: divide-by-a-probability exists to pass, with nothing tuned and no tolerance beyond float.

    Build a transcriptome and a fragment-length distribution, enumerate every placement to get both
    the library histogram and the junction-crossing subset of it, then de-tilt the subset. It must
    come back to the library histogram — shape-exactly, because the de-tilt preserves the total.
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

    # ⚠ and the raw pool is NOT already the answer, or the test would pass with no correction at all
    assert not np.allclose(pool, want, rtol=1e-3)


def test_the_correction_moves_the_pool_MEAN_onto_the_librarys():
    """The same statement in the one number the calibrator actually consumes."""
    transcripts = [[30, 40], [12] * 6, [70], [25, 8, 61], [100, 3]]
    theta = [5.0, 2.0, 9.0, 1.0, 4.0]
    pmf = np.zeros(91)
    pmf[10:81] = 1.0 / 71.0
    library, pool = _enumerate_library(transcripts, theta, pmf)
    pi = crossing_probability(*_csr(transcripts), np.array(theta), len(pmf) - 1)

    def mean(h):
        return float(np.dot(np.arange(h.size), h) / h.sum())

    assert mean(pool) > mean(library)  # the tilt exists, and it is long
    assert mean(detilt_pool(pool, pi)) == pytest.approx(mean(library), rel=1e-12)


def test_the_ratio_form_is_what_makes_a_WRONG_theta_SAFE():
    """⛔ The derivation's worked examples divide by ``A``; production divides by ``A / T``.

    With a deliberately wrong theta — all the weight on the transcript with the shortest exons, which
    has the steepest tilt — the ``A``-only form OVERSHOOTS past the library mean, i.e. it is worse
    than not correcting. The ratio form must not, because a reweighting moves ``A`` and ``T``
    together.
    """
    transcripts = [[30, 40], [6] * 12, [70], [25, 8, 61], [100, 3]]
    theta = [5.0, 2.0, 9.0, 1.0, 4.0]
    pmf = np.zeros(91)
    pmf[10:81] = 1.0 / 71.0
    library, pool = _enumerate_library(transcripts, theta, pmf)

    def mean(h):
        return float(np.dot(np.arange(h.size), h) / h.sum())

    wrong = np.array([0.0, 1.0, 0.0, 0.0, 0.0])  # every copy is the many-short-exon transcript
    lengths, offsets = _csr(transcripts)
    pi_wrong = crossing_probability(lengths, offsets, wrong, len(pmf) - 1)
    _, a_wrong = junction_opportunity(lengths, offsets, wrong, len(pmf) - 1)

    ratio_err = abs(mean(detilt_pool(pool, pi_wrong)) - mean(library))
    a_only_err = abs(mean(detilt_pool(pool, a_wrong / max(a_wrong.max(), 1.0))) - mean(library))
    raw_err = abs(mean(pool) - mean(library))

    assert ratio_err < raw_err, (ratio_err, raw_err)
    assert a_only_err > raw_err, (a_only_err, raw_err)


# ---------------------------------------------------------------------------
# the de-tilt's own contract
# ---------------------------------------------------------------------------


def test_the_detilt_preserves_the_EVIDENCE_COUNT():
    """⚠ It changes the pool's SHAPE, never how much evidence the pool represents.

    ``build_fl_models`` shrinks each pool toward the anchor with a Dirichlet pseudo-count, and the
    strength of that shrink is the pool total. A shape correction that also inflated the total would
    silently weaken the shrinkage as a side effect — two things varied, one of them by accident.
    """
    pi = crossing_probability(*_csr([[40, 60], [10] * 10]), np.array([3.0, 1.0]), 120)
    counts = np.zeros(121)
    counts[20:110] = np.arange(90, dtype=np.float64)
    out = detilt_pool(counts, pi)
    assert out.sum() == pytest.approx(counts.sum())
    assert not np.allclose(out, counts)


def test_the_detilt_cannot_INVENT_mass_where_the_pool_has_none():
    pi = crossing_probability(*_csr([[40, 60]]), np.array([1.0]), 120)
    counts = np.zeros(121)
    counts[50] = 100.0
    out = detilt_pool(counts, pi)
    assert out[50] == pytest.approx(100.0)
    assert not out[np.arange(121) != 50].any()


# ---------------------------------------------------------------------------
# the index adapter
# ---------------------------------------------------------------------------


def test_the_index_adapter_uses_the_REAL_transcripts_and_only_those(mini_index):
    """⛔ ``~is_synthetic``, alone. The manufactured nascent spans are not molecules anyone sequenced,
    and weighting them would put opportunity on exon structures the library does not contain.

    ⚠ NOT ``~is_synthetic & ~is_nrna``: on a real row ``is_nrna`` means "single-exon, so mature is
    nascent", and using it as a realness filter deletes real transcripts (trap 3).
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

    # and it is not the answer you get by weighting everything, or the filter is doing nothing
    assert int((~real).sum()) > 0, "MINI_GTF grew no synthetic rows — this test proves nothing"
    every = []
    for i in range(len(mini_index.t_df)):
        intervals = mini_index.get_exon_intervals(int(i))
        if intervals is not None and len(intervals):
            every.append([int(b - a) for a, b in intervals])
    assert not np.allclose(pi, crossing_probability(*_csr(every), np.ones(len(every)), max_length))


def test_the_index_adapter_agrees_with_the_ENUMERATING_oracle_on_a_real_index(mini_index):
    """One more rung: the adapter, the kernel and a hand enumeration all on the same three
    transcripts, so a wrong exon-length unpacking cannot hide behind the kernel."""
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


def test_build_fl_models_WITHOUT_the_divisor_is_unchanged():
    """``None`` means "no annotation was offered", and then the pool is what it always was.

    That is what makes the correction one thing varied rather than a rewrite of the EB path.
    """
    from types import SimpleNamespace

    from rigel.calibration.fl import build_fl_models
    from rigel.scan_payload import N_FRAGMENT_POOLS, POOL_DNA_INTERGENIC, POOL_RNA_SPLICED

    pools = np.zeros((N_FRAGMENT_POOLS, 51), dtype=np.int64)
    pools[POOL_RNA_SPLICED, 20] = 300
    pools[POOL_RNA_SPLICED, 40] = 100
    pools[POOL_DNA_INTERGENIC, 30] = 500
    anchor = np.zeros(51, dtype=np.int64)
    anchor[25] = 900
    payload = SimpleNamespace(pool_lengths=pools, deposited_lengths=anchor, max_length=50)

    plain = build_fl_models(payload)
    explicit = build_fl_models(payload, junction_opportunity=None)
    np.testing.assert_array_equal(plain.rna_pmf, explicit.rna_pmf)
    np.testing.assert_array_equal(plain.gdna_pmf, explicit.gdna_pmf)

    pi = np.zeros(51)
    pi[1:] = np.linspace(0.1, 0.9, 50)
    tilted = build_fl_models(payload, junction_opportunity=pi)
    assert not np.allclose(tilted.rna_pmf, plain.rna_pmf), "the divisor reached nothing"
    (
        np.testing.assert_array_equal(tilted.gdna_pmf, plain.gdna_pmf),
        "the junction divisor must not touch the gDNA pool",
    )
    assert tilted.n_rna == pytest.approx(plain.n_rna), "the EB evidence weight moved"


def test_EVERY_production_caller_of_build_fl_models_passes_the_divisor():
    """⛔ The divisor is optional so tests without an annotation can still build a model — which
    means forgetting it in production is silent, and the pool goes back to being tilted with nothing
    to show for it. Pin the call sites, the way the anchor's own frame is pinned.

    ⚠ Source-level on purpose: a runtime check would need a full pipeline run per call site, and the
    failure this guards against is somebody adding a fifth caller.
    """
    import ast
    import inspect

    from rigel import pipeline, scan_cache

    for module in (pipeline, scan_cache):
        tree = ast.parse(inspect.getsource(module))
        calls = [
            region
            for region in ast.walk(tree)
            if isinstance(region, ast.Call)
            and isinstance(region.func, ast.Name)
            and region.func.id == "build_fl_models"
        ]
        assert calls, f"{module.__name__} no longer calls build_fl_models — retarget this test"
        for call in calls:
            assert any(k.arg == "junction_opportunity" for k in call.keywords), (
                f"{module.__name__}:{call.lineno} builds FL models without the junction divisor, "
                "so the RNA pool stays tilted long"
            )
