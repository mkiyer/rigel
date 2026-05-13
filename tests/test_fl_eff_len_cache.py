"""Regression tests for ``FragmentLengthModel._build_eff_len_cache``
memoization (Tier-1 perf fix).

The cumulative CDF / first-moment arrays consumed by
``compute_all_transcript_eff_lens`` are expensive to build (one
``np.cumsum`` over the full PMF per call) but constant after
``finalize()``.  These tests pin the caching contract:

* finalized model returns the exact same array object on repeat calls
* caller-visible values are bit-identical to a fresh recompute
* re-finalizing invalidates the cache (mass changes ⇒ new arrays)
* pre-finalize behaviour is unchanged (always recomputes)
"""

from __future__ import annotations

import numpy as np

from rigel.frag_length_model import FragmentLengthModel


def _make_model(seed: int = 7, max_size: int = 200) -> FragmentLengthModel:
    rng = np.random.default_rng(seed)
    counts = rng.integers(0, 100, size=max_size + 1).astype(np.float64)
    return FragmentLengthModel.from_counts(counts, max_size=max_size)


def test_finalized_cache_returns_same_object():
    fl = _make_model()
    cdf1, cmom1 = fl._build_eff_len_cache()
    cdf2, cmom2 = fl._build_eff_len_cache()
    assert cdf1 is cdf2, "CDF cache must be reused after finalize()"
    assert cmom1 is cmom2, "CMOM cache must be reused after finalize()"


def test_cache_matches_fresh_recompute_bit_for_bit():
    fl = _make_model()
    cdf_cached, cmom_cached = fl._build_eff_len_cache()

    # Fresh recompute via the same private helper after a no-op cache invalidation.
    fl._cdf_cache = None
    fl._cmom_cache = None
    cdf_fresh, cmom_fresh = fl._build_eff_len_cache()

    np.testing.assert_array_equal(cdf_cached, cdf_fresh)
    np.testing.assert_array_equal(cmom_cached, cmom_fresh)


def test_compute_all_transcript_eff_lens_unchanged_by_cache():
    fl = _make_model()
    lens = np.array([10, 50, 100, 250, 500], dtype=np.int64)

    eff_first = fl.compute_all_transcript_eff_lens(lens)
    eff_second = fl.compute_all_transcript_eff_lens(lens)
    eff_third = fl.compute_all_transcript_eff_lens(lens)

    np.testing.assert_array_equal(eff_first, eff_second)
    np.testing.assert_array_equal(eff_second, eff_third)


def test_refinalize_invalidates_cache():
    fl = _make_model()
    cdf_before, _ = fl._build_eff_len_cache()
    # Mutate the histogram and re-finalize.
    fl.counts[10] += 500.0
    fl._total_weight = float(fl.counts.sum())
    fl.finalize()
    cdf_after, _ = fl._build_eff_len_cache()
    assert cdf_after is not cdf_before, "finalize() must invalidate the cache"
    # Sanity: arrays differ at/after the mutated bin.
    assert not np.array_equal(cdf_before, cdf_after)


def test_pre_finalize_does_not_cache():
    fl = FragmentLengthModel(max_size=100)
    fl.observe(50)
    assert not fl._finalized
    cdf1, _ = fl._build_eff_len_cache()
    fl.observe(50)  # mutate before finalize
    cdf2, _ = fl._build_eff_len_cache()
    assert fl._cdf_cache is None, "Pre-finalize calls must not populate the cache"
    # Distributions differ ⇒ CDFs differ.
    assert not np.array_equal(cdf1, cdf2)
