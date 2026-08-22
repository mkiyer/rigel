"""`CaptureSampler.partition_array` against brute-force enumeration over every start position.

    TODO item 2   ·   Ledger: "the capture arm is O(transcripts x distinct fragment lengths)"

⛔ **WHY THIS EXISTS.** The capture partition path is 96.6 % of a capture-on simulation (profiled: 123.6 s
of 127.9 s, 1,052,172 `partition` calls) and its cache grew to 18 GB. It is being made fast. The standing
rule is that every divisor must be **unit-tested against brute-force enumeration**, so the optimisation is
gated on equality with a reference that literally loops over start positions and takes the max.

⚠ The oracle here is deliberately naive and slow: an explicit Python loop, no numpy, no shared helper with
the implementation. — a validator that calls the builder's own helper
validates nothing.

The matrix covers what the real panel does NOT, on purpose, because the optimisation must not silently
encode a property of one panel:

* probes of **unequal length** (the chr22 panel is uniformly 120 bp);
* probes whose start-ranges **do** and **do not** interact (measured: 0 % interact at w = 50, 97 % at
  w = 800, so both regimes are live);
* `min_overlap > 1`, which zeroes part of each trapezoid;
* `w > template length`, `w == template length`, probes flush against either boundary;
* a transcript with **no** probes, and a pool where **no** transcript has probes.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.sim.capture import CaptureConfig
from rigel.sim.capture.sampler import CaptureSampler, WeightedInterval


def brute_force_partition(
    intervals: list[WeightedInterval],
    seq_len: int,
    frag_len: int,
    *,
    off_target_weight: float,
    binding_per_base: float,
    min_overlap: int,
) -> float:
    """The definition, spelled out: for every start, the best single probe GROUP, summed.

    Groups are summed within, then maxed across — a probe split across exons contributes the sum of its
    pieces, and a fragment is credited by its best single probe.
    """
    eff_len = seq_len - frag_len + 1
    if eff_len <= 0:
        return 0.0
    total = off_target_weight * eff_len
    for start in range(eff_len):
        per_group: dict[int, float] = {}
        for interval in intervals:
            overlap = min(start + frag_len, interval.end) - max(start, interval.start)
            if overlap < min_overlap or overlap <= 0:
                continue
            per_group[interval.probe_group] = (
                per_group.get(interval.probe_group, 0.0) + overlap * interval.scale
            )
        if per_group:
            total += binding_per_base * max(per_group.values())
    return total


def make_sampler(
    intervals_by_key: dict[int, list[WeightedInterval]], **overrides
) -> CaptureSampler:
    """A sampler whose mRNA intervals are injected directly, bypassing probe-file parsing."""
    config = CaptureConfig(
        probes="unused",
        probe_format="transcript",
        off_target_weight=overrides.get("off_target_weight", 1.0),
        binding_per_base=overrides.get("binding_per_base", 10.0),
        gdna_split_penalty=0.2,
        min_overlap=overrides.get("min_overlap", 1),
    )
    # ⛔⛔ **BUILT THROUGH THE REAL `__init__`, NEVER `__new__`** — only the probe intervals are then
    # injected, which is the one thing this fixture exists to control. It used to bypass the
    # constructor and hand-set four attributes, so any state `__init__` gained afterwards was simply
    # ABSENT here: on 2026-08-19 a new per-width partition memo landed and **242 of this file's tests
    # failed with `AttributeError`** while the sampler worked perfectly in production. That is
    # `TRAPS: a-gate-that-reconstructs` — a fixture that rebuilds its subject tests the rebuild.
    sampler = CaptureSampler(config, transcripts=(), ref_lengths={}, enabled=True)
    sampler._mrna_intervals = intervals_by_key
    return sampler


def iv(start: int, end: int, *, group: int = 0, scale: float = 1.0) -> WeightedInterval:
    return WeightedInterval(start=start, end=end, scale=scale, probe_group=group)


# ── the layouts. Each is a named shape the optimisation could get wrong differently ────────────────
LAYOUTS: dict[str, list[WeightedInterval]] = {
    "single probe, mid-template": [iv(300, 420)],
    "single probe, flush at start": [iv(0, 120)],
    "single probe, flush at end": [iv(880, 1000)],
    "two probes, far apart": [iv(100, 220, group=0), iv(700, 820, group=1)],
    "two probes, adjacent (gap 2)": [iv(100, 220, group=0), iv(222, 342, group=1)],
    "two probes, UNEQUAL length": [iv(100, 160, group=0), iv(300, 700, group=1)],
    "two probes, unequal SCALE": [
        iv(100, 220, group=0, scale=0.2),
        iv(300, 420, group=1, scale=1.0),
    ],
    "three probes, dense": [
        iv(100, 220, group=0),
        iv(240, 360, group=1),
        iv(380, 500, group=2),
    ],
    "one probe SPLIT across two pieces (same group)": [
        iv(100, 160, group=7),
        iv(400, 460, group=7),
    ],
    "no probes": [],
}

#: ⭐ Every layout goes through the batched path, INCLUDING the split probe group. It used to fall back
#: to the per-key loop, which is why the gDNA space — where a probe spanning an intron always splits
#: into several genomic pieces — never got the fast path and stayed at 268 s / 37 GB.
SPLIT_GROUP_LAYOUT = "one probe SPLIT across two pieces (same group)"
BATCHABLE_LAYOUTS = list(LAYOUTS)

FRAGMENT_LENGTHS = [1, 2, 50, 119, 120, 121, 200, 500, 999, 1000, 1001, 2000]
SEQ_LEN = 1000


def assert_batched_path_is_live(sampler: CaptureSampler) -> None:
    """⛔ WITHOUT THIS, THE WHOLE FILE TESTS NOTHING.

    `_flat_probes` returns `None` — falling back to the old per-key loop for the ENTIRE space — as soon
    as one probe group holds more than one interval. The first version of this file put a split-group
    layout in the shared pool, so every case silently ran the fallback: **seven of eight deliberate
    perturbations to the batched code passed**. Assert the path under test is the path being taken.
    """
    assert sampler._flat_probes("mrna") is not None, (
        "the batched path is NOT live for this sampler, so this test is exercising the fallback"
    )


class TestAgainstBruteForce:
    """Every case goes through `partition_array` — the batched path — not the scalar one."""

    @pytest.mark.parametrize("layout_name", BATCHABLE_LAYOUTS)
    @pytest.mark.parametrize("frag_len", FRAGMENT_LENGTHS)
    def test_partition_array_matches_enumeration(self, layout_name: str, frag_len: int) -> None:
        intervals = LAYOUTS[layout_name]
        sampler = make_sampler({0: intervals} if intervals else {})
        if intervals:
            assert_batched_path_is_live(sampler)
        expected = brute_force_partition(
            intervals,
            SEQ_LEN,
            frag_len,
            off_target_weight=1.0,
            binding_per_base=10.0,
            min_overlap=1,
        )
        lengths = np.array([SEQ_LEN], dtype=np.int64)
        actual = float(sampler.partition_array("mrna", [0], lengths, frag_len)[0])
        assert actual == pytest.approx(expected, rel=1e-12, abs=1e-9), (
            f"{layout_name} at w={frag_len}: got {actual}, enumeration says {expected}"
        )

    @pytest.mark.parametrize("layout_name", list(LAYOUTS))
    @pytest.mark.parametrize("frag_len", FRAGMENT_LENGTHS)
    def test_the_scalar_path_also_matches_enumeration(
        self, layout_name: str, frag_len: int
    ) -> None:
        """The per-key fallback is still reachable (split probe groups), so it is gated too."""
        intervals = LAYOUTS[layout_name]
        sampler = make_sampler({0: intervals} if intervals else {})
        expected = brute_force_partition(
            intervals,
            SEQ_LEN,
            frag_len,
            off_target_weight=1.0,
            binding_per_base=10.0,
            min_overlap=1,
        )
        assert sampler.partition("mrna", 0, SEQ_LEN, frag_len) == pytest.approx(
            expected, rel=1e-12, abs=1e-9
        )

    @pytest.mark.parametrize("min_overlap", [1, 5, 60, 121])
    @pytest.mark.parametrize("frag_len", [50, 200, 500])
    def test_min_overlap_gate_matches_enumeration(self, min_overlap: int, frag_len: int) -> None:
        intervals = LAYOUTS["three probes, dense"]
        sampler = make_sampler({0: intervals}, min_overlap=min_overlap)
        assert_batched_path_is_live(sampler)
        expected = brute_force_partition(
            intervals,
            SEQ_LEN,
            frag_len,
            off_target_weight=1.0,
            binding_per_base=10.0,
            min_overlap=min_overlap,
        )
        lengths = np.array([SEQ_LEN], dtype=np.int64)
        actual = float(sampler.partition_array("mrna", [0], lengths, frag_len)[0])
        assert actual == pytest.approx(expected, rel=1e-12, abs=1e-9)

    def test_a_multi_key_pool_matches_enumeration_key_by_key(self) -> None:
        """⚠ The batched path shares ONE buffer across keys. A key writing outside its own slice, or a
        rank colliding between keys, only shows up with several keys of DIFFERENT lengths at once."""
        names = list(BATCHABLE_LAYOUTS)
        intervals_by_key = {i: LAYOUTS[name] for i, name in enumerate(names)}
        sampler = make_sampler(intervals_by_key)
        assert_batched_path_is_live(sampler)
        lengths = np.array([SEQ_LEN + 137 * i for i in range(len(names))], dtype=np.int64)

        for frag_len in (50, 200, 500, 1200):
            actual = sampler.partition_array("mrna", list(range(len(names))), lengths, frag_len)
            for key, name in enumerate(names):
                expected = brute_force_partition(
                    LAYOUTS[name],
                    int(lengths[key]),
                    frag_len,
                    off_target_weight=1.0,
                    binding_per_base=10.0,
                    min_overlap=1,
                )
                assert actual[key] == pytest.approx(expected, rel=1e-12, abs=1e-9), (
                    f"key {key} ({name}) at w={frag_len}"
                )

    def test_a_split_probe_group_is_batched_and_still_summed_within_the_group(self) -> None:
        """A probe split across exons must have its pieces SUMMED before the max across probes.

        ⚠ This is the case the gDNA space is made of — every probe that spans an intron splits when
        projected to genomic coordinates — so it must be on the FAST path, not a fallback.
        """
        intervals = LAYOUTS[SPLIT_GROUP_LAYOUT]
        sampler = make_sampler({0: intervals})
        assert_batched_path_is_live(sampler)
        lengths = np.array([SEQ_LEN], dtype=np.int64)
        for frag_len in (50, 200, 500):
            expected = brute_force_partition(
                intervals,
                SEQ_LEN,
                frag_len,
                off_target_weight=1.0,
                binding_per_base=10.0,
                min_overlap=1,
            )
            assert float(
                sampler.partition_array("mrna", [0], lengths, frag_len)[0]
            ) == pytest.approx(expected, rel=1e-12, abs=1e-9)


class TestPartitionArray:
    """`partition_array` must agree with `partition` element for element — it is the batched form."""

    @pytest.mark.parametrize("frag_len", [50, 200, 500, 1200])
    def test_array_agrees_with_the_scalar_over_a_mixed_pool(self, frag_len: int) -> None:
        names = list(BATCHABLE_LAYOUTS)
        intervals_by_key = {i: LAYOUTS[name] for i, name in enumerate(names)}
        intervals_by_key[len(names)] = []  # a transcript with no probes, in the middle of the pool
        lengths = np.array([SEQ_LEN + 37 * i for i in range(len(intervals_by_key))], dtype=np.int64)
        sampler = make_sampler(intervals_by_key)
        assert_batched_path_is_live(sampler)

        actual = sampler.partition_array("mrna", range(len(lengths)), lengths, frag_len)
        expected = np.array(
            [
                sampler.partition("mrna", key, int(lengths[key]), frag_len)
                for key in range(len(lengths))
            ]
        )
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-9)

    def test_a_key_absent_from_the_call_does_not_leak_into_another(self) -> None:
        """`partition_array` may be called with a SUBSET of keys — the live-abundance filter does
        exactly that. A key left out must not shift the buffer slices of the keys left in."""
        intervals_by_key = {
            0: LAYOUTS["three probes, dense"],
            1: LAYOUTS["two probes, far apart"],
            2: LAYOUTS["single probe, mid-template"],
        }
        sampler = make_sampler(intervals_by_key)
        assert_batched_path_is_live(sampler)
        lengths_all = np.array([SEQ_LEN, SEQ_LEN + 500, SEQ_LEN + 900], dtype=np.int64)
        full = sampler.partition_array("mrna", [0, 1, 2], lengths_all, 200)
        subset = sampler.partition_array("mrna", [0, 2], lengths_all[[0, 2]], 200)
        np.testing.assert_allclose(subset, full[[0, 2]], rtol=1e-12, atol=1e-9)

    def test_a_pool_where_nothing_is_on_panel_is_pure_off_target(self) -> None:
        sampler = make_sampler({})
        lengths = np.array([1000, 500, 250], dtype=np.int64)
        actual = sampler.partition_array("mrna", range(3), lengths, 200)
        np.testing.assert_array_equal(actual, lengths - 200 + 1)

    def test_templates_shorter_than_the_fragment_contribute_nothing(self) -> None:
        sampler = make_sampler({0: LAYOUTS["three probes, dense"]})
        lengths = np.array([150], dtype=np.int64)
        assert sampler.partition_array("mrna", range(1), lengths, 200)[0] == 0.0


class TestNoUnboundedCache:
    """⚠ The 18 GB (RNA) and 38 GB (gDNA) were `_mass_cache`, and it never hit.

    `partition_array` visits each (key, fragment length) pair exactly once, and `sample_starts` is
    driven by a counts dict keyed by exactly that pair — so every entry either path stored was written
    and never read back. The cache is gone; this pins that no new per-call state replaces it.
    """

    def test_no_per_call_state_accumulates_across_fragment_lengths(self) -> None:
        """⚠ **ONE dict is now ALLOWED to hold one entry per WIDTH, and the rule that replaces "no
        growth" is stricter about the thing that actually cost 38 GB** (2026-08-19). `_mass_cache` grew
        per call and was NEVER READ BACK; `_partition_memo` is read back by every later condition
        sharing the capture panel (measured 260 s -> 0.2 s), and it is bounded because a different
        template population CLEARS it. So: every OTHER dict must still not grow, and the memo must be
        bounded by the number of distinct widths — never by the number of calls.
        """
        sampler = make_sampler({i: LAYOUTS["three probes, dense"] for i in range(40)})
        lengths = np.full(40, SEQ_LEN, dtype=np.int64)
        sizes_before = {k: len(v) for k, v in vars(sampler).items() if isinstance(v, dict)}
        widths = list(range(50, 350))
        for frag_len in widths:
            sampler.partition_array("mrna", range(40), lengths, frag_len)
        for name, size in sizes_before.items():
            if name == "_partition_memo":
                continue
            grown = len(getattr(sampler, name)) - size
            assert grown <= 1, (
                f"{name} grew by {grown} entries across 300 fragment lengths. Anything that scales "
                f"with fragment length here is the 38 GB coming back."
            )
        assert len(sampler._partition_memo) == len(widths), (
            "the memo holds exactly one entry per width"
        )

        # ⭐ REPEATING the same sweep must add NOTHING — that is the difference from `_mass_cache`
        for frag_len in widths:
            sampler.partition_array("mrna", range(40), lengths, frag_len)
        assert len(sampler._partition_memo) == len(widths), "a repeat sweep must be pure cache hits"

        # ⛔ and a DIFFERENT template population must not accumulate beside the first
        other = np.full(40, SEQ_LEN + 1, dtype=np.int64)
        for frag_len in widths:
            sampler.partition_array("mrna", range(40), other, frag_len)
        assert len(sampler._partition_memo) == len(widths), (
            "a second population must CLEAR the memo, not accumulate beside it — that is the 38 GB"
        )

    def test_the_memo_returns_what_recomputation_would(self) -> None:
        """⛔ A cache that returns a stale or aliased vector is worse than no cache. Every width is
        compared against a freshly-built sampler's answer, and the returned array must be a COPY (a
        caller mutating its result must not poison the memo)."""
        layout = {i: LAYOUTS["three probes, dense"] for i in range(4)}
        cached_sampler = make_sampler(layout)
        lengths = np.full(4, SEQ_LEN, dtype=np.int64)
        for frag_len in (60, 120, 121, 200):
            first = cached_sampler.partition_array("mrna", range(4), lengths, frag_len)
            first[:] = -12345.0  # a caller mutating its result
            again = cached_sampler.partition_array("mrna", range(4), lengths, frag_len)
            fresh = make_sampler(layout).partition_array("mrna", range(4), lengths, frag_len)
            np.testing.assert_allclose(again, fresh, err_msg=f"memo differs at width {frag_len}")

    def test_sample_starts_does_not_accumulate_state_either(self) -> None:
        sampler = make_sampler({0: LAYOUTS["three probes, dense"]})
        rng = np.random.default_rng(0)
        sizes_before = {k: len(v) for k, v in vars(sampler).items() if isinstance(v, dict)}
        for frag_len in range(50, 350):
            sampler.sample_starts("mrna", 0, SEQ_LEN, frag_len, 3, rng)
        for name, size in sizes_before.items():
            grown = len(getattr(sampler, name)) - size
            assert grown <= 1, f"{name} grew by {grown} entries across 300 fragment lengths"
