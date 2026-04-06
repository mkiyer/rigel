"""Tests for partition scatter functions and partition-native EM equivalence."""

import numpy as np

from rigel.scored_fragments import Locus, ScoredFragments
from rigel.partition import partition_and_free
from rigel.native import (
    build_partition_offsets,
    scatter_candidates_f64,
    scatter_candidates_i32,
    scatter_candidates_u8,
    scatter_units_f64,
    scatter_units_i32,
    scatter_units_u8,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_scored_fragments(n_units, candidates_per_unit, n_transcripts, rng):
    """Create a small ScoredFragments with random data."""
    offsets = np.zeros(n_units + 1, dtype=np.int64)
    for i in range(n_units):
        offsets[i + 1] = offsets[i] + candidates_per_unit[i]
    n_cand = int(offsets[-1])

    return ScoredFragments(
        offsets=offsets,
        t_indices=rng.integers(0, n_transcripts, size=n_cand, dtype=np.int32),
        log_liks=rng.standard_normal(n_cand),
        count_cols=rng.integers(0, 8, size=n_cand, dtype=np.uint8),
        coverage_weights=rng.random(n_cand),
        tx_starts=rng.integers(0, 500, size=n_cand, dtype=np.int32),
        tx_ends=rng.integers(500, 1000, size=n_cand, dtype=np.int32),
        locus_t_indices=rng.integers(0, n_transcripts, size=n_units, dtype=np.int32),
        locus_count_cols=rng.integers(0, 8, size=n_units, dtype=np.uint8),
        is_spliced=rng.choice([True, False], size=n_units),
        gdna_log_liks=rng.standard_normal(n_units),
        genomic_footprints=rng.integers(100, 500, size=n_units, dtype=np.int32),
        frag_ids=np.arange(n_units, dtype=np.int64),
        frag_class=np.zeros(n_units, dtype=np.int8),
        splice_type=np.zeros(n_units, dtype=np.uint8),
        n_units=n_units,
        n_candidates=n_cand,
    )


def _make_loci(n_units, n_transcripts_per_locus, rng):
    """Create loci that partition all units and transcripts."""
    # Simple 2-locus case: split units and transcripts in half
    n_loci = len(n_transcripts_per_locus)

    # Assign units round-robin to loci
    unit_assignments = np.arange(n_units) % n_loci

    # Build transcript indices per locus
    t_offset = 0
    loci = []
    for li in range(n_loci):
        n_t = n_transcripts_per_locus[li]
        t_indices = np.arange(t_offset, t_offset + n_t, dtype=np.int32)
        u_indices = np.where(unit_assignments == li)[0].astype(np.int32)
        loci.append(Locus(
            locus_id=li,
            transcript_indices=t_indices,
            unit_indices=u_indices,
            gdna_span=10000,
            merged_intervals=[("chr1", 0, 10000)],
        ))
        t_offset += n_t

    return loci


# ---------------------------------------------------------------------------
# Test build_partition_offsets
# ---------------------------------------------------------------------------


class TestBuildPartitionOffsets:
    def test_basic(self):
        """Verify partition offsets match manual calculation."""
        # 5 units, 2 loci
        g_offsets = np.array([0, 3, 6, 8, 10, 15], dtype=np.int64)
        locus_units = [
            np.array([0, 2, 4], dtype=np.int32),  # locus 0: units 0,2,4
            np.array([1, 3], dtype=np.int32),       # locus 1: units 1,3
        ]

        result = build_partition_offsets(g_offsets, locus_units, 2)

        assert len(result) == 2

        # Locus 0: unit 0 has 3 cands, unit 2 has 2 cands, unit 4 has 5 cands
        np.testing.assert_array_equal(result[0], [0, 3, 5, 10])

        # Locus 1: unit 1 has 3 cands, unit 3 has 2 cands
        np.testing.assert_array_equal(result[1], [0, 3, 5])

    def test_empty_locus(self):
        """Locus with zero units produces [0] offsets."""
        g_offsets = np.array([0, 5], dtype=np.int64)
        locus_units = [
            np.array([0], dtype=np.int32),
            np.array([], dtype=np.int32),
        ]
        result = build_partition_offsets(g_offsets, locus_units, 2)
        np.testing.assert_array_equal(result[0], [0, 5])
        np.testing.assert_array_equal(result[1], [0])


# ---------------------------------------------------------------------------
# Test scatter functions
# ---------------------------------------------------------------------------


class TestScatterCandidates:
    def test_round_trip_f64(self):
        """Scattered candidate arrays match original segments."""
        rng = np.random.default_rng(42)
        g_offsets = np.array([0, 3, 6, 10], dtype=np.int64)
        data = rng.standard_normal(10)

        locus_units = [
            np.array([0, 2], dtype=np.int32),
            np.array([1], dtype=np.int32),
        ]
        p_offsets = build_partition_offsets(g_offsets, locus_units, 2)

        result = scatter_candidates_f64(data, g_offsets, locus_units, p_offsets, 2)

        # Locus 0: units 0 and 2 → data[0:3] + data[6:10]
        expected_0 = np.concatenate([data[0:3], data[6:10]])
        np.testing.assert_array_equal(result[0], expected_0)

        # Locus 1: unit 1 → data[3:6]
        np.testing.assert_array_equal(result[1], data[3:6])

    def test_round_trip_i32(self):
        """int32 scatter preserves exact values."""
        g_offsets = np.array([0, 2, 5], dtype=np.int64)
        data = np.array([10, 20, 30, 40, 50], dtype=np.int32)

        locus_units = [np.array([1], dtype=np.int32),
                       np.array([0], dtype=np.int32)]
        p_offsets = build_partition_offsets(g_offsets, locus_units, 2)

        result = scatter_candidates_i32(data, g_offsets, locus_units, p_offsets, 2)
        np.testing.assert_array_equal(result[0], [30, 40, 50])
        np.testing.assert_array_equal(result[1], [10, 20])

    def test_round_trip_u8(self):
        """uint8 scatter preserves exact values."""
        g_offsets = np.array([0, 1, 3], dtype=np.int64)
        data = np.array([100, 200, 255], dtype=np.uint8)

        locus_units = [np.array([0, 1], dtype=np.int32)]
        p_offsets = build_partition_offsets(g_offsets, locus_units, 1)

        result = scatter_candidates_u8(data, g_offsets, locus_units, p_offsets, 1)
        np.testing.assert_array_equal(result[0], [100, 200, 255])


class TestScatterUnits:
    def test_round_trip_f64(self):
        """Scattered unit arrays match gathered elements."""
        data = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        locus_units = [
            np.array([0, 3, 4], dtype=np.int32),
            np.array([1, 2], dtype=np.int32),
        ]

        result = scatter_units_f64(data, locus_units, 2)
        np.testing.assert_array_equal(result[0], [1.0, 4.0, 5.0])
        np.testing.assert_array_equal(result[1], [2.0, 3.0])

    def test_round_trip_i32(self):
        data = np.array([10, 20, 30], dtype=np.int32)
        locus_units = [np.array([2, 0], dtype=np.int32)]

        result = scatter_units_i32(data, locus_units, 1)
        np.testing.assert_array_equal(result[0], [30, 10])

    def test_round_trip_u8(self):
        data = np.array([1, 2, 3, 4], dtype=np.uint8)
        locus_units = [
            np.array([0, 2], dtype=np.int32),
            np.array([1, 3], dtype=np.int32),
        ]

        result = scatter_units_u8(data, locus_units, 2)
        np.testing.assert_array_equal(result[0], [1, 3])
        np.testing.assert_array_equal(result[1], [2, 4])


# ---------------------------------------------------------------------------
# Test partition_and_free
# ---------------------------------------------------------------------------


class TestPartitionAndFree:
    def test_data_equivalence(self):
        """Verify partitioned data matches original global CSR slices."""
        rng = np.random.default_rng(123)
        n_units = 20
        n_transcripts = 6
        candidates_per_unit = rng.integers(1, 5, size=n_units)
        em_data = _make_scored_fragments(n_units, candidates_per_unit, n_transcripts, rng)
        loci = _make_loci(n_units, [3, 3], rng)

        # Keep copies of original data for verification
        orig_offsets = em_data.offsets.copy()
        orig_t_indices = em_data.t_indices.copy()
        orig_log_liks = em_data.log_liks.copy()
        orig_coverage_weights = em_data.coverage_weights.copy()
        orig_tx_starts = em_data.tx_starts.copy()
        orig_tx_ends = em_data.tx_ends.copy()
        orig_count_cols = em_data.count_cols.copy()
        orig_is_spliced = em_data.is_spliced.copy()
        orig_gdna_log_liks = em_data.gdna_log_liks.copy()
        orig_genomic_footprints = em_data.genomic_footprints.copy()
        orig_locus_t_indices = em_data.locus_t_indices.copy()
        orig_locus_count_cols = em_data.locus_count_cols.copy()

        partitions = partition_and_free(em_data, loci)

        # em_data arrays should be None
        assert em_data.offsets is None
        assert em_data.t_indices is None
        assert em_data.frag_ids is None

        assert len(partitions) == 2

        for li, locus in enumerate(loci):
            part = partitions[li]
            assert part.locus_id == li
            assert part.n_units == len(locus.unit_indices)

            # Verify per-unit arrays
            for k, u in enumerate(locus.unit_indices):
                assert part.is_spliced[k] == orig_is_spliced[u]
                np.testing.assert_equal(
                    part.gdna_log_liks[k], orig_gdna_log_liks[u])
                assert part.genomic_footprints[k] == orig_genomic_footprints[u]
                assert part.locus_t_indices[k] == orig_locus_t_indices[u]
                assert part.locus_count_cols[k] == orig_locus_count_cols[u]

            # Verify per-candidate arrays
            cand_pos = 0
            for k, u in enumerate(locus.unit_indices):
                g_start = orig_offsets[u]
                g_end = orig_offsets[u + 1]
                n_cand = g_end - g_start

                np.testing.assert_array_equal(
                    part.t_indices[cand_pos:cand_pos + n_cand],
                    orig_t_indices[g_start:g_end])
                np.testing.assert_array_equal(
                    part.log_liks[cand_pos:cand_pos + n_cand],
                    orig_log_liks[g_start:g_end])
                np.testing.assert_array_equal(
                    part.coverage_weights[cand_pos:cand_pos + n_cand],
                    orig_coverage_weights[g_start:g_end])
                np.testing.assert_array_equal(
                    part.tx_starts[cand_pos:cand_pos + n_cand],
                    orig_tx_starts[g_start:g_end])
                np.testing.assert_array_equal(
                    part.tx_ends[cand_pos:cand_pos + n_cand],
                    orig_tx_ends[g_start:g_end])
                np.testing.assert_array_equal(
                    part.count_cols[cand_pos:cand_pos + n_cand],
                    orig_count_cols[g_start:g_end])
                cand_pos += n_cand

            assert cand_pos == part.n_candidates

    def test_total_elements_preserved(self):
        """Total units and candidates across partitions equals original."""
        rng = np.random.default_rng(456)
        n_units = 30
        n_transcripts = 8
        candidates_per_unit = rng.integers(1, 6, size=n_units)
        em_data = _make_scored_fragments(n_units, candidates_per_unit, n_transcripts, rng)
        orig_n_units = em_data.n_units
        orig_n_cand = em_data.n_candidates
        loci = _make_loci(n_units, [3, 3, 2], rng)

        partitions = partition_and_free(em_data, loci)

        total_units = sum(p.n_units for p in partitions.values())
        total_cand = sum(p.n_candidates for p in partitions.values())
        assert total_units == orig_n_units
        assert total_cand == orig_n_cand


class TestNoDuplicateCandidatesPerUnit:
    """Verify that partition data preserves the scorer's no-duplicate invariant.

    The scorer emits at most one candidate per transcript per unit (fragment).
    The extract function in the EM solver relies on this invariant — it sorts
    candidates by local component index without deduplication. If duplicates
    existed, the equivalence-class builder would silently produce wrong results.
    """

    def test_unique_t_indices_per_unit(self):
        """Each unit in partitioned data has unique global transcript indices."""
        rng = np.random.default_rng(999)
        n_units = 50
        n_transcripts = 10

        # Build ScoredFragments with guaranteed unique candidates per unit:
        # each unit gets a random subset of transcripts (no repeats).
        offsets = np.zeros(n_units + 1, dtype=np.int64)
        all_t_indices = []
        for u in range(n_units):
            width = rng.integers(1, min(6, n_transcripts + 1))
            chosen = rng.choice(n_transcripts, size=width, replace=False)
            all_t_indices.extend(chosen)
            offsets[u + 1] = offsets[u] + width
        n_cand = int(offsets[-1])
        t_indices = np.array(all_t_indices, dtype=np.int32)

        em_data = ScoredFragments(
            offsets=offsets,
            t_indices=t_indices,
            log_liks=rng.standard_normal(n_cand),
            count_cols=rng.integers(0, 4, size=n_cand, dtype=np.uint8),
            coverage_weights=np.ones(n_cand),
            tx_starts=np.zeros(n_cand, dtype=np.int32),
            tx_ends=np.full(n_cand, 500, dtype=np.int32),
            locus_t_indices=rng.integers(0, n_transcripts, size=n_units, dtype=np.int32),
            locus_count_cols=np.zeros(n_units, dtype=np.uint8),
            is_spliced=rng.choice([True, False], size=n_units),
            gdna_log_liks=rng.standard_normal(n_units),
            genomic_footprints=rng.integers(100, 500, size=n_units, dtype=np.int32),
            frag_ids=np.arange(n_units, dtype=np.int64),
            frag_class=np.zeros(n_units, dtype=np.int8),
            splice_type=np.zeros(n_units, dtype=np.uint8),
            n_units=n_units,
            n_candidates=n_cand,
        )

        loci = _make_loci(n_units, [5, 5], rng)
        partitions = partition_and_free(em_data, loci)

        for li, part in partitions.items():
            for u in range(part.n_units):
                start = part.offsets[u]
                end = part.offsets[u + 1]
                unit_t = part.t_indices[start:end]
                assert len(unit_t) == len(np.unique(unit_t)), (
                    f"Locus {li}, unit {u}: duplicate transcript indices "
                    f"in partition data: {unit_t}"
                )
