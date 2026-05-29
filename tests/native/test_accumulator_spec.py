"""
Fractional accumulator specification tests.

Each test case has two halves:

1. `test_<case>_reference` — runs against the pure-Python reference
   in `_accumulator_reference.py`. Passes today; locks the spec
   so any change to the reference (or to expected values) requires
   updating this file.

2. `test_<case>_native` — runs against the native C++ implementation
   that lands in Phase 2. Marked `pytest.mark.xfail` until the
   binding is exposed via `rigel.native`. After Phase 3, the xfail
   marker is removed and the tests must pass byte-for-byte against
   the reference.

The 10 cases cover:
- T1  contained, single block, single region
- T2  contained, multi-block (spliced), single region
- T3  two blocks in adjacent regions (boundary mass + flux split)
- T4  two blocks in non-adjacent regions (user-verified §4.5.1)
- T5  three blocks all-adjacent (user-verified §4.5.2)
- T6  fully-spans-region (single block straddles one region)
- T7  mass conservation over 1000 random contained fragments
- T8  flux dedup on adjacent-region two-block fragment
- T10 negative strand attribution
- T11 spliced flag attribution across channels

FL histograms are NOT accumulated (see audit_phase1.md decision
#6 — FL is used downstream in EM only). No FL tests.
"""

from __future__ import annotations

import numpy as np
import pytest

from ._accumulator_reference import (
    Accumulator,
    channel_idx,
)

# Native is not yet implemented; this import is guarded so the
# reference tests can run regardless.
try:
    from rigel.native import Accumulator as NativeAccumulator  # type: ignore
    HAS_NATIVE = True
except (ImportError, AttributeError):
    NativeAccumulator = None  # type: ignore
    HAS_NATIVE = False


XFAIL_NATIVE = pytest.mark.xfail(
    not HAS_NATIVE,
    reason="fractional accumulator native binding lands in Phase 2/3",
    strict=True,
)


# --- Common partition fixtures -----------------------------------------------


def partition_exon_intron_exon() -> Accumulator:
    """Three regions: exon1 (1000,2000), intron (2000,5000), exon2 (5000,6000).
    Four boundaries at positions 1000, 2000, 5000, 6000.
    """
    edges = np.array([1000, 2000, 5000, 6000], dtype=np.int64)
    return Accumulator(region_edges=edges)


def partition_three_adjacent_exons() -> Accumulator:
    """Three contiguous regions: (0,100), (100,200), (200,400).
    Four boundaries at 0, 100, 200, 400.
    """
    edges = np.array([0, 100, 200, 400], dtype=np.int64)
    return Accumulator(region_edges=edges)


# --- T1: contained single block ---------------------------------------------


def test_t1_contained_single_block_reference():
    acc = partition_exon_intron_exon()
    acc.deposit(blocks=[(1100, 1200)], spliced=False, strand_pos=True)
    ch = channel_idx(spliced=False, strand_pos=True)
    assert acc.regions[0].contained[ch] == 1
    assert acc.regions[0].contained.sum() == 1
    assert acc.regions[1].contained.sum() == 0
    assert acc.regions[2].contained.sum() == 0
    for b in acc.boundaries:
        assert b.mass_left.sum() == 0.0
        assert b.mass_right.sum() == 0.0
        assert b.flux.sum() == 0
    assert acc.total_mass_deposited() == pytest.approx(1.0)


@XFAIL_NATIVE
def test_t1_contained_single_block_native():
    raise AssertionError("native impl not landed")


# --- T2: contained, multi-block spliced -------------------------------------


def test_t2_contained_multi_block_spliced_reference():
    """Two blocks both inside exon2 (5000,6000): (5100,5200) and (5400,5500).
    Spliced (CIGAR N between them, but both blocks live in the same region
    because the region is wider than the splice). Strand +.
    """
    acc = partition_exon_intron_exon()
    acc.deposit(
        blocks=[(5100, 5200), (5400, 5500)],
        spliced=True,
        strand_pos=True,
    )
    ch = channel_idx(spliced=True, strand_pos=True)
    assert acc.regions[2].contained[ch] == 1
    assert acc.regions[2].contained.sum() == 1
    assert acc.regions[0].contained.sum() == 0
    assert acc.regions[1].contained.sum() == 0
    for b in acc.boundaries:
        assert b.flux.sum() == 0
        assert b.mass_left.sum() == 0.0
        assert b.mass_right.sum() == 0.0
    assert acc.total_mass_deposited() == pytest.approx(1.0)


@XFAIL_NATIVE
def test_t2_contained_multi_block_spliced_native():
    raise AssertionError("native impl not landed")


# --- T3: two blocks adjacent regions ----------------------------------------


def test_t3_two_block_adjacent_regions_reference():
    """Blocks in regions 0 and 1 (adjacent). Boundary at edges[1]=100.
    Block1 in R0 = (50,100) len=50; Block2 in R1 = (100,180) len=80.
    L = 130. Spliced+, strand+.
    The boundary between R0 and R1 (boundary index 1) receives:
        mass_left  += 50/130
        mass_right += 80/130
        flux       += 1
    Other boundaries are untouched.
    """
    acc = partition_three_adjacent_exons()
    acc.deposit(
        blocks=[(50, 100), (100, 180)],
        spliced=True,
        strand_pos=True,
    )
    ch = channel_idx(spliced=True, strand_pos=True)
    L = 130.0
    b = acc.boundaries[1]
    assert b.mass_left[ch] == pytest.approx(50.0 / L, abs=1e-6)
    assert b.mass_right[ch] == pytest.approx(80.0 / L, abs=1e-6)
    assert b.flux[ch] == 1
    for i in (0, 2, 3):
        assert acc.boundaries[i].mass_left.sum() == 0.0
        assert acc.boundaries[i].mass_right.sum() == 0.0
        assert acc.boundaries[i].flux.sum() == 0
    for r in acc.regions:
        assert r.contained.sum() == 0
    assert acc.total_mass_deposited() == pytest.approx(1.0)


@XFAIL_NATIVE
def test_t3_two_block_adjacent_regions_native():
    raise AssertionError("native impl not landed")


# --- T4: two blocks non-adjacent regions (user-verified §4.5.1) -------------


def test_t4_two_block_non_adjacent_regions_reference():
    """User-verified walkthrough from 00_design.md §4.5.1.

    Partition: exon1 (1000,2000) = R0, intron (2000,5000) = R1, exon2 (5000,6000) = R2.
    Boundaries: B0@1000, B1@2000, B2@5000, B3@6000.

    Fragment: B1=(1800,1950) in R0 (ℓ=150), B2=(5050,5950) in R2 (ℓ=900).
    L = 1050. Strand +, spliced.

    Expected:
      - B1 (right boundary of R0).mass_left.spl_pos  += 150/1050
      - B1 (right boundary of R0).flux.spl_pos       += 1
      - B2 (left boundary of R2).mass_right.spl_pos += 900/1050
      - B2 (left boundary of R2).flux.spl_pos        += 1
      - B1.mass_right untouched (no block in R1)
      - B2.mass_left  untouched (no block in R1)
      - R1 (intronic gap) untouched
      - Terminal boundaries B0, B3 untouched
    """
    acc = partition_exon_intron_exon()
    acc.deposit(
        blocks=[(1800, 1950), (5050, 5950)],
        spliced=True,
        strand_pos=True,
    )
    ch = channel_idx(spliced=True, strand_pos=True)
    L = 1050.0

    b1 = acc.boundaries[1]  # at position 2000
    b2 = acc.boundaries[2]  # at position 5000
    assert b1.mass_left[ch] == pytest.approx(150.0 / L, abs=1e-6)
    assert b1.mass_right[ch] == 0.0
    assert b1.flux[ch] == 1
    assert b2.mass_right[ch] == pytest.approx(900.0 / L, abs=1e-6)
    assert b2.mass_left[ch] == 0.0
    assert b2.flux[ch] == 1

    assert acc.boundaries[0].mass_left.sum() == 0.0
    assert acc.boundaries[0].mass_right.sum() == 0.0
    assert acc.boundaries[0].flux.sum() == 0
    assert acc.boundaries[3].mass_left.sum() == 0.0
    assert acc.boundaries[3].mass_right.sum() == 0.0
    assert acc.boundaries[3].flux.sum() == 0

    for r in acc.regions:
        assert r.contained.sum() == 0

    assert acc.total_mass_deposited() == pytest.approx(1.0, abs=1e-6)


@XFAIL_NATIVE
def test_t4_two_block_non_adjacent_regions_native():
    raise AssertionError("native impl not landed")


# --- T5: three blocks all adjacent (user-verified §4.5.2) -------------------


def test_t5_three_block_all_adjacent_reference():
    """Three adjacent exonic regions (0,100), (100,200), (200,400).

    Blocks: B1 in R0 = (0,100) ℓ=100; B2 in R1 = (100,180) ℓ=80; B3 in R2 = (200,320) ℓ=120.
    L = 300. Strand +, spliced.

    Junction (B1, B2) at boundary B1 (position 100):
      - B1.mass_left.spl_pos  += 100/300
      - B1.mass_right.spl_pos += 80/300
      - B1.flux.spl_pos       += 1
    Junction (B2, B3) at boundary B2 (position 200):
      - B2.mass_left.spl_pos  += 80/300
      - B2.mass_right.spl_pos += 120/300
      - B2.flux.spl_pos       += 1

    Block-2 contributes mass to two different boundaries — see
    00_design.md §6 for the per-block-side mass conservation note.
    """
    acc = partition_three_adjacent_exons()
    acc.deposit(
        blocks=[(0, 100), (100, 180), (200, 320)],
        spliced=True,
        strand_pos=True,
    )
    ch = channel_idx(spliced=True, strand_pos=True)
    L = 300.0
    b1 = acc.boundaries[1]
    b2 = acc.boundaries[2]
    assert b1.mass_left[ch] == pytest.approx(100.0 / L, abs=1e-6)
    assert b1.mass_right[ch] == pytest.approx(80.0 / L, abs=1e-6)
    assert b1.flux[ch] == 1
    assert b2.mass_left[ch] == pytest.approx(80.0 / L, abs=1e-6)
    assert b2.mass_right[ch] == pytest.approx(120.0 / L, abs=1e-6)
    assert b2.flux[ch] == 1
    assert acc.boundaries[0].mass_left.sum() == 0.0
    assert acc.boundaries[3].mass_right.sum() == 0.0
    for r in acc.regions:
        assert r.contained.sum() == 0
    # Middle block (length 80) contributes to both B1.mass_right and
    # B2.mass_left, so total = 1 + 80/300. See 00_design.md §6.
    assert acc.total_mass_deposited() == pytest.approx(1.0 + 80.0 / L, abs=1e-6)


@XFAIL_NATIVE
def test_t5_three_block_all_adjacent_native():
    raise AssertionError("native impl not landed")


# --- T6: fully spans region -------------------------------------------------


def test_t6_fully_spans_region_reference():
    """Single block straddles an entire region.

    Partition: (0,100), (100,200), (200,300). Single block (50, 250) of length 200.
    Decomposed into per-region slices:
      slice in R0: (50, 100)   ℓ = 50
      slice in R1: (100, 200)  ℓ = 100
      slice in R2: (200, 250)  ℓ = 50
    L = 200. Strand +, unspliced.

    Two implicit junction-events:
    (R0, R1) at boundary B1 (position 100):
      B1.mass_left.unspl_pos  += 50/200 = 0.25
      B1.mass_right.unspl_pos += 100/200 = 0.50
      B1.flux.unspl_pos       += 1
    (R1, R2) at boundary B2 (position 200):
      B2.mass_left.unspl_pos  += 100/200 = 0.50
      B2.mass_right.unspl_pos += 50/200 = 0.25
      B2.flux.unspl_pos       += 1

    Total deposited = 0.25 + 0.50 + 0.50 + 0.25 = 1.5.
    Total > 1.0 because the middle slice contributes to both sides
    — see 00_design.md §6 + audit_phase1.md §4 for the new/legacy
    divergence note.
    """
    edges = np.array([0, 100, 200, 300], dtype=np.int64)
    acc = Accumulator(region_edges=edges)
    acc.deposit(
        blocks=[(50, 250)],
        spliced=False,
        strand_pos=True,
    )
    ch = channel_idx(spliced=False, strand_pos=True)
    L = 200.0
    b1 = acc.boundaries[1]
    b2 = acc.boundaries[2]
    assert b1.mass_left[ch] == pytest.approx(50.0 / L, abs=1e-6)
    assert b1.mass_right[ch] == pytest.approx(100.0 / L, abs=1e-6)
    assert b1.flux[ch] == 1
    assert b2.mass_left[ch] == pytest.approx(100.0 / L, abs=1e-6)
    assert b2.mass_right[ch] == pytest.approx(50.0 / L, abs=1e-6)
    assert b2.flux[ch] == 1
    for r in acc.regions:
        assert r.contained.sum() == 0
    assert acc.total_mass_deposited() == pytest.approx(1.5, abs=1e-6)


@XFAIL_NATIVE
def test_t6_fully_spans_region_native():
    raise AssertionError("native impl not landed")


# --- T7: mass conservation over random contained fragments ------------------


def test_t7_mass_conservation_random_reference():
    """1000 random single-block contained fragments. Each is contained → +1
    per fragment. Total mass deposited must equal exactly 1000.
    """
    edges = np.array([0, 1000, 2000, 3000, 4000], dtype=np.int64)
    acc = Accumulator(region_edges=edges)
    rng = np.random.default_rng(42)
    n = 1000
    for _ in range(n):
        r = rng.integers(0, 4)
        region_start = int(edges[r])
        region_end = int(edges[r + 1])
        length = int(rng.integers(50, 200))
        max_start = region_end - length
        start = int(rng.integers(region_start, max_start))
        end = start + length
        acc.deposit(
            blocks=[(start, end)],
            spliced=bool(rng.integers(0, 2)),
            strand_pos=bool(rng.integers(0, 2)),
        )
    assert acc.total_mass_deposited() == pytest.approx(float(n), abs=1e-3)


@XFAIL_NATIVE
def test_t7_mass_conservation_random_native():
    raise AssertionError("native impl not landed")


# --- T8: flux dedup on adjacent-region two-block fragment -------------------


def test_t8_flux_dedup_adjacent_regions_reference():
    """Two-block adjacent-region fragment hits exactly one boundary.
    That boundary gets flux += 1, not +2.
    """
    acc = partition_three_adjacent_exons()
    acc.deposit(
        blocks=[(50, 100), (100, 180)],
        spliced=False,
        strand_pos=True,
    )
    ch = channel_idx(spliced=False, strand_pos=True)
    assert acc.boundaries[1].flux[ch] == 1
    for i in (0, 2, 3):
        assert acc.boundaries[i].flux.sum() == 0


@XFAIL_NATIVE
def test_t8_flux_dedup_adjacent_regions_native():
    raise AssertionError("native impl not landed")


# --- T10: negative strand attribution ---------------------------------------


def test_t10_negative_strand_attribution_reference():
    """Same as T3 but strand=-. Verifies the strand channel routing."""
    acc = partition_three_adjacent_exons()
    acc.deposit(
        blocks=[(50, 100), (100, 180)],
        spliced=True,
        strand_pos=False,
    )
    ch_neg = channel_idx(spliced=True, strand_pos=False)
    ch_pos = channel_idx(spliced=True, strand_pos=True)
    b = acc.boundaries[1]
    assert b.flux[ch_neg] == 1
    assert b.flux[ch_pos] == 0
    assert b.mass_left[ch_neg] > 0.0
    assert b.mass_left[ch_pos] == 0.0
    assert b.mass_right[ch_neg] > 0.0
    assert b.mass_right[ch_pos] == 0.0


@XFAIL_NATIVE
def test_t10_negative_strand_attribution_native():
    raise AssertionError("native impl not landed")


# --- T11: spliced flag attribution across channels --------------------------


def test_t11_spliced_flag_attribution_reference():
    """One contained spliced fragment + one contained unspliced fragment in
    the same region; assert each goes to the correct channel.
    """
    acc = partition_three_adjacent_exons()
    acc.deposit(blocks=[(110, 190)], spliced=True, strand_pos=True)
    acc.deposit(blocks=[(110, 190)], spliced=False, strand_pos=True)
    ch_spl_pos = channel_idx(spliced=True, strand_pos=True)
    ch_unspl_pos = channel_idx(spliced=False, strand_pos=True)
    assert acc.regions[1].contained[ch_spl_pos] == 1
    assert acc.regions[1].contained[ch_unspl_pos] == 1
    ch_spl_neg = channel_idx(spliced=True, strand_pos=False)
    ch_unspl_neg = channel_idx(spliced=False, strand_pos=False)
    assert acc.regions[1].contained[ch_spl_neg] == 0
    assert acc.regions[1].contained[ch_unspl_neg] == 0


@XFAIL_NATIVE
def test_t11_spliced_flag_attribution_native():
    raise AssertionError("native impl not landed")
