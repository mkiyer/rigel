"""
Fractional accumulator specification tests.

Each test case has two halves:

1. ``test_<case>_reference`` — runs against the pure-Python reference in
   ``_accumulator_reference.py``; locks the spec.
2. ``test_<case>_native`` — runs against the native C++ implementation and
   must match the reference byte-for-byte.

Boundary flux is **per side** (``flux_left`` / ``flux_right``): a contiguous
crossing credits both sides of its one boundary; a spliced intron-skip credits
one side of each flanking boundary (no false exon-intron flux). The deposit
``primary`` bit is the channel-0 selector — genome '+' for unspliced, SENSE
for spliced (the scanner orients spliced reads by the splice motif).

FL histograms are NOT accumulated (audit_phase1.md decision #6).
"""

from __future__ import annotations

import numpy as np
import pytest

from ._accumulator_reference import (
    Accumulator,
    channel_idx,
)

try:
    from rigel.native import Accumulator as NativeAccumulator  # type: ignore

    HAS_NATIVE = True
except (ImportError, AttributeError):
    NativeAccumulator = None  # type: ignore
    HAS_NATIVE = False


XFAIL_NATIVE = pytest.mark.xfail(
    not HAS_NATIVE,
    reason="fractional accumulator native binding unavailable",
    strict=True,
)


# --- Common partition fixtures -----------------------------------------------


def partition_exon_intron_exon():
    """exon1 (1000,2000), intron (2000,5000), exon2 (5000,6000)."""
    return np.array([1000, 2000, 5000, 6000], dtype=np.int64)


def partition_three_adjacent_exons():
    """Three contiguous regions: (0,100), (100,200), (200,400)."""
    return np.array([0, 100, 200, 400], dtype=np.int64)


def _both(make, edges):
    return make(boundary_positions=edges)


# --- T1: contained single block ---------------------------------------------


def _check_t1(acc):
    acc.deposit(blocks=[(1100, 1200)], spliced=False, primary=True)
    ch = channel_idx(spliced=False, primary=True)
    assert acc.regions[0].contained[ch] == 1
    assert acc.regions[0].contained.sum() == 1
    assert acc.regions[1].contained.sum() == 0
    assert acc.regions[2].contained.sum() == 0
    for b in acc.boundaries:
        assert b.mass_left.sum() == 0.0
        assert b.mass_right.sum() == 0.0
        assert b.flux_left.sum() == 0
        assert b.flux_right.sum() == 0
    assert acc.total_mass_deposited() == pytest.approx(1.0)


def test_t1_contained_single_block_reference():
    _check_t1(_both(Accumulator, partition_exon_intron_exon()))


@XFAIL_NATIVE
def test_t1_contained_single_block_native():
    _check_t1(_both(NativeAccumulator, partition_exon_intron_exon()))


# --- T2: contained, multi-block spliced -------------------------------------


def _check_t2(acc):
    acc.deposit(blocks=[(5100, 5200), (5400, 5500)], spliced=True, primary=True)
    ch = channel_idx(spliced=True, primary=True)
    assert acc.regions[2].contained[ch] == 1
    assert acc.regions[2].contained.sum() == 1
    assert acc.regions[0].contained.sum() == 0
    assert acc.regions[1].contained.sum() == 0
    for b in acc.boundaries:
        assert b.flux_left.sum() == 0
        assert b.flux_right.sum() == 0
        assert b.mass_left.sum() == 0.0
        assert b.mass_right.sum() == 0.0
    assert acc.total_mass_deposited() == pytest.approx(1.0)


def test_t2_contained_multi_block_spliced_reference():
    _check_t2(_both(Accumulator, partition_exon_intron_exon()))


@XFAIL_NATIVE
def test_t2_contained_multi_block_spliced_native():
    _check_t2(_both(NativeAccumulator, partition_exon_intron_exon()))


# --- T3: two blocks adjacent regions (contiguous crossing) ------------------


def _check_t3(acc):
    # Block1 R0 = (50,100) len 50; Block2 R1 = (100,180) len 80. L=130.
    # Contiguous crossing of boundary 1 → BOTH sides credited.
    acc.deposit(blocks=[(50, 100), (100, 180)], spliced=True, primary=True)
    ch = channel_idx(spliced=True, primary=True)
    L = 130.0
    b = acc.boundaries[1]
    assert b.mass_left[ch] == pytest.approx(50.0 / L, abs=1e-6)
    assert b.mass_right[ch] == pytest.approx(80.0 / L, abs=1e-6)
    assert b.flux_left[ch] == 1
    assert b.flux_right[ch] == 1
    for i in (0, 2, 3):
        assert acc.boundaries[i].mass_left.sum() == 0.0
        assert acc.boundaries[i].mass_right.sum() == 0.0
        assert acc.boundaries[i].flux_left.sum() == 0
        assert acc.boundaries[i].flux_right.sum() == 0
    for r in acc.regions:
        assert r.contained.sum() == 0
    assert acc.total_mass_deposited() == pytest.approx(1.0)


def test_t3_two_block_adjacent_regions_reference():
    _check_t3(_both(Accumulator, partition_three_adjacent_exons()))


@XFAIL_NATIVE
def test_t3_two_block_adjacent_regions_native():
    _check_t3(_both(NativeAccumulator, partition_three_adjacent_exons()))


# --- T4: two blocks non-adjacent (intron skip) — PER-SIDE flux --------------


def _check_t4(acc):
    # B1=(1800,1950) in R0 (150); B2=(5050,5950) in R2 (900). R1 (intron)
    # skipped. L=1050. The exon1→intron boundary (B1) is credited only on its
    # LEFT side; the intron→exon2 boundary (B2) only on its RIGHT side — the
    # intron-facing sides stay zero (no false flux).
    acc.deposit(blocks=[(1800, 1950), (5050, 5950)], spliced=True, primary=True)
    ch = channel_idx(spliced=True, primary=True)
    L = 1050.0
    b1 = acc.boundaries[1]  # position 2000 (exon1 → intron)
    b2 = acc.boundaries[2]  # position 5000 (intron → exon2)
    assert b1.mass_left[ch] == pytest.approx(150.0 / L, abs=1e-6)
    assert b1.mass_right[ch] == 0.0
    assert b1.flux_left[ch] == 1
    assert b1.flux_right[ch] == 0  # intron-facing side: NO false flux
    assert b2.mass_right[ch] == pytest.approx(900.0 / L, abs=1e-6)
    assert b2.mass_left[ch] == 0.0
    assert b2.flux_right[ch] == 1
    assert b2.flux_left[ch] == 0  # intron-facing side: NO false flux

    for i in (0, 3):
        assert acc.boundaries[i].mass_left.sum() == 0.0
        assert acc.boundaries[i].mass_right.sum() == 0.0
        assert acc.boundaries[i].flux_left.sum() == 0
        assert acc.boundaries[i].flux_right.sum() == 0
    for r in acc.regions:
        assert r.contained.sum() == 0
    assert acc.total_mass_deposited() == pytest.approx(1.0, abs=1e-6)


def test_t4_two_block_non_adjacent_regions_reference():
    _check_t4(_both(Accumulator, partition_exon_intron_exon()))


@XFAIL_NATIVE
def test_t4_two_block_non_adjacent_regions_native():
    _check_t4(_both(NativeAccumulator, partition_exon_intron_exon()))


# --- T5: three blocks all adjacent (contiguous) -----------------------------


def _check_t5(acc):
    # B1 R0=(0,100); B2 R1=(100,180); B3 R2=(200,320). L=300. R1 is the
    # ENCOMPASSED interior region (crossed at both boundaries) → its slice mass
    # (80/L) splits 50/50: 40/L to B1.mass_right and 40/L to B2.mass_left. The
    # end slices (R0, R2) keep full mass. Fragment mass conserves to 1.0.
    acc.deposit(blocks=[(0, 100), (100, 180), (200, 320)], spliced=True, primary=True)
    ch = channel_idx(spliced=True, primary=True)
    L = 300.0
    b1 = acc.boundaries[1]
    b2 = acc.boundaries[2]
    assert b1.mass_left[ch] == pytest.approx(100.0 / L, abs=1e-6)
    assert b1.mass_right[ch] == pytest.approx(40.0 / L, abs=1e-6)
    assert b1.flux_left[ch] == 1
    assert b1.flux_right[ch] == 1
    assert b2.mass_left[ch] == pytest.approx(40.0 / L, abs=1e-6)
    assert b2.mass_right[ch] == pytest.approx(120.0 / L, abs=1e-6)
    assert b2.flux_left[ch] == 1
    assert b2.flux_right[ch] == 1
    assert acc.boundaries[0].mass_left.sum() == 0.0
    assert acc.boundaries[3].mass_right.sum() == 0.0
    for r in acc.regions:
        assert r.contained.sum() == 0
    assert acc.total_mass_deposited() == pytest.approx(1.0, abs=1e-6)


def test_t5_three_block_all_adjacent_reference():
    _check_t5(_both(Accumulator, partition_three_adjacent_exons()))


@XFAIL_NATIVE
def test_t5_three_block_all_adjacent_native():
    _check_t5(_both(NativeAccumulator, partition_three_adjacent_exons()))


# --- T6: fully spans region (single block straddles regions, contiguous) ----


def _check_t6(acc):
    # Single block (50,250) over partition [0,100,200,300]: R0=(0,100) gets
    # 50bp, R1=(100,200) is fully ENCOMPASSED (100bp slice), R2=(200,300) gets
    # 50bp. L=200. R1's slice mass (100/L) splits 50/50 → 50/L to B1.mass_right
    # and 50/L to B2.mass_left; the end slices keep full mass. Fragment mass
    # conserves to 1.0 (was the 1.5 double-count bug, PR 5.5).
    acc.deposit(blocks=[(50, 250)], spliced=False, primary=True)
    ch = channel_idx(spliced=False, primary=True)
    L = 200.0
    b1 = acc.boundaries[1]
    b2 = acc.boundaries[2]
    assert b1.mass_left[ch] == pytest.approx(50.0 / L, abs=1e-6)
    assert b1.mass_right[ch] == pytest.approx(50.0 / L, abs=1e-6)
    assert b1.flux_left[ch] == 1
    assert b1.flux_right[ch] == 1
    assert b2.mass_left[ch] == pytest.approx(50.0 / L, abs=1e-6)
    assert b2.mass_right[ch] == pytest.approx(50.0 / L, abs=1e-6)
    assert b2.flux_left[ch] == 1
    assert b2.flux_right[ch] == 1
    for r in acc.regions:
        assert r.contained.sum() == 0
    assert acc.total_mass_deposited() == pytest.approx(1.0, abs=1e-6)


def test_t6_fully_spans_region_reference():
    _check_t6(_both(Accumulator, np.array([0, 100, 200, 300], dtype=np.int64)))


@XFAIL_NATIVE
def test_t6_fully_spans_region_native():
    _check_t6(_both(NativeAccumulator, np.array([0, 100, 200, 300], dtype=np.int64)))


# --- T6b: mass conservation for region-spanning fragments (PR 5.5) ----------


def _check_t6b(acc):
    # Every fragment — contained, simple crossing, or spanning whole regions —
    # deposits total mass exactly 1.0. Partition [0,100,..,500] → 5×100bp
    # regions; fragments span 1..5 regions (with 0..3 encompassed interiors).
    frags = [
        [(10, 50)],  # contained in R0 (M=1)
        [(50, 150)],  # crosses R0→R1 (M=2, no interior)
        [(50, 250)],  # spans R0, R1(enc), R2 (M=3, 1 interior)
        [(50, 450)],  # spans R0..R4 (M=5, 3 interiors)
        [(150, 350)],  # spans R1, R2(enc), R3 (M=3, 1 interior)
    ]
    for blocks in frags:
        acc.deposit(blocks=blocks, spliced=False, primary=True)
    assert acc.total_mass_deposited() == pytest.approx(float(len(frags)), abs=1e-5)


def test_t6b_spanning_mass_conservation_reference():
    edges = np.array([0, 100, 200, 300, 400, 500], dtype=np.int64)
    _check_t6b(_both(Accumulator, edges))


@XFAIL_NATIVE
def test_t6b_spanning_mass_conservation_native():
    edges = np.array([0, 100, 200, 300, 400, 500], dtype=np.int64)
    _check_t6b(_both(NativeAccumulator, edges))


# --- T7: mass conservation over random contained fragments ------------------


def _check_t7(acc, edges):
    rng = np.random.default_rng(42)
    n = 1000
    for _ in range(n):
        r = rng.integers(0, 4)
        region_start = int(edges[r])
        region_end = int(edges[r + 1])
        length = int(rng.integers(50, 200))
        start = int(rng.integers(region_start, region_end - length))
        acc.deposit(
            blocks=[(start, start + length)],
            spliced=bool(rng.integers(0, 2)),
            primary=bool(rng.integers(0, 2)),
        )
    assert acc.total_mass_deposited() == pytest.approx(float(n), abs=1e-3)


def test_t7_mass_conservation_random_reference():
    edges = np.array([0, 1000, 2000, 3000, 4000], dtype=np.int64)
    _check_t7(_both(Accumulator, edges), edges)


@XFAIL_NATIVE
def test_t7_mass_conservation_random_native():
    edges = np.array([0, 1000, 2000, 3000, 4000], dtype=np.int64)
    _check_t7(_both(NativeAccumulator, edges), edges)


# --- T8: per-side flux on contiguous crossing (each side at most once) ------


def _check_t8(acc):
    acc.deposit(blocks=[(50, 100), (100, 180)], spliced=False, primary=True)
    ch = channel_idx(spliced=False, primary=True)
    # Contiguous crossing of boundary 1 → each side credited exactly once.
    assert acc.boundaries[1].flux_left[ch] == 1
    assert acc.boundaries[1].flux_right[ch] == 1
    for i in (0, 2, 3):
        assert acc.boundaries[i].flux_left.sum() == 0
        assert acc.boundaries[i].flux_right.sum() == 0


def test_t8_per_side_flux_adjacent_regions_reference():
    _check_t8(_both(Accumulator, partition_three_adjacent_exons()))


@XFAIL_NATIVE
def test_t8_per_side_flux_adjacent_regions_native():
    _check_t8(_both(NativeAccumulator, partition_three_adjacent_exons()))


# --- T10: negative-/secondary-channel attribution ---------------------------


def _check_t10(acc):
    acc.deposit(blocks=[(50, 100), (100, 180)], spliced=True, primary=False)
    ch_sec = channel_idx(spliced=True, primary=False)
    ch_pri = channel_idx(spliced=True, primary=True)
    b = acc.boundaries[1]
    assert b.flux_left[ch_sec] == 1
    assert b.flux_right[ch_sec] == 1
    assert b.flux_left[ch_pri] == 0
    assert b.flux_right[ch_pri] == 0
    assert b.mass_left[ch_sec] > 0.0
    assert b.mass_left[ch_pri] == 0.0
    assert b.mass_right[ch_sec] > 0.0
    assert b.mass_right[ch_pri] == 0.0


def test_t10_secondary_channel_attribution_reference():
    _check_t10(_both(Accumulator, partition_three_adjacent_exons()))


@XFAIL_NATIVE
def test_t10_secondary_channel_attribution_native():
    _check_t10(_both(NativeAccumulator, partition_three_adjacent_exons()))


# --- T11: spliced flag attribution across channels --------------------------


def _check_t11(acc):
    acc.deposit(blocks=[(110, 190)], spliced=True, primary=True)
    acc.deposit(blocks=[(110, 190)], spliced=False, primary=True)
    assert acc.regions[1].contained[channel_idx(spliced=True, primary=True)] == 1
    assert acc.regions[1].contained[channel_idx(spliced=False, primary=True)] == 1
    assert acc.regions[1].contained[channel_idx(spliced=True, primary=False)] == 0
    assert acc.regions[1].contained[channel_idx(spliced=False, primary=False)] == 0


def test_t11_spliced_flag_attribution_reference():
    _check_t11(_both(Accumulator, partition_three_adjacent_exons()))


@XFAIL_NATIVE
def test_t11_spliced_flag_attribution_native():
    _check_t11(_both(NativeAccumulator, partition_three_adjacent_exons()))


# --- T12: native matches reference byte-for-byte on an intron-skip ----------


@pytest.mark.skipif(not HAS_NATIVE, reason="native binding unavailable")
def test_native_matches_reference_intron_skip():
    edges = partition_exon_intron_exon()
    ref = Accumulator(boundary_positions=edges)
    nat = NativeAccumulator(boundary_positions=edges)
    for acc in (ref, nat):
        acc.deposit(blocks=[(1800, 1950), (5050, 5950)], spliced=True, primary=True)
        acc.deposit(blocks=[(1100, 1300)], spliced=False, primary=False)
    for i in range(len(ref.boundaries)):
        np.testing.assert_allclose(ref.boundaries[i].mass_left, nat.boundaries[i].mass_left)
        np.testing.assert_allclose(ref.boundaries[i].mass_right, nat.boundaries[i].mass_right)
        np.testing.assert_array_equal(ref.boundaries[i].flux_left, nat.boundaries[i].flux_left)
        np.testing.assert_array_equal(ref.boundaries[i].flux_right, nat.boundaries[i].flux_right)
    for i in range(len(ref.regions)):
        np.testing.assert_array_equal(ref.regions[i].contained, nat.regions[i].contained)


# --- FL pools (PR 4c): gDNA fragment-length histograms -----------------------
#
# Partition [0,100,200,400] → regions (0,100), (100,200), (200,400) with
# region_types [exon=2, intron=1, intergenic=0]. FL pool index = type*2 +
# (1 if boundary else 0); gDNA pool = intergenic+intronic, both compartments.


def _fl_deposits(acc):
    # contained in region 2 (intergenic): footprint 50 → INTERGENIC_CONTAINED.
    acc.deposit(blocks=[(210, 260)], spliced=False, primary=True)
    # crossing region 1→2 (intron→intergenic): footprint 100, 50/50 each side.
    acc.deposit(blocks=[(150, 250)], spliced=False, primary=True)
    # spliced fragment → excluded from FL pools entirely.
    acc.deposit(blocks=[(10, 50), (110, 150)], spliced=True, primary=True)


def _check_fl(acc):
    fl = np.asarray(acc.fl_pool_mass)
    assert fl.shape == (6, 1001)
    np.testing.assert_allclose(fl[0, 50], 1.0)  # INTERGENIC_CONTAINED, bin 50
    np.testing.assert_allclose(fl[1, 100], 0.5)  # INTERGENIC_BOUNDARY, bin 100
    np.testing.assert_allclose(fl[3, 100], 0.5)  # INTRONIC_BOUNDARY, bin 100
    np.testing.assert_allclose(fl.sum(), 2.0)  # 2 unspliced frags; spliced excluded


def test_fl_pools_reference():
    edges = partition_three_adjacent_exons()
    types = np.array([2, 1, 0], dtype=np.uint8)
    acc = Accumulator(boundary_positions=edges, region_types=types, max_fl=1000)
    _fl_deposits(acc)
    _check_fl(acc)


@XFAIL_NATIVE
def test_fl_pools_native_matches_reference():
    edges = partition_three_adjacent_exons()
    types = np.array([2, 1, 0], dtype=np.uint8)
    ref = Accumulator(boundary_positions=edges, region_types=types, max_fl=1000)
    nat = NativeAccumulator(boundary_positions=edges, region_types=types, max_fl=1000)
    _fl_deposits(ref)
    _fl_deposits(nat)
    np.testing.assert_allclose(np.asarray(nat.fl_pool_mass), ref.fl_pool_mass)
    _check_fl(nat)


@XFAIL_NATIVE
def test_fl_pools_disabled_without_region_types():
    edges = partition_three_adjacent_exons()
    nat = NativeAccumulator(boundary_positions=edges)  # no region_types → FL off
    nat.deposit(blocks=[(210, 260)], spliced=False, primary=True)
    assert np.asarray(nat.fl_pool_mass).shape == (6, 0)


# --- FL footprint = SPAN, not covered length (regression) --------------------
#
# A fragment whose blocks have an inter-block GAP (paired mates with insert >
# read1+read2 covered bases) must bin its FL at the genomic SPAN
# (max end − min start), NOT the covered length (Σ block lengths). Both blocks
# (210,260) and (310,360) lie in region 2 (intergenic [200,400)): covered =
# 50+50 = 100, span = 360−210 = 150. Pre-fix the pool binned at the covered
# length, collapsing long gDNA fragments to a spike at ~2×read_length and
# leaking typical-length gDNA to RNA in the scorer (which queries the pmf at the
# genomic footprint = span). The gap case was previously untested.


def _fl_gap_deposit(acc):
    acc.deposit(blocks=[(210, 260), (310, 360)], spliced=False, primary=True)


def _check_fl_gap(acc):
    fl = np.asarray(acc.fl_pool_mass)
    np.testing.assert_allclose(fl[0, 150], 1.0)  # INTERGENIC_CONTAINED at SPAN 150
    assert fl[0, 100] == 0.0  # NOT the covered length 100 (the pre-fix bug bin)
    np.testing.assert_allclose(fl.sum(), 1.0)


def test_fl_footprint_is_span_reference():
    edges = partition_three_adjacent_exons()
    types = np.array([2, 1, 0], dtype=np.uint8)
    acc = Accumulator(boundary_positions=edges, region_types=types, max_fl=1000)
    _fl_gap_deposit(acc)
    _check_fl_gap(acc)


@XFAIL_NATIVE
def test_fl_footprint_is_span_native():
    edges = partition_three_adjacent_exons()
    types = np.array([2, 1, 0], dtype=np.uint8)
    ref = Accumulator(boundary_positions=edges, region_types=types, max_fl=1000)
    nat = NativeAccumulator(boundary_positions=edges, region_types=types, max_fl=1000)
    _fl_gap_deposit(ref)
    _fl_gap_deposit(nat)
    np.testing.assert_allclose(np.asarray(nat.fl_pool_mass), ref.fl_pool_mass)
    _check_fl_gap(nat)


# --- Junction strand: per-boundary genomic strand of the splice junction --------
# A spliced crossing records its motif (genomic) strand on every boundary its
# spliced mass touches; unspliced/contained fragments leave it 0. Partition:
# exon1 (1000,2000) | intron (2000,5000) | exon2 (5000,6000) → boundaries at
# indices 0..3 (positions 1000,2000,5000,6000). A spliced crossing exon1→exon2
# deposits at the donor boundary (idx 1, pos 2000) and acceptor (idx 2, pos 5000).


def _js(acc) -> list[int]:
    return [int(acc.boundaries[i].junction_strand) for i in range(len(acc.boundaries))]


def _check_junction_strand_pos(acc):
    acc.deposit(blocks=[(1800, 1950), (5050, 5950)], spliced=True, primary=True, strand=1)
    assert _js(acc) == [0, 1, 1, 0]  # POS on the two junction boundaries only


def _check_junction_strand_neg(acc):
    acc.deposit(blocks=[(1800, 1950), (5050, 5950)], spliced=True, primary=False, strand=2)
    assert _js(acc) == [0, 2, 2, 0]  # NEG; the SENSE/ANTISENSE (primary) channel is independent


def _check_junction_strand_unspliced_noop(acc):
    # an unspliced contiguous crossing of boundary idx 1 records mass but NO junction strand
    acc.deposit(blocks=[(1900, 2100)], spliced=False, primary=True, strand=1)
    assert _js(acc) == [0, 0, 0, 0]
    assert acc.boundaries[1].mass_left.sum() + acc.boundaries[1].mass_right.sum() > 0.0


def test_junction_strand_pos_reference():
    _check_junction_strand_pos(_both(Accumulator, partition_exon_intron_exon()))


@XFAIL_NATIVE
def test_junction_strand_pos_native():
    _check_junction_strand_pos(_both(NativeAccumulator, partition_exon_intron_exon()))


def test_junction_strand_neg_reference():
    _check_junction_strand_neg(_both(Accumulator, partition_exon_intron_exon()))


@XFAIL_NATIVE
def test_junction_strand_neg_native():
    _check_junction_strand_neg(_both(NativeAccumulator, partition_exon_intron_exon()))


def test_junction_strand_unspliced_noop_reference():
    _check_junction_strand_unspliced_noop(_both(Accumulator, partition_exon_intron_exon()))


@XFAIL_NATIVE
def test_junction_strand_unspliced_noop_native():
    _check_junction_strand_unspliced_noop(_both(NativeAccumulator, partition_exon_intron_exon()))


# ---------------------------------------------------------------------------
# Worker-merge DETERMINISM — the control for accumulator v5's test A9.
# ---------------------------------------------------------------------------
#
# v5 §6 claims that moving every channel to integer / fixed-point storage makes the accumulator
# "bit-exact at any thread count", because integer addition is associative — and that this REMOVES
# the known ~2.6 % cross-process nondeterminism rather than merely working around it (memory
# `calibrate_cross_process_nondeterminism`). Test A9 is that claim.
#
# ⚠ Before this, `Accumulator::merge_from` was defined in C++ (accumulator.cpp:234) and called only
# from `BamScanner::scan` — it was NOT bound to Python, so A9 had no infrastructure and the v5 claim
# had no *control*: nothing recorded what today's accumulator actually does. These tests are that
# control. They must keep passing through the v5 rewrite, and the float assertion below is expected
# to STRENGTHEN from "close" to "exact" when the fractional mass channels disappear.


class TestWorkerMergeDeterminism:
    """Deposit one corpus across K workers in K different orders; merge; compare to a single pass."""

    @staticmethod
    def _corpus(rng, n, hi):
        out = []
        for _ in range(n):
            s = int(rng.integers(0, hi - 20))
            w = int(rng.integers(20, 260))
            spliced = bool(rng.integers(0, 2))
            blocks = [(s, min(s + w, hi))]
            if spliced and s + w + 300 < hi:  # a spliced fragment: two blocks over an intron
                blocks = [(s, s + w), (s + w + 200, min(s + w + 200 + w, hi))]
            out.append((blocks, spliced, bool(rng.integers(0, 2)), int(rng.integers(0, 3))))
        return out

    @staticmethod
    def _cuts():
        return np.array([0, 500, 550, 560, 561, 1561, 2000], dtype=np.int64)

    def _run(self, frags, n_workers, rng=None):
        """Shard `frags` across n_workers (in shuffled order when rng is given), merge, return arrays."""
        cuts = self._cuts()
        order = list(range(len(frags)))
        if rng is not None:
            rng.shuffle(order)
        accs = [NativeAccumulator(boundary_positions=cuts) for _ in range(n_workers)]
        for i, k in enumerate(order):
            blocks, spliced, primary, strand = frags[k]
            accs[i % n_workers].deposit(blocks, spliced=spliced, primary=primary, strand=strand)
        merged = accs[0]
        for a in accs[1:]:
            merged.merge_from(a)
        n = merged._native
        return (
            np.asarray(n.regions_contained).copy(),
            np.asarray(n.boundaries_flux_left).copy(),
            np.asarray(n.boundaries_flux_right).copy(),
            np.asarray(n.boundaries_mass_left).copy(),
            np.asarray(n.boundaries_mass_right).copy(),
        )

    def test_merge_equals_single_pass(self):
        """K workers + merge == one accumulator seeing every fragment. The merge itself is correct."""
        frags = self._corpus(np.random.default_rng(20260729), 800, 2000)
        one = self._run(frags, 1)
        many = self._run(frags, 4)
        np.testing.assert_array_equal(one[0], many[0])  # contained counts — INTEGER, exact
        np.testing.assert_array_equal(one[1], many[1])  # flux left    — INTEGER, exact
        np.testing.assert_array_equal(one[2], many[2])  # flux right   — INTEGER, exact
        np.testing.assert_allclose(one[3], many[3], rtol=1e-6)  # mass — float32, see below
        np.testing.assert_allclose(one[4], many[4], rtol=1e-6)

    @pytest.mark.parametrize("n_workers", [1, 2, 4, 8])
    def test_integer_channels_are_bit_identical_at_any_worker_count(self, n_workers):
        """⭐ The INTEGER channels are already exactly worker-count- and order-independent.

        Integer addition is associative and commutative, so no shard split and no merge order can move
        them. This is the property v5 §6 extends to EVERY channel by replacing the float masses with
        fixed-point `uint64` sums — at which point this test's float sibling below becomes exact too.
        """
        frags = self._corpus(np.random.default_rng(7), 800, 2000)
        ref = self._run(frags, 1)
        got = self._run(frags, n_workers, rng=np.random.default_rng(100 + n_workers))
        np.testing.assert_array_equal(ref[0], got[0])
        np.testing.assert_array_equal(ref[1], got[1])
        np.testing.assert_array_equal(ref[2], got[2])

    def test_float_mass_channels_are_NOT_bit_identical_across_shardings(self):
        """⚠ The CONTROL, and the reason v5 §6 is worth doing: the float channels are order-dependent.

        `mass_left`/`mass_right` accumulate `float32 += share` over a data-dependent worker partition,
        and float addition is not associative — so a different sharding gives a different sum in the
        low bits. That is the documented ~2.6 % cross-process nondeterminism at its source.

        This test asserts the CURRENT behaviour: close, but not bit-identical. **When v5's fixed-point
        `uint64` recip lands, this test must be inverted to assert exact equality** — that inversion is
        the v5 §6 / test-A9 deliverable, and if it cannot be inverted the claim is false.
        """
        frags = self._corpus(np.random.default_rng(11), 4000, 2000)
        ref = self._run(frags, 1)
        got = self._run(frags, 8, rng=np.random.default_rng(999))
        np.testing.assert_allclose(ref[3], got[3], rtol=1e-5)
        np.testing.assert_allclose(ref[4], got[4], rtol=1e-5)
        exact = np.array_equal(ref[3], got[3]) and np.array_equal(ref[4], got[4])
        if exact:  # not a failure — it means the corpus was too small to expose the reassociation
            pytest.skip(
                "float masses happened to be bit-identical; corpus too small to reassociate"
            )
