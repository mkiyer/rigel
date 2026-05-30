"""
Pure-Python reference implementation of fractional accumulator.

This module is the canonical interpreter of the fractional accumulator spec
(see docs/accumulator/00_design.md). The native C++ implementation
landing in Phase 2 must reproduce this module's behavior on every
synthetic fragment in tests/native/test_accumulator_spec.py.

Decisions locked in docs/accumulator/audit_phase1.md:
- Per-fragment $L = \\sum_k \\ell_k$; each block-region overlap deposits
  overlap_bp / L.
- Single per-fragment splice flag (not per-junction).
- Single per-fragment strand; ambiguous fragments are dropped upstream
  and never reach `deposit`.
- Region storage: uint32 contained counts only. Contained fragment
  increments by exactly +1 regardless of block count.
- Boundary storage: float32 mass_left/mass_right per channel + uint32
  flux_left/flux_right per channel (per-side). A contiguous crossing
  credits both sides of its boundary; a spliced intron-skip credits one
  side of each flanking boundary.
- **FL histograms are not accumulated.** FL is used downstream by the
  EM scorer (see `rigel.scoring.FragmentScorer` and
  `rigel.frag_length_model.FragmentLengthModel`), not by calibration.
  Spliced fragments are FL-ambiguous (compatible with many isoforms);
  per-region FL accumulation is deferred indefinitely.
- Fully-spans-region: decompose straddling block into per-region slices
  and apply the §4.3 adjacent-pair rule to consecutive slices.

This file is test-only — do not import from src/rigel/.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# Channel encoding: (spliced, primary) -> channel index in 0..3
# 0: unspl_primary, 1: unspl_!primary, 2: spl_primary, 3: spl_!primary
# where primary = genome '+' (unspliced) or SENSE (spliced, align==motif).
def channel_idx(spliced: bool, primary: bool) -> int:
    return (2 if spliced else 0) + (0 if primary else 1)


@dataclass
class Region:
    """One region on a reference. Channels are uint32 counts."""

    contained: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=np.uint32))


@dataclass
class Boundary:
    """One boundary on a reference.

    mass_left[c] = sum of (block_len / L) for blocks lying in the
                   region to the LEFT of this boundary, channel c.
    mass_right[c] = same, for blocks to the RIGHT.
    flux_left[c] / flux_right[c] = per-side count of fragment-events
              touching the LEFT / RIGHT side of this boundary, channel c.
              The left region's slice credits flux_left; the right region's
              slice credits flux_right. Slices are monotonic so each side is
              credited at most once per fragment.
    """

    mass_left: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=np.float32))
    mass_right: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=np.float32))
    flux_left: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=np.uint32))
    flux_right: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=np.uint32))


@dataclass
class Accumulator:
    """Per-reference partition: N regions, N+1 boundaries.

    Region i is bordered by boundaries i (left) and i+1 (right).
    boundaries[0] and boundaries[N] are terminal.

    `boundary_positions` is a sorted int64 array of length N+1 giving the
    boundary positions in genomic coordinates:
        regions[i]    = [boundary_positions[i], boundary_positions[i+1])
        boundaries[i] is positioned at boundary_positions[i]
    Regions are contiguous on the reference; intronic gaps are
    represented as their own regions.
    """

    boundary_positions: np.ndarray  # int64[N+1], strictly increasing
    regions: list[Region] = field(default_factory=list)
    boundaries: list[Boundary] = field(default_factory=list)

    def __post_init__(self) -> None:
        n = len(self.boundary_positions) - 1
        if n < 0:
            raise ValueError("boundary_positions must have length >= 1")
        if not self.regions:
            self.regions = [Region() for _ in range(n)]
        if not self.boundaries:
            self.boundaries = [Boundary() for _ in range(n + 1)]

    # ---- Geometry helpers -------------------------------------------------

    @property
    def n_regions(self) -> int:
        return len(self.regions)

    def region_of_pos(self, pos: int) -> int:
        """Return the region index containing `pos`, or -1 if outside."""
        b = self.boundary_positions
        if pos < b[0] or pos >= b[-1]:
            return -1
        idx = int(np.searchsorted(b, pos, side="right")) - 1
        return idx

    def left_boundary_of(self, region_idx: int) -> int:
        return region_idx

    def right_boundary_of(self, region_idx: int) -> int:
        return region_idx + 1

    # ---- Core deposit -----------------------------------------------------

    def deposit(
        self,
        blocks: list[tuple[int, int]],
        spliced: bool,
        primary: bool,
    ) -> None:
        """Deposit one fragment's evidence.

        blocks: list of (start, end) aligned blocks in genomic order.
                Sum of (end-start) is L. Soft-clipped bp must already
                be excluded.
        spliced: single per-fragment splice classification.
        primary: channel-0 selector — genome '+' for unspliced, SENSE
                 (align_strand == sj_strand) for spliced.
        """
        if not blocks:
            return
        ch = channel_idx(spliced, primary)

        # 1. Expand each block into per-region slices. A block may
        #    straddle one or more region boundaries (the "fully spans
        #    region" case); the resulting slices behave like adjacent
        #    blocks with no inter-block gap.
        slices: list[tuple[int, int, int]] = []  # (region_idx, start, end)
        for blk_start, blk_end in blocks:
            if blk_end <= blk_start:
                continue
            s = max(int(blk_start), int(self.boundary_positions[0]))
            e = min(int(blk_end), int(self.boundary_positions[-1]))
            if e <= s:
                continue
            cur = s
            r = self.region_of_pos(cur)
            while cur < e and r != -1 and r < self.n_regions:
                region_end = int(self.boundary_positions[r + 1])
                slice_end = min(e, region_end)
                slices.append((r, cur, slice_end))
                cur = slice_end
                r += 1

        if not slices:
            return

        L = sum(end - start for _, start, end in slices)
        if L <= 0:
            return
        inv_L = 1.0 / float(L)

        # 3. Single-region (all slices in same region) → contained.
        regions_touched = {sl[0] for sl in slices}
        if len(regions_touched) == 1:
            r = slices[0][0]
            self.regions[r].contained[ch] += np.uint32(1)
            return

        # 4. Crossing path. For each consecutive slice pair, deposit on the
        #    boundary between them per §4.3 of 00_design.md, with PER-SIDE
        #    flux: the left region's slice credits the LEFT side of b_out
        #    (matching its mass_left); the right region's slice credits the
        #    RIGHT side of b_in (matching its mass_right). Contiguous crossing
        #    (b_out == b_in) credits both sides; a spliced jump credits one
        #    side of each flanking boundary. Slices are monotonic → each side
        #    is credited at most once per fragment (no dedup needed).
        for i in range(len(slices) - 1):
            r_left, _, _ = slices[i]
            r_right, _, _ = slices[i + 1]
            len_left = slices[i][2] - slices[i][1]
            len_right = slices[i + 1][2] - slices[i + 1][1]

            b_out = self.right_boundary_of(r_left)
            b_in = self.left_boundary_of(r_right)

            self.boundaries[b_out].mass_left[ch] += np.float32(len_left * inv_L)
            self.boundaries[b_out].flux_left[ch] += np.uint32(1)
            self.boundaries[b_in].mass_right[ch] += np.float32(len_right * inv_L)
            self.boundaries[b_in].flux_right[ch] += np.uint32(1)

    # ---- Invariants -------------------------------------------------------

    def total_mass_deposited(self) -> float:
        """Sum of all mass + counts. For per-block-side mass conservation
        details see docs/accumulator/00_design.md §6.
        """
        total = 0.0
        for r in self.regions:
            total += float(r.contained.sum())
        for b in self.boundaries:
            total += float(b.mass_left.sum()) + float(b.mass_right.sum())
        return total
