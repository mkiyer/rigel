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


# gDNA FL pools (PR 4c): 3 region-types {0=intergenic, 1=intronic, 2=exonic}
# × 2 compartments {0=contained, 1=boundary}. Matches the C++ fl_pool_idx.
N_FL_POOLS = 6


def fl_pool_idx(region_type: int, boundary: bool) -> int:
    return int(region_type) * 2 + (1 if boundary else 0)


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
    #: Splice-junction genomic strand (Strand: POS=1/NEG=2, 0=no junction). Set
    #: from the motif strand on a spliced crossing — the mature-RNA anchor's strand.
    junction_strand: int = 0


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
    # gDNA FL pools (PR 4c): enabled iff region_types (len N) + max_fl > 0.
    region_types: "np.ndarray | None" = None
    max_fl: int = 0
    fl_pool_mass: "np.ndarray | None" = None

    def __post_init__(self) -> None:
        n = len(self.boundary_positions) - 1
        if n < 0:
            raise ValueError("boundary_positions must have length >= 1")
        if not self.regions:
            self.regions = [Region() for _ in range(n)]
        if not self.boundaries:
            self.boundaries = [Boundary() for _ in range(n + 1)]
        if self.region_types is not None and self.max_fl > 0 and len(self.region_types) == n:
            self.fl_pool_mass = np.zeros((N_FL_POOLS, self.max_fl + 1), dtype=np.float64)
        else:
            self.region_types = None
            self.max_fl = 0
            self.fl_pool_mass = None

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
        strand: int = 0,
    ) -> None:
        """Deposit one fragment's evidence.

        blocks: list of (start, end) aligned blocks in genomic order.
                Sum of (end-start) is L. Soft-clipped bp must already
                be excluded.
        spliced: single per-fragment splice classification.
        primary: channel-0 selector — genome '+' for unspliced, SENSE
                 (align_strand == sj_strand) for spliced.
        strand: genomic strand (Strand: POS=1/NEG=2, 0=none). For a SPLICED
                crossing it is the junction motif strand, recorded on every
                boundary the spliced mass touches; ignored otherwise.
        """
        if not blocks:
            return
        ch = channel_idx(spliced, primary)
        js = int(strand) if spliced else 0  # junction strand only for spliced

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

        # gDNA FL pooling (PR 4c): bin at the fragment's genomic SPAN (template
        # footprint = max block end − min block start), NOT the covered length L
        # (= Σ slice lengths). For paired mates with an inter-mate gap the covered
        # length saturates at the read-length sum, collapsing the gDNA FL pmf to a
        # spike at 2×readlen; the scorer queries the pmf at the genomic footprint
        # (span), so the pool must bin at the span. L stays the mass denominator.
        fl_on = (not spliced) and self.fl_pool_mass is not None
        if fl_on:
            fl_span = max(int(be) for _, be in blocks) - min(int(bs) for bs, _ in blocks)
            fl_bin = min(int(fl_span), self.max_fl)
        else:
            fl_bin = 0

        # 3. Single-region (all slices in same region) → contained.
        regions_touched = {sl[0] for sl in slices}
        if len(regions_touched) == 1:
            r = slices[0][0]
            self.regions[r].contained[ch] += np.uint32(1)
            if fl_on:
                pool = fl_pool_idx(int(self.region_types[r]), boundary=False)
                self.fl_pool_mass[pool, fl_bin] += 1.0
            return

        # 4. Crossing path. Distribute each slice's mass across the boundaries it
        #    crosses, conserving fragment mass (§4.3 of 00_design.md). A slice
        #    crosses its LEFT boundary iff it is not the first slice, and its
        #    RIGHT boundary iff it is not the last. A region the fragment
        #    *encompasses* — a fully-traversed interior slice that overlaps BOTH
        #    its boundaries — therefore splits its mass 50/50: half to the RIGHT
        #    side of its left boundary (mass_right) and half to the LEFT side of
        #    its right boundary (mass_left). End slices keep full mass on their
        #    single crossed side. Flux is the integer per-side crossing-event
        #    count (NOT split): the left region's slice credits flux_left of its
        #    right boundary; the right region's slice credits flux_right of its
        #    left boundary. A contiguous crossing credits both sides of its
        #    shared boundary; a spliced jump credits one side of each flanking
        #    boundary (intron-facing sides stay zero). Slices are monotonic →
        #    each side is credited at most once per fragment.
        n_slices = len(slices)
        for i, (r, start, end) in enumerate(slices):
            crosses_left = i > 0
            crosses_right = i < n_slices - 1
            n_cross = (1 if crosses_left else 0) + (1 if crosses_right else 0)
            if n_cross == 0:
                continue  # defensive: single-region fragments handled above
            share = (end - start) * inv_L / n_cross
            if crosses_right:
                b_out = self.right_boundary_of(r)
                self.boundaries[b_out].mass_left[ch] += np.float32(share)
                self.boundaries[b_out].flux_left[ch] += np.uint32(1)
                if js != 0:  # spliced ⇒ this boundary is a junction of strand `js`
                    self.boundaries[b_out].junction_strand = js
            if crosses_left:
                b_in = self.left_boundary_of(r)
                self.boundaries[b_in].mass_right[ch] += np.float32(share)
                self.boundaries[b_in].flux_right[ch] += np.uint32(1)
                if js != 0:
                    self.boundaries[b_in].junction_strand = js

        # gDNA FL (crossing): each slice's fractional mass → the BOUNDARY pool of
        # its region-type, at FL bin = min(footprint, max_fl).
        if fl_on:
            for r, start, end in slices:
                pool = fl_pool_idx(int(self.region_types[r]), boundary=True)
                self.fl_pool_mass[pool, fl_bin] += (end - start) * inv_L

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
