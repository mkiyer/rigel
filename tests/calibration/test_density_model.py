"""Gates for `density_model.count_observable_masks` — WHICH counts are gDNA counts.

The module's only surviving public function, and the gDNA strand-overdispersion fit's seed selector. The
rule is signature arithmetic and nothing else:

* a REGION is count-observable iff it carries **no exon bit** — an exon's contained count holds unspliced
  mature RNA, so it is not a gDNA count;
* a contiguous BOUNDARY is count-observable iff its two flanks **share no exon bit**. Sharing is the whole
  test: one exon-strand continuing across the boundary is what lets unspliced mature RNA cross it. Two
  exons on OPPOSITE strands share no bit, so nothing continues across and the crossing count is gDNA.

⚠ The per-region gDNA DENSITY this module used to impute (flank anchoring, run-fill, global baseline) was
deleted on 2026-08-30 with `run_fill`: it existed only to weight the strand fit's seeds, and the away-half
moment needs no weight. Its tests went with it; these gates cover what remains.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.density_model import count_observable_masks
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
)


def _masks(signatures, ref_id=None):
    sig = np.asarray(signatures, dtype=np.uint8)
    ref = (
        np.zeros(sig.shape[0], dtype=np.int64)
        if ref_id is None
        else np.asarray(ref_id, dtype=np.int64)
    )
    return count_observable_masks(sig, ref)


def test_an_exon_region_is_not_count_observable_and_intron_and_intergenic_are():
    """The region rule, on all four kinds: intergenic (0) and intron are gDNA counts; either exon strand,
    or an exon overlapping an intron, is not."""
    region, _ = _masks(
        [
            0,
            BIT_INTRON_POS,
            BIT_INTRON_NEG,
            BIT_EXON_POS,
            BIT_EXON_NEG,
            BIT_EXON_POS | BIT_INTRON_NEG,
        ]
    )
    np.testing.assert_array_equal(region, [True, True, True, False, False, False])


def test_a_boundary_is_observable_unless_an_EXON_STRAND_CONTINUES_across_it():
    """The boundary rule is SHARING, not presence. exon+|exon+ shares a bit — mature RNA crosses, not
    observable. exon+|intron+ shares none — the exon ends there, so the crossing count is gDNA."""
    _, boundary = _masks([BIT_EXON_POS, BIT_EXON_POS, BIT_INTRON_POS, 0])
    np.testing.assert_array_equal(boundary, [False, True, True])


def test_two_exons_on_OPPOSITE_strands_share_no_bit_so_the_boundary_IS_observable():
    """⚠ The tell that the rule is `sig[lo] & sig[hi] & EXON_BITS`, not `either flank is an exon`: no single
    exon-strand continues across an exon+|exon− junction, so no unspliced mature RNA crosses it."""
    _, boundary = _masks([BIT_EXON_POS, BIT_EXON_NEG])
    np.testing.assert_array_equal(boundary, [True])


def test_masks_are_on_their_own_axes_and_no_boundary_straddles_a_reference():
    """N regions, and E = N − (number of references) boundaries: two single-region references own ZERO
    boundaries between them, so nothing can leak across a reference edge."""
    region, boundary = _masks([BIT_INTRON_POS, BIT_INTRON_POS, BIT_EXON_POS], ref_id=[0, 1, 1])
    assert region.shape == (3,)
    assert boundary.shape == (1,)  # only the (1, 2) pair inside reference 1
    np.testing.assert_array_equal(boundary, [True])  # intron+ | exon+ shares no bit
