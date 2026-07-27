"""Guard: spliced boundary mass is ONE-SIDED (exon side only) — and must stay that way.

Spliced (mature) RNA splices AT an exon↔intron boundary: there is spliced mass on the EXON side and
ZERO on the intron-facing side. Unspliced reads (gDNA + unspliced nascent) contiguously span a boundary
and so sit on BOTH sides. Spliced and unspliced live in SEPARATE accumulator channels (c2/c3 vs c0/c1)
on the same boundary object.

The whole calibration relies on this one-sidedness being correct and on the channel separation:

  * the gDNA-component effective length (``priors.assemble_priors``, the **pooled-seam** design — pools
    ``mass_gdna_right[r] + mass_gdna_left[r+1]``) is valid ONLY because gDNA is the UNSPLICED channel and
    is genuinely two-sided. Pooling the two sides of a seam is correct for gDNA; it would be WRONG for
    spliced/mature RNA, which is one-sided (a symmetric pool would mis-locate it / silently rely on a side
    that is structurally zero).
  * ``splice_junction.region_splice_gdna_frac`` reads spliced flux from the predicate-named EXON side
    only and *expects* the intron-facing side to be zero.

This test pins that contract at the substrate level (where the calibration consumers read it), so that a
future change — e.g. "fix" the accumulator to deposit spliced two-sided, or pool spliced across a seam the
way gDNA is pooled — fails loudly instead of silently corrupting the splice-junction debias. The
accumulator-level twin is ``tests/native/test_accumulator_spec.py`` T4. Design:
``docs/calibration/archive/splice_junction_node_architecture.md`` (§1, §5 Phase 0).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.scan_payload import AccumulatorPayload

# Channel layout (the separation the guard protects): c0 unspliced+, c1 unspliced−, c2 spliced sense,
# c3 spliced antisense.
_C_UNSPL_POS, _C_UNSPL_NEG, _C_SPL_SENSE, _C_SPL_ANTI = 0, 1, 2, 3


def _exon_intron_exon_payload() -> tuple[AccumulatorPayload, RegionArrays]:
    """A +strand exon→intron→exon locus with a spliced read across the junction.

    Regions r0 (+exon, donor), r1 (+intron), r2 (+exon, acceptor); boundaries b0(term), b1 (donor seam
    r0|r1), b2 (acceptor seam r1|r2), b3(term). At each seam: gDNA (unspliced, c0) is two-sided; the
    spliced mature mass (c2) is on the EXON side ONLY (the intron-facing side is zero), as the accumulator
    deposits it (T4). ``mass_left`` is the side LEFT of the boundary; ``mass_right`` the side RIGHT.
    """
    region_contained = np.array(
        [[4, 0, 2, 0], [1, 0, 0, 0], [4, 0, 2, 0]],
        dtype=np.uint32,  # r1 intron: sparse gDNA, no spliced
    )
    flux_left = np.zeros((4, 4), dtype=np.uint32)
    flux_right = np.zeros((4, 4), dtype=np.uint32)
    mass_left = np.zeros((4, 4), dtype=np.float32)
    mass_right = np.zeros((4, 4), dtype=np.float32)

    # b1 = donor seam (exon r0 on the LEFT, intron r1 on the RIGHT):
    #   gDNA two-sided; spliced on the EXON (left) side only, intron (right) side ZERO.
    mass_left[1] = [2.0, 0.0, 3.0, 0.0]  # exon r0 side: gDNA 2 + spliced 3
    mass_right[1] = [2.0, 0.0, 0.0, 0.0]  # intron r1 side: gDNA 2, spliced 0  ← the invariant
    flux_left[1] = [2, 0, 1, 0]
    flux_right[1] = [2, 0, 0, 0]

    # b2 = acceptor seam (intron r1 on the LEFT, exon r2 on the RIGHT):
    #   gDNA two-sided; spliced on the EXON (right) side only, intron (left) side ZERO.
    mass_left[2] = [2.0, 0.0, 0.0, 0.0]  # intron r1 side: gDNA 2, spliced 0  ← the invariant
    mass_right[2] = [2.0, 0.0, 3.0, 0.0]  # exon r2 side: gDNA 2 + spliced 3
    flux_left[2] = [2, 0, 0, 0]
    flux_right[2] = [2, 0, 1, 0]

    payload = AccumulatorPayload(
        boundaries=np.array([0, 100, 200, 300], dtype=np.int64),
        ref_pos_offsets=np.array([0, 4], dtype=np.int64),
        ref_region_offsets=np.array([0, 3], dtype=np.int64),
        ref_boundary_offsets=np.array([0, 4], dtype=np.int64),
        region_contained=region_contained,
        boundary_mass_left=mass_left,
        boundary_mass_right=mass_right,
        boundary_flux_left=flux_left,
        boundary_flux_right=flux_right,
        n_refs=1,
    )
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(3, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * 3, dtype="string"),
            "start": np.array([0, 100, 200], dtype=np.int64),
            "end": np.array([100, 200, 300], dtype=np.int64),
            "length": np.array([100, 100, 100], dtype=np.int64),
            "signature": np.array([BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS], dtype=np.uint8),
        }
    )
    return payload, RegionArrays.from_region_df(region_df, {"chr1": 0})


def test_spliced_boundary_mass_is_one_sided_at_the_intron():
    """The intron region's boundary views carry ZERO spliced mass; the exon sides carry it all."""
    payload, ra = _exon_intron_exon_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)
    R0, R1_INTRON, R2 = 0, 1, 2

    # The INTRON region's two boundary views (both facing exons) must have NO spliced mass — mature RNA
    # never crosses into an intron. This is the contract the splice-junction consumers rely on.
    assert sub.left.mass_spliced[R1_INTRON] == 0.0
    assert sub.right.mass_spliced[R1_INTRON] == 0.0
    assert sub.left.n_spliced[R1_INTRON] == 0
    assert sub.right.n_spliced[R1_INTRON] == 0

    # The EXON regions carry the spliced mass on their intron-facing views (donor: r0's right view;
    # acceptor: r2's left view) — one-sided, as deposited.
    assert sub.right.mass_spliced[R0] == 3.0  # donor seam, exon side
    assert sub.left.mass_spliced[R2] == 3.0  # acceptor seam, exon side


def test_seam_pooling_is_valid_for_gdna_but_spliced_stays_one_sided():
    """gDNA (unspliced) pools symmetrically across a seam; spliced does NOT (it is one-sided)."""
    payload, ra = _exon_intron_exon_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)
    R0, R1_INTRON, R2 = 0, 1, 2

    # The pooled-seam gDNA the eff-len uses = right side of left region + left side of right region.
    # gDNA is two-sided, so BOTH halves contribute (2 + 2 = 4) — pooling is meaningful.
    donor_pool_gdna = sub.right.mass_unspliced[R0] + sub.left.mass_unspliced[R1_INTRON]
    accept_pool_gdna = sub.right.mass_unspliced[R1_INTRON] + sub.left.mass_unspliced[R2]
    assert donor_pool_gdna == 4.0
    assert accept_pool_gdna == 4.0

    # If spliced were pooled the SAME way across a seam, only ONE side ever contributes — the other is
    # structurally zero. This is the footgun the guard names: a symmetric pool of spliced mass is invalid
    # (it does not double a junction's mass; it just rides on the lone exon side). Consumers must read the
    # exon side explicitly, not pool both sides.
    donor_pool_spliced = sub.right.mass_spliced[R0] + sub.left.mass_spliced[R1_INTRON]
    assert donor_pool_spliced == 3.0  # = exon side only; intron side contributed 0
    assert sub.left.mass_spliced[R1_INTRON] == 0.0  # the zero half is structural, not incidental


def test_gdna_channel_excludes_spliced_so_pooled_seam_efflen_is_immune():
    """gDNA mass (the pooled-seam eff-len input) is the UNSPLICED channel; spliced never enters it."""
    payload, ra = _exon_intron_exon_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)

    # The substrate keeps unspliced (c0/c1 → mass_unspliced) and spliced (c2/c3 → mass_spliced) in
    # SEPARATE reductions. The pooled-seam gDNA eff-len pools mass_unspliced; a regression that folded the
    # spliced channel into the unspliced one would change these whole-locus sums. Pin both, independently:
    total_spliced = (
        sub.contained.mass_spliced.sum()
        + sub.left.mass_spliced.sum()
        + sub.right.mass_spliced.sum()
    )
    total_unspliced = (
        sub.contained.mass_unspliced.sum()
        + sub.left.mass_unspliced.sum()
        + sub.right.mass_unspliced.sum()
    )
    # spliced = contained [2,0,2]=4 + left side [0,0,3]=3 + right side [3,0,0]=3  (the two junction halves)
    assert total_spliced == 4.0 + 3.0 + 3.0
    # unspliced (gDNA) = contained [4,1,4]=9 + left side [0,2,2]=4 + right side [2,2,0]=4 — independent of
    # the spliced mass entirely (no overlap, no double-count between the channels).
    assert total_unspliced == 9.0 + 4.0 + 4.0
