"""Phase A of the region-prior redesign (§8): the
`mrna_active_*` classification carried onto the chain by `build_region_statics`, alongside the existing
`free_*` (nascent-active / RNA-crossing) masks. Behaviour-neutral — no solver change — so this only
asserts the new masks are correct across every region type and all four boundary types.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain
from rigel.calibration.region_geometry import build_region_statics
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
    transcript_strand_class,
)


def _substrate(n_regions, n_boundaries):
    """A zero-count substrate of the right shape — the classifier under test is signature-only, so the
    counts are irrelevant and only ``boundary_spliced`` is read (for ``spliced_count``)."""

    def view(n):
        return SimpleNamespace(count=np.zeros((n, 2)))

    return SimpleNamespace(
        region_contained=view(n_regions), boundary_unspliced=view(n_boundaries), boundary_spliced=view(n_boundaries)
    )


def _build_statics(region_sigs):
    """Build a single-reference chain over ``region_sigs`` (genomic order) and return
    ``(chain, statics)``.

    ⭐ ``k`` regions own ``k − 1`` interior boundaries and there are **no terminal slots** — so the old
    "boundary 0 is the reference-start sink" case does not exist to be classified any more. An boundary
    always has a region on both sides, which is what lets `build_region_statics` read its flanks straight
    off the chain's adjacency instead of a ``left_region``/``right_region`` array with ``-1`` holes.
    """
    sig = np.asarray(region_sigs, dtype=np.int64)
    n_reg = sig.shape[0]
    region_arrays = SimpleNamespace(
        strand_class=transcript_strand_class(sig),
        signature=sig,
        region_size_bp=np.full(n_reg, 1000.0),
    )
    n_boundary = max(n_reg - 1, 0)
    chain = build_region_chain(np.array([0, n_reg]), np.array([0, n_boundary]))
    return chain, build_region_statics(chain, region_arrays)


def test_classifier_covers_region_and_boundary_types():
    # N0 exon+ | N1 exon+ | N2 intron+ | N3 intergenic | N4 ambig-exon | N5 ambig-exon
    # boundaries (5): E0(N0|N1) E1(N1|N2) E2(N2|N3) E3(N3|N4) E4(N4|N5) — no terminals
    sigs = [
        BIT_EXON_POS,
        BIT_EXON_POS,
        BIT_INTRON_POS,
        0,
        BIT_EXON_POS | BIT_EXON_NEG,
        BIT_EXON_POS | BIT_EXON_NEG,
    ]
    chain, st = _build_statics(sigs)

    # masks are bool and full-length; the whole-chain invariant mrna_active ⇒ free (nascent) holds.
    for m in (st.free_pos, st.free_neg, st.mrna_active_pos, st.mrna_active_neg):
        assert m.dtype == bool and m.shape[0] == st.n_slots
    assert np.all(~st.mrna_active_pos | st.free_pos)  # mature ⇒ nascent-active (+)
    assert np.all(~st.mrna_active_neg | st.free_neg)  # (−)

    kind, ref = np.asarray(chain.kind), np.asarray(chain.obj_idx)
    reg = np.where(kind == REGION)[0]  # N0..N5 (genomic order)
    bnd = np.where(kind == BOUNDARY)[0]  # E0..E4
    np.testing.assert_array_equal(ref[reg], np.arange(6))  # confirm genomic ordering
    np.testing.assert_array_equal(ref[bnd], np.arange(5))

    def state(i):
        return (
            bool(st.free_pos[i]),
            bool(st.free_neg[i]),
            bool(st.mrna_active_pos[i]),
            bool(st.mrna_active_neg[i]),
        )

    # --- regions (free_pos, free_neg, mrna_pos, mrna_neg) ---
    assert state(reg[0]) == (True, False, True, False)  # exon+  : mature-capable +
    assert state(reg[2]) == (
        True,
        False,
        False,
        False,
    )  # intron+: NASCENT-ONLY + (free but not mature)
    assert state(reg[3]) == (False, False, False, False)  # intergenic: gDNA sink
    assert state(reg[4]) == (True, True, True, True)  # ambig-exon: mature-capable both strands

    # --- boundaries: the four types. ⚠ Indices shift by one against the predecessor because the
    # reference-start terminal slot no longer exists; E0 is now the first REAL boundary, N0|N1.
    assert state(bnd[0]) == (True, False, True, False)  # exon↔exon+   : MATURE-CAPABLE
    assert state(bnd[1]) == (True, False, False, False)  # exon↔intron+ : NASCENT-ONLY
    assert state(bnd[2]) == (False, False, False, False)  # intron↔intergenic: SINK (no + crossing)
    assert state(bnd[3]) == (False, False, False, False)  # intergenic↔ambig-exon : SINK
    assert state(bnd[4]) == (True, True, True, True)  # ambig↔ambig  : AMBIG, mature both strands
