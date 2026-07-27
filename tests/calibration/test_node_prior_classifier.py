"""Phase A of the node-prior redesign (`docs/calibration/archive/node_prior_design.md` §3, §8): the
`mrna_active_*` classification carried onto the chain by `build_node_statics`, alongside the existing
`free_*` (nascent-active / RNA-crossing) masks. Behaviour-neutral — no solver change — so this only
asserts the new masks are correct across every region type and all four boundary types.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain
from rigel.calibration.node_geometry import build_node_statics
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
    transcript_strand_class,
)


def _cview(n):
    """A zero-count per-region/side substrate view of length ``n`` (counts/masses are irrelevant to the
    signature-only classifier under test)."""
    z = np.zeros(n)
    return SimpleNamespace(
        n_unspliced_pos=z.copy(),
        n_unspliced_neg=z.copy(),
        mass_unspliced=z.copy(),
        mass_spliced=z.copy(),
    )


def _build_statics(region_sigs):
    """Build a single-reference chain over ``region_sigs`` (genomic order) and return
    ``(chain, statics)``. Boundaries are the ``N+1`` seams, flanks resolved from the region signatures."""
    sig = np.asarray(region_sigs, dtype=np.int64)
    n_reg = sig.shape[0]
    region_arrays = SimpleNamespace(
        strand_class=transcript_strand_class(sig),
        signature=sig,
        region_size_bp=np.full(n_reg, 1000.0),
    )
    substrate = SimpleNamespace(contained=_cview(n_reg))
    lr = np.arange(-1, n_reg, dtype=np.int64)  # [-1, 0, 1, ..., n_reg-1]
    rr = np.arange(0, n_reg + 1, dtype=np.int64)
    rr[-1] = -1  # [0, 1, ..., n_reg-1, -1]
    n_bnd = lr.shape[0]
    bsub = SimpleNamespace(
        left_region=lr,
        right_region=rr,
        left=_cview(n_bnd),
        right=_cview(n_bnd),
        junction_strand=np.zeros(n_bnd, dtype=np.int8),
    )
    chain = build_node_chain(np.array([0, n_reg]), np.array([0, n_bnd]))
    return chain, build_node_statics(chain, substrate, bsub, region_arrays)


def test_classifier_covers_region_and_boundary_types():
    # R0 exon+ | R1 exon+ | R2 intron+ | R3 intergenic | R4 ambig-exon | R5 ambig-exon
    # boundaries (7): B0(−1|R0) B1(R0|R1) B2(R1|R2) B3(R2|R3) B4(R3|R4) B5(R4|R5) B6(R5|−1)
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
        assert m.dtype == bool and m.shape[0] == st.n_nodes
    assert np.all(~st.mrna_active_pos | st.free_pos)  # mature ⇒ nascent-active (+)
    assert np.all(~st.mrna_active_neg | st.free_neg)  # (−)

    kind, ref = np.asarray(chain.kind), np.asarray(chain.ref_idx)
    reg = np.where(kind == REGION)[0]  # R0..R5 (genomic order)
    bnd = np.where(kind == BOUNDARY)[0]  # B0..B6
    np.testing.assert_array_equal(ref[reg], np.arange(6))  # confirm genomic ordering
    np.testing.assert_array_equal(ref[bnd], np.arange(7))

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

    # --- boundaries: the four types ---
    assert state(bnd[1]) == (True, False, True, False)  # exon↔exon+   : MATURE-CAPABLE
    assert state(bnd[2]) == (True, False, False, False)  # exon↔intron+ : NASCENT-ONLY
    assert state(bnd[3]) == (False, False, False, False)  # intron↔intergenic: SINK (no + crossing)
    assert state(bnd[0]) == (False, False, False, False)  # terminal ↔ exon : SINK
    assert state(bnd[5]) == (True, True, True, True)  # ambig↔ambig  : AMBIG, mature both strands
