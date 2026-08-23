"""STAGE 0 OF THE FIRST-PASS REDESIGN — the structural claims, falsified clause by clause.

`rigel.calibration.structural_claims` derives, from the chain and the statics alone, the set of slots
whose UNSPLICED population is structurally constrained by the annotation — the first pass's training
and solving substrate. Four classes, each carrying a CLAIM the panel's certified slot truth can test
with no solver:

* ``intergenic``          — no RNA strand admissible (``g1_locked``, both kinds): truth must have
                            ``n_nrna == n_mrna == 0``.
* ``ss_intron_region``    — exactly one strand admissible, no exon membership: no contiguous mature
                            RNA fits inside, so truth must have ``n_mrna == 0``.
* ``ss_intron_boundary``  — exactly one strand continuous, no contiguous exon across: an unspliced
                            crossing has no mature term, so truth must have ``n_mrna == 0``.
* ``solvable_exon``       — single-stranded exonic REGION with at least one flanking BOUNDARY that is
                            a splice site with no contiguous exon across it. The claim is the FLANK's.

⭐ The fixture is DESIGNED so every clause has a positive and a negative case: a clean two-exon
transcript (locus A), a single-exon transcript whose flanks are termini (locus B), an antisense
overlap making every slot AMBIG (locus C), and a retained-intron isoform making mature RNA contiguous
across a donor (locus D). Each test is one gate; the perturbation sweep breaks one clause of the fixed
code at a time and watches exactly that gate fire.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain
from rigel.calibration.region_geometry import build_region_statics, g1_locked
from rigel.calibration.splice_graph import (
    FLAG_TES_NEG,
    FLAG_TES_POS,
    FLAG_TSS_NEG,
    FLAG_TSS_POS,
    build_boundary_flags_array,
    build_region_partition_arrays,
    is_splice_site,
)
from rigel.calibration.structural_claims import build_structural_claims
from rigel.types import Strand

from conftest import build_test_index

#: Locus A (+): tA splices [1000,1500) — both exons flank a splice site with nothing contiguous
#: across, so both are solvable and the intron between them is single-stranded.
#: Locus B (+): tB is single-exon — its flanks are TSS/TES termini, which license nothing.
#: Locus C (±): tCp splices [3800,3900) and tCm splices the overlapping [3780,3920) antisense (offset,
#: because identical coordinates on both strands are biologically impossible and the builder warns), so
#: every locus-C slot is AMBIG — including a pure AMBIG intron with no exon membership, the slot only
#: the strand clause excludes (the case that caught an XOR→OR perturbation the exon bits shielded).
#: Locus D (+): tD1 splices [5000,5300); tD2 retains that intron, so mature RNA is contiguous across
#: both splice sites and the "intron" carries exon membership.
#: Locus F (−): tF is locus A on the MINUS strand — without it, every single-stranded slot is ``+``
#: and ``~free_pos`` masquerades as ``g1_locked`` (a perturbation this fixture was silent on).
#: Locus H (chr3, +): tH2 is single-exon and STARTS exactly at tH1's acceptor, so the boundary at 700
#: is a licensed splice-site flank whose account of the exon is INCOMPLETE — tH2's molecules originate
#: at the flank and never pass through it. The completeness theorem's canonical negative case.
GTF = """\
chr1\ttest\texon\t501\t1000\t.\t+\t.\tgene_id "gA"; transcript_id "tA";
chr1\ttest\texon\t1501\t2000\t.\t+\t.\tgene_id "gA"; transcript_id "tA";
chr1\ttest\texon\t2501\t3000\t.\t+\t.\tgene_id "gB"; transcript_id "tB";
chr1\ttest\texon\t3501\t3800\t.\t+\t.\tgene_id "gC"; transcript_id "tCp";
chr1\ttest\texon\t3901\t4200\t.\t+\t.\tgene_id "gC"; transcript_id "tCp";
chr1\ttest\texon\t3501\t3780\t.\t-\t.\tgene_id "gC"; transcript_id "tCm";
chr1\ttest\texon\t3921\t4200\t.\t-\t.\tgene_id "gC"; transcript_id "tCm";
chr1\ttest\texon\t4501\t5000\t.\t+\t.\tgene_id "gD"; transcript_id "tD1";
chr1\ttest\texon\t5301\t5600\t.\t+\t.\tgene_id "gD"; transcript_id "tD1";
chr1\ttest\texon\t4501\t5600\t.\t+\t.\tgene_id "gD"; transcript_id "tD2";
chr1\ttest\texon\t6501\t7000\t.\t-\t.\tgene_id "gF"; transcript_id "tF";
chr1\ttest\texon\t7301\t7600\t.\t-\t.\tgene_id "gF"; transcript_id "tF";
chr2\ttest\texon\t1\t400\t.\t+\t.\tgene_id "gE"; transcript_id "tE";
chr2\ttest\texon\t701\t1000\t.\t+\t.\tgene_id "gE"; transcript_id "tE";
chr3\ttest\texon\t101\t400\t.\t+\t.\tgene_id "gH"; transcript_id "tH1";
chr3\ttest\texon\t701\t1000\t.\t+\t.\tgene_id "gH"; transcript_id "tH1";
chr3\ttest\texon\t701\t1000\t.\t+\t.\tgene_id "gH"; transcript_id "tH2";
"""

#: The partition the GTF above must produce — a fixture-shape guard, so every hand-enumerated
#: expectation below is anchored to verified geometry rather than to an assumption about the builder.
#: ⭐ Locus E (chr2, +) starts at base 1, so its first exon's LEFT flank is a reference terminal
#: (``chain.left == -1``) — the case that falsifies the flank gather's sentinel.
EXPECTED_BOUNDS = {
    "chr1": [
        500,
        1000,
        1500,
        2000,
        2500,
        3000,
        3500,
        3780,
        3800,
        3900,
        3920,
        4200,
        4500,
        5000,
        5300,
        5600,
        6500,
        7000,
        7300,
        7600,
    ],
    "chr2": [400, 700, 1000],
    "chr3": [100, 400, 700, 1000],
}


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    return build_test_index(
        tmp_path_factory,
        GTF,
        name="structural_claims",
        refs={"chr1": 8000, "chr2": 1500, "chr3": 1500},
    )


@pytest.fixture(scope="module")
def built(index):
    """(chain, statics, claims) on the designed index."""
    ra = RegionArrays.from_index(index)
    per_ref = np.bincount(np.asarray(ra.ref_id, np.int64))
    region_offsets = np.concatenate([[0], np.cumsum(per_ref)])
    boundary_offsets = np.concatenate([[0], np.cumsum(np.maximum(per_ref - 1, 0))])
    chain = build_region_chain(region_offsets, boundary_offsets)
    statics = build_region_statics(chain, ra, build_boundary_flags_array(index))
    return chain, statics, build_structural_claims(chain, statics)


def _region_slots(chain) -> np.ndarray:
    """Slot id of region ``i``, in genomic order — the chain is N0 E0 N1 …, so obj order is slot order."""
    return np.flatnonzero(np.asarray(chain.kind) == REGION)


def _boundary_slots(chain) -> np.ndarray:
    return np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)


def test_fixture_shape_is_the_designed_partition(index):
    """Guard: the enumerations below index regions/boundaries by genomic order, so the partition must
    be exactly the designed one. A failure here means the FIXTURE moved, not the predicate."""
    positions, offsets, _types = build_region_partition_arrays(index)
    for f, name in enumerate(("chr1", "chr2", "chr3")):
        lo, hi = int(offsets[f]), int(offsets[f + 1])
        interior = positions[lo + 1 : hi - 1]
        np.testing.assert_array_equal(interior, EXPECTED_BOUNDS[name], err_msg=name)


def test_hand_enumerated_membership(built):
    """The four masks match the designed loci exactly, slot by slot."""
    chain, _statics, claims = built
    r, b = _region_slots(chain), _boundary_slots(chain)

    want_intergenic_regions = [0, 4, 6, 12, 16, 20, 24, 25, 29]
    want_intergenic_boundaries = [0, 3, 4, 5, 6, 11, 12, 15, 16, 19, 22, 23, 26]
    want = np.zeros(chain.n_slots, bool)
    want[r[want_intergenic_regions]] = True
    want[b[want_intergenic_boundaries]] = True
    np.testing.assert_array_equal(claims.intergenic, want, err_msg="intergenic")

    want = np.zeros(chain.n_slots, bool)
    want[r[[2, 18, 22, 27]]] = True  # locus A's, F's, E's and H's introns, and nothing else
    np.testing.assert_array_equal(claims.ss_intron_region, want, err_msg="ss_intron_region")

    want = np.zeros(chain.n_slots, bool)
    want[b[[1, 2, 17, 18, 20, 21, 24, 25]]] = True  # those loci's donors and acceptors only
    np.testing.assert_array_equal(claims.ss_intron_boundary, want, err_msg="ss_intron_boundary")

    want = np.zeros(chain.n_slots, bool)
    want[r[[1, 3, 17, 19, 21, 23, 26, 28]]] = True  # locus A's, F's, E's and H's exons only
    np.testing.assert_array_equal(claims.solvable_exon, want, err_msg="solvable_exon")


def test_masks_are_disjoint_and_kind_scoped(built):
    """The classes partition the substrate: pairwise disjoint, each on its own axis, and ``claimed``
    is their union."""
    chain, _statics, claims = built
    stack = np.stack(
        [
            claims.intergenic,
            claims.ss_intron_region,
            claims.ss_intron_boundary,
            claims.solvable_exon,
        ]
    )
    assert (stack.sum(0) <= 1).all(), "a slot is claimed by two classes at once"
    is_region = np.asarray(chain.kind) == REGION
    assert not (claims.ss_intron_region & ~is_region).any()
    assert not (claims.solvable_exon & ~is_region).any()
    assert not (claims.ss_intron_boundary & is_region).any()
    np.testing.assert_array_equal(claims.claimed, stack.any(0))


def test_solvable_exon_names_its_licensing_flank(built):
    """Locus A's first exon is licensed by its RIGHT flank only (the donor; its left is the TSS) and
    the second by its LEFT flank only — and every named flank really is a splice-site BOUNDARY with no
    contiguous exon across it, which is the claim the confusion matrix will test."""
    chain, statics, claims = built
    r = _region_slots(chain)
    assert not claims.exon_flank_left[r[1]] and claims.exon_flank_right[r[1]]
    assert claims.exon_flank_left[r[3]] and not claims.exon_flank_right[r[3]]
    # locus F, the same shape on the MINUS strand: the acceptor/donor bits land on the other flags,
    # but the licensing geometry is strand-symmetric.
    assert not claims.exon_flank_left[r[17]] and claims.exon_flank_right[r[17]]
    assert claims.exon_flank_left[r[19]] and not claims.exon_flank_right[r[19]]
    # locus E's first exon: its left flank is a REFERENCE TERMINAL (chain.left == -1), which must not
    # license anything — only the right flank (the donor) does. The sentinel clause's gate.
    assert not claims.exon_flank_left[r[21]] and claims.exon_flank_right[r[21]]
    assert claims.exon_flank_left[r[23]] and not claims.exon_flank_right[r[23]]
    # locus H: the shared exon is licensed ONLY via its left flank (the acceptor tH2's TSS shares).
    assert claims.exon_flank_left[r[28]] and not claims.exon_flank_right[r[28]]

    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    mrna_b = np.asarray(statics.mrna_active_pos) | np.asarray(statics.mrna_active_neg)
    sj = is_splice_site(statics.boundary_flags, Strand.POS) | is_splice_site(
        statics.boundary_flags, Strand.NEG
    )
    for side, flank in ((claims.exon_flank_left, left), (claims.exon_flank_right, right)):
        for s in np.flatnonzero(side):
            f = int(flank[s])
            assert f >= 0 and np.asarray(chain.kind)[f] == BOUNDARY
            assert sj[f] and not mrna_b[f]
    assert (claims.exon_flank_left | claims.exon_flank_right)[claims.solvable_exon].all()
    assert not (claims.exon_flank_left | claims.exon_flank_right)[~claims.solvable_exon].any()


def test_ambig_is_excluded_everywhere(built):
    """Locus C: the antisense overlap makes both strands admissible, and an AMBIG slot has no channel —
    it appears in NO class (the deferred stratum's blindness, excluded from the training substrate).
    ⭐ The shared intron [3800,3900) carries NO exon membership, so only the single-strandedness clause
    stands between it and ``ss_intron_region`` — the slot that catches an XOR→OR regression."""
    chain, _statics, claims = built
    r, b = _region_slots(chain), _boundary_slots(chain)
    locus_c = np.concatenate([r[[7, 8, 9, 10, 11]], b[[7, 8, 9, 10]]])
    assert not claims.claimed[locus_c].any()


def test_a_terminus_flank_does_not_license(built):
    """Locus B: a single-exon transcript's flanks are TSS/TES termini — ``is_splice_site`` is False
    there, and a terminus satisfies "no contiguous exon" without carrying an sj, so it must not
    license the exon."""
    chain, _statics, claims = built
    r = _region_slots(chain)
    assert not claims.claimed[r[5]]


def test_a_contiguous_exon_disqualifies_the_flank(built):
    """Locus D: the retained-intron isoform makes mature RNA contiguous across both splice sites, so
    neither exon of tD1 is solvable, the retained intron is not an intron (it carries exon
    membership), and the two splice-site boundaries carry a mature crossing term."""
    chain, _statics, claims = built
    r, b = _region_slots(chain), _boundary_slots(chain)
    assert not claims.claimed[np.concatenate([r[[13, 14, 15]], b[[13, 14]]])].any()


def test_brute_force_equivalence(built):
    """The vectorized module against an independent per-slot loop over the same statics — the working
    rule's brute-force enumeration, so a vectorization defect cannot hide in broadcasting."""
    chain, statics, claims = built
    kind = np.asarray(chain.kind)
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    fp, fn = np.asarray(statics.free_pos, bool), np.asarray(statics.free_neg, bool)
    mp, mn = (
        np.asarray(statics.mrna_active_pos, bool),
        np.asarray(statics.mrna_active_neg, bool),
    )
    sj = is_splice_site(statics.boundary_flags, Strand.POS) | is_splice_site(
        statics.boundary_flags, Strand.NEG
    )
    flags = np.asarray(statics.boundary_flags, np.uint16)
    low_end = (flags & (FLAG_TSS_POS | FLAG_TES_NEG)) != 0  # a transcript's genomic-LOW end here
    high_end = (flags & (FLAG_TES_POS | FLAG_TSS_NEG)) != 0

    for s in range(chain.n_slots):
        g1 = not fp[s] and not fn[s]
        ss = bool(fp[s]) != bool(fn[s])
        mrna = bool(mp[s]) or bool(mn[s])

        def qualifies(f: int) -> bool:
            if f < 0 or kind[f] != BOUNDARY:
                return False
            return bool(sj[f]) and not (mp[f] or mn[f])

        want_left = ss and kind[s] == REGION and mrna and qualifies(int(left[s]))
        want_right = ss and kind[s] == REGION and mrna and qualifies(int(right[s]))
        assert claims.intergenic[s] == g1, s
        assert claims.ss_intron_region[s] == (ss and kind[s] == REGION and not mrna), s
        assert claims.ss_intron_boundary[s] == (ss and kind[s] == BOUNDARY and not mrna), s
        assert claims.exon_flank_left[s] == want_left, s
        assert claims.exon_flank_right[s] == want_right, s
        assert claims.solvable_exon[s] == (want_left or want_right), s
        # the completeness theorem, re-derived per slot: a molecule missing the flank must END at it
        assert claims.exon_flank_left_complete[s] == (want_left and not low_end[int(left[s])]), s
        assert claims.exon_flank_right_complete[s] == (
            want_right and not high_end[int(right[s])]
        ), s


def test_intergenic_is_exactly_g1_locked(built):
    """The intergenic class is ``g1_locked`` and never a restatement of it — one predicate, one home
    (its docstring records why a second home is how region-only variants survive)."""
    _chain, statics, claims = built
    np.testing.assert_array_equal(claims.intergenic, g1_locked(statics.free_pos, statics.free_neg))


def test_flank_completeness_hand_enumerated(built):
    """The completeness theorem on the designed loci: every licensed flank is COMPLETE — a molecule
    overlapping the exon either passes through the flank or ends exactly at it, and none of these
    flanks carries a facing end bit — EXCEPT locus H's, where tH2's TSS coincides with the acceptor,
    so molecules originate at the flank unseen and the flank's account is INCOMPLETE."""
    chain, _statics, claims = built
    r = _region_slots(chain)
    complete = claims.exon_flank_left_complete | claims.exon_flank_right_complete
    for slot in r[[1, 3, 17, 19, 21, 23, 26]]:
        assert complete[slot], f"slot {slot} must have a complete licensed flank"
    assert not complete[r[28]], "locus H's exon has no complete flank — TSS at the acceptor"
    assert claims.exon_flank_left[r[28]], "…but the flank is still LICENSED (bound-grade)"


def test_completeness_is_a_subset_of_the_licence(built):
    """A completeness bit may only be set where the flank is licensed at all — completeness qualifies
    a claim; it never creates one."""
    _chain, _statics, claims = built
    assert not (claims.exon_flank_left_complete & ~claims.exon_flank_left).any()
    assert not (claims.exon_flank_right_complete & ~claims.exon_flank_right).any()
