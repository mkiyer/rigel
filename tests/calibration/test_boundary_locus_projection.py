"""A locus collects REGIONS **and** BOUNDARIES — no boundary is ever folded into a region's total.

⭐⭐ **THE RULE** (owner, 2026-08-08): *an boundary owns the fragments that cross it; a region owns only the
fragments contained in it; nothing is re-attributed.* A locus's boundaries are therefore the boundaries that
**touch its regions** — every region contributes both of its boundaries, so a locus of ``k`` contiguous regions
carries ``k + 1`` boundaries, its two outer ones included.

⛔ **What this replaces.** ``boundary_owner_regions`` folded each boundary's mass into ONE flank region's total so
the region-overlap projection could reach it — a 0-bp boundary has no extent, and
``_project_regions_to_loci`` divides by ``region_size_bp``. That fold then needed a second heuristic (the
intergenic re-key) to stop a locus's far-LEFT boundary vanishing into the dropped intergenic flank. Both are
gone: an boundary is projected as an boundary.

⭐ **The outer boundaries are unambiguous, and that is structural, not lucky.** A locus is bounded by
intergenic sequence and intergenic regions carry no transcripts, so no two loci can contend for a
boundary boundary. Where two loci *could* contend — adjacent regions in different multi-loci — any fragment
crossing the shared boundary overlaps transcripts in both and would have MERGED them into one multi-locus.
The configuration is therefore unreachable for any boundary that carries mass, which is what
:func:`test_a_contended_boundary_carries_no_mass` states as a measurable claim rather than an assumption.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import _boundary_locus_shares, _region_locus_shares
from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.signature import BIT_EXON_POS
from rigel.locus import Locus, MultiLocus


def _regions(bounds, signature=None, ref_id=None) -> RegionArrays:
    """Regions tiling one (or more) references from the given region_bound positions."""
    bounds = np.asarray(bounds, dtype=np.int64)
    starts, ends = bounds[:-1], bounds[1:]
    n = starts.shape[0]
    rid = np.zeros(n, dtype=np.int32) if ref_id is None else np.asarray(ref_id, np.int32)
    n_refs = int(rid.max()) + 1 if n else 1
    offsets = np.searchsorted(rid, np.arange(n_refs + 1)).astype(np.int32)
    return RegionArrays(
        ref_id=rid,
        start=starts,
        end=ends,
        signature=(
            np.full(n, BIT_EXON_POS, dtype=np.uint8)
            if signature is None
            else np.asarray(signature, np.uint8)
        ),
        strand_class=np.zeros(n, dtype=np.int8),
        region_size_bp=(ends - starts).astype(np.float64),
        ref_offsets=offsets,
        n_refs=n_refs,
    )


def _ml(locus_id, blocks) -> MultiLocus:
    loci = tuple(Locus(ref=str(r), ref_id=r, start=s, end=e) for r, s, e in blocks)
    return MultiLocus(
        multi_locus_id=locus_id,
        transcript_indices=np.array([], dtype=np.int32),
        unit_indices=np.array([], dtype=np.int32),
        gdna_span=sum(e - s for _, s, e in blocks),
        loci=loci,
    )


def _dense(idx, lid, share, n_boundaries, n_loci) -> np.ndarray:
    """The (boundary, locus) share triples as a dense matrix — readable assertions, small fixtures only."""
    out = np.zeros((n_boundaries, n_loci), dtype=np.float64)
    out[np.asarray(idx, np.int64), np.asarray(lid, np.int64)] = share
    return out


# --- the rule -------------------------------------------------------------------------------------


def test_a_locus_of_k_regions_carries_k_plus_1_boundaries():
    """⭐ THE RULE, stated as a count. 5 regions; the locus is the middle 3, flanked by intergenic.

    Its boundaries are the 4 boundaries touching those 3 regions: the two interior ones AND the two OUTER ones.
    ⛔ A form that kept only interior boundaries would give 2, and a form that kept only one outer boundary —
    which is what left-keying without the re-key did — would give 3.
    """
    ra = _regions([0, 100, 200, 300, 400, 500])
    ml = [_ml(0, [(0, 100, 400)])]  # regions 1,2,3
    e, lid, w = _boundary_locus_shares(ra, ml, 1)
    assert sorted(e.tolist()) == [0, 1, 2, 3], "expected the 4 boundaries touching regions 1-3"
    np.testing.assert_allclose(w, 1.0)
    assert set(lid.tolist()) == {0}


def test_an_boundary_between_two_intergenic_regions_belongs_to_no_locus():
    """⛔ The complement, and it is what makes the rule a rule rather than 'keep everything'.

    A fragment crossing a boundary with no locus region on either side overlaps no transcript, so it is a
    candidate nowhere and must load no prior.
    """
    ra = _regions([0, 100, 200, 300, 400])
    ml = [_ml(0, [(0, 300, 400)])]  # region 3 only
    e, lid, _w = _boundary_locus_shares(ra, ml, 1)
    assert sorted(e.tolist()) == [2], "only the boundary touching region 3 may be kept"
    assert set(lid.tolist()) == {0}


def test_a_region_touching_ONE_locus_gives_it_everything():
    """⚠ **The share is an ALLOCATION, not an overlap fraction** — easy to misread, and load-bearing.

    ``_region_locus_shares`` normalises across the loci a region touches, so a region overlapping exactly
    one locus contributes **all** its mass there however small the overlap. That is the shipped
    behaviour and it is what makes the projection conserve: a region's mass is never partially discarded,
    only distributed.
    """
    ra = _regions([0, 100, 200])
    ml = [_ml(0, [(0, 0, 160)])]  # region 1 = [100,200) overlaps by only 60 bp
    r_idx, _lid, r_w = _region_locus_shares(ra, ml, 1)
    got = dict(zip(r_idx.tolist(), r_w.tolist()))
    assert got[0] == pytest.approx(1.0)
    assert got[1] == pytest.approx(1.0), "a 60 % overlap with the ONLY locus still allocates 100 %"


def test_an_boundary_takes_the_MAX_share_of_its_two_flanks():
    """A region genuinely split BETWEEN two loci has fractional shares; its boundaries take the larger.

    ⚠ ``max``, not a sum and not a mean: the rule is *"if a region is part of a locus, its two boundaries are
    part of that locus"*, so a boundary inherits the stronger of its two flanks' memberships.

    region 0 = [0,100) lies wholly in locus 0; region 1 = [100,200) is split 50/50 between loci 0 and 1.
    The single boundary between them is therefore ``max(1.0, 0.5) = 1.0`` in locus 0 and
    ``max(0.0, 0.5) = 0.5`` in locus 1.
    """
    ra = _regions([0, 100, 200])
    ml = [_ml(0, [(0, 0, 150)]), _ml(1, [(0, 150, 200)])]
    r_idx, r_lid, r_w = _region_locus_shares(ra, ml, 2)
    got = dict(zip(zip(r_idx.tolist(), r_lid.tolist()), r_w.tolist()))
    assert got[(0, 0)] == pytest.approx(1.0)
    assert got[(1, 0)] == pytest.approx(0.5) and got[(1, 1)] == pytest.approx(0.5)
    e, lid, w = _boundary_locus_shares(ra, ml, 2)
    assert e.tolist() == [0, 0]
    assert dict(zip(lid.tolist(), w.tolist())) == {0: pytest.approx(1.0), 1: pytest.approx(0.5)}


def test_a_contended_boundary_carries_no_mass():
    """⛔⛔ THE CLAIM THE RULE RESTS ON, AS A MEASUREMENT RATHER THAN AN ASSUMPTION.

    Two adjacent regions in DIFFERENT multi-loci would give their shared boundary to both, so its shares sum
    to 2. The owner's argument is that this is unreachable *for a boundary that carries mass*: any fragment
    crossing it overlaps transcripts in both loci, is a candidate in both, and the union-find would
    already have merged them into ONE multi-locus.

    ⭐ So the assertion is not "it cannot happen" but "where it happens, the mass is zero" — and the
    projection reports such boundaries rather than silently renormalising them.
    """
    ra = _regions([0, 100, 200])
    ml = [_ml(0, [(0, 0, 100)]), _ml(1, [(0, 100, 200)])]  # adjacent, no intergenic between
    e, lid, w = _boundary_locus_shares(ra, ml, 2)
    assert e.tolist() == [0, 0], "the contended boundary reaches both loci"
    assert sorted(lid.tolist()) == [0, 1]
    assert float(w.sum()) == pytest.approx(2.0), "shares sum above 1 — the reportable configuration"


def test_two_references_do_not_share_an_boundary():
    """A boundary exists only between two regions of the SAME reference; the axis must not straddle refs."""
    ra = _regions([0, 100, 200, 0, 100, 200][:4], ref_id=[0, 0, 1])
    lo, hi = boundary_region_indices(np.asarray(ra.ref_id))
    assert lo.tolist() == [0], "one boundary, inside reference 0 only"
    assert hi.tolist() == [1]


def test_empty_loci_returns_empty():
    ra = _regions([0, 100, 200])
    e, lid, w = _boundary_locus_shares(ra, [], 0)
    assert e.size == lid.size == w.size == 0


def test_the_region_projection_is_unchanged_by_the_refactor():
    """⛔ ``_project_regions_to_loci`` must keep its exact behaviour — the region half is NOT changing.

    It is now expressed through ``_region_locus_shares``, so this pins that the split introduced no
    drift: shares normalise across the loci a region touches, and a region touching none is dropped.
    """
    from rigel.calibration.priors import _project_regions_to_loci

    ra = _regions([0, 100, 200, 300])
    ml = [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])]  # region 0 split 50/50, regions 1-2 outside
    out = _project_regions_to_loci(ra, ml, 2, {"m": np.array([10.0, 99.0, 99.0])})
    np.testing.assert_allclose(out["m"], [5.0, 5.0])
