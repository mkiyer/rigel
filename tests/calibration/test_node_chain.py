"""The unified region↔boundary node chain (`calibration.node_chain`).

Locks: the genomic B-R-B-…-R-B interleave, the per-reference adjacency (terminals → -1), and consistency with
the existing region↔boundary index maps (a region's chain-neighbours ARE its flanking boundaries; a boundary's
chain-neighbours ARE its flanking regions).
"""

from __future__ import annotations

import numpy as np
from _synthetic import make_gdna_fl_pmf, make_synthetic_payload

from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain
from rigel.calibration.node_geometry import build_node_geometry
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.calibration.region_arrays import boundary_region_indices, region_boundary_indices


def test_interleave_and_adjacency_two_refs():
    # ref0: 2 regions (r0,r1) → boundaries b0,b1,b2 ; ref1: 1 region (r2) → boundaries b3,b4
    rro = np.array([0, 2, 3])
    rbo = np.array([0, 3, 5])
    ch = build_node_chain(rro, rbo)

    assert ch.n_nodes == 3 + 5  # R + B
    assert ch.n_regions == 3 and ch.n_boundaries == 5
    # genomic order ref0: B0 R0 B1 R1 B2  | ref1: B3 R2 B4
    assert list(ch.kind) == [
        BOUNDARY,
        REGION,
        BOUNDARY,
        REGION,
        BOUNDARY,
        BOUNDARY,
        REGION,
        BOUNDARY,
    ]
    assert list(ch.ref_idx) == [0, 0, 1, 1, 2, 3, 2, 4]
    # adjacency (node ids); reference terminals are -1
    assert list(ch.left) == [-1, 0, 1, 2, 3, -1, 5, 6]
    assert list(ch.right) == [1, 2, 3, 4, -1, 6, 7, -1]
    assert list(ch.order) == list(range(8))


def test_terminals_are_sinks():
    rro = np.array([0, 2, 3])
    rbo = np.array([0, 3, 5])
    ch = build_node_chain(rro, rbo)
    # the first node of each ref has no left, the last no right
    assert ch.left[0] == -1 and ch.left[5] == -1
    assert ch.right[4] == -1 and ch.right[7] == -1
    # both ref terminals are BOUNDARY nodes (the B-R-B-…-B shape)
    assert ch.kind[0] == BOUNDARY and ch.kind[4] == BOUNDARY
    assert ch.kind[5] == BOUNDARY and ch.kind[7] == BOUNDARY


def test_consistent_with_index_maps():
    # A richer topology: ref0 k=3, ref1 k=1, ref2 k=2.
    rro = np.array([0, 3, 4, 6])
    rbo = np.array([0, 4, 6, 9])
    ch = build_node_chain(rro, rbo)
    lb, rb = region_boundary_indices(rro, rbo)  # per region → flanking boundary array-indices
    blr, brr = boundary_region_indices(
        rro, rbo
    )  # per boundary → flanking region array-indices (-1 terminal)

    for node in range(ch.n_nodes):
        L, R = int(ch.left[node]), int(ch.right[node])
        if ch.kind[node] == REGION:
            r = int(ch.ref_idx[node])
            # a region's chain neighbours are BOUNDARY nodes == its flanking boundaries
            assert ch.kind[L] == BOUNDARY and ch.ref_idx[L] == lb[r]
            assert ch.kind[R] == BOUNDARY and ch.ref_idx[R] == rb[r]
        else:
            b = int(ch.ref_idx[node])
            # a boundary's chain neighbours are REGION nodes == its flanking regions (or -1 at a terminal)
            if blr[b] >= 0:
                assert L >= 0 and ch.kind[L] == REGION and ch.ref_idx[L] == blr[b]
            else:
                assert L == -1
            if brr[b] >= 0:
                assert R >= 0 and ch.kind[R] == REGION and ch.ref_idx[R] == brr[b]
            else:
                assert R == -1


def test_total_node_count_invariant():
    # n_nodes == R + B == 2R + n_refs (k+1 boundaries per k regions)
    rro = np.array([0, 3, 4, 6])
    rbo = np.array([0, 4, 6, 9])
    ch = build_node_chain(rro, rbo)
    n_refs = rro.shape[0] - 1
    assert ch.n_nodes == 2 * int(rro[-1]) + n_refs


def test_spliced_flux_plumbing_matches_mass_routing():
    """The integer spliced COUNT must be routed EXACTLY like the spliced MASS: same motif/exon gate, same
    faces, zero on regions. It is the Poisson count a message VARIANCE needs (`1/n`), where the mass — a sum
    of fractional per-fragment shares — would over-state it. Inert until the priority-#3 mature-measurement
    channel consumes it (docs/calibration/archive/boundary_spliced_channel_design.md §4.1), so this test guards the
    plumbing, not a behaviour.
    """
    payload, ra = make_synthetic_payload()
    fl = make_gdna_fl_pmf()
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    g = build_node_geometry(
        chain,
        CalibrationSubstrate.from_payload(payload, ra),
        BoundarySubstrate.from_payload(payload),
        ra,
        fl,
        fl,
    )
    is_reg = np.asarray(chain.kind) == REGION
    for m_name, n_name in (
        ("spliced_pos_left", "spliced_n_pos_left"),
        ("spliced_pos_right", "spliced_n_pos_right"),
        ("spliced_neg_left", "spliced_n_neg_left"),
        ("spliced_neg_right", "spliced_n_neg_right"),
    ):
        mass = np.asarray(getattr(g, m_name), float)
        cnt = np.asarray(getattr(g, n_name), float)
        assert cnt.shape == mass.shape
        # regions carry no spliced mass -> and no spliced count
        assert not np.any(cnt[is_reg]), f"{n_name} nonzero on a REGION node"
        # identical support: the count is nonzero exactly where the mass is (same motif/exon gate)
        assert np.array_equal(cnt > 0, mass > 0), f"{n_name} support != {m_name} support"
        # a count is an integer flux and bounds its own fractional mass from above
        assert np.all(cnt >= mass - 1e-9), (
            f"{n_name} < {m_name} (mass sums shares <= 1 per fragment)"
        )
        assert np.allclose(cnt, np.round(cnt)), f"{n_name} is not integral"


def test_zero_opportunity_nodes_carry_zero_mass():
    """THE GEOMETRY<->ACCUMULATOR CROSS-CHECK: a node face with no opportunity must carry no mass.

    `eff == 0` means "no fragment of any length in the FL pmf can be observed here". The accumulator can
    then have deposited nothing, so `mass == 0` and `mass/eff` is 0/0 — an empty node. If mass > 0 while
    eff == 0, the accumulator's deposit rule and the geometry's divisor DISAGREE, and the disagreement is
    SILENT: `build_node_geometry` floors every eff at `_EPS = 1e-9` (node_geometry.py:209-214), so a
    violation becomes a ~1e9 density rather than an exception.

    This is the test that would have caught the `region_eff_length` off-by-one. That bug computed
    `E[max(0, L-l)]` for a discrete `max(0, L-l+1)` start-position count, so a region exactly as long as the
    shortest fragment got eff = 0 while genuinely containing fragments. On ambig_dense_10mb one 100 bp
    region with 2 fragments read rho = 2e9 and consumed 84% of the NPMLE prior's fitted grid. Nothing tested
    the geometry AGAINST the accumulator — only against its own closed form — so it survived.

    NB the guard is on OPPORTUNITY (eff), not on COUNT. A zero count over a real eff is a genuine
    measurement ("my rate is low"); a zero count over zero eff is not a measurement at all.

    The FL spike is set to the fixture's own region length (100 bp) ON PURPOSE. That is the L == l_min
    boundary the off-by-one got wrong, and it is the ONLY place this invariant can bite: with the default
    50 bp FL the fixture's eff is 51 and the test passes on the broken code too (verified by mutation).
    """
    from rigel.calibration.node_geometry import _EPS

    payload, ra = make_synthetic_payload()
    # region_size_bp == 100 for every fixture region; put the shortest fragment exactly there.
    fl = make_gdna_fl_pmf(mean=100, max_size=200)
    assert int(np.asarray(ra.region_size_bp)[0]) == int(np.nonzero(fl)[0].min()), (
        "this test only bites when a region sits exactly at the shortest fragment length"
    )
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    sub = CalibrationSubstrate.from_payload(payload, ra)
    bsub = BoundarySubstrate.from_payload(payload)
    g = build_node_geometry(chain, sub, bsub, ra, fl, fl)

    at_floor = _EPS * 1.001
    for face in ("left", "right"):
        mass = np.asarray(getattr(g, f"mass_{face}"), float)
        for comp in ("gdna", "rna"):
            eff = np.asarray(getattr(g, f"eff_{comp}_{face}"), float)
            bad = (eff <= at_floor) & (mass > 0.0)
            assert not np.any(bad), (
                f"eff_{comp}_{face}: {int(bad.sum())} node(s) have mass > 0 with NO opportunity "
                f"(eff at the _EPS floor) -> mass/eff would read ~{(mass[bad] / _EPS).max():.3g}. "
                f"The accumulator deposited where the geometry says nothing can be observed."
            )
        eff_s = np.asarray(getattr(g, f"eff_spl_{face}"), float)
        spl = np.asarray(getattr(g, f"spliced_pos_{face}"), float) + np.asarray(
            getattr(g, f"spliced_neg_{face}"), float
        )
        bad = (eff_s <= at_floor) & (spl > 0.0)
        assert not np.any(bad), (
            f"eff_spl_{face}: spliced mass deposited with no spliced opportunity"
        )


def test_region_eff_length_zero_iff_shorter_than_every_fragment():
    """`eff == 0` must mean exactly one thing: the region is shorter than the SHORTEST fragment in the pmf.

    Pins the boundary case the off-by-one got wrong. L == l_min contains a fragment at exactly one start
    position, so eff == 1 — NOT 0. Only L < l_min is a true zero-opportunity node.
    """
    from rigel.calibration.effective_length import region_eff_length

    pmf = np.zeros(301)
    pmf[100:201] = 1.0  # uniform FL over [100, 200]; shortest = 100
    L = np.array([50.0, 99.0, 100.0, 101.0, 150.0, 1000.0])
    eff = region_eff_length(L, pmf)
    assert eff[0] == 0.0 and eff[1] == 0.0, "L < l_min must be zero opportunity"
    assert eff[2] > 0.0, "L == l_min contains a fragment at ONE position -> eff must be > 0, not 0"
    assert np.all(eff[2:] > 0.0)
    assert np.all(np.diff(eff) >= 0.0), "eff must be non-decreasing in L"
