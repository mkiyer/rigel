"""A reverse splice face preserves the component densities under either length gap."""

from typing import NamedTuple

import numpy as np
import pytest

from rigel.calibration import splice_graph as SG
from rigel.native import transfer_prepare
from _transfer_harness import NONE, Prepared, _rule

GAPS = [(250, 78), (78, 250), (150, 150)]
#: (TSS, TES, DONOR, ACCEPTOR) by strand; a DONOR bit marks an intron's genomic low end on either strand
BITS = {
    0: (SG.FLAG_TSS_POS, SG.FLAG_TES_POS, SG.FLAG_DONOR_POS, SG.FLAG_ACCEPTOR_POS),
    1: (SG.FLAG_TSS_NEG, SG.FLAG_TES_NEG, SG.FLAG_DONOR_NEG, SG.FLAG_ACCEPTOR_NEG),
}


@pytest.mark.parametrize("gdna_length,rna_length", GAPS)
@pytest.mark.parametrize("junction_rate", [0.0, 2.0])
def test_reverse_splice_recovers_the_boundary_density_mixture(
    gdna_length, rna_length, junction_rate
):
    """Expected counts are independent sums of density times each origin's admissible placements.

    The message is sharp and the boundary deep, so its mode must recover the known mixture.
    The zero-junction case also tests the pure change of count frame.
    """
    case = _density_case(gdna_length, rna_length, junction_rate)
    _check_map(case, 0, 1)


@pytest.mark.parametrize("gdna_length,rna_length", GAPS)
@pytest.mark.parametrize("source,destination", [(1, 2), (2, 1)])
def test_shared_density_changes_count_frame_at_the_intron_face(
    gdna_length, rna_length, source, destination
):
    _check_map(_density_case(gdna_length, rna_length, 0.0), source, destination)


@pytest.mark.parametrize("gdna_length,rna_length", GAPS)
@pytest.mark.parametrize("junction_rate", [0.0, 2.0])
@pytest.mark.parametrize("source,destination", [(0, 1), (1, 0), (2, 1), (1, 2)])
def test_alternative_splice_faces_preserve_each_component_opportunity(
    gdna_length, rna_length, junction_rate, source, destination
):
    case = _density_case(gdna_length, rna_length, junction_rate, kind="alternative")
    _check_map(case, source, destination)
    assert case.prepared.faces.at(source, destination).width == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("gdna_length,rna_length", GAPS)
@pytest.mark.parametrize("junction_rate", [0.0, 2.0])
@pytest.mark.parametrize("source,destination", [(0, 1), (1, 0)])
def test_terminus_outside_face_preserves_each_component_opportunity(
    gdna_length, rna_length, junction_rate, source, destination
):
    case = _density_case(gdna_length, rna_length, junction_rate, kind="terminus")
    _check_map(case, source, destination)


@pytest.mark.parametrize("kind", ["intron", "alternative", "terminus", "terminus_junction"])
@pytest.mark.parametrize("mirror", [False, True])
@pytest.mark.parametrize("strand", [0, 1])
@pytest.mark.parametrize("junction_opportunity", [0.4, 3.0])
def test_each_face_reads_its_own_junction_route_in_every_orientation(
    kind, mirror, strand, junction_opportunity
):
    """A junction's RNA rate is its count over ITS OWN opportunity — not the boundary's ``Er`` — read
    from the route on the junction's exon side: lo or hi, on either strand. A terminus boundary adds
    its own spliced crossing, and its junction's route when that leaves from the outside exon."""
    for gdna_length, rna_length in GAPS:
        case = _density_case(
            gdna_length,
            rna_length,
            2.0,
            kind=kind,
            mirror=mirror,
            strand=strand,
            junction_opportunity=junction_opportunity,
            crossing_rate=1.5,
        )
        far = [(case.far, 1), (1, case.far)] if kind == "alternative" else []
        for source, destination in [(case.near, 1), (1, case.near), *far]:
            _check_map(case, source, destination)


def test_a_zero_rna_opportunity_builds_no_splice_face_and_claims_no_flux():
    """An exon no RNA fragment fits in: no splice face, and its junction flux claims no RNA level."""
    case = _density_case(250, 78, 2.0, rna_opportunity_scale=(0.0, 1.0, 1.0))
    faces = case.prepared.faces
    assert faces.kind_at(1, 0) == NONE and faces.kind_at(0, 1) == NONE
    assert case.prepared.lanes["pos"].flux_at(0, 1) is None


def test_a_vanishing_rna_opportunity_tends_to_no_claim():
    """Down to a vanishing opportunity the splice face still delivers the truth, while the flux level
    flattens towards the zero opportunity's no-claim."""
    entropies = []
    for scale in (1e-3, 1e-6, 1e-9, 1e-12):
        case = _density_case(250, 78, 2.0, rna_opportunity_scale=(scale, 1.0, 1.0))
        if scale >= 1e-6:
            _check_map(case, 1, 0)
        row = case.prepared.lanes["pos"].flux_at(0, 1)
        w = np.exp(row - row.max())
        w /= w.sum()
        entropies.append(float(-(w[w > 0] * np.log(w[w > 0])).sum()))
    assert np.all(np.diff(entropies) > 0)
    assert np.log(case.lam.size) - entropies[-1] < 0.01


class _Case(NamedTuple):
    prepared: Prepared
    lam: np.ndarray
    truth_fg: np.ndarray
    near: int
    far: int


def _density_case(
    gdna_length,
    rna_length,
    junction_rate,
    kind="intron",
    mirror=False,
    strand=0,
    junction_opportunity=1.0,
    crossing_rate=0.0,
    rna_opportunity_scale=(1.0, 1.0, 1.0),
):
    """Three slots — region, boundary, region — with known gDNA, unspliced and spliced densities.

    ``near`` is the exon on the junction's side (left, or right when ``mirror``), ``far`` the other
    region: an intron (``kind="intron"``) or an exon. The junction's count is its rate times ITS OWN
    opportunity, ``junction_opportunity`` times the boundary's ``Er``; its route carries the rate.
    """
    near, far = (2, 0) if mirror else (0, 2)
    tss, tes, donor, acceptor = BITS[strand]
    body_right, body_left = (tss, tes) if strand == 0 else (tes, tss)
    lengths = np.zeros(3)
    lengths[near], lengths[far] = 400, 10000
    ag = lengths - gdna_length + 1.0
    ar = lengths - rna_length + 1.0
    ag[1], ar[1] = gdna_length - 1.0, rna_length - 1.0
    ar = ar * np.asarray(rna_opportunity_scale)
    is_exon = np.zeros(3, bool)
    is_exon[near] = True
    is_exon[far] = kind != "intron"
    # the region away from the near exon lies past the site, and the near exon outside the terminus
    splice_bit = acceptor if mirror else donor
    outside_bit = body_left if mirror else body_right
    flags = {
        "intron": splice_bit,
        "alternative": splice_bit,
        "terminus": outside_bit,
        "terminus_junction": outside_bit | splice_bit,
    }[kind]
    route_rate = 0.0 if kind == "terminus" else junction_rate
    crossing = {"terminus": junction_rate, "terminus_junction": crossing_rate}.get(kind, 0.0)
    rna = np.ones(3)
    rna[near] += route_rate + crossing
    depth = 10000.0
    count = (ag + rna * ar) * depth
    truth_fg = ag / (ag + rna * ar)
    zero = np.zeros((3, 2))
    sj, route, spliced = zero.copy(), zero.copy(), zero.copy()
    sj[1, strand] = route_rate * junction_opportunity * ar[1] * depth
    route[1, strand] = route_rate * depth
    spliced[1, strand] = crossing * ar[1] * depth
    fp = np.full(3, strand == 0)
    fn = ~fp
    sense = 0.5 * truth_fg + 0.99 * (1 - truth_fg)
    cnt = np.column_stack([count * sense, count * (1 - sense)])
    left, right = np.array([-1, 0, 1]), np.array([1, 2, -1])
    lam = np.linspace(-12, 12, 2401)
    tables = transfer_prepare(
        lam=lam,
        is_boundary=np.array([False, True, False]),
        is_exon=is_exon,
        free_pos=fp,
        free_neg=fn,
        exon_pos=is_exon & fp,
        exon_neg=is_exon & fn,
        left=left,
        right=right,
        flags=np.array([0, flags, 0], np.uint16),
        cnt=np.ascontiguousarray(cnt if strand == 0 else cnt[:, ::-1]),
        spliced=spliced,
        sj_count=sj,
        sj_count_lo=zero if mirror else sj,
        sj_count_hi=sj if mirror else zero,
        route_rate_lo=zero if mirror else route,
        route_rate_hi=route if mirror else zero,
        eff_gdna=ag,
        eff_rna=ar,
        has_own_composition=np.ones(3, bool),
        has_strand=True,
        kappa=0.99,
        rho_gdna=depth,
        rho_rna=depth,
        split_live=True,
    )
    return _Case(Prepared(tables, left, right, fp, fn), lam, truth_fg, near, far)


def _check_map(case, source, destination):
    lam, truth_fg = case.lam, case.truth_fg
    source_logodds = np.log(truth_fg[source] / (1 - truth_fg[source]))
    source_row = -0.5 * ((lam - source_logodds) / 0.1) ** 2
    delivered = _rule(case.prepared, source, destination, own=source_row)
    got_logodds = float(lam[np.argmax(delivered)])
    want_logodds = float(np.log(truth_fg[destination] / (1 - truth_fg[destination])))
    assert got_logodds == pytest.approx(want_logodds, abs=2 * (lam[1] - lam[0]))
