"""``capture_eff_length``: a transcript's effective length under capture is its conserved share of every
object it deposits on, each at that object's capture efficiency.

The geometry half: ``transcript_objects`` must map every transcript to every piece its exons overlap, on any
partition it is handed (the coarse partition built here by hand still has exon boundaries in a region's
interior, which the shipped partition no longer produces — losing that fixture would retire the guard); each
object's share must be what the deposit rule gives it, per object, through the reference accumulator; and a
transcript's shares must sum to its fl-marginal length exactly. The pricing half: every efficiency 1 (a field
with no reference) returns the FL-marginal lengths bit-identically; a piece prices at its region, a
contiguous cut at its boundary, a junction at its two splice sites' boundaries less the intron beside them;
a transcript of exons shorter than a fragment is read through its cuts; the ruler reads the result's
efficiencies and the geometry, never the counts.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
import pytest
from _index_builder import build_test_index

from rigel.calibration.capture_eff_length import (
    transcript_capture_eff_lengths,
    transcript_objects,
)
from rigel.calibration.effective_length import contained_eff_length
from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig


def _pmf(mean=200.0, sd=50.0, lo=100, hi=300) -> np.ndarray:
    w = np.arange(hi + 1, dtype=np.float64)
    p = np.exp(-0.5 * ((w - mean) / sd) ** 2)
    p[(w < lo) | (w > hi)] = 0.0
    return p / p.sum()


PMF = _pmf()
MU = float(np.dot(np.arange(PMF.shape[0]), PMF))


def _fl(lengths, pmf=PMF) -> np.ndarray:
    """The fl-marginal length ``Σ_w f(w)(L − w + 1)⁺`` per transcript."""
    w = np.arange(pmf.shape[0], dtype=np.float64)
    L = np.asarray(lengths, dtype=np.float64)[:, None]
    return (pmf[None, :] * np.maximum(L - w[None, :] + 1.0, 0.0)).sum(1)


def _flank_mean(ra: RegionArrays, c) -> np.ndarray:
    """A boundary's efficiency on a field where capture adds over bases and both flanks are longer than a
    fragment: its crossing fragments lie half on each side, so it reads its two flanks' mean."""
    lo, hi = boundary_region_indices(np.asarray(ra.ref_id))
    c = np.asarray(c, dtype=np.float64)
    return 0.5 * (c[lo] + c[hi])


def _cal(ra: RegionArrays, efficiency, reference: float | None, boundary=None) -> CalibrationResult:
    """THE fixture: a result whose capture efficiencies are stated outright — a region's, and a boundary's
    (its flanks' mean unless given) — with counts that carry nothing; the supports are the geometry the
    ruler reads to find an intron piece too short to contain a fragment. ``reference`` ``None`` is a field
    with no enriched mode, where every efficiency must be 1."""
    n = int(ra.n_regions)
    lo, _hi = boundary_region_indices(np.asarray(ra.ref_id))
    ne = lo.shape[0]
    z = np.zeros(n)
    ez = np.zeros(ne)
    c = np.asarray(efficiency, dtype=np.float64)
    return CalibrationResult(
        count_gdna_region=z.copy(),
        count_rna_region=z.copy(),
        count_gdna_boundary=ez.copy(),
        count_rna_boundary=ez.copy(),
        count_rna_spliced_boundary=ez.copy(),
        boundary_mass_per_crossing=np.ones(ne),
        count_rna_sj=np.zeros(0),
        boundary_spliced_mass_per_crossing=np.ones(ne),
        sj_mass_per_crossing=np.ones(0),
        gdna_region_eff_len=contained_eff_length(
            np.asarray(ra.region_size_bp, dtype=np.float64), PMF
        ),
        gdna_boundary_eff_len=np.full(ne, MU - 1.0),
        gdna_boundary_conserved_len=np.full(ne, MU - 1.0),
        rna_region_eff_len=np.asarray(ra.region_size_bp, dtype=np.float64),
        rna_boundary_eff_len=np.full(ne, MU - 1.0),
        gdna_frac_region=z.copy(),
        rna_pos_frac_region=z.copy(),
        rna_neg_frac_region=z.copy(),
        gdna_frac_boundary=ez.copy(),
        rna_pos_frac_boundary=ez.copy(),
        rna_neg_frac_boundary=ez.copy(),
        gdna_density_global=0.01,
        gdna_reference_density=reference,
        gdna_reference_members=0 if reference is None else 1,
        gdna_capture_efficiency_region=c,
        gdna_capture_efficiency_boundary=(
            _flank_mean(ra, c) if boundary is None else np.asarray(boundary, dtype=np.float64)
        ),
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        n_boundaries=ne,
        n_sj=0,
        config=CalibrationConfig(),
    )


def _exon_mask(ra: RegionArrays, a: int, b: int) -> np.ndarray:
    s, e = np.asarray(ra.start), np.asarray(ra.end)
    return (e > a) & (s < b)


def _tidx(idx, tid: str) -> int:
    tdf = idx.t_df
    return int(tdf.loc[tdf["t_id"] == tid, "t_index"].iloc[0])


# Two same-strand transcripts sharing gene g1. t1's first exon [150, 300) is a sub-interval of t0's
# [100, 300); both halves carry the identical EXON_POS signature, so the coarse partition MERGES them
# into a single region [100, 300) whose interior holds t1's exon start (150).
_MISALIGNED_GTF = """\
chr1\ttest\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t501\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t151\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\ttest\texon\t501\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
"""

# six 500 bp exons + a single-exon nascent parent covering the whole 1000..6500 span
_MULTIEXON_GTF = (
    "".join(
        f'chr1\ttest\texon\t{s + 1}\t{s + 500}\t.\t+\t.\tgene_id "gm"; transcript_id "mrna";\n'
        for s in range(1000, 6500, 1000)
    )
    + 'chr1\ttest\texon\t1001\t6500\t.\t+\t.\tgene_id "gm"; transcript_id "nasc";\n'
)

# ten exons of 40 bp at 1,040 bp pitch: every piece shorter than the shortest fragment (100 bp)
_TINY_GTF = "".join(
    f'chr1\ttest\texon\t{1001 + i * 1040}\t{1040 + i * 1040}\t.\t+\t.\tgene_id "gt"; transcript_id "tiny";\n'
    for i in range(10)
)

_TWO_REF_GTF = "".join(
    f'chrA\ttest\texon\t{s + 1}\t{s + 400}\t.\t+\t.\tgene_id "ga"; transcript_id "ta";\n'
    for s in (200, 800)
) + "".join(
    f'chrB\ttest\texon\t{s + 1}\t{s + 400}\t.\t+\t.\tgene_id "gb"; transcript_id "tb";\n'
    for s in (200, 800)
)


@pytest.fixture(scope="module")
def misaligned_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _MISALIGNED_GTF, genome_size=1000, name="misaligned")


@pytest.fixture(scope="module")
def multiexon_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _MULTIEXON_GTF, genome_size=7000, name="multiexon")


@pytest.fixture(scope="module")
def tiny_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _TINY_GTF, genome_size=14000, name="tiny")


@pytest.fixture(scope="module")
def two_ref_index(tmp_path_factory):
    return build_test_index(
        tmp_path_factory, _TWO_REF_GTF, name="tworef", refs={"chrA": 2000, "chrB": 2000}
    )


def _coarsened(idx) -> RegionArrays:
    """The index's regions with adjacent equal-signature neighbours MERGED — a deliberately coarse
    partition whose regions contain interior exon boundaries."""
    n = idx.regions_df
    ref = n["ref_name"].astype(str).to_numpy()
    sig = n["signature"].to_numpy(np.uint8)
    start = n["start"].to_numpy(np.int64)
    end = n["end"].to_numpy(np.int64)
    keep = np.r_[True, (ref[1:] != ref[:-1]) | (sig[1:] != sig[:-1]) | (start[1:] != end[:-1])]
    i = np.flatnonzero(keep)
    last = np.r_[i[1:], ref.size] - 1
    merged = pd.DataFrame(
        {"ref_name": ref[i], "start": start[i], "end": end[last], "signature": sig[i]}
    )
    assert len(merged) < len(n), "the fixture must actually coarsen, or it tests nothing"
    return RegionArrays.from_frame(merged, idx.ref_name_to_id)


# a transcript cut into 1–48 bp pieces by two others' exon ends, with two junctions, for the per-object gate
_PIECES_GTF = (
    'chr1\ttest\texon\t101\t160\t.\t+\t.\tgene_id "gp"; transcript_id "t0";\n'
    'chr1\ttest\texon\t301\t340\t.\t+\t.\tgene_id "gp"; transcript_id "t0";\n'
    'chr1\ttest\texon\t501\t560\t.\t+\t.\tgene_id "gp"; transcript_id "t0";\n'
    'chr1\ttest\texon\t102\t150\t.\t+\t.\tgene_id "gp"; transcript_id "t1";\n'
    'chr1\ttest\texon\t306\t338\t.\t+\t.\tgene_id "gp"; transcript_id "t1";\n'
    'chr1\ttest\texon\t104\t106\t.\t+\t.\tgene_id "gp"; transcript_id "t2";\n'
    'chr1\ttest\texon\t511\t512\t.\t+\t.\tgene_id "gp"; transcript_id "t2";\n'
)

# two 500 bp exons around a 50 bp intron — shorter than the shortest fragment (100 bp)
_SHORT_INTRON_GTF = "".join(
    f'chr1\ttest\texon\t{s + 1}\t{s + 500}\t.\t+\t.\tgene_id "gs"; transcript_id "si";\n'
    for s in (1000, 1550)
)


# two spliced transcripts whose introns hold another transcript's exon, so each intron is three regions: ``tj``'s
# intron pieces beside its junction are 50 bp (shorter than every fragment), ``tk``'s are 400 bp
_MULTI_REGION_INTRON_GTF = (
    'chr1\ttest\texon\t1001\t1500\t.\t+\t.\tgene_id "gj"; transcript_id "tj";\n'
    'chr1\ttest\texon\t3001\t3500\t.\t+\t.\tgene_id "gj"; transcript_id "tj";\n'
    'chr1\ttest\texon\t1551\t2950\t.\t+\t.\tgene_id "gj"; transcript_id "tx";\n'
    'chr1\ttest\texon\t5001\t5500\t.\t+\t.\tgene_id "gk"; transcript_id "tk";\n'
    'chr1\ttest\texon\t7001\t7500\t.\t+\t.\tgene_id "gk"; transcript_id "tk";\n'
    'chr1\ttest\texon\t5901\t6600\t.\t+\t.\tgene_id "gk"; transcript_id "ty";\n'
)


@pytest.fixture(scope="module")
def multi_region_intron_index(tmp_path_factory):
    return build_test_index(
        tmp_path_factory, _MULTI_REGION_INTRON_GTF, genome_size=9000, name="multiintron"
    )


@pytest.fixture(scope="module")
def pieces_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _PIECES_GTF, genome_size=1000, name="pieces")


@pytest.fixture(scope="module")
def short_intron_index(tmp_path_factory):
    return build_test_index(
        tmp_path_factory, _SHORT_INTRON_GTF, genome_size=3000, name="shortintron"
    )


# --- the geometry: every transcript, every object, and the shares partition the start count ----------


@pytest.mark.parametrize("partition", ["coarse", "live"])
def test_every_transcript_maps_to_the_pieces_its_exons_overlap(misaligned_index, partition):
    """No transcript (mature or nascent entity) is dropped, and the piece covering the exon start that
    sits INTERIOR to a merged region is mapped — the off-by-one that skipped it on the coarse partition
    and the live partition, where that start is a region interface, both."""
    idx = misaligned_index
    ra = _coarsened(idx) if partition == "coarse" else RegionArrays.from_index(idx)
    obj = transcript_objects(idx, ra, PMF)
    assert set(int(x) for x in obj.t) == set(range(len(idx.t_df))), "a transcript was dropped"
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    covering = int(np.flatnonzero((starts <= 150) & (ends > 150))[0])
    tdf = idx.t_df
    for ti in range(len(tdf)):
        a, b = int(tdf["start"].iloc[ti]), int(tdf["end"].iloc[ti])
        if a <= 150 < b:
            assert covering in set(int(p) for p in obj.piece[obj.t == ti])


@pytest.mark.parametrize("partition", ["coarse", "live"])
def test_a_transcripts_shares_sum_to_its_fl_marginal_length(misaligned_index, partition):
    """Σ contained + Σ cut shares over a transcript's objects is Σ_w f(w)(L − w + 1)⁺ — on the coarse
    partition, whose pieces spill beyond the exons and must be clipped to them, and on the live one."""
    idx = misaligned_index
    ra = _coarsened(idx) if partition == "coarse" else RegionArrays.from_index(idx)
    obj = transcript_objects(idx, ra, PMF)
    total = np.zeros(len(idx.t_df))
    np.add.at(total, obj.t, obj.contained)
    np.add.at(total, obj.t[obj.cut_row], obj.cut_share)
    np.testing.assert_allclose(total, _fl(idx.t_df["length"].to_numpy()), rtol=1e-9)


def test_the_pieces_are_on_the_transcripts_own_reference(two_ref_index):
    """A second reference numbers its regions after the first's: every piece of a transcript must lie on
    the transcript's reference, or a region index would silently read another chromosome's efficiency."""
    idx = two_ref_index
    ra = RegionArrays.from_index(idx)
    obj = transcript_objects(idx, ra, PMF)
    ref_of_t = idx.t_df["ref"].astype(str).map(idx.ref_name_to_id).to_numpy()
    np.testing.assert_array_equal(np.asarray(ra.ref_id)[obj.piece], ref_of_t[obj.t])


def test_each_share_is_what_the_deposit_rule_gives_the_object(pieces_index):
    """THE GATE, per object, on the index path. Every placement of ``t0`` — cut into 1–48 bp pieces by its
    siblings' exon ends, spliced twice — through the reference accumulator on the index's own partition:
    each piece's contained units are its contained share, each boundary's and each junction's mass its cut
    share, to 1e-12. PERTURBATION: pricing a crossing start once per cut it crosses (the incidence count)
    overstates every cut beside a piece shorter than a fragment."""
    from native._accumulator_reference import Accumulator, DepositOutcome, Partition

    from rigel.types import Strand

    idx = pieces_index
    ra = RegionArrays.from_index(idx)
    w = np.arange(41, dtype=np.float64)
    pmf = np.where(w >= 2.0, 1.0 + 0.5 * np.sin(w), 0.0)
    pmf /= pmf.sum()
    t0 = _tidx(idx, "t0")
    obj = transcript_objects(idx, ra, pmf)
    rows = np.flatnonzero(obj.t == t0)
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    blocks = [(100, 160), (300, 340), (500, 560)]
    bounds = np.r_[starts, ends[-1]]
    sj = [(0, blocks[k][1], blocks[k + 1][0], int(Strand.POS)) for k in range(2)]
    part = Partition.from_region_bounds([bounds], sj=sj)
    offs = np.cumsum([0, 60, 40, 60])
    region = np.zeros(part.n_regions)
    boundary = np.zeros(part.n_boundaries)
    junction = np.zeros(part.n_sj)
    for width in np.flatnonzero(pmf):
        acc = Accumulator(part, max_fragment_length=60)
        for s in range(0, 160 - int(width) + 1):
            k0 = int(np.searchsorted(offs, s, side="right")) - 1
            k1 = int(np.searchsorted(offs, s + width - 1, side="right")) - 1
            introns = tuple((blocks[k][1], blocks[k + 1][0]) for k in range(k0, k1))
            out = acc.deposit(
                0,
                blocks[k0][0] + s - int(offs[k0]),
                blocks[k1][0] + s + int(width) - int(offs[k1]),
                observed_introns=introns,
                align_strand=Strand.POS,
                sj_strand=Strand.POS if introns else Strand.NONE,
            )
            assert out is DepositOutcome.DEPOSITED
        region += pmf[width] * acc.tally.region_contained_count.sum(1)
        boundary += pmf[width] * (
            acc.tally.boundary_unspliced_mass + acc.tally.boundary_spliced_mass
        )
        junction += pmf[width] * acc.tally.sj_mass.sum(1)
    assert np.min(ends[obj.piece[rows]] - starts[obj.piece[rows]]) <= 2, (
        "the fixture must cut tiny pieces"
    )
    np.testing.assert_allclose(region[obj.piece[rows]], obj.contained[rows], rtol=0, atol=1e-12)
    cuts = np.flatnonzero(obj.t[obj.cut_row] == t0)
    left = obj.piece[obj.cut_row[cuts]]
    j = obj.is_junction[cuts]
    assert int(j.sum()) == 2
    np.testing.assert_allclose(boundary[left[~j]], obj.cut_share[cuts[~j]], rtol=0, atol=1e-12)
    np.testing.assert_allclose(junction, obj.cut_share[cuts[j]], rtol=0, atol=1e-12)


# --- the pricing: efficiencies 1 return fl; each object at its own efficiency; the junction ---------


def test_no_reference_returns_the_fl_marginal_lengths_bit_identically(multiexon_index):
    """Capture off, or no gDNA: every efficiency is 1 and the ruler is ``fl`` EXACTLY, not within float
    noise — a contraction from rounding is a systematic bias, not a tolerance."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    fl = np.linspace(900.0, 2000.0, len(idx.t_df))
    eff = transcript_capture_eff_lengths(_cal(ra, np.ones(ra.n_regions), None), ra, idx, fl, PMF)
    np.testing.assert_array_equal(eff, fl)


def test_efficiencies_of_one_everywhere_under_a_reference_contract_nothing(multiexon_index):
    """With a reference and every object at it, the factor is 1 to floating point on every transcript,
    spliced and unspliced alike: a junction between two exons at 1, its introns at 1, prices at 1."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, np.ones(ra.n_regions), 1.0), ra, idx, fl, PMF)
    np.testing.assert_allclose(eff, fl, rtol=1e-12)


def test_a_field_contracts_and_never_expands(multiexon_index):
    """On a field where capture adds over bases — every boundary at its flanks' mean — a junction prices
    at its two exons' mean and no transcript reads longer than its fl-marginal length."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500)] = 1.0  # the first exon captured, everything else depleted
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, PMF)
    assert np.all(eff <= fl + 1e-9)
    assert np.any(eff < fl - 1e-6)


def test_the_length_is_the_share_weighted_mean_of_the_objects_efficiencies(multiexon_index):
    """The claim, computed independently: the six-exon mRNA's 500 bp exons each hold ``501 − E[w]`` of
    contained share and each junction between them ``E[w − 1]`` (no fragment reaches a second cut); with
    exons 1–3 at 1, 4–6 at 0.2 and the introns at 0.001, each junction prices at its two exons' mean."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    level = [1.0, 1.0, 1.0, 0.2, 0.2, 0.2]
    for k, s in enumerate(range(1000, 6500, 1000)):
        c[_exon_mask(ra, s, s + 500)] = level[k]
    m = _tidx(idx, "mrna")
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, PMF)
    contained, cut = 501.0 - MU, MU - 1.0
    junctions = [0.5 * (a + b) for a, b in zip(level[:-1], level[1:])]
    expected = fl[m] * (contained * sum(level) + cut * sum(junctions)) / (6 * contained + 5 * cut)
    assert eff[m] == pytest.approx(expected, rel=1e-9)


def test_a_junction_adds_its_two_sides_it_does_not_average_them(multiexon_index):
    """THE FALSIFICATION TEST for the junction price. A fragment across a junction between two captured
    exons is all captured exon; the gDNA fragments across its low and high boundaries are half exon,
    half intron. With every exon at 1, every intron at 0.001 and every boundary at its flanks' mean
    (0.5005), the mRNA's junctions price at ``0.5005 + 0.5005 − 0.001 = 1`` and the mRNA reads its full
    length. PERTURBATION: pricing a junction at the mean of its two boundaries (0.5005) reads it near half —
    the 1.68× under-price the gDNA-only prices measured. (Its adjacent pieces' mean is right on this uniform
    field and wrong where a piece beside the junction holds no fragment: the tiny-exon and short-intron
    tests below catch that one.)"""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    for s in range(1000, 6500, 1000):
        c[_exon_mask(ra, s, s + 500)] = 1.0
    m = _tidx(idx, "mrna")
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, PMF)
    assert eff[m] == pytest.approx(fl[m], rel=1e-12)


def test_an_intron_piece_too_short_to_contain_a_fragment_reads_its_far_boundary(short_intron_index):
    """A 50 bp intron holds no gDNA fragment against a 100–300 bp pmf, so its own efficiency is the
    population's and says nothing about it: the junction reads the boundary on the intron's far side
    instead — both splice sites' boundaries here, so the junction is their mean — and the intron's own
    efficiency moves nothing. PERTURBATION: reading the intron region's own efficiency makes the length
    follow it."""
    idx = short_intron_index
    ra = RegionArrays.from_index(idx)
    si = _tidx(idx, "si")
    intron = _exon_mask(ra, 1500, 1550)
    assert intron.sum() == 1
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500) | _exon_mask(ra, 1550, 2050)] = 1.0
    cb = _flank_mean(ra, c)
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0, boundary=cb), ra, idx, fl, PMF)
    c2 = c.copy()
    c2[intron] = 0.9
    eff2 = transcript_capture_eff_lengths(_cal(ra, c2, 1.0, boundary=cb), ra, idx, fl, PMF)
    assert eff2[si] == eff[si]
    contained, cut = 501.0 - MU, MU - 1.0
    expected = fl[si] * (2 * contained + cut * 0.5005) / (2 * contained + cut)
    assert eff[si] == pytest.approx(expected, rel=1e-9)


def test_a_junction_reads_the_objects_beside_it_on_a_multi_region_intron(multi_region_intron_index):
    """Every region and every boundary at its own efficiency, so no two objects agree and a junction read
    through a wrong index cannot land on the right number. Each transcript is two 500 bp exons, so each holds
    ``501 − E[w]`` of contained share per exon and ``E[w − 1]`` at its one junction; the junction's objects
    are found by COORDINATE here, never by the module's index arithmetic: its low boundary is where the exon
    below ends, its high boundary where the exon above starts, and the intron pieces beside them are the
    regions touching those two positions — read at their own efficiency when they hold a fragment (``tk``,
    400 bp), at the boundary on their far side when they cannot (``tj``, 50 bp). PERTURBATION: reading the
    exons' other edges, the near side of a short intron piece, or one intron piece for both sides each moves
    the number."""
    idx = multi_region_intron_index
    ra = RegionArrays.from_index(idx)
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    lo, _hi = boundary_region_indices(np.asarray(ra.ref_id))
    n, ne = int(ra.n_regions), int(lo.shape[0])
    c = 0.05 + 0.9 * ((np.arange(n) * 0.618034) % 1.0)
    cb = 0.3 + 0.4 * ((np.arange(ne) * 0.414214) % 1.0)

    def boundary_at(pos):
        return int(np.flatnonzero(ends[lo] == pos)[0])

    def region_from(pos):
        return int(np.flatnonzero(starts == pos)[0])

    def region_to(pos):
        return int(np.flatnonzero(ends == pos)[0])

    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0, boundary=cb), ra, idx, fl, PMF)
    contained, cut = 501.0 - MU, MU - 1.0
    for tid, (e1, e2), (intron_lo, intron_hi) in (
        ("tk", (5000, 7000), (c[region_from(5500)], c[region_to(7000)])),
        ("tj", (1000, 3000), (cb[boundary_at(1550)], cb[boundary_at(2950)])),
    ):
        t = _tidx(idx, tid)
        assert ra.region_size_bp[region_from(e1 + 500)] in (50.0, 400.0)
        c_junction = cb[boundary_at(e1 + 500)] + cb[boundary_at(e2)] - 0.5 * (intron_lo + intron_hi)
        assert c_junction > 0.0
        c_exons = c[region_from(e1)] + c[region_from(e2)]
        expected = fl[t] * (contained * c_exons + cut * c_junction) / (2 * contained + cut)
        assert eff[t] == pytest.approx(expected, rel=1e-12), tid


def test_a_zero_length_fragment_places_nowhere(multiexon_index):
    """A length model with mass at ``w = 0`` (a smoothed real library's once did) deposits no fragment
    there: the lengths equal those of the same pmf with that mass removed, for every transcript."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500)] = 1.0
    fl = _fl(idx.t_df["length"].to_numpy())
    with_zero = PMF.copy()
    with_zero[0] = 0.05
    np.testing.assert_allclose(
        transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, with_zero),
        transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, PMF),
        rtol=1e-12,
    )


def test_a_transcript_of_exons_shorter_than_a_fragment_is_read_through_its_cuts(tiny_index):
    """Ten 40 bp exons against a 100–300 bp pmf: no piece can contain a fragment, so the transcript's whole
    length sits on its junctions. Its exons' own efficiencies carry no weight — they hold no share — and
    with every junction priced at ½ (splice-site boundaries at ½, introns at ½) the factor is exactly ½.
    PERTURBATION: an object set that drops the cuts reads this transcript as nothing."""
    idx = tiny_index
    ra = RegionArrays.from_index(idx)
    t = _tidx(idx, "tiny")
    fl = _fl(idx.t_df["length"].to_numpy())
    exons = np.zeros(ra.n_regions, dtype=bool)
    for i in range(10):
        exons |= _exon_mask(ra, 1000 + i * 1040, 1040 + i * 1040)
    half = np.full(ra.n_regions, 0.5)
    ne = int(boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0])
    eff = transcript_capture_eff_lengths(
        _cal(ra, half, 1.0, boundary=np.full(ne, 0.5)), ra, idx, fl, PMF
    )
    assert eff[t] == pytest.approx(0.5 * fl[t], rel=1e-12)
    loud = half.copy()
    loud[exons] = 0.9
    eff2 = transcript_capture_eff_lengths(
        _cal(ra, loud, 1.0, boundary=np.full(ne, 0.5)), ra, idx, fl, PMF
    )
    assert eff2[t] == eff[t]


def test_the_ruler_reads_the_efficiencies_and_nothing_from_the_counts(multiexon_index):
    """Two results with the same efficiencies and wildly different counts give the same lengths: the
    evidence was weighed upstream (`capture_efficiency`), and the ruler is geometry."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.3)
    c[_exon_mask(ra, 1000, 1500)] = 1.0
    fl = _fl(idx.t_df["length"].to_numpy())
    cal = _cal(ra, c, 1.0)
    loud = dataclasses.replace(
        cal,
        count_gdna_region=np.full(ra.n_regions, 1e6),
        count_gdna_boundary=np.full(cal.n_boundaries, 1e6),
        count_rna_region=np.full(ra.n_regions, 1e6),
    )
    np.testing.assert_array_equal(
        transcript_capture_eff_lengths(cal, ra, idx, fl, PMF),
        transcript_capture_eff_lengths(loud, ra, idx, fl, PMF),
    )


def test_no_nascent_mature_inversion_under_capture(multiexon_index):
    """On a field where capture adds over bases, a nascent parent holds its spliced child's captured exon
    and more depleted bases besides, so ``eff(nascent) ≥ eff(mature)``; and the mature genuinely
    contracts here."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500)] = 1.0
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, PMF)
    mrna, nasc = _tidx(idx, "mrna"), _tidx(idx, "nasc")
    assert eff[nasc] >= eff[mrna] - 1e-9
    assert eff[mrna] < fl[mrna] - 1e-6
