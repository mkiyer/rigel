"""``capture_eff_length``: a transcript's effective length under capture is its conserved share of every object it
deposits on, each at that object's capture efficiency.

The geometry half: ``transcript_objects`` maps every transcript to every piece its exons overlap on any partition it is
handed (the hand-built coarse partition keeps exon boundaries inside a region, which the shipped partition does not
produce); each object's share is what the deposit rule gives it, through the reference accumulator; and a transcript's
shares sum to its fl-marginal length exactly. The pricing half: with every efficiency 1 (no landscape, or one level) the
lengths are the fl-marginal ones bit-identically; a piece prices at its region, a contiguous cut at its boundary, a
junction at its two splice sites' boundaries less the intron beside them, bounded by the larger of the junction cap and
the best object it touches; a placement that leaves a cut's exons across a junction at that junction's crossing level
where it is higher; the cap is the saturated level, raised as far as each junction's spliced fragments support
against the unspliced RNA of its own exons; one factor on every efficiency is one factor on every length; a transcript of exons
shorter than a fragment is read through its cuts; and the ruler reads the result's efficiencies and the geometry, never
the counts.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
import pytest
from _index_builder import build_test_index
from _synthetic import contained_support, normal_pmf

from rigel.calibration.capture_eff_length import (
    cut_efficiencies,
    fit_junction_cap,
    TranscriptFrame,
    object_shares,
    transcript_pieces,
    rna_captured_fraction,
    transcript_capture_eff_lengths,
)
from rigel.calibration.effective_length import contained_eff_length
from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.result import CalibrationResult
from rigel.calibration.shared_length_law import LengthTables
from rigel.calibration.signature import RegionType, coarse_type_array


def transcript_objects(index, region_arrays, pmf):
    """Every transcript's pieces and cuts with its conserved share of each, on one pmf."""
    return object_shares(transcript_pieces(index, region_arrays), pmf)


PMF = normal_pmf()
MU = float(np.dot(np.arange(PMF.shape[0]), PMF))


def _one_law(pmf) -> LengthTables:
    """The EM's tables when both origins and the uncaptured law are one pmf: the efficiencies read as given."""
    return LengthTables(rna=pmf, gdna=pmf, anchor=pmf, decision=None)


def _lengths(cal, ra, idx, fl, tables, junction_cap: float | None = None) -> np.ndarray:
    """The ruler at a stated junction cap: by default the result's saturated level, the cap of a library whose spliced
    counts raise nothing (``fit_junction_cap``'s starting point)."""
    cap = cal.gdna_capture_saturation if junction_cap is None else junction_cap
    return transcript_capture_eff_lengths(cal, TranscriptFrame(idx, ra, tables), fl, cap)


def _flank_mean(ra: RegionArrays, c) -> np.ndarray:
    """A boundary's efficiency on a field where capture adds over bases and both flanks are longer than a fragment: its
    crossing fragments lie half on each side, so it reads its two flanks' mean."""
    lo, hi = boundary_region_indices(np.asarray(ra.ref_id))
    c = np.asarray(c, dtype=np.float64)
    return 0.5 * (c[lo] + c[hi])


def _cal(
    ra: RegionArrays,
    efficiency,
    bulk: float | None,
    boundary=None,
    saturation: float = 1.0,
    rna_region=None,
    sj_mass=None,
) -> CalibrationResult:
    """A result whose capture efficiencies are stated outright (a region's, and a boundary's, its flanks' mean unless
    given), with counts that carry nothing unless given (``rna_region``, each region's unspliced RNA; ``sj_mass``, each
    junction's conserved spliced mass, on the sj axis); the supports are the geometry the ruler reads to find an intron
    piece too short to contain a fragment. ``bulk`` ``None`` is a library with no landscape, where every efficiency
    must be 1. ``saturation`` is the landscape's saturated level in the efficiencies' unit; at its default 1, with every
    efficiency at most 1, the unit is the fully captured level and the junction's bound is 1."""
    n = int(ra.n_regions)
    lo, _hi = boundary_region_indices(np.asarray(ra.ref_id))
    ne = lo.shape[0]
    z = np.zeros(n)
    ez = np.zeros(ne)
    c = np.asarray(efficiency, dtype=np.float64)
    sj = np.zeros(0) if sj_mass is None else np.asarray(sj_mass, dtype=np.float64)
    return CalibrationResult(
        count_gdna_region=z.copy(),
        count_rna_region=z.copy() if rna_region is None else np.asarray(rna_region, np.float64),
        count_gdna_boundary=ez.copy(),
        count_rna_boundary=ez.copy(),
        count_rna_spliced_boundary=ez.copy(),
        boundary_mass_per_crossing=np.ones(ne),
        count_rna_sj=sj,
        boundary_spliced_mass_per_crossing=np.ones(ne),
        sj_mass_per_crossing=np.ones(sj.shape[0]),
        # not read here: the ruler and the cap take each junction's regions from the index
        sj_src_region=np.zeros(sj.shape[0], dtype=np.int64),
        sj_dst_region=np.zeros(sj.shape[0], dtype=np.int64),
        gdna_region_eff_len=contained_eff_length(
            np.asarray(ra.region_size_bp, dtype=np.float64), PMF
        ),
        gdna_boundary_eff_len=np.full(ne, MU - 1.0),
        gdna_boundary_conserved_len=np.full(ne, MU - 1.0),
        gdna_frac_region=z.copy(),
        rna_pos_frac_region=z.copy(),
        rna_neg_frac_region=z.copy(),
        gdna_frac_boundary=ez.copy(),
        rna_pos_frac_boundary=ez.copy(),
        rna_neg_frac_boundary=ez.copy(),
        gdna_density_global=0.01,
        gdna_bulk_density=bulk,
        gdna_capture_weight=1.0,
        gdna_capture_saturation=1.0 if bulk is None else saturation,
        gdna_capture_efficiency_region=c,
        gdna_capture_efficiency_boundary=(
            _flank_mean(ra, c) if boundary is None else np.asarray(boundary, dtype=np.float64)
        ),
        rna_sense_frac=0.9,
        n_regions=n,
        n_boundaries=ne,
        n_sj=sj.shape[0],
    )


def _exon_mask(ra: RegionArrays, a: int, b: int) -> np.ndarray:
    s, e = np.asarray(ra.start), np.asarray(ra.end)
    return (e > a) & (s < b)


def _tidx(idx, tid: str) -> int:
    tdf = idx.t_df
    return int(tdf.loc[tdf["t_id"] == tid, "t_index"].iloc[0])


# Two same-strand transcripts sharing gene g1. t1's first exon [150, 300) is a sub-interval of t0's [100, 300); both
# halves carry the identical EXON_POS signature, so the coarse partition merges them into a single region [100, 300)
# whose interior holds t1's exon start (150).
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
    """The index's regions with adjacent equal-signature neighbours merged: a coarse partition whose regions contain
    interior exon boundaries."""
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

# two 500 bp exons around a 50 bp intron, shorter than the shortest fragment (100 bp)
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
    """No transcript (mature or nascent entity) is dropped, and the piece covering an exon start interior to a merged
    region is mapped, on the coarse partition and on the live one, where that start is a region interface."""
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
    np.testing.assert_allclose(
        total, contained_support(idx.t_df["length"].to_numpy(), PMF), rtol=1e-9
    )


def test_the_pieces_are_on_the_transcripts_own_reference(two_ref_index):
    """A second reference numbers its regions after the first's: every piece of a transcript must lie on the
    transcript's reference, or a region index would read another chromosome's efficiency."""
    idx = two_ref_index
    ra = RegionArrays.from_index(idx)
    obj = transcript_objects(idx, ra, PMF)
    ref_of_t = idx.t_df["ref"].astype(str).map(idx.ref_name_to_id).to_numpy()
    np.testing.assert_array_equal(np.asarray(ra.ref_id)[obj.piece], ref_of_t[obj.t])


def test_each_share_is_what_the_deposit_rule_gives_the_object(pieces_index):
    """Per object, on the index path: every placement of ``t0`` (cut into 1–48 bp pieces by its siblings' exon ends,
    spliced twice) through the reference accumulator on the index's own partition. Each piece's contained units are its
    contained share, and each boundary's and junction's mass its cut share, to 1e-12. Pricing a crossing start once per
    cut it crosses (the incidence count) would overstate every cut beside a piece shorter than a fragment.
    """
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


def test_no_landscape_returns_the_fl_marginal_lengths_bit_identically(multiexon_index):
    """No landscape: every efficiency is 1 and the ruler is ``fl`` exactly, not within float noise, since a contraction
    from rounding would be a systematic bias."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    fl = np.linspace(900.0, 2000.0, len(idx.t_df))
    eff = _lengths(_cal(ra, np.ones(ra.n_regions), None), ra, idx, fl, _one_law(PMF))
    np.testing.assert_array_equal(eff, fl)


@pytest.mark.parametrize("tables", ["one law", "gDNA's own law apart", "every table apart"])
def test_efficiencies_of_one_everywhere_contract_nothing_bit_identically(multiexon_index, tables):
    """A landscape of one level reads every object at the bulk: the factor is exactly 1 on every transcript, spliced and
    unspliced (a junction between two exons at 1, its introns at 1, prices at exactly 1), so a library reading no
    capture has lengths ``fl`` bit for bit, whatever the tables. An efficiency of 1 carries nothing to the scorer's gDNA
    table, though its law differs from gDNA's uncaptured one and from RNA's."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    g = normal_pmf(mean=160.0, sd=40.0, lo=60, hi=300)
    t_g = normal_pmf(mean=230.0, sd=30.0, lo=150, hi=320)
    length_tables = {
        "one law": _one_law(PMF),
        "gDNA's own law apart": LengthTables(PMF, PMF, g, None),
        "every table apart": LengthTables(PMF, t_g, g, None),
    }[tables]
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(_cal(ra, np.ones(ra.n_regions), 1.0), ra, idx, fl, length_tables)
    np.testing.assert_array_equal(eff, fl)


def test_a_field_contracts_and_never_expands(multiexon_index):
    """In the unit of the fully captured level (saturation 1, no object above it), on a field where capture adds over
    bases (every boundary at its flanks' mean), a junction prices at its two exons' mean and no transcript reads longer
    than its fl-marginal length."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500)] = 1.0  # the first exon captured, everything else depleted
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(_cal(ra, c, 1.0), ra, idx, fl, _one_law(PMF))
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
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(_cal(ra, c, 1.0), ra, idx, fl, _one_law(PMF))
    contained, cut = 501.0 - MU, MU - 1.0
    junctions = [0.5 * (a + b) for a, b in zip(level[:-1], level[1:])]
    expected = fl[m] * (contained * sum(level) + cut * sum(junctions)) / (6 * contained + 5 * cut)
    assert eff[m] == pytest.approx(expected, rel=1e-9)


def test_a_junction_adds_its_two_sides_it_does_not_average_them(multiexon_index):
    """A fragment across a junction between two captured exons is all captured exon, while the gDNA fragments across its
    low and high boundaries are half exon, half intron. With every exon at 1, every intron at 0.001 and every boundary
    at its flanks' mean (0.5005), the mRNA's junctions price at ``0.5005 + 0.5005 − 0.001 = 1`` and the mRNA reads its
    full length. Pricing a junction at the mean of its two boundaries (0.5005) would read it near half. (The mean of its
    adjacent pieces is right on this uniform field but wrong where a piece beside the junction holds no fragment; the
    tiny-exon and short-intron tests below cover that.)
    """
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    for s in range(1000, 6500, 1000):
        c[_exon_mask(ra, s, s + 500)] = 1.0
    m = _tidx(idx, "mrna")
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(_cal(ra, c, 1.0), ra, idx, fl, _one_law(PMF))
    assert eff[m] == pytest.approx(fl[m], rel=1e-12)


def test_a_junction_is_captured_no_more_than_a_fully_captured_piece(multiexon_index):
    """Capture saturates: a fragment binds its best single probe part, so a fragment across a junction between two
    captured exons is captured no more than one wholly inside either. Where probes tile the exons to their edges, a gDNA
    fragment crossing an exon's boundary still binds most of a probe, so the boundary reads most of the exon's level
    rather than half, and the plain sum ``c_lo + c_hi − ½(c_intron,lo + c_intron,hi)`` would price the junction well
    above 1. With every exon at 1, every intron at 0.001 and every boundary at 0.835, the mRNA's junctions price at 1,
    it reads its full fl-marginal length, and no transcript reads longer than its own. The uncapped sum would read the
    mRNA long; a cap at the boundaries' own level would read it short.
    """
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    for s in range(1000, 6500, 1000):
        c[_exon_mask(ra, s, s + 500)] = 1.0
    ne = int(boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0])
    m = _tidx(idx, "mrna")
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(_cal(ra, c, 1.0, boundary=np.full(ne, 0.835)), ra, idx, fl, _one_law(PMF))
    assert eff[m] == pytest.approx(fl[m], rel=1e-12)
    assert np.all(eff <= fl * (1.0 + 1e-12))


def test_a_junction_the_introns_out_capture_prices_at_zero_never_below(multiexon_index):
    """Conservation takes the intron pieces' capture off the two boundaries', and where the introns are captured and
    the exons' edges are not it would go negative: every exon and boundary at 1 and every intron at 10 gives
    ``1 + 1 − 10``. A yield is not negative, so the junction prices at 0 and the mRNA reads its exons alone.
    """
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 10.0)
    for s in range(1000, 6500, 1000):
        c[_exon_mask(ra, s, s + 500)] = 1.0
    ne = int(boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0])
    m = _tidx(idx, "mrna")
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(
        _cal(ra, c, 1e-3, boundary=np.ones(ne), saturation=10.0), ra, idx, fl, _one_law(PMF)
    )
    assert eff[m] == pytest.approx(_mrna_length(fl[m], [1.0] * 6, [0.0] * 5), rel=1e-12)


def _mrna_length(fl_m: float, exon_levels, junction_prices) -> float:
    """The six-exon mRNA's length from its objects: each 500 bp exon holds ``501 − E[w]`` of contained share and
    each junction ``E[w − 1]`` (no fragment reaches a second cut)."""
    contained, cut = 501.0 - MU, MU - 1.0
    return (
        fl_m
        * (contained * sum(exon_levels) + cut * sum(junction_prices))
        / (len(exon_levels) * contained + len(junction_prices) * cut)
    )


def test_a_junction_may_out_capture_every_object_it_touches_up_to_saturation(multiexon_index):
    """In the bulk unit a probe spanning a junction binds a spliced RNA fragment as one part, where gDNA sees two, so
    each exon side reads part of the probe's level and the junction more than either. Saturation 100, every exon and
    boundary at 40, the introns at the bulk: the junction prices at ``40 + 40 − ½(1 + 1) = 79``, above every object it
    touches and below saturation. A bound at the best touched object would price it at 40.
    """
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.ones(ra.n_regions)
    for s in range(1000, 6500, 1000):
        c[_exon_mask(ra, s, s + 500)] = 40.0
    ne = int(boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0])
    m = _tidx(idx, "mrna")
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(
        _cal(ra, c, 1e-3, boundary=np.full(ne, 40.0), saturation=100.0), ra, idx, fl, _one_law(PMF)
    )
    assert eff[m] == pytest.approx(_mrna_length(fl[m], [40.0] * 6, [79.0] * 5), rel=1e-12)
    assert eff[m] > _mrna_length(fl[m], [40.0] * 6, [40.0] * 5) * (1.0 + 1e-6)


def test_a_junction_beside_an_object_above_saturation_may_reach_that_objects_level(multiexon_index):
    """Capture saturates per level: an object of a stronger probe reads above the library's saturated level, and a
    junction beside it is captured as that object is, not cut back to the library's level. Saturation 100, every exon
    and boundary at 300, the introns at the bulk: the sum ``599`` is bounded by the best touched object, 300, and the
    mRNA reads every object at 300. A bound at saturation alone would price the junctions at 100, under the exons beside
    them.
    """
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.ones(ra.n_regions)
    for s in range(1000, 6500, 1000):
        c[_exon_mask(ra, s, s + 500)] = 300.0
    ne = int(boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0])
    m = _tidx(idx, "mrna")
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(
        _cal(ra, c, 1e-3, boundary=np.full(ne, 300.0), saturation=100.0), ra, idx, fl, _one_law(PMF)
    )
    assert eff[m] == pytest.approx(300.0 * fl[m], rel=1e-12)


def test_one_factor_on_every_efficiency_is_one_factor_on_every_length(multiexon_index):
    """The unit is free: scaling every efficiency and the saturated level by one factor scales every transcript's length
    by exactly that factor (spliced, unspliced and nascent alike), so the EM, which reads every component in the same
    unit, cannot tell the units apart."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    rng = np.random.default_rng(5)
    c = rng.uniform(0.5, 400.0, ra.n_regions)
    ne = int(boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0])
    b = rng.uniform(0.5, 400.0, ne)
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    one = _lengths(_cal(ra, c, 1e-3, boundary=b, saturation=120.0), ra, idx, fl, _one_law(PMF))
    scaled = _lengths(
        _cal(ra, 7.25 * c, 1e-3, boundary=7.25 * b, saturation=7.25 * 120.0),
        ra,
        idx,
        fl,
        _one_law(PMF),
    )
    np.testing.assert_allclose(scaled, 7.25 * one, rtol=1e-12)


def test_an_intron_piece_too_short_to_contain_a_fragment_reads_its_far_boundary(short_intron_index):
    """A 50 bp intron holds no gDNA fragment against a 100–300 bp pmf, so its own efficiency is the population's and
    says nothing about it: the junction reads the boundary on the intron's far side instead (both splice sites'
    boundaries here, so their mean), and the intron's own efficiency moves nothing.
    """
    idx = short_intron_index
    ra = RegionArrays.from_index(idx)
    si = _tidx(idx, "si")
    intron = _exon_mask(ra, 1500, 1550)
    assert intron.sum() == 1
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500) | _exon_mask(ra, 1550, 2050)] = 1.0
    cb = _flank_mean(ra, c)
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(_cal(ra, c, 1.0, boundary=cb), ra, idx, fl, _one_law(PMF))
    c2 = c.copy()
    c2[intron] = 0.9
    eff2 = _lengths(_cal(ra, c2, 1.0, boundary=cb), ra, idx, fl, _one_law(PMF))
    assert eff2[si] == eff[si]
    contained, cut = 501.0 - MU, MU - 1.0
    expected = fl[si] * (2 * contained + cut * 0.5005) / (2 * contained + cut)
    assert eff[si] == pytest.approx(expected, rel=1e-9)


def test_a_junction_reads_the_objects_beside_it_on_a_multi_region_intron(multi_region_intron_index):
    """Every region and every boundary at its own efficiency, so no two objects agree and a junction read through a
    wrong index cannot land on the right number. Each transcript is two 500 bp exons, holding ``501 − E[w]`` of
    contained share per exon and ``E[w − 1]`` at its one junction. The junction's objects are found here by coordinate,
    never by the module's index arithmetic: its low boundary is where the exon below ends, its high boundary where the
    exon above starts, and the intron pieces beside them are the regions touching those two positions, read at their own
    efficiency when they hold a fragment (``tk``, 400 bp) and at the boundary on their far side when they cannot
    (``tj``, 50 bp).
    """
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

    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(_cal(ra, c, 1.0, boundary=cb), ra, idx, fl, _one_law(PMF))
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
    """A length model with mass at ``w = 0`` deposits no fragment there: the lengths equal those of the same pmf with
    that mass removed, for every transcript."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500)] = 1.0
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    with_zero = PMF.copy()
    with_zero[0] = 0.05
    np.testing.assert_allclose(
        _lengths(_cal(ra, c, 1.0), ra, idx, fl, _one_law(with_zero)),
        _lengths(_cal(ra, c, 1.0), ra, idx, fl, _one_law(PMF)),
        rtol=1e-12,
    )


def test_a_transcript_of_exons_shorter_than_a_fragment_is_read_through_its_cuts(tiny_index):
    """Ten 40 bp exons against a 100–300 bp pmf: no piece can contain a fragment, so the transcript's whole length sits
    on its junctions. Its exons' own efficiencies carry no weight (they hold no share), and with every junction priced
    at ½ (splice-site boundaries at ½, introns at ½) the factor is exactly ½. An object set that dropped the cuts would
    read this transcript as nothing.
    """
    idx = tiny_index
    ra = RegionArrays.from_index(idx)
    t = _tidx(idx, "tiny")
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    exons = np.zeros(ra.n_regions, dtype=bool)
    for i in range(10):
        exons |= _exon_mask(ra, 1000 + i * 1040, 1040 + i * 1040)
    half = np.full(ra.n_regions, 0.5)
    ne = int(boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0])
    eff = _lengths(_cal(ra, half, 1.0, boundary=np.full(ne, 0.5)), ra, idx, fl, _one_law(PMF))
    assert eff[t] == pytest.approx(0.5 * fl[t], rel=1e-12)
    loud = half.copy()
    loud[exons] = 0.9
    eff2 = _lengths(_cal(ra, loud, 1.0, boundary=np.full(ne, 0.5)), ra, idx, fl, _one_law(PMF))
    assert eff2[t] == eff[t]


def test_the_ruler_reads_the_efficiencies_and_nothing_from_the_counts(multiexon_index):
    """Two results with the same efficiencies and wildly different counts give the same lengths: the
    evidence was weighed upstream (`capture_efficiency`), and the ruler is geometry."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.3)
    c[_exon_mask(ra, 1000, 1500)] = 1.0
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    cal = _cal(ra, c, 1.0)
    loud = dataclasses.replace(
        cal,
        count_gdna_region=np.full(ra.n_regions, 1e6),
        count_gdna_boundary=np.full(cal.n_boundaries, 1e6),
        count_rna_region=np.full(ra.n_regions, 1e6),
    )
    np.testing.assert_array_equal(
        _lengths(cal, ra, idx, fl, _one_law(PMF)),
        _lengths(loud, ra, idx, fl, _one_law(PMF)),
    )


def test_no_nascent_mature_inversion_under_capture(multiexon_index):
    """On a field where capture adds over bases, a nascent parent holds its spliced child's captured exon plus more
    depleted bases, so ``eff(nascent) ≥ eff(mature)``; and the mature length contracts here."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500)] = 1.0
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(_cal(ra, c, 1.0), ra, idx, fl, _one_law(PMF))
    mrna, nasc = _tidx(idx, "mrna"), _tidx(idx, "nasc")
    assert eff[nasc] >= eff[mrna] - 1e-9
    assert eff[mrna] < fl[mrna] - 1e-6


# --- the efficiencies are gDNA's, carried to the EM's two tables ------------------------------------------------------


def _random_field(ra, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(0.5, 3.0, ra.n_regions)


def test_at_one_table_rna_reads_its_own_length_plus_the_excess_on_the_uncaptured_law(
    multiexon_index,
):
    """Both tables one law and gDNA's uncaptured law another: every transcript's length is its own fl-marginal length
    plus the excess gDNA's reading carries, ``fl + Σ_o (c_o − 1) share_o(g)``, whatever the RNA law the shares are
    weighted by. Capture adds to the uncaptured yield, and only the addition is gDNA's measurement."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    g = normal_pmf(mean=160.0, sd=40.0, lo=60, hi=300)
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(
        _cal(ra, np.full(ra.n_regions, 3.0), 1.0), ra, idx, fl, LengthTables(PMF, PMF, g, None)
    )
    on_g = transcript_objects(idx, ra, g)
    expected = fl.copy()
    np.add.at(expected, on_g.t, 2.0 * on_g.contained)
    np.add.at(expected, on_g.t[on_g.cut_row], 2.0 * on_g.cut_share)
    np.testing.assert_allclose(eff, expected, rtol=1e-12)


def test_a_reading_below_one_carries_nothing_and_never_goes_below_zero(multiexon_index):
    """Only an excess over 1 is carried. A depleted object read on gDNA's short law, carried to a scorer's table of
    longer fragments that places fewer starts on it, would carry an excess below −1 and enter below zero; it enters as
    read instead, so every length is ``fl`` times the one reading."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    g = normal_pmf(mean=160.0, sd=40.0, lo=60, hi=300)
    t_g = normal_pmf(mean=400.0, sd=30.0, lo=320, hi=480)
    c = 1e-3
    on_g, on_t = transcript_objects(idx, ra, g), transcript_objects(idx, ra, t_g)
    reach = on_t.contained > 0.0
    assert np.any(on_g.contained[reach] / on_t.contained[reach] > 1.0 / (1.0 - c))
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    eff = _lengths(
        _cal(ra, np.full(ra.n_regions, c), 1.0), ra, idx, fl, LengthTables(PMF, t_g, g, None)
    )
    np.testing.assert_allclose(eff, fl * c, rtol=1e-12)


def test_where_gdnas_table_is_its_uncaptured_law_the_efficiencies_are_read_as_given(
    multiexon_index,
):
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    g = normal_pmf(mean=160.0, sd=40.0, lo=60, hi=300)
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    cal = _cal(ra, _random_field(ra, 1), 1.0)
    np.testing.assert_array_equal(
        _lengths(cal, ra, idx, fl, LengthTables(PMF, g, g, None)),
        _lengths(cal, ra, idx, fl, _one_law(PMF)),
    )


def test_each_efficiency_is_carried_by_the_anchors_share_over_the_tables(multiexon_index):
    """``ĉ_o = 1 + (c_o − 1) share_o(g) / share_o(T_G)`` where ``c_o > 1``, weighted by RNA's shares: only an excess
    over 1 is carried. A reading of 1 or below enters as read, and so does an object ``T_G`` cannot reach (here every
    500 bp piece, against a law that only places longer fragments)."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    g = normal_pmf(mean=160.0, sd=40.0, lo=60, hi=300)
    t_g = normal_pmf(mean=560.0, sd=30.0, lo=510, hi=620)
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    cal = _cal(ra, _random_field(ra, 2), 1.0)
    eff = _lengths(cal, ra, idx, fl, LengthTables(PMF, t_g, g, None))
    on_r, on_t, on_g = (transcript_objects(idx, ra, p) for p in (PMF, t_g, g))
    assert np.any(on_t.contained == 0.0) and np.any(on_t.cut_share > 0.0)
    c = np.asarray(cal.gdna_capture_efficiency_region)
    c_cut = cut_efficiencies(
        on_r,
        ra,
        c,
        np.asarray(cal.gdna_capture_efficiency_boundary),
        np.asarray(cal.gdna_region_eff_len),
        1.0,
    )
    assert np.any(c[on_r.piece] < 1.0) and np.any(c[on_r.piece] > 1.0)
    piece = np.where(
        (on_t.contained > 0) & (c[on_r.piece] > 1.0),
        1.0
        + (c[on_r.piece] - 1.0)
        * on_g.contained
        / np.where(on_t.contained > 0, on_t.contained, 1.0),
        c[on_r.piece],
    )
    cut = np.where(
        (on_t.cut_share > 0) & (c_cut > 1.0),
        1.0 + (c_cut - 1.0) * on_g.cut_share / np.where(on_t.cut_share > 0, on_t.cut_share, 1.0),
        c_cut,
    )
    num, den = np.zeros(fl.size), np.zeros(fl.size)
    np.add.at(num, on_r.t, on_r.contained * piece)
    np.add.at(den, on_r.t, on_r.contained)
    np.add.at(num, on_r.t[on_r.cut_row], on_r.cut_share * cut)
    np.add.at(den, on_r.t[on_r.cut_row], on_r.cut_share)
    np.testing.assert_allclose(eff, fl * num / den, rtol=1e-12)
    # and it is not the efficiencies as read
    assert not np.allclose(eff, _lengths(cal, ra, idx, fl, _one_law(PMF)))


# --- kappa: the fraction of RNA's fragments capture added --------------------------------------------------------------


def test_rna_captured_fraction_weights_each_exons_captured_share_by_its_rna(multiexon_index):
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    exon = coarse_type_array(ra.signature) == int(RegionType.EXON)
    c = np.where(exon, 1.0, 5.0)
    rna = np.zeros(ra.n_regions)
    e = np.flatnonzero(exon)
    assert e.size >= 3
    c[e[:3]] = (2.0, 4.0, 0.8)
    rna[e[:3]] = (1.0, 3.0, 6.0)
    rna[~exon] = 100.0  # RNA outside exons does not count, whatever its region's efficiency
    cal = dataclasses.replace(_cal(ra, c, 1.0), count_rna_region=rna)
    expected = (1.0 * 0.5 + 3.0 * 0.75 + 6.0 * 0.0) / 10.0
    assert rna_captured_fraction(cal, ra) == pytest.approx(expected, rel=1e-15)


def test_rna_captured_fraction_is_zero_with_no_landscape_or_no_enrichment(multiexon_index):
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    rna = np.full(ra.n_regions, 10.0)
    assert (
        rna_captured_fraction(
            dataclasses.replace(_cal(ra, np.ones(ra.n_regions), None), count_rna_region=rna), ra
        )
        == 0.0
    )
    assert (
        rna_captured_fraction(
            dataclasses.replace(_cal(ra, np.full(ra.n_regions, 0.9), 1.0), count_rna_region=rna), ra
        )
        == 0.0
    )
    assert (
        rna_captured_fraction(
            dataclasses.replace(
                _cal(ra, np.full(ra.n_regions, 3.0), 1.0), count_rna_region=0.0 * rna
            ),
            ra,
        )
        == 0.0
    )


# --- the observable genome: a template's length is over the starts the mask leaves -----------------------------

# a spliced transcript with N inside its first exon and over its last exon's last 10 bases, its gene's nascent span N
# inside the intron too, and a single-exon transcript N end to end
_MASKED_GTF = (
    'chr1\ttest\texon\t1001\t1500\t.\t+\t.\tgene_id "gm"; transcript_id "tm";\n'
    'chr1\ttest\texon\t3001\t3600\t.\t+\t.\tgene_id "gm"; transcript_id "tm";\n'
    'chr1\ttest\texon\t4001\t4300\t.\t+\t.\tgene_id "gz"; transcript_id "tz";\n'
)
_MASKED_RUNS = {"chr1": [(1200, 1230), (2000, 2100), (3590, 3600), (4000, 4300)]}


@pytest.fixture(scope="module")
def masked_index(tmp_path_factory):
    return build_test_index(
        tmp_path_factory, _MASKED_GTF, genome_size=6000, name="masked", unknown=_MASKED_RUNS
    )


def _observable_starts(idx, t: int, pmf=PMF) -> float:
    """Brute force: every start of every length along template ``t``'s exons, walked through the exon structure into
    genomic bases and refused if its fragment touches a masked base, weighed by the pmf."""
    row = idx.t_df.iloc[t]
    exons = (
        np.array([[int(row["start"]), int(row["end"])]])
        if bool(row["is_synthetic"])
        else idx.get_exon_intervals(t)
    )
    genomic = np.concatenate([np.arange(a, b) for a, b in exons])
    masked = np.zeros(genomic.size, dtype=bool)
    for a, b in _MASKED_RUNS["chr1"]:
        masked |= (genomic >= a) & (genomic < b)
    touched = np.concatenate([[0], np.cumsum(masked)])
    total = 0.0
    for w in range(1, min(genomic.size, pmf.size - 1) + 1):
        s = np.arange(genomic.size - w + 1)
        total += pmf[w] * np.count_nonzero(touched[s + w] == touched[s])
    return total


def test_the_fixture_masks_what_it_says(masked_index):
    """Every template holds masked bases, the single-exon one end to end; the spliced transcript's nascent span is
    there to be read."""
    idx = masked_index
    assert idx.t_df["is_synthetic"].sum() == 1
    assert len(idx.unobservable_df) == len(_MASKED_RUNS["chr1"])


def test_a_templates_length_counts_only_the_starts_the_mask_leaves(masked_index):
    """Through the pipeline's own geometry step: each template's effective length is the pmf-weighted count of its
    starts whose fragment touches no masked base; a wholly masked template keeps the floor of one base, as does a
    transcript shorter than every fragment.
    """
    from rigel.frag_length_model import FragmentLengthModel
    from rigel.pipeline import transcript_geometry

    idx = masked_index
    ra = RegionArrays.from_index(idx)
    fl = FragmentLengthModel.from_pmf(PMF)
    _, est = transcript_geometry(idx, fl, _cal(ra, np.ones(ra.n_regions), None), ra, _one_law(PMF))
    for t in range(len(idx.t_df)):
        assert est.effective_lengths[t] == pytest.approx(
            max(_observable_starts(idx, t), 1.0), rel=1e-12
        ), idx.t_df["t_id"].iat[t]
    tz = _tidx(idx, "tz")
    assert _observable_starts(idx, tz) == 0.0 and est.effective_lengths[tz] == 1.0
    # the EM's yield too, before the estimator's own floor on its output: no landscape hands it over verbatim
    from rigel.calibration.capture_eff_length import observable_fl_eff_lengths

    assert observable_fl_eff_lengths(TranscriptFrame(idx, ra, _one_law(PMF)), fl)[tz] == 1.0


def test_a_masked_templates_shares_partition_its_observable_starts(masked_index):
    """Σ contained + Σ cut shares over a template's objects is its observable start count: no piece holds a masked
    base and no cut joins two pieces a masked base separates, so a fragment across the mask deposits nowhere."""
    idx = masked_index
    obj = transcript_objects(idx, RegionArrays.from_index(idx), PMF)
    total = np.zeros(len(idx.t_df))
    np.add.at(total, obj.t, obj.contained)
    np.add.at(total, obj.t[obj.cut_row], obj.cut_share)
    for t in range(len(idx.t_df)):
        assert total[t] == pytest.approx(_observable_starts(idx, t), rel=1e-9, abs=1e-12), idx.t_df[
            "t_id"
        ].iat[t]


def test_with_nothing_masked_the_lengths_and_shares_are_todays_bit_for_bit(multiexon_index):
    """An index whose genome holds no unknown base gives every template the closed-form length and the partition's
    pieces and cuts, to the bit."""
    from rigel.frag_length_model import FragmentLengthModel
    from rigel.pipeline import transcript_geometry

    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    assert len(idx.unobservable_df) == 0
    fl = FragmentLengthModel.from_pmf(PMF)
    _, est = transcript_geometry(idx, fl, _cal(ra, np.ones(ra.n_regions), None), ra, _one_law(PMF))
    lengths = idx.t_df["length"].to_numpy(np.int64)
    np.testing.assert_array_equal(
        est.effective_lengths, np.maximum(fl.compute_all_transcript_eff_lens(lengths), 1.0)
    )
    obj = transcript_objects(idx, ra, PMF)
    total = np.zeros(len(idx.t_df))
    np.add.at(total, obj.t, obj.contained)
    np.add.at(total, obj.t[obj.cut_row], obj.cut_share)
    np.testing.assert_allclose(total, contained_support(lengths, PMF), rtol=1e-12)


# --- the place weight: a multimapper's place weighs its candidates by the deposit masses this length is built from --


def _path_placements(exons, w):
    """Every placement of a ``w``-base fragment on the path through ``exons`` (genomic ``[s, e)``, ascending): its
    genomic start and end, the introns it crosses, and the genomic blocks it covers."""
    cum = np.concatenate([[0], np.cumsum([e - s for s, e in exons])])
    for u in range(int(cum[-1]) - w + 1):
        k0 = int(np.searchsorted(cum, u, side="right")) - 1
        k1 = int(np.searchsorted(cum, u + w - 1, side="right")) - 1
        start = exons[k0][0] + u - int(cum[k0])
        end = exons[k1][0] + u + w - int(cum[k1])
        blocks = [(max(s, start), min(e, end)) for s, e in exons[k0 : k1 + 1]]
        yield start, end, [(exons[k][1], exons[k + 1][0]) for k in range(k0, k1)], blocks


@pytest.mark.parametrize("junction_cap", [4.0, 9.0])
@pytest.mark.parametrize("fixture", ["multiexon_index", "masked_index"])
def test_a_templates_place_weights_carry_its_lengths_excess_on_the_uncaptured_law(
    fixture, junction_cap, request
):
    """``Σ_x g(ℓ_x) (w(x) − 1) = eff_em_t − fl_t`` over every observable placement of a template on its own path,
    junctions included, on a field at or above 1. ``w(x) = Σ_o m_o(x) c_o`` is the native deposit's masses at the
    reader's efficiencies (a junction at its price), the weight a multimapper's place carries (``place_log_weights``),
    and ``g`` is the law the efficiencies are read on (``tables.anchor``), not the scorer's table, here another law
    (``T_R = T_G ≠ g``). So the weight and the length the EM divides a transcript by carry one excess, in the frame
    calibration measures capture in; the uncaptured part of the length is RNA's own on its own table, which a place's
    weight, common to every candidate there, does not carry. A placement touching a masked base is never observed and
    weighs nothing. Both are handed one junction cap: at a cap other than the result's saturated level (9 against 4) the
    identity holds only if the weight and the length read the cap they were handed.
    """
    from rigel.calibration.capture_eff_length import (
        observable_fl_eff_lengths,
        place_efficiencies,
        place_log_weights,
    )
    from rigel.calibration.splice_graph import build_region_partition_arrays, build_sj_arrays
    from rigel.frag_length_model import FragmentLengthModel
    from rigel._bam_impl.testing import Accumulator
    from rigel.types import Strand

    idx = request.getfixturevalue(fixture)
    ra = RegionArrays.from_index(idx)
    rng = np.random.default_rng(5)
    n_boundaries = boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0]
    # regions at or above 1 and boundaries well above them, so every junction prices above 1 too
    cal = _cal(
        ra,
        rng.uniform(1.0, 2.0, ra.n_regions),
        1.0,
        boundary=rng.uniform(4.0, 6.0, n_boundaries),
        saturation=4.0,
    )
    g = np.zeros(PMF.shape[0])
    g[[120, 190, 260]] = [0.4, 0.35, 0.25]
    fl = observable_fl_eff_lengths(
        TranscriptFrame(idx, ra, _one_law(PMF)), FragmentLengthModel.from_pmf(PMF)
    )
    eff = _lengths(cal, ra, idx, fl, LengthTables(PMF, PMF, g, None), junction_cap)

    bounds, _ref_offsets, types = build_region_partition_arrays(idx)
    acc = Accumulator(
        region_bounds=np.ascontiguousarray(bounds, dtype=np.int64),
        region_types=np.ascontiguousarray(types, dtype=np.uint8),
        max_length=PMF.shape[0],
        ref=0,
    )
    sj = build_sj_arrays(idx)
    acc.set_sj(
        np.ascontiguousarray(sj.offsets, dtype=np.int32),
        np.ascontiguousarray(sj.boundary_right, dtype=np.int32),
        np.ascontiguousarray(sj.strand, dtype=np.int8),
    )
    c_region, c_boundary, c_sj, c_sj_crossing, sj_low_boundary, sj_high_boundary = (
        place_efficiencies(cal, ra, idx, junction_cap)
    )
    assert min(c_region.min(), c_boundary.min(), c_sj.min()) >= 1.0
    masked = [(int(a), int(b)) for _ref, a, b in np.asarray(ra.unobservable)]

    checked = 0
    for t in range(len(idx.t_df)):
        row = idx.t_df.iloc[t]
        exons = (
            [(int(row["start"]), int(row["end"]))]
            if bool(row["is_synthetic"])
            else [(int(a), int(b)) for a, b in idx.get_exon_intervals(t)]
        )
        offsets, kinds, ids, masses, f = [0], [], [], [], []
        for w in np.flatnonzero(g):
            for start, end, introns, blocks in _path_placements(exons, int(w)):
                if any(a < mb and ma < b for a, b in blocks for ma, mb in masked):
                    continue
                k, i, m = acc.placement_masses(
                    start=start, end=end, observed_introns=tuple(introns), sj_strand=int(Strand.POS)
                )
                kinds += list(k)
                ids += list(i)
                masses += list(m)
                offsets.append(len(kinds))
                f.append(g[w])
        if not f:
            continue  # wholly masked: no placement is observed
        log_w = np.asarray(
            place_log_weights(
                np.asarray(offsets, dtype=np.int64),
                np.asarray(kinds, dtype=np.int8),
                np.asarray(ids, dtype=np.int64),
                np.asarray(masses, dtype=np.float64),
                intergenic=np.zeros(len(f), dtype=bool),
                c_region=c_region,
                c_boundary=c_boundary,
                c_sj=c_sj,
                c_sj_crossing=c_sj_crossing,
                sj_low_boundary=sj_low_boundary,
                sj_high_boundary=sj_high_boundary,
            ),
            dtype=np.float64,
        )
        # ``log w`` is float32, as the scorer reads it: half an ulp of each moves each weight by at most that share
        w = np.exp(log_w)
        tolerance = 0.5 * float(np.finfo(np.float32).eps) * float(np.max(np.abs(log_w)))
        assert float(np.sum(np.asarray(f) * (w - 1.0))) == pytest.approx(
            eff[t] - fl[t], abs=max(tolerance, 1e-12) * float(np.sum(np.asarray(f) * w))
        ), row["t_id"]
        checked += 1
    assert checked >= 2


# --- the best part: a fragment that leaves a cut's exons across a junction is priced at that junction where it is higher


def _best_part_reference(idx, ra, cal, pmf, junction_cap, t: int) -> tuple[float, float, int]:
    """``(best part, own prices, placements across two junctions, junctions whose best object is passed over)`` for
    template ``t``, placement by placement.

    Every placement on the template's path goes through the native deposit rule; each mass is priced at its object's
    own efficiency (a junction at its rule price), and then, on a placement that crosses junctions, at the larger of
    that and the crossing level of the nearest junction the placement crosses on each side of the object. A junction's
    crossing level is written out here from the field: its two splice sites' boundaries less half the intron beside
    each, up to the cap, with no best-object term. The sums are the template's lengths under the two readings, on the
    one law ``pmf``. Plain Python over genomic coordinates: no array of the ruler's is read."""
    from rigel._bam_impl.testing import Accumulator
    from rigel.calibration.capture_eff_length import place_efficiencies
    from rigel.calibration.splice_graph import (
        build_region_partition_arrays,
        build_sj_arrays,
        build_sj_geometry_arrays,
    )
    from rigel.native import PLACE_BOUNDARY, PLACE_REGION, PLACE_SJ
    from rigel.types import Strand

    bounds, _ref_offsets, types = build_region_partition_arrays(idx)
    acc = Accumulator(
        region_bounds=np.ascontiguousarray(bounds, dtype=np.int64),
        region_types=np.ascontiguousarray(types, dtype=np.uint8),
        max_length=pmf.shape[0],
        ref=0,
    )
    sj = build_sj_arrays(idx)
    acc.set_sj(
        np.ascontiguousarray(sj.offsets, dtype=np.int32),
        np.ascontiguousarray(sj.boundary_right, dtype=np.int32),
        np.ascontiguousarray(sj.strand, dtype=np.int8),
    )
    c_region, c_boundary, c_sj = place_efficiencies(cal, ra, idx, junction_cap)[:3]
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    lo, _hi = boundary_region_indices(np.asarray(ra.ref_id))
    boundary_at = {int(ends[r]): e for e, r in enumerate(lo)}
    geo = build_sj_geometry_arrays(idx)
    intron_of = {
        (int(ends[s]), int(starts[d])): j
        for j, (s, d, strand) in enumerate(
            zip(geo.src_region, geo.dst_region, geo.strand, strict=True)
        )
        if int(strand) == int(Strand.POS)
    }
    exons = [(int(a), int(b)) for a, b in idx.get_exon_intervals(t)]

    def crossing_level(donor: int, acceptor: int) -> float:
        intron = np.flatnonzero((starts >= donor) & (ends <= acceptor))
        assert intron.shape[0] == 1 and cal.gdna_region_eff_len[intron[0]] > 0.0, (
            "the fixture's intron must be one region that holds a fragment"
        )
        summed = (
            c_boundary[boundary_at[donor]] + c_boundary[boundary_at[acceptor]] - c_region[intron[0]]
        )
        return float(min(max(summed, 0.0), junction_cap))

    crossing = {
        (exons[k][1], exons[k + 1][0]): crossing_level(exons[k][1], exons[k + 1][0])
        for k in range(len(exons) - 1)
    }
    passed_over = sum(float(c_sj[intron_of[i]]) > x for i, x in crossing.items())
    best_part = own_prices = 0.0
    across_two = 0
    for w in np.flatnonzero(pmf):
        for start, end, introns, _blocks in _path_placements(exons, int(w)):
            kinds, ids, masses = acc.placement_masses(
                start=start, end=end, observed_introns=tuple(introns), sj_strand=int(Strand.POS)
            )
            across_two += len(introns) >= 2
            crossed = [(d, a, crossing[(d, a)]) for d, a in introns]
            for kind, i, m in zip(kinds, ids, masses, strict=True):
                if kind == PLACE_REGION:
                    own, span = float(c_region[i]), None
                elif kind == PLACE_BOUNDARY:
                    own, span = float(c_boundary[i]), (int(ends[lo[i]]), int(ends[lo[i]]))
                else:
                    assert kind == PLACE_SJ
                    own = float(c_sj[i])
                    span = (int(ends[geo.src_region[i]]), int(starts[geo.dst_region[i]]))
                price = own
                if span is not None and crossed:
                    below = [(a, p) for d, a, p in crossed if a <= span[0]]
                    above = [(d, p) for d, a, p in crossed if d >= span[1]]
                    if below:
                        price = max(price, max(below)[1])
                    if above:
                        price = max(price, min(above)[1])
                best_part += pmf[w] * m * price
                own_prices += pmf[w] * m * own
    return best_part, own_prices, across_two, passed_over


def test_a_fragment_that_leaves_a_cuts_exons_across_a_junction_binds_its_best_part(pieces_index):
    """Where gDNA's fragment at a cut runs on into the intron, RNA's runs across the junction into the next exon, and
    gDNA has no such placement to read it by. Such a placement is RNA's crossing of that junction too, and a fragment
    binds its best probe part, so its deposit at the cut (a boundary inside the exon, or the junction before it) is
    priced at that junction's crossing level where it is higher; a placement that stays inside the cut's exons keeps
    the cut's own price. The crossing level is conservation's sum up to the cap: the best object beside a junction,
    which can price the junction itself above the cap, is one gDNA object's reading and is not passed on. ``t0``'s
    exons are 60, 40 and 60 bp, cut into pieces by its siblings, under fragments of 30, 70 and 110 bp: every placement
    is priced one by one from the native deposit, and the ruler's closed form must equal the sum, under a cap below the
    boundaries' readings (4) and one above them (9). At the higher cap, pricing each cut at its own efficiency alone
    reads this template about a tenth short."""
    idx = pieces_index
    ra = RegionArrays.from_index(idx)
    pmf = np.zeros(121)
    pmf[[30, 70, 110]] = [0.3, 0.4, 0.3]
    rng = np.random.default_rng(11)
    n_boundaries = boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0]
    # boundaries well above the regions, so the junctions price above the boundaries inside the exons beside them
    cal = _cal(
        ra,
        rng.uniform(1.0, 2.0, ra.n_regions),
        1.0,
        boundary=rng.uniform(2.0, 6.0, n_boundaries),
        saturation=4.0,
    )
    t0 = _tidx(idx, "t0")
    fl = contained_eff_length(idx.t_df["length"].to_numpy(np.float64), pmf)
    for cap, margin in ((4.0, 1.01), (9.0, 1.05)):
        best_part, own_prices, across_two, passed_over = _best_part_reference(
            idx, ra, cal, pmf, cap, t0
        )
        assert across_two > 0, "the fixture must place a fragment across two junctions"
        assert best_part > margin * own_prices, "the fixture must separate the two readings"
        if cap == 4.0:
            assert passed_over > 0, (
                "a junction must be priced at its best object, above its crossing level"
            )
        eff = _lengths(cal, ra, idx, fl, _one_law(pmf), cap)[t0]
        assert eff == pytest.approx(best_part, rel=1e-12)


def test_a_cuts_placements_are_split_on_the_law_its_excess_is_carried_from(pieces_index):
    """An efficiency is gDNA's yield on its uncaptured law ``g``, and only its excess is carried to the scorer's table
    (``_carried``). So the placements that leave a cut's exons are counted on ``g`` too: with RNA's table and gDNA's
    one law and ``g`` another, ``eff_em_t − fl_t = Σ_x g(ℓ_x) (w(x) − 1)`` over every placement of the template, each
    placement's deposits at the best part, both priced one by one from the field and weighed as a multimapper's place
    is. Splitting on RNA's table instead would carry the excess of placements ``g`` does not hold."""
    idx = pieces_index
    ra = RegionArrays.from_index(idx)
    table = np.zeros(121)
    table[[30, 70, 110]] = [0.3, 0.4, 0.3]
    g = np.zeros(121)
    g[[25, 55, 95]] = [0.2, 0.3, 0.5]
    rng = np.random.default_rng(11)
    n_boundaries = boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0]
    cal = _cal(
        ra,
        rng.uniform(1.0, 2.0, ra.n_regions),
        1.0,
        boundary=rng.uniform(2.0, 6.0, n_boundaries),
        saturation=4.0,
    )
    t0 = _tidx(idx, "t0")
    lengths = idx.t_df["length"].to_numpy(np.float64)
    fl = contained_eff_length(lengths, table)
    placements = float(contained_eff_length(lengths, g)[t0])
    # a cap under the boundaries' readings, where a junction is priced at its best object, and one above them
    for cap, margin in ((4.0, 1.01), (9.0, 1.05)):
        eff = _lengths(cal, ra, idx, fl, LengthTables(table, table, g, None), cap)[t0]
        best_part, own_prices, across_two, passed_over = _best_part_reference(
            idx, ra, cal, g, cap, t0
        )
        assert across_two > 0 and best_part > margin * own_prices
        assert passed_over > 0 or cap == 9.0
        assert eff - fl[t0] == pytest.approx(best_part - placements, rel=1e-12)
        # and the place weight is the same rule: Σ_x g(ℓ_x) (w(x) − 1) over the template's placements
        weights, law = _place_weights(idx, ra, cal, g, cap, t0)
        tolerance = 0.5 * float(np.finfo(np.float32).eps) * float(np.max(np.abs(np.log(weights))))
        assert float(np.sum(law * (weights - 1.0))) == pytest.approx(
            eff - fl[t0], abs=max(tolerance, 1e-12) * float(np.sum(law * weights))
        )


def _place_weights(idx, ra, cal, pmf, junction_cap, t: int) -> tuple[np.ndarray, np.ndarray]:
    """``(w, f)``: every placement of template ``t`` on its own path weighed as a multimapper's place is
    (``place_log_weights`` on the native deposit's masses), with the law's mass at its length."""
    from rigel._bam_impl.testing import Accumulator
    from rigel.calibration.capture_eff_length import place_efficiencies, place_log_weights
    from rigel.calibration.splice_graph import build_region_partition_arrays, build_sj_arrays
    from rigel.types import Strand

    bounds, _ref_offsets, types = build_region_partition_arrays(idx)
    acc = Accumulator(
        region_bounds=np.ascontiguousarray(bounds, dtype=np.int64),
        region_types=np.ascontiguousarray(types, dtype=np.uint8),
        max_length=pmf.shape[0],
        ref=0,
    )
    sj = build_sj_arrays(idx)
    acc.set_sj(
        np.ascontiguousarray(sj.offsets, dtype=np.int32),
        np.ascontiguousarray(sj.boundary_right, dtype=np.int32),
        np.ascontiguousarray(sj.strand, dtype=np.int8),
    )
    c_region, c_boundary, c_sj, c_sj_crossing, sj_low_boundary, sj_high_boundary = (
        place_efficiencies(cal, ra, idx, junction_cap)
    )
    exons = [(int(a), int(b)) for a, b in idx.get_exon_intervals(t)]
    offsets, kinds, ids, masses, f = [0], [], [], [], []
    for w in np.flatnonzero(pmf):
        for start, end, introns, _blocks in _path_placements(exons, int(w)):
            k, i, m = acc.placement_masses(
                start=start, end=end, observed_introns=tuple(introns), sj_strand=int(Strand.POS)
            )
            kinds += list(k)
            ids += list(i)
            masses += list(m)
            offsets.append(len(kinds))
            f.append(pmf[w])
    log_w = place_log_weights(
        np.asarray(offsets, dtype=np.int64),
        np.asarray(kinds, dtype=np.int8),
        np.asarray(ids, dtype=np.int64),
        np.asarray(masses, dtype=np.float64),
        intergenic=np.zeros(len(f), dtype=bool),
        c_region=c_region,
        c_boundary=c_boundary,
        c_sj=c_sj,
        c_sj_crossing=c_sj_crossing,
        sj_low_boundary=sj_low_boundary,
        sj_high_boundary=sj_high_boundary,
    )
    return np.exp(np.asarray(log_w, dtype=np.float64)), np.asarray(f)


def test_a_template_with_no_junction_prices_every_cut_at_its_own_boundary(multiexon_index):
    """A synthetic span and a single-exon transcript are gDNA's own placements: no fragment of theirs crosses a
    junction, so every cut keeps its boundary's efficiency and the length is the plain share-weighted sum, to the
    bit, whatever the junctions beside them price at."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    rng = np.random.default_rng(3)
    n_boundaries = boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0]
    cal = _cal(
        ra,
        rng.uniform(1.0, 2.0, ra.n_regions),
        1.0,
        boundary=rng.uniform(2.0, 6.0, n_boundaries),
        saturation=4.0,
    )
    fl = contained_eff_length(idx.t_df["length"].to_numpy(np.float64), PMF)
    obj = transcript_objects(idx, ra, PMF)
    spliced = np.zeros(len(idx.t_df), dtype=bool)
    spliced[obj.t[obj.cut_row[obj.is_junction]]] = True
    cut = np.zeros(len(idx.t_df), dtype=bool)
    cut[obj.t[obj.cut_row]] = True
    unspliced = np.flatnonzero(cut & ~spliced)
    assert unspliced.size and spliced.any(), (
        "the fixture must hold a cut template with no junction beside a spliced one"
    )
    eff = _lengths(cal, ra, idx, fl, _one_law(PMF), 9.0)
    c_region = np.asarray(cal.gdna_capture_efficiency_region)
    c_boundary = np.asarray(cal.gdna_capture_efficiency_boundary)
    lo, _hi = boundary_region_indices(np.asarray(ra.ref_id))
    boundary_after = {int(r): e for e, r in enumerate(lo)}
    for t in unspliced:
        rows = np.flatnonzero(obj.t == t)
        cuts = np.flatnonzero(obj.t[obj.cut_row] == t)
        assert not obj.is_junction[cuts].any()
        # the ruler's own order of additions: the pieces, then the cuts
        num = den = 0.0
        for share, c in zip(obj.contained[rows], c_region[obj.piece[rows]], strict=True):
            num += share * c
            den += share
        own = c_boundary[[boundary_after[int(p)] for p in obj.piece[obj.cut_row[cuts]]]]
        for share, c in zip(obj.cut_share[cuts], own, strict=True):
            num += share * c
            den += share
        assert eff[t] == fl[t] * (num / den)


# --- the junction cap: the saturated level, raised as far as RNA's own spliced split supports -------------------------

# gene ga: one transcript of three 500 bp exons, so the exon regions beside each of its two junctions are its own.
# gene gb: two isoforms sharing their first exon, which is therefore neither junction's own; each second exon is.
_CAP_GTF = "".join(
    f'chr1\ttest\texon\t{s + 1}\t{s + 500}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{tid}";\n'
    for gene, tid, starts in (
        ("ga", "a1", (1000, 2500, 4000)),
        ("gb", "b1", (6000, 7500)),
        ("gb", "b2", (6000, 9000)),
    )
    for s in starts
)

# two transcripts on opposite strands with the same two exons: one intron, one junction per strand
_ANTISENSE_GTF = "".join(
    f'chr1\ttest\texon\t{s + 1}\t{s + 500}\t.\t{strand}\t.\tgene_id "g{name}"; transcript_id "t{name}";\n'
    for strand, name in (("+", "p"), ("-", "n"))
    for s in (1000, 2500)
)


@pytest.fixture(scope="module")
def cap_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _CAP_GTF, genome_size=11000, name="junctioncap")


@pytest.fixture(scope="module")
def antisense_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _ANTISENSE_GTF, genome_size=5000, name="antisense")


def _exon_regions(ra: RegionArrays) -> np.ndarray:
    return coarse_type_array(ra.signature) == int(RegionType.EXON)


def _split_at(idx, ra, c, price, tables=None, rate=None):
    """``(rna_region, sj_mass)``: the unspliced and spliced RNA of a library in which every annotated transcript's
    fragments split between its objects exactly as the ruler's terms say with every junction priced at ``price``. A
    transcript at rate ``ρ`` leaves ``ρ · (its contained share) · ĉ`` in an exon region at efficiency ``c`` and
    ``ρ · (its junction share) · p̂`` on a junction, ``ĉ = 1 + (c − 1)·k`` above 1 with ``k`` the object's share on
    gDNA's uncaptured law over its share on gDNA's table, as the ruler carries an excess, and ``c`` as read at 1 or
    below. ``rate`` is ``ρ`` by transcript id, 1 where not given."""
    from rigel.calibration.splice_graph import build_sj_geometry_arrays

    tables = _one_law(PMF) if tables is None else tables
    pieces = transcript_pieces(idx, ra)
    obj = object_shares(pieces, tables.rna)
    on_anchor, on_table = object_shares(pieces, tables.anchor), object_shares(pieces, tables.gdna)

    def carried(level, k):
        # the ruler's rule: only an excess over 1 is carried; a level of 1 or below enters as read
        return np.where(level > 1.0, 1.0 + (level - 1.0) * k, level)

    piece_term = obj.contained * carried(c[obj.piece], on_anchor.contained / on_table.contained)
    cut_term = obj.cut_share * carried(np.float64(price), on_anchor.cut_share / on_table.cut_share)
    geo = build_sj_geometry_arrays(idx)
    rna = np.zeros(ra.n_regions)
    spliced = np.zeros(np.asarray(geo.src_region).shape[0])
    tdf = idx.t_df
    for t in np.flatnonzero(~tdf["is_synthetic"].to_numpy(dtype=bool)):
        rho = (rate or {}).get(tdf["t_id"].iat[t], 1.0)
        rows = np.flatnonzero(obj.t == t)
        np.add.at(rna, obj.piece[rows], rho * piece_term[rows])
        for k in np.flatnonzero(obj.is_junction & (obj.t[obj.cut_row] == t)):
            below, above = obj.piece[obj.cut_row[k]], obj.piece[obj.cut_row[k] + 1]
            row = np.flatnonzero(
                (np.asarray(geo.src_region) == below)
                & (np.asarray(geo.dst_region) == above)
                & (np.asarray(geo.strand) == tdf["strand"].iat[t])
            )
            assert row.shape[0] == 1, "the fixture's junction is not one row of the sj axis"
            spliced[row[0]] += rho * cut_term[k]
    return rna, spliced


def _captured(ra, exon_level: float, intron_level: float = 1.0) -> np.ndarray:
    return np.where(_exon_regions(ra), exon_level, intron_level).astype(np.float64)


def _n_boundaries(ra) -> int:
    return int(boundary_region_indices(np.asarray(ra.ref_id))[0].shape[0])


# every exon and boundary at 400 and the introns at the bulk: conservation prices a junction at 400 + 400 − 1 = 799
# and the best object beside it is 400, so a cap between them is what the junction reads
_LEVEL, _SUMMED = 400.0, 799.0
_RATES = {"a1": 50.0, "b1": 8.0, "b2": 300.0}


def _cap_library(idx, ra, price, saturation=100.0, rna_scale=None, tables=None):
    c = _captured(ra, _LEVEL)
    rna, spliced = _split_at(idx, ra, c, price, tables, _RATES)
    if rna_scale is not None:
        rna = rna * rna_scale
    return _cal(
        ra,
        c,
        1e-3,
        boundary=np.full(_n_boundaries(ra), _LEVEL),
        saturation=saturation,
        rna_region=rna,
        sj_mass=spliced,
    )


def _grid_step(lo: float, hi: float) -> float:
    """The relative spacing of neighbouring levels on the cap's search grid from ``lo`` to ``hi``."""
    from rigel.config import CONSTANTS

    return (hi / lo) ** (1.0 / (CONSTANTS.calibration.junction_cap_grid_points - 1)) - 1.0


def test_the_junction_cap_is_the_level_rnas_spliced_split_supports(cap_index):
    """The saturated level is 100, and every junction's spliced mass sits against its exons' unspliced RNA as a
    junction price of 600 says: the cap is 600, to the search grid's spacing, whatever each transcript's abundance. The
    saturated level alone would price every junction at the best object beside it, 400.
    """
    idx = cap_index
    ra = RegionArrays.from_index(idx)
    cap = fit_junction_cap(_cap_library(idx, ra, 600.0), TranscriptFrame(idx, ra, _one_law(PMF)))
    assert cap == pytest.approx(600.0, rel=_grid_step(100.0, _SUMMED))
    cal = _cap_library(idx, ra, 600.0)
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    a1 = _tidx(idx, "a1")
    # the transcript's length at the fitted cap is its exons at 400 and its junctions at the cap
    assert _lengths(cal, ra, idx, fl, _one_law(PMF), cap)[a1] == pytest.approx(
        _mrna_length(fl[a1], [_LEVEL] * 3, [cap] * 2), rel=1e-12
    )
    assert _lengths(cal, ra, idx, fl, _one_law(PMF))[a1] == pytest.approx(
        _mrna_length(fl[a1], [_LEVEL] * 3, [_LEVEL] * 2), rel=1e-12
    )


def test_the_junction_cap_reads_the_split_on_the_multipliers_the_ruler_applies(cap_index):
    """Where gDNA's table is not its uncaptured law the ruler carries every efficiency by the anchor's share over the
    table's, so what it multiplies an efficiency by differs object by object. A library split by those multipliers at a
    price of 600 reads a cap of 600; read on the uncarried shares the same library gives another level.
    """
    idx = cap_index
    ra = RegionArrays.from_index(idx)
    g = np.zeros(PMF.shape[0])
    g[[120, 190, 260]] = [0.4, 0.35, 0.25]
    tables = LengthTables(PMF, PMF, g, None)
    cal = _cap_library(idx, ra, 600.0, tables=tables)
    step = _grid_step(100.0, _SUMMED)
    assert fit_junction_cap(cal, TranscriptFrame(idx, ra, tables)) == pytest.approx(600.0, rel=step)
    assert (
        abs(fit_junction_cap(cal, TranscriptFrame(idx, ra, _one_law(PMF))) / 600.0 - 1.0)
        > 2.0 * step
    )


def test_the_junction_cap_reads_the_split_on_the_rulers_own_terms_near_one(cap_index):
    """Near the bulk the carry's form matters: an exon at 3 enters as ``1 + 2k`` and a junction priced 4 as
    ``1 + 3k``, not ``3k`` and ``4k``. A library split by the ruler's terms at a price of 4, with gDNA's law apart from
    both tables, reads a cap of 4, to the search grid's spacing."""
    idx = cap_index
    ra = RegionArrays.from_index(idx)
    g = np.zeros(PMF.shape[0])
    g[[120, 190, 260]] = [0.4, 0.35, 0.25]
    tables = LengthTables(PMF, PMF, g, None)
    c = _captured(ra, 3.0)
    rna, spliced = _split_at(idx, ra, c, 4.0, tables, _RATES)
    cal = _cal(
        ra,
        c,
        1e-3,
        boundary=np.full(_n_boundaries(ra), 3.0),
        saturation=1.5,
        rna_region=rna,
        sj_mass=spliced,
    )
    assert fit_junction_cap(cal, TranscriptFrame(idx, ra, tables)) == pytest.approx(
        4.0, rel=_grid_step(1.5, 5.0)
    )


def test_the_spliced_counts_raise_the_junction_cap_and_never_lower_it(cap_index):
    """Every exon and boundary at 60, the saturated level 100: a junction is priced ``min(119, max(cap, 60))``, so a cap
    between 60 and 119 is what it reads. The spliced counts say 70, under the saturated level, and the cap stays the
    saturated level exactly: gDNA's law, through which the split's opportunity is read, can be wrong on aligned input,
    and a cap the counts could lower would then fall to the best object beside every junction. Counts that say 110
    raise it.
    """
    idx = cap_index
    ra = RegionArrays.from_index(idx)
    c = _captured(ra, 60.0)

    def cap_at(price):
        rna, spliced = _split_at(idx, ra, c, price, rate=_RATES)
        cal = _cal(
            ra,
            c,
            1e-3,
            boundary=np.full(_n_boundaries(ra), 60.0),
            saturation=100.0,
            rna_region=rna,
            sj_mass=spliced,
        )
        return fit_junction_cap(cal, TranscriptFrame(idx, ra, _one_law(PMF)))

    assert cap_at(70.0) == 100.0
    assert cap_at(110.0) == pytest.approx(110.0, rel=_grid_step(100.0, 119.0))


def test_a_cap_the_counts_cannot_tell_apart_stays_the_saturated_level(cap_index):
    """Every exon and boundary at 400, the saturated level 100, and spliced counts that say 150: each junction already
    reads the best object beside it, 400, at every cap up to 400, so the likelihood is flat from 100 to 400 and falls
    above. Of the caps the counts cannot tell apart the fit keeps the first, the saturated level exactly: a junction
    outside the fit, with weaker objects beside it, would read any larger one.
    """
    idx = cap_index
    ra = RegionArrays.from_index(idx)
    assert (
        fit_junction_cap(_cap_library(idx, ra, 150.0), TranscriptFrame(idx, ra, _one_law(PMF)))
        == 100.0
    )


def test_the_junction_cap_is_the_saturated_level_where_no_capture_is_read(cap_index):
    """A library whose landscape has no level above the bulk has a saturated level of exactly 1, and the cap is 1 to the
    bit. The reader may still read an intron under 1, and then conservation prices a junction above 1 (here
    ``1 + 1 − 0.2 = 1.8``) and the spliced counts ask for 1.5; with no capture read there is none for a junction to
    exceed, and a cap fitted there would move a library that was never captured. With no landscape at all the cap is 1.
    """
    idx = cap_index
    ra = RegionArrays.from_index(idx)
    c = _captured(ra, 1.0, intron_level=0.2)
    rna, spliced = _split_at(idx, ra, c, 1.5, rate=_RATES)

    def cap_at(saturation):
        cal = _cal(
            ra,
            c,
            1e-3,
            boundary=np.ones(_n_boundaries(ra)),
            saturation=saturation,
            rna_region=rna,
            sj_mass=spliced,
        )
        return fit_junction_cap(cal, TranscriptFrame(idx, ra, _one_law(PMF)))

    assert cap_at(1.0) == 1.0
    # the same counts under a library that does read capture (a saturated level a hair above 1) are fitted
    assert cap_at(1.0 + 1e-9) == pytest.approx(1.5, rel=_grid_step(1.0, 1.8))
    none = _cal(ra, np.ones(ra.n_regions), None, rna_region=rna, sj_mass=spliced)
    assert fit_junction_cap(none, TranscriptFrame(idx, ra, _one_law(PMF))) == 1.0


def test_a_junctions_evidence_is_the_exon_regions_only_its_own_transcripts_cover(cap_index):
    """The two isoforms of ``gb`` share their first exon, so that region's unspliced RNA is a mixture of both and belongs
    to neither junction: a thousandfold change in it moves the cap by nothing. Each isoform's second exon is its
    junction's own, and the abundant isoform's moves the cap. The fit holds ``gb`` alone (``ga`` carries no fragment).
    """
    idx = cap_index
    ra = RegionArrays.from_index(idx)
    c = _captured(ra, _LEVEL)
    rna, spliced = _split_at(idx, ra, c, 600.0, rate={"a1": 0.0, "b1": 8.0, "b2": 300.0})
    start, end = np.asarray(ra.start), np.asarray(ra.end)
    shared = np.flatnonzero((start == 6000) & (end == 6500))
    own = np.flatnonzero((start == 9000) & (end == 9500))
    assert shared.shape[0] == 1 and own.shape[0] == 1 and rna[shared[0]] > 0.0

    def cap_with(region=None, factor=1.0):
        r = rna.copy()
        if region is not None:
            r[region] *= factor
        cal = _cal(
            ra,
            c,
            1e-3,
            boundary=np.full(_n_boundaries(ra), _LEVEL),
            saturation=100.0,
            rna_region=r,
            sj_mass=spliced,
        )
        return fit_junction_cap(cal, TranscriptFrame(idx, ra, _one_law(PMF)))

    base = cap_with()
    assert base == pytest.approx(600.0, rel=_grid_step(100.0, _SUMMED))
    assert cap_with(shared[0], 1000.0) == base
    assert cap_with(shared[0], 0.0) == base
    assert cap_with(own[0], 4.0) < base * 0.9


@pytest.fixture(scope="module")
def cap_masked_index(tmp_path_factory):
    """The cap fixture with an unobservable run inside ``b2``'s own second exon, which the mask cuts into two parts."""
    return build_test_index(
        tmp_path_factory,
        _CAP_GTF,
        genome_size=11000,
        name="junctioncap_masked",
        unknown={"chr1": [(9200, 9230)]},
    )


def test_a_junctions_own_exon_keeps_its_evidence_where_the_mask_splits_it(cap_masked_index):
    """The mask cuts ``b2``'s own second exon into two observable parts of one region and one transcript. The footprint
    counts each transcript once per region and gives the region the sum of its parts' shares, so with ``b2`` the only
    isoform carrying fragments its junction is still fitted: the cap is the price its split says, and a fourfold change
    in that exon's unspliced RNA moves it. Counted part by part, one transcript would cover the exon twice, the exon
    would belong to no junction, and nothing would be left to fit."""
    idx = cap_masked_index
    ra = RegionArrays.from_index(idx)
    start, end = np.asarray(ra.start), np.asarray(ra.end)
    own = np.flatnonzero((start == 9000) & (end == 9500))
    assert own.shape[0] == 1
    pieces = transcript_pieces(idx, ra)
    b2 = _tidx(idx, "b2")
    assert np.count_nonzero((pieces.t == b2) & (pieces.piece == own[0])) == 2
    c = _captured(ra, _LEVEL)
    rna, spliced = _split_at(idx, ra, c, 600.0, rate={"a1": 0.0, "b1": 0.0, "b2": 300.0})

    def cap_with(factor=1.0):
        r = rna.copy()
        r[own[0]] *= factor
        cal = _cal(
            ra,
            c,
            1e-3,
            boundary=np.full(_n_boundaries(ra), _LEVEL),
            saturation=100.0,
            rna_region=r,
            sj_mass=spliced,
        )
        return fit_junction_cap(cal, TranscriptFrame(idx, ra, _one_law(PMF)))

    base = cap_with()
    assert base == pytest.approx(600.0, rel=_grid_step(100.0, _SUMMED))
    assert cap_with(4.0) < base * 0.9


# gene gc: ``c1`` is two 500 bp exons around a 40 bp one, so nearly every fragment across its first junction leaves
# the short exon across the second; ``c2`` covers the last exon too, so that exon is no junction's own
_LEAVING_GTF = (
    'chr1\ttest\texon\t1001\t1500\t.\t+\t.\tgene_id "gc"; transcript_id "c1";\n'
    'chr1\ttest\texon\t2001\t2040\t.\t+\t.\tgene_id "gc"; transcript_id "c1";\n'
    'chr1\ttest\texon\t3001\t3500\t.\t+\t.\tgene_id "gc"; transcript_id "c1";\n'
    'chr1\ttest\texon\t3001\t3500\t.\t+\t.\tgene_id "gc"; transcript_id "c2";\n'
)


@pytest.fixture(scope="module")
def leaving_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _LEAVING_GTF, genome_size=5000, name="leaving")


def test_the_junction_cap_reads_a_junction_through_the_placements_that_leave_across_the_next(
    leaving_index,
):
    """The fit and the lengths state one model. ``c1``'s first junction has an exon of its own and a price the cap
    cannot move (its two splice sites read 50, so conservation holds it at 99, under the saturated level); its second
    has no exon of its own, and is itself priced at the best object beside it, a splice site reading 700. A fragment
    across the first that leaves the 40 bp exon across the second takes the second's crossing level, its sum up to
    the cap, so the first junction's spliced mass against its own exon's unspliced RNA reads the cap: a library drawn
    at the ruler's own prices under a cap of 600 gives 600 back. A fit that read each junction at its own price would
    hold no junction and return the saturated level, 100; one that passed the second junction's 700 on would find
    nothing for the cap to move below 700."""
    from rigel.calibration.splice_graph import build_sj_geometry_arrays

    idx = leaving_index
    ra = RegionArrays.from_index(idx)
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    c = _captured(ra, _LEVEL)
    lo, _hi = boundary_region_indices(np.asarray(ra.ref_id))
    boundary = np.full(_n_boundaries(ra), _LEVEL)
    boundary[np.isin(ends[lo], (1500, 2000))] = 50.0
    boundary[ends[lo] == 3000] = 700.0
    tables = _one_law(PMF)

    def library(cap):
        """The unspliced RNA of each exon region and each junction's spliced mass, with ``c1`` at rate 40 and every
        object at the ruler's own price under ``cap``."""
        gdna_contained = contained_eff_length(np.asarray(ra.region_size_bp, dtype=np.float64), PMF)
        obj = transcript_objects(idx, ra, PMF)
        price = cut_efficiencies(obj, ra, c, boundary, gdna_contained, cap)
        geo = build_sj_geometry_arrays(idx)
        rna = np.zeros(ra.n_regions)
        spliced = np.zeros(np.asarray(geo.src_region).shape[0])
        t = _tidx(idx, "c1")
        rows = np.flatnonzero(obj.t == t)
        np.add.at(rna, obj.piece[rows], 40.0 * obj.contained[rows] * c[obj.piece[rows]])
        for k in np.flatnonzero(obj.is_junction & (obj.t[obj.cut_row] == t)):
            row = np.flatnonzero(
                (np.asarray(geo.src_region) == obj.piece[obj.cut_row[k]])
                & (np.asarray(geo.dst_region) == obj.piece[obj.cut_row[k] + 1])
            )
            assert row.shape[0] == 1
            spliced[row[0]] += 40.0 * obj.cut_share[k] * price[k]
        return rna, spliced, obj, price

    rna, spliced, obj, price = library(600.0)
    junctions = np.flatnonzero(obj.is_junction & (obj.t[obj.cut_row] == _tidx(idx, "c1")))
    first = junctions[np.argmin(starts[obj.piece[obj.cut_row[junctions] + 1]])]
    # the fixture's premise: the first junction's own price is conservation's 99 and most of its share leaves; the
    # second is priced at its best object, above the cap
    second = junctions[junctions != first][0]
    assert obj.leave_hi_share[first] > 0.5 * obj.cut_share[first]
    assert 99.0 < price[first] < 600.0
    assert price[second] == pytest.approx(700.0)
    cal = _cal(ra, c, 1e-3, boundary=boundary, saturation=100.0, rna_region=rna, sj_mass=spliced)
    # the grid runs to the second junction's sum, its two splice sites less the intron: 400 + 700 − 1
    assert fit_junction_cap(cal, TranscriptFrame(idx, ra, tables)) == pytest.approx(
        600.0, rel=_grid_step(100.0, 1099.0)
    )


def test_a_junction_with_no_exon_region_of_its_own_is_not_fitted(multiexon_index):
    """In the multi-exon fixture an annotated single-exon transcript covers every exon of the spliced one, so no exon
    region is any junction's own and nothing anchors a split: the cap is the saturated level exactly, though the
    spliced counts are those of a junction price of 600.
    """
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = _captured(ra, _LEVEL)
    rna, spliced = _split_at(idx, ra, c, 600.0)
    assert spliced.sum() > 0.0 and rna.sum() > 0.0
    cal = _cal(
        ra,
        c,
        1e-3,
        boundary=np.full(_n_boundaries(ra), _LEVEL),
        saturation=100.0,
        rna_region=rna,
        sj_mass=spliced,
    )
    assert fit_junction_cap(cal, TranscriptFrame(idx, ra, _one_law(PMF))) == 100.0


@pytest.mark.filterwarnings("ignore:.*strand-coincident splice junction")
@pytest.mark.parametrize("abundant", ["tp", "tn"])
def test_a_junction_is_its_two_regions_and_its_strand(antisense_index, abundant):
    """Two transcripts on opposite strands share both exons, so their one intron is two junctions, one per strand (an
    annotation the index warns about and keeps), and each exon region is covered by a transcript that does not use the
    other strand's junction: neither junction has a region of its own, and the cap is the saturated level. Read by its
    two regions alone, the two would be one junction with both transcripts as its users and both exons as its own, and
    with one strand thirty times the other its spliced mass would read a cap near 580.
    """
    idx = antisense_index
    ra = RegionArrays.from_index(idx)
    c = _captured(ra, _LEVEL)
    rna, spliced = _split_at(idx, ra, c, 600.0, rate={"tp": 1.0, "tn": 1.0, abundant: 30.0})
    assert spliced.shape[0] == 2 and np.all(spliced > 0.0)
    cal = _cal(
        ra,
        c,
        1e-3,
        boundary=np.full(_n_boundaries(ra), _LEVEL),
        saturation=100.0,
        rna_region=rna,
        sj_mass=spliced,
    )
    assert fit_junction_cap(cal, TranscriptFrame(idx, ra, _one_law(PMF))) == 100.0


def test_the_lengths_and_the_place_weights_read_the_cap_they_are_handed(cap_index):
    """The cap is the pipeline's to fit once and hand to both of its readers; neither reads the result's saturated level
    in its place. Saturated level 100, exons and boundaries at 400, cap 600: every junction's place weight is 600, and
    ``a1``'s length holds its junctions at 600. Read at the saturated level both would hold 400.
    """
    from rigel.calibration.capture_eff_length import place_efficiencies

    idx = cap_index
    ra = RegionArrays.from_index(idx)
    c = _captured(ra, _LEVEL)
    n_sj = _split_at(idx, ra, c, 600.0)[1].shape[0]
    cal = _cal(
        ra,
        c,
        1e-3,
        boundary=np.full(_n_boundaries(ra), _LEVEL),
        saturation=100.0,
        sj_mass=np.zeros(n_sj),
    )
    c_sj = place_efficiencies(cal, ra, idx, 600.0)[2]
    np.testing.assert_array_equal(c_sj, np.full(n_sj, 600.0))
    fl = contained_support(idx.t_df["length"].to_numpy(), PMF)
    a1 = _tidx(idx, "a1")
    assert transcript_capture_eff_lengths(cal, TranscriptFrame(idx, ra, _one_law(PMF)), fl, 600.0)[
        a1
    ] == pytest.approx(_mrna_length(fl[a1], [_LEVEL] * 3, [600.0] * 2), rel=1e-12)


def test_a_run_builds_the_pieces_once_and_each_tables_shares_once(tmp_path_factory, monkeypatch):
    """The junction cap's fit, the observable lengths and the capture-contracted lengths read one frame: a run builds
    the transcripts' pieces once and each length table's shares once. On a whole genome the pieces are ten million rows
    and a table's shares a second's work, so a reader that rebuilt them would cost seconds and move no number."""
    from _prior_toy import build_toy

    import rigel.calibration.capture_eff_length as CEL
    from rigel.config import PipelineConfig
    from rigel.pipeline import run_pipeline

    toy = build_toy(tmp_path_factory)
    built = {"pieces": 0, "shares": []}
    pieces, shares = CEL.transcript_pieces, CEL.object_shares

    def count_pieces(index, region_arrays):
        built["pieces"] += 1
        return pieces(index, region_arrays)

    def count_shares(frame_pieces, pmf, *rest):
        built["shares"].append(np.asarray(pmf, dtype=np.float64).tobytes())
        return shares(frame_pieces, pmf, *rest)

    monkeypatch.setattr(CEL, "transcript_pieces", count_pieces)
    monkeypatch.setattr(CEL, "object_shares", count_shares)
    run_pipeline(str(toy.bam_path), toy.index, PipelineConfig())
    assert built["pieces"] == 1
    assert built["shares"], "the toy never priced a transcript: no shares were built"
    assert len(built["shares"]) == len(set(built["shares"]))


def test_the_pipeline_fits_the_cap_once_and_hands_it_to_both_readers(tmp_path_factory, monkeypatch):
    """One run of the pipeline on a small contaminated library: the cap the fit returns is the one the transcripts'
    lengths take, the one a multimapper's place weights take, and the one the result reports. A reader handed the
    result's saturated level, or a second fit, would price a junction differently in the length and in the weight.
    """
    from _prior_toy import build_toy

    import rigel.calibration.capture_eff_length as CEL
    from rigel.config import PipelineConfig
    from rigel.pipeline import run_pipeline

    toy = build_toy(tmp_path_factory)
    fitted, handed = [], {}

    def fit(calibration, frame):
        fitted.append(calibration.gdna_capture_saturation)
        return 7.25

    ruler, place = CEL.transcript_capture_eff_lengths, CEL.place_efficiencies

    def spy_ruler(calibration, frame, fl_eff_lengths, junction_cap):
        handed["lengths"] = junction_cap
        return ruler(calibration, frame, fl_eff_lengths, junction_cap)

    def spy_place(calibration, region_arrays, index, junction_cap):
        handed["place weights"] = junction_cap
        return place(calibration, region_arrays, index, junction_cap)

    monkeypatch.setattr(CEL, "fit_junction_cap", fit)
    monkeypatch.setattr(CEL, "transcript_capture_eff_lengths", spy_ruler)
    monkeypatch.setattr(CEL, "place_efficiencies", spy_place)
    result = run_pipeline(str(toy.bam_path), toy.index, PipelineConfig())
    assert len(fitted) == 1
    assert handed == {"lengths": 7.25, "place weights": 7.25}
    assert result.junction_cap == 7.25
