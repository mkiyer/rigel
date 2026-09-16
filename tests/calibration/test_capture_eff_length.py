"""``capture_eff_length``: a transcript's effective length under capture is its own bases at their
pieces' efficiencies, weighted by the fragment-length end taper.

The geometry half: ``transcript_piece_lengths`` must map every transcript to every piece its exons
overlap, on any partition it is handed (the coarse partition built here by hand still has exon
boundaries in a region's interior, which the shipped partition no longer produces — losing that fixture
would retire the guard), and the taper-weighted base counts of a transcript's pieces must sum to its
fl-marginal length exactly, so the per-base frame is a partition of the start count and never a second
definition of it. The contraction half: every efficiency 1 (a field with no reference) returns the
FL-marginal lengths bit-identically; a field contracts and never expands; a piece shorter than a
fragment carries its bases' weight and nothing else; the ruler reads the result's efficiencies and
re-derives nothing from the counts; a nascent parent never reads shorter than its spliced child.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
import pytest
from _index_builder import build_test_index

from rigel.calibration.capture_eff_length import (
    transcript_capture_eff_lengths,
    transcript_piece_lengths,
)
from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig


def _pmf(mean=200.0, sd=50.0, lo=100, hi=300) -> np.ndarray:
    w = np.arange(hi + 1, dtype=np.float64)
    p = np.exp(-0.5 * ((w - mean) / sd) ** 2)
    p[(w < lo) | (w > hi)] = 0.0
    return p / p.sum()


PMF = _pmf()


def _fl(lengths, pmf=PMF) -> np.ndarray:
    """The fl-marginal length ``Σ_w f(w)(L − w + 1)⁺`` per transcript."""
    w = np.arange(pmf.shape[0], dtype=np.float64)
    L = np.asarray(lengths, dtype=np.float64)[:, None]
    return (pmf[None, :] * np.maximum(L - w[None, :] + 1.0, 0.0)).sum(1)


def _cal(region_arrays: RegionArrays, efficiency, reference: float | None) -> CalibrationResult:
    """THE fixture: a result whose per-piece capture efficiencies are stated outright, with counts and
    supports that carry nothing — the ruler must read the efficiencies and nothing else. ``reference``
    ``None`` is a field with no enriched mode, where the efficiencies must all be 1."""
    n = int(region_arrays.n_regions)
    lo, _hi = boundary_region_indices(np.asarray(region_arrays.ref_id))
    ne = lo.shape[0]
    z = np.zeros(n)
    ez = np.zeros(ne)
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
        gdna_region_eff_len=np.asarray(region_arrays.region_size_bp, dtype=np.float64),
        gdna_boundary_eff_len=np.full(ne, 180.0),
        rna_region_eff_len=np.asarray(region_arrays.region_size_bp, dtype=np.float64),
        rna_boundary_eff_len=np.full(ne, 180.0),
        gdna_frac_region=z.copy(),
        rna_pos_frac_region=z.copy(),
        rna_neg_frac_region=z.copy(),
        gdna_frac_boundary=ez.copy(),
        rna_pos_frac_boundary=ez.copy(),
        rna_neg_frac_boundary=ez.copy(),
        gdna_density_global=0.01,
        gdna_reference_density=reference,
        gdna_reference_members=0 if reference is None else 1,
        gdna_capture_efficiency_region=np.asarray(efficiency, dtype=np.float64),
        gdna_capture_efficiency_boundary=np.ones(ne),
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


# --- the geometry: every transcript, every piece, and the taper partitions the start count -----------


@pytest.mark.parametrize("partition", ["coarse", "live"])
def test_every_transcript_maps_to_the_pieces_its_exons_overlap(misaligned_index, partition):
    """No transcript (mature or nascent entity) is dropped, and the piece covering the exon start that
    sits INTERIOR to a merged region is mapped — the off-by-one that skipped it on the coarse partition
    and the live partition, where that start is a region interface, both."""
    idx = misaligned_index
    ra = _coarsened(idx) if partition == "coarse" else RegionArrays.from_index(idx)
    t, piece, _ltau = transcript_piece_lengths(idx, ra, PMF)
    assert set(int(x) for x in t) == set(range(len(idx.t_df))), "a transcript was dropped"
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    covering = int(np.flatnonzero((starts <= 150) & (ends > 150))[0])
    tdf = idx.t_df
    for ti in range(len(tdf)):
        a, b = int(tdf["start"].iloc[ti]), int(tdf["end"].iloc[ti])
        if a <= 150 < b:
            assert covering in set(int(p) for p in piece[t == ti])


@pytest.mark.parametrize("partition", ["coarse", "live"])
def test_a_transcripts_taper_weighted_pieces_sum_to_its_fl_marginal_length(
    misaligned_index, partition
):
    """Σ_p ℓ_p^τ over a transcript's pieces is Σ_w f(w)(L − w + 1)⁺ — on the coarse partition, whose
    pieces spill beyond the exons and must be clipped to them, and on the live one."""
    idx = misaligned_index
    ra = _coarsened(idx) if partition == "coarse" else RegionArrays.from_index(idx)
    t, _piece, ltau = transcript_piece_lengths(idx, ra, PMF)
    total = np.zeros(len(idx.t_df))
    np.add.at(total, t, ltau)
    np.testing.assert_allclose(total, _fl(idx.t_df["length"].to_numpy()), rtol=1e-9)


def test_the_pieces_are_on_the_transcripts_own_reference(two_ref_index):
    """A second reference numbers its regions after the first's: every piece of a transcript must lie on
    the transcript's reference, or a region index would silently read another chromosome's efficiency."""
    idx = two_ref_index
    ra = RegionArrays.from_index(idx)
    t, piece, _ = transcript_piece_lengths(idx, ra, PMF)
    ref_of_t = idx.t_df["ref"].astype(str).map(idx.ref_name_to_id).to_numpy()
    np.testing.assert_array_equal(np.asarray(ra.ref_id)[piece], ref_of_t[t])


# --- the contraction: efficiencies 1 return fl, a field contracts, a tiny exon carries its bases -----


def test_no_reference_returns_the_fl_marginal_lengths_bit_identically(multiexon_index):
    """Capture off, or no gDNA: every efficiency is 1 and the ruler is ``fl`` EXACTLY, not within float
    noise — a contraction from rounding is a systematic bias, not a tolerance."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    fl = np.linspace(900.0, 2000.0, len(idx.t_df))
    eff = transcript_capture_eff_lengths(_cal(ra, np.ones(ra.n_regions), None), ra, idx, fl, PMF)
    np.testing.assert_array_equal(eff, fl)


def test_efficiencies_of_one_everywhere_under_a_reference_contract_nothing(multiexon_index):
    """With a reference and every piece at it, the factor is 1 to floating point on every transcript,
    spliced and unspliced alike: the taper's numerator and denominator are the same sum."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, np.ones(ra.n_regions), 1.0), ra, idx, fl, PMF)
    np.testing.assert_allclose(eff, fl, rtol=1e-12)


def test_a_field_contracts_and_never_expands(multiexon_index):
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500)] = 1.0  # the first exon captured, everything else depleted
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, PMF)
    assert np.all(eff <= fl + 1e-9)
    assert np.any(eff < fl - 1e-6)


def test_the_length_is_the_taper_weighted_mean_of_the_pieces_efficiencies(multiexon_index):
    """The claim, computed independently: the six-exon mRNA with exons 1–3 at efficiency 1 and 4–6 at
    0.2 reads ``Σ_p ℓ_p^τ c̃_p / Σ_p ℓ_p^τ`` with ``ℓ^τ`` from the taper on the 3,000 bp cDNA."""
    from rigel.calibration.effective_length import base_taper

    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    for k, s in enumerate(range(1000, 6500, 1000)):
        c[_exon_mask(ra, s, s + 500)] = 1.0 if k < 3 else 0.2
    m = _tidx(idx, "mrna")
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, PMF)
    tap = base_taper(PMF)
    L = int(idx.t_df["length"].iloc[m])
    first_half = tap.interval_sums(np.array([0]), np.array([1500]), L)[0]
    second_half = tap.interval_sums(np.array([1500]), np.array([3000]), L)[0]
    expected = fl[m] * (first_half * 1.0 + second_half * 0.2) / (first_half + second_half)
    assert eff[m] == pytest.approx(expected, rel=1e-9)


def test_a_transcript_of_exons_shorter_than_a_fragment_carries_its_bases_weight(tiny_index):
    """Ten 40 bp exons against a 100–300 bp pmf: no piece can contain a fragment, so an object-set
    ruler on contained supports read this transcript as 0 and then its floor. Its length is its bases:
    at efficiency ½ on every exon the factor is exactly ½, and with half the exons at 1 and half at
    0.001 it is the taper-weighted base mean — the tapered ends carry less than the middle exons."""
    from rigel.calibration.effective_length import base_taper

    idx = tiny_index
    ra = RegionArrays.from_index(idx)
    t = _tidx(idx, "tiny")
    fl = _fl(idx.t_df["length"].to_numpy())
    exons = [_exon_mask(ra, 1000 + i * 1040, 1040 + i * 1040) for i in range(10)]
    half = np.full(ra.n_regions, 0.001)
    for m in exons:
        half[m] = 0.5
    eff = transcript_capture_eff_lengths(_cal(ra, half, 1.0), ra, idx, fl, PMF)
    assert eff[t] == pytest.approx(0.5 * fl[t], rel=1e-12)
    mixed = np.full(ra.n_regions, 0.001)
    for i, m in enumerate(exons):
        mixed[m] = 1.0 if i < 5 else 0.001
    eff = transcript_capture_eff_lengths(_cal(ra, mixed, 1.0), ra, idx, fl, PMF)
    tap = base_taper(PMF)
    L = 400
    w_first = tap.interval_sums(np.array([0]), np.array([200]), L)[0]
    w_last = tap.interval_sums(np.array([200]), np.array([400]), L)[0]
    expected = fl[t] * (w_first * 1.0 + w_last * 0.001) / (w_first + w_last)
    assert eff[t] == pytest.approx(expected, rel=1e-9)
    assert 0.3 * fl[t] < eff[t] < 0.7 * fl[t]


def test_the_ruler_reads_the_efficiencies_and_nothing_from_the_counts(multiexon_index):
    """Two results with the same efficiencies and wildly different counts and supports give the same
    lengths: the evidence was weighed upstream (`capture_efficiency`), and the ruler is geometry."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.3)
    c[_exon_mask(ra, 1000, 1500)] = 1.0
    fl = _fl(idx.t_df["length"].to_numpy())
    cal = _cal(ra, c, 1.0)
    loud = dataclasses.replace(
        cal,
        count_gdna_region=np.full(ra.n_regions, 1e6),
        gdna_region_eff_len=np.full(ra.n_regions, 1e-9),
        count_gdna_boundary=np.full(cal.n_boundaries, 1e6),
    )
    np.testing.assert_array_equal(
        transcript_capture_eff_lengths(cal, ra, idx, fl, PMF),
        transcript_capture_eff_lengths(loud, ra, idx, fl, PMF),
    )


def test_no_nascent_mature_inversion_under_capture(multiexon_index):
    """A nascent parent's bases contain its spliced child's, and its taper at every exonic base is at
    least the child's (a base interior to the span is interior to the child at most), so
    ``eff(nascent) ≥ eff(mature)`` for any field; and the mature genuinely contracts here."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    c = np.full(ra.n_regions, 0.001)
    c[_exon_mask(ra, 1000, 1500)] = 1.0
    fl = _fl(idx.t_df["length"].to_numpy())
    eff = transcript_capture_eff_lengths(_cal(ra, c, 1.0), ra, idx, fl, PMF)
    mrna, nasc = _tidx(idx, "mrna"), _tidx(idx, "nasc")
    assert eff[nasc] >= eff[mrna] - 1e-9
    assert eff[mrna] < fl[mrna] - 1e-6
