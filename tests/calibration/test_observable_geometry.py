"""The observable genome over the partition (`region_arrays`), one base at a time, and the fixed point every reading
over it keeps.

A start counts only where its fragment touches no base of the index's unobservable mask (`rigel.observable`). The
geometry is each region cut at the mask into pieces, each boundary's bases to the nearest masked base on either side,
and each reference's maximal observable runs, gated here against a walk over the bases; the closed forms read over it
are gated against every placement beside their plain forms (``test_effective_length.py``, ``test_gdna_opportunity.py``).
With nothing masked each reading is its plain form's number to the bit, so an index with an empty mask — the test
chromosome, every synthetic genome — moves no number.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pysam
import pytest
from _index_builder import MINI_GTF, build_test_index

from rigel.calibration.calibrate import _gdna_boundary_conserved_len
from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    contained_eff_length,
    crossing_eff_length,
    observable_contained_eff_length,
    observable_crossing_eff_length,
    observable_cut_shares,
)
from rigel.calibration.gdna_opportunity import (
    contained_opportunity,
    crossing_opportunity,
    observable_gdna_opportunity,
    total_opportunity,
)
from rigel.calibration.signature import RegionType, coarse_type_array
from rigel.calibration.region_arrays import (
    RegionArrays,
    boundary_reaches,
    boundary_region_indices,
    observable_length,
    observable_pieces,
    observable_stretches,
)
from rigel.index import TranscriptIndex

#: Two references cut into regions, and a mask over them: a run at a reference's start, one inside a region, one across a
#: boundary, one filling a region (so touching a boundary from each side), one touching a boundary from its left, two
#: three bases apart, and a reference with none.
REFS = {
    "chr1": [0, 120, 130, 160, 175, 240, 300, 420],
    "chr2": [0, 100, 150, 260],
}
MASK = {
    "chr1": [(0, 6), (125, 128), (155, 165), (175, 240), (290, 300), (330, 333), (336, 340)],
    "chr2": [],
}


def _frame(refs=REFS):
    rows = []
    for ref, bounds in refs.items():
        for s, e in zip(bounds[:-1], bounds[1:]):
            rows.append((ref, s, e, 0))
    return pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])


def _mask_frame(mask=MASK):
    rows = [(ref, s, e) for ref, runs in mask.items() for s, e in runs]
    return pd.DataFrame(rows, columns=["ref", "start", "end"])


def _arrays(refs=REFS, mask=MASK):
    ids = {ref: i for i, ref in enumerate(refs)}
    return RegionArrays.from_frame(_frame(refs), ids, None if mask is None else _mask_frame(mask))


def _masked(ref, mask=MASK):
    flags = np.zeros(REFS[ref][-1], dtype=bool)
    for s, e in mask.get(ref, []):
        flags[s:e] = True
    return flags


def _runs(flags):
    runs, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        if not f and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(flags)))
    return runs


# ── the geometry ──────────────────────────────────────────────────────────────────────────────────


def test_the_pieces_are_the_regions_cut_at_the_mask():
    ra = _arrays()
    piece_region, piece_len = observable_pieces(ra)
    expected_region, expected_len = [], []
    r = 0
    for ref, bounds in REFS.items():
        masked = _masked(ref)
        for s, e in zip(bounds[:-1], bounds[1:]):
            for a, b in _runs(~masked[s:e]):
                expected_region.append(r)
                expected_len.append(b - a)
            r += 1
    assert piece_region.tolist() == expected_region
    assert piece_len.tolist() == expected_len
    # the region [175, 240) is masked whole and has no piece; the observable bases are the rest
    assert 4 not in piece_region.tolist()
    lengths = np.concatenate([np.diff(bounds) for bounds in REFS.values()])
    np.testing.assert_array_equal(
        observable_length(ra),
        [
            sum(1 for p in range(s, e) if not _masked(ref)[p])
            for ref, bounds in REFS.items()
            for s, e in zip(bounds[:-1], bounds[1:])
        ],
    )
    assert observable_length(ra).sum() < lengths.sum()


def test_the_reaches_are_the_bases_to_the_nearest_masked_base():
    ra = _arrays()
    reach_lo, reach_hi = boundary_reaches(ra)
    expected_lo, expected_hi = [], []
    for ref, bounds in REFS.items():
        masked = _masked(ref)
        for x in bounds[1:-1]:
            left = [x - 1 - p for p in range(x - 1, -1, -1) if masked[p]]
            right = [p - x for p in range(x, len(masked)) if masked[p]]
            expected_lo.append(float(left[0]) if left else np.inf)
            expected_hi.append(float(right[0]) if right else np.inf)
    np.testing.assert_array_equal(reach_lo, expected_lo)
    np.testing.assert_array_equal(reach_hi, expected_hi)


def test_the_stretches_are_the_references_cut_at_the_mask():
    got = observable_stretches(_arrays()).tolist()
    expected = [b - a for ref in REFS for a, b in _runs(~_masked(ref))]
    assert got == expected


def test_with_nothing_masked_the_geometry_is_the_partition():
    ra = _arrays(mask=None)
    piece_region, piece_len = observable_pieces(ra)
    np.testing.assert_array_equal(piece_region, np.arange(ra.n_regions))
    np.testing.assert_array_equal(piece_len, ra.end - ra.start)
    reach_lo, reach_hi = boundary_reaches(ra)
    assert np.all(np.isinf(reach_lo)) and np.all(np.isinf(reach_hi))
    assert observable_stretches(ra).tolist() == [b[-1] for b in REFS.values()]


# ── with nothing masked, today's numbers ──────────────────────────────────────────────────────────


def _pmf():
    """Uneven mass on widths 2..40, so a fragment reaches past short pieces and across several boundaries."""
    w = np.arange(41, dtype=np.float64)
    p = np.where(w >= 2.0, 1.0 + 0.5 * np.cos(w), 0.0)
    return p / p.sum()


_TWO_REFS = {"chr1": 2000, "chrB": 1500}


@pytest.fixture(scope="module")
def unmasked_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, MINI_GTF, name="unmasked", refs=_TWO_REFS)


def test_an_index_of_known_sequence_has_nothing_masked(unmasked_index):
    assert len(unmasked_index.unobservable_df) == 0
    assert RegionArrays.from_index(unmasked_index).unobservable.shape == (0, 3)


def test_with_nothing_masked_every_function_is_its_plain_forms_number_to_the_bit(unmasked_index):
    """The fixed point: on an index with nothing masked each observable reading returns exactly what its plain form
    returns, so the switch from one to the other moves no number anywhere a genome has no unobservable base."""
    ra = RegionArrays.from_index(unmasked_index)
    pmf = _pmf()
    np.testing.assert_array_equal(
        observable_contained_eff_length(ra, pmf), contained_eff_length(ra.region_size_bp, pmf)
    )
    n_boundaries = boundary_region_indices(ra.ref_id)[0].shape[0]
    unbounded = float(
        crossing_eff_length(pmf, np.full(1, UNBOUNDED_REACH), np.full(1, UNBOUNDED_REACH))[0]
    )
    np.testing.assert_array_equal(
        observable_crossing_eff_length(ra, pmf), np.full(n_boundaries, unbounded)
    )
    np.testing.assert_array_equal(
        observable_cut_shares(ra, pmf), _gdna_boundary_conserved_len(ra, pmf)
    )
    # the pools over whole regions: contained by each region's type, crossing by its flanks' sorted pair
    width = pmf.shape[0] - 1
    types = coarse_type_array(ra.signature)
    length = (ra.end - ra.start).astype(np.int64)
    lo, hi = boundary_region_indices(ra.ref_id)
    pairs = np.sort(np.stack([types[lo], types[hi]]), axis=0)
    intron_exon = (pairs[0] == RegionType.INTRON) & (pairs[1] == RegionType.EXON)
    intergenic_exon = (pairs[0] == RegionType.INTERGENIC) & (pairs[1] == RegionType.EXON)
    whole = (
        contained_opportunity(length[types == RegionType.INTERGENIC], width),
        contained_opportunity(length[types == RegionType.INTRON], width),
        crossing_opportunity(length[lo][intron_exon], length[hi][intron_exon], width),
        crossing_opportunity(length[lo][intergenic_exon], length[hi][intergenic_exon], width),
        total_opportunity(np.array([2000, 1500]), width),
    )
    observable = observable_gdna_opportunity(ra, width)
    for got, want in zip(observable.pools + (observable.total,), whole):
        np.testing.assert_array_equal(got, want)


def test_the_frame_carries_the_index_mask(tmp_path_factory):
    """An index of a genome with unknown runs hands its mask to the frame, in the frame's reference ids."""
    base = tmp_path_factory.mktemp("masked_frame")
    seq = ("ACGT" * 600)[:2000]
    seq = seq[:50] + "N" * 30 + seq[80:1990] + "N" * 10
    fasta = base / "genome.fa"
    fasta.write_text(">chrB\n" + "ACGT" * 25 + "\n>chr1\n" + seq + "\n")
    pysam.faidx(str(fasta))
    gtf = base / "test.gtf"
    gtf.write_text(MINI_GTF)
    TranscriptIndex.build(fasta, gtf, base / "idx", write_tsv=False)
    index = TranscriptIndex.load(base / "idx")
    ra = RegionArrays.from_index(index)
    chr1 = index.ref_name_to_id["chr1"]
    assert ra.unobservable.tolist() == [[chr1, 50, 80], [chr1, 1990, 2000]]
