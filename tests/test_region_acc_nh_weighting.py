"""Region accumulator NH-weighting correctness tests.

Verifies that multimapper fragments credit each of their ``NH`` hits
with weight ``1/NH`` instead of full weight (old "first-hit wins"
behaviour).  Summed across all hits of one molecule, a multimapper
contributes exactly the same total region weight (1.0) as a unique
mapper — matching the ``1/NH`` on-target crediting convention of
``mappable_effective_length`` so that the ratio estimator
``λ̂ = Σ aᵢ / Σ Eᵢ`` is unbiased.

See docs/calibration/KNOWN_BUGS.md BUG #2.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pysam
import pytest

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer


# ---------------------------------------------------------------------------
# BAM construction helpers (tailored to the mini_index fixture).
#
# mini_index has a 2000bp chr1 with these calibration regions:
#     id 1: [  99, 200)  — g1 exon 1 (tx_pos=True, exon_pos=True)
#     id 3: [ 299, 400)  — g1 exon 2
#     id 5: [ 499, 600)  — g1 exon 3
#     id 7: [ 999,1100)  — g2 exon 1 (tx_neg=True, exon_neg=True)
#     id 9: [1199,1300)  — g2 exon 2
# ---------------------------------------------------------------------------

_REF_NAME = "chr1"
_REF_LEN = 2000
_HEADER = {
    "HD": {"VN": "1.6", "SO": "queryname"},
    "SQ": [{"SN": _REF_NAME, "LN": _REF_LEN}],
}


def _make_pair(
    qname: str,
    frag_start: int,
    *,
    length: int = 50,
    frag_len: int = 150,
    is_secondary: bool = False,
    nh: int = 1,
) -> list[pysam.AlignedSegment]:
    """Create a properly paired FR R1+R2 pair on the + strand.

    Convention (dUTP/TruSeq stranded, sense from + transcript):
      R2 forward at 5' end of fragment (lower genomic pos);
      R1 reverse at 3' end of fragment (higher genomic pos).

    Fragment covers ``[frag_start, frag_start + frag_len)``.
    R2 covers [frag_start, frag_start+length).
    R1 covers [frag_start+frag_len-length, frag_start+frag_len).
    Both reads are unspliced (``length`` M CIGAR).
    """
    r2_pos = frag_start
    r1_pos = frag_start + frag_len - length
    seq = "A" * length
    qual = pysam.qualitystring_to_array("I" * length)

    # R1: reverse strand, at 3' end of fragment.
    r1 = pysam.AlignedSegment()
    r1.query_name = qname
    r1.reference_id = 0
    r1.reference_start = r1_pos
    r1.mapping_quality = 60
    flag_r1 = 0x1 | 0x2 | 0x40 | 0x10 | 0x20  # paired, proper, R1, rev, mate-rev? no
    # Actually: R1 reverse (0x10), mate (R2) not reverse so don't set 0x20.
    flag_r1 = 0x1 | 0x2 | 0x40 | 0x10
    if is_secondary:
        flag_r1 |= 0x100
    r1.flag = flag_r1
    r1.cigar = [(0, length)]  # M
    r1.query_sequence = seq
    r1.query_qualities = qual
    r1.next_reference_id = 0
    r1.next_reference_start = r2_pos
    r1.template_length = -frag_len
    r1.set_tag("NH", int(nh), value_type="i")

    # R2: forward strand, at 5' end of fragment.
    r2 = pysam.AlignedSegment()
    r2.query_name = qname
    r2.reference_id = 0
    r2.reference_start = r2_pos
    r2.mapping_quality = 60
    # R2 forward (no 0x10); mate (R1) is reverse → set 0x20.
    flag_r2 = 0x1 | 0x2 | 0x80 | 0x20
    if is_secondary:
        flag_r2 |= 0x100
    r2.flag = flag_r2
    r2.cigar = [(0, length)]
    r2.query_sequence = seq
    r2.query_qualities = qual
    r2.next_reference_id = 0
    r2.next_reference_start = r1_pos
    r2.template_length = frag_len
    r2.set_tag("NH", int(nh), value_type="i")

    return [r1, r2]


def _write_bam(path: Path, reads: list[pysam.AlignedSegment]) -> str:
    bam_path = str(path / "test.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=_HEADER) as out:
        for r in reads:
            out.write(r)
    sorted_path = str(path / "test.nsorted.bam")
    pysam.sort("-n", "-o", sorted_path, bam_path)
    return sorted_path


def _scan_counts(bam_path: str, index) -> np.ndarray:
    """Return the region_counts matrix [n_regions × 4] as float64."""
    scan_config = BamScanConfig(
        sj_strand_tag="auto",
        n_scan_threads=1,
        include_multimap=True,
    )
    _, _, _, _, region_counts, _ = scan_and_buffer(bam_path, index, scan_config)
    assert region_counts is not None
    arr = np.zeros((len(region_counts), 4), dtype=np.float64)
    arr[:, 0] = region_counts["n_unspliced_pos"].values
    arr[:, 1] = region_counts["n_unspliced_neg"].values
    arr[:, 2] = region_counts["n_spliced_pos"].values
    arr[:, 3] = region_counts["n_spliced_neg"].values
    return arr


# ---------------------------------------------------------------------------
# Unique-mapper baseline
# ---------------------------------------------------------------------------


# Accumulator is stored float32 until converted; tests tolerate that.
_TOL = 1e-5

# The _make_pair helper produces the standard dUTP R1-antisense
# paired-end layout (R1 reverse, R2 forward).  The C++ accumulator
# assigns this to the unspliced_neg column (col=1) because the
# fragment's genomic orientation is derived from R1's mapping.  We
# check col=1 throughout these tests — a unique fragment still
# contributes exactly 1.0 to its region regardless of column, and the
# NH-weighting invariant (Σ counts = # molecules) is column-agnostic.
_COL = 1  # unspliced_neg


class TestUniqueMapperBaseline:
    """A unique (NH=1) fragment contributes exactly 1.0 to its region."""

    def test_unique_unspliced_pos_contributes_one(self, mini_index, tmp_path):
        # Fragment [110, 110+80) = [110, 190), entirely inside region 1 [99,200).
        reads = _make_pair("u1", frag_start=110, length=50, frag_len=80, nh=1)
        bam = _write_bam(tmp_path, reads)

        counts = _scan_counts(bam, mini_index)

        total = counts.sum()
        assert total == pytest.approx(1.0, abs=_TOL), (
            f"A unique fragment should sum to exactly 1.0 across all regions, "
            f"got {total}"
        )
        assert counts[1, _COL] == pytest.approx(1.0, abs=_TOL)
        # Weight lands in exactly one column; other unspliced column and
        # both spliced columns must be zero.
        other_col = 0 if _COL == 1 else 1
        assert counts[:, other_col].sum() == 0.0
        assert counts[:, 2].sum() == 0.0
        assert counts[:, 3].sum() == 0.0


# ---------------------------------------------------------------------------
# Multimapper 1/NH weighting
# ---------------------------------------------------------------------------


class TestMultimapperNHWeighting:
    """NH=k multimapper contributes total weight 1.0 (1/k per hit)."""

    def test_nh3_multimapper_three_regions(self, mini_index, tmp_path):
        """3 hits in 3 distinct exonic regions (1, 3, 5) → each += 1/3."""
        qname = "mm3"
        reads = []
        reads += _make_pair(qname, frag_start=110, frag_len=80, nh=3, is_secondary=False)
        reads += _make_pair(qname, frag_start=310, frag_len=80, nh=3, is_secondary=True)
        reads += _make_pair(qname, frag_start=510, frag_len=80, nh=3, is_secondary=True)
        bam = _write_bam(tmp_path, reads)

        counts = _scan_counts(bam, mini_index)

        total = counts.sum()
        assert total == pytest.approx(1.0, abs=_TOL), (
            f"NH=3 multimapper should contribute 1.0 total across regions, "
            f"got {total}.  (Pre-fix behaviour would give 1.0 from first-hit "
            f"only OR 3.0 if the gate was removed without 1/NH weighting.)"
        )

        expected = 1.0 / 3.0
        for rid in (1, 3, 5):
            assert counts[rid, _COL] == pytest.approx(expected, abs=_TOL), (
                f"Region {rid} should be 1/3 = {expected}, "
                f"got {counts[rid, _COL]}"
            )

    def test_nh2_multimapper_two_regions(self, mini_index, tmp_path):
        """NH=2 multimapper → 0.5 per region."""
        qname = "mm2"
        reads = []
        reads += _make_pair(qname, frag_start=110, frag_len=80, nh=2, is_secondary=False)
        reads += _make_pair(qname, frag_start=510, frag_len=80, nh=2, is_secondary=True)
        bam = _write_bam(tmp_path, reads)

        counts = _scan_counts(bam, mini_index)

        total = counts.sum()
        assert total == pytest.approx(1.0, abs=_TOL)
        assert counts[1, _COL] == pytest.approx(0.5, abs=_TOL)
        assert counts[5, _COL] == pytest.approx(0.5, abs=_TOL)


# ---------------------------------------------------------------------------
# Mixed unique + multimapper
# ---------------------------------------------------------------------------


class TestMixedUniqueAndMultimapper:
    """A unique fragment plus an NH=3 multimapper → total weight 2.0."""

    def test_mixed_sum_matches_molecule_count(self, mini_index, tmp_path):
        reads = []
        reads += _make_pair("u1", frag_start=110, frag_len=80, nh=1)
        reads += _make_pair("mm3", frag_start=120, frag_len=70, nh=3, is_secondary=False)
        reads += _make_pair("mm3", frag_start=310, frag_len=80, nh=3, is_secondary=True)
        reads += _make_pair("mm3", frag_start=510, frag_len=80, nh=3, is_secondary=True)
        bam = _write_bam(tmp_path, reads)

        counts = _scan_counts(bam, mini_index)

        total = counts.sum()
        assert total == pytest.approx(2.0, abs=_TOL)

        # Region 1 should see unique(1.0) + mm(1/3) = 4/3.
        assert counts[1, _COL] == pytest.approx(1.0 + 1.0 / 3.0, abs=_TOL)
        assert counts[3, _COL] == pytest.approx(1.0 / 3.0, abs=_TOL)
        assert counts[5, _COL] == pytest.approx(1.0 / 3.0, abs=_TOL)


# ---------------------------------------------------------------------------
# Intergenic multimapper (hits that don't resolve to transcripts)
# ---------------------------------------------------------------------------


class TestIntergenicNHWeighting:
    """NH=2 multimapper both hits intergenic → 1/2 each, total 1.0."""

    def test_nh2_both_intergenic(self, mini_index, tmp_path):
        qname = "mm_ig"
        reads = []
        reads += _make_pair(qname, frag_start=610, frag_len=80, nh=2, is_secondary=False)
        reads += _make_pair(qname, frag_start=1310, frag_len=80, nh=2, is_secondary=True)
        bam = _write_bam(tmp_path, reads)

        counts = _scan_counts(bam, mini_index)

        total = counts.sum()
        assert total == pytest.approx(1.0, abs=_TOL), (
            f"Intergenic NH=2 multimapper must still sum to 1.0; got {total}"
        )
        assert counts[6, _COL] == pytest.approx(0.5, abs=_TOL)
        assert counts[10, _COL] == pytest.approx(0.5, abs=_TOL)
