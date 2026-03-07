"""Unit tests for fragment-weight and transcript-position scoring functions.

Tests ``genomic_to_transcript_pos`` and ``compute_fragment_weight``
from ``rigel.scoring``.
"""

import numpy as np
import pytest

from rigel.scoring import (
    compute_fragment_weight,
    genomic_to_transcript_pos,
)
from rigel.types import Strand

POS = int(Strand.POS)
NEG = int(Strand.NEG)


# =====================================================================
# genomic_to_transcript_pos
# =====================================================================


def _exons(intervals):
    """Create an (n, 2) int32 array from a list of (start, end) tuples."""
    return np.array(intervals, dtype=np.int32)


class TestGenomicToTranscriptPos:
    """Tests for genomic_to_transcript_pos."""

    def test_single_exon_start(self):
        """Position at exon start → offset 0 (POS strand)."""
        exons = _exons([(100, 200)])
        assert genomic_to_transcript_pos(100, exons, POS, 100) == 0

    def test_single_exon_middle(self):
        """Position inside a single exon."""
        exons = _exons([(100, 200)])
        assert genomic_to_transcript_pos(150, exons, POS, 100) == 50

    def test_single_exon_end(self):
        """Position at exon end → maps to transcript_length."""
        exons = _exons([(100, 200)])
        assert genomic_to_transcript_pos(200, exons, POS, 100) == 100

    def test_multi_exon_second_exon(self):
        """Position in second exon of a two-exon transcript."""
        # Exon1: [100, 150) len=50, Exon2: [300, 400) len=100, total=150
        exons = _exons([(100, 150), (300, 400)])
        # Position 310 is in exon2: offset = 50 + (310 - 300) = 60
        assert genomic_to_transcript_pos(310, exons, POS, 150) == 60

    def test_intron_maps_to_boundary(self):
        """Position in intron maps to preceding exon boundary."""
        exons = _exons([(100, 150), (300, 400)])
        # Position 200 is in intron: breaks at exon2 start (200 < 300)
        # offset from exon1 = 50 (fully consumed)
        assert genomic_to_transcript_pos(200, exons, POS, 150) == 50

    def test_before_first_exon(self):
        """Position before transcript start → offset 0."""
        exons = _exons([(100, 200)])
        assert genomic_to_transcript_pos(50, exons, POS, 100) == 0

    def test_past_last_exon(self):
        """Position past last exon → transcript_length."""
        exons = _exons([(100, 200)])
        assert genomic_to_transcript_pos(500, exons, POS, 100) == 100

    def test_neg_strand_flips(self):
        """NEG strand flips offset: pos 100 (offset 0 on POS) → L on NEG."""
        exons = _exons([(100, 200)])
        assert genomic_to_transcript_pos(100, exons, NEG, 100) == 100

    def test_neg_strand_middle(self):
        """NEG strand middle position."""
        exons = _exons([(100, 200)])
        # POS: offset = 50 → NEG: 100 - 50 = 50
        assert genomic_to_transcript_pos(150, exons, NEG, 100) == 50

    def test_neg_strand_end(self):
        """NEG strand: genomic end → transcript offset 0."""
        exons = _exons([(100, 200)])
        assert genomic_to_transcript_pos(200, exons, NEG, 100) == 0

    def test_three_exons(self):
        """Three-exon transcript: position in third exon."""
        # Exon1: [0, 100) len=100, Exon2: [200, 250) len=50,
        # Exon3: [400, 500) len=100, total=250
        exons = _exons([(0, 100), (200, 250), (400, 500)])
        # Position 420 in exon3: offset = 100 + 50 + (420 - 400) = 170
        assert genomic_to_transcript_pos(420, exons, POS, 250) == 170


# =====================================================================
# compute_fragment_weight
# =====================================================================


class TestComputeFragmentWeight:
    """Tests for compute_fragment_weight (C++ kernel)."""

    def test_plateau_returns_one(self):
        """Fragment entirely in plateau region → weight 1.0."""
        # L=1000, f=100, w=100. Plateau is [100, 900).
        # Fragment at [200, 300) is entirely in plateau.
        assert compute_fragment_weight(200, 300, 1000) == pytest.approx(1.0)

    def test_edge_returns_gt_one(self):
        """Fragment near 5' edge → weight > 1.0."""
        # L=1000, f=50, w=50. Left ramp is [0, 50).
        # Fragment at [0, 50) is entirely in left ramp.
        w = compute_fragment_weight(0, 50, 1000)
        assert w > 1.0

    def test_right_edge_returns_gt_one(self):
        """Fragment near 3' edge → weight > 1.0."""
        # L=1000, f=50, w=50. Right ramp is [950, 1000].
        w = compute_fragment_weight(950, 1000, 1000)
        assert w > 1.0

    def test_symmetry(self):
        """Left and right edges should give same weight for same distance."""
        L = 1000
        w_left = compute_fragment_weight(0, 100, L)
        w_right = compute_fragment_weight(900, 1000, L)
        assert w_left == pytest.approx(w_right, rel=1e-10)

    def test_full_span_fragment(self):
        """Fragment spans entire transcript."""
        # f = L = 500, w = min(500, 250) = 250
        w = compute_fragment_weight(0, 500, 500)
        assert w >= 1.0

    def test_short_transcript(self):
        """Fragment equals short transcript → weight >= 1.0."""
        w = compute_fragment_weight(0, 100, 100)
        assert w >= 1.0

    def test_zero_length_fragment(self):
        """Zero-length fragment → fallback 1.0."""
        assert compute_fragment_weight(50, 50, 1000) == 1.0

    def test_zero_transcript_length(self):
        """Zero transcript length → fallback 1.0."""
        assert compute_fragment_weight(0, 100, 0) == 1.0

    def test_negative_fragment(self):
        """Inverted coordinates → fallback 1.0."""
        assert compute_fragment_weight(100, 50, 1000) == 1.0

    def test_weight_increases_toward_edge(self):
        """Weight should increase as fragment moves from plateau to edge."""
        L = 2000
        w_center = compute_fragment_weight(500, 600, L)
        w_near_edge = compute_fragment_weight(10, 110, L)
        assert w_center < w_near_edge

    def test_half_span_fragment(self):
        """Fragment spanning exactly half the transcript."""
        # f = L/2 = 500, w = 500. Plateau is [500, 500) — empty.
        # Entire [0, 500) is in the left ramp.
        w = compute_fragment_weight(0, 500, 1000)
        assert w >= 1.0


class TestNativeScoringImport:
    """Verify the native C++ scoring module loads and matches Python."""

    def test_native_module_available(self):
        """_scoring_impl should be importable after build."""
        from rigel._scoring_impl import compute_fragment_weight as native_fw
        assert callable(native_fw)

    @pytest.mark.parametrize("fs,fe,tl", [
        (200, 300, 1000),   # plateau
        (0, 50, 1000),      # left edge
        (950, 1000, 1000),  # right edge
        (0, 500, 500),      # full span
        (0, 100, 100),      # short transcript
        (50, 50, 1000),     # zero-length
        (0, 100, 0),        # zero transcript
        (100, 50, 1000),    # inverted
        (0, 500, 1000),     # half span
        (10, 110, 2000),    # near edge
    ])
    def test_native_matches_python(self, fs, fe, tl):
        """C++ kernel must match Python fallback for all cases."""
        from rigel._scoring_impl import compute_fragment_weight as native_fw
        py_result = compute_fragment_weight(fs, fe, tl)
        native_result = native_fw(fs, fe, tl)
        assert native_result == pytest.approx(py_result, rel=1e-12)
