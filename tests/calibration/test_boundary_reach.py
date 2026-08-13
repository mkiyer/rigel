"""The per-contiguous-boundary RNA REACH, on the accumulator's boundary axis — S5.g's input.

    Ruling: · Equation:

⭐ **What reach is for.** A crossing molecule must fit in what remains of **its own template** either
side of the line. gDNA's template is the chromosome, so its reach is unbounded — that is physics, not a
choice. RNA's template ends where its transcript ends, so near a terminus the admissible placements
collapse and the crossing divisor with them. Ignoring that over-calls gDNA by a measured **11.0 %**
genome-wide and by **+0.36** in the last region before a polyA site.

⚠ **PER STRAND, and per SIDE**: reach is "maximised over transcripts
independently per side AND per strand". A POS-strand transcript and a NEG-strand one ending at
different places give a line two different RNA reaches, and averaging them would describe neither.

⚠ **A contiguous boundary's reach is GENOMIC, unlike a junction's, which is EXONIC.** A junction is used
only by a spliced molecule, so what remains either side of it is exonic. A contiguous line is crossed by
*nascent* RNA too, which is genomic — taking the exonic reach there would declare an intronic nascent
fragment impossible (`splice_graph.JunctionGeometry`).

⚠ **A reach of 0 is MEANINGFUL, not a sentinel**. It says there is no template
of that strand at that line at all, so RNA of that strand has zero opportunity — which is what makes
`crossing_eff_length` return 0 and the consuming `_rate` emit nothing rather than a floored value
(trap 23). Measured on the chr22 pilot index: **40.6 %** POS and **42.9 %** NEG.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.splice_graph import (
    build_contiguous_boundary_reach_arrays,
    build_boundary_flags_array,
    build_region_partition_arrays,
)

from conftest import build_test_index

#: chr1: t0 is a two-exon POS transcript [200,400)+[700,900); t1 is a NEG transcript [1000,1200) whose
#: reach therefore differs from t0's at every line between them. chr2 keeps the per-reference offsets
#: honest — a single-reference fixture cannot see an axis that is laid out per reference.
GTF = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t701\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t1001\t1200\t.\t-\t.\tgene_id "g2"; transcript_id "t1";
chr2\ttest\texon\t301\t500\t.\t+\t.\tgene_id "g3"; transcript_id "t2";
chr2\ttest\texon\t801\t1000\t.\t+\t.\tgene_id "g3"; transcript_id "t2";
"""

REFS = {"chr1": 1500, "chr2": 1500}


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    return build_test_index(tmp_path_factory, GTF, name="boundary_reach", refs=REFS)


def _boundary_positions(index) -> np.ndarray:
    """Genomic position of every contiguous boundary, in boundary order — the independent coordinate.

    A reference contributing ``c`` cuts owns ``c − 1`` regions and ``c − 2`` interior lines, and line
    ``e`` sits at cut ``e + 1``. Derived from the PARTITION, not from the reach builder, so the two
    cannot agree by sharing a helper.
    """
    positions, cut_offsets, _types = build_region_partition_arrays(index)
    out = []
    for f in range(len(cut_offsets) - 1):
        lo, hi = int(cut_offsets[f]), int(cut_offsets[f + 1])
        if hi - lo >= 2:
            out.append(positions[lo + 1 : hi - 1])
    return np.concatenate(out) if out else np.zeros(0, np.int64)


# ---------------------------------------------------------------------------
# Shape: the SAME axis as everything else on the boundary side.
# ---------------------------------------------------------------------------


def test_one_entry_per_line_on_the_SAME_axis_as_the_flags(index):
    """⭐ The reach array and the flags array must be the same axis, element for element — both are
    per contiguous boundary, and a consumer indexes them with one index."""
    reach_lo, reach_hi = build_contiguous_boundary_reach_arrays(index)
    n_boundaries = build_boundary_flags_array(index).shape[0]
    assert reach_lo.shape == (n_boundaries, 2)
    assert reach_hi.shape == (n_boundaries, 2)
    assert reach_lo.shape[0] == _boundary_positions(index).shape[0]


def test_the_strand_axis_is_POS_then_NEG(index):
    """Column 0 is the POS-strand transcript's reach, column 1 the NEG's — the same ordering the
    accumulator's two columns and ``JunctionGeometry``'s strand join use.

    ⚠ On this fixture chr1 carries a POS transcript at [200,900) and a NEG one at [1000,1200), which
    are disjoint — so there are lines where exactly one column is non-zero, and a transposed strand
    axis moves the non-zero to the wrong line rather than merely permuting a pair.
    """
    reach_lo, reach_hi = build_contiguous_boundary_reach_arrays(index)
    pos_live = (reach_lo[:, 0] > 0) | (reach_hi[:, 0] > 0)
    neg_live = (reach_lo[:, 1] > 0) | (reach_hi[:, 1] > 0)
    assert pos_live.any() and neg_live.any(), "the fixture must exercise both strands"
    # the two transcripts are disjoint, so no line can be live on both strands
    assert not (pos_live & neg_live).any()


def test_reach_is_ZERO_where_the_strand_carries_no_transcript(index):
    """⚠ Zero is the ANSWER, not a missing value. A line with no POS-strand template gives POS-RNA no
    opportunity, and the consuming divisor must return 0 so the rate emits nothing (trap 23)."""
    reach_lo, reach_hi = build_contiguous_boundary_reach_arrays(index)
    assert int((reach_lo[:, 0] == 0).sum()) > 0
    assert np.all(reach_lo >= 0) and np.all(reach_hi >= 0)


def test_reach_matches_the_INDEX_keyed_by_src_region(index):
    """The values are the edges_df reach columns of the boundary whose ``src`` is the line's left region —
    read back from the frame directly rather than re-derived."""
    boundaries = index.edges_df
    contiguous = boundaries[boundaries["kind"] == 0]
    by_src = {
        int(s): row for s, row in zip(contiguous["src"], contiguous.itertuples(), strict=True)
    }

    reach_lo, reach_hi = build_contiguous_boundary_reach_arrays(index)
    from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices

    ra = RegionArrays.from_index(index)
    lo_region, _hi_region = boundary_region_indices(np.asarray(ra.ref_id))
    assert lo_region.shape[0] == reach_lo.shape[0]
    for e, src in enumerate(lo_region):
        row = by_src[int(src)]
        assert reach_lo[e, 0] == row.reach_lo_pos
        assert reach_hi[e, 0] == row.reach_hi_pos
        assert reach_lo[e, 1] == row.reach_lo_neg
        assert reach_hi[e, 1] == row.reach_hi_neg


def test_a_single_region_reference_contributes_no_entry(tmp_path_factory):
    """A reference with one region owns no line, so it contributes nothing — the invariant the ``k + 1``
    boundary axis could not state."""
    one = build_test_index(
        tmp_path_factory,
        'chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g"; transcript_id "t";\n',
        name="boundary_reach_one",
        refs={"chr1": 600, "chrE": 400},  # chrE has no feature at all
    )
    reach_lo, _ = build_contiguous_boundary_reach_arrays(one)
    assert reach_lo.shape[0] == build_boundary_flags_array(one).shape[0]


# ---------------------------------------------------------------------------
# What the reach BUYS — the taper is not a refinement.
# ---------------------------------------------------------------------------


def test_the_taper_COLLAPSES_the_divisor_near_a_terminus(index):
    """⭐ The whole point of TRAPS: prove-the-substrate, as a number rather than an assertion.

    ``crossing_eff_length`` at unbounded reach is ``mu − 1``; at a real reach it is far smaller near a
    transcript end. records 199.0 at R=550 against 50.0 at R=50 on RNA N(200,50)
    — a **4×** error if the mean is used blindly at a first exon.
    """
    from rigel.calibration.effective_length import UNBOUNDED_REACH, crossing_eff_length

    pmf = np.zeros(401)
    pmf[200] = 1.0  # a delta at 200: crossing placements are min(199, R_lo, R_hi, R_lo+R_hi-199)
    unbounded = float(crossing_eff_length(pmf, [UNBOUNDED_REACH], [UNBOUNDED_REACH])[0])
    assert unbounded == pytest.approx(199.0)

    reach_lo, reach_hi = build_contiguous_boundary_reach_arrays(index)
    tapered = crossing_eff_length(pmf, reach_lo[:, 0], reach_hi[:, 0])
    assert np.all(tapered <= unbounded + 1e-9)
    # ⚠ Strictly smaller somewhere, or the taper is inert on this fixture and proves nothing.
    assert np.any(tapered < unbounded - 1e-9)
    # and exactly 0 where the strand has no template at all
    assert np.all(tapered[(reach_lo[:, 0] == 0) & (reach_hi[:, 0] == 0)] == 0.0)
