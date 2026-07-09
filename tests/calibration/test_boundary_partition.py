"""Boundary-partition geometry: one boundary row per region interface (the calibration boundary nodes).

For a reference with ``k`` regions there are ``k+1`` boundaries at positions ``[r0.start, r0.end, r1.end, …,
r(k-1).end]``, ref-major in genomic order — aligning by construction with the accumulator's boundary indexing.
(The per-boundary annotation flags this partition used to carry were removed in the message-precision collapse;
the solver reads junction strand from the accumulator splice motif — see INDEX_FORMAT_VERSION 7.)
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.regions import build_boundary_partition, build_region_partition
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _tx(t_id: str, strand: Strand, exons: list[tuple[int, int]]) -> Transcript:
    return Transcript(ref="chr1", strand=strand, exons=[Interval(s, e) for s, e in exons], t_id=t_id)


def test_boundary_count_matches_regions_plus_one():
    ta = _tx("TA", Strand.POS, [(1000, 5000), (20000, 24000)])
    tb = _tx("TB", Strand.POS, [(1000, 9000), (16000, 24000)])
    ref_lengths = {"chr1": 30000}
    rdf = build_region_partition([ta, tb], ref_lengths)
    bdf = build_boundary_partition(rdf, ref_lengths)
    assert len(bdf) == len(rdf) + 1
    # positions are exactly the region interfaces [start…, last end]
    starts = rdf["start"].to_numpy(np.int64)
    ends = rdf["end"].to_numpy(np.int64)
    expected = np.append(starts, ends[-1])
    assert np.array_equal(bdf["position"].to_numpy(np.int64), expected)
