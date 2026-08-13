"""``build_transcript_path`` — a transcript as an ordered walk over REGION / BOUNDARY / SPLICE JUNCTION.

⭐⭐ **THE INCLUDE/EXCLUDE RULE IS THE WHOLE CONTENT OF THE FUNCTION**, so every gate here is one clause
of it, stated as a property rather than as a recorded output:

* a multi-exon transcript takes its EXONIC regions and not its intronic ones;
* a single-exon transcript and a synthetic span take every region they cover, introns included;
* a boundary INTERIOR to an exon is crossed contiguously and is included — even when it exists only
  because a signature changed there;
* a transcript's OUTER boundaries are excluded, because the molecule ends at them;
* a splice donor/acceptor position appears as a SPLICE JUNCTION step, never as a BOUNDARY step;
* ⭐ the steps run in TRANSCRIPTION order, so a minus-strand path descends in genomic coordinate.

⛔ **The sj join goes through INTRON COORDINATES and never through the flanking region pair.** The
pair is unique on the shipped partition only because every exon endpoint is forced to be a region bound;
on a coarsened partition it collides, and it carries no strand.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.splice_graph import (
    STEP_BOUNDARY,
    STEP_REGION,
    STEP_SPLICE_SJ,
    build_transcript_path,
)
from rigel.types import IntervalType, Strand


class _Regions:
    """One reference's regions, tiling it — the shape ``RegionArrays`` presents."""

    def __init__(self, bounds):
        self.start = np.asarray(bounds[:-1], dtype=np.int64)
        self.end = np.asarray(bounds[1:], dtype=np.int64)
        self.n_regions = self.start.size
        self.n_refs = 1
        self.ref_offsets = np.array([0, self.n_regions], dtype=np.int64)
        self.ref_id = np.zeros(self.n_regions, dtype=np.int64)


class _Index:
    """The surface ``build_transcript_path`` reads: intervals.feather, t_df, regions_df, edges_df."""

    def __init__(self, tmp_path, bounds, exons_by_t, strands, synthetic=(), spans=None):
        self.index_dir = str(tmp_path)
        self.ref_name_to_id = {"chr1": 0}
        self.ref_names = ["chr1"]
        n_t = len(strands)
        self.num_transcripts = n_t

        rows = [
            {
                "t_index": int(t),
                "ref": "chr1",
                "start": int(a),
                "end": int(b),
                "interval_type": int(IntervalType.EXON),
            }
            for t, exons in exons_by_t.items()
            for a, b in exons
        ]
        pd.DataFrame(rows).to_feather(f"{tmp_path}/intervals.feather")

        spans = spans or {}
        self.t_df = pd.DataFrame(
            {
                "t_index": np.arange(n_t, dtype=np.int64),
                "ref": ["chr1"] * n_t,
                "start": [spans.get(i, (0, 0))[0] for i in range(n_t)],
                "end": [spans.get(i, (0, 0))[1] for i in range(n_t)],
                "strand": np.asarray(strands, dtype=np.int64),
                "is_synthetic": [i in set(synthetic) for i in range(n_t)],
            }
        )

        n_r = len(bounds) - 1
        self.regions_df = pd.DataFrame(
            {
                "node_id": np.arange(n_r, dtype=np.int64),
                "ref_name": ["chr1"] * n_r,
                "start": np.asarray(bounds[:-1], dtype=np.int64),
                "end": np.asarray(bounds[1:], dtype=np.int64),
            }
        )
        # sj boundaries: one per adjacent exon pair, keyed by the regions the intron separates
        jr = []
        for t, exons in exons_by_t.items():
            for (_, e0), (s1, _) in zip(exons[:-1], exons[1:], strict=True):
                src = int(np.flatnonzero(self.regions_df["end"].to_numpy() == e0)[0])
                dst = int(np.flatnonzero(self.regions_df["start"].to_numpy() == s1)[0])
                jr.append((src, dst, int(strands[t])))
        jr = sorted(set(jr))
        self.edges_df = pd.DataFrame(
            {
                "edge_id": np.arange(len(jr), dtype=np.int64),
                "src": [r[0] for r in jr],
                "dst": [r[1] for r in jr],
                "kind": np.ones(len(jr), dtype=np.uint8),  # EDGE_KIND_SJ
                "strand": np.asarray([r[2] for r in jr], dtype=np.int8),
            }
        )


@pytest.fixture
def _patched_sj(monkeypatch):
    """The fixture index carries sj boundaries directly, so the CSR builder is stubbed to the
    identity — this file gates the PATH, not `build_sj_arrays` (which has its own gates and
    was verified 13,482/13,482 against `sj.feather` on the real index)."""
    import rigel.calibration.splice_graph as SG

    class _JA:
        def __init__(self, n):
            self.edge_row = np.arange(n, dtype=np.int64)

    monkeypatch.setattr(SG, "build_sj_arrays", lambda idx: _JA(len(idx.edges_df)))
    return SG


def _kinds(path, t):
    k, o = path.steps(t)
    return list(zip(k.tolist(), o.tolist(), strict=True))


# ──────────────────────────────────────────────────────────────────────────────
# The include / exclude rule
# ──────────────────────────────────────────────────────────────────────────────


def test_a_MULTI_EXON_transcript_takes_its_exons_and_a_sj_between_them(
    tmp_path, _patched_sj
):
    """⭐ The canonical path: region, splice junction, region. ⛔ The INTRON's region is absent — a mature
    molecule has no intronic bases — and the transcript's two OUTER boundaries are absent, because the
    molecule ends at them."""
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]  # pre | exon1 | INTRON | exon2 | post
    idx = _Index(tmp_path, bounds, {0: [(1_000, 2_000), (9_000, 10_000)]}, strands=[Strand.POS])
    path = build_transcript_path(idx, _Regions(bounds))

    assert _kinds(path, 0) == [
        (STEP_REGION, 1),
        (STEP_SPLICE_SJ, 0),
        (STEP_REGION, 3),
    ]


def test_a_boundary_INTERIOR_to_an_exon_is_crossed_and_included(tmp_path, _patched_sj):
    """⭐ A boundary exists wherever the partition region_bound — including where only a SIGNATURE changed (an
    antisense feature overlapping on the other strand). The molecule crosses it contiguously, so it is
    part of the path. ⛔ Excluding it would drop the only object between two halves of one exon."""
    bounds = [0, 1_000, 1_500, 2_000, 3_000]  # exon [1000,2000) is split at 1500
    idx = _Index(tmp_path, bounds, {0: [(1_000, 2_000)]}, strands=[Strand.POS])
    path = build_transcript_path(idx, _Regions(bounds))

    # region 1, the boundary between regions 1 and 2, region 2 — and nothing else
    assert _kinds(path, 0) == [(STEP_REGION, 1), (STEP_BOUNDARY, 1), (STEP_REGION, 2)]


def test_a_SINGLE_EXON_span_keeps_every_region_it_covers_INTRONS_INCLUDED(
    tmp_path, _patched_sj
):
    """⭐ The synthetic shadow span the index manufactures for a gene: one interval, so every region
    under it is crossed contiguously — including the ones that are introns of its mature sibling."""
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]
    idx = _Index(
        tmp_path,
        bounds,
        {0: [(1_000, 2_000), (9_000, 10_000)]},
        strands=[Strand.POS, Strand.POS],
        synthetic=(1,),
        spans={1: (1_000, 10_000)},
    )
    path = build_transcript_path(idx, _Regions(bounds))

    assert _kinds(path, 1) == [
        (STEP_REGION, 1),
        (STEP_BOUNDARY, 1),
        (STEP_REGION, 2),  # the INTRON — present for the span, absent for its mature sibling
        (STEP_BOUNDARY, 2),
        (STEP_REGION, 3),
    ]
    # ...and the mature sibling still does not see it
    assert (STEP_REGION, 2) not in _kinds(path, 0)


def test_a_transcript_with_no_exons_contributes_an_EMPTY_path(tmp_path, _patched_sj):
    bounds = [0, 1_000, 2_000]
    idx = _Index(tmp_path, bounds, {0: [(0, 1_000)]}, strands=[Strand.POS, Strand.POS])
    path = build_transcript_path(idx, _Regions(bounds))
    assert path.offsets[2] - path.offsets[1] == 0


# ──────────────────────────────────────────────────────────────────────────────
# ⭐⭐ Transcription order
# ──────────────────────────────────────────────────────────────────────────────


def test_a_MINUS_strand_path_runs_in_TRANSCRIPTION_order(tmp_path, _patched_sj):
    """⭐⭐ **THE PROPERTY A GENOMIC-ORDER IMPLEMENTATION GETS SILENTLY WRONG.** A minus-strand molecule
    is transcribed from its HIGH coordinate down, so its path descends. ⚠ Invisible to a consumer that
    treats the path as a set or averages symmetrically over it — and wrong for every consumer that reads
    it as a sequence, which is the whole point of a path.

    The plus- and minus-strand transcripts here occupy the SAME regions, so only the ORDER can differ.
    """
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]
    exons = [(1_000, 2_000), (9_000, 10_000)]
    idx = _Index(
        tmp_path, bounds, {0: exons, 1: exons}, strands=[Strand.POS, Strand.NEG]
    )
    path = build_transcript_path(idx, _Regions(bounds))

    fwd = _kinds(path, 0)
    rev = _kinds(path, 1)
    assert [k for k, _ in rev] == [k for k, _ in reversed(fwd)], (
        "the minus-strand path was emitted in genomic order"
    )
    assert [o for k, o in rev if k == STEP_REGION] == [3, 1]
    assert [o for k, o in fwd if k == STEP_REGION] == [1, 3]
    assert rev[0] == (STEP_REGION, 3), "a minus-strand path must START at its highest region"
    # ⭐ and the two sj ids DIFFER, because strand is part of the sj key — two transcripts on
    # opposite strands splicing identical coordinates are two slots, which is the redundancy the key
    # keeps deliberately for error checking.
    assert [o for k, o in fwd if k == STEP_SPLICE_SJ] != [
        o for k, o in rev if k == STEP_SPLICE_SJ
    ]


# ──────────────────────────────────────────────────────────────────────────────
# The sj join
# ──────────────────────────────────────────────────────────────────────────────


def test_two_transcripts_sharing_one_intron_resolve_to_the_SAME_sj_id(
    tmp_path, _patched_sj
):
    """⭐ A sj is a property of the genome, not of a transcript: two isoforms splicing the same
    intron must land on one slot, because the accumulator tallied their fragments there together."""
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]
    idx = _Index(
        tmp_path,
        bounds,
        {0: [(1_000, 2_000), (9_000, 10_000)], 1: [(1_000, 2_000), (9_000, 11_000)]},
        strands=[Strand.POS, Strand.POS],
    )
    path = build_transcript_path(idx, _Regions(bounds))
    sj0 = [o for k, o in _kinds(path, 0) if k == STEP_SPLICE_SJ]
    sj1 = [o for k, o in _kinds(path, 1) if k == STEP_SPLICE_SJ]
    assert sj0 == sj1 == [0]


def test_an_UNRESOLVABLE_intron_RAISES_rather_than_dropping_a_step(tmp_path, _patched_sj):
    """⛔⛔ **THE SILENT FAILURE THIS GUARD EXISTS FOR.** The sj key is derived from REGION
    boundaries, which equal the intron's coordinates only because the partition region_bounds at every exon
    endpoint — measured 0 of 45,609 violations on the shipped index, but ASSUMED by the derivation.

    If it ever broke, the affected transcript would simply lose a step and its path would still read as
    a shorter, perfectly well-formed walk. ⭐ So the miss is loud: a shorter path is indistinguishable
    from a correct one, and that is exactly the class of wrong answer worth refusing.
    """
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]
    idx = _Index(tmp_path, bounds, {0: [(1_000, 2_000), (9_000, 10_000)]}, strands=[Strand.POS])
    idx.edges_df = idx.edges_df.iloc[:0]  # the sj axis no longer carries this intron
    with pytest.raises(ValueError, match="resolved to no splice-junction slot"):
        build_transcript_path(idx, _Regions(bounds))


# ──────────────────────────────────────────────────────────────────────────────
# Contract
# ──────────────────────────────────────────────────────────────────────────────


def test_the_CSR_is_well_formed(tmp_path, _patched_sj):
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]
    idx = _Index(
        tmp_path, bounds, {0: [(1_000, 2_000), (9_000, 10_000)]}, strands=[Strand.POS, Strand.NEG]
    )
    path = build_transcript_path(idx, _Regions(bounds))
    assert path.offsets.shape == (3,)
    assert path.offsets[0] == 0 and path.offsets[-1] == path.kind.shape[0]
    assert bool(np.all(np.diff(path.offsets) >= 0))
    assert path.kind.dtype == np.int8 and path.obj_id.dtype == np.int64
    assert set(np.unique(path.kind).tolist()) <= {
        STEP_REGION,
        STEP_BOUNDARY,
        STEP_SPLICE_SJ,
    }
