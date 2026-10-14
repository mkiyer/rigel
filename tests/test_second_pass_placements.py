"""``second_pass.score_held_fragments`` on fragments held on their PLACEMENTS: the placement term decides.

A two-contig scan with the same exon geometry on both, UNIQUE contained depth at the first exon of each
contig in a known ratio, and multimapping fragments (``NH = 2``, one hit per contig) contained in both.
The multimapper's two placements tie on every existing factor — the same length, the same unspliced
path, no contested intron, no motif — so only the placement abundance ``A_p`` separates them, and it
must read the unique traffic in exactly the ratio the fixture put there; where the fixture puts none,
the record is undecided and the draw is even; where it puts some on one contig only, the zero stays
hard. Specification: ``second_pass.placement_abundance`` and ``combine_factors``.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.fl import build_fl_models
from rigel.calibration.splice_graph import build_region_partition_arrays, build_sj_arrays
from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.second_pass import choose_hypotheses, drain, score_held_fragments

GENOME = 6_000
_M = 0
N_MULTI = 200


def _gtf() -> str:
    rows = []
    for ref in ("chr1", "chr2"):
        rows += [
            f'{ref}\tt\texon\t1001\t1600\t.\t+\t.\tgene_id "gW_{ref}"; transcript_id "tW_{ref}";\n',
            f'{ref}\tt\texon\t2001\t2200\t.\t+\t.\tgene_id "gW_{ref}"; transcript_id "tW_{ref}";\n',
        ]
    return "".join(rows)


def _scan(tmp_path_factory, label: str, unique_depth: tuple[int, int]):
    """Scan a fixture with ``unique_depth[f]`` uniquely mapped contained pairs in contig ``f``'s first exon and
    ``N_MULTI`` multimapping pairs contained in both; return what the scorer reads."""
    import pysam

    base = tmp_path_factory.mktemp(f"placements_{label}")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    fasta.write_text(
        "".join(
            f">{ref}\n" + "\n".join(["N" * 80] * (GENOME // 80)) + "\n" for ref in ("chr1", "chr2")
        )
    )
    pysam.faidx(str(fasta))
    gtf.write_text(_gtf())
    TranscriptIndex.build(str(fasta), str(gtf), str(base / "idx"), write_tsv=False)
    index = TranscriptIndex.load(str(base / "idx"))

    def read(qname, ref_id, pos, mate_pos, is_r1, nh=1, hi=None, secondary=False):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.reference_id = ref_id
        a.reference_start = pos
        a.mapping_quality = 60
        a.flag = 0x1 | 0x2 | (0x40 | 0x20 if is_r1 else 0x80 | 0x10) | (0x100 if secondary else 0)
        a.cigar = [(_M, 100)]
        a.query_sequence = "A" * 100
        a.query_qualities = pysam.qualitystring_to_array("I" * 100)
        a.next_reference_id = ref_id
        a.next_reference_start = mate_pos
        a.set_tags([("NH", nh, "i")] + ([("HI", hi, "i")] if hi is not None else []))
        return a

    reads = []
    for ref_id, depth in enumerate(unique_depth):
        for k in range(depth):
            reads += [
                read(f"u_{ref_id}_{k}", ref_id, 1100, 1300, True),
                read(f"u_{ref_id}_{k}", ref_id, 1300, 1100, False),
            ]
    for k in range(N_MULTI):
        for ref_id, hi in ((0, 1), (1, 2)):
            reads += [
                read(f"mm_{k}", ref_id, 1120, 1320, True, nh=2, hi=hi, secondary=hi == 2),
                read(f"mm_{k}", ref_id, 1320, 1120, False, nh=2, hi=hi, secondary=hi == 2),
            ]
    bam = str(base / "p.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "chr1", "LN": GENOME}, {"SN": "chr2", "LN": GENOME}],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-n", "-o", bam, bam)
    _s, strand_model, _b, payload = scan_and_buffer(
        bam, index, BamScanConfig(sj_strand_tag="XS", total_threads=1)
    )
    _, _, region_types = build_region_partition_arrays(index)
    return payload, region_types, build_sj_arrays(index), strand_model


def _scored(scanned):
    payload, region_types, sj, strand_model = scanned
    scores = score_held_fragments(
        payload,
        fl_models=build_fl_models(payload),
        rna_sense_frac=strand_model.p_r1_sense,
        region_types=region_types,
        sj=sj,
    )
    d = payload.deferred
    multi = np.flatnonzero(np.diff(d.placement_offsets) > 1)
    runs = d.record_hypothesis_offsets
    return payload, scores, multi, runs, region_types, sj


@pytest.fixture(scope="module")
def three_to_one(tmp_path_factory):
    return _scan(tmp_path_factory, "3to1", (30, 10))


@pytest.fixture(scope="module")
def nothing_unique(tmp_path_factory):
    return _scan(tmp_path_factory, "none", (0, 0))


@pytest.fixture(scope="module")
def one_side_only(tmp_path_factory):
    return _scan(tmp_path_factory, "one", (12, 0))


def test_the_fixture_holds_every_multimapper_on_its_two_placements(three_to_one):
    payload = three_to_one[0]
    d = payload.deferred
    assert payload.qc.deferred_multiple_placements == N_MULTI
    assert int((np.diff(d.placement_offsets) == 2).sum()) == N_MULTI
    assert payload.qc.deferred_undetermined_gap == 0, (
        "no gap is open here: the placement is the only question"
    )
    assert set(d.ref.tolist()) == {0, 1}


def test_unique_contained_depth_in_a_3_to_1_ratio_scores_the_placements_3_to_1(three_to_one):
    payload, scores, multi, runs, _rt, _sj = _scored(three_to_one)
    assert multi.size == N_MULTI
    for i in multi:
        r0, r1 = int(runs[i]), int(runs[i + 1])
        assert r1 - r0 == 2
        a = scores.terms.placement_abundance[r0:r1]
        assert a[0] == pytest.approx(3.0 * a[1], rel=1e-12), (
            "the same geometry, three times the unique traffic"
        )
        np.testing.assert_allclose(scores.score[r0:r1], [0.75, 0.25], rtol=1e-12, atol=0.0)
        # the other factors tie and say nothing: same L, same path, no contested intron, no motif
        assert scores.terms.length[r0] == scores.terms.length[r1 - 1]
        assert scores.terms.density[r0] == scores.terms.density[r1 - 1] == 0.0
        assert scores.terms.strand[r0] == scores.terms.strand[r1 - 1]
    assert scores.n_undecided == 0


def test_the_draw_follows_the_placement_scores_and_the_drain_deposits_where_it_drew(three_to_one):
    """Over 200 records at 0.75 the binomial's own four standard deviations bound the drawn share, and the
    drained contained counts on the two exons move by exactly the records drawn to each."""
    payload, scores, multi, runs, region_types, sj = _scored(three_to_one)
    choices = choose_hypotheses(scores, payload, seed=7)
    chr1 = int((choices[multi] == 0).sum())
    sd = np.sqrt(N_MULTI * 0.75 * 0.25)
    assert abs(chr1 - 0.75 * N_MULTI) < 4.0 * sd, f"{chr1} of {N_MULTI} drawn to chr1 against 0.75"
    drained = drain(payload, choices, region_types=region_types, sj=sj)
    d = payload.deferred

    def exon(ref: int) -> int:  # the first exon region of each contig
        return int(payload.ref_region_offsets[ref]) + 1

    before = np.asarray(payload.region_contained_count, dtype=np.int64).sum(axis=1)
    after = np.asarray(drained.region_contained_count, dtype=np.int64).sum(axis=1)
    assert after[exon(0)] - before[exon(0)] == chr1
    assert after[exon(1)] - before[exon(1)] == N_MULTI - chr1
    assert drained.drain.offered_multimapper == N_MULTI
    assert drained.qc.deferred_multiple_placements == 0 and drained.deferred.n_fragments == 0
    assert d.n_fragments == N_MULTI


def test_no_unique_evidence_anywhere_is_undecided_and_splits_evenly(nothing_unique):
    payload, scores, multi, runs, _rt, _sj = _scored(nothing_unique)
    assert multi.size == N_MULTI
    for i in multi:
        r0, r1 = int(runs[i]), int(runs[i + 1])
        assert (
            scores.terms.placement_abundance[r0] == scores.terms.placement_abundance[r1 - 1] == 0.0
        )
        np.testing.assert_array_equal(scores.score[r0:r1], [0.5, 0.5])
    assert scores.n_undecided == N_MULTI


def test_unique_evidence_on_one_side_only_is_a_hard_zero_on_the_other(one_side_only):
    payload, scores, multi, runs, region_types, sj = _scored(one_side_only)
    for i in multi:
        r0, r1 = int(runs[i]), int(runs[i + 1])
        np.testing.assert_array_equal(scores.score[r0:r1], [1.0, 0.0])
    choices = choose_hypotheses(scores, payload, seed=3)
    assert np.all(choices[multi] == 0)
    drained = drain(payload, choices, region_types=region_types, sj=sj)
    exon2 = int(payload.ref_region_offsets[1]) + 1
    assert int(np.asarray(drained.region_contained_count)[exon2].sum()) == 0, (
        "nothing reaches the silent copy"
    )
