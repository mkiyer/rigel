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
import pandas as pd
import pytest
from _index_builder import write_fasta

from rigel.calibration.fl import build_fl_models
from rigel.calibration.region_arrays import RegionArrays
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


def _scan(
    tmp_path_factory,
    label: str,
    unique_depth: tuple[int, int],
    unknown=None,
    *,
    background: int = 0,
    multi_nm: tuple[int, int] = (0, 0),
    unique_nm: int = 0,
):
    """Scan a fixture with ``unique_depth[f]`` uniquely mapped contained pairs in contig ``f``'s first exon and
    ``N_MULTI`` multimapping pairs contained in both, the genome holding ``N`` over *unknown*'s runs; return what
    the scorer reads. ``background`` unique pairs per contig lie in its last intergenic region, the library's gDNA
    bulk; a multimapper's read 1 carries ``multi_nm[f]`` mismatches on contig ``f``, and every unique read
    ``unique_nm``, which sets the library's mismatch rate."""
    import pysam

    base = tmp_path_factory.mktemp(f"placements_{label}")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    write_fasta(fasta, {"chr1": GENOME, "chr2": GENOME}, unknown)
    gtf.write_text(_gtf())
    TranscriptIndex.build(str(fasta), str(gtf), str(base / "idx"), write_tsv=False)
    index = TranscriptIndex.load(str(base / "idx"))

    def read(qname, ref_id, pos, mate_pos, is_r1, nh=1, hi=None, secondary=False, nm=0):
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
        a.set_tags(
            [("NH", nh, "i"), ("NM", nm, "i")] + ([("HI", hi, "i")] if hi is not None else [])
        )
        return a

    reads = []
    for ref_id, depth in enumerate(unique_depth):
        for k in range(depth):
            reads += [
                read(f"u_{ref_id}_{k}", ref_id, 1100, 1300, True, nm=unique_nm),
                read(f"u_{ref_id}_{k}", ref_id, 1300, 1100, False, nm=unique_nm),
            ]
    for ref_id in (0, 1):
        for k in range(background):
            pos = 3000 + 10 * k
            reads += [
                read(f"b_{ref_id}_{k}", ref_id, pos, pos + 200, True, nm=unique_nm),
                read(f"b_{ref_id}_{k}", ref_id, pos + 200, pos, False, nm=unique_nm),
            ]
    for k in range(N_MULTI):
        for ref_id, hi in ((0, 1), (1, 2)):
            reads += [
                read(
                    f"mm_{k}",
                    ref_id,
                    1120,
                    1320,
                    True,
                    nh=2,
                    hi=hi,
                    secondary=hi == 2,
                    nm=multi_nm[ref_id],
                ),
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
    return (
        payload,
        region_types,
        build_sj_arrays(index),
        strand_model,
        RegionArrays.from_index(index),
    )


def _scored(scanned):
    payload, region_types, sj, strand_model, regions = scanned
    scores = score_held_fragments(
        payload,
        fl_models=build_fl_models(payload),
        rna_sense_frac=strand_model.p_r1_sense,
        regions=regions,
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


@pytest.fixture(scope="module")
def one_side_only_over_a_bulk(tmp_path_factory):
    return _scan(tmp_path_factory, "one_bulk", (12, 0), background=100)


#: A library whose unique reads carry one mismatch a mate in 100 bases (rate 0.01), and a multimapper whose read 1
#: matches chr1's copy and mismatches chr2's once; the unique depth is the same on both copies.
@pytest.fixture(scope="module")
def one_mismatch_apart(tmp_path_factory):
    return _scan(tmp_path_factory, "nm", (10, 10), multi_nm=(0, 1), unique_nm=1)


@pytest.fixture(scope="module")
def one_mismatch_in_a_clean_library(tmp_path_factory):
    return _scan(tmp_path_factory, "nm0", (10, 10), multi_nm=(0, 1), unique_nm=0)


#: chr2's copy of the first exon ([1000, 1600)) holds an N run beyond every read, so the index masks it: of the
#: library's 300 bp fragments the copy's observable support is 151 starts, chr1's whole copy 301.
_MASKED_COPY = {"chr2": [(1450, 1600)]}


@pytest.fixture(scope="module")
def one_density_one_copy_masked(tmp_path_factory):
    return _scan(tmp_path_factory, "masked", (301, 151), unknown=_MASKED_COPY)


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


def _bulk_by_hand(scanned) -> float:
    """The pass-one bulk counted directly: the unique contained fragments over the intergenic regions' observable
    contained support, in the unit a contained placement's term reads."""
    from rigel.calibration.effective_length import observable_contained_eff_length
    from rigel.calibration.signature import coarse_type_array

    payload, _rt, _sj, _sm, regions = scanned
    support = observable_contained_eff_length(regions, build_fl_models(payload).global_pmf)
    intergenic = coarse_type_array(regions.signature) == 0
    counts = np.asarray(payload.region_contained_count, dtype=np.float64).sum(axis=1)
    return float(counts[intergenic].sum() / support[intergenic].sum())


def test_a_copy_with_no_unique_traffic_reads_the_bulk_beside_one_that_has_some(
    one_side_only_over_a_bulk,
):
    """gDNA's floor: a copy every fragment of which multimaps reads no unique traffic, but gDNA is one template at a
    uniform rate before capture and capture only enriches, so a place gDNA can be holds at least the bulk. The silent
    copy reads the library's bulk, not a hard zero, and the draw weighs the copies by their terms.

    PERTURBATION: without the floor the silent copy reads 0 and is never drawn; with a sum in place of the bound the
    copy WITH traffic reads its own density plus the bulk.
    """
    from rigel.calibration.effective_length import observable_contained_eff_length

    payload, scores, multi, runs, _rt, _sj = _scored(one_side_only_over_a_bulk)
    bulk = _bulk_by_hand(one_side_only_over_a_bulk)
    regions = one_side_only_over_a_bulk[4]
    support = observable_contained_eff_length(regions, build_fl_models(payload).global_pmf)
    exon = int(payload.ref_region_offsets[0]) + 1  # chr1's first exon, where its unique pairs sit
    own = float(np.asarray(payload.region_contained_count[exon], dtype=np.float64).sum()) / float(
        support[exon]
    )
    assert own > bulk > 0.0
    for i in multi:
        r0, r1 = int(runs[i]), int(runs[i + 1])
        a = scores.terms.placement_abundance[r0:r1]
        assert a[0] == pytest.approx(own, rel=1e-12), "chr1's copy keeps its own unique density"
        assert a[1] == pytest.approx(bulk, rel=1e-12), "chr2's silent copy reads the bulk"
        np.testing.assert_allclose(scores.score[r0:r1], a / a.sum(), rtol=1e-12, atol=0.0)


def test_two_copies_one_mismatch_apart_are_drawn_at_the_mismatch_price(one_mismatch_apart):
    """The draw reads a multimapper's sequence at the EM's own prices: at the library's mismatch rate 0.01 a copy
    with one more mismatch is weighed by ``ε / (3(1 − ε))`` against the one it matches, everything else tied.

    PERTURBATION: a draw that reads no sequence splits these records evenly.
    """
    from rigel.scoring import sequence_log_penalties

    payload, scores, multi, runs, _rt, _sj = _scored(one_mismatch_apart)
    mismatch, _clip = sequence_log_penalties(0.01)
    assert payload.mismatch_rate == pytest.approx(0.01, rel=1e-12)
    assert multi.size == N_MULTI
    for i in multi:
        r0, r1 = int(runs[i]), int(runs[i + 1])
        np.testing.assert_array_equal(scores.terms.sequence[r0:r1], [0.0, mismatch])
        np.testing.assert_allclose(
            scores.score[r0:r1],
            np.array([1.0, np.exp(mismatch)]) / (1.0 + np.exp(mismatch)),
            rtol=1e-12,
        )


def test_in_a_library_with_no_mismatch_a_mismatched_copy_is_never_drawn(
    one_mismatch_in_a_clean_library,
):
    """At a mismatch rate of zero a mismatch is impossible, so the mismatched copy takes no part, as in the EM."""
    payload, scores, multi, runs, region_types, sj = _scored(one_mismatch_in_a_clean_library)
    assert payload.mismatch_rate == 0.0
    for i in multi:
        r0, r1 = int(runs[i]), int(runs[i + 1])
        np.testing.assert_array_equal(scores.score[r0:r1], [1.0, 0.0])
    choices = choose_hypotheses(scores, payload, seed=5)
    assert np.all(choices[multi] == 0)


def test_one_density_over_observable_bases_scores_the_placements_evenly(
    one_density_one_copy_masked,
):
    """Unique depth 301 on chr1's copy and 151 on chr2's, whose last 150 bases are N: one density per observable
    start on both, so the placement term ties and the draw is even. Read over the whole exon the masked copy's
    density is understated by 151/301 and the draw weighs it low, 0.67 against 0.33."""
    payload, scores, multi, runs, _rt, _sj = _scored(one_density_one_copy_masked)
    assert multi.size == N_MULTI
    assert one_density_one_copy_masked[4].unobservable.tolist() == [[1, 1450, 1600]]
    for i in multi:
        r0, r1 = int(runs[i]), int(runs[i + 1])
        a = scores.terms.placement_abundance[r0:r1]
        assert a[0] == pytest.approx(a[1], rel=1e-12), (
            "one density on both copies' observable bases"
        )
        np.testing.assert_allclose(scores.score[r0:r1], [0.5, 0.5], rtol=1e-12, atol=0.0)
    assert scores.n_undecided == N_MULTI, "a tie on every record, so every draw is a fair coin"


def test_a_frame_of_another_partition_is_refused(three_to_one):
    """The frame is the partition the scan deposited into, or the placement term would read another region's
    support: one with another region count is refused before anything is scored."""
    payload, _region_types, sj, strand_model, _regions = three_to_one
    other = RegionArrays.from_frame(
        pd.DataFrame({"ref_name": ["chr1"], "start": [0], "end": [GENOME], "signature": [0]}),
        {"chr1": 0, "chr2": 1},
    )
    with pytest.raises(ValueError, match="different partitions"):
        score_held_fragments(
            payload,
            fl_models=build_fl_models(payload),
            rna_sense_frac=strand_model.p_r1_sense,
            regions=other,
            sj=sj,
        )
