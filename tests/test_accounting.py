"""Every read pair has exactly one outcome and every record exactly one status, both counted by the scan.

A read pair (a template) is a run of records sharing a read name, filtered records included. Its outcome, decided in
this order: excluded (every record unmapped, QC-failed or a skipped duplicate); withheld, by reason (`off_index` >
`chimeric` > `unannotated_splice` > `unpaired`, and `unpaired` by the furthest check any hit reached, `missing_mate` <
`reference` < `strand` < `outward` < `length`); excluded by the multimap option; or placed, unique or multimapping by
its usable hits. A record is filtered, used (its alignment is a location of a usable hit), unused (in a placed
template, in no usable hit) or unplaced (every eligible record of a template with any other outcome). So

    read pairs = excluded + withheld + excluded by the option + placed unique + placed multimapping
    records = filtered + used + unused + unplaced

Each flag counter counts every record carrying its flag; `proper_pair` counts the primary records flagged as a proper
pair. The toy reference is `test_hit_verdict.py`'s: tK (+) [1000,1200) [2100,2400) on chr1, T1 and T2 near 10 kb; the
index holds chr1 only, so chr2 is off the index.
"""

from __future__ import annotations

import random

import pysam
import pytest

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from test_hit_verdict import (
    GENOME,
    R1_FWD,
    R1_REV,
    R2_FWD,
    R2_REV,
    SECONDARY,
    SUPPLEMENTARY,
    _bam,
    _index,
    _pair,
    _record,
)

UNMAPPED, MATE_UNMAPPED, QC_FAIL, DUPLICATE = 0x4, 0x8, 0x200, 0x400


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    return _index(tmp_path_factory.mktemp("accounting"))


def _unmapped(qname, read1, pos, extra=0):
    """A record the aligner left unmapped, placed beside its mate as aligners place one."""
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.flag = 0x1 | UNMAPPED | (0x40 if read1 else 0x80) | extra
    a.reference_id = 0
    a.reference_start = pos
    a.next_reference_id = 0
    a.next_reference_start = pos
    a.query_sequence = "A" * 50
    a.query_qualities = pysam.qualitystring_to_array("I" * 50)
    return a


def _counters_bam_records():
    """Sixteen read names, each a case; the comment gives its outcome and its records' statuses."""
    nh2 = {"nh": 2}
    rec = []
    # a: placed unique, genic (inside tK's first exon): 2 used
    rec += _pair("a_unique", 1050, "50M", 1100, "50M")
    # b: placed multimapping, two usable hits (6000/6200, 8000/8200), intergenic; the cross pairs are rejected,
    # (6000, 8200) for length and (8000, 6200) for facing outward: 4 used
    rec += _pair("b_multi", 6000, "50M", 6200, "50M", **nh2)
    rec += _pair("b_multi", 8000, "50M", 8200, "50M", secondary=True, **nh2)
    # q: placed multimapping; R1's alignment at 6500 is written twice, as STAR writes a mate beside each partner, and
    # pairs with R2 at 6700 and at 6600: two usable hits, and the repeated record is used with the alignment: 4 used
    rec += _pair("q_repeat", 6500, "50M", 6700, "50M", **nh2)
    rec += _pair("q_repeat", 6500, "50M", 6600, "50M", secondary=True, **nh2)
    # c: placed unique; R2's secondary at 9500 pairs only with R1 at 7000, too long: 2 used, 1 unused
    rec += _pair("c_rejected", 7000, "50M", 7200, "50M", **nh2)
    rec.append(_record("c_rejected", R2_REV | SECONDARY, 9500, "50M", 7000, nh=2))
    # d: placed unique; R1's secondary on chr2 has no R2 there, so it is in no candidate hit: 2 used, 1 unused
    rec += _pair("d_nohit", 4000, "50M", 4200, "50M", **nh2)
    rec.append(_record("d_nohit", R1_FWD | SECONDARY, 100, "50M", 4200, nh=2, ref=1))
    # e: placed unique; a supplementary record of R1 that names its primary in SA: 2 used, 1 unused
    rec += _pair("e_supp", 3000, "50M", 3200, "50M")
    rec.append(
        _record("e_supp", R1_FWD | SUPPLEMENTARY, 15000, "25M25S", 3200, sa="chr1,3001,+,50M,60,0;")
    )
    # f: withheld unpaired, missing mate: R2 unmapped (filtered): 1 unplaced, 1 filtered
    rec.append(_record("f_missing", (R1_FWD & ~0x2) | MATE_UNMAPPED, 12000, "50M", 12000))
    rec.append(_unmapped("f_missing", read1=False, pos=12000))
    # g: withheld unpaired, reference: R2 on chr2: 2 unplaced
    rec.append(_record("g_reference", R1_FWD, 13000, "50M", 500))
    rec.append(_record("g_reference", R2_REV, 500, "50M", 13000, ref=1))
    # h: withheld unpaired, strand: both forward: 2 unplaced
    rec.append(_record("h_strand", R1_FWD, 14000, "50M", 14200))
    rec.append(_record("h_strand", R2_FWD, 14200, "50M", 14000))
    # i: withheld unpaired, outward: R1 reverse at 16000, R2 forward at 16300: 2 unplaced
    rec.append(_record("i_outward", R1_REV, 16000, "50M", 16300))
    rec.append(_record("i_outward", R2_FWD, 16300, "50M", 16000))
    # j: withheld unpaired, length: 1,550 bases along the genome, no transcript between: 2 unplaced
    rec += _pair("j_length", 17000, "50M", 18500, "50M")
    # k: withheld, an unannotated splice: [1200,2101) is one base from tK's intron: 2 unplaced
    rec += _pair("k_novel", 1150, "50M901N50M", 2200, "50M", xs="+")
    # l: withheld, off the index: both mates on chr2: 2 unplaced
    rec += _pair("l_off_index", 1000, "50M", 1200, "50M", ref=1)
    # m: withheld, chimeric: R1's record carries SA: 2 unplaced
    rec.append(_record("m_chimeric", R1_FWD, 19000, "50M", 19200, sa="chr1,5001,+,30M20S,60,0;"))
    rec.append(_record("m_chimeric", R2_REV, 19200, "50M", 19000))
    # n: excluded, both mates unmapped, one also QC-failed: 2 filtered
    rec.append(_unmapped("n_unmapped", read1=True, pos=5500, extra=MATE_UNMAPPED | QC_FAIL))
    rec.append(_unmapped("n_unmapped", read1=False, pos=5500, extra=MATE_UNMAPPED))
    # o: withheld unpaired, missing mate: both primaries are skipped duplicates, so its only eligible record is a
    # supplementary one, which is no location: 2 filtered, 1 unplaced
    rec.append(_record("o_nohit", R1_FWD | DUPLICATE, 5600, "50M", 5800))
    rec.append(_record("o_nohit", R2_REV | DUPLICATE, 5800, "50M", 5600))
    rec.append(
        _record(
            "o_nohit", R1_FWD | SUPPLEMENTARY, 15500, "25M25S", 5800, sa="chr1,5601,+,50M,60,0;"
        )
    )
    return rec


#: The counters of the counters BAM, derived by hand from the comments above.
EXPECTED = {
    "total": 40,
    "filtered": 5,
    "qc_fail": 1,
    "unmapped": 3,
    "duplicate": 2,
    "secondary": 6,
    "supplementary": 2,
    "proper_pair": 28,
    "n_read_pairs": 16,
    "n_excluded": 1,
    "n_excluded_by_option": 0,
    "n_placed_unique": 4,
    "n_placed_multimapping": 2,
    "n_fragments": 6,
    "n_intergenic": 5,
    "n_withheld_unpaired": 6,
    "n_withheld_unpaired_missing_mate": 2,
    "n_withheld_unpaired_reference": 1,
    "n_withheld_unpaired_strand": 1,
    "n_withheld_unpaired_outward": 1,
    "n_withheld_unpaired_length": 1,
    "n_withheld_unannotated_splice": 1,
    "n_withheld_off_index": 1,
    "n_withheld_chimeric": 1,
    "n_records_used": 16,
    "n_records_unused": 3,
    "n_records_unplaced": 16,
    "n_hits_rejected_missing_mate": 1,
    "n_hits_rejected_reference": 1,
    "n_hits_rejected_strand": 1,
    "n_hits_rejected_outward": 2,
    "n_hits_rejected_length": 3,
    "n_hits_rejected_contradictory": 0,
    "n_hits_rejected_novel_junction": 1,
    "n_hits_rejected_unsupported": 0,
    "n_hits_rejected_off_index": 1,
    "n_hits_rejected_chimeric": 1,
    "n_nh_above_reported": 0,
}

#: The same BAM under `--no-include-multimap`: the option excludes a read pair with more than one usable hit, b and q,
#: and their 8 records are unplaced. c and d, which the aligner calls multimapping (NH 2 and a secondary record), have
#: one usable hit each and are placed.
EXPECTED_WITHOUT_MULTIMAP = {
    **EXPECTED,
    "n_excluded_by_option": 2,
    "n_placed_multimapping": 0,
    "n_fragments": 4,
    "n_intergenic": 3,
    "n_records_used": 8,
    "n_records_unplaced": 24,
}


def _bam_as_given(base, records, name="given.bam"):
    """The records in the order given, each name's records together or not: no sort."""
    bam = str(base / name)
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": f"chr{i + 1}", "LN": GENOME} for i in range(2)],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in records:
            out.write(r)
    return bam


def _counters(stats, keys):
    return {k: getattr(stats, k) for k in keys}


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("include_multimap", [True, False], ids=["multimap", "no-multimap"])
def test_every_counter_on_the_counters_bam(index, tmp_path, threads, include_multimap):
    expected = EXPECTED if include_multimap else EXPECTED_WITHOUT_MULTIMAP
    stats, _strand, _buffer, _payload = scan_and_buffer(
        _bam(tmp_path, _counters_bam_records()),
        index,
        BamScanConfig(total_threads=threads, include_multimap=include_multimap),
    )
    assert _counters(stats, expected) == expected


def test_the_counters_do_not_depend_on_the_order_of_a_names_records(index, tmp_path):
    """Each name's records shuffled, the names kept together: a template is its records, in any order."""
    by_name: dict[str, list] = {}
    for r in _counters_bam_records():
        by_name.setdefault(r.query_name, []).append(r)
    rng = random.Random(7)
    records = []
    for name in sorted(by_name):
        group = by_name[name]
        rng.shuffle(group)
        records += group
    stats, _strand, _buffer, _payload = scan_and_buffer(
        _bam_as_given(tmp_path, records), index, BamScanConfig(total_threads=1)
    )
    assert _counters(stats, EXPECTED) == EXPECTED


def test_a_filtered_record_of_another_name_ends_a_read_pair(index, tmp_path):
    """A read pair is a run of records sharing a name, filtered records included: a filtered record of another name
    between a pair's two records splits it, in the scan as in the annotation pass, into two read pairs of one mate
    each, both withheld for a missing mate, with the filtered name excluded between them."""
    a1, a2 = _pair("a", 5000, "50M", 5200, "50M")
    records = [
        a1,
        _unmapped("b", True, 5500, MATE_UNMAPPED),
        _unmapped("b", False, 5500, MATE_UNMAPPED),
        a2,
    ]
    stats, _strand, _buffer, _payload = scan_and_buffer(
        _bam_as_given(tmp_path, records), index, BamScanConfig(total_threads=1)
    )
    assert (stats.n_read_pairs, stats.n_excluded, stats.n_fragments) == (3, 1, 0)
    assert stats.n_withheld_unpaired == stats.n_withheld_unpaired_missing_mate == 2
    assert (stats.filtered, stats.n_records_unplaced) == (2, 2)


@pytest.mark.parametrize("order", ["outward first", "length first"])
def test_an_unpaired_read_pair_is_withheld_for_the_furthest_check_any_hit_reached(
    index, tmp_path, order
):
    """R1 at 2700 pairs with R2 at 2600, their 5' ends facing away, and with R2's other alignment at 3900, 1,250
    bases along the genome: the pair rule reached the length check, its last, so the read pair is withheld for its
    length, in whichever order the records come."""
    r1 = _record("two_reasons", R1_FWD, 2700, "50M", 2600)
    outward = _record("two_reasons", R2_REV, 2600, "50M", 2700)
    too_long = _record("two_reasons", R2_REV | SECONDARY, 3900, "50M", 2700)
    records = [r1, outward, too_long] if order == "outward first" else [r1, too_long, outward]
    stats, _strand, _buffer, _payload = scan_and_buffer(
        _bam_as_given(tmp_path, records), index, BamScanConfig(total_threads=1)
    )
    assert (stats.n_hits_rejected_outward, stats.n_hits_rejected_length) == (1, 1)
    assert stats.n_withheld_unpaired == stats.n_withheld_unpaired_length == 1


@pytest.mark.parametrize(
    "field",
    [
        "n_excluded",
        "n_withheld_chimeric",
        "n_placed_unique",
        "n_records_unused",
        "n_withheld_unpaired_length",
    ],
)
def test_counts_that_do_not_reconcile_stop_the_run(index, tmp_path, field):
    """The outcomes sum to the read pairs, the statuses to the records, the unpaired reasons to the unpaired, the placed
    to the fragments, by construction: one counter off, and the scan's check stops the run."""
    from rigel.pipeline import _check_scan_accounting

    stats, _strand, _buffer, _payload = scan_and_buffer(
        _bam(tmp_path, _counters_bam_records()), index, BamScanConfig(total_threads=1)
    )
    _check_scan_accounting(stats)
    setattr(stats, field, getattr(stats, field) + 1)
    with pytest.raises(RuntimeError, match="do not reconcile"):
        _check_scan_accounting(stats)


def test_the_scan_checks_its_counts(index, tmp_path, monkeypatch):
    """The scan runs the reconciliation check on its own counts before anything reads them."""
    import rigel.pipeline as pipeline

    checked = []
    monkeypatch.setattr(pipeline, "_check_scan_accounting", checked.append)
    stats, _strand, _buffer, _payload = scan_and_buffer(
        _bam(tmp_path, _counters_bam_records()), index, BamScanConfig(total_threads=1)
    )
    assert checked == [stats]


# ── the annotated BAM stamps what the scan decided ─────────────────────────────────────────────────────────────────


def _annotated(tmp_path, index, records, **scan):
    """Each record's ZF, ZC, ZH and ZS in a run's annotated BAM, keyed by (name, flag, reference, position)."""
    from rigel.config import EMConfig, PipelineConfig
    from rigel.pipeline import run_pipeline

    out = tmp_path / "annotated.bam"
    config = PipelineConfig(
        scan=BamScanConfig(total_threads=1, **scan),
        em=EMConfig(seed=1, assignment_mode="fractional"),
        annotated_bam_path=out,
    )
    run_pipeline(_bam(tmp_path, records), index, config=config)
    with pysam.AlignmentFile(str(out), "rb") as f:
        return {
            (r.query_name, r.flag, r.reference_name, r.reference_start): {
                t: r.get_tag(t) for t in ("ZF", "ZC", "ZH", "ZS") if r.has_tag(t)
            }
            for r in f.fetch(until_eof=True)
        }


def _key(name, pos, ref="chr1", pick=lambda r: True):
    """The annotated BAM's key of the counters BAM's record at ``pos`` in read pair ``name``."""
    (r,) = [
        r
        for r in _counters_bam_records()
        if r.query_name == name
        and r.reference_start == pos
        and f"chr{r.reference_id + 1}" == ref
        and pick(r)
    ]
    return (name, r.flag, ref, pos)


def test_an_unused_alignment_is_stamped_unused(index, tmp_path):
    """A placed read pair's records outside its usable alignments (a rejected hit's, one that pairs with nothing, a
    supplementary record) carry the fragment's outcome, ZC ``unused_alignment``, ZH 0 and ZS ``unknown``; its used
    records carry its class and their own alignment's splice type."""
    tags = _annotated(tmp_path, index, _counters_bam_records())
    unused = {
        "c_rejected": _key("c_rejected", 9500),
        "d_nohit": _key("d_nohit", 100, ref="chr2"),
        "e_supp": _key("e_supp", 15000),
    }
    used = {
        "c_rejected": _key("c_rejected", 7000),
        "d_nohit": _key("d_nohit", 4000),
        "e_supp": _key("e_supp", 3000),
    }
    for name, key in unused.items():
        assert tags[key] == {
            "ZF": tags[used[name]]["ZF"],
            "ZC": "unused_alignment",
            "ZH": 0,
            "ZS": "unknown",
        }
    placed = {"a_unique", "b_multi", "c_rejected", "d_nohit", "e_supp", "q_repeat"}
    for key, t in tags.items():
        if key[0] in placed and key not in unused.values():
            assert t["ZC"] != "unused_alignment" and t["ZS"] == "unspliced", key


def test_a_read_pair_the_multimap_option_drops_has_no_winning_alignment(index, tmp_path):
    """The option drops b and q, each with two usable hits; c and d, one usable hit each, are placed."""
    tags = _annotated(tmp_path, index, _counters_bam_records(), include_multimap=False)
    dropped = [t for key, t in tags.items() if key[0] in {"b_multi", "q_repeat"}]
    assert len(dropped) == 8
    assert all((t["ZF"], t["ZH"]) == (0x80, 0) for t in dropped)
    placed = [t for key, t in tags.items() if key[0] in {"c_rejected", "d_nohit"}]
    assert len(placed) == 6 and all(t["ZF"] == 0x25 for t in placed)


def test_each_alignment_is_stamped_with_its_own_splice_type(index, tmp_path):
    """A multimapping read pair's placements can splice differently: one across T1's junction and one unspliced, or
    one across tK's intron implicitly and one unspliced. Each record carries its own placement's splice type."""
    records = [
        *_pair("s_annot", 10070, "40M100N10M", 10300, "50M", nh=2, xs="+"),
        *_pair("s_annot", 13000, "50M", 13200, "50M", nh=2, secondary=True),
        *_pair("s_implicit", 1150, "50M", 2150, "50M", nh=2),
        *_pair("s_implicit", 14000, "50M", 14200, "50M", nh=2, secondary=True),
    ]
    tags = _annotated(tmp_path, index, records)

    def splice(name, positions):
        return {t["ZS"] for (n, _f, _ref, pos), t in tags.items() if n == name and pos in positions}

    assert splice("s_annot", {10070, 10300}) == {"spliced_annot"}
    assert splice("s_annot", {13000, 13200}) == {"unspliced"}
    assert splice("s_implicit", {1150, 2150}) == {"spliced_implicit"}
    assert splice("s_implicit", {14000, 14200}) == {"unspliced"}


# ── the annotation pass sees what the scan saw, read pair by read pair ───────────────────────────────────────────


def _two_placed_pairs(first, second):
    """Two read pairs, each placed uniquely with both records used, named ``first`` then ``second``, in this order."""
    return [*_pair(first, 5000, "50M", 5200, "50M"), *_pair(second, 7000, "50M", 7200, "50M")]


def _run_with_annotation(tmp_path, index, bam, before=None):
    from rigel.config import EMConfig, PipelineConfig
    from rigel.pipeline import run_pipeline

    out = tmp_path / "annotated.bam"
    config = PipelineConfig(
        scan=BamScanConfig(total_threads=1),
        em=EMConfig(seed=1, assignment_mode="fractional"),
        annotated_bam_path=out,
    )
    run_pipeline(bam, index, config=config, before_annotated_bam=before)
    return out


def test_the_annotation_pass_sees_what_the_scan_saw(index, tmp_path):
    """The counters BAM, and a BAM a filtered record of another name splits a name in: both passes close each read
    pair at the same record, judge it alike, and the pass's digest equals the scan's."""
    assert _run_with_annotation(tmp_path, index, _bam(tmp_path, _counters_bam_records())).exists()
    a1, a2 = _pair("a", 5000, "50M", 5200, "50M")
    split = [
        a1,
        _unmapped("b", True, 5500, MATE_UNMAPPED),
        _unmapped("b", False, 5500, MATE_UNMAPPED),
        a2,
    ]
    assert _run_with_annotation(
        tmp_path, index, _bam_as_given(tmp_path, split, "split.bam")
    ).exists()


def test_an_annotation_pass_that_judges_otherwise_stops_the_run(index, tmp_path):
    """The pass re-judges every read pair; one told another fragment-length limit than the scan's judges otherwise,
    and the run stops, the annotated BAM deleted, rather than stamp what the scan did not decide."""
    try:
        with pytest.raises(RuntimeError, match="annotation pass"):
            _run_with_annotation(
                tmp_path,
                index,
                _bam(tmp_path, _counters_bam_records()),
                before=lambda _result: index.resolver.set_max_fragment_length(100),
            )
    finally:
        index.resolver.set_max_fragment_length(BamScanConfig().max_frag_length)
    assert not (tmp_path / "annotated.bam").exists()


def test_an_annotation_pass_given_other_read_pairs_stops_the_run(index, tmp_path):
    """Two read pairs alike in every count, their names swapped in the BAM the pass reads: the same outcomes, records
    and statuses, but other read pairs, which only their names tell apart, so the run stops."""
    import shutil

    bam = _bam_as_given(tmp_path, _two_placed_pairs("x1", "x2"), "pairs.bam")
    swapped = _bam_as_given(tmp_path, _two_placed_pairs("x2", "x1"), "swapped.bam")
    with pytest.raises(RuntimeError, match="annotation pass"):
        _run_with_annotation(
            tmp_path, index, bam, before=lambda _result: shutil.copyfile(swapped, bam)
        )
