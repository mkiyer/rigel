"""A pair of mates is a hit only when the two can be one molecule with a valid fragment length.

Every accepted alignment needs a fragment length, which the EM scores by. A hit's two mates must lie on one reference,
on opposite strands, facing inward once each 5' soft clip is walked back, with both mates' aligned bases inside the
molecule; its length, the molecule's bases with every observed skip removed, must lie in 1..max_frag_length along the
genome, or along one annotated transcript that holds every junction the pair observes, that transcript's other introns
in the unaligned stretches cut by the resolver's K rule. A hit with one mate gives no length at all. A hit that fails
is dropped and counted by its reason; a template left with no hit is withheld, not counted intergenic. Measured on the STAR panel (47.9 M templates, per-base
truth), the rule rejects no correct pair (`ISSUES: star-hi-pairs-against-rigels-own-pairing`).

The toy reference, 0-based half-open:
- tK (+): exons [1000,1200) [2100,2400), the intron [1200,2100) for the tolerance case;
- T1 (+): [10050,10110) [10210,10400); T2 (+): [10250,10300) [11300,11800); T3 (+, the second index only):
  [10050,10110) [10210,10300) [11300,11800). T1 holds the hybrid pair's observed junction [10110,10210), T2 the intron
  [10300,11300) that would shorten it, and only T3 holds both.
"""

from __future__ import annotations

import pysam
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import run_pipeline, scan_and_buffer

GENOME = 20_000
R1_FWD, R2_REV = 0x1 | 0x2 | 0x40 | 0x20, 0x1 | 0x2 | 0x80 | 0x10
SECONDARY = 0x100


def _gtf(with_t3: bool) -> str:
    rows = [
        ("gK", "tK", [(1000, 1200), (2100, 2400)]),
        ("gH", "T1", [(10050, 10110), (10210, 10400)]),
        ("gH", "T2", [(10250, 10300), (11300, 11800)]),
    ]
    if with_t3:
        rows.append(("gH", "T3", [(10050, 10110), (10210, 10300), (11300, 11800)]))
    return "".join(
        f'chr1\tt\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "{g}"; transcript_id "{t}";\n'
        for g, t, exons in rows
        for s, e in exons
    )


def _index(base, with_t3=False):
    fasta, gtf = base / "g.fa", base / "a.gtf"
    fasta.write_text(">chr1\n" + "\n".join(["N" * 80] * (GENOME // 80)) + "\n")
    pysam.faidx(str(fasta))
    gtf.write_text(_gtf(with_t3))
    TranscriptIndex.build(str(fasta), str(gtf), str(base / "idx"), write_tsv=False)
    return TranscriptIndex.load(str(base / "idx"))


def _record(qname, flag, pos, cigar, mate_pos, nh=1, hi=None, xs=None):
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.flag = flag
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigarstring = cigar
    n = a.query_length or sum(int(x) for x in __import__("re").findall(r"(\d+)[MIS=X]", cigar))
    a.query_sequence = "A" * n
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.next_reference_id = 0
    a.next_reference_start = mate_pos
    tags = [("NH", nh, "i"), ("NM", 0, "i")] + ([("HI", hi, "i")] if hi is not None else [])
    if xs is not None:
        tags.append(("XS", xs, "A"))
    a.set_tags(tags)
    return a


def _pair(qname, f_pos, f_cigar, r_pos, r_cigar, nh=1, hi=None, secondary=False, xs=None):
    extra = SECONDARY if secondary else 0
    return [
        _record(qname, R1_FWD | extra, f_pos, f_cigar, r_pos, nh, hi, xs),
        _record(qname, R2_REV | extra, r_pos, r_cigar, f_pos, nh, hi, xs),
    ]


def _bam(base, records, name="p.bam"):
    bam = str(base / name)
    header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": GENOME}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in records:
            out.write(r)
    pysam.sort("-n", "-o", bam, bam)
    return bam


def _scan(base, index, records, **scan):
    return _scan_full(base, index, records, **scan)[0]


def _scan_full(base, index, records, **scan):
    """The scan's stats and its accumulator payload."""
    stats, _strand, _buffer, payload = scan_and_buffer(
        _bam(base, records), index, BamScanConfig(total_threads=1, **scan)
    )
    return stats, payload


REASONS = ("missing_mate", "reference", "strand", "outward", "length")


def _rejected(stats):
    return {r: getattr(stats, f"n_pairs_rejected_{r}") for r in REASONS}


CONTROL = ("control", 1010, "100M", 1100, "100M")  # inward, L 190 along the genome


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    return _index(tmp_path_factory.mktemp("pairs"))


@pytest.fixture(scope="module")
def index_t3(tmp_path_factory):
    return _index(tmp_path_factory.mktemp("pairs_t3"), with_t3=True)


def test_a_valid_pair_is_a_fragment(index, tmp_path):
    stats = _scan(tmp_path, index, _pair(*CONTROL))
    assert stats.n_fragments == 1
    assert _rejected(stats) == dict.fromkeys(REASONS, 0)
    assert stats.n_withheld_no_valid_pair == 0


def test_a_pair_facing_outward_is_withheld_not_intergenic(index, tmp_path):
    """F forward at 5300, R reverse ending at 5200: their 5' ends point away from each other."""
    stats = _scan(tmp_path, index, _pair("outward", 5300, "100M", 5100, "100M"))
    assert stats.n_fragments == 0, "an unclipped pair facing outward is no molecule"
    assert stats.n_intergenic_unspliced == 0
    assert _rejected(stats)["outward"] == 1
    assert stats.n_withheld_no_valid_pair == 1


def test_a_5_prime_soft_clip_is_walked_back(index, tmp_path):
    """F `10S90M` at 1020 and R at [1012,1112): on their aligned bases R starts before F's, but F's 5' end, its clip
    walked back, is at 1010, so the mates face inward over [1010,1112) and L is 102."""
    stats = _scan(tmp_path, index, _pair("clipped", 1020, "10S90M", 1012, "100M"))
    assert stats.n_fragments == 1 and sum(_rejected(stats).values()) == 0


def test_a_pair_too_long_along_every_path_is_withheld(index, tmp_path):
    """Mates 2,100 bases apart over no annotated intron: L 2,100 against max_frag_length 1,000."""
    stats = _scan(tmp_path, index, _pair("far", 6000, "100M", 8000, "100M"))
    assert stats.n_fragments == 0
    assert _rejected(stats)["length"] == 1
    assert stats.n_withheld_no_valid_pair == 1


def test_the_tolerance_lets_a_mate_run_k_bases_into_an_intron(index, tmp_path):
    """F [1103,1203) runs 3 bases into tK's intron [1200,2100), R at [2150,2250): 1,147 along the genome, but tK's
    intron is cut whole at K 3, giving L 247. At K 2 the 3 bases contradict tK, which then offers no path."""
    records = _pair("tolerance", 1103, "100M", 2150, "100M")
    at_k3 = _scan(tmp_path, index, records, intron_overhang_tolerance=3)
    assert at_k3.n_fragments == 1 and sum(_rejected(at_k3).values()) == 0
    at_k2 = _scan(tmp_path, index, records, intron_overhang_tolerance=2)
    assert at_k2.n_fragments == 0 and _rejected(at_k2)["length"] == 1


def test_one_transcript_witnesses_the_whole_path(index, index_t3, tmp_path):
    """F `10M100N10M` at 10100 (T1's junction) and R `10M` at 11700: 1,510 along the genome and along T1. T2's
    intron [10300,11300) would give 510, but T2 does not hold the observed junction, so no hybrid path is formed
    (targeted review, finding 1). With T3, which holds both, the pair is a hit at L 510."""
    records = _pair("hybrid", 10100, "10M100N10M", 11700, "10M", xs="+")
    without = _scan(tmp_path, index, records)
    assert without.n_fragments == 0 and _rejected(without)["length"] == 1
    with_t3 = _scan(tmp_path, index_t3, records)
    assert with_t3.n_fragments == 1 and sum(_rejected(with_t3).values()) == 0


def test_an_invalid_hit_is_dropped_and_its_template_kept(index, tmp_path):
    """HI 1 is the control pair, HI 2 a pair facing outward: one hit dropped, one fragment kept."""
    records = _pair("mixed", 1010, "100M", 1100, "100M", nh=2, hi=1) + _pair(
        "mixed", 5300, "100M", 5100, "100M", nh=2, hi=2, secondary=True
    )
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 1
    assert _rejected(stats)["outward"] == 1
    assert stats.n_withheld_no_valid_pair == 0


def _unmapped_mate(qname, mate_pos):
    """R2 unmapped, placed beside its mate as aligners write it; the scan reads no unmapped record."""
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.flag = 0x1 | 0x80 | 0x4
    a.reference_id = 0
    a.reference_start = mate_pos
    a.next_reference_id = 0
    a.next_reference_start = mate_pos
    a.query_sequence = "A" * 100
    a.query_qualities = pysam.qualitystring_to_array("I" * 100)
    return a


def test_a_mate_whose_partner_is_unmapped_is_no_hit(index, tmp_path):
    """R1 at [1010,1110) and R2 unmapped: one mate gives no fragment length, so the template is withheld rather than
    read as a 100-bp fragment."""
    records = [_record("lone", R1_FWD | 0x8, 1010, "100M", 1010), _unmapped_mate("lone", 1010)]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 0, "a lone mate is no molecule"
    assert _rejected(stats)["missing_mate"] == 1
    assert stats.n_withheld_no_valid_pair == 1


def test_an_hi_group_missing_a_mate_is_dropped_and_its_template_kept(index, tmp_path):
    """HI 1 holds both mates; HI 2 holds R1 alone, a half-filled group. The half group is no hit, so calibration is
    offered one placement, not a read-length second one."""
    records = _pair("half", 1010, "100M", 1100, "100M", nh=2, hi=1) + [
        _record("half", R1_FWD | SECONDARY, 1500, "100M", 1600, nh=2, hi=2)
    ]
    stats, payload = _scan_full(tmp_path, index, records)
    assert payload.qc.deferred_multiple_placements == 0, "the one-mate group became a placement"
    assert stats.n_fragments == 1
    assert _rejected(stats)["missing_mate"] == 1
    assert stats.n_withheld_no_valid_pair == 0


def test_a_lone_secondary_without_hi_is_no_hit(index, tmp_path):
    """Without HI, a secondary R1 whose mate is unmapped becomes a one-mate hit; it is dropped, and the primary pair
    stays, calibration's one placement."""
    records = _pair("nohi", 1010, "100M", 1100, "100M", nh=2) + [
        _record("nohi", R1_FWD | SECONDARY | 0x8, 1500, "100M", 1500, nh=2)
    ]
    stats, payload = _scan_full(tmp_path, index, records)
    assert payload.qc.deferred_multiple_placements == 0, "the lone secondary became a placement"
    assert stats.n_fragments == 1
    assert _rejected(stats)["missing_mate"] >= 1
    assert stats.n_withheld_no_valid_pair == 0


# ── multiplicity: a template's distinct valid hits ──────────────────────────────────────────────────


def test_a_template_whose_other_hit_is_invalid_is_unique(index, tmp_path):
    """NH 2, HI 1 a valid pair and HI 2 one facing outward: one hit is left, so the template is a unique mapper, not
    a multimapper whose second place does not exist."""
    records = _pair("one_left", 1010, "100M", 1100, "100M", nh=2, hi=1) + _pair(
        "one_left", 5300, "100M", 5100, "100M", nh=2, hi=2, secondary=True
    )
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 1
    assert stats.n_multimapper_groups == 0, "NH, not the hits left, still set the multiplicity"


def test_two_distinct_valid_hits_are_a_multimapper(index, tmp_path):
    records = _pair("two", 1010, "100M", 1100, "100M", nh=2, hi=1) + _pair(
        "two", 1500, "100M", 1600, "100M", nh=2, hi=2, secondary=True
    )
    stats = _scan(tmp_path, index, records)
    assert stats.n_multimapper_groups == 1


def test_a_spliced_and_a_contiguous_reading_stay_two_hits(index, tmp_path):
    """R1 alike in both groups; R2 `100M` at 1100 in HI 1, `50M300N50M` from 1100 in HI 2. One start, two readings: the
    competition between them is the point (owner, 2026-10-18), so nothing may merge them into one hit."""
    records = _pair("compete", 1010, "100M", 1100, "100M", nh=2, hi=1) + _pair(
        "compete", 1010, "100M", 1100, "50M300N50M", nh=2, hi=2, secondary=True, xs="+"
    )
    stats = _scan(tmp_path, index, records)
    assert stats.n_multimapper_groups == 1


def test_an_nh_above_the_reported_groups_is_counted(index, tmp_path):
    """NH 3 over two HI groups: a truncated output, counted; the two reported hits are the template's multiplicity."""
    records = _pair("trunc", 1010, "100M", 1100, "100M", nh=3, hi=1) + _pair(
        "trunc", 1500, "100M", 1600, "100M", nh=3, hi=2, secondary=True
    )
    stats = _scan(tmp_path, index, records)
    assert stats.n_nh_above_reported == 1
    assert stats.n_multimapper_groups == 1


def test_the_annotated_bam_labels_a_withheld_template_unresolved(index, tmp_path):
    """A withheld template's records are written unresolved (ZF 0), class `unpaired`, not as intergenic gDNA."""
    records = [
        r for k in range(200) for r in _pair(f"ok{k:03d}", 1010 + k, "100M", 1100 + k, "100M")
    ]
    records += _pair("outward", 5300, "100M", 5100, "100M")
    out = tmp_path / "annotated.bam"
    config = PipelineConfig(
        em=EMConfig(seed=7), scan=BamScanConfig(total_threads=1), annotated_bam_path=out
    )
    run_pipeline(_bam(tmp_path, records), index, config=config)
    with pysam.AlignmentFile(str(out), "rb") as bam:
        tags = {(r.query_name, r.is_read1): (r.get_tag("ZF"), r.get_tag("ZC")) for r in bam}
    assert tags[("outward", True)] == tags[("outward", False)] == (0, "unpaired")
    assert tags[("ok000", True)][0] & 0x01, "a placed template stays resolved"
