"""A hit is usable only when its two mates can be one molecule with a valid fragment length, and, when it is spliced,
one annotated transcript explains it: overlaps both mates and holds every junction it observes.

Every accepted alignment needs a fragment length, which the EM scores by. A hit's two mates must lie on one reference,
on opposite strands, their 5' ends facing each other once each 5' soft clip is walked back. The molecule spans both
mates, so a mate may run past the other's 5' end, as a read trimmed at its 5' end starts inside its molecule (owner,
2026-10-20). Its length, the molecule's bases with every observed skip removed, must lie in 1..max_frag_length along the
genome, or along one annotated transcript that overlaps both mates and holds every junction the hit observes, that
transcript's other introns in the unaligned stretches cut by the resolver's K rule. A transcript that misses a mate
never shortens the path between them (owner, 2026-10-19: overlap, not containment). A hit with one mate gives no
length at all. Measured on the STAR panel (per-base truth), the pair rule rejects no correct pair
(`ISSUES: star-hi-pairs-against-rigels-own-pairing`).

The annotated-path rule: a spliced hit is usable only if one annotated transcript overlaps both mates and holds every
junction the hit observes; Rigel does not assemble transcripts. Otherwise it is contradictory (two junctions overlap
unequal, or their motifs name both strands), a novel junction (one no transcript holds), or unsupported (each held, but
no one transcript holds them all over both mates). The molecule is judged first, so the rule withholds only molecules,
each with a fragment length.

A template's hits come from geometry alone, for any aligner: every R1 place with every R2 place on one reference with
opposite strands, the aligner's own pairing (STAR's HI) not read. One alignment written in several records (STAR writes a
mate's record beside each of its partners) is one place. A hit on a reference the index does not hold is dropped, since
Rigel can place nothing there (owner, 2026-10-20). So is a hit with a mate aligned in pieces, a chimeric alignment:
its record carries SA, or its 5' end is hard-clipped (owner, 2026-10-20); a supplementary record joins no hit. A hit
that fails is dropped and counted by its reason; a template left with no hit is withheld, not counted intergenic: off
the index when one of its hits lies on a reference outside the index, else chimeric when one of its hits is, else an
unannotated splice when the annotated-path rule withheld one of its hits, else unpaired.

The toy reference, 0-based half-open:
- tK (+): exons [1000,1200) [2100,2400), the intron [1200,2100) for the tolerance case;
- T1 (+): [10050,10110) [10210,10400); T2 (+): [10250,10300) [11300,11800); T3 (+, the second index only):
  [10050,10110) [10210,10300) [11300,11800). T1 holds the hybrid pair's observed junction [10110,10210), T2 the intron
  [10300,11300) that would shorten it, and only T3 holds both.
"""

from __future__ import annotations

import pandas as pd
import pysam
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import run_pipeline, scan_and_buffer

GENOME = 20_000
R1_FWD, R2_REV = 0x1 | 0x2 | 0x40 | 0x20, 0x1 | 0x2 | 0x80 | 0x10
R1_REV, R2_FWD = R1_FWD ^ 0x30, R2_REV ^ 0x30  # each strand bit and its mate's flipped
SECONDARY = 0x100
SUPPLEMENTARY = 0x800


TOY = [
    ("gK", "tK", [(1000, 1200), (2100, 2400)]),
    ("gH", "T1", [(10050, 10110), (10210, 10400)]),
    ("gH", "T2", [(10250, 10300), (11300, 11800)]),
]
T3 = ("gH", "T3", [(10050, 10110), (10210, 10300), (11300, 11800)])


def _gtf(rows) -> str:
    return "".join(
        f'chr1\tt\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "{g}"; transcript_id "{t}";\n'
        for g, t, exons in rows
        for s, e in exons
    )


def _index(base, rows=TOY, n_refs=1):
    fasta, gtf = base / "g.fa", base / "a.gtf"
    fasta.write_text(
        "".join(
            f">chr{i + 1}\n" + "\n".join(["N" * 80] * (GENOME // 80)) + "\n" for i in range(n_refs)
        )
    )
    pysam.faidx(str(fasta))
    gtf.write_text(_gtf(rows))
    TranscriptIndex.build(str(fasta), str(gtf), str(base / "idx"), write_tsv=False)
    return TranscriptIndex.load(str(base / "idx"))


def _record(
    qname, flag, pos, cigar, mate_pos, nh=1, hi=None, xs=None, ts=None, ref=0, nm=0, sa=None
):
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.flag = flag
    a.reference_id = ref
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigarstring = cigar
    n = a.query_length or sum(int(x) for x in __import__("re").findall(r"(\d+)[MIS=X]", cigar))
    a.query_sequence = "A" * n
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.next_reference_id = ref
    a.next_reference_start = mate_pos
    tags = [(tag, v, "i") for tag, v in (("NH", nh), ("NM", nm), ("HI", hi)) if v is not None]
    if xs is not None:
        tags.append(("XS", xs, "A"))
    if ts is not None:
        tags.append(("ts", ts, "A"))
    if sa is not None:
        tags.append(("SA", sa, "Z"))
    a.set_tags(tags)
    return a


def _pair(qname, f_pos, f_cigar, r_pos, r_cigar, nh=1, hi=None, secondary=False, xs=None, ref=0):
    extra = SECONDARY if secondary else 0
    return [
        _record(qname, R1_FWD | extra, f_pos, f_cigar, r_pos, nh, hi, xs, ref=ref),
        _record(qname, R2_REV | extra, r_pos, r_cigar, f_pos, nh, hi, xs, ref=ref),
    ]


def _bam(base, records, name="p.bam"):
    """The records, name-sorted, under chr1 and as many further references (chr2, ...) as they name; the index has
    chr1 only."""
    bam = str(base / name)
    n_refs = 1 + max(r.reference_id for r in records)
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": f"chr{i + 1}", "LN": GENOME} for i in range(n_refs)],
    }
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


def _scan_buffered(base, index, records, column):
    """The scan's stats, and one column of every fragment it buffered."""
    stats, _strand, buffer, _payload = scan_and_buffer(
        _bam(base, records), index, BamScanConfig(total_threads=1)
    )
    return stats, [
        int(x) for chunk in buffer.iter_chunks_consuming() for x in getattr(chunk, column)
    ]


REASONS = (
    "missing_mate",
    "reference",
    "strand",
    "outward",
    "length",
    "contradictory",
    "novel_junction",
    "unsupported",
    "off_index",
    "chimeric",
)


def _rejected(stats):
    return {r: getattr(stats, f"n_hits_rejected_{r}") for r in REASONS}


CONTROL = ("control", 1010, "100M", 1100, "100M")  # inward, L 190 along the genome


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    return _index(tmp_path_factory.mktemp("pairs"))


@pytest.fixture(scope="module")
def index_t3(tmp_path_factory):
    return _index(tmp_path_factory.mktemp("pairs_t3"), [*TOY, T3])


@pytest.fixture(scope="module")
def index_two_refs(tmp_path_factory):
    """The toy index over chr1 and chr2; every transcript on chr1."""
    return _index(tmp_path_factory.mktemp("pairs_two_refs"), n_refs=2)


def test_a_valid_pair_is_a_fragment(index, tmp_path):
    stats = _scan(tmp_path, index, _pair(*CONTROL))
    assert stats.n_fragments == 1
    assert _rejected(stats) == dict.fromkeys(REASONS, 0)
    assert stats.n_withheld_unpaired == 0


def test_a_pair_facing_outward_is_withheld_not_intergenic(index, tmp_path):
    """F forward at 5300, R reverse ending at 5200: their 5' ends point away from each other."""
    stats = _scan(tmp_path, index, _pair("outward", 5300, "100M", 5100, "100M"))
    assert stats.n_fragments == 0, "an unclipped pair facing outward is no molecule"
    assert stats.n_intergenic == 0
    assert _rejected(stats)["outward"] == 1
    assert stats.n_withheld_unpaired == 1


def test_a_5_prime_soft_clip_is_walked_back(index, tmp_path):
    """F `10S90M` at 1020 and R at [1012,1112): on their aligned bases R starts before F's, but F's 5' end, its clip
    walked back, is at 1010, so the mates face inward over [1010,1112) and L is 102."""
    stats = _scan(tmp_path, index, _pair("clipped", 1020, "10S90M", 1012, "100M"))
    assert stats.n_fragments == 1 and sum(_rejected(stats).values()) == 0


OVERRUN = (
    "overrun",
    1014,
    "42M",
    1010,
    "42M",
)  # 5' ends 1014 and 1052; each mate runs 4 past the other's


def test_mates_that_run_past_each_others_5_prime_ends_are_one_molecule(index, tmp_path):
    """F `42M` at [1014,1056) and R `42M` at [1010,1052): their 5' ends face each other, and each mate's aligned bases
    run 4 past the other's 5' end, as when both reads were trimmed at their 5' ends. One molecule spans both mates,
    [1010,1056): its footprint is 46."""
    stats, footprint = _scan_buffered(tmp_path, index, _pair(*OVERRUN), "genomic_footprint")
    assert stats.n_fragments == 1 and sum(_rejected(stats).values()) == 0, (
        "mates running past each other's 5' ends were read as facing away"
    )
    assert footprint == [46]


def test_mates_whose_5_prime_ends_meet_face_away(index, tmp_path):
    """F at [1052,1094) and R at [1010,1052): the forward mate starts where the reverse one ends, so their 5' ends do
    not face each other, and the mates are no molecule however they abut."""
    stats = _scan(tmp_path, index, _pair("abut", 1052, "42M", 1010, "42M"))
    assert stats.n_fragments == 0 and _rejected(stats)["outward"] == 1


def test_a_molecule_past_both_5_prime_ends_is_measured_whole(index, tmp_path):
    """The pair above spans 46 bases, though its 5' ends lie 38 apart: at max_frag_length 45 it is too long, and at 46
    it is a molecule."""
    too_long = _scan(tmp_path, index, _pair(*OVERRUN), max_frag_length=45)
    assert too_long.n_fragments == 0 and _rejected(too_long)["length"] == 1
    assert _scan(tmp_path, index, _pair(*OVERRUN), max_frag_length=46).n_fragments == 1


def test_spliced_mates_that_run_past_each_other_are_one_molecule(index, tmp_path):
    """F `46M900N50M` at 1154 and R `50M900N46M` at 1150 both cross tK's intron [1200,2100), and each runs 4 past the
    other's 5' end: one molecule over [1150,2150), L 100 with the intron removed once."""
    records = _pair("overrun_spliced", 1154, "46M900N50M", 1150, "50M900N46M", xs="+")
    stats = _scan(tmp_path, index, records, max_frag_length=100)
    assert stats.n_fragments == 1 and sum(_rejected(stats).values()) == 0


NESTED = [
    ("gT", "T", [(1000, 1100), (3000, 3060), (4000, 4100)]),
    ("gU", "U", [(3100, 3400)]),
]


@pytest.fixture(scope="module")
def index_nested(tmp_path_factory):
    """T over [1000,4100), and U's one exon [3100,3400) inside T's second intron [3060,4000)."""
    return _index(tmp_path_factory.mktemp("pairs_nested"), NESTED)


def test_a_pair_the_rule_accepts_is_quantified_whole(index_nested, tmp_path):
    """F `50M1900N50M` at 1050 crosses T's first intron [1100,3000), and R lies at [3150,3250) in U's exon, inside T's
    second intron: 300 bases along the genome once the junction is removed, and T holds the junction and overlaps both
    mates. The pair is one molecule, T's RNA spliced at its first intron and not at its second, and the EM takes it
    whole, though its blocks' transcripts do not meet (T's exons and U's) and its outer span is 2,200 bases."""
    records = _pair("nested", 1050, "50M1900N50M", 3150, "100M", xs="+")
    stats = _scan(tmp_path, index_nested, records)
    assert stats.n_fragments == 1 and sum(_rejected(stats).values()) == 0
    assert stats.n_with_exon == 1, "a pair the pair rule accepts was set aside as a chimera"


def test_a_pair_too_long_along_every_path_is_withheld(index, tmp_path):
    """Mates 2,100 bases apart over no annotated intron: L 2,100 against max_frag_length 1,000."""
    stats = _scan(tmp_path, index, _pair("far", 6000, "100M", 8000, "100M"))
    assert stats.n_fragments == 0
    assert _rejected(stats)["length"] == 1
    assert stats.n_withheld_unpaired == 1


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


def test_a_transcript_between_the_mates_is_no_witness(index, tmp_path):
    """F at [800,900) and R at [2500,2600), both outside tK [1000,2400): 1,800 along the genome. Cutting tK's intron
    [1200,2100), which lies between them, would give 900, but tK overlaps neither mate, so it offers no path: a
    witness overlaps both mates (owner, 2026-10-19). Rejected for length."""
    stats = _scan(tmp_path, index, _pair("between", 800, "100M", 2500, "100M"))
    assert stats.n_fragments == 0, "a transcript that misses both mates witnessed their path"
    assert _rejected(stats)["length"] == 1
    assert stats.n_withheld_unpaired == 1


def test_a_long_spliced_pair_whose_holders_miss_a_mate_is_too_long(index_t3, tmp_path):
    """F `10M100N10M` at 10100 crosses [10110,10210), which T1 and T3 hold; R `10M` at [11850,11860) lies past both
    transcripts' ends: 1,660 along the genome. T3's intron [10300,11300) between the mates would cut it to 660, but T3
    does not overlap R, so it offers no path. Rejected for length, which is judged before the junctions, so not
    withheld as unsupported."""
    stats = _scan(
        tmp_path, index_t3, _pair("long_spliced", 10100, "10M100N10M", 11850, "10M", xs="+")
    )
    assert stats.n_fragments == 0, "a holder that misses a mate witnessed the path"
    assert _rejected(stats)["length"] == 1 and _rejected(stats)["unsupported"] == 0


def test_an_invalid_hit_is_dropped_and_its_template_kept(index, tmp_path):
    """R1 at [1010,1110); R2 at [1100,1200), the control pair, and at [900,1000), where the two face away from each
    other: that hit is dropped, the template kept as one fragment, a unique mapper, not a multimapper whose second place
    does not exist."""
    records = [
        _record("mixed", R1_FWD, 1010, "100M", 1100, nh=2),
        _record("mixed", R2_REV, 1100, "100M", 1010, nh=2),
        _record("mixed", R2_REV | SECONDARY, 900, "100M", 1010, nh=2),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 1
    assert _rejected(stats)["outward"] == 1
    assert stats.n_withheld_unpaired == 0
    assert stats.n_multimapper_groups == 0, "an invalid place still counted toward the multiplicity"


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
    assert stats.n_withheld_unpaired == 1


def test_a_place_with_no_partner_is_no_hit(index, tmp_path):
    """A secondary R1 flagged mate-unmapped is one more R1 place: the search tries it with R2 at [1100,1200), where the
    two face away from each other, so it is no hit, and the primary pair stays, calibration's one placement. The search
    reads geometry, not the aligner's mate flags."""
    records = _pair("nohi", 1010, "100M", 1100, "100M", nh=2) + [
        _record("nohi", R1_FWD | SECONDARY | 0x8, 1500, "100M", 1500, nh=2)
    ]
    stats, payload = _scan_full(tmp_path, index, records)
    assert payload.qc.deferred_multiple_placements == 0, "the lone secondary became a placement"
    assert stats.n_fragments == 1
    assert _rejected(stats)["outward"] == 1 and _rejected(stats)["missing_mate"] == 0
    assert stats.n_withheld_unpaired == 0


# ── pairing: every R1 place with every R2 place, by geometry alone ──────────────────────────────────────────────────


def test_the_search_finds_a_partner_outside_the_gene(index, tmp_path):
    """R1 at [2250,2350) in tK's last exon; R2 at [1500,1600) in tK's intron and at [2400,2500), just past tK's end,
    the true partner (L 250). A pairing by shared transcripts matches R1 with the intron place through gK's nascent
    span and never tries the true partner; the search tries both, and keeps the one molecule."""
    records = [
        _record("outside", R1_FWD, 2250, "100M", 1500, nh=2),
        _record("outside", R2_REV, 1500, "100M", 2250, nh=2),
        _record("outside", R2_REV | SECONDARY, 2400, "100M", 2250, nh=2),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 1, "the true pair was never tried"
    assert _rejected(stats)["outward"] == 1 and stats.n_multimapper_groups == 0


def test_a_record_repeated_for_each_pair_is_one_location(index, tmp_path):
    """As STAR writes a template whose R1 has one place and R2 two: R1's record twice, once beside each R2 place. The
    two R1 records are one alignment, so one place, and the template has two hits, not four."""
    records = [
        _record("repeat", R1_FWD, 1010, "100M", 1100, nh=2),
        _record("repeat", R2_REV, 1100, "100M", 1010, nh=2),
        _record("repeat", R1_FWD | SECONDARY, 1010, "100M", 1050, nh=2),
        _record("repeat", R2_REV | SECONDARY, 1050, "100M", 1010, nh=2),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_multimapper_groups == 1
    assert stats.n_multimapper_alignments == 2, "a repeated record paired again"


def test_a_location_paired_once_is_tried_with_every_partner(index, tmp_path):
    """R1 at [2250,2350) in tK's last exon; R2 at [2300,2400) in that exon (L 150) and at [2400,2500) past tK's end
    (L 250): two molecules, so a multimapper. A pairing by shared transcripts matches R1 with the exon place and leaves
    the other as a one-mate hit."""
    records = [
        _record("twice", R1_FWD, 2250, "100M", 2300, nh=2),
        _record("twice", R2_REV, 2300, "100M", 2250, nh=2),
        _record("twice", R2_REV | SECONDARY, 2400, "100M", 2250, nh=2),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_multimapper_groups == 1 and stats.n_multimapper_alignments == 2
    assert _rejected(stats)["missing_mate"] == 0


def test_a_same_strand_combination_is_no_candidate(index, tmp_path):
    """R1 forward at [1010,1110) and reverse at [3010,3110), R2 reverse at [1100,1200): the two reverse places are no
    pair to judge, so the template is one fragment with no rejection counted."""
    records = [
        _record("strands", R1_FWD, 1010, "100M", 1100, nh=2),
        _record("strands", (R1_FWD | SECONDARY | 0x10) & ~0x20, 3010, "100M", 1100, nh=2),
        _record("strands", R2_REV, 1100, "100M", 1010, nh=2),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 1
    assert _rejected(stats) == dict.fromkeys(REASONS, 0), "a same-strand combination was judged"


def test_mates_on_two_references_are_judged_once(index_two_refs, tmp_path):
    """R1 at two places on chr1 and R2 on chr2, both references in the index: no pair is a candidate, so the first
    place of each is judged, once, and the template is withheld; its reason is counted once, not per combination."""
    records = [
        _record("trans", R1_FWD, 1010, "100M", 1100, nh=2),
        _record("trans", R1_FWD | SECONDARY, 5010, "100M", 1100, nh=2),
        _record("trans", R2_REV, 1100, "100M", 1010, nh=2, ref=1),
    ]
    stats = _scan(tmp_path, index_two_refs, records)
    assert stats.n_fragments == 0
    assert _rejected(stats)["reference"] == 1 and stats.n_withheld_unpaired == 1


def test_an_alias_is_written_with_its_alignments_hit(index, tmp_path):
    """The repeated R1 record of the two-place template, in the annotated BAM: written with the hit of the alignment it
    repeats, not swept up at the end of its template as an intergenic pass-through."""
    records = [
        r for k in range(200) for r in _pair(f"ok{k:03d}", 1010 + k, "100M", 1100 + k, "100M")
    ]
    records += [
        _record("repeat", R1_FWD, 1010, "100M", 1100, nh=2),
        _record("repeat", R2_REV, 1100, "100M", 1010, nh=2),
        _record("repeat", R1_FWD | SECONDARY, 1010, "100M", 1050, nh=2),
        _record("repeat", R2_REV | SECONDARY, 1050, "100M", 1010, nh=2),
    ]
    out = tmp_path / "annotated.bam"
    config = PipelineConfig(
        em=EMConfig(seed=7), scan=BamScanConfig(total_threads=1), annotated_bam_path=out
    )
    run_pipeline(_bam(tmp_path, records), index, config=config)
    with pysam.AlignmentFile(str(out), "rb") as bam:
        r1 = [
            (r.get_tag("ZF"), r.get_tag("ZT"), r.get_tag("ZC"))
            for r in bam
            if r.query_name == "repeat" and r.is_read1
        ]
    assert len(r1) == 2
    assert r1[0] == r1[1], "the alias was not written with its alignment's hit"
    assert r1[0][0] & 0x01 and r1[0][0] != 0x25


# ── references outside the index: Rigel can place nothing there ─────────────────────────────────────────────────


def test_a_pair_off_the_index_is_withheld_not_counted_gdna(index, tmp_path):
    """Both mates on chr2, a reference the BAM names and the index does not hold (a spike-in, alternate or decoy
    contig). It has no partition, no transcript and no length along one, so the hit is dropped and the template
    withheld. It is no intergenic fragment, which `summary.json` counts as gDNA."""
    stats = _scan(tmp_path, index, _pair("spike", 1010, "100M", 1100, "100M", ref=1))
    assert stats.n_fragments == 0 and stats.n_intergenic == 0, (
        "a pair off the index was counted intergenic"
    )
    assert _rejected(stats)["off_index"] == 1
    assert stats.n_withheld_off_index == 1 and stats.n_withheld_unpaired == 0


def test_mates_on_two_references_off_the_index_are_two_references(index, tmp_path):
    """R1 on chr2 and R2 on chr3, neither in the index: mates on two references, as on any two, not one fragment merged
    by coordinate."""
    records = [
        _record("trans_off", R1_FWD, 1010, "100M", 1100, ref=1),
        _record("trans_off", R2_REV, 1100, "100M", 1010, ref=2),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 0 and stats.n_intergenic == 0, (
        "the two references were merged into one fragment"
    )
    assert _rejected(stats)["reference"] == 1 and stats.n_withheld_unpaired == 1


def test_a_place_off_the_index_is_no_place_of_the_template(index, tmp_path):
    """One pair on chr1 (L 190), and the same pair again on chr2, off the index: that place is dropped and counted, and
    the template is placed on chr1 alone, a unique mapper."""
    records = _pair("decoy", 1010, "100M", 1100, "100M", nh=2) + _pair(
        "decoy", 1010, "100M", 1100, "100M", nh=2, secondary=True, ref=1
    )
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 1 and stats.n_intergenic == 0
    assert stats.n_multimapper_groups == 0, (
        "the place off the index was one of the template's places"
    )
    assert _rejected(stats)["off_index"] == 1


def test_the_annotated_bam_labels_a_template_off_the_index(index, tmp_path):
    """A template off the index is written unresolved (ZF 0), class `off_index`, not as intergenic gDNA (ZF 0x25)."""
    records = [
        r for k in range(200) for r in _pair(f"ok{k:03d}", 1010 + k, "100M", 1100 + k, "100M")
    ]
    records += _pair("spike", 1010, "100M", 1100, "100M", ref=1)
    out = tmp_path / "annotated.bam"
    config = PipelineConfig(
        em=EMConfig(seed=7), scan=BamScanConfig(total_threads=1), annotated_bam_path=out
    )
    run_pipeline(_bam(tmp_path, records), index, config=config)
    with pysam.AlignmentFile(str(out), "rb") as bam:
        tags = {(r.query_name, r.is_read1): (r.get_tag("ZF"), r.get_tag("ZC")) for r in bam}
    assert tags[("spike", True)] == tags[("spike", False)] == (0, "off_index")


# ── chimeric alignments: withheld ──────────────────────────────────────────────────────────────────────────────


def _spanning_chimera(qname="chimera"):
    """A spanning chimera as STAR writes it (`--chimOutType WithinBAM HardClip`): R1's first 77 bases at 1010, written
    as the supplementary record, and its last 23 at 1175 with R2 at 1187, the representative pair. The SA tags name
    each other. Read as a pair, the representative's walked clip makes it a molecule of L 189 along the genome."""
    return [
        _record(
            qname, R1_FWD | SUPPLEMENTARY, 1010, "77M23H", 1187, sa="chr1,1176,+,77S23M,255,0;"
        ),
        _record(qname, R1_FWD, 1175, "77S23M", 1187, sa="chr1,1011,+,77M23H,255,0;"),
        _record(qname, R2_REV, 1187, "100M", 1175),
    ]


def test_a_spanning_chimera_is_withheld(index, tmp_path):
    """STAR writes a chimeric read in place of its linear alignment, and Rigel does not quantify chimeras: the template
    is withheld as chimeric. It is no fragment, chimeric or intergenic, and its pieces are not joined into one."""
    stats = _scan(tmp_path, index, _spanning_chimera())
    assert stats.n_fragments == 0 and stats.n_intergenic == 0
    assert _rejected(stats)["chimeric"] == 1 and stats.n_withheld_chimeric == 1


def test_a_chimeric_mate_keeps_a_self_contained_alignment(index, tmp_path):
    """The chimera above, and R1 also aligned whole at 1010 (a secondary record with no SA): the chimeric hit is dropped
    and the template is placed by the whole alignment, a unique mapper."""
    records = _spanning_chimera("alt") + [
        _record("alt", R1_FWD | SECONDARY, 1010, "100M", 1187, nh=2)
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 1 and stats.n_multimapper_groups == 0
    assert _rejected(stats)["chimeric"] == 1


def test_a_secondary_chimeras_pieces_join_no_other_alignment(index, tmp_path):
    """R1 whole at 1010 with R2 at 1187 (L 277), and R1 also aligned as a secondary chimera: 60 bases at 5000 and its
    supplementary record at 15000. The pieces were glued onto R1's primary location; now they join nothing, so the
    fragment spans the pair alone."""
    records = [
        _record("sec", R1_FWD, 1010, "100M", 1187, nh=2),
        _record("sec", R2_REV, 1187, "100M", 1010, nh=2),
        _record(
            "sec", R1_FWD | SECONDARY, 5000, "60M40S", 1187, nh=2, sa="chr1,15001,+,60H40M,0,0;"
        ),
        _record(
            "sec",
            R1_FWD | SECONDARY | SUPPLEMENTARY,
            15000,
            "60H40M",
            1187,
            nh=2,
            sa="chr1,5001,+,60M40S,0,0;",
        ),
    ]
    stats, footprint = _scan_buffered(tmp_path, index, records, "genomic_footprint")
    assert footprint == [277], (
        "the secondary chimera's supplementary record joined the primary location's fragment"
    )
    assert stats.n_multimapper_groups == 0


def test_a_chain_without_sa_withholds_its_mate(index, tmp_path):
    """R1's supplementary record names no partner, and no R1 record carries SA: which alignment it belongs to is
    ambiguous, so every R1 location is withheld with it."""
    records = [
        _record("nosa", R1_FWD, 1010, "60M40S", 1187),
        _record("nosa", R1_FWD | SUPPLEMENTARY, 5000, "60H40M", 1187),
        _record("nosa", R2_REV, 1187, "100M", 1010),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 0 and stats.n_withheld_chimeric == 1


@pytest.mark.parametrize(
    "flag, cigar, pos, usable",
    [
        pytest.param(R1_FWD, "10H90M", 1020, False, id="forward-5prime"),
        pytest.param(R1_FWD, "90M10H", 1010, True, id="forward-3prime"),
        pytest.param(R1_REV, "90M10H", 1197, False, id="reverse-5prime"),
        pytest.param(R1_REV, "10H90M", 1197, True, id="reverse-3prime"),
    ],
)
def test_a_hard_clipped_5prime_end_is_withheld(index, tmp_path, flag, cigar, pos, usable):
    """A hard clip removes the bases from the record, and at a read's 5' end they are the molecule's end: the mate's
    alignment is not wholly in its record, so it is read as chimeric. A hard clip at the 3' end ends no molecule.
    The mate is R2 at 1100 (forward when R1 is reverse) or 1187."""
    mate_flag = R2_REV if flag == R1_FWD else R2_FWD
    mate_pos = 1187 if flag == R1_FWD else 1100
    records = [
        _record("hard", flag, pos, cigar, mate_pos),
        _record("hard", mate_flag, mate_pos, "100M", pos),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == (1 if usable else 0)
    assert stats.n_withheld_chimeric == (0 if usable else 1)


def test_the_annotated_bam_labels_a_chimeric_template(index, tmp_path):
    """A withheld chimeric template's records, its supplementary record too, are written unresolved (ZF 0), class
    `chimeric`. A placed template's supplementary record is written with the same mate's primary location, never swept
    up as intergenic gDNA (ZF 0x25)."""
    records = [
        r for k in range(200) for r in _pair(f"ok{k:03d}", 1010 + k, "100M", 1100 + k, "100M")
    ]
    records += _spanning_chimera()
    records += _spanning_chimera("alt") + [
        _record("alt", R1_FWD | SECONDARY, 1010, "100M", 1187, nh=2)
    ]
    out = tmp_path / "annotated.bam"
    config = PipelineConfig(
        em=EMConfig(seed=7), scan=BamScanConfig(total_threads=1), annotated_bam_path=out
    )
    run_pipeline(_bam(tmp_path, records), index, config=config)
    with pysam.AlignmentFile(str(out), "rb") as bam:
        rows = [(r.query_name, r.is_supplementary, r.get_tag("ZF"), r.get_tag("ZC")) for r in bam]
    chimera = [(zf, zc) for q, _supp, zf, zc in rows if q == "chimera"]
    assert len(chimera) == 3 and set(chimera) == {(0, "chimeric")}
    alt_supp = [zf for q, supp, zf, _zc in rows if q == "alt" and supp]
    assert len(alt_supp) == 1 and alt_supp[0] != 0x25 and alt_supp[0] & 0x01


# ── a record's integer tags: NM and NH ─────────────────────────────────────────────────────────────────────────


def _scan_nm(base, index, records):
    """The scan's stats, and the NM of each fragment it buffered."""
    stats, _strand, buffer, _payload = scan_and_buffer(
        _bam(base, records), index, BamScanConfig(total_threads=1)
    )
    return stats, [int(x) for chunk in buffer.iter_chunks_consuming() for x in chunk.nm]


def test_a_record_without_nm_is_counted_and_read_as_zero(index, tmp_path):
    """A record without NM is counted, and its edit distance read as 0: a disclosed approximation, which can favor an
    alignment of unknown NM over one whose NM is known. STAR writes NM only when asked for it."""
    records = [
        _record("no_nm", R1_FWD, 1010, "100M", 1100, nm=3),
        _record("no_nm", R2_REV, 1100, "100M", 1010, nm=None),
    ]
    stats, nm = _scan_nm(tmp_path, index, records)
    assert nm == [3]
    assert stats.n_nm_tag_missing == 1 and stats.n_nm_tag_invalid == 0


def test_tags_of_another_type_are_counted_and_read_as_absent(index, tmp_path):
    """NM and NH are integers. One held as a string is an input defect: counted, and read as absent (NM 0, NH 1)."""
    records = _pair("typed", 1010, "100M", 1100, "100M")
    records[0].set_tag("NM", "2", "Z")
    records[1].set_tag("NH", "1", "Z")
    stats, nm = _scan_nm(tmp_path, index, records)
    assert nm == [0] and stats.unique == 1
    assert stats.n_nm_tag_invalid == 1 and stats.n_nh_tag_invalid == 1
    assert stats.n_nm_tag_missing == 0 and stats.n_nh_tag_missing == 0


def test_a_record_without_nh_is_counted(index, tmp_path):
    """An aligner that writes no NH (bwa) marks its other alignments secondary instead: each record is counted, and its
    NH reads 1."""
    records = [
        _record("no_nh", R1_FWD, 1010, "100M", 1100, nh=None),
        _record("no_nh", R2_REV, 1100, "100M", 1010, nh=None),
    ]
    stats, _nm = _scan_nm(tmp_path, index, records)
    assert stats.unique == 1
    assert stats.n_nh_tag_missing == 2


def test_records_of_one_alignment_that_disagree_on_nm_keep_the_larger(index, tmp_path):
    """R1's alignment written twice, as STAR writes it beside each partner, with NM 1 and NM 4: an input defect, counted
    once for the alignment, which keeps the larger NM. R2 at 1100 and at 1600 (NM 0) make two hits."""
    records = [
        _record("nm_conflict", R1_FWD, 1010, "100M", 1100, nh=2, nm=1),
        _record("nm_conflict", R2_REV, 1100, "100M", 1010, nh=2),
        _record("nm_conflict", R1_FWD | SECONDARY, 1010, "100M", 1600, nh=2, nm=4),
        _record("nm_conflict", R2_REV | SECONDARY, 1600, "100M", 1010, nh=2),
    ]
    stats, nm = _scan_nm(tmp_path, index, records)
    assert nm == [4, 4], "the alignment did not keep the larger NM"
    assert stats.n_nm_conflicting == 1


# ── multiplicity: a template's distinct valid hits ──────────────────────────────────────────────────


def test_two_distinct_valid_hits_are_a_multimapper(index, tmp_path):
    """R1 at [1010,1110); R2 at [1100,1200) (L 190) and at [1600,1700) (L 690): two molecules."""
    records = [
        _record("two", R1_FWD, 1010, "100M", 1100, nh=2),
        _record("two", R2_REV, 1100, "100M", 1010, nh=2),
        _record("two", R2_REV | SECONDARY, 1600, "100M", 1010, nh=2),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_multimapper_groups == 1 and stats.n_multimapper_alignments == 2


def test_a_spliced_and_a_contiguous_reading_stay_two_hits(index, tmp_path):
    """R1's record twice, as STAR writes it beside each partner, one place; R2 from 1150, contiguous (`100M`, 50 bases
    into tK's intron) and spliced across tK's intron [1200,2100) (`50M900N50M`). One start, two readings of L 240: the
    competition between them is the point (owner, 2026-10-18), so nothing may merge them into one hit."""
    records = [
        _record("compete", R1_FWD, 1010, "100M", 1150, nh=2),
        _record("compete", R2_REV, 1150, "100M", 1010, nh=2),
        _record("compete", R1_FWD | SECONDARY, 1010, "100M", 1150, nh=2),
        _record("compete", R2_REV | SECONDARY, 1150, "50M900N50M", 1010, nh=2, xs="+"),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_multimapper_groups == 1 and stats.n_multimapper_alignments == 2


def test_an_nh_above_the_records_reported_is_counted(index, tmp_path):
    """NH 3, and each mate reported in two records: a truncated output, counted. The places reported are the template's
    multiplicity: R1 at 1010 pairs with R2 at 1100 and 1600, R1 at 1500 with R2 at 1600."""
    records = _pair("trunc", 1010, "100M", 1100, "100M", nh=3) + _pair(
        "trunc", 1500, "100M", 1600, "100M", nh=3, secondary=True
    )
    stats = _scan(tmp_path, index, records)
    assert stats.n_nh_above_reported == 1
    assert stats.n_multimapper_groups == 1 and stats.n_multimapper_alignments == 3


def test_the_aligners_pairing_is_not_read(index, tmp_path):
    """STAR's HI groups pair R1 at 1010 with R2 at 10290 (too long) and R1 at 10220 with R2 at 1100 (facing away); the
    true molecules are the other two combinations, in tK (L 190) and in T1 (L 170). The search reads geometry, not HI,
    so it keeps both: a multimapper, where pairing by HI would withhold the template."""
    records = [
        _record("hi", R1_FWD, 1010, "100M", 10290, nh=2, hi=1),
        _record("hi", R2_REV, 10290, "100M", 1010, nh=2, hi=1),
        _record("hi", R1_FWD | SECONDARY, 10220, "100M", 1100, nh=2, hi=2),
        _record("hi", R2_REV | SECONDARY, 1100, "100M", 10220, nh=2, hi=2),
    ]
    stats = _scan(tmp_path, index, records)
    assert stats.n_fragments == 1, "the pairing followed HI"
    assert stats.n_multimapper_groups == 1 and stats.n_multimapper_alignments == 2
    assert _rejected(stats)["length"] == 1 and _rejected(stats)["outward"] == 1


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


# ── the annotated-path rule: one annotated transcript holds every junction a hit observes ─────────────────────────────


def _withheld(stats):
    return {
        "unpaired": stats.n_withheld_unpaired,
        "unannotated_splice": stats.n_withheld_unannotated_splice,
    }


def test_a_novel_junction_is_withheld_and_a_contiguous_reading_kept(index, tmp_path):
    """HI 1 splices F across [1060,1160), inside gK, where no transcript has an intron; HI 2 reads the same pair
    contiguous. The spliced hit is dropped and the contiguous one kept: one place, so a unique mapper."""
    records = _pair("novel", 1010, "50M100N50M", 1300, "100M", nh=2, hi=1, xs="+") + _pair(
        "novel", 1010, "100M", 1300, "100M", nh=2, hi=2, secondary=True
    )
    stats = _scan(tmp_path, index, records)
    assert _rejected(stats)["novel_junction"] == 1
    assert stats.n_fragments == 1 and stats.n_multimapper_groups == 0
    assert _withheld(stats) == {"unpaired": 0, "unannotated_splice": 0}


def test_an_annotated_and_a_novel_junction_are_a_novel_junction(index, tmp_path):
    """F crosses T1's intron [10110,10210), R crosses [10240,10340), which no transcript has. One matching junction
    once made the hit annotated-spliced; now it is withheld, an unannotated splice."""
    stats = _scan(tmp_path, index, _pair("half", 10100, "10M100N10M", 10220, "20M100N20M", xs="+"))
    assert _rejected(stats)["novel_junction"] == 1
    assert stats.n_fragments == 0
    assert _withheld(stats) == {"unpaired": 0, "unannotated_splice": 1}


def test_annotated_junctions_no_one_transcript_holds_are_unsupported(index, index_t3, tmp_path):
    """F crosses T1's intron [10110,10210), R T2's [10300,11300): each is annotated, and no transcript holds both
    until T3 does."""
    records = _pair("combo", 10100, "10M100N10M", 10290, "10M1000N10M", xs="+")
    without = _scan(tmp_path, index, records)
    assert _rejected(without)["unsupported"] == 1 and without.n_fragments == 0
    assert _withheld(without)["unannotated_splice"] == 1
    with_t3 = _scan(tmp_path, index_t3, records)
    assert sum(_rejected(with_t3).values()) == 0 and with_t3.n_fragments == 1


def test_a_spliced_hit_needs_a_holder_overlapping_both_mates(index, tmp_path):
    """F splices across tK's intron [1200,2100) (`50M900N50M` from 1150); R lies wholly past tK's end, at [2450,2550).
    Short along the genome (L 500), so no witness is needed for its length, but tK, its one holder, does not overlap R:
    no one transcript explains both mates, so the hit is unsupported (owner, 2026-10-19). R at [2350,2450) runs 50
    bases past tK's end and still overlaps it: overlap, not containment, so that hit is usable, its extra bases the
    EM's overhang."""
    past = _scan(tmp_path, index, _pair("past", 1150, "50M900N50M", 2450, "100M", xs="+"))
    assert past.n_fragments == 0, "a holder that misses a mate made the hit usable"
    assert _rejected(past)["unsupported"] == 1
    assert _withheld(past) == {"unpaired": 0, "unannotated_splice": 1}
    overhang = _scan(tmp_path, index, _pair("overhang", 1150, "50M900N50M", 2350, "100M", xs="+"))
    assert sum(_rejected(overhang).values()) == 0 and overhang.n_fragments == 1


def test_mates_whose_junctions_overlap_unequal_contradict(index, tmp_path):
    """F splices at T1's donor 10110, R at 10115, both to 10210: the mates disagree on the donor, so no molecule
    carries both junctions."""
    stats = _scan(tmp_path, index, _pair("donor", 10100, "10M100N10M", 10105, "10M95N15M", xs="+"))
    assert _rejected(stats)["contradictory"] == 1 and stats.n_fragments == 0


def test_motifs_naming_both_strands_contradict(index, tmp_path):
    """T1's junction read on '+' by F and on '-' by R; and, in a second template, a record whose XS and ts
    disagree, whose own motif is ambiguous. Either way the hit names both strands."""
    mates = [
        _record("opposite", R1_FWD, 10100, "10M100N10M", 10105, xs="+"),
        _record("opposite", R2_REV, 10105, "5M100N10M", 10100, xs="-"),
    ]
    stats = _scan(tmp_path, index, mates)
    assert _rejected(stats)["contradictory"] == 1 and stats.n_fragments == 0
    tags = [
        _record("tags", R1_FWD, 10100, "10M100N10M", 10300, xs="+", ts="-"),
        _record("tags", R2_REV, 10300, "10M", 10100),
    ]
    stats = _scan(tmp_path, index, tags)
    assert _rejected(stats)["contradictory"] == 1 and stats.n_fragments == 0


def test_a_junction_is_held_on_its_motifs_strand_or_either_without_one(index, tmp_path):
    """T1's intron [10110,10210) is on '+': read on '-' it is a junction no transcript has; read with no tag it
    matches at either strand."""
    opposite = _scan(tmp_path, index, _pair("minus", 10100, "10M100N10M", 10300, "10M", xs="-"))
    assert _rejected(opposite)["novel_junction"] == 1 and opposite.n_fragments == 0
    untagged = _scan(tmp_path, index, _pair("untagged", 10100, "10M100N10M", 10300, "10M"))
    assert sum(_rejected(untagged).values()) == 0 and untagged.n_fragments == 1


def test_the_rule_judges_a_hit_the_pair_rule_does_not(index, tmp_path):
    """F spliced across [1060,1160) inside gK, where no transcript has an intron, and its mate with no aligned base
    (`100S`): the pair rule leaves such a mate as it was, but no transcript holds the junction, so the hit is withheld
    rather than counted intergenic RNA."""
    records = [
        _record("noalign", R1_FWD, 1010, "50M100N50M", 1300, xs="+"),
        _record("noalign", R2_REV, 1300, "100S", 1010),
    ]
    stats = _scan(tmp_path, index, records)
    assert _rejected(stats)["novel_junction"] == 1
    assert stats.n_fragments == 0 and stats.n_intergenic == 0
    assert _withheld(stats)["unannotated_splice"] == 1


def test_a_spliced_pair_outside_every_gene_is_withheld(index, tmp_path):
    stats = _scan(tmp_path, index, _pair("intergenic", 15000, "50M100N50M", 15200, "100M", xs="+"))
    assert _rejected(stats)["novel_junction"] == 1
    assert stats.n_fragments == 0 and stats.n_intergenic == 0
    assert _withheld(stats) == {"unpaired": 0, "unannotated_splice": 1}


def test_the_length_is_judged_before_the_junctions(index, tmp_path):
    """A pair 2,000 bases long with a novel junction is no molecule: rejected for its length, so the rule withholds
    only molecules, each with a fragment length."""
    stats = _scan(tmp_path, index, _pair("long", 15000, "50M100N50M", 17000, "100M", xs="+"))
    assert _rejected(stats)["length"] == 1 and _rejected(stats)["novel_junction"] == 0
    assert _withheld(stats) == {"unpaired": 1, "unannotated_splice": 0}


def test_a_template_with_a_molecule_the_rule_withheld_is_an_unannotated_splice(index, tmp_path):
    """HI 1 is too long, HI 2 a molecule with a novel junction: the template is an unannotated splice, not
    unpaired."""
    records = _pair("both", 15000, "50M100N50M", 17000, "100M", nh=2, hi=1, xs="+") + _pair(
        "both", 15000, "50M100N50M", 15200, "100M", nh=2, hi=2, secondary=True, xs="+"
    )
    stats = _scan(tmp_path, index, records)
    assert _rejected(stats)["length"] == 1 and _rejected(stats)["novel_junction"] == 1
    assert _withheld(stats) == {"unpaired": 0, "unannotated_splice": 1}


@pytest.fixture(scope="module")
def index_blacklist(tmp_path_factory):
    """The toy index with a splice blacklist of two junctions no transcript has, each at anchors up to 50."""
    from _index_builder import rebuild_with_splice_blacklist

    base = tmp_path_factory.mktemp("pairs_blacklist")
    _index(base)
    rebuild_with_splice_blacklist(
        base / "idx",
        pd.DataFrame(
            {
                "ref": ["chr1", "chr1"],
                "start": [15050, 10240],
                "end": [15150, 10340],
                "max_anchor_left": [50, 50],
                "max_anchor_right": [50, 50],
            }
        ),
    )
    return TranscriptIndex.load(str(base / "idx"))


def test_the_rule_reads_the_junctions_the_blacklist_keeps(index_blacklist, tmp_path):
    """The blacklist erases a junction before the rule reads it: a pair whose only junction is blacklisted is
    unspliced, and a pair with a blacklisted junction and T1's is decided by T1's."""
    alone = _scan(
        tmp_path, index_blacklist, _pair("alone", 15000, "50M100N50M", 15200, "100M", xs="+")
    )
    assert sum(_rejected(alone).values()) == 0 and alone.n_fragments == 1
    kept = _scan(
        tmp_path, index_blacklist, _pair("kept", 10100, "10M100N10M", 10220, "20M100N20M", xs="+")
    )
    assert sum(_rejected(kept).values()) == 0 and kept.n_fragments == 1


def test_the_annotated_bam_labels_an_unannotated_splice(index, tmp_path):
    """A template the rule withheld is written unresolved (ZF 0), class `unannotated_splice`."""
    records = [
        r for k in range(200) for r in _pair(f"ok{k:03d}", 1010 + k, "100M", 1100 + k, "100M")
    ]
    records += _pair("intergenic", 15000, "50M100N50M", 15200, "100M", xs="+")
    out = tmp_path / "annotated.bam"
    config = PipelineConfig(
        em=EMConfig(seed=7), scan=BamScanConfig(total_threads=1), annotated_bam_path=out
    )
    run_pipeline(_bam(tmp_path, records), index, config=config)
    with pysam.AlignmentFile(str(out), "rb") as bam:
        tags = {(r.query_name, r.is_read1): (r.get_tag("ZF"), r.get_tag("ZC")) for r in bam}
    assert tags[("intergenic", True)] == tags[("intergenic", False)] == (0, "unannotated_splice")
