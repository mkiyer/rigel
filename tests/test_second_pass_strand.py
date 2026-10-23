"""The second pass's strand factor for a placement whose junction motif was sequenced.

A held record's placements are compared by the product of their factors, so every placement's strand factor must
be the probability of its observed orientation under its own reading. A placement with a sequenced motif is RNA on
that motif's strand under every hypothesis, so its factor is ``p`` or ``1 - p``. It was left at 1 on the grounds that
it is constant within the placement, which holds; but across a multimapper's placements it then outweighed an
unspliced hit's ½ even on the wrong strand, where its true factor is ``1 - p``. Specification:
``second_pass.strand_terms``.
"""

from __future__ import annotations

import numpy as np
import pytest
from _index_builder import write_fasta

from rigel.calibration.fl import build_fl_models
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.splice_graph import build_sj_arrays
from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.second_pass import score_held_fragments, strand_terms
from rigel.types import Strand

GENOME = 6_000
N_SENSE = (
    50  # uniquely mapped spliced fragments that set the library's strand: R1 on the motif's strand
)
N_ANTISENSE = 2  # and a few on the other strand, so the library is stranded but not perfectly
N_MULTI = 20
_M, _N = 0, 3


def _gtf() -> str:
    rows = []
    for ref in ("chr1", "chr2"):
        rows += [
            f'{ref}\tt\texon\t1001\t1600\t.\t+\t.\tgene_id "gW_{ref}"; transcript_id "tW_{ref}";\n',
            f'{ref}\tt\texon\t2001\t2200\t.\t+\t.\tgene_id "gW_{ref}"; transcript_id "tW_{ref}";\n',
        ]
    return "".join(rows)


@pytest.fixture(scope="module")
def scanned(tmp_path_factory):
    """``N_SENSE`` unique fragments spliced across chr1's annotated junction with R1 on the motif's strand, and
    ``N_MULTI`` multimappers whose first hit is spliced across that junction with R1 on the OTHER strand and whose
    second hit is unspliced, contained in chr2's first exon."""
    import pysam

    base = tmp_path_factory.mktemp("strand_placements")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    write_fasta(fasta, {"chr1": GENOME, "chr2": GENOME})
    gtf.write_text(_gtf())
    TranscriptIndex.build(str(fasta), str(gtf), str(base / "idx"), write_tsv=False)
    index = TranscriptIndex.load(str(base / "idx"))

    def read(
        qname, ref_id, pos, mate_pos, is_r1, reverse, cigar, nh=1, hi=None, secondary=False, xs=None
    ):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.reference_id = ref_id
        a.reference_start = pos
        a.mapping_quality = 60
        a.flag = (
            0x1
            | 0x2
            | (0x40 if is_r1 else 0x80)
            | (0x10 if reverse else 0)
            | (0x20 if not reverse else 0)
            | (0x100 if secondary else 0)
        )
        a.cigar = cigar
        n = sum(length for op, length in cigar if op == _M)
        a.query_sequence = "A" * n
        a.query_qualities = pysam.qualitystring_to_array("I" * n)
        a.next_reference_id = ref_id
        a.next_reference_start = mate_pos
        tags = [("NH", nh, "i")] + ([("HI", hi, "i")] if hi is not None else [])
        if xs is not None:
            tags.append(("XS", xs, "A"))
        a.set_tags(tags)
        return a

    spliced = [(_M, 50), (_N, 400), (_M, 50)]  # [1550, 1600) | intron [1600, 2000) | [2000, 2050)
    reads = []
    for k in range(N_SENSE):  # R1 forward and spliced on a + motif: sense
        reads += [
            read(f"s_{k}", 0, 1550, 2080, True, False, spliced, xs="+"),
            read(f"s_{k}", 0, 2080, 1550, False, True, [(_M, 100)]),
        ]
    for k in range(N_ANTISENSE):  # R2 forward and spliced on a + motif, R1 reverse: antisense
        reads += [
            read(f"a_{k}", 0, 2080, 1550, True, True, [(_M, 100)]),
            read(f"a_{k}", 0, 1550, 2080, False, False, spliced, xs="+"),
        ]
    for k in range(
        N_MULTI
    ):  # hit 1: R2 forward and spliced on a + motif, R1 reverse; hit 2: unspliced on chr2
        reads += [
            read(f"mm_{k}", 0, 2080, 1550, True, True, [(_M, 100)], nh=2, hi=1),
            read(f"mm_{k}", 0, 1550, 2080, False, False, spliced, nh=2, hi=1, xs="+"),
            read(f"mm_{k}", 1, 1320, 1120, True, True, [(_M, 100)], nh=2, hi=2, secondary=True),
            read(f"mm_{k}", 1, 1120, 1320, False, False, [(_M, 100)], nh=2, hi=2, secondary=True),
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
    return payload, RegionArrays.from_index(index), build_sj_arrays(index), strand_model


def test_a_wrong_strand_spliced_placement_carries_one_minus_p_beside_an_unspliced_one(scanned):
    payload, regions, sj, strand_model = scanned
    p = strand_model.p_r1_sense
    assert 0.9 < p < 1.0, (
        "the fixture's library must be stranded, and not perfectly, for 1 - p to be a probability"
    )
    scores = score_held_fragments(
        payload,
        fl_models=build_fl_models(payload),
        rna_sense_frac=p,
        regions=regions,
        sj=sj,
    )
    d = payload.deferred
    multi = np.flatnonzero(np.diff(d.placement_offsets) == 2)
    assert multi.size == N_MULTI, "every multimapper must be held on both of its placements"
    expected_spliced, expected_genomic = None, None
    for i in multi:
        placements = d.placements_of(int(i))
        spliced_p = [q for q in placements if int(d.sj_strand[q]) != int(Strand.NONE)]
        assert len(spliced_p) == 1, "one hit has a sequenced motif, the other none"
        q = spliced_p[0]
        align = int(d.align_strand[q])
        expected_spliced, expected_genomic = strand_terms(
            align=align, implied_strand=int(d.sj_strand[q]), rna_sense_frac=p
        )
        hyp = slice(int(d.hypothesis_offsets[q]), int(d.hypothesis_offsets[q + 1]))
        np.testing.assert_allclose(scores.terms.strand[hyp], expected_spliced, rtol=1e-12, atol=0.0)
    assert expected_spliced == pytest.approx(1.0 - p), (
        "R1 on the strand opposite the motif: the wrong strand"
    )
    assert expected_spliced < expected_genomic, (
        "on the wrong strand the spliced reading must lose to gDNA's ½"
    )
