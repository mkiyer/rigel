"""The second pass's strand factor for a placement whose junction motif was sequenced.

A held record's placements are compared by the product of their factors, so every placement's strand factor must
be the probability of its observed orientation under its own reading. A placement with a sequenced motif is RNA on
that motif's strand under every hypothesis, so its factor is ``p`` or ``1 - p``, not 1: a constant is harmless within
one placement, but across a multimapper's placements a 1 outweighs an unspliced hit's ½ even on the wrong strand.
Specification: ``second_pass.strand_terms``.
"""

from __future__ import annotations

import numpy as np
import pytest
from _bam_builder import R1_FWD, R1_REV, R2_FWD, R2_REV, SECONDARY, record, write_bam
from _index_builder import build_index
from _length_models import payload_fl_models

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.splice_graph import build_sj_arrays
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.second_pass import score_held_fragments, strand_terms
from rigel.types import Strand

GENOME = 6_000
N_SENSE = (
    50  # uniquely mapped spliced fragments that set the library's strand: R1 on the motif's strand
)
N_ANTISENSE = 2  # and a few on the other strand, so the library is stranded but not perfectly
N_MULTI = 20
N_BACKGROUND = 40  # unique unspliced pairs 400 bases long in chr2's intergenic tail: the all-fragment law is not RNA's
_M, _N = 0, 3


def _gtf() -> str:
    rows = []
    for ref in ("chr1", "chr2"):
        rows += [
            f'{ref}\tt\texon\t1001\t1600\t.\t+\t.\tgene_id "gW_{ref}"; transcript_id "tW_{ref}";\n',
            f'{ref}\tt\texon\t2001\t2200\t.\t+\t.\tgene_id "gW_{ref}"; transcript_id "tW_{ref}";\n',
        ]
    return "".join(rows)


def _scan(tmp_path_factory, label: str, multi_xs):
    """``N_SENSE`` unique fragments spliced across chr1's annotated junction with R1 on the motif's strand, and
    ``N_MULTI`` multimappers whose first hit is spliced across that junction with R1 on the other strand, its motif
    tag *multi_xs* (``None``: the aligner wrote none), and whose second hit is unspliced, contained in chr2's first
    exon."""
    base = tmp_path_factory.mktemp(f"strand_placements_{label}")
    index = build_index(base, _gtf(), {"chr1": GENOME, "chr2": GENOME})

    spliced = [(_M, 50), (_N, 400), (_M, 50)]  # [1550, 1600) | intron [1600, 2000) | [2000, 2050)
    reads = []
    for k in range(N_SENSE):  # R1 forward and spliced on a + motif: sense
        reads += [
            record(f"s_{k}", R1_FWD, 1550, spliced, 2080, xs="+"),
            record(f"s_{k}", R2_REV, 2080, [(_M, 100)], 1550),
        ]
    for k in range(N_ANTISENSE):  # R2 forward and spliced on a + motif, R1 reverse: antisense
        reads += [
            record(f"a_{k}", R1_REV, 2080, [(_M, 100)], 1550),
            record(f"a_{k}", R2_FWD, 1550, spliced, 2080, xs="+"),
        ]
    for k in range(N_BACKGROUND):
        pos = 3000 + 10 * k
        reads += [
            record(f"b_{k}", R1_FWD, pos, [(_M, 100)], pos + 300, ref=1),
            record(f"b_{k}", R2_REV, pos + 300, [(_M, 100)], pos, ref=1),
        ]
    # hit 1: R2 forward and spliced on a + motif, R1 reverse; hit 2: unspliced on chr2
    for k in range(N_MULTI):
        reads += [
            record(f"mm_{k}", R1_REV, 2080, [(_M, 100)], 1550, nh=2, hi=1),
            record(f"mm_{k}", R2_FWD, 1550, spliced, 2080, nh=2, hi=1, xs=multi_xs),
            record(f"mm_{k}", R1_REV | SECONDARY, 1320, [(_M, 100)], 1120, ref=1, nh=2, hi=2),
            record(f"mm_{k}", R2_FWD | SECONDARY, 1120, [(_M, 100)], 1320, ref=1, nh=2, hi=2),
        ]
    bam = write_bam(base / "p.bam", reads, {"chr1": GENOME, "chr2": GENOME})
    _s, strand_model, _b, payload = scan_and_buffer(
        bam, index, BamScanConfig(sj_strand_tag="XS", total_threads=1)
    )
    return payload, RegionArrays.from_index(index), build_sj_arrays(index), strand_model


@pytest.fixture(scope="module")
def scanned(tmp_path_factory):
    return _scan(tmp_path_factory, "tagged", "+")


@pytest.fixture(scope="module")
def scanned_untagged(tmp_path_factory):
    return _scan(tmp_path_factory, "untagged", None)


def _spliced_placements(scanned):
    """Score the fixture; yield each multimapper's spliced placement with its hypotheses' slice."""
    payload, regions, sj, strand_model = scanned
    fl = payload_fl_models(payload)
    scores = score_held_fragments(
        payload, fl_models=fl, rna_sense_frac=strand_model.p_r1_sense, regions=regions, sj=sj
    )
    d = payload.deferred
    multi = np.flatnonzero(np.diff(d.placement_offsets) == 2)
    assert multi.size == N_MULTI
    out = []
    for i in multi:
        q = [p for p in d.placements_of(int(i)) if d.observed_introns_of(p).shape[0]]
        assert len(q) == 1, "one hit is spliced across the annotated junction, the other contiguous"
        q = q[0]
        out.append((q, slice(int(d.hypothesis_offsets[q]), int(d.hypothesis_offsets[q + 1]))))
    return scores, fl, d, strand_model.p_r1_sense, out


def test_an_untagged_annotated_splice_reads_its_junctions_strand(scanned_untagged):
    """The deposit books a placement spliced when an observed intron resolves to an annotated junction, whatever the
    aligner tagged; so the draw reads it as RNA on the strand of the transcripts holding the junction, ``1 - p`` for
    R1 on the other strand, never gDNA's ½.

    Perturbation: reading the motif tag alone gives ½ here.
    """
    scores, _fl, d, p, spliced = _spliced_placements(scanned_untagged)
    for q, hyp in spliced:
        assert int(d.sj_strand[q]) == int(Strand.NONE), "the aligner wrote no motif on this hit"
        np.testing.assert_allclose(scores.terms.strand[hyp], 1.0 - p, rtol=1e-12, atol=0.0)


def test_an_observed_annotated_splice_reads_the_rna_length_law(scanned):
    """A placement whose observed intron resolves to an annotated junction is certified RNA under every one of its
    paths, so each takes the RNA length law; the unconditional anchor is for a path whose component is unknown.

    Perturbation: choosing the law from the implied introns alone reads the anchor here.
    """
    scores, fl, d, _p, spliced = _spliced_placements(scanned)
    assert not np.array_equal(fl.rna_pmf, fl.global_pmf), (
        "the two laws must differ for this gate to see"
    )
    for _q, hyp in spliced:
        L = scores.terms.length[hyp]
        np.testing.assert_allclose(
            scores.terms.length_likelihood[hyp], fl.rna_pmf[L], rtol=1e-12, atol=0.0
        )


def test_a_spliced_placement_never_takes_gdnas_bound(scanned):
    """gDNA's floor holds a placement at least at the library's bulk only where gDNA can be. A placement the deposit
    books spliced holds no gDNA, so it keeps its own junction traffic even where that is none: chr2's junction,
    which no unique fragment crossed, reads 0 beside a bulk of 1, and chr1's reads its own traffic.

    Perturbation: flooring every placement reads chr2's junction at the bulk.
    """
    from rigel.calibration.effective_length import observable_contained_eff_length
    from rigel.second_pass import PlacementPath, placement_abundance

    payload, regions, sj, _strand_model = scanned
    d = payload.deferred
    support = observable_contained_eff_length(regions, payload_fl_models(payload).global_pmf)
    chr1_sj = int(payload.ref_sj_offsets[0])
    chr2_sj = int(payload.ref_sj_offsets[1])
    assert (
        float(payload.sj_inv_length_sum[chr2_sj]) == 0.0 < float(payload.sj_inv_length_sum[chr1_sj])
    )
    none = np.zeros(0, np.int64)
    for junction, expected in (
        (chr2_sj, 0.0),
        (chr1_sj, float(payload.sj_inv_length_sum[chr1_sj])),
    ):
        path = PlacementPath(sj=np.array([junction], np.int64), boundaries=none, regions=none)
        assert placement_abundance(payload, path, support, bulk=1.0) == expected
    unspliced = PlacementPath(sj=none, boundaries=none, regions=none)
    assert placement_abundance(payload, unspliced, support, bulk=1e9) == 1e9, (
        "a placement no annotated junction cuts reads at least the bulk"
    )
    assert d.n_placements > 0


def test_a_wrong_strand_spliced_placement_carries_1_minus_p_beside_an_unspliced_one(scanned):
    payload, regions, sj, strand_model = scanned
    p = strand_model.p_r1_sense
    assert 0.9 < p < 1.0, (
        "the fixture's library must be stranded, and not perfectly, for 1 - p to be a probability"
    )
    scores = score_held_fragments(
        payload,
        fl_models=payload_fl_models(payload),
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
