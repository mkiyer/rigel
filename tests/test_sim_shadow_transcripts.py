"""SHADOW TRANSCRIPTS — unannotated transcription the simulator draws from and the index never sees
(owner design, 2026-08-29).

The gate on the MERGE: shadows are appended to the index's transcript list for SIMULATION with
``t_index`` continuing after the index's rows, no nascent, and are REFUSED if the index already knows
their id — a shadow the tool can see is not a shadow. Their fragments are named like any RNA fragment, so
the oracle files them as ``mrna`` without a read-name change (``rigel.sim.read_name.parse_origin``).
"""

from __future__ import annotations

import pytest

from rigel.sim.read_name import parse_origin
from rigel.sim.whole_genome import merge_shadow_transcripts
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _tx(t_id, ref, exons, strand=Strand.POS, t_index=-1):
    t = Transcript(
        ref=ref,
        strand=strand,
        exons=[Interval(a, b) for a, b in exons],
        t_id=t_id,
        g_id="g_" + t_id,
    )
    t.length = t.compute_length()
    t.t_index = t_index
    return t


def _gtf(tmp_path, rows):
    p = tmp_path / "shadow.gtf"
    p.write_text(
        "".join(
            f'{ref}\trigel_test\texon\t{a}\t{b}\t.\t{st}\t.\tgene_id "g_{tid}"; transcript_id "{tid}";\n'
            for tid, ref, st, (a, b) in rows
        )
    )
    return p


def test_shadows_are_appended_after_the_index_rows_with_no_nascent(tmp_path):
    index_tx = [
        _tx("A", "chr1", [(100, 200)], t_index=0),
        _tx("B", "chr1", [(500, 900)], t_index=1),
    ]
    gtf = _gtf(tmp_path, [("S1", "blank", "+", (1000, 2000)), ("S2", "blank", "-", (3000, 3500))])
    merged = merge_shadow_transcripts(list(index_tx), gtf)
    assert [t.t_id for t in merged] == ["A", "B", "S1", "S2"]
    s1, s2 = merged[2], merged[3]
    assert (s1.t_index, s2.t_index) == (2, 3)  # continue after the index's rows
    assert (
        not s1.is_nrna
        and not s1.is_synthetic
        and s1.nrna_t_index == -1
        and s1.nrna_abundance == 0.0
    )
    assert s2.strand == Strand.NEG


def test_a_shadow_the_index_already_knows_is_refused(tmp_path):
    """PERTURBATION of the premise: name a shadow after an annotated transcript and the merge must refuse."""
    index_tx = [_tx("A", "chr1", [(100, 200)], t_index=0)]
    gtf = _gtf(tmp_path, [("A", "blank", "+", (1000, 2000))])
    with pytest.raises(ValueError, match="KNOWN to the index"):
        merge_shadow_transcripts(index_tx, gtf)


def test_a_shadow_fragment_is_filed_as_mrna_by_the_oracle():
    """No read-name change is needed: ``{t_id}:{start}-{end}:{strand}:{index}`` parses as RNA of that
    transcript, so the certified truth shows RNA where the annotation says there is none."""
    o = parse_origin("shadow_U1_SE2kb_plus:100050-100260:+:7")
    assert o.kind == "mrna"
    assert o.transcript_id == "shadow_U1_SE2kb_plus"
