"""An implicit splice deposits when its path is determined, and defers when it is not.

    Rule: ``docs/ACCUMULATOR_DESIGN.md`` §9.1 (owner ruling, 2026-07-29)

A `SPLICE_IMPLICIT` fragment has an annotated intron inside its unsequenced mate gap. The splice motif was
never read, so `sj_strand` comes from the transcript that implied it — which is only legitimate when the
candidates agree about *which* intron is there. If they do not, `L` is undetermined and nothing about the
fragment can be tallied.

⚠ THIS MODULE EXISTS BECAUSE THE FIRST MEASUREMENT LOOKED LIKE A BUG. On three real cfRNA libraries the
rule deferred **100 %** of implicit fragments and deposited none, which is either a true property of real
candidate sets or a mistake in the unanimity test. Both arms are pinned here on scenarios small enough to
reason about, so the real-data number can be read as evidence rather than guessed at.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel import scan_payload
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario


SEED = 4242

#: Exons far enough apart, and reads short enough, that a fragment spanning the junction puts BOTH mates
#: wholly inside an exon and the intron wholly inside the gap — which is what SPLICE_IMPLICIT means.
SIM = ReadSimConfig(
    frag_mean=260,
    frag_std=20,
    frag_min=220,
    frag_max=320,
    read_length=100,
    strand_specificity=1.0,
    seed=SEED,
)


def _scan_qc(scenario, transcripts, monkeypatch) -> dict:
    scenario.add_gene("g1", "+", transcripts)
    result = scenario.build_oracle(n_fragments=1500, sim_config=SIM)
    monkeypatch.setattr(
        scan_payload.AccumulatorPayload,
        "from_scan_result",
        classmethod(lambda cls, scan_result: scan_result["calibration"]),
    )
    _, _, _, _, payload = scan_and_buffer(
        str(result.bam_path), result.index, BamScanConfig(sj_strand_tag="auto")
    )
    return dict(payload["qc"]), payload


@pytest.fixture
def scenario(tmp_path):
    sc = Scenario("implicit", genome_length=8000, seed=SEED, work_dir=tmp_path / "implicit")
    yield sc
    sc.cleanup()


def test_ONE_candidate_transcript_DEPOSITS_with_the_strand_inferred_from_it(scenario, monkeypatch):
    """The unanimous case. One isoform, so the implied intron is not in doubt and the fragment tallies.

    ⛔ If this ever goes to zero deposits, the unanimity test has become unsatisfiable and the whole
    implicit population is silently deferred — which is exactly the failure the real-data census looked
    like at first.
    """
    qc, payload = _scan_qc(
        scenario,
        [{"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100}],
        monkeypatch,
    )
    assert qc["sj_implicit_fragments"] > 0, (
        "no fragment was classified as an implicit splice at all — the scenario stopped producing them, "
        "so this test would pass vacuously"
    )
    assert qc["dropped_ambiguous_path"] == 0, (
        "a single candidate transcript cannot disagree with itself, so nothing here may be deferred"
    )
    # It deposited: its junction was credited, and it is barred from the pure-RNA length pool.
    assert int(np.asarray(payload["sj_count"]).sum()) > 0


def test_TWO_candidates_implying_DIFFERENT_introns_are_DEFERRED(scenario, monkeypatch):
    """The ambiguous case: two isoforms whose introns start at different places inside the same gap.

    Their implied ``L`` differs by 100 bp, so there is no answer to deposit — picking either is picking a
    fragment length at random.
    """
    qc, _ = _scan_qc(
        scenario,
        [
            {"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 1100), (3000, 3200)], "abundance": 100},
        ],
        monkeypatch,
    )
    assert qc["dropped_ambiguous_path"] > 0, (
        "two isoforms imply introns [1200,3000) and [1100,3000) in the same gap — a 100 bp difference in L "
        "— so those fragments have no determined path and must not deposit"
    )


def test_A_RETAINED_INTRON_ISOFORM_ALONE_IS_ENOUGH_TO_DEFER(scenario, monkeypatch):
    """⚠ The case that makes the rule strict on real data, and the reason the census reads 100 %.

    ``t2`` covers the whole locus as one exon, so for it the mate gap holds no intron at all: the fragment
    is unspliced with an ``L`` that includes the gap. That is a different hypothesis from ``t1``'s spliced
    ``L``, so "implies nothing here" has to count as a distinct answer — and GENCODE is full of
    retained-intron isoforms.
    """
    qc, _ = _scan_qc(
        scenario,
        [
            {"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 3200)], "abundance": 100},
        ],
        monkeypatch,
    )
    assert qc["dropped_ambiguous_path"] > 0, (
        "one candidate implies an intron in the gap and the other implies none, which are two different "
        "fragment lengths for one molecule"
    )
