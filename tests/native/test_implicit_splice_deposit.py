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


def _scan(scenario, transcripts):
    scenario.add_gene("g1", "+", transcripts)
    result = scenario.build_oracle(n_fragments=1500, sim_config=SIM)
    _, _, _, payload = scan_and_buffer(
        str(result.bam_path), result.index, BamScanConfig(sj_strand_tag="auto")
    )
    return payload


@pytest.fixture
def scenario(tmp_path):
    sc = Scenario("implicit", genome_length=8000, seed=SEED, work_dir=tmp_path / "implicit")
    yield sc
    sc.cleanup()


def test_ONE_candidate_transcript_DEPOSITS_with_the_strand_inferred_from_it(scenario):
    """The unanimous case. One isoform, so the implied intron is not in doubt and the fragment tallies.

    ⛔ If this ever goes to zero deposits, the unanimity test has become unsatisfiable and the whole
    implicit population is silently deferred — which is exactly the failure the real-data census looked
    like at first.
    """
    payload = _scan(
        scenario, [{"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100}]
    )
    # ⚠ Non-vacuity, and it is read off the UMBRELLA CENSUS rather than off a splice label. The census
    # counts every fragment whose gap needed resolving, which is exactly this scenario's population — while
    # `splice_type` is the scanner's record of what it SAW and is a different axis. ⛔ If this ever goes to
    # zero the scenario stopped producing gap introns and everything below passes vacuously.
    assert payload.gap_resolution.gap_resolved_spliced > 0, (
        "no fragment had an intron in its unsequenced gap at all — the scenario stopped producing them, "
        "so this test would pass vacuously"
    )
    assert payload.qc.deferred_undetermined_gap == 0, (
        "a single candidate transcript cannot disagree with itself, so nothing here may be deferred"
    )
    assert payload.deferred.n_fragments == 0
    # It deposited: its junction was credited, and it is barred from the pure-RNA length pool.
    assert int(payload.sj_count.sum()) > 0


def test_TWO_candidates_implying_DIFFERENT_introns_are_DEFERRED(scenario):
    """The ambiguous case: two isoforms whose introns start at different places inside the same gap.

    Their implied ``L`` differs by 100 bp, so there is no answer to deposit — picking either is picking a
    fragment length at random.
    """
    payload = _scan(
        scenario,
        [
            {"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 1100), (3000, 3200)], "abundance": 100},
        ],
    )
    assert payload.qc.deferred_undetermined_gap > 0, (
        "two isoforms imply introns [1200,3000) and [1100,3000) in the same gap — a 100 bp difference in L "
        "— so those fragments have no determined path and must not deposit"
    )
    # ⭐ And they are HELD, not dropped: the bank carries as many fragments as the counter claims, with
    # both paths on each, which is what makes them recoverable in the second pass.
    assert payload.deferred.n_fragments == payload.qc.deferred_undetermined_gap
    assert int(np.diff(payload.deferred.hypothesis_offsets).min()) >= 2


def test_A_RETAINED_INTRON_ISOFORM_ALONE_IS_ENOUGH_TO_DEFER(scenario):
    """⚠ The case that makes the rule strict on real data, and the reason the census reads 100 %.

    ``t2`` covers the whole locus as one exon, so for it the mate gap holds no intron at all: the fragment
    is unspliced with an ``L`` that includes the gap. That is a different hypothesis from ``t1``'s spliced
    ``L``, so "implies nothing here" has to count as a distinct answer — and GENCODE is full of
    retained-intron isoforms.

    ⚠ **The locus is TIGHT on purpose** — a 300 bp intron, so the unspliced hypothesis's ``L`` is about
    560 bp and stays under ``max_fragment_length``. Spread the exons 1800 bp apart instead and the
    unspliced hypothesis is ruled out **by length** and the fragment deposits; that is the next test, and
    keeping the two apart is what stops either passing for the other's reason.
    """
    payload = _scan(
        scenario,
        [
            {"t_id": "t1", "exons": [(1000, 1200), (1500, 1700)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 1700)], "abundance": 100},
        ],
    )
    assert payload.qc.deferred_undetermined_gap > 0, (
        "one candidate implies an intron in the gap and the other implies none, which are two different "
        "fragment lengths for one molecule"
    )
    # ⛔ THE SUBCLASS NAMES THE QUESTION, and here it is the composition one: the unspliced path against
    # one spliced path is "RNA or gDNA", one bit. `t2` covering the locus as a single exon is what puts the
    # unspliced path in the set, and GENCODE is full of such retained-intron isoforms.
    assert payload.gap_resolution.gap_deferred_rna_or_gdna > 0
    assert payload.deferred.n_fragments == payload.qc.deferred_undetermined_gap
    # ⭐ Every held record carries BOTH paths, one of them empty — the unspliced hypothesis needs no flag,
    # because cutting nothing IS the statement that the gap is real template.
    empty = [
        h
        for h in range(payload.deferred.n_hypotheses)
        if payload.deferred.hypothesis_introns_of(h).size == 0
    ]
    assert len(empty) == payload.deferred.n_fragments, (
        "exactly one hypothesis per held fragment must be the empty (unspliced) path"
    )


def test_A_SPAN_OVER_THE_LIMIT_RULES_OUT_the_retained_intron_hypothesis(scenario):
    """⭐ The other side of the same rule, on a real scan — and it is **not a separate rule**.

    The owner's *"if the genomic span exceeds ``max_fragment_length``, assume it is RNA"* is the ordinary
    hypothesis filter applied to the unspliced path, whose ``L`` **is** that span. Spread ``t1``'s exons
    1800 bp apart and a junction-spanning fragment's unspliced ``L`` is ~2100 bp against a limit of 1000, so
    the retained-intron explanation is deleted and the spliced path stands alone and deposits.

    ⛔ **This is the concern ``SPEC_GAP_PATHS.md`` §8 C3 names, and it is why this test exists.** The filter
    is *not* purely a cost gate: it changes CLASSIFICATION, so the same annotation defers here and resolves
    there depending only on how far apart the exons sit. Measuring that from both sides is the difference
    between a documented consequence and an assumption.
    """
    payload = _scan(
        scenario,
        [
            {"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 3200)], "abundance": 100},
        ],
    )
    assert payload.gap_resolution.gap_resolved_spliced > 0, (
        "no fragment had a gap intron at all, so the span rule is not being exercised — this would pass "
        "vacuously on a scenario that stopped producing mate gaps"
    )
    assert payload.qc.deferred_undetermined_gap == 0, (
        "the unspliced hypothesis's L IS the genomic span, ~2100 bp here against a 1000 bp limit, so it "
        "must be filtered and the fragment's path determined"
    )
    assert payload.deferred.n_fragments == 0
    assert payload.qc.dropped_too_long == 0, (
        "and the SPLICED path is well under the limit, so nothing may be dropped as too long — a fragment "
        "rejected here would mean the filter deleted the wrong hypothesis"
    )
