"""⭐ P4 — the second pass is WIRED IN, ahead of calibration, and the gDNA pools do not move.

     (where the drain sits) and its P4
    Measurement: — the tail and the anchor, scored on the pilot against truth

§2's structural claim is that the drain runs **between the scan and calibration**, which is what lets
calibration run exactly once on the complete tally. Two things follow that are worth gating rather than
assuming: that calibration really does see the drained tally, and that the drain's own seed is its own.

⭐ **THE gDNA CONTROL IS EXACT, AND IT IS STRUCTURAL.** Not "moves by less than something" — the four
gDNA length pools move by **exactly zero**, and there is a derivation:

* a held fragment has ≥ 2 hypotheses, so its gap contains ≥ 1 annotated intron, whose endpoints are cuts;
* if the drain picks ∅ the molecule crosses **both** those lines, making it a multi-line crossing — and
   gives a multi-line crossing **no pool**, because it is a gDNA/RNA mixture;
* if the drain picks a spliced path the fragment used an annotated junction, so it is ``RNA_SPLICED``.

So the drain can only ever touch ``RNA_SPLICED``, and measured on all 8 pilot conditions it does: 100 % of
the ``pool_lengths`` delta lands there and the four gDNA rows are byte-identical.

⚠ **This is NOT TRAPS: pure-and-length-censored.6's gDNA control, and the difference matters.** There the fix could not reach a fragment
with no introns, so *any* movement was a bug. Here the drain reaches fragments it is supposed to reach —
on a gdna100 library a real gDNA fragment whose mate gap spans an annotated intron is genuinely ambiguous
and was genuinely held. What this gate says is that a drained fragment never lands in a pool that is
supposed to be **pure**, whichever hypothesis won.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import _drain_side_buffer, scan_and_buffer
from rigel.scan_payload import (
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTERGENIC_EXON,
    POOL_DNA_INTRON_EXON,
    POOL_DNA_INTRONIC,
    POOL_RNA_SPLICED,
)

SEED = 4242

#: The four pools that are pure gDNA by construction. ⭐ Named from the schema rather than by index, so a
#: reordering of the pool axis cannot silently point this gate at the RNA row.
GDNA_POOLS = (
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTRONIC,
    POOL_DNA_INTRON_EXON,
    POOL_DNA_INTERGENIC_EXON,
)


@pytest.fixture(scope="module")
def scenario():
    """A scenario with gDNA and multi-isoform genes, so fragments are genuinely HELD."""
    from rigel.sim import ReadSimConfig, Scenario

    sc = Scenario("p4", genome_length=20_000, seed=SEED)
    # ⚠ Two isoforms differing in ONE internal intron is what makes a mate gap ambiguous, and it is the
    # only geometry that fills the side buffer at all.
    sc.add_gene(
        "gA",
        "+",
        [
            {"t_id": "tA1", "exons": [(1000, 1600), (2000, 2600), (3000, 3400)], "abundance": 80},
            {"t_id": "tA2", "exons": [(1000, 1700), (1900, 2600), (3000, 3400)], "abundance": 60},
        ],
    )
    sc.add_gene(
        "gB",
        "-",
        [
            {"t_id": "tB1", "exons": [(8000, 8600), (9000, 9600)], "abundance": 50},
            {"t_id": "tB2", "exons": [(8000, 8700), (8900, 9600)], "abundance": 40},
        ],
    )
    sim = ReadSimConfig(
        frag_mean=260,
        frag_std=70,
        frag_min=90,
        frag_max=700,
        read_length=100,
        strand_specificity=0.99,
        seed=SEED,
    )
    result = sc.build_oracle(n_fragments=4000, sim_config=sim, gdna_fraction=0.25)
    yield result
    sc.cleanup()


def _config(**overrides) -> PipelineConfig:
    base = dict(
        em=EMConfig(seed=SEED, assignment_mode="map", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    base.update(overrides)
    return PipelineConfig(**base)


@pytest.fixture(scope="module")
def scanned(scenario):
    config = _config()
    stats, strand_models, _buffer, payload = scan_and_buffer(
        str(scenario.bam_path), scenario.index, config.scan
    )
    del stats
    return scenario.index, strand_models, payload, config


def test_the_fixture_actually_HOLDS_fragments(scanned):
    """⚠ Non-vacuity, first. With an empty side buffer the drain is a documented no-op and every gate below
    would compare a payload with itself."""
    _index, _strand, payload, _cfg = scanned
    assert payload.deferred.n_fragments > 0, (
        "the scenario must produce ambiguous mate gaps or there is nothing for P4 to gate; two isoforms "
        "differing in one internal intron is the geometry that does it"
    )
    assert payload.drain is None, "a freshly scanned payload has not been drained"


def test_the_gDNA_LENGTH_POOLS_DO_NOT_MOVE(scanned):
    """⭐ **THE P4 CONTROL, and it is EXACT.** The four pure-gDNA pools are byte-identical after the drain.

    ⛔ Not a tolerance. A held fragment's gap holds an annotated intron whose endpoints are cuts, so a
    chosen ∅ crosses two lines — a multi-line crossing, which deliberately gives
    **no pool** because it is a gDNA/RNA mixture — and a chosen spliced path used an annotated junction, so
    it is ``RNA_SPLICED``. Either way a drained fragment cannot enter a pool that is supposed to be pure.

    ⚠ **An impure pool is worse than a missing one** (§8), and these four pools are what the gDNA
    fragment-length model is fitted from. A drained RNA fragment landing in ``DNA_INTRONIC`` would fit the
    gDNA length model partly from RNA, and nothing downstream would look wrong.
    """
    index, strand_models, payload, cfg = scanned
    drained = _drain_side_buffer(payload, index, strand_models, seed=cfg.second_pass_seed)
    before = np.asarray(payload.pool_lengths, dtype=np.int64)
    after = np.asarray(drained.pool_lengths, dtype=np.int64)

    for pool in GDNA_POOLS:
        assert np.array_equal(after[pool], before[pool]), (
            f"pool {pool} is pure gDNA by construction and the drain moved it by "
            f"{int(after[pool].sum() - before[pool].sum())} fragments. A drained fragment can only be a "
            f"multi-line crossing (no pool) or an annotated splice (RNA_SPLICED)."
        )
    # ⭐ And the delta is not zero everywhere — otherwise the gate above passes because nothing happened.
    assert after[POOL_RNA_SPLICED].sum() > before[POOL_RNA_SPLICED].sum(), (
        "the drain must have deposited spliced fragments into RNA_SPLICED, or this fixture drained nothing "
        "and the control above is vacuous"
    )


def test_CALIBRATION_sees_the_DRAINED_tally(scenario):
    """⭐ §2's structural claim: the drain runs BEFORE calibration, so calibration runs once on the whole
    tally. Gated through ``run_pipeline`` rather than by reading the wiring.

    ⚠ The observable is the fragment-length model, because that is what the drained mass actually changes —
    the held fragments are the LONG ones, so the anchor an undrained calibration would see is biased short.
    Measured on the pilot: the anchor's mean error against truth goes from −1.6 % to **+0.00 %**.
    """
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.pipeline import run_pipeline

    config = _config()
    result = run_pipeline(scenario.bam_path, scenario.index, config=config)

    _stats, strand_models, _buffer, pass_one = scan_and_buffer(
        str(scenario.bam_path), scenario.index, config.scan
    )
    drained = _drain_side_buffer(
        pass_one, scenario.index, strand_models, seed=config.second_pass_seed
    )
    assert drained.drain is not None and drained.drain.offered > 0

    # The eff lengths the pipeline shipped must match the DRAINED payload's RNA pmf, not pass one's.
    from rigel.frag_length_model import FragmentLengthModel

    exonic = scenario.index.t_df["length"].values.astype(np.int64)
    shipped = np.asarray(result.estimator.effective_lengths, dtype=np.float64)
    for label, payload in (("drained", drained), ("pass one", pass_one)):
        # ⭐ TRAPS: divide-by-a-probability's de-tilt is part of production's route to the RNA pmf, so it is part of this
        # reproduction; it is the same array for both arms, so it cannot decide which one matches.
        fl = build_fl_models(
            payload,
            junction_opportunity=crossing_probability_from_index(
                scenario.index, int(payload.max_length)
            ),
        )
        expected = FragmentLengthModel.from_pmf(
            fl.rna_pmf, fl.max_size
        ).compute_all_transcript_eff_lens(exonic)
        matches = np.allclose(shipped, expected, rtol=0, atol=0)
        if label == "drained":
            assert matches, (
                "calibration's effective lengths do not come from the DRAINED tally; the drain is either "
                "not wired in or runs after the fragment-length models are built"
            )
        else:
            assert not matches, (
                "pass one's tally and the drained tally give the SAME effective lengths, so this gate "
                "cannot tell which one calibration used — the fixture drained too little to be decisive"
            )


def test_the_DRAIN_SEED_is_its_own_and_is_NOT_the_EM_seed(scanned):
    """⛔ Two independent RNG consumers must not share one field. If the drain read ``em.seed``, an EM A/B
    would silently re-draw every held fragment — so the arms would differ in the tally as well as in the
    estimator, and the comparison would be meaningless.

    ⭐ Also §5.2's reproducibility, at the config level: one seed, one answer.
    """
    index, strand_models, payload, _cfg = scanned

    def drained(second_pass_seed: int):
        return _drain_side_buffer(payload, index, strand_models, seed=second_pass_seed)

    first, again, other = drained(1), drained(1), drained(2)
    assert np.array_equal(first.deposited_lengths, again.deposited_lengths), (
        "one drain seed must give one answer"
    )
    assert not np.array_equal(first.deposited_lengths, other.deposited_lengths), (
        "a different drain seed must actually change the draw, or the seed is not reaching it"
    )
    # ⭐ The EM's seed is a different field and cannot reach the drain: PipelineConfig carries both, and
    # only `second_pass_seed` is passed to `_drain_side_buffer` in `run_pipeline`.
    assert _config(em=EMConfig(seed=999)).second_pass_seed == PipelineConfig().second_pass_seed, (
        "changing the EM seed must not change the drain's seed"
    )


def test_an_EMPTY_side_buffer_is_a_documented_no_op(scanned):
    """⚠ A library with no annotated intron in any mate gap is a real state, not an error. The payload must
    come back untouched and still report ``drain is None``, so "pass one only" keeps one spelling."""
    import dataclasses

    from rigel.scan_payload import DeferredFragments, GapCensus

    index, strand_models, payload, _cfg = scanned
    empty = dataclasses.replace(
        payload,
        deferred=DeferredFragments.empty(),
        gap_resolution=GapCensus.zeros(),
        qc=dataclasses.replace(payload.qc, deferred_undetermined_gap=0),
    )
    out = _drain_side_buffer(empty, index, strand_models, seed=0)
    assert out is empty, "with nothing held the drain must return the payload it was given"
    assert out.drain is None
