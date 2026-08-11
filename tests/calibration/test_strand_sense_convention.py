"""⭐ WHAT `rna_sense_frac` MEANS — and why 0.0101 on a "0.99 stranded" library is CORRECT.

    Retires: rank 6 · Corrects: + S5.f-addendum,
    fact 17 and §0 TRAPS: opposite-tilts-must-not-pool · Unblocks: (the strand term)

⛔ **THIS WAS FILED TWICE AS A SIGN BUG. IT IS NOT ONE.** The record said the fitted κ was `1 − truth`
and that "only the exported scalar is mis-labelled". Both statements are wrong, and the reason is a
collision between **two different quantities that were both being called strand specificity**:

===========================================  =====================================================
``ReadSimConfig.strand_specificity``          *"probability an RNA fragment preserves correct read
                                              orientation … an R1↔R2 swap with probability 1 − ss"*
                                              — protocol FIDELITY, and direction-agnostic.
``StrandModel.p_r1_sense``                    ``P(align_strand == the junction's strand)`` — DIRECTIONAL,
                                              and its own docstring already says it: *"High (≈0.95) for
                                              R1-sense libraries (KAPA Stranded). Low (≈0.05) for
                                              R1-antisense libraries (Illumina TruSeq dUTP)."*
===========================================  =====================================================

⭐ **For an R1-antisense protocol these are complements**, so comparing one against the other reads as a
sign error and is not one. The simulator emits R1-antisense (the gDNA path is explicit about it:
``r2_seqs, r1_seqs = _batch_extract_reads(...)`` — the first extracted read becomes **R2**), which is the
most common real protocol. A dUTP library at 99 % fidelity genuinely has a sense fraction near 0.01.

⭐ **And the tool already carries the direction-agnostic quantity, under the right name.**
``StrandModel.strand_specificity`` is ``max(p_r1_sense, p_r1_antisense)`` and it recovers the simulated
parameter directly. That is what this module pins, because its absence is what let the same non-defect be
filed twice.

⚠ **The 166× measurement stands and now has an explanation.** -addendum forced κ to the
nominal 0.99 and a zero-gDNA library read ``f_gdna = 0.4992`` against the fitted value's ``0.0030``. That
is not "the mirror cancels" — it is ``0.0101`` being the **right answer** and ``0.99`` being a different
quantity substituted for it.
"""

from __future__ import annotations

import collections
import dataclasses
import pathlib
import tempfile

import pysam
import pytest

from rigel.calibration.strand_balance import fit_strand_balance
from rigel.config import BamScanConfig, PipelineConfig
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.sim.read_name import parse_origin


SEED = 11

#: ⚠ Two genes, one per strand. With a single-strand locus ``p_r1_sense`` is still well defined but a
#: convention error that swapped the comparison's operands would be invisible — both operands would move
#: together. Opposite-strand genes are what make the comparison discriminating.
GENES = (
    ("g1", "+", [{"t_id": "t1", "exons": [(1000, 1400), (2000, 2400)], "abundance": 100}]),
    ("g2", "-", [{"t_id": "t2", "exons": [(5000, 5400), (6000, 6400)], "abundance": 100}]),
)

#: ⚠ Sampling tolerance on 4,000 fragments, not a tuned constant: the binomial standard error at
#: p ≈ 0.5 is ~0.008, so 0.03 is a shade under 4 sigma. It is loose enough never to flake and far tighter
#: than the 0.5-scale error a genuine convention flip would produce.
TOLERANCE = 0.03


def _strand_model(strand_specificity: float, *, r1_sense: bool = False):
    tmp = pathlib.Path(tempfile.mkdtemp())
    scenario = Scenario("sense", genome_length=9000, seed=SEED, work_dir=tmp / "s")
    for gene, strand, transcripts in GENES:
        scenario.add_gene(gene, strand, transcripts)
    result = scenario.build_oracle(
        n_fragments=4000,
        gdna_fraction=0.0,  # every fragment is mature RNA of known strand
        sim_config=ReadSimConfig(
            frag_mean=220,
            frag_std=40,
            frag_min=100,
            frag_max=400,
            read_length=100,
            strand_specificity=strand_specificity,
            r1_sense=r1_sense,
            seed=SEED,
        ),
    )
    _stats, strand_models, _buffer, _payload = scan_and_buffer(
        str(result.bam_path), result.index, BamScanConfig(sj_strand_tag="auto", total_threads=1)
    )
    scenario.cleanup()
    return strand_models.exonic_spliced


@pytest.mark.parametrize("simulated", [1.0, 0.75, 0.5])
def test_strand_specificity_RECOVERS_the_simulated_parameter(simulated):
    """⭐ The quantity that matches the simulator's own knob, and it is already exposed.

    ``StrandModel.strand_specificity = max(p_r1_sense, p_r1_antisense)`` is direction-agnostic, which is
    what ``ReadSimConfig.strand_specificity`` also is. Measured 2026-08-11: 1.00 → 1.0000,
    0.75 → 0.7494, 0.50 → 0.5034.

    ⚠ These read 0.7701 / 0.5020 until then, and the same pair was restated in ``TRAPS.md`` — two homes
    for one measured number, and both had drifted. The prose there now points here instead.

    ⛔ **This is the comparison the record should always have made.** Comparing the simulated parameter
    against ``rna_sense_frac`` instead is what produced a phantom sign bug, twice.
    """
    model = _strand_model(simulated)
    assert model.n_observations > 0, "no spliced strand observations; the fixture proves nothing"
    assert abs(model.strand_specificity - simulated) < TOLERANCE, (
        f"simulated strand specificity {simulated} but the model recovered "
        f"{model.strand_specificity:.4f}. This is the direction-AGNOSTIC quantity and it is the one that "
        f"must match the simulator's knob."
    )


def test_the_simulator_emits_an_R1_ANTISENSE_library_so_the_SENSE_fraction_is_LOW():
    """⛔ The fact that was mistaken for a sign error, pinned so it cannot be mistaken again.

    ``p_r1_sense`` is ``P(align_strand == the junction's strand)``. The simulator emits R1-antisense
    (dUTP-style), so a perfectly stranded library has a sense fraction of **~0**, not ~1. That is the
    protocol, not a flip: ``StrandModel``'s own docstring gives ≈0.05 for TruSeq dUTP and ≈0.95 for KAPA.

    ⚠ So ``rna_sense_frac ≈ 0.01`` on a "0.99 stranded" simulated library is the correct reading, and
    forcing it to 0.99 substitutes a different quantity — which -addendum measured as
    **166× worse** on a zero-gDNA library.
    """
    model = _strand_model(1.0)
    assert model.p_r1_sense < TOLERANCE, (
        f"p_r1_sense is {model.p_r1_sense:.4f}. The simulator emits R1-antisense reads, so a perfectly "
        f"stranded library must have a sense fraction near ZERO. A value near 1 would mean the emission "
        f"orientation changed, and every strand number in the docs would need re-deriving."
    )
    assert model.p_r1_antisense > 1.0 - TOLERANCE
    assert not model.read1_sense, "read1_sense must report an R1-antisense protocol as such"


def test_rna_sense_frac_IS_p_r1_sense_and_is_therefore_ALSO_low():
    """⭐ `rna_sense_frac` is the Beta posterior mean of exactly ``p_r1_sense`` — the same quantity, the
    same direction, the same convention. It is **not** mis-labelled and it needs no sign flip.

    ⭐ This is what unblocks: the second pass needs
    ``P(align_strand agrees | RNA)`` to score an unspliced fragment's competing strand hypotheses, and
    that is precisely this number, already correct.
    """
    model = _strand_model(1.0)
    balance = fit_strand_balance(model)
    # The Beta(1,1) prior pulls the MLE toward 0.5 by 1/(n+2); on 4,000 observations that is negligible.
    assert abs(balance.rna_sense_frac - model.p_r1_sense) < 0.01, (
        "rna_sense_frac must be the posterior mean of p_r1_sense; a divergence means a second strand "
        "convention has appeared between the model and the balance fit"
    )
    assert balance.rna_sense_frac < TOLERANCE


@pytest.mark.parametrize(
    "r1_sense, want_p_r1_sense",
    [(False, 0.0), (True, 1.0)],
    ids=["R1-antisense (TruSeq dUTP)", "R1-sense (KAPA Stranded)"],
)
def test_the_suite_can_produce_EITHER_protocol_direction(r1_sense, want_p_r1_sense):
    """⭐⭐ **BOTH protocols, and the coverage gap is CLOSED.** ``ReadSimConfig.r1_sense`` is the
    protocol's DIRECTION; the engine's base emission is R1-antisense (dUTP) and ``r1_sense=True`` emits
    R1-sense (KAPA Stranded).

    ⛔ **Parametrised over BOTH rather than testing the new one alone.** An R1-sense-only test passes if
    the direction is hard-wired the other way, which is the defect it is supposed to detect; running the
    pair is what makes it a test of the SWITCH rather than of one setting.

    ⛔ ``strand_specificity`` is the FIDELITY about whichever direction is targeted and must not be
    reused as the direction. ``strand_specificity=0.0`` also emits a perfectly R1-sense library — and
    re-creates the two-quantities-one-name collision this module exists to prevent, since it would make
    ``test_strand_specificity_RECOVERS_the_simulated_parameter`` read 1.0 for a knob set to 0.0.
    """
    model = _strand_model(1.0, r1_sense=r1_sense)
    assert model.n_observations > 0, "no spliced strand observations; the fixture proves nothing"
    assert bool(model.read1_sense) is r1_sense, (
        f"r1_sense={r1_sense} but p_r1_sense is {model.p_r1_sense:.4f}; the protocol DIRECTION did not "
        f"reach the strand model"
    )
    assert abs(model.p_r1_sense - want_p_r1_sense) < TOLERANCE
    # ⭐ The direction moved; the FIDELITY did not. Both protocols are perfectly stranded here.
    assert abs(model.strand_specificity - 1.0) < TOLERANCE


def test_the_two_protocols_are_EXACT_MIRRORS_at_an_IMPERFECT_fidelity():
    """⭐⭐⭐ **THE CLAIM THAT PINS THE SWITCH, and perfect fidelity cannot make it.**

    At ``strand_specificity = 1.0`` the two protocols read 0.0 and 1.0 — but so would any pair of
    hard-wired opposites, so that comparison cannot tell a real switch from two separate code paths.
    At an IMPERFECT fidelity the mirror is a much narrower target::

        p_r1_sense       0.1989   <->   0.8011      sums to exactly 1
        strand_specificity 0.8011  ==   0.8011      the FIDELITY is direction-agnostic
        n_observations     1317    ==   1317        same fragments, one RNG stream

    ⭐ The implementation earns this: an R1-sense library flips exactly the fragments the R1-antisense
    protocol would have KEPT, so the two are per-fragment mirrors drawn from one stream rather than two
    independent simulations that happen to look opposite.
    """
    anti = _strand_model(0.8, r1_sense=False)
    sense = _strand_model(0.8, r1_sense=True)

    assert anti.n_observations == sense.n_observations > 0, (
        "the two protocols must be the SAME fragments differently labelled; a differing observation "
        "count means they are two simulations, not one mirrored"
    )
    assert abs((anti.p_r1_sense + sense.p_r1_sense) - 1.0) < TOLERANCE, (
        f"p_r1_sense must mirror about ½: {anti.p_r1_sense:.4f} + {sense.p_r1_sense:.4f}"
    )
    assert abs(anti.strand_specificity - sense.strand_specificity) < TOLERANCE, (
        f"strand_specificity is the FIDELITY and must not move with the direction: "
        f"{anti.strand_specificity:.4f} vs {sense.strand_specificity:.4f}"
    )
    # ⛔ And it must be the fidelity that was ASKED for, not merely a matched pair of wrong numbers.
    assert abs(anti.strand_specificity - 0.8) < TOLERANCE
    assert not anti.read1_sense and sense.read1_sense


def _deconvolve(*, r1_sense: bool):
    """One full pipeline run with REAL gDNA present. Returns ``(true_f_gdna, result)``.

    ⛔ The gDNA is supplied as a ``GDNAConfig``, not as ``gdna_fraction``. On this scenario
    ``gdna_fraction=0.35`` silently produces **ZERO** gDNA reads — measured — which would make the
    comparison below agree perfectly for the one reason that proves nothing
    (`TRAPS: could-the-arm-have-fired`). The true fraction is therefore counted off the oracle BAM's own
    read names and asserted, rather than assumed from the knob.
    """
    tmp = pathlib.Path(tempfile.mkdtemp())
    scenario = Scenario("proto", genome_length=9000, seed=SEED, work_dir=tmp / "s")
    for gene, strand, transcripts in GENES:
        scenario.add_gene(gene, strand, transcripts)
    result = scenario.build_oracle(
        n_fragments=4000,
        gdna_config=GDNAConfig(
            abundance=200, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000
        ),
        sim_config=ReadSimConfig(
            frag_mean=220,
            frag_std=40,
            frag_min=100,
            frag_max=400,
            read_length=100,
            strand_specificity=1.0,
            r1_sense=r1_sense,
            seed=SEED,
        ),
    )
    origins = collections.Counter()
    with pysam.AlignmentFile(str(result.bam_path)) as handle:
        for record in handle:
            origins[parse_origin(record.query_name).kind] += 1
    true_f_gdna = origins["gdna"] / sum(origins.values())

    config = PipelineConfig()
    config = dataclasses.replace(config, em=dataclasses.replace(config.em, seed=0))
    calibration = run_pipeline(str(result.bam_path), result.index, config).calibration
    scenario.cleanup()
    return true_f_gdna, calibration


def test_the_DECONVOLUTION_recovers_the_SAME_BIOLOGY_under_EITHER_protocol():
    """⭐⭐⭐ **THE DELIVERABLE CLAIM: a protocol is a labelling convention, so the biology recovered
    from it must not depend on which one was used.**

    The two libraries are the same fragments with R1 and R2 exchanged. The gDNA/RNA separation is a
    statement about molecules, so it must come out the same — while ``rna_sense_frac``, which IS the
    protocol, must mirror. A pipeline that recovered a different gDNA fraction under KAPA than under
    TruSeq would be reading the convention as biology.

    ⛔ **The arm can fire.** Measured on this fixture: the true gDNA fraction is **0.937**, so there is
    a great deal to get wrong — and the earlier version of this scenario had ZERO gDNA, where both
    protocols agree at ``f_gdna ≈ 0`` for a reason that has nothing to do with the claim.

    ⚠ Not bit-identical, and it should not be: exchanging R1 and R2 changes which mate carries which
    end, so the scan sees a different record order. The tolerance is on the RECOVERED BIOLOGY.
    """
    true_anti, anti = _deconvolve(r1_sense=False)
    true_sense, sense = _deconvolve(r1_sense=True)

    assert true_anti == true_sense > 0.5, (
        f"the two runs must simulate the SAME molecules ({true_anti:.4f} vs {true_sense:.4f}), and "
        f"there must be substantial gDNA or this comparison could not have differed"
    )

    def f_gdna(calibration):
        g = calibration.library_gdna_fragments
        r = calibration.library_rna_fragments
        return g / (g + r)

    # ⭐ The protocol is READ, and it mirrors.
    assert anti.rna_sense_frac < TOLERANCE < 1.0 - TOLERANCE < sense.rna_sense_frac, (
        f"rna_sense_frac must follow the protocol: antisense {anti.rna_sense_frac:.4f}, "
        f"sense {sense.rna_sense_frac:.4f}"
    )
    assert abs((anti.rna_sense_frac + sense.rna_sense_frac) - 1.0) < TOLERANCE

    # ⭐⭐ The BIOLOGY is not. Both recover the true gDNA fraction, and they agree with each other.
    assert abs(f_gdna(anti) - true_anti) < 0.05, (
        f"R1-antisense recovered f_gdna {f_gdna(anti):.4f} against a truth of {true_anti:.4f}"
    )
    assert abs(f_gdna(sense) - true_sense) < 0.05, (
        f"R1-sense recovered f_gdna {f_gdna(sense):.4f} against a truth of {true_sense:.4f}"
    )
    assert abs(f_gdna(anti) - f_gdna(sense)) < 0.05, (
        f"the two protocols disagree about the LIBRARY's composition — {f_gdna(anti):.4f} vs "
        f"{f_gdna(sense):.4f} — so the tool is reading a labelling convention as biology"
    )


def _r1_orientation(*, r1_sense: bool, strand_specificity: float = 0.8):
    """``{qname: R1 is_reverse}`` for one protocol, on a fixed RNG stream."""
    tmp = pathlib.Path(tempfile.mkdtemp())
    scenario = Scenario("mirror", genome_length=9000, seed=SEED, work_dir=tmp / "s")
    for gene, strand, transcripts in GENES:
        scenario.add_gene(gene, strand, transcripts)
    result = scenario.build_oracle(
        n_fragments=2000,
        gdna_fraction=0.0,
        sim_config=ReadSimConfig(
            frag_mean=220,
            frag_std=40,
            frag_min=100,
            frag_max=400,
            read_length=100,
            strand_specificity=strand_specificity,
            r1_sense=r1_sense,
            seed=SEED,
        ),
    )
    out = {}
    with pysam.AlignmentFile(str(result.bam_path)) as handle:
        for record in handle:
            if record.is_read1:
                out[record.query_name] = bool(record.is_reverse)
    scenario.cleanup()
    return out


def test_the_two_protocols_are_a_PER_FRAGMENT_mirror_not_merely_a_statistical_one():
    """⛔⛔ **THE CLAIM THE STATISTICAL GATES CANNOT MAKE, and it was measured slipping through.**

    ``test_the_two_protocols_are_EXACT_MIRRORS_at_an_IMPERFECT_fidelity`` compares totals —
    ``p_r1_sense`` mirroring, the fidelity and the observation count matching. All three survive an
    implementation that draws the direction as a *different fidelity* rather than inverting the mask::

        shipped     flip = (u >= ss);  if r1_sense: flip = ~flip      complementary sets
        perturbed   flip = (u >= (1 - ss) if r1_sense else ss)        same SIZE, different FRAGMENTS

    Both flip 80 % of fragments at ``ss = 0.8``, so every statistical gate passes — measured: 9/9 green
    under the perturbation. Only a per-fragment comparison separates them.

    ⭐ The shipped rule inverts the mask, so fragment *i* is flipped in exactly one of the two libraries
    and **every** fragment's R1 comes out on the opposite strand. Measured: 2000 of 2000, 100.00 %.
    """
    anti = _r1_orientation(r1_sense=False)
    sense = _r1_orientation(r1_sense=True)
    shared = set(anti) & set(sense)
    assert len(shared) > 500, (
        f"only {len(shared)} fragments are common to both libraries; they must be the SAME molecules "
        f"on one RNG stream or this comparison means nothing"
    )
    same = [q for q in shared if anti[q] == sense[q]]
    assert not same, (
        f"{len(same)} of {len(shared)} fragments have R1 on the SAME strand under both protocols. The "
        f"two must be an exact per-fragment mirror — an R1-sense library flips precisely the fragments "
        f"the R1-antisense protocol kept. A merely statistical mirror means the direction was "
        f"implemented as a second fidelity, which is the two-quantities-one-name collision this module "
        f"exists to prevent."
    )
