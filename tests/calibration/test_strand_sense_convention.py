"""What the library sense mean ``rna_sense_frac`` means, and the fit that produces it.

Two different quantities are called strand specificity. ``ReadSimConfig.strand_specificity`` is the probability an RNA
fragment preserves correct read orientation: protocol fidelity, direction-agnostic. ``StrandModel.p_r1_sense`` is
``P(align_strand == the sj's strand)``, which is directional. The simulator emits R1-antisense, the most common real
protocol, so a dUTP library at 99 % fidelity has a sense fraction near 0.01, and comparing the two would look like a
sign error without being one. The first half pins the direction-agnostic quantity, the protocol switch, and the
per-fragment mirror between the two directions; the second tests ``fit_strand_balance``, the Beta(1,1)-smoothed
posterior-predictive sense mean.
"""

from __future__ import annotations

import collections
import dataclasses
from types import SimpleNamespace

import numpy as np
import pysam
import pytest

from rigel.calibration.strand_balance import StrandBalance, fit_strand_balance
from rigel.config import BamScanConfig, PipelineConfig
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.sim.read_name import parse_origin


SEED = 11

#: Two genes, one per strand. With a single-strand locus ``p_r1_sense`` is still defined, but a convention error that
#: swapped the comparison's operands would be invisible, since both operands would move together.
GENES = (
    ("g1", "+", [{"t_id": "t1", "exons": [(1000, 1400), (2000, 2400)], "abundance": 100}]),
    ("g2", "-", [{"t_id": "t2", "exons": [(5000, 5400), (6000, 6400)], "abundance": 100}]),
)

#: Sampling tolerance on 4,000 fragments: the binomial standard error at p ≈ 0.5 is ~0.008, so 0.03 is just under 4
#: sigma, loose enough never to flake and far tighter than the 0.5-scale error a convention flip would produce.
TOLERANCE = 0.03


def _strand_model(strand_specificity: float, *, r1_sense: bool = False):
    # the scenario owns its working directory, so `cleanup()` below removes it
    scenario = Scenario("sense", genome_length=9000, seed=SEED)
    for gene, strand, transcripts in GENES:
        scenario.add_gene(gene, strand, transcripts)
    result = scenario.build_oracle(
        n_fragments=4000,  # no gDNA: every fragment is mature RNA of known strand
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
def test_strand_specificity_recovers_the_simulated_parameter(simulated):
    """``StrandModel.strand_specificity = max(p_r1_sense, p_r1_antisense)`` is direction-agnostic, like
    ``ReadSimConfig.strand_specificity``, so the fit recovers the simulator's knob across the fidelity range. This is
    the right comparison; comparing the knob against ``rna_sense_frac`` would not be.
    """
    model = _strand_model(simulated)
    assert model.n_observations > 0, "no spliced strand observations; the fixture proves nothing"
    assert abs(model.strand_specificity - simulated) < TOLERANCE, (
        f"simulated strand specificity {simulated} but the model recovered "
        f"{model.strand_specificity:.4f}. This is the direction-AGNOSTIC quantity and it is the one that "
        f"must match the simulator's knob."
    )


def test_the_simulator_emits_an_r1_antisense_library_so_the_sense_fraction_is_low():
    """``p_r1_sense`` is ``P(align_strand == the sj's strand)``. The simulator emits R1-antisense (dUTP-style), so a
    perfectly stranded library has a sense fraction of ~0, not ~1; ``StrandModel``'s own docstring gives ≈0.05 for
    TruSeq dUTP and ≈0.95 for KAPA. So ``rna_sense_frac ≈ 0.01`` on a "0.99 stranded" simulated library is correct, and
    forcing it to 0.99 would read a gDNA fraction near a half on a library with no gDNA.
    """
    model = _strand_model(1.0)
    assert model.p_r1_sense < TOLERANCE, (
        f"p_r1_sense is {model.p_r1_sense:.4f}. The simulator emits R1-antisense reads, so a perfectly "
        f"stranded library must have a sense fraction near ZERO. A value near 1 would mean the emission "
        f"orientation changed, and every strand number in the docs would need re-deriving."
    )
    assert model.p_r1_antisense > 1.0 - TOLERANCE
    assert not model.read1_sense, "read1_sense must report an R1-antisense protocol as such"


def test_rna_sense_frac_is_p_r1_sense_and_is_therefore_also_low():
    """`rna_sense_frac` is the Beta posterior mean of exactly ``p_r1_sense``: the same quantity, direction and
    convention, needing no sign flip. The second pass scores an unspliced fragment's competing strand hypotheses with
    ``P(align_strand agrees | RNA)``, which is this number.
    """
    model = _strand_model(1.0)
    balance = fit_strand_balance(model)
    # The Beta(1,1) prior pulls the MLE toward 0.5 by 1/(n+2), negligible on 4,000 observations.
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
def test_the_suite_can_produce_either_protocol_direction(r1_sense, want_p_r1_sense):
    """``ReadSimConfig.r1_sense`` is the protocol's direction: the engine's base emission is R1-antisense (dUTP) and
    ``r1_sense=True`` emits R1-sense (KAPA Stranded). Both are parametrised, since an R1-sense-only test would pass if
    the direction were hard-wired the other way; the pair tests the switch. ``strand_specificity`` is the fidelity about
    whichever direction is targeted and must not double as the direction: ``strand_specificity=0.0`` also emits a
    perfectly R1-sense library, but would read 1.0 in the recovery test for a knob set to 0.0.
    """
    model = _strand_model(1.0, r1_sense=r1_sense)
    assert model.n_observations > 0, "no spliced strand observations; the fixture proves nothing"
    assert bool(model.read1_sense) is r1_sense, (
        f"r1_sense={r1_sense} but p_r1_sense is {model.p_r1_sense:.4f}; the protocol DIRECTION did not "
        f"reach the strand model"
    )
    assert abs(model.p_r1_sense - want_p_r1_sense) < TOLERANCE
    # The direction moved; the fidelity did not. Both protocols are perfectly stranded here.
    assert abs(model.strand_specificity - 1.0) < TOLERANCE


def test_the_two_protocols_are_exact_mirrors_at_an_imperfect_fidelity():
    """At ``strand_specificity = 1.0`` the two protocols read 0.0 and 1.0, as any pair of hard-wired opposites would. At
    an imperfect fidelity the mirror is a much narrower target::

        p_r1_sense       0.1989   <->   0.8011      sums to exactly 1
        strand_specificity 0.8011  ==   0.8011      the fidelity is direction-agnostic
        n_observations     1317    ==   1317        same fragments, one RNG stream

    An R1-sense library flips exactly the fragments the R1-antisense protocol would have kept, so the two are
    per-fragment mirrors drawn from one stream.
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
    # It must be the fidelity that was asked for, not a matched pair of wrong numbers.
    assert abs(anti.strand_specificity - 0.8) < TOLERANCE
    assert not anti.read1_sense and sense.read1_sense


def _deconvolve(*, r1_sense: bool):
    """One full pipeline run with real gDNA present. Returns ``(true_f_gdna, result)``. The gDNA is supplied as a
    ``GDNAConfig``, not as ``gdna_fraction``, which on this scenario produces zero gDNA reads and would make the
    comparison agree trivially. The true fraction is counted off the oracle BAM's read names and asserted.
    """
    # the scenario owns its working directory, so `cleanup()` below removes it
    scenario = Scenario("proto", genome_length=9000, seed=SEED)
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


def test_the_deconvolution_recovers_the_same_biology_under_either_protocol():
    """A protocol is a labelling convention, so the biology recovered must not depend on which one was used. The two
    libraries are the same fragments with R1 and R2 exchanged: the gDNA/RNA separation must come out the same while
    ``rna_sense_frac``, which is the protocol, mirrors. The true gDNA fraction here is high, so there is much to get
    wrong. The results are not bit-identical: exchanging R1 and R2 changes which mate carries which end, so the scan
    sees a different record order; the tolerance is on the recovered biology.
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

    # The protocol is read, and it mirrors.
    assert anti.rna_sense_frac < TOLERANCE < 1.0 - TOLERANCE < sense.rna_sense_frac, (
        f"rna_sense_frac must follow the protocol: antisense {anti.rna_sense_frac:.4f}, "
        f"sense {sense.rna_sense_frac:.4f}"
    )
    assert abs((anti.rna_sense_frac + sense.rna_sense_frac) - 1.0) < TOLERANCE

    # The biology is not: both recover the true gDNA fraction, and they agree with each other.
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
    # the scenario owns its working directory, so `cleanup()` below removes it
    scenario = Scenario("mirror", genome_length=9000, seed=SEED)
    for gene, strand, transcripts in GENES:
        scenario.add_gene(gene, strand, transcripts)
    result = scenario.build_oracle(
        n_fragments=2000,
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


def test_the_two_protocols_are_a_per_fragment_mirror_not_merely_a_statistical_one():
    """The mirror is per fragment, not only statistical. The totals test above (``p_r1_sense`` mirroring, fidelity and
    observation count matching) would also pass an implementation that drew the direction as a different fidelity
    instead of inverting the mask::

        shipped     flip = (u >= ss);  if r1_sense: flip = ~flip      complementary sets
        alternative flip = (u >= (1 - ss) if r1_sense else ss)        same size, different fragments

    Both flip 80 % of fragments at ``ss = 0.8``. The shipped rule inverts the mask, so fragment *i* is flipped in
    exactly one of the two libraries and every fragment's R1 comes out on the opposite strand.
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


# ── the posterior-predictive sense mean itself, off a minimal stand-in model ──────────────────


def _strand_model_stub(p_r1_sense: float, n_observations: int):
    """Minimal StrandModels stand-in: the posterior-predictive needs only these two."""
    return SimpleNamespace(p_r1_sense=p_r1_sense, n_observations=n_observations)


def test_posterior_mean():
    # n_obs=10, p_r1_sense=0.8 -> n_same=8 -> Beta(9, 3): kappa = 9/12 = 0.75.
    sb = fit_strand_balance(_strand_model_stub(0.8, 10))
    np.testing.assert_allclose(sb.rna_sense_frac, 9.0 / 12.0)
    assert sb.n_observations == 10
    assert not sb.fallback_used


def test_dense_converges_to_the_mle():
    # Abundant spliced reads: the Laplace prior washes out and kappa -> p.
    sb = fit_strand_balance(_strand_model_stub(0.8, 200_000))
    np.testing.assert_allclose(sb.rna_sense_frac, 0.8, atol=1e-4)


def test_kappa_never_at_the_bound():
    # Even a "perfectly stranded" point estimate is pulled off 0/1 by the Beta(1,1) prior.
    for p in (0.0, 1.0):
        sb = fit_strand_balance(_strand_model_stub(p, 10))
        assert 0.0 < sb.rna_sense_frac < 1.0
    np.testing.assert_allclose(
        fit_strand_balance(_strand_model_stub(1.0, 10)).rna_sense_frac, 11.0 / 12.0
    )


def test_sparse_is_pulled_toward_one_half():
    # The prior dominates when there is almost no evidence, and less as evidence arrives.
    sb1 = fit_strand_balance(_strand_model_stub(1.0, 1))  # n_same=1 -> Beta(2, 1) -> 2/3
    sb5 = fit_strand_balance(_strand_model_stub(1.0, 5))  # -> 6/7
    np.testing.assert_allclose(sb1.rna_sense_frac, 2.0 / 3.0)
    np.testing.assert_allclose(sb5.rna_sense_frac, 6.0 / 7.0)
    assert sb1.rna_sense_frac < sb5.rna_sense_frac


def test_zero_spliced_is_symmetric_fallback():
    # No spliced reads -> Beta(1,1): kappa=0.5 (channel neutral), fallback flagged.
    sb = fit_strand_balance(_strand_model_stub(0.5, 0))
    np.testing.assert_allclose(sb.rna_sense_frac, 0.5)
    assert sb.fallback_used and sb.n_observations == 0


def test_returns_strand_balance_type():
    sb = fit_strand_balance(_strand_model_stub(0.8, 20))
    assert isinstance(sb, StrandBalance)
    assert 0.0 < sb.rna_sense_frac < 1.0
