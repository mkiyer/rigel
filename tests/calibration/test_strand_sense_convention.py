"""⭐ WHAT `rna_sense_frac` MEANS — and why 0.0101 on a "0.99 stranded" library is CORRECT.

    Retires: `TODO.md` rank 6   ·   Corrects: `LEDGER.md` S5.f + S5.f-addendum, `CARRY_FORWARD.md` §1
    fact 17 and §0 C4   ·   Unblocks: `SPEC_SECOND_PASS.md` §3.3 (the strand term)

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

⚠ **The 166× measurement stands and now has an explanation.** `LEDGER.md` S5.f-addendum forced κ to the
nominal 0.99 and a zero-gDNA library read ``f_gdna = 0.4992`` against the fitted value's ``0.0030``. That
is not "the mirror cancels" — it is ``0.0101`` being the **right answer** and ``0.99`` being a different
quantity substituted for it.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest

from rigel.calibration.strand_balance import fit_strand_balance
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario


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


def _strand_model(strand_specificity: float):
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
    what ``ReadSimConfig.strand_specificity`` also is. Measured: 1.00 → 1.0000, 0.75 → 0.7701,
    0.50 → 0.5020.

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
    forcing it to 0.99 substitutes a different quantity — which `LEDGER.md` S5.f-addendum measured as
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

    ⭐ This is what unblocks `SPEC_SECOND_PASS.md` §3.3: the second pass needs
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


@pytest.mark.xfail(
    reason=(
        "⛔ COVERAGE GAP, filed not fixed: the simulator hard-codes an R1-ANTISENSE emission, so no "
        "simulated condition ever exercises the R1-sense branch (KAPA-style). `strand_specificity` is a "
        "swap probability about that fixed orientation, never a choice of orientation. Real R1-sense "
        "libraries exist and nothing in the suite covers them. See `TODO.md`."
    ),
    strict=True,
)
def test_the_suite_can_produce_an_R1_SENSE_library():
    """⚠ Strict xfail: the day the simulator gains an orientation switch, this fails and is deleted."""
    model = _strand_model(1.0)
    assert model.read1_sense
