"""The per-transcript RNA prior WEIGHT — stage 6's estimator, gated as properties.

⭐⭐ **THE ONE DERIVATION EVERYTHING RESTS ON**, and the first gate below is it:

    density(r) = mass_rna_region[r] / ceff(r) = Σ_{t ∋ r} A_t / E_t

A transcript ``t`` with ``A_t`` fragments over effective length ``E_t`` puts ``A_t·ceff(r)/E_t`` contained
fragments in region ``r``, so the REGION's own opportunity cancels and what is left is a sum of
per-transcript rates. That is why densities at objects of wildly different sizes are commensurate at all,
why the owner's ``min``-along-the-path is a min OF something, and why the multiplier that turns the soft
min back into a count is the transcript's OWN effective length.

⛔⛔ **AND THE SOFT MIN IS OPPORTUNITY-WEIGHTED, WHICH IS A DELIBERATE DEVIATION FROM THE SPEC AS
WRITTEN.** `TRAPS: a-mean-of-ratios-inherits-the-partition` was measured while designing this prior: a
boundary appears wherever a signature changes, so how many objects a path contains is an artefact of the
annotation. An unweighted mean inherits it; a weighted one does not.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "design"))

import transcript_weights as TW  # noqa: E402
from test_transcript_path import _Index, _Regions  # noqa: E402

from rigel.calibration.effective_length import contained_eff_length  # noqa: E402
from rigel.calibration.result import CalibrationResult  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.types import Strand  # noqa: E402

#: ⭐ A DEGENERATE pmf — every fragment exactly 200 bp — so every effective length is exact arithmetic
#: (``L − 199``) and a gate can state the answer rather than record it.
FL = 200
PMF = np.zeros(FL + 1, dtype=np.float64)
PMF[FL] = 1.0


@pytest.fixture
def _patched_junctions(monkeypatch):
    """⚠ A pytest fixture does not cross modules, so this is declared here rather than imported —
    the same stub ``test_transcript_path`` uses, for the same reason: the fixture index carries its
    junction boundaries directly, and these gates are about the WEIGHT, not about the CSR builder."""
    import rigel.calibration.splice_graph as SG

    class _JA:
        def __init__(self, n):
            self.edge_row = np.arange(n, dtype=np.int64)

    monkeypatch.setattr(SG, "build_junction_edge_arrays", lambda idx: _JA(len(idx.edges_df)))
    return SG


def _calibration(
    n_regions, n_boundaries, n_junctions, *, region_mass, region_opp, boundary_mass=None, boundary_opp=None
) -> CalibrationResult:
    """A real ``CalibrationResult`` — not a stand-in, so a schema change reaches these gates."""
    z_n, z_e, z_j = (np.zeros(n) for n in (n_regions, n_boundaries, n_junctions))
    return CalibrationResult(
        mass_gdna_region=z_n.copy(),
        mass_rna_region=np.asarray(region_mass, dtype=np.float64),
        mass_gdna_boundary=z_e.copy(),
        mass_rna_boundary=z_e.copy() if boundary_mass is None else np.asarray(boundary_mass, dtype=np.float64),
        mass_rna_spliced_boundary=z_e.copy(),
        boundary_mass_per_crossing=np.ones(n_boundaries),
        mass_rna_junction=z_j.copy(),
        boundary_spliced_mass_per_crossing=np.ones(n_boundaries),
        junction_mass_per_crossing=np.ones(n_junctions),
        gdna_region_eff_len=np.ones(n_regions),
        gdna_boundary_eff_len=np.ones(n_boundaries),
        rna_region_eff_len=np.asarray(region_opp, dtype=np.float64),
        rna_boundary_eff_len=np.ones(n_boundaries)
        if boundary_opp is None
        else np.asarray(boundary_opp, dtype=np.float64),
        gdna_frac_region=z_n.copy(),
        rna_pos_frac_region=np.ones(n_regions),
        rna_neg_frac_region=z_n.copy(),
        gdna_frac_boundary=z_e.copy(),
        rna_pos_frac_boundary=np.ones(n_boundaries),
        rna_neg_frac_boundary=z_e.copy(),
        gdna_density_global=0.0,
        rna_sense_frac=1.0,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n_regions,
        n_boundaries=n_boundaries,
        n_junctions=n_junctions,
        config=CalibrationConfig(),
    )


# ──────────────────────────────────────────────────────────────────────────────
# The derivation.
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("mode", TW.MODES)
def test_a_transcript_ALONE_on_its_path_recovers_its_OWN_abundance(
    tmp_path, _patched_junctions, mode
):
    """⭐⭐ **THE DERIVATION, END TO END.** One transcript, one 1,000 bp exon cut by the partition into
    two 500 bp regions plus the boundary between them — three objects of two different kinds and three
    different opportunities. Load each with exactly the mass ``A_t·ceff(o)/E_t`` the derivation predicts
    and the recovered weight must be ``A_t`` itself, for EVERY rung of the dial: when a transcript is
    alone its path is flat, so arithmetic, geometric, harmonic and min must agree.

    ⛔ The gate is that the three objects' opportunities are **301 / 301 / 199** — deliberately unequal
    and spanning both axes — so an implementation that forgot to divide by opportunity, or divided by the
    wrong axis's, cannot land on 8,010 by accident.
    """
    bounds = [0, 1_000, 1_500, 2_000, 3_000]
    idx = _Index(tmp_path, bounds, {0: [(1_000, 2_000)]}, strands=[Strand.POS])
    regions = _Regions(bounds)

    e_t = float(contained_eff_length(np.array([1_000.0]), PMF)[0])  # 801
    ceff_region = float(contained_eff_length(np.array([500.0]), PMF)[0])  # 301
    crossing_opp = 199.0  # a 200 bp fragment crossing a line: w − 1 admissible offsets
    abundance = 8_010.0
    rate = abundance / e_t  # 10.0 fragments per admissible start

    region_mass = np.zeros(regions.n_regions)
    region_mass[1] = region_mass[2] = rate * ceff_region
    region_opp = np.full(regions.n_regions, ceff_region)
    n_boundaries = regions.n_regions - 1
    boundary_mass = np.zeros(n_boundaries)
    boundary_mass[1] = rate * crossing_opp
    cal = _calibration(
        regions.n_regions,
        n_boundaries,
        1,
        region_mass=region_mass,
        region_opp=region_opp,
        boundary_mass=boundary_mass,
        boundary_opp=np.full(n_boundaries, crossing_opp),
    )

    w = TW.build_weights(cal, regions, idx, PMF, mode=mode, opportunity="total")
    assert w[0] == pytest.approx(abundance, rel=1e-9)


def test_build_weights_ITSELF_weights_by_opportunity_not_by_object_count(
    tmp_path, _patched_junctions
):
    """⛔⛔ **THE HOLE THIS GATE WAS ADDED TO CLOSE.** The re-partition gates below call
    ``repartition_invariance`` → ``_power_mean`` directly, so replacing ``build_weights``' own
    ``opp[keep]`` with ``ones`` — the exact defect `TRAPS: a-mean-of-ratios-inherits-the-partition`
    names, and the whole reason this module deviates from the spec as written — fired **zero** gates.
    Found by perturbation, which is the only thing that could have found it.

    ⭐ A 990 bp region beside a 300 bp sliver at a tenth of its density: opportunity-weighted the sliver
    is worth 11 % of the answer, unweighted it is worth half, and the two land **2.7× apart**.
    """
    bounds = [0, 1_000, 1_990, 2_290, 3_000]
    idx = _Index(tmp_path, bounds, {0: [(1_000, 2_290)]}, strands=[Strand.POS])
    regions = _Regions(bounds)
    o_big = float(contained_eff_length(np.array([990.0]), PMF)[0])  # 791
    o_small = float(contained_eff_length(np.array([300.0]), PMF)[0])  # 101

    region_mass = np.zeros(regions.n_regions)
    region_mass[1], region_mass[2] = 3.0 * o_big, 0.3 * o_small
    region_opp = np.ones(regions.n_regions)
    region_opp[1], region_opp[2] = o_big, o_small
    n_e = regions.n_regions - 1
    # ⚠ the boundary between them carries no mass, so it drops out and the two regions are the path
    cal = _calibration(regions.n_regions, n_e, 1, region_mass=region_mass, region_opp=region_opp)

    total_opp = float(contained_eff_length(np.array([1_290.0]), PMF)[0])
    weighted = (o_big + o_small) / (o_big / 3.0 + o_small / 0.3)
    unweighted = 2.0 / (1.0 / 3.0 + 1.0 / 0.3)
    got = TW.build_weights(cal, regions, idx, PMF, mode="harmonic", opportunity="total")[0]

    assert got == pytest.approx(weighted * total_opp, rel=1e-9)
    assert weighted / unweighted > 2.5, "the fixture cannot separate the two forms"
    assert got != pytest.approx(unweighted * total_opp, rel=1e-3)


def test_build_weights_SERVES_the_opportunity_it_was_ASKED_for(tmp_path, _patched_junctions):
    """⛔⛔ **THE SECOND HOLE PERTURBATION FOUND.** Making ``opportunity="total"`` silently return the
    UNSPLICED array fired zero gates: the quantities are compared in their own test, and the derivation
    gate uses a SINGLE-exon transcript, where the two are equal by construction.

    ⭐ So this one is deliberately MULTI-exon, where the ratio is a stated number: two 1,000 bp exons
    give ``ceff(2000) = 1801`` total against ``2·ceff(1000) = 1602`` unspliced.
    """
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]
    idx = _Index(tmp_path, bounds, {0: [(1_000, 2_000), (9_000, 10_000)]}, strands=[Strand.POS])
    regions = _Regions(bounds)
    region_mass = np.zeros(regions.n_regions)
    region_mass[1], region_mass[3] = 500.0, 300.0
    cal = _calibration(
        regions.n_regions,
        regions.n_regions - 1,
        1,
        region_mass=region_mass,
        region_opp=np.full(regions.n_regions, 801.0),
    )

    kw = dict(mode="harmonic")
    w_total = TW.build_weights(cal, regions, idx, PMF, opportunity="total", **kw)[0]
    w_full = TW.build_weights(cal, regions, idx, PMF, opportunity="full", **kw)[0]
    assert w_total > 0.0, "the arm could not have fired"
    assert w_total / w_full == pytest.approx(1_801.0 / 1_602.0, rel=1e-9)


def test_the_two_opportunities_are_DIFFERENT_quantities(tmp_path, _patched_junctions):
    """⛔ ``total`` sums the exon lengths BEFORE ``contained_eff_length`` and ``full`` after, so they
    differ by exactly the intron-crossing a spliced molecule may do and an unspliced one may not.

    ⚠ Without this the two arms could silently be one arm (`TRAPS: could-the-arm-have-fired`): on a
    SINGLE-exon transcript they are equal by construction, which is most of the annotation by count.
    """
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]
    idx = _Index(
        tmp_path,
        bounds,
        {0: [(1_000, 2_000), (9_000, 10_000)], 1: [(1_000, 2_000)]},
        strands=[Strand.POS, Strand.POS],
    )
    unspliced, total = TW.transcript_opportunities(idx, PMF)
    # two 1,000 bp exons: 2·(1000−199) unspliced against (2000−199) total
    assert unspliced[0] == pytest.approx(2 * 801.0)
    assert total[0] == pytest.approx(1_801.0)
    # ⭐ and they COINCIDE on the single-exon transcript, which is why the multi-exon one is the fixture
    assert unspliced[1] == pytest.approx(total[1])


# ──────────────────────────────────────────────────────────────────────────────
# The refusals — what the weight must NOT do.
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("mode", TW.MODES)
@pytest.mark.parametrize("opportunity", TW.OPPORTUNITIES)
def test_a_path_that_is_ENTIRELY_ZERO_gets_weight_EXACTLY_ZERO(
    tmp_path, _patched_junctions, mode, opportunity
):
    """⛔⛔ **THE LOAD-BEARING CASE, NOT AN BOUNDARY CASE.** 4,579 of 8,750 annotated transcripts are silent
    at ``g00``, and a zero weight is EXACTLY absorbing — so this is the whole of the prior's ability to
    say a transcript is absent. Every one of `ROADMAP.md` §4.1's eleven refused mechanisms was a rule
    for resolving doubt that lifted an evidence-free object OFF zero, and the refused first attempt's
    per-object ``+½`` took false-positive mass 18.6 M → 41.6 M by reviving exactly this population.

    ⭐ It must hold on every rung and every multiplier: a floor introduced in any one of them is a floor.
    """
    bounds = [0, 1_000, 2_000, 3_000]
    idx = _Index(tmp_path, bounds, {0: [(1_000, 2_000)]}, strands=[Strand.POS])
    regions = _Regions(bounds)
    cal = _calibration(
        regions.n_regions,
        regions.n_regions - 1,
        1,
        region_mass=np.zeros(regions.n_regions),
        region_opp=np.full(regions.n_regions, 301.0),
    )
    w = TW.build_weights(cal, regions, idx, PMF, mode=mode, opportunity=opportunity)
    assert w[0] == 0.0


@pytest.mark.parametrize("mode", TW.MODES)
def test_a_SPLICE_JUNCTION_step_carries_NO_weight(tmp_path, _patched_junctions, mode):
    """⛔ The prior is UNSPLICED fragments only — a spliced fragment has no gDNA candidate in the EM and
    never enters the split the prior arbitrates — so a junction, which is spliced by construction, must
    not reach the weight.

    ⭐ Falsified by PERTURBATION rather than by inspection: the junction banks are loaded with a mass
    three orders of magnitude above everything else, and the weight must not move by one ulp. Asserting
    the code skips ``STEP_SPLICE_JUNCTION`` would only restate the implementation.
    """
    bounds = [0, 1_000, 2_000, 9_000, 10_000, 11_000]
    idx = _Index(tmp_path, bounds, {0: [(1_000, 2_000), (9_000, 10_000)]}, strands=[Strand.POS])
    regions = _Regions(bounds)
    n_e = regions.n_regions - 1
    region_mass = np.zeros(regions.n_regions)
    region_mass[1], region_mass[3] = 100.0, 60.0
    kw = dict(region_mass=region_mass, region_opp=np.full(regions.n_regions, 801.0))

    quiet = _calibration(regions.n_regions, n_e, 1, **kw)
    loud = _calibration(regions.n_regions, n_e, 1, **kw)
    object.__setattr__(loud, "mass_rna_junction", np.array([1e6]))
    object.__setattr__(loud, "junction_mass_per_crossing", np.array([1e3]))

    a = TW.build_weights(quiet, regions, idx, PMF, mode=mode, opportunity="total")
    b = TW.build_weights(loud, regions, idx, PMF, mode=mode, opportunity="total")
    np.testing.assert_array_equal(a, b)
    assert a[0] > 0.0, "the arm could not have fired — the quiet weight was already zero"


# ──────────────────────────────────────────────────────────────────────────────
# `TRAPS: a-mean-of-ratios-inherits-the-partition`.
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("mode", TW.MODES)
def test_the_OPPORTUNITY_WEIGHTED_soft_min_survives_a_re_partition(mode):
    """⭐⭐ The property the trap says to test directly: subdividing a 10 bp sliver into 2 and then 10
    pieces of the SAME density changes no data, so the estimate must not move."""
    got = TW.repartition_invariance(mode, weighted=True)
    assert max(got) - min(got) < 1e-9, f"{mode} drifted under a re-partition: {got}"


@pytest.mark.parametrize("mode", ["arithmetic", "geometric", "harmonic"])
def test_and_the_UNWEIGHTED_one_does_NOT_which_is_why_the_weighting_is_there(mode):
    """⛔ The other half, and the one that makes the gate above mean something. Without it a constant
    estimator would pass and the opportunity weighting would look like decoration.

    ⚠ ``min`` is excluded because it is invariant either way — it ignores weights entirely — which is
    also why it is the noisiest rung of the dial rather than the safest.
    """
    got = TW.repartition_invariance(mode, weighted=False)
    assert max(got) / min(got) > 1.5, f"{mode} did not drift unweighted: {got}"


def test_the_unweighted_drift_is_toward_the_SLIVER_not_merely_different():
    """⭐ Direction, not just magnitude: subdividing a low sliver drags an unweighted mean DOWN
    monotonically, because each new piece gets a full vote it has no data to justify."""
    got = TW.repartition_invariance("harmonic", weighted=False)
    assert got[0] > got[1] > got[2]
