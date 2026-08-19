"""⭐⭐⭐ THE SPLICE-FLUX REFRAME — which population does an BOUNDARY compare against?

**The defect.** The reframe ``r = rho_tot(dst)/rho_tot(src)`` is a COMPOSITION imputation
(`EQUATIONS.md` §3.5), so its numerator and denominator must be totals over the SAME component set. At an
BOUNDARY the accumulator holds two disjoint populations: the molecules that cross it CONTIGUOUSLY, and the
sj flux — the molecules that SPLICE there. The predecessor put the whole sj flux into one
total per slot and used it on every hop, in both directions and in both twins. Measured against oracle
truth on a two-exon toy, that inflated the intron-facing side of the two ``intron|exon`` BOUNDARIES by
**1.28×** and **1.43×** where the truth ratio is reproduced to 3 % once the term is dropped there.

**The derivation, in one sentence.** A molecule that splices at a position has its body in the exon on
exactly ONE side of it — the genomically LOW side if that position is the sj's low end, the HIGH
side if it is the sj's high end — and it never enters the other flank at all. So the total an BOUNDARY
presents to its low neighbour must count the flux of sj that START there, the total it presents to
its high neighbour the flux of sj that END there, and one number per slot cannot express either.

**Three consequences, each gated below.**

1. ⭐ The decision is **per sj**, not per boundary — an BOUNDARY can be the low end of one sj and the
   high end of another, and then both flank totals are nonzero with *different* fluxes. So the split has
   to happen where the flux is gathered (`build_region_geometry`), not where it is consumed.
2. ⭐⭐ It is keyed on the **SIDE only — DIRECTION DOES NOT ENTER.** A hop joins two adjacent slots
   ``(k, k+1)``; whichever is the source, the pair is the same pair, so ``r = rho_lo[k+1]/rho_hi[k]``
   always. Travelling low→high (mature departing — a **splice-out**) versus high→low (mature arriving — a
   **splice-in**, `DESIGN.md` §0) changes only which is numerator. ⛔ But it is NOT one array per
   direction: within ONE forward pass an BOUNDARY at a sj's low end is the DESTINATION of a hop from
   its low flank (flux INCLUDED) and the SOURCE of the next hop into its high flank (EXCLUDED).
3. ⛔⛔ It must be **GENOMIC, never donor/acceptor.** ``FLAG_DONOR_s`` marks the genomic-LOW end of an
   ``s``-strand intron on BOTH strands (`splice_graph._terminus_and_splice_events`: ``don_bit`` ←
   ``intron_start``), so on ``−`` it sits at the transcript's biological ACCEPTOR. A predicate written in
   transcript terms flips sign with the strand; ``test_a_NEG_sj_splits_IDENTICALLY_to_a_POS_one`` is
   the gate that catches it, and it is why the fields are named ``_lo``/``_hi``.

⛔ **A rule keyed on the coarse region type provably cannot do this**, which is why none is attempted:
on `toy_harness --spec splice_both_strands` three regions are simultaneously intron and exon and
``coarse_type_array`` reports **every** gene-body region as ``exon``, so the string ``intron|exon`` never
appears on the rung where the question is hardest.

⚠⚠ **PERTURBATIONS — TRAPS: self-checking-validator's second half, all EIGHT applied and each one's gate named. Two of
them fired NOTHING against the first version of this file, and both holes were exactly the shape TRAPS: name-the-observable-per-site
predicts** — a change made in three places, gated in one.

| # | perturbation | fires |
|---|---|---|
| 1 | put the whole flux in BOTH banks (the shipped behaviour) | 7 of 12: the sum identity, the placement, both-ends, the frame pair, and all three consumer gates |
| 2 | swap ``_lo``/``_hi`` at the build | 5: the placement gate, ⭐ `−`-strand, both-ends, the frame pair, the solver's publication |
| 3 | key the split on the sj's STRAND, not its genomic end | ⭐⭐⭐ **the `−`-strand gate ALONE.** This is the sign error `EQUATIONS.md` §3.5b warns about, and exactly one gate in the file can see it — which is why it is written as an equality between the two strand arms rather than as a property of either |
| 4 | pair the FORWARD relay's arrays the wrong way round | ⭐⭐ the relay mirror gate |
| 5 | pair the BACKWARD relay's arrays the wrong way round | ⭐⭐ the relay mirror gate |
| 6 | pair the COMBINE's arrays the wrong way round | the combine gate, via ``_capture['_rescale']['r']`` |
| 7 | pair ``_flank_dom``'s LEFT (or RIGHT) lift the wrong way round | ⭐ the boundary-pair gate — and NOTHING else |
| 8 | drop the split, unspliced-only everywhere | 5: the placement, both-ends, the frame pair, and both the relay and combine gates. ⚠ This is the version `EQUATIONS.md` §3.6c records as already measured WORSE, because at an BOUNDARY→EXON step the exon genuinely contains the spliced population |

⛔⛔ **THE TWO HOLES, recorded because they are the lesson.** (a) The first relay gate reconstructed the
delivered gDNA level and compared candidates — and on an unstranded fixture the composition licence is
withheld, so the level crosses UNSCALED, the reconstruction was vacuous, and perturbation 4 passed all
eleven gates. The fix was to stop reconstructing and use a SYMMETRY the relay cannot fake: on a
palindromic chain its two passes must be mirror images. (b) ``_flank_dom`` is a THIRD consumer of the
frame and perturbation 7 fired nothing, because a palindromic fixture exchanges the two flanks of every
boundary pair and the pooled fit is symmetric in them — so that gate needs its own deliberately ASYMMETRIC
chain. ⭐ Between them: three consumers, three gates, three different fixtures, and a count of gates would
have said "twelve" either way.
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from rigel.calibration.messages.relay import RelayPolicy
from rigel.calibration.sweep import solve_chain

from rigel.calibration.region_chain import BOUNDARY, REGION
from rigel.calibration.region_geometry import (
    RegionGeometry,
    init_beliefs,
    region_gdna_geometry,
    region_total_density,
)
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS
from rigel.types import Strand

from _synthetic import make_chain_parts


#: ⚠ These gates exercise HEADPOLICY's operators, so the policy is named EXPLICITLY. ``solve_chain``
#: defaults to ``SilentPolicy``, which sends nothing — every assertion below would then be vacuous, which
#: is TRAPS: could-the-arm-have-fired exactly ("check the arm COULD have changed something").
region_sweep = functools.partial(solve_chain, policy=RelayPolicy())

_N_GRID = 41


def _delta_pmf(length):
    pmf = np.zeros(int(length) + 1)
    pmf[int(length)] = 1.0
    return pmf


#: A two-exon gene: exon / intron / exon, flanked by intergenic. ⭐ The sj runs region 1 → region 3, so
#: the BOUNDARY between regions 1 and 2 is its genomic-LOW end and the one between 2 and 3 its genomic-HIGH end
#: — and the exon of that sj's transcript is on the LOW side of the first and the HIGH side of the
#: second, which is the whole asymmetry under test.
_EXON = np.uint8(BIT_EXON_POS)
_INTRON = np.uint8(BIT_INTRON_POS)
_GENE = [np.uint8(0), _EXON, _INTRON, _EXON, np.uint8(0)]


def _gene_parts(*, strand=Strand.POS, flux=400.0, sj=None):
    """The exon/intron/exon chain, with one sj spanning the intron unless told otherwise."""
    return make_chain_parts(
        _GENE,
        region_size_bp=1000.0,
        region_pos=300.0,
        boundary_pos=120.0,
        sj=[(1, 3, strand, 1e9, 1e9, flux)] if sj is None else sj,
        gdna_fl=_delta_pmf(50),
        rna_fl=_delta_pmf(80),
    )


def _boundary_slots(chain):
    return np.flatnonzero(np.asarray(chain.kind) == BOUNDARY)


def _low_high_boundaries(chain):
    """``(low_end_slot, high_end_slot)`` for the fixture's sj — the two BOUNDARIES beside the intron."""
    kind = np.asarray(chain.kind)
    # chain is N E N E N E N E N; the intron is region index 2, i.e. slot 4
    intron_slot = 4
    assert kind[intron_slot] == REGION
    return intron_slot - 1, intron_slot + 1


# ---------------------------------------------------------------------------
# 1. the split at the SOURCE — build_region_geometry
# ---------------------------------------------------------------------------


def test_the_split_ACCOUNTS_FOR_THE_WHOLE_FLUX_and_nothing_more():
    """⛔ The sum identity, both banks. A sj that fell between the two halves would silently
    disappear from every total, and one counted twice would inflate both."""
    g = _gene_parts().geometry
    np.testing.assert_allclose(
        np.asarray(g.sj_count_lo) + np.asarray(g.sj_count_hi),
        np.asarray(g.sj_count),
    )
    np.testing.assert_allclose(
        np.asarray(g.eff_sj_lo) + np.asarray(g.eff_sj_hi),
        np.asarray(g.eff_sj),
    )
    assert float(np.asarray(g.sj_count).sum()) > 0.0, "the fixture must carry flux"


def test_a_sj_flux_lands_in_LO_at_its_low_end_and_HI_at_its_high_end():
    """⭐ The placement itself, which is the derivation's whole content: at the sj's genomic-LOW end
    the exon is on the LOW side, so the flux belongs to the LOW-flank total and to nothing else."""
    parts = _gene_parts(flux=400.0)
    g, chain = parts.geometry, parts.chain
    lo_boundary, hi_boundary = _low_high_boundaries(chain)
    jc_lo = np.asarray(g.sj_count_lo)[:, 0]
    jc_hi = np.asarray(g.sj_count_hi)[:, 0]
    assert jc_lo[lo_boundary] == pytest.approx(400.0)
    assert jc_hi[lo_boundary] == 0.0
    assert jc_lo[hi_boundary] == 0.0
    assert jc_hi[hi_boundary] == pytest.approx(400.0)
    # and nowhere else — in particular not at any REGION, which stores only CONTAINED fragments
    region_slots = np.asarray(chain.kind) == REGION
    assert float(jc_lo[region_slots].sum()) == 0.0
    assert float(jc_hi[region_slots].sum()) == 0.0
    others = [s for s in _boundary_slots(chain) if s not in (lo_boundary, hi_boundary)]
    assert float(jc_lo[others].sum() + jc_hi[others].sum()) == 0.0


def test_a_NEG_sj_splits_IDENTICALLY_to_a_POS_one():
    """⭐⭐⭐ THE SIGN GATE, and it is the reason this file exists.

    The index flags the genomic-LOW end of a ``−`` intron ``FLAG_DONOR_NEG`` and its genomic-HIGH end
    ``FLAG_ACCEPTOR_NEG`` — but a ``−`` transcript is read right-to-left, so its **biological** splice
    donor is at the genomic-HIGH end. A predicate written in transcript terms therefore places the flux on
    the opposite flank on ``−`` from ``+``, and every downstream ratio is inverted at exactly the objects
    the reframe is about.

    ⭐ The correct behaviour is that NOTHING about the placement changes with the strand: only which
    TRANSCRIPT-strand column the flux is filed under. Asserted as an exact equality between the two arms.
    """
    pos = _gene_parts(strand=Strand.POS, flux=400.0).geometry
    neg = _gene_parts(strand=Strand.NEG, flux=400.0).geometry
    # the flux moves column (genome strand -> transcript strand keying) and NOT slot
    np.testing.assert_allclose(
        np.asarray(pos.sj_count_lo)[:, 0], np.asarray(neg.sj_count_lo)[:, 1]
    )
    np.testing.assert_allclose(
        np.asarray(pos.sj_count_hi)[:, 0], np.asarray(neg.sj_count_hi)[:, 1]
    )
    assert float(np.asarray(neg.sj_count_lo)[:, 0].sum()) == 0.0
    assert float(np.asarray(neg.sj_count_hi)[:, 0].sum()) == 0.0
    # ⛔ and the placement is NOT symmetric, so the equality above has teeth: swapping the two banks
    # would be a different array, not the same one.
    lo_boundary, hi_boundary = _low_high_boundaries(_gene_parts(strand=Strand.NEG).chain)
    assert np.asarray(neg.sj_count_lo)[lo_boundary, 1] > 0.0
    assert np.asarray(neg.sj_count_lo)[hi_boundary, 1] == 0.0


def test_ONE_BOUNDARY_can_be_the_LOW_end_of_one_sj_and_the_HIGH_end_of_another():
    """⭐⭐ The case that makes the split necessary rather than merely tidy — and the case a
    per-BOUNDARY rule cannot represent at all.

    Two sj MEETING at one boundary — ``A = region0 → region2`` and ``B = region1 → region3``, i.e. A's intron
    spans region 1 and B's spans region 2, so A ENDS exactly where B BEGINS. On a 5-region chain
    (``N0 E1 N2 E3 N4 E5 N6 E7 N8``, region ``k`` at slot ``2k``) A's high end is ``left(region2) = slot 3``
    and B's low end is ``right(region1) = slot 3``: **one BOUNDARY, both roles.** Both flank totals are then
    nonzero and they carry DIFFERENT fluxes with DIFFERENT divisors, which no per-BOUNDARY number can hold.

    ⚠ On `spliced_exons` and `splice_both_strands` every sj-bearing BOUNDARY has exactly one
    attachment, so neither of those rungs exercises this — it is gated here instead of assumed.
    """
    parts = make_chain_parts(
        [_EXON, _EXON, _EXON, _EXON, _EXON],
        region_size_bp=1000.0,
        region_pos=300.0,
        boundary_pos=120.0,
        sj=[
            (0, 2, Strand.POS, 1e9, 1e9, 500.0),
            (1, 3, Strand.POS, 40.0, 40.0, 90.0),
        ],
        gdna_fl=_delta_pmf(50),
        rna_fl=_delta_pmf(80),
    )
    g = parts.geometry
    jc_lo, jc_hi = np.asarray(g.sj_count_lo)[:, 0], np.asarray(g.sj_count_hi)[:, 0]
    assert jc_lo[1] == pytest.approx(500.0) and jc_hi[1] == 0.0  # A's low end
    assert jc_lo[5] == pytest.approx(0.0) and jc_hi[5] == pytest.approx(90.0)  # B's high end
    # ⭐⭐ THE SHARED BOUNDARY: A's high end and B's low end at once, each on its own bank
    assert jc_hi[3] == pytest.approx(500.0), "A's flux belongs to the HIGH flank here"
    assert jc_lo[3] == pytest.approx(90.0), "B's flux belongs to the LOW flank here"
    # ⭐ and the divisors go with them — the two sj have different reaches, so a bank that pooled
    # them would show the same number on both sides of the shared boundary.
    ej_lo, ej_hi = np.asarray(g.eff_sj_lo)[:, 0], np.asarray(g.eff_sj_hi)[:, 0]
    assert ej_hi[3] == pytest.approx(ej_lo[1])
    assert ej_lo[3] == pytest.approx(ej_hi[5])
    assert ej_lo[3] != pytest.approx(ej_hi[3]), "the two sj must have distinct divisors"
    # ⛔ and the two flank TOTALS at that boundary differ by more than rounding, which is the consequence
    rho_lo, rho_hi = region_total_density(g, np.full(int(parts.chain.n_slots), 0.4))
    assert not np.isclose(rho_lo[3], rho_hi[3])


def test_the_two_FLANK_TOTALS_differ_by_exactly_that_flanks_own_flux():
    """`region_total_density` returns ``(rho_lo, rho_hi)``: the unspliced total plus, on each side, only
    that side's flux. ⭐ At a REGION both banks are 0, so the pair collapses to one number and every
    sj-free chain is byte-identical to the predecessor."""
    parts = _gene_parts(flux=400.0)
    g, chain = parts.geometry, parts.chain
    f_g = np.full(int(chain.n_slots), 0.4)
    rho_lo, rho_hi = region_total_density(g, f_g)
    mass, eff_g = region_gdna_geometry(g)
    eff_r = np.asarray(g.eff_rna, float)
    unspliced = np.asarray(mass, float) * (
        0.4 / np.where(np.asarray(eff_g, float) > 0, eff_g, np.inf)
        + 0.6 / np.where(eff_r > 0, eff_r, np.inf)
    )
    lo_boundary, hi_boundary = _low_high_boundaries(chain)
    ej = np.asarray(g.eff_sj_lo)[lo_boundary, 0]
    assert ej > 0
    assert rho_lo[lo_boundary] - unspliced[lo_boundary] == pytest.approx(400.0 / ej)
    assert rho_hi[lo_boundary] == pytest.approx(unspliced[lo_boundary])
    assert rho_lo[hi_boundary] == pytest.approx(unspliced[hi_boundary])
    assert rho_hi[hi_boundary] - unspliced[hi_boundary] == pytest.approx(400.0 / ej)
    region_slots = np.asarray(chain.kind) == REGION
    np.testing.assert_allclose(rho_lo[region_slots], rho_hi[region_slots])


def test_a_chain_with_NO_sj_has_ONE_frame_and_is_the_falsification_arm():
    """⛔ The ``--arms base noop`` of this change: with no sj anywhere the two flank totals are
    identically equal, so the pair cannot be measuring the split — it is measuring the rebuild. If this
    fails, every other gate in the file is reading an artefact."""
    parts = _gene_parts(sj=[])
    rho_lo, rho_hi = region_total_density(parts.geometry, np.full(int(parts.chain.n_slots), 0.4))
    np.testing.assert_array_equal(rho_lo, rho_hi)


def test_the_flank_split_is_NOT_A_FACE():
    """⚠ `test_region_geometry.test_NO_FIELD_NAMES_A_FACE` bans a ``_left``/``_right`` pair, and the reason
    is that a 0-bp boundary's own measurement is one set of numbers seen identically from both sides. ⭐ These
    fields do not reintroduce that: they are ONE measurement — the sj axis's flux — partitioned by
    which sj attaches where, and the partition is a property of the SJ, not of the boundary's own
    counting. The unspliced banks stay single, which is what this asserts."""
    fields = set(RegionGeometry.__dataclass_fields__)
    for dead in ("unspliced_count_lo", "unspliced_count_hi", "eff_gdna_lo", "eff_rna_hi"):
        assert dead not in fields
    for name in fields:
        assert not name.endswith("_left") and not name.endswith("_right")


# ---------------------------------------------------------------------------
# 2. the split at the CONSUMER — the relay and the combine, which are twins
# ---------------------------------------------------------------------------


def _sweep(parts, *, kappa=0.5, n_obs=200_000):
    cap: dict = {}
    region_sweep(
        parts.chain,
        parts.statics,
        parts.geometry,
        init_beliefs(
            parts.chain, parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=_N_GRID
        ),
        parts.region_arrays,
        rna_sense_frac=kappa,
        n_rna_obs=n_obs,
        n_gdna_obs=n_obs,
        n_grid=_N_GRID,
        _capture=cap,
    )
    return cap


def test_the_SOLVER_publishes_the_pair_and_they_differ_only_at_the_sj_BOUNDARIES():
    """The frame the whole sweep runs on. ⚠ Published as two arrays precisely because an instrument
    reconstructing a hop's ``r`` must pair them by ROLE, and one array cannot express that."""
    cap = _sweep(_gene_parts(flux=400.0))
    st = cap["_uni_static"]
    rho_lo, rho_hi = np.asarray(st["rho_lo"], float), np.asarray(st["rho_hi"], float)
    lo_boundary, hi_boundary = _low_high_boundaries(_gene_parts().chain)
    differ = set(np.flatnonzero(~np.isclose(rho_lo, rho_hi)).tolist())
    assert differ == {lo_boundary, hi_boundary}, (
        f"only the two sj-bearing BOUNDARIES may differ between flanks; got {sorted(differ)}"
    )
    assert rho_lo[lo_boundary] > rho_hi[lo_boundary], "the LOW flank of a sj's low end holds the flux"
    assert rho_hi[hi_boundary] > rho_lo[hi_boundary], "the HIGH flank of its high end holds it"


def test_THE_RELAYS_TWO_PASSES_ARE_MIRROR_IMAGES_ON_A_PALINDROMIC_CHAIN():
    """⭐⭐⭐ TRAPS: name-the-observable-per-site's gate, read off the RELAY's OWN output — a combine-only gate once let a
    relay-only deletion pass the entire calibration suite — and it is exact rather than approximate.

    **The construction.** The fixture is a PALINDROME: signatures ``[·, exon, intron, exon, ·]``, uniform
    region sizes and counts, and one sj ``region1 → region3`` with equal reaches. So reflecting the chain
    left-to-right maps it onto itself, slot ``k`` ↔ slot ``n−1−k``, and the sj's genomic-LOW end maps
    onto its genomic-HIGH end.

    **The claim.** The forward relay reads its destination's LOW-flank total and its source's HIGH-flank
    one; the backward relay reads the mirror. Under the reflection those two prescriptions are the SAME
    prescription, so on a palindromic chain::

        fwd[k]  ==  bwd[n − 1 − k]        exactly, at every slot

    ⛔ Pair the forward arrays the other way round and the equality breaks at exactly the two
    sj-bearing BOUNDARIES — the forward pass then takes the flux-free total where its mirror takes the
    flux-bearing one. ⚠ It does NOT break for a one-total-per-slot implementation, which is symmetric by
    construction; that case is caught by the placement and frame gates above. The two halves of the change
    therefore have DIFFERENT gates, which is the point of TRAPS: name-the-observable-per-site rather than a gap.

    ⭐ Both arrays compared here are the relay's, before the combine has touched anything.
    """
    parts = _gene_parts(flux=400.0)
    st = _sweep(parts, kappa=0.99)["_uni_static"]
    n = int(parts.chain.n_slots)
    lo_boundary, hi_boundary = _low_high_boundaries(parts.chain)
    assert lo_boundary == n - 1 - hi_boundary, (
        "the fixture must be a palindrome for this gate to mean anything"
    )
    # the frames themselves must mirror, or the chain is not palindromic and the gate is vacuous
    np.testing.assert_allclose(
        np.asarray(st["rho_lo"], float), np.asarray(st["rho_hi"], float)[::-1]
    )
    for a, b in (("fwd_g", "bwd_g"), ("fwd_p", "bwd_p"), ("fwd_pg", "bwd_pg")):
        np.testing.assert_allclose(
            np.asarray(st[a], float),
            np.asarray(st[b], float)[::-1],
            rtol=1e-9,
            atol=0.0,
            err_msg=f"{a} is not the mirror of {b}: the relay's two passes disagree about which "
            "flank total belongs to which role",
        )
    # ⛔ and the mirror is a real constraint here, not an identity: the two flank totals genuinely
    # differ at the sj BOUNDARIES, so a mis-paired pass has somewhere to go wrong.
    assert not np.isclose(
        np.asarray(st["rho_lo"], float)[lo_boundary], np.asarray(st["rho_hi"], float)[lo_boundary]
    )


def test_THE_BOUNDARY_PAIR_LIFT_pairs_the_frames_the_same_way_ON_AN_ASYMMETRIC_CHAIN():
    """⭐ The THIRD consumer of the frame, and it needed its own gate for the reason TRAPS: name-the-observable-per-site gives.

    ``_boundary_pair`` fits ``splice_in_premise_logvar`` from the two flanking BOUNDARIES' fluxes lifted into the
    exon's frame, and the lift is the same reframe ``r`` — so it takes the same flank pair with the same
    role pairing. ⛔ Swapping ITS pairing fires none of the gates above: on the palindromic fixture the
    swap just exchanges the two flanks of every pair, and the pooled fit is symmetric in them, so it is
    genuinely inert there. The chain here is deliberately ASYMMETRIC (unequal region sizes and counts) so
    the two lifts differ, and the gate reads the published per-pair log gap ``d`` directly.

    ⚠ What this protects is a library-level VARIANCE, not a location, so its accuracy weight is
    second-order (TRAPS: a-variance-cannot-fix-a-bias) — but an unobserved term is an unobserved term.
    """
    parts = make_chain_parts(
        [np.uint8(0), _EXON, _INTRON, _EXON, np.uint8(0)],
        region_size_bp=[900.0, 400.0, 5000.0, 1100.0, 2500.0],
        region_pos=[40.0, 260.0, 90.0, 610.0, 55.0],
        boundary_pos=[30.0, 140.0, 75.0, 210.0],
        sj=[(1, 3, Strand.POS, 1e9, 1e9, 400.0)],
        gdna_fl=_delta_pmf(50),
        rna_fl=_delta_pmf(80),
    )
    cap = _sweep(parts, kappa=0.99)
    st = cap["_uni_static"]
    rho_lo, rho_hi = np.asarray(st["rho_lo"], float), np.asarray(st["rho_hi"], float)
    spl = np.asarray(st["spl_p"], float)
    left, right = np.asarray(parts.chain.left, np.int64), np.asarray(parts.chain.right, np.int64)
    is_boundary = np.asarray(parts.chain.kind) == BOUNDARY
    n = int(parts.chain.n_slots)

    def _side(nbr, num, den):
        out = np.zeros(n)
        for i in range(n):
            nb = int(nbr[i])
            if nb < 0 or not is_boundary[nb]:
                continue
            r = num[i] / den[nb] if (den[nb] > 0.0 and num[i] > 0.0) else 1.0
            out[i] = spl[nb] * r
        return out

    fl = _side(left, rho_lo, rho_hi)
    fr = _side(right, rho_hi, rho_lo)
    glv = [g for g in cap["_glv"] if int(g["strand"]) == 0]
    assert glv, "the + strand's boundary-pair fit must be published"
    ok = np.asarray(glv[0]["ok"], bool)
    assert ok.any(), "the fixture must produce at least one live boundary pair"
    expected = np.log(np.maximum(fl, 1e-300)) - np.log(np.maximum(fr, 1e-300))
    np.testing.assert_allclose(np.asarray(glv[0]["d"], float)[ok], expected[ok], rtol=1e-9)
    # ⛔ and the swapped pairing is a DIFFERENT number here, unlike on the palindrome
    swapped = np.log(np.maximum(_side(left, rho_hi, rho_lo), 1e-300)) - np.log(
        np.maximum(_side(right, rho_lo, rho_hi), 1e-300)
    )
    assert not np.allclose(expected[ok], swapped[ok]), (
        "the chain must be asymmetric enough to distinguish the two pairings"
    )


def test_THE_COMBINE_pairs_the_frames_the_same_way():
    """The other twin. ``_capture['_rescale']`` publishes each hop's ``r`` as ``_transport`` computed it, and
    the left-hand message's source is the genomic-LOW neighbour — so the same pairing must appear.
    ⚠ Two entries are published, one per direction; the left-hand one is the first."""
    parts = _gene_parts(flux=400.0)
    cap = _sweep(parts)
    st = cap["_uni_static"]
    rho_lo = np.asarray(st["rho_lo"], float)
    rho_hi = np.asarray(st["rho_hi"], float)
    pins = cap["_rescale"]
    assert len(pins) >= 2, "both directions must publish"
    lo_boundary, hi_boundary = _low_high_boundaries(parts.chain)
    left = np.asarray(parts.chain.left, np.int64)
    right = np.asarray(parts.chain.right, np.int64)
    # the LEFT-hand message: src is the low neighbour, dst presents rho_lo
    for pin, nbr, num, den in (
        (pins[0], left, rho_lo, rho_hi),
        (pins[1], right, rho_hi, rho_lo),
    ):
        r = np.asarray(pin["r"], float)
        for slot in (lo_boundary, hi_boundary):
            s = int(nbr[slot])
            if s < 0 or den[s] <= 0.0 or num[slot] <= 0.0:
                continue
            assert r[slot] == pytest.approx(num[slot] / den[s]), (
                f"slot {slot}: the combine used {r[slot]:.6g} where the role pairing gives "
                f"{num[slot] / den[s]:.6g}"
            )
            # ⛔ and the swapped pairing must be a different number, or the assertion is vacuous
            other = den[slot] / num[s] if num[s] > 0 else np.nan
            assert not np.isclose(num[slot] / den[s], other)


def test_a_sj_free_chain_SOLVES_IDENTICALLY_under_both_pairings():
    """⛔ The solver-level falsification arm. With no sj the two frames are the same array, so the
    role pairing cannot matter — if a sj-free chain's answer depends on it, the plumbing is reading
    something other than the flux split."""
    parts = _gene_parts(sj=[])
    cap = _sweep(parts)
    st = cap["_uni_static"]
    np.testing.assert_array_equal(np.asarray(st["rho_lo"], float), np.asarray(st["rho_hi"], float))
    pins = cap["_rescale"]
    rho = np.asarray(st["rho_lo"], float)
    left = np.asarray(parts.chain.left, np.int64)
    r = np.asarray(pins[0]["r"], float)
    for slot in range(int(parts.chain.n_slots)):
        s = int(left[slot])
        if s < 0 or rho[s] <= 0.0 or rho[slot] <= 0.0:
            continue
        assert r[slot] == pytest.approx(rho[slot] / rho[s])
