"""Exhaustive proof of the splice-junction eligibility predicate (Phase 4-mean, Test A).

The predicate ``splice_junction_eligibility`` is pure 4-bit signature logic. This module pins it down
two ways: (1) an explicit, hand-derived truth table over the eligible pairs, asserted against the full
16×16 grid; (2) structural invariants that hold *independently* of the implementation — strand-swap
symmetry, L/R reflection, the mixed-state (retained-intron) veto, and the matched-set property. Together
they leave no pair unchecked. The companion ``test_splice_junction_realized.py`` proves the same
classification with real transcript geometry and planted read origins (Test B).
"""

from __future__ import annotations

import math

import pytest

from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    N_SIGNATURES,
    TS_NEG,
    TS_POS,
)
from rigel.calibration.splice_junction import (
    SpliceJunction,
    boundary_gdna_fraction,
    splice_junction_eligibility,
)

# Readable signature aliases.
Ep, En, Ip, In = BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG
POS, NEG = (TS_POS,), (TS_NEG,)
BOTH = tuple(sorted((TS_POS, TS_NEG)))

# Hand-derived truth table: the 30 eligible (sig_L, sig_R) -> (exon_side, strands). Everything else None.
# Each entry was reasoned from "a clean exon↔intron flip on a strand, with the imputed region's exon
# strands matching exactly the strands that splice here." See docs/calibration/count_mean_bias_design.md.
EXPECTED_ELIGIBLE: dict[tuple[int, int], tuple[str, tuple[int, ...]]] = {
    # --- 6 clean single-/both-strand junctions (no antisense overlap) ---
    (Ep, Ip): ("L", POS),
    (Ip, Ep): ("R", POS),
    (En, In): ("L", NEG),
    (In, En): ("R", NEG),
    (Ep | En, Ip | In): ("L", BOTH),  # matched AMBIG (predicate-eligible; unrealizable, see Test B)
    (Ip | In, Ep | En): ("R", BOTH),
    # --- 24 antisense-overlap cases that survive the matched-set rule ---
    # + junction, antisense feature on the non-imputed bits is a TSS/intron-end, not a junction:
    (Ep, En | Ip): ("L", POS),
    (Ep, In | Ip): ("L", POS),
    (Ep, En | Ip | In): ("L", POS),
    (Ep | In, Ip): ("L", POS),
    (Ep | In, Ip | In): ("L", POS),
    (Ep | In, En | Ip | In): ("L", POS),
    (Ip, Ep | In): ("R", POS),
    (En | Ip, Ep): ("R", POS),
    (Ip | In, Ep): ("R", POS),
    (Ip | In, Ep | In): ("R", POS),
    (En | Ip | In, Ep): ("R", POS),
    (En | Ip | In, Ep | In): ("R", POS),
    # - junction mirror:
    (En, Ep | In): ("L", NEG),
    (En, Ip | In): ("L", NEG),
    (En, Ep | Ip | In): ("L", NEG),
    (En | Ip, In): ("L", NEG),
    (En | Ip, Ip | In): ("L", NEG),
    (En | Ip, Ep | Ip | In): ("L", NEG),
    (In, En | Ip): ("R", NEG),
    (Ep | In, En): ("R", NEG),
    (Ip | In, En): ("R", NEG),
    (Ip | In, En | Ip): ("R", NEG),
    (Ep | Ip | In, En): ("R", NEG),
    (Ep | Ip | In, En | Ip): ("R", NEG),
}


def _swap_strands(sig: int) -> int:
    """Map a signature's + features to - and vice versa (the strand mirror)."""
    out = 0
    if sig & BIT_EXON_POS:
        out |= BIT_EXON_NEG
    if sig & BIT_EXON_NEG:
        out |= BIT_EXON_POS
    if sig & BIT_INTRON_POS:
        out |= BIT_INTRON_NEG
    if sig & BIT_INTRON_NEG:
        out |= BIT_INTRON_POS
    return out


def _swap_strand_class(strands: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(-s for s in strands))


# ---------------------------------------------------------------------------
# (1) Exhaustive truth-table assertion
# ---------------------------------------------------------------------------


def test_eligible_pairs_match_truth_table_exhaustively():
    for sig_l in range(N_SIGNATURES):
        for sig_r in range(N_SIGNATURES):
            got = splice_junction_eligibility(sig_l, sig_r)
            expected = EXPECTED_ELIGIBLE.get((sig_l, sig_r))
            if expected is None:
                assert got is None, f"({sig_l:#x},{sig_r:#x}) expected ineligible, got {got}"
            else:
                side, strands = expected
                assert got == SpliceJunction(exon_side=side, strands=strands), (
                    f"({sig_l:#x},{sig_r:#x}) expected {expected}, got {got}"
                )


def test_exactly_30_eligible_by_signature_logic():
    n = sum(
        splice_junction_eligibility(a, b) is not None
        for a in range(N_SIGNATURES)
        for b in range(N_SIGNATURES)
    )
    assert n == 30


def test_clean_single_strand_junctions():
    # The four textbook cases, reasoned by hand.
    assert splice_junction_eligibility(Ep, Ip) == SpliceJunction("L", POS)  # exon+ → intron+
    assert splice_junction_eligibility(Ip, Ep) == SpliceJunction("R", POS)  # intron+ → exon+
    assert splice_junction_eligibility(En, In) == SpliceJunction("L", NEG)
    assert splice_junction_eligibility(In, En) == SpliceJunction("R", NEG)


# ---------------------------------------------------------------------------
# (2) Implementation-independent structural invariants over the full grid
# ---------------------------------------------------------------------------


def test_distinct_signatures_required_is_not_assumed():
    # The partition guarantees neighbours differ, but the predicate must still be defined on the
    # diagonal: identical signatures are never a junction (no exon↔intron flip).
    for s in range(N_SIGNATURES):
        assert splice_junction_eligibility(s, s) is None


def test_strand_swap_symmetry():
    # Mirroring +/- on both regions mirrors the result's strand set; the exon side is unchanged.
    for a in range(N_SIGNATURES):
        for b in range(N_SIGNATURES):
            base = splice_junction_eligibility(a, b)
            swapped = splice_junction_eligibility(_swap_strands(a), _swap_strands(b))
            if base is None:
                assert swapped is None
            else:
                assert swapped == SpliceJunction(base.exon_side, _swap_strand_class(base.strands))


def test_left_right_reflection():
    # Reflecting the boundary (swap L/R regions) flips the exon side, same strand set.
    flip = {"L": "R", "R": "L"}
    for a in range(N_SIGNATURES):
        for b in range(N_SIGNATURES):
            base = splice_junction_eligibility(a, b)
            reflected = splice_junction_eligibility(b, a)
            if base is None:
                assert reflected is None
            else:
                assert reflected == SpliceJunction(flip[base.exon_side], base.strands)


def test_mixed_state_strand_is_never_in_the_matched_set():
    # A region that is both exon and intron on a strand (retained intron) is not a clean junction
    # on that strand, so that strand can never appear in an eligible result that uses it as a side.
    for a in range(N_SIGNATURES):
        for b in range(N_SIGNATURES):
            res = splice_junction_eligibility(a, b)
            if res is None:
                continue
            region = a if res.exon_side == "L" else b
            other = b if res.exon_side == "L" else a
            for strand, ebit, ibit in ((TS_POS, Ep, Ip), (TS_NEG, En, In)):
                if strand in res.strands:
                    # exon side: pure exon on this strand; intron side: pure intron on this strand.
                    assert region & ebit and not region & ibit
                    assert other & ibit and not other & ebit


def test_matched_set_property_holds_for_every_eligible_pair():
    # The defining rule: eligible ⇒ the imputed region's exon strands == the matched strand set.
    for a in range(N_SIGNATURES):
        for b in range(N_SIGNATURES):
            res = splice_junction_eligibility(a, b)
            if res is None:
                continue
            region = a if res.exon_side == "L" else b
            exon_strands = set()
            if region & BIT_EXON_POS:
                exon_strands.add(TS_POS)
            if region & BIT_EXON_NEG:
                exon_strands.add(TS_NEG)
            assert exon_strands == set(res.strands)


def test_ambig_exon_with_single_strand_junction_is_ineligible():
    # E+·E- exon imputed from a +-only junction would over-attribute gDNA (mature- unreferenced).
    assert splice_junction_eligibility(Ep | En, Ip) is None
    assert splice_junction_eligibility(Ip, Ep | En) is None


def test_split_exon_sides_are_ineligible():
    # + junction's exon side = L, - junction's exon side = R → no single region to impute.
    assert splice_junction_eligibility(Ep | In, Ip | En) is None


# ---------------------------------------------------------------------------
# boundary_gdna_fraction arithmetic
# ---------------------------------------------------------------------------


def test_two_term_fraction_equal_eff_lengths():
    # With equal eff lengths, the 2-term fraction is just U/(U+S).
    f = boundary_gdna_fraction(
        unspliced_gdna=30.0, unspliced_rna=0.0, spliced=10.0, eff_gdna=1.0, eff_rna=1.0
    )
    assert f == pytest.approx(30.0 / 40.0)


def test_fl_normalization_changes_the_fraction():
    # Longer gDNA fragments cross more readily per molecule → dividing by the larger eff length lowers
    # the gDNA density relative to RNA. The fraction must reflect the eff-length ratio, not raw counts.
    f = boundary_gdna_fraction(
        unspliced_gdna=30.0, unspliced_rna=0.0, spliced=10.0, eff_gdna=2.0, eff_rna=1.0
    )
    rho_g, rho_r = 30.0 / 2.0, 10.0 / 1.0
    assert f == pytest.approx(rho_g / (rho_g + rho_r))


def test_three_term_moves_nascent_to_the_rna_side():
    # Same total unspliced, but split: 20 gDNA + 10 nascent-RNA. The 3-term fraction is lower than the
    # 2-term (which would call all 30 unspliced gDNA).
    two_term = boundary_gdna_fraction(
        unspliced_gdna=30.0, unspliced_rna=0.0, spliced=10.0, eff_gdna=1.0, eff_rna=1.0
    )
    three_term = boundary_gdna_fraction(
        unspliced_gdna=20.0, unspliced_rna=10.0, spliced=10.0, eff_gdna=1.0, eff_rna=1.0
    )
    assert three_term < two_term
    assert three_term == pytest.approx(20.0 / (20.0 + 10.0 + 10.0))


def test_no_crossing_evidence_returns_nan():
    f = boundary_gdna_fraction(
        unspliced_gdna=0.0, unspliced_rna=0.0, spliced=0.0, eff_gdna=1.0, eff_rna=1.0
    )
    assert math.isnan(f)
