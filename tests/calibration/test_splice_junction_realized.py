"""Realized-scenario proof of splice-junction eligibility (Phase 4-mean, Test B).

Test A proves the predicate's *logic* over the abstract 16×16 grid. This module proves it against
**actual transcript geometry**: it constructs overlapping and antisense transcripts in many
combinations, builds the real region partition (``build_region_partition``), and reads off which
``(sig_L, sig_R)`` boundary pairs geometry can actually produce — the **realizable** set. Two things
fall out and are asserted:

1. **Realizability respects splice-motif asymmetry.** A genomic position is a splice donor *or* an
   acceptor, never both, so a boundary cannot be a clean junction on both strands at once. The two
   both-strand pairs the signature predicate marks eligible (``E+·E- │ I+·I-`` and its reflection) are
   therefore **never realized** — proving the predicate is permitted to stay pure signature logic
   because the partition never feeds it those pairs.
2. **The eligibility rationale is real, not conservative.** At an eligible boundary the crossing-spliced
   reads are a *clean* mature reference (only the region's exon strands splice there); at an
   antisense-contaminated boundary they are not, so the same gDNA yields a *wrong* fraction — which is
   exactly why the predicate vetoes it and the caller falls back.
"""

from __future__ import annotations

from rigel.calibration.regions import build_region_partition
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS
from rigel.calibration.splice_junction import boundary_gdna_fraction, splice_junction_eligibility
from rigel.transcript import Interval, Transcript
from rigel.types import Strand

Ep, En, Ip, In = BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG
REF = "chr1"
REF_LEN = 4000


def _tx(strand: Strand, exons: list[tuple[int, int]], t_id: str) -> Transcript:
    return Transcript(ref=REF, strand=strand, exons=[Interval(a, b) for a, b in exons], t_id=t_id)


def _boundary_pairs(transcripts: list[Transcript]) -> list[tuple[int, int]]:
    """All adjacent-region (sig_L, sig_R) pairs from the realized partition of these transcripts."""
    df = build_region_partition(transcripts, {REF: REF_LEN})
    sigs = df["signature"].to_numpy()
    return [(int(sigs[i]), int(sigs[i + 1])) for i in range(len(sigs) - 1)]


_POS3 = [(300, 700), (1100, 1500), (1900, 2300)]  # a + 3-exon gene (introns 700-1100, 1500-1900);
_PERIOD = 800  # exon(400) + intron(400): offsets that are multiples align the genes exon-for-exon.


def _realizable_pairs() -> set[tuple[int, int]]:
    """Enumerate transcript combinations (single, multi-exon, overlapping, antisense) and collect every
    boundary signature pair geometry produces — for *generic* (non-coincident) opposite-strand layouts.

    The antisense sweep skips offsets that are multiples of the gene period, which would place two
    opposite-strand genes at *identical* exon/intron coordinates — the biologically degenerate case
    (distinct genes never share every splice site; motif asymmetry forbids it). That coincident case is
    exercised separately in ``test_coincident_opposite_strand_genes_*``.
    """
    realized: set[tuple[int, int]] = set()

    # Single genes (clean, no overlap) on each strand: single- and multi-exon.
    realized.update(_boundary_pairs([_tx(Strand.POS, _POS3, "p3")]))
    realized.update(_boundary_pairs([_tx(Strand.NEG, _POS3, "n3")]))
    realized.update(_boundary_pairs([_tx(Strand.POS, [(300, 700)], "p1")]))
    realized.update(_boundary_pairs([_tx(Strand.NEG, [(300, 700)], "n1")]))

    # Retained-intron isoforms (force E_s & I_s mixed-state regions) on each strand.
    realized.update(
        _boundary_pairs([_tx(Strand.POS, _POS3, "p3"), _tx(Strand.POS, [(300, 2300)], "pRI")])
    )
    realized.update(
        _boundary_pairs([_tx(Strand.NEG, _POS3, "n3"), _tx(Strand.NEG, [(300, 2300)], "nRI")])
    )

    # + and - multi-exon genes overlapping at a fine sweep of offsets → every antisense combination.
    for off in range(-900, 901, 50):
        if off % _PERIOD == 0:
            continue  # skip identical-coordinate (coincident-junction) degenerate layouts
        neg3 = [(a + off, b + off) for a, b in _POS3]
        if any(a < 0 or b > REF_LEN for a, b in neg3):
            continue
        realized.update(_boundary_pairs([_tx(Strand.POS, _POS3, "p3"), _tx(Strand.NEG, neg3, "n3")]))
    return realized


REALIZABLE = _realizable_pairs()


# ---------------------------------------------------------------------------
# Realizability
# ---------------------------------------------------------------------------


def test_clean_single_strand_junctions_are_realizable_and_eligible():
    for pair in [(Ep, Ip), (Ip, Ep), (En, In), (In, En)]:
        assert pair in REALIZABLE, f"{pair} should be producible by a multi-exon gene"
        assert splice_junction_eligibility(*pair) is not None


def test_no_both_strand_junction_in_generic_layouts():
    # With distinct (non-coincident) opposite-strand splice sites, a boundary is a junction on at most
    # one strand — so the both-strand pairs E+·E-│I+·I- never arise and no realized pair is both-strand
    # eligible. This is the realizability content of the splice-asymmetry constraint, for any
    # biologically plausible annotation.
    assert (Ep | En, Ip | In) not in REALIZABLE
    assert (Ip | In, Ep | En) not in REALIZABLE
    for a, b in REALIZABLE:
        res = splice_junction_eligibility(a, b)
        assert res is None or len(res.strands) == 1, (
            f"realized pair ({a:#x},{b:#x}) unexpectedly both-strand eligible: {res}"
        )


def test_coincident_opposite_strand_genes_are_handled_correctly_if_they_occur():
    # The coordinate-based region partition cannot see splice *motifs*, so it WILL emit E+·E-│I+·I- if
    # two opposite-strand genes are given identical coordinates (motif asymmetry forbids this in real
    # sequence, but calibration only has coordinates). Document + prove the predicate is correct anyway:
    # both strands splice and the region carries both exon bodies, so the combined spliced reference
    # matches the combined exon body → matched-set eligible. No special-casing needed.
    pairs = _boundary_pairs([_tx(Strand.POS, _POS3, "p3"), _tx(Strand.NEG, _POS3, "n3")])
    assert (Ep | En, Ip | In) in pairs  # identical coords DO produce the both-strand boundary
    res = splice_junction_eligibility(Ep | En, Ip | In)
    assert res is not None and len(res.strands) == 2 and res.exon_side == "L"


def test_antisense_overlap_actually_occurs():
    # AMBIG exon (E+·E-) and mixed-state (E·I on one strand) regions must arise from the overlaps,
    # else the sweep proved nothing about antisense.
    seen_sigs = {s for pair in REALIZABLE for s in pair}
    assert (Ep | En) in seen_sigs, "overlapping +/- exons should produce an AMBIG-exon region"
    assert any((s & Ep and s & Ip) or (s & En and s & In) for s in seen_sigs), (
        "retained-intron isoforms should produce a mixed exon·intron region"
    )


def test_antisense_eligible_cases_are_reachable():
    # The matched-set rule admits antisense-overlap junctions (e.g. a + junction whose antisense bit is
    # a TSS, not a junction). At least one such non-clean eligible pair must be realizable.
    clean = {(Ep, Ip), (Ip, Ep), (En, In), (In, En)}
    eligible = {p for p in REALIZABLE if splice_junction_eligibility(*p) is not None}
    assert eligible - clean, "no antisense-overlap eligible boundary was realized"


def test_every_realized_eligible_pair_obeys_the_matched_set_rule():
    for a, b in REALIZABLE:
        res = splice_junction_eligibility(a, b)
        if res is None:
            continue
        region = a if res.exon_side == "L" else b
        exon_strands = {s for s, bit in ((1, Ep), (-1, En)) if region & bit}
        assert exon_strands == set(res.strands)


# ---------------------------------------------------------------------------
# The eligibility rationale, on planted read origins
# ---------------------------------------------------------------------------


def test_eligible_boundary_gives_a_clean_mature_reference():
    # A clean + junction: crossing-unspliced = gDNA; crossing-spliced = mature+ ONLY. The 2-term
    # fraction recovers the planted not-mature share.
    gdna, mature_pos = 30.0, 10.0
    f = boundary_gdna_fraction(
        unspliced_gdna=gdna, unspliced_rna=0.0, spliced=mature_pos, eff_gdna=1.0, eff_rna=1.0
    )
    assert abs(f - gdna / (gdna + mature_pos)) < 1e-12


def test_antisense_contamination_would_bias_the_fraction():
    # Were we (wrongly) to use an antisense-contaminated boundary to impute a + exon, the crossing-
    # spliced would include mature- from the antisense gene, inflating the denominator and
    # UNDER-attributing gDNA. Same gDNA, two spliced references → materially different fractions.
    # This is the harm the predicate's veto prevents.
    gdna, mature_pos, mature_neg = 30.0, 10.0, 40.0
    clean = boundary_gdna_fraction(
        unspliced_gdna=gdna, unspliced_rna=0.0, spliced=mature_pos, eff_gdna=1.0, eff_rna=1.0
    )
    contaminated = boundary_gdna_fraction(
        unspliced_gdna=gdna,
        unspliced_rna=0.0,
        spliced=mature_pos + mature_neg,
        eff_gdna=1.0,
        eff_rna=1.0,
    )
    assert contaminated < clean - 0.2  # a large, real bias — not a rounding artifact


def test_three_term_recovers_pure_gdna_when_nascent_present():
    # With nascent present, the 2-term lumps it into gDNA (over-call); the 3-term, given the strand-
    # cleaned split, moves nascent to the RNA side and recovers the pure gDNA fraction.
    gdna, nascent, mature = 30.0, 15.0, 10.0
    two_term = boundary_gdna_fraction(
        unspliced_gdna=gdna + nascent, unspliced_rna=0.0, spliced=mature, eff_gdna=1.0, eff_rna=1.0
    )
    three_term = boundary_gdna_fraction(
        unspliced_gdna=gdna, unspliced_rna=nascent, spliced=mature, eff_gdna=1.0, eff_rna=1.0
    )
    assert abs(three_term - gdna / (gdna + nascent + mature)) < 1e-12
    assert two_term > three_term  # 2-term over-attributes by the nascent share
