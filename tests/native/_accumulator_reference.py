"""rigel accumulator — the executable SPECIFICATION.

The accumulator is the per-fragment tally built during the single-pass BAM scan; everything downstream
reads only its output. This module is the reference implementation, and the native accumulator is
required to reproduce it **byte for byte**. Where this file and a document disagree, this file wins.


    Spec matrix: ``test_accumulator_spec.py``

THE MODEL
    The genome is a graph. **Regions** are half-open intervals tiling each reference, numbered in genomic
    order. The 0-bp **boundaries** between adjacent regions are **contiguous boundaries**; a **sj boundary** is a
    directed donor→acceptor link taken from the annotation. A **fragment is a path** — its aligned
    blocks, joined across mate gaps and broken by introns::

        ref   |······ n0 ······|·· n1 ··|········ n2 ········|
        region_bounds  0               100      200                  600
        boundaries                  1        2
        path        [====== block ======]                       crosses boundary 1
        path        [= block =]~~intron~~[==== block ====]       crosses nothing; uses a sj

    Regions count fragments **contained** (they fit inside); boundaries count fragments **crossing**.
    ⭐ **Each population stores only the channels something reads**, and they differ: the two banks the
    deconvolution consumes carry ``count`` = Sum 1 (integer) and a reciprocal-opportunity sum
    ``Sum 1/A(w)`` (float64 — there is no fixed point anywhere, and no ``length_sum``, deleted
    2026-08-13); the certified-RNA banks carry fewer, because nothing deconvolves a fragment that is
    already known to be RNA.

    ⚠ A ``spanning`` population — one segment covering a region whole — used to exist and was **deleted
    on evidence**: measured on the panel it reached **0 starved regions** the bounding boundaries did not
    already reach off capture, and 141 regions / 822 fragments (0.008 %) under it. Its mass is not lost,
    since a spanning fragment crosses both of the region's boundaries and is deposited there.
    ⛔ One consequence, and it is structural: **no spliced fragment now touches the region axis at all**
    (a spliced fragment can never be *contained* — both endpoints of an annotated intron are region_bounds).

WHAT EACH OBJECT'S NUMBERS MEAN
    With ``A(w)`` the number of admissible start positions — ``(ell − w + 1)₊`` contained in a region,
    ``w − 1`` crossing a 0-bp boundary — and a component of start-density ``rho`` and length
    distribution ``f``, the deposit is ``1/A(w)`` and::

        E[count]      =  rho * E_f[A(w)]
        E[sum 1/A]    =  rho * P_f(A > 0)          <- the cancellation is CONDITIONAL on its support
                      =  rho                        at a BOUNDARY: P(w >= 2) = 1 for any real library
                      =  rho * P_f(w <= ell)        at a REGION — NOT rho: a fragment with w > ell
                                                    deposits NOTHING, and P(w <= ell) is a
                                                    per-component pmf functional
                                                    (TRAPS: a-cancellation-is-conditional-on-its-support)

    ⚠ **This is why the two rules carry two names and neither is called ``density``.**
    ``inv_length_sum`` (boundary/sj) IS a density, exactly; ``inv_opportunity_sum`` (region) is a
    density SHAPE truncated by its support. One name for the two rules is how the truncation stayed
    invisible (TRAPS: two-masks-one-name).

    ⚠ A ``length_sum`` channel used to be described here as "what separates two components that share a
    mean length" — that argument was RETRACTED (at equal means its row is proportional to ``(1, 1)``
    like the others, so it separates nothing) and the banks were deleted 2026-08-13;
    ``scan_payload.py``'s docstring carries the retraction.

⚠ NO PARTITIONING. Every crossed boundary receives the FULL weight. The chance that a length-``L`` fragment
crosses a given boundary is proportional to ``L − 1`` and the deposit is ``1/(L − 1)``, so the two cancel and
every fragment length contributes equally to each boundary. Dividing by the number of boundaries crossed destroys
that cancellation and makes the answer depend on region spacing — measured up to **3.6× low**. An boundary count
is not a share of a conserved total; there is no total. Density is intensive, not extensive.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field

import numpy as np

from rigel.types import Strand


__all__ = [
    "UNSPLICED_ONLY",
    "N_STRAND_COLUMNS",
    "STRAND_COLUMNS",
    "Accumulator",
    "DeferredFragment",
    "DepositOutcome",
    "FragmentPool",
    "GapResolution",
    "GapHypothesis",
    "Partition",
    "Tally",
]

#: ⭐ THE STRAND CONVENTION, and it is the same one everywhere in this file.
#:
#: Every count and density array is ``[n, N_STRAND_COLUMNS]``, and the two columns **are** the two genome
#: strands — ``Strand.POS`` and ``Strand.NEG`` — following the index's own ``_pos``/``_neg`` column naming
#: (``reach_lo_pos``, ``reach_lo_neg``, …). Nothing is stored transcript-relative.
#:
#: **Sense / antisense is derived, never stored.** A sj boundary carries its own genomic strand, so a
#: consumer that wants transcript-relative counts computes ``sense = (fragment strand == sj
#: strand)``. Storing it instead would put two conventions into one schema — which is exactly the defect
#: this replaced.
#:
#: ⚠ ``Strand`` has **four** values, not two: OR semantics make ``POS | NEG == AMBIGUOUS``, and ``NONE``
#: means no strand at all. Only POS and NEG name a column, so the mapping is a dict and a fragment whose
#: strand is neither raises ``KeyError`` here rather than being quietly filed under one of them. The
#: earlier version of this file mapped "anything that is not POS" to the minus column, which credited such
#: a fragment to the *wrong strand* instead of rejecting it.
#:
#: ⭐ **The column is indexed by ``align_strand`` — where the read ALIGNED — never by ``sj_strand``.**
#: The two are independent: ``align_strand`` is where the read sat on the genome, ``sj_strand`` is which way an
#: intron was spliced (from its ``GT..AG`` motif). Only a spliced read has the second, and it never selects
#: a column. Mixing them is what put one array into two conventions.
STRAND_COLUMNS: dict[int, int] = {Strand.POS: 0, Strand.NEG: 1}
N_STRAND_COLUMNS = len(STRAND_COLUMNS)


#: ⭐⭐⭐ ONE NUMERIC CONVENTION: **a COUNT is an integer, a FRACTION is float64.** There is no fixed
#: point anywhere in the tally, and no scale constant to decode.
#:
#: ⛔ **The predecessor accumulated every fraction as ``round(2^32 / placements)`` in uint64**, on the
#: argument that integer addition is associative and therefore bit-identical across worker counts.
#: That argument was sound and the price it quoted was not: the ~2.6 % it cited was measured on a
#: **float32** accumulator (~3.7e-7 per cell), and float64 is ~1e-16 — 3.4e9x finer, reaching the
#: deliverable at ~1e-11, five orders below ``EMConfig.convergence_delta``.
#:
#: ⭐⭐ **And the fixed point was LESS accurate, measured against exact rational arithmetic.** On the
#: reciprocal-opportunity theorem (every admissible placement deposits ``1/A``, so each length
#: contributes exactly one density unit) the two representations miss the integer answer by::
#:
#:      region_len 151    fixed 7.0e-10    float64 5.8e-15      120,000x better
#:      region_len 400    fixed 1.7e-08    float64 1.0e-13      170,000x better
#:      region_len 1000   fixed 2.0e-07    float64 2.8e-13      714,000x better
#:
#: ⚠ The exactness the old gates asserted was a property of their FIXTURES, not of the representation:
#: ``1/2 + 1/3 + 1/6`` lands back on ``2^32`` because two rounding errors cancel, while ``1/3 + 1/3 +
#: 1/3`` is one quantum short — and float64 is exact on both.
#:
#: ⛔ **What is genuinely given up is bit-identity across worker counts**, since float addition is not
#: associative. Owner ruling 2026-08-10: one convention, and this is it.


class DepositOutcome(enum.Enum):
    """Why a fragment did or did not deposit. Every rejection is counted, never silent."""

    DEPOSITED = "deposited"
    TOO_LONG = "dropped_too_long"
    EMPTY = "dropped_empty"
    #: The fragment has no single genome strand — NONE, or AMBIGUOUS (``POS | NEG``, which
    #: ``build_fragment`` produces when the mates agree in reference orientation). The column axis IS the
    #: genome strand, so there is no column to credit. Required by design §10.3, which lists
    #: strand-undefined fragments among the denominators the accumulator must emit.
    STRAND_UNDEFINED = "dropped_strand_undefined"
    #: ⭐ **Two or more hypotheses survived, so the path is not determined** — and therefore neither is
    #: ``L``, either quantum, the pool bin, or the set of boundaries the fragment crosses. It deposits on
    #: nothing and goes to the **deferred queue**, where the second pass resolves it with the fragment-length
    # distribution and the transcript abundances.,
    #:
    #: ⚠ **Not "dropped".** The fragment is retained in full and the conservation identity is
    #: ``deposited + deferred + dropped_* == offered``. The qc key says ``deferred`` for that reason: a
    #: name that said ``dropped`` for a population that is kept is how a recoverable loss gets read as a
    #: permanent one.
    DEFERRED = "deferred_undetermined_gap"


@dataclass(frozen=True, slots=True)
class GapHypothesis:
    """ONE hypothesis about what a fragment's **unsequenced** gaps contain.

    A fragment's mate gap may hold no intron, one, or several, and which it is cannot be observed — the
    bases are not there. Each candidate transcript determines exactly one answer (its own introns lying
    inside the gaps), so the hypotheses are finite and small, and two transcripts implying the same
    introns are ONE hypothesis.

    ⭐ **The empty path is the GENOMIC hypothesis.** Cutting nothing means the gap is real template, i.e.
    the molecule is gDNA — or nascent RNA, which is the same unspliced span. That is why the accumulator
    needs no separate "could this be gDNA?" flag and why the nascent shadow transcript is not a candidate:
    it *is* this hypothesis, and it is always in the set unless something rules it out.

    ⚠ ``introns`` are the IMPLIED ones only. Introns the CIGAR actually stated are region_bound under **every**
    hypothesis and are passed to :meth:`Accumulator.deposit` separately, because they are not in doubt.
    """

    #: The implied introns as ``(start, end)`` pairs. Empty ⇒ the genomic hypothesis.
    introns: tuple[tuple[int, int], ...] = ()
    #: The strand the supporting transcripts imply. ⚠ An INFERENCE, and it is only used when no motif was
    #: sequenced anywhere on the fragment — an observed motif always wins.
    sj_strand: int = Strand.NONE
    #: Which candidate transcripts imply this path. The second pass weights hypotheses by their
    #: abundance, so it needs the ids; the first pass never reads them.
    supporting_t_inds: tuple[int, ...] = ()

    @property
    def is_unspliced(self) -> bool:
        """True for the empty path — region_bound nothing, so the gap is template."""
        return not self.introns


#: The one hypothesis every fragment has: region_bound nothing beyond what was sequenced. A fragment with no
#: unsequenced gap, or with no annotated intron in the gap it has, has exactly this and deposits with no
#: arbitration at all — ⭐ the degenerate case is the general case, not a branch.
UNSPLICED_ONLY: tuple[GapHypothesis, ...] = (GapHypothesis(),)


class GapResolution(enum.Enum):
    """⭐ The umbrella, and its subclasses — owner ruling, 2026-08-01.

    Every fragment for which the enumeration produced **at least one non-genomic hypothesis** is counted
    here: its ``L`` depends on whether a gap intron is region_bound, so it is a fragment "needing further
    partitioning". The subclasses are exhaustive and mutually exclusive, so the counts close:

        Sum(GapResolution) == the umbrella total
        the three DEFERRED_* == qc["deferred_undetermined_gap"]

    ⛔ **This is its own axis and NOT a `splice_type`.** The umbrella region_bounds ACROSS the splice census: a
    certified-RNA ``SPLICED_ANNOT`` fragment with an intron in its mate gap needs resolving exactly as
    much as an ``UNSPLICED`` one does. Putting these values on ``splice_type`` would need two labels per
    fragment and would break TRAPS: pure-and-length-censored.0's property that the splice census sums to the library.

    ⚠ These classify the ARBITRATION, not the deposit. A ``RESOLVED_*`` fragment can still be rejected
    afterwards as ``TOO_LONG`` — that is a different question and it has its own counter.

    ⛔ **THERE IS NO ``RESOLVED_UNSPLICED``, AND IT IS NOT AN OMISSION.** An earlier version of this enum
    had one, described as "the genomic hypothesis survived alone because every spliced path was longer than
    ``max_fragment_length``". That is impossible. A spliced hypothesis REGION_BOUNDS bases the genomic one keeps, so
    ``L_spliced <= L_genomic`` always; the one filter is ``L <= max_fragment_length``; therefore if the
    genomic path survives the filter **every** spliced path survives it too, and the survivor set can never
    be exactly ``{genomic}`` while a spliced path was offered — which is the condition for being in this
    census at all. The class could not be entered by any fragment.

    ⭐ It is deleted rather than left at zero, and the property it depended on is pinned directly by
    ``test_gap_hypothesis_arbitration.test_the_GENOMIC_hypothesis_is_ALWAYS_the_LONGEST`` — the *reason*,
    not the consequence, so a future filter that broke the ordering would fail there instead of silently
    making a deleted class necessary again.
    """

    #: One hypothesis survived. It region_bounds something — see the class docstring for why it cannot be the
    #: genomic one — so the gap intron is real and ``L`` excludes it.
    RESOLVED_SPLICED = "gap_resolved_spliced"
    #: ⛔ The genomic path against exactly one spliced path. The open question is **RNA or gDNA** — one
    #: bit, and it is the composition question calibration exists to answer.
    DEFERRED_RNA_OR_GDNA = "gap_deferred_rna_or_gdna"
    #: ⛔ Two or more spliced paths and no genomic one: the molecule is certified RNA (gDNA cannot be
    #: spliced) and the open question is purely **which structure**.
    DEFERRED_WHICH_INTRONS = "gap_deferred_which_introns"
    #: ⛔ Both questions at once.
    DEFERRED_BOTH = "gap_deferred_both"


@dataclass(frozen=True, slots=True)
class DeferredFragment:
    """A fragment held for the second pass, stored WHOLE.

    ⭐ The fragment is stored, never its consequences. Object ids are large, derived, and would have to be
    kept consistent with the partition; the fragment is small and replays exactly. The drain re-enters
    :meth:`Accumulator.deposit` with the chosen hypothesis, so there is no second deposit path, no
    duplicated crossing logic, and byte-identity with the native accumulator is preserved for free.
    """

    ref: int
    start: int
    end: int
    align_strand: int
    sj_strand: int
    observed_introns: tuple[tuple[int, int], ...]
    hypotheses: tuple[GapHypothesis, ...]


class FragmentPool(enum.IntEnum):
    """Fragment-length pools, chosen so that each is **pure by construction**.

    Purity is what removes the circularity: the length models are fitted from populations known to be
    one component, so nothing is ever estimated from the fragments it will later explain.

    There is deliberately **no pool** for an exonic contained fragment or a multi-boundary crossing — those
    are gDNA/RNA mixtures, and an impure pool is worse than a missing one.
    """

    DNA_INTERGENIC = 0  # contained in an intergenic region
    DNA_INTRONIC = 1  # contained in an intronic region
    DNA_INTRON_EXON = 2  # crossing exactly one boundary, flanks {intron, exon}
    DNA_INTERGENIC_EXON = 3  # crossing exactly one boundary, flanks {intergenic, exon}
    RNA_SPLICED = 4  # using an annotated sj, splice OBSERVED


#: Coarse region types, as ``signature.coarse_type_array`` emits them.
_TYPE_INTERGENIC, _TYPE_INTRON, _TYPE_EXON = 0, 1, 2

#: Flank-type pair -> splash pool. The key is sorted, so it is order-insensitive.
_SPLASH_POOL = {
    (_TYPE_INTRON, _TYPE_EXON): FragmentPool.DNA_INTRON_EXON,
    (_TYPE_INTERGENIC, _TYPE_EXON): FragmentPool.DNA_INTERGENIC_EXON,
}

#: Contained region type -> pure gDNA pool. An exonic region is a mixture and is absent by design.
_CONTAINED_POOL = {
    _TYPE_INTERGENIC: FragmentPool.DNA_INTERGENIC,
    _TYPE_INTRON: FragmentPool.DNA_INTRONIC,
}


# ---------------------------------------------------------------------------
# the partition: the region_bound axis, the region types, and the sj CSR
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Partition:
    """Everything the deposit addresses. Three axes, off by one from each other per reference.

    A reference contributing ``c`` region_bound positions owns ``c − 1`` regions and ``c − 2`` contiguous boundaries
    (its interior boundaries). A reference with no regions contributes no region_bounds at all::

        region_bounds    0        100       200       600        c = 4
        regions   [  n0  ][   n1   ][   n2   ]            c - 1 = 3
        boundaries            boundary 1    boundary 2               c - 2 = 2

    SpliceJunctions are a CSR keyed by the **donor region_bound index**, which is the index the deposit has already
    computed while locating the boundaries its path crosses. It is cheap because every annotated intron has
    both endpoints as region_bounds, so "is this intron annotated?" reduces to "are both endpoints region_bounds, and is
    the pair registered?" — and if the start is not a region_bound, the table is never consulted.
    """

    region_bounds: np.ndarray  # int64[n_region_bounds] — flat, reference-major, ascending within a reference
    ref_region_bound_offsets: np.ndarray  # int64[n_refs + 1] — CSR over region_bounds
    region_types: np.ndarray  # uint8[n_regions] — 0 intergenic / 1 intron / 2 exon
    ref_region_offsets: np.ndarray  # int64[n_refs + 1]
    ref_boundary_offsets: np.ndarray  # int64[n_refs + 1]
    sj_offsets: np.ndarray  # int64[n_region_bounds + 1] — CSR over the donor region_bound index
    sj_boundary_right: np.ndarray  # int64[n_sj] — flat region_bound index of the intron's high end
    sj_strand: np.ndarray  # int8[n_sj]

    @property
    def n_refs(self) -> int:
        return int(self.ref_region_bound_offsets.shape[0]) - 1

    @property
    def n_regions(self) -> int:
        return int(self.ref_region_offsets[-1])

    @property
    def n_boundaries(self) -> int:
        return int(self.ref_boundary_offsets[-1])

    @property
    def n_sj(self) -> int:
        return int(self.sj_boundary_right.shape[0])

    @classmethod
    def from_region_bounds(cls, region_bounds_per_ref, region_types=None, sj=()) -> "Partition":
        """Build from per-reference region_bound lists.

        ``sj`` are ``(ref, intron_start, intron_end, sj_strand)``; both endpoints must be region_bounds on
        that reference. SpliceJunction ids are assigned by sorting on ``(donor region_bound, acceptor region_bound, sj_strand)``,
        so they are a deterministic function of the partition alone.
        """
        region_bounds_per_ref = [np.asarray(c, dtype=np.int64) for c in region_bounds_per_ref]
        n_refs = len(region_bounds_per_ref)
        for r, c in enumerate(region_bounds_per_ref):
            if c.size and (c.size < 2 or np.any(np.diff(c) <= 0)):
                raise ValueError(f"reference {r}: region_bounds must strictly increase, got {c.tolist()}")

        region_bound_offsets = np.zeros(n_refs + 1, np.int64)
        region_offsets = np.zeros(n_refs + 1, np.int64)
        boundary_offsets = np.zeros(n_refs + 1, np.int64)
        for r, c in enumerate(region_bounds_per_ref):
            region_bound_offsets[r + 1] = region_bound_offsets[r] + c.size
            region_offsets[r + 1] = region_offsets[r] + max(c.size - 1, 0)
            boundary_offsets[r + 1] = boundary_offsets[r] + max(c.size - 2, 0)

        region_bounds = np.concatenate(region_bounds_per_ref) if n_refs else np.zeros(0, np.int64)
        if region_types is None:
            types = np.zeros(int(region_offsets[-1]), np.uint8)
        else:
            types = np.concatenate([np.asarray(t, np.uint8) for t in region_types])
        if types.shape[0] != int(region_offsets[-1]):
            raise ValueError(
                f"region_types has {types.shape[0]} entries; the partition has "
                f"{int(region_offsets[-1])} regions"
            )

        left_boundaries, right_boundaries, sj_strands = [], [], []
        for ref, intron_start, intron_end, sj_strand in sj:
            donor = _exact_region_bound(region_bounds, region_bound_offsets, ref, intron_start)
            acceptor = _exact_region_bound(region_bounds, region_bound_offsets, ref, intron_end)
            if donor < 0 or acceptor < 0:
                raise ValueError(
                    f"sj [{intron_start}, {intron_end}) on reference {ref} has an endpoint that "
                    f"is not a region_bound. Every annotated intron endpoint is a region_bound by construction, so this "
                    f"is a partition/annotation mismatch, not an unannotated sj."
                )
            left_boundaries.append(donor)
            right_boundaries.append(acceptor)
            sj_strands.append(int(sj_strand))
        boundary_left = np.asarray(left_boundaries, np.int64)
        boundary_right = np.asarray(right_boundaries, np.int64)
        sj_strand = np.asarray(sj_strands, np.int8)
        order = np.lexsort((sj_strand, boundary_right, boundary_left))
        boundary_left, boundary_right, sj_strand = boundary_left[order], boundary_right[order], sj_strand[order]

        n_region_bounds = int(region_bound_offsets[-1])
        sj_offsets = np.zeros(n_region_bounds + 1, np.int64)
        np.cumsum(np.bincount(boundary_left, minlength=n_region_bounds), out=sj_offsets[1:])
        return cls(
            region_bounds=region_bounds,
            ref_region_bound_offsets=region_bound_offsets,
            region_types=types,
            ref_region_offsets=region_offsets,
            ref_boundary_offsets=boundary_offsets,
            sj_offsets=sj_offsets,
            sj_boundary_right=boundary_right,
            sj_strand=sj_strand,
        )


def _exact_region_bound(region_bounds, region_bound_offsets, ref: int, position: int) -> int:
    """The flat region_bound index of ``position`` on ``ref``, or -1 if it is not a region_bound there."""
    first, last = int(region_bound_offsets[ref]), int(region_bound_offsets[ref + 1])
    if last <= first:
        return -1
    k = first + int(np.searchsorted(region_bounds[first:last], position))
    return k if k < last and int(region_bounds[k]) == position else -1


def _normalise_introns(introns, start: int, end: int) -> tuple[list[tuple[int, int]], int]:
    """The introns as a sorted, DISJOINT set clipped to ``[start, end)``, and how many were absorbed.

    ⚠ The caller is contracted to hand over a clean set, but a real BAM does not always allow it: the
    splice-motif strand tag is read once per RECORD (``XS`` for STAR, ``ts`` for minimap2 — the scanner
    auto-detects), so a pair whose mates disagree about an intron's acceptor produces two overlapping
    introns for one molecule. Measured on MO_3021: one read group in 875,670
    (``chr8:138206290-138206943``, three mutually overlapping introns).

    Normalising here — rather than trusting the caller — is what lets ``L`` be *defined* as the total
    length of the path's segments instead of computed by a second, independent formula. Two formulas for
    one quantity is how the same number ends up different in two places; with an overlap they disagreed
    by up to 1.5×, and a wide overlap drove ``L`` negative so a good fragment was silently discarded.
    """
    absorbed = 0
    merged: list[tuple[int, int]] = []
    for raw_start, raw_end in sorted((int(s), int(e)) for s, e in introns):
        intron_start, intron_end = max(raw_start, start), min(raw_end, end)
        if intron_end <= intron_start:  # zero-length, or entirely outside the fragment
            absorbed += 1
            continue
        # Overlapping OR ABUTTING introns are both malformed and both merge. Abutting means a
        # zero-length exon between them, which is physically impossible — there is no molecular
        # difference between a transcript with one and without it — so a single molecule can never
        # legitimately use two abutting introns. The index cannot produce the pair either: a
        # zero-length exon is dropped when the exon arrays are built, which fuses its two flanking
        # introns into one. So an abutting pair in a fragment is an alignment artifact.
        if merged and intron_start <= merged[-1][1]:
            previous_start, previous_end = merged[-1]
            merged[-1] = (previous_start, max(previous_end, intron_end))
            absorbed += 1
            continue
        merged.append((intron_start, intron_end))
    return merged, absorbed


def _segments(start: int, end: int, introns) -> list[tuple[int, int]]:
    """The path's contiguous genomic segments: ``[start, end)`` with the introns region_bound out.

    ``introns`` must already be sorted and disjoint (:func:`_normalise_introns`).
    """
    out, cursor = [], start
    for intron_start, intron_end in introns:
        if intron_start > cursor:
            out.append((cursor, intron_start))
        cursor = intron_end
    if end > cursor:
        out.append((cursor, end))
    return out


# ---------------------------------------------------------------------------
# the tally
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class Tally:
    """The accumulator's output — per-object sums over the fragments that landed on each object.

    ``count`` = ``Sum 1`` carries the statistical power (a Beta-Binomial needs an integer); the
    reciprocal-opportunity sums carry the level under TWO rules with TWO names —
    ``inv_length_sum`` = ``Sum 1/(w−1)`` (boundary/sj: an exact model-free density,
    ``E = rho·P(w≥2) = rho``) and ``inv_opportunity_sum`` = ``Sum 1/(ell−w+1)`` (region:
    ``E = rho·P(w≤ell)``, truncated by its own support —
    TRAPS: a-cancellation-is-conditional-on-its-support); ``mass`` is the conserved fragment count, and
    sums to ONE per fragment.

    ⭐⭐ **THE POPULATIONS DO NOT ALL CARRY THE SAME CHANNELS, AND THAT IS THE DESIGN.** A channel is
    stored where something reads it and nowhere else::

        region_contained    count  inv_opportunity_sum          the mixture, on the region axis
        boundary_unspliced    count  inv_length_sum        mass   the mixture, on the boundary axis
        boundary_spliced      count                        mass   certified RNA — nothing deconvolves it
        sj          count  inv_length_sum               inv_length_sum is LIVE in second_pass

    ⛔ Six banks were removed on that rule (three ``region_spanning_*``, the two spliced-boundary moments and
    ``sj_length_sum``). ⚠ A future channel is added the same way — because a named consumer needs it,
    not because the shape would be tidier.

    The COUNT banks are integers, so a per-worker merge of them is bit-identical at any thread count;
    the reciprocal and mass banks are float64 (one numeric convention, 2026-08-11), reproducible within
    a derived tolerance rather than by byte.
    """

    region_contained_count: np.ndarray  # uint32[n_regions, 2]
    #: ⭐ float64[n_regions] — ONE column. See :meth:`Tally.zeros`: the length moments are
    #: strand-AGNOSTIC, and every consumer summed the two columns before using them.
    region_contained_inv_opportunity_sum: np.ndarray
    #: uint32[n_regions, 2] — the path's FIRST covered base, by align strand; one per accepted
    #: fragment, so ``Σ (both columns) == qc.deposited`` — THE ledger invariant. ⭐ Opportunity ``ℓ``
    #: for every fragment length; wall-blind only at the template's DOWNSTREAM end (the consumer
    #: side-selects against ``region_end_count``, whose wall is the opposite one).
    region_start_count: np.ndarray
    #: uint32[n_regions, 2] — the MIRROR: the path's LAST covered base, by align strand. Its own
    #: ledger: ``Σ == qc.deposited``. Wall-blind only at the template's UPSTREAM end.
    region_end_count: np.ndarray
    #: uint32[n_regions, 2] — regions STRICTLY SPANNED: every base covered by ONE segment and
    #: neither path endpoint inside. ⚠ Opportunity ``(w − ℓ − 1)₊`` — a pmf functional per component
    #: BY DESIGN (the length-solve channel's input; TRAPS: a-cancellation-is-conditional-on-its-support
    #: applies to it as it does to `contained`). A region a spliced path JUMPS over is not covered and
    #: books nothing. A CONTAINED fragment books start and end in one region and never span, so
    #: ``contained ≤ min(start, end)`` per region, and ``span ≡ 0`` wherever ``ℓ ≥ w_max − 1``.
    region_span_count: np.ndarray
    boundary_unspliced_count: np.ndarray  # uint32[n_boundaries, 2]
    boundary_unspliced_inv_length_sum: np.ndarray  # float64[n_boundaries] — ONE column
    #: ⭐⭐ float64[n_boundaries] — **THE CONSERVED MASS**. See :meth:`Accumulator.deposit`.
    #: A COUNT and a MASS are two different deposits and one number cannot be both: ``count`` is
    #: extensive and discrete (a Beta-Binomial needs integers) and is ``+1`` on every boundary a fragment
    #: crosses, so a fragment contributes ``max(K, 1)`` of them; this sums to ONE per fragment.
    #:
    #: ⛔ **ONE COLUMN, NOT TWO, AND THAT IS DELIBERATE.** Every other boundary bank is ``[n_boundaries, 2]``
    #: because the two genome strands are read separately — ``strand_deconv`` reads
    #: ``boundary_unspliced.count`` per column. Nothing reads a mass per strand: the mass exists to convert
    #: an object-incidence total into a fragment count, and that question has no strand in it. A column
    #: nothing reads is half the bank wasted by construction, which is the defect the surviving
    #: ``*_length_sum`` banks already carry and this one declines to repeat.
    boundary_unspliced_mass: np.ndarray
    #: uint32[n_boundaries, 2] — certified RNA: gDNA cannot be spliced. ⭐ **COUNT AND MASS ONLY.** Its two
    #: length moments were removed: nothing deconvolves a certified-RNA crossing, so they had no
    #: consumer, and `pool_lengths`' RNA_SPLICED row already carries that population's lengths.
    boundary_spliced_count: np.ndarray
    #: ⭐ uint64[n_boundaries] — the same rule, routed by the same ``spliced`` flag.
    #:
    #: ⛔ **A PARTIAL BY CONSTRUCTION, AND NOT A CONSERVATION LEDGER.** A spliced fragment's blocks that
    #: contain no interior boundary deposit nothing here — their accounting is on the sj axis — so
    #: this sums to ``crossed_block_len / L`` per fragment, never to 1. That is correct: it is a
    #: per-BOUNDARY certified-RNA term, exactly commensurate with the unspliced mass at the same boundary
    #: (both are "the share of this fragment's bases adjacent to this boundary"), which is what makes the
    #: two safe to compare there. ⛔ It is NOT "the number of spliced fragments at this boundary".
    #:
    #: ⭐ It exists so that ``mass`` is not the ONE channel that ignores the spliced/unspliced split.
    #: Every boundary channel is selected by one tuple at deposit time; a spliced fragment's mass landing
    #: in the unspliced bank would put certified RNA into the competition the prior arbitrates.
    boundary_spliced_mass: np.ndarray
    #: uint32[n_sj, 2] — ⛔⛔ **BOTH GENOME-STRAND COLUMNS ARE RETAINED, AND NOT BECAUSE ANYTHING READS
    #: THEM YET** (owner ruling, 2026-08-08). A sj is stranded by its genomic splicing MOTIF, so
    #: the strand of the *fragments* on it looks redundant, and every consumer today sums the two.
    #:
    #: ⭐ **The reason is aligner-artifact detection.** Aligners emit false-positive ``N`` CIGAR ops from
    #: plain genomic DNA. ``rigel.splice_blacklist`` catches those the sister tool ``alignable`` has
    #: enumerated by coordinate — an a-priori list, and far from complete. The EMPIRICAL detector is
    #: this column: in a stranded library a real sj inherits the global strand specificity, while
    #: an artifact deposits on BOTH strands and deviates from it. ⚠ Unstranded data cannot use it
    #: (κ = ½ leaves nothing to deviate from), which is a property of the detector, not a reason to drop
    #: the column.
    #:
    #: ⛔ The discriminating information lives ONLY in the split — a clean sj and an artifactual
    #: one carry the same total. Gated by
    #: ``test_the_sj_STRAND_SPLIT_IS_RETAINED_FOR_ALIGNER_ARTIFACT_DETECTION``.
    sj_count: np.ndarray
    #: uint64[n_sj] — ⭐ LIVE: ``second_pass.py`` scores a held fragment's sj evidence with it.
    #: ⚠ ``sj_length_sum`` is gone for the same reason the spliced boundary moments are.
    sj_inv_length_sum: np.ndarray
    #: uint64[n_sj] — ⭐⭐⭐ **THE CONSERVED MASS'S THIRD AXIS, AND IT IS WHAT MAKES A LIBRARY FRAGMENT
    #: COUNT COMPUTABLE AT ALL.** A spliced fragment's block that contains no interior boundary deposits
    #: nothing on either boundary bank, and it is not ``contained`` either — its path spans a sj, so
    #: it lies in no single region. Before this bank such a fragment existed on the incidence axis
    #: (``sj_count``) and on no conserved one, so no sum over conserved banks could count it.
    #:
    #: ⭐ **Measured on the origin-split oracle, ladder g50 capture_off: 1,222,375 of 4,830,713 RNA
    #: fragments — 25.3 % — are in that population**, against 0 of 4,997,761 gDNA fragments, because
    #: gDNA cannot splice. The gDNA side's conserved count was therefore already exact (1.000x deposited)
    #: while RNA's read 0.747x.
    #:
    #: ⛔ **THE RULE ADDS A BOUNDARY CLASS; IT DOES NOT RE-APPORTION AN EXISTING ONE.** A block that
    #: crosses at least one boundary is unchanged — its bases go to boundaries exactly as before — so
    #: ``boundary_spliced_mass`` and ``boundary_unspliced_mass`` are byte-identical to what they were, and the
    #: commensurability ``boundary_spliced_mass`` documents ("the share of this fragment's bases adjacent to
    #: this boundary", directly comparable with the unspliced mass at the same boundary) survives. Only the
    #: blocks that previously disposed of nothing are affected, and they now give their whole
    #: ``block_len / L`` to the annotated sj bounding them, shared equally.
    #:
    #: ⚠ A boundary counts only where the intron RESOLVED to an annotated sj. A block bounded
    #: solely by unannotated introns still has nowhere conserved to send its bases — the same residual
    #: the unspliced rule has always had — so the identity is exact over deposited, ANNOTATED fragments,
    #: which is one qualifier weaker than before rather than none.
    #:
    #: ⭐⭐⭐ **float64[n_sj, 2] — TWO COLUMNS, and the only mass that has any.** `accumulator.h`'s
    #: one-value ruling was reversed here and only here (owner, 2026-08-12): an ARTIFACTUAL sj
    #: accumulates SYMMETRICALLY on both strands like gDNA, so the existing strand model can detect one
    #: given a per-strand observable, and a per-strand mass is also what makes artifact filtering
    #: single-pass rather than tally-filter-re-accumulate. ⛔ The columns are `sj_count`'s columns.
    sj_mass: np.ndarray
    pool_lengths: np.ndarray  # int64[5, max_fragment_length + 1] — binned at L, once per fragment
    #: uint32[max_fragment_length + 1] — ⭐ **TRAPS: a-purity-filter-is-a-length-filter: EVERY deposited fragment, binned at its own L, with no
    # purity condition.** The five pure pools above are deliberately CONDITIONED
    #: §8: an impure pool is worse than a missing one), so they cannot serve as the unconditional anchor an
    #: empirical-Bayes shrinkage needs — which is why that anchor was taken from the SCANNER, which
    # measures length by two other rules over another population. This
    #: row removes that reason: anchor and pools become one measurement of one quantity.
    #:
    #: ⚠ **It is "unconditional GIVEN DEPOSIT", and the name says so.** It excludes what the accumulator
    #: rejects — over the length limit, ambiguous path, strand-undefined, empty — every one of which is
    #: counted in ``qc``. That is exactly the population the pools are drawn from, which is what makes it
    #: the right anchor and not merely a convenient one.
    #:
    #: ⭐ Its invariant is the same externally-checkable form as ``region_start_count``'s, and it is a
    #: DIFFERENT statement: ``Σ deposited_lengths == Σ region_start_count == qc.deposited``. The first says
    #: every fragment was binned by length, the second that every fragment was located in space.
    deposited_lengths: np.ndarray
    qc: dict[str, int] = field(default_factory=dict)
    #: ⭐ The deferred queue: fragments whose path is not determined, held WHOLE for the second pass.
    #: ⚠ Not an array, so the parity gate reaches it through :meth:`Tally.deferred_arrays` — the flattening
    #: is specified here so the two languages still compare one representation.
    deferred: list[DeferredFragment] = field(default_factory=list)
    #: ⭐ The umbrella census — one entry per :class:`GapResolution`, incremented at ONE site.
    gap_resolution: dict[str, int] = field(default_factory=dict)

    @classmethod
    def zeros(cls, n_regions: int, n_boundaries: int, n_sj: int, max_length: int) -> "Tally":
        def counts(rows):
            return np.zeros((rows, N_STRAND_COLUMNS), np.uint32)

        def fraction(rows):
            # ⭐ ONE column. The length moments are strand-AGNOSTIC: which strand a read aligned to says
            # nothing about whether the molecule was gDNA or RNA, and every consumer summed the two.
            return np.zeros(rows, np.float64)

        return cls(
            region_contained_count=counts(n_regions),
            region_contained_inv_opportunity_sum=fraction(n_regions),
            region_start_count=counts(n_regions),
            region_end_count=counts(n_regions),
            region_span_count=counts(n_regions),
            boundary_unspliced_count=counts(n_boundaries),
            boundary_unspliced_inv_length_sum=fraction(n_boundaries),
            boundary_unspliced_mass=fraction(n_boundaries),
            boundary_spliced_count=counts(n_boundaries),
            boundary_spliced_mass=fraction(n_boundaries),
            sj_count=counts(n_sj),
            sj_inv_length_sum=fraction(n_sj),
            # ⭐ TWO COLUMNS, unlike every other mass here — see the field.
            sj_mass=np.zeros((n_sj, N_STRAND_COLUMNS), np.float64),
            pool_lengths=np.zeros((len(FragmentPool), max_length + 1), np.int64),
            deposited_lengths=np.zeros(max_length + 1, np.uint32),
            qc={outcome.value: 0 for outcome in DepositOutcome}
            | {
                "unannotated_introns": 0,
                "contradictory_sj_strand": 0,
                "introns_absorbed": 0,
            },
            deferred=[],
            gap_resolution={cls.value: 0 for cls in GapResolution},
        )

    def deferred_canonical(self) -> list["DeferredFragment"]:
        """The deferred queue in its ONE canonical order — sorted on each record's own content.

        ⭐ **Factored out because the DRAIN needs the same order.** The bank is the one tally whose order
        is observable, and the second pass indexes into it: a choice vector is meaningless unless the
        producer of the scores and the consumer of the draw walk the queue identically. Defining the order
        once, here, is what makes ``choices[i]`` refer to the same fragment in both languages.

        Mirrors ``Accumulator::deferred_canonical`` in the C++.
        """
        return sorted(
            self.deferred,
            key=lambda f: (
                f.ref,
                f.start,
                f.end,
                f.align_strand,
                f.sj_strand,
                f.observed_introns,
                tuple((p.introns, p.sj_strand, p.supporting_t_inds) for p in f.hypotheses),
            ),
        )

    def deferred_arrays(self) -> dict[str, np.ndarray]:
        """The deferred queue flattened to the CSR the payload carries.

        ⭐ Specified HERE, so the readable list above and the native arrays are one representation with
        one definition rather than two implementations that have to be argued equal. The parity gate
        compares this.

        Two nested variable-length levels — fragments hold hypotheses, hypotheses hold introns — so there
        are two offset arrays. Offsets are cumulative and always start at 0, so ``n`` is
        ``len(offsets) - 1`` and an empty deferred queue is ``[0]``, never ``[]``.

        ⭐ **SORTED, and that is what makes it bit-identical at any worker count.** Every other bank is a
        sum of integers, and integer addition is associative, so a per-worker merge is exact whatever
        order the chunks arrived in. The deferred queue is a **list**, and a list has an order — so concatenating
        per-worker deferred queues would give a different byte sequence at 1, 2, 4 and 8 workers even though the
        contents are identical. Sorting on the record's own content is the canonical form: two records
        that tie on this key are identical records, so their relative order cannot be observed.
        """
        frag_fields = ("ref", "start", "end", "align_strand", "sj_strand")
        out: dict[str, list[int]] = {name: [] for name in frag_fields}
        ordered = self.deferred_canonical()
        out |= {
            "observed_intron_offsets": [0],
            "observed_introns": [],
            "hypothesis_offsets": [0],
            "hypothesis_sj_strand": [],
            "hypothesis_intron_offsets": [0],
            "hypothesis_introns": [],
            "hypothesis_t_offsets": [0],
            "hypothesis_t": [],
        }
        for frag in ordered:
            for name in frag_fields:
                out[name].append(getattr(frag, name))
            for start, end in frag.observed_introns:
                out["observed_introns"] += [start, end]
            out["observed_intron_offsets"].append(len(out["observed_introns"]) // 2)
            for path in frag.hypotheses:
                out["hypothesis_sj_strand"].append(int(path.sj_strand))
                for start, end in path.introns:
                    out["hypothesis_introns"] += [start, end]
                out["hypothesis_intron_offsets"].append(len(out["hypothesis_introns"]) // 2)
                out["hypothesis_t"] += list(path.supporting_t_inds)
                out["hypothesis_t_offsets"].append(len(out["hypothesis_t"]))
            out["hypothesis_offsets"].append(len(out["hypothesis_sj_strand"]))
        return {name: np.asarray(values, dtype=np.int64) for name, values in out.items()}


# ---------------------------------------------------------------------------
# the deposit
# ---------------------------------------------------------------------------


class Accumulator:
    """Deposit fragments into a :class:`Tally` against a :class:`Partition`.

    One public method. A fragment arrives as its genomic extent on ONE reference plus the introns inside
    it; mates are expected to be merged into that extent already, because the unsequenced gap between
    them is part of the molecule and must count toward ``L``.
    """

    __slots__ = ("partition", "max_fragment_length", "tally")

    def __init__(self, partition: Partition, max_fragment_length: int = 1000) -> None:
        self.partition = partition
        self.max_fragment_length = int(max_fragment_length)
        self.tally = Tally.zeros(
            partition.n_regions, partition.n_boundaries, partition.n_sj, self.max_fragment_length
        )

    @property
    def n_boundaries(self) -> int:
        return self.partition.n_boundaries

    def deposit(
        self,
        ref: int,
        start: int,
        end: int,
        observed_introns=(),
        align_strand: int = Strand.POS,
        sj_strand: int = Strand.NONE,
        hypotheses: tuple[GapHypothesis, ...] = UNSPLICED_ONLY,
    ) -> DepositOutcome:
        """Deposit one fragment, **or deferred queue it**; return which, and why.

        ``[start, end)`` is the full genomic extent — leftmost block start to rightmost block end, mate
        gap included. ``introns`` are the introns the CIGAR actually stated, as ``(start, end)`` pairs:
        they are region_bound under **every** hypothesis, because they are not in doubt.

        ⭐ **THE ACCUMULATOR IS THE ARBITER** (owner ruling, 2026-08-01). ``hypotheses`` is the set of
        competing answers about the fragment's *unsequenced* gaps — see :class:`GapHypothesis`, where the empty
        path is the genomic one. This method filters the set, and then:

        * **exactly one survives** → deposit it. The path is determined and nothing is in doubt;
        * **two or more survive** → :attr:`DepositOutcome.DEFERRED`, and the fragment goes WHOLE
          into :attr:`Tally.deferred` for the second pass, which has the fragment-length distribution and
          the transcript abundances to discriminate with.

        ⛔ The caller used to decide this and pass a bool, and the previous version of this docstring said
        the accumulator *could not* decide it because only the caller had the candidate list. It could not
        decide it because the caller **collapsed the answer before handing it over**. Given the set, the
        decision belongs here — where the outcome was already being reported. The rule is expected to keep
        changing as the second pass is built, which is the other reason it lives in exactly one place.

        ⭐ **The one filter is ``max_fragment_length``, and it is not a new rule.** Short-read chemistry
        does not sequence molecules past it — the same statement that makes ``TOO_LONG`` a rejection — so
        a hypothesis implying a longer ``L`` is not a molecule this library contains. Applying it to the
        *genomic* hypothesis is exactly the rule "a fragment whose genomic span exceeds the limit must be
        RNA": the genomic path's ``L`` **is** that span. ⚠ Unless the filter would empty the set, in which
        case the survivors stand and the ordinary ``TOO_LONG`` rejection counts them, as it does today.

        ⚠ The filter changes CLASSIFICATION, not just cost: hypotheses at 400 and 1200 are determined with
        it and ambiguous without. It is the second pass's likelihood applied early, at a length where the
        fragment-length distribution has no mass — defensible, and its false-positive rate is a
        measurement rather than an assumption.

        ⭐ **TWO STRANDS, AND THEY ARE INDEPENDENT.** Every read has the first; only a splice has the second.

        ``align_strand``
            The genomic strand the read **aligned** to, ``+`` or ``−``. Every read has one. It selects the
            array column and nothing else.
        ``sj_strand``
            A splice junction's strand, read from the aligner's **genomic-motif** tag (``XS`` for STAR,
            ``ts`` for minimap2 — auto-detected). A ``GT..AG`` intron is on ``+``, its reverse complement
            ``CT..AC`` on ``−``. Only spliced reads carry one, and it is used *only* to resolve an intron
            against the annotation (see :meth:`_sj_edge_id`).

        ⚠ The two are **unrelated**: an aligned strand says where the read sat, a splice strand says which
        way the intron was spliced, and neither constrains the other. Comparing them yields *sense* versus
        *antisense*, which is a derived quantity a consumer may compute and this accumulator never stores.
        ⭐ The shipped code collapses that comparison into one bool named ``primary``
        (``bam_scanner.cpp:1493``), which is how a dUTP library ended up with 0.6 % of spliced fragments in
        the column labelled *sense*.

        ⚠ ``sj_strand`` here is **observed**, per fragment. ``Partition.sj_strand`` is **annotated**, per
        sj boundary. One quantity, two sources; :meth:`_sj_edge_id` compares them.
        """
        # ⚠ Checked FIRST, and before any geometry: the strand is a property of the fragment alone, and a
        # fragment with no single genome strand has no column in any bank. The scanner's gate at
        # `bam_scanner.cpp:1474-1480` did this before the old accumulator ever saw the fragment; doing it
        # here is what lets the loss be COUNTED instead of vanishing.
        column = STRAND_COLUMNS.get(align_strand)
        if column is None:
            return self._reject(DepositOutcome.STRAND_UNDEFINED)

        p = self.partition
        first_region_bound, last_region_bound = int(p.ref_region_bound_offsets[ref]), int(p.ref_region_bound_offsets[ref + 1])
        if last_region_bound - first_region_bound < 2:
            return self._reject(DepositOutcome.EMPTY)
        region_bounds = p.region_bounds[first_region_bound:last_region_bound]

        # clip to the reference; L is the CLIPPED length, so the placement count stays consistent
        start, end = max(int(start), int(region_bounds[0])), min(int(end), int(region_bounds[-1]))
        if end <= start:
            return self._reject(DepositOutcome.EMPTY)

        # ── arbitration: which hypotheses survive, and is exactly one left? ───────────────────────
        #
        # ⚠ AFTER the strand and the clip, and the order is part of the contract. Every fragment must
        # count exactly ONCE and a fragment can fail several ways. A fragment with no genome strand is
        # not recoverable by the second pass — that pass resolves which PATH, not which strand the read
        # aligned to — so the strand rejection must win over the deferred queue. And the clip has to come first
        # because a hypothesis is filtered on its `L`, which is measured after clipping.
        scored = [
            (hypothesis, *self._hypothesis_length(start, end, observed_introns, hypothesis))
            for hypothesis in hypotheses
        ]
        survivors = [row for row in scored if row[1] <= self.max_fragment_length] or scored
        self._record_gap_resolution(hypotheses, survivors)
        if len(survivors) > 1:
            self.tally.deferred.append(
                DeferredFragment(
                    ref=int(ref),
                    start=start,
                    end=end,
                    align_strand=int(align_strand),
                    sj_strand=int(sj_strand),
                    observed_introns=tuple((int(s), int(e)) for s, e in observed_introns),
                    hypotheses=tuple(hypotheses),
                )
            )
            return self._reject(DepositOutcome.DEFERRED)

        # ⚠ `region_bound_introns` and not `introns`: these are the introns actually removed from the molecule —
        # the observed ones UNIONED with the surviving hypothesis's implied ones, normalised and clipped.
        # Naming them apart from `observed_introns` is what stops the two being confused downstream.
        hypothesis, length, region_bound_introns, absorbed = survivors[0]
        segments = _segments(start, end, region_bound_introns)
        if length <= 0:
            return self._reject(DepositOutcome.EMPTY)
        if length > self.max_fragment_length:
            return self._reject(DepositOutcome.TOO_LONG)
        self.tally.qc["introns_absorbed"] += absorbed
        # ⭐ An implied strand is used ONLY when no motif was sequenced anywhere on this fragment. An
        # observed motif is evidence; an implied strand is an inference from the annotation, and mixing
        # an inference into an observation is how `primary` went wrong.
        if sj_strand == Strand.NONE:
            sj_strand = hypothesis.sj_strand

        # ── which annotated sj does this path use? ─────────────────────────────────────────
        # ⚠ A contradictory motif strand disqualifies the whole fragment's splices, so it is checked once
        # here rather than per intron: `sj_strand` is the OR of the per-record strand-tag values, so
        # AMBIGUOUS means the mates disagreed about the same molecule. That is contradictory EVIDENCE, not
        # missing evidence, and it is counted on its own denominator — folding it into
        # `unannotated_introns` would poison the one metric whose job is measuring annotation coverage.
        # ⭐ Resolved PER INTRON POSITION, with -1 where the annotation has no such sj, because
        # the conserved mass needs to know WHICH of a block's two ends is a sj boundary. A
        # filtered list cannot answer that — dropping the unannotated entries destroys the alignment
        # between intron `i` and the gap between blocks `i` and `i+1`.
        if sj_strand == Strand.AMBIGUOUS:
            sj_id_at_gap: list[int] = [-1] * len(region_bound_introns)
            self.tally.qc["contradictory_sj_strand"] += 1
        else:
            sj_id_at_gap = [
                self._sj_edge_id(ref, intron_start, intron_end, sj_strand)
                for intron_start, intron_end in region_bound_introns
            ]
            self.tally.qc["unannotated_introns"] += sum(1 for jid in sj_id_at_gap if jid < 0)
        sj_ids = [jid for jid in sj_id_at_gap if jid >= 0]
        spliced = bool(sj_ids)

        # ⚠ The path's own first and last COVERED base, not the fragment's extent. A leading or trailing
        # intron means the molecule does not begin at `start`: with introns [(150,480)] over [150,500) the
        # path is [480,500) and lives in a different region entirely. Using the extent would attribute the
        # fragment — and the start-count invariant — to a region it never touches.
        first_base, last_base = segments[0][0], segments[-1][1] - 1

        t = self.tally
        region_base, boundary_base = int(p.ref_region_offsets[ref]), int(p.ref_boundary_offsets[ref])
        r_first = self._local_region(region_bounds, first_base)
        r_last = self._local_region(region_bounds, last_base)
        t.region_start_count[region_base + r_first, column] += 1
        t.region_end_count[region_base + r_last, column] += 1
        # ⭐ TRAPS: a-purity-filter-is-a-length-filter: the unconditional length histogram, incremented HERE — beside the start count and the
        # DEPOSITED counter — so all three describe one population by construction rather than by
        # agreement. ``length`` is already clipped to the reference and gated by the length limit above.
        t.deposited_lengths[length] += 1
        t.qc[DepositOutcome.DEPOSITED.value] += 1

        # ── crossings, per contiguous SEGMENT of the path ─────────────────────────────────────────
        # A boundary is crossed iff it lies strictly inside a segment, so per segment the crossed boundaries are a
        # contiguous index range and no container is needed. A region is SPANNED iff ONE segment crosses
        # both of its boundaries — not merely "both boundaries crossed", which would count a region the fragment
        # JUMPS OVER, whose two boundaries are touched by the two flanking segments from opposite sides.
        # ⭐ TWO channels on the spliced bank, FOUR on the unspliced one, and the asymmetry is the
        # design rather than an omission. A spliced crossing is certified RNA: nothing deconvolves it,
        # so its length moments have no consumer and are not stored. The unspliced bank is the mixture,
        # and the length banks are built on exactly its moments.
        boundary_count = t.boundary_spliced_count if spliced else t.boundary_unspliced_count
        boundary_mass = t.boundary_spliced_mass if spliced else t.boundary_unspliced_mass
        # ⚠ 0 at L == 1: a length-1 molecule cannot cross a 0-bp boundary, and 1/(L-1) would divide by zero.
        inv_boundary = 1.0 / (length - 1) if length >= 2 else 0.0
        n_crossed, sole_boundary = 0, -1
        for block, (seg_start, seg_end) in enumerate(segments):
            first = int(np.searchsorted(region_bounds, seg_start, side="right"))
            last = int(np.searchsorted(region_bounds, seg_end, side="left"))
            for boundary in range(first, last):
                boundary_count[boundary_base + boundary - 1, column] += 1
                if not spliced:
                    t.boundary_unspliced_inv_length_sum[boundary_base + boundary - 1] += inv_boundary
            # ── STRICT SPAN: regions this segment covers whole, excluding the path's own endpoint
            # regions. ⚠ `side="left"`/`- 1` rather than the crossing loop's bounds: a NON-FIRST
            # segment starts exactly ON a region_bound (its acceptor), and the region beginning there
            # IS fully covered — the crossing window would miss it (TRAPS: could-the-arm-have-fired's
            # shape, caught at design time by the three-intron worked example).
            f2 = int(np.searchsorted(region_bounds, seg_start, side="left"))
            l2 = int(np.searchsorted(region_bounds, seg_end, side="right")) - 1
            for r_loc in range(f2, l2):
                if r_loc == r_first or r_loc == r_last:
                    continue
                t.region_span_count[region_base + r_loc, column] += 1
            # ── ⭐⭐⭐ THE CONSERVED MASS, per SLICE, over ONE BOUNDARY SET ─────────────────────────
            # The crossed region_bounds split this block into `last - first + 1` slices. Each slice's
            # `slice_len / length` is shared EQUALLY between the objects that bound it, so every bounded
            # slice disposes of exactly its own bases::
            #
            #     Sum over the fragment  =  Sum slice_len / length  =  1
            #
            # ⭐⭐ **A SJ IS A BOUNDARY EXACTLY LIKE A BOUNDARY**, and that is the whole rule. A block's
            # interior boundaries are the boundaries it crosses; its two ENDS are boundaries too whenever the
            # intron there resolved to an annotated sj. So a fragment's 1.0 is shared across every
            # object it crosses — boundaries and sj together — rather than boundaries first and sj
            # only with what is left over.
            #
            # ⛔ The predecessor gave a boundary-crossing block's bases entirely to boundaries, so a sj whose
            # two flanking blocks both crossed a boundary received NOTHING while `sj_count` credited it. That
            # still conserved — the total was 1.0 — but it is not a sharing, and it left `q_sj = 0` on
            # 35 of 8,436 crossed sj on ladder g50 capture_off.
            #
            # ⭐ Coverage-weighted, NOT `1/K`. Both conserve; only this one says WHERE the fragment sat,
            # and only this one is expressible per base — which is how the two are told apart at all.
            # ⚠ An UNSPLICED path has no sj boundaries, so this reduces to the previous rule
            # exactly and `boundary_unspliced_mass` is byte-identical. A single block with no boundary and no
            # annotated sj is bounded by nothing and deposits nothing: for a one-block path that is
            # the CONTAINED case, already whole in `region_contained_count`; for a multi-block one it is an
            # unannotated intron's block, whose bases have nowhere conserved to go.
            left_sj = sj_id_at_gap[block - 1] if block > 0 else -1
            right_sj = sj_id_at_gap[block] if block < len(sj_id_at_gap) else -1
            n_slices = last - first + 1
            for i in range(n_slices):
                lo = seg_start if i == 0 else int(region_bounds[first + i - 1])
                hi = seg_end if i == n_slices - 1 else int(region_bounds[first + i])
                left_boundary = first + i - 1 if i > 0 else -1
                right_boundary = first + i if i < n_slices - 1 else -1
                # The block's own ends are sj boundaries; its interior ends are boundaries. A slice
                # therefore has at most one boundary of each kind on each side, never both.
                left_jid = left_sj if left_boundary < 0 else -1
                right_jid = right_sj if right_boundary < 0 else -1
                n_bounds = (
                    int(left_boundary >= 0)
                    + int(right_boundary >= 0)
                    + int(left_jid >= 0)
                    + int(right_jid >= 0)
                )
                if n_bounds == 0:
                    continue
                share = (hi - lo) / (length * n_bounds)
                for boundary in (left_boundary, right_boundary):
                    if boundary >= 0:
                        boundary_mass[boundary_base + boundary - 1] += share
                for jid in (left_jid, right_jid):
                    if jid >= 0:
                        # ⭐ `column` — the SAME column `sj_count` is deposited at below, so
                        # `mass[c]/count[c]` is a per-strand mean and not a cross-population ratio.
                        t.sj_mass[jid, column] += share
            if last > first:
                sole_boundary = first if (n_crossed == 0 and last - first == 1) else -1
                n_crossed += last - first

        for jid in sj_ids:
            t.sj_count[jid, column] += 1
            # ⚠ `inv_length_sum` only. `sj_length_sum` was removed: nothing read it, and `pool_lengths`
            # already carries the spliced population's length distribution (RNA_SPLICED) unconditioned.
            t.sj_inv_length_sum[jid] += inv_boundary

        # ── contained: the WHOLE path lies inside ONE region ────────────────────────────────────────
        # ⚠ Not merely "crossed no boundary". An unannotated intron can swallow every boundary between two
        # blocks, leaving a fragment that crosses nothing yet straddles two regions — crediting that as
        # contained would place its whole length in a region it only partly overlaps. Such a fragment is
        # neither contained nor crossing: it deposits on no object but is still counted (start_count),
        # so the loss is visible rather than silent.
        first_region = self._local_region(region_bounds, first_base)
        contained_region = -1
        if not sj_ids and first_region == self._local_region(region_bounds, last_base):
            contained_region = region_base + first_region
            # ⭐⭐ THE RECIPROCAL-OPPORTUNITY DEPOSIT. A length-`w` fragment contained in a region of
            # length `ell` had `ell − w + 1` admissible start positions, so `1/(ell − w + 1)` cancels
            # the opportunity ON ITS OWN SUPPORT: `E[SUM] = rho * P(w <= ell)`, NOT `rho` — a fragment
            # with `w > ell` deposits NOTHING here, and `P(w <= ell)` is a per-component pmf functional
            # (TRAPS: a-cancellation-is-conditional-on-its-support). ⛔ The predecessor deposited
            # `1/L`, which does not cancel `(ell − w + 1)` at all — measured, that channel read 25.67
            # density units for short fragments and 1.60 for long ones at the same true density.
            # ⛔ The BOUNDARY rule `1/(L−1)` is NOT this rule's `ell -> 0` limit (that limit is 0 for
            # every `w >= 2`); it is a DIFFERENT relation — crossing a designated point, `A = w − 1` at
            # every `ell` — whose support factor `P(w >= 2)` is 1 for any real library.
            # ⚠ `A >= 1` here is the support restated: the fragment IS contained, so `w <= ell` — which
            # is exactly why `E[SUM] = rho * P(w <= ell)` and not `rho`.
            region_len = int(region_bounds[first_region + 1]) - int(region_bounds[first_region])
            t.region_contained_count[contained_region, column] += 1
            t.region_contained_inv_opportunity_sum[contained_region] += 1.0 / (
                region_len - length + 1
            )

        pool = self._pool(spliced, contained_region, sole_boundary, region_base)
        if pool is not None:
            t.pool_lengths[pool, length] += 1
        return DepositOutcome.DEPOSITED

    # -- the drain --------------------------------------------------------------------------------

    def drain(self, choices) -> dict[str, int]:
        """⭐ **THE SECOND PASS'S DRAIN.** Replay each held fragment with ONE chosen hypothesis.

        ``choices[i]`` is an index into hypothesis set of the ``i``-th record of
        :meth:`Tally.deferred_canonical` — the queue's one canonical order, which is why that order is
        defined in exactly one place.

        ⭐ **There is no second tally path.** Each record re-enters :meth:`deposit` with its chosen
        hypothesis **alone**: a set of size one, so the arbitration is degenerate and the fragment either
        deposits or is rejected by the ordinary rules. Every crossing rule, the quantum, the pools and
        ``L`` itself are reached through the same code that ran in pass one, so byte-identity with the
        native accumulator is preserved for free rather than argued.

        ⛔ **The drain's outcomes do NOT go into ``gap_resolution``, and that is structural rather than
        stylistic.** The census has no ``gap_resolved_unspliced`` class because pass-one arbitration can
        never produce one — a spliced path always region_bounds bases the genomic path keeps, so the genomic path
        can never be the sole survivor. ⭐ But the DRAIN can *choose* it. Folding drain outcomes into that
        census would require resurrecting the class S1 deleted for a proven reason.

        ⚠ **And left alone means RESTORED, not merely not-added-to.** ``deposit`` reaches
        ``_record_gap_resolution`` on every fragment, and for a set of size one that lands on
        ``RESOLVED_SPLICED`` — while a chosen ∅ hits the ``all(h.is_unspliced)`` early return and is not
        counted at all. So a naive drain would grow the census by however many draws happened to pick a
        spliced path: a census that depends on the RNG. This method snapshots it and puts it back, and the
        information the census cannot express lives on the drain's own axis as
        ``chose_genomic`` / ``chose_spliced``.

        ⭐ **After the drain the tally describes the FINAL state**: nothing is held, so the bank is empty
        *and* ``deferred_undetermined_gap`` is 0 and the three ``gap_deferred_*`` are 0. That keeps
        "the counter and the fragments it counts are the same population" absolute — the payload refuses a
        bank that disagrees with its counter at the door, and a drained payload must not need an exception.
        Pass one's numbers are returned here, and the payload keeps them in its ``drain`` block.

        Returns the drain's own counters plus ``census_before``. The conservation statement is
        ``deposited + dropped_* == offered``.
        """
        held = self.tally.deferred_canonical()
        if len(choices) != len(held):
            raise ValueError(
                f"drain got {len(choices)} choices for {len(held)} held fragments. One choice per held "
                f"record, in `deferred_canonical` order — a length mismatch means the producer of the "
                f"scores and the consumer of the draw walked different queues."
            )
        counters = {
            "offered": len(held),
            "deposited": 0,
            "dropped_too_long": 0,
            "dropped_empty": 0,
            "dropped_strand_undefined": 0,
            "chose_genomic": 0,
            "chose_spliced": 0,
            "census_before": dict(self.tally.gap_resolution),
        }
        census_before = dict(self.tally.gap_resolution)
        # ⚠ Emptied FIRST. `deposit` appends to this list when it defers, and a set of size one cannot
        # defer — but draining into a live queue would make that assumption invisible instead of checked,
        # and the assertion below is what checks it.
        self.tally.deferred = []
        for fragment, choice in zip(held, choices):
            path = fragment.hypotheses[int(choice)]
            counters["chose_genomic" if path.is_unspliced else "chose_spliced"] += 1
            outcome = self.deposit(
                fragment.ref,
                fragment.start,
                fragment.end,
                observed_introns=fragment.observed_introns,
                align_strand=fragment.align_strand,
                sj_strand=fragment.sj_strand,
                hypotheses=(path,),
            )
            if outcome is DepositOutcome.DEFERRED:
                raise AssertionError(
                    "a hypothesis set of size one deferred, which is unreachable: arbitration defers "
                    "only when two or more hypotheses survive."
                )
            counters[outcome.value] += 1
        if (
            counters["deposited"] + sum(counters[k] for k in counters if k.startswith("dropped_"))
            != counters["offered"]
        ):
            raise AssertionError(
                "the drain lost a fragment; deposited + dropped_* must equal offered"
            )

        # ⭐ The census restored to pass one's, and the held population now zero everywhere it is counted.
        self.tally.gap_resolution.update(census_before)
        for key in (
            GapResolution.DEFERRED_RNA_OR_GDNA,
            GapResolution.DEFERRED_WHICH_INTRONS,
            GapResolution.DEFERRED_BOTH,
        ):
            self.tally.gap_resolution[key.value] = 0
        self.tally.qc[DepositOutcome.DEFERRED.value] = 0
        return counters

    # -- helpers ----------------------------------------------------------------------------------

    def _reject(self, outcome: DepositOutcome) -> DepositOutcome:
        self.tally.qc[outcome.value] += 1
        return outcome

    @staticmethod
    def _hypothesis_length(start, end, observed, path: GapHypothesis):
        """``(L, normalised introns, absorbed)`` for one hypothesis.

        ⭐ ONE definition of ``L``, and one code path to it, whether the hypothesis wins or is only being
        scored for the filter. The observed introns and the implied ones are UNIONED — they are disjoint
        by construction (the enumeration never searches a gap the CIGAR already explained) so this is a
        union and not a merge, and ``_normalise_introns`` sorts and clips it either way.
        """
        introns, absorbed = _normalise_introns(tuple(observed) + path.introns, start, end)
        return sum(b - a for a, b in _segments(start, end, introns)), introns, absorbed

    def _record_gap_resolution(self, hypotheses, survivors) -> None:
        """⭐ The umbrella census — the ONE site, so the subclasses cannot drift out of sync.

        Counted only when the enumeration had something to arbitrate: a fragment whose only hypothesis
        is the unspliced one never had a question to answer, and counting it would quietly make the
        umbrella's denominator the whole library.
        """
        if all(h.is_unspliced for h in hypotheses):
            return
        n_spliced = sum(1 for h, *_ in survivors if not h.is_unspliced)
        unspliced_survives = any(h.is_unspliced for h, *_ in survivors)
        if len(survivors) == 1:
            # ⛔ The sole survivor is necessarily the SPLICED one, and there is no branch for the other
            # case because it cannot occur: the genomic path is always the longest, so it can never be the
            # last one left. See :class:`GapResolution`.
            resolution = GapResolution.RESOLVED_SPLICED
        elif not unspliced_survives:
            resolution = (
                GapResolution.DEFERRED_WHICH_INTRONS
            )  # certified RNA; only the structure is open
        elif n_spliced == 1:
            resolution = GapResolution.DEFERRED_RNA_OR_GDNA  # one bit: was anything spliced at all?
        else:
            resolution = GapResolution.DEFERRED_BOTH
        self.tally.gap_resolution[resolution.value] += 1

    @staticmethod
    def _local_region(region_bounds: np.ndarray, position: int) -> int:
        """Index within this reference of the region containing ``position``."""
        return min(max(int(np.searchsorted(region_bounds, position, side="right")) - 1, 0), region_bounds.size - 2)

    def _sj_edge_id(self, ref: int, intron_start: int, intron_end: int, sj_strand: int) -> int:
        """The sj-boundary id for this intron, or -1 if it is not annotated.

        One to three iterations at human scale: 70.4 % of region_bounds are not a donor at all, and over those
        that are, the mean fan-out is 1.31.
        """
        p = self.partition
        donor = _exact_region_bound(p.region_bounds, p.ref_region_bound_offsets, ref, intron_start)
        if donor < 0:
            return -1
        acceptor = _exact_region_bound(p.region_bounds, p.ref_region_bound_offsets, ref, intron_end)
        if acceptor < 0:
            return -1
        for k in range(int(p.sj_offsets[donor]), int(p.sj_offsets[donor + 1])):
            if int(p.sj_boundary_right[k]) != acceptor:
                continue
            # ⚠ The strand filter applies only when the motif strand is DEFINITE. AMBIGUOUS means "no
            # strand information", not "a strand that matches nothing" — it used to be treated as the
            # latter, which silently demoted an annotated spliced fragment to an unspliced deposit,
            # credited `unannotated_introns` (poisoning the one metric that measures annotation coverage),
            # and dropped it from the pure RNA pool the length model is fitted from.
            #
            # The filter exists only to separate two sj sharing a coordinate pair. ⭐ Measured: 0 of
            # 404,168 human sj coordinates are annotated on both strands, so it can only ever lose a
            # match, never disambiguate one. The resolver already reads a non-definite strand as "either"
            # (`sj_lookup_into` falls back to the union of POS and NEG), so this agrees with it.
            if sj_strand in STRAND_COLUMNS and int(p.sj_strand[k]) != int(sj_strand):
                continue
            return k
        return -1

    def _pool(self, spliced, contained_region, sole_boundary, region_base):
        """The one length pool this fragment belongs to, or ``None``.

        Priority, so that every pool stays pure: a splice is unambiguously RNA; a contained fragment is
        typed by its region; a single-boundary crossing is a "splash" read typed by its two flank types.
        Anything else — an exonic contained fragment, a multi-boundary crossing — is a mixture and enters
        nothing.

        ⭐ **DETERMINACY, NOT PROVENANCE**. There used to be an
        ``sj_implicit`` condition here barring a fragment whose splice was inferred rather than
        sequenced, on the grounds that a length partly inferred from the annotation is a product of the
        model the pool is used to fit. It is gone, and so is the flag: a fragment only reaches this boundary
        when **exactly one hypothesis survived**, so its ``L`` is not in doubt at all.

        ⚠ Measured before deleting it, because the two criteria disagree and the disagreement is large:
        on the chr22 pilot the pool reads **+0.67 % mean / +2.40 % sd** against truth under determinacy
        and **−9.58 % / −22.46 %** under provenance. Excluding inferred lengths preferentially excludes
        fragments whose mates sit far apart — **a purity filter on a length pool is a length filter**.
        """
        types = self.partition.region_types
        if spliced:
            return FragmentPool.RNA_SPLICED
        if contained_region >= 0:
            return _CONTAINED_POOL.get(int(types[contained_region]))
        if sole_boundary >= 0:
            flanks = (int(types[region_base + sole_boundary - 1]), int(types[region_base + sole_boundary]))
            return _SPLASH_POOL.get(tuple(sorted(flanks)))
        return None
