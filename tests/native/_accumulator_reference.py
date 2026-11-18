"""The rigel accumulator: the executable specification.

The reference implementation of the per-fragment tally built during the single-pass BAM scan. The native
accumulator must reproduce it byte for byte. It defines the deposit rule end to end (the region, boundary and
sj axes, the hypothesis arbitration over a fragment's unsequenced gaps, the deferred queue and its canonical
order, the conserved mass, the five structural length pools and the QC ledger) and the drain that replays a
held fragment with one chosen hypothesis. Tests: ``test_accumulator_spec.py`` (the deposit rule),
``test_accumulator_native_parity.py`` (the C++ against this file), ``test_gap_hypothesis_arbitration.py``,
``test_accumulator_drain.py`` and ``test_conserved_mass.py``.

THE MODEL
    The genome is a graph. Regions are half-open intervals tiling each reference, numbered in genomic
    order. The 0-bp gaps between adjacent regions are contiguous boundaries; a sj boundary is a
    directed donor→acceptor link taken from the annotation. A fragment is a path: its aligned blocks,
    joined across mate gaps and broken by introns::

        ref   |······ n0 ······|·· n1 ··|········ n2 ········|
        region_bounds  0               100      200                  600
        boundaries                  1        2
        path        [====== block ======]                       crosses boundary 1
        path        [= block =]~~intron~~[==== block ====]       crosses nothing; uses a sj

    Regions count fragments contained (they fit inside); boundaries count fragments crossing. Each
    population stores only the channels something reads: the two banks the deconvolution consumes carry
    ``count`` = Sum 1 (integer), and the boundary one also a reciprocal-opportunity sum ``Sum 1/A(w)``
    (float64); the certified-RNA banks carry fewer, because nothing deconvolves a fragment already known
    to be RNA. No spliced fragment touches the region axis, since both endpoints of an annotated intron
    are region bounds and a spliced fragment can never be contained.

WHAT EACH OBJECT'S NUMBERS MEAN
    With ``A(w)`` the number of admissible start positions (``w − 1`` crossing a 0-bp boundary) and a
    component of start density ``rho`` and length distribution ``f``, the deposit is ``1/A(w)`` and::

        E[count]      =  rho * E_f[A(w)]
        E[sum 1/A]    =  rho * P_f(A > 0)          <- the cancellation holds only on its support
                      =  rho                        at a boundary: P(w >= 2) = 1 for any real library

NO PARTITIONING
    Every crossed boundary receives the full weight. The chance that a length-``L`` fragment crosses a
    given boundary is proportional to ``L − 1`` and the deposit is ``1/(L − 1)``, so the two cancel and
    every fragment length contributes equally to each boundary. Dividing by the number of boundaries
    crossed would break that cancellation and make the answer depend on region spacing. A boundary count
    is not a share of a conserved total: density is intensive, not extensive.
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
    "Placement",
    "DepositOutcome",
    "FragmentPool",
    "GapResolution",
    "GapHypothesis",
    "Partition",
    "Tally",
]

#: The strand convention, the same everywhere in this file.
#:
#: Every count and density array is ``[n, N_STRAND_COLUMNS]``, and the two columns are the two genome strands,
#: ``Strand.POS`` and ``Strand.NEG``, following the index's ``_pos``/``_neg`` column naming (``reach_lo_pos``,
#: ``reach_lo_neg``, …). Nothing is stored transcript-relative. Sense/antisense is derived, never stored: a sj boundary
#: carries its own genomic strand, so a consumer computes ``sense = (fragment strand == sj strand)``.
#:
#: ``Strand`` has four values: OR semantics make ``POS | NEG == AMBIGUOUS``, and ``NONE`` means no strand. Only POS and
#: NEG name a column, so the mapping is a dict and any other strand raises ``KeyError`` rather than being filed under
#: the minus column.
#:
#: The column is indexed by ``align_strand`` (where the read aligned), never by ``sj_strand`` (which way an intron was
#: spliced, from its ``GT..AG`` motif). Only a spliced read has the second, and it never selects a column.
STRAND_COLUMNS: dict[int, int] = {Strand.POS: 0, Strand.NEG: 1}
N_STRAND_COLUMNS = len(STRAND_COLUMNS)


#: One numeric convention: a count is an integer, a fraction is float64. There is no fixed point anywhere in the tally.
#: Float addition is not associative, so fractions are reproducible across worker counts only within a derived
#: tolerance; in exchange float64 is accurate where a fixed point in ``round(2^32 / placements)`` would miss the
#: reciprocal-opportunity theorem's exact answer by orders of magnitude more.


class DepositOutcome(enum.Enum):
    """Why a fragment did or did not deposit. Every rejection is counted, never silent."""

    DEPOSITED = "deposited"
    TOO_LONG = "dropped_too_long"
    EMPTY = "dropped_empty"
    #: The fragment has no single genome strand: NONE, or AMBIGUOUS (``POS | NEG``, which ``build_fragment`` produces
    #: when the mates agree in reference orientation). The column axis is the genome strand, so there is no column to
    #: credit; strand-undefined fragments are one of the denominators the accumulator emits.
    STRAND_UNDEFINED = "dropped_strand_undefined"
    #: Two or more hypotheses survived, so the path is not determined, and neither is ``L``, its ``1/A(w)`` deposit, the
    #: pool bin, or the set of boundaries crossed. The fragment deposits on nothing and goes to the deferred queue,
    #: where the second pass resolves it with the fragment-length distribution and the pass-one tally. It is retained in
    #: full, not dropped: the conservation identity is ``deposited + deferred + dropped_* == offered``.
    DEFERRED = "deferred_undetermined_gap"
    #: Two or more placements of one multimapping fragment admit a surviving (placement, hypothesis), so which placement
    #: the molecule took is undetermined. Held whole, every placement and every hypothesis, for the second pass, which
    #: weighs each placement by the unique traffic at its own objects. The identity is ``deposited +
    #: deferred_undetermined_gap + deferred_multiple_placements + dropped_* == offered``.
    DEFERRED_PLACEMENTS = "deferred_multiple_placements"


@dataclass(frozen=True, slots=True)
class GapHypothesis:
    """One hypothesis about what a fragment's unsequenced gaps contain.

    A fragment's mate gap may hold no intron, one, or several, and which cannot be observed. Each candidate transcript
    determines one answer (its own introns inside the gaps), so the hypotheses are few, and two transcripts implying the
    same introns are one hypothesis.

    The empty path is the genomic hypothesis: cutting nothing means the gap is real template, so the molecule is gDNA,
    or RNA unspliced there, the same span. So the accumulator needs no separate "could this be gDNA?" flag, and a
    synthetic span is not a separate candidate: it is this hypothesis, always in the set unless something rules it out.

    ``introns`` are the implied ones only. Introns the CIGAR stated are cut under every hypothesis and are passed to
    :meth:`Accumulator.deposit` separately.
    """

    #: The implied introns as ``(start, end)`` pairs. Empty means the genomic hypothesis.
    introns: tuple[tuple[int, int], ...] = ()
    #: The strand of the transcripts that imply this path. An inference, used only when no motif was sequenced anywhere
    #: on the fragment; an observed motif always wins.
    sj_strand: int = Strand.NONE

    @property
    def is_unspliced(self) -> bool:
        """True for the empty path: it cuts nothing, so the gap is template."""
        return not self.introns


#: The one hypothesis every fragment has: cut nothing beyond what was sequenced. A fragment with no unsequenced gap, or
#: no annotated intron in it, has exactly this and deposits without arbitration.
UNSPLICED_ONLY: tuple[GapHypothesis, ...] = (GapHypothesis(),)


class GapResolution(enum.Enum):
    """The arbitration census: an umbrella class and its subclasses.

    Every fragment for which the enumeration produced at least one non-genomic hypothesis is counted here, since its
    ``L`` depends on whether a gap intron is cut. The subclasses are exhaustive and mutually exclusive, so the counts
    close::

        Sum(GapResolution) == the umbrella total
        the three DEFERRED_* == qc["deferred_undetermined_gap"]

    This is its own axis, not a ``splice_type``: the umbrella cuts across the splice census (a certified-RNA
    ``SPLICED_ANNOT`` fragment with an intron in its mate gap needs resolving as much as an ``UNSPLICED`` one), and the
    splice census must keep summing to the library. The classes describe the arbitration, not the deposit; a
    ``RESOLVED_*`` fragment can still be rejected afterwards as ``TOO_LONG``.

    There is no ``RESOLVED_UNSPLICED`` because no fragment could enter it. A spliced hypothesis cuts bases the genomic
    one keeps, so ``L_spliced <= L_genomic``; the one filter is ``L <= max_fragment_length``; so if the genomic path
    survives every spliced path survives too, and the survivors are never exactly ``{genomic}`` while a spliced path was
    offered.
    """

    #: One hypothesis survived. It cuts something (see the class docstring for why it cannot be the genomic one), so the
    #: gap intron is real and ``L`` excludes it.
    RESOLVED_SPLICED = "gap_resolved_spliced"
    #: The genomic path against exactly one spliced path. The open question is RNA or gDNA, the composition question
    #: calibration answers.
    DEFERRED_RNA_OR_GDNA = "gap_deferred_rna_or_gdna"
    #: Two or more spliced paths and no genomic one: the molecule is certified RNA (gDNA cannot be
    #: spliced) and the open question is purely which structure.
    DEFERRED_WHICH_INTRONS = "gap_deferred_which_introns"
    #: Both questions at once.
    DEFERRED_BOTH = "gap_deferred_both"


@dataclass(frozen=True, slots=True)
class Placement:
    """One alignment of a fragment on one reference, with every explanation of its unsequenced gaps.

    A uniquely mapped fragment has one; a multimapping fragment has one per resolved alignment, each with its own
    reference, extent, strands, observed introns and gap hypotheses. This is what the scanner offers
    (:meth:`Accumulator.offer`) and what a held record stores, clipped to its reference.
    """

    ref: int
    start: int
    end: int
    align_strand: int
    sj_strand: int
    observed_introns: tuple[tuple[int, int], ...] = ()
    hypotheses: tuple[GapHypothesis, ...] = UNSPLICED_ONLY
    #: The hit's sequence over both mates: its edits and its soft-clipped read bases. Never deposited; the second
    #: pass's draw prices a multimapper's placements by it.
    nm: int = 0
    clip: int = 0

    @property
    def key(self):
        """The placement's own content, in the canonical order's comparison form (Python tuple order,
        prefix rule included). Two placements equal here deposit identically, whatever their sequence."""
        return (
            self.ref,
            self.start,
            self.end,
            self.align_strand,
            self.sj_strand,
            self.observed_introns,
            tuple((p.introns, p.sj_strand) for p in self.hypotheses),
        )


@dataclass(frozen=True, slots=True)
class DeferredFragment:
    """A fragment held for the second pass, stored whole: every placement, every hypothesis.

    The fragment is stored, never its consequences: object ids are large, derived and tied to the partition, while the
    fragment is small and replays exactly. The drain re-enters :meth:`Accumulator.deposit` at the chosen placement with
    the chosen hypothesis, so there is no second deposit path and byte-identity with the native accumulator follows. A
    gap-held unique mapper is a record with one placement; a multimapper holds several, in their own canonical order, so
    the record does not depend on the order the aligner listed the hits in.
    """

    placements: tuple[Placement, ...]

    @property
    def n_placements(self) -> int:
        return len(self.placements)

    @property
    def key(self):
        """The record's content first, then its placements' sequences: a record's sequence must never move it
        past a record its content orders, and records whose content ties are ordered by it, so the order is
        total."""
        return tuple(p.key for p in self.placements), tuple((p.nm, p.clip) for p in self.placements)


class FragmentPool(enum.IntEnum):
    """Fragment-length pools, each defined by structure.

    RNA_SPLICED is certified RNA, since gDNA cannot splice. The four DNA_* pools target gDNA but are not pure (RNA
    inside introns and unannotated transcription reach them too), and ``rigel.calibration.fl`` replaces purity with a
    two-pool contrast. There is no pool for an exonic contained fragment or a multi-boundary crossing, which are
    gDNA/RNA mixtures by structure.
    """

    DNA_INTERGENIC = 0  # contained in an intergenic region
    DNA_INTRONIC = 1  # contained in an intronic region
    DNA_INTRON_EXON = 2  # crossing exactly one boundary, flanks {intron, exon}
    DNA_INTERGENIC_EXON = 3  # crossing exactly one boundary, flanks {intergenic, exon}
    RNA_SPLICED = 4  # using an annotated sj on a determined path, sequenced or implied


#: Coarse region types, as ``signature.coarse_type_array`` emits them.
_TYPE_INTERGENIC, _TYPE_INTRON, _TYPE_EXON = 0, 1, 2

#: Flank-type pair -> splash pool. The key is sorted, so it is order-insensitive.
_SPLASH_POOL = {
    (_TYPE_INTRON, _TYPE_EXON): FragmentPool.DNA_INTRON_EXON,
    (_TYPE_INTERGENIC, _TYPE_EXON): FragmentPool.DNA_INTERGENIC_EXON,
}

#: Contained region type -> its gDNA pool. An exonic region is a mixture and is absent by design.
_CONTAINED_POOL = {
    _TYPE_INTERGENIC: FragmentPool.DNA_INTERGENIC,
    _TYPE_INTRON: FragmentPool.DNA_INTRONIC,
}


# ---------------------------------------------------------------------------
# the partition: the region-bound axis, the region types, and the sj CSR
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Partition:
    """Everything the deposit addresses: three axes, off by one from each other per reference.

    A reference contributing ``c`` region-bound positions owns ``c − 1`` regions and ``c − 2`` contiguous (interior)
    boundaries. A reference with no regions contributes no region bounds::

        region_bounds    0        100       200       600        c = 4
        regions   [  n0  ][   n1   ][   n2   ]            c - 1 = 3
        boundaries            boundary 1    boundary 2               c - 2 = 2

    Splice junctions are a CSR keyed by the donor region-bound index, which the deposit has already computed while
    locating the boundaries its path crosses. Every annotated intron has both endpoints as region bounds, so "is this
    intron annotated?" reduces to "are both endpoints region bounds, and is the pair registered?", and the table is not
    consulted when the start is not a bound.
    """

    region_bounds: (
        np.ndarray
    )  # int64[n_region_bounds]: flat, reference-major, ascending within a reference
    ref_region_bound_offsets: np.ndarray  # int64[n_refs + 1]: CSR over region_bounds
    region_types: np.ndarray  # uint8[n_regions]: 0 intergenic / 1 intron / 2 exon
    ref_region_offsets: np.ndarray  # int64[n_refs + 1]
    ref_boundary_offsets: np.ndarray  # int64[n_refs + 1]
    sj_offsets: np.ndarray  # int64[n_region_bounds + 1]: CSR over the donor region_bound index
    sj_boundary_right: np.ndarray  # int64[n_sj]: flat region_bound index of the intron's high end
    sj_strand: np.ndarray  # int8[n_sj]

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
        """Build from per-reference region-bound lists.

        ``sj`` are ``(ref, intron_start, intron_end, sj_strand)``; both endpoints must be region bounds
        on that reference. Splice-junction ids are assigned by sorting on
        ``(donor bound, acceptor bound, sj_strand)``, so they are a deterministic function of the
        partition alone.
        """
        region_bounds_per_ref = [np.asarray(c, dtype=np.int64) for c in region_bounds_per_ref]
        n_refs = len(region_bounds_per_ref)
        for r, c in enumerate(region_bounds_per_ref):
            if c.size and (c.size < 2 or np.any(np.diff(c) <= 0)):
                raise ValueError(
                    f"reference {r}: region_bounds must strictly increase, got {c.tolist()}"
                )

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
        boundary_left, boundary_right, sj_strand = (
            boundary_left[order],
            boundary_right[order],
            sj_strand[order],
        )

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
    """The flat region-bound index of ``position`` on ``ref``, or -1 if it is not a bound there."""
    first, last = int(region_bound_offsets[ref]), int(region_bound_offsets[ref + 1])
    if last <= first:
        return -1
    k = first + int(np.searchsorted(region_bounds[first:last], position))
    return k if k < last and int(region_bounds[k]) == position else -1


def _normalise_introns(introns, start: int, end: int) -> list[tuple[int, int]]:
    """The introns as a sorted, disjoint set clipped to ``[start, end)``.

    A real BAM does not always give a clean set: the splice-motif strand tag is read once per record (``XS`` for STAR,
    ``ts`` for minimap2, auto-detected), so a pair whose mates disagree about an intron's acceptor yields two
    overlapping introns for one molecule. Normalising here lets ``L`` be defined as the total length of the path's
    segments, with no second formula that could disagree (``span − Σ introns`` goes wrong on overlaps and negative on a
    wide one).
    """
    merged: list[tuple[int, int]] = []
    for raw_start, raw_end in sorted((int(s), int(e)) for s, e in introns):
        intron_start, intron_end = max(raw_start, start), min(raw_end, end)
        if intron_end <= intron_start:  # zero-length, or entirely outside the fragment
            continue
        # Overlapping and abutting introns are both malformed and both merge. Abutting introns imply a zero-length exon
        # between them, which is physically impossible (a transcript with one is molecularly identical to one without),
        # so a molecule cannot use both. The index cannot produce the pair either: a zero-length exon is dropped when
        # the exon arrays are built, fusing its flanking introns. So an abutting pair is an alignment artifact.
        if merged and intron_start <= merged[-1][1]:
            previous_start, previous_end = merged[-1]
            merged[-1] = (previous_start, max(previous_end, intron_end))
            continue
        merged.append((intron_start, intron_end))
    return merged


def _segments(start: int, end: int, introns) -> list[tuple[int, int]]:
    """The path's contiguous genomic segments: ``[start, end)`` with the introns cut out.

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
    """The accumulator's output: per-object sums over the fragments that landed on each object.

    ``count`` = ``Sum 1`` carries the statistical power (a Beta-Binomial needs an integer); the reciprocal-opportunity
    sum ``inv_length_sum`` = ``Sum 1/(w−1)`` carries the level (at a boundary or sj, an exact model-free density, ``E =
    rho·P(w≥2) = rho``); ``mass`` is the conserved fragment count and sums to one per fragment.

    A channel is stored where something reads it and nowhere else::

        region_contained    count                                the mixture, on the region axis
        boundary_unspliced    count  inv_length_sum        mass   the mixture, on the boundary axis
        boundary_spliced      count                        mass   certified RNA; nothing deconvolves it
        sj          count  inv_length_sum               inv_length_sum is read by second_pass

    The count banks are integers, so their per-worker merge is bit-identical at any thread count; the reciprocal and
    mass banks are float64, reproducible within a derived tolerance.
    """

    region_contained_count: np.ndarray  # uint32[n_regions, 2]
    #: uint32[n_regions, 2]: the path's first covered base, by align strand; one per accepted fragment, so ``Σ (both
    #: columns) == qc.deposited``, the ledger invariant. Opportunity ``ℓ`` for every fragment length; wall-blind only at
    #: the template's downstream end.
    region_start_count: np.ndarray
    boundary_unspliced_count: np.ndarray  # uint32[n_boundaries, 2]
    boundary_unspliced_inv_length_sum: np.ndarray  # float64[n_boundaries]: one column
    #: float64[n_boundaries]: the conserved mass (see :meth:`Accumulator.deposit`). A count and a mass are different
    #: deposits: ``count`` is discrete and ``+1`` on every boundary a fragment crosses, so a fragment contributes
    #: ``max(K, 1)`` of them, while the mass sums to one per fragment.
    #:
    #: One column: the other boundary banks keep both genome strands because the strand fit's boundary seeds read
    #: ``boundary_unspliced.count`` per column, but nothing reads a mass per strand; the mass converts an
    #: object-incidence total into a fragment count, a question with no strand in it.
    boundary_unspliced_mass: np.ndarray
    #: uint32[n_boundaries, 2]: certified RNA, since gDNA cannot be spliced. Count and mass only: nothing deconvolves a
    #: certified-RNA crossing, and `pool_lengths`' RNA_SPLICED row already carries that population's lengths.
    boundary_spliced_count: np.ndarray
    #: float64[n_boundaries]: the same rule, routed by the same ``spliced`` flag. A partial by construction, not a
    #: conservation ledger: a per-boundary certified-RNA term, not "the number of spliced fragments at this boundary".
    #: It keeps a spliced fragment's mass out of the unspliced bank, where certified RNA would enter the competition the
    #: prior arbitrates.
    boundary_spliced_mass: np.ndarray
    #: uint32[n_sj, 2]: both genome-strand columns are kept for aligner-artifact detection. A sj is stranded by its
    #: splicing motif, so the strand of the fragments on it looks redundant, and every consumer sums the two.
    #:
    #: Aligners emit false-positive ``N`` CIGAR ops from plain genomic DNA. In a stranded library a real sj inherits
    #: the global strand specificity while an artifact deposits on both strands, so this split is an empirical
    #: detector; a clean sj and an artifactual one carry the same total. Unstranded data cannot use it (κ = ½ leaves
    #: nothing to deviate from).
    sj_count: np.ndarray
    #: float64[n_sj]: ``second_pass.py`` scores a held fragment's sj evidence with it.
    sj_inv_length_sum: np.ndarray
    #: float64[n_sj, 2]: the conserved mass's third axis, which makes a library fragment count computable. A spliced
    #: fragment's block containing no interior boundary deposits nothing on either boundary bank and is not
    #: ``contained`` (its path spans a sj). Without this bank such a fragment exists on the incidence axis
    #: (``sj_count``) and on no conserved one, and that population is all RNA, since gDNA cannot splice.
    #:
    #: A sj bounds a block exactly as a contiguous boundary does (:meth:`Accumulator.deposit`), so a slice bounded by
    #: both shares its bases between them. A sj counts only where the intron resolved to an annotated sj; a block
    #: bounded solely by unannotated introns has nowhere conserved to send its bases, so the identity is exact over
    #: deposited, annotated fragments.
    #:
    #: Two columns, the only mass with any: an artifactual sj accumulates symmetrically on both strands like gDNA, so a
    #: per-strand mass lets artifact filtering run in a single pass. The columns are `sj_count`'s columns.
    sj_mass: np.ndarray
    pool_lengths: np.ndarray  # int64[5, max_fragment_length + 1]: binned at L, once per fragment
    #: uint32[max_fragment_length + 1]: every deposited fragment, binned at its own L, with no purity condition. The
    #: five pools above are conditioned on structure, so they cannot serve as the unconditional anchor an
    #: empirical-Bayes shrinkage needs; this row is that anchor, measured over the same population by the same rule.
    #:
    #: It is unconditional given deposit: it excludes what the accumulator rejects (over the length limit, ambiguous
    #: path, strand-undefined, empty), each counted in ``qc``, which is exactly the population the pools are drawn from.
    #: Its invariant ``Σ deposited_lengths == Σ region_start_count == qc.deposited`` says every fragment was binned by
    #: length, as the start count says every fragment was located in space.
    deposited_lengths: np.ndarray
    qc: dict[str, int] = field(default_factory=dict)
    #: The deferred queue: fragments whose path is not determined, held whole for the second pass. Not an array, so the
    #: parity test reaches it through :meth:`Tally.deferred_arrays`, whose flattening is specified here.
    deferred: list[DeferredFragment] = field(default_factory=list)
    #: The umbrella census: one entry per :class:`GapResolution`, incremented at one site.
    gap_resolution: dict[str, int] = field(default_factory=dict)

    @classmethod
    def zeros(cls, n_regions: int, n_boundaries: int, n_sj: int, max_length: int) -> "Tally":
        def counts(rows):
            return np.zeros((rows, N_STRAND_COLUMNS), np.uint32)

        def fraction(rows):
            # One column. The length moments are strand-agnostic: the strand a read aligned to says nothing about
            # whether the molecule was gDNA or RNA, and every consumer sums the two.
            return np.zeros(rows, np.float64)

        return cls(
            region_contained_count=counts(n_regions),
            region_start_count=counts(n_regions),
            boundary_unspliced_count=counts(n_boundaries),
            boundary_unspliced_inv_length_sum=fraction(n_boundaries),
            boundary_unspliced_mass=fraction(n_boundaries),
            boundary_spliced_count=counts(n_boundaries),
            boundary_spliced_mass=fraction(n_boundaries),
            sj_count=counts(n_sj),
            sj_inv_length_sum=fraction(n_sj),
            # Two columns, unlike every other mass here; see the field.
            sj_mass=np.zeros((n_sj, N_STRAND_COLUMNS), np.float64),
            pool_lengths=np.zeros((len(FragmentPool), max_length + 1), np.int64),
            deposited_lengths=np.zeros(max_length + 1, np.uint32),
            qc={outcome.value: 0 for outcome in DepositOutcome},
            deferred=[],
            gap_resolution={cls.value: 0 for cls in GapResolution},
        )

    def deferred_canonical(self) -> list["DeferredFragment"]:
        """The deferred queue in its one canonical order, sorted on each record's own content.

        The drain needs the same order: the second pass indexes into the queue, so a choice vector means something only
        if the producer of the scores and the consumer of the draw walk it identically. Defining the order once makes
        ``choices[i]`` refer to the same fragment in both languages. Mirrors ``Accumulator::deferred_canonical`` in the
        C++.
        """
        return sorted(self.deferred, key=lambda f: f.key)

    def deferred_arrays(self) -> dict[str, np.ndarray]:
        """The deferred queue flattened to the CSR the payload carries.

        Specified here so the readable list and the native arrays have one definition; the parity test compares this.
        Three nested variable-length levels (fragments hold placements, placements hold observed introns and hypotheses,
        hypotheses hold implied introns) give four offset arrays. Offsets are cumulative and start at 0, so
        ``n_fragments`` is ``len(placement_offsets) - 1`` and an empty queue is ``[0]``, never ``[]``. A gap-held unique
        mapper is a record with one placement.

        Sorted, which makes it bit-identical at any worker count. Every other bank is a sum, merged exactly (uint32) or
        within a derived tolerance (float64) in any chunk order; the queue is a list, so concatenated per-worker queues
        would differ in order at 1, 2, 4 and 8 workers with identical contents. Two records that tie on the sort key are
        identical, so their relative order is unobservable.
        """
        placement_fields = ("ref", "start", "end", "align_strand", "sj_strand", "nm", "clip")
        out: dict[str, list[int]] = {name: [] for name in placement_fields}
        ordered = self.deferred_canonical()
        out |= {
            "placement_offsets": [0],
            "observed_intron_offsets": [0],
            "observed_introns": [],
            "hypothesis_offsets": [0],
            "hypothesis_sj_strand": [],
            "hypothesis_intron_offsets": [0],
            "hypothesis_introns": [],
        }
        for frag in ordered:
            for placement in frag.placements:
                for name in placement_fields:
                    out[name].append(getattr(placement, name))
                for start, end in placement.observed_introns:
                    out["observed_introns"] += [start, end]
                out["observed_intron_offsets"].append(len(out["observed_introns"]) // 2)
                for path in placement.hypotheses:
                    out["hypothesis_sj_strand"].append(int(path.sj_strand))
                    for start, end in path.introns:
                        out["hypothesis_introns"] += [start, end]
                    out["hypothesis_intron_offsets"].append(len(out["hypothesis_introns"]) // 2)
                out["hypothesis_offsets"].append(len(out["hypothesis_sj_strand"]))
            out["placement_offsets"].append(len(out["ref"]))
        return {name: np.asarray(values, dtype=np.int64) for name, values in out.items()}


# ---------------------------------------------------------------------------
# the deposit
# ---------------------------------------------------------------------------


class Accumulator:
    """Deposit fragments into a :class:`Tally` against a :class:`Partition`.

    A fragment arrives as its genomic extent on one reference plus the introns inside it. Mates must already be merged
    into that extent, because the unsequenced gap between them is part of the molecule and counts toward ``L``.
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
        nm: int = 0,
        clip: int = 0,
    ) -> DepositOutcome:
        """Deposit one fragment, or hold it for the second pass; return which, and why.

        ``[start, end)`` is the full genomic extent, leftmost block start to rightmost block end, mate gap included.
        ``introns`` are the introns the CIGAR stated, as ``(start, end)`` pairs, cut under every hypothesis.

        ``hypotheses`` is the set of competing answers about the fragment's unsequenced gaps (see
        :class:`GapHypothesis`; the empty path is the genomic one). This method filters the set, then:

        * exactly one survives → deposit it; the path is determined;
        * two or more survive → :attr:`DepositOutcome.DEFERRED`, and the fragment goes whole into :attr:`Tally.deferred`
          for the second pass, which has the fragment-length distribution and the pass-one tally to discriminate with.

        The decision lives here, in one place, so the caller hands over every answer uncollapsed.

        The one filter is ``max_fragment_length``, the same statement that makes ``TOO_LONG`` a rejection: short-read
        chemistry does not sequence molecules past it, so a hypothesis implying a longer ``L`` is not a molecule this
        library contains. Applied to the genomic hypothesis it is the rule "a fragment whose genomic span exceeds the
        limit must be RNA", since the genomic path's ``L`` is that span. If the filter would empty the set, the
        survivors stand and the ordinary ``TOO_LONG`` rejection counts them. The filter changes classification, not just
        cost (hypotheses at 400 and 1200 are determined with it and ambiguous without): it is the second pass's
        likelihood applied early, where the fragment-length distribution has no mass.

        Two independent strands. Every read has the first; only a splice has the second.

        ``align_strand``
            The genomic strand the read aligned to, ``+`` or ``−``. It selects the array column and nothing else.
        ``sj_strand``
            A splice junction's strand, read from the aligner's genomic-motif tag (``XS`` for STAR, ``ts`` for minimap2, auto-detected). A ``GT..AG`` intron is on ``+``, its reverse complement ``CT..AC`` on ``−``. It is used only to resolve an intron against the annotation (see :meth:`_sj_edge_id`).

        Neither constrains the other. Comparing them gives sense versus antisense, which a consumer may derive and this
        accumulator never stores; storing it as one bool would put a dUTP library's spliced fragments in the sense
        column. ``sj_strand`` here is observed, per fragment; ``Partition.sj_strand`` is annotated, per sj boundary, and
        :meth:`_sj_edge_id` compares them.
        """
        # Checked first, before any geometry: a fragment with no single genome strand has no column in any bank.
        # Rejecting it here rather than upstream lets the loss be counted.
        admitted = self._admit(ref, start, end, align_strand)
        if isinstance(admitted, DepositOutcome):
            return self._reject(admitted)
        column, region_bounds, start, end = admitted

        # ── arbitration: which hypotheses survive, and is exactly one left? ───────────────────────
        #
        # After the strand and the clip, and the order is part of the contract: every fragment counts exactly once, and
        # a fragment can fail several ways. The second pass resolves the path, not the strand, so the strand rejection
        # must win over the deferral. The clip comes first because a hypothesis is filtered on its `L`, measured after
        # clipping.
        _scored, survivors = self._arbitrate(start, end, observed_introns, hypotheses)
        self._record_gap_resolution(hypotheses, survivors)
        if len(survivors) > 1:
            self.tally.deferred.append(
                DeferredFragment(
                    placements=(
                        self._placement(
                            ref,
                            start,
                            end,
                            align_strand,
                            sj_strand,
                            observed_introns,
                            hypotheses,
                            nm,
                            clip,
                        ),
                    )
                )
            )
            return self._reject(DepositOutcome.DEFERRED)
        return self._deposit_survivor(
            ref, column, region_bounds, start, end, sj_strand, survivors[0]
        )

    def offer(self, placements) -> DepositOutcome:
        """Offer one fragment with every placement it aligned to; deposit it, hold it, or reject it.

        A uniquely mapped fragment is the one-placement case and goes through :meth:`deposit` unchanged, counters
        included. A multimapping fragment arrives with one :class:`Placement` per resolved alignment. The rule:

        * per placement, the strand check and the clip exclude that placement alone (a placement with no genome strand
          has no column; one that clips to nothing is not on its reference). If none remains the fragment is rejected
          once: ``STRAND_UNDEFINED`` when every exclusion was the strand, ``EMPTY`` otherwise;
        * identical placements (same reference, clipped extent, strands, observed introns and hypotheses, which pairing
          secondary alignments can produce twice) are one placement, as two transcripts implying one path are one
          hypothesis;
        * exactly one placement remains → :meth:`deposit` at it with all its hypotheses: its one open question, if any,
          is its gap, held, counted and censused as a unique mapper's is;
        * several remain → the survivors are the union over placements of the (placement, hypothesis) pairs whose ``L``
          is within the fragment-length limit (the one filter, applied to the union because the molecule took exactly
          one element of it); if the union is empty every pair stands and the ordinary ``TOO_LONG`` rejection counts the
          chosen one at the drain. One survivor deposits at its placement with that hypothesis alone; two or more hold
          the fragment whole as :attr:`DepositOutcome.DEFERRED_PLACEMENTS`.

        A multi-placement fragment leaves the umbrella gap census alone, resolved or held: its classes describe how the
        gap at one placement was resolved, which has no single answer here, and the record carries its placements so a
        census is derivable from the bank.
        """
        admitted = []
        strand_rejected = empty_rejected = 0
        for placement in placements:
            result = self._admit(
                placement.ref, placement.start, placement.end, placement.align_strand
            )
            if isinstance(result, DepositOutcome):
                if result is DepositOutcome.STRAND_UNDEFINED:
                    strand_rejected += 1
                else:
                    empty_rejected += 1
                continue
            _column, _region_bounds, start, end = result
            clipped = self._placement(
                placement.ref,
                start,
                end,
                placement.align_strand,
                placement.sj_strand,
                placement.observed_introns,
                placement.hypotheses,
                placement.nm,
                placement.clip,
            )
            admitted.append(clipped)
        # Identical placements are one: of a group, the first in content-then-sequence order stays (the fewest
        # mismatches, then the fewest clipped bases). The deposit is the group's either way.
        admitted.sort(key=lambda pl: (pl.key, (pl.nm, pl.clip)))
        admitted = [pl for i, pl in enumerate(admitted) if i == 0 or pl.key != admitted[i - 1].key]
        if not admitted:
            return self._reject(
                DepositOutcome.STRAND_UNDEFINED
                if strand_rejected and not empty_rejected
                else DepositOutcome.EMPTY
            )
        if len(admitted) == 1:
            one = admitted[0]
            return self.deposit(
                one.ref,
                one.start,
                one.end,
                observed_introns=one.observed_introns,
                align_strand=one.align_strand,
                sj_strand=one.sj_strand,
                hypotheses=one.hypotheses,
                nm=one.nm,
                clip=one.clip,
            )
        # Within the record the placements sit in their canonical order (the sort above), so the record does not depend
        # on the order the aligner listed the hits in.
        union = [
            (pl, row)
            for pl in admitted
            for row in self._arbitrate(pl.start, pl.end, pl.observed_introns, pl.hypotheses)[0]
            if row[1] <= self.max_fragment_length
        ]
        if not union:
            union = [
                (pl, row)
                for pl in admitted
                for row in self._arbitrate(pl.start, pl.end, pl.observed_introns, pl.hypotheses)[0]
            ]
        if len(union) == 1:
            pl, (hypothesis, *_rest) = union[0]
            census_before = dict(self.tally.gap_resolution)
            outcome = self.deposit(
                pl.ref,
                pl.start,
                pl.end,
                observed_introns=pl.observed_introns,
                align_strand=pl.align_strand,
                sj_strand=pl.sj_strand,
                hypotheses=(hypothesis,),
                nm=pl.nm,
                clip=pl.clip,
            )
            self.tally.gap_resolution.update(census_before)
            return outcome
        self.tally.deferred.append(DeferredFragment(placements=tuple(admitted)))
        return self._reject(DepositOutcome.DEFERRED_PLACEMENTS)

    def _admit(self, ref: int, start: int, end: int, align_strand: int):
        """The strand check and the clip, counting nothing: ``(column, region_bounds, start, end)`` with the clipped
        extent, or the :class:`DepositOutcome` that rejects the placement."""
        column = STRAND_COLUMNS.get(align_strand)
        if column is None:
            return DepositOutcome.STRAND_UNDEFINED
        p = self.partition
        first_region_bound, last_region_bound = (
            int(p.ref_region_bound_offsets[ref]),
            int(p.ref_region_bound_offsets[ref + 1]),
        )
        if last_region_bound - first_region_bound < 2:
            return DepositOutcome.EMPTY
        region_bounds = p.region_bounds[first_region_bound:last_region_bound]
        # clip to the reference; L is the clipped length, so the placement count stays consistent
        start, end = max(int(start), int(region_bounds[0])), min(int(end), int(region_bounds[-1]))
        if end <= start:
            return DepositOutcome.EMPTY
        return column, region_bounds, start, end

    def _arbitrate(self, start, end, observed_introns, hypotheses):
        """Every hypothesis scored ``(hypothesis, L, cut introns)``, and the survivors of the one
        filter ``L <= max_fragment_length`` — all of them when the filter would empty the set."""
        scored = [
            (hypothesis, *self._hypothesis_length(start, end, observed_introns, hypothesis))
            for hypothesis in hypotheses
        ]
        survivors = [row for row in scored if row[1] <= self.max_fragment_length] or scored
        return scored, survivors

    @staticmethod
    def _placement(
        ref, start, end, align_strand, sj_strand, observed_introns, hypotheses, nm=0, clip=0
    ) -> Placement:
        return Placement(
            ref=int(ref),
            start=int(start),
            end=int(end),
            align_strand=int(align_strand),
            sj_strand=int(sj_strand),
            observed_introns=tuple((int(s), int(e)) for s, e in observed_introns),
            hypotheses=tuple(hypotheses),
            nm=int(nm),
            clip=int(clip),
        )

    def _deposit_survivor(self, ref, column, region_bounds, start, end, sj_strand, survivor):
        """The write-out of the one surviving hypothesis: ``L``, the pools, every channel of every object the path
        touches."""
        # `cut_introns`, not `introns`: the introns removed from the molecule, the observed ones unioned with the
        # surviving hypothesis's implied ones, normalised and clipped.
        hypothesis, length, cut_introns = survivor
        segments = _segments(start, end, cut_introns)
        p = self.partition
        if length <= 0:
            return self._reject(DepositOutcome.EMPTY)
        if length > self.max_fragment_length:
            return self._reject(DepositOutcome.TOO_LONG)
        # An implied strand is used only when no motif was sequenced anywhere on this fragment: an observed motif is
        # evidence, an implied strand an inference from the annotation.
        if sj_strand == Strand.NONE:
            sj_strand = hypothesis.sj_strand

        # ── which annotated sj does this path use? ─────────────────────────────────────────
        # A contradictory motif strand disqualifies all the fragment's splices, so it is checked once here: `sj_strand`
        # is the OR of the per-record tag values, so AMBIGUOUS means the mates disagreed about the same molecule, and no
        # splice is trusted. Resolved per intron position, with -1 where the annotation has no such sj, because the
        # conserved mass needs to know which of a block's two ends is a sj boundary; a filtered list would lose the
        # alignment between intron `i` and the gap between blocks `i` and `i+1`.
        if sj_strand == Strand.AMBIGUOUS:
            sj_id_at_gap: list[int] = [-1] * len(cut_introns)
        else:
            sj_id_at_gap = [
                self._sj_edge_id(ref, intron_start, intron_end, sj_strand)
                for intron_start, intron_end in cut_introns
            ]
        sj_ids = [jid for jid in sj_id_at_gap if jid >= 0]
        spliced = bool(sj_ids)

        # The path's own first and last covered base, not the fragment's extent. With introns [(150,480)] over [150,500)
        # the path is [480,500) and lives in a different region; using the extent would attribute the fragment, and the
        # start count, to a region it never touches.
        first_base, last_base = segments[0][0], segments[-1][1] - 1

        t = self.tally
        region_base, boundary_base = (
            int(p.ref_region_offsets[ref]),
            int(p.ref_boundary_offsets[ref]),
        )
        r_first = self._local_region(region_bounds, first_base)
        t.region_start_count[region_base + r_first, column] += 1
        # The unconditional length histogram, incremented here beside the start count and the deposited counter, so all
        # three describe one population by construction. ``length`` is already clipped to the reference and gated by the
        # length limit.
        t.deposited_lengths[length] += 1
        t.qc[DepositOutcome.DEPOSITED.value] += 1

        # ── crossings, per contiguous segment of the path ─────────────────────────────────────────
        # A boundary is crossed iff it lies strictly inside a segment, so per segment the crossed boundaries are a
        # contiguous index range. The spliced bank has two channels and the unspliced one four: a spliced crossing is
        # certified RNA and nothing deconvolves it, so its length moments have no consumer; the unspliced bank is the
        # mixture the length banks are built on.
        boundary_count = t.boundary_spliced_count if spliced else t.boundary_unspliced_count
        boundary_mass = t.boundary_spliced_mass if spliced else t.boundary_unspliced_mass
        # 0 at L == 1: a length-1 molecule cannot cross a 0-bp boundary, and 1/(L-1) would divide by 0.
        inv_boundary = 1.0 / (length - 1) if length >= 2 else 0.0
        n_crossed, sole_boundary = 0, -1
        for block, (seg_start, seg_end) in enumerate(segments):
            first = int(np.searchsorted(region_bounds, seg_start, side="right"))
            last = int(np.searchsorted(region_bounds, seg_end, side="left"))
            for boundary in range(first, last):
                boundary_count[boundary_base + boundary - 1, column] += 1
                if not spliced:
                    t.boundary_unspliced_inv_length_sum[boundary_base + boundary - 1] += (
                        inv_boundary
                    )
            # ── the conserved mass, per slice, over one boundary set ───────────────────────────────
            # The crossed bounds split this block into `last - first + 1` slices. Each slice's
            # `slice_len / length` is shared equally between the objects that bound it, so every bounded
            # slice disposes of exactly its own bases::
            #
            #     Sum over the fragment  =  Sum slice_len / length  =  1
            #
            # A sj bounds a block exactly as a contiguous boundary does. A block's interior boundaries are
            # the ones it crosses; its two ends are boundaries too whenever the intron there resolved to an
            # annotated sj. So a fragment's 1.0 is shared across every object it crosses, boundaries and sj
            # together. Giving a boundary-crossing block's bases entirely to boundaries would still conserve,
            # but would leave a sj whose two flanking blocks both crossed a boundary with nothing while
            # `sj_count` credits it.
            #
            # Coverage-weighted, not `1/K`: both conserve, but only this one says where the fragment sat and
            # is expressible per base. An unspliced path has no sj boundaries, so this reduces to the
            # boundary-only rule. A single block bounded by nothing deposits nothing: for a one-block path that
            # is the contained case, already whole in `region_contained_count`; for a multi-block one it is an
            # unannotated intron's block, whose bases have nowhere conserved to go.
            left_sj = sj_id_at_gap[block - 1] if block > 0 else -1
            right_sj = sj_id_at_gap[block] if block < len(sj_id_at_gap) else -1
            n_slices = last - first + 1
            for i in range(n_slices):
                lo = seg_start if i == 0 else int(region_bounds[first + i - 1])
                hi = seg_end if i == n_slices - 1 else int(region_bounds[first + i])
                left_boundary = first + i - 1 if i > 0 else -1
                right_boundary = first + i if i < n_slices - 1 else -1
                # The block's own ends are sj boundaries; its interior ends are boundaries. A slice has at most one
                # boundary of each kind on each side.
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
                        # `column`: the same column `sj_count` is deposited at below, so `mass[c]/count[c]` is a
                        # per-strand mean, not a cross-population ratio.
                        t.sj_mass[jid, column] += share
            if last > first:
                sole_boundary = first if (n_crossed == 0 and last - first == 1) else -1
                n_crossed += last - first

        for jid in sj_ids:
            t.sj_count[jid, column] += 1
            # `inv_length_sum` only: nothing reads a length sum here, and `pool_lengths` already carries the spliced
            # population's length distribution (RNA_SPLICED).
            t.sj_inv_length_sum[jid] += inv_boundary

        # ── contained: the whole path lies inside one region ────────────────────────────────────────
        # Not merely "crossed no boundary". An unannotated intron can swallow every boundary between two blocks, leaving
        # a fragment that crosses nothing yet straddles two regions; crediting it as contained would place its whole
        # length in a region it only partly overlaps. Such a fragment deposits on no object but is still counted
        # (start_count), so the loss is visible.
        first_region = self._local_region(region_bounds, first_base)
        contained_region = -1
        if not sj_ids and first_region == self._local_region(region_bounds, last_base):
            contained_region = region_base + first_region
            t.region_contained_count[contained_region, column] += 1

        pool = self._pool(spliced, contained_region, sole_boundary, region_base)
        if pool is not None:
            t.pool_lengths[pool, length] += 1
        return DepositOutcome.DEPOSITED

    # -- the drain --------------------------------------------------------------------------------

    def drain(self, choices) -> dict[str, int]:
        """The second pass's drain: replay each held fragment with one chosen hypothesis.

        ``choices[i]`` indexes the ``i``-th record's hypotheses over its placements in order
        (:meth:`Tally.deferred_canonical`, the queue's one canonical order), so it names a placement and a hypothesis at
        once.

        Each record re-enters :meth:`deposit` with its chosen hypothesis alone, a set of size one, so the arbitration is
        degenerate and the fragment deposits or is rejected by the ordinary rules. Every crossing rule, the pools and
        ``L`` are reached through the code that ran in pass one, so byte-identity with the native accumulator follows.

        The drain's outcomes stay out of ``gap_resolution``. The census has no ``gap_resolved_unspliced`` class because
        pass-one arbitration cannot produce one, but the drain can choose ∅. And ``deposit`` reaches
        ``_record_gap_resolution`` on every fragment: a size-one set lands on ``RESOLVED_SPLICED``, while a chosen ∅
        hits the ``all(h.is_unspliced)`` early return, so an unguarded drain would move the census with the RNG. This
        method snapshots the census and restores it; what the census cannot express lives on the drain's own axis as
        ``chose_no_implied_intron`` / ``chose_implied_intron``.

        After the drain the tally describes the final state: nothing is held, so the bank is empty,
        ``deferred_undetermined_gap`` is 0 and the three ``gap_deferred_*`` are 0. The counter and the fragments it
        counts stay one population (the payload refuses a bank that disagrees with its counter). Pass one's numbers are
        returned, and the payload keeps them in its ``drain`` block.

        Returns the drain's own counters plus ``census_before``. The conservation statement is ``deposited + dropped_*
        == offered``.
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
            "offered_multimapper": 0,
            "deposited": 0,
            "dropped_too_long": 0,
            "dropped_empty": 0,
            "dropped_strand_undefined": 0,
            "chose_no_implied_intron": 0,
            "chose_implied_intron": 0,
            "census_before": dict(self.tally.gap_resolution),
        }
        census_before = dict(self.tally.gap_resolution)
        # Emptied first. `deposit` appends to this list when it defers, and a set of size one cannot defer; draining
        # into an empty queue lets the assertion below check that.
        self.tally.deferred = []
        for fragment, choice in zip(held, choices):
            # The choice is a local index into the record's run of hypotheses over its placements, in the order the
            # flattening emits them, so it names a placement and a hypothesis at once.
            placement, path = self._locate(fragment, int(choice))
            chose = "chose_no_implied_intron" if path.is_unspliced else "chose_implied_intron"
            counters[chose] += 1
            counters["offered_multimapper"] += int(fragment.n_placements > 1)
            outcome = self.deposit(
                placement.ref,
                placement.start,
                placement.end,
                observed_introns=placement.observed_introns,
                align_strand=placement.align_strand,
                sj_strand=placement.sj_strand,
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

        # The census restored to pass one's, and the held population now zero everywhere it is counted.
        self.tally.gap_resolution.update(census_before)
        for key in (
            GapResolution.DEFERRED_RNA_OR_GDNA,
            GapResolution.DEFERRED_WHICH_INTRONS,
            GapResolution.DEFERRED_BOTH,
        ):
            self.tally.gap_resolution[key.value] = 0
        self.tally.qc[DepositOutcome.DEFERRED.value] = 0
        self.tally.qc[DepositOutcome.DEFERRED_PLACEMENTS.value] = 0
        return counters

    @staticmethod
    def _locate(fragment: DeferredFragment, choice: int) -> tuple[Placement, GapHypothesis]:
        """The (placement, hypothesis) a record-local hypothesis index names."""
        remaining = choice
        for placement in fragment.placements:
            if remaining < len(placement.hypotheses):
                return placement, placement.hypotheses[remaining]
            remaining -= len(placement.hypotheses)
        raise IndexError(
            f"choice {choice} is outside the record's {sum(len(p.hypotheses) for p in fragment.placements)} hypotheses"
        )

    # -- helpers ----------------------------------------------------------------------------------

    def _reject(self, outcome: DepositOutcome) -> DepositOutcome:
        self.tally.qc[outcome.value] += 1
        return outcome

    @staticmethod
    def _hypothesis_length(start, end, observed, path: GapHypothesis):
        """``(L, normalised introns)`` for one hypothesis.

        One definition of ``L`` and one code path to it, whether the hypothesis wins or is only scored for the filter.
        The observed and implied introns are disjoint by construction (the enumeration never searches a gap the CIGAR
        already explained), so this is a union, and ``_normalise_introns`` sorts and clips it.
        """
        introns = _normalise_introns(tuple(observed) + path.introns, start, end)
        return sum(b - a for a, b in _segments(start, end, introns)), introns

    def _record_gap_resolution(self, hypotheses, survivors) -> None:
        """The umbrella census, incremented at this one site so the subclasses cannot drift. Counted only when the
        enumeration had something to arbitrate: a fragment whose only hypothesis is the unspliced one had no question to
        answer, and counting it would make the umbrella's denominator the whole library.
        """
        if all(h.is_unspliced for h in hypotheses):
            return
        n_spliced = sum(1 for h, *_ in survivors if not h.is_unspliced)
        unspliced_survives = any(h.is_unspliced for h, *_ in survivors)
        if len(survivors) == 1:
            # The sole survivor is necessarily the spliced one: the genomic path is always the longest, so it can never
            # be the last one left (see :class:`GapResolution`).
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
        return min(
            max(int(np.searchsorted(region_bounds, position, side="right")) - 1, 0),
            region_bounds.size - 2,
        )

    def _sj_edge_id(self, ref: int, intron_start: int, intron_end: int, sj_strand: int) -> int:
        """The sj-boundary id for this intron, or -1 if it is not annotated.

        One to three iterations at human scale: most region bounds are not a donor at all, and the mean
        fan-out over those that are is near one.
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
            # The strand filter applies only when the motif strand is definite (POS or NEG); any other value matches on
            # coordinates alone, as the resolver's `sj_lookup_into` does (it falls back to the union of POS and NEG).
            # `deposit` never passes AMBIGUOUS: it rejects such a fragment's splices before the lookup. The filter only
            # separates two sj sharing a coordinate pair, and no human sj coordinate pair is annotated on both strands.
            if sj_strand in STRAND_COLUMNS and int(p.sj_strand[k]) != int(sj_strand):
                continue
            return k
        return -1

    def _pool(self, spliced, contained_region, sole_boundary, region_base):
        """The one length pool this fragment belongs to, or ``None``.

        Priority, by structure: a splice is certified RNA; a contained fragment is typed by its region; a
        single-boundary crossing is a "splash" read typed by its two flank types. Anything else (an exonic contained
        fragment, a multi-boundary crossing) is a mixture and enters nothing.

        Membership is determinacy, not provenance: a fragment reaches this point only when exactly one hypothesis
        survived, so its ``L`` is not in doubt, whether its splice was sequenced or inferred. Barring inferred splices
        would preferentially exclude fragments whose mates sit far apart and bias the pool's own length distribution.
        """
        types = self.partition.region_types
        if spliced:
            return FragmentPool.RNA_SPLICED
        if contained_region >= 0:
            return _CONTAINED_POOL.get(int(types[contained_region]))
        if sole_boundary >= 0:
            flanks = (
                int(types[region_base + sole_boundary - 1]),
                int(types[region_base + sole_boundary]),
            )
            return _SPLASH_POOL.get(tuple(sorted(flanks)))
        return None
