"""rigel.scan_payload — the BAM scan's accumulator tally, as a typed Python object.

The C++ scanner returns ``result["calibration"]`` as a dict of flat ndarrays (``BamScanner::build_result``); this
module validates the schema and reshapes the two-column banks, and everything downstream reads only this object. The
field names are the specification's ``Tally`` names (``tests/native/_accumulator_reference.py``) and the C++ keys,
character for character.

Three axes, off by one from each other per reference::

    region_bounds    0        100       200       600        c = 4 region_bounds on this reference
    regions   [  n0  ][   n1   ][   n2   ]            c - 1 = 3 regions
    boundaries            boundary 1    boundary 2               c - 2 = 2 contiguous boundaries

A reference with fewer than two region bounds owns no region and no boundary. Splice-junction (sj) boundaries are their
own axis, sliced by ``ref_sj_offsets`` and concatenated in reference order, so a junction's id is its slot.

A count is an integer and a fraction is float64, so no bank needs decoding::

    count           Sum 1                    integer   exact at any worker count
    inv_length_sum  Sum 1/A(w)               float64   A = w−1 crossing a 0-bp boundary or sj
    mass            Sum slice_len/(L·bounds) float64   the conserved fragment count

The float banks agree across worker counts to ~1e-15, not bit for bit, since float addition is not associative. At a
boundary the opportunity ``w−1`` and the deposit ``1/(w−1)`` cancel, so ``inv_length_sum`` is an exact, model-free
density. The trailing ``2`` on every bank is the genome strand, ``Strand.POS`` then ``Strand.NEG``. This object holds
views over buffers the C++ owns and keeps them alive; a cast that copied would double peak memory.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Any

import numpy as np

from .stats import library_mismatch_rate


#: The two columns of every bank are the two genome strands. Of `Strand`'s four values only POS and NEG name a
#: column; the accumulator rejects a fragment carrying neither.
N_STRAND_COLUMNS = 2

#: Five fragment-length pools, each defined by structure, in the order of `rigel::accumulator::FragmentPool`.
N_FRAGMENT_POOLS = 5

#: The pool axis, named: the accumulator's enum, kept with the schema so that no consumer holds its own copy.
#: The two DNA_* contained pools are gDNA-dominated, not pure (`calibration.fl` deconvolves them by a two-pool
#: contrast). RNA_SPLICED used an annotated sj on the fragment's one surviving path, sequenced or implied; gDNA cannot
#: be spliced. The two *_EXON crossing pools hold the gDNA that capture moves out of containment and enter the gDNA
#: model through their own opportunity. An exonic contained fragment and a multi-boundary crossing have no pool:
#: they are gDNA/RNA mixtures by structure.
POOL_DNA_INTERGENIC = 0  # contained in an intergenic region — gDNA-dominated
POOL_DNA_INTRONIC = 1  # contained in an intronic region — gDNA-dominated
POOL_DNA_INTRON_EXON = 2  # crossing one boundary, flanks {intron, exon} — on-target gDNA
POOL_DNA_INTERGENIC_EXON = 3  # crossing one boundary, {intergenic, exon} — on-target gDNA
POOL_RNA_SPLICED = 4  # used an annotated sj on its one surviving path — certified RNA


@dataclass(frozen=True, slots=True)
class ScanQC:
    """The accumulator's QC counters, named as the specification's ``Tally.qc`` keys.

    Typed so that a misspelled key fails here instead of reading as zero far away.
    """

    deposited: int  # fragments that reached an object; == sum(region_start_count)
    dropped_too_long: int  # L above --max-fragment-length
    dropped_empty: int  # no path left after clipping to the reference
    dropped_strand_undefined: int  # align_strand named no column, so there was none to credit
    #: More than one hypothesis survived, so the gap is undetermined; held whole for the second pass, not dropped.
    deferred_undetermined_gap: int
    #: More than one placement admits a surviving (placement, hypothesis), so the placement is undetermined; held
    #: whole. ``deposited + deferred_undetermined_gap + deferred_multiple_placements + dropped_* == offered``.
    deferred_multiple_placements: int

    @classmethod
    def from_dict(cls, qc: dict[str, Any]) -> "ScanQC":
        expected = {field.name for field in dataclasses.fields(cls)}
        missing = expected - set(qc)
        if missing:
            raise ValueError(f"the scan's qc block is missing {sorted(missing)}.")
        return cls(**{name: int(qc[name]) for name in expected})


@dataclass(frozen=True, slots=True)
class GapCensus:
    """How each fragment whose gap needed resolving was resolved, named as the ``Tally.gap_resolution`` keys.

    A fragment is counted when the enumeration produced at least one spliced hypothesis for it, since its ``L`` then
    depends on whether a gap intron is cut. The classes are exhaustive and mutually exclusive, and the three
    ``gap_deferred_*`` sum to ``qc.deferred_undetermined_gap``. The census cuts across the splice census rather than
    extending ``splice_type``: a spliced fragment with an intron in its mate gap needs resolving as much as an
    unspliced one, and the splice census must still sum to the library.

    There is no ``gap_resolved_unspliced`` because it is unreachable. A spliced hypothesis cuts bases the unspliced
    one keeps, so ``L_spliced <= L_unspliced``, and the one arbitration filter is ``L <= max_fragment_length``: if the
    unspliced path survives, every spliced path does too.
    """

    #: One hypothesis survived, and it cuts the gap intron. Classifies the arbitration, not the deposit: the fragment
    #: can still be rejected afterwards as too long.
    gap_resolved_spliced: int
    #: The unspliced path against exactly one spliced path: the open question is RNA or gDNA.
    gap_deferred_rna_or_gdna: int
    #: Two or more spliced paths and no unspliced one: certified RNA, and the open question is which structure.
    gap_deferred_which_introns: int
    #: Both questions at once.
    gap_deferred_both: int

    @classmethod
    def zeros(cls) -> "GapCensus":
        """No fragment had a gap to resolve (a real state: no annotated intron in any mate gap)."""
        return cls(**{field.name: 0 for field in dataclasses.fields(cls)})

    @classmethod
    def from_dict(cls, census: dict[str, Any]) -> "GapCensus":
        expected = {field.name for field in dataclasses.fields(cls)}
        missing = expected - set(census)
        if missing:
            raise ValueError(
                f"the scan's gap_resolution block is missing {sorted(missing)}. The subclasses are "
                f"exhaustive by construction, so a missing one is a partition that does not close."
            )
        return cls(**{name: int(census[name]) for name in expected})

    @property
    def deferred(self) -> int:
        """The three ``gap_deferred_*`` summed — which must equal ``qc.deferred_undetermined_gap``."""
        return (
            self.gap_deferred_rna_or_gdna + self.gap_deferred_which_introns + self.gap_deferred_both
        )


@dataclass(frozen=True, slots=True)
class DrainQC:
    """What the second pass's drain did: the audit trail of how the final tally was reached.

    After the drain nothing is held: the bank is empty and the held counters are 0, so a counter and the fragments it
    counts stay one population. Pass one's numbers live here instead.

    The drain keeps its own census rather than extending :class:`GapCensus`: it can choose the genomic path, which
    pass-one arbitration never leaves as the sole survivor. ``chose_no_implied_intron`` / ``chose_implied_intron``
    record that composition.
    """

    #: Pass one's ``qc.deferred_undetermined_gap + qc.deferred_multiple_placements``, the conservation's denominator.
    offered: int
    #: Of those, the records held on their placements (multimappers); the rest were held on their gap.
    offered_multimapper: int
    deposited: int
    dropped_too_long: int
    dropped_empty: int
    dropped_strand_undefined: int
    #: The held records whose chosen hypothesis implies no gap intron, against those whose chosen hypothesis implies
    #: one. Sums to ``offered``, not ``deposited``: a chosen hypothesis can still be rejected.
    chose_no_implied_intron: int
    chose_implied_intron: int
    #: Pass one's arbitration census, verbatim, so the per-class before and after stay readable.
    census_before: GapCensus

    @property
    def dropped(self) -> int:
        return self.dropped_too_long + self.dropped_empty + self.dropped_strand_undefined

    @property
    def conserved(self) -> bool:
        """Every offered fragment either deposited or was rejected, exactly once."""
        return self.deposited + self.dropped == self.offered

    def __post_init__(self) -> None:
        if not self.conserved:
            raise ValueError(
                f"the drain offered {self.offered} fragments but accounts for "
                f"{self.deposited + self.dropped} (deposited {self.deposited} + dropped {self.dropped}); "
                f"a drained fragment that is neither deposited nor rejected has been lost."
            )
        n_chosen = self.chose_no_implied_intron + self.chose_implied_intron
        if n_chosen != self.offered:
            raise ValueError(
                f"the drain chose {n_chosen} hypotheses for "
                f"{self.offered} offered fragments; exactly one hypothesis wins each whole fragment."
            )
        if not 0 <= self.offered_multimapper <= self.offered:
            raise ValueError(
                f"the drain reports {self.offered_multimapper} multimapping records among {self.offered} "
                f"offered; the multimappers are a subset of the held population."
            )


#: Every two-column bank, with its axis and dtype. One table for `from_scan_result`'s validation and the drain's
#: delta, so a new channel cannot reach one and miss the other.
BANK_AXES: tuple[tuple[str, str, Any], ...] = (
    ("region_contained_count", "region", np.uint32),
    # the start bank: the path's first covered base, by align strand; sums over both columns to qc.deposited
    ("region_start_count", "region", np.uint32),
    ("boundary_unspliced_count", "boundary", np.uint32),
    ("boundary_spliced_count", "boundary", np.uint32),
    ("sj_count", "sj", np.uint32),
    # The one mass with a strand axis. An artifactual sj accumulates on both strands like gDNA, so a per-strand mass
    # lets the strand model detect one (`SpliceJunction::mass`). Its columns are `sj_count`'s, so `mass[c] / count[c]`
    # is a per-strand mean.
    ("sj_mass", "sj", np.float64),
)


#: The single-column additive banks, all float64 fractions, with the axis and dtype `from_scan_result` validates.
#: Kept apart from `BANK_AXES`, whose contract is two columns: merging them would give these a strand column in every
#: loop that reads that table.
SINGLE_COLUMN_AXES: tuple[tuple[str, str, Any], ...] = (
    ("boundary_unspliced_inv_length_sum", "boundary", np.float64),
    ("sj_inv_length_sum", "sj", np.float64),
    ("boundary_unspliced_mass", "boundary", np.float64),
    ("boundary_spliced_mass", "boundary", np.float64),
)

#: Every additive channel with its axis (``"library"``: library-wide, not per reference), derived from the two tables
#: above so the drain adds every one. Only `region_start_count` and `deposited_lengths` have a checkable sum
#: (`qc.deposited`); a skipped `pool_lengths` would go short silently.
ADDITIVE_AXES: tuple[tuple[str, str], ...] = (
    *((name, axis) for name, axis, _dtype in BANK_AXES),
    *((name, axis) for name, axis, _dtype in SINGLE_COLUMN_AXES),
    ("pool_lengths", "library"),
    ("deposited_lengths", "library"),
)


#: The deferred bank's per-placement arrays, in the order :meth:`Tally.deferred_arrays` emits them. The last two are
#: the hit's sequence, which the canonical sort reads only when a record's content otherwise ties.
DEFERRED_RECORD_FIELDS = ("ref", "start", "end", "align_strand", "sj_strand", "nm", "clip")


@dataclass(frozen=True, slots=True)
class DeferredFragments:
    """The side buffer: fragments held whole for the second pass.

    A fragment is held when its unsequenced mate gap has more than one surviving explanation (which introns it holds
    cannot be observed, and deciding needs a fragment-length law pass one has not yet built), or when more than one of
    a multimapper's placements admits a surviving hypothesis. The fragment itself is stored, not its consequences: it
    is small and replays exactly through ``Accumulator::deposit``.

    Three nested levels: fragments hold placements; placements hold observed introns and hypotheses; hypotheses hold
    implied introns. A placement is one alignment on one reference with its own extent, strands, observed introns and
    gap hypotheses; a gap-held unique mapper has one. Offsets are cumulative and start at 0, so an empty level is
    ``[0]``.

    This is the one bank whose order is observable, so the C++ export sorts it on each record's own content before it
    crosses the ABI; otherwise the same BAM would give different bytes at different worker counts. Records that tie on
    that key are identical, so no tie-break is needed. Every array is ``int64``, the strand columns included.
    """

    placement_offsets: np.ndarray  # int64[n_fragments + 1] — fragments → placements
    ref: np.ndarray  # int64[n_placements] — the reference the drain replays onto
    start: np.ndarray  # int64[n_placements] — the clipped extent, which the drain replays
    end: np.ndarray  # int64[n_placements]
    align_strand: np.ndarray  # int64[n_placements]
    sj_strand: np.ndarray  # int64[n_placements] — the observed motif strand, if any
    observed_intron_offsets: np.ndarray  # int64[n_placements + 1] — in pairs, into observed_introns
    observed_introns: (
        np.ndarray
    )  # int64[2 * n_observed] — flat (start, end); cut under every hypothesis
    hypothesis_offsets: np.ndarray  # int64[n_placements + 1] — placements → hypotheses
    hypothesis_sj_strand: (
        np.ndarray
    )  # int64[n_hypotheses] — inferred; an observed motif always wins
    hypothesis_intron_offsets: np.ndarray  # int64[n_hypotheses + 1] — in pairs
    hypothesis_introns: np.ndarray  # int64[2 * n_implied] — the implied introns; empty ⇒ unspliced
    #: int64[n_placements] — the hit's edits and soft-clipped read bases over both mates. Never deposited: the draw
    #: prices a multimapper's placements by them, as the EM prices its hits.
    nm: np.ndarray
    clip: np.ndarray

    @property
    def n_fragments(self) -> int:
        return int(self.placement_offsets.shape[0]) - 1

    @property
    def n_placements(self) -> int:
        return int(self.ref.shape[0])

    @property
    def n_hypotheses(self) -> int:
        return int(self.hypothesis_sj_strand.shape[0])

    @property
    def record_hypothesis_offsets(self) -> np.ndarray:
        """``int64[n_fragments + 1]``: each record's run of the flat hypothesis axis over its placements; a
        choice is a local index into this run."""
        return self.hypothesis_offsets[self.placement_offsets]

    def placements_of(self, i: int) -> range:
        """The placements of fragment ``i``, as indices into the per-placement arrays."""
        return range(int(self.placement_offsets[i]), int(self.placement_offsets[i + 1]))

    def placement_of_hypothesis(self, h: int) -> int:
        """Which placement flat hypothesis ``h`` belongs to."""
        return int(np.searchsorted(self.hypothesis_offsets, h, side="right")) - 1

    @classmethod
    def empty(cls) -> "DeferredFragments":
        """Nothing was deferred: every offset array is ``[0]`` and every other array is empty."""
        return cls.from_dict(
            {
                field.name: np.asarray(
                    [0] if field.name.endswith("_offsets") else [], dtype=np.int64
                )
                for field in dataclasses.fields(cls)
            }
        )

    @classmethod
    def from_dict(cls, deferred: dict[str, Any]) -> "DeferredFragments":
        """Validate and adopt the flat arrays the C++ emits.

        The nested CSRs are checked rather than trusted: an offset array of the right length can still be
        inconsistent, and a truncated bank would score a fragment against another fragment's hypotheses.
        """
        arrays = _int64_columns(deferred, [field.name for field in dataclasses.fields(cls)])
        n = int(arrays["placement_offsets"].shape[0]) - 1
        if n < 0:
            raise ValueError(
                "deferred['placement_offsets'] is empty; offsets are cumulative and start at 0, so an "
                "empty bank is [0] and never []"
            )
        n_placements = int(arrays["placement_offsets"][-1])
        for name in DEFERRED_RECORD_FIELDS:
            if arrays[name].shape != (n_placements,):
                raise ValueError(
                    f"deferred[{name!r}] has {arrays[name].shape[0]} entries but placement_offsets "
                    f"describes {n_placements} placements"
                )
        for name in ("observed_intron_offsets", "hypothesis_offsets"):
            if arrays[name].shape[0] != n_placements + 1:
                raise ValueError(
                    f"deferred[{name!r}] has {arrays[name].shape[0]} entries but there are "
                    f"{n_placements} placements, so it must have {n_placements + 1}"
                )
        for offsets_name, values_name, stride in (
            ("placement_offsets", "ref", 1),
            ("observed_intron_offsets", "observed_introns", 2),
            ("hypothesis_offsets", "hypothesis_sj_strand", 1),
            ("hypothesis_intron_offsets", "hypothesis_introns", 2),
        ):
            _check_csr(arrays, offsets_name, values_name, stride)
        n_hypotheses = int(arrays["hypothesis_offsets"][-1])
        if arrays["hypothesis_intron_offsets"].shape[0] != n_hypotheses + 1:
            raise ValueError(
                f"deferred['hypothesis_intron_offsets'] has "
                f"{arrays['hypothesis_intron_offsets'].shape[0]} entries but there are "
                f"{n_hypotheses} hypotheses, so it must have {n_hypotheses + 1}"
            )
        # A record is held because two or more (placement, hypothesis) pairs survived; fewer means the bank lost some.
        runs = np.diff(arrays["hypothesis_offsets"][arrays["placement_offsets"]])
        if n and int(runs.min()) < 2:
            bad = int(np.argmin(runs))
            raise ValueError(
                f"deferred fragment {bad} carries {int(runs[bad])} hypotheses over its placements. A "
                f"fragment is deferred only when two or more survived, so every record must hold at least two."
            )
        return cls(**arrays)

    def observed_introns_of(self, p: int) -> np.ndarray:
        """Placement ``p``'s observed introns as a ``[k, 2]`` view; cut under every hypothesis."""
        lo, hi = int(self.observed_intron_offsets[p]), int(self.observed_intron_offsets[p + 1])
        return self.observed_introns[2 * lo : 2 * hi].reshape(hi - lo, 2)

    def hypothesis_introns_of(self, h: int) -> np.ndarray:
        """Hypothesis ``h``'s implied introns as a ``[k, 2]`` view; empty ⇒ the unspliced hypothesis."""
        lo, hi = int(self.hypothesis_intron_offsets[h]), int(self.hypothesis_intron_offsets[h + 1])
        return self.hypothesis_introns[2 * lo : 2 * hi].reshape(hi - lo, 2)


def _int64_columns(deferred: dict[str, Any], names: list[str]) -> dict[str, np.ndarray]:
    """The deferred bank's arrays, each present, ``int64`` and one-dimensional."""
    missing = set(names) - set(deferred)
    if missing:
        raise ValueError(
            f"the scan's deferred bank is missing {sorted(missing)}. Every array is part of one "
            f"record, so a missing one is a fragment the second pass cannot replay."
        )
    arrays: dict[str, np.ndarray] = {}
    for name in names:
        array = np.ascontiguousarray(deferred[name])
        if array.dtype != np.int64:
            raise ValueError(
                f"deferred[{name!r}] has dtype {array.dtype}, expected int64. One bank, one dtype — "
                f"a narrowed or widened array compares equal by value and would hide the change."
            )
        if array.ndim != 1:
            raise ValueError(f"deferred[{name!r}] has shape {array.shape}, expected one dimension")
        arrays[name] = array
    return arrays


def _check_csr(
    arrays: dict[str, np.ndarray], offsets_name: str, values_name: str, stride: int
) -> None:
    """One CSR level of the deferred bank: offsets that start at 0, never decrease and end at the values' count."""
    offsets = arrays[offsets_name]
    if offsets.shape[0] == 0 or int(offsets[0]) != 0:
        raise ValueError(
            f"deferred[{offsets_name!r}] must start at 0; offsets are cumulative and an empty "
            f"level is [0]"
        )
    if np.any(np.diff(offsets) < 0):
        bad = int(np.argmax(np.diff(offsets) < 0))
        raise ValueError(
            f"deferred[{offsets_name!r}] decreases at {bad}; a CSR offset array cannot go backwards"
        )
    if int(offsets[-1]) * stride != arrays[values_name].shape[0]:
        raise ValueError(
            f"deferred[{offsets_name!r}] ends at {int(offsets[-1])} but "
            f"{values_name} has {arrays[values_name].shape[0]} entries "
            f"({int(offsets[-1])} x {stride} expected)"
        )


@dataclass(frozen=True, slots=True)
class AccumulatorPayload:
    """One BAM scan's tally. Views over C++-owned buffers; this object is the keep-alive."""

    # -- the partition, echoed back so a consumer can locate every object without reloading the index --
    region_bounds: (
        np.ndarray
    )  # int64[n_region_bounds] — flat, reference-major, ascending per reference
    ref_region_bound_offsets: np.ndarray  # int64[n_refs + 1] — CSR over region_bounds
    ref_region_offsets: np.ndarray  # int64[n_refs + 1]
    ref_boundary_offsets: np.ndarray  # int64[n_refs + 1] — contiguous boundaries
    ref_sj_offsets: np.ndarray  # int64[n_refs + 1] — sj boundaries

    # -- regions: two banks, each two genome-strand columns; a contained path books both --
    region_contained_count: (
        np.ndarray
    )  # uint32[n_regions, 2] — the whole path lies inside the region
    #: uint32[n_regions, 2] — the path's first covered base, by align strand; sums over both columns to
    #: qc.deposited. Opportunity ℓ for every fragment length; wall-blind only at the template's downstream end.
    region_start_count: np.ndarray

    # -- contiguous boundaries: the 0-bp boundary between two adjacent regions --
    boundary_unspliced_count: np.ndarray  # uint32[n_boundaries, 2] — the mixture being deconvolved
    boundary_unspliced_inv_length_sum: (
        np.ndarray
    )  # float64[n_boundaries] — one column, strand-agnostic
    #: float64[n_boundaries] — the conserved mass. ``boundary_unspliced_count`` adds 1 on every boundary a fragment
    #: crosses; this sums to one per fragment, so a consumer can turn an object total into a fragment count. One
    #: column: the question it answers has no strand.
    boundary_unspliced_mass: np.ndarray
    #: uint32[n_boundaries, 2] — certified RNA, since gDNA cannot be spliced. Count and mass only: nothing
    #: deconvolves a certified-RNA crossing, so its length moments are not kept.
    boundary_spliced_count: np.ndarray
    #: float64[n_boundaries] — the mass rule, routed by the same spliced flag. A per-boundary partial, not a ledger
    #: and not "the number of spliced fragments here".
    boundary_spliced_mass: np.ndarray

    # -- sj boundaries: one exact donor->acceptor jump. Pure RNA by construction --
    #: uint32[n_sj, 2] — both genome-strand columns are kept although consumers sum them, because the split detects
    #: aligner artifacts (false-positive ``N`` operations on genomic DNA) that ``rigel.splice_blacklist`` does not
    #: list. In a stranded library a real sj inherits the library's strand specificity, while an artifact deposits
    #: on both strands; a clean and an artifactual sj carry the same total, so the information is only in the split.
    sj_count: np.ndarray
    #: float64[n_sj] — read by ``second_pass`` to score a held fragment's sj evidence.
    sj_inv_length_sum: np.ndarray
    #: float64[n_sj, 2] — the conserved mass's third axis. A spliced fragment whose block holds no interior boundary
    #: deposits on neither boundary bank and lies in no single region, so without this it would be on no conserved
    #: axis and the RNA side's fragment count would fall short. A slice bounded by a sj and a boundary shares its
    #: bases between them.
    sj_mass: np.ndarray

    # -- the fragment-length pools, binned at L, once per fragment --
    pool_lengths: np.ndarray  # int64[5, max_length + 1]

    #: uint32[max_length + 1] — every deposited fragment, binned at its own L, with no purity condition. The five
    #: pools are conditioned, so none is unconditional; this is the anchor `calibration.fl` builds the length laws
    #: against (gDNA's pool shrinks toward it, and an empty RNA pool returns it). Unconditional given deposit: it
    #: excludes what the accumulator rejects, which is exactly the population the pools are drawn from.
    deposited_lengths: np.ndarray

    #: The side buffer: exactly the fragments the two ``qc.deferred_*`` counters count, so
    #: ``deposited + deferred_* + dropped_* == offered``.
    deferred: DeferredFragments
    #: How each gap was resolved — its own axis, cutting across the splice census.
    gap_resolution: GapCensus

    qc: ScanQC
    max_length: int  # the fragment-length limit applied to L, and the pool-histogram width
    n_refs: int
    #: The library's unique unspliced hits: their edits and aligned read bases (M, I, =, X), so a cached scan drains
    #: at its own mismatch rate.
    mismatches: int
    aligned_bases: int

    #: Provenance over regions and boundaries: the sj axis is meaningless against a different sj CSR, and an index
    #: rebuild can rewrite ``edges.feather`` while ``regions.feather``, all `index.partition_hash` covers, stays
    #: byte-identical. `None` when the scanner ran without an index to hash.
    graph_hash: str | None = None

    #: Set once the second pass has drained this tally, ``None`` while the side buffer is held. A drained payload has
    #: the same shape (an empty bank, zero held counters), so this is the only way to tell the two apart; pass one's
    #: numbers are in it (:class:`DrainQC`).
    drain: DrainQC | None = None

    # -- derived, never stored ------------------------------------------------------------------------

    @property
    def mismatch_rate(self) -> float | None:
        """The library's mismatch rate, which the draw prices a multimapper's sequence at."""
        return library_mismatch_rate(self.mismatches, self.aligned_bases)

    @property
    def n_regions(self) -> int:
        return int(self.region_start_count.shape[0])

    @property
    def n_boundaries(self) -> int:
        return int(self.boundary_unspliced_count.shape[0])

    @property
    def n_sj(self) -> int:
        return int(self.sj_count.shape[0])

    def with_drain(self, delta: dict[str, np.ndarray], drain: DrainQC) -> "AccumulatorPayload":
        """Pass one plus the drain's delta: the tally calibration consumes.

        ``delta`` holds one globally shaped array per additive channel, from replaying the held fragments through
        ``Accumulator.deposit``; this method does the arithmetic and the bookkeeping. Both payloads exist afterwards,
        so the drain's contribution to every channel is a subtraction. The result's bank is empty, its held counters
        are 0, and ``drain`` records what pass one held.
        """
        if self.drain is not None:
            raise ValueError(
                "this payload has already been drained. The drain consumes the side buffer, so a second "
                "one would deposit nothing and silently double the bookkeeping."
            )
        missing = {name for name, _axis in ADDITIVE_AXES} - set(delta)
        if missing:
            raise ValueError(
                f"the drain delta is missing {sorted(missing)}; every additive channel must be present, "
                f"because a silently absent one reads as a tally that simply saw fewer fragments."
            )

        totals: dict[str, np.ndarray] = {}
        for name, _axis in ADDITIVE_AXES:
            before = getattr(self, name)
            added = np.ascontiguousarray(delta[name], dtype=before.dtype)
            if added.shape != before.shape:
                raise ValueError(
                    f"the drain delta for {name!r} has shape {added.shape}, expected {before.shape}"
                )
            total = before + added
            # Both terms are non-negative, so a sum below either one has wrapped; a uint32 count can wrap on a deep
            # library, and would then read as a plausible small number.
            if np.any(total < before) or np.any(total < added):
                raise ValueError(
                    f"adding the drain's delta to {name!r} overflowed {before.dtype}; the tally cannot "
                    f"represent pass one and the drain together."
                )
            totals[name] = total

        return dataclasses.replace(
            self,
            **totals,
            # The held counters go to 0 with the bank; pass one's count is kept as `drain.offered`.
            qc=dataclasses.replace(
                self.qc,
                deposited=self.qc.deposited + drain.deposited,
                dropped_too_long=self.qc.dropped_too_long + drain.dropped_too_long,
                dropped_empty=self.qc.dropped_empty + drain.dropped_empty,
                dropped_strand_undefined=(
                    self.qc.dropped_strand_undefined + drain.dropped_strand_undefined
                ),
                deferred_undetermined_gap=0,
                deferred_multiple_placements=0,
            ),
            # `gap_resolved_spliced` stays pass one's: the census classifies pass one's arbitration (:class:`DrainQC`).
            gap_resolution=dataclasses.replace(
                self.gap_resolution,
                gap_deferred_rna_or_gdna=0,
                gap_deferred_which_introns=0,
                gap_deferred_both=0,
            ),
            deferred=DeferredFragments.empty(),
            drain=drain,
        )

    @classmethod
    def from_scan_result(
        cls, scan_result: dict[str, Any], graph_hash: str | None = None
    ) -> "AccumulatorPayload":
        """The scan's ``result["calibration"]`` dict, validated: every bank's dtype, shape and axis, the partition's
        offsets against its region-bound axis, and the ledgers that must close (the start bank and the length
        histogram each book every deposited fragment once; the deferred bank holds exactly the fragments the held
        counters count)."""
        cal = scan_result["calibration"]
        n_refs = int(cal["n_refs"])
        n_cols = int(cal["n_strand_columns"])
        if n_cols != N_STRAND_COLUMNS:
            raise ValueError(
                f"the scan reports n_strand_columns={n_cols}, expected {N_STRAND_COLUMNS}. The two "
                f"columns are the two genome strands; any other count means the C++ and this schema "
                f"disagree about the channel axis."
            )
        max_length = int(cal["max_length"])
        region_bounds, offsets = _partition(cal, n_refs)
        rows_on = {
            "region": int(offsets["ref_region_offsets"][-1]),
            "boundary": int(offsets["ref_boundary_offsets"][-1]),
            "sj": int(offsets["ref_sj_offsets"][-1]),
        }
        banks: dict[str, np.ndarray] = {}
        for name, axis, dtype in BANK_AXES:
            banks[name] = _bank(cal, name, rows_on[axis], dtype)
        for name, axis, dtype in SINGLE_COLUMN_AXES:
            banks[name] = _single_column_bank(cal, name, rows_on[axis], dtype)

        qc = ScanQC.from_dict(cal["qc"])
        # The ledger: every deposited fragment books its first covered base in region_start_count, once.
        n_booked = int(banks["region_start_count"].sum())
        if n_booked != qc.deposited:
            raise ValueError(
                f"region_start_count sums to {n_booked} but {qc.deposited} fragments were deposited; "
                f"the ledger books every deposited fragment exactly once in region_start_count."
            )
        pool_lengths, deposited_lengths = _length_banks(cal, max_length, qc)
        deferred, gap_resolution = _held_fragments(cal, qc)
        return cls(
            region_bounds=region_bounds,
            **offsets,
            **banks,
            pool_lengths=pool_lengths,
            deposited_lengths=deposited_lengths,
            deferred=deferred,
            gap_resolution=gap_resolution,
            qc=qc,
            max_length=max_length,
            n_refs=n_refs,
            mismatches=int(cal["unique_unspliced_mismatches"]),
            aligned_bases=int(cal["unique_unspliced_bases"]),
            graph_hash=graph_hash,
        )


def _partition(cal: dict[str, Any], n_refs: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """The region-bound axis and the four per-reference offset arrays. The region and boundary offsets are re-derived
    from the region-bound axis rather than trusted: every consumer slices by them, and an offset array of the right
    length can still be inconsistent."""
    offsets = {
        name: _offsets(cal, name, n_refs)
        for name in (
            "ref_region_bound_offsets",
            "ref_region_offsets",
            "ref_boundary_offsets",
            "ref_sj_offsets",
        )
    }
    region_bounds = np.ascontiguousarray(cal["region_bounds"], dtype=np.int64)
    if int(offsets["ref_region_bound_offsets"][-1]) != region_bounds.shape[0]:
        raise ValueError(
            f"ref_region_bound_offsets ends at {int(offsets['ref_region_bound_offsets'][-1])} but region_bounds has "
            f"{region_bounds.shape[0]} entries"
        )
    per_ref_region_bounds = np.diff(offsets["ref_region_bound_offsets"])
    for name, per_ref in (
        ("ref_region_offsets", np.maximum(per_ref_region_bounds - 1, 0)),
        ("ref_boundary_offsets", np.maximum(per_ref_region_bounds - 2, 0)),
    ):
        expected = np.zeros(n_refs + 1, np.int64)
        np.cumsum(per_ref, out=expected[1:])
        if not np.array_equal(offsets[name], expected):
            bad = int(np.argmax(offsets[name] != expected))
            raise ValueError(
                f"{name}[{bad}] is {int(offsets[name][bad])} but the region_bound axis implies "
                f"{int(expected[bad])}. A reference contributing c region_bounds owns c-1 regions and c-2 "
                f"interior boundaries, and none at all below two region_bounds."
            )
    return region_bounds, offsets


def _length_banks(
    cal: dict[str, Any], max_length: int, qc: ScanQC
) -> tuple[np.ndarray, np.ndarray]:
    """The five pools' length histograms, ``[5, max_length + 1]``, and the unconditional one, which must bin every
    deposited fragment exactly once: it anchors every length law."""
    pool_lengths = _exact_dtype(cal, "pool_lengths", np.int64)
    if pool_lengths.size != N_FRAGMENT_POOLS * (max_length + 1):
        raise ValueError(
            f"pool_lengths has {pool_lengths.size} entries, expected "
            f"{N_FRAGMENT_POOLS} x (max_length + 1) = {N_FRAGMENT_POOLS * (max_length + 1)}"
        )
    deposited_lengths = _exact_dtype(cal, "deposited_lengths", np.uint32)
    if deposited_lengths.shape != (max_length + 1,):
        raise ValueError(
            f"deposited_lengths has shape {deposited_lengths.shape}, expected ({max_length + 1},)"
        )
    n_binned = int(deposited_lengths.sum())
    if n_binned != qc.deposited:
        raise ValueError(
            f"deposited_lengths sums to {n_binned} but {qc.deposited} fragments were "
            "deposited; the unconditional length histogram must bin every one of them exactly once."
        )
    return pool_lengths.reshape(N_FRAGMENT_POOLS, max_length + 1), deposited_lengths


def _held_fragments(cal: dict[str, Any], qc: ScanQC) -> tuple[DeferredFragments, GapCensus]:
    """The side buffer and the gap census. The bank must hold exactly the fragments the held counters count (a
    truncated cache would arrive short), and the census's deferred classes must sum to its counter."""
    deferred = DeferredFragments.from_dict(cal["deferred"])
    held = qc.deferred_undetermined_gap + qc.deferred_multiple_placements
    if deferred.n_fragments != held:
        raise ValueError(
            f"the deferred bank holds {deferred.n_fragments} fragments but "
            f"qc.deferred_undetermined_gap + qc.deferred_multiple_placements is {held}; the counters "
            f"and the fragments they count must be the same population, or the second pass silently "
            f"drains a different one."
        )
    multi = np.diff(deferred.placement_offsets) > 1
    if int(multi.sum()) != qc.deferred_multiple_placements:
        raise ValueError(
            f"the deferred bank holds {int(multi.sum())} multi-placement records but "
            f"qc.deferred_multiple_placements is {qc.deferred_multiple_placements}."
        )
    gap_resolution = GapCensus.from_dict(cal["gap_resolution"])
    if gap_resolution.deferred != qc.deferred_undetermined_gap:
        raise ValueError(
            f"the gap census's three deferred_* sum to {gap_resolution.deferred} but "
            f"qc.deferred_undetermined_gap is {qc.deferred_undetermined_gap}; the subclasses are "
            f"exhaustive by construction, so this is a partition that does not close."
        )
    return deferred, gap_resolution


def _offsets(cal: dict[str, Any], name: str, n_refs: int) -> np.ndarray:
    array = np.ascontiguousarray(cal[name], dtype=np.int64)
    if array.shape != (n_refs + 1,):
        raise ValueError(f"{name} has shape {array.shape}, expected ({n_refs + 1},)")
    return array


def _exact_dtype(cal: dict[str, Any], name: str, dtype: type) -> np.ndarray:
    """``cal[name]`` as a contiguous array, its dtype asserted rather than coerced.

    A coercion would widen or narrow silently, and either compares equal by value, so a schema drift would pass every
    value check. Without a dtype ``ascontiguousarray`` is a no-op on the C++ arrays, which stay views.
    """
    array = np.ascontiguousarray(cal[name])
    if array.dtype != dtype:
        raise ValueError(
            f"{name} has dtype {array.dtype}, expected {np.dtype(dtype)}. A widened or narrowed array "
            f"compares equal by value and would hide the change."
        )
    return array


def _bank(cal: dict[str, Any], name: str, rows: int, dtype: type) -> np.ndarray:
    """One two-column bank, reshaped from the flat array the C++ emits; its dtype is asserted by
    :func:`_exact_dtype`."""
    array = _exact_dtype(cal, name, dtype)
    if array.size != rows * N_STRAND_COLUMNS:
        raise ValueError(
            f"{name} has {array.size} entries, which does not divide into {rows} rows x "
            f"{N_STRAND_COLUMNS} strand columns"
        )
    return array.reshape(rows, N_STRAND_COLUMNS)


def _single_column_bank(cal: dict[str, Any], name: str, rows: int, dtype: type) -> np.ndarray:
    """One single-column bank, a density or a mass, with :func:`_bank`'s dtype contract. Kept separate because a
    two-column bank asserts the strand axis is meaningful for its channel and this one asserts it is not."""
    array = _exact_dtype(cal, name, dtype)
    if array.shape != (rows,):
        raise ValueError(
            f"{name} has shape {array.shape}, expected ({rows},). This bank has ONE value per object, "
            f"not one per strand column — a {(rows, N_STRAND_COLUMNS)} array here means the emitter "
            f"gave it a strand axis it does not have."
        )
    return array


__all__ = [
    "N_FRAGMENT_POOLS",
    "N_STRAND_COLUMNS",
    "POOL_DNA_INTERGENIC",
    "POOL_DNA_INTERGENIC_EXON",
    "POOL_DNA_INTRONIC",
    "POOL_DNA_INTRON_EXON",
    "POOL_RNA_SPLICED",
    "AccumulatorPayload",
    "DeferredFragments",
    "GapCensus",
    "ScanQC",
]
