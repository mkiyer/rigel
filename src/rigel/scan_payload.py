"""rigel.scan_payload — the BAM scan's accumulator tally, as a typed Python object.

    Spec: ``tests/native/_accumulator_reference.py``

The C++ scanner returns ``result["calibration"]`` as a dict of **flat** ndarrays
(``BamScanner::build_result``); this module reshapes the two-column banks and validates the schema.
Everything downstream reads only this object.

⭐ **THE FIELD NAMES ARE THE SPECIFICATION'S ``Tally`` FIELD NAMES, CHARACTER FOR CHARACTER**, and so are
the C++ payload keys. One quantity, one name, in all three places — so there is no mapping table anywhere
that could drift. `tests/test_accumulator_payload.py` checks the schema against ``Tally`` itself rather
than against a written-out list, for the same reason.

THE AXES — three of them, off by one from each other per reference::

    cuts    0        100       200       600        c = 4 cuts on this reference
    regions   [  n0  ][   n1   ][   n2   ]            c - 1 = 3 regions
    lines            line 1    line 2               c - 2 = 2 contiguous edges

A reference contributing ``c`` cuts owns ``c − 1`` regions and ``c − 2`` interior lines; one contributing
none owns neither, which is legal. Junction edges are their own axis, sliced by ``ref_sj_offsets``; the
flat slot order is the per-reference banks concatenated in reference order, which is what lets a
junction-edge id simply BE its slot.

WHAT THE NUMBERS MEAN. ⭐⭐ **ONE NUMERIC CONVENTION: a COUNT is an integer, a FRACTION is float64.**
There is no fixed point and no scale constant, so nothing anywhere decodes a bank::

    count           Sum 1                    integer   — exact, and reproduces across worker counts
    inv_length_sum  Sum 1/placements         float64   placements = L at a region, L−1 at a 0-bp line
    mass            Sum slice_len/(L·bounds) float64   — the conserved fragment count

⛔ **float64 is not a concession, it is the more accurate choice here** (measured 2026-08-11). Against
exact rational arithmetic on the reciprocal-opportunity theorem the fixed point it replaced missed the
answer by 7.0e-10 … 2.0e-07 while float64 misses by 5.8e-15 … 2.8e-13 — 1e5–7e5× closer.
⚠ **What it costs is bit-identity across worker counts**, since float addition is not associative: the
integer banks still reproduce exactly, the float ones agree to ~1e-15. Owner ruling 2026-08-11 — tests
validate the float banks within a DERIVED tolerance. `TRAPS: integer-channels-reproduce`.

⚠ **``inv_length_sum`` is NOT called ``density`` on purpose.** It is an exact, model-free density at an
edge — the opportunity ``L−1`` and the deposit ``1/(L−1)`` cancel identically — and it is *not* a density
at a region, where the opportunity is ``(region − L + 1)₊`` and nothing cancels. One word for two concepts is
the defect this naming avoids.

⛔⛔ **RETRACTED, AND DELETED WITH ITS BANKS (2026-08-13).** This paragraph read: *"``length_sum`` exists
because the other two are blind to one real case. At an edge the count row is ``(mu_g − 1, mu_r − 1)``
and the inv-length row is ``(1, 1)``, so the determinant is ``mu_g − mu_r``: when gDNA and RNA share a
mean fragment length the pair carries zero information about the split, at any depth. ``length_sum`` is
an independent tilt and removes that blind spot."*
⭐ **The premise is right and the conclusion does not follow.** At ``mu_g = mu_r = mu`` the third row is
``(mu, mu)`` — proportional to ``(1, 1)``, exactly like the other two — so the 3x2 system is still rank
one and ``length_sum`` removes nothing. It is an independent tilt only when the means already differ,
which is precisely when the first two rows are already sufficient. `TRAPS: equal-lengths-carry-no-composition`
says the same thing from the other side, and `ROADMAP.md` §1.4 closed the whole θ-independent-length
programme by measurement on 2026-08-10.
⚠ **It survived because nothing consumed it**: the banks reached ``PopulationView.length_sum`` and
stopped, so no test could disagree with the claim. A justification with no consumer is a claim with no
gate — `TRAPS: a-guard-outlives-its-divisor`, in prose form.

The trailing ``2`` on every bank is the **genome strand** — ``Strand.POS`` then ``Strand.NEG``, without
exception. Sense/antisense is transcript-relative, derived by the consumer from the junction's own
strand, and never stored. ⭐ ``sj_mass`` joined them on 2026-08-13 and is the only mass that has: see
``JunctionEdge::mass`` for the premise that changed and why it did not reach the other two masses.

⚠ **OWNERSHIP: this object holds VIEWS, and it is the keep-alive.** ``np.ascontiguousarray(x, dtype=D)`` is
a **no-op** when the array already has dtype ``D``, so nothing here copies — the buffers belong to
capsules owned by the C++ side. Someone "adding a cast for safety" would silently double peak memory
against a 1.04 M-region partition. Do not.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Any

import numpy as np


#: The two columns of every bank ARE the two genome strands. ⚠ `Strand` has four values — OR semantics
#: make POS|NEG == AMBIGUOUS, and NONE means no strand — so only two of them name a column, and a fragment
#: carrying neither is rejected by the accumulator rather than filed under one.
N_STRAND_COLUMNS = 2

#: Five fragment-length pools, each pure by construction (design §8). The order is
#: `rigel::accumulator::FragmentPool` and `_accumulator_reference.FragmentPool`.
N_FRAGMENT_POOLS = 5

#: The pool axis, named. ⚠ **These live here, with the schema, and not in the consumer that indexes by
#: them** — they are the accumulator's own enum in a third language, and a consumer holding a private
#: copy is how three files end up disagreeing about which row is which. A disagreement here would fit the
#: gDNA length model from the RNA pool and nothing downstream would look wrong;
#: `tests/calibration/test_fl.py` pins them against the executable specification's enum itself.
#:
#: Purity is the point (design §8): the two DNA_* contained pools are ~99 % gDNA on real data, and
#: RNA_SPLICED used an ANNOTATED junction with the splice OBSERVED — gDNA cannot be spliced. ⭐ The two
#: *_EXON "splash" pools are the only ON-TARGET gDNA population, so they are named rather than folded
#: into the gDNA model: on-target gDNA runs ~42 bp shorter than off-target (§8.2), and the shipped model
#: read a gDNA mean of 146.05 against the pure intergenic pool's 88.0 precisely by pooling them in.
#: There is deliberately NO pool for an exonic contained fragment or a multi-line crossing — those are
#: gDNA/RNA mixtures, and an impure pool is worse than a missing one.
POOL_DNA_INTERGENIC = 0  # contained in an intergenic region — pure gDNA
POOL_DNA_INTRONIC = 1  # contained in an intronic region — pure gDNA
POOL_DNA_INTRON_EXON = 2  # crossing one line, flanks {intron, exon} — on-target gDNA
POOL_DNA_INTERGENIC_EXON = 3  # crossing one line, {intergenic, exon} — on-target gDNA
POOL_RNA_SPLICED = 4  # used an annotated junction, splice OBSERVED — pure RNA


@dataclass(frozen=True, slots=True)
class ScanQC:
    """The denominators design §10.3 requires the accumulator to emit.

    Not optional and not derivable afterwards: **every conservation statement downstream must be able to
    name what it excluded.** Typed rather than a dict so that a misspelled denominator fails at this
    boundary instead of silently reading as zero somewhere far away.

    The field names are the specification's own ``Tally.qc`` keys.
    """

    deposited: int  # fragments that reached an object; == sum(region_start_count)
    dropped_too_long: int  # L above --max-fragment-length
    dropped_empty: int  # no path left after clipping to the reference
    dropped_strand_undefined: int  # align_strand named no column, so there was none to credit
    #: ⭐ >1 surviving hypothesis, so the fragment's gap is undetermined. ⚠ NOT dropped — it is held
    #: WHOLE for the second pass, and the identity is `deposited + deferred + dropped_* == offered`.
    deferred_undetermined_gap: int
    unannotated_introns: int  # observed introns with no annotated junction
    contradictory_sj_strand: int  # the mates' motif tags disagreed; no splice trusted
    introns_absorbed: int  # overlapping or abutting introns merged away

    @classmethod
    def from_dict(cls, qc: dict[str, Any]) -> "ScanQC":
        expected = {field.name for field in dataclasses.fields(cls)}
        missing = expected - set(qc)
        if missing:
            raise ValueError(
                f"the scan's qc block is missing {sorted(missing)}. Every one of these is a reported "
                f"denominator, so a missing key is a statement that cannot name what it excluded."
            )
        return cls(**{name: int(qc[name]) for name in expected})


@dataclass(frozen=True, slots=True)
class GapCensus:
    """⭐ The umbrella census: how each fragment whose gap needed resolving was resolved.

    Every fragment for which the enumeration produced at least one **spliced** hypothesis is counted here —
    its ``L`` depends on whether a gap intron is cut. The subclasses are exhaustive and mutually exclusive,
    and the three ``gap_deferred_*`` are exactly ``qc.deferred_undetermined_gap``.

    ⛔ **Its own axis, and NOT a splice type.** The umbrella cuts ACROSS the splice census: a certified-RNA
    ``SPLICED_ANNOT`` fragment with an intron in its mate gap needs resolving exactly as much as an
    ``UNSPLICED`` one does, so putting these on ``splice_type`` would need two labels per fragment and would
    break TRAPS: pure-and-length-censored's property that the splice census sums to the library.

    ⛔ **There is no ``gap_resolved_unspliced``, and that is not an omission.** A spliced hypothesis cuts
    bases the unspliced one keeps, so ``L_spliced <= L_unspliced`` always, and the one arbitration filter is
    ``L <= max_fragment_length``. If the unspliced path survives the filter then every spliced path does
    too — so the unspliced path can never be the sole survivor, and the class it would name is unreachable.
    The ordering is pinned by ``test_gap_hypothesis_arbitration.py``.

    The field names are the specification's own ``Tally.gap_resolution`` keys.
    """

    #: One hypothesis survived, and it necessarily cuts something: the gap intron is real and ``L``
    #: excludes it. ⚠ Classifies the ARBITRATION, not the deposit — such a fragment can still be rejected
    #: afterwards as ``TOO_LONG``, which is a different question with its own counter.
    gap_resolved_spliced: int
    #: ⛔ The unspliced path against exactly one spliced path. The open question is **RNA or gDNA** — one
    #: bit, and it is the composition question calibration exists to answer.
    gap_deferred_rna_or_gdna: int
    #: ⛔ Two or more spliced paths and no unspliced one: gDNA cannot be spliced, so the molecule is
    #: certified RNA and the open question is purely **which structure**.
    gap_deferred_which_introns: int
    #: ⛔ Both questions at once.
    gap_deferred_both: int

    @classmethod
    def zeros(cls) -> "GapCensus":
        """No fragment had a gap to resolve. ⚠ The ONE spelling of that, so a hand-built payload cannot
        invent a second — and it is a real state: a library with no annotated intron in any mate gap."""
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
    """⭐ **What the second pass's DRAIN did** — the audit trail of how the final tally was reached.

     The drain replays each held fragment with one chosen hypothesis, so
    after it **nothing is held**: the payload's bank is empty, ``qc.deferred_undetermined_gap`` is 0 and
    the three ``gap_deferred_*`` are 0. ⭐ That is deliberate — it keeps *"the counter and the fragments it
    counts are the same population"* absolute, so a drained payload needs no exception at the door. Pass
    one's numbers do not vanish; they live here.

    ⛔ **The drain has its OWN axis and does not extend `GapCensus`.** That census has no
    ``gap_resolved_unspliced`` class because pass-one arbitration cannot produce one — the genomic path is
    always the longest, so it can never be the sole survivor. ⭐ The drain, however, *chooses*, and it can
    choose ∅. ``chose_genomic`` / ``chose_spliced`` is the composition the drain assigned and is exactly
    what the census could not have recorded without resurrecting a class S1 deleted for a proven reason.
    """

    #: Pass one's ``qc.deferred_undetermined_gap`` — ⭐ the denominator the conservation is checked against.
    offered: int
    deposited: int
    dropped_too_long: int
    dropped_empty: int
    dropped_strand_undefined: int
    #: The drain's own census: ∅ against a spliced path. ⚠ Sums to ``offered``, not to ``deposited`` — a
    #: chosen hypothesis can still be rejected afterwards, which is a different question.
    chose_genomic: int
    chose_spliced: int
    # Pass one's arbitration census, kept verbatim so the per-class before/after
    #: is still readable off a drained payload.
    census_before: GapCensus

    @property
    def dropped(self) -> int:
        return self.dropped_too_long + self.dropped_empty + self.dropped_strand_undefined

    @property
    def conserved(self) -> bool:
        """⭐ §6.2's identity: every offered fragment either deposited or was rejected, exactly once."""
        return self.deposited + self.dropped == self.offered

    def __post_init__(self) -> None:
        if not self.conserved:
            raise ValueError(
                f"the drain offered {self.offered} fragments but accounts for "
                f"{self.deposited + self.dropped} (deposited {self.deposited} + dropped {self.dropped}); "
                f"a drained fragment that is neither deposited nor rejected has been lost."
            )
        if self.chose_genomic + self.chose_spliced != self.offered:
            raise ValueError(
                f"the drain chose {self.chose_genomic + self.chose_spliced} hypotheses for "
                f"{self.offered} offered fragments; exactly one hypothesis wins each whole fragment."
            )

    @classmethod
    def from_dict(cls, drain: dict[str, Any]) -> "DrainQC":
        expected = {field.name for field in dataclasses.fields(cls)} - {"census_before"}
        missing = expected - set(drain)
        if missing:
            raise ValueError(f"the drain block is missing {sorted(missing)}")
        census = drain["census_before"]
        return cls(
            **{name: int(drain[name]) for name in expected},
            census_before=census if isinstance(census, GapCensus) else GapCensus.from_dict(census),
        )


#: ⭐ Every two-column bank, with the axis it is indexed on and its dtype. ⚠ ONE table: `from_scan_result`
#: validates against it and the drain adds a per-reference delta into it, so a new channel cannot reach one
#: and miss the other.
BANK_AXES: tuple[tuple[str, str, Any], ...] = (
    ("region_contained_count", "region", np.uint32),
    ("edge_unspliced_count", "edge", np.uint32),
    ("edge_spliced_count", "edge", np.uint32),
    ("sj_count", "sj", np.uint32),
    # ⭐⭐ THE ONE MASS WITH A STRAND, and the only non-integer row in this table. `accumulator.h`'s
    # one-value ruling was reversed on the SJ axis alone (owner, 2026-08-12) because its premise
    # changed: an artifactual junction accumulates SYMMETRICALLY on both strands like gDNA, so the
    # existing strand model can detect one — given a per-strand observable. `JunctionEdge::mass` carries
    # the full reversal. ⛔ The columns are `sj_count`'s columns, so `mass[c]/count[c]` is a per-strand
    # mean. ⚠ The ruling still STANDS for `edge_unspliced_mass`, which has no such consumer.
    ("sj_mass", "sj", np.float64),
)


#: ⭐ The single-column additive banks, with the axis and dtype `from_scan_result` validates them on.
#: ⚠ Kept beside `BANK_AXES` rather than folded into it: that table's contract is "every TWO-column
#: bank", and widening it to mean "every bank" would silently give the mass a strand column in every
#: loop that reads it.
SINGLE_COLUMN_AXES: tuple[tuple[str, str, Any], ...] = (
    # ⭐⭐⭐ ONE NUMERIC CONVENTION: a COUNT is an integer, a FRACTION is float64. Every row here is a
    # fraction, so every row is float64; there is no fixed point and no scale constant to decode.
    # ⚠ The two `*_length_sum` banks were the integer exception and are GONE (2026-08-13) — see the
    # module docstring for why their stated justification did not survive measurement.
    ("region_contained_inv_opportunity_sum", "region", np.float64),
    ("edge_unspliced_inv_length_sum", "edge", np.float64),
    ("sj_inv_length_sum", "sj", np.float64),
    ("edge_unspliced_mass", "edge", np.float64),
    ("edge_spliced_mass", "edge", np.float64),
)

#: ⭐ **EVERY additive array channel, with the axis it is indexed on** — the two-column banks plus the
#: three that are not banks. ``"library"`` means the axis is library-wide rather than per reference.
#:
#: ⚠ Derived from `BANK_AXES` rather than restated, because the drain must add every additive channel and
#: miss none. `region_start_count` and `deposited_lengths` are the two externally-checkable invariants (each
#: sums to `qc.deposited`), so a drain that skipped either would be caught — but `pool_lengths` would just
#: go quietly short, and nothing downstream would look wrong.
ADDITIVE_AXES: tuple[tuple[str, str], ...] = (
    *((name, axis) for name, axis, _dtype in BANK_AXES),
    *((name, axis) for name, axis, _dtype in SINGLE_COLUMN_AXES),
    ("region_start_count", "region"),
    ("pool_lengths", "library"),
    ("deposited_lengths", "library"),
)


#: The deferred bank's arrays, in the order :meth:`Tally.deferred_arrays` emits them. ⚠ Read off the
#: dataclass below rather than written out twice; the tuple exists because the validation needs an order.
_DEFERRED_RECORD_FIELDS = ("ref", "start", "end", "align_strand", "sj_strand")


@dataclass(frozen=True, slots=True)
class DeferredFragments:
    """⭐ **THE SIDE BUFFER.** Fragments whose gap has more than one surviving explanation, held WHOLE.

    A fragment's unsequenced mate gap may hold no intron, one, or several, and *which* cannot be observed —
    the bases are not there. Deciding needs a fragment-length distribution that does not exist until the
    first pass is over, so the first pass does not guess: it enumerates, and when more than one hypothesis
    survives it holds the fragment here.

    ⭐ **The FRAGMENT is stored, never its consequences.** Object ids are large, derived, and would have to
    be kept consistent with the partition; the fragment is small and replays exactly. The second pass
    re-enters ``Accumulator::deposit`` with the chosen hypothesis, so there is one tally path and
    byte-identity with the specification is preserved for free.

    Two nested variable-length levels — fragments hold hypotheses, hypotheses hold introns — so there are
    two offset arrays at each level. Offsets are cumulative and always start at 0, so an empty bank is
    ``[0]`` and never ``[]``, and ``n_fragments`` is ``len(offsets) - 1``.

    ⛔ **THIS IS THE ONE BANK WHOSE ORDER IS OBSERVABLE.** Every other is a sum of integers and integer
    addition is associative, so a per-worker merge is exact whatever order the chunks arrived in. This is a
    list, so the C++ export **sorts it on the record's own content** before it crosses the ABI — otherwise
    the same BAM would give a different byte sequence at 1, 2, 4 and 8 workers with identical contents. Two
    records that tie on that key are identical records, so no tie-break is needed or possible.

    ⚠ Every array is ``int64``, including the two strand columns, which are ``int32`` everywhere else in the
    scanner. One dtype for one bank: the parity gate compares dtypes, and a narrowing conversion at the
    boundary compares equal by value.
    """

    ref: np.ndarray  # int64[n] — which reference; the drain replays onto THAT cut axis
    start: np.ndarray  # int64[n] — the CLIPPED extent, because that is what the drain must replay
    end: np.ndarray  # int64[n]
    align_strand: np.ndarray  # int64[n]
    sj_strand: np.ndarray  # int64[n] — the OBSERVED motif strand, if any
    observed_intron_offsets: np.ndarray  # int64[n + 1] — in PAIRS, into observed_introns
    observed_introns: (
        np.ndarray
    )  # int64[2 * n_observed] — flat (start, end); cut under EVERY hypothesis
    hypothesis_offsets: np.ndarray  # int64[n + 1] — into the per-hypothesis arrays below
    hypothesis_sj_strand: (
        np.ndarray
    )  # int64[n_hypotheses] — INFERRED; an observed motif always wins
    hypothesis_intron_offsets: np.ndarray  # int64[n_hypotheses + 1] — in PAIRS
    hypothesis_introns: np.ndarray  # int64[2 * n_implied] — the IMPLIED introns; empty ⇒ unspliced
    hypothesis_t_offsets: np.ndarray  # int64[n_hypotheses + 1]
    hypothesis_t: np.ndarray  # int64[n_supporting] — which transcripts imply this path

    @property
    def n_fragments(self) -> int:
        return int(self.hypothesis_offsets.shape[0]) - 1

    @property
    def n_hypotheses(self) -> int:
        return int(self.hypothesis_sj_strand.shape[0])

    @classmethod
    def empty(cls) -> "DeferredFragments":
        """Nothing was deferred. ⚠ Offsets are cumulative and start at 0, so every offset array is ``[0]``
        and never ``[]`` — spelled once here so a hand-built payload cannot get that boundary wrong."""
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

        ⚠ The two nested CSRs are re-derived rather than trusted, exactly as ``ref_region_offsets`` is. An
        offset array of the right LENGTH can still be inconsistent, and the second pass indexes every one
        of these — a truncated bank would score a fragment against another fragment's hypotheses, which is
        a wrong answer that looks entirely plausible.
        """
        names = [field.name for field in dataclasses.fields(cls)]
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
                raise ValueError(
                    f"deferred[{name!r}] has shape {array.shape}, expected one dimension"
                )
            arrays[name] = array

        n = int(arrays["hypothesis_offsets"].shape[0]) - 1
        if n < 0:
            raise ValueError(
                "deferred['hypothesis_offsets'] is empty; offsets are cumulative and start at 0, so an "
                "empty bank is [0] and never []"
            )
        for name in _DEFERRED_RECORD_FIELDS:
            if arrays[name].shape != (n,):
                raise ValueError(
                    f"deferred[{name!r}] has {arrays[name].shape[0]} entries but the offsets describe "
                    f"{n} fragments"
                )
        for offsets_name, values_name, stride in (
            ("observed_intron_offsets", "observed_introns", 2),
            ("hypothesis_offsets", "hypothesis_sj_strand", 1),
            ("hypothesis_intron_offsets", "hypothesis_introns", 2),
            ("hypothesis_t_offsets", "hypothesis_t", 1),
        ):
            offsets = arrays[offsets_name]
            if offsets.shape[0] == 0 or int(offsets[0]) != 0:
                raise ValueError(
                    f"deferred[{offsets_name!r}] must start at 0; offsets are cumulative and an empty "
                    f"level is [0]"
                )
            if np.any(np.diff(offsets) < 0):
                bad = int(np.argmax(np.diff(offsets) < 0))
                raise ValueError(
                    f"deferred[{offsets_name!r}] decreases at {bad}; a CSR offset array cannot go "
                    f"backwards"
                )
            if int(offsets[-1]) * stride != arrays[values_name].shape[0]:
                raise ValueError(
                    f"deferred[{offsets_name!r}] ends at {int(offsets[-1])} but "
                    f"{values_name} has {arrays[values_name].shape[0]} entries "
                    f"({int(offsets[-1])} x {stride} expected)"
                )
        n_hypotheses = int(arrays["hypothesis_offsets"][-1])
        for name in ("hypothesis_intron_offsets", "hypothesis_t_offsets"):
            if arrays[name].shape[0] != n_hypotheses + 1:
                raise ValueError(
                    f"deferred[{name!r}] has {arrays[name].shape[0]} entries but there are "
                    f"{n_hypotheses} hypotheses, so it must have {n_hypotheses + 1}"
                )
        # ⭐ A fragment is deferred BECAUSE two or more hypotheses survived, so a record carrying fewer
        # than two is a bank that lost them — and the second pass would then "choose" from a set of one and
        # deposit an answer nothing supported.
        runs = np.diff(arrays["hypothesis_offsets"])
        if n and int(runs.min()) < 2:
            bad = int(np.argmin(runs))
            raise ValueError(
                f"deferred fragment {bad} carries {int(runs[bad])} hypotheses. A fragment is deferred "
                f"only when two or more survived, so every record must hold at least two."
            )
        return cls(**arrays)

    def observed_introns_of(self, i: int) -> np.ndarray:
        """Fragment ``i``'s observed introns as an ``[k, 2]`` view. Cut under **every** hypothesis."""
        lo, hi = int(self.observed_intron_offsets[i]), int(self.observed_intron_offsets[i + 1])
        return self.observed_introns[2 * lo : 2 * hi].reshape(hi - lo, 2)

    def hypothesis_introns_of(self, h: int) -> np.ndarray:
        """Hypothesis ``h``'s IMPLIED introns as an ``[k, 2]`` view. Empty ⇒ the unspliced hypothesis."""
        lo, hi = int(self.hypothesis_intron_offsets[h]), int(self.hypothesis_intron_offsets[h + 1])
        return self.hypothesis_introns[2 * lo : 2 * hi].reshape(hi - lo, 2)

    def supporting_t_of(self, h: int) -> np.ndarray:
        """Which candidate transcripts imply hypothesis ``h``. Empty for the unspliced one."""
        lo, hi = int(self.hypothesis_t_offsets[h]), int(self.hypothesis_t_offsets[h + 1])
        return self.hypothesis_t[lo:hi]


@dataclass(frozen=True, slots=True)
class AccumulatorPayload:
    """One BAM scan's tally. Views over C++-owned buffers; this object is the keep-alive."""

    # -- the partition, echoed back so a consumer can locate every object without reloading the index --
    cut_positions: np.ndarray  # int64[n_cuts] — flat, reference-major, ascending within a reference
    ref_cut_offsets: np.ndarray  # int64[n_refs + 1] — CSR over cut_positions
    ref_region_offsets: np.ndarray  # int64[n_refs + 1]
    ref_edge_offsets: np.ndarray  # int64[n_refs + 1] — contiguous edges
    ref_sj_offsets: np.ndarray  # int64[n_refs + 1] — junction edges

    # -- regions: two disjoint populations, each two genome-strand columns --
    region_contained_count: np.ndarray  # uint32[n_regions, 2] — the whole path lies inside the region
    #: ⭐ uint64[n_regions] — ONE column. The length moments are strand-AGNOSTIC: which strand a read
    #: aligned to says nothing about whether the molecule was gDNA or RNA, and every consumer summed
    #: the two columns. ⛔ The COUNTS keep both — the strand model is a Beta-Binomial over them.
    region_contained_inv_opportunity_sum: np.ndarray
    region_start_count: np.ndarray  # uint32[n_regions] — THE invariant; sums to qc.deposited

    # -- contiguous edges: the 0-bp line between two adjacent regions --
    edge_unspliced_count: np.ndarray  # uint32[n_edges, 2] — the mixture being deconvolved
    edge_unspliced_inv_length_sum: np.ndarray  # uint64[n_edges] — ONE column, strand-agnostic
    #: ⭐⭐ uint64[n_edges] — **THE CONSERVED MASS**, fixed point at ``INV_LENGTH_SCALE``. A COUNT and a
    #: MASS are two different deposits and one number cannot be both: ``edge_unspliced_count`` is ``+1``
    #: on every line a fragment crosses, so a fragment books ``max(K, 1)`` of them; this sums to ONE per
    #: fragment. That is what lets a consumer turn an object-incidence total into a FRAGMENT COUNT
    #: without manufacturing one from a density. ⛔ ONE column, not two — nothing reads a mass per
    #: strand, and the question it answers has no strand in it.
    edge_unspliced_mass: np.ndarray
    #: uint32[n_edges, 2] — certified RNA: gDNA cannot be spliced. ⭐ COUNT AND MASS ONLY: nothing
    #: deconvolves a certified-RNA crossing, so its two length moments had no consumer and are gone.
    edge_spliced_count: np.ndarray
    #: ⭐ uint64[n_edges] — the same rule, routed by the same ``spliced`` flag, so ``mass`` is not the
    #: one channel that ignores the split. ⛔ A PARTIAL, never a conservation ledger: it sums to
    #: ``crossed_block_len / L`` per fragment. A per-LINE certified-RNA term, commensurate with the
    #: unspliced mass at the same line — NOT "the number of spliced fragments here".
    edge_spliced_mass: np.ndarray

    # -- junction edges: one exact donor->acceptor jump. Pure RNA by construction --
    #: uint32[n_sj, 2] — ⛔⛔ **BOTH GENOME-STRAND COLUMNS ARE RETAINED, AND NOT BECAUSE ANYTHING READS
    #: THEM YET** (owner ruling, 2026-08-08). A junction is stranded by its genomic splicing MOTIF, so
    #: the strand of the *fragments* on it looks redundant, and every consumer today sums the two.
    #:
    #: ⭐ **The reason is aligner-artifact detection.** Aligners emit false-positive ``N`` CIGAR ops from
    #: plain genomic DNA. ``rigel.splice_blacklist`` catches those the sister tool ``alignable`` has
    #: enumerated by coordinate — an a-priori list, and far from complete. The EMPIRICAL detector is
    #: this column: in a stranded library a real junction inherits the global strand specificity, while
    #: an artifact deposits on BOTH strands and deviates from it. ⚠ Unstranded data cannot use it
    #: (κ = ½ leaves nothing to deviate from), which is a property of the detector, not a reason to drop
    #: the column.
    #:
    #: ⛔ The discriminating information lives ONLY in the split — a clean junction and an artifactual
    #: one carry the same total. Gated by
    #: ``test_the_junction_STRAND_SPLIT_IS_RETAINED_FOR_ALIGNER_ARTIFACT_DETECTION``.
    sj_count: np.ndarray
    #: uint64[n_sj] — ⭐ LIVE in ``second_pass``, which scores a held fragment's junction evidence
    #: with it. ⚠ ``sj_length_sum`` is gone for the same reason the spliced edge moments are.
    sj_inv_length_sum: np.ndarray
    #: uint64[n_sj] — ⭐⭐⭐ **THE CONSERVED MASS'S THIRD AXIS, and what makes a LIBRARY FRAGMENT COUNT
    #: COMPUTABLE.** A spliced fragment's block containing no interior line deposits on neither edge
    #: bank, and is not ``contained`` either — its path spans a junction, so it lies in no single region.
    #: Such a fragment existed on the incidence axis (``sj_count``) and on no conserved one.
    #:
    #: ⭐ Measured on the origin-split oracle at ladder g50 capture_off: **1,222,375 of 4,830,713 RNA
    #: fragments (25.3 %)** are in that population, against **0 of 4,997,761** gDNA fragments — gDNA
    #: cannot splice, so its conserved count was already exact while RNA's read 0.747x deposited.
    #:
    #: ⛔ **It ADDS a boundary class rather than re-apportioning one**: a block that crossed a line is
    #: untouched, so ``edge_unspliced_mass`` and ``edge_spliced_mass`` are byte-identical to what they
    #: were. Gates: ``tests/native/test_conserved_mass.py`` claim 5.
    sj_mass: np.ndarray

    # -- the fragment-length pools, binned at L, once per fragment --
    pool_lengths: np.ndarray  # int64[5, max_length + 1]

    #: uint32[max_length + 1] — ⭐ **TRAPS: a-purity-filter-is-a-length-filter: EVERY deposited fragment, binned at its own L, no purity
    #: condition.** The five pools above are deliberately CONDITIONED (an impure pool is worse than a
    #: missing one), so none of them is an unconditional anchor — which is why the empirical-Bayes
    #: shrinkage in `calibration.fl` took its anchor from the SCANNER, which measures length by two
    # other rules over another population. This is that anchor, in
    #: the accumulator's own frame.
    #: ⚠ "Unconditional GIVEN DEPOSIT": it excludes what the accumulator rejects (too long, ambiguous
    #: path, strand-undefined, empty), each counted in ``qc``. That is exactly the population the pools
    #: are drawn from, which is what makes it the right anchor rather than merely a convenient one.
    deposited_lengths: np.ndarray

    #: ⭐ **The side buffer** — fragments whose gap has more than one surviving explanation, held WHOLE for
    #: the second pass. ⚠ NOT a loss: ``deposited + deferred + dropped_* == offered``, and this bank holds
    #: the fragments ``qc.deferred_undetermined_gap`` counts.
    deferred: DeferredFragments
    #: ⭐ How each gap was resolved — its own axis, cutting across the splice census.
    gap_resolution: GapCensus

    qc: ScanQC
    max_length: int  # the fragment-length limit applied to L, and the pool-histogram width
    n_refs: int

    #: ⭐ Provenance, and it must cover **regions AND edges**. The payload is edge-keyed by construction —
    #: its junction axis is meaningless against a different junction CSR — and `index.partition_hash`
    #: deliberately covers `regions.feather` only. A 2026-07-29 flag fix rewrote every `edges.feather` while
    #: leaving every `regions.feather` byte-identical, so a regions-only key would have verified CLEAN against
    #: a stale cache. `None` when the scanner was driven without an index to hash.
    graph_hash: str | None = None

    #: ⭐ **Set once the second pass has DRAINED this tally**, and ``None`` while the side buffer is still
    # held. ⚠ It is the only way to tell the two states apart, because a
    #: drained payload is deliberately indistinguishable in shape: its bank is empty and its held counters
    #: are 0 precisely so that *"the counter and the fragments it counts are the same population"* needs no
    #: exception. Pass one's numbers live in here — see :class:`DrainQC`.
    drain: DrainQC | None = None

    # -- derived, never stored ------------------------------------------------------------------------

    @property
    def n_strand_columns(self) -> int:
        return N_STRAND_COLUMNS

    @property
    def n_regions(self) -> int:
        return int(self.region_start_count.shape[0])

    @property
    def n_edges(self) -> int:
        return int(self.edge_unspliced_count.shape[0])

    @property
    def n_sj(self) -> int:
        return int(self.sj_count.shape[0])

    def with_drain(self, delta: dict[str, np.ndarray], drain: DrainQC) -> "AccumulatorPayload":
        """⭐ **Pass one plus the drain's delta** — the tally calibration actually consumes.

         ``delta`` holds one globally-shaped array per additive channel, as
        produced by replaying the held fragments through :meth:`Accumulator.deposit`; this method is only
        the arithmetic and the bookkeeping.

        ⭐ **The delta arrives as its own object rather than being accumulated in place**, which is what
        makes the drain's contribution *observable*: both payloads exist, so every channel's before/after
        is a subtraction rather than a rerun — which is the reason to prefer
        this shape, and it is also why the drain needs no new C++ — every channel is already exported.

        After this: the bank is empty, the held counters are 0, and ``drain`` says what pass one held.
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
            # ⛔ Both terms are non-negative integers, so the sum can only fail to be >= either term by
            # WRAPPING. Counts are uint32 and a real library can approach that; a silent wrap would read
            # as a plausible small number in the one bank every density is computed from.
            if np.any(total < before) or np.any(total < added):
                raise ValueError(
                    f"adding the drain's delta to {name!r} overflowed {before.dtype}; the tally cannot "
                    f"represent pass one and the drain together."
                )
            totals[name] = total

        return dataclasses.replace(
            self,
            **totals,
            # ⚠ `deferred_undetermined_gap` goes to 0 with the bank, not kept at pass one's value: the two
            # must describe one population, and pass one's count is preserved as `drain.offered`.
            qc=dataclasses.replace(
                self.qc,
                deposited=self.qc.deposited + drain.deposited,
                dropped_too_long=self.qc.dropped_too_long + drain.dropped_too_long,
                dropped_empty=self.qc.dropped_empty + drain.dropped_empty,
                dropped_strand_undefined=(
                    self.qc.dropped_strand_undefined + drain.dropped_strand_undefined
                ),
                deferred_undetermined_gap=0,
            ),
            # ⛔ `gap_resolved_spliced` is pass one's and is NOT extended. The census classifies pass one's
            # ARBITRATION, and it has no class for a chosen ∅ — see :class:`DrainQC`.
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
        cal = scan_result.get("calibration")
        if cal is None:
            raise ValueError(
                "scan_result['calibration'] is None; BamScanner.set_regions was not called"
            )

        n_refs = int(cal["n_refs"])
        n_cols = int(cal["n_strand_columns"])
        if n_cols != N_STRAND_COLUMNS:
            raise ValueError(
                f"the scan reports n_strand_columns={n_cols}, expected {N_STRAND_COLUMNS}. The two "
                f"columns are the two genome strands; any other count means the C++ and this schema "
                f"disagree about the channel axis."
            )
        max_length = int(cal["max_length"])

        offsets = {
            name: _offsets(cal, name, n_refs)
            for name in (
                "ref_cut_offsets",
                "ref_region_offsets",
                "ref_edge_offsets",
                "ref_sj_offsets",
            )
        }
        cut_positions = np.ascontiguousarray(cal["cut_positions"], dtype=np.int64)
        if int(offsets["ref_cut_offsets"][-1]) != cut_positions.shape[0]:
            raise ValueError(
                f"ref_cut_offsets ends at {int(offsets['ref_cut_offsets'][-1])} but cut_positions has "
                f"{cut_positions.shape[0]} entries"
            )

        # ⚠ Re-derived from the cut axis rather than trusted. An offset array of the right LENGTH can
        # still be inconsistent, and every consumer slices by these — a per-reference offset that drifts
        # is the defect class that once dropped 476,719 of 476,732 fragments while every golden passed.
        per_ref_cuts = np.diff(offsets["ref_cut_offsets"])
        for name, per_ref in (
            ("ref_region_offsets", np.maximum(per_ref_cuts - 1, 0)),
            ("ref_edge_offsets", np.maximum(per_ref_cuts - 2, 0)),
        ):
            expected = np.zeros(n_refs + 1, np.int64)
            np.cumsum(np.where(per_ref_cuts > 0, per_ref, 0), out=expected[1:])
            if not np.array_equal(offsets[name], expected):
                bad = int(np.argmax(offsets[name] != expected))
                raise ValueError(
                    f"{name}[{bad}] is {int(offsets[name][bad])} but the cut axis implies "
                    f"{int(expected[bad])}. A reference contributing c cuts owns c-1 regions and c-2 "
                    f"interior lines, and none at all below two cuts."
                )

        n_regions = int(offsets["ref_region_offsets"][-1])
        n_edges = int(offsets["ref_edge_offsets"][-1])
        n_sj = int(offsets["ref_sj_offsets"][-1])

        rows_on = {"region": n_regions, "edge": n_edges, "sj": n_sj}
        banks: dict[str, np.ndarray] = {}
        for name, axis, dtype in BANK_AXES:
            banks[name] = _bank(cal, name, rows_on[axis], dtype)
        for name, axis, dtype in SINGLE_COLUMN_AXES:
            banks[name] = _single_column_bank(cal, name, rows_on[axis], dtype)

        region_start_count = np.ascontiguousarray(cal["region_start_count"], dtype=np.uint32)
        if region_start_count.shape != (n_regions,):
            raise ValueError(
                f"region_start_count has shape {region_start_count.shape}, expected ({n_regions},)"
            )
        pool_lengths = np.ascontiguousarray(cal["pool_lengths"], dtype=np.int64)
        if pool_lengths.size != N_FRAGMENT_POOLS * (max_length + 1):
            raise ValueError(
                f"pool_lengths has {pool_lengths.size} entries, expected "
                f"{N_FRAGMENT_POOLS} x (max_length + 1) = {N_FRAGMENT_POOLS * (max_length + 1)}"
            )

        deposited_lengths = np.ascontiguousarray(cal["deposited_lengths"], dtype=np.uint32)
        if deposited_lengths.shape != (max_length + 1,):
            raise ValueError(
                f"deposited_lengths has shape {deposited_lengths.shape}, expected ({max_length + 1},)"
            )
        # ⭐ THE TRAPS: a-purity-filter-is-a-length-filter INVARIANT, checked at the door. Same externally-checkable form as
        # ``sum(region_start_count) == deposited`` and a DIFFERENT statement: that one says every fragment
        # was located in space, this one that every fragment was binned by length. A histogram that is
        # the anchor for every FL model in the tool must not be allowed in one fragment short.
        n_binned = int(deposited_lengths.sum())
        if n_binned != int(cal["qc"]["deposited"]):
            raise ValueError(
                f"deposited_lengths sums to {n_binned} but {int(cal['qc']['deposited'])} fragments were "
                "deposited; the unconditional length histogram must bin every one of them exactly once."
            )

        qc = ScanQC.from_dict(cal["qc"])
        deferred = DeferredFragments.from_dict(cal["deferred"])
        # ⭐ THE CONSERVATION HALF, refused at the door. `qc.deferred_undetermined_gap` says how many
        # fragments were held; this bank is supposed to BE them. The identity
        # `deposited + deferred + dropped_* == offered` is worth nothing if the deferred term is a number
        # with no fragments behind it — and a cache can be truncated or partially written, which is
        # precisely how a bank would arrive short of the counter that describes it.
        if deferred.n_fragments != qc.deferred_undetermined_gap:
            raise ValueError(
                f"the deferred bank holds {deferred.n_fragments} fragments but "
                f"qc.deferred_undetermined_gap is {qc.deferred_undetermined_gap}; the counter and the "
                f"fragments it counts must be the same population, or the second pass silently drains a "
                f"different one."
            )
        gap_resolution = GapCensus.from_dict(cal["gap_resolution"])
        if gap_resolution.deferred != qc.deferred_undetermined_gap:
            raise ValueError(
                f"the gap census's three deferred_* sum to {gap_resolution.deferred} but "
                f"qc.deferred_undetermined_gap is {qc.deferred_undetermined_gap}; the subclasses are "
                f"exhaustive by construction, so this is a partition that does not close."
            )
        return cls(
            cut_positions=cut_positions,
            **offsets,
            **banks,
            region_start_count=region_start_count,
            pool_lengths=pool_lengths.reshape(N_FRAGMENT_POOLS, max_length + 1),
            deposited_lengths=deposited_lengths,
            deferred=deferred,
            gap_resolution=gap_resolution,
            qc=qc,
            max_length=max_length,
            n_refs=n_refs,
            graph_hash=graph_hash,
        )


def _offsets(cal: dict[str, Any], name: str, n_refs: int) -> np.ndarray:
    array = np.ascontiguousarray(cal[name], dtype=np.int64)
    if array.shape != (n_refs + 1,):
        raise ValueError(f"{name} has shape {array.shape}, expected ({n_refs + 1},)")
    return array


def _bank(cal: dict[str, Any], name: str, rows: int, dtype: type) -> np.ndarray:
    """One two-column bank, reshaped from the flat array the C++ emits.

    ⚠ The dtype is asserted rather than coerced. ``ascontiguousarray`` would happily widen a uint32 count
    to uint64, which compares equal by value — so a schema drift would pass every value check downstream
    while the payload silently stopped matching the specification.
    """
    array = np.ascontiguousarray(cal[name])
    if array.dtype != dtype:
        raise ValueError(
            f"{name} has dtype {array.dtype}, expected {np.dtype(dtype)}. Counts are uint32 and "
            f"densities uint64; a widened array compares equal by value and would hide the change."
        )
    if array.size != rows * N_STRAND_COLUMNS:
        raise ValueError(
            f"{name} has {array.size} entries, which does not divide into {rows} rows x "
            f"{N_STRAND_COLUMNS} strand columns"
        )
    return array.reshape(rows, N_STRAND_COLUMNS)


def _single_column_bank(cal: dict[str, Any], name: str, rows: int, dtype: type) -> np.ndarray:
    """One SINGLE-column bank — the conserved masses. Same dtype contract as :func:`_bank`.

    ⚠ Separate from :func:`_bank` rather than parameterised by a column count, because the two carry
    different claims: a two-column bank asserts the strand axis is meaningful for that channel, and this
    one asserts it is not. A shared function taking ``n_columns=1`` would read as the same statement.
    """
    array = np.ascontiguousarray(cal[name])
    if array.dtype != dtype:
        raise ValueError(
            f"{name} has dtype {array.dtype}, expected {np.dtype(dtype)}. These are uint64 sums (a "
            f"fixed point for the reciprocal ones), and a narrowed array wraps at realistic depth."
        )
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
