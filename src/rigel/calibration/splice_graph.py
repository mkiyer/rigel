"""rigel.calibration.splice_graph — the v8 SPLICE GRAPH: regions + boundaries.

The calibration partition: what the accumulator deposits into, and the structure the solver reads.


**A region** is a genomic interval; regions tile each reference and are numbered in genomic order.
**An boundary** is a transition, always ``src < dst`` so genomic order is a topological order:

* ``CONTIGUOUS(i, i+1)`` — the genomic point ``end(i) == start(i+1)``. Carries the 8 structural
  :ref:`flag bits <flags>`.
* ``JUNCTION(donor, acceptor, strand)`` — one per distinct annotated intron.

**Every distinct exon endpoint is a region_bound, and adjacent equal-signature regions are NOT merged.** The merge is
what the predecessor partition did and what made it blind to transcript termini: ⭐ **53.4 %** of real human
transcript termini (232,451 of 435,291) fall strictly *inside* a merged region and vanish from the partition
entirely. Measured consequence: 752,654 merged regions → **1,043,881** regions (+38.7 %), median 151 bp,
15,687 of length 1, plus **404,168** junction boundaries.

⚠ The often-quoted **59.5 %** is that same statistic computed under the ``~is_synthetic & ~is_nrna``
transcript filter — i.e. with 26,475 real single-exon transcripts excluded, which is the very bug
described below. Under the shipped filter it is **53.4 %**.

⭐ **ONE transcript filter: ``~is_synthetic``.** Every manufactured nRNA span row is ``is_synthetic``,
so that single predicate already excludes all of them — which is what P1G_SCOPE §5's "211/211 signature
false positives are nRNA-span boundaries" actually requires.

⚠ **It was briefly TWO filters, and the second one was a measured bug** (plan TRAPS: specificity-and-sense-are-complements proposed
``~is_synthetic & ~is_nrna`` for the flags and reaches, reasoning that "an nRNA span's ends are not real
transcript termini"). The reasoning is right; the predicate is not. On a **non-synthetic** row
``is_nrna`` does not mean "manufactured span" — it means **this real transcript is single-exon, so its
mature and nascent forms coincide**. Measured on the human annotation: all **26,475** rows that are
``is_nrna & ~is_synthetic`` have ``n_exons == 1`` and **none** is a ``RIGEL_NRNA_*`` row, so the extra
clause deleted the TSS/TES of 26,475 real transcripts — **52,104 distinct terminus positions** — which
is exactly the visibility v8 exists to buy.

.. _flags:

**Structural flags** (contiguous boundaries): ``TSS_s``/``TES_s``/``DONOR_s``/``ACCEPTOR_s`` for ``s ∈ {+,−}``.
They are NOT mutually exclusive — a position that is both a terminus for one transcript and a splice site
for another is exactly the case the 4-bit signature is structurally blind to.

**Reaches** — how many bases of the molecule's own sequence remain either side of an boundary. An RNA
molecule must fit inside its transcript, so this is what makes the crossing divisor taper near a
terminus (measured worth 11.0 % of the mature opportunity genome-wide, and ⭐ `Σ1/(L−1)` reads only
0.108 ρ twenty bases from a transcript end).

⚠ **The two boundary kinds carry DIFFERENT reaches, deliberately.**

* A **CONTIGUOUS** boundary is crossed by gDNA, by nascent RNA and by mature RNA alike. Nascent RNA is an
  ordinary transcript that happens to be single-exon and to span its whole gene, and a genomic distance
  is never shorter than the exonic distance within the same span — so the widest RNA molecule covering
  a position is the nascent one, and the reach is the **genomic** distance to a covering transcript's
  span ends. ⚠ It is therefore **nonzero inside an intron**: that is where nascent RNA lives, and an
  exonic reach would declare zero RNA opportunity across every intron in the genome.
* A **JUNCTION** boundary is used only by a molecule that spliced across it, so what remains either side is
  **exonic**. :func:`_junction_edges` is deliberately left on the exonic reach.

⚠ **Named ``lo``/``hi``, NOT ``donor``/``acceptor``.** Boundaries store ``src < dst``, so ``src`` is genomically
LEFT whatever the strand — but for a NEG-strand junction the biological donor is on the RIGHT. The
divisor is symmetric in the two reaches so no number changes, but ``reach_left``/``reach_right`` name the two sides genomically; a
``reach_donor`` would mislabel roughly
half of the 404,168 junctions. ``lo``/``hi`` means genomically lower/higher and cannot be misread.

⚠ **Per STRAND**, because the mature-crossing gate is per strand: at an AMBIG boundary a ``+`` and a ``−``
transcript have different reaches, and a strand-agnostic maximum would over-state one of them — on exactly
the population that carries 31.6 % of the calibration suite's error mass.

**Maximal over isoforms, independently per side** (decision TRAPS: two-gaussians-one-latent). A genomic position usually belongs to
several transcripts that disagree about where the molecule ends, and the isoform abundances are precisely
what calibration does not know. Measured (plan TRAPS: divide-by-a-probability): against the alternative of taking the best *realisable*
pair from a single isoform, the two agree exactly on **93.9 %** of disagreeing junctions with a mean ratio
of **0.9989**, so the simple independent maximum is used. It is one-sided — it over-states opportunity,
hence under-states ``ρ_mature``, hence over-states ``f_g``.

**A reach of 0 is meaningful, not a sentinel:** no strand-``s`` RNA molecule can occupy that side of the
boundary, so the opportunity is genuinely zero and the trapezoid returns 0 for free. On a contiguous boundary
that now means *outside every transcript span on that strand*; on a junction boundary, *no strand-``s``
molecule splices across it*.
"""

from __future__ import annotations

import dataclasses
import warnings
from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd

from ..transcript import Transcript
from ..types import Strand
from .signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    N_SIGNATURES,
    coarse_type_array,
)

__all__ = [
    "REGION_COLUMNS",
    "REGION_COLUMN_DTYPES",
    "BOUNDARY_COLUMNS",
    "BOUNDARY_COLUMN_DTYPES",
    "EDGE_KIND_CONTIGUOUS",
    "EDGE_KIND_JUNCTION",
    "FLAG_TSS_POS",
    "FLAG_TSS_NEG",
    "FLAG_TES_POS",
    "FLAG_TES_NEG",
    "FLAG_DONOR_POS",
    "FLAG_DONOR_NEG",
    "FLAG_ACCEPTOR_POS",
    "FLAG_ACCEPTOR_NEG",
    "build_splice_graph",
    "build_region_partition_arrays",
    "build_boundary_flags_array",
    "build_junction_edge_arrays",
    "build_contiguous_boundary_reach_arrays",
    "build_junction_geometry_arrays",
    "build_transcript_path",
    "JunctionEdgeArrays",
    "JunctionGeometry",
    "TranscriptPath",
    "STEP_REGION",
    "STEP_BOUNDARY",
    "STEP_SPLICE_JUNCTION",
    "is_terminus",
    "is_splice_site",
    "validate_graph",
    "load_regions",
    "load_boundaries",
]

REGION_COLUMNS = ["node_id", "ref_name", "start", "end", "length", "signature"]
REGION_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "node_id": np.int64,
    "start": np.int64,
    "end": np.int64,
    "length": np.int64,
    "signature": np.uint8,
}

EMPTY_I64 = np.zeros(0, np.int64)

EDGE_KIND_CONTIGUOUS = 0
EDGE_KIND_JUNCTION = 1

BOUNDARY_COLUMNS = [
    "edge_id",
    "src",
    "dst",
    "kind",
    "strand",
    "flags",
    "reach_lo_pos",
    "reach_hi_pos",
    "reach_lo_neg",
    "reach_hi_neg",
]
BOUNDARY_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "edge_id": np.int64,
    "src": np.int64,
    "dst": np.int64,
    "kind": np.uint8,
    "strand": np.int8,
    "flags": np.uint16,
    "reach_lo_pos": np.int32,
    "reach_hi_pos": np.int32,
    "reach_lo_neg": np.int32,
    "reach_hi_neg": np.int32,
}

# Structural bits on a CONTIGUOUS boundary, at its genomic position. Not mutually exclusive.
FLAG_TSS_POS = np.uint16(1 << 0)
FLAG_TSS_NEG = np.uint16(1 << 1)
FLAG_TES_POS = np.uint16(1 << 2)
FLAG_TES_NEG = np.uint16(1 << 3)
FLAG_DONOR_POS = np.uint16(1 << 4)
FLAG_DONOR_NEG = np.uint16(1 << 5)
FLAG_ACCEPTOR_POS = np.uint16(1 << 6)
FLAG_ACCEPTOR_NEG = np.uint16(1 << 7)


# ---------------------------------------------------------------------------
# Flattening: list[Transcript] -> the exon arrays the vectorised build needs
# ---------------------------------------------------------------------------


class _Exons:
    """Flat per-exon arrays for one filter, plus per-exon cumulative exonic offsets.

    ``before[i]`` is the transcript's exonic length strictly before exon ``i``; ``total[i]`` its full
    exonic length. Together they give the reach either side of any position inside that exon without a
    per-transcript loop.
    """

    __slots__ = ("ref", "start", "end", "strand", "tid", "before", "total", "n")

    def __init__(self, transcripts, keep):
        ref, start, end, strand, tid = [], [], [], [], []
        for k, tx in enumerate(transcripts):
            if not keep(tx):
                continue
            for ex in tx.exons:
                if int(ex.end) <= int(ex.start):
                    continue  # v7 skips a zero/negative-length exon silently; match it exactly
                ref.append(str(tx.ref))
                start.append(int(ex.start))
                end.append(int(ex.end))
                strand.append(int(tx.strand))
                tid.append(k)
        self.n = len(start)
        self.ref = np.asarray(ref, dtype=object)
        self.start = np.asarray(start, dtype=np.int64)
        self.end = np.asarray(end, dtype=np.int64)
        self.strand = np.asarray(strand, dtype=np.int8)
        self.tid = np.asarray(tid, dtype=np.int64)
        if self.n == 0:
            self.before = np.zeros(0, np.int64)
            self.total = np.zeros(0, np.int64)
            return
        # exons of one transcript are contiguous in this array and already in genomic order
        order = np.lexsort((self.start, self.tid))
        for a in ("ref", "start", "end", "strand", "tid"):
            setattr(self, a, getattr(self, a)[order])
        length = self.end - self.start
        csum = np.cumsum(length)
        first = np.flatnonzero(np.r_[True, self.tid[1:] != self.tid[:-1]])
        last = np.r_[first[1:], self.n] - 1
        base = np.repeat(np.r_[0, csum[last[:-1]]], np.diff(np.r_[first, self.n]))
        self.before = csum - base - length
        self.total = np.repeat(
            csum[last] - np.r_[0, csum[last[:-1]]], np.diff(np.r_[first, self.n])
        )


def _introns_of(ex: _Exons):
    """Per-transcript intron instances: ``(ref, start, end, strand, exonic_before_start, exonic_total)``."""
    if ex.n < 2:
        z = np.zeros(0, np.int64)
        return ex.ref[:0], z, z, ex.strand[:0], z, z
    same = ex.tid[1:] == ex.tid[:-1]
    lo = ex.end[:-1][same]
    hi = ex.start[1:][same]
    ok = hi > lo
    return (
        ex.ref[:-1][same][ok],
        lo[ok],
        hi[ok],
        ex.strand[:-1][same][ok],
        (ex.before[:-1] + (ex.end - ex.start)[:-1])[same][ok],  # exonic bases before the intron
        ex.total[:-1][same][ok],
    )


# ---------------------------------------------------------------------------
# The builder
# ---------------------------------------------------------------------------


def build_splice_graph(
    transcripts: list[Transcript],
    ref_lengths: Mapping[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build ``(regions_df, edges_df)`` — the v8 splice graph. Fully vectorised; no per-region Python object.

    Deterministic by construction: ``np.unique`` sorts, region ids are assigned by position, boundaries by
    ``(src, kind, dst)``. No dict iteration order, no hashing, no parallel reduction ⇒ byte-identical
    across runs and platforms.
    """
    reflen = {str(k): int(v) for k, v in ref_lengths.items()}
    for name, length in reflen.items():
        if length < 0:
            raise ValueError(f"Reference {name!r} has negative length {length}.")

    ex = _Exons(transcripts, lambda tx: _is_real(tx, reflen))
    for name, s, e in ((str(r), int(a), int(b)) for r, a, b in zip(ex.ref, ex.start, ex.end)):
        if s < 0 or e > reflen[name]:
            raise ValueError(
                f"Transcript interval [{s}, {e}) is outside reference {name!r} length {reflen[name]}."
            )

    by_ref = _ref_slices(ex.ref)
    intr = _introns_of(ex)
    intr_by_ref = _ref_slices(intr[0])

    n_rows, e_rows = [], []
    region_base = 0
    for name, length in reflen.items():
        if length == 0:
            continue
        rows = by_ref.get(name, EMPTY_I64)
        i_rows = intr_by_ref.get(name, EMPTY_I64)
        # ── 1. the region_bound set ────────────────────────────────────────────────────────────────────────
        region_bounds = np.unique(
            np.concatenate([ex.start[rows], ex.end[rows], np.array([0, length], np.int64)])
        )
        starts, ends = region_bounds[:-1], region_bounds[1:]
        n = starts.shape[0]

        # ── 2. the 4-bit signature, by difference arrays over the region index ──────────────────────
        sig = np.zeros(n, np.uint8)
        for bit, iv_s, iv_e in _signature_intervals(ex, rows, intr, i_rows):
            if iv_s.size == 0:
                continue
            d = np.zeros(n + 1, np.int32)
            np.add.at(d, np.searchsorted(region_bounds, iv_s), 1)
            np.add.at(d, np.searchsorted(region_bounds, iv_e), -1)
            sig |= np.where(np.cumsum(d)[:-1] > 0, bit, 0).astype(np.uint8)

        n_rows.append(
            pd.DataFrame(
                {
                    "node_id": np.arange(region_base, region_base + n, dtype=np.int64),
                    "ref_name": name,
                    "start": starts,
                    "end": ends,
                    "length": ends - starts,
                    "signature": sig,
                }
            )
        )

        # ── 3. contiguous boundaries: flags + per-strand reaches at each interior interface ────────────
        pos = region_bounds[1:-1]  # the n-1 interior interfaces
        flags = np.zeros(pos.shape[0], np.uint16)
        if pos.size:
            for arr, bit in _terminus_and_splice_events(ex, rows, intr, i_rows):
                if arr.size:
                    hit = np.searchsorted(pos, arr)
                    ok = (hit < pos.size) & (pos[np.clip(hit, 0, max(pos.size - 1, 0))] == arr)
                    np.bitwise_or.at(flags, hit[ok], bit)
        c_reach = _contiguous_reaches(ex, rows, pos)

        # ── 4. junction boundaries: one per distinct (intron_start, intron_end, strand) ────────────────
        j_src, j_dst, j_strand, j_reach = _junction_edges(region_bounds, intr, i_rows, name)

        src = np.concatenate([np.arange(n - 1, dtype=np.int64), j_src]) + region_base
        dst = np.concatenate([np.arange(1, n, dtype=np.int64), j_dst]) + region_base
        kind = np.concatenate(
            [
                np.full(max(n - 1, 0), EDGE_KIND_CONTIGUOUS, np.uint8),
                np.full(j_src.size, EDGE_KIND_JUNCTION, np.uint8),
            ]
        )
        strand = np.concatenate([np.zeros(max(n - 1, 0), np.int8), j_strand])
        allflags = np.concatenate([flags, np.zeros(j_src.size, np.uint16)])
        reach = {k: np.concatenate([c_reach[k], j_reach[k]]) for k in c_reach}
        e_rows.append(
            pd.DataFrame(
                {
                    "src": src,
                    "dst": dst,
                    "kind": kind,
                    "strand": strand,
                    "flags": allflags,
                    **{k: v.astype(np.int32) for k, v in reach.items()},
                }
            )
        )
        region_base += n

    regions_df = (
        pd.concat(n_rows, ignore_index=True) if n_rows else pd.DataFrame(columns=REGION_COLUMNS)
    )
    edges_df = (
        pd.concat(e_rows, ignore_index=True) if e_rows else pd.DataFrame(columns=BOUNDARY_COLUMNS[1:])
    )
    if len(edges_df):
        # ⚠ sorted by (src, kind, dst, STRAND). The design doc says (src, kind, dst), but that is not a
        # TOTAL order: two strand-coincident junctions (§3.3, G18) share all three and differ only in
        # strand, so the order would be ambiguous and determinism would depend on the input order.
        # GENCODE contains zero such junctions — which is exactly why only a synthetic case can catch it.
        # Adding strand REFINES the documented order, so the "out-boundaries of a region are contiguous" CSR
        # contract is unchanged.
        edges_df = edges_df.sort_values(
            ["src", "kind", "dst", "strand"], kind="stable"
        ).reset_index(drop=True)
    edges_df.insert(0, "edge_id", np.arange(len(edges_df), dtype=np.int64))
    return _coerce(regions_df, REGION_COLUMNS, REGION_COLUMN_DTYPES), _coerce(
        edges_df, BOUNDARY_COLUMNS, BOUNDARY_COLUMN_DTYPES
    )


def _signature_intervals(ev: _Exons, ei, ev_i, ii):
    """The four ``(bit, starts, ends)`` interval sets that define a region's 4-bit signature.

    ONE definition, consumed by the builder (cumulative difference arrays over the region index) and
    by validator **I3b** (midpoint containment). Sharing the interval sets and differing in the
    evaluation is the point: the two can then only agree by both being right.
    """
    return (
        (
            BIT_EXON_POS,
            ev.start[ei][ev.strand[ei] == Strand.POS],
            ev.end[ei][ev.strand[ei] == Strand.POS],
        ),
        (
            BIT_EXON_NEG,
            ev.start[ei][ev.strand[ei] == Strand.NEG],
            ev.end[ei][ev.strand[ei] == Strand.NEG],
        ),
        (
            BIT_INTRON_POS,
            ev_i[1][ii][ev_i[3][ii] == Strand.POS],
            ev_i[2][ii][ev_i[3][ii] == Strand.POS],
        ),
        (
            BIT_INTRON_NEG,
            ev_i[1][ii][ev_i[3][ii] == Strand.NEG],
            ev_i[2][ii][ev_i[3][ii] == Strand.NEG],
        ),
    )


def _is_real(tx, reflen) -> bool:
    """⭐ THE transcript filter — ONE definition, shared by the builder and the validator.

    It was duplicated verbatim in both, in the module whose headline defect was exactly a filter
    divergence: I13 would have been validating the flags against its own private copy of the
    predicate, so the two could drift apart and agree with each other while both were wrong.

    Manufactured nRNA spans are ``is_synthetic``; a **non**-synthetic ``is_nrna`` row is a real
    single-exon transcript (mature ≡ nascent), whose termini are real. See the module docstring.
    """
    return (
        not tx.is_synthetic
        and bool(tx.exons)
        and tx.strand in (Strand.POS, Strand.NEG)
        and str(tx.ref) in reflen
    )


def _ref_slices(ref_arr):
    """reference name -> the row indices for that reference, in the flat per-exon/intron arrays.

    Rows are grouped contiguously by reference by construction, so this is one pass over the run
    boundaries. It replaced a dict-of-lists built by a Python loop over every row (1.0 s of a 3.5 s
    human build, and the same 'reference name -> rows' concept under a second name).
    """
    n = ref_arr.shape[0]
    if n == 0:
        return {}
    region_bound = np.flatnonzero(ref_arr[1:] != ref_arr[:-1]) + 1
    return {
        str(ref_arr[a]): np.arange(a, b, dtype=np.int64)
        for a, b in zip(np.r_[0, region_bound], np.r_[region_bound, n])
    }


#: The 8 structural bits, as (strand, TSS, TES, DONOR, ACCEPTOR). One tuple, so the builder and the
#: validator cannot disagree about which bit means what.
_BIT_SPEC = (
    (Strand.POS, FLAG_TSS_POS, FLAG_TES_POS, FLAG_DONOR_POS, FLAG_ACCEPTOR_POS),
    (Strand.NEG, FLAG_TSS_NEG, FLAG_TES_NEG, FLAG_DONOR_NEG, FLAG_ACCEPTOR_NEG),
)


def _terminus_and_splice_events(ma: _Exons, mj, ma_i, mi_):
    """(positions, flag_bit) for each of the 8 structural bits — the BUILDER's derivation.

    ⚠ Emits an entry only for a class that HAS events, which is right for building (an empty class
    sets no bit) and wrong for validating. :func:`_events_independently` is I13's counterpart: all
    eight classes, always, by a different algorithm.
    """
    out = []
    for st, tss_bit, tes_bit, don_bit, acc_bit in _BIT_SPEC:
        m = mj[ma.strand[mj] == st] if mj.size else mj
        if m.size:
            # a transcript's own outer boundaries: exonic offset 0 and offset == total
            first = m[ma.before[m] == 0]
            last = m[(ma.before[m] + (ma.end - ma.start)[m]) == ma.total[m]]
            lo_boundary, hi_boundary = ma.start[first], ma.end[last]
            # 5' terminus is the genomically-low boundary on +, the high boundary on −
            tss = lo_boundary if st == Strand.POS else hi_boundary
            tes = hi_boundary if st == Strand.POS else lo_boundary
            out.append((np.unique(tss), tss_bit))
            out.append((np.unique(tes), tes_bit))
        k = mi_[ma_i[3][mi_] == st] if mi_.size else mi_
        if k.size:
            out.append((np.unique(ma_i[1][k]), don_bit))  # intron_start
            out.append((np.unique(ma_i[2][k]), acc_bit))  # intron_end
    return out


def _events_independently(ex: _Exons, rows, intr, i_rows):
    """⭐ I13's SECOND OPINION: the same 8 event sets, derived by a DIFFERENT algorithm.

    The same discipline as I3b. The builder finds a transcript's outer exons by cumulative-exonic-
    offset arithmetic (``before == 0`` / ``before + len == total``); this finds them by per-transcript
    ``min(start)`` and ``max(end)`` via ``np.minimum.at`` — no cumulative offsets at all. The two can
    then only agree by both being right.

    ⚠ **All EIGHT entries, always**, even for a class with zero events on this reference. The builder
    emits only non-empty classes, which is correct for building and useless for validating: a
    spurious ``TSS_NEG`` on a reference with no NEG transcript would never be compared. Returned as
    ``(positions, bit)``, positions possibly empty.
    """
    out = []
    if rows.size:
        tid = ex.tid[rows]
        uniq, inv = np.unique(tid, return_inverse=True)
        lo = np.full(uniq.size, np.iinfo(np.int64).max, np.int64)
        hi = np.full(uniq.size, np.iinfo(np.int64).min, np.int64)
        np.minimum.at(lo, inv, ex.start[rows])
        np.maximum.at(hi, inv, ex.end[rows])
        strand_of = np.zeros(uniq.size, np.int8)
        strand_of[inv] = ex.strand[rows]
    for st, tss_bit, tes_bit, don_bit, acc_bit in _BIT_SPEC:
        if rows.size:
            s = strand_of == st
            five, three = (lo[s], hi[s]) if st == Strand.POS else (hi[s], lo[s])
        else:
            five = three = EMPTY_I64
        out.append((np.unique(five), tss_bit))
        out.append((np.unique(three), tes_bit))
        k = i_rows[intr[3][i_rows] == st] if i_rows.size else i_rows
        out.append((np.unique(intr[1][k]) if k.size else EMPTY_I64, don_bit))
        out.append((np.unique(intr[2][k]) if k.size else EMPTY_I64, acc_bit))
    return out


def _contiguous_reaches(ex: _Exons, rows, pos):
    """Per-strand RNA reach either side of each interior interface, over TRANSCRIPT SPANS.

    "Reach" is how much RNA molecule can still exist on each side of the interface — the quantity the
    fragment-placement count is truncated by near a transcript end.

    Nascent RNA is an ordinary transcript that happens to be single-exon and to span its whole gene, and
    a genomic distance is never shorter than the exonic distance inside the same span. So the widest RNA
    molecule covering a position is the nascent one, and the reach is the **genomic** distance to the
    ends of a covering transcript's span — maximised per side and per strand, independently (TRAPS: two-gaussians-one-latent).

    ⚠ In particular the reach inside an intron is **not zero**: that is precisely where nascent RNA
    lives. An exonic reach declares zero RNA opportunity across every intron in the genome, which is
    backwards.

    ⚠ This is the **contiguous-boundary** rule only. A junction boundary is used solely by a molecule that
    spliced across it, so what remains either side of it is exonic — see :func:`_junction_edges`, which
    is deliberately left on the exonic reach.
    """
    keys = ("reach_lo_pos", "reach_hi_pos", "reach_lo_neg", "reach_hi_neg")
    out = {k: np.zeros(pos.shape[0], np.int64) for k in keys}
    if pos.size == 0 or rows.size == 0:
        return out

    # per-transcript genomic span on this reference: [min exon start, max exon end)
    uniq, inv = np.unique(ex.tid[rows], return_inverse=True)
    span_start = np.full(uniq.size, np.iinfo(np.int64).max, np.int64)
    span_end = np.full(uniq.size, np.iinfo(np.int64).min, np.int64)
    np.minimum.at(span_start, inv, ex.start[rows])
    np.maximum.at(span_end, inv, ex.end[rows])
    strand_of = np.zeros(uniq.size, np.int8)
    strand_of[inv] = ex.strand[rows]

    for s, key_lo, key_hi in (
        (Strand.POS, "reach_lo_pos", "reach_hi_pos"),
        (Strand.NEG, "reach_lo_neg", "reach_hi_neg"),
    ):
        keep = strand_of == s
        if not keep.any():
            continue
        lo_boundary, hi_boundary = span_start[keep], span_end[keep]

        # HIGH side: among spans starting at or before p, the one reaching furthest right. If that
        # furthest end is still left of p then no span covers p and the reach is 0.
        n_span = lo_boundary.size
        by_start = np.argsort(lo_boundary, kind="stable")
        furthest_end = np.maximum.accumulate(hi_boundary[by_start])
        i = np.searchsorted(lo_boundary[by_start], pos, side="right")
        covered = i > 0  # some span starts at or before p
        far = furthest_end[np.maximum(i, 1) - 1]  # clamped: the value is discarded where ~covered
        out[key_hi] = np.maximum(np.where(covered, far, pos) - pos, 0)

        # LOW side: among spans ending at or after p, the one starting furthest left. Symmetric.
        by_end = np.argsort(hi_boundary, kind="stable")
        furthest_start = np.minimum.accumulate(lo_boundary[by_end][::-1])[::-1]
        j = np.searchsorted(hi_boundary[by_end], pos, side="left")
        covered = j < n_span  # some span ends at or after p
        near = furthest_start[np.minimum(j, n_span - 1)]  # clamped, as above
        out[key_lo] = np.maximum(pos - np.where(covered, near, pos), 0)
    return out


def _junction_edges(region_bounds, ma_i, mi_, ref_name):
    """Distinct junction boundaries for one reference: ``(src, dst, strand, reach dict)``, region-local ids."""
    empty = np.zeros(0, np.int64)
    keys = {
        k: empty.copy() for k in ("reach_lo_pos", "reach_hi_pos", "reach_lo_neg", "reach_hi_neg")
    }
    if mi_.size == 0:
        return empty, empty, np.zeros(0, np.int8), keys
    a, b, st = ma_i[1][mi_], ma_i[2][mi_], ma_i[3][mi_]
    lo, tot = ma_i[4][mi_], ma_i[5][mi_]
    key = np.stack([a, b, st.astype(np.int64)], axis=1)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    ua, ub, ust = uniq[:, 0], uniq[:, 1], uniq[:, 2].astype(np.int8)
    src = np.searchsorted(region_bounds, ua) - 1  # the region ENDING at intron_start
    dst = np.searchsorted(region_bounds, ub)  # the region STARTING at intron_end
    if not (np.all(region_bounds[src + 1] == ua) and np.all(region_bounds[dst] == ub)):
        raise ValueError(f"ref {ref_name!r}: a junction endpoint is not a region interface (I5)")
    out = {k: np.zeros(uniq.shape[0], np.int64) for k in keys}
    for s, kl, kh in (
        (Strand.POS, "reach_lo_pos", "reach_hi_pos"),
        (Strand.NEG, "reach_lo_neg", "reach_hi_neg"),
    ):
        m = st == s
        if m.any():
            np.maximum.at(out[kl], inv[m], lo[m])
            np.maximum.at(out[kh], inv[m], (tot - lo)[m])
    return src, dst, ust, out


def _coerce(df, columns, dtypes):
    df = df.copy()
    for c in columns:
        if c not in df.columns:
            df[c] = pd.Series(dtype="string" if c == "ref_name" else dtypes.get(c, object))
    df = df.loc[:, columns]
    for c, dt in dtypes.items():
        if df[c].dtype != dt:
            df[c] = df[c].astype(dt)
    if "ref_name" in columns:
        df["ref_name"] = df["ref_name"].astype("string")
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Validation — I1..I13. ⚠ I3b, I4, I11 and I13 need the transcripts and therefore run at
# BUILD ONLY; a load-time call passes none and gets the graph-internal checks (I1/I2/I5-I9/I12).
# ---------------------------------------------------------------------------


def validate_graph(regions_df, edges_df, ref_lengths: Mapping[str, int], transcripts=None) -> None:
    """Assert invariants I1–I13. Raises :class:`ValueError` on the first violation.

    ⚠ **``transcripts`` gates four invariants, not one**: I3b (the signature, recomputed), I4 (the
    interfaces ARE the events), I11 (every transcript walks) and I13 (the flags ARE the events).
    Without it only the graph-internal checks run — I1/I2/I5–I9/I12 — which cannot detect a
    ``signature`` or ``flags`` column that has drifted from the annotation. The build passes
    transcripts; the load path does not, because reconstructing them costs ~3 s at human scale.
    """
    reflen = {str(k): int(v) for k, v in ref_lengths.items()}
    nid = regions_df["node_id"].to_numpy(np.int64)
    start = regions_df["start"].to_numpy(np.int64)
    end = regions_df["end"].to_numpy(np.int64)
    length = regions_df["length"].to_numpy(np.int64)
    sig = regions_df["signature"].to_numpy(np.uint8)
    ref = regions_df["ref_name"].astype(str).to_numpy()

    # I2 — ids are 0..n-1 in row order
    if not np.array_equal(nid, np.arange(nid.size, dtype=np.int64)):
        raise ValueError("I2: region_id is not 0..n-1 in row order. Rebuild the index.")
    # I3a — signature range (I3b, the independent recomputation, needs the transcripts)
    if nid.size and int(sig.max()) >= N_SIGNATURES:
        raise ValueError(f"I3: signature {int(sig.max())} out of range. Rebuild the index.")
    if not np.array_equal(length, end - start):
        raise ValueError("I1: stored region length != end - start. Rebuild the index.")
    if nid.size and int(length.min()) <= 0:
        raise ValueError("I1: a region has non-positive length. Rebuild the index.")

    # I1 — regions tile each reference exactly. ⚠ Reference slices are computed ONCE from the run
    # boundaries; `np.flatnonzero(ref == name)` inside the loop is 286 full-array string comparisons
    # over 1.04 M rows and cost 3.5 s on its own.
    slices = _ref_slices(ref)
    for name, L in reflen.items():
        sl = slices.get(name)
        if sl is None:
            if L > 0:
                raise ValueError(f"I1: reference {name!r} (length {L}) has no regions.")
            continue
        lo, hi = int(sl[0]), int(sl[-1]) + 1
        if start[lo] != 0 or end[hi - 1] != L:
            raise ValueError(
                f"I1: reference {name!r} spans [{start[lo]}, {end[hi - 1]}), expected [0, {L})."
            )
        if hi - lo > 1 and not np.array_equal(end[lo : hi - 1], start[lo + 1 : hi]):
            raise ValueError(f"I1: reference {name!r} has a gap or overlap between regions.")
    if len(slices) != sum(1 for _n, L in reflen.items() if L > 0):
        raise ValueError("I1: region rows are not grouped contiguously by reference.")

    src = edges_df["src"].to_numpy(np.int64)
    dst = edges_df["dst"].to_numpy(np.int64)
    kind = edges_df["kind"].to_numpy(np.uint8)
    estrand = edges_df["strand"].to_numpy(np.int8)

    # I8 — src < dst on EVERY boundary  ⇒ genomic order is a topological order
    if src.size and not np.all(src < dst):
        raise ValueError("I8: an boundary has src >= dst; genomic order is not a topological order.")
    if src.size and not np.all(ref[src] == ref[dst]):
        raise ValueError("I5: an boundary crosses references.")
    # I12 — row order is the id order
    if not np.array_equal(edges_df["edge_id"].to_numpy(np.int64), np.arange(src.size)):
        raise ValueError("I12: edge_id is not the row index.")
    if src.size > 1:
        # (src, kind, dst, strand) packed into one strictly-increasing int64 key. STRAND is part of the
        # key because it is part of the order (see the sort in build_splice_graph): without it two
        # strand-coincident junctions collide and read as a duplicate.
        n_regions_ = max(nid.size, 1)
        key = ((src * 2 + kind.astype(np.int64)) * n_regions_ + dst) * 4 + estrand.astype(np.int64)
        if not np.all(np.diff(key) > 0):
            raise ValueError(
                "I12: boundaries are not sorted by (src, kind, dst, strand), or a duplicate exists."
            )

    # I7 — contiguous boundaries are exactly {(i, i+1)} within each reference
    c = kind == EDGE_KIND_CONTIGUOUS
    if not np.all(dst[c] == src[c] + 1):
        raise ValueError("I7: a CONTIGUOUS boundary is not between adjacent regions.")
    # ⚠ `sum(... (ref == name).any())` here was 286 full-array object-dtype comparisons over 1.04 M
    # rows — 1.68 s of a 1.91 s load-time validate_graph, i.e. 88 % of it — and it recomputed a
    # number `slices` already holds. This is the same trap the I1 block above documents.
    n_expected = nid.size - len(slices)
    if int(c.sum()) != n_expected:
        raise ValueError(f"I7: {int(c.sum())} contiguous boundaries, expected {n_expected}.")
    if np.any(estrand[c] != 0):
        raise ValueError("I7: a CONTIGUOUS boundary carries a strand.")

    # I5/I6 — junctions land on interfaces, one row per distinct (ref, donor, acceptor, strand)
    j = kind == EDGE_KIND_JUNCTION
    if j.any():
        if not np.all(np.isin(estrand[j], [Strand.POS, Strand.NEG])):
            raise ValueError("I6: a JUNCTION boundary has an invalid strand.")
        # ⚠ STRAND-COINCIDENT junctions are BIOLOGICALLY IMPOSSIBLE and are reported, not tolerated
        # silently. Splicing reads a non-palindromic motif: the reverse complement of a GT..AG intron
        # begins CT, so the same interval cannot be a valid intron on both strands. Measured: ZERO in
        # GENCODE. One in a GTF therefore means the ANNOTATION is wrong (a simulator, or a strand
        # column filled in by hand), and silence would let it propagate into every downstream density.
        # It is a WARNING rather than an error on purpose: the graph handles the case correctly (the
        # sort key carries strand, so the two boundaries stay distinct and ordered), and the G18 test builds
        # exactly this input to prove that. Raising would make the guard untestable.
        pair = src[j] * max(nid.size, 1) + dst[j]
        _u, _cnt = np.unique(pair, return_counts=True)
        n_coincident = int((_cnt > 1).sum())
        if n_coincident:
            warnings.warn(
                f"{n_coincident} strand-coincident splice junction(s): the same (donor, acceptor) "
                f"annotated on BOTH strands. This is biologically impossible — splice motifs are "
                f"non-palindromic (a GT..AG intron reverse-complements to CT..AC) — so the source "
                f"annotation is likely wrong. The graph keeps them as distinct boundaries and is correct "
                f"either way.",
                RuntimeWarning,
                stacklevel=2,
            )
        key = np.stack([src[j], dst[j], estrand[j].astype(np.int64)], axis=1)
        if np.unique(key, axis=0).shape[0] != key.shape[0]:
            raise ValueError("I6: a (donor, acceptor, strand) junction appears more than once.")
        # I9 — both endpoints carry that strand's exon bit
        for s, bit in ((Strand.POS, BIT_EXON_POS), (Strand.NEG, BIT_EXON_NEG)):
            m = j & (estrand == s)
            if m.any() and not (np.all(sig[src[m]] & bit) and np.all(sig[dst[m]] & bit)):
                raise ValueError(f"I9: a strand-{int(s)} junction endpoint lacks its exon bit.")

    # I3b — signature recomputed; I4 — interfaces ARE the events; I11 — every transcript walks;
    # I13 — the flags ARE the events
    if transcripts is not None:
        _validate_against_transcripts(transcripts, reflen, regions_df, edges_df)


def _validate_against_transcripts(transcripts, reflen, regions_df, edges_df) -> None:
    """I3b (the signature, recomputed) + I4 (interfaces ARE the events) + I11 (every transcript
    walks) + I13 (the flags ARE the events), vectorised.

    ⚠ Written without a per-transcript Python loop on purpose: the loop form cost 10.4 s on the human
    annotation against §8's 5 s budget, and it is the check most likely to be skipped when it is slow —
    which would be exactly backwards, since I11 is the one that catches real bugs.
    """
    ref = regions_df["ref_name"].astype(str).to_numpy()
    start = regions_df["start"].to_numpy(np.int64)
    end = regions_df["end"].to_numpy(np.int64)
    sig = regions_df["signature"].to_numpy(np.uint8)
    src = edges_df["src"].to_numpy(np.int64)
    dst = edges_df["dst"].to_numpy(np.int64)
    kind = edges_df["kind"].to_numpy(np.uint8)
    estrand = edges_df["strand"].to_numpy(np.int8)

    ex = _Exons(transcripts, lambda tx: _is_real(tx, reflen))
    by_ref = _ref_slices(ex.ref)
    intr = _introns_of(ex)
    intr_by_ref = _ref_slices(intr[0])

    cm = kind == EDGE_KIND_CONTIGUOUS
    csrc = src[cm]  # already ascending: boundaries sort by (src, kind, dst, strand)
    cflags = edges_df["flags"].to_numpy(np.uint16)[cm]

    n_regions = ref.size
    jm = kind == EDGE_KIND_JUNCTION
    jkey = np.sort((src[jm] * n_regions + dst[jm]) * 4 + estrand[jm].astype(np.int64))

    slices = _ref_slices(ref)
    for name, L in reflen.items():
        m = slices.get(name)
        if m is None:
            continue
        base, region_bounds = int(m[0]), np.concatenate([start[m], end[m][-1:]])

        # I4 — the interior interfaces are EXACTLY the annotation events
        rows = by_ref.get(name, EMPTY_I64)
        i_rows = intr_by_ref.get(name, EMPTY_I64)
        events = (
            np.unique(np.concatenate([ex.start[rows], ex.end[rows]])) if rows.size else EMPTY_I64
        )
        events = events[(events > 0) & (events < L)]
        interior = end[m[:-1]] if m.size > 1 else EMPTY_I64
        if not np.array_equal(np.unique(interior), events):
            extra = np.setdiff1d(interior, events).size
            missing = np.setdiff1d(events, interior).size
            raise ValueError(
                f"I4: reference {name!r} interior interfaces are not exactly the annotation events "
                f"({extra} extra, {missing} missing)."
            )

        # ⭐ I3b — recompute the signature from each region's MIDPOINT, by direct interval
        # containment. Deliberately a DIFFERENT algorithm from the builder's cumulative-difference
        # sweep, so the two can only agree by both being right. A region is homogeneous by
        # construction (every interval endpoint is a region_bound), so its midpoint decides it.
        mid = (start[m] + end[m]) // 2
        want = np.zeros(mid.size, np.uint8)
        for bit, iv_s, iv_e in _signature_intervals(ex, rows, intr, i_rows):
            if iv_s.size:
                # covered iff #(start <= mid) > #(end <= mid), for half-open [start, end)
                n_started = np.searchsorted(np.sort(iv_s), mid, side="right")
                n_ended = np.searchsorted(np.sort(iv_e), mid, side="right")
                want |= np.where(n_started > n_ended, bit, 0).astype(np.uint8)
        if not np.array_equal(want, sig[m]):
            k = int(np.flatnonzero(want != sig[m])[0])
            raise ValueError(
                f"I3: reference {name!r} region {int(m[k])} [{int(start[m][k])}, "
                f"{int(end[m][k])}) has signature {int(sig[m][k])}, recomputed {int(want[k])}."
            )

        # ⭐ I13 — each structural bit is set at EXACTLY the positions that generate it.
        #
        # Both directions, per bit: no event without its flag, no flag without its event. It
        # subsumes the design doc's I10 (terminus counts) and it is the check that would have caught
        # the flag-filter bug immediately — the flags were built from `~is_synthetic & ~is_nrna`,
        # which on a NON-synthetic row does not mean "manufactured span" but "single-exon
        # transcript", so the TSS/TES of 26,475 real human transcripts were silently dropped. Nothing
        # else noticed, because the flags had no consumer yet and so no number moved.
        #
        # ⚠ It compares ALL EIGHT bits, including classes with zero events on this reference —
        # `_events_independently` returns an empty position set rather than omitting the entry. The
        # builder's own emitter omits empty classes (correct for building, useless for validating),
        # and with it a spurious TSS_NEG on a reference with no NEG transcript was never compared.
        #
        # ⚠ The converse of I4 does NOT hold bit-wise: an interior interface may carry NO flag at
        # all, when adjacent exons of one transcript are BOOKENDED (a zero-length intron). That is an
        # exon endpoint which is neither a terminus nor a splice site — a real third case, measured
        # ZERO times in GENCODE but constructed by G14.
        if m.size > 1:
            interior = end[m[:-1]]
            at = np.zeros(interior.size, np.uint16)
            e_of_region = np.searchsorted(csrc, m[:-1])
            hit = (e_of_region < csrc.size) & (
                csrc[np.clip(e_of_region, 0, max(csrc.size - 1, 0))] == m[:-1]
            )
            at[hit] = cflags[e_of_region[hit]]
            for arr, bit in _events_independently(ex, rows, intr, i_rows):
                want = np.isin(interior, arr)
                got = (at & bit) != 0
                if not np.array_equal(want, got):
                    k = int(np.flatnonzero(want != got)[0])
                    raise ValueError(
                        f"I13: reference {name!r} position {int(interior[k])} "
                        f"{'is missing' if want[k] else 'wrongly carries'} flag bit {int(bit):#06x}."
                    )

        # I11a — every mature exon boundary lands exactly on a region interface
        if rows.size:
            for arr, what in ((ex.start[rows], "start"), (ex.end[rows], "end")):
                hit = np.searchsorted(region_bounds, arr)
                bad = (hit >= region_bounds.size) | (region_bounds[np.clip(hit, 0, region_bounds.size - 1)] != arr)
                if bad.any():
                    raise ValueError(
                        f"I11: reference {name!r} has {int(bad.sum())} exon {what}s that are not region "
                        f"interfaces (first at {int(arr[bad][0])})."
                    )

        # I11b — every mature intron is realised by a JUNCTION boundary on that strand
        if i_rows.size:
            a, b, st = intr[1][i_rows], intr[2][i_rows], intr[3][i_rows].astype(np.int64)
            js = base + np.searchsorted(region_bounds, a) - 1
            jd = base + np.searchsorted(region_bounds, b)
            want = (js * n_regions + jd) * 4 + st
            if jkey.size == 0:
                miss = np.ones(want.shape, bool)
            else:
                pos = np.searchsorted(jkey, want)
                miss = (pos >= jkey.size) | (jkey[np.clip(pos, 0, jkey.size - 1)] != want)
            if miss.any():
                k = int(np.flatnonzero(miss)[0])
                raise ValueError(
                    f"I11: reference {name!r} has {int(miss.sum())} mature introns with no JUNCTION "
                    f"boundary (first [{int(a[k])}, {int(b[k])}) strand {int(st[k])})."
                )


# ---------------------------------------------------------------------------
# The scanner adapter + load
# ---------------------------------------------------------------------------


def build_region_partition_arrays(index) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flatten the graph into the ``BamScanner.set_regions`` ABI — THE partition the accumulator
    deposits into.

    For each reference in ``index.ref_names`` the region_bound positions are ``[n0.start, n0.end, n1.end, …]``
    over that reference's regions (which tile it contiguously), so ``k`` regions give ``k + 1`` positions.
    References with zero regions contribute none.

    :meth:`RegionArrays.from_index <rigel.calibration.region_arrays.RegionArrays.from_index>` reads
    the same frame; that is what keeps the calibration geometry and the scanner from addressing
    different partitions.

    Returns ``(region_bounds int64[P], ref_pos_offsets int64[n_refs+1], region_types uint8[N])``, where
    ``region_types`` is the coarse type (0=intergenic, 1=intron, 2=exon) aligned 1:1 with the regions —
    the gDNA FL-pool axis.
    """
    regions_df = index.regions_df
    ref_names = index.ref_names
    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in regions_df.groupby("ref_name", sort=False)
    }
    positions: list[np.ndarray] = []
    types: list[np.ndarray] = []
    offsets = np.zeros(len(ref_names) + 1, dtype=np.int64)
    for i, ref in enumerate(ref_names):
        grp = by_ref.get(ref)
        if grp is None or len(grp) == 0:
            positions.append(np.empty(0, dtype=np.int64))
            types.append(np.empty(0, dtype=np.uint8))
            offsets[i + 1] = offsets[i]
            continue
        starts = grp["start"].to_numpy(np.int64, copy=False)
        ends = grp["end"].to_numpy(np.int64, copy=False)
        region_bounds = np.empty(len(starts) + 1, dtype=np.int64)
        region_bounds[:-1] = starts
        region_bounds[-1] = ends[-1]
        positions.append(region_bounds)
        types.append(coarse_type_array(grp["signature"].to_numpy()))
        offsets[i + 1] = offsets[i] + region_bounds.shape[0]
    return (
        np.concatenate(positions) if positions else np.empty(0, dtype=np.int64),
        offsets,
        np.concatenate(types) if types else np.empty(0, dtype=np.uint8),
    )


@dataclass(frozen=True, slots=True)
class JunctionEdgeArrays:
    """The junction boundaries, re-indexed onto **the accumulator's flat region_bound axis** as a CSR.

    The deposit rule has to answer one question per observed intron, in the BAM-scan hot loop: *is this
    intron an annotated junction, and if so which boundary?* This is the lookup table for it.

    ⭐ It is cheap because of a measured property of the graph: **every annotated intron has both of its
    endpoints as partition region_bounds** (100.00 % of 404,168 on the human annotation — forced, since a region_bound is
    placed at every exon endpoint). So "is this intron annotated?" reduces to *"are both endpoints region_bounds,
    and is the pair registered?"* — and finding the region_bound index is the binary search the deposit already
    performs to locate the crossed boundaries. If an intron's start is not a region_bound it is unannotated and the
    table is never consulted.

    Keyed by the **donor** region_bound index, i.e. the flat index into
    ``build_region_partition_arrays(index)[0]`` of the intron's LOW endpoint::

        for k in range(offsets[boundary_left], offsets[boundary_left + 1]):
            if boundary_right[k] == observed_right_boundary:   # 1-3 iterations at human scale
                credit junction boundary k          # <- the SLOT is the id; see below

    ⚠ **The junction-boundary id IS the CSR slot ``k``.** ``edge_row`` is *not* the id — it is the key for
    joining a payload row back to ``index.edges_df``, and using it to index a junction bank is a
    memory-safety bug: the highest junction row is **1,447,755** in a **404,168**-entry bank, so the write
    lands past the end, and even in bounds it permutes every surviving row. The accumulator never receives
    this column.

    ⚠ **The slot ordering is a contract**, because the id is the rank: this function and the reference
    accumulator's ``Partition.from_region_bounds`` must sort identically, or every row permutes and the native
    build's byte-identity gate compares two different labellings. They disagreed once — ``(acceptor,
    donor)`` against ``(strand, acceptor, donor)`` — and it is reachable only through a strand-coincident
    junction pair. Pinned by ``test_the_csr_slot_order_matches_the_reference_accumulator``.

    ``strand`` disambiguates that (biologically improbable, constructible) case of two junctions sharing a
    coordinate pair and differing only in strand; the caller matches on it when the observed motif strand
    is known.
    """

    offsets: np.ndarray  # int64[P + 1] — CSR over the flat region_bound axis, P = region_bounds.size
    boundary_right: np.ndarray  # int64[J] — flat region_bound index of the intron's HIGH endpoint
    edge_row: np.ndarray  # int64[J] — row in index.edges_df. A JOIN KEY, not the junction-boundary id
    strand: np.ndarray  # int8[J]  — the junction's genomic strand (Strand POS/NEG)


def build_junction_edge_arrays(index) -> JunctionEdgeArrays:
    """Build the :class:`JunctionEdgeArrays` CSR for ``index``.

    A junction boundary stores region ids; the accumulator works in region_bound indices. For a reference whose first
    region is ``region_base`` and whose first region_bound is ``region_bound_base``, region ``i`` spans
    ``[region_bounds[region_bound_base + i - region_base], region_bounds[region_bound_base + i - region_base + 1])``, so::

        donor region_bound    = region_bound_base + (src - region_base) + 1      # the intron starts where src ENDS
        acceptor region_bound = region_bound_base + (dst - region_base)          # and ends where dst BEGINS
    """
    regions_df, edges_df = index.regions_df, index.edges_df
    _, region_bound_offsets, _ = build_region_partition_arrays(index)
    n_region_bounds = int(region_bound_offsets[-1])

    # per-reference region_base, in the same reference order the region_bound axis uses
    counts = regions_df.groupby("ref_name", sort=False).size()
    region_base = np.zeros(len(index.ref_names) + 1, dtype=np.int64)
    for i, ref in enumerate(index.ref_names):
        region_base[i + 1] = region_base[i] + int(counts.get(ref, 0))

    is_junction = edges_df["kind"].to_numpy(np.uint8) == EDGE_KIND_JUNCTION
    src = edges_df["src"].to_numpy(np.int64)[is_junction]
    dst = edges_df["dst"].to_numpy(np.int64)[is_junction]
    strand = edges_df["strand"].to_numpy(np.int8)[is_junction]
    edge_row = np.flatnonzero(is_junction).astype(np.int64)

    # which reference each junction belongs to: region ids are contiguous per reference (I2)
    ref_of = np.searchsorted(region_base, src, side="right") - 1
    shift = region_bound_offsets[ref_of] - region_base[ref_of]
    boundary_left = src + shift + 1
    boundary_right = dst + shift

    # ⚠ The sort key includes STRAND, and it must match the reference accumulator's
    # ``Partition.from_region_bounds`` exactly — the junction-boundary id IS the rank in this order, so the two
    # disagreeing would permute every row and break byte-identity. It shows up only on a donor/acceptor
    # pair carrying two strands, which is biologically impossible and therefore only ever reachable from
    # a synthetic stress test.
    order = np.lexsort((strand, boundary_right, boundary_left))
    boundary_left, boundary_right = boundary_left[order], boundary_right[order]
    offsets = np.zeros(n_region_bounds + 1, dtype=np.int64)
    np.cumsum(np.bincount(boundary_left, minlength=n_region_bounds), out=offsets[1:])
    return JunctionEdgeArrays(offsets, boundary_right, edge_row[order], strand[order])


def build_boundary_flags_array(index) -> np.ndarray:
    """The structural flags on **the accumulator's contiguous-boundary axis** — one entry per boundary.

    ⭐ **There is no padding, because there are no terminal slots.** The predecessor
    (``build_boundary_flags_array``) emitted ``k + 1`` entries per reference — the ``k − 1`` interior
    interfaces plus two data-free terminals — purely so every region had an object on each side. A
    contiguous boundary is the boundary BETWEEN two adjacent regions and there is no such boundary before the first or
    after the last, so a reference with ``k`` regions contributes exactly ``k − 1`` entries and the
    off-by-one commentary that used to live here goes with the slots.

    A contiguous boundary IS the interface to the right of its ``src`` region, so the flags are keyed by
    ``src``. Junction boundaries carry no flags — they are not a genomic position — and are excluded.

    Returns ``uint16[E]`` with ``E == ref_boundary_offsets[-1]``, aligned element for element with the
    payload's contiguous-boundary axis.
    """
    regions_df, edges_df = index.regions_df, index.edges_df
    contiguous = edges_df["kind"].to_numpy(np.uint8) == EDGE_KIND_CONTIGUOUS
    by_region = np.zeros(len(regions_df), dtype=np.uint16)
    by_region[edges_df["src"].to_numpy(np.int64)[contiguous]] = edges_df["flags"].to_numpy(np.uint16)[
        contiguous
    ]

    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in regions_df.groupby("ref_name", sort=False)
    }
    out: list[np.ndarray] = []
    for ref in index.ref_names:
        grp = by_ref.get(ref)
        if grp is None or len(grp) == 0:
            continue
        ids = grp.index.to_numpy(np.int64)  # == region_id (I2), and contiguous within a reference
        out.append(by_region[ids[:-1]])  # the k-1 interior boundaries; a 1-region reference contributes none
    return np.concatenate(out) if out else np.zeros(0, dtype=np.uint16)


def build_contiguous_boundary_reach_arrays(index) -> tuple[np.ndarray, np.ndarray]:
    """The RNA **reach** on the accumulator's contiguous-boundary axis — ``(reach_lo, reach_hi)``,
     ``float64[E, 2]``, column 0 the POS-strand transcript's and column 1 the NEG's.

     A crossing molecule must fit in what remains of **its own template** either side of the boundary. gDNA's
     template is the chromosome, so its reach is unbounded — physics, not a choice. RNA's ends where its
     transcript ends, and ignoring that over-calls gDNA by a measured **11.0 %** genome-wide
    This is the array that lets :func:`effective_length.crossing_eff_length`
     taper the RNA divisor.

     ⚠ **PER STRAND and per SIDE**: reach is "maximised over transcripts
     independently per side AND per strand". A POS transcript and a NEG one ending in different places
     give one boundary two different RNA reaches, and a single averaged number describes neither.

     ⚠ **GENOMIC, unlike a junction's EXONIC reach.** A junction is used only by a spliced molecule, so
     what remains either side of it is exonic; a contiguous boundary is also crossed by *nascent* RNA, which
     is genomic. Taking the exonic reach here would declare an intronic nascent fragment impossible
     (:class:`JunctionGeometry`).

     ⚠ **A reach of 0 is the ANSWER, not a missing value** — no template of that strand at that boundary, so
     that strand's RNA has zero opportunity and the divisor is legitimately 0. The consumer must treat 0
     as "emit nothing" rather than flooring it. Measured on the chr22
     pilot index: **40.6 %** of contiguous boundaries have no POS template and **42.9 %** no NEG template.

     ⭐ **Keyed by ``src``, exactly as :func:`build_boundary_flags_array` is**, and laid out per reference in
     ``index.ref_names`` order, so the two arrays are the SAME axis element for element and a consumer
     indexes both with one index. A reference with ``k`` regions contributes ``k − 1`` entries; one with a
     single region contributes none.
    """
    regions_df, edges_df = index.regions_df, index.edges_df
    contiguous = edges_df["kind"].to_numpy(np.uint8) == EDGE_KIND_CONTIGUOUS
    src = edges_df["src"].to_numpy(np.int64)[contiguous]

    # (n_regions, 2) scratch keyed by the region whose RIGHT interface the boundary is, then sliced per
    # reference. Regions with no outgoing contiguous boundary keep 0, which is also the correct reach for a
    # boundary that does not exist — they are dropped by the ``[:-1]`` slice below either way.
    def by_region(column_pos: str, column_neg: str) -> np.ndarray:
        out = np.zeros((len(regions_df), 2), dtype=np.float64)
        out[src, 0] = edges_df[column_pos].to_numpy(np.float64)[contiguous]
        out[src, 1] = edges_df[column_neg].to_numpy(np.float64)[contiguous]
        return out

    lo_by_region = by_region("reach_lo_pos", "reach_lo_neg")
    hi_by_region = by_region("reach_hi_pos", "reach_hi_neg")

    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in regions_df.groupby("ref_name", sort=False)
    }
    lo_out: list[np.ndarray] = []
    hi_out: list[np.ndarray] = []
    for ref in index.ref_names:
        grp = by_ref.get(ref)
        if grp is None or len(grp) == 0:
            continue
        ids = grp.index.to_numpy(np.int64)  # == region_id (I2), contiguous within a reference
        lo_out.append(lo_by_region[ids[:-1]])  # the k-1 interior boundaries
        hi_out.append(hi_by_region[ids[:-1]])
    empty = np.zeros((0, 2), dtype=np.float64)
    return (
        np.concatenate(lo_out) if lo_out else empty,
        np.concatenate(hi_out) if hi_out else empty.copy(),
    )


@dataclass(frozen=True, slots=True)
class JunctionGeometry:
    """The junction boundaries on **the accumulator's junction axis**, in its own slot order.

    A junction boundary is not a chain slot — the graph is a DAG but not a polytree, so a junction must be a
    **factor on its endpoint regions** and never a message channel. This
    is what a consumer needs to place that factor: where the junction attaches, whose transcript it
    belongs to, and how much of its own template remains either side.

    ⭐ **``strand`` is the TRANSCRIPT strand**, and it is the whole reason the accumulator does not store
    sense/antisense. The accumulator's ``sj_count`` columns are the **genome** strand the read aligned
    to; which transcript the molecule came from is a property of the annotation, stated here. "Sense
    derived from a junction's own strand" is exactly this join.

    ⚠ **``reach_lo``/``reach_hi`` are EXONIC**, unlike a contiguous boundary's (which are genomic spans, so
    that nascent RNA inside an intron is not declared impossible). Only a spliced molecule uses a
    junction, so what remains either side of it is exonic — see the module docstring. A reach of 0 is
    meaningful, not a sentinel.

    ⚠ ``lo``/``hi`` are genomically lower/higher, never donor/acceptor: boundaries store ``src < dst``, so on
    a NEG-strand junction the biological donor is on the right.
    """

    src_region: np.ndarray  # int64[J] — global region id; the intron STARTS where this region ends
    dst_region: np.ndarray  # int64[J] — global region id; the intron ENDS where this region begins
    strand: np.ndarray  # int8[J]  — the junction's own genomic strand == the TRANSCRIPT strand
    reach_lo: np.ndarray  # float64[J] — exonic reach on the junction's OWN strand
    reach_hi: np.ndarray  # float64[J]

    @property
    def n_junctions(self) -> int:
        return int(self.src_region.shape[0])


def build_junction_geometry_arrays(index) -> JunctionGeometry:
    """Build :class:`JunctionGeometry` for ``index``, in the accumulator's junction slot order.

    ⭐ **The slot order is not recomputed here.** It is read off
    :func:`build_junction_edge_arrays`, whose ``edge_row`` is already the ``edges_df`` row of each
    junction slot. That ordering is a byte-identity contract against the reference accumulator's
    ``Partition.from_region_bounds``, and two implementations of one ordering is how the same quantity comes to
    differ in two places — so there is one, and this joins through it.
    """
    rows = build_junction_edge_arrays(index).edge_row
    boundaries = index.edges_df
    strand = boundaries["strand"].to_numpy(np.int8)[rows]
    is_pos = strand == np.int8(Strand.POS)

    def reach(column_pos: str, column_neg: str) -> np.ndarray:
        # A junction row carries ONE strand, so only that strand's pair was populated by the builder;
        # selecting on the strand is what stops a NEG junction reading a POS column of zeros and
        # declaring itself opportunity-free.
        return np.where(
            is_pos,
            boundaries[column_pos].to_numpy(np.int64)[rows],
            boundaries[column_neg].to_numpy(np.int64)[rows],
        ).astype(np.float64)

    return JunctionGeometry(
        src_region=boundaries["src"].to_numpy(np.int64)[rows],
        dst_region=boundaries["dst"].to_numpy(np.int64)[rows],
        strand=strand,
        reach_lo=reach("reach_lo_pos", "reach_lo_neg"),
        reach_hi=reach("reach_hi_pos", "reach_hi_neg"),
    )


def is_terminus(flags: np.ndarray, strand: int) -> np.ndarray:
    """``TSS_s or TES_s`` — a transcript of strand ``s`` starts or ends here (graph doc §2.3)."""
    bits = (FLAG_TSS_POS | FLAG_TES_POS) if strand == Strand.POS else (FLAG_TSS_NEG | FLAG_TES_NEG)
    return (np.asarray(flags, dtype=np.uint16) & bits) != 0


def is_splice_site(flags: np.ndarray, strand: int) -> np.ndarray:
    """``DONOR_s or ACCEPTOR_s`` — a strand-``s`` intron begins or ends here.

    ⚠ NOT the complement of :func:`is_terminus`: the two are independent bits, and a position can be
    both — the case a 4-bit signature is structurally blind to. ⭐ Measured over all 1,043,595 human
    contiguous boundaries: **terminus-only 40.70 %, splice-only 58.31 %, BOTH 0.99 %** (10,337 boundaries), and
    **zero** carrying neither. So the both-bits case is real and worth handling, but it is 1 % of
    boundaries — not "the majority", as the design doc and an earlier draft of this docstring claimed.
    """
    bits = (
        (FLAG_DONOR_POS | FLAG_ACCEPTOR_POS)
        if strand == Strand.POS
        else (FLAG_DONOR_NEG | FLAG_ACCEPTOR_NEG)
    )
    return (np.asarray(flags, dtype=np.uint16) & bits) != 0


def load_regions(path) -> pd.DataFrame:
    df = pd.read_feather(str(path))
    missing = set(REGION_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(
            f"regions.feather at {path} is missing {sorted(missing)}. Rebuild the index."
        )
    return _coerce(df, REGION_COLUMNS, REGION_COLUMN_DTYPES)


def load_boundaries(path) -> pd.DataFrame:
    df = pd.read_feather(str(path))
    missing = set(BOUNDARY_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(
            f"edges.feather at {path} is missing {sorted(missing)}. Rebuild the index."
        )
    return _coerce(df, BOUNDARY_COLUMNS, BOUNDARY_COLUMN_DTYPES)


# ══════════════════════════════════════════════════════════════════════════════════════════════════
# THE TRANSCRIPT PATH — a transcript as an ordered walk over REGIONs, BOUNDARYs and SPLICE JUNCTIONs
# ══════════════════════════════════════════════════════════════════════════════════════════════════

#: The three kinds of step a transcript's path takes. ⭐ Deliberately the same ``(kind, obj_id)`` idiom
#: ``RegionChain`` already uses for the solve's slots — one addressing convention for one graph — with a
#: third kind, because a splice junction is neither a position nor an interval.
STEP_REGION = 0
STEP_BOUNDARY = 1
STEP_SPLICE_JUNCTION = 2


@dataclasses.dataclass(frozen=True, slots=True)
class TranscriptPath:
    """A CSR of every transcript's ordered walk through the graph.

    ``offsets[t]:offsets[t + 1]`` is transcript ``t``'s steps, and step ``s`` is
    ``(kind[s], obj_id[s])`` where ``obj_id`` indexes the REGION axis, the BOUNDARY axis or the SPLICE
    JUNCTION axis according to ``kind``.

    ⭐⭐ **THE STEPS ARE IN TRANSCRIPTION ORDER, 5' to 3'** — so a minus-strand transcript's steps run
    DESCENDING in genomic coordinate. ⚠ That is the one property a genomic-order implementation gets
    silently wrong: it is invisible to a consumer that treats a path as a SET or averages symmetrically
    over it, and wrong for every consumer that reads it as a sequence.

    ⭐ Verified on the shipped suite index against the independently-checked ``_transcript_region_incidence``
    (2026-08-12): **0 of 15,669 transcripts differ** on the region axis or the boundary axis, the splice
    steps number **45,609** — exactly the annotated (transcript, intron) count — and the ordering holds
    on **7,312/7,312** minus-strand and **8,357/8,357** plus-strand transcripts.
    """

    offsets: np.ndarray  # int64[n_transcripts + 1]
    kind: np.ndarray  # int8[n_steps]
    obj_id: np.ndarray  # int64[n_steps]

    @property
    def n_transcripts(self) -> int:
        return int(self.offsets.shape[0]) - 1

    def steps(self, t: int) -> tuple[np.ndarray, np.ndarray]:
        """``(kind, obj_id)`` for one transcript, in transcription order."""
        lo, hi = int(self.offsets[t]), int(self.offsets[t + 1])
        return self.kind[lo:hi], self.obj_id[lo:hi]


def _splice_junction_by_intron(index) -> dict[tuple, int]:
    """``(ref_name, intron_start, intron_end, strand) -> sj_id``, the ONLY admissible join key.

    ⛔ **NOT the flanking region pair, and not a row order.** The region pair happens to be unique on
    the shipped partition — 13,482 distinct pairs for 13,482 junctions — but only because every exon
    endpoint is forced to be a region region_bound; on a coarsened partition it collides (measured: 638 ambiguous
    pairs hiding 1,552 junctions). And ``sj.feather``'s row order is grouped ALPHABETICALLY by reference
    while the graph axis uses ``index.ref_names`` (FASTA order), which diverge on any genome carrying
    chr1/chr2/chr10.

    ⚠ The ``sj_id`` is a DENSE RANK over the junction boundaries, so it is a within-run join key and never a
    durable identifier: dropping one annotated intron renumbers almost every surviving slot.
    """
    ja = build_junction_edge_arrays(index)
    rows = index.edges_df.iloc[np.asarray(ja.edge_row, dtype=np.int64)]
    src = rows["src"].to_numpy(np.int64)
    dst = rows["dst"].to_numpy(np.int64)
    strand = rows["strand"].to_numpy(np.int64)
    regions = index.regions_df
    n_end = regions["end"].to_numpy(np.int64)
    n_start = regions["start"].to_numpy(np.int64)
    n_ref = regions["ref_name"].to_numpy()
    # ⭐ `src` is the genomically LOWER region on BOTH strands, so the intron is [end[src], start[dst]).
    # Never `donor`/`acceptor`: those are 5'/3' and therefore strand-dependent, and on this index they
    # would name the wrong end of 6,527 of 13,482 junctions.
    return {
        (str(n_ref[s]), int(n_end[s]), int(n_start[d]), int(st)): j
        for j, (s, d, st) in enumerate(zip(src, dst, strand, strict=True))
    }


def build_transcript_path(index, region_arrays) -> TranscriptPath:
    """Every transcript's ordered walk over REGIONs, BOUNDARYs and SPLICE JUNCTIONs.

    ⭐⭐ **WHAT A TRANSCRIPT'S PATH IS, AND WHAT IT DELIBERATELY EXCLUDES** — the include/exclude rule is
    the whole content of this function, and it is derived from what a fragment FROM this transcript can
    physically occupy:

    * **REGION** — every region an exon overlaps. A multi-exon transcript therefore takes its exonic
      regions and NOT its intronic ones: a mature molecule has no intronic bases. A single-exon
      transcript, and a synthetic span the index manufactured, take every region they cover.
    * **BOUNDARY** — a boundary the transcript crosses *contiguously*, i.e. one strictly INTERIOR to a
      single exon. ⭐ An exon-interior boundary that merely marks a signature change (an antisense
      feature overlapping on the other strand) IS crossed and IS included.
    * ⛔ **The transcript's OUTER boundaries are excluded** — its TSS and TES. Nothing from this
      transcript crosses them, because the molecule ends there.
    * ⛔ **A splice donor/acceptor boundary is excluded as a BOUNDARY step** and appears as the
      SPLICE JUNCTION step instead. A molecule at that position either splices (the sj) or reads
      through (a different molecule, and a different transcript's path).
    * **SPLICE JUNCTION** — one per adjacent exon pair, resolved to its exact ``sj_id`` by intron
      coordinates.

    ⚠ **Steps are emitted in TRANSCRIPTION order.** The regions and boundaries of one exon run
    ascending in genomic coordinate; the whole path is then reversed for a minus-strand transcript.

    ⛔ Annotation-only and sample-independent — it reads the index and the partition and nothing else,
    so it is valid for every condition and could be precomputed at index build.
    """
    import os

    from ..types import IntervalType
    from .region_arrays import region_right_boundary

    starts = np.asarray(region_arrays.start, dtype=np.int64)
    ends = np.asarray(region_arrays.end, dtype=np.int64)
    ref_off = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    right_boundary = region_right_boundary(np.asarray(region_arrays.ref_id))
    name_to_id = index.ref_name_to_id
    sj_of_intron = _splice_junction_by_intron(index)

    n_t = int(index.num_transcripts)
    per_t: list[list[tuple[int, int]]] = [[] for _ in range(n_t)]
    # ⛔ Resolved BEFORE the exon walk, not after: the splice-junction join key contains the strand, so
    # a strand read later is a strand read as 0 — which silently resolves no junction at all.
    strand_of = np.zeros(n_t, dtype=np.int64)
    tdf = index.t_df
    if tdf is not None and "strand" in tdf.columns:
        strand_of[tdf["t_index"].to_numpy(np.int64)] = tdf["strand"].to_numpy().astype(np.int64)

    def _exon_steps(ref_name, a: int, b: int):
        """One exon's REGION and interior BOUNDARY steps, interleaved, genomically ascending."""
        rid = name_to_id.get(str(ref_name))
        if rid is None:
            return None
        lo0, hi0 = int(ref_off[rid]), int(ref_off[rid + 1])
        lo = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        hi = lo0 + int(np.searchsorted(starts[lo0:hi0], b, side="left"))
        if hi <= lo:
            return None
        out: list[tuple[int, int]] = [(STEP_REGION, lo)]
        for r in range(lo + 1, hi):
            # the boundary between region r-1 and r is INTERIOR to this exon, so it is crossed
            out.append((STEP_BOUNDARY, int(right_boundary[r - 1])))
            out.append((STEP_REGION, r))
        return out, lo, hi - 1

    iv = pd.read_feather(os.path.join(index.index_dir, "intervals.feather"))
    ex = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]
    ex = ex.sort_values(["t_index", "start"], kind="stable")
    seen: set[int] = set()
    prev_end: dict[int, int] = {}  # t -> the previous exon's genomic END (the intron's low side)
    unresolved: list[tuple[int, tuple]] = []

    for t, ref_name, a, b in zip(ex["t_index"], ex["ref"], ex["start"], ex["end"], strict=True):
        t = int(t)
        seen.add(t)
        res = _exon_steps(ref_name, int(a), int(b))
        if res is None:
            continue
        steps, _lo, _hi = res
        if t in prev_end:
            key = (str(ref_name), prev_end[t], int(a), int(strand_of[t]))
            sj = sj_of_intron.get(key)
            if sj is None:
                unresolved.append((t, key))
            else:
                per_t[t].append((STEP_SPLICE_JUNCTION, int(sj)))
        per_t[t].extend(steps)
        prev_end[t] = int(b)

    # ⛔⛔ A TRANSCRIPT'S ANNOTATED INTRON THAT RESOLVES TO NO SLOT IS A DEFECT, NOT A GAP TO SKIP.
    # The junction key is derived from REGION boundaries (`end[src]`, `start[dst]`), which equal the
    # intron's coordinates only because the partition region_bounds at every exon endpoint — measured 0 of 45,609
    # violations on the shipped index. That invariant is ASSUMED by the derivation, so it is asserted
    # here: were it ever to break, every affected transcript would silently lose a step and its path
    # would read as a shorter, well-formed walk. A silent wrong answer is the one outcome to refuse.
    if unresolved:
        t0, k0 = unresolved[0]
        raise ValueError(
            f"{len(unresolved)} annotated intron(s) resolved to no splice-junction slot; first is "
            f"transcript {t0} intron {k0!r}. The junction axis is keyed on region boundaries, so this "
            f"means an exon endpoint is not a region region_bound — the index and the partition disagree."
        )

    if tdf is not None and "is_synthetic" in tdf.columns:
        syn = tdf[tdf["is_synthetic"].to_numpy(dtype=bool)]
        for t, ref_name, a, b in zip(
            syn["t_index"], syn["ref"], syn["start"], syn["end"], strict=True
        ):
            if int(t) in seen:
                continue
            res = _exon_steps(ref_name, int(a), int(b))
            if res is not None:
                per_t[int(t)].extend(res[0])

    # ⭐ TRANSCRIPTION ORDER. Everything above is genomically ascending; a minus-strand transcript is
    # transcribed from its high coordinate down, so its whole path reverses — steps AND their order.
    for t in np.flatnonzero(strand_of == int(Strand.NEG)):
        per_t[int(t)].reverse()

    counts = np.fromiter((len(p) for p in per_t), dtype=np.int64, count=n_t)
    offsets = np.zeros(n_t + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    total = int(offsets[-1])
    kind = np.empty(total, dtype=np.int8)
    obj_id = np.empty(total, dtype=np.int64)
    i = 0
    for p in per_t:
        for k, o in p:
            kind[i] = k
            obj_id[i] = o
            i += 1
    return TranscriptPath(offsets=offsets, kind=kind, obj_id=obj_id)
