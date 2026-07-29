"""rigel.calibration.splice_graph — the v8 SPLICE GRAPH: nodes + edges.

The calibration partition: what the accumulator deposits into, and the structure the solver reads.

    Design: ``docs/CARRY_FORWARD.md``   ·   Consumer: ``docs/CARRY_FORWARD.md``
    Implementation plan + measured amendments: ``docs/CARRY_FORWARD.md``

**A node** is a genomic interval; nodes tile each reference and are numbered in genomic order.
**An edge** is a transition, always ``src < dst`` so genomic order is a topological order:

* ``CONTIGUOUS(i, i+1)`` — the genomic point ``end(i) == start(i+1)``. Carries the 8 structural
  :ref:`flag bits <flags>`.
* ``JUNCTION(donor, acceptor, strand)`` — one per distinct annotated intron.

**Every distinct exon endpoint is a cut, and adjacent equal-signature nodes are NOT merged.** The merge is
what the predecessor partition did and what made it blind to transcript termini: ⭐ **53.4 %** of real human
transcript termini (232,451 of 435,291) fall strictly *inside* a merged region and vanish from the partition
entirely. Measured consequence: 752,654 merged regions → **1,043,881** nodes (+38.7 %), median 151 bp,
15,687 of length 1, plus **404,168** junction edges.

⚠ The often-quoted **59.5 %** is that same statistic computed under the ``~is_synthetic & ~is_nrna``
transcript filter — i.e. with 26,475 real single-exon transcripts excluded, which is the very bug
described below. Under the shipped filter it is **53.4 %**.

⭐ **ONE transcript filter: ``~is_synthetic``.** Every manufactured nRNA span row is ``is_synthetic``,
so that single predicate already excludes all of them — which is what P1G_SCOPE §5's "211/211 signature
false positives are nRNA-span edges" actually requires.

⚠ **It was briefly TWO filters, and the second one was a measured bug** (plan F1 proposed
``~is_synthetic & ~is_nrna`` for the flags and reaches, reasoning that "an nRNA span's ends are not real
transcript termini"). The reasoning is right; the predicate is not. On a **non-synthetic** row
``is_nrna`` does not mean "manufactured span" — it means **this real transcript is single-exon, so its
mature and nascent forms coincide**. Measured on the human annotation: all **26,475** rows that are
``is_nrna & ~is_synthetic`` have ``n_exons == 1`` and **none** is a ``RIGEL_NRNA_*`` row, so the extra
clause deleted the TSS/TES of 26,475 real transcripts — **52,104 distinct terminus positions** — which
is exactly the visibility v8 exists to buy.

.. _flags:

**Structural flags** (contiguous edges): ``TSS_s``/``TES_s``/``DONOR_s``/``ACCEPTOR_s`` for ``s ∈ {+,−}``.
They are NOT mutually exclusive — a position that is both a terminus for one transcript and a splice site
for another is exactly the case the 4-bit signature is structurally blind to.

**Reaches** — how many bases of the molecule's own sequence remain either side of an edge. A mature
molecule must fit inside its transcript, so this is what makes the crossing divisor taper near a terminus
(v5 §4.2/§4.3; measured worth 11.0 % of the mature opportunity genome-wide).

⚠ **Named ``lo``/``hi``, NOT ``donor``/``acceptor``.** Edges store ``src < dst``, so ``src`` is genomically
LEFT whatever the strand — but for a NEG-strand junction the biological donor is on the RIGHT. The
divisor is symmetric in the two reaches so no number changes, but ``reach_donor`` would mislabel roughly
half of the 404,168 junctions. ``lo``/``hi`` means genomically lower/higher and cannot be misread.

⚠ **Per STRAND**, because the mature-crossing gate is per strand: at an AMBIG seam a ``+`` and a ``−``
transcript have different reaches, and a strand-agnostic maximum would over-state one of them — on exactly
the population that carries 31.6 % of the calibration suite's error mass.

**Maximal over isoforms, independently per side** (decision D2). A genomic position usually belongs to
several transcripts that disagree about where the molecule ends, and the isoform abundances are precisely
what calibration does not know. Measured (plan C3): against the alternative of taking the best *realisable*
pair from a single isoform, the two agree exactly on **93.9 %** of disagreeing junctions with a mean ratio
of **0.9989**, so the simple independent maximum is used. It is one-sided — it over-states opportunity,
hence under-states ``ρ_mature``, hence over-states ``f_g``.

**A reach of 0 is meaningful, not a sentinel:** no strand-``s`` mature molecule crosses this edge, so the
opportunity is genuinely zero and the trapezoid returns 0 for free.
"""

from __future__ import annotations

import warnings
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
    "NODE_COLUMNS",
    "NODE_COLUMN_DTYPES",
    "EDGE_COLUMNS",
    "EDGE_COLUMN_DTYPES",
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
    "build_node_partition_arrays",
    "build_boundary_flags_array",
    "is_terminus",
    "is_splice_site",
    "validate_graph",
    "load_nodes",
    "load_edges",
]

NODE_COLUMNS = ["node_id", "ref_name", "start", "end", "length", "signature"]
NODE_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "node_id": np.int64,
    "start": np.int64,
    "end": np.int64,
    "length": np.int64,
    "signature": np.uint8,
}

EMPTY_I64 = np.zeros(0, np.int64)

EDGE_KIND_CONTIGUOUS = 0
EDGE_KIND_JUNCTION = 1

EDGE_COLUMNS = [
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
EDGE_COLUMN_DTYPES: dict[str, type | np.dtype] = {
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

# Structural bits on a CONTIGUOUS edge, at its genomic position. Not mutually exclusive.
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
    """Build ``(nodes_df, edges_df)`` — the v8 splice graph. Fully vectorised; no per-node Python object.

    Deterministic by construction: ``np.unique`` sorts, node ids are assigned by position, edges by
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
    node_base = 0
    for name, length in reflen.items():
        if length == 0:
            continue
        rows = by_ref.get(name, EMPTY_I64)
        i_rows = intr_by_ref.get(name, EMPTY_I64)
        # ── 1. the cut set ────────────────────────────────────────────────────────────────────────
        cuts = np.unique(
            np.concatenate([ex.start[rows], ex.end[rows], np.array([0, length], np.int64)])
        )
        starts, ends = cuts[:-1], cuts[1:]
        n = starts.shape[0]

        # ── 2. the 4-bit signature, by difference arrays over the node index ──────────────────────
        sig = np.zeros(n, np.uint8)
        for bit, iv_s, iv_e in _signature_intervals(ex, rows, intr, i_rows):
            if iv_s.size == 0:
                continue
            d = np.zeros(n + 1, np.int32)
            np.add.at(d, np.searchsorted(cuts, iv_s), 1)
            np.add.at(d, np.searchsorted(cuts, iv_e), -1)
            sig |= np.where(np.cumsum(d)[:-1] > 0, bit, 0).astype(np.uint8)

        n_rows.append(
            pd.DataFrame(
                {
                    "node_id": np.arange(node_base, node_base + n, dtype=np.int64),
                    "ref_name": name,
                    "start": starts,
                    "end": ends,
                    "length": ends - starts,
                    "signature": sig,
                }
            )
        )

        # ── 3. contiguous edges: flags + per-strand reaches at each interior interface ────────────
        pos = cuts[1:-1]  # the n-1 interior interfaces
        flags = np.zeros(pos.shape[0], np.uint16)
        if pos.size:
            for arr, bit in _terminus_and_splice_events(ex, rows, intr, i_rows):
                if arr.size:
                    hit = np.searchsorted(pos, arr)
                    ok = (hit < pos.size) & (pos[np.clip(hit, 0, max(pos.size - 1, 0))] == arr)
                    np.bitwise_or.at(flags, hit[ok], bit)
        c_reach = _contiguous_reaches(ex, rows, pos)

        # ── 4. junction edges: one per distinct (intron_start, intron_end, strand) ────────────────
        j_src, j_dst, j_strand, j_reach = _junction_edges(cuts, intr, i_rows, name)

        src = np.concatenate([np.arange(n - 1, dtype=np.int64), j_src]) + node_base
        dst = np.concatenate([np.arange(1, n, dtype=np.int64), j_dst]) + node_base
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
        node_base += n

    nodes_df = (
        pd.concat(n_rows, ignore_index=True) if n_rows else pd.DataFrame(columns=NODE_COLUMNS)
    )
    edges_df = (
        pd.concat(e_rows, ignore_index=True) if e_rows else pd.DataFrame(columns=EDGE_COLUMNS[1:])
    )
    if len(edges_df):
        # ⚠ sorted by (src, kind, dst, STRAND). The design doc says (src, kind, dst), but that is not a
        # TOTAL order: two strand-coincident junctions (§3.3, G18) share all three and differ only in
        # strand, so the order would be ambiguous and determinism would depend on the input order.
        # GENCODE contains zero such junctions — which is exactly why only a synthetic case can catch it.
        # Adding strand REFINES the documented order, so the "out-edges of a node are contiguous" CSR
        # contract is unchanged.
        edges_df = edges_df.sort_values(
            ["src", "kind", "dst", "strand"], kind="stable"
        ).reset_index(drop=True)
    edges_df.insert(0, "edge_id", np.arange(len(edges_df), dtype=np.int64))
    return _coerce(nodes_df, NODE_COLUMNS, NODE_COLUMN_DTYPES), _coerce(
        edges_df, EDGE_COLUMNS, EDGE_COLUMN_DTYPES
    )


def _signature_intervals(ev: _Exons, ei, ev_i, ii):
    """The four ``(bit, starts, ends)`` interval sets that define a node's 4-bit signature.

    ONE definition, consumed by the builder (cumulative difference arrays over the node index) and
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
    cut = np.flatnonzero(ref_arr[1:] != ref_arr[:-1]) + 1
    return {
        str(ref_arr[a]): np.arange(a, b, dtype=np.int64)
        for a, b in zip(np.r_[0, cut], np.r_[cut, n])
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
            # a transcript's own outer edges: exonic offset 0 and offset == total
            first = m[ma.before[m] == 0]
            last = m[(ma.before[m] + (ma.end - ma.start)[m]) == ma.total[m]]
            lo_edge, hi_edge = ma.start[first], ma.end[last]
            # 5' terminus is the genomically-low edge on +, the high edge on −
            tss = lo_edge if st == Strand.POS else hi_edge
            tes = hi_edge if st == Strand.POS else lo_edge
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


def _contiguous_reaches(ma: _Exons, mj, pos):
    """Per-strand maximal reach either side of each interior interface, over MATURE transcripts.

    A strand-``s`` mature molecule crosses interface ``p`` contiguously iff one of its exons contains ``p``
    strictly inside. Its reach on the low side is (exonic bases before that exon) + (p − exon.start); the
    high side is the remainder. Maximal over transcripts, independently per side (D2).
    """
    out = {
        k: np.zeros(pos.shape[0], np.int64)
        for k in ("reach_lo_pos", "reach_hi_pos", "reach_lo_neg", "reach_hi_neg")
    }
    if pos.size == 0 or mj.size == 0:
        return out
    lo_i = np.searchsorted(pos, ma.start[mj], side="right")  # first interface > exon.start
    hi_i = np.searchsorted(pos, ma.end[mj], side="left")  # first interface >= exon.end
    cnt = np.maximum(hi_i - lo_i, 0)
    if int(cnt.sum()) == 0:
        return out
    ex = np.repeat(np.arange(mj.size), cnt)  # which exon row
    off = np.arange(int(cnt.sum())) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    slot = np.repeat(lo_i, cnt) + off  # which interface
    p = pos[slot]
    lo = ma.before[mj][ex] + (p - ma.start[mj][ex])
    hi = ma.total[mj][ex] - lo
    st = ma.strand[mj][ex]
    for s, kl, kh in (
        (Strand.POS, "reach_lo_pos", "reach_hi_pos"),
        (Strand.NEG, "reach_lo_neg", "reach_hi_neg"),
    ):
        m = st == s
        if m.any():
            np.maximum.at(out[kl], slot[m], lo[m])
            np.maximum.at(out[kh], slot[m], hi[m])
    return out


def _junction_edges(cuts, ma_i, mi_, ref_name):
    """Distinct junction edges for one reference: ``(src, dst, strand, reach dict)``, node-local ids."""
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
    src = np.searchsorted(cuts, ua) - 1  # the node ENDING at intron_start
    dst = np.searchsorted(cuts, ub)  # the node STARTING at intron_end
    if not (np.all(cuts[src + 1] == ua) and np.all(cuts[dst] == ub)):
        raise ValueError(f"ref {ref_name!r}: a junction endpoint is not a node interface (I5)")
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


def validate_graph(nodes_df, edges_df, ref_lengths: Mapping[str, int], transcripts=None) -> None:
    """Assert invariants I1–I13. Raises :class:`ValueError` on the first violation.

    ⚠ **``transcripts`` gates four invariants, not one**: I3b (the signature, recomputed), I4 (the
    interfaces ARE the events), I11 (every transcript walks) and I13 (the flags ARE the events).
    Without it only the graph-internal checks run — I1/I2/I5–I9/I12 — which cannot detect a
    ``signature`` or ``flags`` column that has drifted from the annotation. The build passes
    transcripts; the load path does not, because reconstructing them costs ~3 s at human scale.
    """
    reflen = {str(k): int(v) for k, v in ref_lengths.items()}
    nid = nodes_df["node_id"].to_numpy(np.int64)
    start = nodes_df["start"].to_numpy(np.int64)
    end = nodes_df["end"].to_numpy(np.int64)
    length = nodes_df["length"].to_numpy(np.int64)
    sig = nodes_df["signature"].to_numpy(np.uint8)
    ref = nodes_df["ref_name"].astype(str).to_numpy()

    # I2 — ids are 0..n-1 in row order
    if not np.array_equal(nid, np.arange(nid.size, dtype=np.int64)):
        raise ValueError("I2: node_id is not 0..n-1 in row order. Rebuild the index.")
    # I3a — signature range (I3b, the independent recomputation, needs the transcripts)
    if nid.size and int(sig.max()) >= N_SIGNATURES:
        raise ValueError(f"I3: signature {int(sig.max())} out of range. Rebuild the index.")
    if not np.array_equal(length, end - start):
        raise ValueError("I1: stored node length != end - start. Rebuild the index.")
    if nid.size and int(length.min()) <= 0:
        raise ValueError("I1: a node has non-positive length. Rebuild the index.")

    # I1 — nodes tile each reference exactly. ⚠ Reference slices are computed ONCE from the run
    # boundaries; `np.flatnonzero(ref == name)` inside the loop is 286 full-array string comparisons
    # over 1.04 M rows and cost 3.5 s on its own.
    slices = _ref_slices(ref)
    for name, L in reflen.items():
        sl = slices.get(name)
        if sl is None:
            if L > 0:
                raise ValueError(f"I1: reference {name!r} (length {L}) has no nodes.")
            continue
        lo, hi = int(sl[0]), int(sl[-1]) + 1
        if start[lo] != 0 or end[hi - 1] != L:
            raise ValueError(
                f"I1: reference {name!r} spans [{start[lo]}, {end[hi - 1]}), expected [0, {L})."
            )
        if hi - lo > 1 and not np.array_equal(end[lo : hi - 1], start[lo + 1 : hi]):
            raise ValueError(f"I1: reference {name!r} has a gap or overlap between nodes.")
    if len(slices) != sum(1 for _n, L in reflen.items() if L > 0):
        raise ValueError("I1: node rows are not grouped contiguously by reference.")

    src = edges_df["src"].to_numpy(np.int64)
    dst = edges_df["dst"].to_numpy(np.int64)
    kind = edges_df["kind"].to_numpy(np.uint8)
    estrand = edges_df["strand"].to_numpy(np.int8)

    # I8 — src < dst on EVERY edge  ⇒ genomic order is a topological order
    if src.size and not np.all(src < dst):
        raise ValueError("I8: an edge has src >= dst; genomic order is not a topological order.")
    if src.size and not np.all(ref[src] == ref[dst]):
        raise ValueError("I5: an edge crosses references.")
    # I12 — row order is the id order
    if not np.array_equal(edges_df["edge_id"].to_numpy(np.int64), np.arange(src.size)):
        raise ValueError("I12: edge_id is not the row index.")
    if src.size > 1:
        # (src, kind, dst, strand) packed into one strictly-increasing int64 key. STRAND is part of the
        # key because it is part of the order (see the sort in build_splice_graph): without it two
        # strand-coincident junctions collide and read as a duplicate.
        n_nodes_ = max(nid.size, 1)
        key = ((src * 2 + kind.astype(np.int64)) * n_nodes_ + dst) * 4 + estrand.astype(np.int64)
        if not np.all(np.diff(key) > 0):
            raise ValueError(
                "I12: edges are not sorted by (src, kind, dst, strand), or a duplicate exists."
            )

    # I7 — contiguous edges are exactly {(i, i+1)} within each reference
    c = kind == EDGE_KIND_CONTIGUOUS
    if not np.all(dst[c] == src[c] + 1):
        raise ValueError("I7: a CONTIGUOUS edge is not between adjacent nodes.")
    # ⚠ `sum(... (ref == name).any())` here was 286 full-array object-dtype comparisons over 1.04 M
    # rows — 1.68 s of a 1.91 s load-time validate_graph, i.e. 88 % of it — and it recomputed a
    # number `slices` already holds. This is the same trap the I1 block above documents.
    n_expected = nid.size - len(slices)
    if int(c.sum()) != n_expected:
        raise ValueError(f"I7: {int(c.sum())} contiguous edges, expected {n_expected}.")
    if np.any(estrand[c] != 0):
        raise ValueError("I7: a CONTIGUOUS edge carries a strand.")

    # I5/I6 — junctions land on interfaces, one row per distinct (ref, donor, acceptor, strand)
    j = kind == EDGE_KIND_JUNCTION
    if j.any():
        if not np.all(np.isin(estrand[j], [Strand.POS, Strand.NEG])):
            raise ValueError("I6: a JUNCTION edge has an invalid strand.")
        # ⚠ STRAND-COINCIDENT junctions are BIOLOGICALLY IMPOSSIBLE and are reported, not tolerated
        # silently. Splicing reads a non-palindromic motif: the reverse complement of a GT..AG intron
        # begins CT, so the same interval cannot be a valid intron on both strands. Measured: ZERO in
        # GENCODE. One in a GTF therefore means the ANNOTATION is wrong (a simulator, or a strand
        # column filled in by hand), and silence would let it propagate into every downstream density.
        # It is a WARNING rather than an error on purpose: the graph handles the case correctly (the
        # sort key carries strand, so the two edges stay distinct and ordered), and the G18 test builds
        # exactly this input to prove that. Raising would make the guard untestable.
        pair = src[j] * max(nid.size, 1) + dst[j]
        _u, _cnt = np.unique(pair, return_counts=True)
        n_coincident = int((_cnt > 1).sum())
        if n_coincident:
            warnings.warn(
                f"{n_coincident} strand-coincident splice junction(s): the same (donor, acceptor) "
                f"annotated on BOTH strands. This is biologically impossible — splice motifs are "
                f"non-palindromic (a GT..AG intron reverse-complements to CT..AC) — so the source "
                f"annotation is likely wrong. The graph keeps them as distinct edges and is correct "
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
        _validate_against_transcripts(transcripts, reflen, nodes_df, edges_df)


def _validate_against_transcripts(transcripts, reflen, nodes_df, edges_df) -> None:
    """I3b (the signature, recomputed) + I4 (interfaces ARE the events) + I11 (every transcript
    walks) + I13 (the flags ARE the events), vectorised.

    ⚠ Written without a per-transcript Python loop on purpose: the loop form cost 10.4 s on the human
    annotation against §8's 5 s budget, and it is the check most likely to be skipped when it is slow —
    which would be exactly backwards, since I11 is the one that catches real bugs.
    """
    ref = nodes_df["ref_name"].astype(str).to_numpy()
    start = nodes_df["start"].to_numpy(np.int64)
    end = nodes_df["end"].to_numpy(np.int64)
    sig = nodes_df["signature"].to_numpy(np.uint8)
    src = edges_df["src"].to_numpy(np.int64)
    dst = edges_df["dst"].to_numpy(np.int64)
    kind = edges_df["kind"].to_numpy(np.uint8)
    estrand = edges_df["strand"].to_numpy(np.int8)

    ex = _Exons(transcripts, lambda tx: _is_real(tx, reflen))
    by_ref = _ref_slices(ex.ref)
    intr = _introns_of(ex)
    intr_by_ref = _ref_slices(intr[0])

    cm = kind == EDGE_KIND_CONTIGUOUS
    csrc = src[cm]  # already ascending: edges sort by (src, kind, dst, strand)
    cflags = edges_df["flags"].to_numpy(np.uint16)[cm]

    n_nodes = ref.size
    jm = kind == EDGE_KIND_JUNCTION
    jkey = np.sort((src[jm] * n_nodes + dst[jm]) * 4 + estrand[jm].astype(np.int64))

    slices = _ref_slices(ref)
    for name, L in reflen.items():
        m = slices.get(name)
        if m is None:
            continue
        base, cuts = int(m[0]), np.concatenate([start[m], end[m][-1:]])

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

        # ⭐ I3b — recompute the signature from each node's MIDPOINT, by direct interval
        # containment. Deliberately a DIFFERENT algorithm from the builder's cumulative-difference
        # sweep, so the two can only agree by both being right. A node is homogeneous by
        # construction (every interval endpoint is a cut), so its midpoint decides it.
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
                f"I3: reference {name!r} node {int(m[k])} [{int(start[m][k])}, "
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
            e_of_node = np.searchsorted(csrc, m[:-1])
            hit = (e_of_node < csrc.size) & (
                csrc[np.clip(e_of_node, 0, max(csrc.size - 1, 0))] == m[:-1]
            )
            at[hit] = cflags[e_of_node[hit]]
            for arr, bit in _events_independently(ex, rows, intr, i_rows):
                want = np.isin(interior, arr)
                got = (at & bit) != 0
                if not np.array_equal(want, got):
                    k = int(np.flatnonzero(want != got)[0])
                    raise ValueError(
                        f"I13: reference {name!r} position {int(interior[k])} "
                        f"{'is missing' if want[k] else 'wrongly carries'} flag bit {int(bit):#06x}."
                    )

        # I11a — every mature exon edge lands exactly on a node interface
        if rows.size:
            for arr, what in ((ex.start[rows], "start"), (ex.end[rows], "end")):
                hit = np.searchsorted(cuts, arr)
                bad = (hit >= cuts.size) | (cuts[np.clip(hit, 0, cuts.size - 1)] != arr)
                if bad.any():
                    raise ValueError(
                        f"I11: reference {name!r} has {int(bad.sum())} exon {what}s that are not node "
                        f"interfaces (first at {int(arr[bad][0])})."
                    )

        # I11b — every mature intron is realised by a JUNCTION edge on that strand
        if i_rows.size:
            a, b, st = intr[1][i_rows], intr[2][i_rows], intr[3][i_rows].astype(np.int64)
            js = base + np.searchsorted(cuts, a) - 1
            jd = base + np.searchsorted(cuts, b)
            want = (js * n_nodes + jd) * 4 + st
            if jkey.size == 0:
                miss = np.ones(want.shape, bool)
            else:
                pos = np.searchsorted(jkey, want)
                miss = (pos >= jkey.size) | (jkey[np.clip(pos, 0, jkey.size - 1)] != want)
            if miss.any():
                k = int(np.flatnonzero(miss)[0])
                raise ValueError(
                    f"I11: reference {name!r} has {int(miss.sum())} mature introns with no JUNCTION "
                    f"edge (first [{int(a[k])}, {int(b[k])}) strand {int(st[k])})."
                )


# ---------------------------------------------------------------------------
# The scanner adapter + load
# ---------------------------------------------------------------------------


def build_node_partition_arrays(index) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flatten the graph into the ``BamScanner.set_regions`` ABI — THE partition the accumulator
    deposits into.

    For each reference in ``index.ref_names`` the cut positions are ``[n0.start, n0.end, n1.end, …]``
    over that reference's nodes (which tile it contiguously), so ``k`` nodes give ``k + 1`` positions.
    References with zero nodes contribute none.

    :meth:`RegionArrays.from_index <rigel.calibration.region_arrays.RegionArrays.from_index>` reads
    the same frame; that is what keeps the calibration geometry and the scanner from addressing
    different partitions.

    Returns ``(cut_positions int64[P], ref_pos_offsets int64[n_refs+1], node_types uint8[N])``, where
    ``node_types`` is the coarse type (0=intergenic, 1=intron, 2=exon) aligned 1:1 with the nodes —
    the gDNA FL-pool axis.
    """
    nodes_df = index.nodes_df
    ref_names = index.ref_names
    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in nodes_df.groupby("ref_name", sort=False)
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
        cuts = np.empty(len(starts) + 1, dtype=np.int64)
        cuts[:-1] = starts
        cuts[-1] = ends[-1]
        positions.append(cuts)
        types.append(coarse_type_array(grp["signature"].to_numpy()))
        offsets[i + 1] = offsets[i] + cuts.shape[0]
    return (
        np.concatenate(positions) if positions else np.empty(0, dtype=np.int64),
        offsets,
        np.concatenate(types) if types else np.empty(0, dtype=np.uint8),
    )


def build_boundary_flags_array(index) -> np.ndarray:
    """The structural flags, re-indexed onto **the accumulator's boundary axis**.

    ⚠ **The two index spaces are off by one per reference, and in opposite directions.** A reference
    with ``k`` nodes has ``k − 1`` contiguous edges (one per interior interface) but ``k + 1``
    accumulator boundary slots (the interfaces *plus* the two reference terminals). So::

        slot 0        -> 0   (reference start terminal: not a transition, no edge exists)
        slot 1..k-1   -> the contiguous edge between nodes (i-1, i)
        slot k        -> 0   (reference end terminal)

    The two zeroed terminals are not padding — they are genuinely flag-less. A terminal boundary
    never carried a deposit either (``00_design.md`` §6 invariant 4), so a zero there is the correct
    answer rather than a missing one.

    Returns ``uint16[B]`` with ``B == build_node_partition_arrays(index)[1][-1]``, aligned element
    for element with the payload's ``ref_boundary_offsets``.
    """
    nodes_df, edges_df = index.nodes_df, index.edges_df
    # A contiguous edge IS the interface to the right of its src node, so key the flags by src.
    # Junction edges carry no flags (they are not a genomic position) and are excluded.
    contiguous = edges_df["kind"].to_numpy(np.uint8) == EDGE_KIND_CONTIGUOUS
    by_node = np.zeros(len(nodes_df), dtype=np.uint16)
    by_node[edges_df["src"].to_numpy(np.int64)[contiguous]] = edges_df["flags"].to_numpy(np.uint16)[
        contiguous
    ]

    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in nodes_df.groupby("ref_name", sort=False)
    }
    out: list[np.ndarray] = []
    for ref in index.ref_names:
        grp = by_ref.get(ref)
        if grp is None or len(grp) == 0:
            continue
        ids = grp.index.to_numpy(np.int64)  # == node_id (I2), and contiguous within a reference
        k = ids.shape[0]
        slots = np.zeros(k + 1, dtype=np.uint16)
        slots[1:k] = by_node[ids[: k - 1]]
        out.append(slots)
    return np.concatenate(out) if out else np.zeros(0, dtype=np.uint16)


def is_terminus(flags: np.ndarray, strand: int) -> np.ndarray:
    """``TSS_s or TES_s`` — a transcript of strand ``s`` starts or ends here (graph doc §2.3)."""
    bits = (FLAG_TSS_POS | FLAG_TES_POS) if strand == Strand.POS else (FLAG_TSS_NEG | FLAG_TES_NEG)
    return (np.asarray(flags, dtype=np.uint16) & bits) != 0


def is_splice_site(flags: np.ndarray, strand: int) -> np.ndarray:
    """``DONOR_s or ACCEPTOR_s`` — a strand-``s`` intron begins or ends here.

    ⚠ NOT the complement of :func:`is_terminus`: the two are independent bits, and a position can be
    both — the case a 4-bit signature is structurally blind to. ⭐ Measured over all 1,043,595 human
    contiguous edges: **terminus-only 40.70 %, splice-only 58.31 %, BOTH 0.99 %** (10,337 edges), and
    **zero** carrying neither. So the both-bits case is real and worth handling, but it is 1 % of
    seams — not "the majority", as the design doc and an earlier draft of this docstring claimed.
    """
    bits = (
        (FLAG_DONOR_POS | FLAG_ACCEPTOR_POS)
        if strand == Strand.POS
        else (FLAG_DONOR_NEG | FLAG_ACCEPTOR_NEG)
    )
    return (np.asarray(flags, dtype=np.uint16) & bits) != 0


def load_nodes(path) -> pd.DataFrame:
    df = pd.read_feather(str(path))
    missing = set(NODE_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(
            f"nodes.feather at {path} is missing {sorted(missing)}. Rebuild the index."
        )
    return _coerce(df, NODE_COLUMNS, NODE_COLUMN_DTYPES)


def load_edges(path) -> pd.DataFrame:
    df = pd.read_feather(str(path))
    missing = set(EDGE_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(
            f"edges.feather at {path} is missing {sorted(missing)}. Rebuild the index."
        )
    return _coerce(df, EDGE_COLUMNS, EDGE_COLUMN_DTYPES)
