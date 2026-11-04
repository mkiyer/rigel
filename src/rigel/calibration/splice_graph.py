"""The splice graph: the calibration partition of regions and boundaries.

Builds the structure the accumulator deposits into and the solver reads, validates it against the annotation it came
from, and provides the array factories the rest of calibration addresses it with (the region partition, the sj arrays
and geometry, the boundary flags). Layer 1, the payload view.

A region is a genomic interval; regions tile each reference in genomic order. A boundary is a transition, always
``src < dst``, so genomic order is a topological order:

* ``CONTIGUOUS(i, i+1)`` — the genomic point ``end(i) == start(i+1)``, carrying the 8 structural flag bits
  ``TSS_s`` / ``TES_s`` / ``DONOR_s`` / ``ACCEPTOR_s`` for ``s`` in ``{+, -}``. The bits are not mutually exclusive:
  one position can be a terminus of one transcript and a splice site of another.
* ``SJ(donor, acceptor, strand)`` — one per distinct annotated intron, carrying its exonic reaches
  (:class:`SpliceJunctionGeometry`); a contiguous boundary's reach columns hold 0.

Every distinct exon endpoint is a region bound, and adjacent equal-signature regions are not merged, because a merged
region hides every transcript terminus that falls strictly inside it. The transcript filter is :func:`_is_real`
(``~is_synthetic``), shared by the builder and the validator; a real single-exon transcript passes it, since its
termini are real.
"""

from __future__ import annotations

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
    "EDGE_KIND_SJ",
    "FLAG_TSS_POS",
    "FLAG_TSS_NEG",
    "FLAG_TES_POS",
    "FLAG_TES_NEG",
    "FLAG_DONOR_POS",
    "FLAG_DONOR_NEG",
    "FLAG_ACCEPTOR_POS",
    "FLAG_ACCEPTOR_NEG",
    "FLAG_TERMINUS",
    "FLAG_JUNCTION",
    "build_splice_graph",
    "build_region_partition_arrays",
    "build_boundary_flags_array",
    "build_sj_arrays",
    "build_sj_geometry_arrays",
    "SpliceJunctionArrays",
    "SpliceJunctionGeometry",
    "validate_graph",
    "load_regions",
    "load_boundaries",
]

REGION_COLUMNS = ["region_id", "ref_name", "start", "end", "length", "signature"]
REGION_COLUMN_DTYPES: dict[str, type | np.dtype] = {
    "region_id": np.int64,
    "start": np.int64,
    "end": np.int64,
    "length": np.int64,
    "signature": np.uint8,
}

EMPTY_I64 = np.zeros(0, np.int64)

EDGE_KIND_CONTIGUOUS = 0
EDGE_KIND_SJ = 1

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
#: the four terminus bits, and the four junction bits, together
FLAG_TERMINUS = FLAG_TSS_POS | FLAG_TSS_NEG | FLAG_TES_POS | FLAG_TES_NEG
FLAG_JUNCTION = FLAG_DONOR_POS | FLAG_DONOR_NEG | FLAG_ACCEPTOR_POS | FLAG_ACCEPTOR_NEG


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
                    continue  # a zero- or negative-length exon covers no base
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
    """Build ``(regions_df, edges_df)`` — the splice graph. Fully vectorised; no per-region Python object.

    Deterministic by construction: ``np.unique`` sorts, region ids are assigned by position,
    boundaries by ``(src, kind, dst)``. No dict iteration order, no hashing and no parallel
    reduction, so the output is byte-identical across runs and platforms.
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
    region_parts, edge_parts = [], []
    region_base = 0
    for name, length in reflen.items():
        if length == 0:
            continue
        rows = by_ref.get(name, EMPTY_I64)
        i_rows = intr_by_ref.get(name, EMPTY_I64)
        regions, edges = _reference_graph(name, length, ex, rows, intr, i_rows, region_base)
        region_parts.append(regions)
        edge_parts.append(edges)
        region_base += len(regions)
    return _assemble(region_parts, edge_parts)


def _reference_graph(name, length, ex: _Exons, rows, intr, i_rows, region_base: int):
    """One reference's ``(regions, edges)`` frames, region ids offset by ``region_base``.

    Every exon endpoint and both reference ends are region bounds. A contiguous boundary sits at each interior
    interface and carries its structural flags; there is one sj boundary per distinct ``(intron_start, intron_end,
    strand)``. A contiguous row's reach columns hold 0: the one reader of a reach takes the sj rows.
    """
    region_bounds = np.unique(
        np.concatenate([ex.start[rows], ex.end[rows], np.array([0, length], np.int64)])
    )
    starts, ends = region_bounds[:-1], region_bounds[1:]
    n = starts.shape[0]
    regions = pd.DataFrame(
        {
            "region_id": np.arange(region_base, region_base + n, dtype=np.int64),
            "ref_name": name,
            "start": starts,
            "end": ends,
            "length": ends - starts,
            "signature": _signatures(region_bounds, ex, rows, intr, i_rows),
        }
    )
    pos = region_bounds[1:-1]  # the n-1 interior interfaces
    flags = _contiguous_flags(pos, ex, rows, intr, i_rows)
    j_src, j_dst, j_strand, j_reach = _sj_edges(region_bounds, intr, i_rows, name)
    n_contiguous = max(n - 1, 0)
    edges = pd.DataFrame(
        {
            "src": np.concatenate([np.arange(n - 1, dtype=np.int64), j_src]) + region_base,
            "dst": np.concatenate([np.arange(1, n, dtype=np.int64), j_dst]) + region_base,
            "kind": np.concatenate(
                [
                    np.full(n_contiguous, EDGE_KIND_CONTIGUOUS, np.uint8),
                    np.full(j_src.size, EDGE_KIND_SJ, np.uint8),
                ]
            ),
            "strand": np.concatenate([np.zeros(n_contiguous, np.int8), j_strand]),
            "flags": np.concatenate([flags, np.zeros(j_src.size, np.uint16)]),
            **{
                k: np.concatenate([np.zeros(pos.shape[0], np.int64), v]).astype(np.int32)
                for k, v in j_reach.items()
            },
        }
    )
    return regions, edges


def _signatures(region_bounds, ex: _Exons, rows, intr, i_rows) -> np.ndarray:
    """Each region's 4-bit signature, by difference arrays over the region index."""
    n = region_bounds.shape[0] - 1
    sig = np.zeros(n, np.uint8)
    for bit, iv_s, iv_e in _signature_intervals(ex, rows, intr, i_rows):
        if iv_s.size == 0:
            continue
        d = np.zeros(n + 1, np.int32)
        np.add.at(d, np.searchsorted(region_bounds, iv_s), 1)
        np.add.at(d, np.searchsorted(region_bounds, iv_e), -1)
        sig |= np.where(np.cumsum(d)[:-1] > 0, bit, 0).astype(np.uint8)
    return sig


def _contiguous_flags(pos, ex: _Exons, rows, intr, i_rows) -> np.ndarray:
    """The structural flags at each interior interface ``pos``: a bit is set where one of its events falls."""
    flags = np.zeros(pos.shape[0], np.uint16)
    if not pos.size:
        return flags
    for arr, bit in _terminus_and_splice_events(ex, rows, intr, i_rows):
        if arr.size:
            hit = np.searchsorted(pos, arr)
            ok = (hit < pos.size) & (pos[np.clip(hit, 0, max(pos.size - 1, 0))] == arr)
            np.bitwise_or.at(flags, hit[ok], bit)
    return flags


def _assemble(region_parts, edge_parts) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The references' frames as one graph: boundaries in ``(src, kind, dst, strand)`` order, ids by row."""
    regions_df = (
        pd.concat(region_parts, ignore_index=True)
        if region_parts
        else pd.DataFrame(columns=REGION_COLUMNS)
    )
    edges_df = (
        pd.concat(edge_parts, ignore_index=True)
        if edge_parts
        else pd.DataFrame(columns=BOUNDARY_COLUMNS[1:])
    )
    if len(edges_df):
        # Strand completes the order: two strand-coincident sj share (src, kind, dst). It only refines the order, so a
        # region's out-boundaries stay contiguous for the CSR.
        edges_df = edges_df.sort_values(
            ["src", "kind", "dst", "strand"], kind="stable"
        ).reset_index(drop=True)
    edges_df.insert(0, "edge_id", np.arange(len(edges_df), dtype=np.int64))
    return _coerce(regions_df, REGION_COLUMNS, REGION_COLUMN_DTYPES), _coerce(
        edges_df, BOUNDARY_COLUMNS, BOUNDARY_COLUMN_DTYPES
    )


def _signature_intervals(ev: _Exons, ei, ev_i, ii):
    """The four ``(bit, starts, ends)`` interval sets that define a region's 4-bit signature.

    Shared by the builder (difference arrays over the region index) and validator I3b (midpoint containment). The two
    evaluate them by different algorithms, so they agree only if both are right.
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
    """The transcript filter, shared by the builder and the validator so the two cannot drift apart.

    Manufactured spans are ``is_synthetic``; a non-synthetic single-exon row is a real transcript with real termini.
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
    boundaries rather than a Python loop over every row.
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
    """``(positions, flag_bit)`` for each structural bit that has events: the builder's derivation.

    A class with no event sets no bit, so it is omitted; :func:`_events_independently` is the validator's counterpart.
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
    """Validator I13's own derivation of the 8 event sets, by a different algorithm from the builder's.

    The builder finds a transcript's outer exons by cumulative exonic offsets; this takes each transcript's
    ``min(start)`` and ``max(end)``. All eight entries are returned, an empty class included, so a spurious flag on a
    reference with no transcript of that strand is still compared. Returns ``(positions, bit)`` pairs.
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


def _sj_edges(region_bounds, ma_i, mi_, ref_name):
    """Distinct sj boundaries for one reference: ``(src, dst, strand, reach dict)``, region-local ids."""
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
    src = np.searchsorted(region_bounds, ua) - 1  # the region ending at intron_start
    dst = np.searchsorted(region_bounds, ub)  # the region starting at intron_end
    if not (np.all(region_bounds[src + 1] == ua) and np.all(region_bounds[dst] == ub)):
        raise ValueError(f"ref {ref_name!r}: a sj endpoint is not a region interface (I5)")
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
# Validation — invariants I1..I13. I3b, I4, I11 and I13 need the transcripts and run only at build; a load-time call
# gets the graph-internal checks (I1/I2/I5-I9/I12).
# ---------------------------------------------------------------------------


def validate_graph(regions_df, edges_df, ref_lengths: Mapping[str, int], transcripts=None) -> None:
    """Assert invariants I1–I13; raises :class:`ValueError` on the first violation.

    ``transcripts`` enables the four annotation checks: I3b (the signature, recomputed), I4 (the interfaces are the
    events), I11 (every transcript walks the graph) and I13 (the flags are the events). Without it only the
    graph-internal checks run, which cannot see a ``signature`` or ``flags`` column that has drifted from the
    annotation. The build passes transcripts; the load path does not, since rebuilding them is slow at human scale.
    """
    reflen = {str(k): int(v) for k, v in ref_lengths.items()}
    ref = regions_df["ref_name"].astype(str).to_numpy()
    # The reference slices come once from the run boundaries; a per-reference string comparison over every row would
    # dominate the validation.
    slices = _ref_slices(ref)
    _check_regions(regions_df, reflen, slices)
    _check_boundaries(regions_df, edges_df, ref, slices)
    if transcripts is not None:
        _validate_against_transcripts(transcripts, reflen, regions_df, edges_df)


def _check_regions(regions_df, reflen: dict[str, int], slices: dict) -> None:
    """I1 (the regions tile each reference exactly), I2 (ids are row order) and I3a (the signature's range)."""
    nid = regions_df["region_id"].to_numpy(np.int64)
    start = regions_df["start"].to_numpy(np.int64)
    end = regions_df["end"].to_numpy(np.int64)
    length = regions_df["length"].to_numpy(np.int64)
    sig = regions_df["signature"].to_numpy(np.uint8)
    if not np.array_equal(nid, np.arange(nid.size, dtype=np.int64)):
        raise ValueError("I2: region_id is not 0..n-1 in row order. Rebuild the index.")
    if nid.size and int(sig.max()) >= N_SIGNATURES:
        raise ValueError(f"I3: signature {int(sig.max())} out of range. Rebuild the index.")
    if not np.array_equal(length, end - start):
        raise ValueError("I1: stored region length != end - start. Rebuild the index.")
    if nid.size and int(length.min()) <= 0:
        raise ValueError("I1: a region has non-positive length. Rebuild the index.")

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


def _check_boundaries(regions_df, edges_df, ref, slices: dict) -> None:
    """I5 (a boundary stays on its reference), I7 (the contiguous boundaries are the adjacent pairs), I8 (``src <
    dst``) and I12 (row order is the sorted id order), then the sj boundaries (:func:`_check_sj_boundaries`)."""
    n_regions = len(regions_df)
    src = edges_df["src"].to_numpy(np.int64)
    dst = edges_df["dst"].to_numpy(np.int64)
    kind = edges_df["kind"].to_numpy(np.uint8)
    estrand = edges_df["strand"].to_numpy(np.int8)

    # I8 — src < dst on every boundary  ⇒ genomic order is a topological order
    if src.size and not np.all(src < dst):
        raise ValueError("I8: a boundary has src >= dst; genomic order is not a topological order.")
    if src.size and not np.all(ref[src] == ref[dst]):
        raise ValueError("I5: a boundary crosses references.")
    if not np.array_equal(edges_df["edge_id"].to_numpy(np.int64), np.arange(src.size)):
        raise ValueError("I12: edge_id is not the row index.")
    if src.size > 1:
        # (src, kind, dst, strand) packed into one int64 key that must strictly increase. Strand is in it because it is
        # in the order; without it two strand-coincident sj would read as a duplicate.
        key = ((src * 2 + kind.astype(np.int64)) * max(n_regions, 1) + dst) * 4 + estrand.astype(
            np.int64
        )
        if not np.all(np.diff(key) > 0):
            raise ValueError(
                "I12: boundaries are not sorted by (src, kind, dst, strand), or a duplicate exists."
            )

    # I7 — contiguous boundaries are exactly {(i, i+1)} within each reference
    c = kind == EDGE_KIND_CONTIGUOUS
    if not np.all(dst[c] == src[c] + 1):
        raise ValueError("I7: a CONTIGUOUS boundary is not between adjacent regions.")
    n_expected = n_regions - len(slices)
    if int(c.sum()) != n_expected:
        raise ValueError(f"I7: {int(c.sum())} contiguous boundaries, expected {n_expected}.")
    if np.any(estrand[c] != 0):
        raise ValueError("I7: a CONTIGUOUS boundary carries a strand.")

    j = kind == EDGE_KIND_SJ
    if j.any():
        sig = regions_df["signature"].to_numpy(np.uint8)
        _check_sj_boundaries(src[j], dst[j], estrand[j], sig, n_regions)


def _check_sj_boundaries(src, dst, strand, sig, n_regions: int) -> None:
    """I5/I6 (the sj land on interfaces, one row per distinct ``(donor, acceptor, strand)``) and I9 (both endpoints
    carry that strand's exon bit)."""
    if not np.all(np.isin(strand, [Strand.POS, Strand.NEG])):
        raise ValueError("I6: a SJ boundary has an invalid strand.")
    # A strand-coincident sj pair is biologically impossible (splice motifs are not palindromic: GT..AG
    # reverse-complements to CT..AC), so one means the annotation is wrong. It warns rather than raises because the
    # graph handles it correctly: strand is in the sort key, so the two boundaries stay distinct and ordered.
    _u, _cnt = np.unique(src * max(n_regions, 1) + dst, return_counts=True)
    n_coincident = int((_cnt > 1).sum())
    if n_coincident:
        warnings.warn(
            f"{n_coincident} strand-coincident splice junction(s): the same (donor, acceptor) "
            f"annotated on BOTH strands. This is biologically impossible — splice motifs are "
            f"non-palindromic (a GT..AG intron reverse-complements to CT..AC) — so the source "
            f"annotation is likely wrong. The graph keeps them as distinct boundaries and is correct "
            f"either way.",
            RuntimeWarning,
            stacklevel=4,
        )
    key = np.stack([src, dst, strand.astype(np.int64)], axis=1)
    if np.unique(key, axis=0).shape[0] != key.shape[0]:
        raise ValueError("I6: a (donor, acceptor, strand) sj appears more than once.")
    for s, bit in ((Strand.POS, BIT_EXON_POS), (Strand.NEG, BIT_EXON_NEG)):
        m = strand == s
        if m.any() and not (np.all(sig[src[m]] & bit) and np.all(sig[dst[m]] & bit)):
            raise ValueError(f"I9: a strand-{int(s)} sj endpoint lacks its exon bit.")


@dataclass(frozen=True, slots=True)
class _Reference:
    """One reference's slice of the graph and of the annotation, for the annotation checks: its name and length, its
    region rows ``m`` (the first is ``base``) with their bounds, and its exon and intron rows."""

    name: str
    length: int
    m: np.ndarray
    base: int
    region_bounds: np.ndarray
    rows: np.ndarray
    i_rows: np.ndarray


def _validate_against_transcripts(transcripts, reflen, regions_df, edges_df) -> None:
    """I3b (the signature, recomputed), I4 (the interfaces are the events), I11 (every transcript walks the graph) and
    I13 (the flags are the events), vectorised so the check stays cheap enough to run on every build."""
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
    jm = kind == EDGE_KIND_SJ
    jkey = np.sort((src[jm] * n_regions + dst[jm]) * 4 + estrand[jm].astype(np.int64))

    slices = _ref_slices(ref)
    for name, L in reflen.items():
        m = slices.get(name)
        if m is None:
            continue
        r = _Reference(
            name=name,
            length=L,
            m=m,
            base=int(m[0]),
            region_bounds=np.concatenate([start[m], end[m][-1:]]),
            rows=by_ref.get(name, EMPTY_I64),
            i_rows=intr_by_ref.get(name, EMPTY_I64),
        )
        _check_interfaces_are_events(r, ex, end)
        _check_signatures_recomputed(r, ex, intr, start, end, sig)
        _check_flags_are_events(r, ex, intr, end, csrc, cflags)
        _check_exon_ends_are_interfaces(r, ex)
        _check_introns_are_sj(r, intr, jkey, n_regions)


def _check_interfaces_are_events(r: _Reference, ex: _Exons, end) -> None:
    """I4 — the interior interfaces are exactly the annotation events."""
    events = (
        np.unique(np.concatenate([ex.start[r.rows], ex.end[r.rows]])) if r.rows.size else EMPTY_I64
    )
    events = events[(events > 0) & (events < r.length)]
    interior = end[r.m[:-1]] if r.m.size > 1 else EMPTY_I64
    if not np.array_equal(np.unique(interior), events):
        extra = np.setdiff1d(interior, events).size
        missing = np.setdiff1d(events, interior).size
        raise ValueError(
            f"I4: reference {r.name!r} interior interfaces are not exactly the annotation events "
            f"({extra} extra, {missing} missing)."
        )


def _check_signatures_recomputed(r: _Reference, ex: _Exons, intr, start, end, sig) -> None:
    """I3b — recompute the signature from each region's midpoint by direct interval containment, a different algorithm
    from the builder's difference sweep. Every interval endpoint is a region bound, so a region is homogeneous and its
    midpoint decides it."""
    m = r.m
    mid = (start[m] + end[m]) // 2
    want = np.zeros(mid.size, np.uint8)
    for bit, iv_s, iv_e in _signature_intervals(ex, r.rows, intr, r.i_rows):
        if iv_s.size:
            # covered iff #(start <= mid) > #(end <= mid), for half-open [start, end)
            n_started = np.searchsorted(np.sort(iv_s), mid, side="right")
            n_ended = np.searchsorted(np.sort(iv_e), mid, side="right")
            want |= np.where(n_started > n_ended, bit, 0).astype(np.uint8)
    if not np.array_equal(want, sig[m]):
        k = int(np.flatnonzero(want != sig[m])[0])
        raise ValueError(
            f"I3: reference {r.name!r} region {int(m[k])} [{int(start[m][k])}, "
            f"{int(end[m][k])}) has signature {int(sig[m][k])}, recomputed {int(want[k])}."
        )


def _check_flags_are_events(r: _Reference, ex: _Exons, intr, end, csrc, cflags) -> None:
    """I13 — each structural bit is set at exactly the positions that generate it, in both directions, for all eight
    bits, including a class with no event on this reference. The converse of I4 does not hold bit-wise: the interface
    between bookended exons of one transcript (a zero-length intron) carries no flag."""
    m = r.m
    if m.size <= 1:
        return
    interior = end[m[:-1]]
    at = np.zeros(interior.size, np.uint16)
    e_of_region = np.searchsorted(csrc, m[:-1])
    hit = (e_of_region < csrc.size) & (
        csrc[np.clip(e_of_region, 0, max(csrc.size - 1, 0))] == m[:-1]
    )
    at[hit] = cflags[e_of_region[hit]]
    for arr, bit in _events_independently(ex, r.rows, intr, r.i_rows):
        want = np.isin(interior, arr)
        got = (at & bit) != 0
        if not np.array_equal(want, got):
            k = int(np.flatnonzero(want != got)[0])
            raise ValueError(
                f"I13: reference {r.name!r} position {int(interior[k])} "
                f"{'is missing' if want[k] else 'wrongly carries'} flag bit {int(bit):#06x}."
            )


def _check_exon_ends_are_interfaces(r: _Reference, ex: _Exons) -> None:
    """I11a — every annotated exon boundary lands exactly on a region interface."""
    if not r.rows.size:
        return
    bounds = r.region_bounds
    for arr, what in ((ex.start[r.rows], "start"), (ex.end[r.rows], "end")):
        hit = np.searchsorted(bounds, arr)
        bad = (hit >= bounds.size) | (bounds[np.clip(hit, 0, bounds.size - 1)] != arr)
        if bad.any():
            raise ValueError(
                f"I11: reference {r.name!r} has {int(bad.sum())} exon {what}s that are not region "
                f"interfaces (first at {int(arr[bad][0])})."
            )


def _check_introns_are_sj(r: _Reference, intr, jkey, n_regions: int) -> None:
    """I11b — every annotated intron is realised by a SJ boundary on that strand."""
    if not r.i_rows.size:
        return
    a, b, st = intr[1][r.i_rows], intr[2][r.i_rows], intr[3][r.i_rows].astype(np.int64)
    js = r.base + np.searchsorted(r.region_bounds, a) - 1
    jd = r.base + np.searchsorted(r.region_bounds, b)
    want = (js * n_regions + jd) * 4 + st
    if jkey.size == 0:
        miss = np.ones(want.shape, bool)
    else:
        pos = np.searchsorted(jkey, want)
        miss = (pos >= jkey.size) | (jkey[np.clip(pos, 0, jkey.size - 1)] != want)
    if miss.any():
        k = int(np.flatnonzero(miss)[0])
        raise ValueError(
            f"I11: reference {r.name!r} has {int(miss.sum())} annotated introns with no SJ "
            f"boundary (first [{int(a[k])}, {int(b[k])}) strand {int(st[k])})."
        )


# ---------------------------------------------------------------------------
# The scanner adapter + load
# ---------------------------------------------------------------------------


def build_region_partition_arrays(index) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flatten the graph into the ``BamScanner.set_regions`` ABI: the partition the accumulator deposits into.

    For each reference in ``index.ref_names`` the bounds are ``[n0.start, n0.end, n1.end, …]`` over its regions, so
    ``k`` regions give ``k + 1`` positions; a reference with no region contributes none.
    :meth:`RegionArrays.from_index <rigel.calibration.region_arrays.RegionArrays.from_index>` reads the same frame,
    so calibration and the scanner address one partition.

    Returns ``(region_bounds int64[P], ref_pos_offsets int64[n_refs+1], region_types uint8[N])``, ``region_types``
    the coarse type (0=intergenic, 1=intron, 2=exon) aligned 1:1 with the regions.
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
class SpliceJunctionArrays:
    """The sj boundaries re-indexed onto the accumulator's flat region-bound axis, as a CSR.

    It is the deposit rule's per-intron lookup: is this observed intron an annotated sj, and which one? Every annotated
    intron's endpoints are region bounds (a bound sits at every exon endpoint), so the lookup is the binary search the
    deposit already does, then a scan of one CSR row, keyed by the flat index into
    ``build_region_partition_arrays(index)[0]`` of the intron's low endpoint::

        for k in range(offsets[boundary_left], offsets[boundary_left + 1]):
            if boundary_right[k] == observed_right_boundary:
                credit sj boundary k          # the CSR slot k is the sj-boundary id

    ``edge_row`` is only the join key back to ``index.edges_df`` and must never index an sj bank: edge rows number
    every boundary, so they run past a bank sized to the sj count and permute it even in bounds. The accumulator never
    receives it. Because the id is the rank, the slot order must match the reference accumulator's
    ``Partition.from_region_bounds`` exactly; ``strand`` is in the sort key and separates two sj that share a
    coordinate pair.

    Gate: ``tests/calibration/test_splice_graph.py::test_the_csr_slot_order_matches_the_reference_accumulator``.
    """

    offsets: (
        np.ndarray
    )  # int64[P + 1] — CSR over the flat region_bound axis, P = region_bounds.size
    boundary_right: np.ndarray  # int64[J] — flat region_bound index of the intron's high endpoint
    edge_row: np.ndarray  # int64[J] — row in index.edges_df: a join key, not the sj-boundary id
    strand: np.ndarray  # int8[J]  — the sj's genomic strand (Strand POS/NEG)


def build_sj_arrays(index) -> SpliceJunctionArrays:
    """Build the :class:`SpliceJunctionArrays` CSR for ``index``.

    An sj boundary stores region ids; the accumulator works in region-bound indices. For a reference whose first region
    is ``region_base`` and whose first bound is ``region_bound_base``::

        left bound  = region_bound_base + (src - region_base) + 1      # the intron starts where src ends
        right bound = region_bound_base + (dst - region_base)          # and ends where dst begins
    """
    regions_df, edges_df = index.regions_df, index.edges_df
    _, region_bound_offsets, _ = build_region_partition_arrays(index)
    n_region_bounds = int(region_bound_offsets[-1])

    # per-reference region_base, in the same reference order the region_bound axis uses
    counts = regions_df.groupby("ref_name", sort=False).size()
    region_base = np.zeros(len(index.ref_names) + 1, dtype=np.int64)
    for i, ref in enumerate(index.ref_names):
        region_base[i + 1] = region_base[i] + int(counts.get(ref, 0))

    is_sj = edges_df["kind"].to_numpy(np.uint8) == EDGE_KIND_SJ
    src = edges_df["src"].to_numpy(np.int64)[is_sj]
    dst = edges_df["dst"].to_numpy(np.int64)[is_sj]
    strand = edges_df["strand"].to_numpy(np.int8)[is_sj]
    edge_row = np.flatnonzero(is_sj).astype(np.int64)

    # which reference each sj belongs to: region ids are contiguous per reference (I2)
    ref_of = np.searchsorted(region_base, src, side="right") - 1
    shift = region_bound_offsets[ref_of] - region_base[ref_of]
    boundary_left = src + shift + 1
    boundary_right = dst + shift

    # The sort key must match the reference accumulator's ``Partition.from_region_bounds`` exactly, since the
    # sj-boundary id is the rank in this order. Strand matters only for a coordinate pair annotated on both strands.
    order = np.lexsort((strand, boundary_right, boundary_left))
    boundary_left, boundary_right = boundary_left[order], boundary_right[order]
    offsets = np.zeros(n_region_bounds + 1, dtype=np.int64)
    np.cumsum(np.bincount(boundary_left, minlength=n_region_bounds), out=offsets[1:])
    return SpliceJunctionArrays(offsets, boundary_right, edge_row[order], strand[order])


def build_boundary_flags_array(index) -> np.ndarray:
    """The structural flags on the accumulator's contiguous-boundary axis, one entry per boundary.

    A contiguous boundary is the interface to the right of its ``src`` region, so the flags are keyed by ``src``; a
    reference with ``k`` regions contributes ``k - 1`` entries, with no padding. Sj boundaries are not a genomic
    position and carry no flags. Returns ``uint16[E]``, ``E == ref_boundary_offsets[-1]``, aligned with the payload's
    contiguous-boundary axis.
    """
    regions_df, edges_df = index.regions_df, index.edges_df
    contiguous = edges_df["kind"].to_numpy(np.uint8) == EDGE_KIND_CONTIGUOUS
    by_region = np.zeros(len(regions_df), dtype=np.uint16)
    by_region[edges_df["src"].to_numpy(np.int64)[contiguous]] = edges_df["flags"].to_numpy(
        np.uint16
    )[contiguous]

    by_ref: dict[str, pd.DataFrame] = {
        ref: grp for ref, grp in regions_df.groupby("ref_name", sort=False)
    }
    out: list[np.ndarray] = []
    for ref in index.ref_names:
        grp = by_ref.get(ref)
        if grp is None or len(grp) == 0:
            continue
        ids = grp.index.to_numpy(np.int64)  # == region_id (I2), and contiguous within a reference
        out.append(
            by_region[ids[:-1]]
        )  # the k-1 interior boundaries; a 1-region reference contributes none
    return np.concatenate(out) if out else np.zeros(0, dtype=np.uint16)


@dataclass(frozen=True, slots=True)
class SpliceJunctionGeometry:
    """The sj boundaries on the accumulator's sj axis, in its slot order.

    An sj is a factor on its two endpoint regions, never a chain slot or message channel, because the graph is a DAG
    but not a polytree; this holds what a consumer needs to place that factor. ``strand`` is the transcript strand,
    from which sense is derived: the accumulator's ``sj_count`` columns are the genome strand the read aligned to.

    A reach is how many bases of the molecule's own sequence remain either side of the sj. Only a molecule that spliced
    across it uses an sj, so the reaches are exonic. They are ``lo`` / ``hi`` (genomically lower / higher), never
    donor / acceptor, because ``src`` is genomically left whatever the strand. The graph stores them per strand (an sj
    row fills only its own strand's pair, so two strand-coincident sj never read each other's) and maximal over
    isoforms independently per side, which over-states opportunity and so errs toward a higher gDNA fraction. A reach
    of 0 means no RNA of that strand splices across the sj: zero opportunity, not a sentinel.
    """

    src_region: np.ndarray  # int64[J] — global region id; the intron starts where this region ends
    dst_region: np.ndarray  # int64[J] — global region id; the intron ends where this region begins
    strand: np.ndarray  # int8[J]  — the sj's own genomic strand, which is the transcript strand
    reach_lo: np.ndarray  # float64[J] — exonic reach on the sj's own strand
    reach_hi: np.ndarray  # float64[J]

    @property
    def n_sj(self) -> int:
        return int(self.src_region.shape[0])


def build_sj_geometry_arrays(index) -> SpliceJunctionGeometry:
    """Build :class:`SpliceJunctionGeometry` for ``index``, in the accumulator's sj slot order.

    The slot order is read off :func:`build_sj_arrays` (its ``edge_row``), never recomputed, so the order that must
    match the reference accumulator has one implementation.
    """
    rows = build_sj_arrays(index).edge_row
    boundaries = index.edges_df
    strand = boundaries["strand"].to_numpy(np.int8)[rows]
    is_pos = strand == np.int8(Strand.POS)

    def reach(column_pos: str, column_neg: str) -> np.ndarray:
        # An sj row fills only its own strand's pair, so select on the strand: a NEG sj reading the POS column would
        # read zeros and look opportunity-free.
        return np.where(
            is_pos,
            boundaries[column_pos].to_numpy(np.int64)[rows],
            boundaries[column_neg].to_numpy(np.int64)[rows],
        ).astype(np.float64)

    return SpliceJunctionGeometry(
        src_region=boundaries["src"].to_numpy(np.int64)[rows],
        dst_region=boundaries["dst"].to_numpy(np.int64)[rows],
        strand=strand,
        reach_lo=reach("reach_lo_pos", "reach_lo_neg"),
        reach_hi=reach("reach_hi_pos", "reach_hi_neg"),
    )


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
