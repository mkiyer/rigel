"""Capture-aware EM effective lengths: a transcript's conserved shares of its objects at their capture efficiencies.

Layer 2 (opportunity). Under hybrid capture only the probe-enriched part of a transcript's footprint is sampled, and
contracting gDNA alone would concentrate it against RNA, so every transcript's EM length is contracted by
calibration's own per-object gDNA capture efficiencies, on the objects the deposit rule uses. A transcript's pieces
are the regions its exons cover (an unspliced entity's single exon is its span); between consecutive pieces is a cut,
a boundary or a junction. With ``S_p`` its contained share of piece ``p`` and ``M_k`` its conserved share at cut ``k``
(``effective_length.conserved_cut_shares``), which together total the fl-marginal length exactly,

    eff_t = fl_t · (Σ_pieces S_p c_p + Σ_cuts M_k c_k) / (Σ_pieces S_p + Σ_cuts M_k).

A piece prices at its region's efficiency, a contiguous cut at its boundary's, a junction by conservation of bases
(:func:`_junction_prices`), up to a cap RNA's own spliced fragments set (:func:`fit_junction_cap`). A placement that
leaves a cut's exons across a junction is priced there at that junction's crossing level where it is higher
(:func:`_best_part`): gDNA has no such placement to read it by. With no landscape,
or every efficiency exactly 1, ``fl_eff_lengths`` comes back bit-identically, whatever the length tables. Efficiencies are yields against the bulk, so a captured transcript's length exceeds ``fl``, and a
factor common to every component cancels in the EM. Pieces are cut at the index's unobservable mask; a template the
mask does not touch is unchanged to the bit.

Gate: ``tests/calibration/test_capture_eff_length.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..config import CONSTANTS
from ..types import Strand
from .effective_length import conserved_cut_shares, contained_eff_length
from .region_arrays import overlapping_region_runs, region_right_boundary
from .signature import RegionType, coarse_type_array

if TYPE_CHECKING:
    from ..index import TranscriptIndex
    from .region_arrays import RegionArrays
    from .result import CalibrationResult
    from .shared_length_law import LengthTables

__all__ = [
    "TranscriptFrame",
    "TranscriptObjects",
    "fit_junction_cap",
    "observable_fl_eff_lengths",
    "place_efficiencies",
    "place_log_weights",
    "rna_captured_fraction",
    "transcript_capture_eff_lengths",
]


@dataclass(frozen=True, slots=True)
class TranscriptObjects:
    """Every transcript's pieces and cuts with its conserved share of each, in its own coordinates.

    Pieces, one row each in transcript then genomic order: ``t`` the transcript, ``piece`` the region,
    ``contained`` its contained share. Cuts, one row per pair of consecutive pieces of a transcript:
    ``cut_row`` the left piece's row (the right one's is ``cut_row + 1``), ``cut_share`` the conserved
    share, ``is_junction`` whether the transcript splices there. A cut's HOME is the exon it lies in, or the two exons
    a junction joins: ``leave_lo`` and ``leave_hi`` are the cuts (rows of the cut arrays, ``-1`` for none) that are the
    junctions ending the home below and above, and ``leave_lo_share``, ``leave_hi_share`` and ``leave_both_share`` the
    parts of ``cut_share`` from placements that leave the home across the one, the other, or both; the rest of the
    share is from placements inside the home. Annotation and pmf only.
    """

    t: np.ndarray
    piece: np.ndarray
    contained: np.ndarray
    cut_row: np.ndarray
    cut_share: np.ndarray
    is_junction: np.ndarray
    leave_lo: np.ndarray
    leave_hi: np.ndarray
    leave_lo_share: np.ndarray
    leave_hi_share: np.ndarray
    leave_both_share: np.ndarray


@dataclass(frozen=True, slots=True)
class TranscriptPieces:
    """The pmf-free half of :class:`TranscriptObjects`: every transcript's pieces, each piece's length, and at each
    cut the bases of the transcript left of it (``through``) and right of it (``rest``), counted to the nearest
    masked base. ``run_t`` and ``run_length`` are each transcript's runs — its stretches no masked base interrupts —
    and ``masked`` the bases the mask took from each transcript. ``home_lo`` and ``home_hi`` are the bases of each
    cut's home either side of it, and ``leave_lo`` and ``leave_hi`` the cuts that are the junctions ending the home
    below and above (:func:`_homes`)."""

    t: np.ndarray
    piece: np.ndarray
    length: np.ndarray
    cut_row: np.ndarray
    through: np.ndarray
    rest: np.ndarray
    is_junction: np.ndarray
    run_t: np.ndarray
    run_length: np.ndarray
    masked: np.ndarray
    home_lo: np.ndarray
    home_hi: np.ndarray
    leave_lo: np.ndarray
    leave_hi: np.ndarray


class TranscriptFrame:
    """Every template's pieces and cuts (:func:`transcript_pieces`) and their shares on the EM's length tables
    (:func:`object_shares`, and :func:`_carry_shares` where capture's excess is carried), each built on first use and
    kept. The junction cap's fit, the observable lengths and the capture-contracted lengths read the same frame; on a
    whole genome it is ten million rows and each table's shares a second's work, so it is built once and dropped once
    the lengths are."""

    def __init__(
        self,
        index: "TranscriptIndex",
        region_arrays: "RegionArrays",
        tables: "LengthTables",
        n_threads: int = 1,
    ) -> None:
        self.index = index
        self.region_arrays = region_arrays
        self.tables = tables
        self.n_threads = n_threads

    @cached_property
    def pieces(self) -> TranscriptPieces:
        return transcript_pieces(self.index, self.region_arrays)

    @cached_property
    def rna(self) -> TranscriptObjects:
        """The shares on RNA's table."""
        return object_shares(self.pieces, self.tables.rna, self.n_threads)

    @cached_property
    def carry(self) -> tuple[TranscriptObjects, TranscriptObjects] | None:
        """``(on_anchor, on_table)``, or ``None`` where nothing is carried (:func:`_carry_shares`)."""
        return _carry_shares(self.pieces, self.rna, self.tables, self.n_threads)


def _observable_parts(
    t: np.ndarray,
    piece: np.ndarray,
    k: np.ndarray,
    ref: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    u: np.ndarray,
    mask: np.ndarray,
):
    """Each piece ``[a, b)`` on reference ``ref`` cut at the mask's runs into its observable parts: a run inside it
    splits it, one over an end trims it, one over the whole of it drops it. Each part keeps its transcript, region
    and exon and its start ``u`` in the transcript's full coordinates, masked bases included. Returns ``(t, piece,
    k, length, u)`` for the parts, in the pieces' order."""
    length = (b - a).astype(np.float64)
    mask = np.asarray(mask, dtype=np.int64)
    if mask.shape[0] == 0 or t.size == 0:
        return t, piece, k, length, u
    # the first masked run ending past each piece's start, on one ascending key: reference, then position (every
    # position is below 2^40)
    span = np.int64(1) << np.int64(40)
    first = np.searchsorted(mask[:, 0] * span + mask[:, 2], ref * span + a, side="right")
    at = np.minimum(first, mask.shape[0] - 1)
    hit = (first < mask.shape[0]) & (mask[at, 0] == ref) & (mask[at, 1] < b)
    if not hit.any():
        return t, piece, k, length, u
    rows, lo, hi = [], [], []
    for i in np.flatnonzero(hit):
        pos, m = int(a[i]), int(first[i])
        while m < mask.shape[0] and mask[m, 0] == ref[i] and mask[m, 1] < b[i]:
            if mask[m, 1] > pos:
                rows.append(i)
                lo.append(pos)
                hi.append(int(mask[m, 1]))
            pos = max(pos, int(mask[m, 2]))
            m += 1
        if pos < b[i]:
            rows.append(i)
            lo.append(pos)
            hi.append(int(b[i]))
    # every piece keeps its place: one part where no run touches it, its parts in order where one does (the loop
    # writes them ascending), none where a run covers it
    rows = np.asarray(rows, dtype=np.int64)
    n_parts = np.where(hit, 0, 1)
    n_parts[hit] = np.bincount(rows, minlength=t.size)[hit]
    row = np.repeat(np.arange(t.size), n_parts)
    slot = np.cumsum(n_parts) - n_parts
    pa, pb = np.empty(row.size, dtype=np.int64), np.empty(row.size, dtype=np.int64)
    keep = np.flatnonzero(~hit)
    pa[slot[keep]], pb[slot[keep]] = a[keep], b[keep]
    split = slot[rows] + (np.arange(rows.size) - np.searchsorted(rows, rows, side="left"))
    pa[split], pb[split] = np.asarray(lo, dtype=np.int64), np.asarray(hi, dtype=np.int64)
    return t[row], piece[row], k[row], (pb - pa).astype(np.float64), u[row] + (pa - a[row])


def transcript_pieces(index: "TranscriptIndex", region_arrays: "RegionArrays") -> TranscriptPieces:
    """The observable parts of the regions every transcript's exons (or an unspliced entity's span) cover, in
    transcript then genomic order, and between consecutive parts a cut — a junction where they belong to different
    exons with a region between them. No part holds a base of the index's unobservable mask
    (``region_arrays.unobservable``) and no cut joins two parts a masked base separates: a fragment touching one is
    never observed. With nothing masked the parts are the pieces. Annotation only: no length law enters."""
    t, rid, a, b = _template_spans(index)
    tk, piece, k, seg_a, seg_b = _region_overlaps(t, rid, a, b, region_arrays)
    length = (seg_b - seg_a).astype(np.float64)
    n_t = int(index.t_df.shape[0])
    whole = _total_by(tk, length, n_t)
    tk, piece, k, length, full_start = _observable_parts(
        tk, piece, k, rid[k], seg_a, seg_b, _start_within(tk, length), region_arrays.unobservable
    )
    return _runs(tk, piece, k, length, full_start, n_t, whole - _total_by(tk, length, n_t))


def _template_spans(index: "TranscriptIndex"):
    """``(t, ref_id, start, end)``: every transcript's exons, a nascent span's one exon its whole span, in transcript
    then genomic order: the index's exon CSR (:meth:`TranscriptIndex.build_exon_csr`). ``ref_id`` is ``-1`` on a
    reference the index does not carry."""
    offsets, starts, ends, _ = index.build_exon_csr()
    t = np.repeat(np.arange(offsets.size - 1, dtype=np.int64), np.diff(offsets))
    ref = index.t_df["ref"].astype(str).to_numpy()[t]
    rid = pd.Series(ref).map(index.ref_name_to_id).fillna(-1).to_numpy(np.int64)
    return t, rid, starts.astype(np.int64), ends.astype(np.int64)


def _region_overlaps(t, rid, a, b, region_arrays: "RegionArrays"):
    """Each span expanded into the regions it covers: ``(t, region, span row, start, end)`` per piece, a piece being a
    span's overlap with one region."""
    starts = np.asarray(region_arrays.start, dtype=np.int64)
    ends = np.asarray(region_arrays.end, dtype=np.int64)
    lo, hi = overlapping_region_runs(
        rid, a, b, starts, ends, np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    )
    width = np.where((rid >= 0) & (hi > lo), hi - lo, 0)
    k = np.repeat(np.arange(t.size), width)
    piece = np.repeat(lo, width) + (
        np.arange(int(width.sum()), dtype=np.int64) - np.repeat(np.cumsum(width) - width, width)
    )
    return t[k], piece, k, np.maximum(starts[piece], a[k]), np.minimum(ends[piece], b[k])


def _start_within(tk: np.ndarray, length: np.ndarray) -> np.ndarray:
    """Each piece's start in its transcript's full coordinates (the transcript's bases before it), masked bases
    included."""
    first = np.r_[True, tk[1:] != tk[:-1]] if tk.size else np.zeros(0, dtype=bool)
    csum = np.cumsum(length)
    return (
        csum
        - np.repeat(csum[first] - length[first], np.diff(np.r_[np.flatnonzero(first), tk.size]))
        - length
    )


def _total_by(tk: np.ndarray, length: np.ndarray, n_t: int) -> np.ndarray:
    """``float64[n_t]`` — each transcript's total over its pieces."""
    out = np.zeros(n_t)
    np.add.at(out, tk, length)
    return out


def _runs(tk, piece, k, length, full_start, n_t: int, masked) -> TranscriptPieces:
    """The pieces' runs — consecutive parts of one transcript with no masked base between them — and the cuts, each
    joining two parts of one run, with the transcript's bases either side of it."""
    full_end = full_start + length
    joined = (
        np.r_[(tk[1:] == tk[:-1]) & (full_start[1:] == full_end[:-1]), False]
        if tk.size
        else np.zeros(0, dtype=bool)
    )
    run_first = np.r_[True, ~joined[:-1]] if tk.size else np.zeros(0, dtype=bool)
    run_id = np.cumsum(run_first) - 1
    run_lo = full_start[run_first]
    run_hi = full_end[np.r_[np.flatnonzero(run_first)[1:] - 1, tk.size - 1]] if tk.size else run_lo
    cut = np.flatnonzero(joined)
    is_junction = (k[cut] != k[cut + 1]) & (piece[cut + 1] > piece[cut] + 1)
    home_lo, home_hi, leave_lo, leave_hi = _homes(cut, is_junction, full_start, full_end)
    return TranscriptPieces(
        t=tk,
        piece=piece,
        length=length,
        cut_row=cut,
        through=full_end[cut] - run_lo[run_id[cut]],
        rest=run_hi[run_id[cut]] - full_end[cut],
        is_junction=is_junction,
        run_t=tk[run_first],
        run_length=run_hi - run_lo,
        masked=masked,
        home_lo=home_lo,
        home_hi=home_hi,
        leave_lo=leave_lo,
        leave_hi=leave_hi,
    )


def _homes(cut: np.ndarray, is_junction: np.ndarray, full_start: np.ndarray, full_end: np.ndarray):
    """``(home_lo, home_hi, leave_lo, leave_hi)`` per cut. A cut's HOME is the exon it lies in (a contiguous cut) or the
    two exons it joins (a junction), an exon being the pieces contiguous cuts join. ``home_lo`` and ``home_hi`` are the
    home's bases below and above the cut, so a placement with more bases than that on a side leaves the home there;
    ``leave_lo`` and ``leave_hi`` are the cuts it leaves across, the junctions ending the home at each end, ``-1``
    where the transcript or a masked base ends it and nothing leaves."""
    n = full_start.shape[0]
    none = np.zeros(0, dtype=np.int64)
    if n == 0:
        return np.zeros(0), np.zeros(0), none, none
    same_exon = np.zeros(n, dtype=bool)  # piece i and piece i + 1 lie in one exon
    same_exon[cut[~is_junction]] = True
    first = np.flatnonzero(np.r_[True, ~same_exon[:-1]])
    last = np.r_[first[1:] - 1, n - 1]
    exon = np.repeat(np.arange(first.shape[0]), last - first + 1)
    cut_after = np.full(n, -1, dtype=np.int64)  # the cut joining piece i to piece i + 1
    cut_after[cut] = np.arange(cut.shape[0])
    # an exon's pieces end where no contiguous cut follows, so a cut at either end of it is a junction
    below = np.where(first > 0, cut_after[np.maximum(first - 1, 0)], -1)
    above = cut_after[last]
    lo_exon, hi_exon = exon[cut], exon[cut + 1]
    return (
        full_end[cut] - full_start[first][lo_exon],
        full_end[last][hi_exon] - full_end[cut],
        below[lo_exon],
        above[hi_exon],
    )


def observable_fl_eff_lengths(frame: TranscriptFrame, rna_fl) -> np.ndarray:
    """Every template's fl-marginal length over the starts the index's mask leaves: ``Σ_runs Σ_w f(w) (run − w + 1)₊``
    over its runs, the stretches no masked base interrupts, since a fragment touching a masked base is never observed.
    A template the mask does not touch keeps the closed form of its whole length, to the bit; every template keeps the
    closed form's floor of one base, as one shorter than every fragment does."""
    lengths = frame.index.t_df["length"].to_numpy(np.int64)
    eff = rna_fl.compute_all_transcript_eff_lens(lengths)
    pieces = frame.pieces
    touched = pieces.masked > 0.0
    if not touched.any():
        return eff
    on = touched[pieces.run_t]
    over_runs = np.zeros(lengths.size)
    np.add.at(
        over_runs,
        pieces.run_t[on],
        rna_fl.compute_all_transcript_eff_lens(
            pieces.run_length[on].astype(np.int64), min_value=0.0
        ),
    )
    eff[touched] = np.maximum(over_runs[touched], CONSTANTS.fragment_length.min_effective_length)
    return eff


def object_shares(
    pieces: TranscriptPieces, pmf: np.ndarray, n_threads: int = 1
) -> TranscriptObjects:
    """The pieces' contained shares and the cuts' conserved shares on one pmf. A zero-length fragment places
    nowhere — the accumulator deposits none — so the pmf's mass at ``w = 0`` is dropped, and the shares partition
    the fl-marginal length the pipeline builds, which excludes it too. Every cut's shares are its own, so the cuts are
    priced a block at a time (``CONSTANTS.resources.share_block_rows``) and the temporaries stay a block's size."""
    pmf = np.asarray(pmf, dtype=np.float64).copy()
    pmf[0] = 0.0
    cut = pieces.cut_row
    n = cut.shape[0]
    whole, low, high, both = (np.empty(n, dtype=np.float64) for _ in range(4))
    rows = CONSTANTS.resources.share_block_rows
    for lo in range(0, n, rows):
        hi = min(lo + rows, n)
        block = cut[lo:hi]
        below, above = pieces.length[block], pieces.length[block + 1]

        def share(reach_lo, reach_hi):
            left, right = conserved_cut_shares(
                pmf, below, above, reach_lo[lo:hi], reach_hi[lo:hi], n_threads=n_threads
            )
            return left + right

        # a placement leaves the cut's home where it holds more bases on a side than the home does: the share at the
        # template's reach less the share at the home's is the part that leaves
        whole[lo:hi] = share(pieces.through, pieces.rest)
        inside = share(pieces.home_lo, pieces.home_hi)
        low[lo:hi] = np.maximum(share(pieces.through, pieces.home_hi) - inside, 0.0)
        high[lo:hi] = np.maximum(share(pieces.home_lo, pieces.rest) - inside, 0.0)
        both[lo:hi] = np.maximum(whole[lo:hi] - inside - low[lo:hi] - high[lo:hi], 0.0)
    return TranscriptObjects(
        t=pieces.t,
        piece=pieces.piece,
        contained=contained_eff_length(pieces.length, pmf),
        cut_row=cut,
        cut_share=whole,
        is_junction=pieces.is_junction,
        leave_lo=pieces.leave_lo,
        leave_hi=pieces.leave_hi,
        leave_lo_share=low,
        leave_hi_share=high,
        leave_both_share=both,
    )


def cut_efficiencies(
    objects: TranscriptObjects,
    region_arrays: "RegionArrays",
    c_region: np.ndarray,
    c_boundary: np.ndarray,
    gdna_contained: np.ndarray,
    junction_cap: float,
    law: TranscriptObjects | None = None,
) -> np.ndarray:
    """Every cut's efficiency. Its own is a contiguous cut's boundary's, a junction's from :func:`_junction_prices` at
    ``junction_cap``; the part of its share from placements that leave its home across a junction is priced at that
    junction's crossing level where it is higher (:func:`_best_part`), the split taken on ``law``, the objects on the
    law the efficiencies are read on (``objects`` by default). ``A`` is the piece below the cut and ``B`` the piece
    above it, in genomic order."""
    right_of = region_right_boundary(np.asarray(region_arrays.ref_id))
    A = objects.piece[objects.cut_row]
    B = objects.piece[objects.cut_row + 1]
    # a cut's lower piece has its upper piece above it on its reference, so a boundary to its right
    below = right_of[A]
    if below.size and int(below.min()) < 0:
        raise ValueError(
            "a cut's lower piece is the last region of its reference: the transcript frame is corrupt"
        )
    out = c_boundary[below]
    crossing = np.zeros(out.shape[0])
    j = np.flatnonzero(objects.is_junction)
    if j.size:
        out = out.copy()
        out[j], crossing[j] = _junction_prices(
            A[j], B[j], right_of, c_region, c_boundary, gdna_contained, junction_cap
        )
    law = objects if law is None else law
    return _best_part(
        out,
        _across(crossing, law.leave_lo),
        _across(crossing, law.leave_hi),
        law.leave_lo_share,
        law.leave_hi_share,
        law.leave_both_share,
        law.cut_share,
    )


def _across(value: np.ndarray, leave: np.ndarray) -> np.ndarray:
    """Per cut, ``value`` at the junction its placements leave its home across (``leave``, a row of the cuts); 0 where
    the home ends with no junction. An efficiency is never negative, so 0 there raises nothing."""
    return np.where(leave >= 0, value[np.maximum(leave, 0)], 0.0)


def _best_part(own, x_lo, x_hi, share_lo, share_hi, share_both, share) -> np.ndarray:
    """Every cut's price from its own (``own``) and the crossing levels (``x_lo``, ``x_hi``) of the junctions its
    placements leave its home across, below and above (:func:`_across`).

    gDNA's fragment at a cut near an exon's end runs on into the intron; RNA's runs across the junction into the next
    exon, and gDNA has no such placement to read it by. That placement is RNA's crossing of the junction too, and a
    fragment binds its best single probe part, so its deposit at the cut is priced at the junction's crossing level
    where that is higher::

        c_cut = own + ( lo·(x_lo − own)₊ + hi·(x_hi − own)₊ + both·(max(x_lo, x_hi) − own)₊ ) / share

    ``lo``, ``hi`` and ``both`` are the parts of the cut's share from placements leaving across the junction below,
    the one above, or both; the part inside the home keeps ``own``. A
    junction's crossing level is conservation's sum up to the junction cap, ``min(summed, cap)``: what RNA's crossing
    captures there, as the library's junctions read it. The best object beside a junction is not passed on. It is one
    gDNA object's own reading, the largest of four, and it belongs to the junction it lies beside. A cut whose
    placements all stay home, a template with no junction, and a field of one level keep ``own`` to the bit."""
    excess = (
        share_lo * np.maximum(x_lo - own, 0.0)
        + share_hi * np.maximum(x_hi - own, 0.0)
        + share_both * np.maximum(np.maximum(x_lo, x_hi) - own, 0.0)
    )
    return own + np.divide(excess, share, out=np.zeros_like(excess), where=share > 0.0)


def _junction_prices(
    A, B, right_of, c_region, c_boundary, gdna_contained, junction_cap
) -> tuple[np.ndarray, np.ndarray]:
    """``(price, crossing)`` at the junction between exon pieces ``A`` (below it, in genomic order) and ``B`` (above
    it). The price is by conservation of bases::

        c_junction = c_lo + c_hi − ½ (c_intron,lo + c_intron,hi)

    gDNA never deposits on a junction and capture is local, so the price comes from the objects within a fragment of
    it. A gDNA fragment across the low boundary (the end of the exon below) holds the junction fragment's low-exon
    bases plus intron, one across the high boundary its high-exon bases plus intron, and averaged over positions each
    carries half a fragment of intron; where capture adds over bases, that is the two boundaries less the intron
    pieces beside them.

    The price is at least 0 and at most the larger of ``junction_cap`` and the best object the junction touches (its
    two pieces and two boundaries). A fragment binds its best single probe part, so it is captured no more than one
    inside a fully captured piece; without the cap, probes tiled to the exons' edges would count the exon's capture
    about twice. The cap lets a probe spanning the junction (one part for spliced RNA, two for gDNA) out-capture every
    gDNA object beside it; the best-object term lets a junction beside a stronger probe reach that object's level.
    The cap is the library's (:func:`fit_junction_cap`); every other term is the junction's own.

    The crossing level is the price without the best-object term, ``min(summed, junction_cap)``: what RNA's crossing
    captures at the junction, which a placement that leaves a cut's exons across it takes with it
    (:func:`_best_part`)."""
    summed, best = _junction_parts(A, B, right_of, c_region, c_boundary, gdna_contained)
    return np.minimum(summed, np.maximum(junction_cap, best)), np.minimum(summed, junction_cap)


def _junction_parts(
    A, B, right_of, c_region, c_boundary, gdna_contained
) -> tuple[np.ndarray, np.ndarray]:
    """``(summed, best)``: conservation's price at the junction between pieces ``A`` and ``B``, floored at 0, and the
    best object the junction touches. :func:`_junction_prices` is ``min(summed, max(cap, best))``."""
    junction_lo, junction_hi = right_of[A], right_of[B - 1]
    intron_lo, intron_hi = A + 1, B - 1
    # an intron piece too short to contain a gDNA fragment reads the boundary on its far side
    c_intron_lo = np.where(
        gdna_contained[intron_lo] > 0.0, c_region[intron_lo], c_boundary[right_of[intron_lo]]
    )
    c_intron_hi = np.where(
        gdna_contained[intron_hi] > 0.0,
        c_region[intron_hi],
        c_boundary[right_of[intron_hi - 1]],
    )
    summed = np.maximum(
        c_boundary[junction_lo] + c_boundary[junction_hi] - 0.5 * (c_intron_lo + c_intron_hi), 0.0
    )
    best = np.maximum.reduce(
        [c_region[A], c_region[B], c_boundary[junction_lo], c_boundary[junction_hi]]
    )
    return summed, best


def fit_junction_cap(calibration: "CalibrationResult", frame: TranscriptFrame) -> float:
    """The level a junction's price is capped at (:func:`_junction_prices`): the saturated level calibration publishes
    (the landscape's, at the capture weight), raised as far as RNA's own spliced fragments support.

    The saturated level is gDNA's, and gDNA never deposits on a junction; RNA does. In an exon region RNA and gDNA are
    the same placements under the same probes, so a transcript's fragments split between a junction and an exon region
    beside it as the ruler's own terms say: its conserved share at the junction times the junction's price, against its
    contained share of the region times the region's efficiency, each as the ruler carries it (:func:`_carried`). With
    ``W`` a junction's conserved spliced mass and ``m`` the unspliced RNA calibration booked in the exon regions beside
    it::

        W | W + m  ~  Binomial( a(C) / (a(C) + D) ),      a(C) = mean_t [ X_t + X'_t·(p_t(C) − 1) ]

    over the annotated transcripts ``t`` using the junction: ``X_t`` is the transcript's share there and ``X'_t`` that
    share times the carry factor (``X' = X`` where nothing is carried, so the term is ``X·p``); ``p_t(C)`` is the price
    the ruler gives the transcript's cut at cap ``C``, the junction rule ``min(summed, max(C, best))`` with the
    placements that leave across the next junction at its crossing level ``min(summed, C)`` where that is higher
    (:func:`_best_part`); ``D`` is the ruler's term of those regions at their efficiencies. The cap is the ``C`` that
    maximises the junctions' joint likelihood on a geometric grid from the saturated level to the largest uncapped
    price the fit reads, the first maximum taken. The counts can raise the cap and never lower it, so where they ask
    for no more the rule is the saturated level's, bit for bit.

    A junction's evidence is its own footprint. An exon region beside it enters only when every annotated transcript
    covering the region uses the junction, so no other junction's and no other isoform's fragments are read; a
    junction is a row of the sj axis, its two regions and its strand. A junction with no such region, with no fragment
    in its footprint, or whose term the cap cannot move (neither its own price nor the crossing level of a junction
    its placements leave across) is not in the fit, and with none in it the cap is the saturated level. So it is at a
    saturated level of exactly 1 (no level above the bulk, or a capture weight of 0): no capture is read, and there is
    none for a junction to exceed.
    """
    from .splice_graph import build_sj_geometry_arrays

    saturation = float(calibration.gdna_capture_saturation)
    if saturation == 1.0:
        return saturation
    index, region_arrays = frame.index, frame.region_arrays
    c_region = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    c_boundary = np.asarray(calibration.gdna_capture_efficiency_boundary, dtype=np.float64)
    obj = frame.rna
    # the ruler's own terms: a piece enters at its share times its region's efficiency as carried, and a junction
    # priced p above 1 at its share plus the carried excess, ``share + slope·(p − 1)`` (:func:`_carried`)
    c_piece = c_region[obj.piece]
    cut_slope = obj.cut_share
    law = obj
    carry = frame.carry
    if carry is not None:
        on_anchor, on_table = carry
        c_piece = _carried(c_piece, on_anchor.contained, on_table.contained)
        cut_slope = obj.cut_share * _carry_factor(on_anchor.cut_share, on_table.cut_share)
        law = on_anchor
    piece_term = obj.contained * c_piece

    # each annotated transcript's junction cuts, on the sj axis
    geo = build_sj_geometry_arrays(index)
    src = np.asarray(geo.src_region, dtype=np.int64)
    dst = np.asarray(geo.dst_region, dtype=np.int64)
    n_regions, n_sj = int(region_arrays.n_regions), src.shape[0]
    synthetic = index.t_df["is_synthetic"].to_numpy(dtype=bool)
    strand = index.t_df["strand"].to_numpy(dtype=np.int64)
    cut_t = obj.t[obj.cut_row]
    sj_of_cut = _sj_rows(
        (src, dst, np.asarray(geo.strand, dtype=np.int64)),
        (obj.piece[obj.cut_row], obj.piece[obj.cut_row + 1], strand[cut_t]),
        n_regions,
    )
    used = obj.is_junction & ~synthetic[cut_t] & (sj_of_cut >= 0)
    users = np.bincount(sj_of_cut[used], minlength=n_sj)

    # the junction's own footprint: an exon region beside it that only its transcripts cover. The mask can cut a
    # transcript's piece of a region into parts, so a transcript counts once per region and its term there is the sum
    # over its parts; every transcript covering a region holds the same observable parts of it, so they agree.
    annotated = ~synthetic[obj.t]
    pair, part_of = np.unique(
        obj.t[annotated] * n_regions + obj.piece[annotated], return_inverse=True
    )
    cover = np.bincount(pair % n_regions, minlength=n_regions)
    region_term = np.zeros(n_regions)
    np.maximum.at(
        region_term, pair % n_regions, np.bincount(part_of, weights=piece_term[annotated])
    )
    rna = np.maximum(np.asarray(calibration.count_rna_region, dtype=np.float64), 0.0)
    anchor = np.zeros(n_sj)
    anchor_term = np.zeros(n_sj)
    for side in (src, dst):
        own = (users > 0) & (cover[side] == users) & (region_term[side] > 0.0)
        anchor += np.where(own, rna[side], 0.0)
        anchor_term += np.where(own, region_term[side], 0.0)
    spliced = np.asarray(calibration.sj_conserved_mass, dtype=np.float64)

    # every junction cut as the ruler prices it: its own rule's two bounds, and the junctions its placements leave
    # across (:func:`_best_part`), as rows of the junction cuts
    cuts = np.flatnonzero(obj.is_junction)
    at = np.full(obj.cut_row.shape[0], -1, dtype=np.int64)
    at[cuts] = np.arange(cuts.shape[0])
    summed, best = _junction_parts(
        obj.piece[obj.cut_row][cuts],
        obj.piece[obj.cut_row + 1][cuts],
        region_right_boundary(np.asarray(region_arrays.ref_id)),
        c_region,
        c_boundary,
        np.asarray(calibration.gdna_region_eff_len, dtype=np.float64),
    )
    leave = [
        np.where(to[cuts] >= 0, at[np.maximum(to[cuts], 0)], -1)
        for to in (law.leave_lo, law.leave_hi)
    ]
    leaving = (law.leave_lo_share[cuts], law.leave_hi_share[cuts], law.leave_both_share[cuts])
    slope, sj, counted = cut_slope[cuts], sj_of_cut[cuts], used[cuts]

    # a junction is in the fit when the cap can move its term: its own price (conservation's sum above the saturated
    # level and the best object beside it), or the crossing level of a junction its placements leave across (a sum
    # above the saturated level)
    moves = summed > np.maximum(saturation, best)
    reads = [moves] + [
        (to >= 0) & (part + leaving[2] > 0.0) & (summed[np.maximum(to, 0)] > saturation)
        for to, part in zip(leave, leaving[:2], strict=True)
    ]
    moved = np.bincount(sj[counted], weights=np.logical_or.reduce(reads)[counted], minlength=n_sj)
    held = np.bincount(sj[counted], weights=slope[counted], minlength=n_sj)
    fit = (moved > 0) & (held > 0.0) & (anchor_term > 0.0) & (spliced + anchor > 0.0)
    in_fit = np.flatnonzero(counted & fit[np.maximum(sj, 0)])
    if in_fit.size == 0:
        return saturation
    # the grid ends at the largest sum the fit's junctions read, their own or one they leave across
    read = np.zeros(cuts.shape[0], dtype=bool)
    read[in_fit[reads[0][in_fit]]] = True
    for to, reaches in zip(leave, reads[1:], strict=True):
        read[to[in_fit[reaches[in_fit]]]] = True
    levels = np.geomspace(
        saturation, float(summed[read].max()), CONSTANTS.calibration.junction_cap_grid_points
    )

    # the fit's own cuts: each one's bounds, the sums of the junctions it leaves across, its shares and its sj
    own_sum, own_best = summed[in_fit], best[in_fit]
    sum_lo, sum_hi = (_across(summed, to)[in_fit] for to in leave)
    parts = [part[in_fit] for part in (*leaving, law.cut_share[cuts])]
    share, slope, sj = obj.cut_share[cuts][in_fit], slope[in_fit], sj[in_fit]
    W, m, D, n = spliced[fit], anchor[fit], anchor_term[fit], users[fit]
    log_likelihood = np.empty(levels.shape[0])
    for i, level in enumerate(levels):
        # the ruler's term at each cut with the cap at this level, the excess over 1 carried, and its mean over the
        # annotated transcripts using the junction
        price = _best_part(
            np.minimum(own_sum, np.maximum(level, own_best)),
            np.minimum(sum_lo, level),
            np.minimum(sum_hi, level),
            *parts,
        )
        term = np.where(price > 1.0, share + slope * (price - 1.0), share * price)
        a = np.bincount(sj, weights=term, minlength=n_sj)[fit] / n
        log_likelihood[i] = (
            W * (np.log(a) - np.log(a + D)) + m * (np.log(D) - np.log(a + D))
        ).sum()
    return float(levels[int(np.argmax(log_likelihood))])


def _sj_rows(sj: tuple, cut: tuple, n_regions: int) -> np.ndarray:
    """The sj row of each cut, ``-1`` where it has none: the row with the cut's two regions and its transcript's
    strand. ``sj`` and ``cut`` are ``(low region, high region, strand)``."""

    def key(lo, hi, strand):
        return (lo * (n_regions + 1) + hi) * len(Strand) + strand

    sj_key = key(*sj)
    if sj_key.shape[0] == 0:
        return np.full(cut[0].shape[0], -1, dtype=np.int64)
    order = np.argsort(sj_key, kind="stable")
    cut_key = key(*cut)
    at = np.minimum(np.searchsorted(sj_key[order], cut_key), sj_key.shape[0] - 1)
    return np.where(sj_key[order][at] == cut_key, order[at], -1)


def place_efficiencies(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    index: "TranscriptIndex",
    junction_cap: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``(c_region, c_boundary, c_sj, c_sj_crossing, sj_low_boundary, sj_high_boundary)``: the efficiencies a
    multimapper's place weight reads (:func:`place_log_weights`), on the accumulator's flat object axes: calibration's
    own at a region and a boundary, and at an annotated junction the price its own rule gives it
    (``_junction_prices`` at the same ``junction_cap``), from its two exon regions, with its crossing level, the part
    of that price a place across it passes on (:func:`_best_part`); and the boundary at each junction's low and high
    splice site, which orders a place's junctions among the boundaries it crosses. A junction reads 1 with no
    landscape, as every object does."""
    from .splice_graph import build_sj_geometry_arrays

    c_region = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    c_boundary = np.asarray(calibration.gdna_capture_efficiency_boundary, dtype=np.float64)
    geo = build_sj_geometry_arrays(index)
    src = np.asarray(geo.src_region, dtype=np.int64)
    dst = np.asarray(geo.dst_region, dtype=np.int64)
    right_of = region_right_boundary(np.asarray(region_arrays.ref_id))
    # every boundary's src lies below its dst (the graph's I8, checked on load), so dst - 1 is a region
    low, high = right_of[src], right_of[dst - 1]
    if calibration.gdna_bulk_density is None or src.size == 0:
        one = np.ones(src.size, dtype=np.float64)
        return c_region, c_boundary, one, one, low, high
    c_sj, c_sj_crossing = _junction_prices(
        src,
        dst,
        right_of,
        c_region,
        c_boundary,
        np.asarray(calibration.gdna_region_eff_len, dtype=np.float64),
        junction_cap,
    )
    return c_region, c_boundary, c_sj, c_sj_crossing, low, high


def place_log_weights(
    object_offsets: np.ndarray,
    object_kind: np.ndarray,
    object_index: np.ndarray,
    object_mass: np.ndarray,
    *,
    intergenic: np.ndarray,
    c_region: np.ndarray,
    c_boundary: np.ndarray,
    c_sj: np.ndarray,
    c_sj_crossing: np.ndarray,
    sj_low_boundary: np.ndarray,
    sj_high_boundary: np.ndarray,
) -> np.ndarray:
    """``log w`` per buffered row: a multimapper's place's capture weight, ``w(x) = Σ_o m_o(x) c_o``, the deposit
    rule's mass of the place at every object it deposits on (CSR over rows, ``rigel.native.PLACE_*`` kinds, flat
    ids) times the object's capture efficiency (:func:`place_efficiencies`), a place across junctions binding its
    best part (:func:`_place_best_part`). Its sum over the gDNA component's
    placements, on the law the efficiencies are read on, is that component's capture-contracted length, and its excess
    over 1 summed over a transcript's is the excess of the transcript's length over its own (where every efficiency is
    at least 1), so every candidate at the place carries it.

    A row with no object (a unique fragment's, on which the weight is common to every candidate) reads 0, and so does
    an intergenic place, because intergenic sequence is uncaptured. float32, as the scorer reads it.
    """
    from ..native import PLACE_BOUNDARY, PLACE_REGION, PLACE_SJ

    offsets = np.asarray(object_offsets, dtype=np.int64)
    n_rows = offsets.size - 1
    kind = np.asarray(object_kind)
    index = np.asarray(object_index, dtype=np.int64)
    mass = np.asarray(object_mass, dtype=np.float64)
    c = np.ones(mass.size, dtype=np.float64)
    for code, eff in ((PLACE_REGION, c_region), (PLACE_BOUNDARY, c_boundary), (PLACE_SJ, c_sj)):
        sel = kind == code
        c[sel] = np.asarray(eff, dtype=np.float64)[index[sel]]
    per_row = np.diff(offsets)
    row = np.repeat(np.arange(n_rows), per_row)
    c = _place_best_part(row, kind, index, c, c_sj_crossing, sj_low_boundary, sj_high_boundary)
    w = np.zeros(n_rows, dtype=np.float64)
    np.add.at(w, row, mass * c)
    out = np.zeros(n_rows, dtype=np.float32)
    weighted = (per_row > 0) & ~np.asarray(intergenic, dtype=bool)
    with np.errstate(divide="ignore"):
        out[weighted] = np.log(w[weighted])
    return out


def _place_best_part(
    row: np.ndarray,
    kind: np.ndarray,
    index: np.ndarray,
    c: np.ndarray,
    c_sj_crossing: np.ndarray,
    sj_low_boundary: np.ndarray,
    sj_high_boundary: np.ndarray,
) -> np.ndarray:
    """Every deposit's efficiency with a place across junctions binding its best part: a boundary the place crosses
    and a junction it uses each read the larger of their own efficiency and the crossing level of the nearest junction
    the place uses on each side of them. It is :func:`_best_part` one placement at a time: the junction ending the
    object's home below is the nearest the place uses below it, and a place that stays home uses none there. ``row``
    is each object's row; a place's objects lie on one reference, where a boundary's index orders them."""
    from ..native import PLACE_BOUNDARY, PLACE_SJ

    is_sj = kind == PLACE_SJ
    j = np.flatnonzero(is_sj)
    if j.size == 0:
        return c
    low_of = np.asarray(sj_low_boundary, dtype=np.int64)
    high_of = np.asarray(sj_high_boundary, dtype=np.int64)
    o = np.flatnonzero(is_sj | (kind == PLACE_BOUNDARY))
    o_low = np.where(is_sj[o], low_of[np.where(is_sj[o], index[o], 0)], index[o])
    o_high = np.where(is_sj[o], high_of[np.where(is_sj[o], index[o], 0)], index[o])
    # the places' junctions, by row then position: one alignment's introns do not overlap, so either end orders them
    j = j[np.lexsort((low_of[index[j]], row[j]))]
    span = int(o_high.max()) + 2
    key_low = row[j] * span + low_of[index[j]]
    key_high = row[j] * span + high_of[index[j]]
    out = c.copy()
    crossing = np.asarray(c_sj_crossing, dtype=np.float64)[index[j]]
    last = j.size - 1
    # below: the last of the row's junctions that ends at or below the object
    at = np.clip(np.searchsorted(key_high, row[o] * span + o_low, side="right") - 1, 0, last)
    there = (row[j[at]] == row[o]) & (high_of[index[j[at]]] <= o_low)
    out[o] = np.maximum(out[o], np.where(there, crossing[at], 0.0))
    # above: the first of the row's junctions that starts at or above the object
    at = np.clip(np.searchsorted(key_low, row[o] * span + o_high, side="left"), 0, last)
    there = (row[j[at]] == row[o]) & (low_of[index[j[at]]] >= o_high)
    out[o] = np.maximum(out[o], np.where(there, crossing[at], 0.0))
    return out


def transcript_capture_eff_lengths(
    calibration: "CalibrationResult",
    frame: TranscriptFrame,
    fl_eff_lengths: np.ndarray,
    junction_cap: float,
) -> np.ndarray:
    """``eff_em_t = fl_t · Σ_o share_o(T_R) ĉ_o / Σ_o share_o(T_R)`` over the transcript's pieces and cuts.

    The efficiencies ``c_o`` are the result's (a piece's region's, a contiguous cut's boundary's, a junction's from the
    objects beside it, capped at ``junction_cap``, :func:`fit_junction_cap`; a cut's placements that leave its exons
    across a junction at that junction's crossing level where it is higher, :func:`_best_part`): gDNA's yields, read on gDNA's
    uncaptured law ``g`` (``tables.anchor``). Capture adds to the uncaptured yield, and the uncaptured part is RNA's own share on
    its own table, so only the excess is carried to the scorer's gDNA table ``T_G`` (``tables.gdna``):
    ``ĉ_o = 1 + (c_o − 1)·share_o(g)/share_o(T_G)`` where ``c_o > 1``, and ``c_o`` as read elsewhere and where ``T_G``
    cannot reach (:func:`_carried`). Every efficiency exactly 1 gives ``fl`` to the bit, whatever the tables; where
    ``T_G`` is ``g`` the efficiencies are used as read. With no landscape ``fl`` is returned unchanged, and a transcript
    with no share on any object (shorter than every fragment) keeps ``fl``.
    """
    fl = np.asarray(fl_eff_lengths, dtype=np.float64)
    if calibration.gdna_bulk_density is None:
        return fl.copy()
    c_region = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    c_boundary = np.asarray(calibration.gdna_capture_efficiency_boundary, dtype=np.float64)
    gdna_contained = np.asarray(calibration.gdna_region_eff_len, dtype=np.float64)
    obj = frame.rna
    c_piece = c_region[obj.piece]
    carry = frame.carry
    # a cut's placements are split on the law its efficiencies are read on, as its excess is carried from there
    c_cut = cut_efficiencies(
        obj,
        frame.region_arrays,
        c_region,
        c_boundary,
        gdna_contained,
        junction_cap,
        obj if carry is None else carry[0],
    )
    if carry is not None:
        on_anchor, on_table = carry
        c_piece = _carried(c_piece, on_anchor.contained, on_table.contained)
        c_cut = _carried(c_cut, on_anchor.cut_share, on_table.cut_share)
    num = np.zeros(fl.shape[0])
    den = np.zeros(fl.shape[0])
    np.add.at(num, obj.t, obj.contained * c_piece)
    np.add.at(den, obj.t, obj.contained)
    cut_t = obj.t[obj.cut_row]
    np.add.at(num, cut_t, obj.cut_share * c_cut)
    np.add.at(den, cut_t, obj.cut_share)
    factor = np.divide(num, den, out=np.ones_like(num), where=den > 0.0)
    return fl * factor


def _carry_shares(
    pieces: TranscriptPieces, obj: TranscriptObjects, tables: "LengthTables", n_threads: int = 1
) -> tuple[TranscriptObjects, TranscriptObjects] | None:
    """``(on_anchor, on_table)``: the objects' shares on gDNA's uncaptured law and on the scorer's gDNA table, by
    which an efficiency is carried to that table (:func:`_carried`); ``None`` where the table is the uncaptured law
    and nothing is carried. ``obj`` is the shares on RNA's table."""
    if np.array_equal(tables.gdna, tables.anchor):
        return None
    on_table = (
        obj
        if np.array_equal(tables.gdna, tables.rna)
        else object_shares(pieces, tables.gdna, n_threads)
    )
    return object_shares(pieces, tables.anchor, n_threads), on_table


def _carry_factor(share_anchor: np.ndarray, share_table: np.ndarray) -> np.ndarray:
    """``share(g) / share(T_G)`` per object, by which capture's excess is carried to the scorer's gDNA table; 1 where
    ``T_G`` has no share there."""
    reach = share_table > 0.0
    return np.where(reach, share_anchor / np.where(reach, share_table, 1.0), 1.0)


def _carried(c: np.ndarray, share_anchor: np.ndarray, share_table: np.ndarray) -> np.ndarray:
    """``ĉ`` per object. Capture adds to the uncaptured yield, and the uncaptured part is RNA's own share on its own
    table, so only an excess over 1 is carried: ``1 + (c − 1)·share(g)/share(T_G)``. A reading of 1 or below carries
    nothing and enters as read, and so does an object ``T_G`` has no share of."""
    carry = (c > 1.0) & (share_table > 0.0)
    return np.where(carry, 1.0 + (c - 1.0) * _carry_factor(share_anchor, share_table), c)


def rna_captured_fraction(calibration: "CalibrationResult", region_arrays: "RegionArrays") -> float:
    """``kappa``: the fraction of RNA's fragments that capture added, ``Σ_exons n_RNA (1 − 1/c)₊ / Σ_exons n_RNA``.

    ``c`` is each exon region's gDNA capture efficiency against the bulk: of an object's fragments at yield ``c`` the
    uncaptured background accounts for ``1/c`` and capture for the rest. RNA's count weights the exons, because the
    question is what capture did to RNA's census. 0 with no landscape, or with no RNA at an exon.
    """
    if calibration.gdna_bulk_density is None:
        return 0.0
    c = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    n = np.clip(np.asarray(calibration.count_rna_region, dtype=np.float64), 0.0, None)
    exon = coarse_type_array(region_arrays.signature) == int(RegionType.EXON)
    n = np.where(exon, n, 0.0)
    if n.sum() <= 0.0:
        return 0.0
    gain = np.where(c > 1.0, 1.0 - 1.0 / np.where(c > 1.0, c, 1.0), 0.0)
    return float((n * gain).sum() / n.sum())
