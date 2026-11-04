"""Capture-aware EM effective lengths: a transcript's conserved shares of its objects at their capture efficiencies.

Layer 2 (opportunity). Under hybrid capture only the probe-enriched part of a transcript's footprint is sampled, and
contracting gDNA alone would concentrate it against RNA, so every transcript's EM length is contracted by
calibration's own per-object gDNA capture efficiencies, on the objects the deposit rule uses. A transcript's pieces
are the regions its exons cover (an unspliced entity's single exon is its span); between consecutive pieces is a cut,
a boundary or a junction. With ``S_p`` its contained share of piece ``p`` and ``M_k`` its conserved share at cut ``k``
(``effective_length.conserved_cut_shares``), which together total the fl-marginal length exactly,

    eff_t = fl_t · (Σ_pieces S_p c_p + Σ_cuts M_k c_k) / (Σ_pieces S_p + Σ_cuts M_k).

A piece prices at its region's efficiency, a contiguous cut at its boundary's, a junction by conservation of bases
(:func:`_junction_prices`). With no landscape, or every efficiency exactly 1, ``fl_eff_lengths`` comes back
bit-identically. Efficiencies are yields against the bulk, so a captured transcript's length exceeds ``fl``, and a
factor common to every component cancels in the EM. Pieces are cut at the index's unobservable mask; a template the
mask does not touch is unchanged to the bit.

Gate: ``tests/calibration/test_capture_eff_length.py``.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..types import IntervalType
from .effective_length import conserved_cut_shares, contained_eff_length
from .region_arrays import overlapping_region_runs, region_right_boundary
from .signature import RegionType, coarse_type_array

if TYPE_CHECKING:
    from ..index import TranscriptIndex
    from .region_arrays import RegionArrays
    from .result import CalibrationResult
    from .shared_length_law import LengthTables

__all__ = [
    "TranscriptObjects",
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
    share, ``is_junction`` whether the transcript splices there. Annotation and pmf only.
    """

    t: np.ndarray
    piece: np.ndarray
    contained: np.ndarray
    cut_row: np.ndarray
    cut_share: np.ndarray
    is_junction: np.ndarray
    n_transcripts: int


@dataclass(frozen=True, slots=True)
class _TranscriptPieces:
    """The pmf-free half of :class:`TranscriptObjects`: every transcript's pieces, each piece's length, and at each
    cut the bases of the transcript left of it (``through``) and right of it (``rest``), counted to the nearest
    masked base. ``run_t`` and ``run_length`` are each transcript's runs — its stretches no masked base interrupts —
    and ``masked`` the bases the mask took from each transcript."""

    t: np.ndarray
    piece: np.ndarray
    length: np.ndarray
    cut_row: np.ndarray
    through: np.ndarray
    rest: np.ndarray
    is_junction: np.ndarray
    n_transcripts: int
    run_t: np.ndarray
    run_length: np.ndarray
    masked: np.ndarray


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
    # the first masked run ending past each piece's start, on one ascending key: reference, then position
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
    keep = np.flatnonzero(~hit)
    row = np.concatenate([keep, np.asarray(rows, dtype=np.int64)])
    pa = np.concatenate([a[keep], np.asarray(lo, dtype=np.int64)])
    pb = np.concatenate([b[keep], np.asarray(hi, dtype=np.int64)])
    order = np.lexsort((pa, row))
    row, pa, pb = row[order], pa[order], pb[order]
    return t[row], piece[row], k[row], (pb - pa).astype(np.float64), u[row] + (pa - a[row])


def _transcript_pieces(
    index: "TranscriptIndex", region_arrays: "RegionArrays"
) -> _TranscriptPieces:
    """The observable parts of the regions every transcript's exons (or an unspliced entity's span) cover, in
    transcript then genomic order, and between consecutive parts a cut — a junction where they belong to different
    exons with a region between them. No part holds a base of the index's unobservable mask
    (``region_arrays.unobservable``) and no cut joins two parts a masked base separates: a fragment touching one is
    never observed. With nothing masked the parts are the pieces. Annotation only: no length law enters."""
    starts = np.asarray(region_arrays.start, dtype=np.int64)
    ends = np.asarray(region_arrays.end, dtype=np.int64)
    ref_off = np.asarray(region_arrays.ref_offsets, dtype=np.int64)

    iv = pd.read_feather(os.path.join(index.index_dir, "intervals.feather"))
    ex = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]
    parts = [
        (
            ex["t_index"].to_numpy(np.int64),
            ex["ref"].astype(str).to_numpy(),
            ex["start"].to_numpy(np.int64),
            ex["end"].to_numpy(np.int64),
        )
    ]
    tdf = index.t_df
    syn = tdf[tdf["is_synthetic"].to_numpy(dtype=bool)]
    syn = syn[~np.isin(syn["t_index"].to_numpy(np.int64), parts[0][0])]
    parts.append(
        (
            syn["t_index"].to_numpy(np.int64),
            syn["ref"].astype(str).to_numpy(),
            syn["start"].to_numpy(np.int64),
            syn["end"].to_numpy(np.int64),
        )
    )
    t = np.concatenate([p[0] for p in parts])
    ref = np.concatenate([p[1] for p in parts])
    a = np.concatenate([p[2] for p in parts])
    b = np.concatenate([p[3] for p in parts])
    order = np.lexsort((a, t))
    t, ref, a, b = t[order], ref[order], a[order], b[order]
    rid = pd.Series(ref).astype(str).map(index.ref_name_to_id).fillna(-1).to_numpy(np.int64)
    lo, hi = overlapping_region_runs(rid, a, b, starts, ends, ref_off)
    width = np.where((rid >= 0) & (hi > lo), hi - lo, 0)

    # expand each exon into the regions it covers
    k = np.repeat(np.arange(t.size), width)
    piece = np.repeat(lo, width) + (
        np.arange(int(width.sum()), dtype=np.int64) - np.repeat(np.cumsum(width) - width, width)
    )
    tk = t[k]
    seg_a = np.maximum(starts[piece], a[k])
    seg_b = np.minimum(ends[piece], b[k])
    length = (seg_b - seg_a).astype(np.float64)

    # each piece's start in the transcript's full coordinates, masked bases included
    n_t = int(tdf.shape[0])
    whole = np.zeros(n_t)
    np.add.at(whole, tk, length)
    first = np.r_[True, tk[1:] != tk[:-1]] if tk.size else np.zeros(0, dtype=bool)
    csum = np.cumsum(length)
    full_start = (
        csum
        - np.repeat(csum[first] - length[first], np.diff(np.r_[np.flatnonzero(first), tk.size]))
        - length
    )
    tk, piece, k, length, full_start = _observable_parts(
        tk, piece, k, rid[k], seg_a, seg_b, full_start, region_arrays.unobservable
    )

    # the runs: consecutive parts of one transcript with no masked base between them; a cut joins two parts of one run
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
    kept = np.zeros(n_t)
    np.add.at(kept, tk, length)
    return _TranscriptPieces(
        t=tk,
        piece=piece,
        length=length,
        cut_row=cut,
        through=full_end[cut] - run_lo[run_id[cut]],
        rest=run_hi[run_id[cut]] - full_end[cut],
        is_junction=(k[cut] != k[cut + 1]) & (piece[cut + 1] > piece[cut] + 1),
        n_transcripts=n_t,
        run_t=tk[run_first],
        run_length=run_hi - run_lo,
        masked=whole - kept,
    )


def observable_fl_eff_lengths(
    index: "TranscriptIndex", region_arrays: "RegionArrays", rna_fl
) -> np.ndarray:
    """Every template's fl-marginal length over the starts the index's mask leaves: ``Σ_runs Σ_w f(w) (run − w + 1)₊``
    over its runs, the stretches no masked base interrupts, since a fragment touching a masked base is never observed.
    A template the mask does not touch keeps the closed form of its whole length, to the bit; every template keeps the
    closed form's floor of one base, as one shorter than every fragment does."""
    lengths = index.t_df["length"].to_numpy(np.int64)
    eff = rna_fl.compute_all_transcript_eff_lens(lengths)
    pieces = _transcript_pieces(index, region_arrays)
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
    eff[touched] = np.maximum(over_runs[touched], 1.0)
    return eff


def _objects(pieces: _TranscriptPieces, pmf: np.ndarray) -> TranscriptObjects:
    """The pieces' contained shares and the cuts' conserved shares on one pmf. A zero-length fragment places
    nowhere — the accumulator deposits none — so the pmf's mass at ``w = 0`` is dropped, and the shares partition
    the fl-marginal length the pipeline builds, which excludes it too."""
    pmf = np.asarray(pmf, dtype=np.float64).copy()
    pmf[0] = 0.0
    cut = pieces.cut_row
    left, right = conserved_cut_shares(
        pmf, pieces.length[cut], pieces.length[cut + 1], pieces.through, pieces.rest
    )
    return TranscriptObjects(
        t=pieces.t,
        piece=pieces.piece,
        contained=contained_eff_length(pieces.length, pmf),
        cut_row=cut,
        cut_share=left + right,
        is_junction=pieces.is_junction,
        n_transcripts=pieces.n_transcripts,
    )


def _cut_efficiencies(
    objects: TranscriptObjects,
    region_arrays: "RegionArrays",
    c_region: np.ndarray,
    c_boundary: np.ndarray,
    gdna_contained: np.ndarray,
    saturation: float,
) -> np.ndarray:
    """Every cut's efficiency: a contiguous cut's boundary's, a junction's from :func:`_junction_prices`. ``A`` is the
    piece below the cut and ``B`` the piece above it, in genomic order."""
    right_of = region_right_boundary(np.asarray(region_arrays.ref_id))
    A = objects.piece[objects.cut_row]
    B = objects.piece[objects.cut_row + 1]
    out = c_boundary[np.maximum(right_of[A], 0)]
    j = np.flatnonzero(objects.is_junction)
    if j.size:
        out = out.copy()
        out[j] = _junction_prices(
            A[j], B[j], right_of, c_region, c_boundary, gdna_contained, saturation
        )
    return out


def _junction_prices(
    A, B, right_of, c_region, c_boundary, gdna_contained, saturation
) -> np.ndarray:
    """The junction between exon pieces ``A`` (below it, in genomic order) and ``B`` (above it), by conservation of
    bases::

        c_junction = c_lo + c_hi − ½ (c_intron,lo + c_intron,hi)

    gDNA never deposits on a junction and capture is local, so the price comes from the objects within a fragment of
    it. A gDNA fragment across the low boundary (the end of the exon below) holds the junction fragment's low-exon
    bases plus intron, one across the high boundary its high-exon bases plus intron, and averaged over positions each
    carries half a fragment of intron; where capture adds over bases, that is the two boundaries less the intron
    pieces beside them.

    The price is at least 0 and at most the larger of the landscape's saturated level and the best object the junction
    touches (its two pieces and two boundaries). A fragment binds its best single probe part, so it is captured no more
    than one inside a fully captured piece; without the cap, probes tiled to the exons' edges would count the exon's
    capture about twice. The saturation term lets a probe spanning the junction (one part for spliced RNA, two for
    gDNA) out-capture every gDNA object beside it; the best-object term lets a junction beside a stronger probe reach
    that object's level."""
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
    bound = np.maximum(
        saturation,
        np.maximum.reduce(
            [c_region[A], c_region[B], c_boundary[junction_lo], c_boundary[junction_hi]]
        ),
    )
    return np.clip(
        c_boundary[junction_lo] + c_boundary[junction_hi] - 0.5 * (c_intron_lo + c_intron_hi),
        0.0,
        bound,
    )


def place_efficiencies(
    calibration: "CalibrationResult", region_arrays: "RegionArrays", index: "TranscriptIndex"
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(c_region, c_boundary, c_sj)``: the efficiencies a multimapper's place weight reads
    (:func:`place_log_weights`), on the accumulator's flat object axes: calibration's own at a region and a boundary,
    and at an annotated junction the price the transcripts' lengths give it (``_junction_prices``), from its two exon
    regions. A junction reads 1 with no landscape, as every object does."""
    from .splice_graph import build_sj_geometry_arrays

    c_region = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    c_boundary = np.asarray(calibration.gdna_capture_efficiency_boundary, dtype=np.float64)
    geo = build_sj_geometry_arrays(index)
    src = np.asarray(geo.src_region, dtype=np.int64)
    if calibration.gdna_bulk_density is None or src.size == 0:
        return c_region, c_boundary, np.ones(src.size, dtype=np.float64)
    c_sj = _junction_prices(
        src,
        np.asarray(geo.dst_region, dtype=np.int64),
        region_right_boundary(np.asarray(region_arrays.ref_id)),
        c_region,
        c_boundary,
        np.asarray(calibration.gdna_region_eff_len, dtype=np.float64),
        calibration.gdna_capture_saturation,
    )
    return c_region, c_boundary, c_sj


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
) -> np.ndarray:
    """``log w`` per buffered row: a multimapper's place's capture weight, ``w(x) = Σ_o m_o(x) c_o``, the deposit
    rule's mass of the place at every object it deposits on (CSR over rows, ``rigel.native.PLACE_*`` kinds, flat
    ids) times the object's capture efficiency (:func:`place_efficiencies`). Its sum over a component's placements,
    on the law the efficiencies are read on, is the component's capture-contracted length, so every candidate at the
    place carries it.

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
    w = np.zeros(n_rows, dtype=np.float64)
    np.add.at(w, np.repeat(np.arange(n_rows), per_row), mass * c)
    out = np.zeros(n_rows, dtype=np.float32)
    weighted = (per_row > 0) & ~np.asarray(intergenic, dtype=bool)
    with np.errstate(divide="ignore"):
        out[weighted] = np.log(w[weighted])
    return out


def transcript_capture_eff_lengths(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    index: "TranscriptIndex",
    fl_eff_lengths: np.ndarray,
    tables: "LengthTables",
) -> np.ndarray:
    """``eff_em_t = fl_t · Σ_o share_o(T_R) ĉ_o / Σ_o share_o(T_R)`` over the transcript's pieces and cuts.

    The efficiencies ``c_o`` are the result's (a piece's region's, a contiguous cut's boundary's, a junction's from the
    objects beside it): gDNA's yields, read on gDNA's uncaptured law ``g`` (``tables.anchor``). Each is carried to the
    scorer's gDNA table ``T_G`` (``tables.gdna``) as ``ĉ_o = c_o · share_o(g) / share_o(T_G)``. Where the two tables
    are one law the per-fragment capture weight cancels and this is ``Σ_o c_o share_o(g)``, the null model exactly;
    where ``T_G`` is ``g`` the efficiencies are used as read. An object ``T_G`` cannot reach keeps ``c_o``. With no
    landscape ``fl`` is returned unchanged, and a transcript with no share on any object (shorter than every fragment)
    keeps ``fl``.
    """
    fl = np.asarray(fl_eff_lengths, dtype=np.float64)
    if calibration.gdna_bulk_density is None:
        return fl.copy()
    c_region = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    c_boundary = np.asarray(calibration.gdna_capture_efficiency_boundary, dtype=np.float64)
    gdna_contained = np.asarray(calibration.gdna_region_eff_len, dtype=np.float64)
    pieces = _transcript_pieces(index, region_arrays)
    obj = _objects(pieces, tables.rna)
    c_piece = c_region[obj.piece]
    c_cut = _cut_efficiencies(
        obj,
        region_arrays,
        c_region,
        c_boundary,
        gdna_contained,
        calibration.gdna_capture_saturation,
    )
    if not np.array_equal(tables.gdna, tables.anchor):
        on_anchor = _objects(pieces, tables.anchor)
        on_table = obj if np.array_equal(tables.gdna, tables.rna) else _objects(pieces, tables.gdna)
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


def _carried(c: np.ndarray, share_anchor: np.ndarray, share_table: np.ndarray) -> np.ndarray:
    """``c · share(g) / share(T_G)`` per object; ``c`` where ``T_G`` has no share there."""
    reach = share_table > 0.0
    return np.where(reach, c * share_anchor / np.where(reach, share_table, 1.0), c)


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
