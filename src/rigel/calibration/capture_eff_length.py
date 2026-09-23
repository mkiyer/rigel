"""Capture-aware EM effective lengths — a transcript's conserved shares of its objects at their capture
efficiencies.

Under hybrid capture a transcript's usable length is not its full length: only the probe-enriched part
of its footprint is sampled, and contracting the gDNA component alone would leave it artificially
concentrated against the RNA components. So EVERY transcript's EM effective length is contracted by the
calibration's own efficiencies, read off the deconvolved gDNA against the fully captured level
(`capture_efficiency`; `CalibrationResult.gdna_capture_efficiency_region` / ``_boundary``) — on the same
objects the deposit rule puts the transcript's fragments on (gate:
``tests/calibration/test_capture_eff_length.py``).

THE LENGTH IS A SUM OVER OBJECTS, CONSERVED. A transcript's pieces are the regions its exons cover (an
unspliced entity's single exon is its span, introns included), in its own coordinates; between two
consecutive pieces is a CUT — the boundary between them, or a junction where the transcript splices. Its
placements put their unit on those objects exactly as the accumulator does: a contained placement on its
piece, a crossing one split over the cuts it crosses (`effective_length.conserved_cut_shares`). So

    eff_t = fl_t · (Σ_pieces S_p c_p + Σ_cuts M_k c_k) / (Σ_pieces S_p + Σ_cuts M_k),

``S_p`` the transcript's contained share of piece ``p`` and ``M_k`` its conserved share at cut ``k``; the
shares total the fl-marginal length exactly, however many cuts one fragment crosses, and the contraction
scales the fl-marginal ``fl_eff_lengths`` the pipeline built. A piece prices at its region's efficiency
and a contiguous cut at its boundary's: there the transcript's fragments are gDNA's.

A JUNCTION IS PRICED BY CONSERVATION OF BASES. gDNA never deposits on a junction, and capture is local, so
a junction's efficiency comes from the objects within a fragment of it. In genomic order, a fragment across
a junction holds bases of the exon below it and of the exon above it; a gDNA fragment across the junction's
LOW boundary (the end of the exon below) holds the same low-exon bases and intron bases, one across its HIGH
boundary (the start of the exon above) intron bases and the same high-exon bases, and averaged over the
crossing positions each carries half a fragment of intron. Where capture adds over bases,

    c_junction = c_lo + c_hi − ½ (c_intron,lo + c_intron,hi),

the junction's two boundaries less the intron pieces beside them — a piece too short to contain a fragment
has no contained count, and reads the boundary on its far side instead. The sum is never below 0, and it may
exceed 1: the two sides' exon capture is added, not averaged, and each is a boundary's unclipped half-exon
level. A field with no reference — every efficiency 1 — returns ``fl_eff_lengths``
bit-identically; only the captured case contracts.
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

if TYPE_CHECKING:
    from ..index import TranscriptIndex
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

__all__ = ["TranscriptObjects", "transcript_capture_eff_lengths", "transcript_objects"]


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


def transcript_objects(
    index: "TranscriptIndex", region_arrays: "RegionArrays", rna_fl_pmf: np.ndarray
) -> TranscriptObjects:
    """Every transcript's :class:`TranscriptObjects` on the RNA pmf: the regions its exons (or an unspliced
    entity's span) cover, and between consecutive ones a cut — a junction where they belong to different
    exons with a region between them. A zero-length fragment places nowhere — the accumulator deposits
    none — so the pmf's mass at ``w = 0`` is dropped, and the shares partition the fl-marginal length the
    pipeline builds, which excludes it too."""
    pmf = np.asarray(rna_fl_pmf, dtype=np.float64).copy()
    pmf[0] = 0.0
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
    if tdf is not None and "is_synthetic" in tdf.columns:
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
    length = (np.minimum(ends[piece], b[k]) - np.maximum(starts[piece], a[k])).astype(np.float64)

    # transcript coordinates: the bases left of each cut and right of it
    n_t = int(tdf.shape[0])
    total = np.zeros(n_t)
    np.add.at(total, tk, length)
    first = np.r_[True, tk[1:] != tk[:-1]] if tk.size else np.zeros(0, dtype=bool)
    csum = np.cumsum(length)
    through = csum - np.repeat(
        csum[first] - length[first], np.diff(np.r_[np.flatnonzero(first), tk.size])
    )
    cut = np.flatnonzero(np.r_[tk[1:] == tk[:-1], False]) if tk.size else np.zeros(0, np.int64)
    left, right = conserved_cut_shares(
        pmf, length[cut], length[cut + 1], through[cut], total[tk[cut]] - through[cut]
    )
    return TranscriptObjects(
        t=tk,
        piece=piece,
        contained=contained_eff_length(length, pmf),
        cut_row=cut,
        cut_share=left + right,
        is_junction=(k[cut] != k[cut + 1]) & (piece[cut + 1] > piece[cut] + 1),
        n_transcripts=n_t,
    )


def _cut_efficiencies(
    objects: TranscriptObjects,
    region_arrays: "RegionArrays",
    c_region: np.ndarray,
    c_boundary: np.ndarray,
    gdna_contained: np.ndarray,
) -> np.ndarray:
    """Every cut's efficiency: a contiguous cut's boundary's; a junction's by conservation of bases,
    ``c_lo + c_hi − ½(c_intron,lo + c_intron,hi)``, never below 0 — ``A`` the piece below the junction and
    ``B`` the piece above it, in genomic order."""
    right_of = region_right_boundary(np.asarray(region_arrays.ref_id))
    A = objects.piece[objects.cut_row]
    B = objects.piece[objects.cut_row + 1]
    out = c_boundary[np.maximum(right_of[A], 0)]
    j = np.flatnonzero(objects.is_junction)
    if j.size:
        junction_lo, junction_hi = right_of[A[j]], right_of[B[j] - 1]
        intron_lo, intron_hi = A[j] + 1, B[j] - 1
        # an intron piece too short to contain a gDNA fragment reads the boundary on its far side
        c_intron_lo = np.where(
            gdna_contained[intron_lo] > 0.0, c_region[intron_lo], c_boundary[right_of[intron_lo]]
        )
        c_intron_hi = np.where(
            gdna_contained[intron_hi] > 0.0,
            c_region[intron_hi],
            c_boundary[right_of[intron_hi - 1]],
        )
        out = out.copy()
        out[j] = np.maximum(
            c_boundary[junction_lo] + c_boundary[junction_hi] - 0.5 * (c_intron_lo + c_intron_hi),
            0.0,
        )
    return out


def transcript_capture_eff_lengths(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    index: "TranscriptIndex",
    fl_eff_lengths: np.ndarray,
    rna_fl_pmf: np.ndarray,
) -> np.ndarray:
    """``eff_em_t = fl_t · Σ_o share_o c_o / Σ_o share_o`` over the transcript's pieces and cuts (gate:
    ``tests/calibration/test_capture_eff_length.py``).

    The efficiencies are the result's: a piece's region's, a contiguous cut's boundary's, a junction's from
    the objects beside it. With no reference (capture off, or no gDNA) every efficiency is 1 and ``fl`` is
    returned verbatim; a transcript with no share on any object (shorter than every fragment) keeps ``fl``.
    A piece's and a boundary's efficiency lie in ``[0, 1]``; a junction's adds its two sides' exon capture
    and can exceed 1, so ``eff_em`` can exceed ``fl`` where well-captured junctions carry the transcript.
    """
    fl = np.asarray(fl_eff_lengths, dtype=np.float64)
    if calibration.gdna_reference_density is None:
        return fl.copy()
    c_region = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    c_boundary = np.asarray(calibration.gdna_capture_efficiency_boundary, dtype=np.float64)
    gdna_contained = np.asarray(calibration.gdna_region_eff_len, dtype=np.float64)
    obj = transcript_objects(index, region_arrays, rna_fl_pmf)
    c_cut = _cut_efficiencies(obj, region_arrays, c_region, c_boundary, gdna_contained)
    num = np.zeros(fl.shape[0])
    den = np.zeros(fl.shape[0])
    np.add.at(num, obj.t, obj.contained * c_region[obj.piece])
    np.add.at(den, obj.t, obj.contained)
    cut_t = obj.t[obj.cut_row]
    np.add.at(num, cut_t, obj.cut_share * c_cut)
    np.add.at(den, cut_t, obj.cut_share)
    factor = np.divide(num, den, out=np.ones_like(num), where=den > 0.0)
    return fl * factor
