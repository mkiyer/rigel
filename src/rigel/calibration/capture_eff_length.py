"""Capture-aware EM effective lengths — a transcript's own bases at their pieces' capture efficiencies.

Under hybrid capture a transcript's usable length is not its full length: only the probe-enriched part
of its footprint is sampled, and contracting the gDNA component alone would leave it artificially
concentrated against the RNA components. So EVERY transcript's EM effective length is contracted by the
same per-piece efficiencies — the calibration's own, read off the deconvolved gDNA against the fully
captured level (`capture_efficiency`; `CalibrationResult.gdna_capture_efficiency_region`) — over the
transcript's own bases (gate: ``tests/calibration/test_capture_eff_length.py``).

THE LENGTH IS A SUM OVER BASES. A transcript's effective length under capture is the sum over its start
positions of the efficiency of the fragment starting there, and a fragment's efficiency is the mean
per-base efficiency over the bases it covers, so

    eff_t = Σ_x c̃(x) · τ(x),      τ(x) = E_f[(starts whose fragment covers x) / w],

over the transcript's own bases ``x`` in transcript coordinates, with ``τ`` the fragment-length end taper
(`effective_length.BaseTaper`; ``Σ_x τ(x)`` is the fl-marginal length exactly). With ``c̃`` constant on a
piece — a region of the partition the transcript's exons overlap —

    factor_t = Σ_p ℓ_p^τ · c̃_p / Σ_p ℓ_p^τ,      eff_t = fl_t · factor_t,

``ℓ_p^τ`` the taper-weighted count of the transcript's bases in piece ``p``. No boundary object, no
junction object and no contained support enters the length: a junction-spanning start is counted
through the bases it covers, a piece shorter than a fragment carries its bases' weight, and ``span = fl``
holds for every structure by construction. The region set is every piece the transcript's exons overlap
(a spliced mRNA drops its introns; an unspliced entity's single exon is its span, introns included). The
contraction scales the FL-marginal ``fl_eff_lengths`` the pipeline built, so a field with no reference —
every efficiency 1 — returns them bit-identically; only the captured case contracts.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..types import IntervalType
from .effective_length import base_taper
from .region_arrays import overlapping_region_runs

if TYPE_CHECKING:
    from ..index import TranscriptIndex
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

__all__ = ["transcript_capture_eff_lengths", "transcript_piece_lengths"]


def transcript_piece_lengths(
    index: "TranscriptIndex", region_arrays: "RegionArrays", rna_fl_pmf: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(t, piece, ltau)`` — for every transcript, the pieces its exons (or an unspliced entity's span)
    overlap and the taper-weighted count ``ℓ_p^τ`` of its bases in each, in transcript coordinates.

    Exons are taken in genomic order per transcript; the taper is symmetric in the two ends, so the
    strand does not matter. Annotation and pmf only — sample-independent, so a consumer computes it once
    per index and pmf.
    """
    starts = np.asarray(region_arrays.start, dtype=np.int64)
    ends = np.asarray(region_arrays.end, dtype=np.int64)
    ref_off = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    taper = base_taper(rna_fl_pmf)

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
    ok = (rid >= 0) & (hi > lo)

    # transcript coordinates: cumulative exon lengths per transcript, restarting at each transcript
    exlen = b - a
    n_t = int(tdf.shape[0])
    L_t = np.zeros(n_t, dtype=np.int64)
    np.add.at(L_t, t, exlen)
    off = np.cumsum(exlen) - exlen
    first = np.r_[True, t[1:] != t[:-1]]
    off = off - np.repeat(off[first], np.diff(np.r_[np.flatnonzero(first), t.size]))

    # expand each exon into its pieces
    width = np.where(ok, hi - lo, 0)
    ramp = np.arange(int(width.sum()), dtype=np.int64) - np.repeat(np.cumsum(width) - width, width)
    k = np.repeat(np.arange(t.size), width)
    piece = np.repeat(lo, width) + ramp
    x0 = np.maximum(starts[piece], a[k]) - a[k] + off[k]
    x1 = np.minimum(ends[piece], b[k]) - a[k] + off[k]
    tk = t[k]
    return tk, piece, taper.interval_sums(x0, x1, L_t[tk])


def transcript_capture_eff_lengths(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    index: "TranscriptIndex",
    fl_eff_lengths: np.ndarray,
    rna_fl_pmf: np.ndarray,
) -> np.ndarray:
    """``eff_em_t = fl_t · Σ_p ℓ_p^τ c̃_p / Σ_p ℓ_p^τ`` — each transcript's FL-marginal length times the
    taper-weighted mean of its pieces' capture efficiencies (gate:
    ``tests/calibration/test_capture_eff_length.py``).

    The efficiencies are the result's, ``gdna_capture_efficiency_region``: the posterior mean of each
    piece's clipped gDNA density against the located enriched mode of the fitted landscape, 1 everywhere
    when the field carries no reference (capture off, or no gDNA), so that case returns ``fl`` verbatim.
    One O(incidence) pass does every transcript at once; ``fl · factor ≤ fl`` since every efficiency is
    at most 1, and a global reference shared by every transcript keeps ``eff(unspliced) ≥ eff(spliced)``
    for a parent and its child by construction.
    """
    fl = np.asarray(fl_eff_lengths, dtype=np.float64)
    c = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    if calibration.gdna_reference_density is None:
        return fl.copy()
    t, piece, ltau = transcript_piece_lengths(index, region_arrays, rna_fl_pmf)
    num = np.zeros(fl.shape[0])
    den = np.zeros(fl.shape[0])
    np.add.at(num, t, ltau * c[piece])
    np.add.at(den, t, ltau)
    factor = np.where(den > 0.0, num / np.maximum(den, 1e-300), 1.0)
    return np.minimum(fl * factor, fl)
