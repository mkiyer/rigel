"""Capture-aware EM effective lengths — the gDNA enrichment contraction, applied to every component.

Under hybrid capture a transcript's usable length is not its full length: only the probe-enriched part
of its footprint is sampled. Contracting the gDNA component alone leaves it artificially concentrated,
so it out-competes the RNA components for the enriched reads. This module therefore contracts EVERY
transcript's EM effective length by the SAME per-region gDNA enrichment, over the transcript's own
region set.

gDNA is the right readout of that enrichment because it is source-uniform within a locus, so its
per-region density is the probe pattern and carries none of the expression dynamic range that would
poison a coverage-based readout. The region set is the regions the transcript's exons overlap: exon
regions for a spliced mRNA, since introns are gaps, and the full span for an unspliced component,
whose single full-span exon covers introns too.

The contraction is applied in FL-marginal units — it scales the existing ``effective_lengths`` by the
enrichment ratio over the transcript's regions — so a uniform enrichment reduces exactly to the input
length and a capture-off run is bit-identical; only the captured case contracts. It introduces no new
readout, reusing the calibration's per-region gDNA mass and the same divisors the calibrator itself
used.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..types import IntervalType
from .region_arrays import region_right_boundary

if TYPE_CHECKING:
    from ..index import TranscriptIndex
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

__all__ = ["transcript_capture_eff_lengths"]


def _transcript_region_incidence(
    index: "TranscriptIndex", region_arrays: "RegionArrays"
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-transcript membership — the regions, boundaries and splice junctions a component crosses.

    Returns ``(inc_t_reg, inc_reg, inc_t_bnd, inc_bnd, inc_t_junc, inc_junc_left, inc_junc_right)``:
    region incidence ``(t, r)``; interior-boundary incidence ``(t, e)`` where ``e`` is a CONTIGUOUS
    BOUNDARY INDEX, since a boundary is a first-class object on its own axis and a consumer indexes
    the per-boundary arrays directly; and splice-sj incidence ``(t, r_left, r_right)``, one per
    adjacent exon pair of a multi-exon mRNA, where ``r_left`` is the previous exon's last region and
    ``r_right`` this exon's first region. The intron between them carries no gDNA and is no
    genomic-adjacent boundary, so the sj's crossing mass is imputed by
    ``transcript_capture_eff_lengths`` from the two flanking exon densities. That stitches the
    spliced transcript into one contiguous ruler, so its ``span_full`` equals its FL-marginal
    length; dropping the sj instead under-states the mature footprint and the resulting
    ``fl/span_full`` inflation lifts a spliced mRNA's EM effective length above its unspliced
    parent's, which is impossible.

    A component's effective length is the IPR over exactly the regions it occupies contiguously (the
    transcript-structure gate, from first principles):

    * A **region** is in the set if an exon (mRNA) / the full span (nRNA) overlaps it.
    * A **boundary** is in the set iff the component crosses it *without a splice* — i.e. it lies STRICTLY
      INTERIOR to a single exon / span. For an exon ``[a, b)`` spanning region range ``[lo, hi)`` the interior
      boundaries are ``r ∈ [lo, hi-1)``: their genomic positions ``end[r] = start[r+1]`` all satisfy ``a < · < b``
      (``lo`` is the first region with ``end > a``; ``hi-1`` the last with ``start < b``). Boundaries at the exon
      BOUNDARIES (splice donor/acceptor, or the transcript's outer ends) sit at ``a`` or ``b`` ⇒ index ``< lo`` or
      ``≥ hi-1`` ⇒ excluded automatically. Introns lie in no exon's range ⇒ their regions and boundaries are
      excluded. So a MULTI-exon mRNA drops its introns + splice-junction boundaries but KEEPS an
      exon-interior boundary that merely marks a signature change (an antisense feature overlapping on the
      other strand — crossed contiguously); a SINGLE-exon mRNA / nRNA keeps every interior region (introns
      included, for nRNA); the outer boundaries are never interior ⇒ excluded (they belong to gDNA, which
      spans the chromosome). Annotation-only (sample-independent) — could be precomputed at index build.
    """
    starts = np.asarray(region_arrays.start, dtype=np.int64)
    ends = np.asarray(region_arrays.end, dtype=np.int64)
    ref_off = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    name_to_id = index.ref_name_to_id
    r_t: list[np.ndarray] = []
    r_r: list[np.ndarray] = []
    b_t: list[np.ndarray] = []
    b_r: list[np.ndarray] = []
    j_t: list[int] = []  # splice-junction boundary: transcript
    j_l: list[int] = []  # left-flank region (the previous exon's LAST region)
    j_r: list[int] = []  # right-flank region (this exon's FIRST region)
    seen: set[int] = set()
    prev_last: dict[int, int] = {}  # t → last region of its previous (genomically earlier) exon

    def _add(t: int, ref_name: object, a: int, b: int) -> "tuple[int, int] | None":
        rid = name_to_id.get(str(ref_name))
        if rid is None:
            return None
        lo0, hi0 = int(ref_off[rid]), int(ref_off[rid + 1])
        # regions overlapping [a, b): first with end_r > a (contains/after a) through last with start_r < b.
        lo = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        hi = lo0 + int(np.searchsorted(starts[lo0:hi0], b, side="left"))
        if hi > lo:
            r_r.append(np.arange(lo, hi, dtype=np.int64))
            r_t.append(np.full(hi - lo, int(t), dtype=np.int64))
            if (
                hi - 1 > lo
            ):  # interior boundaries r ∈ [lo, hi-1): boundaries crossed contiguously (no splice)
                b_r.append(np.arange(lo, hi - 1, dtype=np.int64))
                b_t.append(np.full(hi - 1 - lo, int(t), dtype=np.int64))
            return lo, hi
        return None

    iv = pd.read_feather(os.path.join(index.index_dir, "intervals.feather"))
    ex = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]
    # genomic order per transcript so consecutive rows of one transcript are ADJACENT exons — the pairs
    # whose SPLICE JUNCTION must be stitched (the intron between them carries no gDNA ⇒ it is not a
    # genomic-adjacent boundary; its crossing mass is imputed downstream from the flanking EXON densities).
    ex = ex.sort_values(["t_index", "start"], kind="stable")
    for t, ref_name, a, b in zip(ex["t_index"], ex["ref"], ex["start"], ex["end"], strict=True):
        t = int(t)
        res = _add(t, ref_name, int(a), int(b))
        seen.add(t)
        if res is not None:
            lo, hi = res
            if t in prev_last:  # exon→exon sj to this transcript's previous exon
                j_t.append(t)
                j_l.append(prev_last[t])
                j_r.append(lo)
            prev_last[t] = hi - 1

    tdf = index.t_df
    if tdf is not None and "is_synthetic" in tdf.columns:
        syn = tdf[tdf["is_synthetic"].to_numpy(dtype=bool)]
        for t, ref_name, a, b in zip(
            syn["t_index"], syn["ref"], syn["start"], syn["end"], strict=True
        ):
            if int(t) in seen:
                continue
            _add(int(t), ref_name, int(a), int(b))  # single-exon spans (nRNA) → no splice junctions

    e = np.empty(0, dtype=np.int64)
    # The boundary axis is emitted as a BOUNDARY index, not a left-region index. A boundary is a
    # first-class object with its own axis, so a consumer indexes the per-boundary arrays directly;
    # returning the left region instead forces every caller through a region-shaped copy, which
    # reads as an attribution of the boundary's mass to a region and is not one.
    right_boundary = region_right_boundary(np.asarray(region_arrays.ref_id))
    b_boundaries = right_boundary[np.concatenate(b_r)] if b_r else e
    return (
        np.concatenate(r_t) if r_t else e,
        np.concatenate(r_r) if r_r else e,
        np.concatenate(b_t) if b_t else e,
        b_boundaries,
        np.asarray(j_t, dtype=np.int64) if j_t else e,
        np.asarray(j_l, dtype=np.int64) if j_l else e,
        np.asarray(j_r, dtype=np.int64) if j_r else e,
    )


# KDE smoothing params (log-density space): a fixed bandwidth + a peak-prominence floor. Standard KDE
# smoothing, not tuned to a target; validated across the 16-scenario suite (kde_mode_scan.py). Could be made
# data-driven (Silverman) later if needed.
_KDE_BW = 0.4
_KDE_PROM = 0.05


def _global_reference_density(mass: np.ndarray, support: np.ndarray) -> "float | None":
    """The global enriched-mode gDNA reference density for the eff-length contraction.

    The rightmost significant peak of the **mass-weighted** log-density KDE over the per-region gDNA
    densities ``ρ = mass/support`` — the fully-captured gDNA level, detected from the data with no assumption
    about probe locations. Mass-weighting is the
    key: a small captured panel is a tiny COUNT bump but the dominant MASS peak (enriched regions carry ~100×
    the mass), so its enriched mode is detectable. Unimodal (capture-off / no enrichment) ⇒ the single mode
    ⇒ every region lands at ``w = 1`` ⇒ no contraction. The result is SNAPPED to a real region density so a
    uniform field returns its density EXACTLY (factor 1, capture-off bit-identical). Returns ``None`` if
    there is too little gDNA to detect a reference (⇒ no contraction)."""
    m = np.asarray(mass, dtype=np.float64)
    s = np.maximum(np.asarray(support, dtype=np.float64), 1e-9)
    rho = m / s
    ok = np.isfinite(rho) & (rho > 1e-12) & (m > 0.0)
    if int(ok.sum()) < 5:
        return None
    x = np.log(rho[ok])
    wt = m[ok]
    grid = np.linspace(float(x.min()) - 1.0, float(x.max()) + 1.0, 512)
    wn = wt / wt.sum()
    d = (grid[:, None] - x[None, :]) / _KDE_BW
    km = (wn[None, :] * np.exp(-0.5 * d * d)).sum(1)  # mass-weighted log-density KDE
    pk = np.where((km[1:-1] >= km[:-2]) & (km[1:-1] > km[2:]))[0] + 1
    if pk.size == 0:
        mode = grid[int(np.argmax(km))]
    elif pk.size == 1:
        mode = grid[int(pk[0])]
    else:  # rightmost peak with height ≥ _KDE_PROM of the tallest (a real mode, not a tail wiggle)
        h = km[pk]
        sig = pk[h >= _KDE_PROM * float(h.max())]
        mode = grid[int((sig if sig.size else pk)[-1])]
    # snap to the nearest ACTUAL region density: exact ρ under a uniform field (⇒ factor 1), a real density
    # under capture — no grid-quantization contraction is fabricated.
    return float(rho[ok][int(np.argmin(np.abs(x - mode)))])


def transcript_capture_eff_lengths(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    index: "TranscriptIndex",
    fl_eff_lengths: np.ndarray,
) -> np.ndarray:
    """Capture-contract each transcript's EM effective length by the per-region gDNA-enrichment density,
    against a single GLOBAL reference density ``ρ_ref`` (the fully-captured level; ``_global_reference_density``).

    ``eff_em_t = fl_t · factor_t``, ``factor_t = [Σ_{n∈t} S_n·min(ρ_n/ρ_ref, 1)] / [Σ_{n∈t} S_n]`` — the
    enrichment-weighted fraction of the transcript's footprint that survives at the reference density, over
    exactly the regions it occupies CONTIGUOUSLY (``_transcript_region_incidence``), differing ONLY in the region
    set:

    * a per-region CONTAINED region at effective support ``S_r = E[max(0, L_r − ℓ)]`` (mass ``m_r``);
    * a per-interior-BOUNDARY crossing object at support ``S_e = gdna_boundary_eff_len[e] = E_f[w-1]``
      (mass ``m_e = mass_gdna_boundary[e]``), for boundaries the transcript crosses without a splice,
      i.e. interior to an exon;
    * a per-SPLICE-SJ crossing object (multi-exon mRNA), same crossing support ``S_j`` but with its
      mass imputed from the two flanking exon densities, ``m_j = 0.5*(rho_left + rho_right)*S_j``.
      The intron between the exons holds no gDNA, so the sj's enrichment is that of the exonic
      sequence a spliced fragment covers — neither zero, which dropping it would imply, nor full
      length, which the FL marginal implies.

    gDNA, a contiguous genomic interval, takes ALL regions; an unspliced component keeps every
    interior region, introns included; a spliced mRNA takes its exon regions plus its interior and
    splice-junction boundaries, dropping only the introns. Keeping the sj boundaries makes a spliced
    mRNA's ``span_full`` equal its FL-marginal length. Without them the ``fl/span_full`` ratio
    exceeds 1, growing with exon count, and inflates a spliced mRNA's effective length above its
    unspliced parent's, which is impossible: the parent's genomic region set strictly contains the
    child's.

    A single O(incidence) pass (``np.add.at``) does every transcript at once. Properties:

    * uniform gDNA (capture off, unimodal density) gives ``rho_ref = rho``, so every region has
      ``min(rho/rho_ref, 1) = 1``, the factor is 1 and ``eff_em == fl``, bit-identical to the
      FL-marginal length. That holds for a noise-free uniform field; on Poisson-noisy capture-off
      data the mass-weighted mode can sit slightly above the median and manufacture a small
      spurious contraction, and there is no unimodality guard against it;
    * no detectable gDNA gives ``rho_ref = None`` and a factor of 1, i.e. no contraction;
    * concentrated gDNA (capture) leaves depleted regions with ``rho_n`` far below ``rho_ref``, so
      ``min(m_n/rho_ref, S_n)`` falls well below ``S_n`` and the effective length contracts to the
      enriched footprint.

    The reference density is a single GLOBAL number shared by every transcript. That makes
    ``eff(unspliced) >= eff(spliced)`` hold by construction, with no inversion, and it is stable
    because gDNA barely varies across loci, unlike RNA. A per-transcript reference instead
    contracts on within-transcript density variation including noise, which fires even with no gDNA
    present.
    """
    fl = np.asarray(fl_eff_lengths, dtype=np.float64)
    n_t = fl.shape[0]

    # per-region CONTAINED object (mass, effective support) and per-interior-BOUNDARY crossing object. The
    # boundary between region r and r+1 is keyed to r — the SAME objects the gDNA component uses
    # (priors._gdna_region_arrays).
    contained_m = np.asarray(calibration.mass_gdna_region, dtype=np.float64)
    contained_S = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9)
    contained_ev = contained_m + np.asarray(calibration.mass_rna_region, dtype=np.float64)
    # The per-BOUNDARY crossing objects, on their own axis. `inc_bnd` is a BOUNDARY index, so these are
    # indexed directly — no region-shaped copy, and nothing that reads as an attribution to a region.
    boundary_m = np.asarray(calibration.mass_gdna_boundary, dtype=np.float64)
    boundary_S = np.maximum(np.asarray(calibration.gdna_boundary_eff_len, dtype=np.float64), 0.0)
    # A SPLICE sj is not a contiguous boundary, so it has no entry on the boundary axis — but it is
    # still a crossing, and gDNA's crossing divisor is the same everywhere (UNBOUNDED_REACH both
    # sides gives mu_g - 1). Take it from the boundary supports, not a length model recomputed here: one
    # definition, and it cannot drift from the one the calibrator divided by.
    crossing_S = float(boundary_S[boundary_S > 0.0][0]) if np.any(boundary_S > 0.0) else 0.0

    rt, rr, bt, br, jt, jl, jr = _transcript_region_incidence(index, region_arrays)
    # GLOBAL reference density ρ_ref = the enriched mode of the MASS-WEIGHTED region-density KDE — the
    # fully-captured gDNA level detected from the data (no probe assumptions), SHARED across all transcripts
    # so eff(nascent) ≥ eff(mature) by construction. Unimodal (capture-off / no enrichment) ⇒ single mode ⇒
    # every region weighs 1 and there is no contraction. A per-transcript reference would instead
    # contract on within-transcript density variation including noise, firing even with no gDNA.
    rho_ref = _global_reference_density(contained_m, contained_S)
    if rho_ref is None or rho_ref <= 0.0:
        return fl.copy()  # no detectable gDNA reference ⇒ no contraction
    inv = 1.0 / rho_ref
    # Per-transcript enrichment-weighted length num = Σ_n min(m_n/ρ_ref, S_n) = Σ_n S_n·min(ρ_n/ρ_ref, 1),
    # the uniform-case length span_full = Σ_n S_n, and the contained evidence (multimapper shrinkage), over
    # the region set (regions + interior boundaries + splice-junction boundaries). factor = num/span_full ∈ (0, 1].
    num = np.zeros(n_t)
    span_full = np.zeros(n_t)
    c_ev = np.zeros(n_t)
    if rt.size:
        np.add.at(num, rt, np.minimum(contained_m[rr] * inv, contained_S[rr]))
        np.add.at(span_full, rt, contained_S[rr])
        np.add.at(c_ev, rt, contained_ev[rr])
    if bt.size:
        np.add.at(num, bt, np.minimum(boundary_m[br] * inv, boundary_S[br]))
        np.add.at(span_full, bt, boundary_S[br])
    if jt.size:
        # SPLICE-SJ boundaries (multi-exon mRNA). The intron between the two exons carries no gDNA, so
        # the sj crossing is NOT a genomic-adjacent boundary — its mass is IMPUTED from the two flanking
        # EXON densities ρ = m/S (the exonic sequence a sj-spanning fragment actually covers), at the
        # same crossing support every genomic boundary uses. Stitching these in makes span_full == fl
        # for a spliced mRNA, so the sj-dropped fl/span_full inflation — which lifts a spliced
        # transcript's eff_em above its unspliced parent's, an impossible inversion — vanishes.
        # Under uniform gDNA m_j = ρ·S_j like every other region,
        # so factor stays EXACTLY 1 (capture-off bit-identical); under capture the sj contributes at
        # its flanking-exon enrichment, not the fabricated full-length weight.
        rho_l = contained_m[jl] / contained_S[jl]
        rho_r = contained_m[jr] / contained_S[jr]
        # The sj boundary's SUPPORT is the gDNA crossing effective length: one number, the same one
        # every contiguous boundary uses, taken from the boundary supports rather than re-derived
        # here, so it cannot drift from the divisor the calibrator applied. The `0.5·(rho_l + rho_r)`
        # below is a genuine AVERAGE OF DENSITIES — the sj's imputed density is the mean of its two
        # flanks — and is unrelated to the support.
        s_j = np.full(jt.shape[0], crossing_S, dtype=np.float64)
        m_j = 0.5 * (rho_l + rho_r) * s_j
        np.add.at(num, jt, np.minimum(m_j * inv, s_j))
        np.add.at(span_full, jt, s_j)

    with np.errstate(divide="ignore", invalid="ignore"):
        # factor = Σ min(m_n/ρ_ref, S_n) / Σ S_n ∈ (0, 1] (num ≤ span_full since min(·, S_n) ≤ S_n). Under
        # uniform gDNA every region sits at ρ_ref ⇒ num = span_full ⇒ factor 1 (capture-off bit-identical);
        # under capture depleted regions contribute min(m_n/ρ_ref, S_n) ≪ S_n ⇒ contracts to the enriched
        # footprint. ONE global ρ_ref for every transcript ⇒ eff(unspliced) ≥ eff(spliced), no inversion.
        factor = np.where(span_full > 1e-9, num / np.maximum(span_full, 1e-9), 1.0)
        # multimapper-blindness shrinkage: shrink the contraction toward 1 (no contraction) on sparse
        # CONTAINED evidence (the accumulator is unique-mapper-fed), smoothly (w = C/(C+1), magic-free).
        w = c_ev / (c_ev + 1.0)
        factor = w * factor + (1.0 - w)
    return np.minimum(fl * factor, fl)
