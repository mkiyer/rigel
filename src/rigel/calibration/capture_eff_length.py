"""Capture-aware EM effective lengths — extend the gDNA IPR contraction to the mRNA/nRNA components.

**The bug it fixes.** The hybrid-capture inverse-participation-ratio (IPR) effective-length contraction
was applied to the **gDNA component only** (``priors.assemble_priors``); mRNA and nRNA components used the
capture-blind full FL-marginal length. Under capture the gDNA component is artificially concentrated and
out-competes nRNA (and mature) for the capture-enriched reads — the nascent→gDNA siphon (nRNA effective
length measured 5–12× the gDNA component's on the siphon loci).

**The fix.** Contract *every* transcript's EM effective length by the **same per-region gDNA enrichment**,
over the transcript's own region set. gDNA is source-uniform within a locus, so its per-region density is
the probe-enrichment pattern — carrying none of the ~10⁴-fold expression dynamic range that would poison a
coverage-based readout. The region set is the regions the transcript's exons overlap: **exon regions for
mRNA** (introns are gaps), the **full span for nRNA** (its single full-span exon covers introns too).

Contracted in **FL-marginal units** (scale the existing ``effective_lengths`` by the enrichment IPR over
its regions), so a uniform enrichment reduces *exactly* to the input length ⇒ capture-off is bit-identical;
only the captured case contracts. No new readout — it reuses the calibration's per-region gDNA mass and the
same IPR shape as the gDNA component (the standard 1-pseudocount convention + ``1e-9`` numerical floors used
throughout calibration; not tuned constants). The density-correct node model (effective-support
divisors, averaged-side-density pooled seams, transport-free) is documented in
``docs/calibration/effective_length_redesign_plan.md`` §8 (the *why* is in
``capture_effective_length_design.md``).
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..types import IntervalType

if TYPE_CHECKING:
    from ..index import TranscriptIndex
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

__all__ = ["transcript_capture_eff_lengths"]


def _transcript_node_incidence(
    index: "TranscriptIndex", region_arrays: "RegionArrays"
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
]:
    """Per-transcript **node** membership — the regions, boundaries, AND splice junctions a component crosses.

    Returns ``(inc_t_reg, inc_reg, inc_t_bnd, inc_bnd, inc_t_junc, inc_junc_left, inc_junc_right)``: region
    incidence ``(t, r)``; interior-boundary incidence ``(t, r)`` where boundary ``r`` is the seam between
    regions ``r`` and ``r+1`` (left-region keyed, matching ``priors._gdna_region_node_arrays``); and
    SPLICE-JUNCTION incidence ``(t, r_left, r_right)`` — one per adjacent exon pair of a multi-exon mRNA,
    where ``r_left`` is the previous exon's last region and ``r_right`` this exon's first region. The intron
    between them carries no gDNA (no genomic-adjacent seam), so the junction's crossing mass is IMPUTED by
    ``transcript_capture_eff_lengths`` from the two flanking exon densities — stitching the spliced
    transcript into one contiguous ruler so its ``span_full`` equals its FL-marginal length (else the
    junction-dropped ``span_full`` under-states the mature footprint and the ``fl/span_full`` inflation
    lifts a spliced mRNA's EM eff-length ABOVE its nascent parent's — the physically impossible inversion).

    A component's effective length is the IPR over exactly the nodes it occupies contiguously (the
    transcript-structure gate, from first principles):

    * A **region** is in the set if an exon (mRNA) / the full span (nRNA) overlaps it.
    * A **boundary** is in the set iff the component crosses it *without a splice* — i.e. it lies STRICTLY
      INTERIOR to a single exon / span. For an exon ``[a, b)`` spanning region range ``[lo, hi)`` the interior
      seams are ``r ∈ [lo, hi-1)``: their genomic positions ``end[r] = start[r+1]`` all satisfy ``a < · < b``
      (``lo`` is the first region with ``end > a``; ``hi-1`` the last with ``start < b``). Seams at the exon
      EDGES (splice donor/acceptor, or the transcript's outer ends) sit at ``a`` or ``b`` ⇒ index ``< lo`` or
      ``≥ hi-1`` ⇒ excluded automatically. Introns lie in no exon's range ⇒ their regions and boundaries are
      excluded. So a MULTI-exon mRNA drops its introns + splice-junction boundaries but KEEPS an
      exon-interior boundary that merely marks a signature change (an antisense feature overlapping on the
      other strand — crossed contiguously); a SINGLE-exon mRNA / nRNA keeps every interior node (introns
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
    j_t: list[int] = []  # splice-junction seam: transcript
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
            if hi - 1 > lo:  # interior seams r ∈ [lo, hi-1): boundaries crossed contiguously (no splice)
                b_r.append(np.arange(lo, hi - 1, dtype=np.int64))
                b_t.append(np.full(hi - 1 - lo, int(t), dtype=np.int64))
            return lo, hi
        return None

    iv = pd.read_feather(os.path.join(index.index_dir, "intervals.feather"))
    ex = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]
    # genomic order per transcript so consecutive rows of one transcript are ADJACENT exons — the pairs
    # whose SPLICE JUNCTION must be stitched (the intron between them carries no gDNA ⇒ it is not a
    # genomic-adjacent seam; its crossing mass is imputed downstream from the flanking EXON densities).
    ex = ex.sort_values(["t_index", "start"], kind="stable")
    for t, ref_name, a, b in zip(ex["t_index"], ex["ref"], ex["start"], ex["end"], strict=True):
        t = int(t)
        res = _add(t, ref_name, int(a), int(b))
        seen.add(t)
        if res is not None:
            lo, hi = res
            if t in prev_last:  # exon→exon junction to this transcript's previous exon
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
    return (
        np.concatenate(r_t) if r_t else e,
        np.concatenate(r_r) if r_r else e,
        np.concatenate(b_t) if b_t else e,
        np.concatenate(b_r) if b_r else e,
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
    about probe locations (``docs/calibration/enriched_mode_reference_density.md``). Mass-weighting is the
    key: a small captured panel is a tiny COUNT bump but the dominant MASS peak (enriched nodes carry ~100×
    the mass), so its enriched mode is detectable. Unimodal (capture-off / no enrichment) ⇒ the single mode
    ⇒ every node lands at ``w = 1`` ⇒ no contraction. The result is SNAPPED to a real node density so a
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
    # snap to the nearest ACTUAL node density: exact ρ under a uniform field (⇒ factor 1), a real density
    # under capture — no grid-quantization contraction is fabricated.
    return float(rho[ok][int(np.argmin(np.abs(x - mode)))])


def _pooled_seam_arrays(calibration, region_arrays):
    """The left-keyed POOLED-SEAM node arrays ``(seam_mass, seam_support)`` — THE one seam node model both
    the transcript contraction (``transcript_capture_eff_lengths``) and the gDNA component
    (``priors._gdna_region_node_arrays``) share, so the two contract on an identical node basis.

    Seam ``r`` is the boundary between region ``r`` and ``r+1`` (genomically adjacent, same reference):
    ``mass = mass_gdna_right[r] + mass_gdna_left[r+1]`` (the two halves POOLED) and ``support =
    ½·(gdna_boundary_len[r] + gdna_boundary_len[r+1])`` (the AVERAGED per-side density lengths, the
    deposition-faithful divisor). Zero at terminal / cross-reference boundaries. The gDNA path re-keys some
    seams to their right flank (intergenic outer boundaries); the transcript path takes them as-is."""
    right = np.asarray(calibration.mass_gdna_right, dtype=np.float64)
    left = np.asarray(calibration.mass_gdna_left, dtype=np.float64)
    side_len = np.asarray(calibration.gdna_boundary_len, dtype=np.float64)
    ref_id = np.asarray(region_arrays.ref_id)
    n = right.shape[0]
    seam_m = np.zeros(n, dtype=np.float64)  # seam r = boundary between region r and r+1 (left-keyed)
    seam_S = np.zeros(n, dtype=np.float64)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]  # internal seam: genomically adjacent, same reference
        seam_m[:-1] = np.where(same, right[:-1] + left[1:], 0.0)
        seam_S[:-1] = np.where(same, 0.5 * (side_len[:-1] + side_len[1:]), 0.0)
    return seam_m, seam_S


def transcript_capture_eff_lengths(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    index: "TranscriptIndex",
    fl_eff_lengths: np.ndarray,
) -> np.ndarray:
    """Capture-contract each transcript's EM effective length by the per-node gDNA-enrichment density,
    against a single GLOBAL reference density ``ρ_ref`` (the fully-captured level; ``_global_reference_density``).

    ``eff_em_t = fl_t · factor_t``, ``factor_t = [Σ_{n∈t} S_n·min(ρ_n/ρ_ref, 1)] / [Σ_{n∈t} S_n]`` — the
    enrichment-weighted fraction of the transcript's footprint that survives at the reference density, over
    exactly the nodes it occupies CONTIGUOUSLY (``_transcript_node_incidence``), differing ONLY in the node
    set:

    * a per-region CONTAINED node at effective support ``S_r = E[max(0, L_r − ℓ)]`` (mass ``m_r``);
    * a per-interior-boundary POOLED SEAM node at averaged per-side support ``S_s = ½·(E[min(ℓ,L_r)] +
      E[min(ℓ,L_{r+1})])`` (mass ``m_s = right[r] + left[r+1]``) — for boundaries the transcript
      crosses without a splice (interior to an exon);
    * a per-SPLICE-JUNCTION seam node (multi-exon mRNA), same crossing support ``S_j`` but with its mass
      IMPUTED from the two flanking exon densities ``m_j = ½·(ρ_left + ρ_right)·S_j`` — the intron between
      the exons holds no gDNA, so the junction's enrichment is that of the exonic sequence a spliced
      fragment covers, not zero (dropping it) nor full length (the FL-marginal's implicit weight).

    gDNA (contiguous genomic interval) takes ALL nodes; nRNA (single unspliced exon) keeps every interior
    node (introns included); a spliced mRNA takes its exon regions + its interior + splice-junction seams
    (dropping only the introns). Keeping the junction seams makes a spliced mRNA's ``span_full`` equal its
    FL-marginal length — WITHOUT them the ``fl/span_full`` ratio exceeds 1 (growing with exon count) and
    inflates a spliced mRNA's eff-length ABOVE its nascent parent's, the physically impossible inversion
    (a nascent's genomic node set strictly contains its mature child's).

    A single O(incidence) pass (``np.add.at``) does every transcript at once. Properties:

    * **uniform gDNA** (capture off, unimodal density) ⇒ ρ_ref = ρ ⇒ every node ``min(ρ/ρ_ref, 1) = 1``
      ⇒ factor 1 ⇒ ``eff_em = fl`` — bit-identical to the FL-marginal length. NOTE this holds for a
      *noise-free* uniform field; on Poisson-noisy capture-off data the mass-weighted mode can sit slightly
      above the median and manufacture a small spurious contraction. The deferred unimodal / gDNA-abundance
      guard (require a genuine bimodal separation, or cap ρ_ref at a high mass-weighted-density quantile for
      outlier resistance) would neutralise both this and a single high-mass region dominating ρ_ref
      genome-wide; see docs/calibration/enriched_mode_reference_density.md;
    * **no detectable gDNA** ⇒ ``ρ_ref = None`` ⇒ factor 1 (no contraction);
    * **concentrated gDNA** (capture) ⇒ depleted nodes have ρ_n ≪ ρ_ref ⇒ ``min(m_n/ρ_ref, S_n) ≪ S_n``
      ⇒ the eff-len contracts to the enriched footprint.

    A single GLOBAL ``ρ_ref`` (shared by every transcript) makes ``eff(nascent) ≥ eff(mature)`` hold by
    construction — no inversion — and is stable because gDNA barely varies across loci (a few-fold even in
    cancer), unlike RNA. This REPLACES the former per-transcript ``ρ* = G_c/E_c``, which contracted on
    within-transcript density variation including noise — it fired even with NO gDNA and drove the nascent
    siphon. Full derivation + 16-scenario validation: ``docs/calibration/enriched_mode_reference_density.md``
    and ``efflen_shared_reference_fix_plan.md``.
    """
    fl = np.asarray(fl_eff_lengths, dtype=np.float64)
    n_t = fl.shape[0]

    # per-region CONTAINED node (mass, effective support) and per-interior-seam POOLED node. The seam
    # between region r and r+1 (left-region keyed) pools both boundary halves at the averaged per-side
    # density support — the SAME node the gDNA component uses (priors._gdna_region_node_arrays).
    contained_m = np.asarray(calibration.mass_gdna_contained, dtype=np.float64)
    contained_S = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9)
    contained_ev = contained_m + np.asarray(calibration.mass_rna_contained, dtype=np.float64)
    side_len = np.asarray(calibration.gdna_boundary_len, dtype=np.float64)  # for the junction-seam support
    seam_m, seam_S = _pooled_seam_arrays(calibration, region_arrays)  # the SHARED seam node model

    rt, rr, bt, br, jt, jl, jr = _transcript_node_incidence(index, region_arrays)
    # GLOBAL reference density ρ_ref = the enriched mode of the MASS-WEIGHTED node-density KDE — the
    # fully-captured gDNA level detected from the data (no probe assumptions), SHARED across all transcripts
    # so eff(nascent) ≥ eff(mature) by construction. Unimodal (capture-off / no enrichment) ⇒ single mode ⇒
    # every node w=1 ⇒ no contraction. Replaces the per-transcript ρ*=G_c/E_c, which contracted on
    # within-transcript density variation incl. noise — it fired even with NO gDNA, driving the nascent
    # siphon. See docs/calibration/enriched_mode_reference_density.md.
    rho_ref = _global_reference_density(contained_m, contained_S)
    if rho_ref is None or rho_ref <= 0.0:
        return fl.copy()  # no detectable gDNA reference ⇒ no contraction
    inv = 1.0 / rho_ref
    # Per-transcript enrichment-weighted length num = Σ_n min(m_n/ρ_ref, S_n) = Σ_n S_n·min(ρ_n/ρ_ref, 1),
    # the uniform-case length span_full = Σ_n S_n, and the contained evidence (multimapper shrinkage), over
    # the node set (regions + interior seams + splice-junction seams). factor = num/span_full ∈ (0, 1].
    num = np.zeros(n_t)
    span_full = np.zeros(n_t)
    c_ev = np.zeros(n_t)
    if rt.size:
        np.add.at(num, rt, np.minimum(contained_m[rr] * inv, contained_S[rr]))
        np.add.at(span_full, rt, contained_S[rr])
        np.add.at(c_ev, rt, contained_ev[rr])
    if bt.size:
        np.add.at(num, bt, np.minimum(seam_m[br] * inv, seam_S[br]))
        np.add.at(span_full, bt, seam_S[br])
    if jt.size:
        # SPLICE-JUNCTION seams (multi-exon mRNA). The intron between the two exons carries no gDNA, so
        # the junction crossing is NOT a genomic-adjacent seam — its mass is IMPUTED from the two flanking
        # EXON densities ρ = m/S (the exonic sequence a junction-spanning fragment actually covers), at the
        # SAME pooled crossing support S_j = ½·(gdna_boundary_len[left] + gdna_boundary_len[right]) the
        # genomic seams use. Stitching these in makes span_full == fl for a spliced mRNA, so the
        # junction-dropped fl/span_full inflation (which lifted eff_em(mature) ABOVE eff_em(nascent) — the
        # physically impossible inversion) vanishes. Under uniform gDNA m_j = ρ·S_j like every other node,
        # so factor stays EXACTLY 1 (capture-off bit-identical); under capture the junction contributes at
        # its flanking-exon enrichment, not the fabricated full-length weight.
        rho_l = contained_m[jl] / contained_S[jl]
        rho_r = contained_m[jr] / contained_S[jr]
        s_j = 0.5 * (side_len[jl] + side_len[jr])
        m_j = 0.5 * (rho_l + rho_r) * s_j
        np.add.at(num, jt, np.minimum(m_j * inv, s_j))
        np.add.at(span_full, jt, s_j)

    with np.errstate(divide="ignore", invalid="ignore"):
        # factor = Σ min(m_n/ρ_ref, S_n) / Σ S_n ∈ (0, 1] (num ≤ span_full since min(·, S_n) ≤ S_n). Under
        # uniform gDNA every node sits at ρ_ref ⇒ num = span_full ⇒ factor 1 (capture-off bit-identical);
        # under capture depleted nodes contribute min(m_n/ρ_ref, S_n) ≪ S_n ⇒ contracts to the enriched
        # footprint. ONE global ρ_ref for every transcript ⇒ eff(nascent) ≥ eff(mature) (no inversion).
        factor = np.where(span_full > 1e-9, num / np.maximum(span_full, 1e-9), 1.0)
        # multimapper-blindness shrinkage: shrink the contraction toward 1 (no contraction) on sparse
        # CONTAINED evidence (the accumulator is unique-mapper-fed), smoothly (w = C/(C+1), magic-free).
        w = c_ev / (c_ev + 1.0)
        factor = w * factor + (1.0 - w)
    return np.minimum(fl * factor, fl)
