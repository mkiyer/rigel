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

from ..types import STRAND_POS, IntervalType

if TYPE_CHECKING:
    from ..index import TranscriptIndex
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

__all__ = ["transcript_capture_eff_lengths", "build_transcript_warm_start"]


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


def _kde_significant_mode(
    mass: np.ndarray, support: np.ndarray, *, rightmost: bool
) -> "float | None":
    """The rightmost (``rightmost=True``) or leftmost significant peak of the **mass-weighted** log-density
    KDE over the per-region gDNA densities ``ρ = mass/support``, snapped to a real node density.

    Mass-weighting is the key: a small captured panel is a tiny COUNT bump but the dominant MASS peak
    (enriched nodes carry ~100× the mass), so its enriched mode is detectable
    (``docs/calibration/enriched_mode_reference_density.md``). The snap makes a uniform field return its
    density EXACTLY (factor 1, capture-off bit-identical). Returns ``None`` if there is too little gDNA to
    detect a mode (``< 5`` valid nodes). The sole selection axis is which END of the significant-peak set to
    pick — shared verbatim by the enriched reference (rightmost) and the depleted floor (leftmost) so an
    honestly unimodal field yields the SAME snapped density from both."""
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
    else:  # significant peaks: height ≥ _KDE_PROM of the tallest (a real mode, not a tail wiggle)
        h = km[pk]
        sig = pk[h >= _KDE_PROM * float(h.max())]
        cand = sig if sig.size else pk
        mode = grid[int(cand[-1] if rightmost else cand[0])]
    # snap to the nearest ACTUAL node density: exact ρ under a uniform field (⇒ factor 1), a real density
    # under capture — no grid-quantization contraction is fabricated.
    return float(rho[ok][int(np.argmin(np.abs(x - mode)))])


def _global_reference_density(mass: np.ndarray, support: np.ndarray) -> "float | None":
    """The global enriched-mode gDNA reference density for the eff-length contraction — the RIGHTMOST
    significant peak of the mass-weighted log-density KDE (:func:`_kde_significant_mode`), i.e. the
    fully-captured gDNA level detected from the data with no assumption about probe locations. Unimodal
    (capture-off / no enrichment) ⇒ the single mode ⇒ every node lands at ``w = 1`` ⇒ no contraction. ``None``
    when there is too little gDNA to detect a reference (⇒ no contraction)."""
    return _kde_significant_mode(mass, support, rightmost=True)


def _global_depleted_density(mass: np.ndarray, support: np.ndarray) -> "float | None":
    """The global depleted (off-target) gDNA density — the LEFTMOST significant peak of the SAME
    mass-weighted log-density KDE as :func:`_global_reference_density`.

    Feeds the warm-start capture floor ``ε_floor = rho_depleted / rho_ref``: the residual fraction of the
    fully-captured level that off-target (unenriched) sequence still carries. When the field is unimodal the
    leftmost and rightmost significant peaks coincide ⇒ ``rho_depleted == rho_ref`` ⇒ ``ε_floor → 1`` (no
    floor). Returns ``None`` under the same ``< 5`` valid-node guard (⇒ no floor)."""
    return _kde_significant_mode(mass, support, rightmost=False)


def _pooled_seam_mass(right: np.ndarray, left: np.ndarray, ref_id: np.ndarray) -> np.ndarray:
    """Left-keyed POOLED seam mass: seam ``r`` (the boundary between genomically adjacent same-reference
    regions ``r`` and ``r+1``) pools both halves ``right[r] + left[r+1]``; ``0`` at terminal / cross-reference
    boundaries. The mass half of :func:`_pooled_seam_arrays`, shared with the RNA warm-start seam so the
    gDNA and RNA crossing densities sit on ONE node basis."""
    n = right.shape[0]
    out = np.zeros(n, dtype=np.float64)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]  # internal seam: genomically adjacent, same reference
        out[:-1] = np.where(same, right[:-1] + left[1:], 0.0)
    return out


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
    seam_m = _pooled_seam_mass(right, left, ref_id)  # left-keyed pooled seam mass (the shared node model)
    seam_S = np.zeros(n, dtype=np.float64)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]  # internal seam: genomically adjacent, same reference
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


def _capture_correction(rho_ref: "float | None", rho_depleted: "float | None"):
    """Return ``correct(rho_rna, rho_gdna, ev) → capture-corrected density`` (all array-valued).

    Under capture, off-probe RNA is depleted; the local gDNA enrichment ``ε = rho_gdna / rho_ref`` (vs the
    enriched mode) is an endogenous readout of that node's capture efficiency, so ``rho_rna / ε`` lifts the
    reading back to its true level. Two guards against over-correction (dividing by a tiny ε), both
    magic-free: ``ε`` is floored at the DEPLETED mode ``rho_depleted / rho_ref`` (the off-target background —
    unimodal ⇒ floor ``= 1`` ⇒ no correction), and the lift is shrunk toward none by the per-node evidence
    ``w = ev/(ev+1)`` (the same 1-pseudocount form used across calibration). With no detectable reference the
    correction is the identity — the warm start degrades to the raw stranded-density bottleneck."""
    if rho_ref is None or rho_ref <= 0.0:
        return lambda rho_rna, rho_gdna, ev: rho_rna
    inv_ref = 1.0 / rho_ref
    eps_floor = rho_depleted / rho_ref if (rho_depleted is not None and rho_depleted > 0.0) else 1.0

    def correct(rho_rna, rho_gdna, ev):
        eps = np.clip(rho_gdna * inv_ref, eps_floor, 1.0)
        w = ev / (ev + 1.0)
        return rho_rna * (1.0 + w * (1.0 / eps - 1.0))

    return correct


def build_transcript_warm_start(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    index: "TranscriptIndex",
    effective_lengths_em: "np.ndarray | None",
) -> np.ndarray:
    """A calibration-informed per-transcript EM warm-start seed (float64[n_transcripts], counts).

    A transcript's density can be no higher than its scarcest node, so its abundance CEILING is the
    ``min`` over its nodes of the capture-corrected stranded RNA density (``rho ≈ true abundance θ``); the
    seed is that ceiling × the transcript's capture-contracted EM effective length (``θ·E`` = the expected
    OBSERVED count the EM accumulates). The ``min`` is the nascent gate: a nascent shadow bottlenecks on its
    intron / exon↔intron-crossing nodes (≈0 without real nascent), while its shared exons — high, but really
    the mature parent's — never bind. Three node ROLES (``_transcript_node_incidence``), each reading the
    per-strand density :class:`~rigel.calibration.result.RnaWarmStart` exposes on the transcript's OWN strand:

    * **contained** region (exon or intron) — ``rho_contained_s``;
    * **exon↔intron boundary crossing** (nascent) — ``rho_crossing_s`` over the pooled seam;
    * **splice junction** (mature) — ``rho_spliced`` on the junction's fixed motif strand, at the donor
      (region ``r_left``'s right side) and acceptor (region ``r_right``'s left side). This is the mature
      signal that discriminates exon-sharing isoforms; a junction with no spliced reads on the transcript's
      strand contributes ``0`` and correctly bottlenecks that isoform down.

    Observability (``_transcript_node_incidence`` gives structure; the data decides): a node enters the
    ``min`` only if its evidence mass (gDNA+RNA) ``> 0`` — so a genuinely unobserved node (a micro-region, a
    mappability hole) is EXCLUDED, never punished, while an observed-but-RNA-empty node (a covered intron
    with no nascent) is a real, binding zero. A transcript with NO observable node returns ``NaN`` — the
    signal ``"fall back to the coverage seed"`` (never a spurious zero).

    Requires the Phase-A :attr:`CalibrationResult.rna_warm_start` and the capture-contracted
    ``effective_lengths_em``; returns all-``NaN`` (full fallback) if either is absent.
    """
    ws = calibration.rna_warm_start
    if effective_lengths_em is None or ws is None:
        return np.full(int(index.num_transcripts), np.nan)
    eff = np.asarray(effective_lengths_em, dtype=np.float64)
    n_t = eff.shape[0]
    strand = np.asarray(index.t_to_strand_arr)

    # gDNA density (the capture-correction ε) + observed evidence (the observability gate + shrinkage), per
    # node role, on the SAME node model the gDNA component uses (contained region + pooled seam).
    gdna_S = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9)
    gdna_m = np.asarray(calibration.mass_gdna_contained, dtype=np.float64)
    contained_gdna = gdna_m / gdna_S
    contained_ev = gdna_m + np.asarray(calibration.mass_rna_contained, dtype=np.float64)
    ref_id = np.asarray(region_arrays.ref_id)
    seam_gdna_m, seam_S = _pooled_seam_arrays(calibration, region_arrays)
    seam_rna_m = _pooled_seam_mass(
        np.asarray(calibration.mass_rna_right, dtype=np.float64),
        np.asarray(calibration.mass_rna_left, dtype=np.float64),
        ref_id,
    )
    seam_gdna = seam_gdna_m / np.maximum(seam_S, 1e-9)
    seam_ev = seam_gdna_m + seam_rna_m

    correct = _capture_correction(
        _global_reference_density(gdna_m, gdna_S), _global_depleted_density(gdna_m, gdna_S)
    )
    ceiling = np.full(n_t, np.inf)

    def _pick(rho_pos, rho_neg, t_idx, r_idx):
        """The RNA density on each transcript's OWN strand at its node."""
        return np.where(strand[t_idx] == STRAND_POS, rho_pos[r_idx], rho_neg[r_idx])

    def _bottleneck(t_idx, rho_rna, rho_gdna, ev):
        """Fold observable nodes' corrected densities into the running per-transcript ``min``."""
        obs = ev > 0.0
        if np.any(obs):
            np.minimum.at(ceiling, t_idx[obs], correct(rho_rna[obs], rho_gdna[obs], ev[obs]))

    rt, rr, bt, br, jt, jl, jr = _transcript_node_incidence(index, region_arrays)
    if rt.size:  # CONTAINED region nodes
        _bottleneck(rt, _pick(ws.rho_contained_pos, ws.rho_contained_neg, rt, rr),
                    contained_gdna[rr], contained_ev[rr])
    if bt.size:  # exon↔intron BOUNDARY CROSSING (nascent), pooled seam
        _bottleneck(bt, _pick(ws.rho_crossing_pos, ws.rho_crossing_neg, bt, br),
                    seam_gdna[br], seam_ev[br])
    if jt.size:  # SPLICE JUNCTION (mature): donor at r_left's right, acceptor at r_right's left; motif-gated
        donor = np.where(ws.spliced_strand_right[jl] == strand[jt], ws.rho_spliced_right[jl], 0.0)
        accept = np.where(ws.spliced_strand_left[jr] == strand[jt], ws.rho_spliced_left[jr], 0.0)
        _bottleneck(jt, donor, contained_gdna[jl], contained_ev[jl])
        _bottleneck(jt, accept, contained_gdna[jr], contained_ev[jr])

    warm = np.full(n_t, np.nan)  # NaN ⇒ no observable node ⇒ fall back to the coverage seed
    seen = np.isfinite(ceiling)
    warm[seen] = ceiling[seen] * eff[seen]
    return warm
