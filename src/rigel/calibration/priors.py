"""assemble_priors — bridge from CalibrationResult to the per-locus EM prior (PR 6).

Turns the calibration's per-region deconvolved mass + geometric length into the **two
per-locus Dirichlet scalars** the locus EM consumes — ``rna_prior_count`` and
``gdna_prior_count`` — plus the per-locus gDNA-component effective length (the IPR).

The prior's only job is to split each locus's unspliced fragments between gDNA
and RNA; it does **not** attribute RNA mass to individual transcripts (that is
what the EM is for). See ``docs/acc_caljointmodel/prs/PR06_integrate.md`` §I and
``docs/caljointmodel/04_interface_contract.md`` §6.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..locus import MultiLocus
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

# Numerical floor for the gDNA-component effective length: matches the EM's own
# default (``run_batch_locus_em_partitioned`` floors at 1.0), avoiding a zero
# denominator when the EM normalises the gDNA component's abundance.
_GDNA_EFF_LEN_FLOOR = 1.0


@dataclass(frozen=True, slots=True)
class LocusPriors:
    """Per-locus EM prior scalars (float64[n_loci], indexed by ``multi_locus_id``)."""

    gdna_prior_count: np.ndarray  # gDNA-component Dirichlet pseudocount
    rna_prior_count: np.ndarray  # RNA-group Dirichlet pseudocount (the EM splits it by evidence)
    gdna_eff_len: np.ndarray  # IPR-contracted effective length of the gDNA component


def _project_regions_to_loci(
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
    n_loci: int,
    arrays: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Overlap-weighted projection of per-region arrays to per-locus sums.

    For each region, computes its fractional overlap with each ``MultiLocus``
    block, normalises the shares across the loci it touches, and distributes
    each region array's value by that share. Regions overlapping no locus
    (intergenic) are dropped. Adapted from the pre-burn
    ``adaptive_prior._project_to_loci`` (the overlap math only).
    """
    out = {name: np.zeros(n_loci, dtype=np.float64) for name in arrays}
    if n_loci == 0:
        return out

    # Group locus blocks by reference, sorted ascending by start.
    blocks_by_ref: dict[int, list[tuple[int, int, int]]] = {}
    for ml in multi_loci:
        lid = int(ml.multi_locus_id)
        for blk in ml.loci:
            if blk.end > blk.start:
                blocks_by_ref.setdefault(int(blk.ref_id), []).append(
                    (int(blk.start), int(blk.end), lid)
                )
    for blocks in blocks_by_ref.values():
        blocks.sort()

    starts = region_arrays.start
    ends = region_arrays.end
    ref_offsets = region_arrays.ref_offsets
    for ref_id in range(int(region_arrays.n_refs)):
        blocks = blocks_by_ref.get(ref_id)
        if not blocks:
            continue
        lo, hi = int(ref_offsets[ref_id]), int(ref_offsets[ref_id + 1])
        block_starts = np.fromiter((b[0] for b in blocks), dtype=np.int64, count=len(blocks))
        for r in range(lo, hi):
            r_start = int(starts[r])
            r_end = int(ends[r])
            r_len = r_end - r_start
            if r_len <= 0:
                continue
            cand_hi = int(np.searchsorted(block_starts, r_end, side="left"))
            raw: dict[int, float] = {}
            for b_start, b_end, lid in blocks[:cand_hi]:
                if b_end <= r_start:
                    continue
                overlap = min(b_end, r_end) - max(b_start, r_start)
                if overlap > 0:
                    raw[lid] = raw.get(lid, 0.0) + overlap / r_len
            total = sum(raw.values())
            if total <= 0.0:
                continue
            for lid, raw_share in raw.items():
                share = raw_share / total
                for name, arr in arrays.items():
                    out[name][lid] += share * float(arr[r])
    return out


def _transport_boundary_flux(
    contained: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    length: np.ndarray,
    boundary_capacity: np.ndarray,
    ref_id: np.ndarray,
    max_iter: int = 8,
) -> np.ndarray:
    """Length-bias-free boundary-flux transport of per-region gDNA mass.

    Each internal boundary's pooled gDNA mass is re-attributed to its two sides ∝ ``density_ratio·𝓔``
    — the local gDNA density ratio (``g/L`` over the library average, length-bias-free) × the
    directional boundary effective length ``𝓔(L)=E[min(ℓ,L)]`` (``boundary_capacity``). This moves
    capture smear off the depleted (e.g. intronic) side and onto the probed side that generated it.
    Iterates until the
    total mass moved is sub-count (< 1 fragment-equivalent), capped at ``max_iter``; total
    mass is conserved (ref-edge / cross-ref sides stay in place). See
    docs/futureprs/phase6_boundary_flux_transport_plan.md.
    """
    contained = np.asarray(contained, dtype=np.float64)
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    length = np.maximum(np.asarray(length, dtype=np.float64), 1e-9)
    boundary_capacity = np.asarray(boundary_capacity, dtype=np.float64)
    r = contained.shape[0]
    if r <= 1:
        return contained + left + right
    same = np.asarray(ref_id)[:-1] == np.asarray(ref_id)[1:]  # boundary (i,i+1) internal
    total_len = float(length.sum())
    g = contained + left + right
    prev = g
    for _ in range(max_iter):
        gdna_density_global = g.sum() / total_len if total_len > 0.0 else 0.0
        density_ratio = np.where(
            gdna_density_global > 0.0, (g / length) / max(gdna_density_global, 1e-12), 1.0
        )
        w = density_ratio * boundary_capacity
        pooled = right[:-1] + left[1:]  # boundary (i,i+1) pooled gDNA mass
        denom = w[:-1] + w[1:]
        ok = same & (denom > 0.0)
        share_l = np.where(ok, pooled * w[:-1] / np.maximum(denom, 1e-12), 0.0)
        share_r = np.where(ok, pooled * w[1:] / np.maximum(denom, 1e-12), 0.0)
        out = contained.copy()
        out[:-1] += share_l  # → left region of the boundary
        out[1:] += share_r  # → right region of the boundary
        keep_left = np.ones(r, dtype=bool)
        keep_left[1:] = ~ok  # keep sides whose boundary was not transported
        keep_right = np.ones(r, dtype=bool)
        keep_right[:-1] = ~ok
        out = out + np.where(keep_left, left, 0.0) + np.where(keep_right, right, 0.0)
        moved = float(np.abs(out - prev).sum())
        prev = out
        g = out
        if moved < 1.0:  # sub-count: no further observable change
            break
    return g


def assemble_priors(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
) -> LocusPriors:
    """Build the per-locus EM prior from the acyclic calibration result.

    The gDNA node set (PR-A; ``effective_length_redesign_plan.md`` §7) is the region-contained nodes PLUS
    one **pooled seam node** per internal boundary — *not* the ``contained+left+right`` fold (PR-1):

        region node r:    mass = mass_gdna_contained[r],         length = region_size_bp[r]  (genomic)
        seam node (r,r+1): mass = mass_gdna_right[r]+mass_gdna_left[r+1] (the boundary's two halves POOLED),
                           length = min(gdna_boundary_len[r], gdna_boundary_len[r+1]) ≈ fl_mean  (FL-scale)

    keyed to the left-flank region r. Both conservation laws hold: **mass** — ``Σ contained + Σ seam`` =
    the total deconvolved gDNA (each non-terminal boundary side counted once; terminals are zero);
    **genomic length** — ``Σ region_size_bp`` = the covered genomic span (regions tile; seams add 0
    genomic bp — their FL length is a density normalizer inside the dimensionless IPR only). The region +
    seam masses/participations project to loci by genomic-overlap ``share``::

        gdna_prior_count = Σ share * (contained + seam)                          (deconvolved gDNA count)
        rna_prior_count  = Σ share * rna_region            (UNSPLICED RNA; spliced withheld — see below)
        gdna_eff_len     = (G+1)² / [ Σ share*(contained²/L_g + seam²/L_e) + (2G+1)/span ],  capped at span
                           G = Σ share*(contained+seam),  L_g = region_size_bp,  span = Σ share*region_size_bp

    **Why pool, not split:** entering a seam's two halves ``s_L,s_R`` as separate IPR nodes contributes
    support ``s_L²/L_e + s_R²/L_e ≈ s²/(2 L_e)``; pooling contributes ``(s_L+s_R)²/L_e = s²/L_e`` — by
    Cauchy-Schwarz the split halves the IPR support (halves the contraction), nearly doubling the leak
    (measured 13.5% vs pooled 5.5%). Pooling quarantines the captured intron↔exon crossing mass at FL
    density (``s²/fl_mean``), recovering the gDNA concentration that PR-1's region-dilution (``s²/intron_bp``)
    lost — **without** the old transport (which moved mass) or the FL-aware ``gdna_geom_len`` (which summed
    to ~8% > the genomic span, ≈2× on small exons, violating length conservation). ``gdna_eff_len`` is the
    **inverse participation ratio** of the per-node gDNA density (its reciprocal makes the gDNA component's
    per-position rate the local density, so under capture the support contracts to the probed footprint and
    gDNA competes at its true density), correct by construction — uniform gDNA ⇒ ``eff_len = span``.
    **Laplace-smoothed** by one fragment-equivalent of uniform per-base prior mass (the ``(2G+1)/span``
    term), the canonical add-one prior with no tunable constant. Measured flagship leak:
    transport 7.97% → PR-1 8.98% → pooled-seam **5.52%** (below the old transport baseline).

    The add-one is the canonical Laplace prior — there is **no tunable shrinkage constant**. It
    shrinks the IPR toward the uniform geometric ``span`` exactly in proportion to the gDNA
    evidence: a tiny sparse-but-concentrated mass stays ≈ ``span`` (the uniform prior dominates),
    so the EM cannot amplify it past the calibration's call; abundant gDNA (capture, G ≫ 1)
    recovers the empirical concentration; a locus with no deconvolved gDNA gives ``span`` exactly
    (G=0, the multimap-blind null — no special case). See
    docs/futureprs/phase6_boundary_flux_transport_plan.md.
    """
    if calibration.n_regions != region_arrays.n_regions:
        raise ValueError(
            f"calibration has {calibration.n_regions} regions but region_arrays has "
            f"{region_arrays.n_regions}; they must address the same partition."
        )

    # CONSERVATION-CORRECT POOLED-SEAM node set (PR-A; effective_length_redesign_plan.md §7). The gDNA
    # node set is the region-contained nodes PLUS one **pooled seam node** per internal boundary:
    #   • region node r: mass = mass_gdna_contained[r],  length = region_size_bp[r] (genomic bp);
    #   • seam node (r,r+1): mass = mass_gdna_right[r] + mass_gdna_left[r+1] (the boundary's BOTH halves
    #     POOLED — gDNA is two-sided), length = min(gdna_boundary_len[r], gdna_boundary_len[r+1]) ≈ fl_mean
    #     (the FL-scale density normalizer; NOT genomic — seams contribute 0 to the genomic span). The seam
    #     is keyed to its left-flank region r.
    # Conservation: mass — Σ contained + Σ seam = Σ contained + Σ both sides (each non-terminal side once;
    # terminals are zero) = the total deconvolved gDNA; length — Σ region_size_bp = the covered genomic span
    # (regions tile; seams add 0). The seam mass is POOLED, not split into two side nodes: splitting a
    # seam's mass s into s_L,s_R as separate IPR nodes contributes support s_L²/Le + s_R²/Le ≈ s²/(2Le),
    # whereas pooling contributes (s_L+s_R)²/Le = s²/Le — by Cauchy-Schwarz the split HALVES the IPR
    # support (halves the contraction), nearly doubling the leak (measured 13.5% vs pooled 5.5%). Pooling
    # quarantines the captured intron↔exon crossing mass at FL density (s²/fl_mean), recovering the
    # concentration PR-1's region-dilution lost (s²/intron_bp) — WITHOUT the transport's mass-moving or the
    # gdna_geom_len length-violation. Measured flagship leak: transport 7.97% → PR-1 8.98% → pooled 5.52%.
    length = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    ref_id = np.asarray(region_arrays.ref_id)
    contained = np.asarray(calibration.mass_gdna_contained, dtype=np.float64)
    side_left = np.asarray(calibration.mass_gdna_left, dtype=np.float64)
    side_right = np.asarray(calibration.mass_gdna_right, dtype=np.float64)
    side_len = np.maximum(np.asarray(calibration.gdna_boundary_len, dtype=np.float64), 1e-9)
    pooled = np.zeros_like(contained)  # seam mass keyed to its left-flank region
    pooled_supp = np.zeros_like(contained)  # seam participation s²/Le keyed to its left-flank region
    if contained.shape[0] > 1:
        same = ref_id[:-1] == ref_id[1:]  # internal seams (genomically adjacent, same reference)
        seam_mass = side_right[:-1] + side_left[1:]  # POOL the boundary's two halves into one node
        seam_len = np.minimum(side_len[:-1], side_len[1:])  # FL-scale; the tighter flank cap
        pooled[:-1] = np.where(same, seam_mass, 0.0)
        pooled_supp[:-1] = np.where(same, seam_mass**2 / seam_len, 0.0)
    gdna_region = contained + pooled  # total gDNA per region = contained + its right-seam node
    # RNA prior = the UNSPLICED RNA mass only. Spliced fragments have no gDNA candidate in the EM
    # (gDNA does not splice) → they are guaranteed-RNA and the EM assigns them directly; counting
    # them in rna_prior_count would double-count them and unfairly inflate the RNA side of the
    # gDNA-vs-RNA unspliced split (the prior arbitrates only the unspliced fragments, so a_g+a_r
    # should equal the unspliced competing mass). mass_rna_* is spliced-inclusive (node
    # conservation); subtracting mass_rna_spliced here yields (1−g)·M_unspliced per region.
    rna_region = np.maximum(
        calibration.mass_rna_contained
        + calibration.mass_rna_left
        + calibration.mass_rna_right
        - calibration.mass_rna_spliced,
        0.0,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        # IPR participation = region term (contained²/genomic_len) + seam term (pooled²/Le). The seam term
        # is pre-computed above (pooled_supp); both ride region r's locus overlap share in the projection.
        reg_supp = np.where(length > 0.0, contained**2 / np.maximum(length, 1e-9), 0.0)
        gdna_sq_over_len = reg_supp + pooled_supp

    proj = _project_regions_to_loci(
        region_arrays,
        multi_loci,
        len(multi_loci),
        {
            "gdna": gdna_region,
            "rna": rna_region,
            "support": gdna_sq_over_len,
            "span": length,
            # the CONTAINED (unique-mapper) mass per locus — the calibration-blindness discriminator
            # for the eff-len guard below (calibration's accumulator is fed by unique mappers only).
            "gdna_contained": np.asarray(calibration.mass_gdna_contained, dtype=np.float64),
            "rna_contained": np.asarray(calibration.mass_rna_contained, dtype=np.float64),
        },
    )
    gdna_locus, support_locus, span = proj["gdna"], proj["support"], proj["span"]
    # gDNA effective length = the inverse participation ratio of the per-region gDNA
    # mass, Laplace-smoothed by one fragment-equivalent of uniform (per-base) prior
    # mass. Add c_r = L_r/span pseudo-mass (total 1, spread per base) to each
    # region's mass g_r and form the IPR (G̃)² / Σ(g̃²/L) with G̃ = G + 1. In the
    # projected quantities this is closed-form, since Σ(2·g·c/L) = 2G/span and
    # Σ(c²/L) = 1/span:
    #     eff_len = (G + 1)² / [ Σ(g²/L) + (2G + 1)/span ],   G = Σ g  (= gdna_locus)
    # The add-one is the canonical Laplace prior — there is NO tunable shrinkage constant. It shrinks the
    # support toward the uniform span in proportion to the gDNA evidence: G = 0 → span; abundant gDNA
    # (capture, G ≫ 1) → the empirical concentration (the IPR); capped at span. NOTE the Laplace term
    # (2G+1)/span is at GENOMIC scale while a pooled SEAM support term is at FL scale (≈ fl_mean ≪ span),
    # so the add-one does NOT by itself protect a locus whose entire gDNA signal is FL-scale seam smear —
    # there the seam term can dominate and over-contract. The contained-evidence gate below is what guards
    # that multimap-blind case (NOT the +1 prior).
    with np.errstate(divide="ignore", invalid="ignore"):
        safe_span = np.maximum(span, 1e-9)
        smoothed_support = support_locus + (2.0 * gdna_locus + 1.0) / safe_span
        eff_len = np.minimum((gdna_locus + 1.0) ** 2 / smoothed_support, span)

    # CONTAINED-EVIDENCE GATE (calibration multimapper-blindness). Calibration's accumulator is fed by
    # UNIQUE mappers only (bam_scanner.cpp: deposit_to_accumulator is inside `if (is_unique_mapper)`; the
    # EM buffer separately gets multimappers at 1/NH). A locus dominated by multimappers — e.g. identical
    # paralogs — therefore has ~zero CONTAINED mass; the only gDNA signal it retains is the sparse,
    # asymmetric unique-flank smear the pooled seam imports from adjacent introns at FL density, which
    # over-contracts gdna_eff_len and lets the gDNA component over-claim the (degenerate) RNA budget. Where
    # calibration is structurally blind (contained gDNA + contained RNA ≈ 0) the prior must be honestly
    # uninformative: revert gdna_eff_len → span (no contraction). This fires ONLY at contained == 0 (NOT a
    # multimap-fraction blend, which would erode the capture contraction on real gDNA-bearing duplicated
    # loci); on every locus calibration can see it is the identity, so the capture leak win is preserved by
    # construction. The real fix (multimapper-aware accumulator) is the mappability initiative; this is the
    # principled interim. See docs/calibration/splice_junction_node_architecture.md siblings / project_mappability.
    blind = (proj["gdna_contained"] + proj["rna_contained"]) < 1e-6
    eff_len = np.where(blind, span, eff_len)

    return LocusPriors(
        gdna_prior_count=gdna_locus,
        rna_prior_count=proj["rna"],
        gdna_eff_len=np.maximum(eff_len, _GDNA_EFF_LEN_FLOOR),
    )


__all__ = ["LocusPriors", "assemble_priors"]
