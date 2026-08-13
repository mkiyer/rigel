"""assemble_priors — bridge from CalibrationResult to the per-locus EM prior (PR 6).

Turns the calibration's per-region deconvolved mass + geometric length into the **two
per-locus Dirichlet scalars** the locus EM consumes — ``rna_prior_count`` and
``gdna_prior_count`` — plus the per-locus gDNA-component effective length (the IPR).

The prior's only job is to split each locus's unspliced fragments between gDNA
and RNA; it does **not** attribute RNA mass to individual transcripts (that is
what the EM is for).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .region_arrays import boundary_region_indices
from .signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS

if TYPE_CHECKING:
    from ..locus import MultiLocus
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

# A region with none of these strand/type bits is intergenic — it overlaps no locus and is dropped by
# the per-locus projection, so a boundary whose left flank is such a region must be re-attributed to its
# (locus) right flank or its gDNA is lost (see boundary_owner_regions).
_RNA_SIGNATURE_BITS = BIT_EXON_POS | BIT_EXON_NEG | BIT_INTRON_POS | BIT_INTRON_NEG

# Numerical floor for the gDNA-component effective length: matches the EM's own
# default (``run_batch_locus_em_partitioned`` floors at 1.0), avoiding a zero
# denominator when the EM normalises the gDNA component's abundance.
_GDNA_EFF_LEN_FLOOR = 1.0


# gDNA eff-len reference density ρ* (eff = θ_g/ρ*): the mass-weighted CONTAINED (exon) density G_c/E_c —
# robust across capture on/off, stranded/unstranded, and gDNA level.


@dataclass(frozen=True, slots=True)
class LocusPriors:
    """Per-locus EM prior scalars (float64[n_loci], indexed by ``multi_locus_id``)."""

    gdna_prior_count: np.ndarray  # gDNA-component Dirichlet pseudocount
    rna_prior_count: np.ndarray  # RNA-group Dirichlet pseudocount (the EM splits it by evidence)
    gdna_eff_len: np.ndarray  # IPR-contracted effective length of the gDNA component


def _region_locus_shares(
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
    n_loci: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(region_idx, locus_idx, share)`` — THE region↔locus overlap, computed exactly once.

    For each region, the fractional overlap with each ``MultiLocus`` block, normalised across the loci
    it touches. A region overlapping no locus (intergenic) emits nothing and is thereby dropped. The
    overlap math is the pre-burn ``adaptive_prior._project_to_loci``'s, unchanged.

    ⭐ **Published as triples rather than folded straight into sums**, because the BOUNDARY axis needs the
    same shares (:func:`_boundary_locus_shares`) and a second traversal computing the same predicate is how
    two homes for one rule come about (``TRAPS: a-test-that-redefines``).
    """
    r_idx: list[int] = []
    l_idx: list[int] = []
    weight: list[float] = []
    if n_loci == 0:
        return (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.float64))

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
                r_idx.append(r)
                l_idx.append(lid)
                weight.append(raw_share / total)
    return (
        np.asarray(r_idx, dtype=np.int64),
        np.asarray(l_idx, dtype=np.int64),
        np.asarray(weight, dtype=np.float64),
    )


def _boundary_locus_shares(
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
    n_loci: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(edge_idx, locus_idx, share)`` — **a locus's boundaries are the boundaries that touch its regions.**

    ⭐⭐ **THE RULE** (owner, 2026-08-08): an BOUNDARY owns the fragments that cross it, a REGION owns only the
    fragments contained in it, and nothing is re-attributed between them. Every region contributes both of
    its boundaries, so a locus of ``k`` contiguous regions carries ``k + 1`` boundaries — its two OUTER ones
    included, which is correct because a fragment crossing a locus's outer boundary overlaps the locus
    and is therefore one of its EM candidates.

        share(e, L) = max( share(lo(e), L), share(hi(e), L) )

    ``max`` and not a sum: the rule is *"if a region is part of a locus, its two boundaries are part of that
    locus"*, so a boundary inherits the stronger of its two flanks' memberships rather than accumulating
    them.

    ⛔ **This replaces ``boundary_owner_regions``**, which folded a boundary's mass into ONE flank region's total so
    the region projection could reach it — a 0-bp boundary has no extent and ``_region_locus_shares`` divides
    by the region length. The fold then needed the intergenic re-key to stop a locus's far-LEFT boundary
    vanishing into its dropped intergenic flank. Projecting an boundary AS an boundary removes both.

    ⚠ **Shares can sum above 1 only for a CONTENDED boundary** — adjacent regions in different multi-loci —
    and that boundary carries no mass: any fragment crossing it overlaps transcripts in both loci, so it is
    a candidate in both and the union-find has already merged them into one multi-locus. The
    configuration is therefore unreachable for a boundary with mass, and it is *reported* by
    :func:`contended_boundaries` rather than silently renormalised.
    """
    r_idx, l_idx, w = _region_locus_shares(region_arrays, multi_loci, n_loci)
    lo, hi = boundary_region_indices(np.asarray(region_arrays.ref_id))
    if r_idx.size == 0 or lo.size == 0:
        return (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.float64))

    # per-region CSR over the (locus, share) pairs, so each boundary can gather both of its flanks' pairs
    n_regions = int(region_arrays.n_regions)
    order = np.argsort(r_idx, kind="stable")
    r_sorted, l_sorted, w_sorted = r_idx[order], l_idx[order], w[order]
    region_off = np.searchsorted(r_sorted, np.arange(n_regions + 1)).astype(np.int64)

    parts_e: list[np.ndarray] = []
    parts_p: list[np.ndarray] = []
    for flank in (lo, hi):
        counts = region_off[flank + 1] - region_off[flank]
        if not counts.sum():
            continue
        live = counts > 0
        # expand each boundary once per (locus, share) pair its flank carries
        parts_e.append(np.repeat(np.flatnonzero(live).astype(np.int64), counts[live]))
        starts = region_off[flank][live]
        c = counts[live]
        ramp = np.arange(int(c.sum()), dtype=np.int64) - np.repeat(
            np.cumsum(c) - c, c
        )
        parts_p.append(np.repeat(starts, c) + ramp)
    if not parts_e:
        return (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.float64))

    edge_ids = np.concatenate(parts_e)
    pair_ids = np.concatenate(parts_p)
    locus_ids = l_sorted[pair_ids]
    shares = w_sorted[pair_ids]

    # reduce duplicates — an boundary whose two flanks are both in L appears twice — keeping the MAX
    key = edge_ids * np.int64(n_loci) + locus_ids
    uniq, inv = np.unique(key, return_inverse=True)
    out = np.zeros(uniq.size, dtype=np.float64)
    np.maximum.at(out, inv, shares)
    return (uniq // np.int64(n_loci), uniq % np.int64(n_loci), out)


def contended_boundaries(
    region_arrays: "RegionArrays", multi_loci: "list[MultiLocus]", n_loci: int
) -> np.ndarray:
    """``int64[]`` — boundaries whose locus shares sum above 1, i.e. reached by two multi-loci at once.

    ⛔ **Reported, never renormalised.** The rule in :func:`_boundary_locus_shares` says such a boundary cannot
    carry mass; a caller that wants to *prove* that on real data needs the list, and silently rescaling
    the shares would destroy the evidence. Expected to be empty or mass-free.
    """
    e, _lid, w = _boundary_locus_shares(region_arrays, multi_loci, n_loci)
    if e.size == 0:
        return np.zeros(0, np.int64)
    n_boundaries = int(e.max()) + 1
    total = np.zeros(n_boundaries, dtype=np.float64)
    np.add.at(total, e, w)
    return np.flatnonzero(total > 1.0 + 1e-9).astype(np.int64)


def _project_regions_to_loci(
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
    n_loci: int,
    arrays: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Overlap-weighted projection of per-REGION arrays to per-locus sums.

    ⚠ Regions only. The crossing axis is projected by :func:`_boundary_locus_shares`, because a boundary is a
    first-class object and is not carried by a region.
    """
    out = {name: np.zeros(n_loci, dtype=np.float64) for name in arrays}
    r_idx, l_idx, w = _region_locus_shares(region_arrays, multi_loci, n_loci)
    for name, arr in arrays.items():
        if r_idx.size:
            np.add.at(out[name], l_idx, w * np.asarray(arr, dtype=np.float64)[r_idx])
    return out




def _sum_by_locus(
    idx: np.ndarray, lid: np.ndarray, share: np.ndarray, values: np.ndarray, n_loci: int
) -> np.ndarray:
    """``float64[n_loci]`` — ``Σ share · values[idx]``, the one weighted scatter both axes use."""
    out = np.zeros(n_loci, dtype=np.float64)
    if idx.size:
        np.add.at(out, lid, share * np.asarray(values, dtype=np.float64)[idx])
    return out


def assemble_priors(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
) -> LocusPriors:
    """Build the per-locus EM prior from the calibration result.

    ⭐⭐⭐ **A REGION OWNS THE FRAGMENTS CONTAINED IN IT; AN BOUNDARY OWNS THE FRAGMENTS THAT CROSS IT; NOTHING
    IS RE-ATTRIBUTED** (owner, 2026-08-08). A locus collects both kinds of object — its regions by genomic
    overlap (:func:`_region_locus_shares`) and its boundaries by touching those regions
    (:func:`_boundary_locus_shares`)::

        {gdna,rna}_prior_count = Σ_regions share(r,L)·mass_c_region[r]
                               + Σ_boundaries share(e,L)·mass_c_boundary[e]·q[e]

        gdna_eff_len = clamp( w·elen + (1−w)·span ),          w = C/(C+1)
            elen = Σ_regions share·min(m_r/ρ_ref, S_r) + Σ_boundaries share·min(m_e/ρ_ref, S_e)
            span = Σ_regions share·S_r                 + Σ_boundaries share·S_e

    ⭐ **THE PRIOR IS A CONSERVED FRAGMENT COUNT.** The EM adds these scalars straight to its own soft
    counts (``G = n_gdna + a_g``, ``em_solver.cpp:apply_grouped_prior_update``), where ``n_gdna`` counts
    the gDNA fragments that are candidates in this multi-locus — each ONCE, since a multi-locus is a
    connected component of transcripts linked by shared fragments. The region term is already such a count
    (a contained fragment deposits on exactly one region). Only the CROSSING term is converted, by one
    multiply against ``q = boundary_mass_per_crossing`` — the accumulator's own ``mass / count`` at that
    boundary, which undoes the ``+1``-per-crossed-boundary inflation. ``q`` is a geometry,
    ``[min(w−1,a) + min(w−1,b)] / 2(w−1)`` under a uniform field.

    ⛔ **A locus's OUTER boundaries are included, and that is the point.** A fragment crossing a locus's
    boundary overlaps the locus, so it is one of its EM candidates and must load its prior. ⚠ It follows
    that a *first-base* count of the locus's fragments is NOT this quantity — it drops exactly the
    straddlers — so an oracle built that way reads a one-way excess here that is semantics, not error.

    ⛔ **What this replaced, so it does not come back.** ``boundary_owner_regions`` folded each boundary's mass
    into ONE flank region's total, because ``_region_locus_shares`` divides by the region length and so
    cannot see a 0-bp object. The fold then needed the intergenic RE-KEY to stop a locus's far-left boundary
    vanishing into its dropped intergenic flank. Projecting an boundary AS an boundary removes both, and makes
    this the same operation ``transcript_capture_eff_lengths`` performs over a transcript's own objects.

    ⭐ **The zero-opportunity guard is now structural.** ``min(m/ρ_ref, S)`` is applied PER OBJECT, so an
    object with ``S = 0`` contributes exactly 0 to ``elen`` and 0 to ``span`` with no test and no floor.
    The predecessor summed supports first and therefore needed an explicit
    ``mass_where_there_is_opportunity`` and a ``1e-9`` cap; folding two objects into one ``min()`` also
    UNDER-contracted a captured exon whose boundary ran into a depleted intron. Both are gone.

    **The bedrock invariant — factor 1 under uniform gDNA.** Dividing each object's mass by its EFFECTIVE
    sampling support makes its density ``m/S`` exactly the true ρ under a uniform (unenriched) library,
    because the accumulator deposits ``ρ·E_f[(L−w+1)+]`` of contained mass on a region and ``ρ·E_f[w−1]``
    of crossing mass on a boundary. Every object's ``min(m/ρ_ref, S)`` then returns ``S``, so
    ``gdna_eff_len == span`` exactly: an unenriched library contracts NOTHING. Using the genomic
    ``region_size_bp`` instead would understate short-region density and fabricate a contraction with no
    capture bias present (verified factor 0.878 vs the correct 1.000). Under capture the contraction
    falls below ``span`` toward the probed footprint.

    ⚠ **The RNA prior is the UNSPLICED RNA mass only.** A spliced fragment has no gDNA candidate in the
    EM (gDNA does not splice), so it is assigned directly and counting it here would inflate the RNA side
    of a split that arbitrates only unspliced fragments. ``mass_rna_boundary`` is spliced-inclusive, so
    ``mass_rna_spliced_boundary`` is subtracted. ⛔ The JUNCTION flux is deliberately NOT added, for the same
    reason (owner ruling, 2026-07-30) — a locus whose RNA is fully spliced SHOULD get a near-zero
    ``rna_prior_count``.

    ⚠ **The contraction is SHRUNK toward the uniform span on the contained evidence** ``C``, by
    ``w = C/(C+1)`` — one pseudo-observation, no tunable. Calibration's accumulator is fed by unique
    mappers only, so a multimapper-dominated locus has little contained mass and an unreliable reference
    read; ``C = 0`` ⇒ ``span`` exactly.
    """
    if calibration.n_regions != region_arrays.n_regions:
        raise ValueError(
            f"calibration has {calibration.n_regions} regions but region_arrays has "
            f"{region_arrays.n_regions}; they must address the same partition."
        )
    n_loci = len(multi_loci)
    r_idx, r_lid, r_w = _region_locus_shares(region_arrays, multi_loci, n_loci)
    e_idx, e_lid, e_w = _boundary_locus_shares(region_arrays, multi_loci, n_loci)

    def by_region(values):
        return _sum_by_locus(r_idx, r_lid, r_w, values, n_loci)

    def by_boundary(values):
        return _sum_by_locus(e_idx, e_lid, e_w, values, n_loci)

    # ⭐ THE TWO PSEUDOCOUNTS. The region term is already a fragment count; only the crossing term is
    # converted, by the accumulator's own conserved mass-per-crossing at that boundary.
    q = np.asarray(calibration.boundary_mass_per_crossing, dtype=np.float64)
    gdna_boundary = np.asarray(calibration.mass_gdna_boundary, dtype=np.float64) * q
    rna_boundary = (
        np.maximum(
            np.asarray(calibration.mass_rna_boundary, dtype=np.float64)
            - np.asarray(calibration.mass_rna_spliced_boundary, dtype=np.float64),
            0.0,
        )
        * q
    )
    gdna_locus = np.maximum(by_region(calibration.mass_gdna_region) + by_boundary(gdna_boundary), 0.0)
    rna_locus = np.maximum(by_region(calibration.mass_rna_region) + by_boundary(rna_boundary), 0.0)

    # gDNA effective length: every object contracted against the SHARED global ρ_ref, PER OBJECT, so the
    # gDNA-vs-transcript density comparison sits on one scale. ρ_ref None (no detectable gDNA) ⇒ no
    # contraction. This is `transcript_capture_eff_lengths`' operation over the locus's object set.
    from .capture_eff_length import _global_reference_density

    region_m = np.asarray(calibration.mass_gdna_region, dtype=np.float64)
    region_s = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 0.0)
    boundary_m = np.asarray(calibration.mass_gdna_boundary, dtype=np.float64)
    boundary_s = np.maximum(np.asarray(calibration.gdna_boundary_eff_len, dtype=np.float64), 0.0)
    rho_ref = _global_reference_density(region_m, calibration.gdna_region_eff_len)
    if rho_ref is None or rho_ref <= 0.0:
        region_e, boundary_e = region_s, boundary_s
    else:
        inv = 1.0 / rho_ref
        region_e = np.minimum(region_m * inv, region_s)
        boundary_e = np.minimum(boundary_m * inv, boundary_s)

    span = by_region(region_s) + by_boundary(boundary_s)
    elen = by_region(region_e) + by_boundary(boundary_e)
    contained_ev = np.maximum(
        by_region(
            np.asarray(calibration.mass_gdna_region, dtype=np.float64)
            + np.asarray(calibration.mass_rna_region, dtype=np.float64)
        ),
        0.0,
    )
    w = contained_ev / (contained_ev + 1.0)
    eff_len = w * elen + (1.0 - w) * span

    return LocusPriors(
        gdna_prior_count=gdna_locus,
        rna_prior_count=rna_locus,
        # Clamp into [min(floor, span), span]: the 1 bp floor matches the EM's own eff-len floor but must
        # never exceed the locus's own effective span, or a degenerate sub-basepair span (a microexon-only
        # locus) would return eff_len > span, breaking eff_len ∈ (0, span].
        gdna_eff_len=np.minimum(np.maximum(eff_len, _GDNA_EFF_LEN_FLOOR), np.maximum(span, 1e-9)),
    )


__all__ = ["LocusPriors", "assemble_priors", "contended_boundaries"]
