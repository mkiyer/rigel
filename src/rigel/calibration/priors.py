"""assemble_priors — the bridge from CalibrationResult to the per-locus EM prior.

Turns the calibration's per-object deconvolved mass and geometric length into the two per-locus
Dirichlet scalars the locus EM consumes — ``rna_prior_count`` and ``gdna_prior_count`` — plus the
per-locus gDNA-component effective length (the inverse participation ratio of the deconvolved gDNA
mass over its supports).

The prior's only job is to split each locus's unspliced fragments between gDNA and RNA; it does not
attribute RNA mass to individual transcripts, which is what the EM is for.

Layer: LAYER 7. It reads a finished `CalibrationResult` and never re-solves anything.
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

# A region with none of these strand/type bits is intergenic — it overlaps no locus and is dropped by the
# per-locus projection.
# ⛔ No re-attribution happens here: a boundary is projected AS a boundary, so a boundary beside a
# dropped intergenic region does not have to be re-keyed to its other flank to keep its gDNA.
# The constant's one consumer is ``tests/calibration/test_prior_vs_oracle.py``, which imports it to
# rebuild the same in-locus predicate rather than restate the bits.
_RNA_SIGNATURE_BITS = BIT_EXON_POS | BIT_EXON_NEG | BIT_INTRON_POS | BIT_INTRON_NEG

# Numerical floor for the gDNA-component effective length: matches the EM's own
# default (``run_batch_locus_em_partitioned`` floors at 1.0), avoiding a zero
# denominator when the EM normalises the gDNA component's abundance.
_GDNA_EFF_LEN_FLOOR = 1.0


@dataclass(frozen=True, slots=True)
class LocusPriors:
    """Per-locus EM prior scalars (float64[n_loci], indexed by ``multi_locus_id``)."""

    gdna_prior_count: np.ndarray  # gDNA-component Dirichlet pseudocount
    rna_prior_count: np.ndarray  # RNA-group Dirichlet pseudocount (the EM splits it by evidence)
    gdna_eff_len: np.ndarray  # capture-contracted effective length of the gDNA component


def _region_locus_shares(
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
    n_loci: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(region_idx, locus_idx, share)`` — THE region↔locus overlap, computed exactly once.

    For each region, the fractional overlap with each ``MultiLocus`` block, normalised across the loci
    it touches. A region overlapping no locus (intergenic) emits nothing and is thereby dropped.

    Published as triples rather than folded straight into sums, because the BOUNDARY axis needs the
    same shares (:func:`_boundary_locus_shares`) and a second traversal computing the same predicate is how
    two homes for one rule come about.
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

    starts = np.asarray(region_arrays.start, dtype=np.int64)
    ends = np.asarray(region_arrays.end, dtype=np.int64)
    ref_offsets = region_arrays.ref_offsets
    for ref_id in range(int(region_arrays.n_refs)):
        blocks = blocks_by_ref.get(ref_id)
        if not blocks:
            continue
        lo, hi = int(ref_offsets[ref_id]), int(ref_offsets[ref_id + 1])
        # The regions partition the reference, so each block overlaps one contiguous run of them: the
        # first region ending after the block starts through the last starting before it ends. Visiting
        # each region's overlapping blocks in sorted-block order keeps every float sum in the same order
        # as a scan over all earlier blocks, while touching only the pairs that overlap.
        b_start = np.fromiter((b[0] for b in blocks), dtype=np.int64, count=len(blocks))
        b_end = np.fromiter((b[1] for b in blocks), dtype=np.int64, count=len(blocks))
        first = lo + np.searchsorted(ends[lo:hi], b_start, side="right")
        last = lo + np.searchsorted(starts[lo:hi], b_end, side="left")
        width = np.maximum(last - first, 0)
        pair_block = np.repeat(np.arange(len(blocks), dtype=np.int64), width)
        pair_region = np.repeat(first, width) + (
            np.arange(int(width.sum()), dtype=np.int64) - np.repeat(np.cumsum(width) - width, width)
        )
        order = np.argsort(pair_region, kind="stable")
        pair_region, pair_block = pair_region[order], pair_block[order]
        cut = np.flatnonzero(np.diff(pair_region)) + 1
        for run in np.split(np.arange(pair_region.size), cut):
            if run.size == 0:
                continue
            r = int(pair_region[run[0]])
            r_start = int(starts[r])
            r_end = int(ends[r])
            r_len = r_end - r_start
            if r_len <= 0:
                continue
            raw: dict[int, float] = {}
            for j in pair_block[run]:
                bs, be, lid = blocks[int(j)]
                overlap = min(be, r_end) - max(bs, r_start)
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
    """``(edge_idx, locus_idx, share)`` — a locus's boundaries are the boundaries that touch its regions.

    THE RULE: a BOUNDARY owns the fragments that cross it, a REGION owns only the
    fragments contained in it, and nothing is re-attributed between them. Every region contributes both of
    its boundaries, so a locus of ``k`` contiguous regions carries ``k + 1`` boundaries — its two OUTER ones
    included, which is correct because a fragment crossing a locus's outer boundary overlaps the locus
    and is therefore one of its EM candidates.

        share(e, L) = max( share(lo(e), L), share(hi(e), L) )

    ``max`` and not a sum: the rule is *"if a region is part of a locus, its two boundaries are part of that
    locus"*, so a boundary inherits the stronger of its two flanks' memberships rather than accumulating
    them.

    ⛔ Do not fold a boundary's mass into one flank region's total instead: ``_region_locus_shares``
    divides by the region length and a 0-bp boundary has no extent, so such a fold then needs an
    intergenic re-key to stop a locus's far-left boundary vanishing into its dropped intergenic flank.
    Projecting a boundary AS a boundary removes both.

    Shares can sum above 1 only for a CONTENDED boundary — adjacent regions in different multi-loci —
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
        ramp = np.arange(int(c.sum()), dtype=np.int64) - np.repeat(np.cumsum(c) - c, c)
        parts_p.append(np.repeat(starts, c) + ramp)
    if not parts_e:
        return (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.float64))

    edge_ids = np.concatenate(parts_e)
    pair_ids = np.concatenate(parts_p)
    locus_ids = l_sorted[pair_ids]
    shares = w_sorted[pair_ids]

    # reduce duplicates — a boundary whose two flanks are both in L appears twice — keeping the MAX
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

    Regions only. The crossing axis is projected by :func:`_boundary_locus_shares`, because a boundary is a
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

    A REGION OWNS THE FRAGMENTS CONTAINED IN IT; A BOUNDARY OWNS THE FRAGMENTS THAT CROSS IT; NOTHING
    IS RE-ATTRIBUTED. A locus collects both kinds of object — its regions by genomic
    overlap (:func:`_region_locus_shares`) and its boundaries by touching those regions
    (:func:`_boundary_locus_shares`)::

        {gdna,rna}_prior_count = Σ_regions share(r,L)·mass_c_region[r]
                               + Σ_boundaries share(e,L)·mass_c_boundary[e]·q[e]

        gdna_eff_len = clamp( Σ_regions share·S_r·c̃_r  +  Σ_boundaries share·S_e·c̃_e )

    ``c̃`` is each object's capture efficiency, the calibration's own
    (`CalibrationResult.gdna_capture_efficiency_region` / ``_boundary``; `capture_efficiency`) and ``S``
    its support.

    THE PRIOR IS A CONSERVED FRAGMENT COUNT. The EM adds these scalars straight to its own soft
    counts (``G = n_gdna + a_g``, ``em_solver.cpp:apply_grouped_prior_update``), where ``n_gdna`` counts
    the gDNA fragments that are candidates in this multi-locus — each ONCE, since a multi-locus is a
    connected component of transcripts linked by shared fragments. The region term is already such a count
    (a contained fragment deposits on exactly one region). Only the CROSSING term is converted, by one
    multiply against ``q = boundary_mass_per_crossing`` — the accumulator's own ``mass / count`` at that
    boundary, which undoes the ``+1``-per-crossed-boundary inflation. ``q`` is a geometry,
    ``[min(w−1,a) + min(w−1,b)] / 2(w−1)`` under a uniform field.

    ⛔ A locus's OUTER boundaries are included, and that is the point: a fragment crossing a locus's
    boundary overlaps the locus, so it is one of its EM candidates and must load its prior. It follows
    that a *first-base* count of the locus's fragments is NOT this quantity — it drops exactly the
    straddlers — so an oracle built that way reads a one-way excess here that is semantics, not error.

    THE gDNA EFFECTIVE LENGTH COUNTS THE SAME OBJECTS AS THE COUNT. The prior's count is the calibration's
    gDNA mass on the locus's regions and on its boundaries, so the length the EM divides that count by is
    those same objects' supports at their own efficiencies: every region's contained support ``S_r`` at
    ``c̃_r`` and every boundary's crossing support ``S_e`` at ``c̃_e``. ⛔ Two other forms are refused by
    measurement. A length over the locus's BASES at the regions' efficiencies alone — the transcript
    ruler's form, right for a template's capture — drops the boundary objects whose masses the count
    keeps, and where the calibration's crossing masses sit above their geometry the gDNA component reads
    denser than its objects, over-claims the exonic unspliced fragments and every probed gene under-calls
    (the test chromosome's `g50 ss.99 ON` row: gene-level Σ|Δ| 25,633 → 38,174 against 23,967 here). And
    the boundary support converted by ``q`` to count each crossing start once, as the count counts each
    fragment once, collapses the length where pieces are short: the gDNA component then saturates, the
    EM stops responding to its own gDNA pseudocount (the thermometer's injection gate on a contaminated
    toy goes insensitive) and both capture-OFF strata read 1–2 % worse, for 21,733 on that one row; the
    crossing support enters at its full opportunity, as the count's ``q`` undoes an inflation of MASS
    and not of starts.

    The bedrock invariant — factor 1 under uniform gDNA. With no reference every efficiency is exactly
    1 and ``gdna_eff_len == span == Σ S`` bit-identically: an unenriched library contracts NOTHING and
    reads what it read before the efficiencies existed. Under capture a depleted object contributes its
    support at its efficiency and the length contracts toward the probed footprint.

    The RNA prior is the UNSPLICED RNA mass only. A spliced fragment has no gDNA candidate in the
    EM (gDNA does not splice), so it is assigned directly and counting it here would inflate the RNA side
    of a split that arbitrates only unspliced fragments. ``count_rna_boundary`` is spliced-inclusive, so
    ``count_rna_spliced_boundary`` is subtracted. ⛔ The SJ flux is deliberately NOT added, for the same
    reason — a locus whose RNA is fully spliced SHOULD get a near-zero ``rna_prior_count``.

    No floor and no shrinkage: the efficiencies are posterior means under the population landscape, so a
    locus with little evidence reads the population's own level, never a fabricated 0 and never the
    uncontracted span (the multimapper floor ``C/(C+1)`` this replaces cost the unprobed class 30×;
    `ISSUES: ruler-multimapper-floor-caps-the-correction`).
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

    # THE TWO PSEUDOCOUNTS. The region term is already a fragment count; only the crossing term is
    # converted, by the accumulator's own conserved mass-per-crossing at that boundary.
    q = np.asarray(calibration.boundary_mass_per_crossing, dtype=np.float64)
    gdna_boundary = np.asarray(calibration.count_gdna_boundary, dtype=np.float64) * q
    rna_boundary = (
        np.maximum(
            np.asarray(calibration.count_rna_boundary, dtype=np.float64)
            - np.asarray(calibration.count_rna_spliced_boundary, dtype=np.float64),
            0.0,
        )
        * q
    )
    gdna_locus = np.maximum(
        by_region(calibration.count_gdna_region) + by_boundary(gdna_boundary), 0.0
    )
    rna_locus = np.maximum(by_region(calibration.count_rna_region) + by_boundary(rna_boundary), 0.0)

    # THE gDNA EFFECTIVE LENGTH: the count's own objects at their own supports and efficiencies.
    region_s = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 0.0)
    boundary_s = np.maximum(np.asarray(calibration.gdna_boundary_eff_len, dtype=np.float64), 0.0)
    c_region = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    c_boundary = np.asarray(calibration.gdna_capture_efficiency_boundary, dtype=np.float64)
    span = by_region(region_s) + by_boundary(boundary_s)
    eff_len = by_region(region_s * c_region) + by_boundary(boundary_s * c_boundary)

    return LocusPriors(
        gdna_prior_count=gdna_locus,
        rna_prior_count=rna_locus,
        # Clamp into [min(floor, span), span]: the 1 bp floor matches the EM's own eff-len floor but must
        # never exceed the locus's own uncontracted span, or a degenerate sub-basepair span (a
        # microexon-only locus) would return eff_len > span, breaking eff_len ∈ (0, span].
        gdna_eff_len=np.minimum(np.maximum(eff_len, _GDNA_EFF_LEN_FLOOR), np.maximum(span, 1e-9)),
    )


__all__ = ["LocusPriors", "assemble_priors", "contended_boundaries"]
