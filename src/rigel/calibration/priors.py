"""``assemble_priors``: from a finished ``CalibrationResult`` to the per-locus EM prior inputs (layer 7, assemble).

Calibration gives the locus EM two numbers per locus: ``gdna_count``, its count of the locus's gDNA fragments, and the
gDNA component's effective length, gDNA's conserved shares of the locus's objects at their capture efficiencies
(``Σ S_r c_r + Σ M_e c_e``). The EM forms its two pseudocounts from ``gdna_count`` and its own count of the locus's
fragments (``pipeline.em_pseudocounts``), so calibration's RNA count does not enter and the two stages never have to
agree on which fragments are spliced. Nothing here re-solves calibration or attributes RNA to transcripts.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .region_arrays import boundary_region_indices, overlapping_region_runs

if TYPE_CHECKING:
    from ..locus import MultiLocus
    from .region_arrays import RegionArrays
    from .result import CalibrationResult


@dataclass(frozen=True, slots=True)
class LocusPriors:
    """Per-locus EM prior inputs (float64[n_loci], indexed by ``multi_locus_id``)."""

    gdna_count: np.ndarray  # calibration's count of the locus's gDNA fragments
    gdna_eff_len: np.ndarray  # capture-contracted effective length of the gDNA component


def _region_locus_shares(
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
    n_loci: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(region_idx, locus_idx, share)``: each region's fractional overlap with each ``MultiLocus`` block.

    The overlaps are normalised across the loci a region touches; a region overlapping no locus (intergenic) emits
    nothing. They are returned as triples so the boundary projection (:func:`_boundary_locus_shares`) reuses them
    rather than re-deriving the overlap.
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
    for ref_id in range(int(region_arrays.n_refs)):
        blocks = blocks_by_ref.get(ref_id)
        if not blocks:
            continue
        # Each block overlaps one contiguous run of regions, `[first, last)` (`overlapping_region_runs`). A region's
        # blocks are visited in sorted-block order, so every float sum has a fixed order and only overlapping pairs
        # are touched.
        b_start = np.fromiter((b[0] for b in blocks), dtype=np.int64, count=len(blocks))
        b_end = np.fromiter((b[1] for b in blocks), dtype=np.int64, count=len(blocks))
        first, last = overlapping_region_runs(
            np.full(len(blocks), ref_id, dtype=np.int64),
            b_start,
            b_end,
            starts,
            ends,
            region_arrays.ref_offsets,
        )
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
    *,
    region_shares: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(edge_idx, locus_idx, share)``: a locus's boundaries are the boundaries that touch its regions.

    A boundary owns the fragments that cross it and a region only those contained in it; nothing is re-attributed
    between them. Every region contributes both of its boundaries, so a locus of ``k`` contiguous regions carries
    ``k + 1``, its two outer ones included, because a fragment crossing an outer boundary overlaps the locus and is one
    of its EM candidates::

        share(e, L) = max( share(lo(e), L), share(hi(e), L) )

    ``max``, not a sum: a boundary belongs to a locus if either flank does, with the stronger membership. A boundary
    must be projected as a boundary rather than folded into a flank region's total, because the region shares divide
    by region length and a 0-bp boundary has no extent. Shares sum above 1 only at a contended boundary (adjacent
    regions in different multi-loci), which carries no mass: a fragment crossing it would be a candidate in both loci,
    and the union-find would have merged them. Its shares are therefore not renormalised.

    ``region_shares`` takes :func:`_region_locus_shares`'s triples when the caller already holds them, so the overlap
    is computed once on the pipeline's path.
    """
    r_idx, l_idx, w = (
        _region_locus_shares(region_arrays, multi_loci, n_loci)
        if region_shares is None
        else region_shares
    )
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
    """Build the per-locus EM prior inputs from the calibration result.

    A locus collects its regions by genomic overlap (:func:`_region_locus_shares`) and its boundaries by touching those
    regions (:func:`_boundary_locus_shares`)::

        gdna_count = Σ_regions share(r,L)·gdna_mass_region[r] + Σ_boundaries share(e,L)·gdna_mass_boundary[e]·q[e]

        gdna_eff_len = Σ_regions share·S_r·c̃_r  +  Σ_boundaries share·M_e·c̃_e

    ``c̃`` is each object's calibrated capture efficiency (``CalibrationResult.gdna_capture_efficiency_region`` /
    ``_boundary``), ``S_r`` a region's contained support and ``M_e`` gDNA's conserved share at a boundary
    (``gdna_boundary_conserved_len``).

    The count is a conserved fragment count, because the EM compares it with its own count of the locus's fragments,
    which holds each fragment once. A contained fragment deposits on exactly one region, so the region term already is
    one; the crossing term is converted by ``q = boundary_mass_per_crossing``, the accumulator's own ``mass / count``
    at the boundary, which undoes the one-per-crossed-boundary inflation (under a uniform field
    ``q = [min(w−1,a) + min(w−1,b)] / 2(w−1)``). The outer boundaries are included because their straddling fragments
    are EM candidates of the locus, so a first-base count of the locus's fragments differs from this one by exactly
    those straddlers.

    The length counts what the count counts: gDNA's share of the same objects at their own efficiencies. ``M_e``
    splits a crossing start over the boundaries its fragment crosses as the deposit rule splits the fragment
    (``effective_length.conserved_cut_shares``), so the locus holds each start once, as the count holds each fragment
    once; the crossing support ``E_f[w − 1]`` would count a start at every boundary crossed and, over pieces shorter
    than a fragment, read gDNA at a fraction of the field's density. ``M_e`` must stay gDNA's own geometry, never the
    count's pooled ``q``, which is RNA's where RNA holds most of a boundary's crossings, and the length must keep the
    boundary objects whose masses the count keeps.

    With no landscape, or a one-basin landscape, every efficiency is exactly 1 and ``gdna_eff_len == Σ S_r + Σ M_e``
    bit-identically: an unenriched library contracts nothing. Under capture each share is weighted by its object's
    yield against the bulk, the unit the transcripts' lengths share (``capture_eff_length``), so a common factor
    cancels in the EM. There is no floor and no shrinkage: the efficiencies are posterior means under the landscape,
    so a locus with little evidence reads the population's level.

    No RNA count is produced. The EM's gDNA share is a share of every fragment in the locus, spliced ones included
    (they skip the E-step but enter every M-step), while calibration's RNA count is its unspliced RNA only, so using it
    would lean the EM toward gDNA in proportion to the spliced share. The RNA side is the EM's own count less this
    one, formed where the EM's count exists.
    """
    if calibration.n_regions != region_arrays.n_regions:
        raise ValueError(
            f"calibration has {calibration.n_regions} regions but region_arrays has "
            f"{region_arrays.n_regions}; they must address the same partition."
        )
    n_loci = len(multi_loci)
    r_idx, r_lid, r_w = _region_locus_shares(region_arrays, multi_loci, n_loci)
    e_idx, e_lid, e_w = _boundary_locus_shares(
        region_arrays, multi_loci, n_loci, region_shares=(r_idx, r_lid, r_w)
    )

    def by_region(values):
        return _sum_by_locus(r_idx, r_lid, r_w, values, n_loci)

    def by_boundary(values):
        return _sum_by_locus(e_idx, e_lid, e_w, values, n_loci)

    # the gDNA count: the region term is already a fragment count; the crossing term is converted by the
    # accumulator's own mass-per-crossing at the boundary
    q = np.asarray(calibration.boundary_mass_per_crossing, dtype=np.float64)
    gdna_boundary = np.asarray(calibration.count_gdna_boundary, dtype=np.float64) * q
    gdna_locus = np.maximum(
        by_region(calibration.count_gdna_region) + by_boundary(gdna_boundary), 0.0
    )

    # the gDNA effective length: gDNA's conserved shares of the count's own objects at their efficiencies
    region_s = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 0.0)
    boundary_m = np.maximum(
        np.asarray(calibration.gdna_boundary_conserved_len, dtype=np.float64), 0.0
    )
    c_region = np.asarray(calibration.gdna_capture_efficiency_region, dtype=np.float64)
    c_boundary = np.asarray(calibration.gdna_capture_efficiency_boundary, dtype=np.float64)
    eff_len = by_region(region_s * c_region) + by_boundary(boundary_m * c_boundary)

    return LocusPriors(
        gdna_count=gdna_locus,
        # in the efficiencies' unit, the bulk; no floor, because the EM reads a zero yield as "cannot emit"
        gdna_eff_len=eff_len,
    )


def intergenic_gdna_rate(calibration: "CalibrationResult", region_arrays: "RegionArrays") -> float:
    """``λ``, the gDNA rate at an intergenic place: fragments per start of the gDNA law, in the unit every capture
    efficiency is read in, so it stands beside a locus's gDNA component, whose rate is its count over its
    capture-weighted length.

    Intergenic sequence is uncaptured, so the rate is the bulk, ``gdna_bulk_density``: the uncaptured gDNA density, a
    mode of the landscape. Under capture the pooled intergenic count over its support is not the bulk, because an
    aligner places captured genic fragments on intergenic copies. With no landscape every efficiency is 1, the unit is
    a start, and the rate is the intergenic regions' gDNA count over their contained support.
    """
    if calibration.gdna_bulk_density is not None:
        return float(calibration.gdna_bulk_density)
    from .signature import RegionType, coarse_type_array

    intergenic = coarse_type_array(region_arrays.signature) == int(RegionType.INTERGENIC)
    support = float(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64)[intergenic].sum())
    count = float(np.asarray(calibration.count_gdna_region, dtype=np.float64)[intergenic].sum())
    return count / support if support > 0.0 else 0.0


__all__ = ["LocusPriors", "assemble_priors", "intergenic_gdna_rate"]
