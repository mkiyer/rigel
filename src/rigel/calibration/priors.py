"""assemble_priors — bridge from CalibrationResult to the per-locus EM prior (PR 6).

Turns the calibration's per-region deconvolved mass + exposure into the **two
per-locus Dirichlet scalars** the locus EM consumes — ``alpha_rna_add`` and
``alpha_gdna_add`` — plus the per-locus gDNA-component effective length.

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

    alpha_gdna_add: np.ndarray  # gDNA-component Dirichlet pseudocount
    alpha_rna_add: np.ndarray  # RNA-group Dirichlet pseudocount (the EM splits it by evidence)
    gdna_eff_len: np.ndarray  # exposure-weighted physical length of the gDNA component


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


def assemble_priors(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
    *,
    prior_weight: float,
) -> LocusPriors:
    """Build the per-locus EM prior from the calibration result (PR 6 §I.2).

    Per region ``r``: the deconvolved RNA mass ``D_r = M_d_cont + M_d_left +
    M_d_right`` (D1, no ½) and the exposure-weighted physical length
    ``e_r = ω_r · L_phys_r``. These project to loci by overlap share ``φ``; then
    per locus ``L``::

        alpha_rna_add[L]  = prior_weight · Σ_r φ_{L,r} · D_r
        alpha_gdna_add[L] = prior_weight · ρ_0 · Σ_r φ_{L,r} · e_r     (modeled gDNA mass)
        gdna_eff_len[L]   = max(1, Σ_r φ_{L,r} · e_r)
    """
    if calibration.n_regions != region_arrays.n_regions:
        raise ValueError(
            f"calibration has {calibration.n_regions} regions but region_arrays has "
            f"{region_arrays.n_regions}; they must address the same partition."
        )

    n_loci = len(multi_loci)
    d_region = calibration.mass_d_contained + calibration.mass_d_left + calibration.mass_d_right
    e_region = calibration.omega * region_arrays.region_size_bp

    proj = _project_regions_to_loci(
        region_arrays, multi_loci, n_loci, {"d": d_region, "e": e_region}
    )
    d_locus = proj["d"]
    e_locus = proj["e"]

    w = float(prior_weight)
    return LocusPriors(
        alpha_gdna_add=w * float(calibration.rho_0) * e_locus,
        alpha_rna_add=w * d_locus,
        gdna_eff_len=np.maximum(e_locus, _GDNA_EFF_LEN_FLOOR),
    )


__all__ = ["LocusPriors", "assemble_priors"]
