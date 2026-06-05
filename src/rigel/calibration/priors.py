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


def _transport_boundary_flux(
    contained: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    length: np.ndarray,
    e_cap: np.ndarray,
    ref_id: np.ndarray,
    max_iter: int = 8,
) -> np.ndarray:
    """Length-bias-free boundary-flux transport of per-region gDNA mass.

    Each internal boundary's pooled gDNA mass is re-attributed to its two sides ∝ ``ω·𝓔``
    — exposure (density ``g/L``, length-bias-free) × the directional boundary effective
    length ``𝓔(L)=E[min(ℓ,L)]`` (``e_cap``). This moves capture smear off the unexposed
    (e.g. intronic) side and onto the probed side that generated it. Iterates until the
    total mass moved is sub-count (< 1 fragment-equivalent), capped at ``max_iter``; total
    mass is conserved (ref-edge / cross-ref sides stay in place). See
    docs/futureprs/phase6_boundary_flux_transport_plan.md.
    """
    contained = np.asarray(contained, dtype=np.float64)
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    length = np.maximum(np.asarray(length, dtype=np.float64), 1e-9)
    e_cap = np.asarray(e_cap, dtype=np.float64)
    r = contained.shape[0]
    if r <= 1:
        return contained + left + right
    same = np.asarray(ref_id)[:-1] == np.asarray(ref_id)[1:]  # boundary (i,i+1) internal
    ltot = float(length.sum())
    g = contained + left + right
    prev = g
    for _ in range(max_iter):
        gdna_density_global = g.sum() / ltot if ltot > 0.0 else 0.0
        exposure = np.where(
            gdna_density_global > 0.0, (g / length) / max(gdna_density_global, 1e-12), 1.0
        )
        w = exposure * e_cap
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
    *,
    prior_weight: float,
) -> LocusPriors:
    """Build the per-locus EM prior from the acyclic calibration result.

    Boundary-crossing gDNA is first re-attributed to its origin region by the
    length-bias-free boundary-flux transport (``ω·𝓔``); the transported per-region gDNA
    ``G_r`` and RNA ``D_r`` project to loci by genomic-overlap share ``φ``::

        alpha_gdna_add[L] = prior_weight · Σ_r φ·G_r            (deconvolved gDNA count)
        alpha_rna_add[L]  = prior_weight · Σ_r φ·D_r            (deconvolved RNA count)
        IPR[L]            = (Σ_r φ·G_r)² / Σ_r φ·(G_r²/E_r)     (concentrated support)
        gdna_eff_len[L]   = (1−π)·span + π·IPR,  π = G/(G+κ)    (power-shrunk support)

    ``E_r = gdna_geom_len`` is the FL-aware gDNA support of region ``r`` — the bases a gDNA
    fragment can overlap (contained + both boundary crossings ≈ ``R + L̄``), not the bare
    region length. The IPR is the effective support of the gDNA over those lengths; its
    reciprocal makes the gDNA component's per-position rate equal the local gDNA density, so
    under capture (gDNA piled on exons) the support contracts to the exons and gDNA competes
    at its true density.

    A *concentrated* support, though, is only trustworthy in proportion to the gDNA **count**
    backing it — a tiny sparse mass that merely *looks* concentrated would get a tiny eff-len
    and the EM would amplify it far past the calibration's call. So the eff-len is shrunk from
    the IPR toward the uniform geometric ``span`` by the per-locus power ``π = G/(G+κ)``
    (``κ = config.gdna_eff_len_shrink ≈ 1/φ``): abundant gDNA → ``π→1`` → trust the
    concentration; sparse gDNA → ``π→0`` → uniform support. A locus with no deconvolved gDNA
    falls back to ``span`` (multimap blindness). See
    docs/futureprs/phase6_boundary_flux_transport_plan.md.
    """
    if calibration.n_regions != region_arrays.n_regions:
        raise ValueError(
            f"calibration has {calibration.n_regions} regions but region_arrays has "
            f"{region_arrays.n_regions}; they must address the same partition."
        )

    # Transport uses the physical region length for its length-bias-free density (ω); the
    # IPR uses the FL-aware gDNA support E_r = gdna_geom_len (contained + boundary crossings
    # ≈ R + L̄) — the bases a gDNA fragment can actually overlap, not the bare region length.
    length = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    geom = np.asarray(calibration.gdna_geom_len, dtype=np.float64)
    g_region = _transport_boundary_flux(
        calibration.mass_gdna_contained,
        calibration.mass_gdna_left,
        calibration.mass_gdna_right,
        length,
        calibration.gdna_boundary_len,
        np.asarray(region_arrays.ref_id),
    )
    d_region = (
        calibration.mass_rna_contained + calibration.mass_rna_left + calibration.mass_rna_right
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        g2_over_e = np.where(geom > 0.0, g_region**2 / np.maximum(geom, 1e-9), 0.0)

    proj = _project_regions_to_loci(
        region_arrays,
        multi_loci,
        len(multi_loci),
        {"g": g_region, "d": d_region, "s": g2_over_e, "span": geom},
    )
    big_g, s_g, span = proj["g"], proj["s"], proj["span"]
    with np.errstate(divide="ignore", invalid="ignore"):
        ipr = np.where(s_g > 0.0, big_g**2 / np.maximum(s_g, 1e-30), span)
    # Power-shrink the support from the concentrated IPR toward the uniform geometric span,
    # in proportion to the gDNA count backing the concentration: π = G/(G+κ). Abundant gDNA
    # (capture) → π→1 → trust the IPR; sparse/spurious gDNA → π→0 → uniform span, so the EM
    # cannot amplify a tiny concentrated mass past the calibration's call. κ ≈ 1/φ.
    kappa = float(calibration.config.gdna_eff_len_shrink)
    if kappa > 0.0:
        power = np.where(big_g > 0.0, big_g / (big_g + kappa), 0.0)
        eff_len = (1.0 - power) * span + power * ipr
    else:
        eff_len = ipr

    w = float(prior_weight)
    return LocusPriors(
        alpha_gdna_add=w * big_g,
        alpha_rna_add=w * proj["d"],
        gdna_eff_len=np.maximum(eff_len, _GDNA_EFF_LEN_FLOOR),
    )


__all__ = ["LocusPriors", "assemble_priors"]
