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

    The per-region gDNA is the **conservation-correct node mass over genomic bp**: ``gdna_region`` =
    ``contained + left + right`` (the gDNA overlapping the region's genomic interval) and the per-region
    length is the bare genomic ``region_size_bp``. Both conserve exactly — ``Σ_r gdna_region`` = the total
    deconvolved gDNA mass (CONSERVATION OF MASS) and ``Σ_r region_size_bp`` = the covered genomic span
    (CONSERVATION OF GENOMIC LENGTH, since regions tile). ``gdna_region`` and ``rna_region`` project to
    loci by genomic-overlap ``share``::

        gdna_prior_count = Σ_r share * gdna_region    (deconvolved gDNA count)
        rna_prior_count  = Σ_r share * rna_region     (deconvolved UNSPLICED RNA count;
                                                       spliced mass withheld — see below)
        gdna_eff_len     = (G+1)² / [ Σ share*(gdna²/L) + (2G+1)/span ],  G = Σ share*gdna,
                           L = region_size_bp,  span = Σ share*region_size_bp   (capped at span)

    This REPLACES the earlier boundary-flux transport + the FL-aware ``gdna_geom_len`` (= contained
    support + both boundary supports), which summed to ~8% more than the genomic span (≈2× on small
    exons — it double-counts each seam and adds FL-mass on top of the contained support), violating
    genomic-length conservation and biasing the IPR span. ``gdna_eff_len`` is the **inverse
    participation ratio** of the per-region gDNA *density* ``gdna_region/region_size_bp`` (its reciprocal
    makes the gDNA component's per-position rate equal the local gDNA density, so under capture the
    support contracts to the probed exons and gDNA competes at its true density), correct by
    construction — uniform gDNA (``g_r ∝ region_size_bp``) ⇒ ``eff_len = span`` (the genomic span).
    **Laplace-smoothed** by one fragment-equivalent of uniform per-base prior mass: adding
    ``c_r = L_r/span`` to each ``gdna_region`` and forming ``(G+1)²/Σ(g̃²/L)`` gives the closed
    form above, since ``Σ 2·g·c/L = 2G/span`` and ``Σ c²/L = 1/span``.

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

    # CONSERVATION-CORRECT per-region gDNA (PR-1). Per region: the gDNA mass overlapping its genomic
    # interval — contained + the two boundary sides' portions — over the region's genomic bp. Both
    # quantities conserve EXACTLY:
    #   • mass:   Σ_r (contained + left + right) = the total deconvolved gDNA mass — each crossing
    #             fragment's mass is split once across the two regions it spans (the two sides);
    #   • length: Σ_r region_size_bp = the covered genomic span — regions tile the genome.
    # This REPLACES the boundary-flux transport + the FL-aware gdna_geom_len, which summed to ~8% more
    # than the genomic span (≈2× on small exons: gdna_geom_len = region_eff + both sides double-counts
    # each seam and adds FL-mass on top of the contained support) — a genomic-length-conservation
    # violation that biased the IPR span (the cap + the uniform reference). With genomic-bp lengths the
    # IPR is correct by construction: uniform gDNA (g_r ∝ region_size_bp) ⇒ eff_len = span (the genomic
    # span); only capture-concentrated gDNA contracts. No transport — the length-bias it compensated
    # for is gone (the FL-aware length was the source of that bias).
    length = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    gdna_region = (
        np.asarray(calibration.mass_gdna_contained, dtype=np.float64)
        + np.asarray(calibration.mass_gdna_left, dtype=np.float64)
        + np.asarray(calibration.mass_gdna_right, dtype=np.float64)
    )
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
        gdna_sq_over_len = np.where(
            length > 0.0, gdna_region**2 / np.maximum(length, 1e-9), 0.0
        )

    proj = _project_regions_to_loci(
        region_arrays,
        multi_loci,
        len(multi_loci),
        {
            "gdna": gdna_region,
            "rna": rna_region,
            "support": gdna_sq_over_len,
            "span": length,
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
    # The add-one is the canonical Laplace prior — there is NO tunable shrinkage
    # constant. It shrinks the support toward the uniform span exactly in
    # proportion to the gDNA evidence: G = 0 → span (uniform; the multimap-blind
    # null, no special case); a tiny sparse-but-concentrated mass → still ≈ span
    # (the +1 uniform prior dominates), so the EM cannot amplify it past the
    # calibration's call; abundant gDNA (capture, G ≫ 1) → the empirical
    # concentration (the IPR). Capped at span: the gDNA cannot be more spread than
    # the locus's uniform length.
    with np.errstate(divide="ignore", invalid="ignore"):
        safe_span = np.maximum(span, 1e-9)
        smoothed_support = support_locus + (2.0 * gdna_locus + 1.0) / safe_span
        eff_len = np.minimum((gdna_locus + 1.0) ** 2 / smoothed_support, span)

    return LocusPriors(
        gdna_prior_count=gdna_locus,
        rna_prior_count=proj["rna"],
        gdna_eff_len=np.maximum(eff_len, _GDNA_EFF_LEN_FLOOR),
    )


__all__ = ["LocusPriors", "assemble_priors"]
