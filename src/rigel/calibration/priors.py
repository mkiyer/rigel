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

from .signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS

if TYPE_CHECKING:
    from ..locus import MultiLocus
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

# A region with none of these strand/type bits is intergenic — it overlaps no locus and is dropped by
# the per-locus projection, so a seam whose left flank is such a region must be re-attributed to its
# (locus) right flank or its gDNA is lost (see _gdna_region_node_arrays).
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


def _gdna_region_node_arrays(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-region node arrays ``(gdna_region, support_len, pooled, seam_len)`` for the gDNA region +
    pooled-seam nodes — the SHARED node model with ``transcript_capture_eff_lengths`` (via _pooled_seam_arrays).

    The density-correct, transport-free gDNA node model (``docs/CARRY_FORWARD.md`` §8 —
    the first-principles rev). The node set over a region chain is the per-region CONTAINED node plus
    one POOLED SEAM node per internal boundary (genomically adjacent, same reference), keyed to its
    left-flank region ``r``::

        region node r:     mass m_r = mass_gdna_contained[r],
                           effective support S_r = gdna_region_eff_len[r] = E[max(0, L_r − ℓ)]
        seam node (r,r+1): mass m_s = mass_gdna_right[r] + mass_gdna_left[r+1]   (the two halves POOLED),
                           effective support S_s = gdna_boundary_len[r] + gdna_boundary_len[r+1]  (the SUM of
                           the two flanking per-side density lengths gdna_boundary_len)

    **Why these supports — the bedrock invariant.** Driving the reference accumulator under uniform
    genomic gDNA at density ρ, the expected masses are ``m_r = ρ·E[max(0,L_r−ℓ)]`` (a contained fragment
    must FIT) and ``m_s = ρ·(E[min(ℓ,L_r)] + E[min(ℓ,L_{r+1})])/2`` (each side captures only ``min(ℓ,L_side)``
    of a crossing fragment, so the pooled mass is ρ times the SUM of the two stored side density
    lengths — each of which is already E[min(ℓ,L)]/2, so the sum is ½·(E[min_r]+E[min_{r+1}]). ⚠ Naming
    it the AVERAGE of gdna_boundary_len is what applied the ½ twice (D6, fixed 2026-07-29).
    Dividing each node's mass by these supports makes EVERY node density ``m_n/S_n = ρ``, so the
    enrichment contraction ``min(m_n/ρ_ref, S_n)`` (applied per node in ``assemble_priors`` against the
    shared ρ_ref) returns ``S_n`` EXACTLY (contraction factor 1) under uniform gDNA — an unenriched library
    contracts nothing, EVEN for regions shorter than ``E[ℓ]``.
    Two divisors are WRONG: the genomic ``region_size_bp`` understates short-region density (verified
    factor 0.878 under uniform); and the count crossing length ``E[ℓ]`` over-states the seam support for
    short flanks (a fragment can only deposit ``min(ℓ,L_side)``, not ``ℓ``), under-contracting exon-flank
    seams and inflating the gDNA→RNA leak — the summed side density length is the deposition-faithful
    support.

    **Why POOL the seam, not split it.** The two halves ``s_L,s_R`` are one physical crossing event, so the
    pooled node (one node at the summed support) is the faithful representation; splitting them into two
    separate nodes would double-count the crossing and over-contract.

    Returns ``(gdna_region, support_len, pooled, seam_len)``, each float64[R] keyed to region ``r``::

        gdna_region[r]   = m_r + m_s(r,r+1)   total node mass on region r (contained + pooled seam)
        support_len[r]   = S_r + S_s(r,r+1)   total effective support (Σ S), for the span + mass projection
        pooled[r]        = m_s(r,r+1)         the pooled-seam mass ALONE  (per-node contraction in priors)
        seam_len[r]      = S_s(r,r+1)         the pooled-seam support ALONE

    ``assemble_priors`` contracts PER NODE — ``min(contained/ρ_ref, S_r) + min(pooled/ρ_ref, S_s)`` —
    matching ``transcript_capture_eff_lengths``, NOT over the folded ``gdna_region`` (which would
    under-contract a captured exon whose seam runs into a depleted intron).

    Mass conservation (no mass moved — transport-free): ``Σ gdna_region = Σ contained +
    Σ_{internal} (right[r] + left[r+1])`` — every non-terminal boundary side counted exactly once
    (terminal / cross-reference sides carry zero on real data and are excluded).
    """
    from .capture_eff_length import (
        _pooled_seam_arrays,
    )  # THE shared seam node model (transcript + gDNA)

    contained = np.asarray(calibration.mass_gdna_contained, dtype=np.float64)
    region_eff_len = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9)
    ref_id = np.asarray(region_arrays.ref_id)
    n = contained.shape[0]

    pooled = np.zeros(n, dtype=np.float64)  # pooled seam mass, attributed to a flank region
    seam_len = np.zeros(n, dtype=np.float64)  # seam effective support, attributed to a flank region
    if n > 1:
        seam_mass, seam_support = _pooled_seam_arrays(
            calibration, region_arrays
        )  # left-keyed, length n
        same = ref_id[:-1] == ref_id[1:]  # internal seam: genomically adjacent, same reference
        # ATTRIBUTE each seam to a flank REGION so the locus projection picks it up. Default: the LEFT
        # flank r. BUT a locus's far-LEFT outer boundary is an intergenic→(exon/intron) seam whose left
        # flank is INTERGENIC — a region that overlaps no locus and is dropped by _project_regions_to_loci.
        # Keying to the left flank there SILENTLY LOSES that boundary's crossing gDNA, under-counting the
        # locus gDNA prior AND inflating the eff-length. The far-RIGHT boundary is already kept (its left
        # flank is the locus's last region), so this restores symmetry: attribute the seam to the RIGHT
        # flank whenever the left flank is intergenic (no RNA-signature bits) and the right flank is not.
        sig = np.asarray(region_arrays.signature).astype(np.int64)
        ig = (
            sig & _RNA_SIGNATURE_BITS
        ) == 0  # intergenic: no exon/intron bit ⇒ dropped by the projection
        rekey_right = same & ig[:-1] & ~ig[1:]  # far-left outer boundary: intergenic → locus region
        owner = np.where(rekey_right, np.arange(1, n), np.arange(0, n - 1))
        np.add.at(
            pooled, owner, seam_mass[:-1]
        )  # a first-region node may own its right seam + a rekeyed one
        np.add.at(seam_len, owner, seam_support[:-1])

    gdna_region = contained + pooled
    support_len = region_eff_len + seam_len
    return gdna_region, support_len, pooled, seam_len


def assemble_priors(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
) -> LocusPriors:
    """Build the per-locus EM prior from the calibration result.

    The gDNA node set is the **density-correct, transport-free** region + pooled-seam model
    (``docs/CARRY_FORWARD.md`` §8; see ``_gdna_region_node_arrays``):

        region node r:    mass = mass_gdna_contained[r],   effective support S_r = E[max(0, L_r − ℓ)]
        seam node (r,r+1): mass = mass_gdna_right[r] + mass_gdna_left[r+1] (the two halves POOLED),
                           effective support S_s = gdna_boundary_len[r] + gdna_boundary_len[r+1]  (the SUM of
                           the flanking per-side density lengths gdna_boundary_len)

    keyed to the left-flank region r. The region + seam masses / participations project to loci by
    genomic-overlap ``share``::

        gdna_prior_count = Σ share * (contained + seam)                       (deconvolved gDNA count)
        rna_prior_count  = Σ share * rna_region         (UNSPLICED RNA; spliced withheld — see below)
        gdna_eff_len     = (G+1)² / [ Σ share*(contained²/S_r + seam²/S_s) + (2G+1)/span ], capped at span
                           G = Σ share*(contained+seam),  span = Σ share*(S_r + S_s)  (EFFECTIVE support)

    **The bedrock invariant — factor 1 under uniform gDNA.** Dividing each node's mass by its EFFECTIVE
    sampling support (``S_r`` for contained, ``S_s = ½·(E[min(ℓ,L_r)]+E[min(ℓ,L_{r+1})])`` for the pooled
    crossing — the SUM of the two stored per-side density lengths, NOT ``E[ℓ]`` which over-states it) makes the per-node
    density ``m_n/S_n`` exactly the true ρ under a uniform (unenriched) library — because the accumulator
    deposits ``ρ·E[max(0,L−ℓ)]`` of contained mass and ``ρ·S_s`` of pooled crossing mass.
    The Laplace-smoothed IPR then returns ``span`` EXACTLY (``eff_len = span`` ⇒ contraction factor 1):
    an unenriched library contracts NOTHING. Using the genomic ``region_size_bp`` instead would
    understate short-region density and fabricate a contraction even with no capture bias (verified
    factor 0.878 vs the correct 1.000) — that was the latent defect this redesign removes. Under
    capture (concentrated gDNA) the IPR contracts below ``span`` toward the probed footprint, so the
    gDNA component competes at its true local density. This is **transport-free**: no mass is moved
    (no boundary-flux redistribution — that non-physical heuristic is gone), and the pooled seam
    quarantines the captured intron↔exon crossing mass at crossing density — recovering
    the gDNA concentration a region-dilution divisor would lose, by Cauchy–Schwarz strictly better than
    splitting the seam into two side nodes.

    **Laplace-smoothed** by one fragment-equivalent of uniform support (the ``(2G+1)/span`` term), the
    canonical add-one prior with no tunable constant: ``G = 0`` ⇒ ``span`` exactly; abundant gDNA
    (capture, ``G ≫ 1``) ⇒ the empirical concentration; in between, the IPR is shrunk toward the
    uniform ``span`` in proportion to the gDNA evidence, so the EM cannot amplify a tiny concentrated
    mass past the calibration's call.
    """
    if calibration.n_regions != region_arrays.n_regions:
        raise ValueError(
            f"calibration has {calibration.n_regions} regions but region_arrays has "
            f"{region_arrays.n_regions}; they must address the same partition."
        )

    # Density-correct, transport-free gDNA node model (docs/CARRY_FORWARD.md §8): per-region
    # CONTAINED node (effective support gdna_region_eff_len = E[max(0,L−ℓ)]) + one POOLED SEAM node per
    # internal boundary (support = the SUM of the flanking gdna_boundary_len = ½·Σ E[min(ℓ,L)]),
    # keyed to the left-flank region — the SAME node model _pooled_seam_arrays gives the transcript
    # contraction (EFFECTIVE, not genomic, supports; the factor-1-under-uniform bedrock).
    gdna_region, support_len, pooled, seam_len = _gdna_region_node_arrays(
        calibration, region_arrays
    )

    # SHARED global reference density — the SAME ρ_ref every transcript contracts against, so the
    # gDNA-vs-transcript density comparison sits on ONE scale. The enrichment contraction is applied PER
    # NODE (contained region node + pooled-seam node separately), identical to transcript_capture_eff_lengths:
    #   elen = min(contained/ρ_ref, S_region) + min(pooled/ρ_ref, S_seam).
    # Folding the two into one min() would UNDER-contract a captured exon whose seam runs into a depleted
    # intron (up to ~13% per region under capture, verified). ρ_ref None (no detectable gDNA) ⇒ elen =
    # support (no contraction). See docs/CARRY_FORWARD.md.
    from .capture_eff_length import _global_reference_density

    contained = np.asarray(calibration.mass_gdna_contained, dtype=np.float64)
    region_eff = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9)
    rho_ref = _global_reference_density(contained, calibration.gdna_region_eff_len)
    if rho_ref is None or rho_ref <= 0.0:
        elen = support_len.copy()
    else:
        inv = 1.0 / rho_ref
        elen = np.minimum(contained * inv, region_eff) + np.minimum(pooled * inv, seam_len)

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

    proj = _project_regions_to_loci(
        region_arrays,
        multi_loci,
        len(multi_loci),
        {
            "gdna": gdna_region,
            "rna": rna_region,
            "span": support_len,  # Σ S — the EFFECTIVE support (region_eff_len + summed seams), NOT genomic
            # the CONTAINED (unique-mapper) mass per locus — the calibration-blindness discriminator for the
            # eff-len guard below (calibration's accumulator is fed by unique mappers only).
            "gdna_contained": np.asarray(calibration.mass_gdna_contained, dtype=np.float64),
            "rna_contained": np.asarray(calibration.mass_rna_contained, dtype=np.float64),
            # per-region enrichment-weighted node length (global ρ_ref) → the gDNA component's eff-length.
            "elen": elen,
        },
    )
    gdna_locus, span = proj["gdna"], proj["span"]
    # gDNA EM effective length = the enrichment-weighted length of the locus's gDNA against the SHARED global
    # ρ_ref: eff = Σ_locus share·min(m_n/ρ_ref, S_n) = proj["elen"]. gDNA's node set is ALL the locus's nodes
    # (every region + boundary over its span; the intergenic regions outside are dropped by the projection).
    # Under uniform gDNA every node is at ρ_ref ⇒ elen = support ⇒ eff = span (no contraction). Using the SAME
    # ρ_ref for gDNA AND every transcript puts the gDNA-vs-transcript density comparison on ONE consistent
    # scale (and gives eff(nascent) ≥ eff(mature)). Replaces the former per-locus ρ* = G_c/E_c.
    #
    # CONTAINED-EVIDENCE SHRINKAGE (calibration multimapper-blindness): the accumulator is fed by UNIQUE
    # mappers only, so a multimapper-dominated locus (identical paralogs) has little CONTAINED mass and an
    # unreliable reference read. Shrink the contracted length toward the uniform span on the reliable
    # contained evidence C, smoothly (w = C/(C+1), one pseudo-observation, magic-free). See project_mappability.
    contained_ev = np.maximum(proj["gdna_contained"] + proj["rna_contained"], 0.0)
    w = contained_ev / (contained_ev + 1.0)
    eff_len = w * proj["elen"] + (1.0 - w) * span

    return LocusPriors(
        gdna_prior_count=gdna_locus,
        rna_prior_count=proj["rna"],
        # Clamp into [min(floor, span), span]: the 1 bp floor (matching the EM's own eff-len floor) applies
        # to every real locus, but must never exceed the locus's own effective span — otherwise a degenerate
        # sub-basepair span (e.g. a microexon-only locus, region shorter than a fragment ⇒ E[max(0,L−ℓ)]≈0)
        # would return eff_len > span, breaking eff_len ∈ (0, span]. No effect on real loci (span ≫ 1).
        gdna_eff_len=np.minimum(np.maximum(eff_len, _GDNA_EFF_LEN_FLOOR), np.maximum(span, 1e-9)),
    )


__all__ = ["LocusPriors", "assemble_priors"]
