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
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-region (mass, IPR participation, effective support) for the gDNA region + pooled-seam nodes.

    The density-correct, transport-free gDNA node model (``effective_length_redesign_plan.md`` §8 —
    the first-principles rev). The node set over a region chain is the per-region CONTAINED node plus
    one POOLED SEAM node per internal boundary (genomically adjacent, same reference), keyed to its
    left-flank region ``r``::

        region node r:     mass m_r = mass_gdna_contained[r],
                           effective support S_r = gdna_region_eff_len[r] = E[max(0, L_r − ℓ)]
        seam node (r,r+1): mass m_s = mass_gdna_right[r] + mass_gdna_left[r+1]   (the two halves POOLED),
                           effective support S_s = ½·(E[min(ℓ,L_r)] + E[min(ℓ,L_{r+1})])  (the AVERAGE of
                           the two flanking per-side density lengths gdna_boundary_len)

    **Why these supports — the bedrock invariant.** Driving the reference accumulator under uniform
    genomic gDNA at density ρ, the expected masses are ``m_r = ρ·E[max(0,L_r−ℓ)]`` (a contained fragment
    must FIT) and ``m_s = ρ·(E[min(ℓ,L_r)] + E[min(ℓ,L_{r+1})])/2`` (each side captures only ``min(ℓ,L_side)``
    of a crossing fragment, so the pooled mass is ρ times the AVERAGE of the two side density lengths).
    Dividing each node's mass by these supports makes EVERY node density ``m_n/S_n = ρ``, and the
    Laplace-smoothed IPR ``(G+1)²/[Σ m_n²/S_n + (2G+1)/ΣS_n]`` then returns ``ΣS_n`` EXACTLY (contraction
    factor 1) for any ρ — an unenriched library contracts nothing, EVEN for regions shorter than ``E[ℓ]``.
    Two divisors are WRONG: the genomic ``region_size_bp`` understates short-region density (verified
    factor 0.878 under uniform); and the count crossing length ``E[ℓ]`` over-states the seam support for
    short flanks (a fragment can only deposit ``min(ℓ,L_side)``, not ``ℓ``), under-contracting exon-flank
    seams and inflating the gДНК→RNA leak — the averaged side density length is the deposition-faithful
    support.

    **Why POOL the seam, not split it.** Entering the two halves ``s_L,s_R`` as separate nodes contributes
    ``s_L²/S_s + s_R²/S_s ≤ (s_L+s_R)²/S_s`` (Cauchy–Schwarz) — splitting at most halves the IPR support
    (over-contracting). The halves are also one physical crossing event, so the pooled node is faithful.

    Returns ``(gdna_region, participation, support_len)``, each float64[R] keyed to region ``r``::

        gdna_region[r]   = m_r + m_s(r,r+1)        total node mass on region r
        participation[r] = m_r²/S_r + m_s²/S_s      the Σ m²/S terms
        support_len[r]   = S_r + S_s(r,r+1)         the Σ S effective-support terms

    Mass conservation (no mass moved — transport-free): ``Σ gdna_region = Σ contained +
    Σ_{internal} (right[r] + left[r+1])`` — every non-terminal boundary side counted exactly once
    (terminal / cross-reference sides carry zero on real data and are excluded).
    """
    contained = np.asarray(calibration.mass_gdna_contained, dtype=np.float64)
    side_left = np.asarray(calibration.mass_gdna_left, dtype=np.float64)
    side_right = np.asarray(calibration.mass_gdna_right, dtype=np.float64)
    region_eff_len = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9)
    # per-region per-side density length E[min(ℓ,L_r)] (the mass a crossing fragment deposits on region r)
    side_density_len = np.asarray(calibration.gdna_boundary_len, dtype=np.float64)
    ref_id = np.asarray(region_arrays.ref_id)
    n = contained.shape[0]

    pooled = np.zeros(n, dtype=np.float64)  # pooled seam mass, attributed to a flank region
    seam_len = np.zeros(n, dtype=np.float64)  # seam effective support, attributed to a flank region
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]  # internal seam: genomically adjacent, same reference
        seam_mass = np.where(same, side_right[:-1] + side_left[1:], 0.0)  # POOL the boundary's two halves
        # seam support = AVERAGE of the two flanking per-side density lengths E[min(ℓ,L)] — the
        # deposition-faithful divisor (each side captures min(ℓ,L_side) of a crossing fragment, so the
        # pooled mass is ρ·(E[min_left]+E[min_right])/2 under uniform gDNA). NOT E[ℓ], which over-states
        # the support for short flanks and under-contracts.
        seam_support = np.where(same, 0.5 * (side_density_len[:-1] + side_density_len[1:]), 0.0)
        # ATTRIBUTE each seam to a flank REGION so the locus projection picks it up. Default: the LEFT
        # flank r. BUT a locus's far-LEFT outer boundary is an intergenic→(exon/intron) seam whose left
        # flank is INTERGENIC — a region that overlaps no locus and is dropped by _project_regions_to_loci.
        # Keying to the left flank there SILENTLY LOSES that boundary's pure intergenic-crossing gDNA (both
        # the intergenic side_right and the first-region side_left), under-counting the locus gDNA prior AND
        # inflating the gDNA-component IPR eff-length (a high-density crossing node vanishes from Σm²/S, so
        # eff_len lengthens). The far-RIGHT boundary is already kept (its left flank is the locus's last
        # region), so this restores the symmetry: attribute the seam to the RIGHT flank whenever the left
        # flank is intergenic (no RNA-signature bits) and the right flank is not — otherwise keep the left.
        sig = np.asarray(region_arrays.signature).astype(np.int64)
        ig = (sig & _RNA_SIGNATURE_BITS) == 0  # intergenic: no exon/intron bit ⇒ dropped by the projection
        rekey_right = same & ig[:-1] & ~ig[1:]  # far-left outer boundary: intergenic → locus region
        owner = np.where(rekey_right, np.arange(1, n), np.arange(0, n - 1))
        np.add.at(pooled, owner, seam_mass)  # accumulate: a first-region node may own its right seam + this
        np.add.at(seam_len, owner, seam_support)

    gdna_region = contained + pooled
    safe_seam = np.maximum(
        seam_len, 1e-9
    )  # floor the divisor consistently with the seam_len>0 gate
    with np.errstate(divide="ignore", invalid="ignore"):
        participation = contained**2 / region_eff_len + np.where(
            seam_len > 0.0, pooled**2 / safe_seam, 0.0
        )
    support_len = region_eff_len + seam_len
    return gdna_region, participation, support_len


def assemble_priors(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
) -> LocusPriors:
    """Build the per-locus EM prior from the calibration result.

    The gDNA node set is the **density-correct, transport-free** region + pooled-seam model
    (``effective_length_redesign_plan.md`` §8; see ``_gdna_region_node_arrays``):

        region node r:    mass = mass_gdna_contained[r],   effective support S_r = E[max(0, L_r − ℓ)]
        seam node (r,r+1): mass = mass_gdna_right[r] + mass_gdna_left[r+1] (the two halves POOLED),
                           effective support S_s = ½·(E[min(ℓ,L_r)] + E[min(ℓ,L_{r+1})])  (the average of
                           the flanking per-side density lengths gdna_boundary_len)

    keyed to the left-flank region r. The region + seam masses / participations project to loci by
    genomic-overlap ``share``::

        gdna_prior_count = Σ share * (contained + seam)                       (deconvolved gDNA count)
        rna_prior_count  = Σ share * rna_region         (UNSPLICED RNA; spliced withheld — see below)
        gdna_eff_len     = (G+1)² / [ Σ share*(contained²/S_r + seam²/S_s) + (2G+1)/span ], capped at span
                           G = Σ share*(contained+seam),  span = Σ share*(S_r + S_s)  (EFFECTIVE support)

    **The bedrock invariant — factor 1 under uniform gDNA.** Dividing each node's mass by its EFFECTIVE
    sampling support (``S_r`` for contained, ``S_s = ½·(E[min(ℓ,L_r)]+E[min(ℓ,L_{r+1})])`` for the pooled
    crossing — the AVERAGE per-side density length, NOT ``E[ℓ]`` which over-states it) makes the per-node
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

    # Density-correct, transport-free gDNA node model (effective_length_redesign_plan.md §8): per-region
    # CONTAINED node (effective support gdna_region_eff_len = E[max(0,L−ℓ)]) + one POOLED SEAM node per
    # internal boundary (effective support = ½ the sum of the flanking E[min(ℓ,L)] = averaged gdna_boundary_len),
    # keyed to the left-flank region. The participation Σ m²/S and the span Σ S use these EFFECTIVE supports — NOT
    # region_size_bp, which understates short-region density and fabricates a contraction under uniform
    # gDNA. _gdna_region_node_arrays carries the bedrock factor-1-under-uniform proof. NOTE: the same ρ_ref
    # is shared with the transcript contraction (capture_eff_length), but the two apply it slightly
    # differently — here the contraction min() runs over the FOLDED region+seam node mass/support, whereas
    # capture_eff_length applies min() PER node (contained node and seam node separately). Since
    # min(a+b, Sa+Sb) ≥ min(a,Sa)+min(b,Sb), the gDNA component is marginally LESS contracted than an
    # equivalent transcript over the same nodes. Both are factor-1 under uniform gDNA (no divergence there);
    # the difference appears only where a node is enriched and its seam depleted. TODO(cleanup): unify via a
    # single shared per-node contraction helper so the two paths are byte-identical (docs/calibration/
    # siphon_measurement.md notes this as a follow-up).
    gdna_region, gdna_sq_over_len, support_len = _gdna_region_node_arrays(
        calibration, region_arrays
    )

    # SHARED global reference density — the SAME ρ_ref the transcript contraction uses. gDNA and EVERY
    # transcript must contract against ONE reference, else the gDNA-vs-transcript density comparison sits on
    # inconsistent per-locus scales. Per-region enrichment-weighted node length elen[r] = min(gdna_region[r]/
    # ρ_ref, support_len[r]) (the node's support that survives at the fully-captured density); ρ_ref None (no
    # gDNA) ⇒ elen = support (no contraction). See docs/calibration/enriched_mode_reference_density.md.
    from .capture_eff_length import _global_reference_density

    rho_ref = _global_reference_density(
        calibration.mass_gdna_contained, calibration.gdna_region_eff_len
    )
    if rho_ref is None or rho_ref <= 0.0:
        elen = support_len.copy()
    else:
        elen = np.minimum(gdna_region / rho_ref, support_len)

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
            "span": support_len,  # Σ S — the EFFECTIVE support (region_eff_len + averaged seams), NOT genomic
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
        gdna_eff_len=np.minimum(
            np.maximum(eff_len, _GDNA_EFF_LEN_FLOOR), np.maximum(span, 1e-9)
        ),
    )


__all__ = ["LocusPriors", "assemble_priors"]
