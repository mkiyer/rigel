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

import os
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


# gDNA eff-len reference density ρ* (eff = θ_g/ρ*). Two production-supported modes, env-selectable:
#   "contained" (DEFAULT, production): the mass-weighted CONTAINED (exon) density G_c/E_c — robust across
#      capture on/off, stranded/unstranded, and gDNA level on today's (observed) calibration.
#   "kmeans" (OPTIONAL, the TARGET): the magic-free per-locus ENRICHED-MODE reference (see below). Principled
#      and preferred, but it needs a cleanly BIMODAL node-density landscape: on OBSERVED calibration it is
#      currently worse than "contained" (it hallucinates a split from noise on unimodal/capture-off loci) and
#      only wins once the calibration's per-node gDNA accuracy improves (perfect-calibration study). Keep it
#      here so we can flip production to it after the calibration-accuracy work. See
#      docs/calibration/effective_length_state_and_roadmap.md.
_RHOSTAR_MODE = os.environ.get("RIGEL_RHOSTAR", "contained")


def _per_locus_kmeans_rhostar(
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
    n_loci: int,
    density: np.ndarray,
    weight: np.ndarray,
) -> np.ndarray:
    """Per-locus ENRICHED-MODE reference density via a MAGIC-FREE log-space support-weighted 1-sided k-means.

    Collects each locus's per-region gDNA densities ``ρ_n = density[r]`` (weight = support ``S_n``) for
    regions with positive density, then fits two log-space centroids by iterating to a fixed point with the
    split at their geometric midpoint. ``ρ_ref`` = the support-weighted geomean of the ENRICHED (above-split)
    cluster. NO tunable constant (the split is self-determined; k=2 is the depleted-vs-captured physics),
    degenerate-safe (uniform / single node → that density). Robust to the few-enriched/many-depleted panel
    regime where a fixed quantile collapses into the depleted cluster. (enriched-mode-estimator workflow;
    docs/calibration/effective_length_state_and_roadmap.md.)
    """
    per: list[list[tuple[float, float]]] = [[] for _ in range(n_loci)]
    blocks_by_ref: dict[int, list[tuple[int, int, int]]] = {}
    for ml in multi_loci:
        for blk in ml.loci:
            if blk.end > blk.start:
                blocks_by_ref.setdefault(int(blk.ref_id), []).append(
                    (int(blk.start), int(blk.end), int(ml.multi_locus_id))
                )
    for blocks in blocks_by_ref.values():
        blocks.sort()
    starts, ends, ref_off = region_arrays.start, region_arrays.end, region_arrays.ref_offsets
    for ref_id in range(int(region_arrays.n_refs)):
        blocks = blocks_by_ref.get(ref_id)
        if not blocks:
            continue
        lo, hi = int(ref_off[ref_id]), int(ref_off[ref_id + 1])
        bstarts = np.fromiter((b[0] for b in blocks), dtype=np.int64, count=len(blocks))
        for r in range(lo, hi):
            if density[r] <= 0.0 or weight[r] <= 0.0:
                continue
            rs, re = int(starts[r]), int(ends[r])
            for b_start, b_end, lid in blocks[: int(np.searchsorted(bstarts, re, side="left"))]:
                if b_end > rs and min(b_end, re) - max(b_start, rs) > 0:
                    per[lid].append((float(density[r]), float(weight[r])))
    out = np.zeros(n_loci, dtype=np.float64)
    for lid, pairs in enumerate(per):
        if not pairs:
            continue
        d = np.array([p[0] for p in pairs])
        wts = np.array([p[1] for p in pairs])
        x = np.log(d)
        if np.ptp(x) < 1e-12:  # uniform / single node ⇒ that density (degenerate-safe)
            out[lid] = float(d[0])
            continue
        c_hi, c_lo = float(x.max()), float(x.min())
        for _ in range(100):  # iter cap is numerical convergence, not a modeling constant
            tau = 0.5 * (c_hi + c_lo)  # self-determined inter-centroid geometric midpoint (no free knob)
            up = x >= tau
            whi, wlo = wts[up].sum(), wts[~up].sum()
            n_hi = (wts[up] * x[up]).sum() / whi if whi > 0 else c_hi
            n_lo = (wts[~up] * x[~up]).sum() / wlo if wlo > 0 else c_lo
            if abs(n_hi - c_hi) < 1e-12 and abs(n_lo - c_lo) < 1e-12:
                c_hi, c_lo = n_hi, n_lo
                break
            c_hi, c_lo = n_hi, n_lo
        out[lid] = float(np.exp(c_hi))  # support-weighted geomean of the enriched (above-split) cluster
    return out


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

    **The bedrock invariant — factor 1 under uniform gDNA.** Dividing each node's mass by its
    EFFECTIVE sampling support (``S_r`` for contained, ``S_s = E[ℓ]`` for the pooled crossing) makes the
    per-node density ``m_n/S_n`` exactly the true ρ under a uniform (unenriched) library — because the
    accumulator deposits ``ρ·E[max(0,L−ℓ)]`` of contained mass and ``ρ·E[ℓ]`` of pooled crossing mass.
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
    # gDNA. _gdna_region_node_arrays carries the bedrock factor-1-under-uniform proof; capture_eff_length
    # reuses the SAME helper, so the gDNA component and the transcript contraction share one node model.
    gdna_region, gdna_sq_over_len, support_len = _gdna_region_node_arrays(
        calibration, region_arrays
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
            # CONTAINED-footprint IPR pieces (m_c²/S_c and S_c) — for the exon-competition density ρ* below.
            "c_part": np.asarray(calibration.mass_gdna_contained, dtype=np.float64) ** 2
            / np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9),
            "c_span": np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9),
        },
    )
    gdna_locus, span = proj["gdna"], proj["span"]
    # gDNA EM effective length = the component's TOTAL node mass held at the exon-competition DENSITY ρ* —
    # the per-node effective-length view (each node contributes ℓ_n = m_n/ρ*, and eff = Σ ℓ_n = θ/ρ*). gDNA
    # is a contiguous genomic interval, so its node set is ALL the locus's nodes: every region PLUS every
    # boundary over its span — including the two OUTER boundaries that bookend the locus (exclusively-gDNA
    # crossing fragments, no RNA candidate) but NOT the intergenic regions outside it (dropped by the
    # projection). ρ* is read on the CONTAINED (exon) footprint — where the ambiguous fragments mature RNA
    # competes for actually sit (a boundary-crossing fragment is likelihood-resolved to gDNA/nascent, never
    # spliced mature):
    #
    #   E_c = contained-footprint Laplace IPR over the locus's CONTAINED nodes only:
    #           E_c = (G_c+1)² / [ Σ m_c²/S_c + (2G_c+1)/span_c ], capped at span_c   (G_c = Σ m_c)
    #         ⇒ ρ* = G_c/E_c, the mass-weighted exon density where mature competes.
    #   eff_len = θ_g/ρ* = (θ_g/G_c)·E_c, capped at the total span; G_c→0 (multimap-blind) ⇒ span.
    #
    # WHY read ρ* on the exon footprint (vs an all-node mean): the all-node density G/eff is a mean over the
    # whole footprint; under capture the depleted introns/large-mass moderate seams drag it below the exon
    # density, so gDNA is under-weighted where it competes → the gDNA→mature leak (measured: the all-node
    # arithmetic mean and the Lehmer M3/M2 both leave a larger mature leak than this contained read). Under
    # uniform gDNA E_c = span_c and θ_g/G_c = span/span_c ⇒ eff_len = span (factor 1 — capture-off
    # bit-identical). This is the SAME per-node ℓ_n = m_n/ρ* the transcript components use, so gDNA and a
    # full-span nRNA sharing a locus land at the same eff-len up to gDNA's two extra outer-boundary supports
    # (negligible — verified ~2.5% on a single-transcript locus). Derivation: rhostar-derivation workflow.
    Gc = proj["gdna_contained"]
    Pc = proj["c_part"]
    span_c = np.maximum(proj["c_span"], 1e-9)
    with np.errstate(divide="ignore", invalid="ignore"):
        E_c = np.minimum((Gc + 1.0) ** 2 / (Pc + (2.0 * Gc + 1.0) / span_c), span_c)

    # CONTAINED-EVIDENCE SHRINKAGE (calibration multimapper-blindness). The accumulator is fed by UNIQUE
    # mappers only, so a multimapper-dominated locus (identical paralogs) has little CONTAINED mass and an
    # unreliable E_c. Shrink E_c toward its uniform span_c on the reliable contained evidence C, smoothly
    # (w = C/(C+1), one pseudo-observation, magic-free); the G_c→0 guard then sends eff_len→span. See
    # project_mappability.
    contained_ev = np.maximum(proj["gdna_contained"] + proj["rna_contained"], 0.0)
    w = contained_ev / (contained_ev + 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        E_c = np.exp(w * np.log(np.maximum(E_c, 1e-9)) + (1.0 - w) * np.log(span_c))
        ratio = np.where(Gc > 1e-9, gdna_locus / np.maximum(Gc, 1e-9), 1.0)
        eff_len = np.where(Gc > 1e-9, np.minimum(ratio * E_c, span), span)

    if _RHOSTAR_MODE == "kmeans":
        # OPTIONAL (RIGEL_RHOSTAR=kmeans): override the reference density with the magic-free per-locus
        # ENRICHED-MODE (log-space support-weighted 1-sided k-means). eff = θ_g/ρ*, capped at span, then the
        # SAME contained-evidence shrinkage. The TARGET production reference — see the module note + the
        # roadmap doc; on today's observed calibration "contained" (the default above) is more robust.
        supp = np.maximum(np.asarray(calibration.gdna_region_eff_len, dtype=np.float64), 1e-9)
        dens = np.asarray(calibration.mass_gdna_contained, dtype=np.float64) / supp
        rho_star = _per_locus_kmeans_rhostar(region_arrays, multi_loci, len(multi_loci), dens, supp)
        with np.errstate(divide="ignore", invalid="ignore"):
            have = (gdna_locus > 1e-9) & (rho_star > 1e-9)
            eff_raw = np.where(
                have, np.minimum(gdna_locus / np.where(have, rho_star, 1.0), span), span
            )
            eff_len = np.exp(
                w * np.log(np.maximum(eff_raw, 1e-9)) + (1.0 - w) * np.log(np.maximum(span, 1e-9))
            )

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
