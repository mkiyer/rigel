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

from .region_arrays import edge_node_indices
from .signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS

if TYPE_CHECKING:
    from ..locus import MultiLocus
    from .region_arrays import RegionArrays
    from .result import CalibrationResult

# A node with none of these strand/type bits is intergenic — it overlaps no locus and is dropped by
# the per-locus projection, so a line whose left flank is such a node must be re-attributed to its
# (locus) right flank or its gDNA is lost (see edge_owner_nodes).
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


def edge_owner_nodes(calibration: "CalibrationResult", region_arrays: "RegionArrays") -> np.ndarray:
    """``int64[E]`` — which NODE each contiguous edge is attributed to for the locus projection.

    An edge is a 0-bp line, so it has no genomic extent to overlap a locus with; it must be carried by
    one of its two flanking nodes, which do. Default: the **left** flank ``lo``.

    ⚠ **Except at a locus's far-LEFT outer line.** That is an intergenic→(exon/intron) seam whose left
    flank is INTERGENIC — a node that overlaps no locus and is dropped by
    :func:`_project_regions_to_loci`. Keying to the left flank there SILENTLY LOSES that line's crossing
    gDNA, under-counting the locus gDNA prior AND inflating its eff-length. The far-RIGHT line is
    already kept (its left flank is the locus's last node), so attributing to the RIGHT flank whenever
    the left is intergenic and the right is not restores the symmetry.
    """
    lo, hi = edge_node_indices(np.asarray(region_arrays.ref_id))
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    ig = (
        sig & _RNA_SIGNATURE_BITS
    ) == 0  # intergenic: no exon/intron bit ⇒ dropped by the projection
    return np.where(ig[lo] & ~ig[hi], hi, lo)


def _component_node_arrays(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    mass_node: np.ndarray,
    mass_edge: np.ndarray,
    eff_node: np.ndarray,
    eff_edge: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-node arrays ``(mass_total, support_total, edge_mass, edge_support)`` for ONE component —
    the SHARED node model with ``transcript_capture_eff_lengths``.

    ⭐ **Parameterised by component (P1).** It was ``_gdna_node_arrays`` and reached into the gDNA
    fields directly. ``assemble_priors`` now needs the identical fold for RNA — because converting a
    mass into a density requires THAT COMPONENT'S OWN opportunity as the divisor — and one function
    called twice is the alternative to a second implementation that drifts (`CARRY_FORWARD.md` §3
    trap 27).

    The density-correct, transport-free gDNA node model. The object set over a chain is the per-node
    CONTAINED object plus the per-line CROSSING object, each attributed to a flank node
    (:func:`edge_owner_nodes`) so the locus projection can pick it up::

        node r:  mass m_r = mass_gdna_node[r],  effective support S_r = gdna_node_eff_len[r]
        edge e:  mass m_e = mass_gdna_edge[e],  effective support S_e = gdna_edge_eff_len[e]

    ⭐ **THE POOLING IS GONE, NOT RE-DERIVED (S5.f).** This function used to call
    ``capture_eff_length._pooled_seam_arrays`` to add ``mass_gdna_right[r] + mass_gdna_left[r+1]`` back
    together and sum the two halved per-side lengths ``gdna_boundary_len[r] + gdna_boundary_len[r+1]``.
    The calibrator had split one crossing into two faces and this put it back — a no-op with a history,
    since that exact sum-then-halve pattern hid an *exact factor of 2* for months
    (`CARRY_FORWARD.md` §3 trap 2). ``chain_edge_deconv`` now returns one number per line, so there is
    nothing to pool, no ½ anywhere, and no second implementation to keep in step.

    **Why these supports — the bedrock invariant.** Under uniform genomic gDNA at density ρ the
    accumulator deposits ``m_r = ρ·E_f[(L_r − w + 1)+]`` on a node (a contained fragment must FIT) and
    ``m_e = ρ·E_f[w − 1]`` on a line (a crossing fragment has ``w − 1`` admissible offsets). Those are
    exactly ``gdna_node_eff_len`` and ``gdna_edge_eff_len``, so EVERY object's density ``m/S`` is ρ, the
    enrichment contraction ``min(m/ρ_ref, S)`` returns ``S`` exactly, and an unenriched library
    contracts NOTHING — even for nodes shorter than one fragment.

    ⛔ The genomic ``region_size_bp`` is the wrong node divisor: it ignores the fit-inside constraint,
    so it understates a short node's density and fabricates a contraction with no capture bias present
    (verified factor 0.878 under uniform).

    Returns four ``float64[N]`` keyed to node ``r``::

        gdna_total[r]    = m_r + Σ_{e owned by r} m_e        total mass carried by node r
        support_total[r] = S_r + Σ_{e owned by r} S_e        total effective support (Σ S)
        edge_mass[r]     = Σ_{e owned by r} m_e              the crossing mass ALONE
        edge_support[r]  = Σ_{e owned by r} S_e              the crossing support ALONE

    ``assemble_priors`` contracts PER OBJECT — ``min(m_r/ρ_ref, S_r) + min(edge_mass/ρ_ref, edge_support)``
    — matching ``transcript_capture_eff_lengths``, NOT over the folded ``gdna_total`` (which would
    under-contract a captured exon whose line runs into a depleted intron).

    ⚠ **OBJECT conservation, not FRAGMENT conservation** — ``Σ mass_total = Σ m_r + Σ m_e``, every
    OBJECT counted exactly once. That is not the same as every fragment counted once: a fragment
    deposits on ``max(K, 1)`` objects. Converting to fragments is ``assemble_priors``' job and is what
    the supports returned here are for.
    """
    node_support = np.maximum(np.asarray(eff_node, dtype=np.float64), 0.0)
    node_mass = _mass_where_there_is_opportunity(mass_node, node_support)
    n = node_mass.shape[0]

    edge_mass = np.zeros(n, dtype=np.float64)
    edge_support = np.zeros(n, dtype=np.float64)
    if calibration.n_edges:
        owner = edge_owner_nodes(calibration, region_arrays)
        e_support = np.maximum(np.asarray(eff_edge, dtype=np.float64), 0.0)
        # a node may own its own right line AND a re-keyed one from its left, hence np.add.at
        np.add.at(edge_mass, owner, _mass_where_there_is_opportunity(mass_edge, e_support))
        np.add.at(edge_support, owner, e_support)

    return node_mass + edge_mass, node_support + edge_support, edge_mass, edge_support


def _mass_where_there_is_opportunity(mass: np.ndarray, support: np.ndarray) -> np.ndarray:
    """Drop a component's mass on objects where that component had **zero opportunity**.

    ⛔ **Both sides of the pooled ratio, or neither.** ``rho = Σm / ΣS`` is a Poisson rate — observations
    over exposure. An object with zero exposure that nonetheless carries mass is a contradiction in the
    model, and leaving its mass in the numerator inflates ``rho`` with nothing in the denominator to pay
    for it. It must contribute to neither: it carries no information about that component's density.

    ⚠ **This is not a hypothetical.** ``contained_eff_length`` is exactly 0 wherever a node is shorter
    than the shortest fragment in that component's pmf — measured on the chr22 pilot against its own
    measured pure pools, **21.7 % of nodes for RNA and 18.7 % for gDNA** (`CARRY_FORWARD.md` §1 fact 9
    records 12.4 % genome-wide). The solver can still put mass there, because ``f_g`` is an inference
    and not a fact, so ``mass > 0`` with ``support == 0`` is an ordinary configuration.

    ⭐ Found by perturbation P1e: flooring the divisor to ``1e-9`` instead of testing ``support > 0``
    left every test green, because the only zero-support fixture also had zero mass. The floor would
    have turned a single stray fragment on a 40 bp node into a density of ``1e9`` — `CARRY_FORWARD.md`
    §3 trap 23, which is how a "no data" default of 100 % gDNA once seeded false gDNA into neighbouring
    exons. The same contract as ``node_geometry._rate`` on the solver side.
    """
    m = np.asarray(mass, dtype=np.float64)
    return np.where(np.asarray(support, dtype=np.float64) > 0.0, m, 0.0)


def _density_times_span(mass: np.ndarray, support: np.ndarray, span_bp: np.ndarray) -> np.ndarray:
    """``(Σ mass / Σ support) · span_bp`` — a per-object mass total converted to a FRAGMENT COUNT.

    The three arrays are already per-locus sums. ⛔ Zero support returns **0**, never a floored
    division: an object with no opportunity for a component must emit nothing (`CARRY_FORWARD.md` §3
    trap 23). Negative mass cannot occur (``CalibrationResult`` rejects it) but is clamped anyway,
    because a prior pseudocount below zero is not a prior.
    """
    m = np.asarray(mass, dtype=np.float64)
    s = np.asarray(support, dtype=np.float64)
    rho = np.divide(m, s, out=np.zeros_like(m), where=s > 0.0)
    return np.maximum(rho * np.asarray(span_bp, dtype=np.float64), 0.0)


def assemble_priors(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
) -> LocusPriors:
    """Build the per-locus EM prior from the calibration result.

    Each component's object set is the **density-correct, transport-free** node + line model
    (see :func:`_component_node_arrays`), on that component's OWN opportunity:

        node r:  mass = mass_c_node[r],  effective support S_r = E_c[(L_r − w + 1)+]
        edge e:  mass = mass_c_edge[e],  effective support S_e = E_c[w − 1]

    with each line attributed to a flank node (:func:`edge_owner_nodes`). Masses and supports project
    to loci by genomic-overlap ``share``::

        rho_c            = Σ share*mass_c / Σ share*S_c            <- pooled density, RATIO OF SUMS
        {gdna,rna}_prior_count = rho_c * span_bp                   <- the SAME genomic span for both
        gdna_eff_len     = (G+1)² / [ Σ share*(node²/S_r + edge²/S_e) + (2G+1)/span ], capped at span
                           G = Σ share*(node+edge),  span = Σ share*(S_r + S_e)  (EFFECTIVE support)

    ⭐ **THE PRIOR IS A FRAGMENT COUNT, AND A SUM OF PER-OBJECT MASSES IS NOT ONE.** The EM adds these
    scalars straight to its own fragment counts (``G = n_gdna + a_g``,
    ``em_solver.cpp:apply_grouped_prior_update``), so they must be in fragment units. But the accumulator
    deposits ``+1`` on **every** object a fragment touches — ``max(K, 1)`` of them, ``K`` being the lines
    crossed — so summing per-object masses gives an object-incidence count::

        incidences(w) = max( 1 , (w-1)/s )        for a partition of spacing s

    Counts are conserved exactly where every node is longer than every fragment, and become a
    **length-weighted** count where they are not — and 56.7 % of human nodes are shorter than one 200 bp
    fragment. Because the weight is the fragment's own length it does NOT cancel between two components
    with different mean lengths: measured on the chr22 pilot, gDNA deposits **1.031** incidences per
    fragment against RNA **≈1.17**, so the raw sum under-called gDNA by 13–19 %. Deterministically:
    re-tiling 1200 bp from one node to twelve moved the RNA prior **2.19×** with the library unchanged
    (`tests/calibration/test_prior_units.py`).

    Density is intensive — it is pooled as a ratio of sums and then **integrated over the span**. That is
    the one-line statement of the model, and it is partition-free because ``mass_c/S_c`` is.

    ⚠ **The pooling is support-weighted, and that is an approximation, not an identity.**
    ``Σm/ΣS`` is the support-weighted mean density, so ``rho·span`` is exact only where ``rho_c`` is
    uniform *within* the locus. It is the same pooling `derive.gdna_density_global` and
    `CARRY_FORWARD.md` §2's ``rho_bg = Σg/ΣE`` already use — a ratio of sums, never a mean of ratios —
    and it is a strict improvement on the raw sum, but a locus with a strong internal density gradient
    carries a second-order residual. Under uniform ρ the two components' weightings agree exactly, so
    the g:r **ratio** is exact there regardless.

    ⚠ **``span_bp`` is the GENOMIC span and is the same number for both components.** Any
    component-specific span (``ΣS_c``, say) would re-introduce exactly the tilt this removes — ``ρ_c·ΣS_c``
    *is* the raw sum. Edges contribute mass and support but **zero** ``span_bp``: a contiguous edge is a
    0-bp line, which is correct and not an omission.

    **The bedrock invariant — factor 1 under uniform gDNA.** Dividing each object's mass by its
    EFFECTIVE sampling support makes its density ``m/S`` exactly the true ρ under a uniform (unenriched)
    library, because the accumulator deposits ``ρ·E_f[(L−w+1)+]`` of contained mass on a node and
    ``ρ·E_f[w−1]`` of crossing mass on a line. The Laplace-smoothed IPR then returns ``span`` EXACTLY
    (``eff_len = span`` ⇒ contraction factor 1): an unenriched library contracts NOTHING. Using the
    genomic ``region_size_bp`` instead would understate short-node density and fabricate a contraction
    even with no capture bias (verified factor 0.878 vs the correct 1.000) — that was the latent defect
    this design removes. Under capture (concentrated gDNA) the IPR contracts below ``span`` toward the
    probed footprint, so the gDNA component competes at its true local density. This is
    **transport-free**: no mass is moved (no boundary-flux redistribution — that non-physical heuristic
    is gone), and the line object quarantines the captured intron↔exon crossing mass at crossing
    density, recovering the gDNA concentration a node-dilution divisor would lose.

    ⭐ **There is no seam POOLING any more.** The predecessor split each crossing into two faces and
    then summed them back here; ``chain_edge_deconv`` returns one number per line, so the split, the
    re-pool and the ``gdna_boundary_len`` ½ that travelled with them are all gone (`S5_DESIGN_LOG.md`
    §4; `CARRY_FORWARD.md` §3 trap 2).

    **Laplace-smoothed** by one fragment-equivalent of uniform support (the ``(2G+1)/span`` term), the
    canonical add-one prior with no tunable constant: ``G = 0`` ⇒ ``span`` exactly; abundant gDNA
    (capture, ``G ≫ 1``) ⇒ the empirical concentration; in between, the IPR is shrunk toward the
    uniform ``span`` in proportion to the gDNA evidence, so the EM cannot amplify a tiny concentrated
    mass past the calibration's call.
    """
    if calibration.n_nodes != region_arrays.n_regions:
        raise ValueError(
            f"calibration has {calibration.n_nodes} nodes but region_arrays has "
            f"{region_arrays.n_regions}; they must address the same partition."
        )

    # Density-correct, transport-free object model, PER COMPONENT: the per-node CONTAINED object +
    # the per-line CROSSING object, each attributed to a flank node — the SAME object model
    # transcript_capture_eff_lengths contracts on (EFFECTIVE, not genomic, supports; the
    # factor-1-under-uniform bedrock). ⭐ Called twice: converting a mass to a density needs THAT
    # COMPONENT'S own opportunity, and using one divisor for both is the tilt P1 removes.
    gdna_region, support_len, pooled, seam_len = _component_node_arrays(
        calibration,
        region_arrays,
        calibration.mass_gdna_node,
        calibration.mass_gdna_edge,
        calibration.gdna_node_eff_len,
        calibration.gdna_edge_eff_len,
    )

    # SHARED global reference density — the SAME ρ_ref every transcript contracts against, so the
    # gDNA-vs-transcript density comparison sits on ONE scale. The enrichment contraction is applied PER
    # NODE (contained region node + pooled-seam node separately), identical to transcript_capture_eff_lengths:
    #   elen = min(contained/ρ_ref, S_region) + min(pooled/ρ_ref, S_seam).
    # Folding the two into one min() would UNDER-contract a captured exon whose seam runs into a depleted
    # intron (up to ~13% per region under capture, verified). ρ_ref None (no detectable gDNA) ⇒ elen =
    # support (no contraction). See docs/CARRY_FORWARD.md.
    from .capture_eff_length import _global_reference_density

    contained = np.asarray(calibration.mass_gdna_node, dtype=np.float64)
    region_eff = np.maximum(np.asarray(calibration.gdna_node_eff_len, dtype=np.float64), 1e-9)
    rho_ref = _global_reference_density(contained, calibration.gdna_node_eff_len)
    if rho_ref is None or rho_ref <= 0.0:
        elen = support_len.copy()
    else:
        inv = 1.0 / rho_ref
        elen = np.minimum(contained * inv, region_eff) + np.minimum(pooled * inv, seam_len)

    # RNA prior = the UNSPLICED RNA mass only. Spliced fragments have no gDNA candidate in the EM
    # (gDNA does not splice) → they are guaranteed-RNA and the EM assigns them directly; counting
    # them in rna_prior_count would double-count them and unfairly inflate the RNA side of the
    # gDNA-vs-RNA unspliced split (the prior arbitrates only the unspliced fragments, so a_g+a_r
    # should equal the unspliced competing mass). mass_rna_edge is spliced-inclusive (per-edge
    # conservation); subtracting mass_rna_spliced_edge here leaves (1−g)·unspliced on every line.
    #
    # ⚠ The JUNCTION flux (``mass_rna_junction``) is deliberately NOT added. It is certified RNA in
    # exactly the sense the spliced crossings are withheld for, so feeding it in would load the RNA
    # side of a split that arbitrates only unspliced fragments — and a locus whose RNA is fully spliced
    # SHOULD have a near-zero rna_prior_count, because its unspliced fragments really are gDNA or
    # nascent. The result exports the flux for QC; the prior does not read it (owner ruling, 2026-07-30).
    rna_edge_unspliced = np.maximum(
        np.asarray(calibration.mass_rna_edge, dtype=np.float64)
        - np.asarray(calibration.mass_rna_spliced_edge, dtype=np.float64),
        0.0,
    )
    rna_region, rna_support, _, _ = _component_node_arrays(
        calibration,
        region_arrays,
        calibration.mass_rna_node,
        rna_edge_unspliced,
        calibration.rna_node_eff_len,
        calibration.rna_edge_eff_len,
    )
    rna_region = np.maximum(rna_region, 0.0)

    proj = _project_regions_to_loci(
        region_arrays,
        multi_loci,
        len(multi_loci),
        {
            "gdna": gdna_region,
            "rna": rna_region,
            # the two per-component OPPORTUNITY totals — the divisors that turn a mass into a density
            "gdna_support": support_len,
            "rna_support": rna_support,
            # ⭐ the GENOMIC span: the number of start positions, the SAME for both components. Nodes
            # only — an edge is a 0-bp line and contributes none.
            "span_bp": np.asarray(region_arrays.region_size_bp, dtype=np.float64),
            "span": support_len,  # Σ S — the EFFECTIVE support (region_eff_len + summed seams), NOT genomic
            # the CONTAINED (unique-mapper) mass per locus — the calibration-blindness discriminator for the
            # eff-len guard below (calibration's accumulator is fed by unique mappers only).
            "gdna_contained": np.asarray(calibration.mass_gdna_node, dtype=np.float64),
            "rna_contained": np.asarray(calibration.mass_rna_node, dtype=np.float64),
            # per-region enrichment-weighted node length (global ρ_ref) → the gDNA component's eff-length.
            "elen": elen,
        },
    )
    span = proj["span"]

    # ⭐ MASS -> DENSITY -> FRAGMENT COUNT. The pooled density is a ratio of sums (never a mean of
    # ratios, `CARRY_FORWARD.md` §2), integrated over the locus's GENOMIC span — the same span for both
    # components, so the ratio a_g:a_r is ρ_g:ρ_r and carries no length tilt.
    # ⛔ `where=` and not a floored divisor: a locus with no opportunity for a component must emit
    # NOTHING, never a floored division (`CARRY_FORWARD.md` §3 trap 23 — a "no data" default of 100 %
    # gDNA was actively seeding false gDNA into neighbouring exons).
    span_bp = proj["span_bp"]
    gdna_locus = _density_times_span(proj["gdna"], proj["gdna_support"], span_bp)
    rna_locus = _density_times_span(proj["rna"], proj["rna_support"], span_bp)
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
        rna_prior_count=rna_locus,
        # Clamp into [min(floor, span), span]: the 1 bp floor (matching the EM's own eff-len floor) applies
        # to every real locus, but must never exceed the locus's own effective span — otherwise a degenerate
        # sub-basepair span (e.g. a microexon-only locus, region shorter than a fragment ⇒ E[max(0,L−ℓ)]≈0)
        # would return eff_len > span, breaking eff_len ∈ (0, span]. No effect on real loci (span ≫ 1).
        gdna_eff_len=np.minimum(np.maximum(eff_len, _GDNA_EFF_LEN_FLOOR), np.maximum(span, 1e-9)),
    )


__all__ = ["LocusPriors", "assemble_priors", "edge_owner_nodes"]
