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


def _gdna_node_arrays(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-node arrays ``(gdna_total, support_total, edge_mass, edge_support)`` for the gDNA node set —
    the SHARED node model with ``transcript_capture_eff_lengths``.

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

    Mass conservation (no mass moved — transport-free): ``Σ gdna_total = Σ m_r + Σ m_e``, every object
    counted exactly once.
    """
    node_mass = np.asarray(calibration.mass_gdna_node, dtype=np.float64)
    node_support = np.maximum(np.asarray(calibration.gdna_node_eff_len, dtype=np.float64), 1e-9)
    n = node_mass.shape[0]

    edge_mass = np.zeros(n, dtype=np.float64)
    edge_support = np.zeros(n, dtype=np.float64)
    if calibration.n_edges:
        owner = edge_owner_nodes(calibration, region_arrays)
        # a node may own its own right line AND a re-keyed one from its left, hence np.add.at
        np.add.at(edge_mass, owner, np.asarray(calibration.mass_gdna_edge, dtype=np.float64))
        np.add.at(edge_support, owner, np.asarray(calibration.gdna_edge_eff_len, dtype=np.float64))

    return node_mass + edge_mass, node_support + edge_support, edge_mass, edge_support


def assemble_priors(
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    multi_loci: "list[MultiLocus]",
) -> LocusPriors:
    """Build the per-locus EM prior from the calibration result.

    The gDNA object set is the **density-correct, transport-free** node + line model
    (see :func:`_gdna_node_arrays`):

        node r:  mass = mass_gdna_node[r],  effective support S_r = E_f[(L_r − w + 1)+]
        edge e:  mass = mass_gdna_edge[e],  effective support S_e = E_f[w − 1]

    with each line attributed to a flank node (:func:`edge_owner_nodes`). Masses and supports project
    to loci by genomic-overlap ``share``::

        gdna_prior_count = Σ share * (node + edge)                            (deconvolved gDNA count)
        rna_prior_count  = Σ share * rna_node_total   (UNSPLICED RNA; spliced withheld — see below)
        gdna_eff_len     = (G+1)² / [ Σ share*(node²/S_r + edge²/S_e) + (2G+1)/span ], capped at span
                           G = Σ share*(node+edge),  span = Σ share*(S_r + S_e)  (EFFECTIVE support)

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

    # Density-correct, transport-free gDNA object model: the per-node CONTAINED object (support
    # gdna_node_eff_len) + the per-line CROSSING object (support gdna_edge_eff_len), each attributed to
    # a flank node — the SAME object model transcript_capture_eff_lengths contracts on (EFFECTIVE, not
    # genomic, supports; the factor-1-under-uniform bedrock).
    gdna_region, support_len, pooled, seam_len = _gdna_node_arrays(calibration, region_arrays)

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
    rna_region = np.asarray(calibration.mass_rna_node, dtype=np.float64).copy()
    if calibration.n_edges:
        np.add.at(rna_region, edge_owner_nodes(calibration, region_arrays), rna_edge_unspliced)
    rna_region = np.maximum(rna_region, 0.0)

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
            "gdna_contained": np.asarray(calibration.mass_gdna_node, dtype=np.float64),
            "rna_contained": np.asarray(calibration.mass_rna_node, dtype=np.float64),
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


__all__ = ["LocusPriors", "assemble_priors", "edge_owner_nodes"]
