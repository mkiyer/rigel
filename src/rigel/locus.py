"""rigel.locus — Locus graph construction and gDNA prior aggregation.

Handles everything between "buffer scan complete" and "per-locus EM
loop": connected-component partitioning of transcripts linked by
shared fragments, and Empirical-Bayes-style gDNA Dirichlet priors
aggregated from calibration regions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from .native import connected_components as _cc_native

from .scored_fragments import Locus, ScoredFragments

from .index import TranscriptIndex

if TYPE_CHECKING:
    from .calibration import CalibrationResult


# ---------------------------------------------------------------------------
# Locus builder: C++ union-find connected components
# ---------------------------------------------------------------------------


def build_loci(
    em_data: ScoredFragments,
    index: TranscriptIndex,
) -> list[Locus]:
    """Build loci as connected components of transcripts linked by fragments.

    Uses C++ union-find (disjoint-set with path compression and union
    by rank) for fast component detection.

    Parameters
    ----------
    em_data : ScoredFragments
        Global EM data (mRNA + nRNA candidates + per-unit metadata).
    index : TranscriptIndex
        Reference index.

    Returns
    -------
    list[Locus]
    """
    n_transcripts = index.num_transcripts
    offsets = em_data.offsets
    t_indices = em_data.t_indices
    n_units = em_data.n_units

    if n_units == 0 or len(t_indices) == 0:
        return []

    # C++ union-find returns per-component transcript and unit lists
    # in CSR form (offsets + flat index arrays), already sorted ascending.
    n_comp, comp_t_offsets, comp_t_indices, comp_u_offsets, comp_u_indices = _cc_native(
        offsets,
        t_indices,
        np.int32(n_transcripts),
    )

    # Pre-extract transcript coordinates as numpy arrays.
    t_starts_all = index.t_df["start"].values
    t_ends_all = index.t_df["end"].values

    # Integer ref codes for fast sort/compare within merge loops.
    # The "ref" column is already categorical; use its codes directly
    # instead of an O(N log N) np.unique sort on 457K string objects.
    ref_cat = index.t_df["ref"].cat
    _ref_names = ref_cat.categories.values
    _ref_codes = ref_cat.codes.values

    loci = []
    for lid in range(n_comp):
        t_lo = comp_t_offsets[lid]
        t_hi = comp_t_offsets[lid + 1]
        t_idx = comp_t_indices[t_lo:t_hi].copy()

        # Sort transcript intervals by (ref_code, start), then merge
        # overlapping spans to compute the genomic footprint.
        rc = _ref_codes[t_idx]
        ss = t_starts_all[t_idx]
        ee = t_ends_all[t_idx]
        order = np.lexsort((ss, rc))

        merged = []
        span = 0
        prev_rc = int(rc[order[0]])
        prev_s = int(ss[order[0]])
        prev_e = int(ee[order[0]])
        for k in range(1, len(order)):
            j = order[k]
            rj, sj, ej = int(rc[j]), int(ss[j]), int(ee[j])
            if rj != prev_rc or sj > prev_e:
                merged.append((_ref_names[prev_rc], prev_s, prev_e))
                span += prev_e - prev_s
                prev_rc, prev_s, prev_e = rj, sj, ej
            else:
                if ej > prev_e:
                    prev_e = ej
        merged.append((_ref_names[prev_rc], prev_s, prev_e))
        span += prev_e - prev_s

        loci.append(
            Locus(
                locus_id=lid,
                transcript_indices=t_idx,
                unit_indices=comp_u_indices[comp_u_offsets[lid] : comp_u_offsets[lid + 1]].copy(),
                gdna_span=max(span, 1),
                merged_intervals=merged,
            )
        )

    return loci


# ---------------------------------------------------------------------------
# Calibration γ aggregation per locus
# ---------------------------------------------------------------------------


def compute_locus_priors(
    loci: list[Locus],
    index: TranscriptIndex,
    calibration: "CalibrationResult",
    *,
    c_base: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-locus Dirichlet priors from calibration.

    For each locus we compute a **local gDNA mixing fraction** γ_ℓ
    representing "what fraction of fragments in this locus do we
    expect to be gDNA?".  γ_ℓ is set to

        γ_ℓ = Σ_i  region_e_gdna[i]   /   Σ_i  region_n_total[i]

    where the sum runs over all calibration regions overlapping the
    locus's merged transcript-span intervals.

    * ``region_e_gdna[i]`` is the **expected** gDNA count in region
      *i* computed by calibration: ``min(λ_G × mappable_bp[i],
      n_total[i])``.  λ_G is the global gDNA fragment density (per
      mappable base) fitted by the v4 EM and applies everywhere on
      the genome — intergenic, intronic, exonic alike.  The
      ``min(·, n_total)`` clip is a defensive cap that should
      rarely fire outside very small low-coverage regions.
    * ``region_n_total[i]`` is the **observed** total fragment count
      in region *i* (all four strand/splice columns summed).  This
      is *all* fragments seen in the region, not gDNA-only.

    So γ_ℓ = (expected gDNA frags in locus) / (observed total frags
    in locus).  This is a density-aware local mixing prior: in an
    RNA-rich region n_total is large so γ_ℓ is small; in a
    pure-intergenic region n_total ≈ λ_G·mbp so γ_ℓ ≈ 1.

    The priors returned to the EM solver are:

        α_gDNA = γ_ℓ × c_base
        α_RNA  = (1 − γ_ℓ) × c_base

    Downstream, the C++ EM solver uses them in **two** ways:

    1. As a Dirichlet pseudocount in the M-step (negligible vs the
       per-locus fragment count for all but tiny loci).
    2. As the **warm-start ratio** for θ: at iter 0,
       ``θ_gDNA = (α_gDNA/α_RNA) × θ_RNA_total``.  This is what
       actually matters — it places gDNA at the density-expected
       share before EM iterations begin.

    A mode-aware baseline (+0.5 for VBEM, 0 for MAP) is added by
    the solver to compensate for VBEM sparsification bias.

    Loci with no overlapping calibration regions fall back to the
    genome-wide γ = Σ region_e_gdna / Σ region_n_total.

    Returns
    -------
    alpha_gdna : np.ndarray, shape (n_loci,), float64 ≥ 0
    alpha_rna : np.ndarray, shape (n_loci,), float64 ≥ 0
    """
    n_loci = len(loci)
    alpha_gdna = np.empty(n_loci, dtype=np.float64)
    alpha_rna = np.empty(n_loci, dtype=np.float64)

    region_e_gdna = calibration.region_e_gdna
    region_n = calibration.region_n_total

    region_cr = getattr(index, "region_cr", None)
    if region_cr is None or region_e_gdna is None or region_n is None:
        # Global fallback
        total_e = float(region_e_gdna.sum()) if region_e_gdna is not None else 0.0
        total_n = float(region_n.sum()) if region_n is not None else 0.0
        gamma = total_e / max(total_n, 1.0)
        for li in range(n_loci):
            alpha_gdna[li] = gamma * c_base
            alpha_rna[li] = (1.0 - gamma) * c_base
        return alpha_gdna, alpha_rna

    # Global fallback γ for loci with no overlapping regions
    total_e = float(region_e_gdna.sum())
    total_n = float(region_n.sum())
    fallback_gamma = total_e / max(total_n, 1.0)

    for li, locus in enumerate(loci):
        e_sum = 0.0
        n_sum = 0.0
        has_overlap = False
        for ref, start, end in locus.merged_intervals:
            for _s, _e, rid in region_cr.overlap(ref, start, end):
                e_sum += float(region_e_gdna[rid])
                n_sum += float(region_n[rid])
                has_overlap = True

        if has_overlap and n_sum > 0:
            gamma = e_sum / n_sum
        else:
            gamma = fallback_gamma

        alpha_gdna[li] = gamma * c_base
        alpha_rna[li] = (1.0 - gamma) * c_base

    return alpha_gdna, alpha_rna
