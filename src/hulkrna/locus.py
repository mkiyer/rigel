"""hulkrna.locus — Locus graph construction and EM initialization.

Handles everything between "buffer scan complete" and "per-locus EM
loop": connected-component partitioning, per-locus EM data extraction,
nRNA initialization, and Empirical Bayes gDNA priors.
"""

import logging
import math
from collections import defaultdict

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components as _scipy_cc

from .estimator import (
    AbundanceEstimator,
    EM_PRIOR_EPSILON,
    Locus,
    LocusEMInput,
    ScoredFragments,
    estimate_kappa,
)
from .index import TranscriptIndex
from .strand_model import StrandModels

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Empirical Bayes hyperparameters
# ---------------------------------------------------------------------------
#: Fragments needed for locus-level to dominate.
EB_K_LOCUS = 20.0
#: Fragments needed for chrom-level to dominate.
EB_K_CHROM = 50.0

#: Minimum strand specificity denominator (2*SS − 1).
#: When the denominator is below this threshold, the strand signal is too
#: weak to reliably estimate nRNA or gDNA fractions — we return 0.0.
STRAND_DENOM_MIN: float = 0.2

#: Minimum denominator for the strand weight ``2s − 1`` below which
#: the strand component of the hybrid estimator is zeroed out.
_STRAND_DENOM_EPS: float = 0.01


# ---------------------------------------------------------------------------
# Locus builder: scipy connected_components on fragment→transcript graph
# ---------------------------------------------------------------------------


def build_loci(
    em_data: ScoredFragments,
    index: TranscriptIndex,
) -> list[Locus]:
    """Build loci as connected components of transcripts linked by fragments.

    Uses scipy sparse graph connected_components for vectorized
    component detection.

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
    nt = index.num_transcripts
    nrna_base = em_data.nrna_base_index
    offsets = em_data.offsets
    t_indices = em_data.t_indices
    n_units = em_data.n_units

    if n_units == 0 or len(t_indices) == 0:
        return []

    # Map all candidate indices to transcript space [0, nt).
    all_t = t_indices.copy()
    nrna_mask = all_t >= nrna_base
    all_t[nrna_mask] -= nrna_base

    seg_lengths = np.diff(offsets).astype(np.intp)

    # Build sparse adjacency graph by connecting each candidate's
    # transcript to the first candidate's transcript within each unit.
    nonempty = seg_lengths > 0
    first_idx = offsets[:-1].astype(np.intp)

    # First transcript per unit
    first_t = np.empty(n_units, dtype=np.int32)
    first_t[nonempty] = all_t[first_idx[nonempty]]
    first_t[~nonempty] = -1

    # Expand first_t to candidate level
    first_t_expanded = np.repeat(first_t, seg_lengths)

    # Build edges: (first_t_expanded[i], all_t[i]) for every candidate
    row = first_t_expanded
    col = all_t
    valid = (row != col) & (row >= 0) & (col >= 0) & (row < nt) & (col < nt)
    row = row[valid]
    col = col[valid]

    if len(row) > 0:
        data = np.ones(len(row), dtype=np.int8)
        adj = coo_matrix((data, (row, col)), shape=(nt, nt))
        n_comp, labels = _scipy_cc(adj, directed=False)
    else:
        labels = np.arange(nt, dtype=np.int32)

    # Find active transcripts (those appearing in at least one unit)
    active = np.unique(all_t[all_t >= 0])

    # Group active transcripts by component label
    active_labels = labels[active]
    unique_labels, inverse = np.unique(active_labels, return_inverse=True)

    # Build component_map: label → sorted transcript list
    component_map: dict[int, list[int]] = {}
    for i, t in enumerate(active):
        lbl = int(unique_labels[inverse[i]])
        component_map.setdefault(lbl, []).append(int(t))

    # Assign EM units to components (using first transcript's label)
    unit_labels = np.full(n_units, -1, dtype=np.int32)
    unit_labels[nonempty] = labels[first_t[nonempty]]

    # Build per-component unit lists
    root_to_units: dict[int, list[int]] = {lbl: [] for lbl in component_map}
    for u in range(n_units):
        lbl = int(unit_labels[u])
        if lbl >= 0 and lbl in root_to_units:
            root_to_units[lbl].append(u)

    # Build Locus objects
    loci = []
    for lid, (lbl, t_list) in enumerate(sorted(component_map.items())):
        t_arr = np.array(sorted(t_list), dtype=np.int32)
        u_arr = np.array(sorted(root_to_units.get(lbl, [])), dtype=np.int32)

        loci.append(Locus(
            locus_id=lid,
            transcript_indices=t_arr,
            unit_indices=u_arr,
        ))

    return loci


# ---------------------------------------------------------------------------
# Locus EM data builder
# ---------------------------------------------------------------------------


def build_locus_em_data(
    locus: Locus,
    em_data: ScoredFragments,
    estimator: AbundanceEstimator,
    index: TranscriptIndex,
    mean_frag: float,
    gdna_init: float,
    *,
    _cache: dict | None = None,
) -> LocusEMInput:
    """Extract and renumber global ScoredFragments into a per-locus sub-problem.

    Component layout per locus::

        [0, n_t)           - mRNA (one per local transcript)
        [n_t, 2*n_t)       - nRNA (one per local transcript)
        [2*n_t]            - gDNA (ONE shadow for entire locus)

    Only UNSPLICED units get a gDNA candidate.  Spliced units see
    only mRNA/nRNA components.

    Parameters
    ----------
    locus : Locus
    em_data : ScoredFragments
    estimator : AbundanceEstimator
    index : TranscriptIndex
    mean_frag : float
    gdna_init : float
        Empirical Bayes estimated gDNA count for this locus.
    _cache : dict or None
        Optional pre-extracted arrays to avoid repeated DataFrame
        access.  Expected keys: ``"t_starts"``, ``"t_ends"``,
        ``"t_lengths"``, ``"local_map"`` (reusable scratch buffer).
    """
    t_arr = locus.transcript_indices
    n_t = len(t_arr)
    gdna_idx = 2 * n_t  # single gDNA component index
    n_components = 2 * n_t + 1

    nrna_base = em_data.nrna_base_index
    n_local_units = len(locus.unit_indices)

    # Build global → local mapping array (fast lookup, no dict)
    max_global = max(int(nrna_base) + int(t_arr.max()) + 1,
                     int(t_arr.max()) + 1) if n_t > 0 else 0

    # Reuse pre-allocated local_map buffer if provided
    if _cache is not None and "local_map" in _cache:
        local_map = _cache["local_map"]
        if len(local_map) < max_global:
            # Grow if needed (rare — only when locus has very large indices)
            local_map = np.full(max_global, -1, dtype=np.int32)
            _cache["local_map"] = local_map
        else:
            local_map[:max_global] = -1
    else:
        local_map = np.full(max_global, -1, dtype=np.int32)
    for local_i in range(n_t):
        gt = int(t_arr[local_i])
        local_map[gt] = local_i                     # mRNA
        local_map[nrna_base + gt] = n_t + local_i   # nRNA

    # Gather all candidate ranges for this locus's units at once
    global_offsets = em_data.offsets
    unit_starts = global_offsets[locus.unit_indices].astype(np.intp)
    unit_ends = global_offsets[locus.unit_indices + 1].astype(np.intp)
    seg_lens = unit_ends - unit_starts
    total_cands = int(seg_lens.sum())

    if total_cands > 0:
        # Build flat index array for all candidates in this locus.
        cum_lens = np.empty(n_local_units + 1, dtype=np.intp)
        cum_lens[0] = 0
        np.cumsum(seg_lens, out=cum_lens[1:])

        global_flat_idx = (
            np.repeat(unit_starts, seg_lens)
            + np.arange(total_cands, dtype=np.intp)
            - np.repeat(cum_lens[:-1], seg_lens)
        )

        # Extract all candidate data at once
        all_gidx = em_data.t_indices[global_flat_idx]
        all_ll = em_data.log_liks[global_flat_idx]
        all_cc = em_data.count_cols[global_flat_idx]
        all_cw = em_data.coverage_weights[global_flat_idx]
        all_txs = em_data.tx_starts[global_flat_idx]
        all_txe = em_data.tx_ends[global_flat_idx]

        # Map global indices to local component indices
        safe_gidx = np.clip(all_gidx, 0, max_global - 1)
        all_lidx = local_map[safe_gidx]
        all_lidx[all_gidx < 0] = -1
        all_lidx[all_gidx >= max_global] = -1

        # Which local unit each candidate belongs to
        unit_of_cand = np.repeat(np.arange(n_local_units, dtype=np.int32),
                                 seg_lens)

        # Filter valid candidates
        valid = all_lidx >= 0
        v_lidx = all_lidx[valid]
        v_ll = all_ll[valid]
        v_cc = all_cc[valid]
        v_cw = all_cw[valid]
        v_txs = all_txs[valid]
        v_txe = all_txe[valid]
        v_unit = unit_of_cand[valid]

        # Deduplicate: keep best log_lik per (unit, local_idx).
        compound_key = v_unit.astype(np.int64) * n_components + v_lidx
        order = np.lexsort((-v_ll, compound_key))
        sorted_keys = compound_key[order]

        first_mask = np.empty(len(sorted_keys), dtype=bool)
        first_mask[0] = True
        first_mask[1:] = sorted_keys[1:] != sorted_keys[:-1]

        dedup_lidx = v_lidx[order][first_mask]
        dedup_ll = v_ll[order][first_mask]
        dedup_cc = v_cc[order][first_mask]
        dedup_cw = v_cw[order][first_mask]
        dedup_txs = v_txs[order][first_mask]
        dedup_txe = v_txe[order][first_mask]
        dedup_unit = v_unit[order][first_mask]
    else:
        dedup_lidx = np.empty(0, dtype=np.int32)
        dedup_ll = np.empty(0, dtype=np.float64)
        dedup_cc = np.empty(0, dtype=np.uint8)
        dedup_cw = np.empty(0, dtype=np.float64)
        dedup_txs = np.empty(0, dtype=np.int32)
        dedup_txe = np.empty(0, dtype=np.int32)
        dedup_unit = np.empty(0, dtype=np.int32)

    # Add gDNA candidates for unspliced units
    is_spl = em_data.is_spliced[locus.unit_indices]
    gdna_lls = em_data.gdna_log_liks[locus.unit_indices]
    valid_gdna = (~is_spl) & np.isfinite(gdna_lls)
    n_gdna = int(valid_gdna.sum())

    # Compute true locus span for gDNA per-fragment effective length.
    # Span extends from leftmost transcript/fragment start to rightmost
    # transcript/fragment end, capturing any overhanging fragments.
    _t_starts_all = _cache["t_starts"] if _cache is not None else index.t_df["start"].values
    _t_ends_all = _cache["t_ends"] if _cache is not None else index.t_df["end"].values
    t_starts = _t_starts_all[t_arr]
    t_ends = _t_ends_all[t_arr]
    locus_start = int(t_starts.min())
    locus_end = int(t_ends.max())
    footprints = em_data.genomic_footprints[locus.unit_indices]
    locus_span = float(locus_end - locus_start)

    if n_gdna > 0:
        gdna_units = np.arange(n_local_units, dtype=np.int32)[valid_gdna]
        gdna_lidx_arr = np.full(n_gdna, gdna_idx, dtype=np.int32)
        gdna_ll_arr = gdna_lls[valid_gdna].copy()
        gdna_cc_arr = np.zeros(n_gdna, dtype=np.uint8)

        # Per-fragment effective length correction for gDNA is now
        # handled at EM time via BiasProfile.  Store tx_start=0 and
        # tx_end=footprint so the bias model sees frag_len=footprint
        # on a profile of length=locus_span, yielding
        # eff_len = max(locus_span - footprint + 1, 1).
        gdna_footprints_arr = footprints[valid_gdna].astype(np.int32)

        final_lidx = np.concatenate([dedup_lidx, gdna_lidx_arr])
        final_ll = np.concatenate([dedup_ll, gdna_ll_arr])
        final_cc = np.concatenate([dedup_cc, gdna_cc_arr])
        final_cw = np.concatenate([dedup_cw, np.ones(n_gdna, dtype=np.float64)])
        # gDNA: tx_start=0, tx_end=footprint (bias model uses frag_len)
        final_txs = np.concatenate([dedup_txs, np.zeros(n_gdna, dtype=np.int32)])
        final_txe = np.concatenate([dedup_txe, gdna_footprints_arr])
        final_unit = np.concatenate([dedup_unit, gdna_units])
    else:
        final_lidx = dedup_lidx
        final_ll = dedup_ll
        final_cc = dedup_cc
        final_cw = dedup_cw
        final_txs = dedup_txs
        final_txe = dedup_txe
        final_unit = dedup_unit

    # Sort by unit to reconstruct CSR
    if len(final_unit) > 0:
        sort_order = np.argsort(final_unit, kind='stable')
        final_lidx = final_lidx[sort_order]
        final_ll = final_ll[sort_order]
        final_cc = final_cc[sort_order]
        final_cw = final_cw[sort_order]
        final_txs = final_txs[sort_order]
        final_txe = final_txe[sort_order]
        final_unit = final_unit[sort_order]

        bin_counts = np.bincount(final_unit, minlength=n_local_units)
        local_offsets = np.empty(n_local_units + 1, dtype=np.int64)
        local_offsets[0] = 0
        np.cumsum(bin_counts, out=local_offsets[1:])
    else:
        local_offsets = np.zeros(n_local_units + 1, dtype=np.int64)

    # Per-unit locus_t and locus_ct (trivial gather)
    local_locus_t = em_data.locus_t_indices[locus.unit_indices]
    local_locus_ct = em_data.locus_count_cols[locus.unit_indices]

    # Build local init vectors
    #
    # Only mRNA unambig_counts go into unambig_totals because those
    # fragments are deterministically assigned (NOT sent to the EM),
    # so they are disjoint from em_totals.  nrna_init is derived from
    # the SAME fragments the EM processes (intronic strand excess),
    # so putting it here would double-count in the M-step.
    unambig_totals = np.zeros(n_components, dtype=np.float64)
    unambig_totals[:n_t] = estimator.unambig_counts[t_arr].sum(axis=1)
    # nRNA slot intentionally left at 0 — see nrna double-counting fix.
    # gDNA slot intentionally left at 0 — gdna_init is derived from the
    # SAME unspliced fragments the EM processes (strand-based EB prior),
    # so adding it to unambig_totals would double-count in the M-step
    # (theta = unambig + em_totals + prior).  gdna_init is used only
    # for prior gating below (prior[gdna] = 0 when gdna_init = 0).

    # Build effective lengths — all set to 1.0 because per-fragment
    # effective length correction is now applied at EM time via
    # BiasProfile (see _apply_bias_correction in estimator.py).
    eff_len = np.ones(n_components, dtype=np.float64)

    # Build per-component lengths for uniform bias correction.
    #
    # Component layout:
    #   [0, n_t)       → mRNA: length = exonic transcript length
    #   [n_t, 2*n_t)   → nRNA: length = genomic span (incl introns)
    #   [2*n_t]        → gDNA: length = locus span
    #
    # Phase H optimisation: pass a flat int64 array of component
    # lengths to _apply_bias_correction_uniform directly, avoiding
    # 2*n_t+1 BiasProfile object + np.arange allocations per locus.
    _t_lengths_all = _cache["t_lengths"] if _cache is not None else index.t_df["length"].values
    mrna_lengths = _t_lengths_all[t_arr]
    nrna_lengths = (
        _t_ends_all[t_arr]
        - _t_starts_all[t_arr]
    )
    bias_profiles = np.empty(n_components, dtype=np.int64)
    bias_profiles[:n_t] = mrna_lengths
    bias_profiles[n_t:2*n_t] = nrna_lengths
    bias_profiles[2*n_t] = int(locus_span)

    # Build prior — start at epsilon (numerical floor) for every
    # component; the real prior mass comes from the OVR coverage
    # weighting applied later in run_locus_em().
    prior = np.full(n_components, EM_PRIOR_EPSILON, dtype=np.float64)

    # Zero gDNA prior when there are no unspliced fragments
    if gdna_init == 0.0:
        prior[gdna_idx] = 0.0

    # Zero nRNA prior for single-exon transcripts
    if estimator._transcript_spans is not None:
        if estimator._exonic_lengths is not None:
            t_exon = estimator._exonic_lengths[t_arr]
        else:
            t_exon = estimator._t_eff_len[t_arr] + mean_frag - 1.0
        single_exon = estimator._transcript_spans[t_arr] <= t_exon
        prior[n_t:2*n_t][single_exon] = 0.0

    # Zero nRNA prior when nrna_init is zero for that transcript.
    # We read nrna_init directly (not from unambig_totals, which no
    # longer carries nrna_init after the double-counting fix).
    nrna_init_local = estimator.nrna_init[t_arr]
    prior[n_t:2*n_t][nrna_init_local == 0.0] = 0.0

    # Zero unambig_totals for dead components—prevents the EM M-step
    # (theta_new = unambig_totals + em_totals + prior) from seeding
    # non-zero theta for components whose prior has been zeroed.
    unambig_totals[prior == 0.0] = 0.0

    return LocusEMInput(
        locus=locus,
        offsets=local_offsets,
        t_indices=final_lidx.astype(np.int32),
        log_liks=final_ll.astype(np.float64),
        count_cols=final_cc.astype(np.uint8),
        coverage_weights=final_cw.astype(np.float64),
        tx_starts=final_txs.astype(np.int32),
        tx_ends=final_txe.astype(np.int32),
        locus_t_indices=local_locus_t.astype(np.int32),
        locus_count_cols=local_locus_ct.astype(np.uint8),
        n_transcripts=n_t,
        n_components=n_components,
        local_to_global_t=t_arr.copy(),
        unambig_totals=unambig_totals,
        nrna_init=estimator.nrna_init[t_arr].copy(),
        gdna_init=gdna_init,
        effective_lengths=eff_len,
        prior=prior,
        bias_profiles=bias_profiles,
        nrna_frac_alpha=estimator.nrna_frac_alpha[t_arr].copy(),
        nrna_frac_beta=estimator.nrna_frac_beta[t_arr].copy(),
    )


# ---------------------------------------------------------------------------
# nRNA initialization (strand-corrected intronic evidence)
# ---------------------------------------------------------------------------


def compute_nrna_init(
    transcript_intronic_sense: np.ndarray,
    transcript_intronic_antisense: np.ndarray,
    transcript_spans: np.ndarray,
    exonic_lengths: np.ndarray,
    mean_frag: float,
    strand_models: StrandModels,
) -> np.ndarray:
    """Compute per-transcript nRNA initialization from intronic evidence.

    Model::

        sense_int   = gDNA_int/2 + nRNA_int × SS
        anti_int    = gDNA_int/2 + nRNA_int × (1-SS)

    Exact solution::

        nRNA_int = (sense_int - anti_int) / (2SS - 1)

    When 2SS − 1 ≤ 0.2, returns zeros.
    Single-exon transcripts get nrna_init = 0.
    """
    ss = strand_models.strand_specificity
    denom = 2.0 * ss - 1.0

    if denom <= STRAND_DENOM_MIN:
        nrna_init = np.zeros_like(transcript_intronic_sense)
    else:
        raw = (
            transcript_intronic_sense - transcript_intronic_antisense
        ) / denom
        nrna_init = np.maximum(0.0, raw)

    intronic_span = np.maximum(transcript_spans - exonic_lengths, 0.0)
    nrna_init[intronic_span <= 0] = 0.0

    return nrna_init


# ---------------------------------------------------------------------------
# Empirical Bayes gDNA prior (locus → chromosome → global)
# ---------------------------------------------------------------------------


def compute_gdna_rate_from_strand(
    sense: float,
    antisense: float,
    ss: float,
) -> float:
    """Strand-corrected gDNA rate from sense/antisense counts.

    G = 2(A·SS - S·(1-SS)) / (2SS-1)
    rate = G / (S + A)  clamped to [0, 1]

    Returns 0.0 when evidence is insufficient.
    """
    denom_ss = 2.0 * ss - 1.0
    if denom_ss <= STRAND_DENOM_MIN:
        return 0.0
    total = sense + antisense
    if total == 0:
        return 0.0
    g = 2.0 * (antisense * ss - sense * (1.0 - ss)) / denom_ss
    g = max(g, 0.0)
    rate = g / total
    return min(rate, 1.0)


def compute_gdna_rate_hybrid(
    sense: float,
    antisense: float,
    exonic_bp: float,
    ss: float,
    intergenic_density: float,
) -> tuple[float, float]:
    """Hybrid density + strand gDNA rate estimate.

    Combines strand-based and density-based estimators using
    inverse-variance weighting ``W = (2s − 1)²``.

    * **Strand component** — isolates gDNA from sense/antisense
      imbalance (same formula as :func:`compute_gdna_rate_from_strand`).
    * **Density component** — estimates gDNA from intergenic background
      density scaled by the region's exonic territory.

    For stranded libraries (W ≈ 1), the strand signal dominates.
    For weakly-stranded libraries (W ≈ 0), the density signal provides
    a fallback that is strictly better than the strand-only 0.0.

    Parameters
    ----------
    sense, antisense : float
        Unspliced sense/antisense fragment counts.
    exonic_bp : float
        Total exonic base pairs in this region.
    ss : float
        Library strand specificity ∈ [0.5, 1.0].
    intergenic_density : float
        Background gDNA density (frags / bp) from intergenic regions.
        Pass 0.0 to disable the density component.

    Returns
    -------
    rate : float
        Estimated gDNA fraction ∈ [0, 1].
    evidence : float
        Total fragment count (sense + antisense).
    """
    total = sense + antisense
    if total == 0:
        return 0.0, 0.0

    # --- Strand component ---
    denom_ss = 2.0 * ss - 1.0
    W = denom_ss ** 2 if denom_ss > _STRAND_DENOM_EPS else 0.0

    if W > 0:
        g = 2.0 * (antisense * ss - sense * (1.0 - ss)) / denom_ss
        g = max(g, 0.0)
        strand_rate = min(g / total, 1.0)
    else:
        strand_rate = 0.0

    # --- Density component ---
    if exonic_bp > 0 and intergenic_density > 0:
        expected_gdna = intergenic_density * exonic_bp
        density_rate = min(expected_gdna / total, 1.0)
    else:
        density_rate = 0.0

    # --- Weighted combination ---
    rate = W * strand_rate + (1.0 - W) * density_rate
    return min(max(rate, 0.0), 1.0), total


def compute_eb_gdna_priors(
    loci: list[Locus],
    em_data: ScoredFragments,
    estimator: AbundanceEstimator,
    index: TranscriptIndex,
    strand_models: StrandModels,
    *,
    intergenic_density: float = 0.0,
    kappa_chrom: float | None = None,
    kappa_locus: float | None = None,
    mom_min_evidence_chrom: float = 50.0,
    mom_min_evidence_locus: float = 30.0,
    kappa_min: float = 2.0,
    kappa_max: float = 200.0,
    kappa_fallback: float = 5.0,
    kappa_min_obs: int = 20,
) -> list[float]:
    """Compute empirical Bayes gDNA initialization per locus.

    Hierarchical weighted shrinkage: locus → chromosome → global.

    Uses a hybrid density + strand signal for gDNA rate estimation
    at each hierarchical level (see :func:`compute_gdna_rate_hybrid`).
    When ``kappa_chrom`` or ``kappa_locus`` is ``None``, the shrinkage
    concentration is auto-estimated via Method of Moments using
    :func:`~hulkrna.estimator.estimate_kappa`.

    Parameters
    ----------
    loci : list[Locus]
    em_data : ScoredFragments
    estimator : AbundanceEstimator
    index : TranscriptIndex
    strand_models : StrandModels
    intergenic_density : float
        Background gDNA density (frags / bp) from intergenic regions.
        Pass 0.0 to fall back to strand-only estimation.
    kappa_chrom : float or None
        Shrinkage pseudo-count for chrom → global.  ``None`` (default)
        auto-estimates via Method of Moments.
    kappa_locus : float or None
        Shrinkage pseudo-count for locus → chrom.  ``None`` (default)
        auto-estimates via Method of Moments.
    mom_min_evidence_chrom : float
        Minimum evidence for a chromosome to contribute to MoM κ_chrom.
    mom_min_evidence_locus : float
        Minimum evidence for a locus to contribute to MoM κ_locus.
    kappa_min, kappa_max : float
        Clamps for MoM-estimated κ values.
    kappa_fallback : float
        Fallback κ when too few units pass the evidence filter.
    kappa_min_obs : int
        Minimum number of units for MoM; fewer triggers fallback.

    Returns
    -------
    list[float]
        gDNA init count per locus (same order as ``loci``).
    """
    ss = strand_models.strand_specificity
    t_sense = estimator.transcript_unspliced_sense
    t_anti = estimator.transcript_unspliced_antisense
    t_refs = index.t_df["ref"].values
    t_exonic = index.t_df["length"].values.astype(np.float64)

    # --- Global level (hybrid density+strand) ---
    total_sense = float(t_sense.sum())
    total_anti = float(t_anti.sum())
    total_exonic_bp = float(t_exonic.sum())
    global_rate, global_n = compute_gdna_rate_hybrid(
        total_sense, total_anti, total_exonic_bp, ss, intergenic_density,
    )

    # --- Chromosome level (hybrid density+strand) ---
    chrom_sense: dict[str, float] = defaultdict(float)
    chrom_anti: dict[str, float] = defaultdict(float)
    chrom_exonic: dict[str, float] = defaultdict(float)
    for t_idx in range(len(t_sense)):
        ref = str(t_refs[t_idx])
        chrom_sense[ref] += t_sense[t_idx]
        chrom_anti[ref] += t_anti[t_idx]
        chrom_exonic[ref] += t_exonic[t_idx]

    chrom_rate: dict[str, float] = {}
    chrom_n: dict[str, float] = {}
    for ref in chrom_sense:
        r, n = compute_gdna_rate_hybrid(
            chrom_sense[ref], chrom_anti[ref], chrom_exonic[ref],
            ss, intergenic_density,
        )
        chrom_rate[ref] = r
        chrom_n[ref] = n

    # MoM κ for chrom → global shrinkage
    if kappa_chrom is None:
        chrom_rate_arr = np.array(list(chrom_rate.values()), dtype=np.float64)
        chrom_n_arr = np.array(list(chrom_n.values()), dtype=np.float64)
        k_chrom = estimate_kappa(
            chrom_rate_arr, chrom_n_arr, mom_min_evidence_chrom,
            kappa_min=kappa_min, kappa_max=kappa_max,
            kappa_fallback=kappa_fallback, kappa_min_obs=kappa_min_obs,
        )
        logger.debug(f"gDNA κ_chrom (MoM auto): {k_chrom:.1f}")
    else:
        k_chrom = kappa_chrom

    logger.info(
        f"gDNA EB: global_rate={global_rate:.4f} "
        f"(N={global_n:.0f}), κ_chrom={k_chrom:.1f}, "
        f"n_chroms={len(chrom_rate)}"
    )

    # Shrink chrom toward global
    chrom_shrunk: dict[str, float] = {}
    for ref in chrom_rate:
        n = chrom_n[ref]
        w = n / (n + k_chrom) if (n + k_chrom) > 0 else 0.0
        chrom_shrunk[ref] = w * chrom_rate[ref] + (1.0 - w) * global_rate

    # --- Per-locus level (hybrid density+strand) ---
    locus_rates: list[float] = []
    locus_evidence: list[float] = []
    locus_parents: list[float] = []

    for locus in loci:
        t_arr = locus.transcript_indices
        locus_sense = float(t_sense[t_arr].sum())
        locus_anti = float(t_anti[t_arr].sum())
        locus_exonic_bp = float(t_exonic[t_arr].sum())

        locus_rate, locus_n = compute_gdna_rate_hybrid(
            locus_sense, locus_anti, locus_exonic_bp, ss, intergenic_density,
        )

        # Determine primary chromosome for this locus
        ref_counts: dict[str, float] = defaultdict(float)
        for t_idx in t_arr:
            ref = str(t_refs[int(t_idx)])
            ref_counts[ref] += t_sense[int(t_idx)] + t_anti[int(t_idx)]
        if ref_counts:
            primary_ref = max(ref_counts, key=ref_counts.get)
        elif len(t_arr) > 0:
            primary_ref = str(t_refs[int(t_arr[0])])
        else:
            primary_ref = ""

        parent_rate = chrom_shrunk.get(primary_ref, global_rate)

        locus_rates.append(locus_rate)
        locus_evidence.append(locus_n)
        locus_parents.append(parent_rate)

    # MoM κ for locus → chrom shrinkage
    if kappa_locus is None:
        locus_rate_arr = np.array(locus_rates, dtype=np.float64)
        locus_n_arr = np.array(locus_evidence, dtype=np.float64)
        k_locus = estimate_kappa(
            locus_rate_arr, locus_n_arr, mom_min_evidence_locus,
            kappa_min=kappa_min, kappa_max=kappa_max,
            kappa_fallback=kappa_fallback, kappa_min_obs=kappa_min_obs,
        )
        logger.debug(f"gDNA κ_locus (MoM auto): {k_locus:.1f}")
    else:
        k_locus = kappa_locus

    # --- Diagnostic logging ---
    locus_rate_arr = np.array(locus_rates, dtype=np.float64)
    n_zero_gdna = int(np.sum(locus_rate_arr == 0.0))
    n_high_gdna = int(np.sum(locus_rate_arr >= 0.05))
    median_locus_rate = float(np.median(locus_rate_arr)) if len(locus_rates) > 0 else 0.0
    logger.info(
        f"gDNA EB: κ_locus={k_locus:.1f}, "
        f"median_locus_rate={median_locus_rate:.4f}, "
        f"n_zero={n_zero_gdna}/{len(loci)}, "
        f"n_high(≥0.05)={n_high_gdna}/{len(loci)}"
    )

    # Shrink locus toward parent (chromosome) and compute gdna_init
    gdna_inits: list[float] = []
    for i, locus in enumerate(loci):
        locus_n = locus_evidence[i]
        locus_rate = locus_rates[i]
        parent_rate = locus_parents[i]

        w = locus_n / (locus_n + k_locus) if (locus_n + k_locus) > 0 else 0.0
        shrunk_rate = w * locus_rate + (1.0 - w) * parent_rate

        # gDNA init = shrunk_rate × unspliced fragment count in locus
        n_unspliced = 0
        for u in locus.unit_indices:
            if not em_data.is_spliced[u]:
                n_unspliced += 1

        gdna_init = shrunk_rate * n_unspliced
        gdna_inits.append(max(gdna_init, 0.0))

    return gdna_inits
