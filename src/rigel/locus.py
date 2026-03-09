"""rigel.locus — Locus graph construction and EM initialization.

Handles everything between "buffer scan complete" and "per-locus EM
loop": connected-component partitioning, per-locus EM data extraction,
nRNA initialization, and Empirical Bayes gDNA priors.
"""

import logging
from collections import defaultdict

import numpy as np
from ._em_impl import connected_components as _cc_native

from .config import EMConfig
from .estimator import (
    AbundanceEstimator,
    EM_PRIOR_EPSILON,
    Locus,
    LocusEMInput,
    ScoredFragments,
    _STRAND_DENOM_EPS,
    estimate_kappa,
)
from .index import TranscriptIndex
from .strand_model import StrandModels

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Empirical Bayes hyperparameters
# ---------------------------------------------------------------------------

#: Minimum strand specificity denominator (2*SS − 1).
#: When the denominator is below this threshold, the strand signal is too
#: weak to reliably estimate nRNA or gDNA fractions — we return 0.0.
STRAND_DENOM_MIN: float = 0.2


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
    nt = index.num_transcripts
    nrna_base = em_data.nrna_base_index
    offsets = em_data.offsets
    t_indices = em_data.t_indices
    n_units = em_data.n_units

    if n_units == 0 or len(t_indices) == 0:
        return []

    # Build nRNA → transcript CSR mapping for connected components
    t_to_nrna = index.t_to_nrna_arr
    num_nrna = index.num_nrna
    nrna_counts = np.zeros(num_nrna, dtype=np.int32)
    np.add.at(nrna_counts, t_to_nrna, 1)
    nrna_to_t_offsets = np.zeros(num_nrna + 1, dtype=np.int32)
    np.cumsum(nrna_counts, out=nrna_to_t_offsets[1:])
    nrna_to_t_indices = np.empty(nt, dtype=np.int32)
    cursor = nrna_to_t_offsets[:-1].copy()
    for t in range(nt):
        ni = t_to_nrna[t]
        nrna_to_t_indices[cursor[ni]] = t
        cursor[ni] += 1

    # C++ union-find returns per-component transcript and unit lists
    # in CSR form, already sorted ascending.
    n_comp, ct_off, ct_flat, cu_off, cu_flat = _cc_native(
        offsets, t_indices,
        np.int32(nrna_base), np.int32(nt),
        nrna_to_t_offsets, nrna_to_t_indices,
    )

    loci = []
    for lid in range(n_comp):
        loci.append(Locus(
            locus_id=lid,
            transcript_indices=ct_flat[ct_off[lid]:ct_off[lid + 1]].copy(),
            unit_indices=cu_flat[cu_off[lid]:cu_off[lid + 1]].copy(),
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

        [0, n_t)                - mRNA (one per local transcript)
        [n_t, n_t + n_nrna)     - nRNA (one per unique nRNA span)
        [n_t + n_nrna]          - gDNA (ONE shadow for entire locus)

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
    nrna_base = em_data.nrna_base_index
    n_local_units = len(locus.unit_indices)

    # ------------------------------------------------------------------
    # Build unique local nRNA mapping for this locus
    # ------------------------------------------------------------------
    global_t_to_nrna = index.t_to_nrna_arr  # global transcript → global nRNA
    # Unique global nRNA indices used by this locus's transcripts
    global_nrna_for_locus = global_t_to_nrna[t_arr]  # [n_t]
    unique_global_nrna, inverse = np.unique(global_nrna_for_locus,
                                            return_inverse=True)
    n_nrna = len(unique_global_nrna)
    # inverse[i] gives local nRNA index for local transcript i
    local_t_to_local_nrna = inverse.astype(np.int32)
    local_to_global_nrna = unique_global_nrna.astype(np.int32)

    # Build nRNA→transcript CSR (local indices)
    nrna_to_t_offsets = np.zeros(n_nrna + 1, dtype=np.int32)
    for local_t in range(n_t):
        nrna_to_t_offsets[local_t_to_local_nrna[local_t] + 1] += 1
    np.cumsum(nrna_to_t_offsets, out=nrna_to_t_offsets)
    nrna_to_t_indices = np.empty(n_t, dtype=np.int32)
    fill_pos = nrna_to_t_offsets[:-1].copy()
    for local_t in range(n_t):
        ln = local_t_to_local_nrna[local_t]
        nrna_to_t_indices[fill_pos[ln]] = local_t
        fill_pos[ln] += 1

    gdna_idx = n_t + n_nrna  # single gDNA component index
    n_components = n_t + n_nrna + 1

    # Build global → local mapping array (fast lookup, no dict)
    # mRNA: global_t → local_t
    # nRNA: nrna_base + global_nrna → n_t + local_nrna
    max_mRNA_global = int(t_arr.max()) + 1 if n_t > 0 else 0
    max_nRNA_global = int(nrna_base) + int(unique_global_nrna.max()) + 1 \
        if n_nrna > 0 else 0
    max_global = max(max_mRNA_global, max_nRNA_global, 1)

    # Reuse pre-allocated local_map buffer if provided
    if _cache is not None and "local_map" in _cache:
        local_map = _cache["local_map"]
        if len(local_map) < max_global:
            local_map = np.full(max_global, -1, dtype=np.int32)
            _cache["local_map"] = local_map
        else:
            local_map[:max_global] = -1
    else:
        local_map = np.full(max_global, -1, dtype=np.int32)
    for local_i in range(n_t):
        gt = int(t_arr[local_i])
        local_map[gt] = local_i                     # mRNA
    for local_n in range(n_nrna):
        gn = int(unique_global_nrna[local_n])
        local_map[nrna_base + gn] = n_t + local_n   # nRNA

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
    #   [0, n_t)                → mRNA: length = exonic transcript length
    #   [n_t, n_t + n_nrna)    → nRNA: length = genomic span of the nRNA
    #   [n_t + n_nrna]         → gDNA: length = locus span
    _t_lengths_all = _cache["t_lengths"] if _cache is not None else index.t_df["length"].values
    mrna_lengths = _t_lengths_all[t_arr]
    # nRNA span comes from the nrna_df (global nRNA table)
    nrna_spans = (index.nrna_df["end"].values - index.nrna_df["start"].values)
    nrna_lengths = nrna_spans[unique_global_nrna]
    bias_profiles = np.empty(n_components, dtype=np.int64)
    bias_profiles[:n_t] = mrna_lengths
    bias_profiles[n_t:n_t + n_nrna] = nrna_lengths
    bias_profiles[n_t + n_nrna] = int(locus_span)

    # Build prior — start at epsilon (numerical floor) for every
    # component; the real prior mass comes from the OVR coverage
    # weighting applied later in run_locus_em().
    prior = np.full(n_components, EM_PRIOR_EPSILON, dtype=np.float64)

    # Zero gDNA prior when there are no unspliced fragments
    if gdna_init == 0.0:
        prior[gdna_idx] = 0.0

    # Zero nRNA prior for nRNAs where all sharing transcripts are single-exon
    if estimator._transcript_spans is not None:
        if estimator._exonic_lengths is not None:
            t_exon = estimator._exonic_lengths[t_arr]
        else:
            t_exon = estimator._t_eff_len[t_arr] + mean_frag - 1.0
        single_exon = estimator._transcript_spans[t_arr] <= t_exon
        # For each local nRNA, check if ALL sharing transcripts are single-exon
        for ln in range(n_nrna):
            ts_start = nrna_to_t_offsets[ln]
            ts_end = nrna_to_t_offsets[ln + 1]
            all_single = True
            for j in range(ts_start, ts_end):
                if not single_exon[nrna_to_t_indices[j]]:
                    all_single = False
                    break
            if all_single:
                prior[n_t + ln] = 0.0

    # Zero nRNA prior when nrna_init is zero for that nRNA.
    nrna_init_local = estimator.nrna_init[unique_global_nrna]
    prior[n_t:n_t + n_nrna][nrna_init_local == 0.0] = 0.0

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
        n_nrna=n_nrna,
        n_components=n_components,
        local_to_global_t=t_arr.copy(),
        local_to_global_nrna=local_to_global_nrna,
        local_t_to_local_nrna=local_t_to_local_nrna,
        nrna_to_t_offsets=nrna_to_t_offsets,
        nrna_to_t_indices=nrna_to_t_indices,
        unambig_totals=unambig_totals,
        nrna_init=estimator.nrna_init[unique_global_nrna].copy(),
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
    intronic_sense: np.ndarray,
    intronic_antisense: np.ndarray,
    nrna_spans: np.ndarray,
    nrna_max_exonic: np.ndarray,
    strand_models: StrandModels,
) -> np.ndarray:
    """Compute per-nRNA initialization from intronic evidence.

    Model::

        sense_int   = gDNA_int/2 + nRNA_int × SS
        anti_int    = gDNA_int/2 + nRNA_int × (1-SS)

    Exact solution::

        nRNA_int = (sense_int - anti_int) / (2SS - 1)

    When 2SS − 1 ≤ 0.2, returns zeros.
    nRNAs with no intronic span (max_exonic >= span) get nrna_init = 0.

    Parameters
    ----------
    intronic_sense, intronic_antisense : np.ndarray
        float64[num_nrna] - per-nRNA intronic counts.
    nrna_spans : np.ndarray
        float64[num_nrna] - genomic span (end - start) per nRNA.
    nrna_max_exonic : np.ndarray
        float64[num_nrna] - max exonic length among sharing transcripts.
    strand_models : StrandModels
    """
    strand_spec = strand_models.strand_specificity
    denom = 2.0 * strand_spec - 1.0

    if denom <= STRAND_DENOM_MIN:
        nrna_init = np.zeros_like(intronic_sense)
    else:
        raw = (
            intronic_sense - intronic_antisense
        ) / denom
        nrna_init = np.maximum(0.0, raw)

    intronic_span = np.maximum(nrna_spans - nrna_max_exonic, 0.0)
    nrna_init[intronic_span <= 0] = 0.0

    return nrna_init


# ---------------------------------------------------------------------------
# Empirical Bayes gDNA prior (locus → chromosome → global)
# ---------------------------------------------------------------------------


def compute_gdna_rate_from_strand(
    sense: float,
    antisense: float,
    strand_spec: float,
) -> float:
    """Strand-corrected gDNA rate from sense/antisense counts.

    G = 2(A·SS - S·(1-SS)) / (2SS-1)
    rate = G / (S + A)  clamped to [0, 1]

    Returns 0.0 when evidence is insufficient.
    """
    denom_ss = 2.0 * strand_spec - 1.0
    if denom_ss <= STRAND_DENOM_MIN:
        return 0.0
    total = sense + antisense
    if total == 0:
        return 0.0
    g = 2.0 * (antisense * strand_spec - sense * (1.0 - strand_spec)) / denom_ss
    g = max(g, 0.0)
    rate = g / total
    return min(rate, 1.0)


def compute_gdna_rate_hybrid(
    sense: float,
    antisense: float,
    exonic_bp: float,
    strand_spec: float,
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
    strand_spec : float
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
    denom_ss = 2.0 * strand_spec - 1.0
    inv_var_weight = denom_ss ** 2 if denom_ss > _STRAND_DENOM_EPS else 0.0

    if inv_var_weight > 0:
        g = 2.0 * (antisense * strand_spec - sense * (1.0 - strand_spec)) / denom_ss
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
    rate = inv_var_weight * strand_rate + (1.0 - inv_var_weight) * density_rate
    return min(max(rate, 0.0), 1.0), total


def _compute_chrom_gdna_rates(
    t_sense: np.ndarray,
    t_anti: np.ndarray,
    t_exonic: np.ndarray,
    t_refs: np.ndarray,
    strand_spec: float,
    intergenic_density: float,
    global_rate: float,
    k_chrom: float,
) -> dict[str, float]:
    """Compute per-chromosome gDNA rates shrunk toward global.

    Returns a dict mapping chromosome name to the shrunk gDNA rate.
    """
    chrom_sense: dict[str, float] = defaultdict(float)
    chrom_anti: dict[str, float] = defaultdict(float)
    chrom_exonic: dict[str, float] = defaultdict(float)
    for t_idx in range(len(t_sense)):
        ref = str(t_refs[t_idx])
        chrom_sense[ref] += t_sense[t_idx]
        chrom_anti[ref] += t_anti[t_idx]
        chrom_exonic[ref] += t_exonic[t_idx]

    chrom_shrunk: dict[str, float] = {}
    for ref in chrom_sense:
        rate, evidence = compute_gdna_rate_hybrid(
            chrom_sense[ref], chrom_anti[ref], chrom_exonic[ref],
            strand_spec, intergenic_density,
        )
        shrink_weight = evidence / (evidence + k_chrom) if (evidence + k_chrom) > 0 else 0.0
        chrom_shrunk[ref] = shrink_weight * rate + (1.0 - shrink_weight) * global_rate

    return chrom_shrunk


def _compute_per_locus_gdna_rates(
    loci: list,
    t_sense: np.ndarray,
    t_anti: np.ndarray,
    t_exonic: np.ndarray,
    t_refs: np.ndarray,
    strand_spec: float,
    intergenic_density: float,
    chrom_shrunk: dict[str, float],
    global_rate: float,
) -> tuple[list[float], list[float], list[float]]:
    """Compute per-locus gDNA rates and their parent (chromosome) rates.

    Returns (locus_rates, locus_evidence, locus_parents).
    """
    locus_rates: list[float] = []
    locus_evidence: list[float] = []
    locus_parents: list[float] = []

    for locus in loci:
        t_arr = locus.transcript_indices
        locus_sense = float(t_sense[t_arr].sum())
        locus_anti = float(t_anti[t_arr].sum())
        locus_exonic_bp = float(t_exonic[t_arr].sum())

        locus_rate, locus_n = compute_gdna_rate_hybrid(
            locus_sense, locus_anti, locus_exonic_bp, strand_spec, intergenic_density,
        )

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

        locus_rates.append(locus_rate)
        locus_evidence.append(locus_n)
        locus_parents.append(chrom_shrunk.get(primary_ref, global_rate))

    return locus_rates, locus_evidence, locus_parents


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
    kappa_min: float = EMConfig().nrna_frac_kappa_min,
    kappa_max: float = EMConfig().nrna_frac_kappa_max,
    kappa_fallback: float = EMConfig().nrna_frac_kappa_fallback,
    kappa_min_obs: int = EMConfig().nrna_frac_kappa_min_obs,
) -> list[float]:
    """Compute empirical Bayes gDNA initialization per locus.

    Hierarchical weighted shrinkage: locus → chromosome → global.

    Uses a hybrid density + strand signal for gDNA rate estimation
    at each hierarchical level (see :func:`compute_gdna_rate_hybrid`).
    When ``kappa_chrom`` or ``kappa_locus`` is ``None``, the shrinkage
    concentration is auto-estimated via Method of Moments using
    :func:`~rigel.estimator.estimate_kappa`.

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
    strand_spec = strand_models.strand_specificity
    # Fan out per-nRNA unspliced counts to per-transcript for gDNA EB.
    # Divide by fan count so that summing across transcripts sharing an
    # nRNA gives the original per-nRNA total (no double-counting).
    if estimator._t_to_nrna is not None:
        _t2n = estimator._t_to_nrna
        _fan = np.bincount(_t2n, minlength=estimator.num_nrna).astype(np.float64)
        _fan = np.maximum(_fan, 1.0)
        t_sense = estimator.transcript_unspliced_sense[_t2n] / _fan[_t2n]
        t_anti = estimator.transcript_unspliced_antisense[_t2n] / _fan[_t2n]
    else:
        t_sense = estimator.transcript_unspliced_sense
        t_anti = estimator.transcript_unspliced_antisense
    t_refs = index.t_df["ref"].values
    t_exonic = index.t_df["length"].values.astype(np.float64)

    # --- Global level (hybrid density+strand) ---
    total_sense = float(t_sense.sum())
    total_anti = float(t_anti.sum())
    total_exonic_bp = float(t_exonic.sum())
    global_rate, global_n = compute_gdna_rate_hybrid(
        total_sense, total_anti, total_exonic_bp, strand_spec, intergenic_density,
    )

    # --- Chromosome level: κ estimation ---
    _mom_kw = dict(
        kappa_min=kappa_min, kappa_max=kappa_max,
        kappa_fallback=kappa_fallback, kappa_min_obs=kappa_min_obs,
    )
    if kappa_chrom is None:
        # Need raw chrom rates for MoM estimation before shrinking
        chrom_sense: dict[str, float] = defaultdict(float)
        chrom_anti: dict[str, float] = defaultdict(float)
        chrom_exonic: dict[str, float] = defaultdict(float)
        for t_idx in range(len(t_sense)):
            ref = str(t_refs[t_idx])
            chrom_sense[ref] += t_sense[t_idx]
            chrom_anti[ref] += t_anti[t_idx]
            chrom_exonic[ref] += t_exonic[t_idx]
        chrom_rate_arr = np.empty(len(chrom_sense), dtype=np.float64)
        chrom_n_arr = np.empty(len(chrom_sense), dtype=np.float64)
        for i, ref in enumerate(chrom_sense):
            rate, evidence = compute_gdna_rate_hybrid(
                chrom_sense[ref], chrom_anti[ref], chrom_exonic[ref],
                strand_spec, intergenic_density,
            )
            chrom_rate_arr[i] = rate
            chrom_n_arr[i] = evidence
        k_chrom = estimate_kappa(
            chrom_rate_arr, chrom_n_arr, mom_min_evidence_chrom,
            **_mom_kw,
        )
        logger.debug(f"gDNA κ_chrom (MoM auto): {k_chrom:.1f}")
        n_chroms = len(chrom_sense)
    else:
        k_chrom = kappa_chrom
        n_chroms = 0  # not computed

    logger.info(
        f"gDNA EB: global_rate={global_rate:.4f} "
        f"(N={global_n:.0f}), κ_chrom={k_chrom:.1f}"
    )

    # --- Chromosome level: shrink toward global ---
    chrom_shrunk = _compute_chrom_gdna_rates(
        t_sense, t_anti, t_exonic, t_refs,
        strand_spec, intergenic_density, global_rate, k_chrom,
    )

    # --- Per-locus level ---
    locus_rates, locus_evidence, locus_parents = _compute_per_locus_gdna_rates(
        loci, t_sense, t_anti, t_exonic, t_refs,
        strand_spec, intergenic_density, chrom_shrunk, global_rate,
    )

    # MoM κ for locus → chrom shrinkage
    if kappa_locus is None:
        locus_rate_arr = np.array(locus_rates, dtype=np.float64)
        locus_n_arr = np.array(locus_evidence, dtype=np.float64)
        k_locus = estimate_kappa(
            locus_rate_arr, locus_n_arr, mom_min_evidence_locus,
            **_mom_kw,
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

        shrink_weight = locus_n / (locus_n + k_locus) if (locus_n + k_locus) > 0 else 0.0
        shrunk_rate = shrink_weight * locus_rate + (1.0 - shrink_weight) * parent_rate

        # gDNA init = shrunk_rate × unspliced fragment count in locus
        n_unspliced = 0
        for u in locus.unit_indices:
            if not em_data.is_spliced[u]:
                n_unspliced += 1

        gdna_init = shrunk_rate * n_unspliced
        gdna_inits.append(max(gdna_init, 0.0))

    return gdna_inits
