"""rigel.locus — Locus graph construction and EM initialization.

Handles everything between "buffer scan complete" and "per-locus EM
loop": connected-component partitioning, per-locus EM data extraction,
and Empirical Bayes gDNA priors.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from .native import connected_components as _cc_native

from .scored_fragments import Locus, LocusEMInput, ScoredFragments

from .native import EM_PRIOR_EPSILON
from .index import TranscriptIndex

if TYPE_CHECKING:
    from .calibration import GDNACalibration
    from .estimator import AbundanceEstimator


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
    # to_numpy(dtype=object) converts the Arrow-backed string column
    # to a plain object array, avoiding repeated pyarrow.compute.take
    # overhead in the per-locus loop below.
    t_ref_strs = index.t_df["ref"].to_numpy(dtype=object, na_value="")
    t_starts_all = index.t_df["start"].values
    t_ends_all = index.t_df["end"].values

    # Integer ref codes for fast sort/compare within merge loops
    _ref_names, _ref_codes = np.unique(t_ref_strs, return_inverse=True)

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
# Locus EM data builder
# ---------------------------------------------------------------------------


def build_locus_em_data(
    locus: Locus,
    em_data: ScoredFragments,
    estimator: AbundanceEstimator,
    index: TranscriptIndex,
    gdna_init: float,
    *,
    _cache: dict | None = None,
) -> LocusEMInput:
    """Extract and renumber global ScoredFragments into a per-locus sub-problem.

    Component layout per locus::

        [0, n_t)                - mRNA (one per local transcript)
        [n_t]          - gdna (single gDNA component)

    Only UNSPLICED units get a gDNA candidate.
    Spliced units see only transcript components.

    Parameters
    ----------
    locus : Locus
    em_data : ScoredFragments
    estimator : AbundanceEstimator
    index : TranscriptIndex
    gdna_init : float
        Empirical Bayes estimated gDNA count for this locus.
    _cache : dict or None
        Optional pre-extracted arrays to avoid repeated DataFrame
        access.  Expected keys: ``"t_lengths"``, ``"local_map"``
        (reusable scratch buffer).
    """
    t_arr = locus.transcript_indices
    n_t = len(t_arr)
    n_local_units = len(locus.unit_indices)

    gdna_idx = n_t  # single gDNA component index
    n_components = n_t + 1

    # Build global → local mapping array (fast lookup, no dict)
    # transcript: global_t → local_t [0, n_t)
    max_global = int(t_arr.max()) + 1 if n_t > 0 else 1

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
        local_map[gt] = local_i

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
        unit_of_cand = np.repeat(np.arange(n_local_units, dtype=np.int32), seg_lens)

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

    # Use pre-computed union genomic footprint for gDNA bias profile.
    footprints = em_data.genomic_footprints[locus.unit_indices]
    locus_span = float(locus.gdna_span)

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
        sort_order = np.argsort(final_unit, kind="stable")
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
    # so they are disjoint from em_totals.
    unambig_totals = np.zeros(n_components, dtype=np.float64)
    unambig_totals[:n_t] = estimator.unambig_counts[t_arr].sum(axis=1)
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
    #   [0, n_t)     → transcript: length = exonic transcript length
    #   [n_t]        → gdna: length = locus span
    _t_lengths_all = _cache["t_lengths"] if _cache is not None else index.t_df["length"].values
    mrna_lengths = _t_lengths_all[t_arr]
    bias_profiles = np.empty(n_components, dtype=np.int64)
    bias_profiles[:n_t] = mrna_lengths
    bias_profiles[gdna_idx] = int(locus_span)

    # Build prior — start at epsilon (numerical floor) for every
    # component; the real prior mass comes from the OVR coverage
    # weighting applied later in run_locus_em().
    prior = np.full(n_components, EM_PRIOR_EPSILON, dtype=np.float64)

    # Zero gDNA prior when there are no unspliced fragments
    if gdna_init == 0.0:
        prior[gdna_idx] = 0.0

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
        gdna_init=gdna_init,
        effective_lengths=eff_len,
        prior=prior,
        bias_profiles=bias_profiles,
    )


# ---------------------------------------------------------------------------
# Calibration γ aggregation per locus
# ---------------------------------------------------------------------------


def compute_gdna_locus_gammas(
    loci: list[Locus],
    index: TranscriptIndex,
    calibration: "GDNACalibration",
) -> np.ndarray:
    """Aggregate calibration region posteriors per locus.

    For each locus, finds overlapping calibration regions via
    ``index.region_cr`` (cgranges) and computes a fragment-weighted
    aggregate gDNA posterior:

        γ_locus = Σ(γ_r × n_r) / Σ(n_r)

    Falls back to ``calibration.mixing_proportion`` (global π) when no
    regions overlap a locus.

    Returns
    -------
    locus_gammas : np.ndarray, shape (n_loci,), float64 ∈ [0, 1]
    """
    n_loci = len(loci)
    locus_gammas = np.empty(n_loci, dtype=np.float64)

    region_gamma = calibration.region_posteriors
    region_n = calibration.region_n_total
    fallback = calibration.mixing_proportion

    # Use cgranges for spatial overlap if available
    region_cr = getattr(index, "region_cr", None)
    if region_cr is None or region_gamma is None or region_n is None:
        locus_gammas[:] = fallback
        return locus_gammas

    for li, locus in enumerate(loci):
        wsum = 0.0
        nsum = 0.0
        for ref, start, end in locus.merged_intervals:
            for _s, _e, rid in region_cr.overlap(ref, start, end):
                n_r = float(region_n[rid])
                wsum += float(region_gamma[rid]) * n_r
                nsum += n_r

        locus_gammas[li] = wsum / nsum if nsum > 0.0 else fallback

    return locus_gammas
