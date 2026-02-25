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

from .estimator import AbundanceEstimator, ScanData, Locus, LocusEMInput
from .index import HulkIndex
from .strand_model import StrandModels

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Empirical Bayes hyperparameters
# ---------------------------------------------------------------------------
#: Fragments needed for locus-level to dominate.
EB_K_LOCUS = 20.0
#: Fragments needed for chrom-level to dominate.
EB_K_CHROM = 50.0


# ---------------------------------------------------------------------------
# Locus builder: scipy connected_components on fragment→transcript graph
# ---------------------------------------------------------------------------


def build_loci(
    em_data: ScanData,
    index: HulkIndex,
) -> list[Locus]:
    """Build loci as connected components of transcripts linked by fragments.

    Uses scipy sparse graph connected_components for vectorized
    component detection.

    Parameters
    ----------
    em_data : ScanData
        Global EM data (mRNA + nRNA candidates + per-unit metadata).
    index : HulkIndex
        Reference index.

    Returns
    -------
    list[Locus]
    """
    nt = index.num_transcripts
    nrna_base = em_data.nrna_base_index
    offsets = em_data.offsets
    t_indices = em_data.t_indices
    t_to_g = index.t_to_g_arr
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
        g_set = sorted(set(int(t_to_g[t]) for t in t_arr))
        g_arr = np.array(g_set, dtype=np.int32)
        u_arr = np.array(sorted(root_to_units.get(lbl, [])), dtype=np.int32)

        loci.append(Locus(
            locus_id=lid,
            transcript_indices=t_arr,
            gene_indices=g_arr,
            unit_indices=u_arr,
        ))

    return loci


# ---------------------------------------------------------------------------
# Locus EM data builder
# ---------------------------------------------------------------------------


def build_locus_em_data(
    locus: Locus,
    em_data: ScanData,
    counter: AbundanceEstimator,
    index: HulkIndex,
    mean_frag: float,
    gdna_init: float,
) -> LocusEMInput:
    """Extract and renumber global ScanData into a per-locus sub-problem.

    Component layout per locus::

        [0, n_t)           - mRNA (one per local transcript)
        [n_t, 2*n_t)       - nRNA (one per local transcript)
        [2*n_t]            - gDNA (ONE shadow for entire locus)

    Only UNSPLICED units get a gDNA candidate.  Spliced units see
    only mRNA/nRNA components.

    Parameters
    ----------
    locus : Locus
    em_data : ScanData
    counter : AbundanceEstimator
    index : HulkIndex
    mean_frag : float
    gdna_init : float
        Empirical Bayes estimated gDNA count for this locus.
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
        dedup_unit = v_unit[order][first_mask]
    else:
        dedup_lidx = np.empty(0, dtype=np.int32)
        dedup_ll = np.empty(0, dtype=np.float64)
        dedup_cc = np.empty(0, dtype=np.uint8)
        dedup_unit = np.empty(0, dtype=np.int32)

    # Add gDNA candidates for unspliced units
    is_spl = em_data.is_spliced[locus.unit_indices]
    gdna_lls = em_data.gdna_log_liks[locus.unit_indices]
    valid_gdna = (~is_spl) & np.isfinite(gdna_lls)
    n_gdna = int(valid_gdna.sum())

    # Compute true locus span for gDNA per-fragment effective length.
    # Span extends from leftmost transcript/fragment start to rightmost
    # transcript/fragment end, capturing any overhanging fragments.
    t_starts = index.t_df["start"].values[t_arr]
    t_ends = index.t_df["end"].values[t_arr]
    locus_start = int(t_starts.min())
    locus_end = int(t_ends.max())
    footprints = em_data.genomic_footprints[locus.unit_indices]
    locus_span = float(locus_end - locus_start)

    if n_gdna > 0:
        gdna_units = np.arange(n_local_units, dtype=np.int32)[valid_gdna]
        gdna_lidx_arr = np.full(n_gdna, gdna_idx, dtype=np.int32)
        gdna_ll_arr = gdna_lls[valid_gdna].copy()
        gdna_cc_arr = np.zeros(n_gdna, dtype=np.uint8)

        # Per-fragment effective length correction for gDNA.
        # Subtract log(max(locus_span - footprint + 1, 1)) from each
        # gDNA candidate's log-likelihood.
        gdna_footprints = footprints[valid_gdna].astype(np.float64)
        gdna_per_frag_eff = np.maximum(
            locus_span - gdna_footprints + 1.0, 1.0
        )
        gdna_ll_arr -= np.log(gdna_per_frag_eff)

        final_lidx = np.concatenate([dedup_lidx, gdna_lidx_arr])
        final_ll = np.concatenate([dedup_ll, gdna_ll_arr])
        final_cc = np.concatenate([dedup_cc, gdna_cc_arr])
        final_unit = np.concatenate([dedup_unit, gdna_units])
    else:
        final_lidx = dedup_lidx
        final_ll = dedup_ll
        final_cc = dedup_cc
        final_unit = dedup_unit

    # Sort by unit to reconstruct CSR
    if len(final_unit) > 0:
        sort_order = np.argsort(final_unit, kind='stable')
        final_lidx = final_lidx[sort_order]
        final_ll = final_ll[sort_order]
        final_cc = final_cc[sort_order]
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
    unique_totals = np.zeros(n_components, dtype=np.float64)
    unique_totals[:n_t] = counter.unique_counts[t_arr].sum(axis=1)
    unique_totals[n_t:2*n_t] = counter.nrna_init[t_arr]
    unique_totals[gdna_idx] = gdna_init

    # Build effective lengths — all set to 1.0 because per-fragment
    # effective length correction is now baked into each candidate's
    # log-likelihood (mRNA/nRNA during scan, gDNA above).
    eff_len = np.ones(n_components, dtype=np.float64)

    # Build prior
    prior = np.full(n_components, counter.em_prior, dtype=np.float64)

    # Zero gDNA prior when there are no unspliced fragments
    if gdna_init == 0.0:
        prior[gdna_idx] = 0.0

    # Zero nRNA prior for single-exon transcripts
    if counter._transcript_spans is not None:
        if counter._exonic_lengths is not None:
            t_exon = counter._exonic_lengths[t_arr]
        else:
            t_exon = counter._t_eff_len[t_arr] + mean_frag - 1.0
        single_exon = counter._transcript_spans[t_arr] <= t_exon
        prior[n_t:2*n_t][single_exon] = 0.0

    # Zero nRNA prior when nrna_init is zero for that transcript
    nrna_init_local = unique_totals[n_t:2*n_t]
    prior[n_t:2*n_t][nrna_init_local == 0.0] = 0.0

    return LocusEMInput(
        locus=locus,
        offsets=local_offsets,
        t_indices=final_lidx.astype(np.int32),
        log_liks=final_ll.astype(np.float64),
        count_cols=final_cc.astype(np.uint8),
        locus_t_indices=local_locus_t.astype(np.int32),
        locus_count_cols=local_locus_ct.astype(np.uint8),
        n_transcripts=n_t,
        n_components=n_components,
        local_to_global_t=t_arr.copy(),
        unique_totals=unique_totals,
        nrna_init=counter.nrna_init[t_arr].copy(),
        gdna_init=gdna_init,
        effective_lengths=eff_len,
        prior=prior,
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
    ss = strand_models.exonic_spliced.strand_specificity
    denom = 2.0 * ss - 1.0

    if denom <= 0.2:
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
    if denom_ss <= 0.2:
        return 0.0
    total = sense + antisense
    if total == 0:
        return 0.0
    g = 2.0 * (antisense * ss - sense * (1.0 - ss)) / denom_ss
    g = max(g, 0.0)
    rate = g / total
    return min(rate, 1.0)


def compute_eb_gdna_priors(
    loci: list[Locus],
    em_data: ScanData,
    counter: AbundanceEstimator,
    index: HulkIndex,
    strand_models: StrandModels,
    *,
    k_locus: float = EB_K_LOCUS,
    k_chrom: float = EB_K_CHROM,
) -> list[float]:
    """Compute empirical Bayes gDNA initialization per locus.

    Hierarchical weighted shrinkage: locus → chromosome → global.

    Parameters
    ----------
    loci : list[Locus]
    em_data : ScanData
    counter : AbundanceEstimator
    index : HulkIndex
    strand_models : StrandModels
    k_locus : float
        Shrinkage strength for locus → chrom (default 20).
    k_chrom : float
        Shrinkage strength for chrom → global (default 50).

    Returns
    -------
    list[float]
        gDNA init count per locus (same order as ``loci``).
    """
    ss = strand_models.exonic_spliced.strand_specificity
    g_sense = counter.gene_sense_all
    g_anti = counter.gene_antisense_all
    g_refs = index.g_df["ref"].values

    # --- Global level ---
    total_sense = float(g_sense.sum())
    total_anti = float(g_anti.sum())
    global_rate = compute_gdna_rate_from_strand(
        total_sense, total_anti, ss,
    )
    global_n = total_sense + total_anti

    # --- Chromosome level ---
    chrom_sense: dict[str, float] = defaultdict(float)
    chrom_anti: dict[str, float] = defaultdict(float)
    for g_idx in range(len(g_sense)):
        ref = str(g_refs[g_idx])
        chrom_sense[ref] += g_sense[g_idx]
        chrom_anti[ref] += g_anti[g_idx]

    chrom_rate: dict[str, float] = {}
    chrom_n: dict[str, float] = {}
    for ref in chrom_sense:
        cs = chrom_sense[ref]
        ca = chrom_anti[ref]
        chrom_rate[ref] = compute_gdna_rate_from_strand(cs, ca, ss)
        chrom_n[ref] = cs + ca

    # Shrink chrom toward global
    chrom_shrunk: dict[str, float] = {}
    for ref in chrom_rate:
        n = chrom_n[ref]
        w = n / (n + k_chrom) if (n + k_chrom) > 0 else 0.0
        chrom_shrunk[ref] = w * chrom_rate[ref] + (1.0 - w) * global_rate

    # --- Per-locus level ---
    gdna_inits: list[float] = []
    for locus in loci:
        g_arr = locus.gene_indices
        locus_sense = float(g_sense[g_arr].sum())
        locus_anti = float(g_anti[g_arr].sum())
        locus_n = locus_sense + locus_anti

        locus_rate = compute_gdna_rate_from_strand(
            locus_sense, locus_anti, ss,
        )

        # Determine primary chromosome for this locus
        ref_counts: dict[str, float] = defaultdict(float)
        for g_idx in g_arr:
            ref = str(g_refs[int(g_idx)])
            ref_counts[ref] += g_sense[int(g_idx)] + g_anti[int(g_idx)]
        if ref_counts:
            primary_ref = max(ref_counts, key=ref_counts.get)
        else:
            if len(locus.transcript_indices) > 0:
                primary_ref = str(
                    index.t_df["ref"].values[
                        int(locus.transcript_indices[0])
                    ]
                )
            else:
                primary_ref = ""

        parent_rate = chrom_shrunk.get(primary_ref, global_rate)

        # Shrink locus toward parent (chromosome)
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
