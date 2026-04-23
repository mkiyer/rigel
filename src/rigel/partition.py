"""rigel.partition — Array-by-array scatter from global CSR to per-locus partitions.

Scatters the monolithic ``ScoredFragments`` global CSR into per-locus
``LocusPartition`` objects, freeing each global array immediately after
scatter to bound peak memory.
"""

import numpy as np

from .native import (
    build_partition_offsets,
    scatter_candidates_f64,
    scatter_candidates_i32,
    scatter_candidates_u8,
    scatter_units_f64,
    scatter_units_i32,
    scatter_units_i64,
    scatter_units_u8,
)
from .scored_fragments import LocusPartition, ScoredFragments


def partition_and_free(
    em_data: ScoredFragments,
    loci: list,
) -> dict[int, LocusPartition]:
    """Scatter global CSR into per-locus partitions, freeing each global
    array immediately after scatter.

    Parameters
    ----------
    em_data : ScoredFragments
        Global CSR arrays. All data arrays are set to None after scatter.
    loci : list[Locus]
        Locus objects with ``unit_indices`` arrays.

    Returns
    -------
    dict[int, LocusPartition]
        Mapping from locus index to LocusPartition.
    """
    n_loci = len(loci)
    locus_units = [locus.unit_indices for locus in loci]

    # ---- Build per-locus CSR offsets ----
    offsets_list = build_partition_offsets(
        em_data.offsets, locus_units, n_loci
    )

    # ---- Scatter per-candidate arrays (largest first) ----
    # g_offsets must stay alive during this phase.
    CAND_ARRAYS = [
        ("log_liks", scatter_candidates_f64),
        ("coverage_weights", scatter_candidates_f64),
        ("t_indices", scatter_candidates_i32),
        ("tx_starts", scatter_candidates_i32),
        ("tx_ends", scatter_candidates_i32),
        ("count_cols", scatter_candidates_u8),
    ]
    cand_results = {}
    for attr, scatter_fn in CAND_ARRAYS:
        global_arr = getattr(em_data, attr)
        cand_results[attr] = scatter_fn(
            global_arr, em_data.offsets, locus_units, offsets_list, n_loci
        )
        setattr(em_data, attr, None)
        del global_arr

    # g_offsets no longer needed
    em_data.offsets = None

    # ---- Scatter per-unit arrays ----
    UNIT_ARRAYS = [
        ("gdna_log_liks", scatter_units_f64),
        ("locus_t_indices", scatter_units_i32),
        ("locus_count_cols", scatter_units_u8),
        ("is_spliced", scatter_units_u8),
        ("frag_ids", scatter_units_i64),
        ("frag_class", scatter_units_u8),
        ("splice_type", scatter_units_u8),
    ]
    unit_results = {}
    for attr, scatter_fn in UNIT_ARRAYS:
        global_arr = getattr(em_data, attr)
        if global_arr.dtype == np.bool_:
            global_arr = global_arr.view(np.uint8)
        elif global_arr.dtype == np.int8:
            global_arr = global_arr.view(np.uint8)
        unit_results[attr] = scatter_fn(global_arr, locus_units, n_loci)
        setattr(em_data, attr, None)
        del global_arr

    # ---- Assemble LocusPartition objects ----
    partitions = {}
    for li in range(n_loci):
        partitions[li] = LocusPartition(
            locus_id=li,
            n_units=len(offsets_list[li]) - 1,
            n_candidates=int(offsets_list[li][-1]),
            offsets=offsets_list[li],
            t_indices=cand_results["t_indices"][li],
            log_liks=cand_results["log_liks"][li],
            count_cols=cand_results["count_cols"][li],
            coverage_weights=cand_results["coverage_weights"][li],
            tx_starts=cand_results["tx_starts"][li],
            tx_ends=cand_results["tx_ends"][li],
            is_spliced=unit_results["is_spliced"][li],
            gdna_log_liks=unit_results["gdna_log_liks"][li],
            locus_t_indices=unit_results["locus_t_indices"][li],
            locus_count_cols=unit_results["locus_count_cols"][li],
            frag_ids=unit_results["frag_ids"][li],
            frag_class=unit_results["frag_class"][li],
            splice_type=unit_results["splice_type"][li],
        )

    return partitions
