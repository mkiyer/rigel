"""rigel.locus_partition — Per-locus CSR scatter for the per-locus EM solver.

Scatters the monolithic ``ScoredFragments`` global CSR into per-locus
``LocusPartition`` objects, freeing each global array immediately after
scatter to bound peak memory.

Unrelated to region/annotation partitioning. The name reflects what the
module actually does: it builds the input arrays the per-locus EM
requires (one ``LocusPartition`` per connected-component locus).
"""

import numpy as np

from .native import (
    build_partition_offsets,
    scatter_candidates_f32,
    scatter_candidates_f64,
    scatter_candidates_i32,
    scatter_candidates_u8,
    scatter_units_f32,
    scatter_units_f64,
    scatter_units_i32,
    scatter_units_i64,
    scatter_units_u8,
)
from .scored_fragments import LocusPartition, ScoredFragments


def _float_candidate_scatter(arr: np.ndarray):
    if arr.dtype == np.float32:
        return scatter_candidates_f32
    if arr.dtype == np.float64:
        return scatter_candidates_f64
    raise TypeError(f"Expected float32 or float64 candidate payload, got {arr.dtype}")


def _float_unit_scatter(arr: np.ndarray):
    if arr.dtype == np.float32:
        return scatter_units_f32
    if arr.dtype == np.float64:
        return scatter_units_f64
    raise TypeError(f"Expected float32 or float64 unit payload, got {arr.dtype}")


# Dtype-driven selectors: a table entry whose scatter_fn is one of these is
# resolved to its concrete scatter by calling it with the global array (no
# by-name special-casing — a new float array routes correctly automatically).
_FLOAT_SELECTORS = (_float_candidate_scatter, _float_unit_scatter)


def partition_and_free(
    em_data: ScoredFragments,
    multi_loci: list,
) -> dict[int, LocusPartition]:
    """Scatter global CSR into per-locus partitions, freeing each global
    array immediately after scatter.

    Parameters
    ----------
    em_data : ScoredFragments
        Global CSR arrays. All data arrays are set to None after scatter.
    multi_loci : list[MultiLocus]
        MultiLocus objects with ``unit_indices`` arrays.

    Returns
    -------
    dict[int, LocusPartition]
        Mapping from multi-locus index to LocusPartition.
    """
    n_loci = len(multi_loci)
    locus_units = [locus.unit_indices for locus in multi_loci]

    # ---- Build per-locus CSR offsets ----
    offsets_list = build_partition_offsets(em_data.offsets, locus_units, n_loci)

    # ---- Scatter per-candidate arrays (largest first) ----
    # g_offsets must stay alive during this phase.
    CAND_ARRAYS = [
        ("log_liks", _float_candidate_scatter),
        ("coverage_weights", _float_candidate_scatter),
        ("t_indices", scatter_candidates_i32),
        ("count_cols", scatter_candidates_u8),
    ]
    cand_results = {}
    for attr, scatter_fn in CAND_ARRAYS:
        global_arr = getattr(em_data, attr)
        if scatter_fn in _FLOAT_SELECTORS:
            scatter_fn = scatter_fn(global_arr)
        cand_results[attr] = scatter_fn(
            global_arr, em_data.offsets, locus_units, offsets_list, n_loci
        )
        setattr(em_data, attr, None)
        del global_arr

    # g_offsets no longer needed
    em_data.offsets = None

    # ---- Scatter per-unit arrays ----
    UNIT_ARRAYS = [
        ("gdna_log_liks", _float_unit_scatter),
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
        elif scatter_fn in _FLOAT_SELECTORS:
            scatter_fn = scatter_fn(global_arr)
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
            is_spliced=unit_results["is_spliced"][li],
            gdna_log_liks=unit_results["gdna_log_liks"][li],
            locus_t_indices=unit_results["locus_t_indices"][li],
            locus_count_cols=unit_results["locus_count_cols"][li],
            frag_ids=unit_results["frag_ids"][li],
            frag_class=unit_results["frag_class"][li],
            splice_type=unit_results["splice_type"][li],
        )

    return partitions
