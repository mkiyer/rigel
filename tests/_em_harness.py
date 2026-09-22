"""The locus-EM harness the estimator tests drive.

`_make_locus_em_data` builds the per-locus arrays the native EM consumes from a compact description of
units and their count columns; `_run_and_assign` runs one EM and returns the mRNA / nRNA / gDNA totals it
assigned. A plain harness rather than a fixture, so a test composes it with whatever index and priors it
needs — and a module of its own, because importing `conftest` by name is ambiguous once a sub-directory
has one.
"""

from __future__ import annotations

import numpy as np

from rigel.locus import Locus, MultiLocus
from rigel.scored_fragments import ScoredFragments
from rigel.splice import SpliceStrandCol

# Default column for UNSPLICED_SENSE
_UNSPLICED_SENSE = int(SpliceStrandCol.UNSPLICED_SENSE)


def _make_locus_em_data(
    t_indices_per_unit,
    log_liks_per_unit=None,
    count_cols_per_unit=None,
    num_transcripts=None,
    rc=None,
    include_nrna=False,
    include_gdna=False,
    nrna_log_lik=-2.0,
    gdna_log_lik=0.0,
    gdna_prior_count=1.0,
):
    """Build (ScoredFragments, [Locus], gdna_prior_count) for batch EM tests.

    Returns a tuple (em_data, loci, gdna_prior_count_arr) that ``_run_and_assign``
    partitions and hands to ``AbundanceEstimator.run_batch_locus_em_partitioned``.

    Parameters
    ----------
    t_indices_per_unit : list[list[int]]
        GLOBAL mRNA transcript indices per unit.
    rc : AbundanceEstimator or None
        If provided, used for unambig_counts.
    include_nrna : bool
        If True, units are marked as unspliced (``is_spliced=False``).
        Nothing is appended for nascent RNA: synthetic nRNA transcripts are
        ordinary transcripts, listed among a unit's candidates like any other.
    include_gdna : bool
        If True, units are unspliced with finite ``gdna_log_liks``, so the
        locus EM appends a gDNA candidate to each, and ``gdna_prior_count``
        is passed through (otherwise the prior count is 0).
    """
    if num_transcripts is None:
        all_t = [t for unit in t_indices_per_unit for t in unit]
        num_transcripts = (max(all_t) + 1) if all_t else 1

    n_t = num_transcripts
    n_units = len(t_indices_per_unit)

    # Build the transcript-candidate CSR (no gDNA — the locus EM appends that per unit)
    offsets = [0]
    flat_t = []
    flat_lk = []
    flat_cc = []

    for u, t_list in enumerate(t_indices_per_unit):
        for j, t_idx in enumerate(t_list):
            flat_t.append(t_idx)
            flat_lk.append(log_liks_per_unit[u][j] if log_liks_per_unit else 0.0)
            flat_cc.append(count_cols_per_unit[u][j] if count_cols_per_unit else _UNSPLICED_SENSE)
        offsets.append(len(flat_t))

    n_candidates = len(flat_t)

    # Per-unit locus tracking
    locus_t = np.full(n_units, -1, dtype=np.int32)
    locus_cc = np.zeros(n_units, dtype=np.uint8)
    for u, t_list in enumerate(t_indices_per_unit):
        if t_list:
            locus_t[u] = t_list[0]
            locus_cc[u] = count_cols_per_unit[u][0] if count_cols_per_unit else _UNSPLICED_SENSE

    # is_spliced: True (spliced) → no gDNA candidate.
    # For include_nrna or include_gdna, set False (unspliced).
    if include_nrna or include_gdna:
        is_spliced = np.zeros(n_units, dtype=bool)
    else:
        is_spliced = np.ones(n_units, dtype=bool)

    # gDNA log-likelihoods per unit
    if include_gdna:
        gdna_log_liks = np.full(n_units, gdna_log_lik, dtype=np.float64)
    else:
        gdna_log_liks = np.full(n_units, -np.inf, dtype=np.float64)

    em_data = ScoredFragments(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.array(flat_lk, dtype=np.float64),
        count_cols=np.array(flat_cc, dtype=np.uint8),
        coverage_weights=np.ones(n_candidates, dtype=np.float64),
        locus_t_indices=locus_t,
        locus_count_cols=locus_cc,
        is_spliced=is_spliced,
        gdna_log_liks=gdna_log_liks,
        frag_ids=np.arange(n_units, dtype=np.int64),
        frag_class=np.zeros(n_units, dtype=np.int8),
        splice_type=np.zeros(n_units, dtype=np.uint8),
        n_units=n_units,
        n_candidates=n_candidates,
    )

    locus = MultiLocus(
        multi_locus_id=0,
        transcript_indices=np.arange(n_t, dtype=np.int32),
        unit_indices=np.arange(n_units, dtype=np.int32),
        gdna_span=10000,
        loci=(Locus(ref="chr1", ref_id=0, start=0, end=10000),),
    )
    loci = [locus]

    gdna_prior_count_arr = np.array([gdna_prior_count if include_gdna else 0.0], dtype=np.float64)

    return em_data, loci, gdna_prior_count_arr


def _run_and_assign(rc, em_data, loci=None, gdna_prior_count=None, *, em_iterations=10):
    """Run batch locus EM via the partitioned path. Returns pool_counts dict.

    Accepts either the tuple form (em_data, loci, gdna_prior_count) separately,
    or the tuple returned by ``_make_locus_em_data`` as ``em_data``.
    """
    from rigel.locus_partition import partition_and_free

    # Unpack tuple form from _make_locus_em_data
    if isinstance(em_data, tuple):
        em_data, loci, gdna_prior_count = em_data

    # Partition ScoredFragments into per-locus LocusPartition objects
    partitions = partition_and_free(em_data, loci)

    # Build the 9-tuples and transcript index lists expected by C++
    partition_tuples = [
        (
            p.offsets,
            p.t_indices,
            p.log_liks,
            p.coverage_weights,
            p.count_cols,
            p.is_spliced,
            p.gdna_log_liks,
            p.locus_t_indices,
            p.locus_count_cols,
        )
        for p in [partitions[i] for i in range(len(loci))]
    ]
    locus_t_lists = [loc.transcript_indices for loc in loci]

    total_gdna, _locus_mrna, _locus_gdna = rc.run_batch_locus_em_partitioned(
        partition_tuples,
        locus_t_lists,
        gdna_prior_count,
        em_iterations=em_iterations,
    )
    rc._gdna_em_total += total_gdna

    return {
        "mrna": float(rc.em_counts.sum()),
        "nrna": float(rc.nrna_em_count),
        "gdna": float(total_gdna),
    }
