"""The locus-EM harness the estimator tests drive.

`_make_locus_em_data` builds the per-locus arrays the native EM consumes from a compact description of
units and their count columns; `_run_and_assign` runs one EM and returns the totals it assigned; `em_partition` is
one locus's partition as the EM reads it, its arrays by name. A plain harness
rather than a fixture, so a test composes it with whatever index and priors it needs, and a module of its own,
because importing `conftest` by name is ambiguous once a sub-directory has one.
"""

from __future__ import annotations

import dataclasses
from types import SimpleNamespace

import numpy as np

from rigel.locus import Locus, MultiLocus
from rigel.scored_fragments import ScoredFragments
from rigel.splice import SpliceStrandCol

# Default column for UNSPLICED_SENSE
_UNSPLICED_SENSE = int(SpliceStrandCol.UNSPLICED_SENSE)


def em_partition(
    offsets,
    t_indices,
    log_liks,
    coverage_weights,
    count_cols,
    gdna_excluded,
    gdna_log_liks,
    intergenic_log_liks,
):
    """One locus's partition as ``run_batch_locus_em_partitioned`` reads it: each array by its name, as a
    ``LocusPartition`` carries them (the candidates' CSR, then the units')."""
    return SimpleNamespace(
        offsets=offsets,
        t_indices=t_indices,
        log_liks=log_liks,
        coverage_weights=coverage_weights,
        count_cols=count_cols,
        gdna_excluded=gdna_excluded,
        gdna_log_liks=gdna_log_liks,
        intergenic_log_liks=intergenic_log_liks,
    )


def _make_locus_em_data(
    t_indices_per_unit,
    log_liks_per_unit=None,
    count_cols_per_unit=None,
    num_transcripts=None,
    rc=None,
    include_nascent=False,
    include_gdna=False,
    nascent_log_lik=-2.0,
    gdna_log_lik=0.0,
    gdna_prior_count=1.0,
):
    """Build (ScoredFragments, [Locus], gdna_prior_count) for batch EM tests.

    Returns a tuple (em_data, loci, gdna_prior_count_arr) that ``_run_and_assign``
    partitions and hands to ``AbundanceEstimator.run_batch_locus_em_partitioned``.

    Parameters
    ----------
    t_indices_per_unit : list[list[int]]
        Global transcript indices per unit.
    rc : AbundanceEstimator or None
        If provided, used for unambig_counts.
    include_nascent : bool
        If True, units are marked as unspliced (``gdna_excluded=False``).
        Nothing is appended for nascent RNA: synthetic spans are ordinary
        transcripts, listed among a unit's candidates like any other.
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

    # gdna_excluded: True (spliced) → no gDNA candidate.
    # For include_nascent or include_gdna, set False (unspliced).
    if include_nascent or include_gdna:
        gdna_excluded = np.zeros(n_units, dtype=bool)
    else:
        gdna_excluded = np.ones(n_units, dtype=bool)

    # gDNA log-likelihoods per unit
    if include_gdna:
        gdna_log_liks = np.full(n_units, gdna_log_lik, dtype=np.float32)
    else:
        gdna_log_liks = np.full(n_units, -np.inf, dtype=np.float32)

    em_data = ScoredFragments(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.array(flat_lk, dtype=np.float32),
        count_cols=np.array(flat_cc, dtype=np.uint8),
        coverage_weights=np.ones(n_candidates, dtype=np.float32),
        gdna_excluded=gdna_excluded,
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
        span_bp=10000,
        loci=(Locus(ref="chr1", ref_id=0, start=0, end=10000),),
    )
    # as production builds them (`locus.build_multi_loci`), a scored set with no unit makes no locus
    loci = [locus] if n_units else []

    gdna_prior_count_arr = np.array(
        [gdna_prior_count if include_gdna else 0.0] * len(loci), dtype=np.float64
    )

    return em_data, loci, gdna_prior_count_arr


def _run_and_assign(rc, em_data, loci=None, gdna_prior_count=None, *, em_iterations=10):
    """Run batch locus EM via the partitioned path, with *rc*'s EM budget set to ``em_iterations`` EM steps. Returns
    the totals assigned: ``annotated`` (every transcript's EM count), ``nascent`` (the synthetic spans') and ``gdna``.

    Accepts either the tuple form (em_data, loci, gdna_prior_count) separately,
    or the tuple returned by ``_make_locus_em_data`` as ``em_data``.
    """
    from rigel.locus_partition import partition_and_free

    # Unpack tuple form from _make_locus_em_data
    if isinstance(em_data, tuple):
        em_data, loci, gdna_prior_count = em_data

    # Partition ScoredFragments into per-locus LocusPartition objects
    partitions = partition_and_free(em_data, loci)

    locus_t_lists = [loc.transcript_indices for loc in loci]

    rc.em_config = dataclasses.replace(rc.em_config, iterations=em_iterations)
    total_gdna, _locus_rna, _locus_gdna = rc.run_batch_locus_em_partitioned(
        [partitions[i] for i in range(len(loci))], locus_t_lists, gdna_prior_count
    )
    rc._gdna_em_total += total_gdna

    return {
        "annotated": float(rc.em_counts.sum()),
        "nascent": float(rc.synthetic_count),
        "gdna": float(total_gdna),
    }
