"""rigel.scored_fragments — Data containers for the EM solver.

Pure dataclasses with no logic.  These are shared between scan.py
(producer), locus.py (consumer/builder), and estimator.py (EM driver).

- ``ScoredFragments`` — global CSR arrays linking fragment units to
  candidate transcripts with log-likelihoods.
- ``Locus`` — connected component of transcripts linked by shared
  fragments.
- ``LocusEMInput`` — per-locus sub-problem with locally-renumbered
  components, ready for the C++ EM solver.
"""

from dataclasses import dataclass

import numpy as np


# ======================================================================
# ScoredFragments — pre-computed CSR arrays for vectorized EM
# ======================================================================


@dataclass(slots=True)
class ScoredFragments:
    """Fragment-level candidate data from the BAM scan pass.

    Contains pre-computed CSR (compressed sparse row) arrays linking
    each ambiguous fragment unit to its candidate transcripts and
    their log-likelihoods.  Produced by ``FragmentRouter`` in scan.py
    and consumed during locus EM construction.

    The global ScoredFragments contains transcript candidates only
    (NO gDNA).  gDNA candidates are added per-locus during locus EM
    construction.

    Attributes
    ----------
    offsets : np.ndarray
        int64[n_units + 1] — CSR offsets into flat arrays.
    t_indices : np.ndarray
        int32[n_candidates] — candidate transcript indices.
    log_liks : np.ndarray
        float64[n_candidates] — log(P_strand × P_insert) per candidate.
    count_cols : np.ndarray
        uint8[n_candidates] — internal column index per candidate (0–7).
    coverage_weights : np.ndarray
        float64[n_candidates] — coverage weight per candidate.
    tx_starts : np.ndarray
        int32[n_candidates] — 0-based transcript-space start position.
    tx_ends : np.ndarray
        int32[n_candidates] — 0-based transcript-space end position
        (exclusive).
    locus_t_indices : np.ndarray
        int32[n_units] — best transcript index per unit.
    locus_count_cols : np.ndarray
        uint8[n_units] — count column for the locus transcript.
    is_spliced : np.ndarray
        bool[n_units] — True if this unit is a spliced fragment.
    gdna_log_liks : np.ndarray
        float64[n_units] — pre-computed gDNA log-likelihood per unit.
        -inf for spliced units.
    genomic_footprints : np.ndarray
        int32[n_units] — genomic footprint per unit.
    frag_ids : np.ndarray
        int64[n_units] — buffer frag_id for each EM unit.
    frag_class : np.ndarray
        int8[n_units] — fragment class code per unit.
    splice_type : np.ndarray
        uint8[n_units] — SpliceType enum value per unit.
    n_units : int
        Number of ambiguous units.
    n_candidates : int
        Total number of (unit, candidate) entries.
    """

    offsets: np.ndarray
    t_indices: np.ndarray
    log_liks: np.ndarray
    count_cols: np.ndarray
    coverage_weights: np.ndarray
    tx_starts: np.ndarray
    tx_ends: np.ndarray
    locus_t_indices: np.ndarray
    locus_count_cols: np.ndarray
    is_spliced: np.ndarray
    gdna_log_liks: np.ndarray
    genomic_footprints: np.ndarray
    frag_ids: np.ndarray
    frag_class: np.ndarray
    splice_type: np.ndarray
    n_units: int
    n_candidates: int


# ======================================================================
# Locus — connected component of transcripts linked by shared fragments
# ======================================================================


@dataclass(slots=True)
class Locus:
    """A connected component of transcripts linked by shared fragments.

    Attributes
    ----------
    locus_id : int
        Sequential label (0-based).
    transcript_indices : np.ndarray
        int32 — global transcript indices in this locus.
    unit_indices : np.ndarray
        int32 — EM unit indices (rows in global CSR) belonging to
        this locus.
    gdna_span : int
        Total merged genomic footprint (bp) across all chromosomes.
    merged_intervals : list[tuple[str, int, int]]
        Merged (ref, start, end) intervals for this locus, cached
        during ``build_loci`` to avoid recomputation.
    """

    locus_id: int
    transcript_indices: np.ndarray
    unit_indices: np.ndarray
    gdna_span: int
    merged_intervals: list


# ======================================================================
# LocusEMInput — per-locus sub-problem with local component indices
# ======================================================================


@dataclass(slots=True)
class LocusEMInput:
    """Per-locus EM input: locally-renumbered components and likelihoods.

    Component layout::

        [0, n_t)     — transcript (one per transcript in the locus)
        [n_t]        — gdna (single gDNA component)

    Total components = n_transcripts + 1.

    Only UNSPLICED units have a gDNA candidate.  Spliced units compete
    only among transcript components.
    """

    locus: Locus
    offsets: np.ndarray
    t_indices: np.ndarray
    log_liks: np.ndarray
    count_cols: np.ndarray
    coverage_weights: np.ndarray
    tx_starts: np.ndarray
    tx_ends: np.ndarray
    locus_t_indices: np.ndarray
    locus_count_cols: np.ndarray
    n_transcripts: int
    n_components: int
    local_to_global_t: np.ndarray
    unambig_totals: np.ndarray
    gdna_init: float
    effective_lengths: np.ndarray
    prior: np.ndarray
    bias_profiles: np.ndarray | list | None


# ======================================================================
# LocusPartition — per-locus CSR subset for partitioned EM
# ======================================================================


@dataclass(slots=True)
class LocusPartition:
    """Per-locus CSR subset with contiguous, 0-indexed arrays.

    Produced by ``partition_and_free()`` which scatters the global
    ``ScoredFragments`` CSR into per-locus partitions.  Each partition
    is self-contained: its ``offsets`` array defines a local CSR over
    ``n_units`` rows and ``n_candidates`` total candidate entries.

    Transcript indices (``t_indices``) remain in **global** transcript
    space.  Remapping to local indices is deferred to the C++ extraction
    function, consistent with the existing ``extract_locus_sub_problem``.
    """

    locus_id: int
    n_units: int
    n_candidates: int

    # CSR structure
    offsets: np.ndarray  # int64[n_units + 1]

    # Per-candidate arrays (indexed by offsets)
    t_indices: np.ndarray  # int32 — GLOBAL transcript indices
    log_liks: np.ndarray  # float64
    count_cols: np.ndarray  # uint8
    coverage_weights: np.ndarray  # float64
    tx_starts: np.ndarray  # int32
    tx_ends: np.ndarray  # int32

    # Per-unit arrays
    is_spliced: np.ndarray  # uint8 (bool viewed as uint8 for C++)
    gdna_log_liks: np.ndarray  # float64
    genomic_footprints: np.ndarray  # int32
    locus_t_indices: np.ndarray  # int32
    locus_count_cols: np.ndarray  # uint8

    # Per-unit annotation metadata (not passed to C++ EM)
    frag_ids: np.ndarray  # int64[n_units] — buffer frag_id per unit
    frag_class: np.ndarray  # uint8[n_units] — fragment class code
    splice_type: np.ndarray  # uint8[n_units] — SpliceType enum value
