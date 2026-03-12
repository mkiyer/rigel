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

    The global ScoredFragments contains mRNA + nRNA candidates only
    (NO gDNA).  gDNA candidates are added per-locus during locus EM
    construction.

    Attributes
    ----------
    offsets : np.ndarray
        int64[n_units + 1] — CSR offsets into flat arrays.
    t_indices : np.ndarray
        int32[n_candidates] — candidate transcript or nRNA shadow indices.
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
    nrna_base_index : int
        First nRNA shadow index (= num_transcripts).
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
    nrna_base_index: int


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
    """

    locus_id: int
    transcript_indices: np.ndarray
    unit_indices: np.ndarray


# ======================================================================
# LocusEMInput — per-locus sub-problem with local component indices
# ======================================================================


@dataclass(slots=True)
class LocusEMInput:
    """Per-locus EM input: locally-renumbered components and likelihoods.

    Component layout::

        [0, n_t)                — mRNA (one per transcript in the locus)
        [n_t, n_t + n_nrna)     — nRNA (one per unique nRNA span)
        [n_t + n_nrna]          — gDNA (ONE merged shadow for the locus)

    Total components = n_transcripts + n_nrna + 1.

    Only UNSPLICED units have a gDNA candidate.  Spliced units compete
    only among mRNA/nRNA components.
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
    n_nrna: int
    n_components: int
    local_to_global_t: np.ndarray
    local_to_global_nrna: np.ndarray
    local_t_to_local_nrna: np.ndarray
    nrna_to_t_offsets: np.ndarray
    nrna_to_t_indices: np.ndarray
    unambig_totals: np.ndarray
    nrna_init: np.ndarray
    gdna_init: float
    effective_lengths: np.ndarray
    prior: np.ndarray
    bias_profiles: np.ndarray | list | None
