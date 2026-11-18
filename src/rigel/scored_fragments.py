"""rigel.scored_fragments — the EM's input containers, pure dataclasses.

``ScoredFragments`` is the global CSR linking fragment units to candidate transcripts, built by ``rigel.routing``;
``LocusPartition`` is one locus's slice of it, scattered by ``locus_partition.partition_and_free``.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class ScoredFragments:
    """The EM's candidate data, scored from the fragment buffer by ``FragmentRouter``: CSR arrays linking each
    ambiguous fragment unit to its candidate transcripts and their log-likelihoods. Transcript candidates only; the
    gDNA candidate is added per locus.

    Attributes
    ----------
    offsets : np.ndarray
        int64[n_units + 1] — CSR offsets into flat arrays.
    t_indices : np.ndarray
        int32[n_candidates] — candidate transcript indices.
    log_liks : np.ndarray
        float32[n_candidates] — per candidate, ``log P_strand + log P_insert`` plus the overhang penalty
        (log-penalty × overhang bases); a multimapper's candidate also carries its hit's sequence price (NM and
        soft-clipped bases) and its place weight (``native/scoring.cpp``).
    count_cols : np.ndarray
        uint8[n_candidates] — ``SpliceStrandCol`` index per candidate (0–5).
    coverage_weights : np.ndarray
        float32[n_candidates] — coverage weight per candidate.
    is_spliced : np.ndarray
        bool[n_units] — True when the unit has no gDNA reading, at its locus or intergenic: gDNA can explain
        none of its alignments (each carries a sequenced junction; the scorer's ``gdna_can_explain``), or
        its locus gDNA term fell outside the pruning bound of the unit's best RNA candidate
        (``gdna_competes``) and it has no intergenic place. An unspliced or implicit unit is True only in
        the second case.
    gdna_log_liks : np.ndarray
        float32[n_units] — the unit's gDNA log-likelihood, summed over its places in the locus; -inf where the
        locus's gDNA is no candidate.
    intergenic_log_liks : np.ndarray
        float32[n_units] — the sum of the gDNA terms at a multimapper's places with no transcript candidate,
        which no locus carries: the EM reads it at the intergenic rate, and it is never pruned. -inf where no
        intergenic place survived its gates.
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
    is_spliced: np.ndarray
    gdna_log_liks: np.ndarray
    frag_ids: np.ndarray
    frag_class: np.ndarray
    splice_type: np.ndarray
    n_units: int = 0
    n_candidates: int = 0
    intergenic_log_liks: np.ndarray | None = None


@dataclass(slots=True)
class LocusPartition:
    """One locus's slice of the global CSR, contiguous and 0-indexed: ``offsets`` is a local CSR over ``n_units`` rows
    and ``n_candidates`` entries. ``t_indices`` stay global; the C++ ``extract_locus_sub_problem_from_partition``
    remaps them."""

    locus_id: int
    n_units: int
    n_candidates: int

    # CSR structure
    offsets: np.ndarray  # int64[n_units + 1]

    # Per-candidate arrays (indexed by offsets)
    t_indices: np.ndarray  # int32 — global transcript indices
    log_liks: np.ndarray  # float32
    count_cols: np.ndarray  # uint8
    coverage_weights: np.ndarray  # float32

    # Per-unit arrays
    is_spliced: np.ndarray  # uint8 (bool viewed as uint8 for C++)
    gdna_log_liks: np.ndarray  # float32

    # Per-unit annotation metadata (not passed to C++ EM)
    frag_ids: np.ndarray  # int64[n_units] — buffer frag_id per unit
    frag_class: np.ndarray  # uint8[n_units] — fragment class code
    splice_type: np.ndarray  # uint8[n_units] — SpliceType enum value

    intergenic_log_liks: np.ndarray | None = None  # float32[n_units]
