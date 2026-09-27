"""rigel.scored_fragments — Data containers for the EM solver.

Pure dataclasses with no logic.

- ``ScoredFragments`` — global CSR arrays linking fragment units to
  candidate transcripts with log-likelihoods, built by ``scan.py``.
- ``LocusPartition`` — per-locus CSR subset for the partitioned native
  EM path, scattered from it by ``locus_partition.partition_and_free``.

The locus containers live in :mod:`rigel.locus`.
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
        float32[n_candidates] — log(P_strand × P_insert) per candidate.
    count_cols : np.ndarray
        uint8[n_candidates] — ``SpliceStrandCol`` index per candidate (0–9).
    coverage_weights : np.ndarray
        float32[n_candidates] — coverage weight per candidate.
    is_spliced : np.ndarray
        bool[n_units] — True if gDNA cannot explain this unit: certified RNA (the scorer's
        ``gdna_can_explain``). An artifact or an implicit splice within the maximum fragment length is
        False; a sequenced junction, or an implicit splice beyond that length, is True.
    gdna_log_liks : np.ndarray
        float32[n_units] — pre-computed gDNA log-likelihood per unit.
        -inf where ``is_spliced``.
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
    space; the C++ ``extract_locus_sub_problem_from_partition`` remaps them
    to local indices.
    """

    locus_id: int
    n_units: int
    n_candidates: int

    # CSR structure
    offsets: np.ndarray  # int64[n_units + 1]

    # Per-candidate arrays (indexed by offsets)
    t_indices: np.ndarray  # int32 — GLOBAL transcript indices
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
