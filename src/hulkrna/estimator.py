"""
hulkrna.estimator — Fragment abundance estimation with locus-level EM.

AbundanceEstimator manages per-transcript and per-gene fragment count
accumulation combined with Bayesian abundance estimation via EM.

For uniquely-mapping fragments, assignment is deterministic.  For ambiguous
fragments, transcript abundances are estimated via per-locus EM with a
coverage-weighted One Virtual Read (OVR) prior, then counts are accumulated
from the converged posterior.

Locus-level EM architecture
----------------------------
The transcriptome is partitioned into **loci** — connected components of
transcripts linked by shared fragments.  Each locus gets its own
independent EM instance solving:

    mRNA (per-transcript) + nRNA (per-transcript) + gDNA (one per locus)

The gDNA shadow is a **single genomic component** per locus, built from
the coverage islands of unspliced fragment alignments.  It represents the
genome as a competing source, NOT attributed to genes.  gDNA is a property
of the reference chromosomes; fragments classified as gDNA are reported
in genomic coordinates.

Posterior expected-count assignment
----------------------------------
After EM convergence, ambiguous RNA units are accumulated as expected
counts from the RNA-normalized posterior distribution rather than by
stochastic winner sampling.

EM algorithm
------------
MAP-EM with effective-length normalization and a coverage-weighted
One Virtual Read (OVR) Dirichlet prior.

1. Initialize theta from unique counts + OVR prior.
2. E-step: posterior proportional to likelihood * theta_t / eff_len_t.
3. M-step: theta proportional to unique + em_totals + prior.
4. Repeat until convergence (|delta| < ``_EM_CONVERGENCE_DELTA``).
5. Assign: add posterior-expected counts per ambiguous unit.
"""

import logging
import math as _math
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import digamma as _digamma

from .types import Strand
from .config import EMConfig, TranscriptGeometry
from .annotate import POOL_CODE_MRNA, POOL_CODE_NRNA, POOL_CODE_GDNA
from .categories import (
    ANTISENSE_COLS,
    SpliceType,
    SpliceStrandCol,
    NUM_SPLICE_STRAND_COLS,
    SPLICED_COLS,
)

logger = logging.getLogger(__name__)

# Epsilon added to theta before taking log to prevent log(0) in EM.
_EM_LOG_EPSILON = 1e-300

# Numerical-stability floor used inside locus.py to mark
# eligible prior components.  The actual prior values are
# computed in run_locus_em() from alpha + gamma * OVR.
EM_PRIOR_EPSILON = 1e-10

# Relative convergence tolerance for EM theta updates.
_EM_CONVERGENCE_DELTA = 1e-6

# Upper bound on fragment length for compound-key arithmetic in
# bias correction.  Fragments longer than this are clamped — any
# footprint >= transcript length already yields effective_length=1,
# so clamping is functionally neutral.
_MAX_FRAG_LEN = 1_000_000

# Default mean fragment length used when no observations are available
# (e.g. purely intergenic loci with no measured insert sizes).
_DEFAULT_MEAN_FRAG = 200.0

# Fraction of EM iteration budget allocated to SQUAREM acceleration.
_SQUAREM_BUDGET_DIVISOR = 3

# Number of VBEM/MAP-EM iterations after pruning zero-weight components.
_POST_PRUNE_ITERS = 10


def _apply_bias_correction_uniform(log_liks, t_indices, tx_starts, tx_ends,
                                    profile_lengths):
    """Uniform bias fast-path: pure vectorized NumPy, no Python loops.

    Under uniform bias:
    - effective_length(frag_len) = max(L - frag_len + 1, 1)
    - fragment_weight(start, end) = 1.0 (for valid spans)

    So the correction is simply ``-log(max(L - frag_len + 1, 1))``.

    Parameters
    ----------
    log_liks : np.ndarray
        float64[n_candidates] — modified in-place.
    t_indices : np.ndarray
        int32[n_candidates] — component indices into *profile_lengths*.
    tx_starts, tx_ends : np.ndarray
        int32[n_candidates] — transcript-relative coordinates.
    profile_lengths : np.ndarray
        int64[n_components] — transcript/component lengths.
    """
    frag_lens = np.clip(
        (tx_ends - tx_starts).astype(np.int64), 0, _MAX_FRAG_LEN - 1,
    )
    eff_lens = np.maximum(
        profile_lengths[t_indices] - frag_lens + 1, 1,
    ).astype(np.float64)
    log_liks -= np.log(eff_lens)


def _apply_bias_correction(log_liks, t_indices, tx_starts, tx_ends,
                           bias_profiles):
    """Apply per-candidate positional bias correction **in-place**.

    For each candidate *i* with component ``t_indices[i]`` and
    transcript coordinates ``[tx_starts[i], tx_ends[i])``, computes::

        log W(f|t) − log L̃_t(l_f)

    and adds it to ``log_liks[i]``.

    Under uniform bias this equals ``−log(max(L − l_f + 1, 1))``,
    exactly reproducing the former baked-in effective-length
    correction.

    Parameters
    ----------
    bias_profiles : np.ndarray or list
        If ``np.ndarray`` of int64 component lengths, dispatches
        directly to the vectorized uniform fast-path (no BiasProfile
        objects needed).  If a list of ``BiasProfile`` objects, checks
        for uniform profiles and falls back to per-candidate lookup.
    """
    n = len(t_indices)
    if n == 0:
        return

    # --- Fast path: numpy array of component lengths -----------------
    # This avoids creating BiasProfile objects entirely (Phase H).
    if isinstance(bias_profiles, np.ndarray):
        _apply_bias_correction_uniform(
            log_liks, t_indices, tx_starts, tx_ends, bias_profiles,
        )
        return

    if not bias_profiles:
        return

    # --- Fast path: all uniform BiasProfile objects -------------------
    if all(p.is_uniform for p in bias_profiles):
        profile_lengths = np.array(
            [p.length for p in bias_profiles], dtype=np.int64,
        )
        _apply_bias_correction_uniform(
            log_liks, t_indices, tx_starts, tx_ends, profile_lengths,
        )
        return

    # --- General path: LUT + per-candidate fragment weight -----------
    frag_lens = (tx_ends - tx_starts).astype(np.int64)
    # Clamp to [0, MAX_FLEN) so the compound key arithmetic is safe.
    # Discordant minimap2 pairs can produce genomic footprints > 1 Mb;
    # any footprint ≥ profile length already yields effective_length=1,
    # so clamping is functionally neutral.
    frag_lens = np.clip(frag_lens, 0, _MAX_FRAG_LEN - 1)
    compound = t_indices.astype(np.int64) * _MAX_FRAG_LEN + frag_lens
    unique_compound, inverse = np.unique(compound, return_inverse=True)

    log_eff = np.empty(len(unique_compound), dtype=np.float64)
    for j in range(len(unique_compound)):
        uc = int(unique_compound[j])
        cidx, fl = divmod(uc, _MAX_FRAG_LEN)
        log_eff[j] = _math.log(
            max(bias_profiles[int(cidx)].effective_length(int(fl)), _EM_LOG_EPSILON)
        )

    # Subtract log L̃ for every candidate (vectorised via inverse)
    log_liks -= log_eff[inverse]

    # --- Fragment weight ---------------------------------------------
    # Compute log W(f|t) for each candidate.  Under uniform bias
    # W = 1.0 → log W = 0, so the loop body is a no-op.  For
    # non-uniform profiles (Phase 2+) this contributes.
    for i in range(n):
        w = bias_profiles[int(t_indices[i])].fragment_weight(
            int(tx_starts[i]), int(tx_ends[i]),
        )
        if w != 1.0:
            log_liks[i] += _math.log(max(w, _EM_LOG_EPSILON))


def _build_equiv_classes(offsets, t_indices, log_liks, coverage_weights):
    """Group EM units by candidate component set into equivalence classes.

    Units sharing the exact same ordered set of candidate component
    indices are merged.  Each class stores component indices, a dense
    log-likelihood matrix, and a dense coverage-weight matrix.

    Parameters
    ----------
    offsets : np.ndarray
        int64[n_units + 1] CSR row pointers.
    t_indices : np.ndarray
        int32[n_candidates] candidate component indices.
    log_liks : np.ndarray
        float64[n_candidates] per-candidate log-likelihoods.
    coverage_weights : np.ndarray
        float64[n_candidates] per-candidate coverage weights.

    Returns
    -------
    list[tuple[np.ndarray, np.ndarray, np.ndarray]]
        Each element is ``(comp_indices, log_lik_matrix, wt_matrix)``.
    """
    n_units = len(offsets) - 1
    if n_units == 0 or len(t_indices) == 0:
        return []

    class_map: dict[tuple, list[int]] = {}
    for u in range(n_units):
        start = int(offsets[u])
        end = int(offsets[u + 1])
        if start == end:
            continue
        key = tuple(t_indices[start:end].tolist())
        class_map.setdefault(key, []).append(start)

    result: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for key, start_list in class_map.items():
        comp_idx = np.array(key, dtype=np.int32)
        k = len(comp_idx)
        n = len(start_list)
        ll_matrix = np.empty((n, k), dtype=np.float64)
        wt_matrix = np.empty((n, k), dtype=np.float64)
        for i, s in enumerate(start_list):
            ll_matrix[i, :] = log_liks[s:s + k]
            wt_matrix[i, :] = coverage_weights[s:s + k]
        result.append((comp_idx, ll_matrix, wt_matrix))

    return result


def _em_step(theta, ec_data, log_eff_len, unique_totals, prior, em_totals):
    """One EM iteration: theta -> theta_new, updating *em_totals* in place.

    Returns the normalised theta_new.  ``em_totals`` is overwritten with
    the raw posterior column sums from this E-step so that
    ``unique_totals + em_totals + prior`` gives the un-normalised alpha.
    """
    log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len
    em_totals[:] = 0.0

    for comp_idx, ll_matrix, scratch in ec_data:
        np.add(ll_matrix, log_weights[comp_idx], out=scratch)
        max_r = scratch.max(axis=1, keepdims=True)
        scratch -= max_r
        np.exp(scratch, out=scratch)
        row_sum = scratch.sum(axis=1, keepdims=True)

        bad = (row_sum == 0) | ~np.isfinite(row_sum)
        if bad.any():
            row_sum[bad] = 1.0
            scratch /= row_sum
            scratch[bad.ravel()] = 0.0
        else:
            scratch /= row_sum

        em_totals[comp_idx] += scratch.sum(axis=0)

    theta_new = unique_totals + em_totals + prior
    total = theta_new.sum()
    if total > 0:
        theta_new /= total
    return theta_new


def _vbem_step(alpha, ec_data, log_eff_len, unique_totals, prior, em_totals):
    """One VBEM iteration: alpha -> alpha_new, updating *em_totals* in place.

    Uses ``digamma(alpha_k) - digamma(sum(alpha))`` instead of
    ``log(theta)`` for the E-step weights.  The digamma function maps
    small alpha values to large negative numbers (e.g. digamma(0.01) ≈
    −100.6), which exponentially suppresses low-evidence components —
    much more aggressively than ``log(theta)`` in MAP-EM.

    Returns the updated alpha_new (un-normalised Dirichlet parameters).
    """
    alpha_sum = alpha.sum()
    log_weights = (
        _digamma(np.maximum(alpha, _EM_LOG_EPSILON))
        - _digamma(max(alpha_sum, _EM_LOG_EPSILON))
        - log_eff_len
    )
    em_totals[:] = 0.0

    for comp_idx, ll_matrix, scratch in ec_data:
        np.add(ll_matrix, log_weights[comp_idx], out=scratch)
        max_r = scratch.max(axis=1, keepdims=True)
        scratch -= max_r
        np.exp(scratch, out=scratch)
        row_sum = scratch.sum(axis=1, keepdims=True)

        bad = (row_sum == 0) | ~np.isfinite(row_sum)
        if bad.any():
            row_sum[bad] = 1.0
            scratch /= row_sum
            scratch[bad.ravel()] = 0.0
        else:
            scratch /= row_sum

        em_totals[comp_idx] += scratch.sum(axis=0)

    alpha_new = unique_totals + em_totals + prior
    return alpha_new


# ======================================================================
# ScanData - pre-computed CSR arrays for vectorized EM
# ======================================================================


@dataclass(slots=True)
class ScanData:
    """Fragment-level candidate data from the BAM scan pass.

    Contains pre-computed CSR (compressed sparse row) arrays linking
    each ambiguous fragment unit to its candidate transcripts and
    their log-likelihoods.  Produced by ``_scan_and_build_em_data()``
    and consumed during locus EM construction.

    The global ScanData contains mRNA + nRNA candidates only (NO gDNA).
    gDNA candidates are added per-locus during locus EM construction.

    Attributes
    ----------
    offsets : np.ndarray
        int64[n_units + 1] - CSR offsets into flat arrays.
    t_indices : np.ndarray
        int32[n_candidates] - candidate transcript or nRNA shadow indices.
    log_liks : np.ndarray
        float64[n_candidates] - log(P_strand x P_insert) per candidate.
    count_cols : np.ndarray
        uint8[n_candidates] - internal column index per candidate (0-7).
    coverage_weights : np.ndarray
        float64[n_candidates] - coverage weight per candidate.  mRNA
        candidates get position-based weight >= 1.0 (plateau = 1.0,
        edge > 1.0).  nRNA and gDNA candidates get weight 1.0.
    tx_starts : np.ndarray
        int32[n_candidates] - 0-based transcript-space start position
        of the fragment for each candidate.  For mRNA these are spliced
        transcript coordinates; for nRNA these are genomic-relative
        positions.  Used by the positional bias model.
    tx_ends : np.ndarray
        int32[n_candidates] - 0-based transcript-space end position
        (exclusive) for each candidate.  ``tx_ends - tx_starts`` gives
        the fragment footprint in transcript space.
    locus_t_indices : np.ndarray
        int32[n_units] - best transcript index per unit.
    locus_count_cols : np.ndarray
        uint8[n_units] - count column for the locus transcript.
    is_spliced : np.ndarray
        bool[n_units] - True if this unit is a spliced fragment.
        Unspliced units get a gDNA candidate in the locus EM.
    gdna_log_liks : np.ndarray
        float64[n_units] - pre-computed gDNA log-likelihood per unit.
        -inf for spliced units (gDNA not applicable).
    genomic_footprints : np.ndarray
        int32[n_units] - genomic footprint (aligned extent on genome)
        per unit.  Used for per-fragment effective length correction
        of gDNA candidates in ``build_locus_em_data``.
    frag_ids : np.ndarray
        int64[n_units] - buffer frag_id for each EM unit.
        Used to map EM assignment results back to BAM read-name groups
        for annotated BAM output.
    frag_class : np.ndarray
        int8[n_units] - fragment class code per unit (FRAG_* constants).
    splice_type : np.ndarray
        uint8[n_units] - SpliceType enum value per unit.
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
# Locus - connected component of transcripts linked by shared fragments
# ======================================================================


@dataclass(slots=True)
class Locus:
    """A connected component of transcripts linked by shared fragments.

    Attributes
    ----------
    locus_id : int
        Sequential label (0-based).
    transcript_indices : np.ndarray
        int32 - global transcript indices in this locus.
    gene_indices : np.ndarray
        int32 - global gene indices in this locus.
    unit_indices : np.ndarray
        int32 - EM unit indices (rows in global CSR) belonging to this locus.
    """

    locus_id: int
    transcript_indices: np.ndarray
    gene_indices: np.ndarray
    unit_indices: np.ndarray


# ======================================================================
# LocusEMInput - per-locus EM sub-problem with local component indices
# ======================================================================


@dataclass(slots=True)
class LocusEMInput:
    """Per-locus EM input: locally-renumbered components and likelihoods.

    Contains the data needed to run one independent EM instance for a
    single locus (connected component of transcripts).

    Component layout::

        [0, n_t)       - mRNA (one per transcript in the locus)
        [n_t, 2*n_t)   - nRNA (one per transcript, nascent shadow)
        [2*n_t]        - gDNA (ONE merged shadow for the entire locus)

    Total components = 2 * n_transcripts + 1.

    Only UNSPLICED units have a gDNA candidate.  Spliced units compete
    only among mRNA/nRNA components.

    Attributes
    ----------
    locus : Locus
        The parent locus.
    offsets : np.ndarray
        int64[n_local_units + 1] - CSR row pointers.
    t_indices : np.ndarray
        int32[n_local_candidates] - LOCAL component indices.
    log_liks : np.ndarray
        float64[n_local_candidates] - log-likelihoods.
    count_cols : np.ndarray
        uint8[n_local_candidates] - column indices for count accumulation.
    coverage_weights : np.ndarray
        float64[n_candidates] - coverage weight per candidate.  mRNA
        candidates get position-based weight >= 1.0 (plateau = 1.0,
        edge > 1.0).  nRNA and gDNA candidates get weight 1.0.
    tx_starts : np.ndarray
        int32[n_local_candidates] - transcript-space start positions.
    tx_ends : np.ndarray
        int32[n_local_candidates] - transcript-space end positions.
    locus_t_indices : np.ndarray
        int32[n_local_units] - best transcript index per unit (global).
    locus_count_cols : np.ndarray
        uint8[n_local_units] - count column for the best transcript.
    n_transcripts : int
        Number of transcripts in the locus.
    n_components : int
        2 * n_transcripts + 1 (mRNA + nRNA + 1 gDNA).
    local_to_global_t : np.ndarray
        int32[n_transcripts] - maps local transcript index to global.
    unique_totals : np.ndarray
        float64[n_components] - init counts (unique + shadow inits).
    nrna_init : np.ndarray
        float64[n_transcripts] - per-transcript nRNA init.
    gdna_init : float
        Scalar gDNA init from empirical Bayes estimation.
    effective_lengths : np.ndarray
        float64[n_components] - effective lengths for EM normalization.
    prior : np.ndarray
        float64[n_components] - Prior eligibility mask.  Eligible
        components carry ``EM_PRIOR_EPSILON``; ineligible ones are 0.0.
        The OVR prior (``gamma * coverage_weights``) is added to eligible
        components during ``run_locus_em``.
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
    unique_totals: np.ndarray
    nrna_init: np.ndarray
    gdna_init: float
    effective_lengths: np.ndarray
    prior: np.ndarray
    bias_profiles: np.ndarray | list | None


# ======================================================================
# AbundanceEstimator
# ======================================================================


class AbundanceEstimator:
    """Bayesian fragment abundance estimator with locus-level EM.

    Accumulates per-transcript and per-gene fragment counts from both
    deterministic (unique) and probabilistic (EM) assignment paths.
    Estimates transcript-level mRNA, nRNA, and gDNA abundances.

    Locus-level EM architecture::

        Per locus:
        [0, n_t)       - mRNA transcript components
        [n_t, 2*n_t)   - nRNA shadow per transcript
        [2*n_t]        - gDNA: ONE shadow per locus

    Fragment routing:

    - SPLICED_ANNOT + unique -> deterministic mRNA (no EM)
    - All other -> enter locus EM
      - Spliced fragments: mRNA + nRNA candidates only
      - Unspliced fragments: mRNA + nRNA + gDNA candidates

    Parameters
    ----------
    num_transcripts : int
    num_genes : int
    em_config : EMConfig or None
        EM algorithm configuration.  Defaults to ``EMConfig()``.
    geometry : TranscriptGeometry or None
        Pre-computed transcript/gene geometry.  None for minimal
        construction (e.g. in tests).
    """

    def __init__(
        self,
        num_transcripts: int,
        num_genes: int,
        *,
        em_config: EMConfig | None = None,
        geometry: TranscriptGeometry | None = None,
    ):
        if em_config is None:
            em_config = EMConfig()
        self.em_config = em_config
        self.num_transcripts = num_transcripts
        self.num_genes = num_genes
        self._rng = np.random.default_rng(em_config.seed)

        # nRNA shadow base index (for global CSR component numbering).
        self.nrna_base_index = num_transcripts

        # --- Geometry: from TranscriptGeometry or defaults ---
        if geometry is not None:
            self._t_to_g = np.asarray(geometry.t_to_g, dtype=np.int32)
            self._gene_spans = np.asarray(
                geometry.gene_spans, dtype=np.float64
            )
            self._mean_frag = geometry.mean_frag
            self._intronic_spans = np.asarray(
                geometry.intronic_spans, dtype=np.float64
            )
            self._transcript_spans = np.asarray(
                geometry.transcript_spans, dtype=np.float64
            )
            self._exonic_lengths = np.asarray(
                geometry.exonic_lengths, dtype=np.float64
            )
            self._t_eff_len = np.maximum(
                np.asarray(geometry.effective_lengths, dtype=np.float64),
                1.0,
            )
        else:
            self._t_to_g = None
            self._gene_spans = None
            self._mean_frag = _DEFAULT_MEAN_FRAG
            self._intronic_spans = None
            self._transcript_spans = None
            self._exonic_lengths = None
            self._t_eff_len = np.ones(
                num_transcripts, dtype=np.float64
            )

        self.unique_counts = np.zeros(
            (num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64
        )
        self.em_counts = np.zeros(
            (num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64
        )

        # Per-transcript high-confidence EM counts (posterior >= threshold).
        self.em_high_conf_counts = np.zeros(
            (num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64
        )

        # Per-transcript nRNA shadow initialization and EM counts.
        self.nrna_init = np.zeros(num_transcripts, dtype=np.float64)
        self.nrna_em_counts = np.zeros(num_transcripts, dtype=np.float64)

        # --- gDNA: locus-level, NOT per-gene ---
        # Total gDNA count as assigned by locus EM (sum across all loci)
        self._gdna_em_total = 0.0
        # Per-gene gDNA summary (for gene-level output compatibility).
        # Accumulated by mapping gDNA back to overlapping genes.
        self.gdna_gene_summary = np.zeros(num_genes, dtype=np.float64)

        # Per-locus gDNA results stored as list of dicts
        self.gdna_locus_results: list[dict] = []

        # Per-transcript gDNA locus attribution (for reporting).
        self.gdna_locus_counts = np.zeros(
            (num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64
        )

        # --- Pre-EM strand accumulators (for gDNA init via EB) ---
        # Gene-level: all single-gene UNSPLICED fragments.
        # Accumulated during scan pass for empirical Bayes gDNA prior.
        self.gene_sense_all = np.zeros(num_genes, dtype=np.float64)
        self.gene_antisense_all = np.zeros(num_genes, dtype=np.float64)

        # --- Pre-EM intronic accumulators (for nRNA init) ---
        self.transcript_intronic_sense = np.zeros(
            num_transcripts, dtype=np.float64
        )
        self.transcript_intronic_antisense = np.zeros(
            num_transcripts, dtype=np.float64
        )

        # Per-transcript confidence tracking:
        self._em_posterior_sum: np.ndarray | None = None
        self._em_n_assigned: np.ndarray | None = None

    @property
    def effective_lengths(self) -> np.ndarray:
        """Per-transcript effective lengths (read-only)."""
        return self._t_eff_len

    @property
    def gdna_em_count(self) -> float:
        """Total EM-assigned gDNA count across all loci."""
        return self._gdna_em_total

    @property
    def nrna_em_count(self) -> float:
        """Total EM-assigned nRNA count."""
        return float(self.nrna_em_counts.sum())

    @property
    def gdna_total(self) -> float:
        """Total gDNA count (EM-assigned)."""
        return self._gdna_em_total

    @property
    def gdna_contamination_rate(self) -> float:
        """Fraction of all fragments attributed to gDNA."""
        total_rna = float(self.unique_counts.sum() + self.em_counts.sum())
        total = total_rna + self.gdna_total
        if total == 0:
            return 0.0
        return self.gdna_total / total

    @property
    def t_counts(self) -> np.ndarray:
        """Total transcript counts (unique + EM), shape (N_t, 8)."""
        return self.unique_counts + self.em_counts

    @property
    def t_high_conf_counts(self) -> np.ndarray:
        """High-confidence transcript counts (unique + EM above threshold)."""
        return self.unique_counts + self.em_high_conf_counts

    # ------------------------------------------------------------------
    # Strand classification (internal)
    # ------------------------------------------------------------------

    @staticmethod
    def is_antisense(
        exon_strand,
        gene_strand: int,
        strand_model,
    ) -> bool:
        """Classify whether a fragment is antisense using the trained model."""
        if getattr(strand_model, '_finalized', False):
            p = strand_model.strand_likelihood_int(exon_strand, gene_strand)
        else:
            p = strand_model.strand_likelihood(
                exon_strand, Strand(gene_strand),
            )
        return p < 0.5

    # ------------------------------------------------------------------
    # Unique assignment
    # ------------------------------------------------------------------

    def assign_unique(self, resolved, index, strand_models) -> None:
        """Deterministic assignment for a truly unique fragment."""
        t_inds = resolved.t_inds

        if len(t_inds) == 0:
            return

        t_idx = int(next(iter(t_inds)))
        g_idx = int(index.t_to_g_arr[t_idx])
        gene_strand = int(index.g_to_strand_arr[g_idx])
        sm = strand_models.model_for_category(resolved.splice_type)
        anti = self.is_antisense(resolved.exon_strand, gene_strand, sm)
        col = SpliceStrandCol.from_category(resolved.splice_type, anti)
        self.unique_counts[t_idx, col] += 1.0

    # ------------------------------------------------------------------
    # Per-locus EM solver
    # ------------------------------------------------------------------

    def run_locus_em(
        self,
        locus_em: LocusEMInput,
        *,
        em_iterations: int = 1000,
        em_convergence_delta: float = _EM_CONVERGENCE_DELTA,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Run EM for a single locus sub-problem.

        Uses equivalence-class grouping: units sharing the same set of
        candidate components are merged into dense matrices, replacing
        the flat CSR ``reduceat``/``repeat``/``add.at`` inner loop with
        contiguous matrix operations per class.  The computation is
        mathematically identical (same posteriors, same M-step sums).

        Component layout per locus::

            [0, n_t)       - mRNA transcripts
            [n_t, 2*n_t)   - nRNA shadows
            [2*n_t]        - gDNA (single locus shadow)

        Only unspliced units have a gDNA candidate.  Spliced units
        compete only among mRNA/nRNA.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (theta, alpha) - converged parameters for this locus.
        """
        n_total = locus_em.n_components
        prior = locus_em.prior.copy()
        unique_totals = locus_em.unique_totals.copy()

        offsets = locus_em.offsets
        t_indices = locus_em.t_indices
        log_liks = locus_em.log_liks

        n_units = len(offsets) - 1

        # Apply per-candidate positional bias correction.  For uniform
        # profiles this is identical to the former baked-in
        # -log(max(t_len - frag_len + 1, 1)) correction.  Mutating
        # log_liks in-place is safe because each LocusEMInput is
        # consumed exactly once.
        _apply_bias_correction(
            log_liks, t_indices,
            locus_em.tx_starts, locus_em.tx_ends,
            locus_em.bias_profiles,
        )

        if n_units == 0 or len(t_indices) == 0:
            alpha = unique_totals + prior
            total = alpha.sum()
            theta = alpha / total if total > 0 else alpha
            return theta, alpha

        eff_len = locus_em.effective_lengths
        log_eff_len = np.log(eff_len)

        # Build equivalence classes: group units by candidate set.
        # Pre-allocate scratch arrays per class to avoid per-iteration
        # memory allocation (measured 1.78× speedup on chr17 mega-locus).
        ec_list = _build_equiv_classes(
            offsets, t_indices, log_liks, locus_em.coverage_weights,
        )
        ec_data: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        ec_wt: list[np.ndarray] = []
        for comp_idx, ll_matrix, wt_matrix in ec_list:
            scratch = np.empty_like(ll_matrix)
            ec_data.append((comp_idx, ll_matrix, scratch))
            ec_wt.append(wt_matrix)

        # Coverage-weighted warm start + One Virtual Read (OVR) prior.
        #
        # Warm start: distribute each ambiguous fragment's initial
        # share proportional to coverage weights (the "trapezoid"
        # geometric plausibility), restricted to prior-eligible
        # components.
        #
        # OVR prior: the same coverage distribution, scaled to sum
        # to ``em_prior_gamma`` total virtual reads, becomes the
        # Dirichlet prior.  This replaces the old flat pseudocount.
        # The prior for component k is:
        #
        #     prior_k = alpha + gamma * (coverage_totals_k / N_ambiguous)
        #
        # where alpha = em_prior_alpha (flat Dirichlet pseudocount)
        # and gamma = em_prior_gamma (OVR scale factor).
        #
        # Ineligible components (nRNA without intronic evidence,
        # gDNA without unspliced fragments) keep prior = 0.0 and
        # must not receive warm-start weight — the EM cannot drive
        # self-seeded phantom counts back to zero.
        eligible = prior > 0.0
        theta_init = unique_totals.copy()
        coverage_totals = np.zeros(n_total, dtype=np.float64)
        n_ambiguous = 0

        for i, (comp_idx, ll_matrix, _scratch) in enumerate(ec_data):
            wt_matrix = ec_wt[i]
            n = ll_matrix.shape[0]
            n_ambiguous += n
            has_p = eligible[comp_idx]
            # Mask weights by prior eligibility
            masked_wt = wt_matrix * has_p[np.newaxis, :]
            row_sums = masked_wt.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1.0
            shares = masked_wt / row_sums       # (n, k) normalized
            per_comp = shares.sum(axis=0)        # (k,)
            theta_init[comp_idx] += per_comp
            coverage_totals[comp_idx] += per_comp

        # Build prior: alpha (flat) + gamma-scaled OVR coverage weights.
        alpha = self.em_config.prior_alpha
        gamma = self.em_config.prior_gamma
        if n_ambiguous > 0 and gamma > 0.0:
            ovr = gamma * coverage_totals / n_ambiguous
        else:
            ovr = np.zeros(n_total, dtype=np.float64)
        # prior_k = alpha + OVR_k for eligible, 0 for ineligible.
        prior = np.where(eligible, alpha + ovr, 0.0)

        em_totals = np.zeros(n_total, dtype=np.float64)
        max_sq_iters = max(em_iterations // _SQUAREM_BUDGET_DIVISOR, 1)

        if self.em_config.mode == "vbem":
            # ---- VBEM with SQUAREM acceleration -----------------
            # Fixed-point iteration on Dirichlet parameters alpha.
            # E-step weights: digamma(alpha_k) - digamma(sum(alpha))
            # M-step: alpha_new = prior + unique_totals + em_totals
            alpha_dir = theta_init + prior  # initial Dirichlet params

            for _ in range(max_sq_iters):
                # --- two plain VBEM steps ---
                a1 = _vbem_step(
                    alpha_dir, ec_data, log_eff_len,
                    unique_totals, prior, em_totals,
                )
                a2 = _vbem_step(
                    a1, ec_data, log_eff_len,
                    unique_totals, prior, em_totals,
                )

                # --- SQUAREM extrapolation on alpha ---
                r = a1 - alpha_dir
                v = (a2 - a1) - r

                sv2 = np.dot(v, v)
                if sv2 == 0.0:
                    a_extrap = a2
                else:
                    srv = np.dot(r, v)
                    step = max(-srv / sv2, 1.0)
                    a_extrap = (
                        alpha_dir + 2.0 * step * r + step * step * v
                    )
                    # Alpha must stay positive.
                    np.maximum(a_extrap, _EM_LOG_EPSILON, out=a_extrap)

                # --- stabilisation VBEM step ---
                a_new = _vbem_step(
                    a_extrap, ec_data, log_eff_len,
                    unique_totals, prior, em_totals,
                )

                # Convergence on normalised theta.
                a_sum_old = alpha_dir.sum()
                a_sum_new = a_new.sum()
                delta = np.abs(
                    a_new / a_sum_new - alpha_dir / a_sum_old
                ).sum()
                alpha_dir = a_new

                if delta < em_convergence_delta:
                    break

            total = alpha_dir.sum()
            theta = alpha_dir / total if total > 0 else alpha_dir
            alpha_out = alpha_dir  # Dirichlet params ARE the alpha

        else:
            # ---- MAP-EM with SQUAREM acceleration ---------------
            # (Varadhan & Roland 2008)
            #
            # Each SQUAREM iteration performs 3 EM evaluations:
            #   1. theta_1 = M(theta_0)
            #   2. theta_2 = M(theta_1)
            #   3. theta_s = M(extrapolated) — stabilisation step
            # With SqS3 step length: alpha = -(r'v) / (v'v), >= 1.
            theta = theta_init + prior
            total = theta.sum()
            if total > 0:
                theta /= total

            for _ in range(max_sq_iters):
                # --- two plain EM steps ---
                theta_1 = _em_step(
                    theta, ec_data, log_eff_len,
                    unique_totals, prior, em_totals,
                )
                theta_2 = _em_step(
                    theta_1, ec_data, log_eff_len,
                    unique_totals, prior, em_totals,
                )

                # --- SQUAREM extrapolation ---
                r = theta_1 - theta
                v = (theta_2 - theta_1) - r

                sv2 = np.dot(v, v)
                if sv2 == 0.0:
                    theta_extrap = theta_2
                else:
                    srv = np.dot(r, v)
                    alpha_step = max(-srv / sv2, 1.0)

                    theta_extrap = (
                        theta
                        + 2.0 * alpha_step * r
                        + alpha_step * alpha_step * v
                    )
                    np.maximum(theta_extrap, 0.0, out=theta_extrap)
                    total = theta_extrap.sum()
                    if total > 0.0:
                        theta_extrap /= total
                    else:
                        theta_extrap = theta_2

                # --- stabilisation EM step ---
                theta_new = _em_step(
                    theta_extrap, ec_data, log_eff_len,
                    unique_totals, prior, em_totals,
                )

                delta = np.abs(theta_new - theta).sum()
                theta = theta_new

                if delta < em_convergence_delta:
                    break

            alpha_out = unique_totals + em_totals + prior

        # ----------------------------------------------------------
        # Post-EM pruning: zero components with no unique evidence
        # and low evidence ratio, then re-run EM to redistribute.
        # ----------------------------------------------------------
        if self.em_config.prune_threshold is not None:
            data_count = unique_totals + em_totals
            evidence_ratio = data_count / np.maximum(alpha_out, _EM_LOG_EPSILON)
            prune_mask = (
                (unique_totals == 0)
                & (evidence_ratio < self.em_config.prune_threshold)
            )
            # Never prune gDNA (last component).
            prune_mask[-1] = False

            if prune_mask.any():
                prior[prune_mask] = 0.0

                if self.em_config.mode == "vbem":
                    alpha_out[prune_mask] = _EM_LOG_EPSILON
                    for _ in range(_POST_PRUNE_ITERS):
                        alpha_out = _vbem_step(
                            alpha_out, ec_data, log_eff_len,
                            unique_totals, prior, em_totals,
                        )
                    total = alpha_out.sum()
                    theta = alpha_out / total if total > 0 else alpha_out
                else:
                    theta[prune_mask] = 0.0
                    total = theta.sum()
                    if total > 0:
                        theta /= total
                    for _ in range(_POST_PRUNE_ITERS):
                        theta = _em_step(
                            theta, ec_data, log_eff_len,
                            unique_totals, prior, em_totals,
                        )
                    alpha_out = unique_totals + em_totals + prior

        return theta, alpha_out

    # ------------------------------------------------------------------
    # Per-locus posterior assignment
    # ------------------------------------------------------------------

    def assign_locus_ambiguous(
        self,
        locus_em: LocusEMInput,
        theta: np.ndarray,
        *,
        confidence_threshold: float = 0.95,
        unit_annotations: list | None = None,
    ) -> float:
        """Assign ambiguous units within one locus, scatter to global arrays.

        Each unit is assigned fractionally across all components (mRNA,
        nRNA, gDNA) using the converged posterior probabilities.  This
        eliminates the binary gDNA threshold and prevents cliff-edge
        wipeout behaviour.

        Fully vectorized — no per-unit Python loop.

        Parameters
        ----------
        locus_em : LocusEMInput
        theta : np.ndarray
        confidence_threshold : float
        unit_annotations : list or None
            If provided, a list to which ``(global_unit_idx, best_tid,
            best_gid, pool_code, posterior)`` tuples are appended for
            each unit.  Used by the annotated BAM feature.

        Returns
        -------
        float
            Fractional count of fragments assigned to gDNA in this locus.
        """
        offsets = locus_em.offsets
        t_indices = locus_em.t_indices
        log_liks = locus_em.log_liks
        count_cols_arr = locus_em.count_cols
        n_units = len(offsets) - 1

        if n_units == 0 or len(t_indices) == 0:
            return 0.0

        n_t = locus_em.n_transcripts
        gdna_idx = 2 * n_t  # the single gDNA component index
        local_to_global_t = locus_em.local_to_global_t
        seg_lengths = np.diff(offsets).astype(np.intp)

        # Compute posteriors from converged theta
        eff_len = locus_em.effective_lengths
        log_weights = (
            np.log(theta + _EM_LOG_EPSILON) - np.log(eff_len)
        )
        log_posteriors = log_liks + log_weights[t_indices]

        offsets_int = offsets[:-1].astype(np.intp)
        seg_max = np.maximum.reduceat(log_posteriors, offsets_int)
        log_posteriors -= np.repeat(seg_max, seg_lengths)
        posteriors = np.exp(log_posteriors)
        seg_sum = np.add.reduceat(posteriors, offsets_int)
        # Guard: segments where every candidate is -inf → zero posteriors
        bad_seg = (seg_sum == 0) | ~np.isfinite(seg_sum)
        seg_sum[bad_seg] = 1.0
        posteriors /= np.repeat(seg_sum, seg_lengths)
        if bad_seg.any():
            bad_mask = np.repeat(bad_seg, seg_lengths)
            posteriors[bad_mask] = 0.0

        # --- Classify all candidates at once ---
        mrna_mask = t_indices < n_t
        nrna_mask = (t_indices >= n_t) & (t_indices < gdna_idx)
        gdna_mask = t_indices == gdna_idx

        # --- Vectorized mRNA assignment ---
        if mrna_mask.any():
            mrna_p = posteriors[mrna_mask]
            mrna_local_t = t_indices[mrna_mask]
            mrna_global_t = local_to_global_t[mrna_local_t]
            mrna_cols = count_cols_arr[mrna_mask]

            np.add.at(
                self.em_counts,
                (mrna_global_t, mrna_cols),
                mrna_p,
            )

            # High-confidence: per-unit max of mRNA posteriors
            # Set non-mRNA posteriors to 0, then reduceat for per-unit max
            mrna_p_for_max = np.where(mrna_mask, posteriors, 0.0)
            unit_max_mrna = np.maximum.reduceat(
                mrna_p_for_max, offsets_int,
            )
            conf_units = unit_max_mrna >= confidence_threshold
            if conf_units.any():
                conf_expanded = np.repeat(conf_units, seg_lengths)
                hc_mask = mrna_mask & conf_expanded
                if hc_mask.any():
                    hc_p = posteriors[hc_mask]
                    hc_local_t = t_indices[hc_mask]
                    hc_global_t = local_to_global_t[hc_local_t]
                    hc_cols = count_cols_arr[hc_mask]
                    np.add.at(
                        self.em_high_conf_counts,
                        (hc_global_t, hc_cols),
                        hc_p,
                    )

            # Confidence tracking
            if self._em_posterior_sum is None:
                self._em_posterior_sum = np.zeros(
                    self.num_transcripts, dtype=np.float64,
                )
                self._em_n_assigned = np.zeros(
                    self.num_transcripts, dtype=np.float64,
                )
            np.add.at(
                self._em_posterior_sum,
                mrna_global_t, mrna_p * mrna_p,
            )
            np.add.at(self._em_n_assigned, mrna_global_t, mrna_p)

        # --- Vectorized nRNA assignment ---
        if nrna_mask.any():
            nrna_p = posteriors[nrna_mask]
            nrna_local_t = t_indices[nrna_mask] - n_t
            nrna_global_t = local_to_global_t[nrna_local_t]
            np.add.at(self.nrna_em_counts, nrna_global_t, nrna_p)

        # --- Vectorized gDNA assignment ---
        gdna_count = 0.0
        if gdna_mask.any():
            gdna_p = posteriors[gdna_mask]
            gdna_count = float(gdna_p.sum())

            # Per-unit gDNA sum for locus attribution
            gdna_p_full = np.where(gdna_mask, posteriors, 0.0)
            gdna_per_unit = np.add.reduceat(
                gdna_p_full, offsets_int,
            )
            # Scatter to gdna_locus_counts using per-unit locus_t/locus_ct
            locus_t_arr = locus_em.locus_t_indices
            locus_ct_arr = locus_em.locus_count_cols
            valid_gdna_units = (gdna_per_unit > 0) & (locus_t_arr >= 0)
            if valid_gdna_units.any():
                np.add.at(
                    self.gdna_locus_counts,
                    (locus_t_arr[valid_gdna_units],
                     locus_ct_arr[valid_gdna_units]),
                    gdna_per_unit[valid_gdna_units],
                )

        # --- Per-unit annotation extraction (opt-in) ---
        if unit_annotations is not None:
            global_unit_indices = locus_em.locus.unit_indices
            t_to_g = self._t_to_g
            for u in range(n_units):
                s = int(offsets[u])
                e = int(offsets[u + 1])
                if e <= s:
                    continue
                seg_posteriors = posteriors[s:e]
                seg_t_indices = t_indices[s:e]
                best_local = int(np.argmax(seg_posteriors))
                best_posterior = float(seg_posteriors[best_local])
                best_comp = int(seg_t_indices[best_local])

                if best_comp < n_t:
                    # mRNA
                    g_tid = int(local_to_global_t[best_comp])
                    g_gid = int(t_to_g[g_tid])
                    pool_code = POOL_CODE_MRNA
                elif best_comp < gdna_idx:
                    # nRNA
                    g_tid = int(local_to_global_t[best_comp - n_t])
                    g_gid = int(t_to_g[g_tid])
                    pool_code = POOL_CODE_NRNA
                else:
                    # gDNA
                    g_tid = -1
                    g_gid = -1
                    pool_code = POOL_CODE_GDNA
                unit_annotations.append((
                    int(global_unit_indices[u]),
                    g_tid,
                    g_gid,
                    pool_code,
                    best_posterior,
                    e - s,  # n_candidates
                ))

        return gdna_count

    # ------------------------------------------------------------------
    # Confidence / posterior metrics
    # ------------------------------------------------------------------

    def posterior_mean(self) -> np.ndarray:
        """Per-transcript mean posterior at assignment."""
        n_assigned = self._em_n_assigned
        if n_assigned is None:
            return np.full(self.num_transcripts, np.nan)
        post_sum = self._em_posterior_sum
        with np.errstate(divide="ignore", invalid="ignore"):
            result = np.where(
                n_assigned > 0,
                post_sum / n_assigned,
                np.nan,
            )
        return result

    # ------------------------------------------------------------------
    # Output - primary counts
    # ------------------------------------------------------------------

    def get_counts_df(self, index) -> pd.DataFrame:
        """Primary transcript-level counts with identifiers."""
        total = self.t_counts
        unique = self.unique_counts

        count_total = total.sum(axis=1)
        count_unique = unique.sum(axis=1)
        count_spliced = total[:, list(SPLICED_COLS)].sum(axis=1)
        count_em = self.em_counts.sum(axis=1)
        count_high_conf = self.t_high_conf_counts.sum(axis=1)
        n_gdna = self.gdna_locus_counts.sum(axis=1)
        pmean = self.posterior_mean()

        df = pd.DataFrame({
            "transcript_id": index.t_df["t_id"].values,
            "gene_id": index.t_df["g_id"].values,
            "gene_name": index.t_df["g_name"].values,
            "effective_length": self._t_eff_len,
            "count": count_total,
            "count_unique": count_unique,
            "count_spliced": count_spliced,
            "count_em": count_em,
            "count_high_conf": count_high_conf,
            "n_gdna": n_gdna,
            "posterior_mean": pmean,
        })
        return df

    def get_gene_counts_df(self, index) -> pd.DataFrame:
        """Primary gene-level counts with gDNA summary."""
        t_to_g = index.t_to_g_arr
        total = self.t_counts
        unique = self.unique_counts

        n_genes = index.num_genes
        g_total = np.zeros((n_genes, NUM_SPLICE_STRAND_COLS), dtype=np.float64)
        np.add.at(g_total, t_to_g, total)

        g_unique = np.zeros((n_genes, NUM_SPLICE_STRAND_COLS), dtype=np.float64)
        np.add.at(g_unique, t_to_g, unique)

        count_total = g_total.sum(axis=1)
        count_unique = g_unique.sum(axis=1)
        count_spliced = g_total[:, list(SPLICED_COLS)].sum(axis=1)

        count_em_arr = np.zeros(n_genes, dtype=np.float64)
        np.add.at(count_em_arr, t_to_g, self.em_counts.sum(axis=1))

        count_hc_arr = np.zeros(n_genes, dtype=np.float64)
        np.add.at(count_hc_arr, t_to_g, self.t_high_conf_counts.sum(axis=1))

        # nRNA per gene
        n_nrna_per_t = self.nrna_init + self.nrna_em_counts
        n_nrna = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_nrna, t_to_g, n_nrna_per_t)

        # gDNA per gene (from genomic-coordinate gDNA mapped to genes)
        n_gdna = self.gdna_gene_summary.copy()

        # Antisense unique counts per gene
        t_antisense = unique[:, list(ANTISENSE_COLS)].sum(axis=1)
        n_antisense = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_antisense, t_to_g, t_antisense)

        # Pre-EM strand accumulators
        n_sense_all = self.gene_sense_all
        n_antisense_all = self.gene_antisense_all

        # Intronic accumulators per gene
        n_intronic_sense = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_intronic_sense, t_to_g, self.transcript_intronic_sense)
        n_intronic_antisense = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_intronic_antisense, t_to_g, self.transcript_intronic_antisense)

        # gDNA rate per gene
        denom = count_total + n_gdna
        with np.errstate(divide="ignore", invalid="ignore"):
            gdna_rate = np.where(denom > 0, n_gdna / denom, 0.0)

        # Gene effective length: abundance-weighted mean of transcript
        # effective lengths.  For genes with zero counts, use the
        # unweighted mean of transcript effective lengths (since there
        # are no fragments to weight, the effective lengths already
        # incorporate the fragment length model).
        t_eff = self._t_eff_len
        t_counts_flat = (self.unique_counts + self.em_counts).sum(axis=1)
        g_eff_num = np.zeros(n_genes, dtype=np.float64)
        g_eff_den = np.zeros(n_genes, dtype=np.float64)
        np.add.at(g_eff_num, t_to_g, t_counts_flat * t_eff)
        np.add.at(g_eff_den, t_to_g, t_counts_flat)
        # Fallback for zero-count genes: unweighted mean of transcript eff lens
        g_eff_sum = np.zeros(n_genes, dtype=np.float64)
        g_eff_cnt = np.zeros(n_genes, dtype=np.float64)
        np.add.at(g_eff_sum, t_to_g, t_eff)
        np.add.at(g_eff_cnt, t_to_g, 1.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            g_eff_mean = np.where(
                g_eff_cnt > 0, g_eff_sum / g_eff_cnt, 1.0,
            )
        with np.errstate(divide="ignore", invalid="ignore"):
            g_eff_len = np.where(
                g_eff_den > 0,
                g_eff_num / g_eff_den,
                g_eff_mean,
            )
        g_eff_len = np.maximum(g_eff_len, 1.0)

        df = pd.DataFrame({
            "gene_id": index.g_df["g_id"].values,
            "gene_name": index.g_df["g_name"].values,
            "effective_length": g_eff_len,
            "count": count_total,
            "count_unique": count_unique,
            "count_spliced": count_spliced,
            "count_em": count_em_arr,
            "count_high_conf": count_hc_arr,
            "n_nrna": n_nrna,
            "n_gdna": n_gdna,
            "n_antisense": n_antisense,
            "n_sense_all": n_sense_all,
            "n_antisense_all": n_antisense_all,
            "n_intronic_sense": n_intronic_sense,
            "n_intronic_antisense": n_intronic_antisense,
            "gdna_rate": gdna_rate,
        })
        return df

    # ------------------------------------------------------------------
    # Output - detail (long format QC breakdown)
    # ------------------------------------------------------------------

    def get_detail_df(self, index) -> pd.DataFrame:
        """Detailed counts in long format for QC."""
        t_ids = index.t_df["t_id"].values
        g_ids = index.t_df["g_id"].values

        frames = []
        for source_name, counts in (
            ("unique", self.unique_counts),
            ("em", self.em_counts),
        ):
            cat_counts = np.zeros(
                (self.num_transcripts, len(SpliceType)),
                dtype=np.float64,
            )
            for cat in SpliceType:
                sense_col = SpliceStrandCol.from_category(cat, False)
                anti_col = SpliceStrandCol.from_category(cat, True)
                cat_counts[:, int(cat)] = (
                    counts[:, sense_col] + counts[:, anti_col]
                )

            nz_t, nz_cat = np.nonzero(cat_counts)
            if len(nz_t) == 0:
                continue
            cat_order = [c.name.lower() for c in SpliceType]
            frames.append(
                pd.DataFrame({
                    "transcript_id": t_ids[nz_t],
                    "gene_id": g_ids[nz_t],
                    "category": pd.Categorical(
                        [SpliceType(c).name.lower() for c in nz_cat],
                        categories=cat_order,
                    ),
                    "source": source_name,
                    "count": cat_counts[nz_t, nz_cat],
                })
            )

        if not frames:
            return pd.DataFrame(
                columns=[
                    "transcript_id", "gene_id",
                    "category", "source", "count",
                ]
            )
        return pd.concat(frames, ignore_index=True)
