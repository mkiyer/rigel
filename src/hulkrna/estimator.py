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
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .types import Strand
from .config import EMConfig, TranscriptGeometry
from ._em_impl import run_locus_em_native as _run_locus_em_native
from .annotate import POOL_CODE_MRNA, POOL_CODE_NRNA, POOL_CODE_GDNA
from .splice import (
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


# NOTE: _apply_bias_correction_uniform, _apply_bias_correction,
# _build_equiv_classes, _em_step, _vbem_step have been moved to C++
# in hulkrna._em_impl (em_solver.cpp).  The entire SQUAREM loop now
# executes as a single native call via run_locus_em_native().


# ======================================================================
# ScoredFragments - pre-computed CSR arrays for vectorized EM
# ======================================================================


@dataclass(slots=True)
class ScoredFragments:
    """Fragment-level candidate data from the BAM scan pass.

    Contains pre-computed CSR (compressed sparse row) arrays linking
    each ambiguous fragment unit to its candidate transcripts and
    their log-likelihoods.  Produced by ``_scan_and_build_em_data()``
    and consumed during locus EM construction.

    The global ScoredFragments contains mRNA + nRNA candidates only (NO gDNA).
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
    unit_indices : np.ndarray
        int32 - EM unit indices (rows in global CSR) belonging to this locus.
    """

    locus_id: int
    transcript_indices: np.ndarray
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
    em_config : EMConfig or None
        EM algorithm configuration.  Defaults to ``EMConfig()``.
    geometry : TranscriptGeometry or None
        Pre-computed transcript/gene geometry.  None for minimal
        construction (e.g. in tests).
    """

    def __init__(
        self,
        num_transcripts: int,
        *,
        em_config: EMConfig | None = None,
        geometry: TranscriptGeometry | None = None,
    ):
        if em_config is None:
            em_config = EMConfig()
        self.em_config = em_config
        self.num_transcripts = num_transcripts
        self._rng = np.random.default_rng(em_config.seed)

        # nRNA shadow base index (for global CSR component numbering).
        self.nrna_base_index = num_transcripts

        # --- Geometry: from TranscriptGeometry or defaults ---
        if geometry is not None:
            self._t_to_g = np.asarray(geometry.t_to_g, dtype=np.int32)
            self._mean_frag = geometry.mean_frag
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
            self._mean_frag = _DEFAULT_MEAN_FRAG
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

        # Per-locus gDNA results stored as list of dicts
        self.gdna_locus_results: list[dict] = []

        # Per-transcript gDNA locus attribution (for reporting).
        self.gdna_locus_counts = np.zeros(
            (num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64
        )

        # --- Pre-EM strand accumulators (for gDNA init via EB) ---
        # Transcript-level: all same-strand UNSPLICED fragments.
        # Accumulated during scan pass for empirical Bayes gDNA prior.
        # Gene-level totals are derived via np.add.at(g_arr, t_to_g, t_arr).
        self.transcript_unspliced_sense = np.zeros(
            num_transcripts, dtype=np.float64
        )
        self.transcript_unspliced_antisense = np.zeros(
            num_transcripts, dtype=np.float64
        )

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
        ref_strand: int,
        strand_model,
    ) -> bool:
        """Classify whether a fragment is antisense using the trained model.

        Parameters
        ----------
        exon_strand : int
            Observed exon-block strand of the fragment.
        ref_strand : int
            Strand of the reference transcript (or gene — they are
            identical by definition).
        strand_model
            Trained strand model with ``strand_likelihood*`` methods.
        """
        if getattr(strand_model, '_finalized', False):
            p = strand_model.strand_likelihood_int(exon_strand, ref_strand)
        else:
            p = strand_model.strand_likelihood(
                exon_strand, Strand(ref_strand),
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
        t_strand = int(index.t_to_strand_arr[t_idx])
        sm = strand_models.model_for_category(resolved.splice_type)
        anti = self.is_antisense(resolved.exon_strand, t_strand, sm)
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
        # Bias correction expects an np.ndarray of int64 component
        # lengths (uniform fast path).  The C++ solver handles this
        # internally, so we just need to pass the raw arrays.
        #
        # The entire EM pipeline — bias correction, equivalence-class
        # building, coverage-weighted warm start, OVR prior, SQUAREM
        # loop, and post-prune redistribution — is executed in a single
        # C++ call via run_locus_em_native().
        theta, alpha_out, _em_totals = _run_locus_em_native(
            locus_em.offsets,
            locus_em.t_indices,
            locus_em.log_liks,
            locus_em.coverage_weights,
            locus_em.tx_starts,
            locus_em.tx_ends,
            locus_em.bias_profiles,
            locus_em.unique_totals,
            locus_em.effective_lengths,
            (locus_em.prior > 0).astype(np.float64),
            locus_em.n_components,
            self.em_config.prior_alpha,
            self.em_config.prior_gamma,
            em_iterations,
            em_convergence_delta,
            self.em_config.mode == "vbem",
            self.em_config.prune_threshold if self.em_config.prune_threshold is not None else -1.0,
            _POST_PRUNE_ITERS,
        )
        return np.asarray(theta), np.asarray(alpha_out)

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

        # gDNA per gene (derived from per-transcript gdna_locus_counts)
        n_gdna = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_gdna, t_to_g, self.gdna_locus_counts.sum(axis=1))

        # Antisense unique counts per gene
        t_antisense = unique[:, list(ANTISENSE_COLS)].sum(axis=1)
        n_antisense = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_antisense, t_to_g, t_antisense)

        # Pre-EM strand accumulators (aggregate from transcript-level)
        n_sense_all = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_sense_all, t_to_g, self.transcript_unspliced_sense)
        n_antisense_all = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_antisense_all, t_to_g, self.transcript_unspliced_antisense)

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
