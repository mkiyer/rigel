"""
hulkrna.counter - Read count accumulation with Variational Bayes EM.

ReadCounter accumulates discrete read counts into transcript arrays.
Each fragment contributes exactly one count (1.0) to a single transcript.

For uniquely-mapping fragments, assignment is deterministic.  For ambiguous
fragments, transcript abundances are estimated via Variational Bayes EM
(VBEM) then each fragment is assigned to exactly one transcript by
sampling from the posterior probability distribution.

Posterior sampling assignment
-----------------------------
After VBEM convergence, each ambiguous unit is assigned to a single
candidate by drawing from the posterior distribution.  This preserves
discrete integer counts (one fragment → one transcript) while correctly
recapitulating the posterior probabilities over many fragments.

Unlike MAP (argmax) assignment, sampling avoids winner-take-all
behavior when posteriors are similar (e.g. 100 candidates each at
~1%), and unlike expected-value assignment, it preserves the discrete
nature of the biological counting process.

VBEM algorithm
--------------
Uses Variational Bayes EM with a sparse Dirichlet prior to robustly
estimate transcript abundances.  The digamma function ψ replaces
log in the E-step, which prevents "superset collapse" — the failure
mode where a transcript containing all exons absorbs all shared-region
reads.

1. Initialize α from unique counts + Dirichlet prior α₀.
2. E-step: posterior ∝ likelihood × exp(ψ(α) − log(ℓ_eff)).
3. M-step: α_t = α₀ + unique_t + Σ_f P(t | f).
4. Repeat until convergence.
5. Assign: sample one candidate per unit from posterior → 1 count.

The sparse prior (α₀ ≪ 1) pushes low-evidence components toward zero
via the digamma function's behavior: ψ(x) → −∞ much faster than
log(x) as x → 0, creating a natural sparsity-inducing regularizer.
This is equivalent to the approach used by salmon's VBEM.
"""

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import digamma

from .types import Strand
from .categories import (
    ANTISENSE_COLS,
    CountCategory,
    CountCol,
    NUM_COUNT_COLS,
    SPLICED_COLS,
)

logger = logging.getLogger(__name__)

# Baseline Dirichlet prior for gDNA shadow components.
# Shadows should not be sparsified (they model gDNA contamination,
# not isoform abundance), so they use a non-informative prior.
_SHADOW_PRIOR = 1.0

# Number of standard EM warm-up iterations before switching to VBEM.
# Standard EM uses log(theta) instead of digamma(alpha), which avoids
# the degeneracy where low-alpha components are irrecoverably pushed
# to zero.  The warm-up provides VBEM with a good starting point.
_VBEM_WARMUP_ITERS = 5

# Number of VBEM re-estimation passes.  After the initial VBEM
# converges with fixed models, we re-estimate strand models from the
# converged posteriors and re-converge.  Repeating this allows the
# abundance estimates and model parameters to co-adapt.
_REEST_PASSES = 2

# Epsilon added to theta before taking log to prevent log(0) in EM.
_EM_LOG_EPSILON = 1e-300

# Relative convergence tolerance for VBEM alpha updates.
_VBEM_CONVERGENCE_DELTA = 1e-6

# Bayesian priors used during EM model re-estimation.
_STRAND_PRIOR_N = 10.0   # pseudo-observations for strand
_INSERT_PRIOR_N = 50.0   # pseudo-observations for insert

# Numeric floor used for strand log-probability clamping.
_REEST_LOG_FLOOR = 1e-10


# ======================================================================
# EMData - pre-computed CSR arrays for vectorized EM
# ======================================================================


@dataclass(slots=True)
class EMData:
    """Pre-computed CSR arrays for vectorized EM.

    Attributes
    ----------
    offsets : np.ndarray
        int64[n_units + 1] - CSR offsets into flat arrays.
    t_indices : np.ndarray
        int32[n_candidates] - candidate transcript or shadow indices.
    log_liks : np.ndarray
        float64[n_candidates] - log(P_strand x P_insert) per candidate.
    count_cols : np.ndarray
        uint8[n_candidates] - internal column index per candidate (0-7).
    locus_t_indices : np.ndarray
        int32[n_units] - best transcript index per unit.
    locus_count_cols : np.ndarray
        uint8[n_units] - count column for the locus transcript.
    n_units : int
        Number of ambiguous units.
    n_candidates : int
        Total number of (unit, candidate) entries.
    gdna_base_index : int
        First shadow index.
    cand_strand_dir : np.ndarray or None
        int8[n_candidates] — strand direction per candidate.
        +1 = sense (exon matches gene/POS), −1 = antisense, 0 = ambiguous.
        Used for EM model re-estimation.
    cand_is_shadow : np.ndarray or None
        bool[n_candidates] — True for gDNA shadow candidates.
    cand_insert_sizes : np.ndarray or None
        int32[n_candidates] — insert size per candidate (expanded from
        per-unit; same for all candidates within a unit).
    cand_splice_ll : np.ndarray or None
        float64[n_candidates] — splice penalty log-likelihood per
        candidate (non-zero only for gDNA shadows; fixed during EM).
    init_p_rna_sense : float
        Initial RNA strand specificity (from exonic_spliced model).
    init_p_gdna_sense : float
        Initial gDNA strand balance (from intergenic model).
    init_rna_insert_ll : np.ndarray or None
        float64[max_insert+1] — initial RNA insert log-probabilities.
    init_gdna_insert_ll : np.ndarray or None
        float64[max_insert+1] — initial gDNA insert log-probabilities.
    insert_size_max : int
        Maximum insert size for histogram dimensioning.
    """

    offsets: np.ndarray
    t_indices: np.ndarray
    log_liks: np.ndarray
    count_cols: np.ndarray
    locus_t_indices: np.ndarray
    locus_count_cols: np.ndarray
    n_units: int
    n_candidates: int
    gdna_base_index: int
    cand_strand_dir: np.ndarray | None = None
    cand_is_shadow: np.ndarray | None = None
    cand_insert_sizes: np.ndarray | None = None
    cand_splice_ll: np.ndarray | None = None
    init_p_rna_sense: float = 0.5
    init_p_gdna_sense: float = 0.5
    init_rna_insert_ll: np.ndarray | None = None
    init_gdna_insert_ll: np.ndarray | None = None
    insert_size_max: int = 1000


# ======================================================================
# ReadCounter
# ======================================================================


class ReadCounter:
    """Accumulates read counts into transcript arrays.

    Count arrays have shape (N, 8) where columns correspond to
    4 categories × 2 strands (sense/antisense).

    Uses Variational Bayes EM (VBEM) with a sparse Dirichlet prior
    for robust abundance estimation, matching the approach used by
    salmon.  The digamma function in the E-step prevents superset
    collapse and induces sparsity in the transcript abundance vector.

    Parameters
    ----------
    num_transcripts : int
    num_genes : int
    seed : int or None
        Random seed for reproducible posterior sampling assignment.
    alpha : float
        Dirichlet prior hyperparameter for VBEM (default 0.01).
        Small values (≪ 1) induce sparsity; larger values (≥ 1)
        approach the standard EM behavior.
    shadow_init : np.ndarray or None
        Per-gene shadow initialization, shape (num_genes,).
    """

    def __init__(
        self,
        num_transcripts: int,
        num_genes: int,
        *,
        seed=None,
        alpha: float = 0.01,
        shadow_init: np.ndarray | None = None,
        effective_lengths: np.ndarray | None = None,
        t_to_g: np.ndarray | None = None,
    ):
        self.num_transcripts = num_transcripts
        self.num_genes = num_genes
        self.vbem_prior = alpha
        self._rng = np.random.default_rng(seed)

        # Base index for per-gene shadow components in EM theta vector.
        self.gdna_base_index = num_transcripts

        # Transcript → gene index mapping (for per-gene shadow eff_len)
        self._t_to_g = (
            np.asarray(t_to_g, dtype=np.int32) if t_to_g is not None
            else None
        )

        # --- Effective lengths for EM normalization ---
        # Per-transcript effective length (exonic length − mean_frag + 1,
        # floored to 1).  Used in the EM E-step to convert raw counts
        # (theta) to abundance-per-nucleotide (rho) so that longer
        # transcripts do not automatically attract more shared-exon
        # reads.  When None, all transcripts are assumed to have the
        # same effective length (no normalization).
        if effective_lengths is not None:
            self._t_eff_len = np.maximum(
                np.asarray(effective_lengths, dtype=np.float64), 1.0
            )
        else:
            self._t_eff_len = np.ones(
                num_transcripts, dtype=np.float64
            )

        self.unique_counts = np.zeros(
            (num_transcripts, NUM_COUNT_COLS), dtype=np.float64
        )
        self.em_counts = np.zeros(
            (num_transcripts, NUM_COUNT_COLS), dtype=np.float64
        )

        # Per-transcript high-confidence EM counts (posterior >= threshold).
        self.em_high_conf_counts = np.zeros(
            (num_transcripts, NUM_COUNT_COLS), dtype=np.float64
        )

        # Per-gene gDNA shadow tracking
        if shadow_init is not None:
            self.shadow_init = np.asarray(shadow_init, dtype=np.float64)
        else:
            self.shadow_init = np.zeros(num_genes, dtype=np.float64)

        # Per-gene EM-assigned gDNA counts.
        self.gdna_em_counts = np.zeros(num_genes, dtype=np.float64)

        # Per-transcript locus attribution of gDNA assignments.
        self.gdna_locus_counts = np.zeros(
            (num_transcripts, NUM_COUNT_COLS), dtype=np.float64
        )

        # Converged transcript abundance from EM (set by run_em)
        self._converged_theta: np.ndarray | None = None
        self._converged_alpha: np.ndarray | None = None

        # Per-transcript confidence tracking:
        # _em_posterior_sum[t] = sum of posteriors at which t was chosen
        # _em_n_assigned[t] = number of units assigned to t
        self._em_posterior_sum: np.ndarray | None = None
        self._em_n_assigned: np.ndarray | None = None

    @property
    def effective_lengths(self) -> np.ndarray:
        """Per-transcript effective lengths (read-only)."""
        return self._t_eff_len

    def _build_eff_len_vector(self) -> np.ndarray:
        """Build full effective-length vector for EM (transcripts + shadows).

        For transcript components [0 .. N_t), uses _t_eff_len.
        For shadow components [N_t .. N_t + N_g), uses the maximum
        effective length among transcripts within each gene, so that
        the shadow competes fairly with the gene's transcripts.
        """
        n_total = self.num_transcripts + self.num_genes
        eff_len = np.ones(n_total, dtype=np.float64)
        eff_len[:self.num_transcripts] = self._t_eff_len
        # Per-gene shadow eff_len: max over transcripts in that gene.
        # This is conservative: the shadow "covers" at least as much
        # effective length as the longest transcript.
        if hasattr(self, '_t_to_g') and self._t_to_g is not None:
            for g in range(self.num_genes):
                mask = self._t_to_g == g
                if mask.any():
                    eff_len[self.gdna_base_index + g] = self._t_eff_len[mask].max()
        else:
            # Fallback: use mean transcript eff_len for all shadows
            if self.num_transcripts > 0:
                mean_eff = self._t_eff_len.mean()
                eff_len[self.gdna_base_index:] = mean_eff
        return eff_len

    @property
    def gdna_unique_count(self) -> float:
        """Total shadow_init (empirical Bayes gDNA prior)."""
        return float(self.shadow_init.sum())

    @property
    def gdna_em_count(self) -> float:
        """Total EM-assigned gDNA count."""
        return float(self.gdna_em_counts.sum())

    @property
    def gdna_total(self) -> float:
        """Total gDNA count (shadow_init + EM-assigned)."""
        return self.gdna_unique_count + self.gdna_em_count

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
        """High-confidence transcript counts (unique + EM above threshold).

        Unique counts are always high-confidence.  EM counts are
        included only when the RNA-normalized posterior of the winner
        met or exceeded ``confidence_threshold`` during assignment.
        """
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
        """Classify whether a fragment is antisense using the trained model.

        Uses the strand model decision boundary (p < 0.5).

        Returns
        -------
        bool
            True if antisense (p < 0.5), False if sense (p >= 0.5).
        """
        p = strand_model.strand_likelihood(exon_strand, Strand(gene_strand))
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
        sm = strand_models.model_for_category(resolved.count_cat)
        anti = self.is_antisense(resolved.exon_strand, gene_strand, sm)
        col = CountCol.from_category(resolved.count_cat, anti)
        self.unique_counts[t_idx, col] += 1.0

    # ------------------------------------------------------------------
    # EM iteration (vectorized)
    # ------------------------------------------------------------------

    def run_em(self, em_data: EMData, *, em_iterations: int = 10) -> None:
        """Run Variational Bayes EM with two-pass model re-estimation.

        Uses **VBEM** with a sparse Dirichlet prior, effective length
        normalization, and post-convergence re-estimation of strand
        models.

        The algorithm proceeds in phases:

        1. **Warm-up** — standard EM (5 iterations) to seed α.
        2. **VBEM pass 0** — converge with initial (empirical) models.
        3. **Re-estimate** strand models from converged posteriors.
        4. **VBEM pass 1** — re-converge with updated log-likelihoods.

        This two-pass strategy avoids the instability of mid-loop
        re-estimation: posteriors must be converged before they can
        reliably inform model parameters.

        Initial strand and insert models come from empirical per-category
        observations (trained during the BAM scan).  After the first
        VBEM convergence, strand models are re-estimated from the
        converged posterior-weighted fragment observations, improving
        RNA/gDNA discrimination.

        The digamma function ψ replaces log(θ) in the E-step:

            log q(t|f) = log L(f|t) + ψ(α_t) − log(ℓ_eff_t)

        This prevents "superset collapse" and induces sparsity.
        """
        n_total = self.num_transcripts + self.num_genes
        gdna_base = self.gdna_base_index

        # Per-component prior: sparse for transcripts, baseline for
        # gDNA shadows.
        prior = np.full(n_total, self.vbem_prior, dtype=np.float64)
        prior[gdna_base:] = _SHADOW_PRIOR

        # Initialize alpha from unique counts
        unique_totals = np.zeros(n_total, dtype=np.float64)
        unique_totals[:self.num_transcripts] = self.unique_counts.sum(axis=1)
        unique_totals[gdna_base:gdna_base + self.num_genes] = self.shadow_init

        if em_data.n_units == 0 or em_data.n_candidates == 0:
            alpha = unique_totals + prior
            self._converged_alpha = alpha
            self._converged_theta = alpha / alpha.sum()
            return

        offsets = em_data.offsets
        t_indices = em_data.t_indices
        log_liks = em_data.log_liks  # Modified in-place during re-estimation
        seg_lengths = np.diff(offsets)

        # Pre-compute log effective lengths for the E-step
        eff_len = self._build_eff_len_vector()
        log_eff_len = np.log(eff_len)

        # ── Model re-estimation setup ───────────────────────────────
        # Check if per-candidate metadata is available for model
        # re-estimation.  If not, fall back to fixed-model EM.
        can_reestimate = (
            em_data.cand_strand_dir is not None
            and em_data.cand_is_shadow is not None
            and em_data.cand_insert_sizes is not None
            and em_data.cand_splice_ll is not None
            and em_data.init_rna_insert_ll is not None
            and em_data.init_gdna_insert_ll is not None
        )

        if can_reestimate:
            cand_strand_dir = em_data.cand_strand_dir
            cand_is_shadow = em_data.cand_is_shadow
            cand_insert_sizes = em_data.cand_insert_sizes
            cand_splice_ll = em_data.cand_splice_ll
            max_isize = em_data.insert_size_max

            rna_mask = ~cand_is_shadow
            gdna_mask = cand_is_shadow

            # Pre-compute strand direction masks (invariant across iterations)
            sense_mask = cand_strand_dir == 1
            anti_mask = cand_strand_dir == -1
            ambig_mask = cand_strand_dir == 0

            # Insert size validity mask (insert_size > 0)
            valid_insert = cand_insert_sizes > 0
            # Clamp insert sizes into histogram range
            clamped_isizes = np.clip(cand_insert_sizes, 0, max_isize)

            # Current model parameters (evolve during EM)
            p_rna_sense = em_data.init_p_rna_sense
            p_gdna_sense = em_data.init_p_gdna_sense
            rna_insert_ll = em_data.init_rna_insert_ll.copy()
            gdna_insert_ll = em_data.init_gdna_insert_ll.copy()

            # Insert histogram priors (scaled copies of initial PMFs)
            _init_rna_pmf = np.exp(rna_insert_ll)
            _init_rna_pmf /= _init_rna_pmf.sum()
            _rna_hist_prior = _init_rna_pmf * _INSERT_PRIOR_N

            _init_gdna_pmf = np.exp(gdna_insert_ll)
            _init_gdna_pmf /= _init_gdna_pmf.sum()
            _gdna_hist_prior = _init_gdna_pmf * _INSERT_PRIOR_N

            # Pre-compute the non-strand component of log_liks.
            # We only re-estimate the strand component during EM;
            # re-estimating insert separately for RNA/gDNA creates
            # a self-reinforcing asymmetry where the dominant
            # component gets a denser histogram → better insert LL
            # → becomes more dominant.  Instead, we keep the
            # original insert + splice LL fixed (insert cancels
            # between RNA and gDNA as in the original design).
            _init_strand_ll = np.full(len(log_liks), np.log(0.5))
            _init_strand_ll[rna_mask & sense_mask] = np.log(
                max(p_rna_sense, _REEST_LOG_FLOOR)
            )
            _init_strand_ll[rna_mask & anti_mask] = np.log(
                max(1.0 - p_rna_sense, _REEST_LOG_FLOOR)
            )
            _init_strand_ll[gdna_mask & sense_mask] = np.log(
                max(p_gdna_sense, _REEST_LOG_FLOOR)
            )
            _init_strand_ll[gdna_mask & anti_mask] = np.log(
                max(1.0 - p_gdna_sense, _REEST_LOG_FLOOR)
            )
            # non-strand = initial_log_liks - initial_strand
            _non_strand_ll = log_liks - _init_strand_ll

        def _reestimate_and_rescore(posteriors):
            """Re-estimate strand/insert models and recompute log_liks.

            Uses posterior weights to decompose fragments into RNA and
            gDNA components.  **Strand models** are re-estimated and
            used to update log_liks.  **Insert models** are
            re-estimated for monitoring but NOT used for log_lik
            updates — using separate insert models for RNA and gDNA
            creates a self-reinforcing asymmetry where the dominant
            component gets better insert LL, so we keep insert
            cancellation from the original design.

            Bayesian regularization prevents instability: initial
            model parameters serve as priors (pseudo-observations).
            """
            nonlocal p_rna_sense, p_gdna_sense, rna_insert_ll, gdna_insert_ll

            # ── Strand model re-estimation (Beta prior) ─────────────
            # Only RNA strand specificity is re-estimated.  gDNA is
            # definitionally unstranded (~0.5) — re-estimating it
            # would contaminate the model with RNA signal from
            # fragments that have high gDNA posteriors in early
            # iterations, collapsing RNA/gDNA discrimination.
            rna_sense_wt = (
                posteriors[rna_mask & sense_mask].sum()
                + em_data.init_p_rna_sense * _STRAND_PRIOR_N
            )
            rna_anti_wt = (
                posteriors[rna_mask & anti_mask].sum()
                + (1.0 - em_data.init_p_rna_sense) * _STRAND_PRIOR_N
            )
            p_rna_sense = rna_sense_wt / (rna_sense_wt + rna_anti_wt)
            # p_gdna_sense stays fixed at init value

            # ── Insert model re-estimation (monitoring only) ────────
            # Estimate separate RNA/gDNA insert histograms but do NOT
            # use them for scoring (insert cancels by design).
            rna_hist = _rna_hist_prior.copy()
            gdna_hist = _gdna_hist_prior.copy()

            rna_valid = rna_mask & valid_insert
            gdna_valid = gdna_mask & valid_insert

            if rna_valid.any():
                np.add.at(
                    rna_hist, clamped_isizes[rna_valid],
                    posteriors[rna_valid]
                )
            if gdna_valid.any():
                np.add.at(
                    gdna_hist, clamped_isizes[gdna_valid],
                    posteriors[gdna_valid]
                )

            rna_total = rna_hist.sum()
            rna_insert_ll[:] = (
                np.log(rna_hist + 1.0)
                - np.log(rna_total + max_isize + 1)
            )
            gdna_total = gdna_hist.sum()
            gdna_insert_ll[:] = (
                np.log(gdna_hist + 1.0)
                - np.log(gdna_total + max_isize + 1)
            )

            # ── Recompute log_liks (strand only) ────────────────────
            # Only the strand component is updated; insert + splice
            # are preserved from the initial scoring.
            log_p_sense_rna = np.log(max(p_rna_sense, _REEST_LOG_FLOOR))
            log_p_anti_rna = np.log(max(1.0 - p_rna_sense, _REEST_LOG_FLOOR))
            log_p_sense_gdna = np.log(max(p_gdna_sense, _REEST_LOG_FLOOR))
            log_p_anti_gdna = np.log(max(1.0 - p_gdna_sense, _REEST_LOG_FLOOR))
            log_half = np.log(0.5)

            strand_ll = np.full(len(log_liks), log_half, dtype=np.float64)
            strand_ll[rna_mask & sense_mask] = log_p_sense_rna
            strand_ll[rna_mask & anti_mask] = log_p_anti_rna
            strand_ll[gdna_mask & sense_mask] = log_p_sense_gdna
            strand_ll[gdna_mask & anti_mask] = log_p_anti_gdna

            log_liks[:] = strand_ll + _non_strand_ll

        # ── Phase 1: Standard EM warm-up ────────────────────────────
        theta = (unique_totals + prior)
        theta /= theta.sum()
        n_warmup = min(_VBEM_WARMUP_ITERS, em_iterations)

        em_totals = np.zeros(n_total, dtype=np.float64)
        for it in range(n_warmup):
            log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len
            log_posteriors = log_liks + log_weights[t_indices]

            seg_max = np.maximum.reduceat(log_posteriors, offsets[:-1])
            log_posteriors -= np.repeat(seg_max, seg_lengths)
            posteriors = np.exp(log_posteriors)
            seg_sum = np.add.reduceat(posteriors, offsets[:-1])
            posteriors /= np.repeat(seg_sum, seg_lengths)

            em_totals = np.zeros(n_total, dtype=np.float64)
            np.add.at(em_totals, t_indices, posteriors)

            theta = unique_totals + em_totals + prior
            theta /= theta.sum()

        logger.debug(
            f"  VBEM warm-up: {n_warmup} standard EM iterations"
        )

        # ── Phase 2: VBEM with digamma ─────────────────────────────
        # Two-pass strategy for model re-estimation:
        #   Pass 0: VBEM with initial (empirical) models → converge
        #   Pass 1+: re-estimate strand from converged posteriors,
        #            then re-converge VBEM with updated log_liks
        # This avoids the instability of mid-loop re-estimation
        # (posteriors must be converged before they can reliably
        # inform model parameters).
        alpha = unique_totals + em_totals + prior
        n_vbem = em_iterations - n_warmup
        n_passes = _REEST_PASSES if can_reestimate else 1
        total_vbem_iters = 0

        for reest_pass in range(n_passes):
            for it in range(n_vbem):
                # E-step: posterior ∝ L(f|t) × exp(ψ(α_t) − log(ℓ_eff_t))
                log_weights = digamma(alpha) - log_eff_len
                log_posteriors = log_liks + log_weights[t_indices]

                # Segment-wise log-sum-exp normalization
                seg_max = np.maximum.reduceat(log_posteriors, offsets[:-1])
                log_posteriors -= np.repeat(seg_max, seg_lengths)
                posteriors = np.exp(log_posteriors)
                seg_sum = np.add.reduceat(posteriors, offsets[:-1])
                posteriors /= np.repeat(seg_sum, seg_lengths)

                # M-step: accumulate fractional counts per component
                em_totals = np.zeros(n_total, dtype=np.float64)
                np.add.at(em_totals, t_indices, posteriors)

                # Update Dirichlet hyperparameters
                alpha_new = unique_totals + em_totals + prior

                delta = np.abs(alpha_new - alpha).sum() / max(alpha.sum(), 1.0)
                alpha = alpha_new
                total_vbem_iters += 1

                logger.debug(
                    f"  VBEM pass {reest_pass} iteration {it + 1}: "
                    f"delta={delta:.6f}"
                )

                if delta < _VBEM_CONVERGENCE_DELTA:
                    logger.debug(
                        f"  VBEM pass {reest_pass} converged after "
                        f"{it + 1} iterations"
                    )
                    break

            # After convergence, re-estimate models for next pass
            if can_reestimate and reest_pass < n_passes - 1:
                _reestimate_and_rescore(posteriors)
                logger.debug(
                    f"  Re-estimated strand: p_rna_sense={p_rna_sense:.4f}"
                )

        logger.debug(f"  VBEM total iterations: {total_vbem_iters}")

        self._converged_alpha = alpha
        self._converged_theta = alpha / alpha.sum()

        # Store re-estimated model parameters for downstream use
        if can_reestimate:
            self._em_p_rna_sense = p_rna_sense
            self._em_p_gdna_sense = p_gdna_sense
            logger.info(
                f"  EM re-estimated models: p_rna_sense={p_rna_sense:.4f}, "
                f"p_gdna_sense={p_gdna_sense:.4f}"
            )

        gdna_alpha_sum = alpha[gdna_base:gdna_base + self.num_genes].sum()
        shadow_init_sum = self.shadow_init.sum()
        logger.info(
            f"  gDNA alpha (sum) = {gdna_alpha_sum:.1f}, "
            f"shadow_init (sum) = {shadow_init_sum:.0f} / "
            f"{unique_totals.sum():.0f}"
        )

    # ------------------------------------------------------------------
    # Posterior sampling assignment (one fragment → one transcript)
    # ------------------------------------------------------------------

    def assign_ambiguous(
        self,
        em_data: EMData,
        *,
        gdna_threshold: float = 0.5,
        confidence_threshold: float = 0.95,
    ) -> None:
        """Assign each ambiguous fragment via thresholded posterior sampling.

        For each EM unit the posterior over candidates is computed from
        the converged theta.  The RNA/gDNA decision is deterministic:

        * If ``sum(posterior_RNA) >= gdna_threshold``, the unit is
          classified as RNA and one transcript is sampled from the
          RNA-renormalized posteriors.
        * Otherwise the unit is assigned to its gDNA shadow.

        For RNA assignments whose winning (RNA-normalized) posterior
        meets or exceeds ``confidence_threshold``, the count is also
        added to ``em_high_conf_counts``.

        Must be called after ``run_em()``.

        Parameters
        ----------
        gdna_threshold : float
            Minimum sum of RNA posteriors to classify a unit as RNA.
            0.0 → never assign to gDNA; 1.0 → only shadow-free units
            are RNA.
        confidence_threshold : float
            Minimum RNA-normalized posterior for an assignment to be
            considered high-confidence (default 0.95).
        """
        if em_data.n_units == 0 or em_data.n_candidates == 0:
            return

        alpha = self._converged_alpha
        if alpha is None:
            raise RuntimeError(
                "run_em() must be called before assign_ambiguous()"
            )

        offsets = em_data.offsets
        t_indices = em_data.t_indices
        log_liks = em_data.log_liks
        count_cols_arr = em_data.count_cols
        seg_lengths = np.diff(offsets)

        # Compute final posteriors using converged VBEM alpha
        eff_len = self._build_eff_len_vector()
        log_weights = digamma(alpha) - np.log(eff_len)

        log_posteriors = log_liks + log_weights[t_indices]

        # Segment-wise softmax normalization
        seg_max = np.maximum.reduceat(log_posteriors, offsets[:-1])
        log_posteriors -= np.repeat(seg_max, seg_lengths)
        posteriors = np.exp(log_posteriors)
        seg_sum = np.add.reduceat(posteriors, offsets[:-1])
        posteriors /= np.repeat(seg_sum, seg_lengths)

        # Per-unit thresholded sampling
        gdna_base = self.gdna_base_index
        rng = self._rng
        n_units = em_data.n_units
        winner_idx = np.empty(n_units, dtype=np.int64)
        winner_is_rna = np.empty(n_units, dtype=bool)
        winner_conf = np.empty(n_units, dtype=np.float64)

        for u in range(n_units):
            s, e = int(offsets[u]), int(offsets[u + 1])
            p = posteriors[s:e]
            t_seg = t_indices[s:e]

            rna_mask = t_seg < gdna_base
            has_shadow = (~rna_mask).any()
            p_rna = p[rna_mask]
            p_rna_total = float(p_rna.sum())

            if p_rna_total >= gdna_threshold or not has_shadow:
                # --- Assign to RNA: sample among RNA candidates ---
                rna_local = np.where(rna_mask)[0]
                if len(rna_local) == 1:
                    choice = 0
                    conf = 1.0
                else:
                    rna_p_norm = p_rna / p_rna_total
                    choice = rng.choice(len(rna_local), p=rna_p_norm)
                    conf = float(rna_p_norm[choice])
                winner_idx[u] = s + rna_local[choice]
                winner_is_rna[u] = True
                winner_conf[u] = conf
            else:
                # --- Assign to gDNA shadow ---
                gdna_local = np.where(~rna_mask)[0]
                if len(gdna_local) == 1:
                    choice = 0
                else:
                    p_gdna = p[gdna_local]
                    p_gdna_total = float(p_gdna.sum())
                    if p_gdna_total > 0.0:
                        p_gdna_norm = p_gdna / p_gdna_total
                    else:
                        p_gdna_norm = np.full(
                            len(gdna_local),
                            1.0 / len(gdna_local),
                            dtype=np.float64,
                        )
                    choice = int(rng.choice(len(gdna_local), p=p_gdna_norm))
                winner_idx[u] = s + gdna_local[choice]
                winner_is_rna[u] = False
                winner_conf[u] = 0.0

        winner_t = t_indices[winner_idx]
        winner_col = count_cols_arr[winner_idx]

        # --- RNA assignments: discrete 1.0 each ---
        if winner_is_rna.any():
            rna_sel = np.where(winner_is_rna)[0]
            rna_t = winner_t[rna_sel]
            rna_col = winner_col[rna_sel]
            rna_conf = winner_conf[rna_sel]
            np.add.at(self.em_counts, (rna_t, rna_col), 1.0)

            # High-confidence assignments
            hc_mask = rna_conf >= confidence_threshold
            if hc_mask.any():
                hc_sel = rna_sel[hc_mask]
                np.add.at(
                    self.em_high_conf_counts,
                    (winner_t[hc_sel], winner_col[hc_sel]),
                    1.0,
                )

            # Track confidence per transcript
            if self._em_posterior_sum is None:
                self._em_posterior_sum = np.zeros(
                    self.num_transcripts, dtype=np.float64
                )
                self._em_n_assigned = np.zeros(
                    self.num_transcripts, dtype=np.int64
                )
            np.add.at(self._em_posterior_sum, rna_t, rna_conf)
            np.add.at(self._em_n_assigned, rna_t, 1)

        # --- gDNA assignments: discrete 1.0 each ---
        gdna_mask_final = ~winner_is_rna
        if gdna_mask_final.any():
            gdna_sel = np.where(gdna_mask_final)[0]
            gene_indices = winner_t[gdna_sel] - gdna_base
            np.add.at(self.gdna_em_counts, gene_indices, 1.0)

            # Locus attribution
            locus_t = em_data.locus_t_indices[gdna_sel]
            locus_ct = em_data.locus_count_cols[gdna_sel]
            valid = locus_t >= 0
            if valid.any():
                np.add.at(
                    self.gdna_locus_counts,
                    (locus_t[valid], locus_ct[valid]),
                    1.0,
                )

    # ------------------------------------------------------------------
    # Confidence / posterior metrics
    # ------------------------------------------------------------------

    def posterior_mean(self) -> np.ndarray:
        """Per-transcript mean posterior at assignment.

        For each transcript, this is the average posterior probability
        at which a unit was assigned to it.  Returns NaN for
        transcripts with no EM assignments.

        Returns
        -------
        np.ndarray
            float64[num_transcripts]
        """
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
        """Primary transcript-level counts with identifiers.

        Columns: transcript_id, gene_id, gene_name, count,
        count_unique, count_spliced, count_em, count_high_conf,
        n_gdna, posterior_mean.
        """
        total = self.t_counts  # (N_t, 8)
        unique = self.unique_counts  # (N_t, 8)

        # Total count = sum across all 8 columns
        count_total = total.sum(axis=1)
        count_unique = unique.sum(axis=1)

        # Spliced = SPLICED_ANNOT + SPLICED_UNANNOT (both strands)
        count_spliced = total[:, list(SPLICED_COLS)].sum(axis=1)

        # EM count
        count_em = self.em_counts.sum(axis=1)

        # High-confidence count (unique + EM above threshold)
        count_high_conf = self.t_high_conf_counts.sum(axis=1)

        # Per-transcript gDNA locus attribution
        n_gdna = self.gdna_locus_counts.sum(axis=1)

        # Posterior mean
        pmean = self.posterior_mean()

        df = pd.DataFrame({
            "transcript_id": index.t_df["t_id"].values,
            "gene_id": index.t_df["g_id"].values,
            "gene_name": index.t_df["g_name"].values,
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
        """Primary gene-level counts with identifiers and gDNA summary.

        Columns: gene_id, gene_name, count, count_unique,
        count_spliced, count_em, count_high_conf, n_gdna,
        n_antisense, gdna_rate.
        """
        t_to_g = index.t_to_g_arr
        total = self.t_counts
        unique = self.unique_counts

        n_genes = index.num_genes
        g_total = np.zeros((n_genes, NUM_COUNT_COLS), dtype=np.float64)
        np.add.at(g_total, t_to_g, total)

        g_unique = np.zeros((n_genes, NUM_COUNT_COLS), dtype=np.float64)
        np.add.at(g_unique, t_to_g, unique)

        count_total = g_total.sum(axis=1)
        count_unique = g_unique.sum(axis=1)

        count_spliced = g_total[:, list(SPLICED_COLS)].sum(axis=1)
        count_em_arr = np.zeros(n_genes, dtype=np.float64)
        np.add.at(count_em_arr, t_to_g, self.em_counts.sum(axis=1))

        # High-confidence counts aggregated to gene level
        count_hc_arr = np.zeros(n_genes, dtype=np.float64)
        np.add.at(
            count_hc_arr, t_to_g, self.t_high_conf_counts.sum(axis=1)
        )

        n_gdna = self.gdna_em_counts.copy()

        # Antisense unique counts per gene
        t_antisense = unique[:, list(ANTISENSE_COLS)].sum(axis=1)
        n_antisense = np.zeros(n_genes, dtype=np.float64)
        np.add.at(n_antisense, t_to_g, t_antisense)

        # gDNA rate per gene
        denom = count_total + n_gdna
        with np.errstate(divide="ignore", invalid="ignore"):
            gdna_rate = np.where(denom > 0, n_gdna / denom, 0.0)

        df = pd.DataFrame({
            "gene_id": index.g_df["g_id"].values,
            "gene_name": index.g_df["g_name"].values,
            "count": count_total,
            "count_unique": count_unique,
            "count_spliced": count_spliced,
            "count_em": count_em_arr,
            "count_high_conf": count_hc_arr,
            "n_gdna": n_gdna,
            "n_antisense": n_antisense,
            "gdna_rate": gdna_rate,
        })
        return df

    # ------------------------------------------------------------------
    # Output - detail (long format QC breakdown)
    # ------------------------------------------------------------------

    def get_detail_df(self, index) -> pd.DataFrame:
        """Detailed counts in long format for QC.

        Each row is a nonzero (transcript, category, source) combination.
        Columns: transcript_id, gene_id, category, source, count.
        """
        t_ids = index.t_df["t_id"].values
        g_ids = index.t_df["g_id"].values

        frames = []
        for source_name, counts in (
            ("unique", self.unique_counts),
            ("em", self.em_counts),
        ):
            # Sum sense + antisense per category -> (N_t, 4)
            cat_counts = np.zeros(
                (self.num_transcripts, len(CountCategory)),
                dtype=np.float64,
            )
            for cat in CountCategory:
                sense_col = CountCol.from_category(cat, False)
                anti_col = CountCol.from_category(cat, True)
                cat_counts[:, int(cat)] = (
                    counts[:, sense_col] + counts[:, anti_col]
                )

            nz_t, nz_cat = np.nonzero(cat_counts)
            if len(nz_t) == 0:
                continue
            cat_order = [c.name.lower() for c in CountCategory]
            frames.append(
                pd.DataFrame({
                    "transcript_id": t_ids[nz_t],
                    "gene_id": g_ids[nz_t],
                    "category": pd.Categorical(
                        [CountCategory(c).name.lower() for c in nz_cat],
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

    # ------------------------------------------------------------------
    # Output - gDNA summary dict
    # ------------------------------------------------------------------

    def gdna_summary(self) -> dict:
        """Summary dict of gDNA contamination estimates."""
        return {
            "gdna_unique_count": float(self.gdna_unique_count),
            "gdna_em_count": float(self.gdna_em_count),
            "gdna_total": float(self.gdna_total),
            "gdna_contamination_rate": float(
                round(self.gdna_contamination_rate, 6)
            ),
            "rna_unique_total": float(self.unique_counts.sum()),
            "rna_em_total": float(self.em_counts.sum()),
        }
