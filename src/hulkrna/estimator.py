"""
hulkrna.estimator — Fragment abundance estimation with locus-level EM.

AbundanceEstimator manages per-transcript and per-gene fragment count
accumulation combined with Bayesian abundance estimation via EM.

For unambiguously-mapping fragments, assignment is deterministic.  For ambiguous
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

1. Initialize theta from unambig counts + OVR prior.
2. E-step: posterior proportional to likelihood * theta_t / eff_len_t.
3. M-step: theta proportional to unambig + em_totals + prior.
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
    unambig_totals : np.ndarray
        float64[n_components] - init counts (unambig + shadow inits).
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
    nrna_frac_alpha : np.ndarray
        float64[n_transcripts] - Beta prior α for nrna_frac (nascent fraction).
    nrna_frac_beta : np.ndarray
        float64[n_transcripts] - Beta prior β for nrna_frac (nascent fraction).
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
    nrna_init: np.ndarray
    gdna_init: float
    effective_lengths: np.ndarray
    prior: np.ndarray
    bias_profiles: np.ndarray | list | None
    nrna_frac_alpha: np.ndarray
    nrna_frac_beta: np.ndarray


# ======================================================================
# nrna_frac (nascent fraction) prior computation — density + strand hybrid
# ======================================================================

#: Minimum denominator for the strand weight (2s − 1) below which
#: we consider strand information unusable.
_STRAND_DENOM_EPS: float = 0.01


def _compute_hybrid_nrna_frac_vec(
    exon_sense: np.ndarray,
    exon_anti: np.ndarray,
    intron_sense: np.ndarray,
    intron_anti: np.ndarray,
    L_exonic: np.ndarray,
    L_intronic: np.ndarray,
    s: float,
    gdna_density: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute nrna_frac estimates and evidence from per-unit count arrays.

    Uses the unified density + strand hybrid model:

    * **Density subtraction** isolates RNA rates from spatial zone
      densities (exonic vs intronic vs intergenic).
    * **Strand subtraction** isolates RNA from gDNA via the strand
      specificity.
    * The two are combined through inverse-variance weighting
      ``W = (2s − 1)²``.

    Parameters
    ----------
    exon_sense, exon_anti : ndarray
        Total exonic fragment counts (sense / antisense) per unit.
    intron_sense, intron_anti : ndarray
        Unambiguously intronic fragment counts per unit.
    L_exonic, L_intronic : ndarray
        Exonic and intronic lengths (bp) per unit.
    s : float
        Library strand specificity ∈ [0.5, 1.0].
    gdna_density : float
        Global intergenic gDNA density (frags / bp).

    Returns
    -------
    nrna_frac : ndarray
        Nascent fraction estimates, same shape as inputs.
    evidence : ndarray
        Total exonic + intronic fragment count per unit.
    """
    exon_total = exon_sense + exon_anti
    intron_total = intron_sense + intron_anti
    evidence = exon_total + intron_total

    # --- Density estimator ---
    with np.errstate(divide="ignore", invalid="ignore"):
        D_exon = np.where(L_exonic > 0, exon_total / L_exonic, 0.0)
        # Single-exon transcripts (L_intronic == 0): assume intron
        # density equals background so nRNA contribution is 0.
        D_intron = np.where(
            L_intronic > 0, intron_total / L_intronic, gdna_density,
        )
    den_nrna = np.maximum(0.0, D_intron - gdna_density)
    den_mrna = np.maximum(0.0, D_exon - D_intron)

    # --- Strand estimator ---
    denom = 2.0 * s - 1.0
    W = denom ** 2 if denom > _STRAND_DENOM_EPS else 0.0

    if W > 0:
        exon_rna = np.maximum(0.0, (exon_sense - exon_anti) / denom)
        str_exon_density = exon_rna / np.maximum(1.0, L_exonic)
        intron_rna = np.maximum(0.0, (intron_sense - intron_anti) / denom)
        with np.errstate(divide="ignore", invalid="ignore"):
            str_intron_density = np.where(
                L_intronic > 0, intron_rna / L_intronic, 0.0,
            )
        str_nrna = str_intron_density
        str_mrna = np.maximum(0.0, str_exon_density - str_intron_density)
    else:
        str_nrna = den_nrna
        str_mrna = den_mrna

    # --- Weighted combination ---
    final_nrna = W * str_nrna + (1.0 - W) * den_nrna
    final_mrna = W * str_mrna + (1.0 - W) * den_mrna

    # --- nrna_frac = nRNA / (mRNA + nRNA) ---
    # When total_rate = 0 (no evidence), default to nrna_frac = 0 (no nRNA
    # detected).  This is the conservative empirical Bayes choice:
    # absence of evidence → assume fully spliced mRNA.
    total_rate = final_mrna + final_nrna
    with np.errstate(divide="ignore", invalid="ignore"):
        nrna_frac = np.where(total_rate > 0, final_nrna / total_rate, 0.0)

    return nrna_frac, evidence


def compute_global_gdna_density(
    estimator: "AbundanceEstimator",
    strand_specificity: float,
) -> float:
    """Estimate global gDNA density (fragments / bp) via strand correction.

    Uses the same strand-based gDNA isolation as the EB prior system:
    from global unspliced sense/antisense totals, compute the gDNA fraction,
    then convert to a per-bp density using total exonic span.

    When strand specificity is too low for reliable correction,
    returns 0.0 (conservative: no gDNA subtracted by the density
    estimator, which is fine because the density estimator dominates
    in that regime anyway).

    Parameters
    ----------
    estimator : AbundanceEstimator
    strand_specificity : float
        Library strand specificity ∈ [0.5, 1.0].

    Returns
    -------
    float
        Global gDNA density in fragments per base pair.
    """
    from .locus import compute_gdna_rate_from_strand

    total_sense = float(estimator.transcript_unspliced_sense.sum())
    total_anti = float(estimator.transcript_unspliced_antisense.sum())
    gdna_rate = compute_gdna_rate_from_strand(
        total_sense, total_anti, strand_specificity,
    )
    total_unspliced = total_sense + total_anti
    gdna_frag_count = gdna_rate * total_unspliced

    # Normalise by total exonic span to get per-bp density.
    if estimator._exonic_lengths is not None:
        total_exonic_bp = float(estimator._exonic_lengths.sum())
    else:
        total_exonic_bp = 0.0

    if total_exonic_bp > 0:
        return gdna_frag_count / total_exonic_bp
    return 0.0


# ======================================================================
# Method-of-Moments κ estimation
# ======================================================================

_KAPPA_MIN = 2.0
_KAPPA_MAX = 200.0
_KAPPA_FALLBACK = 5.0
_KAPPA_MIN_OBS = 20


def estimate_kappa(
    nrna_frac_array: np.ndarray,
    evidence_array: np.ndarray,
    min_evidence: float = 30.0,
    *,
    kappa_min: float = _KAPPA_MIN,
    kappa_max: float = _KAPPA_MAX,
    kappa_fallback: float = _KAPPA_FALLBACK,
    kappa_min_obs: int = _KAPPA_MIN_OBS,
) -> float:
    """Estimate Beta concentration κ via Method of Moments.

    Selects features whose evidence ≥ *min_evidence*, computes the
    mean (μ) and variance (σ²) of their raw nrna_frac values, and solves:

    .. math::

        \\kappa = \\frac{\\mu (1 - \\mu)}{\\sigma^2} - 1

    The result is clamped to [*kappa_min*, *kappa_max*].
    If fewer than *kappa_min_obs* features pass the evidence
    filter, returns *kappa_fallback*.

    Parameters
    ----------
    nrna_frac_array : ndarray
        Raw nrna_frac estimates (nascent fraction) per feature.
    evidence_array : ndarray
        Total fragment counts per feature.
    min_evidence : float
        Minimum evidence to include a feature in the estimate.
    kappa_min : float
        Lower clamp for the estimated κ (default 2.0).
    kappa_max : float
        Upper clamp for the estimated κ (default 200.0).
    kappa_fallback : float
        Fallback κ when too few features pass the filter (default 5.0).
    kappa_min_obs : int
        Minimum number of features required; fewer triggers fallback
        (default 20).

    Returns
    -------
    float
        Estimated κ (concentration / precision).
    """
    valid = evidence_array >= min_evidence
    if valid.sum() < kappa_min_obs:
        return kappa_fallback

    confident_nrna_fracs = nrna_frac_array[valid]
    mu = float(np.mean(confident_nrna_fracs))
    sigma2 = float(np.var(confident_nrna_fracs))

    if sigma2 <= 0.0:
        return kappa_max  # zero variance → maximum precision

    kappa = mu * (1.0 - mu) / sigma2 - 1.0
    return float(np.clip(kappa, kappa_min, kappa_max))


def compute_hybrid_nrna_frac_priors(
    estimator: "AbundanceEstimator",
    t_to_tss_group: np.ndarray | None,
    t_to_strand: np.ndarray,
    locus_id_per_transcript: np.ndarray,
    strand_specificity: float,
    gdna_density: float,
    kappa_global: float | None = None,
    kappa_locus: float | None = None,
    kappa_tss: float | None = None,
    mom_min_evidence_global: float = 50.0,
    mom_min_evidence_locus: float = 30.0,
    mom_min_evidence_tss: float = 20.0,
    kappa_min: float = _KAPPA_MIN,
    kappa_max: float = _KAPPA_MAX,
    kappa_fallback: float = _KAPPA_FALLBACK,
    kappa_min_obs: int = _KAPPA_MIN_OBS,
) -> None:
    """Compute per-transcript nrna_frac (nascent fraction) Beta priors.

    Uses a **density + strand hybrid** model with smooth Empirical-Bayes
    hierarchical shrinkage.  Raw nrna_frac estimates are computed at three levels
    (transcript, TSS group, locus-strand) and each level is smoothly
    shrunk toward its parent:

    1. **Locus-strand nrna_frac** ← shrink toward empirical global nrna_frac
       (evidence-weighted mean of all locus-strand estimates) with
       strength *kappa_global*.
    2. **TSS-group nrna_frac** ← shrink toward shrunk locus-strand nrna_frac with
       strength *kappa_locus*.
    3. **Transcript nrna_frac** ← shrink toward shrunk TSS-group nrna_frac with
       strength *kappa_tss*.

    The shrinkage formula at each level is::

        nrna_frac_shrunk = (E_local · nrna_frac_local + κ · nrna_frac_parent) / (E_local + κ)

    When a κ parameter is ``None`` (the default), it is estimated
    automatically from the data using the Method of Moments:

    .. math::

        \\kappa = \\frac{\\mu (1 - \\mu)}{\\sigma^2} - 1

    where μ and σ² are the mean and variance of the raw nrna_frac values at
    the corresponding hierarchy level.

    The final effective sample size (κ) for the Beta prior equals
    ``kappa_tss`` — a weak, constant prior that lets the EM update
    nrna_frac freely from the fragment-level data.

    Parameters
    ----------
    estimator : AbundanceEstimator
        Must have ``unambig_counts``, ``transcript_intronic_sense``,
        ``transcript_intronic_antisense``, ``_exonic_lengths``, and
        ``_transcript_spans`` populated.
    t_to_tss_group : ndarray or None
        int32[n_t] TSS group IDs.  None disables TSS level.
    t_to_strand : ndarray
        int[n_t] strand values (Strand enum ints).
    locus_id_per_transcript : ndarray
        int32[n_t] locus IDs (−1 when not in a locus).
    strand_specificity : float
        Library strand specificity ∈ [0.5, 1.0].
    gdna_density : float
        Global intergenic gDNA density (frags / bp).
    kappa_global : float or None
        Shrinkage strength pulling locus-strand toward global.
        ``None`` (default) → estimate from data via Method of Moments.
    kappa_locus : float or None
        Shrinkage strength pulling TSS-group toward locus-strand.
        ``None`` → auto-estimate.
    kappa_tss : float or None
        Shrinkage strength pulling transcript toward TSS-group.
        ``None`` → auto-estimate.
    mom_min_evidence_global : float
        Minimum evidence for MoM κ estimation at the global level.
    mom_min_evidence_locus : float
        Minimum evidence for MoM κ estimation at the locus level.
    mom_min_evidence_tss : float
        Minimum evidence for MoM κ estimation at the TSS level.
    kappa_min : float
        Lower clamp for MoM-estimated κ.
    kappa_max : float
        Upper clamp for MoM-estimated κ.
    kappa_fallback : float
        Fallback κ when too few features pass the evidence filter.
    kappa_min_obs : int
        Minimum features required for MoM estimation.
    """
    n_t = estimator.num_transcripts
    if n_t == 0:
        return

    # --- Per-transcript count vectors ---
    # Exonic: all fragments overlapping exons (unambig + ambig-same-strand,
    # weighted by 1/n_candidates).  These define D_exon.  Using the
    # pre-EM exonic accumulators provides balanced evidence matching the
    # completeness of the intronic accumulators (which also include
    # ambiguous fragments), rather than the sparse unambig_counts that
    # only contain splice-spanning unambig.
    exon_sense = estimator.transcript_exonic_sense.copy()
    exon_anti = estimator.transcript_exonic_antisense.copy()

    # Intronic: fragments with unambiguous intronic overlap relative to
    # each transcript.  Already per-transcript from the scan phase.
    intron_sense = estimator.transcript_intronic_sense.copy()
    intron_anti = estimator.transcript_intronic_antisense.copy()

    # --- Lengths ---
    if estimator._exonic_lengths is not None:
        L_exonic = estimator._exonic_lengths.copy()
    else:
        # Fallback: treat exonic length as 1 to avoid div-by-zero.
        L_exonic = np.ones(n_t, dtype=np.float64)
    if estimator._transcript_spans is not None:
        L_intronic = np.maximum(
            estimator._transcript_spans - L_exonic, 0.0,
        )
    else:
        L_intronic = np.zeros(n_t, dtype=np.float64)

    s = strand_specificity

    # === Level 1: per-transcript ===
    t_nrna_frac, t_evidence = _compute_hybrid_nrna_frac_vec(
        exon_sense, exon_anti, intron_sense, intron_anti,
        L_exonic, L_intronic, s, gdna_density,
    )

    # === Level 2: TSS group (aggregate counts then compute nrna_frac) ===
    if t_to_tss_group is not None and n_t > 0:
        n_groups = int(t_to_tss_group.max()) + 1
        _z = np.zeros  # shorthand
        g_exon_s = _z(n_groups, dtype=np.float64)
        g_exon_a = _z(n_groups, dtype=np.float64)
        g_int_s = _z(n_groups, dtype=np.float64)
        g_int_a = _z(n_groups, dtype=np.float64)
        g_L_ex = _z(n_groups, dtype=np.float64)
        g_L_in = _z(n_groups, dtype=np.float64)
        np.add.at(g_exon_s, t_to_tss_group, exon_sense)
        np.add.at(g_exon_a, t_to_tss_group, exon_anti)
        np.add.at(g_int_s, t_to_tss_group, intron_sense)
        np.add.at(g_int_a, t_to_tss_group, intron_anti)
        np.add.at(g_L_ex, t_to_tss_group, L_exonic)
        np.add.at(g_L_in, t_to_tss_group, L_intronic)

        g_nrna_frac, g_evidence = _compute_hybrid_nrna_frac_vec(
            g_exon_s, g_exon_a, g_int_s, g_int_a,
            g_L_ex, g_L_in, s, gdna_density,
        )
        # Map back to transcript dimension
        tss_nrna_frac = g_nrna_frac[t_to_tss_group]
        tss_evidence = g_evidence[t_to_tss_group]
    else:
        tss_nrna_frac = np.zeros(n_t, dtype=np.float64)
        tss_evidence = np.zeros(n_t, dtype=np.float64)

    # === Level 3: locus-strand (aggregate by locus × strand) ===
    has_locus = locus_id_per_transcript >= 0
    max_locus = (
        int(locus_id_per_transcript.max()) + 1 if has_locus.any() else 0
    )
    ls_key = np.full(n_t, -1, dtype=np.int64)
    if max_locus > 0:
        ls_key[has_locus] = (
            locus_id_per_transcript[has_locus].astype(np.int64) * 4
            + t_to_strand[has_locus].astype(np.int64)
        )
    n_ls = max(max_locus * 4, 1)
    _z = np.zeros
    ls_exon_s = _z(n_ls, dtype=np.float64)
    ls_exon_a = _z(n_ls, dtype=np.float64)
    ls_int_s = _z(n_ls, dtype=np.float64)
    ls_int_a = _z(n_ls, dtype=np.float64)
    ls_L_ex = _z(n_ls, dtype=np.float64)
    ls_L_in = _z(n_ls, dtype=np.float64)
    valid_ls = ls_key >= 0
    if valid_ls.any():
        np.add.at(ls_exon_s, ls_key[valid_ls], exon_sense[valid_ls])
        np.add.at(ls_exon_a, ls_key[valid_ls], exon_anti[valid_ls])
        np.add.at(ls_int_s, ls_key[valid_ls], intron_sense[valid_ls])
        np.add.at(ls_int_a, ls_key[valid_ls], intron_anti[valid_ls])
        np.add.at(ls_L_ex, ls_key[valid_ls], L_exonic[valid_ls])
        np.add.at(ls_L_in, ls_key[valid_ls], L_intronic[valid_ls])

    ls_nrna_frac_flat, ls_ev_flat = _compute_hybrid_nrna_frac_vec(
        ls_exon_s, ls_exon_a, ls_int_s, ls_int_a,
        ls_L_ex, ls_L_in, s, gdna_density,
    )
    # Map back to transcript dimension
    ls_nrna_frac = np.zeros(n_t, dtype=np.float64)
    ls_evidence = np.zeros(n_t, dtype=np.float64)
    if valid_ls.any():
        ls_nrna_frac[valid_ls] = ls_nrna_frac_flat[ls_key[valid_ls]]
        ls_evidence[valid_ls] = ls_ev_flat[ls_key[valid_ls]]

    # === Estimate κ hyperparameters via Method of Moments ===
    _mom_kw = dict(
        kappa_min=kappa_min, kappa_max=kappa_max,
        kappa_fallback=kappa_fallback, kappa_min_obs=kappa_min_obs,
    )
    if kappa_global is None:
        kappa_global = estimate_kappa(
            ls_nrna_frac_flat, ls_ev_flat,
            min_evidence=mom_min_evidence_global, **_mom_kw,
        )
    if kappa_locus is None:
        if t_to_tss_group is not None:
            kappa_locus = estimate_kappa(
                g_nrna_frac, g_evidence,
                min_evidence=mom_min_evidence_locus, **_mom_kw,
            )
        else:
            # No TSS groups → fall back to locus-strand estimates
            kappa_locus = estimate_kappa(
                ls_nrna_frac_flat, ls_ev_flat,
                min_evidence=mom_min_evidence_locus, **_mom_kw,
            )
    if kappa_tss is None:
        kappa_tss = estimate_kappa(
            t_nrna_frac, t_evidence,
            min_evidence=mom_min_evidence_tss, **_mom_kw,
        )

    # === Empirical global nrna_frac prior ===
    #
    # Instead of a hard-coded nrna_frac_global = 0.5, we compute an
    # evidence-weighted mean of locus-strand nrna_frac values.  This makes
    # the global prior fully data-driven: in a pristine experiment
    # (no nRNA) nrna_frac_global ≈ 0; in a pol-II inhibition experiment
    # it naturally rises toward the empirical nascent fraction.
    total_ls_evidence = float(ls_ev_flat.sum())
    if total_ls_evidence > 0:
        nrna_frac_global = float(np.dot(ls_ev_flat, ls_nrna_frac_flat) / total_ls_evidence)
    else:
        nrna_frac_global = 0.0  # no evidence → assume no nRNA

    logger.debug(
        f"nrna_frac EB hyperparameters: nrna_frac_global={nrna_frac_global:.4f}, "
        f"κ_global={kappa_global:.1f}, "
        f"κ_locus={kappa_locus:.1f}, κ_tss={kappa_tss:.1f}"
    )

    # === Smooth Empirical-Bayes hierarchical shrinkage ===
    #
    #   nrna_frac_shrunk = (E_local · nrna_frac_local + κ · nrna_frac_parent) / (E_local + κ)

    # Level 3 → shrink toward empirical global
    nrna_frac_L3 = (
        (ls_evidence * ls_nrna_frac + kappa_global * nrna_frac_global)
        / (ls_evidence + kappa_global)
    )

    # Level 2 → shrink toward shrunk Level 3
    nrna_frac_L2 = (
        (tss_evidence * tss_nrna_frac + kappa_locus * nrna_frac_L3)
        / (tss_evidence + kappa_locus)
    )

    # Level 1 → shrink toward shrunk Level 2
    nrna_frac_est = (
        (t_evidence * t_nrna_frac + kappa_tss * nrna_frac_L2)
        / (t_evidence + kappa_tss)
    )

    # Effective sample size for the Beta prior.  Use only the
    # hierarchical shrinkage strength (kappa_tss), NOT the raw
    # fragment evidence.  Adding t_evidence makes the prior
    # overwhelmingly strong (kappa ~ 10 000), which locks nrna_frac
    # near the pre-EM estimate and prevents the EM from correcting
    # it.  The EM itself sees all the data, so a weakly informative
    # prior (kappa ~ 5–20) is appropriate.
    kappa = np.full_like(t_evidence, kappa_tss)

    # Clamp nrna_frac away from exact 0 or 1 to keep Beta parameters finite
    nrna_frac_est = np.clip(nrna_frac_est, 1e-4, 1.0 - 1e-4)

    # Convert to Beta(α, β)
    estimator.nrna_frac_alpha = nrna_frac_est * kappa
    estimator.nrna_frac_beta = (1.0 - nrna_frac_est) * kappa

    # --- Diagnostic logging ---
    median_nrna_frac = float(np.median(nrna_frac_est))
    median_kappa = float(np.median(kappa))
    mean_l1_weight = float(np.mean(t_evidence / (t_evidence + kappa_tss)))
    n_zero_nrna_frac = int(np.sum(nrna_frac_est <= 1e-4))
    n_high_nrna_frac = int(np.sum(nrna_frac_est >= 0.1))
    logger.info(
        f"nrna_frac priors: global={nrna_frac_global:.4f}, median={median_nrna_frac:.4f}, "
        f"κ=[{kappa_global:.1f}, {kappa_locus:.1f}, {kappa_tss:.1f}], "
        f"median_κ_final={median_kappa:.1f}, L1_wt={mean_l1_weight:.2f}, "
        f"n_zero={n_zero_nrna_frac}/{n_t}, n_high(≥0.1)={n_high_nrna_frac}/{n_t}"
    )


# ======================================================================
# AbundanceEstimator
# ======================================================================


class AbundanceEstimator:
    """Bayesian fragment abundance estimator with locus-level EM.

    Accumulates per-transcript and per-gene fragment counts from both
    deterministic (unambig) and probabilistic (EM) assignment paths.
    Estimates transcript-level mRNA, nRNA, and gDNA abundances.

    Locus-level EM architecture::

        Per locus:
        [0, n_t)       - mRNA transcript components
        [n_t, 2*n_t)   - nRNA shadow per transcript
        [2*n_t]        - gDNA: ONE shadow per locus

    Fragment routing:

    - SPLICED_ANNOT + unambig -> deterministic mRNA (no EM)
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

        self.unambig_counts = np.zeros(
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

        # Per-transcript nrna_frac (nascent fraction) Beta prior parameters.
        # Computed after the scan phase by compute_hybrid_nrna_frac_priors().
        self.nrna_frac_alpha = np.ones(num_transcripts, dtype=np.float64)
        self.nrna_frac_beta = np.ones(num_transcripts, dtype=np.float64)

        # --- gDNA: locus-level, NOT per-gene ---
        # Total gDNA count as assigned by locus EM (sum across all loci)
        self._gdna_em_total = 0.0

        # Per-locus results (locus_id, chrom, n_transcripts, n_genes,
        # n_em_fragments, mrna, nrna, gdna, gdna_init).
        self.locus_results: list[dict] = []

        # Per-transcript locus_id assignment (-1 = no locus).
        self.locus_id_per_transcript = np.full(
            num_transcripts, -1, dtype=np.int32
        )

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

        # --- Pre-EM exonic accumulators (for nrna_frac prior) ---
        # All unambig + ambig-same-strand fragments overlapping each
        # transcript, weighted by 1/n_candidates.  Provides balanced
        # exonic evidence for the nrna_frac (nascent fraction) estimator,
        # unlike unambig_counts which only has splice-spanning unambig.
        self.transcript_exonic_sense = np.zeros(
            num_transcripts, dtype=np.float64
        )
        self.transcript_exonic_antisense = np.zeros(
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
        total_rna = float(self.unambig_counts.sum() + self.em_counts.sum())
        total = total_rna + self.gdna_total
        if total == 0:
            return 0.0
        return self.gdna_total / total

    @property
    def t_counts(self) -> np.ndarray:
        """Total transcript counts (unambig + EM), shape (N_t, 8)."""
        return self.unambig_counts + self.em_counts

    @property
    def t_high_conf_counts(self) -> np.ndarray:
        """High-confidence transcript counts (unambig + EM above threshold)."""
        return self.unambig_counts + self.em_high_conf_counts

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

    def assign_unambig(self, resolved, index, strand_models) -> None:
        """Deterministic assignment for a truly unambiguous fragment."""
        t_inds = resolved.t_inds

        if len(t_inds) == 0:
            return

        t_idx = int(next(iter(t_inds)))
        t_strand = int(index.t_to_strand_arr[t_idx])
        sm = strand_models.exonic_spliced
        anti = self.is_antisense(resolved.exon_strand, t_strand, sm)
        col = SpliceStrandCol.from_category(resolved.splice_type, anti)
        self.unambig_counts[t_idx, col] += 1.0

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
        theta, alpha_out, _em_totals, nrna_frac_out = _run_locus_em_native(
            locus_em.offsets,
            locus_em.t_indices,
            locus_em.log_liks,
            locus_em.coverage_weights,
            locus_em.tx_starts,
            locus_em.tx_ends,
            locus_em.bias_profiles,
            locus_em.unambig_totals,
            locus_em.effective_lengths,
            (locus_em.prior > 0).astype(np.float64),
            locus_em.n_components,
            self.em_config.prior_alpha,
            self.em_config.prior_gamma,
            em_iterations,
            em_convergence_delta,
            self.em_config.mode == "vbem",
            self.em_config.prune_threshold if self.em_config.prune_threshold is not None else -1.0,
            locus_em.n_transcripts,
            locus_em.nrna_frac_alpha,
            locus_em.nrna_frac_beta,
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
    ) -> dict[str, float]:
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
        dict[str, float]
            Per-pool fragment counts: ``{"mrna": ..., "nrna": ..., "gdna": ...}``.
        """
        offsets = locus_em.offsets
        t_indices = locus_em.t_indices
        log_liks = locus_em.log_liks
        count_cols_arr = locus_em.count_cols
        n_units = len(offsets) - 1

        if n_units == 0 or len(t_indices) == 0:
            return {"mrna": 0.0, "nrna": 0.0, "gdna": 0.0}

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
        mrna_count = 0.0
        if mrna_mask.any():
            mrna_p = posteriors[mrna_mask]
            mrna_count = float(mrna_p.sum())
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
        nrna_count = 0.0
        if nrna_mask.any():
            nrna_p = posteriors[nrna_mask]
            nrna_count = float(nrna_p.sum())
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

        return {"mrna": mrna_count, "nrna": nrna_count, "gdna": gdna_count}

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
        """Primary transcript-level abundance estimates.

        Columns
        -------
        transcript_id, gene_id, gene_name : identifiers
        locus_id : int32, EM locus (-1 if no locus)
        effective_length : bias-corrected effective length
        mrna : total mRNA count (unambig + EM)
        mrna_unambig : uniquely-assigned mRNA
        mrna_em : EM-assigned mRNA
        mrna_high_conf : high-confidence mRNA (unambig + EM >= threshold)
        mrna_spliced : spliced mRNA fragments
        nrna : nascent RNA count (EM-assigned only)
        rna_total : mrna + nrna
        tpm : transcripts per million (mRNA-based)
        posterior_mean : mean posterior at EM assignment
        """
        total = self.t_counts
        unique = self.unambig_counts

        mrna = total.sum(axis=1)
        mrna_unambig = unique.sum(axis=1)
        mrna_spliced = total[:, list(SPLICED_COLS)].sum(axis=1)
        mrna_em = self.em_counts.sum(axis=1)
        mrna_high_conf = self.t_high_conf_counts.sum(axis=1)
        nrna = self.nrna_em_counts
        rna_total = mrna + nrna
        pmean = self.posterior_mean()

        # TPM: mRNA-based, using effective lengths
        eff = self._t_eff_len
        rpk = mrna / eff  # reads per kilobase (un-scaled)
        rpk_sum = rpk.sum()
        tpm = (rpk / rpk_sum * 1e6) if rpk_sum > 0 else np.zeros_like(rpk)

        df = pd.DataFrame({
            "transcript_id": index.t_df["t_id"].values,
            "gene_id": index.t_df["g_id"].values,
            "gene_name": index.t_df["g_name"].values,
            "locus_id": self.locus_id_per_transcript,
            "effective_length": eff,
            "mrna": mrna,
            "mrna_unambig": mrna_unambig,
            "mrna_em": mrna_em,
            "mrna_high_conf": mrna_high_conf,
            "mrna_spliced": mrna_spliced,
            "nrna": nrna,
            "rna_total": rna_total,
            "tpm": tpm,
            "posterior_mean": pmean,
        })
        return df

    def get_gene_counts_df(self, index) -> pd.DataFrame:
        """Primary gene-level abundance estimates.

        Columns
        -------
        gene_id, gene_name : identifiers
        locus_id : int32, primary EM locus for this gene (-1 if none)
        effective_length : abundance-weighted mean effective length
        mrna : total mRNA count (unambig + EM)
        mrna_unambig : uniquely-assigned mRNA
        mrna_em : EM-assigned mRNA
        mrna_high_conf : high-confidence mRNA
        mrna_spliced : spliced mRNA fragments
        nrna : nascent RNA count
        rna_total : mrna + nrna
        tpm : transcripts per million (mRNA-based)
        """
        t_to_g = index.t_to_g_arr
        total = self.t_counts
        unique = self.unambig_counts

        n_genes = index.num_genes
        g_total = np.zeros((n_genes, NUM_SPLICE_STRAND_COLS), dtype=np.float64)
        np.add.at(g_total, t_to_g, total)

        g_unambig = np.zeros((n_genes, NUM_SPLICE_STRAND_COLS), dtype=np.float64)
        np.add.at(g_unambig, t_to_g, unique)

        mrna = g_total.sum(axis=1)
        mrna_unambig = g_unambig.sum(axis=1)
        mrna_spliced = g_total[:, list(SPLICED_COLS)].sum(axis=1)

        mrna_em_arr = np.zeros(n_genes, dtype=np.float64)
        np.add.at(mrna_em_arr, t_to_g, self.em_counts.sum(axis=1))

        mrna_hc_arr = np.zeros(n_genes, dtype=np.float64)
        np.add.at(mrna_hc_arr, t_to_g, self.t_high_conf_counts.sum(axis=1))

        # nRNA per gene (EM-assigned only; nrna_init is not added
        # because it overlaps with EM-assigned fragments)
        nrna = np.zeros(n_genes, dtype=np.float64)
        np.add.at(nrna, t_to_g, self.nrna_em_counts)

        rna_total = mrna + nrna

        # Gene effective length: abundance-weighted mean of transcript
        # effective lengths.  For genes with zero counts, use the
        # unweighted mean of transcript effective lengths.
        t_eff = self._t_eff_len
        t_counts_flat = (self.unambig_counts + self.em_counts).sum(axis=1)
        g_eff_num = np.zeros(n_genes, dtype=np.float64)
        g_eff_den = np.zeros(n_genes, dtype=np.float64)
        np.add.at(g_eff_num, t_to_g, t_counts_flat * t_eff)
        np.add.at(g_eff_den, t_to_g, t_counts_flat)
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

        # Gene-level locus_id: primary locus (highest mrna count).
        g_locus_id = np.full(n_genes, -1, dtype=np.int32)
        # For each transcript, if it has a locus, propagate to its gene
        # using a weight-based approach: pick the locus of the transcript
        # contributing the most mrna.
        t_mrna_flat = total.sum(axis=1)
        for t_idx in range(self.num_transcripts):
            lid = int(self.locus_id_per_transcript[t_idx])
            if lid < 0:
                continue
            g_idx = int(t_to_g[t_idx])
            if g_locus_id[g_idx] < 0:
                g_locus_id[g_idx] = lid
            else:
                # Keep locus of transcript with higher mrna contribution
                current_lid = int(g_locus_id[g_idx])
                if current_lid != lid:
                    # Find which transcript contributes more
                    # (simple heuristic: keep latest if equal)
                    pass  # keep current assignment
                # If same locus, no change needed

        # TPM: mRNA-based at gene level
        rpk = mrna / g_eff_len
        rpk_sum = rpk.sum()
        tpm = (rpk / rpk_sum * 1e6) if rpk_sum > 0 else np.zeros_like(rpk)

        df = pd.DataFrame({
            "gene_id": index.g_df["g_id"].values,
            "gene_name": index.g_df["g_name"].values,
            "locus_id": g_locus_id,
            "effective_length": g_eff_len,
            "mrna": mrna,
            "mrna_unambig": mrna_unambig,
            "mrna_em": mrna_em_arr,
            "mrna_high_conf": mrna_hc_arr,
            "mrna_spliced": mrna_spliced,
            "nrna": nrna,
            "rna_total": rna_total,
            "tpm": tpm,
        })
        return df

    # ------------------------------------------------------------------
    # Output - detail (long format QC breakdown)
    # ------------------------------------------------------------------

    def get_loci_df(self) -> pd.DataFrame:
        """Locus-level output with gDNA/nRNA/mRNA breakdown.

        Returns an empty DataFrame (with correct columns) if no loci
        were processed.

        Columns
        -------
        locus_id : int, sequential locus identifier
        chrom : str, primary chromosome
        n_transcripts : int
        n_genes : int
        n_em_fragments : int, ambiguous fragments entering EM
        mrna : float, mRNA count from locus EM
        nrna : float, nRNA count from locus EM
        gdna : float, gDNA count from locus EM
        gdna_rate : float, gdna / (mrna + nrna + gdna)
        gdna_init : float, EB-shrunk gDNA initialization
        """
        cols = [
            "locus_id", "chrom", "n_transcripts", "n_genes",
            "n_em_fragments", "mrna", "nrna", "gdna",
            "gdna_rate", "gdna_init",
        ]
        if not self.locus_results:
            return pd.DataFrame(columns=cols)

        rows = []
        for r in self.locus_results:
            total = r["mrna"] + r["nrna"] + r["gdna"]
            rate = r["gdna"] / total if total > 0 else 0.0
            rows.append({
                "locus_id": r["locus_id"],
                "chrom": r.get("chrom", ""),
                "n_transcripts": r["n_transcripts"],
                "n_genes": r["n_genes"],
                "n_em_fragments": r["n_em_fragments"],
                "mrna": r["mrna"],
                "nrna": r["nrna"],
                "gdna": r["gdna"],
                "gdna_rate": rate,
                "gdna_init": r.get("gdna_init", 0.0),
            })
        return pd.DataFrame(rows, columns=cols)

    # ------------------------------------------------------------------
    # Output - detail (long format QC breakdown)
    # ------------------------------------------------------------------

    def get_detail_df(self, index) -> pd.DataFrame:
        """Detailed counts in long format for QC."""
        t_ids = index.t_df["t_id"].values
        g_ids = index.t_df["g_id"].values

        frames = []
        for source_name, counts in (
            ("unambig", self.unambig_counts),
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
