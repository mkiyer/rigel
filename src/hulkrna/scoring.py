"""hulkrna.scoring — Pure scoring functions and pre-computed scoring context.

Every function in this module takes explicit arguments (no closures,
no captured state).  The ``ScoringContext`` dataclass holds all
pre-computed parameters that were previously captured as locals inside
``_scan_and_build_em_data``'s nested closures.

Design contract:
    - All functions are stateless pure functions.
    - ``ScoringContext`` maps directly to a C struct.
    - Every function is a candidate for Cython / C acceleration.
"""

import math
from dataclasses import dataclass

import numpy as np

from .categories import SpliceType
from .frag_length_model import _TAIL_DECAY_LP
from .types import Strand

# ---------------------------------------------------------------------------
# Constants (single source of truth for all scoring/penalty values)
# ---------------------------------------------------------------------------

#: Floor value for log-safe clamping to avoid log(0).
LOG_SAFE_FLOOR = 1e-10

#: Pre-computed log(0.5) — used for uninformative strand log-probabilities.
LOG_HALF = math.log(0.5)

#: Default gDNA splice penalty for SPLICED_UNANNOT fragments.
DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT = 0.01

#: Default overhang alpha: each base of overhang reduces probability by 100×.
DEFAULT_OVERHANG_ALPHA = 0.01

#: Default mismatch alpha: each edit-distance mismatch (NM tag) reduces
#: probability by 10×.
DEFAULT_MISMATCH_ALPHA = 0.1

# Int constants for hot-path scoring (avoid enum construction overhead).
STRAND_POS = int(Strand.POS)         # 1
STRAND_NEG = int(Strand.NEG)         # 2
SPLICE_UNSPLICED = int(SpliceType.UNSPLICED)           # 0
SPLICE_UNANNOT = int(SpliceType.SPLICED_UNANNOT)       # 1
SPLICE_ANNOT = int(SpliceType.SPLICED_ANNOT)           # 2

#: Default gDNA splice penalties per SpliceType (int keys for fast lookup).
GDNA_SPLICE_PENALTIES = {
    SPLICE_UNSPLICED: 1.0,
    SPLICE_UNANNOT: DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    SPLICE_ANNOT: 0.0,  # never used — spliced fragments don't see gDNA
}


def overhang_alpha_to_log_penalty(alpha: float) -> float:
    """Convert user-facing alpha to internal log-penalty.

    Parameters
    ----------
    alpha : float
        Per-base overhang penalty in [0, 1].
        - 0.0 → hard binary gate (any overhang → −∞)
        - 0.01 → aggressive (−4.605 per base)
        - 1.0 → no penalty (off)

    Returns
    -------
    float
        ``log(alpha)`` or ``−inf`` when alpha ≤ 0.
    """
    if alpha <= 0.0:
        return -np.inf
    return np.log(alpha)


#: Default overhang log-penalty.
DEFAULT_OVERHANG_LOG_PENALTY = overhang_alpha_to_log_penalty(
    DEFAULT_OVERHANG_ALPHA
)

#: Default mismatch log-penalty.
DEFAULT_MISMATCH_LOG_PENALTY = overhang_alpha_to_log_penalty(
    DEFAULT_MISMATCH_ALPHA
)


# ---------------------------------------------------------------------------
# ScoringContext — replaces 15+ closure-captured locals
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ScoringContext:
    """Pre-computed scoring parameters, built once per pipeline run.

    Replaces the 15+ captured locals in the old closure-based design.
    Every scoring function receives this as its first argument.
    Maps directly to a C struct for future acceleration.
    """

    # Strand model: RNA
    log_p_sense: float
    log_p_antisense: float
    anti_flag: bool              # True when protocol is RF (p_sense < 0.5)

    # Strand model: intergenic (gDNA)
    ig_p: float                  # intergenic strand raw probability
    log_ig_p: float              # log(ig_p)

    # Penalty parameters
    overhang_log_penalty: float
    mismatch_log_penalty: float
    gdna_splice_penalties: dict  # int-keyed

    # Fragment-length LUT (pre-finalized)
    fl_log_prob: np.ndarray | None  # numpy array or None
    fl_max_size: int
    fl_tail_base: float

    # Index arrays (borrowed references, never copied)
    t_strand_arr: np.ndarray     # int8[n_transcripts]
    g_strand_arr: np.ndarray     # int8[n_genes]
    t_to_g: np.ndarray           # int32[n_transcripts]
    nrna_base: int               # offset for nRNA indices in CSR

    @staticmethod
    def from_models(
        strand_models,
        frag_length_models,
        index,
        counter,
        *,
        overhang_log_penalty: float | None = None,
        mismatch_log_penalty: float | None = None,
        gdna_splice_penalties: dict | None = None,
    ) -> "ScoringContext":
        """Build a ScoringContext from trained models and index.

        Parameters
        ----------
        strand_models : StrandModels
        frag_length_models : FragmentLengthModels
        index : HulkIndex
        counter : AbundanceEstimator
        overhang_log_penalty : float or None
        mismatch_log_penalty : float or None
        gdna_splice_penalties : dict or None
        """
        rna_sm = strand_models._best_rna_model()
        p_sense = (
            rna_sm._cached_p_sense
            if rna_sm._finalized
            else rna_sm.p_r1_sense
        )
        p_antisense = 1.0 - p_sense

        ig_sm = strand_models.intergenic
        ig_p = (
            ig_sm._cached_p_sense
            if ig_sm._finalized
            else ig_sm.p_r1_sense
        )

        fl_model = frag_length_models.global_model
        fl_log_prob = fl_model._log_prob  # numpy array or None
        fl_max_size = fl_model.max_size
        fl_tail_base: float = getattr(fl_model, "_tail_base", 0.0)

        return ScoringContext(
            log_p_sense=math.log(max(p_sense, LOG_SAFE_FLOOR)),
            log_p_antisense=math.log(max(p_antisense, LOG_SAFE_FLOOR)),
            anti_flag=p_sense < 0.5,
            ig_p=ig_p,
            log_ig_p=math.log(max(ig_p, LOG_SAFE_FLOOR)),
            overhang_log_penalty=(
                overhang_log_penalty
                if overhang_log_penalty is not None
                else DEFAULT_OVERHANG_LOG_PENALTY
            ),
            mismatch_log_penalty=(
                mismatch_log_penalty
                if mismatch_log_penalty is not None
                else DEFAULT_MISMATCH_LOG_PENALTY
            ),
            gdna_splice_penalties=(
                gdna_splice_penalties or GDNA_SPLICE_PENALTIES
            ),
            fl_log_prob=fl_log_prob,
            fl_max_size=fl_max_size,
            fl_tail_base=fl_tail_base,
            t_strand_arr=index.t_to_strand_arr,
            g_strand_arr=index.g_to_strand_arr,
            t_to_g=index.t_to_g_arr,
            nrna_base=counter.nrna_base_index,
        )


# ---------------------------------------------------------------------------
# Pure scoring functions (C-translatable)
# ---------------------------------------------------------------------------


def frag_len_log_lik(ctx: ScoringContext, flen: int) -> float:
    """Fragment-length log-likelihood (0.0 when unavailable)."""
    if flen <= 0:
        return 0.0
    if ctx.fl_log_prob is not None:
        if flen <= ctx.fl_max_size:
            return float(ctx.fl_log_prob[flen])
        # Exponential tail decay beyond max_size
        return ctx.fl_tail_base + (flen - ctx.fl_max_size) * _TAIL_DECAY_LP
    return 0.0


def strand_log_probs(
    ctx: ScoringContext, exon_strand: int,
) -> tuple[float, float]:
    """Return (log_same, log_diff) strand log-probabilities.

    When exon_strand is informative (POS or NEG), returns the
    sense/antisense log-probs from the RNA strand model.
    Otherwise, returns LOG_HALF for both.
    """
    if exon_strand == 1 or exon_strand == 2:
        return ctx.log_p_sense, ctx.log_p_antisense
    return LOG_HALF, LOG_HALF


def score_mrna_candidate(
    ctx: ScoringContext,
    exon_strand: int,
    t_strand: int,
    splice_type: int,
    frag_len_ll: float,
    overhang_bp: int,
    nm: int,
) -> tuple[float, int]:
    """Compute data log-likelihood and count column for one mRNA candidate.

    Returns (log_lik, count_col_idx).
    """
    log_same, log_diff = strand_log_probs(ctx, exon_strand)
    log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0

    if exon_strand == t_strand:
        log_strand = log_same
        anti = ctx.anti_flag
    else:
        log_strand = log_diff
        anti = not ctx.anti_flag

    log_lik = (
        log_strand
        + frag_len_ll
        + overhang_bp * ctx.overhang_log_penalty
        + log_nm
    )
    count_col_idx = splice_type * 2 + int(anti)
    return log_lik, count_col_idx


def score_nrna_candidate(
    ctx: ScoringContext,
    exon_strand: int,
    gene_strand: int,
    frag_len_ll: float,
    overhang_bp: int,
    nm: int,
) -> float:
    """Compute data log-likelihood for a nascent RNA (nRNA) candidate.

    Returns log_lik (float).
    """
    log_same, log_diff = strand_log_probs(ctx, exon_strand)
    log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0

    if exon_strand == gene_strand:
        log_strand = log_same
    else:
        log_strand = log_diff

    return (
        log_strand
        + frag_len_ll
        + overhang_bp * ctx.overhang_log_penalty
        + log_nm
    )


def score_gdna_candidate(
    ctx: ScoringContext,
    exon_strand: int,
    splice_type: int,
    genomic_footprint: int,
) -> float:
    """Compute data log-likelihood for the gDNA pseudo-component.

    Returns log_lik (float).
    """
    splice_pen = ctx.gdna_splice_penalties.get(splice_type, 1.0)

    if exon_strand == 1 or exon_strand == 2:
        p_ig = ctx.ig_p if exon_strand == 1 else (1.0 - ctx.ig_p)
    else:
        p_ig = 0.5

    log_strand = math.log(max(p_ig, LOG_SAFE_FLOOR))
    log_insert = frag_len_log_lik(ctx, genomic_footprint)

    return (
        log_strand
        + log_insert
        + math.log(max(splice_pen, LOG_SAFE_FLOOR))
    )


# ---------------------------------------------------------------------------
# Standalone scoring (for unit tests / external callers)
# ---------------------------------------------------------------------------
# These match the old _score_gdna_candidate / _score_candidate signatures
# so that tests can call them without a ScoringContext.


def score_gdna_standalone(
    exon_strand: int,
    splice_type: int,
    frag_length: int,
    strand_models,
    frag_length_models,
    gdna_splice_penalties: dict | None = None,
) -> float:
    """Compute gDNA log-likelihood using model objects directly.

    Drop-in replacement for the old ``_score_gdna_candidate`` function.
    """
    penalties = gdna_splice_penalties or GDNA_SPLICE_PENALTIES
    splice_pen = penalties.get(splice_type, 1.0)

    p_strand = strand_models.intergenic.strand_likelihood_int(
        exon_strand, STRAND_POS,
    )

    log_p_insert = (
        frag_length_models.global_model.log_likelihood(frag_length)
        if frag_length > 0
        else 0.0
    )

    return (
        math.log(max(p_strand, LOG_SAFE_FLOOR))
        + log_p_insert
        + math.log(max(splice_pen, LOG_SAFE_FLOOR))
    )
