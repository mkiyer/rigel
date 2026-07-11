"""rigel.scoring — scoring parameter assembly for the native scorer.

The fragment-scoring math itself lives in the C++ extension
(``src/rigel/native/scoring.cpp``).  This module only holds
``overhang_alpha_to_log_penalty()`` and the ``FragmentScorer`` builder,
which gathers the trained strand / fragment-length models and index
geometry into a single parameter bundle and hands them to a
``NativeFragmentScorer``.
"""

import math
from dataclasses import dataclass

import numpy as np

from .splice import SPLICE_UNSPLICED

# ---------------------------------------------------------------------------
# Constants (single source of truth for all scoring/penalty values)
# ---------------------------------------------------------------------------

#: Floor value for log-safe clamping to avoid log(0).
LOG_SAFE_FLOOR = 1e-10

#: Pre-computed log(0.5) — used for uninformative strand log-probabilities.
LOG_HALF = math.log(0.5)

#: Default overhang alpha: each base of overhang reduces probability by 10×.
DEFAULT_OVERHANG_ALPHA = 0.1

#: Default mismatch alpha: each edit-distance mismatch (NM tag) reduces
#: probability by 10×.
DEFAULT_MISMATCH_ALPHA = 0.1


#: Default gDNA splice penalties per SpliceType (int keys for fast lookup). Only the unspliced entry is
#: consumed — spliced fragments are structurally incompatible with the (genomic) gDNA component, so they
#: never reach a gDNA splice penalty.
GDNA_SPLICE_PENALTIES = {
    SPLICE_UNSPLICED: 1.0,
}


def overhang_alpha_to_log_penalty(alpha: float) -> float:
    """Convert user-facing alpha to internal log-penalty.

    Parameters
    ----------
    alpha : float
        Per-base overhang penalty in [0, 1].
        - 0.0 → hard binary gate (any overhang → −∞)
        - 0.01 → aggressive (−4.605 per base)
        - 0.1 → default (−2.303 per base)
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
DEFAULT_OVERHANG_LOG_PENALTY = overhang_alpha_to_log_penalty(DEFAULT_OVERHANG_ALPHA)

#: Default mismatch log-penalty.
DEFAULT_MISMATCH_LOG_PENALTY = overhang_alpha_to_log_penalty(DEFAULT_MISMATCH_ALPHA)


# ---------------------------------------------------------------------------
# FragmentScorer — parameter bundle for the native scorer
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FragmentScorer:
    """Pre-computed scoring parameters, built once per pipeline run.

    Bundles the trained strand / fragment-length models, penalty
    parameters, and index geometry.  ``from_models`` assembles these into
    a ``NativeFragmentScorer`` (the C++ scorer in ``scoring.cpp``) that
    does the actual per-fragment scoring.
    """

    # Strand model: RNA
    log_p_sense: float
    log_p_antisense: float
    r1_antisense: bool  # True when R1-antisense protocol (p_sense < 0.5)

    # Penalty parameters
    overhang_log_penalty: float
    mismatch_log_penalty: float
    gdna_splice_penalties: dict  # int-keyed

    # Fragment-length LUT — RNA model (pre-finalized)
    fl_log_prob: np.ndarray | None  # numpy array or None
    fl_max_size: int
    fl_tail_base: float

    # Fragment-length LUT — gDNA model (pre-finalized)
    gdna_fl_log_prob: np.ndarray | None
    gdna_fl_max_size: int
    gdna_fl_tail_base: float

    # Index arrays (borrowed references, never copied)
    t_strand_arr: np.ndarray  # int8[n_transcripts]
    t_to_g: np.ndarray  # int32[n_transcripts]

    # Transcript geometry used by RNA scoring.
    t_length_arr: np.ndarray  # int32[n_transcripts] — spliced exonic length

    @staticmethod
    def from_models(
        strand_models,
        rna_fl,
        gdna_fl,
        index,
        *,
        overhang_log_penalty: float | None = None,
        mismatch_log_penalty: float | None = None,
        gdna_splice_penalties: dict | None = None,
        pruning_min_posterior: float = 1e-4,
    ) -> "FragmentScorer":
        """Build a FragmentScorer from trained models and index.

        Parameters
        ----------
        strand_models : StrandModels
        rna_fl : FragmentLengthModel
            Finalised RNA fragment-length scoring model.
        gdna_fl : FragmentLengthModel
            Finalised gDNA fragment-length scoring model.
        index : TranscriptIndex
        overhang_log_penalty : float or None
        mismatch_log_penalty : float or None
        gdna_splice_penalties : dict or None
        pruning_min_posterior : float
            Minimum posterior threshold for candidate pruning.
            Lower values are more conservative (keep more candidates).
        """
        rna_sm = strand_models.exonic_spliced
        p_sense = rna_sm._cached_p_sense
        p_antisense = rna_sm._cached_p_antisense

        fl_log_prob = rna_fl._log_prob  # numpy array or None
        fl_max_size = rna_fl.max_size
        fl_tail_base: float = getattr(rna_fl, "_tail_base", 0.0)

        gdna_fl_log_prob = gdna_fl._log_prob
        gdna_fl_max_size = gdna_fl.max_size
        gdna_fl_tail_base: float = getattr(gdna_fl, "_tail_base", 0.0)
        # Per-transcript geometry arrays for RNA scoring.
        t_length_arr = index.t_df["length"].values.astype(np.int32)

        # Build per-transcript exon CSR arrays for coverage-weight
        # and genomic→transcript coordinate mapping. This replaces
        # a 457K-iteration Python loop with a single vectorized call.
        exon_offsets, exon_starts, exon_ends, exon_cumsum = index.build_exon_csr()

        ctx = FragmentScorer(
            log_p_sense=math.log(max(p_sense, LOG_SAFE_FLOOR)),
            log_p_antisense=math.log(max(p_antisense, LOG_SAFE_FLOOR)),
            r1_antisense=p_sense < 0.5,
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
            gdna_splice_penalties=(gdna_splice_penalties or GDNA_SPLICE_PENALTIES),
            fl_log_prob=fl_log_prob,
            fl_max_size=fl_max_size,
            fl_tail_base=fl_tail_base,
            gdna_fl_log_prob=gdna_fl_log_prob,
            gdna_fl_max_size=gdna_fl_max_size,
            gdna_fl_tail_base=gdna_fl_tail_base,
            t_strand_arr=index.t_to_strand_arr,
            t_to_g=index.t_to_g_arr,
            t_length_arr=t_length_arr,
        )

        # Build native C++ scoring context for hot-path acceleration.
        # Cast arrays to exact dtypes expected by the C++ nanobind binding
        # to tolerate callers (including test mocks) that provide int64 etc.
        from .native import NativeFragmentScorer

        # Pool-separated likelihood pruning: compute Δ = -log(ε)
        _eps = max(pruning_min_posterior, 1e-300)
        max_ll_delta = -math.log(_eps) if _eps < 1.0 else 0.0

        # Build is_nrna array for pool separation.
        # Only synthetic (rigel-generated) nRNA spans belong to the nRNA pool.
        # Annotated single-exon transcripts are valid mRNA and must be scored
        # as mRNA even though they can serve as nRNA components.
        t_is_nrna_arr = index.t_df["is_synthetic"].values.astype(np.uint8)

        native_ctx = NativeFragmentScorer(
            log_p_sense=float(ctx.log_p_sense),
            log_p_antisense=float(ctx.log_p_antisense),
            r1_antisense=bool(ctx.r1_antisense),
            overhang_log_penalty=float(ctx.overhang_log_penalty),
            mismatch_log_penalty=float(ctx.mismatch_log_penalty),
            fl_log_prob=ctx.fl_log_prob,
            fl_max_size=int(ctx.fl_max_size),
            fl_tail_base=float(ctx.fl_tail_base),
            gdna_fl_log_prob=ctx.gdna_fl_log_prob,
            gdna_fl_max_size=int(ctx.gdna_fl_max_size),
            gdna_fl_tail_base=float(ctx.gdna_fl_tail_base),
            t_strand_arr=np.ascontiguousarray(ctx.t_strand_arr, dtype=np.int8),
            t_length_arr=np.ascontiguousarray(ctx.t_length_arr, dtype=np.int32),
            exon_offsets=np.ascontiguousarray(exon_offsets, dtype=np.int32),
            exon_starts=np.ascontiguousarray(exon_starts, dtype=np.int32),
            exon_ends=np.ascontiguousarray(exon_ends, dtype=np.int32),
            exon_cumsum=np.ascontiguousarray(exon_cumsum, dtype=np.int32),
            t_is_nrna_arr=np.ascontiguousarray(t_is_nrna_arr),
            pruning_max_ll_delta=max_ll_delta,
        )
        object.__setattr__(ctx, "_native_ctx", native_ctx)

        return ctx
