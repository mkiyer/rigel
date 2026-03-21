"""rigel.scoring — Pure scoring functions and pre-computed scoring context.

Every function in this module takes explicit arguments (no closures,
no captured state).  The ``FragmentScorer`` dataclass holds all
pre-computed parameters that were previously captured as locals inside
``_scan_and_build_em_data``'s nested closures.

Design contract:
    - All functions are stateless pure functions.
    - ``FragmentScorer`` maps directly to a C struct.
    - Every function is a candidate for Cython / C acceleration.
"""

import math
from dataclasses import dataclass

import numpy as np

from .splice import SPLICE_UNSPLICED, SPLICE_UNANNOT, SPLICE_ANNOT

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
DEFAULT_OVERHANG_LOG_PENALTY = overhang_alpha_to_log_penalty(DEFAULT_OVERHANG_ALPHA)

#: Default mismatch log-penalty.
DEFAULT_MISMATCH_LOG_PENALTY = overhang_alpha_to_log_penalty(DEFAULT_MISMATCH_ALPHA)


# ---------------------------------------------------------------------------
# FragmentScorer — replaces 15+ closure-captured locals
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FragmentScorer:
    """Pre-computed scoring parameters, built once per pipeline run.

    Replaces the 15+ captured locals in the old closure-based design.
    Every scoring function receives this as its first argument.
    Maps directly to a C struct for future acceleration.
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
    nrna_base: int  # offset for nRNA indices in CSR

    # Per-transcript lengths for per-fragment effective length correction
    t_length_arr: np.ndarray  # int32[n_transcripts] — spliced exonic length
    t_span_arr: np.ndarray  # int32[n_transcripts] — genomic span (incl introns)
    t_start_arr: np.ndarray  # int32[n_transcripts] — genomic start coordinate

    # Pre-computed exon positions for fast bisect-based position mapping.
    # Dict mapping t_idx → (exon_starts: tuple[int], exon_ends: tuple[int],
    #                        cumsum_before: tuple[int]).
    # Python tuples of ints for zero-overhead bisect access.
    _t_exon_data: dict | None = None

    @staticmethod
    def from_models(
        strand_models,
        frag_length_models,
        index,
        estimator,
        *,
        overhang_log_penalty: float | None = None,
        mismatch_log_penalty: float | None = None,
        gdna_splice_penalties: dict | None = None,
    ) -> "FragmentScorer":
        """Build a FragmentScorer from trained models and index.

        Parameters
        ----------
        strand_models : StrandModels
        frag_length_models : FragmentLengthModels
        index : TranscriptIndex
        estimator : AbundanceEstimator
        overhang_log_penalty : float or None
        mismatch_log_penalty : float or None
        gdna_splice_penalties : dict or None
        """
        rna_sm = strand_models.exonic_spliced
        p_sense = rna_sm._cached_p_sense
        p_antisense = rna_sm._cached_p_antisense

        rna_fl = frag_length_models.rna_model
        fl_log_prob = rna_fl._log_prob  # numpy array or None
        fl_max_size = rna_fl.max_size
        fl_tail_base: float = getattr(rna_fl, "_tail_base", 0.0)

        gdna_fl = frag_length_models.gdna_model
        gdna_fl_log_prob = gdna_fl._log_prob
        gdna_fl_max_size = gdna_fl.max_size
        gdna_fl_tail_base: float = getattr(gdna_fl, "_tail_base", 0.0)

        # Per-transcript length arrays for per-fragment effective length
        t_length_arr = index.t_df["length"].values.astype(np.int32)
        t_span_arr = (index.t_df["end"].values - index.t_df["start"].values).astype(np.int32)
        t_start_arr = index.t_df["start"].values.astype(np.int32)

        # Pre-compute per-transcript exon positions as Python tuples
        # for fast bisect-based genomic → transcript coordinate mapping.
        t_exon_data: dict = {}
        for t_idx in range(index.num_transcripts):
            exon_ivs = index.get_exon_intervals(t_idx)
            if exon_ivs is not None and len(exon_ivs) > 0:
                starts = tuple(exon_ivs[:, 0].tolist())
                ends = tuple(exon_ivs[:, 1].tolist())
                lengths = exon_ivs[:, 1] - exon_ivs[:, 0]
                cumsum_before = tuple(np.concatenate(([0], np.cumsum(lengths[:-1]))).tolist())
                t_exon_data[t_idx] = (starts, ends, cumsum_before)

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
            nrna_base=estimator.nrna_base_index,
            t_length_arr=t_length_arr,
            t_span_arr=t_span_arr,
            t_start_arr=t_start_arr,
            _t_exon_data=t_exon_data,
        )

        # Build native C++ scoring context for hot-path acceleration.
        # Cast arrays to exact dtypes expected by the C++ nanobind binding
        # to tolerate callers (including test mocks) that provide int64 etc.
        from .native import NativeFragmentScorer

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
            t_span_arr=np.ascontiguousarray(ctx.t_span_arr, dtype=np.int32),
            t_start_arr=np.ascontiguousarray(ctx.t_start_arr, dtype=np.int32),
            nrna_base=int(ctx.nrna_base),
            t_exon_data=ctx._t_exon_data,
            t_to_nrna_arr=np.ascontiguousarray(index.t_to_nrna_arr, dtype=np.int32),
            nrna_span_arr=np.ascontiguousarray(
                (index.nrna_df["end"].values - index.nrna_df["start"].values).astype(np.int32)
            ),
            nrna_start_arr=np.ascontiguousarray(index.nrna_df["start"].values.astype(np.int32)),
        )
        object.__setattr__(ctx, "_native_ctx", native_ctx)

        return ctx
