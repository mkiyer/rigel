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

import bisect
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

    # Per-transcript lengths for per-fragment effective length correction
    t_length_arr: np.ndarray     # int32[n_transcripts] — spliced exonic length
    t_span_arr: np.ndarray       # int32[n_transcripts] — genomic span (incl introns)
    t_start_arr: np.ndarray      # int32[n_transcripts] — genomic start coordinate

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

        # Per-transcript length arrays for per-fragment effective length
        t_length_arr = index.t_df["length"].values.astype(np.int32)
        t_span_arr = (
            index.t_df["end"].values - index.t_df["start"].values
        ).astype(np.int32)
        t_start_arr = index.t_df["start"].values.astype(np.int32)

        # Pre-compute per-transcript exon positions as Python tuples
        # for fast bisect-based genomic → transcript coordinate mapping.
        t_exon_data: dict = {}
        for t_idx in range(index.num_transcripts):
            exon_ivs = index.get_exon_intervals(t_idx)
            if exon_ivs is not None and len(exon_ivs) > 0:
                starts = tuple(int(x) for x in exon_ivs[:, 0])
                ends = tuple(int(x) for x in exon_ivs[:, 1])
                cs = 0
                cb: list[int] = []
                for s, e in zip(starts, ends):
                    cb.append(cs)
                    cs += e - s
                t_exon_data[t_idx] = (starts, ends, tuple(cb))

        ctx = ScoringContext(
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
            t_length_arr=t_length_arr,
            t_span_arr=t_span_arr,
            t_start_arr=t_start_arr,
            _t_exon_data=t_exon_data,
        )

        # Build native C++ scoring context for hot-path acceleration.
        # Cast arrays to exact dtypes expected by the C++ nanobind binding
        # to tolerate callers (including test mocks) that provide int64 etc.
        from hulkrna._scoring_impl import NativeScoringContext
        native_ctx = NativeScoringContext(
            log_p_sense=float(ctx.log_p_sense),
            log_p_antisense=float(ctx.log_p_antisense),
            anti_flag=bool(ctx.anti_flag),
            overhang_log_penalty=float(ctx.overhang_log_penalty),
            mismatch_log_penalty=float(ctx.mismatch_log_penalty),
            fl_log_prob=ctx.fl_log_prob,
            fl_max_size=int(ctx.fl_max_size),
            fl_tail_base=float(ctx.fl_tail_base),
            t_strand_arr=np.ascontiguousarray(
                ctx.t_strand_arr, dtype=np.int8),
            t_length_arr=np.ascontiguousarray(
                ctx.t_length_arr, dtype=np.int32),
            t_span_arr=np.ascontiguousarray(
                ctx.t_span_arr, dtype=np.int32),
            t_start_arr=np.ascontiguousarray(
                ctx.t_start_arr, dtype=np.int32),
            nrna_base=int(ctx.nrna_base),
            t_exon_data=ctx._t_exon_data,
        )
        object.__setattr__(ctx, '_native_ctx', native_ctx)

        return ctx


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


# ---------------------------------------------------------------------------
# Standalone scoring (for unit tests / external callers)
# ---------------------------------------------------------------------------


def genomic_to_transcript_pos(
    genomic_pos: int,
    exon_intervals: np.ndarray,
    strand: int,
    transcript_length: int,
) -> int:
    """Map a genomic position to a transcript-relative 5'→3' coordinate.

    Walks through sorted exon intervals ``(start, end)`` accumulating
    spliced offset until the exon containing *genomic_pos* is found.
    Positions falling in introns are mapped to the preceding exon
    boundary.  Negative-strand transcripts are flipped so that
    offset 0 is the 5' end.

    Parameters
    ----------
    genomic_pos : int
        Genomic coordinate to map.
    exon_intervals : np.ndarray
        ``(n_exons, 2)`` int32 array of ``[start, end)`` intervals,
        sorted by genomic start position.
    strand : int
        Transcript strand (``Strand.POS`` or ``Strand.NEG``).
    transcript_length : int
        Total spliced exonic length of the transcript.

    Returns
    -------
    int
        Transcript-relative position in ``[0, transcript_length]``.
    """
    offset = 0
    n = len(exon_intervals)
    for i in range(n):
        ex_start = int(exon_intervals[i, 0])
        ex_end = int(exon_intervals[i, 1])
        if genomic_pos < ex_start:
            # Position is upstream of (or in an intron before) this exon
            break
        if genomic_pos < ex_end:
            # Position is inside this exon
            offset += genomic_pos - ex_start
            break
        offset += ex_end - ex_start
    else:
        # Position is past the last exon
        offset = transcript_length

    offset = max(0, min(offset, transcript_length))
    if strand == STRAND_NEG:
        offset = transcript_length - offset
    return offset


# ---------------------------------------------------------------------------
# compute_fragment_weight: C++ kernel
# ---------------------------------------------------------------------------

from hulkrna._scoring_impl import compute_fragment_weight


def score_gdna_standalone(
    exon_strand: int,
    splice_type: int,
    frag_length: int,
    strand_models,
    frag_length_models,
    gdna_splice_penalties: dict | None = None,
) -> float:
    """Compute gDNA log-likelihood using model objects directly.

    Intended for unit tests and external callers that do not construct
    a ``ScoringContext``.
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


def genomic_to_transcript_pos_bisect(
    genomic_pos: int,
    exon_starts: tuple[int, ...],
    exon_ends: tuple[int, ...],
    cumsum_before: tuple[int, ...],
    strand: int,
    transcript_length: int,
) -> int:
    """Map genomic position to transcript coordinate using binary search.

    Like ``genomic_to_transcript_pos`` but uses pre-computed cumulative
    sums and ``bisect`` for O(log n) exon lookup instead of O(n) linear
    scan with numpy int() conversions.

    Parameters
    ----------
    genomic_pos : int
        Genomic coordinate to map.
    exon_starts : tuple[int, ...]
        Sorted exon start coordinates (Python ints).
    exon_ends : tuple[int, ...]
        Sorted exon end coordinates (Python ints).
    cumsum_before : tuple[int, ...]
        Cumulative spliced length before each exon.
    strand : int
        Transcript strand (``Strand.POS`` or ``Strand.NEG``).
    transcript_length : int
        Total spliced exonic length.

    Returns
    -------
    int
        Transcript-relative position in ``[0, transcript_length]``.
    """
    # Find rightmost exon whose start <= genomic_pos
    ei = bisect.bisect_right(exon_starts, genomic_pos) - 1
    if ei < 0:
        # Before first exon
        offset = 0
    elif genomic_pos >= exon_ends[ei]:
        # In intron after exon ei (or past last exon)
        offset = cumsum_before[ei] + (exon_ends[ei] - exon_starts[ei])
    else:
        # Inside exon ei
        offset = cumsum_before[ei] + (genomic_pos - exon_starts[ei])
    if offset < 0:
        offset = 0
    elif offset > transcript_length:
        offset = transcript_length
    if strand == STRAND_NEG:
        offset = transcript_length - offset
    return offset
