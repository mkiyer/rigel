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
            t_length_arr=t_length_arr,
            t_span_arr=t_span_arr,
            t_start_arr=t_start_arr,
            _t_exon_data=t_exon_data,
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
    t_strand: int,
    frag_len_ll: float,
    overhang_bp: int,
    nm: int,
) -> float:
    """Compute data log-likelihood for a nascent RNA (nRNA) candidate.

    nRNA is a transcript-level phenomenon: the strand comparison uses
    the transcript strand (not gene strand).

    Returns log_lik (float).
    """
    log_same, log_diff = strand_log_probs(ctx, exon_strand)
    log_nm = nm * ctx.mismatch_log_penalty if nm > 0 else 0.0

    if exon_strand == t_strand:
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


def compute_coverage_weight(
    frag_start: int,
    frag_end: int,
    transcript_length: int,
) -> float:
    """Inverse coverage-capacity weight for a fragment on a transcript.

    Uses the trapezoid coverage model: under uniform random sampling,
    the expected coverage at transcript-relative position *x* on a
    transcript of spliced length *L* for a fragment of length *f* is
    ``min(x, w, L-x)`` where ``w = min(f, L/2)``.

    A fragment in the plateau region gets weight 1.0.  A fragment near
    a transcript edge gets weight > 1.0, reflecting that it is less
    likely to be generated yet was still observed — stronger evidence
    for that transcript.

    Parameters
    ----------
    frag_start : int
        Fragment 5' position in transcript-relative coordinates ``[0, L)``.
    frag_end : int
        Fragment 3' position in transcript-relative coordinates ``(0, L]``.
    transcript_length : int
        Spliced exonic length of the transcript.

    Returns
    -------
    float
        Weight >= 1.0.  Plateau → 1.0; edge → > 1.0.
    """
    f_len = frag_end - frag_start
    if f_len <= 0 or transcript_length <= 0:
        return 1.0

    L = float(transcript_length)
    w = min(f_len, L / 2.0)

    total_area = 0.0

    # Left ramp: c(x) = x for x in [0, w)
    s = max(0.0, min(float(frag_start), w))
    e = max(0.0, min(float(frag_end), w))
    if s < e:
        total_area += (e - s) * (s + e) / 2.0

    # Plateau: c(x) = w for x in [w, L - w)
    s = max(w, min(float(frag_start), L - w))
    e = max(w, min(float(frag_end), L - w))
    if s < e:
        total_area += (e - s) * w

    # Right ramp: c(x) = L - x for x in [L - w, L]
    s = max(L - w, min(float(frag_start), L))
    e = max(L - w, min(float(frag_end), L))
    if s < e:
        total_area += (e - s) * ((L - s) + (L - e)) / 2.0

    area_per_base = total_area / f_len
    area_per_base = max(area_per_base, 1.0)

    return w / area_per_base


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


# ---------------------------------------------------------------------------
# Vectorized batch functions (Fix A, B of P2 optimization)
# ---------------------------------------------------------------------------


def compute_coverage_weight_batch(
    frag_starts: np.ndarray,
    frag_ends: np.ndarray,
    transcript_lengths: np.ndarray,
) -> np.ndarray:
    """Vectorized coverage weight for arrays of candidates.

    Equivalent to calling ``compute_coverage_weight`` per element,
    but 10-50× faster by eliminating per-call Python overhead and
    using NumPy vectorized arithmetic.

    Parameters
    ----------
    frag_starts, frag_ends : np.ndarray
        Fragment positions in transcript-relative coordinates.
    transcript_lengths : np.ndarray
        Spliced exonic lengths of the transcripts.

    Returns
    -------
    np.ndarray
        Float64 array of weights (>= 1.0).
    """
    fs = np.asarray(frag_starts, dtype=np.float64)
    fe = np.asarray(frag_ends, dtype=np.float64)
    L = np.asarray(transcript_lengths, dtype=np.float64)

    f_len = fe - fs
    n = len(fs)
    result = np.ones(n, dtype=np.float64)
    valid = (f_len > 0) & (L > 0)
    if not np.any(valid):
        return result

    # Work on valid subset only
    fsv = fs[valid]
    fev = fe[valid]
    flv = f_len[valid]
    Lv = L[valid]
    w = np.minimum(flv, Lv / 2.0)

    total_area = np.zeros(int(np.sum(valid)), dtype=np.float64)

    # Left ramp: c(x) = x for x in [0, w)
    s = np.clip(fsv, 0.0, w)
    e = np.clip(fev, 0.0, w)
    mask = s < e
    total_area += np.where(mask, (e - s) * (s + e) / 2.0, 0.0)

    # Plateau: c(x) = w for x in [w, L-w)
    Lmw = Lv - w
    s = np.clip(fsv, w, Lmw)
    e = np.clip(fev, w, Lmw)
    mask = s < e
    total_area += np.where(mask, (e - s) * w, 0.0)

    # Right ramp: c(x) = L-x for x in [L-w, L]
    s = np.clip(fsv, Lmw, Lv)
    e = np.clip(fev, Lmw, Lv)
    mask = s < e
    total_area += np.where(mask, (e - s) * ((Lv - s) + (Lv - e)) / 2.0, 0.0)

    area_per_base = total_area / flv
    area_per_base = np.maximum(area_per_base, 1.0)
    result[valid] = w / area_per_base

    return result


def frag_len_log_lik_batch(
    ctx: "ScoringContext",
    frag_lengths: np.ndarray,
) -> np.ndarray:
    """Vectorized fragment-length log-likelihood via LUT lookup.

    Equivalent to calling ``frag_len_log_lik`` per element.

    Parameters
    ----------
    ctx : ScoringContext
        Pre-computed context with fragment-length LUT.
    frag_lengths : np.ndarray
        Int32 array of fragment lengths (negative values → 0.0).

    Returns
    -------
    np.ndarray
        Float64 array of log-likelihoods.
    """
    result = np.zeros(len(frag_lengths), dtype=np.float64)
    if ctx.fl_log_prob is None:
        return result
    pos = frag_lengths > 0
    in_range = pos & (frag_lengths <= ctx.fl_max_size)
    if np.any(in_range):
        result[in_range] = ctx.fl_log_prob[frag_lengths[in_range]]
    tail = pos & (frag_lengths > ctx.fl_max_size)
    if np.any(tail):
        result[tail] = (
            ctx.fl_tail_base
            + (frag_lengths[tail].astype(np.float64) - ctx.fl_max_size)
            * _TAIL_DECAY_LP
        )
    return result


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
