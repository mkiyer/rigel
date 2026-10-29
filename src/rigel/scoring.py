"""rigel.scoring — the native fragment scorer, built from one run's models.

The fragment-scoring math lives in the C++ extension (``src/rigel/native/scoring.cpp``). :func:`build_scorer` gathers
what it reads — the trained strand and fragment-length models, the penalties and the index's transcript geometry —
into one ``NativeFragmentScorer``, which ``rigel.routing.FragmentRouter`` drives over the buffer. A hit's sequence is
priced at the library's own mismatch rate (:func:`sequence_log_penalties`).
"""

import math

import numpy as np

from .config import CONSTANTS, FragmentScoringConfig


#: A base is one of four.
_N_BASES = 4


def sequence_log_penalties(mismatch_rate: float | None) -> tuple[float, float]:
    """A hit's two sequence prices per base, each against an aligned base that matches: ``(mismatch, clip)``.

    At the library's mismatch rate ε an aligned base reads its reference base with probability 1 − ε and each of the
    three others with ε/3, while a soft-clipped base is placed nowhere and reads at the composition of the four bases,
    ¼. So a mismatch costs log[ε / (3(1 − ε))] and a clipped base log[1 / (4(1 − ε))]
    (``tests/test_hit_sequence.py`` enumerates the emissions). A library with no mismatch makes one impossible (−inf);
    with no unique unspliced base to read ε on, a hit's sequence is not priced.
    """
    if mismatch_rate is None:
        return 0.0, 0.0
    eps = float(mismatch_rate)
    if not 0.0 <= eps < (_N_BASES - 1) / _N_BASES:
        raise ValueError(
            f"a mismatch rate of {eps:g} is no library's: at or past {(_N_BASES - 1) / _N_BASES:g} a mismatch is "
            "at least as likely as a match"
        )
    mismatch = math.log(eps / ((_N_BASES - 1) * (1.0 - eps))) if eps > 0.0 else -math.inf
    return mismatch, -math.log(_N_BASES * (1.0 - eps))


def build_scorer(
    strand_models,
    rna_fl,
    gdna_fl,
    index,
    *,
    overhang_log_penalty: float = FragmentScoringConfig().overhang_log_penalty,
    mismatch_log_penalty: float = 0.0,
    clip_log_penalty: float = 0.0,
    pruning_min_posterior: float = FragmentScoringConfig.pruning_min_posterior,
):
    """The run's ``NativeFragmentScorer``.

    Parameters
    ----------
    strand_models : StrandModels
    rna_fl, gdna_fl : FragmentLengthModel
        The finalised RNA and gDNA fragment-length scoring models.
    index : TranscriptIndex
    overhang_log_penalty : float
        Log-penalty per base of overhang (``FragmentScoringConfig``'s).
    mismatch_log_penalty, clip_log_penalty : float
        A multimapper's hit's sequence prices, per NM mismatch and per soft-clipped base
        (:func:`sequence_log_penalties`); 0, the default, prices no sequence.
    pruning_min_posterior : float
        A candidate whose likelihood falls below this fraction of its pool's best is pruned.
    """
    from .native import NativeFragmentScorer

    p_sense = strand_models.exonic_spliced.p_r1_sense
    p_antisense = strand_models.exonic_spliced.p_r1_antisense
    floor = CONSTANTS.scoring.strand_probability_floor
    # Pool-separated likelihood pruning: Δ = -log(threshold)
    threshold = max(pruning_min_posterior, CONSTANTS.scoring.pruning_floor)
    # The coverage weight's exon CSR, built once by the index and cached there
    exon_offsets, exon_starts, exon_ends, exon_cumsum = index.build_exon_csr()
    t_df = index.t_df
    return NativeFragmentScorer(
        log_p_sense=math.log(max(p_sense, floor)),
        log_p_antisense=math.log(max(p_antisense, floor)),
        r1_antisense=bool(p_sense < 0.5),
        overhang_log_penalty=float(overhang_log_penalty),
        mismatch_log_penalty=float(mismatch_log_penalty),
        clip_log_penalty=float(clip_log_penalty),
        fl_log_prob=rna_fl._log_prob,
        fl_max_size=int(rna_fl.max_size),
        fl_tail_base=float(getattr(rna_fl, "_tail_base", 0.0)),
        gdna_fl_log_prob=gdna_fl._log_prob,
        gdna_fl_max_size=int(gdna_fl.max_size),
        gdna_fl_tail_base=float(getattr(gdna_fl, "_tail_base", 0.0)),
        t_strand_arr=np.ascontiguousarray(index.t_to_strand_arr, dtype=np.int8),
        t_length_arr=np.ascontiguousarray(t_df["length"].values, dtype=np.int32),
        exon_offsets=np.ascontiguousarray(exon_offsets, dtype=np.int32),
        exon_starts=np.ascontiguousarray(exon_starts, dtype=np.int32),
        exon_ends=np.ascontiguousarray(exon_ends, dtype=np.int32),
        exon_cumsum=np.ascontiguousarray(exon_cumsum, dtype=np.int32),
        # Only rigel's synthetic spans form the span pool; an annotated single-exon transcript is pruned with the
        # annotated ones, even one covering a multi-exon transcript's span
        t_is_synthetic_arr=np.ascontiguousarray(t_df["is_synthetic"].values, dtype=np.uint8),
        # Each transcript's reference, as a code: a multimapper's genic places are told apart on it
        t_ref_arr=np.ascontiguousarray(
            np.unique(t_df["ref"].astype(str).to_numpy(), return_inverse=True)[1], dtype=np.int32
        ),
        pruning_max_ll_delta=-math.log(threshold) if threshold < 1.0 else 0.0,
    )
