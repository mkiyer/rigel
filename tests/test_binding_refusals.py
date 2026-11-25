"""Inputs the native bindings refuse once, at the boundary, so the loops behind them read unguarded: the scorer's
count matrix (a deterministic count of the wrong shape was skipped silently) and the annotation writer's arrays (their
length was a separate argument, never checked against them)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.native import (
    BamAnnotationWriter,
    FragmentResolver,
    NativeFragmentScorer,
    StreamingScorer,
)


def _scorer(n_transcripts: int = 2) -> NativeFragmentScorer:
    one = np.ones(n_transcripts, np.int32)
    return NativeFragmentScorer(
        log_p_sense=0.0,
        log_p_antisense=0.0,
        r1_antisense=False,
        overhang_log_penalty=0.0,
        mismatch_log_penalty=0.0,
        clip_log_penalty=0.0,
        fl_log_prob=None,
        gdna_fl_log_prob=None,
        t_strand_arr=np.ones(n_transcripts, np.int8),
        t_length_arr=100 * one,
        exon_offsets=np.arange(n_transcripts + 1, dtype=np.int32),
        exon_starts=np.zeros(n_transcripts, np.int32),
        exon_ends=100 * one,
        exon_cumsum=np.zeros(n_transcripts, np.int32),
        t_is_synthetic_arr=np.zeros(n_transcripts, np.uint8),
        t_ref_arr=np.zeros(n_transcripts, np.int32),
        pruning_max_ll_delta=10.0,
    )


def test_the_scorer_takes_a_count_matrix_of_one_row_per_transcript_and_six_columns():
    StreamingScorer(_scorer(), np.zeros((2, 6)))


@pytest.mark.parametrize("shape", [(2, 5), (3, 6), (1, 6)])
def test_the_scorer_refuses_a_count_matrix_of_another_shape(shape):
    with pytest.raises(ValueError, match="unambig_counts is not"):
        StreamingScorer(_scorer(), np.zeros(shape))


def test_the_writer_refuses_annotation_arrays_of_unequal_length(tmp_path):
    writer = BamAnnotationWriter(FragmentResolver(), "auto", True, True, 1000, 0)
    n = 3
    arrays = [
        np.zeros(n, np.int64),
        np.zeros(n, np.int32),
        np.zeros(n, np.int32),
        np.zeros(n, np.uint8),
        np.zeros(n - 1, np.float32),
        np.zeros(n, np.int8),
        np.zeros(n, np.int32),
        np.zeros(n, np.int32),
    ]
    with pytest.raises(ValueError, match="not one length"):
        writer.write(str(tmp_path / "in.bam"), str(tmp_path / "out.bam"), *arrays, [], [], [])
