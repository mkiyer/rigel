"""`quant_accuracy.py`'s TPM rows carry no transcript that has no effective length, on either side.

A transcript with no effective length has no rate, so it can carry no TPM. A shadow transcript is simulated but not
annotated: it is in the truth table and absent from the tool's, so its merged length is 0. Floored to 1, its whole
fragment count became a rate, which dominated the truth TPM's denominator and shrank every other row's ``tpm_true``.
Its fragments are not hidden: they are reported in a column of their own.
"""

from __future__ import annotations

import pandas as pd
import pytest

import quant_accuracy as qa


def _annotated():
    """Three annotated transcripts the tool scores perfectly: its TPM is the truth's own."""
    truth = pd.DataFrame(
        {"transcript_id": ["a", "b", "c"], "observed_mrna_fragments": [100.0, 50.0, 0.0]}
    )
    eff = [1000.0, 500.0, 800.0]
    rate = [100.0 / 1000.0, 50.0 / 500.0, 0.0]
    tpm = [r / sum(rate) * 1e6 for r in rate]
    quant = pd.DataFrame(
        {
            "transcript_id": ["a", "b", "c"],
            "count": [100.0, 50.0, 0.0],
            "count_em": [0.0, 0.0, 0.0],
            "tpm": tpm,
            "effective_length": eff,
        }
    )
    return truth, quant


def test_a_perfect_estimate_scores_zero_tpm_error():
    truth, quant = _annotated()
    assert qa.score_transcripts(quant, truth)["tpm_abs_err"] == pytest.approx(0.0, abs=1e-6)


def test_a_shadow_transcript_does_not_move_the_annotated_rows_tpm():
    truth, quant = _annotated()
    shadow = pd.concat(
        [truth, pd.DataFrame({"transcript_id": ["shadow"], "observed_mrna_fragments": [50.0]})],
        ignore_index=True,
    )
    row = qa.score_transcripts(quant, shadow)
    assert row["tpm_abs_err"] == pytest.approx(0.0, abs=1e-6), row["tpm_abs_err"]
    assert row["no_length_fragments"] == 50.0
    assert row["n_no_length"] == 1


def test_a_tool_row_with_zero_length_is_out_of_both_vectors():
    """A transcript the tool holds at length 0, with fragments on both sides: out of the estimate's TPM vector and its
    denominator as well as the truth's, so the other rows still compare exactly."""
    truth, quant = _annotated()
    truth = pd.concat(
        [truth, pd.DataFrame({"transcript_id": ["z"], "observed_mrna_fragments": [30.0]})],
        ignore_index=True,
    )
    rate = [100.0 / 1000.0, 50.0 / 500.0, 0.0, 30.0]  # the tool floors a zero length at 1
    tpm = [r / sum(rate) * 1e6 for r in rate]
    quant = pd.concat(
        [
            quant,
            pd.DataFrame(
                {
                    "transcript_id": ["z"],
                    "count": [30.0],
                    "count_em": [0.0],
                    "tpm": [0.0],
                    "effective_length": [0.0],
                }
            ),
        ],
        ignore_index=True,
    ).assign(tpm=tpm)
    row = qa.score_transcripts(quant, truth)
    assert row["tpm_abs_err"] == pytest.approx(0.0, abs=1e-6), row["tpm_abs_err"]
    assert row["no_length_fragments"] == 30.0
    # the length-free count rows still see it
    assert row["count_true"] == 180.0


def test_the_count_rows_still_carry_the_shadow_fragments():
    truth, quant = _annotated()
    shadow = pd.concat(
        [truth, pd.DataFrame({"transcript_id": ["shadow"], "observed_mrna_fragments": [50.0]})],
        ignore_index=True,
    )
    row = qa.score_transcripts(quant, shadow)
    assert row["count_true"] == 200.0 and row["fn_mass"] == 50.0
