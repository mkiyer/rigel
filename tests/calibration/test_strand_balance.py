"""The strand protocol decision: whether calibration reads a library as stranded.

`strand_discriminability` decides by the Bayes factor of a free sense fraction κ against κ = ½ on the spliced strand
counts `fit_strand_balance` fits, never by a sampling band, so a balanced deep library reads unstranded and the channel
is dead at κ = ½ exactly. The pipeline's unstranded warning is that same decision, not a second test.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from rigel.calibration.strand_balance import fit_strand_balance, strand_discriminability
from rigel.pipeline import warn_if_calibration_strand_unidentifiable
from rigel.strand_model import StrandModel, StrandModels


def test_the_strand_channel_is_a_protocol_decision_not_a_sampling_band():
    """The strand channel is live by the Bayes factor of a free κ against κ = ½ exactly, not by a sampling band: an
    unbiased estimate of (κ−½)² floored at zero would read a balanced unstranded library live about a third of the time.
    A deep unstranded library fitting κ̂ ≈ 1σ from ½ reads dead, the stranded case live by a wide margin, and at κ = ½
    exactly the channel is dead at every N.
    """
    kappa, n = 0.500298, 3_090_137.0
    assert strand_discriminability(kappa, n) == 0.0
    assert 4.0 * ((kappa - 0.5) ** 2 - 0.25 / n) > 0.0, (
        "the replaced form's verdict on the same fit"
    )
    assert strand_discriminability(0.009883, 3_093_241.0) == pytest.approx(
        4.0 * (0.009883 - 0.5) ** 2
    )
    for n in (10.0, 1e3, 1e6):
        assert strand_discriminability(0.5, n) == 0.0


def test_the_occam_penalty_grows_with_the_spliced_sample():
    """ln BF₁₀ ≈ ½·[z² − ln(2N/π)]: the excursion a sampling fluctuation must clear to read as a protocol grows with the
    sample (3σ is live at N = 100 and dead at N = 10⁶, 4σ live at 10⁶ and dead at 10⁸), so no multiple of σ is the
    threshold; the free parameter's Occam penalty is."""

    def at(z, n):
        return strand_discriminability(0.5 + z * np.sqrt(0.25 / n), n)

    assert at(3.0, 100.0) > 0.0 and at(3.0, 1e6) == 0.0
    assert at(4.0, 1e6) > 0.0 and at(4.0, 1e8) == 0.0


def test_no_spliced_observations_is_no_protocol_decision():
    """With no spliced fragment the two hypotheses have equal marginal likelihood, so the channel is dead
    whatever κ a caller claims (`calibrate` raises before this on a real library)."""
    assert strand_discriminability(0.99, 0.0) == 0.0
    assert strand_discriminability(0.5, 0.0) == 0.0


def test_the_pipeline_warns_exactly_when_calibration_reads_the_library_unstranded(caplog):
    """The pipeline's warning is calibration's own strand decision (`strand_discriminability` on `fit_strand_balance`'s
    fit and its spliced count), not a second test: at N = 10⁵ a balanced spliced 2×2 warns, a sense fraction of 0.99 or
    0.52 is silent, and a 3σ excursion, which a 99 % normal-approximation band would read as stranded, warns because the
    Bayes factor reads it unstranded. No spliced observation is silent (`calibrate` raises there)."""

    def warns(n_same, n_opp):
        models = StrandModels(exonic_spliced=StrandModel(pos_pos=n_same, pos_neg=n_opp))
        caplog.clear()
        with caplog.at_level(logging.WARNING, logger="rigel.pipeline"):
            warn_if_calibration_strand_unidentifiable(models)
        return "reads the library as unstranded" in caplog.text

    n = 100_000
    assert warns(n // 2, n // 2)
    assert not warns(99_000, 1_000)
    assert not warns(52_000, 48_000)
    n_same = 50_474  # posterior-mean sense fraction 0.504740, 3.0σ above ½
    band = fit_strand_balance(
        StrandModels(exonic_spliced=StrandModel(pos_pos=n_same, pos_neg=n - n_same))
    )
    assert 2.9 < (band.rna_sense_frac - 0.5) / np.sqrt(0.25 / n) < 3.1
    assert strand_discriminability(band.rna_sense_frac, band.n_observations) == 0.0
    assert warns(n_same, n - n_same)
    assert not warns(0, 0)
