"""The library RNA sense fraction κ, ``rna_sense_frac``, and the strand protocol decision (layer 4).

κ is the sense rate of the genuine junctions over annotated spliced unique mappers: the ``StrandModel``'s
:attr:`~rigel.strand_model.StrandModel.p_r1_sense` reads the per-sj table with splice artifacts and reversed junctions
removed (:func:`rigel.junction_fit.genuine_sense_fraction`), and ``rna_sense_frac = (n_same + 1) / (n_obs + 2)`` here
is that fit's κ. It strand-cleans the count density and parameterises the per-region strand likelihood. With no
spliced read the ``Beta(1, 1)`` prior gives κ = ½ and ``fallback_used`` is set; calibration reads such a library as
unstranded, with a dead strand channel (:func:`strand_discriminability`), never as an error.

κ's own posterior width (:meth:`StrandModel.posterior_variance`) is not a strand overdispersion, the spread of single
junctions about κ: the strand likelihood is binomial and has none.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.special import betaln

if TYPE_CHECKING:
    from ..strand_model import StrandModels


@dataclass(frozen=True, slots=True)
class StrandBalance:
    """Fitted RNA strand mean: the posterior-mean sense fraction κ + its provenance."""

    rna_sense_frac: float  # RNA sense mean (posterior mean), in (0, 1) — used by the deconv
    n_observations: int  # spliced strand observations (fragments) backing the posterior
    fallback_used: bool  # True when there are no spliced observations (rna_sense_frac → 0.5)

    @property
    def live(self) -> bool:
        """Whether the strand channel is live: its :func:`strand_discriminability` is positive."""
        return strand_discriminability(self.rna_sense_frac, self.n_observations) > 0.0

    @property
    def channel_kappa(self) -> float:
        """κ as a reader of the strand channel takes it: the fitted κ on a live channel, ½ on a dead one."""
        return self.rna_sense_frac if self.live else 0.5


def fit_strand_balance(strand_model: "StrandModels") -> StrandBalance:
    """RNA strand mean: κ's posterior is ``Beta(n_same + 1, n_opp + 1)`` from the ``StrandModel``'s spliced 2×2
    counts, and ``rna_sense_frac`` is its mean. With ``n_obs == 0`` that is ½ and ``fallback_used`` is set.
    """
    n_obs = float(strand_model.n_observations)
    p_sense = float(
        strand_model.p_r1_sense
    )  # the genuine junctions' sense count over n_obs (0.5 when n_obs == 0)
    n_same = p_sense * n_obs
    n_opp = n_obs - n_same
    a = n_same + 1.0  # rna_sense_frac posterior Beta(a, b) — Laplace prior Beta(1, 1)
    b = n_opp + 1.0
    rna_sense_frac = a / (a + b)  # posterior mean = (n_same + 1) / (n_obs + 2)
    return StrandBalance(
        rna_sense_frac=rna_sense_frac,
        n_observations=int(round(n_obs)),
        fallback_used=(n_obs == 0.0),
    )


def strand_discriminability(kappa, n_rna_obs) -> float:
    """The library's strand discriminability: ``4(κ−½)²`` if the protocol preserves strand, else exactly 0.

    It is the library-level factor of the strand Fisher information (the kernel's ``strand_evidence``) and the
    message layer's test of whether the strand split is a witness (``messages.ChainView.strand_live``).

    The decision is a Bayes factor on the spliced sense counts (``n_same`` of ``N``, from :func:`fit_strand_balance`).
    H0, unstranded: read 1's strand is independent of the transcript's, so κ = ½ exactly and the split is
    ``Binomial(N, ½)``. H1: κ free under the fit's ``Beta(1, 1)`` prior. Then

        ln BF₁₀ = N·ln 2 + ln B(a, b),   a = n_same + 1 = κ̂·(N + 2),   b = n_opp + 1 = (1 − κ̂)·(N + 2),

    and the channel is live iff ``BF₁₀ > 1`` (equal prior odds, so no threshold constant). For large ``N``,
    ``ln BF₁₀ ≈ ½·[z² − ln(2N/π)]`` with ``z = (κ̂ − ½)/√(¼/N)``: the free parameter's Occam penalty grows with the
    sample, so a sampling fluctuation of an unstranded library is not read as a protocol. gDNA must stay out of this
    test: its strand mean is ½ by symmetry, and a gDNA term would switch off a gDNA-free library's channel. With no
    spliced observation the two hypotheses tie and the channel is dead.

    Gate: ``tests/calibration/test_strand_balance.py``.
    """
    n = float(n_rna_obs)
    if not n > 0.0:
        return 0.0
    a = kappa * (n + 2.0)
    b = (1.0 - kappa) * (n + 2.0)
    ln_bf = n * np.log(2.0) + float(betaln(a, b))
    return 4.0 * (kappa - 0.5) ** 2 if ln_bf > 0.0 else 0.0


__all__ = ["StrandBalance", "fit_strand_balance", "strand_discriminability"]
