"""The library RNA sense fraction κ, ``rna_sense_frac`` (layer 4).

κ is the sense rate of the genuine junctions over annotated spliced unique mappers: the ``StrandModel``'s
:attr:`~rigel.strand_model.StrandModel.p_r1_sense` reads the per-sj table with splice artifacts and reversed junctions
removed (:func:`rigel.junction_fit.genuine_sense_fraction`), and ``rna_sense_frac = (n_same + 1) / (n_obs + 2)`` here
is that fit's κ. It strand-cleans the count density and parameterises the per-region strand likelihood. With no
spliced read the ``Beta(1, 1)`` prior gives κ = ½ and ``fallback_used`` is set; calibration reads such a library as
unstranded, with a dead strand channel (``region_init.strand_discriminability``), never as an error.

κ's own posterior width (:meth:`StrandModel.posterior_variance`) is not a strand overdispersion, the spread of single
junctions about κ: the strand likelihood is binomial and has none.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..strand_model import StrandModels


@dataclass(frozen=True, slots=True)
class StrandBalance:
    """Fitted RNA strand mean: the posterior-mean sense fraction κ + its provenance."""

    rna_sense_frac: float  # RNA sense mean (posterior mean), in (0, 1) — used by the deconv
    n_observations: int  # spliced strand observations (fragments) backing the posterior
    fallback_used: bool  # True when there are no spliced observations (rna_sense_frac → 0.5)


def fit_strand_balance(strand_model: "StrandModels") -> StrandBalance:
    """RNA strand mean: κ's posterior is ``Beta(n_same + 1, n_opp + 1)`` from the ``StrandModel``'s spliced 2×2
    counts, and ``rna_sense_frac`` is its mean. With ``n_obs == 0`` that is ½ and ``fallback_used`` is set.
    """
    n_obs = float(strand_model.n_observations)
    p_sense = float(strand_model.p_r1_sense)  # MLE n_same / n_obs (0.5 when n_obs == 0)
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


__all__ = ["StrandBalance", "fit_strand_balance"]
