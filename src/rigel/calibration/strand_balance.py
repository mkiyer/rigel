"""RNA strand-balance model: the library RNA sense mean ``rna_sense_frac`` (κ).

Fitted from the posterior-predictive of the library sense rate over annotated spliced unique mappers —
the 2×2 contingency in the live ``StrandModel``, itself the marginal of the per-sj strand table. The
sense-rate posterior is ``Beta(n_same + 1, n_opp + 1)`` and ``rna_sense_frac`` is its mean,
``(n_same + 1) / (n_obs + 2)``. It strand-cleans the count density and parameterises the per-region
strand likelihood.

Zero spliced reads give ``Beta(1, 1)`` and hence κ = 0.5, which ``fallback_used`` announces;
calibration reads such a library as unstranded — its strand channel is dead
(``region_init.strand_discriminability``) — never as an error.

⛔ κ's own posterior width (:meth:`StrandModel.posterior_variance`) is not a strand overdispersion — the
spread of individual sj about κ, a different axis, which is 0 by policy.
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
    """RNA strand mean: the posterior mean of the spliced sense rate.

    The RNA strand Beta-Binomial is the posterior-predictive ``BB(n, n_same + 1, n_opp + 1)`` from the
    ``StrandModel``'s 2×2 spliced-fragment counts, and ``rna_sense_frac`` is its posterior mean. With
    ``n_obs == 0`` the ``Beta(1, 1)`` prior gives ``rna_sense_frac = 0.5`` (the strand channel adds no pull)
    and ``fallback_used`` is set.
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
