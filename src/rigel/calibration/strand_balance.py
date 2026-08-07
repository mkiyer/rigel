"""RNA strand-balance model: the library RNA sense mean ``rna_sense_frac``.

Fit from the **posterior-predictive** of the library sense rate over annotated spliced unique mappers (the
2×2 contingency in the live ``StrandModel``, itself the marginal of the per-junction SJ strand table). The
sense-rate posterior is ``Beta(n_same + 1, n_opp + 1)`` and ``rna_sense_frac`` is its **mean**,
``(n_same + 1) / (n_obs + 2)``: it strand-cleans the count density and parameterises the per-node strand
likelihood. Zero spliced reads ⇒ ``Beta(1,1)`` ⇒ 0.5, and ``calibrate`` raises ``CalibrationStrandError``
(a real RNA-seq library always has spliced reads).

⚠ **This module used to also return that posterior's own width, ``1/(n_obs+3)``, under the name
``rna_strand_overdispersion``.** It had no consumer anywhere in ``src/``, and the name collided with the
deconv's genuine RNA strand Beta-Binomial overdispersion — a different quantity on a different axis (the
spread of junctions ABOUT the mean, fitted in :func:`gdna_strand.fit_rna_strand_from_sj_table` and applied
in ψ's ``simplex_logodds._mixture_strand_loglik``, versus the width of the mean itself). **Deleted 2026-07-28**, which removes the
collision without needing a second name. Restoring it is one line if a consumer ever appears; the strand
model's own :meth:`StrandModel.posterior_variance` already reports the same power question.
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
    """RNA strand mean = posterior-predictive of the spliced sense rate (PR 9).

    The RNA strand Beta-Binomial is the posterior-predictive ``BB(n, n_same + 1, n_opp + 1)`` from the
    ``StrandModel``'s 2×2 spliced-fragment counts, and ``rna_sense_frac`` is its posterior mean. With
    ``n_obs == 0`` the ``Beta(1, 1)`` prior gives ``rna_sense_frac = 0.5`` (the strand channel adds no pull).
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
