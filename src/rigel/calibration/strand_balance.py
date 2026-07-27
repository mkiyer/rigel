"""RNA strand-balance model: rna_sense_frac (used by the deconv) + rna_strand_overdispersion (QC-only).

Fit from the **posterior-predictive** of the library sense rate over annotated spliced unique mappers (the 2×2
contingency in the live ``StrandModel``). The sense-rate posterior is ``Beta(n_same + 1, n_opp + 1)``:

* ``rna_sense_frac`` = the posterior **mean** ``(n_same + 1) / (n_obs + 2)`` — the LIVE output: it strand-cleans
  the count density and parameterises the per-node strand likelihood. Zero spliced reads ⇒ ``Beta(1,1)`` ⇒ 0.5,
  and ``calibrate`` raises ``CalibrationStrandError`` (a real RNA-seq library always has spliced reads).
* ``rna_strand_overdispersion`` = the posterior-predictive overdispersion ``1 / (n_obs + 3)`` — a **QC-only**
  strand-model statistical-power diagnostic; NOT fed into the deconv. (Distinct from the deconv's *biological*
  between-boundary RNA overdispersion, fitted separately in :mod:`gdna_strand.fit_rna_strand_from_substrate`.)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..strand_model import StrandModels


@dataclass(frozen=True, slots=True)
class StrandBalance:
    """Fitted RNA strand model: posterior-mean rna_sense_frac (live) + rna_strand_overdispersion (QC)."""

    rna_sense_frac: float  # RNA sense mean (posterior mean), in (0, 1) — used by the deconv
    rna_strand_overdispersion: float  # 1/(n_obs+3); QC strand-power diagnostic, not in the deconv
    n_observations: int  # spliced strand observations (fragments) backing the posterior
    fallback_used: bool  # True when there are no spliced observations (rna_sense_frac → 0.5)


def fit_strand_balance(strand_model: "StrandModels") -> StrandBalance:
    """RNA strand model = posterior-predictive of the spliced sense rate (PR 9).

    The RNA strand Beta-Binomial is the posterior-predictive
    ``BB(n, n_same + 1, n_opp + 1)`` from the ``StrandModel``'s 2×2
    spliced-fragment counts; ``rna_sense_frac`` is its posterior mean and
    ``rna_strand_overdispersion = 1 / (n_obs + 3)`` its overdispersion. With ``n_obs == 0`` the
    ``Beta(1, 1)`` prior gives ``rna_sense_frac = 0.5`` (the strand channel adds no pull).
    """
    n_obs = float(strand_model.n_observations)
    p_sense = float(strand_model.p_r1_sense)  # MLE n_same / n_obs (0.5 when n_obs == 0)
    n_same = p_sense * n_obs
    n_opp = n_obs - n_same
    a = n_same + 1.0  # rna_sense_frac posterior Beta(a, b) — Laplace prior Beta(1, 1)
    b = n_opp + 1.0
    rna_sense_frac = a / (a + b)  # posterior mean = (n_same + 1) / (n_obs + 2)
    rna_strand_overdispersion = 1.0 / (
        a + b + 1.0
    )  # BB overdispersion of Beta(a, b) = 1 / (n_obs + 3)
    return StrandBalance(
        rna_sense_frac=rna_sense_frac,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_observations=int(round(n_obs)),
        fallback_used=(n_obs == 0.0),
    )


__all__ = ["StrandBalance", "fit_strand_balance"]
