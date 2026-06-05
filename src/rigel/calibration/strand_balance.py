"""RNA strand-balance model (D2): rna_sense_frac (mean) + rna_strand_overdispersion (overdispersion).

The E-step's strand log-Bayes-factor (PR 4) compares an observed sense count
against an RNA Beta-Binomial ``BB(n, rna_sense_frac, rna_strand_overdispersion)`` and a gDNA
``BB(n, 0.5, gdna_strand_overdispersion)``. This module supplies the RNA parameters from the
**posterior-predictive** of the library sense rate (PR 9).

``rna_sense_frac`` and ``rna_strand_overdispersion`` are a single object: the rna_sense_frac **posterior**
``Beta(n_same + 1, n_opp + 1)`` over annotated spliced unique mappers (an
uncontaminated, per-fragment measure of library strand specificity — the 2×2
contingency carried by the live ``StrandModel``). The RNA strand BB is the
posterior-predictive ``BB(n, n_same + 1, n_opp + 1)``, so:

* ``rna_sense_frac`` = the posterior **mean** ``(n_same + 1) / (n_obs + 2)``;
* ``rna_strand_overdispersion`` = the posterior-predictive **overdispersion** ``1 / (n_obs + 3)`` —
  a pure function of the spliced-read count (statistical power).

The strand discriminator then degrades gracefully with the *evidence* rather than
via an arbitrary floor: with abundant spliced reads ``rna_strand_overdispersion → 0`` and the BB
collapses to the **Binomial** (the correct sharp limit); with few reads
``rna_strand_overdispersion`` grows so a single off-strand read no longer scores as gDNA; with **zero**
spliced reads ``Beta(1, 1)`` gives ``rna_sense_frac = 0.5`` (symmetric), so balanced gDNA
still routes to gDNA and the channel adds no spurious pull. No method-of-moments
fit, no overdispersion floor, no minimum-observation fallback (PR 9 replaced all
three — and, by dropping the per-region substrate pool, removed the boundary
double-count that biased the old MoM).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..strand_model import StrandModels


@dataclass(frozen=True, slots=True)
class StrandBalance:
    """Fitted RNA strand model: posterior-mean rna_sense_frac + Beta-Binomial overdispersion rna_strand_overdispersion."""

    rna_sense_frac: float  # RNA sense mean (posterior mean), in (0, 1)
    rna_strand_overdispersion: float  # RNA strand BB overdispersion = 1 / (n_obs + 3), in (0, 1)
    n_observations: int  # spliced strand observations (fragments) backing the posterior
    fallback_used: bool  # True when there are no spliced observations (rna_sense_frac → 0.5)
    fallback_reason: str


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
        fallback_reason="no spliced observations" if n_obs == 0.0 else "",
    )


__all__ = ["StrandBalance", "fit_strand_balance"]
