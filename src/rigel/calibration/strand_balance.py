"""RNA strand-balance model: rna_sense_frac (used by the deconv) + rna_strand_overdispersion (QC).

Fit from the **posterior-predictive** of the library sense rate over annotated spliced unique
mappers — an uncontaminated, per-fragment measure of strand specificity (the 2x2 contingency
carried by the live ``StrandModel``). The sense-rate posterior is ``Beta(n_same + 1, n_opp + 1)``:

* ``rna_sense_frac`` = the posterior **mean** ``(n_same + 1) / (n_obs + 2)``. This is the live
  output: it strand-cleans the count density and parameterises the per-node strand likelihood
  (``strand_likelihood.strand_loglik``). With **zero** spliced reads ``Beta(1, 1)`` gives
  ``rna_sense_frac = 0.5`` (symmetric) and ``calibrate`` raises ``CalibrationStrandError`` —
  a real RNA-seq library always carries spliced reads.
* ``rna_strand_overdispersion`` = the posterior-predictive **overdispersion** ``1 / (n_obs + 3)``,
  a pure function of the spliced-read count, i.e. a **diagnostic of the strand model's
  statistical power**. It is **QC-only and NOT fed into the deconv** — it is a *thin-count
  statistical-power* quantity, distinct from the biological between-region RNA strand
  overdispersion the deconv actually uses.

**Do not confuse this QC quantity with the deconv's RNA strand overdispersion.** The deconv's RNA
Beta-Binomial overdispersion is a separate, *biological*, between-boundary quantity fitted from
boundary-side spliced counts in :mod:`gdna_strand` (``fit_rna_strand_from_substrate``) and applied
symmetrically with the gDNA overdispersion in :mod:`strand_likelihood` (see docs/em_strand/05). The
PR-9-era claim that "the deconv uses the Binomial strand limit" is **superseded**: that earlier
measurement (negligible / slightly worse silent-gene FP) concerned wiring *this 1/(n_obs+3)
posterior-predictive widening* into the deconv while gDNA was also Binomial — a different quantity
in a different regime. ``rna_sense_frac`` (the mean) is still the only field of this class the
deconv reads.

(PR 9 replaced the old method-of-moments fit / overdispersion floor / minimum-observation
fallback for *this* quantity, and by dropping the per-region substrate pool removed the boundary
double-count that biased the old MoM.)
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
