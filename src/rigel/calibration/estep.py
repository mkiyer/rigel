"""E-step — per-region soft allocation of mass into gDNA vs RNA (doc 03 §3).

Runs identically on each of the substrate's three views (contained, left,
right). It combines a **count** log-Bayes-factor (NB gDNA-only vs gDNA+RNA
mixture) and a **strand** log-Bayes-factor (Beta-Binomial gDNA vs RNA) into a
per-region gDNA mixing proportion ``π_g``, then soft-allocates the unspliced
mass; the spliced mass is deterministic RNA (gDNA splice-artifacts are removed
upstream by the ``alignable`` splice blacklist, not modelled here; PR 7).

**D1 (power from flux, density from mass).** The log-Bayes-factors are driven by
the integer **flux** (``n_unspliced``, oriented ``k_sense``); that ``π_g``
allocates the fractional **mass**. For the contained view flux == mass; for the
boundary views the integer crossing flux drives ``π_g`` and the fractional
per-side mass is split.

**D7 (strand-ambiguous regions).** A region with transcripts on both strands has
no valid sense split, so its strand log-Bayes-factor is omitted (set to 0) and
``π_g`` comes from the count channel + prior. Intergenic (TS_NONE) regions keep
the strand term (gDNA is unstranded → an arbitrary sense is neutral).

The count channel is **silent where the mixture's RNA mean ``μ_d`` is 0** (the
two hypotheses coincide) — the case on the very first pass (PR 4a), before any
RNA mass is estimated. It activates once the PR 5 outer loop feeds ``M_d`` back.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.special import expit, logit
from scipy.stats import betabinom, nbinom

from .signature import TS_AMBIG, TS_NEG

if TYPE_CHECKING:
    from .substrate import SubstrateView

# Keep π_g strictly inside (0, 1) so allocated masses never land exactly on
# 0 or n. Spec: doc 03 §8 `_PI_CLIP` — a machine-scale guard, not a tunable.
_PI_CLIP = 1.0e-6


@dataclass(frozen=True, slots=True)
class Allocation:
    """Per-region soft-allocated masses from one E-step view (all float64[R])."""

    m_g: np.ndarray  # gDNA mass: the unspliced split alone (spliced is RNA, PR 7)
    m_d: np.ndarray  # RNA mass: unspliced split + all spliced
    m_g_unspl: np.ndarray  # unspliced gDNA mass alone (for the M-step, PR 5)
    m_d_unspl: np.ndarray  # unspliced RNA mass alone (count-channel μ_d + M-step, PR 5)
    pi_g: np.ndarray  # gDNA mixing proportion used
    k_sense: np.ndarray  # int64[R] — oriented sense flux (BB-strand M-step, PR 5)


def _sense_flux(view: "SubstrateView", ts_class: np.ndarray) -> np.ndarray:
    """Orient the unspliced genome-strand flux to transcript sense (doc 03 §3.2).

    TS_POS / TS_NONE → genome '+'; TS_NEG → genome '−'. NONE is intergenic, where
    gDNA is unstranded so the choice is neutral. The library's sense convention
    is already inside κ_rna (measured with the same align-matches-transcript
    convention), so no extra flip is applied. AMBIG regions skip the strand term
    (masked in :func:`_llr_strand`), so their value here is unused.
    """
    return np.where(ts_class == TS_NEG, view.n_unspliced_neg, view.n_unspliced_pos)


def _llr_count(
    n_u: np.ndarray, mu_g: np.ndarray, m_d_prev: np.ndarray, inv_dispersion: float
) -> np.ndarray:
    """NB log-Bayes-factor: gDNA-only mean μ_g vs mixture μ_g + μ_d (doc 03 §3.1).

    Zero where μ_d (``m_d_prev``) is 0 or μ_g is 0: the gDNA-only and mixture
    hypotheses coincide, so the count channel carries no information. On the
    PR 4a single pass ``m_d_prev`` is all-zero and this is identically 0.
    """
    live = (m_d_prev > 0.0) & (mu_g > 0.0)
    if not np.any(live):
        return np.zeros(n_u.shape, dtype=np.float64)
    mu_mix = mu_g + m_d_prev
    p_g = inv_dispersion / (inv_dispersion + mu_g)
    p_mix = inv_dispersion / (inv_dispersion + mu_mix)
    ll_g = nbinom.logpmf(n_u, inv_dispersion, p_g)
    ll_mix = nbinom.logpmf(n_u, inv_dispersion, p_mix)
    return np.where(live, ll_g - ll_mix, 0.0)


def _llr_strand(
    k_sense: np.ndarray,
    n_u: np.ndarray,
    ts_class: np.ndarray,
    *,
    kappa_rna: float,
    rho_r_bb: float,
    rho_d_bb: float,
) -> np.ndarray:
    """BB log-Bayes-factor: gDNA BB(0.5, ρ_d) vs RNA BB(κ_rna, ρ_r) (doc 03 §3.2).

    Evaluated only where ``n_u > 0`` and the region's transcript strand is
    defined (not AMBIG — D7); 0 elsewhere. The gDNA Beta-Binomial is symmetric
    about 0.5 (``κ_d`` fixed), so it always has support on ``[0, n]`` and the
    ratio never produces ``inf − inf``.
    """
    out = np.zeros(n_u.shape, dtype=np.float64)
    valid = (n_u > 0) & (ts_class != TS_AMBIG)
    if not np.any(valid):
        return out
    a_d = 0.5 * (1.0 - rho_d_bb) / rho_d_bb
    a_r = kappa_rna * (1.0 - rho_r_bb) / rho_r_bb
    b_r = (1.0 - kappa_rna) * (1.0 - rho_r_bb) / rho_r_bb
    kv = k_sense[valid]
    nv = n_u[valid]
    ll_g = betabinom.logpmf(kv, nv, a_d, a_d)
    ll_d = betabinom.logpmf(kv, nv, a_r, b_r)
    out[valid] = ll_g - ll_d
    return out


def estep_view(
    view: "SubstrateView",
    ts_class: np.ndarray,
    *,
    omega: np.ndarray,
    rho_0: float,
    L_eff: np.ndarray,
    exposure_dispersion: float,
    kappa_rna: float,
    rho_r_bb: float,
    rho_d_bb: float,
    pi_g_prior: np.ndarray,
    m_d_unspl_prev: np.ndarray,
) -> Allocation:
    """Soft-allocate one view's mass into gDNA / RNA (doc 03 §3.3–3.6)."""
    inv_dispersion = 1.0 / exposure_dispersion
    mu_g = omega * rho_0 * L_eff
    n_u = view.n_unspliced  # int64 flux — drives the log-Bayes-factors (D1)
    k_sense = _sense_flux(view, ts_class)

    llr_count = _llr_count(n_u.astype(np.float64), mu_g, m_d_unspl_prev, inv_dispersion)
    llr_strand = _llr_strand(
        k_sense, n_u, ts_class, kappa_rna=kappa_rna, rho_r_bb=rho_r_bb, rho_d_bb=rho_d_bb
    )
    pi_g = expit(logit(pi_g_prior) + llr_count + llr_strand)
    pi_g = np.clip(pi_g, _PI_CLIP, 1.0 - _PI_CLIP)

    # π_g (from flux) allocates the fractional unspliced mass (density from mass).
    m_g_unspl = pi_g * view.mass_unspliced
    m_d_unspl = view.mass_unspliced - m_g_unspl
    # Spliced fragments are deterministic RNA (PR 7): no gDNA splice-artifact
    # carve-out — artifacts are removed upstream by the alignable splice
    # blacklist. gDNA mass is the unspliced split alone.
    m_g = m_g_unspl
    m_d = m_d_unspl + view.mass_spliced
    return Allocation(
        m_g=m_g, m_d=m_d, m_g_unspl=m_g_unspl, m_d_unspl=m_d_unspl, pi_g=pi_g, k_sense=k_sense
    )


__all__ = ["Allocation", "estep_view"]
